#include "seprocessor.h"

SingleEndProcessor::SingleEndProcessor(Options* opt, BwtFmiDB * tbwtfmiDB){
    mOptions = opt;
    mProduceFinished = false;
    mFinishedThreads = 0;
    mFilter = new Filter(opt);
    mOutStream = NULL;
    mZipFile = NULL;
    mUmiProcessor = new UmiProcessor(opt);
    mLeftWriter =  NULL;
    mFailedWriter = NULL;

    mDuplicate = NULL;
    if(mOptions->duplicate.enabled) {
        mDuplicate = new Duplicate(mOptions);
    }
    this->tbwtfmiDB = tbwtfmiDB;
    fileoutname.clear();
    sortedKOFreqTupleVector.clear();
    sortedPathwayFreqTupleVector.clear();
    rarefaction_map.clear();
    sortedOrgKOFreqVec.clear();
    preOrgKOMMap.clear();
}

SingleEndProcessor::~SingleEndProcessor() {
    delete mFilter;
    if(mDuplicate) {
        delete mDuplicate;
        mDuplicate = NULL;
    }
    if (mUmiProcessor) {
        delete mUmiProcessor;
        mUmiProcessor = NULL;
    }
}

void SingleEndProcessor::initOutput() {
    if(!mOptions->failedOut.empty())
        mFailedWriter = new WriterThread(mOptions, mOptions->failedOut);
    if(mOptions->out1.empty())
        return;
    mLeftWriter = new WriterThread(mOptions, mOptions->out1);
}

void SingleEndProcessor::closeOutput() {
    if(mLeftWriter) {
        delete mLeftWriter;
        mLeftWriter = NULL;
    }
    if(mFailedWriter) {
        delete mFailedWriter;
        mFailedWriter = NULL;
    }
}

void SingleEndProcessor::initConfig(ThreadConfig* config) {
    if(mOptions->out1.empty())
        return;

    if(mOptions->split.enabled) {
        config->initWriterForSplit();
    }
}

bool SingleEndProcessor::process(){
    if(!mOptions->split.enabled)
        initOutput();

    initPackRepository();
    std::thread producer(std::bind(&SingleEndProcessor::producerTask, this));

    //TODO: get the correct cycles
    int cycle = 151;
    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    TransSearcher** transSearchers = new TransSearcher*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        configs[t] = new ThreadConfig(mOptions, t, false);
        initConfig(configs[t]);
        transSearchers[t] = new TransSearcher(tbwtfmiDB, mOptions);
    }
    std::thread** threads = new thread*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        threads[t] = new std::thread(std::bind(&SingleEndProcessor::consumerTask, this, configs[t], transSearchers[t]));
    }

    std::thread* leftWriterThread = NULL;
    std::thread* failedWriterThread = NULL;
    if(mLeftWriter)
        leftWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mLeftWriter));
    if(mFailedWriter)
        failedWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mFailedWriter));

    producer.join();
    for(int t=0; t<mOptions->thread; t++){
        threads[t]->join();
    }

    if(!mOptions->split.enabled) {
        if(leftWriterThread)
            leftWriterThread->join();
        if(failedWriterThread)
            failedWriterThread->join();
    }

    if(mOptions->verbose)
        loginfo("start to generate reports\n");

    // merge stats and read filter results
    vector<Stats*> preStats;
    vector<Stats*> postStats;
    vector<FilterResult*> filterResults;
    for(int t=0; t<mOptions->thread; t++){
        preStats.push_back(configs[t]->getPreStats1());
        postStats.push_back(configs[t]->getPostStats1());
        filterResults.push_back(configs[t]->getFilterResult());
    }
    Stats* finalPreStats = Stats::merge(preStats);
    Stats* finalPostStats = Stats::merge(postStats);
    FilterResult* finalFilterResult = FilterResult::merge(filterResults);
    
    mOptions->mHomoSearchOptions.totalOrigReads = finalPreStats->getReads();

    // read filter results to the first thread's
    for(int t=1; t<mOptions->thread; t++){
        preStats.push_back(configs[t]->getPreStats1());
        postStats.push_back(configs[t]->getPostStats1());
    }
    int* dupHist = NULL;
    double* dupMeanTlen = NULL;
    double* dupMeanGC = NULL;
    double dupRate = 0.0;
    if(mOptions->duplicate.enabled) {
        dupHist = new int[mOptions->duplicate.histSize];
        memset(dupHist, 0, sizeof(int) * mOptions->duplicate.histSize);
        dupMeanGC = new double[mOptions->duplicate.histSize];
        memset(dupMeanGC, 0, sizeof(double) * mOptions->duplicate.histSize);
        dupRate = mDuplicate->statAll(dupHist, dupMeanGC, mOptions->duplicate.histSize);
        cerr << endl;
        cerr << "Duplication rate (may be overestimated since this is SE data): " << dupRate * 100.0 << "%" << endl;
    }
    

    S2FReport();
    S2FReportTuple mS2FReportTuple = std::make_tuple(sortedKOFreqTupleVector, rarefaction_map, sortedPathwayFreqTupleVector, sortedOrgKOFreqVec);
    sortedKOFreqTupleVector.clear();
    rarefaction_map.clear();
    sortedPathwayFreqTupleVector.clear();
    sortedOrgKOFreqVec.clear();
    // make HTML report
    HtmlReporter hr(mOptions);
    hr.setDupHist(dupHist, dupMeanGC, dupRate);
    hr.report(mS2FReportTuple, finalFilterResult, finalPreStats, finalPostStats);

    // clean up
    for(int t=0; t<mOptions->thread; t++){
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
        delete transSearchers[t];
        transSearchers[t] = NULL;
    }

    delete finalPreStats;
    delete finalPostStats;
    delete finalFilterResult;

    if(mOptions->duplicate.enabled) {
        delete[] dupHist;
        delete[] dupMeanGC;
    }

    delete[] threads;
    delete[] configs;
    delete[] transSearchers;

    if(leftWriterThread)
        delete leftWriterThread;
    if(failedWriterThread)
        delete failedWriterThread;

    if(!mOptions->split.enabled)
        closeOutput();

    return true;
}

bool SingleEndProcessor::processSingleEnd(ReadPack* pack, ThreadConfig* config, TransSearcher * transSearcher){
    string outstr;
    string failedOut;
    int readPassed = 0;
    std::string KOTag = "";
    std::multimap<std::string, std::pair<std::string, double> >  preOrgKOAbunMMap;
    
    for(int p=0;p<pack->count;p++){

        // original read1
        Read* or1 = pack->data[p];

        // stats the original read before trimming
        config->getPreStats1()->statRead(or1);

        // handling the duplication profiling
        if(mDuplicate)
            mDuplicate->statRead(or1);

        // filter by index
        if(mOptions->indexFilter.enabled && mFilter->filterByIndex(or1)) {
            delete or1;
            continue;
        }
        
        // umi processing
        if(mOptions->umi.enabled)
            mUmiProcessor->process(or1);

        int frontTrimmed = 0;
        // trim in head and tail, and apply quality cut in sliding window
        Read* r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1, frontTrimmed);

        if(r1 != NULL) {
            if(mOptions->polyGTrim.enabled)
                PolyX::trimPolyG(r1, config->getFilterResult(), mOptions->polyGTrim.minLen);
        }

        if(r1 != NULL && mOptions->adapter.enabled){
            bool trimmed = false;
            if(mOptions->adapter.hasSeqR1)
                trimmed = AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.sequence, false);
            bool incTrimmedCounter = !trimmed;
            if(mOptions->adapter.hasFasta) {
                AdapterTrimmer::trimByMultiSequences(r1, config->getFilterResult(), mOptions->adapter.seqsInFasta, false, incTrimmedCounter);
            }

            if (mOptions->adapter.polyA) {
                AdapterTrimmer::trimPolyA(r1, config->getFilterResult(), false, incTrimmedCounter);
            }
        }

        if(r1 != NULL) {
            if(mOptions->polyXTrim.enabled)
                PolyX::trimPolyX(r1, config->getFilterResult(), mOptions->polyXTrim.minLen);
        }

        if(r1 != NULL) {
            if( mOptions->trim.maxLen1 > 0 && mOptions->trim.maxLen1 < r1->length())
                r1->resize(mOptions->trim.maxLen1);
        }

        int result = mFilter->passFilter(r1);
        
        config->addFilterResult(result, 1);

        if (r1 != NULL && result == PASS_FILTER) {
            KOTag.clear();
            transSearcher->transSearch(r1, KOTag, preOrgKOAbunMMap);
            if (KOTag.length() > 0 && mLeftWriter) {
                outstr += r1->toStringWithTag(KOTag);
                KOTag.clear();
            }

            // stats the read after filtering
            config->getPostStats1()->statRead(r1);
            readPassed++;
        } else if (mFailedWriter) {
            failedOut += or1->toStringWithTag(FAILED_TYPES[result]);
        }
        

        delete or1;
        // if no trimming applied, r1 should be identical to or1
        if(r1 != or1 && r1 != NULL)
            delete r1;
    }

    if(mOptions->verbose){
        logMtx.lock();
        auto rCount = mOptions->transSearch.tmpReadKOPairVec.size();
        auto kCount = mOptions->transSearch.KOSet.size();
        logMtx.unlock();
        std::string str = "Mapped " + std::to_string(rCount) + " reads to " + std::to_string(kCount) + " KOs";
        loginfo(str);
    }
    
    //for species hit
    if (mOptions->mHomoSearchOptions.profiling && preOrgKOAbunMMap.size() > 0) {
        std::vector<std::string> tmpUniqOrgVec;
        for (auto it = preOrgKOAbunMMap.begin(), end = preOrgKOAbunMMap.end();
                it != end; it = preOrgKOAbunMMap.upper_bound(it->first)) {
            tmpUniqOrgVec.push_back(it->first);
        }//get unique species

        std::multimap<std::string, double> tmpKOAbunMMap;
        std::multimap<std::string, std::string> tmpOrgKOMMap;
        
        double hits_num;
        for (auto & it : tmpUniqOrgVec) { 
            auto itr = preOrgKOAbunMMap.equal_range(it);
            for(auto itt = itr.first; itt != itr.second; ++ itt){
                tmpKOAbunMMap.insert(std::make_pair(itt->second.first, itt->second.second));
            }
            
            for (auto itr = tmpKOAbunMMap.begin(), end = tmpKOAbunMMap.end();
                    itr != end; itr = tmpKOAbunMMap.upper_bound(itr->first)) {
                auto range = tmpKOAbunMMap.equal_range(itr->first);
                hits_num = std::accumulate(range.first, range.second, 0.000f, [](double a, std::pair<std::string, double>b) {
                    return a + b.second;
                });
                if (hits_num >= 1) {
                    tmpOrgKOMMap.insert(std::make_pair(it, itr->first));
                }
            }// get  kos -> abundance
            tmpKOAbunMMap.clear();
        }
        preOrgKOAbunMMap.clear();
        tmpUniqOrgVec.clear();
        
        mSpecMtx.lock();
        preOrgKOMMap.insert(tmpOrgKOMMap.begin(), tmpOrgKOMMap.end());
        mSpecMtx.unlock();
        tmpOrgKOMMap.clear();
    }
        
    // if splitting output, then no lock is need since different threads write different files
    if(!mOptions->split.enabled)
        mOutputMtx.lock();
    if(mOptions->outputToSTDOUT) {
        fwrite(outstr.c_str(), 1, outstr.length(), stdout);
    } else if(mOptions->split.enabled) {
        // split output by each worker thread
        if(!mOptions->out1.empty())
            config->getWriter1()->writeString(outstr);
    } 

    if(mLeftWriter) {
        char* ldata = new char[outstr.size()];
        memcpy(ldata, outstr.c_str(), outstr.size());
        mLeftWriter->input(ldata, outstr.size());
    }
    if(mFailedWriter && !failedOut.empty()) {
        // write failed data
        char* fdata = new char[failedOut.size()];
        memcpy(fdata, failedOut.c_str(), failedOut.size());
        mFailedWriter->input(fdata, failedOut.size());
    }
    if(!mOptions->split.enabled)
        mOutputMtx.unlock();

    if(mOptions->split.byFileLines)
        config->markProcessed(readPassed);
    else
        config->markProcessed(pack->count);

    delete pack->data;
    delete pack;

    return true;
}

void SingleEndProcessor::initPackRepository() {
    mRepo.packBuffer = new ReadPack*[PACK_NUM_LIMIT];
    memset(mRepo.packBuffer, 0, sizeof(ReadPack*)*PACK_NUM_LIMIT);
    mRepo.writePos = 0;
    mRepo.readPos = 0;
    //mRepo.readCounter = 0;
    
}

void SingleEndProcessor::destroyPackRepository() {
    delete mRepo.packBuffer;
    mRepo.packBuffer = NULL;
}

void SingleEndProcessor::producePack(ReadPack* pack){
    //std::unique_lock<std::mutex> lock(mRepo.mtx);
    /*while(((mRepo.writePos + 1) % PACK_NUM_LIMIT)
        == mRepo.readPos) {
        //mRepo.repoNotFull.wait(lock);
    }*/

    mRepo.packBuffer[mRepo.writePos] = pack;
    mRepo.writePos++;

    /*if (mRepo.writePos == PACK_NUM_LIMIT)
        mRepo.writePos = 0;*/

    //mRepo.repoNotEmpty.notify_all();
    //lock.unlock();
}

void SingleEndProcessor::consumePack(ThreadConfig* config, TransSearcher * transSearcher){
    ReadPack* data;
    //std::unique_lock<std::mutex> lock(mRepo.mtx);
    // buffer is empty, just wait here.
    /*while(mRepo.writePos % PACK_NUM_LIMIT == mRepo.readPos % PACK_NUM_LIMIT) {
        if(mProduceFinished){
            //lock.unlock();
            return;
        }
        //mRepo.repoNotEmpty.wait(lock);
    }*/

    mInputMtx.lock();
    while(mRepo.writePos <= mRepo.readPos) {
        usleep(1000);
        if(mProduceFinished) {
            mInputMtx.unlock();
            return;
        }
    }
    data = mRepo.packBuffer[mRepo.readPos];
    mRepo.readPos++;

    /*if (mRepo.readPos >= PACK_NUM_LIMIT)
        mRepo.readPos = 0;*/
    mInputMtx.unlock();

    //lock.unlock();
    //mRepo.repoNotFull.notify_all();

    processSingleEnd(data, config, transSearcher);

}

void SingleEndProcessor::producerTask()
{
    if(mOptions->verbose)
        loginfo("start to load data");
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    bool splitSizeReEvaluated = false;
    Read** data = new Read*[PACK_SIZE];
    memset(data, 0, sizeof(Read*)*PACK_SIZE);
    FastqReader reader(mOptions->in1, true, mOptions->phred64);
    int count=0;
    bool needToBreak = false;
    while(true){
        Read* read = reader.read();
        // TODO: put needToBreak here is just a WAR for resolve some unidentified dead lock issue 
        if(!read || needToBreak){
            // the last pack
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            data = NULL;
            if(read) {
                delete read;
                read = NULL;
            }
            break;
        }
        data[count] = read;
        count++;
        // configured to process only first N reads
        if(mOptions->readsToProcess >0 && count + readNum >= mOptions->readsToProcess) {
            needToBreak = true;
        }
        if(mOptions->verbose && count + readNum >= lastReported + 1000000) {
            lastReported = count + readNum;
            string msg = "loaded " + to_string((lastReported/1000000)) + "M reads";
            loginfo(msg);
        }
        // a full pack
        if(count == PACK_SIZE || needToBreak){
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            //re-initialize data for next pack
            data = new Read*[PACK_SIZE];
            memset(data, 0, sizeof(Read*)*PACK_SIZE);
            // if the consumer is far behind this producer, sleep and wait to limit memory usage
            while(mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT){
                //cerr<<"sleep"<<endl;
                slept++;
                usleep(1000);
            }
            readNum += count;
            // if the writer threads are far behind this producer, sleep and wait
            // check this only when necessary
            if(readNum % (PACK_SIZE * PACK_IN_MEM_LIMIT) == 0 && mLeftWriter) {
                while(mLeftWriter->bufferLength() > PACK_IN_MEM_LIMIT) {
                    slept++;
                    usleep(1000);
                }
            }
            // reset count to 0
            count = 0;
            // re-evaluate split size
            // TODO: following codes are commented since it may cause threading related conflicts in some systems
            /*if(mOptions->split.needEvaluation && !splitSizeReEvaluated && readNum >= mOptions->split.size) {
                splitSizeReEvaluated = true;
                // greater than the initial evaluation
                if(readNum >= 1024*1024) {
                    size_t bytesRead;
                    size_t bytesTotal;
                    reader.getBytes(bytesRead, bytesTotal);
                    mOptions->split.size *=  (double)bytesTotal / ((double)bytesRead * (double) mOptions->split.number);
                    if(mOptions->split.size <= 0)
                        mOptions->split.size = 1;
                }
            }*/
        }
    }

    //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
    mProduceFinished = true;
    if(mOptions->verbose)
        loginfo("all reads loaded, start to monitor thread status");
    //lock.unlock();

    // if the last data initialized is not used, free it
    if(data != NULL)
        delete[] data;
}

void SingleEndProcessor::consumerTask(ThreadConfig* config, TransSearcher * transSearcher)
{
    while(true) {
        if(config->canBeStopped()){
            mFinishedThreads++;
            break;
        }
        while(mRepo.writePos <= mRepo.readPos) {
            if(mProduceFinished)
                break;
            usleep(1000);
        }
        //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
        if(mProduceFinished && mRepo.writePos == mRepo.readPos){
            mFinishedThreads++;
            if(mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " data processing completed";
                loginfo(msg);
            }
            //lock.unlock();
            break;
        }
        if(mProduceFinished){
            if(mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " is processing the " + to_string(mRepo.readPos) + " / " + to_string(mRepo.writePos) + " pack";
                loginfo(msg);
            }
            consumePack(config, transSearcher);
            //lock.unlock();
        } else {
            //lock.unlock();
            consumePack(config, transSearcher);
        }
    }

    if(mFinishedThreads == mOptions->thread) {
        if(mLeftWriter)
            mLeftWriter->setInputCompleted();
        if(mFailedWriter)
            mFailedWriter->setInputCompleted();
    }

    if(mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " finished";
        loginfo(msg);
    }
}

void SingleEndProcessor::writeTask(WriterThread* config)
{
    while(true) {
        if(config->isCompleted()){
            // last check for possible threading related issue
            config->output();
            break;
        }
        config->output();
    }

    if(mOptions->verbose) {
        string msg = config->getFilename() + " writer finished";
        loginfo(msg);
    }
}

void SingleEndProcessor::S2FReport(){
    //for reads KO map;
    if(mOptions->mHomoSearchOptions.profiling && mOptions->mHomoSearchOptions.prefix.size() != 0){
        fileoutname.clear();
        fileoutname = mOptions->mHomoSearchOptions.prefix + "_read_KO_map.txt";
        std::ofstream * fout = new std::ofstream();
        fout->open(fileoutname.c_str(), std::ofstream::out);
        if (!fout->is_open()) error_exit("Can not open abundance file: " + fileoutname);
        if (mOptions->verbose) loginfo("Starting to write read KO mapping table");
        *fout << "Read_id" << "\t" << "KO" << "\n";
        for(auto & it : mOptions->transSearch.tmpReadKOPairVec){
            //std::cout << it.first << " " << it.second << "\n";
            *fout << trimName(it.first) << "\t" << it.second << "\n";
        }
        fout->flush();
        fout->close();
        if(fout) delete fout;
        if (mOptions->verbose) loginfo("Finish to write read KO mapping table");
    }
    
    //for KO abundance file;
    mOptions->transSearch.transSearchMappedReads = 0;
    mOptions->transSearch.transSearchMappedReads = mOptions->transSearch.tmpReadKOPairVec.size();
    std::unordered_map<std::string, int> tmpKOFreqUMap;
    for(auto & it : mOptions->transSearch.tmpReadKOPairVec){
        tmpKOFreqUMap[it.second]++;
    }
    
    if(mOptions->samples.size() > 0){
        mOptions->transSearch.sampleKOAbunUMap.insert(tmpKOFreqUMap.begin(), tmpKOFreqUMap.end());
    }
    
    auto tmpSortedKOFreqVec = sortUMapToVector(tmpKOFreqUMap);
    if (mOptions->mHomoSearchOptions.prefix.size() != 0) {
        fileoutname.clear();
        fileoutname = mOptions->mHomoSearchOptions.prefix + "_KO_abundance.txt";
        
        sortedKOFreqTupleVector.clear();
        std::tuple <std::string, int, std::string> tmpKOTuple;
        
        std::ofstream * fout = new std::ofstream();
        fout->open(fileoutname.c_str(), std::ofstream::out);
        if (!fout->is_open()) error_exit("Can not open abundance file: " + fileoutname);
        if (mOptions->verbose) loginfo("Starting to write KO abundance table");
        *fout << "KO_id" << "\t" << "Reads_count" << "\t" << "KO_name" << "\n";
        for (auto & it : tmpSortedKOFreqVec) {
            auto KOName = mOptions->mHomoSearchOptions.ko_fullname_map.find(it.first);
            if(KOName != mOptions->mHomoSearchOptions.ko_fullname_map.end()){
                tmpKOTuple = make_tuple(it.first, it.second, KOName->second);
                sortedKOFreqTupleVector.push_back(tmpKOTuple);
                *fout << it.first << "\t" << it.second << "\t" << KOName->second << "\n";
            }
        }
        fout->flush();
        fout->close();
        if (fout) delete fout;
        if (mOptions->verbose) loginfo("Finish to write KO abundance table");
    }

    tmpSortedKOFreqVec.clear();
    mOptions->transSearch.tmpReadKOPairVec.clear();
    
     //for pathway KO map;
    if (mOptions->mHomoSearchOptions.profiling && tmpKOFreqUMap.size() > 0) {
        fileoutname.clear();
        fileoutname = mOptions->mHomoSearchOptions.prefix + "_pathway_hits.txt";
        std::ofstream * fout = new std::ofstream();
        fout->open(fileoutname.c_str(), std::ofstream::out);
        if (!fout->is_open()) error_exit("Can not open pathway hits file: " + fileoutname);
        if (mOptions->verbose) loginfo("Starting to write pathway hits table");
        *fout << "Pathway_ID" << "\t" << "Pathway_Name" << "\t" << "KO_ID" << "\t" << "KO_count" << "\t" << "KO_Name" << "\n";
        std::vector<std::string> tmpPathwayVec; //for report;
        std::string pathwayID;
        std::string pathwayName;
        for (auto & it : mOptions->mHomoSearchOptions.pathway_ko_multimap) {
            auto itr = tmpKOFreqUMap.find(it.second); //get the KO abundance;
            if (itr != tmpKOFreqUMap.end()) {
                tmpPathwayVec.push_back(it.first); //for report;
                std::string::size_type pos = it.first.find_first_of(":");
                pathwayID = it.first.substr(0, pos);
                pathwayName = it.first.substr(pos + 1);
                auto itt = mOptions->mHomoSearchOptions.ko_fullname_map.find(itr->first); //get KO full name;
                if (itt != mOptions->mHomoSearchOptions.ko_fullname_map.end()) {//get full KO name;
                    *fout << pathwayID << "\t" << pathwayName << "\t" << itr->first << "\t" << itr->second << "\t" << itt->second << "\n";
                }
                pathwayID.clear();
                pathwayName.clear();
            }
        }
        fout->flush();
        fout->close();
        if (fout) delete fout;
        if (mOptions->verbose) loginfo("Finish to write pathway hits table");

        //get the freq pathway map;
        std::unordered_map<std::string, int> tmpPathwayMap;
        for (auto & it : tmpPathwayVec) {
            tmpPathwayMap[it]++;
        }
        tmpPathwayVec.clear();
        
        std::vector<std::tuple<std::string, double, int, int> > tmpPathwayDoubleIntVec;
        for(auto & it : tmpPathwayMap){
            auto itr = mOptions->mHomoSearchOptions.pathway_ko_stats_umap.find(it.first);
            if(itr != mOptions->mHomoSearchOptions.pathway_ko_stats_umap.end()){
                double perc = getPercetageInt(it.second, itr->second);
                auto itt = std::make_tuple(it.first, perc, it.second, itr->second);
                tmpPathwayDoubleIntVec.push_back(itt);
            }
        }
        tmpPathwayMap.clear();
        
        //sort the freq pathway map;
        auto tmpPathwayFreqSortedVec = sortTupleVector(tmpPathwayDoubleIntVec);
        tmpPathwayDoubleIntVec.clear();
        
        std::tuple<std::string, double, std::string, int, int> tmpPathwayTuple;
        sortedPathwayFreqTupleVector.clear();
        for (auto & it : tmpPathwayFreqSortedVec) {
            std::string pathwayIDName = get<0>(it);
            std::string::size_type pos = pathwayIDName.find_first_of(":");
            pathwayID = pathwayIDName.substr(0, pos);
            pathwayName = pathwayIDName.substr(pos + 1);
            tmpPathwayTuple = make_tuple(pathwayID, get<1>(it), pathwayName, get<2>(it), get<3>(it));
            sortedPathwayFreqTupleVector.push_back(tmpPathwayTuple);
            pathwayID.clear();
            pathwayName.clear();
        }
        tmpPathwayFreqSortedVec.clear();
    }

//rarefaction curve;
    if (mOptions->mHomoSearchOptions.profiling) {
        std::vector<std::string> reshuff_vec;
        for (auto & it : tmpKOFreqUMap) {
            for (int i = 0; i < it.second; i++) {
                reshuff_vec.push_back(it.first);
            }
        }
        auto future_rarefaction = std::async(std::launch::async,
                [](std::vector<std::string> reshuff_vec,
                long total_reads_html_report) {
                    std::random_shuffle(reshuff_vec.begin(), reshuff_vec.end());
                    int total = reshuff_vec.size();
                    double ratio = total_reads_html_report / total;
                    int step = 100;
                    int step_size = floor(total / step);
                    auto first = reshuff_vec.begin();
                    std::map<int, int> rarefaction_map_tmp;
                    rarefaction_map_tmp.insert(pair<int, int>(0, 0));
                    for (int i = 1; i < step; i++) {
                        auto last = reshuff_vec.begin() + step_size * i;
                                std::vector<std::string> rarefaction_vec(first, last);
                                std::sort(rarefaction_vec.begin(), rarefaction_vec.end());
                                int unic = std::unique(rarefaction_vec.begin(), rarefaction_vec.end()) - rarefaction_vec.begin();
                                rarefaction_map_tmp[(int) round((step_size * i) * ratio)] = unic;
                                rarefaction_vec.clear();
                    }
                    std::sort(reshuff_vec.begin(), reshuff_vec.end());
                    int unic = std::unique(reshuff_vec.begin(), reshuff_vec.end()) - reshuff_vec.begin();
                    rarefaction_map_tmp.insert(pair<int, int>(total_reads_html_report, unic));
                    reshuff_vec.clear();
                    return (rarefaction_map_tmp);
                    rarefaction_map_tmp.clear();
                }, reshuff_vec, mOptions->mHomoSearchOptions.totalOrigReads);
                
       future_rarefaction.wait();
       rarefaction_map = future_rarefaction.get();
       reshuff_vec.clear();
    }
    tmpKOFreqUMap.clear();

    //for species hit
    if (mOptions->mHomoSearchOptions.profiling && preOrgKOMMap.size() > 0) {
        std::vector<std::string> tmpUniqOrgVec;
        std::unordered_set<std::string> KOUSet;
        std::unordered_map<std::string, int> tmpSpecHItMap;
        for(auto it = preOrgKOMMap.begin(), end = preOrgKOMMap.end();
                it != end; it = preOrgKOMMap.upper_bound(it->first)){
            auto itr = preOrgKOMMap.equal_range(it->first);
            KOUSet.clear();
            for(auto itt = itr.first; itt != itr.second; ++ itt){
                KOUSet.insert(itt->second);
            }
            tmpSpecHItMap.insert(std::make_pair(it->first, KOUSet.size()));
        }
        preOrgKOMMap.clear();
        tmpUniqOrgVec.clear();
        KOUSet.clear();

        sortedOrgKOFreqVec = sortUMapToVector(tmpSpecHItMap);
        tmpSpecHItMap.clear();
        fileoutname.clear();
        fileoutname = mOptions->mHomoSearchOptions.prefix + "_species_hits.txt";
        std::ofstream * fout = new std::ofstream();
        fout->open(fileoutname.c_str(), std::ofstream::out);
        if (!fout->is_open()) error_exit("Can not open species hits file: " + fileoutname);
        if (mOptions->verbose) loginfo("Starting to write species hits table");
        *fout << "Species_Name" << "\t" << "Number_of_KOs" << "\n";
        for (auto & it : sortedOrgKOFreqVec) {
            *fout << it.first << "\t" << it.second << "\n";
        }
        fout->flush();
        fout->close();
        if (fout) delete fout;
        if (mOptions->verbose) loginfo("Finish to write species hits table");
    }
}
