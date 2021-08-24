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
    mReadsKOWriter = NULL;
    
    mDuplicate = NULL;
    if(mOptions->duplicate.enabled) {
        mDuplicate = new Duplicate(mOptions);
    }
    this->tbwtfmiDB = tbwtfmiDB;
    fileoutname.clear();
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

    if (mOptions->outputReadsKOMap && !mOptions->outReadsKOMap.empty()) {
        mReadsKOWriter = new WriterThread(mOptions, mOptions->outReadsKOMap);
    }
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

    if (mReadsKOWriter) {
        delete mReadsKOWriter;
        mReadsKOWriter = NULL;
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
    for(int t=0; t<mOptions->thread; t++){
        configs[t] = new ThreadConfig(mOptions, tbwtfmiDB, t, false);
        initConfig(configs[t]);
    }
    std::thread** threads = new thread*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        threads[t] = new std::thread(std::bind(&SingleEndProcessor::consumerTask, this, configs[t]));
    }

    std::thread* leftWriterThread = NULL;
    std::thread* failedWriterThread = NULL;
    std::thread* readsKOMapWriterThread = NULL;
    if(mLeftWriter)
        leftWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mLeftWriter));
    if(mFailedWriter)
        failedWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mFailedWriter));
    if (mReadsKOWriter)
        readsKOMapWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mReadsKOWriter));

    producer.join();
    for(int t=0; t<mOptions->thread; t++){
        threads[t]->join();
    }

    if(!mOptions->split.enabled) {
        if(leftWriterThread)
            leftWriterThread->join();
        if(failedWriterThread)
            failedWriterThread->join();
        if (readsKOMapWriterThread)
            readsKOMapWriterThread->join();
    }

    if(mOptions->verbose)
        loginfo("start to generate reports\n");

    // merge stats and read filter results
    vector<Stats*> preStats;
    vector<Stats*> postStats;
    vector<FilterResult*> filterResults;
    vector< std::unordered_map<std::string, uint32 > > totalKoFreqVecResults;
    totalKoFreqVecResults.reserve(mOptions->thread);
    vector< std::unordered_map<std::string, std::unordered_map<std::string, double> > > totalOrgKOFreqVecResults;
    totalOrgKOFreqVecResults.reserve(mOptions->thread);
    for(int t=0; t<mOptions->thread; t++){
        preStats.push_back(configs[t]->getPreStats1());
        postStats.push_back(configs[t]->getPostStats1());
        filterResults.push_back(configs[t]->getFilterResult());
        totalKoFreqVecResults.push_back(configs[t]->getTransSearcher()->getSubKoFreqUMap());
        totalOrgKOFreqVecResults.push_back(configs[t]->getTransSearcher()->getSubOrgKOAbunUMap());
    }
    Stats* finalPreStats = Stats::merge(preStats);
    Stats* finalPostStats = Stats::merge(postStats);
    FilterResult* finalFilterResult = FilterResult::merge(filterResults);
    
    mOptions->mHomoSearchOptions.nTotalReads = finalPreStats->getReads(); //change to both reads??????
    mOptions->mHomoSearchOptions.nCleanReads = finalPostStats->getReads();
    prepareResults(totalKoFreqVecResults, totalOrgKOFreqVecResults);

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
    
    
    // make JSON report
    JsonReporter jr(mOptions);
    jr.setDupHist(dupHist, dupMeanGC, dupRate);
    jr.report(finalFilterResult, finalPreStats, finalPostStats);
    
    // make HTML report
    HtmlReporter hr(mOptions);
    hr.setDupHist(dupHist, dupMeanGC, dupRate);
    hr.report(finalFilterResult, finalPreStats, finalPostStats);

    // clean up
    for(int t=0; t<mOptions->thread; t++){
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
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

    if(leftWriterThread)
        delete leftWriterThread;
    if(failedWriterThread)
        delete failedWriterThread;
    if (readsKOMapWriterThread)
        delete readsKOMapWriterThread;

    if(!mOptions->split.enabled)
        closeOutput();

    return true;
}

bool SingleEndProcessor::processSingleEnd(ReadPack* pack, ThreadConfig* config){
    string outstr;
    string failedOut;
    int readPassed = 0;
    std::string outReadsKOMapStr = "";
    std::string koTag = "";
    int mappedReads = 0;
    
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
            koTag = config->getTransSearcher()->transSearch(r1);

            if (koTag.length() > 0) {
                mappedReads++;
                if (mLeftWriter) {
                    outstr += r1->toStringWithTag(koTag);
                }
                if (mReadsKOWriter) {
                    outReadsKOMapStr += trimName(r1->mName) + "\t" + koTag + "\n";
                }
                koTag.clear();
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
    

    if (mOptions->verbose) {
        logMtx.lock();
        mOptions->transSearch.nTransMappedKOReads += mappedReads;
        auto rCount = mOptions->transSearch.nTransMappedKOReads;
        auto kCount = mOptions->transSearch.koUSet.size();
        logMtx.unlock();
        std::string str = "Mapped " + std::to_string(rCount) + " reads to " + std::to_string(kCount) + " KOs";
        loginfo(str);
    }

    mappedReads = 0;
        
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
    
    if (mReadsKOWriter && !outReadsKOMapStr.empty()) {
        char* tdata = new char[outReadsKOMapStr.size()];
        memcpy(tdata, outReadsKOMapStr.c_str(), outReadsKOMapStr.size());
        mReadsKOWriter->input(tdata, outReadsKOMapStr.size());
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

void SingleEndProcessor::consumePack(ThreadConfig* config){
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

    processSingleEnd(data, config);

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

void SingleEndProcessor::consumerTask(ThreadConfig* config)
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
            consumePack(config);
            //lock.unlock();
        } else {
            //lock.unlock();
            consumePack(config);
        }
    }

    if(mFinishedThreads == mOptions->thread) {
        if(mLeftWriter)
            mLeftWriter->setInputCompleted();
        if(mFailedWriter)
            mFailedWriter->setInputCompleted();
        if (mReadsKOWriter)
            mReadsKOWriter->setInputCompleted();
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

void SingleEndProcessor::prepareResults(std::vector< std::unordered_map<std::string, uint32 > > & totalKoFreqVecResults,
        std::vector< std::unordered_map<std::string, std::unordered_map<std::string, double> > > & totalOrgKOFreqVecResults) {

    if (mOptions->mHomoSearchOptions.prefix.size() == 0) {
        error_exit("sample prefix is not specific, quit now!");
    }

    //1. merge KO freq map;
    std::unordered_map<std::string, uint32 > totalKoFreqUMapResults;
    mOptions->transSearch.nTransMappedKOReads = 0;
    for (auto & it : totalKoFreqVecResults) {
        for (auto & itr : it) {
            totalKoFreqUMapResults[itr.first] += itr.second;
            mOptions->transSearch.nTransMappedKOReads += itr.second;
        }
    }
    totalKoFreqVecResults.clear();
    totalKoFreqVecResults.shrink_to_fit();

    if (mOptions->samples.size() > 0) {
        int sampleId = mOptions->getWorkingSampleId(mOptions->mHomoSearchOptions.prefix); // get the working sample id;
        mOptions->samples.at(sampleId).totalKoFreqUMapResults = totalKoFreqUMapResults;
    }

    mOptions->transSearch.nTransMappedKOs = totalKoFreqUMapResults.size();

    auto tmpSortedKOFreqVec = sortUMapToVector(totalKoFreqUMapResults);

    std::tuple <std::string, uint32, std::string> tmpKOTuple;
    fileoutname.clear();
    fileoutname = mOptions->mHomoSearchOptions.prefix + "_KO_abundance.txt";
    std::ofstream * fout = new std::ofstream();
    fout->open(fileoutname.c_str(), std::ofstream::out);
    if (!fout->is_open()) error_exit("Can not open abundance file: " + fileoutname);
    if (mOptions->verbose) loginfo("Starting to write KO abundance table");
    *fout << "#KO_id" << "\t" << "Reads_count" << "\t" << "KO_name" << "\n";
    for (auto & it : tmpSortedKOFreqVec) {
        auto KOName = mOptions->mHomoSearchOptions.ko_fullname_map.find(it.first);
        if (KOName != mOptions->mHomoSearchOptions.ko_fullname_map.end()) {
            tmpKOTuple = make_tuple(it.first, it.second, KOName->second);
            mOptions->transSearch.sortedKOFreqTupleVector.push_back(tmpKOTuple);
            *fout << it.first << "\t" << it.second << "\t" << KOName->second << "\n";
            mOptions->transSearch.nTransMappedKOReads += it.second;
        }
    }

    fout->flush();
    fout->close();
    if (fout) delete fout;
    if (mOptions->verbose) loginfo("Finish to write KO abundance table");
    tmpSortedKOFreqVec.clear();
    tmpSortedKOFreqVec.shrink_to_fit();

    //2. rarefiction curve;
    if (mOptions->mHomoSearchOptions.profiling) {
        std::vector<std::string> reshuff_vec;
        reshuff_vec.reserve(mOptions->transSearch.nTransMappedKOReads);
        for (auto & it : totalKoFreqUMapResults) {
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
                    int step = 50;
                    int step_size = floor(total / step);
                    auto first = reshuff_vec.begin();
                    std::map<long, int> rarefaction_map_tmp;
                    rarefaction_map_tmp.insert(pair<long, int>(0, 0));
                    for (int i = 1; i < step; i++) {
                        auto last = reshuff_vec.begin() + step_size * i;
                                std::vector<std::string> rarefaction_vec(first, last);
                                std::sort(rarefaction_vec.begin(), rarefaction_vec.end());
                                int unic = std::unique(rarefaction_vec.begin(), rarefaction_vec.end()) - rarefaction_vec.begin();
                                rarefaction_map_tmp[(long) round((step_size * i) * ratio)] = unic;
                                rarefaction_vec.clear();
                    }
                    std::sort(reshuff_vec.begin(), reshuff_vec.end());
                    int unic = std::unique(reshuff_vec.begin(), reshuff_vec.end()) - reshuff_vec.begin();
                    rarefaction_map_tmp.insert(pair<long, int>(total_reads_html_report, unic));
                    reshuff_vec.clear();
                    return (rarefaction_map_tmp);
                    rarefaction_map_tmp.clear();
                }, reshuff_vec, mOptions->mHomoSearchOptions.nTotalReads);

        future_rarefaction.wait();
        mOptions->transSearch.rarefactionMap = future_rarefaction.get();
        reshuff_vec.clear();
        reshuff_vec.shrink_to_fit();
    }

    //3. pathway
    if (mOptions->mHomoSearchOptions.profiling) {
        fileoutname.clear();
        fileoutname = mOptions->mHomoSearchOptions.prefix + "_pathway_hits.txt";
        std::ofstream * fout = new std::ofstream();
        fout->open(fileoutname.c_str(), std::ofstream::out);
        if (!fout->is_open()) error_exit("Can not open pathway hits file: " + fileoutname);
        if (mOptions->verbose) loginfo("Starting to write pathway hits table");
        *fout << "Pathway_ID" << "\t" << "Pathway_Name" << "\t" << "KO_ID" << "\t" << "KO_count" << "\t" << "KO_Name" << "\n";

        std::vector<std::string> tmpPathwayVec; //for report;
        tmpPathwayVec.reserve(totalKoFreqUMapResults.size());
        std::string pathwayID;
        std::string pathwayName;
        for (auto & it : mOptions->mHomoSearchOptions.pathway_ko_multimap) {
            auto itr = totalKoFreqUMapResults.find(it.second); //get the KO abundance;
            if (itr != totalKoFreqUMapResults.end()) {
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

        totalKoFreqUMapResults.clear();

        //get the freq pathway map;
        std::unordered_map<std::string, int> tmpPathwayMap;
        for (auto & it : tmpPathwayVec) {
            tmpPathwayMap[it]++;
        }

        if (mOptions->samples.size() > 0) {
            int sampleId = mOptions->getWorkingSampleId(mOptions->mHomoSearchOptions.prefix); // get the working sample id;
            mOptions->samples.at(sampleId).totalPathwayMap = tmpPathwayMap;
        }
        
        mOptions->transSearch.nMappedPathways = tmpPathwayMap.size();

        tmpPathwayVec.clear();
        tmpPathwayVec.shrink_to_fit();

        std::vector<std::tuple<std::string, double, int, int> > tmpPathwayDoubleIntVec;
        for (auto & it : tmpPathwayMap) {
            auto itr = mOptions->mHomoSearchOptions.pathway_ko_stats_umap.find(it.first);
            if (itr != mOptions->mHomoSearchOptions.pathway_ko_stats_umap.end()) {
                double perc = getPercentage(it.second, itr->second);
                auto itt = std::make_tuple(it.first, perc, it.second, itr->second);
                tmpPathwayDoubleIntVec.push_back(itt);
            }
        }
        tmpPathwayMap.clear();

        //sort the freq pathway map;
        auto tmpPathwayFreqSortedVec = sortTupleVector(tmpPathwayDoubleIntVec);
        tmpPathwayDoubleIntVec.clear();
        tmpPathwayDoubleIntVec.shrink_to_fit();

        std::tuple<std::string, double, std::string, int, int> tmpPathwayTuple;

        for (auto & it : tmpPathwayFreqSortedVec) {
            std::string pathwayIDName = get<0>(it);
            std::string::size_type pos = pathwayIDName.find_first_of(":");
            pathwayID = pathwayIDName.substr(0, pos);
            pathwayName = pathwayIDName.substr(pos + 1);
            tmpPathwayTuple = make_tuple(pathwayID, get<1>(it), pathwayName, get<2>(it), get<3>(it));
            mOptions->transSearch.sortedPathwayFreqTupleVector.push_back(tmpPathwayTuple); //for report in html;
            pathwayID.clear();
            pathwayName.clear();
        }
        tmpPathwayFreqSortedVec.clear();
        tmpPathwayFreqSortedVec.shrink_to_fit();
    }

    totalKoFreqUMapResults.clear();

    //4.for species;
    if (mOptions->mHomoSearchOptions.profiling) {
        std::unordered_map<std::string, int> orgKOUMap;
        std::multimap<std::string, std::unordered_map<std::string, double> > tmpOrgKOFeqUMap;
        for (auto & it : totalOrgKOFreqVecResults) {
            tmpOrgKOFeqUMap.insert(it.begin(), it.end());
        }

        for (auto it = tmpOrgKOFeqUMap.begin(); it != tmpOrgKOFeqUMap.end(); it = tmpOrgKOFeqUMap.upper_bound(it->first)) {
            auto org = it->first;
            auto orgKO = tmpOrgKOFeqUMap.equal_range(org);

            std::multimap<std::string, double> tmpKOFreqMMap;
            for (auto & itt = orgKO.first; itt != orgKO.second; itt++) {
                tmpKOFreqMMap.insert(itt->second.begin(), itt->second.end());
            }

            std::unordered_map<std::string, double> tmpKOFreqMap;
            for (auto itk = tmpKOFreqMMap.begin(); itk != tmpKOFreqMMap.end(); itk = tmpKOFreqMMap.upper_bound(itk->first)) {
                auto ko = itk->first;
                auto koFreq = tmpKOFreqMMap.equal_range(ko);
                for (auto & itko = koFreq.first; itko != koFreq.second; itko++) {
                    tmpKOFreqMap[itko->first] += itko->second;
                }
            }
            tmpKOFreqMMap.clear();
            int nKOs = 0;
            for (auto & ko : tmpKOFreqMap) {
                if (ko.second > 1) {
                    nKOs++;
                }
            }
            tmpKOFreqMap.clear();
            orgKOUMap[org] = nKOs;
        }
        tmpOrgKOFeqUMap.clear();

        if (mOptions->samples.size() > 0) {
            int sampleId = mOptions->getWorkingSampleId(mOptions->mHomoSearchOptions.prefix); // get the working sample id;
            mOptions->samples.at(sampleId).totalOrgKOUMap = orgKOUMap;
        }
        
        auto sortedOrgKOFreqVec = sortUMapToVector(orgKOUMap);
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
        mOptions->transSearch.sortedOrgFreqVec = sortedOrgKOFreqVec;
        mOptions->transSearch.nMappedOrgs = sortedOrgKOFreqVec.size();
        sortedOrgKOFreqVec.clear();
        sortedOrgKOFreqVec.shrink_to_fit();
    }

    time_t t_finished = time(NULL);
    mOptions->transSearch.endTime = t_finished;
    mOptions->transSearch.timeLapse = difftime(mOptions->transSearch.endTime, mOptions->transSearch.startTime);
    if (mOptions->samples.size() > 0) {
        int sampleId = mOptions->getWorkingSampleId(mOptions->mHomoSearchOptions.prefix); // get the working sample id;
        mOptions->samples.at(sampleId).totalRawReads = mOptions->mHomoSearchOptions.nTotalReads;
        mOptions->samples.at(sampleId).totalCleanReads = mOptions->mHomoSearchOptions.nCleanReads;
        mOptions->samples.at(sampleId).cleanReadsRate = double(mOptions->samples.at(sampleId).totalCleanReads * 100) /  double(mOptions->samples.at(sampleId).totalRawReads);
        
        mOptions->samples.at(sampleId).nKODb = mOptions->transSearch.nKODB;
        mOptions->samples.at(sampleId).nKO = mOptions->transSearch.nTransMappedKOs;
        mOptions->samples.at(sampleId).koRate = double(mOptions->samples.at(sampleId).nKO * 100) / double(mOptions->samples.at(sampleId).nKODb);
        
        mOptions->samples.at(sampleId).transSearchMappedKOReads = mOptions->transSearch.nTransMappedKOReads;
        mOptions->samples.at(sampleId).mappedKOReadsRate = double(mOptions->samples.at(sampleId).transSearchMappedKOReads * 100) / double(mOptions->samples.at(sampleId).totalRawReads);
        
        mOptions->samples.at(sampleId).nPathwaysDb = mOptions->transSearch.nPathwaysDB;
        mOptions->samples.at(sampleId).nMappedPathways = mOptions->transSearch.nMappedPathways;
        mOptions->samples.at(sampleId).nOrgsDB = mOptions->transSearch.nOrgsDB;
        mOptions->samples.at(sampleId).nMappedOrgs = mOptions->transSearch.nMappedOrgs;
        mOptions->samples.at(sampleId).startTime = mOptions->transSearch.startTime;
        mOptions->samples.at(sampleId).endTime = mOptions->transSearch.endTime;
        mOptions->samples.at(sampleId).timeLapse = mOptions->transSearch.timeLapse;
        mOptions->samples.at(sampleId).rarefactionMap = mOptions->transSearch.rarefactionMap;
    }
}
