#include "seprocessor.h"

SingleEndProcessor::SingleEndProcessor(Options* & opt, BwtFmiDB * tbwtfmiDB){
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
    if (mFilter) {
        delete mFilter; 
        mFilter = NULL;
    }
    
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

    if (mOptions->outputReadsAnnoMap && !mOptions->outReadsKOMap.empty()) {
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
    std::thread** threads = new thread*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        configs[t] = new ThreadConfig(mOptions, tbwtfmiDB, t, false);
        initConfig(configs[t]);
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

    destroyPackRepository();
    
    if(mOptions->verbose){
        mOptions->longlog ? loginfolong("start to generate reports\n") : loginfo("start to generate reports\n");
    }

    // merge stats and read filter results
    vector<Stats*> preStats;
    vector<Stats*> postStats;
    vector<FilterResult*> filterResults;
    vector<std::map<const uint32 *, uint32> > totalIdFreqVecResults;
    totalIdFreqVecResults.reserve(mOptions->thread);
    for(int t=0; t<mOptions->thread; t++){
        preStats.push_back(configs[t]->getPreStats1());
        postStats.push_back(configs[t]->getPostStats1());
        filterResults.push_back(configs[t]->getFilterResult());
        totalIdFreqVecResults.push_back(configs[t]->getTransSearcher()->getIdFreqSubMap());
    }
    Stats* finalPreStats = Stats::merge(preStats);
    Stats* finalPostStats = Stats::merge(postStats);
    FilterResult* finalFilterResult = FilterResult::merge(filterResults);
    
    mOptions->mHomoSearchOptions.nTotalReads = finalPreStats->getReads(); //change to both reads??????
    mOptions->mHomoSearchOptions.nCleanReads = finalPostStats->getReads();
    mOptions->transSearch.totalIdFreqUMapResults = TransSearcher::merge(totalIdFreqVecResults);
    mOptions->transSearch.nTransMappedIds = mOptions->transSearch.totalIdFreqUMapResults.size();
    totalIdFreqVecResults.clear();
    
    prepareResults();

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
    string * outstr = new string();
    string failedOut;
    int readPassed = 0;
    std::string * outReadsKOMapStr = new string();
    //std::string koTag = "";
    uint32* orthId = NULL;
    std::set<uint32 *> idSet;
    
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
        
                // fix MGI
        if(mOptions->fixMGI) {
            or1->fixMGI();
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
            orthId = NULL;
            
            config->getTransSearcher()->transSearch(r1, orthId);

            if (orthId != NULL) {
                if(mOptions->verbose){
                    idSet.insert(orthId);
                }
                mappedReads++;
                if (mLeftWriter) {
                    *outstr += r1->toStringWithTag(orthId);
                }
                if (mReadsKOWriter) {
                    *outReadsKOMapStr += trimName(r1->mName) + "\t" + "s2f_" + paddingOs(std::to_string(*orthId)) + "\n";
                }
                orthId = NULL;
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
        mOptions->transSearch.nTransMappedIdReads += mappedReads;
        logMtx.lock();
        auto rCount = long(mOptions->transSearch.nTransMappedIdReads);
        mOptions->transSearch.idUSet.insert(idSet.begin(), idSet.end());
        auto iCount = mOptions->transSearch.idUSet.size();
//        auto kCount = mOptions->transSearch.koUSet.size();
//        auto gCount = mOptions->transSearch.goUSet.size();
        logMtx.unlock(); 
        if (mOptions->longlog) {
            std::string str = "Mapped " + std::to_string(rCount) + " reads to " + std::to_string(iCount) + " s2f ids!";
            loginfolong(str);
        } else {
            std::string str = "Mapped \033[1;31m" + std::to_string(rCount) + "\033[0m reads to \033[1;32m" + std::to_string(iCount) + "\033[0m s2f ids!";
            loginfo(str, false);
        }
    }

    mappedReads = 0;
    idSet.clear();
    // if splitting output, then no lock is need since different threads write different files
    if(!mOptions->split.enabled)
        mOutputMtx.lock();
    if(mOptions->outputToSTDOUT) {
        fwrite(outstr->c_str(), 1, outstr->length(), stdout);
    } else if(mOptions->split.enabled) {
        // split output by each worker thread
        if(!mOptions->out1.empty())
            config->getWriter1()->writeString(*outstr);
    } 

    if(mLeftWriter) {
        char* ldata = new char[outstr->size()];
        memcpy(ldata, outstr->c_str(), outstr->size());
        mLeftWriter->input(ldata, outstr->size());
    }
    if(mFailedWriter && !failedOut.empty()) {
        // write failed data
        char* fdata = new char[failedOut.size()];
        memcpy(fdata, failedOut.c_str(), failedOut.size());
        mFailedWriter->input(fdata, failedOut.size());
    }
    
    if (mReadsKOWriter && !outReadsKOMapStr->empty()) {
        char* tdata = new char[outReadsKOMapStr->size()];
        memcpy(tdata, outReadsKOMapStr->c_str(), outReadsKOMapStr->size());
        mReadsKOWriter->input(tdata, outReadsKOMapStr->size());
    }
    
    if(!mOptions->split.enabled)
        mOutputMtx.unlock();

    if(mOptions->split.byFileLines)
        config->markProcessed(readPassed);
    else
        config->markProcessed(pack->count);

    if(outstr){
        delete outstr;
        outstr = NULL;
    }
    
    if(outReadsKOMapStr){
        delete outReadsKOMapStr;
        outReadsKOMapStr = NULL;
    }
    
    delete[] pack->data;
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
    if(mRepo.packBuffer) delete[] mRepo.packBuffer; mRepo.packBuffer = NULL;
}

void SingleEndProcessor::producePack(ReadPack* pack){
    mRepo.packBuffer[mRepo.writePos] = pack;
    mRepo.writePos++;
}

void SingleEndProcessor::consumePack(ThreadConfig* config){
    ReadPack* data;
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
    mInputMtx.unlock();
    processSingleEnd(data, config);
}

void SingleEndProcessor::producerTask(){
    if(mOptions->verbose){
        mOptions->longlog ? loginfolong("start to load data") : loginfo("start to load data");
    }
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    bool splitSizeReEvaluated = false;
    Read** data = new Read*[PACK_SIZE];
    memset(data, 0, sizeof(Read*)*PACK_SIZE);
    FastqReader reader(mOptions->in1, true, mOptions->phred64, mOptions->fastqBufferSize);
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
            string msg = "\nloaded " + to_string((lastReported/1000000)) + "M reads";
            mOptions->longlog ? loginfolong(msg) : loginfo(msg);
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
        }
    }

    //std::unique_lock<std::mutex> lock(mRepo.readCounterMtx);
    mProduceFinished = true;
    if(mOptions->verbose){
        mOptions->longlog ? loginfolong("all reads loaded, start to monitor thread status") : loginfo("all reads loaded, start to monitor thread status");
    }
    //lock.unlock();

    // if the last data initialized is not used, free it
    if (data != NULL) {
        for (int i = 0; i < PACK_SIZE; ++i) {
            if (data[i] != NULL) {
                delete data[i];
                data[i] = NULL;
            }
        }
        delete[] data;
        data = NULL;
    }
}

void SingleEndProcessor::consumerTask(ThreadConfig* config){
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
                //loginfo(msg, false);
            }
            break;
        }
        if(mProduceFinished){
            if(mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " is processing the " + to_string(mRepo.readPos) + " / " + to_string(mRepo.writePos) + " pack";
                //loginfo(msg, false);
            }
            consumePack(config);
        } else {
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
        string msg = "\nthread " + to_string(config->getThreadId() + 1) + " finished";
        mOptions->longlog ? loginfolong(msg) : loginfo(msg);
    }
}

void SingleEndProcessor::writeTask(WriterThread* config){
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
        mOptions->longlog ? loginfolong(msg) : loginfo(msg);
    }
}

void SingleEndProcessor::prepareResults() {

    if (mOptions->mHomoSearchOptions.prefix.size() == 0) {
        error_exit("sample prefix is not specific, quit now!");
    }

    //s2f_id abundance
    if(!mOptions->transSearch.totalIdFreqUMapResults.empty()){
        //sort the map;

        fileoutname.clear();
        fileoutname = mOptions->mHomoSearchOptions.prefix + "_s2fid_abundance.txt";
        std::ofstream* fout = new std::ofstream();
        fout->open(fileoutname.c_str(), std::ofstream::out);
        if(!fout->is_open()) error_exit("Can not open abundance file: " + fileoutname);
        if (mOptions->verbose) {
            mOptions->longlog ? loginfolong("Starting to write gene abundance table") : loginfo("Starting to write gene abundance table");
        }
        *fout << "#s2f_id\tReads_cout\tannotation\tcore_ortho_freq\n";
        if(mOptions->transSearch.nTransMappedIdReads != 0) mOptions->transSearch.nTransMappedIdReads = 0;
        for(const auto & it : mOptions->transSearch.totalIdFreqUMapResults){
            mOptions->transSearch.nTransMappedIdReads += it.second;
            auto itt = mOptions->mHomoSearchOptions.fullDbMap.find(it.first);
            if(itt != mOptions->mHomoSearchOptions.fullDbMap.end()){
                *fout << "s2f_" << std::setfill('0') << std::setw(10) << *(it.first) << "\t" <<  it.second << "\t" << itt->second.ko << "|" << itt->second.go << "|" << itt->second.symbol << "|" << itt->second.gene << "\t" << itt->second.coreOrthoPer << "\n";
                if(itt->second.coreOrthoPer >= 0.90) mOptions->transSearch.nMappedCoreOrthos++;
            } else {
                *fout << "s2f_U" << "\t" << *(it.first) << "\t" <<  it.second << "\tU|U|U|U\t0\n";
            }
        }
        
        fout->flush();
        fout->close();
        if(fout != NULL){
            delete fout;
            fout = NULL;
        }
        
        if (mOptions->verbose) {
            mOptions->longlog ? loginfolong("Finish to write s2f id abundance table") : loginfo("Finish to write s2f id abundance table");
        }
        
        //2 rarefaction curve;        
        if (mOptions->mHomoSearchOptions.profiling && mOptions->transSearch.nTransMappedIdReads > 0) {  
            std::vector<uint32> reshuff_vec;
            reshuff_vec.reserve(mOptions->transSearch.nTransMappedIdReads);
            for(const auto & it : mOptions->transSearch.totalIdFreqUMapResults){
                for(int i = 0; i < it.second; i++){
                    reshuff_vec.push_back(*(it.first));
                }
            }
            auto future_rarefaction = std::async(std::launch::async,
                    [](std::vector<uint32> reshuff_vec,
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
                                    std::vector<uint32> rarefaction_vec(first, last);
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
                    }, reshuff_vec, mOptions->mHomoSearchOptions.nTotalReads);

            future_rarefaction.wait();
            mOptions->transSearch.rarefactionIdMap = future_rarefaction.get();
            reshuff_vec.clear();
            reshuff_vec.shrink_to_fit();
        }
        
        if (mOptions->samples.size() > 0) {
            int sampleId = mOptions->getWorkingSampleId(mOptions->mHomoSearchOptions.prefix); // get the working sample id;
            mOptions->samples.at(sampleId).totalIdFreqUMapResults = mOptions->transSearch.totalIdFreqUMapResults;
        }
    }

    time_t t_finished = time(NULL);
    mOptions->transSearch.endTime = t_finished;
    mOptions->transSearch.timeLapse = difftime(mOptions->transSearch.endTime, mOptions->transSearch.startTime);
    if (mOptions->samples.size() > 0) {
        int sampleId = mOptions->getWorkingSampleId(mOptions->mHomoSearchOptions.prefix); // get the working sample id;
        mOptions->samples.at(sampleId).totalRawReads = mOptions->mHomoSearchOptions.nTotalReads;
        mOptions->samples.at(sampleId).totalCleanReads = mOptions->mHomoSearchOptions.nCleanReads;
        mOptions->samples.at(sampleId).cleanReadsRate = double(mOptions->samples.at(sampleId).totalCleanReads * 100) / double(mOptions->samples.at(sampleId).totalRawReads);
        mOptions->samples.at(sampleId).startTime = mOptions->transSearch.startTime;
        mOptions->samples.at(sampleId).endTime = mOptions->transSearch.endTime;
        mOptions->samples.at(sampleId).timeLapse = mOptions->transSearch.timeLapse;
        mOptions->samples.at(sampleId).transSearchMappedIdReads = mOptions->transSearch.nTransMappedIdReads;
        mOptions->samples.at(sampleId).mappedIdReadsRate = double(mOptions->samples.at(sampleId).transSearchMappedIdReads * 100) / double(mOptions->samples.at(sampleId).totalRawReads);
        mOptions->samples.at(sampleId).nId = mOptions->transSearch.nTransMappedIds;
        mOptions->samples.at(sampleId).nIdDb = mOptions->transSearch.nIdDB;
        mOptions->samples.at(sampleId).idRate = double(mOptions->samples.at(sampleId).nId * 100) / double(mOptions->samples.at(sampleId).nIdDb);
        mOptions->samples.at(sampleId).coreOrthosDb = mOptions->transSearch.coreOrthosDb;
        mOptions->samples.at(sampleId).nMappedCoreOrthos = mOptions->transSearch.nMappedCoreOrthos;
        mOptions->samples.at(sampleId).nMappedCoreOrthoRate = double(mOptions->samples.at(sampleId).nMappedCoreOrthos * 100) / double(mOptions->transSearch.coreOrthosDb);
        mOptions->samples.at(sampleId).rarefactionIdMap = mOptions->transSearch.rarefactionIdMap;
    }
        
}
