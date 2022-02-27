#include "peprocessor.h"

PairEndProcessor::PairEndProcessor(Options* & opt, BwtFmiDB * & tbwtfmiDB) {
    mOptions = opt;
    this->tbwtfmiDB = tbwtfmiDB;
    mProduceFinished = false;
    mFinishedThreads = 0;
    mFilter = new Filter(opt);
    mOutStream1 = NULL;
    mZipFile1 = NULL;
    mOutStream2 = NULL;
    mZipFile2 = NULL;
    mUmiProcessor = new UmiProcessor(opt);

    int isizeBufLen = mOptions->insertSizeMax + 1;
    mInsertSizeHist = new atomic_long[isizeBufLen];
    memset(mInsertSizeHist, 0, sizeof (atomic_long)*isizeBufLen);
    mLeftWriter = NULL;
    mRightWriter = NULL;
    mUnpairedLeftWriter = NULL;
    mUnpairedRightWriter = NULL;
    mMergedWriter = NULL;
    mFailedWriter = NULL;
    mReadsKOWriter = NULL;

    mDuplicate = NULL;
    if (mOptions->duplicate.enabled) {
        mDuplicate = new Duplicate(mOptions);
    }
    fileoutname.clear();
}

PairEndProcessor::~PairEndProcessor() {
    if (mFilter){
        delete mFilter;
        mFilter = NULL;
    }
    
    if (mInsertSizeHist) {
        delete mInsertSizeHist; 
        mInsertSizeHist = NULL;
    }
    
    if (mDuplicate) {
        delete mDuplicate;
        mDuplicate = NULL;
    }
    
    if (mUmiProcessor) {
        delete mUmiProcessor;
        mUmiProcessor = NULL;
    }
    
    destroyPackRepository();
}

void PairEndProcessor::initOutput() {
    if (!mOptions->unpaired1.empty())
        mUnpairedLeftWriter = new WriterThread(mOptions, mOptions->unpaired1);

    if (!mOptions->unpaired2.empty() && mOptions->unpaired2 != mOptions->unpaired1)
        mUnpairedRightWriter = new WriterThread(mOptions, mOptions->unpaired2);

    if (mOptions->merge.enabled) {
        if (!mOptions->merge.out.empty())
            mMergedWriter = new WriterThread(mOptions, mOptions->merge.out);
    }

    if (!mOptions->failedOut.empty())
        mFailedWriter = new WriterThread(mOptions, mOptions->failedOut);

    if (mOptions->out1.empty())
        return;

    mLeftWriter = new WriterThread(mOptions, mOptions->out1);
    if (!mOptions->out2.empty())
        mRightWriter = new WriterThread(mOptions, mOptions->out2);

    if (mOptions->outputReadsAnnoMap && !mOptions->outReadsKOMap.empty()) {
        mReadsKOWriter = new WriterThread(mOptions, mOptions->outReadsKOMap);
    }
}

void PairEndProcessor::closeOutput() {
    if (mLeftWriter) {
        delete mLeftWriter;
        mLeftWriter = NULL;
    }
    if (mRightWriter) {
        delete mRightWriter;
        mRightWriter = NULL;
    }
    if (mMergedWriter) {
        delete mMergedWriter;
        mMergedWriter = NULL;
    }
    if (mFailedWriter) {
        delete mFailedWriter;
        mFailedWriter = NULL;
    }
    if (mUnpairedLeftWriter) {
        delete mUnpairedLeftWriter;
        mLeftWriter = NULL;
    }
    if (mUnpairedRightWriter) {
        delete mUnpairedRightWriter;
        mRightWriter = NULL;
    }

    if (mReadsKOWriter) {
        delete mReadsKOWriter;
        mReadsKOWriter = NULL;
    }
}

void PairEndProcessor::initConfig(ThreadConfig* config) {
    if (mOptions->out1.empty())
        return;
    if (mOptions->split.enabled) {
        config->initWriterForSplit();
    }
}

bool PairEndProcessor::process() {
    if (!mOptions->split.enabled)
        initOutput();

    initPackRepository();
    std::thread producer(std::bind(&PairEndProcessor::producerTask, this));
    //TODO: get the correct cycles
    int cycle = 151;
    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    std::thread** threads = new thread*[mOptions->thread];
    for (int t = 0; t < mOptions->thread; t++) {
        configs[t] = new ThreadConfig(mOptions, tbwtfmiDB, t, true);
        initConfig(configs[t]);
        threads[t] = new std::thread(std::bind(&PairEndProcessor::consumerTask, this, configs[t]));
    }

    std::thread* leftWriterThread = NULL;
    std::thread* rightWriterThread = NULL;
    std::thread* unpairedLeftWriterThread = NULL;
    std::thread* unpairedRightWriterThread = NULL;
    std::thread* mergedWriterThread = NULL;
    std::thread* failedWriterThread = NULL;
    std::thread* readsKOMapWriterThread = NULL;
    if (mLeftWriter)
        leftWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mLeftWriter));
    if (mRightWriter)
        rightWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mRightWriter));
    if (mUnpairedLeftWriter)
        unpairedLeftWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mUnpairedLeftWriter));
    if (mUnpairedRightWriter)
        unpairedRightWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mUnpairedRightWriter));
    if (mMergedWriter)
        mergedWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mMergedWriter));
    if (mFailedWriter)
        failedWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mFailedWriter));
    if (mReadsKOWriter)
        readsKOMapWriterThread = new std::thread(std::bind(&PairEndProcessor::writeTask, this, mReadsKOWriter));

    producer.join();
    for (int t = 0; t < mOptions->thread; t++) {
        threads[t]->join();
    }

    if (!mOptions->split.enabled) {
        if (leftWriterThread)
            leftWriterThread->join();
        if (rightWriterThread)
            rightWriterThread->join();
        if (unpairedLeftWriterThread)
            unpairedLeftWriterThread->join();
        if (unpairedRightWriterThread)
            unpairedRightWriterThread->join();
        if (mergedWriterThread)
            mergedWriterThread->join();
        if (failedWriterThread)
            failedWriterThread->join();
        if (readsKOMapWriterThread)
            readsKOMapWriterThread->join();
    }

    if (mOptions->verbose){
        mOptions->longlog ? loginfolong("start to generate reports\n") : loginfo("start to generate reports\n");
    }

    // merge stats and filter results
    vector<Stats*> preStats1;
    vector<Stats*> postStats1;
    vector<Stats*> preStats2;
    vector<Stats*> postStats2;
    vector<FilterResult*> filterResults;
//    vector< std::unordered_map<std::string, uint32 > > totalKoFreqVecResults;
//    totalKoFreqVecResults.reserve(mOptions->thread);
//    vector< std::unordered_map<std::string, std::unordered_map<std::string, double> > > totalOrgKOFreqVecResults;
//    totalOrgKOFreqVecResults.reserve(mOptions->thread);
//    vector<std::unordered_map<std::string, uint32> > totalGoFreqVecResults;
//    totalGoFreqVecResults.reserve(mOptions->thread);

    vector<std::map<const uint32 *, uint32> > totalIdFreqVecResults;
    totalIdFreqVecResults.reserve(mOptions->thread);
    for (int t = 0; t < mOptions->thread; t++) {
        preStats1.push_back(configs[t]->getPreStats1());
        postStats1.push_back(configs[t]->getPostStats1());
        preStats2.push_back(configs[t]->getPreStats2());
        postStats2.push_back(configs[t]->getPostStats2());
        filterResults.push_back(configs[t]->getFilterResult());
//        totalKoFreqVecResults.push_back(configs[t]->getTransSearcher()->getSubKoFreqUMap());
//        totalOrgKOFreqVecResults.push_back(configs[t]->getTransSearcher()->getSubOrgKOAbunUMap());
//        totalGoFreqVecResults.push_back(configs[t]->getTransSearcher()->getSubGoFreqUMap());
        totalIdFreqVecResults.push_back(configs[t]->getTransSearcher()->getIdFreqSubMap());
    }
    Stats* finalPreStats1 = Stats::merge(preStats1);
    Stats* finalPostStats1 = Stats::merge(postStats1);
    Stats* finalPreStats2 = Stats::merge(preStats2);
    Stats* finalPostStats2 = Stats::merge(postStats2);
    FilterResult* finalFilterResult = FilterResult::merge(filterResults);
    mOptions->mHomoSearchOptions.nTotalReads = finalPreStats1->getReads(); //change to both reads??????
    mOptions->mHomoSearchOptions.nCleanReads = finalPostStats1->getReads();
    
    mOptions->transSearch.totalIdFreqUMapResults = TransSearcher::merge(totalIdFreqVecResults);
    mOptions->transSearch.nTransMappedIds = mOptions->transSearch.totalIdFreqUMapResults.size();
    totalIdFreqVecResults.clear();
    
    prepareResults();
    //prepareResults(totalKoFreqVecResults, totalOrgKOFreqVecResults, totalGoFreqVecResults, totalIdFreqVecResults);

    int* dupHist = NULL;
    double* dupMeanTlen = NULL;
    double* dupMeanGC = NULL;
    double dupRate = 0.0;
    if (mOptions->duplicate.enabled) {
        dupHist = new int[mOptions->duplicate.histSize];
        memset(dupHist, 0, sizeof (int) * mOptions->duplicate.histSize);
        dupMeanGC = new double[mOptions->duplicate.histSize];
        memset(dupMeanGC, 0, sizeof (double) * mOptions->duplicate.histSize);
        dupRate = mDuplicate->statAll(dupHist, dupMeanGC, mOptions->duplicate.histSize);
        cerr << endl;
        cerr << "Duplication rate: " << dupRate * 100.0 << "%" << endl;
    }

    // insert size distribution
    int peakInsertSize = getPeakInsertSize();
    cerr << endl;
    cerr << "Insert size peak (evaluated by paired-end reads): " << peakInsertSize << endl;

    if (mOptions->merge.enabled) {
        //        cerr << endl;
        //        cerr << "Read pairs merged: " << finalFilterResult->mMergedPairs << endl;
        if (finalPostStats1->getReads() > 0) {
            double postMergedPercent = 100.0 * finalFilterResult->mMergedPairs / finalPostStats1->getReads();
            double preMergedPercent = 100.0 * finalFilterResult->mMergedPairs / finalPreStats1->getReads();
            //            cerr << "% of original read pairs: " << preMergedPercent << "%" << endl;
            //            cerr << "% in reads after filtering: " << postMergedPercent << "%" << endl;
        }
        cerr << endl;
    }

    JsonReporter jr(mOptions);
    jr.setDupHist(dupHist, dupMeanGC, dupRate);
    jr.setInsertHist(mInsertSizeHist, peakInsertSize);
    jr.report(finalFilterResult, finalPreStats1, finalPostStats1, finalPreStats2, finalPostStats2);

    // make HTML report
    HtmlReporter hr(mOptions);
    hr.setDupHist(dupHist, dupMeanGC, dupRate);
    hr.setInsertHist(mInsertSizeHist, peakInsertSize);
    hr.report(finalFilterResult, finalPreStats1, finalPostStats1, finalPreStats2, finalPostStats2);

    // clean up
    for (int t = 0; t < mOptions->thread; t++) {
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
    }

    delete finalPreStats1;
    delete finalPostStats1;
    delete finalPreStats2;
    delete finalPostStats2;
    delete finalFilterResult;

    if (mOptions->duplicate.enabled) {
        delete[] dupHist;
        delete[] dupMeanGC;
    }

    delete[] threads;
    delete[] configs;

    if (leftWriterThread)
        delete leftWriterThread;
    if (rightWriterThread)
        delete rightWriterThread;
    if (unpairedLeftWriterThread)
        delete unpairedLeftWriterThread;
    if (unpairedRightWriterThread)
        delete unpairedRightWriterThread;
    if (mergedWriterThread)
        delete mergedWriterThread;
    if (failedWriterThread)
        delete failedWriterThread;
    if (readsKOMapWriterThread)
        delete readsKOMapWriterThread;

    if (!mOptions->split.enabled)
        closeOutput();

    return true;
}

int PairEndProcessor::getPeakInsertSize() {
    int peak = 0;
    long maxCount = -1;
    for (int i = 0; i < mOptions->insertSizeMax; i++) {
        if (mInsertSizeHist[i] > maxCount) {
            peak = i;
            maxCount = mInsertSizeHist[i];
        }
    }
    return peak;
}

bool PairEndProcessor::processPairEnd(ReadPairPack* pack, ThreadConfig* config) {
    string * outstr1 = new string();
    string * outstr2 = new string();
    string unpairedOut1;
    string unpairedOut2;
    string singleOutput;
    string mergedOutput;
    string failedOut;
    std::string * outReadsKOMapStr = new string();
    //std::string koTag = "";
    uint32 * orthId = NULL;
    int mappedReads = 0;
    std::set<uint32 *> idSet;

    int readPassed = 0;
    int mergedCount = 0;
    for (int p = 0; p < pack->count; p++) {
        ReadPair* pair = pack->data[p];
        Read* or1 = pair->mLeft;
        Read* or2 = pair->mRight;

        int lowQualNum1 = 0;
        int nBaseNum1 = 0;
        int lowQualNum2 = 0;
        int nBaseNum2 = 0;

        // stats the original read before trimming
        config->getPreStats1()->statRead(or1);
        config->getPreStats2()->statRead(or2);

        // handling the duplication profiling
        if (mDuplicate)
            mDuplicate->statPair(or1, or2);

        // filter by index
        if (mOptions->indexFilter.enabled && mFilter->filterByIndex(or1, or2)) {
            delete pair;
            continue;
        }
        
        if(mOptions->fixMGI) {
            or1->fixMGI();
            or2->fixMGI();
        }
        
        // umi processing
        if (mOptions->umi.enabled)
            mUmiProcessor->process(or1, or2);

        // trim in head and tail, and apply quality cut in sliding window
        int frontTrimmed1 = 0;
        int frontTrimmed2 = 0;
        Read* r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1, frontTrimmed1);
        Read* r2 = mFilter->trimAndCut(or2, mOptions->trim.front2, mOptions->trim.tail2, frontTrimmed2);

        if (r1 != NULL && r2 != NULL) {
            if (mOptions->polyGTrim.enabled)
                PolyX::trimPolyG(r1, r2, config->getFilterResult(), mOptions->polyGTrim.minLen);
        }
        bool isizeEvaluated = false;
        if (r1 != NULL && r2 != NULL && (mOptions->adapter.enabled || mOptions->correction.enabled)) {
            OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire, mOptions->overlapDiffPercentLimit / 100.0);
            // we only use thread 0 to evaluae ISIZE
            if (config->getThreadId() == 0) {
                statInsertSize(r1, r2, ov, frontTrimmed1, frontTrimmed2);
                isizeEvaluated = true;
            }
            if (mOptions->correction.enabled) {
                BaseCorrector::correctByOverlapAnalysis(r1, r2, config->getFilterResult(), ov);
            }
            if (mOptions->adapter.enabled) {
                bool trimmed = AdapterTrimmer::trimByOverlapAnalysis(r1, r2, config->getFilterResult(), ov, frontTrimmed1, frontTrimmed2);
                bool trimmed1 = trimmed;
                bool trimmed2 = trimmed;
                if (!trimmed) {
                    if (mOptions->adapter.hasSeqR1)
                        trimmed1 = AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.sequence, false);
                    if (mOptions->adapter.hasSeqR2)
                        trimmed2 = AdapterTrimmer::trimBySequence(r2, config->getFilterResult(), mOptions->adapter.sequenceR2, true);
                }
                if (mOptions->adapter.hasFasta) {
                    AdapterTrimmer::trimByMultiSequences(r1, config->getFilterResult(), mOptions->adapter.seqsInFasta, false, !trimmed1);
                    AdapterTrimmer::trimByMultiSequences(r2, config->getFilterResult(), mOptions->adapter.seqsInFasta, true, !trimmed2);
                }

                if (mOptions->adapter.polyA) {
                    AdapterTrimmer::trimPolyA(r1, config->getFilterResult(), false, !trimmed1);
                    AdapterTrimmer::trimPolyA(r2, config->getFilterResult(), true, !trimmed2);
                }
            }
        }

        if (config->getThreadId() == 0 && !isizeEvaluated && r1 != NULL && r2 != NULL) {
            OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire, mOptions->overlapDiffPercentLimit / 100.0);
            statInsertSize(r1, r2, ov, frontTrimmed1, frontTrimmed2);
            isizeEvaluated = true;
        }

        if (r1 != NULL && r2 != NULL) {
            if (mOptions->polyXTrim.enabled)
                PolyX::trimPolyX(r1, r2, config->getFilterResult(), mOptions->polyXTrim.minLen);
        }

        if (r1 != NULL && r2 != NULL) {
            if (mOptions->trim.maxLen1 > 0 && mOptions->trim.maxLen1 < r1->length())
                r1->resize(mOptions->trim.maxLen1);
            if (mOptions->trim.maxLen2 > 0 && mOptions->trim.maxLen2 < r2->length())
                r2->resize(mOptions->trim.maxLen2);
        }

        Read* merged = NULL;
        // merging mode
        bool mergeProcessed = false;

        if (!mergeProcessed) {

            int result1 = mFilter->passFilter(r1);
            int result2 = mFilter->passFilter(r2);

            config->addFilterResult(max(result1, result2), 2);

            if (r1 != NULL && result1 == PASS_FILTER && r2 != NULL && result2 == PASS_FILTER) {

                if(orthId != NULL) orthId = NULL;
                if (mOptions->outputToSTDOUT && !mOptions->merge.enabled) {
                    //singleOutput += r1->toString() + r2->toString();
                } else {
                    OverlapResult ov = OverlapAnalysis::analyze(r1, r2, mOptions->overlapDiffLimit, mOptions->overlapRequire, mOptions->overlapDiffPercentLimit / 100.0);
                    if (ov.overlapped) {
                        merged = OverlapAnalysis::merge(r1, r2, ov);
                        int result = mFilter->passFilter(merged);
                        if (result == PASS_FILTER) {
                            config->getTransSearcher()->transSearch(merged, orthId);
                        } else {
                            config->getTransSearcher()->transSearch(r1, r2, orthId);
                        }
                        delete merged;
                    } else {
                        config->getTransSearcher()->transSearch(r1, r2, orthId);
                    }
                }
                
                if (orthId != NULL) {
                    if(mOptions->verbose){
                       idSet.insert(orthId); 
                    }
                    mappedReads++;
                    if (mLeftWriter && mRightWriter) {
                        *outstr1 += r1->toStringWithTag(orthId);
                        *outstr2 += r2->toStringWithTag(orthId);
                    }
                    if (mReadsKOWriter) {
                        *outReadsKOMapStr += trimName(r1->mName) + "\t" + "s2f_" + std::to_string(*orthId) + "\n";
                    }
                }
                // stats the read after filtering
                if (!mOptions->merge.enabled) {
                    config->getPostStats1()->statRead(r1);
                    config->getPostStats2()->statRead(r2);
                }
                readPassed++;
            }
        }

        delete pair;
        // if no trimming applied, r1 should be identical to or1
        if (r1 != or1 && r1 != NULL)
            delete r1;
        // if no trimming applied, r1 should be identical to or1
        if (r2 != or2 && r2 != NULL)
            delete r2;
    }
    
    if (mOptions->verbose) {
        mOptions->transSearch.nTransMappedIdReads += mappedReads;
        logMtx.lock();
        auto rCount = long(mOptions->transSearch.nTransMappedIdReads);
//        auto kCount = mOptions->transSearch.koUSet.size();
//        auto gCount = mOptions->transSearch.goUSet.size();
        mOptions->transSearch.idUSet.insert(idSet.begin(), idSet.end());
        auto iCount = mOptions->transSearch.idUSet.size();
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
    if (!mOptions->split.enabled)
        mOutputMtx.lock();
    if (mOptions->outputToSTDOUT) {
        // STDOUT output
        // if it's merging mode, write the merged reads to STDOUT
        // otherwise write interleaved single output
        if (mOptions->merge.enabled)
            fwrite(mergedOutput.c_str(), 1, mergedOutput.length(), stdout);
        else
            fwrite(singleOutput.c_str(), 1, singleOutput.length(), stdout);
    } else if (mOptions->split.enabled) {
        // split output by each worker thread
        if (!mOptions->out1.empty())
            config->getWriter1()->writeString(*outstr1);
        if (!mOptions->out2.empty())
            config->getWriter2()->writeString(*outstr2);
    }

    if (mMergedWriter && !mergedOutput.empty()) {
        // write merged data
        char* mdata = new char[mergedOutput.size()];
        memcpy(mdata, mergedOutput.c_str(), mergedOutput.size());
        mMergedWriter->input(mdata, mergedOutput.size());
    }

    if (mFailedWriter && !failedOut.empty()) {
        // write failed data
        char* fdata = new char[failedOut.size()];
        memcpy(fdata, failedOut.c_str(), failedOut.size());
        mFailedWriter->input(fdata, failedOut.size());
    }

    // normal output by left/right writer thread
    if (mRightWriter && mLeftWriter && (!outstr1->empty() || !outstr2->empty())) {
        // write PE
        char* ldata = new char[outstr1->size()];
        memcpy(ldata, outstr1->c_str(), outstr1->size());
        mLeftWriter->input(ldata, outstr1->size());

        char* rdata = new char[outstr2->size()];
        memcpy(rdata, outstr2->c_str(), outstr2->size());
        mRightWriter->input(rdata, outstr2->size());
    } else if (mLeftWriter && !singleOutput.empty()) {
        // write singleOutput
        char* ldata = new char[singleOutput.size()];
        memcpy(ldata, singleOutput.c_str(), singleOutput.size());
        mLeftWriter->input(ldata, singleOutput.size());
    }
    // output unpaired reads
    if (!unpairedOut1.empty() || !unpairedOut2.empty()) {
        if (mUnpairedLeftWriter && mUnpairedRightWriter) {
            // write PE
            char* unpairedData1 = new char[unpairedOut1.size()];
            memcpy(unpairedData1, unpairedOut1.c_str(), unpairedOut1.size());
            mUnpairedLeftWriter->input(unpairedData1, unpairedOut1.size());

            char* unpairedData2 = new char[unpairedOut2.size()];
            memcpy(unpairedData2, unpairedOut2.c_str(), unpairedOut2.size());
            mUnpairedRightWriter->input(unpairedData2, unpairedOut2.size());
        } else if (mUnpairedLeftWriter) {
            char* unpairedData = new char[unpairedOut1.size() + unpairedOut2.size() ];
            memcpy(unpairedData, unpairedOut1.c_str(), unpairedOut1.size());
            memcpy(unpairedData + unpairedOut1.size(), unpairedOut2.c_str(), unpairedOut2.size());
            mUnpairedLeftWriter->input(unpairedData, unpairedOut1.size() + unpairedOut2.size());
        }
    }

    if (mReadsKOWriter && !outReadsKOMapStr->empty()) {
        char* tdata = new char[outReadsKOMapStr->size()];
        memcpy(tdata, outReadsKOMapStr->c_str(), outReadsKOMapStr->size());
        mReadsKOWriter->input(tdata, outReadsKOMapStr->size());
    }

    if (!mOptions->split.enabled)
        mOutputMtx.unlock();

    if (mOptions->split.byFileLines)
        config->markProcessed(readPassed);
    else
        config->markProcessed(pack->count);

    if (mOptions->merge.enabled) {
        config->addMergedPairs(mergedCount);
    }
    
    if(outstr1){
        delete outstr1;
        outstr1 = NULL;
    }
    
    if(outstr2){
        delete outstr2;
        outstr2 = NULL;
    }
    
    if(outReadsKOMapStr){
        delete outReadsKOMapStr;
        outReadsKOMapStr = NULL;
    }

    delete pack->data;
    delete pack;

    return true;
}

void PairEndProcessor::statInsertSize(Read* r1, Read* r2, OverlapResult& ov, int frontTrimmed1, int frontTrimmed2) {
    int isize = mOptions->insertSizeMax;
    if (ov.overlapped) {
        if (ov.offset > 0)
            isize = r1->length() + r2->length() - ov.overlap_len + frontTrimmed1 + frontTrimmed2;
        else
            isize = ov.overlap_len + frontTrimmed1 + frontTrimmed2;
    }

    if (isize > mOptions->insertSizeMax)
        isize = mOptions->insertSizeMax;

    mInsertSizeHist[isize]++;
}

bool PairEndProcessor::processRead(Read* r, ReadPair* originalPair, bool reversed) {
    // do something here
    return true;
}

void PairEndProcessor::initPackRepository() {
    mRepo.packBuffer = new ReadPairPack*[PACK_NUM_LIMIT];
    memset(mRepo.packBuffer, 0, sizeof (ReadPairPack*) * PACK_NUM_LIMIT);
    mRepo.writePos = 0;
    mRepo.readPos = 0;

}

void PairEndProcessor::destroyPackRepository() {
    if(mRepo.packBuffer) delete mRepo.packBuffer; mRepo.packBuffer = NULL;
}

void PairEndProcessor::producePack(ReadPairPack* pack) {
    mRepo.packBuffer[mRepo.writePos] = pack;
    mRepo.writePos++;
}

void PairEndProcessor::consumePack(ThreadConfig* config) {
    ReadPairPack* data;
    mInputMtx.lock();
    while (mRepo.writePos <= mRepo.readPos) {
        usleep(1000);
        if (mProduceFinished) {
            mInputMtx.unlock();
            return;
        }
    }
    data = mRepo.packBuffer[mRepo.readPos];
    mRepo.readPos++;
    mInputMtx.unlock();
    processPairEnd(data, config);
}

void PairEndProcessor::producerTask() {
    if (mOptions->verbose){
        mOptions->longlog ? loginfolong("start to load data") : loginfo("start to load data");
    }
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    bool splitSizeReEvaluated = false;
    ReadPair** data = new ReadPair*[PACK_SIZE];
    memset(data, 0, sizeof (ReadPair*) * PACK_SIZE);
    FastqReaderPair reader(mOptions->in1, mOptions->in2, true, mOptions->phred64, mOptions->interleavedInput, mOptions->fastqBufferSize);
    int count = 0;
    bool needToBreak = false;
    while (true) {
        ReadPair* read = reader.read();
        // TODO: put needToBreak here is just a WAR for resolve some unidentified dead lock issue 
        if (!read || needToBreak) {
            // the last pack
            ReadPairPack* pack = new ReadPairPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            data = NULL;
            if (read) {
                delete read;
                read = NULL;
            }
            break;
        }
        data[count] = read;
        count++;
        // configured to process only first N reads
        if (mOptions->readsToProcess > 0 && count + readNum >= mOptions->readsToProcess) {
            needToBreak = true;
        }
        if (mOptions->verbose && count + readNum >= lastReported + 1000000) {
            lastReported = count + readNum;
            string msg = "\nloaded " + to_string((lastReported / 1000000)) + "M read pairs";
            mOptions->longlog ? loginfolong(msg) : loginfo(msg);
        }
        // a full pack
        if (count == PACK_SIZE || needToBreak) {
            ReadPairPack* pack = new ReadPairPack;
            pack->data = data;
            pack->count = count;           
            producePack(pack);
            //re-initialize data for next pack
            data = new ReadPair*[PACK_SIZE];
            memset(data, 0, sizeof (ReadPair*) * PACK_SIZE);
            // if the consumer is far behind this producer, sleep and wait to limit memory usage
            while (mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT) {
                slept++;
                usleep(1000);
            }
            readNum += count;
            // if the writer threads are far behind this producer, sleep and wait
            // check this only when necessary
            if (readNum % (PACK_SIZE * PACK_IN_MEM_LIMIT) == 0 && mLeftWriter) {
                while ((mLeftWriter && mLeftWriter->bufferLength() > PACK_IN_MEM_LIMIT) || (mRightWriter && mRightWriter->bufferLength() > PACK_IN_MEM_LIMIT)) {
                    slept++;
                    usleep(1000);
                }
            }
            // reset count to 0
            count = 0;
        }
    }

    mProduceFinished = true;
    if (mOptions->verbose){
        mOptions->longlog ? loginfolong("all reads loaded, start to monitor thread status") : loginfo("all reads loaded, start to monitor thread status");
    }

    // if the last data initialized is not used, free it
    if (data != NULL)
        delete[] data;
}

void PairEndProcessor::consumerTask(ThreadConfig* config) {
    while (true) {            
        if (config->canBeStopped()) {
            mFinishedThreads++;
            break;
        }
        while (mRepo.writePos <= mRepo.readPos) {
            if (mProduceFinished)
                break;
            usleep(1000);
        }

        if (mProduceFinished && mRepo.writePos == mRepo.readPos) {
            mFinishedThreads++;
            if (mOptions->verbose) {
                string msg = "\nthread " + to_string(config->getThreadId() + 1) + " data processing completed";
                 //mOptions->longlog ? loginfolong(msg) : loginfo(msg);
            }
            break;
        }
        if (mProduceFinished) {
            if (mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " is processing the " + to_string(mRepo.readPos) + " / " + to_string(mRepo.writePos) + " pack";
                // mOptions->longlog ? loginfolong(msg) : loginfo(msg);
            }
            consumePack(config);
            //lock.unlock();
        } else {
            //lock.unlock();
            consumePack(config);
        }
    }

    if (mFinishedThreads == mOptions->thread) {
        if (mLeftWriter)
            mLeftWriter->setInputCompleted();
        if (mRightWriter)
            mRightWriter->setInputCompleted();
        if (mUnpairedLeftWriter)
            mUnpairedLeftWriter->setInputCompleted();
        if (mUnpairedRightWriter)
            mUnpairedRightWriter->setInputCompleted();
        if (mMergedWriter)
            mMergedWriter->setInputCompleted();
        if (mFailedWriter)
            mFailedWriter->setInputCompleted();
        if (mReadsKOWriter)
            mReadsKOWriter->setInputCompleted();
    }

    if (mOptions->verbose) {
        string msg = "\nthread " + to_string(config->getThreadId() + 1) + " finished";
        mOptions->longlog ? loginfolong(msg) : loginfo(msg);
    }
}

void PairEndProcessor::writeTask(WriterThread* config) {
    while (true) {
        if (config->isCompleted()) {
            // last check for possible threading related issue
            config->output();
            break;
        }
        config->output();
    }

    if (mOptions->verbose) {
        string msg = config->getFilename() + " writer finished";
        mOptions->longlog ? loginfolong(msg) : loginfo(msg);
    }
}

void PairEndProcessor::prepareResults() {

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
        mOptions->transSearch.nMappedCoreOrthos = 0;
        for(const auto & it : mOptions->transSearch.totalIdFreqUMapResults){
            mOptions->transSearch.nTransMappedIdReads += it.second;
            auto itt = mOptions->mHomoSearchOptions.fullDbMap.find(it.first);
            if(itt != mOptions->mHomoSearchOptions.fullDbMap.end()){
                *fout << "s2f_" << *(it.first) << "\t" <<  it.second << "\t" << itt->second.ko << "|" << itt->second.go << "|" << itt->second.symbol << "|" << itt->second.gene << "\t" << itt->second.coreOrthoPer << "\n";
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
            std::string msg = "Finish to write s2f id abundance table!";
             mOptions->longlog ? loginfolong(msg) : loginfo(msg);
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