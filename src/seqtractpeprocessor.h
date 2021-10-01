#ifndef SEQTRACTPEPROCESSOR_H
#define SEQTRACTPEPROCESSOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <condition_variable>
#include <mutex>
#include <thread>
#include "options.h"
#include "fastqreader.h"
#include "threadsconfig2.h"
#include "read.h"
#include "util.h"


using namespace std;

struct ReadPack{
    Read** data;
    int count;
};
typedef struct ReadPack ReadPack;

struct ReadRepository{
    ReadPack** packBuffer;
    size_t readPos;
    size_t writePos;
    size_t readCounter;
    std::mutex mtx;
    std::mutex readCounterMtx;
    std::condition_variable repoNotFull;
    std::condition_variable repoNotEmpty;
};
typedef struct ReadPairRepository ReadPairRepository;


class SeqTractPeProcessor {
public:
    SeqTractPeProcessor(Options * opt);
    ~SeqTractPeProcessor();
    bool process();
    
private:
    bool processReads(ReadPack* pack);
    void initPackRepository();
    void destroyPackRepository();
    void producePack(ReadPack* pack);
    void consumePack();
    void producerTask();
    void consumerTask();
    void writeTask(ThreadsConfig2* config);
    
private:
    Options* mOptions;
    ReadRepository mRepo;
    bool mProduceFinished;
    ThreadsConfig2** mConfigs;
    int mSampleSize;
    std::unordered_set<std::string> featureUSet;
};

#endif /* SEQTRACTPEPROCESSOR_H */

