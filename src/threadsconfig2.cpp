#include "threadsconfig2.h"
#include "util.h"
#include <memory.h>
#include <unistd.h>

ThreadsConfig2::ThreadsConfig2(Options* opt, int threadId) {
    mOptions = opt;
    mThreadId = threadId;
    mWriter1 = NULL;
    mInputCounter = 0;
    mOutputCounter = 0;
    mInputCompleted = false;
    mRingBuffer = new char*[PACK_NUM_LIMIT];
    memset(mRingBuffer, 0, sizeof(char*) * PACK_NUM_LIMIT);
    mRingBufferSizes = new size_t[PACK_NUM_LIMIT];
    memset(mRingBufferSizes, 0, sizeof(char*) * PACK_NUM_LIMIT);
    string fullpath = "";
    if(mThreadId < mOptions->mSeqExtractions.targetGenesSubVec.size()){
        auto filename = mOptions->mSeqExtractions.targetGenesSubVec[mThreadId] + mOptions->mSeqExtractions.suffix;
        fullpath = joinpath(mOptions->mSeqExtractions.outputDir, filename);
    } else {
        fullpath = mOptions->mSeqExtractions.undeterminedFileNameOut;
    }
    initWriter(fullpath);
}

ThreadsConfig2::~ThreadsConfig2() {
    cleanup();
    if (mRingBuffer) {
        delete mRingBuffer;
        mRingBuffer = NULL;
    }
    if (mRingBufferSizes) {
        delete mRingBufferSizes;
        mRingBufferSizes = NULL;
    }
}

void ThreadsConfig2::initWriter(string filename1) {
    deleteWriter();
    mWriter1 = new Writer(filename1, mOptions->compression);
}

void ThreadsConfig2::initWriter(ofstream* stream) {
    deleteWriter();
    mWriter1 = new Writer(stream);
}

void ThreadsConfig2::initWriter(gzFile gzfile) {
    deleteWriter();
    mWriter1 = new Writer(gzfile);
}

void ThreadsConfig2::deleteWriter() {
    if(mWriter1 != NULL){
        delete mWriter1;
        mWriter1 = NULL;
    }
}

void ThreadsConfig2::cleanup(){
    deleteWriter();
}

bool ThreadsConfig2::isCompleted(){
    return mInputCompleted && (mOutputCounter == mInputCounter);
}

void ThreadsConfig2::setInputCompleted(){
    mInputCompleted = true;
}

void ThreadsConfig2::input(char* data, size_t size){
    long target = mInputCounter % PACK_NUM_LIMIT;
    mRingBuffer[target] = data;
    mRingBufferSizes[target] = size;
    mInputCounter++;
}

void ThreadsConfig2::output(){
    if(mOutputCounter >= mInputCounter){
        usleep(100);
    }
    while(mOutputCounter < mInputCounter){
        long target = mOutputCounter % PACK_NUM_LIMIT;
        mWriter1->write(mRingBuffer[target], mRingBufferSizes[target]);
        delete mRingBuffer[target];
        mRingBuffer[target] = NULL;
        mOutputCounter++;
    }
}