#ifndef THREADSCONFIG2_H
#define THREADSCONFIG2_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <atomic>
#include "writer.h"
#include "options.h"

using namespace std;

class ThreadsConfig2 {
public:
    ThreadsConfig2(Options* opt, int threadId);
    ~ThreadsConfig2();
    
    void initWriter(string filename1);
    void initWriter(ofstream* stream);
    void initWriter(gzFile gzfile);
    
    int getThreadId(){return mThreadId;};
    void cleanup();
    
    bool isCompleted();
    void output();
    void input(char* data, size_t size);
    void setInputCompleted();
    
private:
    void deleteWriter();
    
private:
    Writer* mWriter1;
    Options* mOptions;
    
    int mThreadId;
    bool mInputCompleted;
    atomic_long mInputCounter;
    atomic_long mOutputCounter;
    char** mRingBuffer;
    size_t* mRingBufferSizes;
};

#endif /* THREADSCONFIG2_H */

