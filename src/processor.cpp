#include "processor.h"
#include "peprocessor.h"
#include "seprocessor.h"

Processor::Processor(Options* opt){
    mOptions = opt;
}

Processor::~Processor(){
}

bool Processor::process(BwtFmiDB * tbwtfmiDB) {
    if(mOptions->isPaired()) {
        PairEndProcessor p(mOptions, tbwtfmiDB);
        p.process();
    } else {
        SingleEndProcessor p(mOptions, tbwtfmiDB);
        p.process();
    }

    return true;
}