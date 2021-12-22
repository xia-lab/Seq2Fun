#ifndef BWTFMIDB_H
#define BWTFMIDB_H

#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include "util.h"
#include "options.h"

#include "include/ncbi-blast+/algo/blast/core/blast_seg.h"
#include "include/ncbi-blast+/algo/blast/core/blast_filter.h"
#include "include/ncbi-blast+/algo/blast/core/blast_encoding.h"

extern "C" {
#include "bwt/fmi.h"
#include "bwt/bwt.h"
#include "bwt/sequence.h"
}
using namespace std;

class BwtFmiDB {
public:
    BwtFmiDB(Options * & opt);
    ~BwtFmiDB();

    //for trans search
    BWT * tbwt;
    FMI * tfmi;
    AlphabetStruct * tastruct;
    SegParameters * tblast_seg_params;
    double tdb_length;
    bool Transsearch;
    
private:
    void init();
    
private:
    Options * mOptions;
};

#endif /* BWTFMIDB_H */

