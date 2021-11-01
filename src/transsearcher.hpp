#ifndef TRANSSEARCHER_HPP
#define TRANSSEARCHER_HPP

#include <stdint.h>
#include <assert.h>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <list>
#include <cmath>
#include <algorithm>
#include <mutex>
#include <iostream>
#include <sstream>
#include <vector>
#include <iterator>
#include <string>
#include <cstring>
#include <climits>
#include <map>
#include <utility>
#include <functional>
#include <locale>
#include <stdio.h>
#include <cmath>

#include "util.h"
#include "algo/blast/core/blast_seg.h"
#include "algo/blast/core/blast_filter.h"
#include "algo/blast/core/blast_encoding.h"
#include "read.h"
#include "options.h"
#include "fragment.h"
#include "options.h"
#include "bwtfmiDB.h"
#include "common.h"

extern "C" {
#include "bwt/bwt.h"
}

const double LN_2 = 0.6931471805;
const double LAMBDA = 0.3176;
const double LN_K = -2.009915479;

class TransSearcher {
protected:
    uint8_t codon_to_int(const char* codon);
    uint8_t revcomp_codon_to_int(const char* codon);

    uint8_t nuc2int[256];
    uint8_t compnuc2int[256];
    char codon2aa[256];
    uint8_t aa2int[256];

    std::map<char, std::vector<char>> blosum_subst;
    int8_t blosum62diag[20];
    int8_t b62[20][20];

    std::string translations[6];
    std::multimap<unsigned int, Fragment *, std::greater<unsigned int>> fragments;
    std::vector<SI *> best_matches_SI;
    std::vector<SI *> longest_matches_SI;
    std::vector<std::string> best_matches;
    std::vector<std::string> longest_fragments;
    
    unsigned int best_match_score = 0;
    
    std::string outputStr = "";
    std::string extraoutput = "";
    std::string extraoutputGO = "";
    std::string extraoutputId = "";
    double query_len;
    uint32_t read_count = 0;

    void clearFragments();
    unsigned int calcScore(const std::string &);
    unsigned int calcScore(const std::string &, int);
    unsigned int calcScore(const std::string &, size_t, size_t, int);
    void addAllMismatchVariantsAtPosSI(const Fragment *, unsigned int, size_t, SI *); // used in Greedy mode
    Fragment * getNextFragment(unsigned int);
    void eval_match_scores(SI *si, Fragment *);
    void getAllFragmentsBits(const std::string & line);
    void getLongestFragmentsBits(const std::string & line);
    void flush_output();
    void preProcess();
    void doProcess();
    void postProcess();

protected:
    void classify_length();
    void classify_greedyblosum();
    void ids_from_SI(SI *);
    void ids_from_SI_recursive(SI *);
    std::set<char *> match_ids;
    
    geneKoGoComb tmpGKG;
    std::vector<std::string> tmpKOVec;
    std::set<std::string> orgSet;
    std::unordered_set<std::string> koUSet;
//    std::pair<std::string, std::string> tmpReadKOPair;
    std::unordered_multimap<std::string, std::string> tmpKOProteinUMMap;
    std::unordered_multimap<std::string, std::pair<std::string, double> > tmpOrgKOAbunUMMap;
    std::multimap<std::string, double> tmpKOFreqMMap; //for org;
    std::unordered_map<std::string, double> tmpKOFreqUMap; //for org;
    std::unordered_map<std::string, uint32 > subKoFreqUMap;
    std::unordered_multimap<std::string, std::unordered_map<std::string, double> > subOrgKOAbunUMMap;
    std::unordered_map<std::string, std::unordered_map<std::string, double> > priOrgKOAbunUMap;
    std::unordered_map<std::string, std::vector<std::string> > subPathwayUMap;
    //std::multimap<std::string, std::pair<std::string, double> >  tmpOrgKOAbunMMap;
    
    std::vector<std::string> tmpGOVec;
    std::unordered_map<std::string, uint32 > subGoFreqUMap;
    std::unordered_set<std::string> goUSet;
    std::vector<std::string> tmpGoVec;
    
    std::vector<std::string> tmpIdVec;
    std::unordered_set<std::string> idUSet;
    std::unordered_map<std::string, uint32 > subIdFreqUMap;
    
    Options * mOptions;   
    BwtFmiDB * tbwtfmiDB;
    
public:
    TransSearcher(Options * opt, BwtFmiDB * tbwtfmiDB);
    std::string transSearch(Read * item);
    std::string transSearch(Read * item1, Read * item2);
    inline std::unordered_map<std::string, uint32 > getSubKoFreqUMap(){return subKoFreqUMap;};
    inline std::unordered_map<std::string, uint32 > getSubGoFreqUMap(){return subGoFreqUMap;};
    inline std::unordered_map<std::string, uint32 > getSubIdFreqUMap(){return subIdFreqUMap;};
    std::unordered_map<std::string, std::unordered_map<std::string, double> > getSubOrgKOAbunUMap();
    std::unordered_map<std::string, std::vector<std::string> > getSubPathwayUMap();
};


#endif /* TRANSSEARCHER_HPP */
