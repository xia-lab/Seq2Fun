#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>
#include <unordered_set>
#include <tuple>
#include <atomic>

#include "common.h"

using namespace std;

#define UMI_LOC_NONE 0
#define UMI_LOC_INDEX1 1
#define UMI_LOC_INDEX2 2
#define UMI_LOC_READ1 3
#define UMI_LOC_READ2 4
#define UMI_LOC_PER_INDEX 5
#define UMI_LOC_PER_READ 6

class MergeOptions {
public:

    MergeOptions() {
        enabled = false;
        includeUnmerged = false;
        out = "";
    }
public:
    bool enabled;
    bool includeUnmerged;
    string out;
};

class DuplicationOptions {
public:

    DuplicationOptions() {
        enabled = true;
        keylen = 12;
        histSize = 32;
    }
public:
    bool enabled;
    int keylen;
    int histSize;
};

class IndexFilterOptions {
public:

    IndexFilterOptions() {
        enabled = false;
        threshold = 0;
    }
public:
    vector<string> blacklist1;
    vector<string> blacklist2;
    bool enabled;
    int threshold;
};

class LowComplexityFilterOptions {
public:

    LowComplexityFilterOptions() {
        enabled = false;
        threshold = 0.3;
    }
public:
    bool enabled;
    double threshold;
};

class OverrepresentedSequenceAnasysOptions {
public:

    OverrepresentedSequenceAnasysOptions() {
        enabled = false;
        sampling = 20;
    }
public:
    bool enabled;
    int sampling;
};

class PolyGTrimmerOptions {
public:

    PolyGTrimmerOptions() {
        enabled = false;
        minLen = 10;
    }
public:
    bool enabled;
    int minLen;
};

class PolyXTrimmerOptions {
public:

    PolyXTrimmerOptions() {
        enabled = false;
        minLen = 10;
    }
public:
    bool enabled;
    int minLen;
};

class UMIOptions {
public:

    UMIOptions() {
        enabled = false;
        location = UMI_LOC_NONE;
        length = 0;
        skip = 0;
    }
public:
    bool enabled;
    int location;
    int length;
    int skip;
    string prefix;
    string separator;
};

class CorrectionOptions {
public:

    CorrectionOptions() {
        enabled = true;
    }
public:
    bool enabled;
};

class QualityCutOptions {
public:

    QualityCutOptions() {
        enabledFront = false;
        enabledTail = false;
        enabledRight = false;
        windowSizeShared = 4;
        qualityShared = 20;
        windowSizeFront = windowSizeShared;
        qualityFront = qualityShared;
        windowSizeTail = windowSizeShared;
        qualityTail = qualityShared;
        windowSizeRight = windowSizeShared;
        qualityRight = qualityShared;
    }
public:
    // enable 5' cutting by quality
    bool enabledFront;
    // enable 3' cutting by quality
    bool enabledTail;
    // enable agressive cutting mode
    bool enabledRight;
    // the sliding window size
    int windowSizeShared;
    // the mean quality requirement
    int qualityShared;
    // the sliding window size for cutting by quality in 5'
    int windowSizeFront;
    // the mean quality requirement for cutting by quality in 5'
    int qualityFront;
    // the sliding window size for cutting by quality in 3'
    int windowSizeTail;
    // the mean quality requirement for cutting by quality in 3'
    int qualityTail;
    // the sliding window size for cutting by quality in aggressive mode
    int windowSizeRight;
    // the mean quality requirement for cutting by quality in aggressive mode
    int qualityRight;
};

class SplitOptions {
public:

    SplitOptions() {
        enabled = false;
        needEvaluation = false;
        number = 0;
        size = 0;
        digits = 4;
        byFileNumber = false;
        byFileLines = false;
    }
public:
    bool enabled;
    // number of files
    int number;
    // lines of each file
    long size;
    // digits number of file name prefix, for example 0001 means 4 digits
    int digits;
    // need evaluation?
    bool needEvaluation;
    bool byFileNumber;
    bool byFileLines;
};

class AdapterOptions {
public:

    AdapterOptions() {
        enabled = true;
        polyA = true;
        hasSeqR1 = false;
        hasSeqR2 = false;
        detectAdapterForPE = false;
    }
public:
    bool enabled;
    bool polyA;
    string sequence;
    string sequenceR2;
    string detectedAdapter1;
    string detectedAdapter2;
    vector<string> seqsInFasta;
    string fastaFile;
    bool hasSeqR1;
    bool hasSeqR2;
    bool hasFasta;
    bool detectAdapterForPE;
};

class TrimmingOptions {
public:

    TrimmingOptions() {
        front1 = 0;
        tail1 = 0;
        front2 = 0;
        tail2 = 0;
        maxLen1 = 0;
        maxLen2 = 0;
    }
public:
    // trimming first cycles for read1
    int front1;
    // trimming last cycles for read1
    int tail1;
    // trimming first cycles for read2
    int front2;
    // trimming last cycles for read2
    int tail2;
    // max length of read1
    int maxLen1;
    // max length of read2
    int maxLen2;
};

class QualityFilteringOptions {
public:

    QualityFilteringOptions() {
        enabled = true;
        // '0' = Q15
        qualifiedQual = '0';
        unqualifiedPercentLimit = 40;
        nBaseLimit = 5;
    }
public:
    // quality filter enabled
    bool enabled;
    // if a base's quality phred score < qualifiedPhred, then it's considered as a low_qual_base
    char qualifiedQual;
    // if low_qual_base_num > lowQualLimit, then discard this read
    int unqualifiedPercentLimit;
    // if n_base_number > nBaseLimit, then discard this read
    int nBaseLimit;
    // if average qual score < avgQualReq, then discard this read
    int avgQualReq;
};

class ReadLengthFilteringOptions {
public:

    ReadLengthFilteringOptions() {
        enabled = false;
        requiredLength = 15;
        maxLength = 0;
    }
public:
    // length filter enabled
    bool enabled;
    // if read_length < requiredLength, then this read is discard
    int requiredLength;
    // length limit, 0 for no limitation
    int maxLength;
};

enum Mode {
    tMEM,
    tGREEDY
};

enum CodonTable {
    codontable1,
    codontable2,
    codontable3,
    codontable4,
    codontable5,
    codontable6,
    codontable9,
    codontable10,
    codontable12,
    codontable13,
    codontable14,
    codontable16,
    codontable21,
    codontable22,
    codontable24,
    codontable26,
    codontable27,
    codontable29,
    codontable30,
    codontable31,
    codontable33
};

class Sample{
public:
    Sample(){
        prefix = "";
        path = "";
        in1 = "";
        in2 = "";
        feature = "";

        nKO = 0;
        nKODb = 0;
        koRate = 0;
        
        nGO = 0;
        
        transSearchMappedKOReads = 0;
        transSearchMappedKOReads = 0;
        mappedKOReadsRate = 0;
        mappedGOReadsRate = 0;
        
        totalCleanReads = 0;
        cleanReadsRate = 0;
        totalRawReads = 0;
        timeLapse = 0;
        tMode = tGREEDY;
        tCodonTable = codontable1;
        nPathwaysDb = 0;
        nMappedPathways = 0;
        nOrgsDB = 0;
        nMappedOrgs = 0;
        totalKoFreqUMapResults.clear();
        totalGoFreqUMapResults.clear();
        rarefactionMap.clear();
        totalPathwayMap.clear();
        totalOrgKOUMap.clear();
        totoalReadsQualityVec.clear();
    }
    
public:
    string prefix;
    string path;
    string in1;
    string in2;
    string feature;
    int nKO;
    int nKODb;
    double koRate;
    int nGO;
    long transSearchMappedKOReads;
    long transSearchMappedGOReads;
    double mappedKOReadsRate;
    double mappedGOReadsRate;
    long totalCleanReads;
    double cleanReadsRate;
    long totalRawReads;
    int nPathwaysDb;
    int nMappedPathways;
    int nOrgsDB;
    int nMappedOrgs;
    time_t startTime;
    time_t endTime;
    long timeLapse;
    Mode tMode;
    CodonTable tCodonTable;
    std::unordered_map<std::string, uint32 > totalKoFreqUMapResults;
    std::unordered_map<std::string, uint32 > totalGoFreqUMapResults;
    std::map<long, int> rarefactionMap;
    std::unordered_map<std::string, int> totalPathwayMap;
    std::unordered_map<std::string, int> totalOrgKOUMap;
    std::vector<std::tuple<string, std::string, int> > totoalReadsQualityVec;
};

class TransSearchOptions {
public:
    TransSearchOptions() {
        mode = tGREEDY;
        codonTable = codontable1;
        SEG = true;
        useEvalue = false;
        minEvalue = 0.01;
        minAAFragLength = 0;

        misMatches = 2;
        minScore = 80;
        seedLength = 7;
        allFragments = false;

        size_t max_matches_SI = 10000;
        size_t max_match_ids = 10000;

        transSearchFinished = false;
        sampleKOAbunUMap.clear();
        koUSet.clear();
        goUSet.clear();
        nTransMappedKOReads = 0;
        nTransMappedGOReads = 0;
        nKODB = 0;
        nGODB = 0;
        nTransMappedKOs = 0;
        nTransMappedGOs = 0;
        sortedKOFreqTupleVector.clear();
        rarefactionMap.clear();
        sortedPathwayFreqTupleVector.clear();
        nPathwaysDB = 0;
        sortedOrgFreqVec.clear();
        nOrgsDB = 0;
        nMappedOrgs = 0;
        timeLapse = 0;
    }

    void reset2Default() {
        transSearchFinished = false;
        sampleKOAbunUMap.clear();
        koUSet.clear();
        goUSet.clear();
        nTransMappedKOReads = 0;
        nTransMappedGOReads = 0;
        nTransMappedKOs = 0;
        nTransMappedGOs = 0;
        sortedKOFreqTupleVector.clear();
        rarefactionMap.clear();
        sortedPathwayFreqTupleVector.clear();
        sortedOrgFreqVec.clear();
        nMappedOrgs = 0;
        timeLapse = 0;
    }
    
public:
    string tfmi;
    string tmode;
    string tCodonTable;
    bool SEG;
    bool useEvalue;
    double minEvalue;
    unsigned int minAAFragLength;
    unsigned int misMatches;
    unsigned int minScore;
    unsigned int seedLength;
    unsigned int maxTransLength;
    bool allFragments;

    size_t max_matches_SI;
    size_t max_match_ids;

    time_t startTime;
    time_t endTime;
    long timeLapse;
    Mode mode;
    CodonTable codonTable;
    std::unordered_set<std::string> koUSet;
    std::unordered_set<std::string> goUSet;
    atomic_long nTransMappedKOReads;
    long nTransMappedGOReads;
    unsigned int nKODB;
    unsigned int nGODB;
    unsigned int nTransMappedKOs;
    unsigned int nTransMappedGOs;
    int nMappedPathways;
    int nPathwaysDB;
    int nOrgsDB;
    int nMappedOrgs;
    std::vector< std::tuple <std::string, uint32, std::string> > sortedKOFreqTupleVector;
    std::vector<std::tuple<std::string, double, std::string, int, int> > sortedPathwayFreqTupleVector;
    std::map<long, int> rarefactionMap;
    std::vector<std::pair<std::string, int> > sortedOrgFreqVec;
    std::unordered_map<std::string, int> sampleKOAbunUMap;
    bool transSearchFinished;
};

class geneKoGoComb{
public:
    std::string ko;
    std::string go;
    std::string spec;
    std::string getKo(){ return ko;};
    std::string getGo(){ return go;};
    std::string getSpec(){ return spec;};
    
public:
    geneKoGoComb(){
        ko = "UNASSIGNED";
        go = "UNASSIGNED";
        spec = "UNASSIGNED";
    }
};

class HomoSearchOptions {
public:
    HomoSearchOptions() {
        profiling = false;
        genemap = "";
        prefix = "";
        sampleTable = "";
        pathway = "";
        genefa = "";
        nCleanReads = 0;
        nTotalReads = 0;
        commandStr = "";
    }

    void reset2Default() {
        nCleanReads = 0;
        nTotalReads = 0;
        commandStr = "";
    }

public:
    std::string genemap;
    std::string prefix;
    std::string sampleTable;
    std::string pathway;
    std::string genefa;
    bool profiling;
    long nCleanReads;
    long nTotalReads;

    std::map<std::string, geneKoGoComb> fullDbMap;
    //std::map<std::string, std::string> db_map;
    //std::map<std::string, std::string> org_map;
    std::multimap<std::string, std::string> pathway_ko_multimap;
    std::unordered_map<std::string, int> pathway_ko_stats_umap;
    std::map<std::string, std::string> ko_fullname_map;
    std::map<std::string, std::string> geneGoMap;

    std::ifstream filein;
    std::string fileName;
    std::string commandStr;
};

class SeqExtractions{
public:
    std::string sampleMappedTableStr;
    std::string targetGeneTableStr;
    std::string outputDir;
    int compression;
    std::unordered_set<std::string> numFeaturesProcessedUSet;
    std::string undeterminedFileName;
    std::string undeterminedFileNameIn;
    std::string undeterminedFileNameOut;
    
    std::vector<std::string> targetGenesVec;
    std::vector<std::string> targetGenesSubVec;
    std::vector<std::string> samplesVecF;
    std::vector<std::string> samplesVecR;
    std::vector<std::string> samplesVec;
    
    bool paired;
    std::string suffix;
    long numTotalReads;
    int numFeatures;
    
public:
    SeqExtractions(){
        sampleMappedTableStr = "";
        targetGeneTableStr = "";
        outputDir = "";
        compression = 6;
        targetGenesVec.clear();
        samplesVecF.clear();
        samplesVecR.clear();
        samplesVec.clear();
        undeterminedFileName = "";
        undeterminedFileNameIn = "";
        undeterminedFileNameOut = "";
        paired = false;
        suffix = "";
        numTotalReads = 0;
        numFeatures = 0;
        numFeaturesProcessedUSet.clear();
    }
};

class Options {
public:
    Options();
    void init();
    bool isPaired();
    bool validate();
    bool adapterCuttingEnabled();
    bool polyXTrimmingEnabled();
    string getAdapter1();
    string getAdapter2();
    void initIndexFiltering(string blacklistFile1, string blacklistFile2, int threshold = 0);
    vector<string> makeListFromFileByLine(string filename);
    bool shallDetectAdapter(bool isR2 = false);
    void loadFastaAdapters();
    void readDB();
    void parseSampleTable();
    void mkSelectedPathwayDB();
    void readSampleExtraction();

public:
    // file name of read1 input
    string in1;
    // file name of read2 input
    string in2;
    // file name of read1 output
    string out1;
    // file name of read2 output
    string out2;
    // enable output mapped reads
    bool outputMappedCleanReads;
    bool outputReadsKOMap;
    string outReadsKOMap;

    //seq2fun file dir;
    string seq2funProgPath;
    //seq2fun dir;
    string seq2funDir;
    
    string internalDBDir;
    //for internal database

    bool screenout;

    // file name of unpaired read1 output
    string unpaired1;
    // file name of unpaired read2 output
    string unpaired2;
    // file name of failed reads output
    string failedOut;
    // json file
    string jsonFile;
    // html file
    string htmlFile;
    // html report title
    string reportTitle;
    // compression level
    int compression;
    // the input file is using phred64 quality scoring
    bool phred64;
    // do not rewrite existing files
    bool dontOverwrite;
    // read STDIN
    bool inputFromSTDIN;
    // write STDOUT
    bool outputToSTDOUT;
    // the input R1 file is interleaved
    bool interleavedInput;
    // only process first N reads
    int readsToProcess;
    // fix the MGI ID tailing issue
    bool fixMGI;
    // worker thread number
    int thread;
    // fastq reads buffer size
    size_t fastqBufferSize;
    // trimming options
    TrimmingOptions trim;
    // quality filtering options
    QualityFilteringOptions qualfilter;
    // length filtering options
    ReadLengthFilteringOptions lengthFilter;
    // adapter options
    AdapterOptions adapter;
    // multiple file splitting options
    SplitOptions split;
    // options for quality cutting
    QualityCutOptions qualityCut;
    // options for base correction
    CorrectionOptions correction;
    // options for UMI
    UMIOptions umi;
    // 3' end polyG trimming, default for Illumina NextSeq/NovaSeq
    PolyGTrimmerOptions polyGTrim;
    // 3' end polyX trimming
    PolyXTrimmerOptions polyXTrim;
    // for overrepresentation analysis
    OverrepresentedSequenceAnasysOptions overRepAnalysis;
    map<string, long> overRepSeqs1;
    map<string, long> overRepSeqs2;
    int seqLen1;
    int seqLen2;
    // low complexity filtering
    LowComplexityFilterOptions complexityFilter;
    // black lists for filtering by index
    IndexFilterOptions indexFilter;
    // options for duplication profiling
    DuplicationOptions duplicate;
    // max value of insert size
    int insertSizeMax;
    // overlap analysis threshold
    int overlapRequire;
    int overlapDiffLimit;
    int overlapDiffPercentLimit;
    // output debug information
    bool verbose;
    bool debug;
    // merge options
    MergeOptions merge;
    TransSearchOptions transSearch;
    HomoSearchOptions mHomoSearchOptions;
    //BwtFmiDB tbwtfmiDB;
    std::vector<Sample> samples;
    int getWorkingSampleId(string & samplePrefix); 
    SeqExtractions mSeqExtractions;
};

#endif
