#include "options.h"
#include "util.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include "fastareader.h"

Options::Options(){
    in1 = "";
    in2 = "";
    out1 = "";
    out2 = "";
    unpaired1 = "";
    unpaired2 = "";
    failedOut = "";
    screenout = false;
    reportTitle = "Seq2Fun report";
    thread = 2;
    compression = 2;
    phred64 = false;
    dontOverwrite = false;
    inputFromSTDIN = false;
    outputToSTDOUT = false;
    readsToProcess = 0;
    interleavedInput = false;
    insertSizeMax = 512;
    overlapRequire = 30;
    overlapDiffLimit = 5;
    overlapDiffPercentLimit = 20;
    verbose = false;
    debug = false;
    longlog = false;
    outputMappedCleanReads = false;
    outputReadsAnnoMap = false;
    outReadsKOMap = "";
    seqLen1 = 151;
    seqLen2 = 151;
    fixMGI = false;
    fastqBufferSize = 1<<20;
    samples.clear();
}

Options::~Options(){
//    for(const auto it : mHomoSearchOptions.idDbMap){
//        if(it.second){
//            delete it.second;
//        }
//    }
//    mHomoSearchOptions.idDbMap.clear();
    
//    if(mHomoSearchOptions != NULL) {
//        delete mHomoSearchOptions; 
//        mHomoSearchOptions = NULL;
//    }
}

void Options::init() {
}

bool Options::isPaired() {
    return in2.length() > 0 || interleavedInput;
}

bool Options::adapterCuttingEnabled() {
    if(adapter.enabled){
        if(isPaired() || !adapter.sequence.empty())
            return true;
    }
    return false;
}

bool Options::polyXTrimmingEnabled() {
    return polyXTrim.enabled;
}

void Options::loadFastaAdapters() {
    if(adapter.fastaFile.empty()) {
        adapter.hasFasta = false;
        return;
    }

    check_file_valid(adapter.fastaFile);

    FastaReader reader(adapter.fastaFile);
    reader.readAll();

    map<string, string> contigs = reader.contigs();
    map<string, string>::iterator iter;
    for(iter = contigs.begin(); iter != contigs.end(); iter++) {
        if(iter->second.length()>=6) {
            adapter.seqsInFasta.push_back(iter->second);
        }
        else {
            cerr << "skip too short adapter sequence in " <<  adapter.fastaFile << " (6bp required): " << iter->second << endl;
        }
    }

    if(adapter.seqsInFasta.size() > 0) {
        adapter.hasFasta = true;
    } else {
        adapter.hasFasta = false;
    }
}

bool Options::validate() {
    if(in1.empty()) {
        if(!in2.empty())
            error_exit("read2 input is specified by <in2>, but read1 input is not specified by <in1>");
        if(inputFromSTDIN)
            in1 = "/dev/stdin";
        else
            error_exit("read1 input should be specified by --in1, or enable --stdin if you want to read STDIN");
    } else {
        check_file_valid(in1);
    }

    if(!in2.empty()) {
        check_file_valid(in2);
    }
    
//    if (fmi.empty()) {
//        error_exit("protein database is not provided");
//    } else {
//        check_file_valid(fmi);
//    }
//
//    if (dbk.empty()) {
//        error_exit("protein KO database is not provided");
//    } else {
//        check_file_valid(dbk);
//    }
    
//    if (abd.length() <= 0){
//        error_exit("basename of output files is not provided, please specify the basename");
//    }
    

//    if(abd.length() > 0 & screenout == true){
//        error_exit("you must specify either output basename or screenout");
//    }
//    
//    if(screenout == true & profiling == true){
//        cerr << "You specified both --screenout and --profiling mode, profiling mode will be turn off" << endl;
//        profiling = false;
//    }
    
    if(merge.enabled) {
        if(split.enabled) {
            error_exit("splitting mode cannot work with merging mode");
        }
        if(in2.empty() && !interleavedInput) {
            error_exit("read2 input should be specified by --in2 for merging mode");
        }
        // enable correction if it's not enabled
        if(!correction.enabled)
            correction.enabled = true;
        if(merge.out.empty() && !outputToSTDOUT && !out1.empty() && out2.empty()) {
            cerr << "You specified --out1, but haven't specified --merged_out in merging mode. Using --out1 to store the merged reads" << endl << endl;
            merge.out = out1;
            out1 = "";
        }
        if(merge.includeUnmerged) {
            if(!out1.empty()) {
                cerr << "You specified --include_unmerged in merging mode. Ignoring argument --out1 = " << out1 << endl;
                out1 = "";
            }
            if(!out2.empty()) {
                cerr << "You specified --include_unmerged in merging mode. Ignoring argument --out2 = " << out2 << endl;
                out2 = "";
            }
            if(!unpaired1.empty()) {
                cerr << "You specified --include_unmerged in merging mode. Ignoring argument --unpaired1 = " << unpaired1 << endl;
                unpaired1 = "";
            }
            if(!unpaired2.empty()) {
                cerr << "You specified --include_unmerged in merging mode. Ignoring argument --unpaired1 = " << unpaired2 << endl;
                unpaired2 = "";
            }
        }
        if(merge.out.empty() && !outputToSTDOUT) {
            error_exit("In merging mode, you should either specify --merged_out or enable --stdout");
        }
        if(!merge.out.empty()) {
            if(merge.out == out1)
                error_exit("--merged_out and --out1 shouldn't have same file name");
            if(merge.out == out2)
                error_exit("--merged_out and --out2 shouldn't have same file name");
            if(merge.out == unpaired1)
                error_exit("--merged_out and --unpaired1 shouldn't have same file name");
            if(merge.out == unpaired2)
                error_exit("--merged_out and --unpaired2 shouldn't have same file name");
        }
    } else {
        // not in merging mode
        if(!merge.out.empty()) {
            cerr << "You haven't enabled merging mode (-m/--merge), ignoring argument --merged_out = " << merge.out << endl;
            merge.out = "";
        }
    }

    // if output to STDOUT, then...
    if(outputToSTDOUT) {
        if(split.enabled) {
            error_exit("splitting mode cannot work with stdout mode");
        }
        cerr << "Streaming uncompressed ";
        if(merge.enabled)
            cerr << "merged";
        else if(isPaired())
            cerr << "interleaved";
        cerr << " reads to STDOUT..." << endl;
        if(isPaired() && !merge.enabled)
            cerr << "Enable interleaved output mode for paired-end input." << endl;
        cerr << endl;
    }

    if(in2.empty() && !interleavedInput && !out2.empty()) {
        error_exit("read2 output is specified (--out2), but neighter read2 input is not specified (--in2), nor read1 is interleaved.");
    }

    if(!in2.empty() || interleavedInput) {
        if(!out1.empty() && out2.empty()) {
            error_exit("paired-end input, read1 output should be specified together with read2 output (--out2 needed) ");
        }
        if(out1.empty() && !out2.empty()) {
            if(!merge.enabled)
                error_exit("paired-end input, read1 output should be specified (--out1 needed) together with read2 output ");
        }
    }

    if(!in2.empty() && interleavedInput) {
        error_exit("<in2> is not allowed when <in1> is specified as interleaved mode by (--interleaved_in)");
    }

    if(!out1.empty()) {
        //check_file_writable(out1);
        if(out1 == out2) {
            error_exit("read1 output (--out1) and read1 output (--out2) should be different");
        }
        if(dontOverwrite && file_exists(out1)) {
            error_exit(out1 + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
    }
    if(!out2.empty()) {
        //check_file_writable(out2);
        if(dontOverwrite && file_exists(out2)) {
            error_exit(out2 + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
    }
    if(!isPaired()) {
        if(!unpaired1.empty()) {
            cerr << "Not paired-end mode. Ignoring argument --unpaired1 = " << unpaired1 << endl;
            unpaired1 = "";
        }
        if(!unpaired2.empty()) {
            cerr << "Not paired-end mode. Ignoring argument --unpaired2 = " << unpaired2 << endl;
            unpaired2 = "";
        }
    }
    if(split.enabled) {
        if(!unpaired1.empty()) {
            cerr << "Outputing unpaired reads is not supported in splitting mode. Ignoring argument --unpaired1 = " << unpaired1 << endl;
            unpaired1 = "";
        }
        if(!unpaired2.empty()) {
            cerr << "Outputing unpaired reads is not supported in splitting mode. Ignoring argument --unpaired2 = " << unpaired2 << endl;
            unpaired2 = "";
        }
    }
    if(!unpaired1.empty()) {
        if(dontOverwrite && file_exists(unpaired1)) {
            error_exit(unpaired1 + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
        if(unpaired1 == out1)
            error_exit("--unpaired1 and --out1 shouldn't have same file name");
        if(unpaired1 == out2)
            error_exit("--unpaired1 and --out2 shouldn't have same file name");
    }
    if(!unpaired2.empty()) {
        if(dontOverwrite && file_exists(unpaired2)) {
            error_exit(unpaired2 + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
        if(unpaired2 == out1)
            error_exit("--unpaired2 and --out1 shouldn't have same file name");
        if(unpaired2 == out2)
            error_exit("--unpaired2 and --out2 shouldn't have same file name");
    }
    if(!failedOut.empty()) {
        if(dontOverwrite && file_exists(failedOut)) {
            error_exit(failedOut + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
        if(failedOut == out1)
            error_exit("--failed_out and --out1 shouldn't have same file name");
        if(failedOut == out2)
            error_exit("--failed_out and --out2 shouldn't have same file name");
        if(failedOut == unpaired1)
            error_exit("--failed_out and --unpaired1 shouldn't have same file name");
        if(failedOut == unpaired2)
            error_exit("--failed_out and --unpaired2 shouldn't have same file name");
        if(failedOut == merge.out)
            error_exit("--failed_out and --merged_out shouldn't have same file name");
    }

    if(dontOverwrite) {
        if(file_exists(jsonFile)) {
            error_exit(jsonFile + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
        if(file_exists(htmlFile)) {
            error_exit(htmlFile + " already exists and you have set to not rewrite output files by --dont_overwrite");
        }
    }

    if(compression < 1 || compression > 9)
        error_exit("compression level (--compression) should be between 1 ~ 9, 1 for fastest, 9 for smallest");

    if(readsToProcess < 0)
        error_exit("the number of reads to process (--reads_to_process) cannot be negative");

    if(thread < 1) {
        thread = 1;
    } else if(thread > 32) {
        cerr << "WARNING: Seq2Fun uses up to 32 threads although you specified " << thread << endl;
        thread = 32;
    }

    if(trim.front1 < 0 || trim.front1 > 30)
        error_exit("trim_front1 (--trim_front1) should be 0 ~ 30, suggest 0 ~ 4");

    if(trim.tail1 < 0 || trim.tail1 > 100)
        error_exit("trim_tail1 (--trim_tail1) should be 0 ~ 100, suggest 0 ~ 4");

    if(trim.front2 < 0 || trim.front2 > 30)
        error_exit("trim_front2 (--trim_front2) should be 0 ~ 30, suggest 0 ~ 4");

    if(trim.tail2 < 0 || trim.tail2 > 100)
        error_exit("trim_tail2 (--trim_tail2) should be 0 ~ 100, suggest 0 ~ 4");

    if(qualfilter.qualifiedQual - 33 < 0 || qualfilter.qualifiedQual - 33 > 93)
        error_exit("qualitified phred (--qualified_quality_phred) should be 0 ~ 93, suggest 10 ~ 20");

    if(qualfilter.avgQualReq < 0 || qualfilter.avgQualReq  > 93)
        error_exit("average quality score requirement (--average_qual) should be 0 ~ 93, suggest 20 ~ 30");

    if(qualfilter.unqualifiedPercentLimit < 0 || qualfilter.unqualifiedPercentLimit > 100)
        error_exit("unqualified percent limit (--unqualified_percent_limit) should be 0 ~ 100, suggest 20 ~ 60");

    if(qualfilter.nBaseLimit < 0 || qualfilter.nBaseLimit > 50)
        error_exit("N base limit (--n_base_limit) should be 0 ~ 50, suggest 3 ~ 10");

    if(lengthFilter.requiredLength < 0 )
        error_exit("length requirement (--length_required) should be >0, suggest 15 ~ 100");

    if(overlapDiffPercentLimit < 0 || overlapDiffPercentLimit > 100)
        error_exit("the maximum percentage of mismatched bases to detect overlapped region (--overlap_diff_percent_limit) should be 0 ~ 100, suggest 20 ~ 60");

    if(split.enabled ) {
        if(split.digits < 0 || split.digits > 10)
            error_exit("you have enabled splitting output to multiple files, the digits number of file name prefix (--split_prefix_digits) should be 0 ~ 10.");
        
        if(split.byFileNumber) {
            if(split.number < 2 || split.number >= 1000)
                error_exit("you have enabled splitting output by file number, the number of files (--split) should be 2 ~ 999.");
            // thread number cannot be more than the number of file to split
            if(thread > split.number)
                thread = split.number;
        }

        if(split.byFileLines) {
            if(split.size < 1000/4)
                error_exit("you have enabled splitting output by file lines, the file lines (--split_by_lines) should be >= 1000.");
        }
    }

    if(qualityCut.enabledFront || qualityCut.enabledTail || qualityCut.enabledRight) {
        if(qualityCut.windowSizeShared < 1 || qualityCut.windowSizeShared > 1000)
            error_exit("the sliding window size for cutting by quality (--cut_window_size) should be between 1~1000.");
        if(qualityCut.qualityShared < 1 || qualityCut.qualityShared > 30)
            error_exit("the mean quality requirement for cutting by quality (--cut_mean_quality) should be 1 ~ 30, suggest 15 ~ 20.");
        if(qualityCut.windowSizeFront < 1 || qualityCut.windowSizeFront > 1000)
            error_exit("the sliding window size for cutting by quality (--cut_front_window_size) should be between 1~1000.");
        if(qualityCut.qualityFront < 1 || qualityCut.qualityFront > 30)
            error_exit("the mean quality requirement for cutting by quality (--cut_front_mean_quality) should be 1 ~ 30, suggest 15 ~ 20.");
        if(qualityCut.windowSizeTail < 1 || qualityCut.windowSizeTail > 1000)
            error_exit("the sliding window size for cutting by quality (--cut_tail_window_size) should be between 1~1000.");
        if(qualityCut.qualityTail < 1 || qualityCut.qualityTail > 30)
            error_exit("the mean quality requirement for cutting by quality (--cut_tail_mean_quality) should be 1 ~ 30, suggest 13 ~ 20.");
        if(qualityCut.windowSizeRight < 1 || qualityCut.windowSizeRight > 1000)
            error_exit("the sliding window size for cutting by quality (--cut_right_window_size) should be between 1~1000.");
        if(qualityCut.qualityRight < 1 || qualityCut.qualityRight > 30)
            error_exit("the mean quality requirement for cutting by quality (--cut_right_mean_quality) should be 1 ~ 30, suggest 15 ~ 20.");
    }

    if(adapter.sequence!="auto" && !adapter.sequence.empty()) {
        // validate adapter sequence for single end adapter trimming
        if(adapter.sequence.length() <= 3)
            error_exit("the sequence of <adapter_sequence> should be longer than 3");

        // validate bases
        for(int i=0; i<adapter.sequence.length(); i++) {
            char c = adapter.sequence[i];
            if(c!='A' && c!='T' && c!='C' && c!='G') {
                error_exit("the adapter <adapter_sequence> can only have bases in {A, T, C, G}, but the given sequence is: " + adapter.sequence);
            }
        }

        adapter.hasSeqR1 = true;
    }

    if(adapter.sequenceR2!="auto" && !adapter.sequenceR2.empty()) {
        // validate adapter sequenceR2 for single end adapter trimming
        if(adapter.sequenceR2.length() <= 3)
            error_exit("the sequence of <adapter_sequence_r2> should be longer than 3");

        // validate bases
        for(int i=0; i<adapter.sequenceR2.length(); i++) {
            char c = adapter.sequenceR2[i];
            if(c!='A' && c!='T' && c!='C' && c!='G') {
                error_exit("the adapter <adapter_sequence_r2> can only have bases in {A, T, C, G}, but the given sequenceR2 is: " + adapter.sequenceR2);
            }
        }

        adapter.hasSeqR2 = true;
    }

    if(correction.enabled && !isPaired()) {
        cerr << "WARNING: base correction is only appliable for paired end data, ignoring -c/--correction" << endl;
        correction.enabled = false;
    }

    if(umi.enabled) {
        if(umi.location == UMI_LOC_READ1 || umi.location == UMI_LOC_READ2 || umi.location == UMI_LOC_PER_READ) {
            if(umi.length<1 || umi.length>100)
                error_exit("UMI length should be 1~100");
            if(umi.skip<0 || umi.skip>100)
                error_exit("The base number to skip after UMI <umi_skip> should be 0~100");
        }else {
            if(umi.skip>0)
                error_exit("Only if the UMI location is in read1/read2/per_read, you can skip bases after UMI");
            if(umi.length>0)
                error_exit("Only if the UMI location is in read1/read2/per_read, you can set the UMI length");
        }
        if(!umi.prefix.empty()) {
            if(umi.prefix.length() >= 10)
                error_exit("UMI prefix should be shorter than 10");
            for(int i=0; i<umi.prefix.length(); i++) {
                char c = umi.prefix[i];
                if( !(c>='A' && c<='Z') && !(c>='a' && c<='z') && !(c>='0' && c<='9')) {
                    error_exit("UMI prefix can only have characters and numbers, but the given is: " + umi.prefix);
                }
            }
        }
        if(!umi.separator.empty()) {
            if(umi.separator.length()>10)
                error_exit("UMI separator cannot be longer than 10 base pairs");
            // validate bases
            for(int i=0; i<umi.separator.length(); i++) {
                char c = umi.separator[i];
                if(c!='A' && c!='T' && c!='C' && c!='G') {
                    error_exit("UMI separator can only have bases in {A, T, C, G}, but the given sequence is: " + umi.separator);
                }
            }
        }

    }

    if(overRepAnalysis.sampling < 1 || overRepAnalysis.sampling > 10000)
        error_exit("overrepresentation_sampling should be 1~10000");

    return true;
}

bool Options::shallDetectAdapter(bool isR2) {
    if(!adapter.enabled)
        return false;

    if(isR2) {
        return isPaired() && adapter.detectAdapterForPE && adapter.sequenceR2 == "auto";
    } else {
        if(isPaired())
            return adapter.detectAdapterForPE && adapter.sequence == "auto";
        else
            return adapter.sequence == "auto";
    }
}

void Options::initIndexFiltering(string blacklistFile1, string blacklistFile2, int threshold) {
    if(blacklistFile1.empty() && blacklistFile2.empty())
        return;

    if(!blacklistFile1.empty()){
        check_file_valid(blacklistFile1);
        indexFilter.blacklist1 = makeListFromFileByLine(blacklistFile1);
    }

    if(!blacklistFile2.empty()){
        check_file_valid(blacklistFile2);
        indexFilter.blacklist2 = makeListFromFileByLine(blacklistFile2);
    }

    if(indexFilter.blacklist1.empty() && indexFilter.blacklist2.empty())
        return;

    indexFilter.enabled = true;
    indexFilter.threshold = threshold;
}

vector<string> Options::makeListFromFileByLine(string filename) {
    vector<string> ret;
    ifstream file;
    file.open(filename.c_str(), ifstream::in);
    const int maxLine = 1000000;
    char line[maxLine];
    cerr << "filter by index, loading " << filename << endl;
    while(file.getline(line, maxLine)){
        // trim \n, \r or \r\n in the tail
        int readed = strlen(line);
        if(readed >=2 ){
            if(line[readed-1] == '\n' || line[readed-1] == '\r'){
                line[readed-1] = '\0';
                if(line[readed-2] == '\r')
                    line[readed-2] = '\0';
            }
        }
        string linestr(line);
        for(int i=0; i<linestr.length(); i++) {
            if(linestr[i] != 'A' && linestr[i] != 'T' && linestr[i] != 'C' && linestr[i] != 'G') {
                error_exit("processing " + filename + ", each line should be one barcode, which can only contain A/T/C/G");
            }
        }
        cerr << linestr << endl;
        ret.push_back(linestr);
    }
    cerr << endl;
    return ret;
}

string Options::getAdapter1(){
    if(adapter.sequence == "" || adapter.sequence == "auto")
        return "unspecified";
    else
        return adapter.sequence;
}

string Options::getAdapter2(){
    if(adapter.sequenceR2 == "" || adapter.sequenceR2 == "auto")
        return "unspecified";
    else
        return adapter.sequenceR2;
}

void Options::readDB() {
    //read gene ko species map;
    if(mHomoSearchOptions.genemap.empty())  error_exit("Gene KO GO species map file is empty : " + mHomoSearchOptions.genemap);
    if (verbose) {
        std::string msg = "Reading gene GO KO species map from file " + mHomoSearchOptions.genemap;
        longlog ? loginfolong(msg) : loginfo(msg);
    }
    std::set<std::string> KUSet;
    std::set<std::string> orgUSet;
    std::set<std::string> GOUSet;

    mHomoSearchOptions.filein.open(mHomoSearchOptions.genemap.c_str());
    if (!mHomoSearchOptions.filein.is_open()) error_exit("Can not open gene KO GO species map file : " + mHomoSearchOptions.genemap);
    const int maxLine = 10000;
    char line[maxLine];
    int readed = 0;
    vector<string> splVec;
    splVec.reserve(8);
    geneKoGoComb gkg;
    std::multimap<const uint32, const geneKoGoComb> fullidDbMMap;
    transSearch.orthIdSet.clear();
    std::map<const std::string, const uint32> tmpIdMap;
    while (mHomoSearchOptions.filein.getline(line, maxLine)) {
        readed = strlen(line);
        if (readed >= 2) {
            if (line[readed - 1] == '\n' || line[readed - 1] == '\r') {
                line[readed - 1] = '\0';
                if (line[readed - 2] == '\r') {
                    line[readed - 2] = '\0';
                }
            }
            string lineStr(line);
            splVec.clear();
            splitStr(lineStr, splVec, "\t");
            if (splVec.size() == 8) {
                uint32 id = (uint32) stoi(removeStr(splVec[1], "s2f_"));
                transSearch.orthIdSet.insert(id);
                tmpIdMap.insert(std::make_pair(splVec[0], id));
                auto ko = splVec[2];
                gkg.ko = ko;
                if (ko != "U") {
                    KUSet.insert(ko);
                }
                
                auto go = splVec[3];
                gkg.go = go;
                if (go != "U") {
                    GOUSet.insert(go);
                }
                
                gkg.symbol = splVec[4];
                gkg.gene = splVec[5];
                gkg.spec = splVec[6];
                orgUSet.insert(splVec[6]);
                gkg.coreOrthoPer = (float) std::stof(splVec[7]);
                fullidDbMMap.insert(std::make_pair(id, gkg));
            }
        }
    }

    if (transSearch.orthIdSet.empty()) error_exit("No ortholog id was detected in map file : " + mHomoSearchOptions.genemap);
    
    mHomoSearchOptions.filein.close();
    mHomoSearchOptions.filein.clear();
    transSearch.nKODB = KUSet.size();
    transSearch.nGODB = GOUSet.size();
    transSearch.nOrgsDB = orgUSet.size();
    transSearch.nIdDB = transSearch.orthIdSet.size();
    KUSet.clear();
    GOUSet.clear();
    orgUSet.clear();

    for (const auto & it : tmpIdMap) {
        auto itr = transSearch.orthIdSet.find(it.second);
        if (itr != transSearch.orthIdSet.end()) {
            mHomoSearchOptions.idDbMap.insert(std::make_pair(it.first, &(*itr)));
        }
    }
    tmpIdMap.clear();

    transSearch.coreOrthosDb = 0;
    for (const auto & it : transSearch.orthIdSet) {
        auto itr = fullidDbMMap.equal_range(it);
        gkg.ko = itr.first->second.ko;
        gkg.go = itr.first->second.go;
        gkg.symbol = itr.first->second.symbol;
        gkg.gene = itr.first->second.gene;
        gkg.coreOrthoPer = itr.first->second.coreOrthoPer;
        if(itr.first->second.coreOrthoPer >= 0.90) transSearch.coreOrthosDb++;
        mHomoSearchOptions.fullDbMap.insert(std::make_pair(&it, gkg));
    }
    fullidDbMMap.clear();
    
    //read ko full name file;
    mHomoSearchOptions.fileName.clear();
    mHomoSearchOptions.fileName = internalDBDir + "ko_fullname.txt";
    if(!mHomoSearchOptions.fileName.empty()) {
        if (verbose) {
            std::string msg = "Reading KO full name from file: " + mHomoSearchOptions.fileName;
            longlog ? loginfolong(msg) : loginfo(msg);
        }
        mHomoSearchOptions.filein.open(mHomoSearchOptions.fileName.c_str());
        if(!mHomoSearchOptions.filein.is_open()) error_exit("Can not open KO full name map file : " + mHomoSearchOptions.fileName);
        for(std::string s1, s2;
                std::getline(mHomoSearchOptions.filein, s1, '\t') && std::getline(mHomoSearchOptions.filein, s2);){
             mHomoSearchOptions.ko_fullname_map.insert(std::make_pair(s1, s2));
        }
        mHomoSearchOptions.filein.close();
        mHomoSearchOptions.filein.clear();
    } else {
        error_exit("KO full name file is empty: " + mHomoSearchOptions.fileName);
    }
    
    //read pathway ko map;
    mHomoSearchOptions.fileName.clear();
    mHomoSearchOptions.fileName = internalDBDir + "pathway_ko.txt";
    if (!mHomoSearchOptions.fileName.empty()) {
        if (verbose) {
            std::string msg = "Reading KO pathway map from file: " + mHomoSearchOptions.fileName;
            longlog ? loginfolong(msg) : loginfo(msg);
        }
        mHomoSearchOptions.filein.open(mHomoSearchOptions.fileName.c_str());
        if (!mHomoSearchOptions.filein.is_open()) error_exit("Can not open KO pathway map file : " + mHomoSearchOptions.fileName);
        std::unordered_set<std::string> pathwayUSet;
        std::string s1, s2;
        while (mHomoSearchOptions.filein >> s1 >> s2) {
            pathwayUSet.insert(s1);
            mHomoSearchOptions.pathway_ko_multimap.insert(std::make_pair(s1, s2));
        }
        transSearch.nPathwaysDB = pathwayUSet.size();
        pathwayUSet.clear();
        mHomoSearchOptions.filein.close();
        mHomoSearchOptions.filein.clear();
    } else {
        error_exit("KO pathway map file is empty: " + mHomoSearchOptions.fileName);
    }
    
    mHomoSearchOptions.fileName.clear();
    mHomoSearchOptions.fileName = internalDBDir + "pathway_ko_stats.txt";
    if (!mHomoSearchOptions.fileName.empty()) {
        if (verbose) {
            std::string msg = "Reading KO pathway stats map from file: " + mHomoSearchOptions.fileName;
            longlog ? loginfolong(msg) : loginfo(msg);
        }
        mHomoSearchOptions.filein.open(mHomoSearchOptions.fileName.c_str());
        if (!mHomoSearchOptions.filein.is_open()) error_exit("Can not open KO pathway map stats file : " + mHomoSearchOptions.fileName);
        std::string s1;
        int s2;
        while (mHomoSearchOptions.filein >> s1 >> s2) {
            mHomoSearchOptions.pathway_ko_stats_umap.insert(std::make_pair(s1, s2));
        }
        mHomoSearchOptions.filein.close();
        mHomoSearchOptions.filein.clear();
    } else {
        error_exit("KO pathway map stats file is empty: " + mHomoSearchOptions.fileName);
    }
    
}

void Options::parseSampleTable() {
    if (!mHomoSearchOptions.sampleTable.empty()) {
        check_file_valid(mHomoSearchOptions.sampleTable);
        if (verbose) {
            std::string msg = "Reading sample table from file " + mHomoSearchOptions.sampleTable;
            longlog ? loginfolong(msg) : loginfo(msg);
        }
    } else {
        error_exit("sample table file should be specified by --sampletable");
    }
    
    ifstream file;
    file.open(mHomoSearchOptions.sampleTable.c_str(), ifstream::in);
    const int maxLine = 1000;
    char line[maxLine];
    while(file.getline(line, maxLine)){
        // trim \n, \r or \r\n in the tail
        int readed = strlen(line);
        if(readed >=2 ){
            if(line[readed-1] == '\n' || line[readed-1] == '\r'){
                line[readed-1] = '\0';
                if(line[readed-2] == '\r')
                    line[readed-2] = '\0';
            }
        }
        string linestr(line);

        // comment or header line

        vector<string> splitted;
        splitStr(linestr, splitted, "\t");
        // a valid line need 4 columns: name, left, center, right

        Sample s;
        s.prefix = trimStr(splitted[0]);
        s.path = dirname(s.prefix);
        s.in1 = trimStr(splitted[1]);
        if (splitted.size() == 3) {
            s.feature = trimStr(splitted[2]);
        } else if (splitted.size() == 4) {
            s.in2 = trimStr(splitted[2]);
            s.feature = trimStr(splitted[3]);
        } else {
            error_exit("sample table must be 3 columns with sample name, forward reads and sample feature or "
                    "4 columns with sample name, forward reads, reverse reads and sample feature");
        }

        this->samples.push_back(s);
    }
}

void Options::mkSelectedPathwayDB(){
//    mHomoSearchOptions.filein.open(mHomoSearchOptions.pathway.c_str());
//    if(!mHomoSearchOptions.filein.is_open()) error_exit("Can not open file: " + mHomoSearchOptions.pathway);
//    std::vector<std::string> pathwayVec;
//    std::string pathwayId;
//    while(mHomoSearchOptions.filein >> pathwayId){
//        pathwayVec.push_back(pathwayId);
//    }
//    mHomoSearchOptions.filein.close();
//    mHomoSearchOptions.filein.clear();
//    
//    if(verbose){
//        std::string msg = to_string(pathwayVec.size()) + " pathways have been provided";
//        loginfo(msg);
//    }
//    
//    //extract pathways and kos;
//    std::multimap<std::string, std::string> selected_pathway_ko_multimap;
//    std::set<std::string> uniqKOSet;
//    for(auto & it : pathwayVec){
//        auto pathwayKO = mHomoSearchOptions.pathway_ko_multimap.equal_range(it); 
//        for(auto itt = pathwayKO.first; itt != pathwayKO.second; ++ itt){
//            uniqKOSet.insert(itt->second);
//            selected_pathway_ko_multimap.insert(std::make_pair(itt->first, itt->second));
//        }
//    }
//    
//    pathwayVec.clear();
//    mHomoSearchOptions.pathway_ko_multimap.clear();
//    mHomoSearchOptions.pathway_ko_multimap.insert(selected_pathway_ko_multimap.begin(), selected_pathway_ko_multimap.end());
//    selected_pathway_ko_multimap.clear();
//    
//    if(verbose){
//        std::string msg = "pathway KO map size is " + to_string(mHomoSearchOptions.pathway_ko_multimap.size());
//        loginfo(msg);
//    }
//    
//    //extract protein IDs;
//    //swap key values
//    std::multimap<std::string, std::string> KOProteinMMap;
//    for(auto & it : mHomoSearchOptions.db_map){
//        KOProteinMMap.insert(std::make_pair(it.second, it.first));
//    }
//    mHomoSearchOptions.db_map.clear();
//    
//    //update new protein ko map;
//    std::vector<std::string> selectedProteinIDVec;
//    for(auto & it : uniqKOSet){
//        auto KOProtein = KOProteinMMap.equal_range(it);
//        for(auto & itr = KOProtein.first; itr != KOProtein.second; ++ itr){
//            selectedProteinIDVec.push_back(itr->second);
//            mHomoSearchOptions.db_map.insert(std::make_pair(itr->second, itr->first));
//        }
//    }
//    
//    KOProteinMMap.clear();
//
//    if (verbose) {
//        std::string msg = "Protein KO map size is " + to_string(mHomoSearchOptions.db_map.size());
//        loginfo(msg);
//    }
//    
//    //read protein.fasta file;
//    FastaReader pathwayProteinFasta(mHomoSearchOptions.genefa);
//    pathwayProteinFasta.readAll();
//    
//    auto proteinSeqMap = pathwayProteinFasta.contigs();
//
//    if (mkdir("selected_pathway_database", 0777) == -1) error_exit("Can not create directory selected_pathway_database or you must remove it first");
//    std::string foutProtein = "selected_pathway_database/selected_pathway_protein_aas.pep.fasta";
//    
//    std::ofstream * fout = new std::ofstream();
//    fout->open(foutProtein.c_str(), std::ofstream::out);
//    if(!fout->is_open()) error_exit("Can not open file: " + foutProtein);
//    for(auto & it : selectedProteinIDVec){
//        auto itr = proteinSeqMap.find(it);
//        if(itr != proteinSeqMap.end()){
//            *fout << ">" << itr->first << "\n" << itr->second << "\n";
//        }
//    }
//    fout->close();
//    if(fout) delete fout;
//    selectedProteinIDVec.clear();
//
//    std::string bwtfm_cmd = seq2funDir + "/bin/mkbwt";
//    std::string fmi_cmd = seq2funDir + "/bin/mkfmi";
//    std::string p_name = "selected_pathway_database/selected_pathway_protein_aas.pep.fasta";
//    std::string f_name = "selected_pathway_database/proteins";
//
//    std::stringstream comandss;
//    comandss << bwtfm_cmd << " -n " << thread << " -a ACDEFGHIKLMNPQRSTVWY " << "-o " << f_name << " " << p_name << "\n";
//    std::cout << comandss.str() << std::endl;
//    system(comandss.str().c_str());
//
//    comandss.str();
//
//    comandss << fmi_cmd << " " << f_name << "\n";
//    system(comandss.str().c_str());
//    comandss.str();
//    transSearch.tfmi = f_name + ".fmi";
}

int Options::getWorkingSampleId(string & samplePrefix){
    int j = 0;
    for(int i = 0; i < samples.size(); i++){
        if(samples.at(i).prefix == samplePrefix){
            j = i;
        }
    }
    return j;
}

void Options::readSampleExtraction(){
    if(mSeqExtractions.sampleMappedTableStr.empty()){
        error_exit("sampleMappedTable file should be specified by --sampleMappedTable");
    } else {
        check_file_valid(mSeqExtractions.sampleMappedTableStr);
        if (verbose) {
            std::string msg = "Reading sampleMappedTable for seq extraction from file " + mSeqExtractions.sampleMappedTableStr;
            longlog ? loginfolong(msg) : loginfo(msg);
        }
    }
    
    ifstream file;
    file.open(mSeqExtractions.sampleMappedTableStr.c_str(), ifstream::in);
    if(!file.is_open()) error_exit("can not open sampleMappedTable, please check it!");
    const int maxLine = 1000;
    char line[maxLine];
    
    mSeqExtractions.samplesVecF.clear();
    mSeqExtractions.samplesVecR.clear();
    std::set<std::string> samUset;
    std::vector<std::string> tmpVec;
    std::string strF, strR;
    while(file.getline(line, maxLine)){
        string lineStr = trimEnd(line);
        tmpVec.clear();
        tmpVec = split2(lineStr, '\t');
        strF = "";
        strR = "";
        if(tmpVec.size() == 1){
            strF = tmpVec.front();
            if(strF.length() > 0){
                if(samUset.size() > 0 && samUset.find(strF) != samUset.end()){
                    error_exit("You have duplicated values in sample file, please remove the duplicated value: " + strF);
                }
                samUset.insert(strF);
                mSeqExtractions.samplesVecF.emplace_back(strF);
                mSeqExtractions.paired = false;
            } else {
                error_exit("You have incorrect sample names!");
            }
        } else if(tmpVec.size() == 2){
            strF = tmpVec.front();
            if (strF.length() > 0) {
                if (samUset.size() > 0 && samUset.find(strF) != samUset.end()) {
                    error_exit("You have duplicated values in sample file, please remove the duplicated value: " + strF);
                }
                samUset.insert(strF);
                mSeqExtractions.samplesVecF.emplace_back(strF);
            } else {
                error_exit("You have incorrect sample names!");
            }
            
            strR = tmpVec.back();
            if (strR.length() > 0) {
                if (samUset.size() > 0 && samUset.find(strR) != samUset.end()) {
                    error_exit("You have duplicated values in sample file, please remove the duplicated value: " + strR);
                }
                samUset.insert(strR);
                mSeqExtractions.samplesVecR.emplace_back(strR);
                mSeqExtractions.paired = true;
            } else {
                error_exit("You have incorrect sample names!");
            }
        } else {
            error_exit("Your sample table is not correct: " + mSeqExtractions.sampleMappedTableStr);
        }
    }
    file.close();
    
    if(mSeqExtractions.samplesVecF.empty()){
        error_exit("Your sample table is not correct: " + mSeqExtractions.sampleMappedTableStr);
    }
    
    if(mSeqExtractions.paired && mSeqExtractions.samplesVecR.empty()){
        error_exit("Your sample table is not correct: " + mSeqExtractions.sampleMappedTableStr);
    }
    
    if(mSeqExtractions.targetGeneTableStr.empty()){
        error_exit("targetGeneTable file must be specified by --targetGeneTable");
    } else {
        check_file_valid(mSeqExtractions.targetGeneTableStr);
        if(verbose){
            std::string msg = "Reading targetGeneTable for seq extraction from file " + mSeqExtractions.targetGeneTableStr;
            longlog ? loginfolong(msg) : loginfo(msg); 
        }
    }
    
    
    {ifstream file;
    file.open(mSeqExtractions.targetGeneTableStr.c_str(), ifstream::in);
    if(!file.is_open()) error_exit("can not open targetGeneTable, please check it!");
    if (verbose) {
        std::string msg = "Reading gene table for seq extraction from file " + mSeqExtractions.targetGeneTableStr;
        longlog ? loginfolong(msg) : loginfo(msg);
    }
    mSeqExtractions.targetGenesVec.clear();
    std::set<std::string> geneSet;
    int i = 0;
    while(file.getline(line, maxLine)){
        string lineStr = trimEnd(line);
        if(lineStr.length() > 0){
            if(!geneSet.empty() && geneSet.find(lineStr) == geneSet.end()){
                error_exit("You have duplicated values gene table, please remove the duplicated value: " + lineStr);
            }
            mSeqExtractions.targetGenesVec.emplace_back(lineStr);
            if(starts_with(lineStr, "s2f_")){
                s2fid4Strct = true;
            }
//            } else if(starts_with(lineStr, "K")){
//                s2fid4Strct = false;
//            } else {
//                s2fid4Strct = true;
//            }
        }
    }
    file.close();
    }
    
    if (!file_exists(mSeqExtractions.outputDir)) {
        mkdir(mSeqExtractions.outputDir.c_str(), 0777);
    }

    if (file_exists(mSeqExtractions.outputDir) && !is_directory(mSeqExtractions.outputDir)) {
        error_exit(mSeqExtractions.outputDir + " is a file, not a directory");
    }

    if (!file_exists(mSeqExtractions.outputDir) || !is_directory(mSeqExtractions.outputDir)) {
        error_exit(mSeqExtractions.outputDir + " is not a directory, or cannot be created");
    }
}
