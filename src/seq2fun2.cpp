#include <stdio.h>
#include <time.h>
#include <sstream>
#include <map>
#include <fstream>
#include <cmath>
#include <thread>
#include <vector>
#include <unordered_map>
#include <utility> 

#include "fastqreader.h"
#include "unittest.h"
#include "cmdline.h"
#include "util.h"
#include "options.h"
#include "processor.h"
#include "evaluator.h"
#include "bwtfmiDB.h"
#include "htmlreporterall.h"

string command;
mutex logmtx;

int main(int argc, char* argv[]) {
    // display version info if no argument is given
    time_t t_begin = time(NULL);

    if (argc == 1) {
        cerr << "Seq2Fun2: high-throughput functional profiling of RNA-seq data for non-model organisms" << endl << "version " << SEQ2FUNR_VER << endl;
    }
    if (argc == 2 && strcmp(argv[1], "test") == 0) {
        UnitTest tester;
        tester.run();
        return 0;
    }
    if (argc == 2 && (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0)) {
        cerr << "Seq2Fun " << SEQ2FUNR_VER << endl;
        return 0;
    }

    cmdline::parser cmd;
    // input/output
    cmd.add<string>("sampletable", 's', "(recommended) sample table must consist of 3 columns (sample prefix name (sample01), forward reads name (sample01_R1.fq.gz), group info (control) for single-reads or 4 columns (sample prefix name (sample01), forward reads (sample01_R1.fq.gz), reverse reads (sample01_R2.fq.gz), group info (control) for paired-end reads. The columns must be separated by tab", false, "");
    cmd.add<string>("in1", 'i', "read1 input file name", false, "");
    cmd.add<string>("in2", 'I', "read2 input file name", false, "");
    cmd.add<string>("prefix", 'X', "prefix name for output files, eg: sample01", false, "");
    cmd.add("outputMappedCleanReads", 0, "enable output mapped clean reads into fastq.gz files, by default is false, using --outputMappedCleanReads to enable it");
    cmd.add("outputReadsKOMap", 0, "enable output mapped clean reads-KO map into .gz files, by default is false, using --outputReadsKOMap to enable it");

    // Homology search;
    cmd.add<string>("genemap", 'D', "gene/protein KO species map", false, "");
    cmd.add("profiling", 0, "profiling mode, by default is false, using --profiling to enable it");

    // translated search
    cmd.add<string>("tfmi", 'd', "fmi index of Protein database", false, "");
    cmd.add<string>("mode", 'K', "searching mode either tGREEDY or tMEM (maximum exactly match). By default greedy", false, "tGREEDY");
    cmd.add<int>("mismatch", 'E', "number of mismatched amino acid in sequence comparison with protein database with default value 2", false, 2);
    cmd.add<int>("minscore", 'j', "minimum matching score of amino acid sequence in comparison with protein database with default value 80", false, 80);
    cmd.add<int>("minlength", 'J', "minimum matching length of amino acid sequence in comparison with protein database with default value 19, for GREEDY and 13 for MEM model", false, 0);
    cmd.add<int>("maxtranslength", 'm', "maximum cutoff of translated peptides, it must be no less than minlength, with default 60", false, 60);
    cmd.add("allFragments", 0, "enable this function will force Seq2Fun to use all the translated AA fragments with length > minlength. This will slightly help to classify reads contain the true stop codon and start codon; This could have limited impact on the accuracy for comparative study and enable this function will slow down the Seq2Fun. by default is false, using --allFragments to enable it");
    cmd.add<string>("codontable", 0, "select the codon table (same as blastx in NCBI), we provide 20 codon tables from 'https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG31'. By default is the codontable1 (Standard Code)", false, "codontable1");
    cmd.add<string>("dbDir", 0, "dir for internal database such as ko_fullname.txt", false, "");
    
    //selected pathways
    cmd.add<string>("pathway", 'Z', "list of selected pathways for target pathways analysis", false, "");
    cmd.add<string>("genefa", 'z', "the gene/protein sequences fasta file for retrieving proteins in selected pathways to construct database", false, "");

    // threading
    cmd.add<int>("thread", 'w', "worker thread number, default is 2", false, 2);
    cmd.add("verbose", 'V', "enable verbose");
    cmd.add("debug", 0, "enable debug");
    
    cmd.add<int>("reads_buffer", 0, "specify reads buffer size (MB) for each file.", false, 1);
    cmd.add("fix_mgi_id", 0, "the MGI FASTQ ID format is not compatible with many BAM operation tools, enable this option to fix it.");

    cmd.add("phred64", '6', "indicate the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)");
    cmd.add<int>("reads_to_process", 0, "specify how many reads/pairs to be processed. Default 0 means process all reads.", false, 0);

    // adapter
    cmd.add("disable_adapter_trimming", 'A', "adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled");
    cmd.add<string>("adapter_sequence", 'a', "the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped.", false, "auto");
    cmd.add<string>("adapter_sequence_r2", 0, "the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence>", false, "auto");
    cmd.add<string>("adapter_fasta", 0, "specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file", false, "");
    cmd.add("detect_adapter_for_pe", 0, "by default, the auto-detection for adapter is for SE data input only, turn on this option to enable it for PE data.");
    //polyA tail
    cmd.add("no_trim_polyA", 0, "by default, ployA tail will be trimmed. If this option is specified, polyA trimming is disabled");

    // trimming
    cmd.add<int>("trim_front1", 'f', "trimming how many bases in front for read1, default is 0", false, 0);
    cmd.add<int>("trim_tail1", 't', "trimming how many bases in tail for read1, default is 0", false, 0);
    cmd.add<int>("max_len1", 'b', "if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation", false, 0);
    cmd.add<int>("trim_front2", 'F', "trimming how many bases in front for read2. If it's not specified, it will follow read1's settings", false, 0);
    cmd.add<int>("trim_tail2", 'T', "trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings", false, 0);
    cmd.add<int>("max_len2", 'B', "if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2. Default 0 means no limitation. If it's not specified, it will follow read1's settings", false, 0);

    // polyG tail trimming
    cmd.add("trim_poly_g", 'g', "force polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data");
    cmd.add<int>("poly_g_min_len", 0, "the minimum length to detect polyG in the read tail. 10 by default.", false, 10);
    cmd.add("disable_trim_poly_g", 'G', "disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data");

    // polyX tail trimming
    cmd.add("trim_poly_x", 'x', "enable polyX trimming in 3' ends.");
    cmd.add<int>("poly_x_min_len", 0, "the minimum length to detect polyX in the read tail. 10 by default.", false, 10);

    // cutting by quality
    cmd.add("cut_front", '5', "move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.");
    cmd.add("cut_tail", '3', "move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.");
    cmd.add("cut_right", 'r', "move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.");
    cmd.add<int>("cut_window_size", 'W', "the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4", false, 4);
    cmd.add<int>("cut_mean_quality", 'M', "the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20)", false, 20);
    cmd.add<int>("cut_front_window_size", 0, "the window size option of cut_front, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_front_mean_quality", 0, "the mean quality requirement option for cut_front, default to cut_mean_quality if not specified", false, 20);
    cmd.add<int>("cut_tail_window_size", 0, "the window size option of cut_tail, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_tail_mean_quality", 0, "the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified", false, 20);
    cmd.add<int>("cut_right_window_size", 0, "the window size option of cut_right, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_right_mean_quality", 0, "the mean quality requirement option for cut_right, default to cut_mean_quality if not specified", false, 20);

    // quality filtering
    cmd.add("disable_quality_filtering", 'Q', "quality filtering is enabled by default. If this option is specified, quality filtering is disabled");
    cmd.add<int>("qualified_quality_phred", 'q', "the quality value that a base is qualified. Default 20 means phred quality >=Q15 is qualified.", false, 20);
    cmd.add<int>("unqualified_percent_limit", 'u', "how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%", false, 40);
    cmd.add<int>("n_base_limit", 'n', "if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5", false, 5);
    cmd.add<int>("average_qual", 'e', "if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement", false, 0);

    // length filtering
    cmd.add("disable_length_filtering", 'L', "length filtering is enabled by default. If this option is specified, length filtering is disabled");
    cmd.add<int>("length_required", 'l', "reads shorter than length_required will be discarded, default is 60.", false, 60);
    cmd.add<int>("length_limit", 0, "reads longer than length_limit will be discarded, default 0 means no limitation.", false, 0);

    // low complexity filtering
    cmd.add("no_low_complexity_filter", 0, "disable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).");
    cmd.add<int>("complexity_threshold", 'Y', "the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required.", false, 30);

    // filter by indexes
    cmd.add<string>("filter_by_index1", 0, "specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line", false, "");
    cmd.add<string>("filter_by_index2", 0, "specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line", false, "");
    cmd.add<int>("filter_by_index_threshold", 0, "the allowed difference of index barcode for index filtering, default 0 means completely identical.", false, 0);

    // base correction in overlapped regions of paired end data
    cmd.add("disable_correction", 'c', "disenable base correction in overlapped regions (only for PE data), default is enabled");
    cmd.add<int>("overlap_len_require", 'v', "the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default.", false, 30);
    cmd.add<int>("overlap_diff_limit", 0, "the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default.", false, 5);
    cmd.add<int>("overlap_diff_percent_limit", 0, "the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%.", false, 20);

    // umi
    cmd.add("umi", 'U', "enable unique molecular identifier (UMI) preprocessing");
    cmd.add<string>("umi_loc", 0, "specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none", false, "");
    cmd.add<int>("umi_len", 0, "if the UMI is in read1/read2, its length should be provided", false, 0);
    cmd.add<string>("umi_prefix", 0, "if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default", false, "");
    cmd.add<int>("umi_skip", 0, "if the UMI is in read1/read2, Seq2Fun can skip several bases following UMI, default is 0", false, 0);

    // overrepresented sequence analysis
    cmd.add("overrepresentation_analysis", 'p', "enable overrepresented sequence analysis.");
    cmd.add<int>("overrepresentation_sampling", 'P', "one in (--overrepresentation_sampling) reads will be computed for overrepresentation analysis (1~10000), smaller is slower, default is 20.", false, 20);


    // deprecated options
    cmd.add("cut_by_quality5", 0, "DEPRECATED, use --cut_front instead.");
    cmd.add("cut_by_quality3", 0, "DEPRECATED, use --cut_tail instead.");
    cmd.add("cut_by_quality_aggressive", 0, "DEPRECATED, use --cut_right instead.");
    cmd.add("discard_unmerged", 0, "DEPRECATED, no effect now, see the introduction for merging.");

    cmd.parse_check(argc, argv);

    if (argc == 1) {
        cerr << cmd.usage() << endl;
        return 0;
    }

    Options opt;

    opt.seq2funProgPath = string(argv[0]);
    opt.seq2funDir = removeStr(opt.seq2funProgPath, "bin/seq2fun2");
    opt.internalDBDir = cmd.get<string>("dbDir") == "" ? opt.seq2funDir + "database" : cmd.get<string>("dbDir");
    opt.internalDBDir = checkDirEnd(opt.internalDBDir);

    // threading
    opt.thread = cmd.get<int>("thread");
    int n_t = std::thread::hardware_concurrency();
    opt.thread = std::min(std::min(opt.thread, 32), n_t);

    opt.compression = 4;
    opt.readsToProcess = cmd.get<int>("reads_to_process");
    if(cmd.get<int>("reads_buffer") < 1) {
        error_exit("reads_buffer should be greater or equal to 1MB.");
    }
    opt.phred64 = cmd.exist("phred64");
    opt.verbose = cmd.exist("verbose");
    opt.debug = cmd.exist("debug");
    opt.fixMGI = cmd.exist("fix_mgi_id");

    // adapter cutting
    opt.adapter.enabled = !cmd.exist("disable_adapter_trimming");
    opt.adapter.detectAdapterForPE = cmd.exist("detect_adapter_for_pe");
    opt.adapter.sequence = cmd.get<string>("adapter_sequence");
    opt.adapter.sequenceR2 = cmd.get<string>("adapter_sequence_r2");
    opt.adapter.fastaFile = cmd.get<string>("adapter_fasta");
    if (opt.adapter.sequenceR2 == "auto" && !opt.adapter.detectAdapterForPE && opt.adapter.sequence != "auto") {
        opt.adapter.sequenceR2 = opt.adapter.sequence;
    }
    if (!opt.adapter.fastaFile.empty()) {
        opt.loadFastaAdapters();
    }
    //polyA tail trimming
    opt.adapter.polyA = !cmd.exist("no_trim_polyA");

    // trimming
    opt.trim.front1 = cmd.get<int>("trim_front1");
    opt.trim.tail1 = cmd.get<int>("trim_tail1");
    opt.trim.maxLen1 = cmd.get<int>("max_len1");
    // read2 settings follows read1 if it's not specified
    if (cmd.exist("trim_front2"))
        opt.trim.front2 = cmd.get<int>("trim_front2");
    else
        opt.trim.front2 = opt.trim.front1;
    if (cmd.exist("trim_tail2"))
        opt.trim.tail2 = cmd.get<int>("trim_tail2");
    else
        opt.trim.tail2 = opt.trim.tail1;
    if (cmd.exist("max_len2"))
        opt.trim.maxLen2 = cmd.get<int>("max_len2");
    else
        opt.trim.maxLen2 = opt.trim.maxLen1;

    // polyG tail trimming
    if (cmd.exist("trim_poly_g") && cmd.exist("disable_trim_poly_g")) {
        error_exit("You cannot enabled both trim_poly_g and disable_trim_poly_g");
    } else if (cmd.exist("trim_poly_g")) {
        opt.polyGTrim.enabled = true;
    } else if (cmd.exist("disable_trim_poly_g")) {
        opt.polyGTrim.enabled = false;
    }
    opt.polyGTrim.minLen = cmd.get<int>("poly_g_min_len");

    // polyX tail trimming
    if (cmd.exist("trim_poly_x")) {
        opt.polyXTrim.enabled = true;
    }
    opt.polyXTrim.minLen = cmd.get<int>("poly_x_min_len");


    // sliding window cutting by quality
    opt.qualityCut.enabledFront = cmd.exist("cut_front");
    // back compatible with old versions
    if (!opt.qualityCut.enabledFront) {
        opt.qualityCut.enabledFront = cmd.exist("cut_by_quality5");
        if (opt.qualityCut.enabledFront)
            cerr << "WARNING: cut_by_quality5 is deprecated, please use cut_front instead." << endl;
    }
    opt.qualityCut.enabledTail = cmd.exist("cut_tail");
    // back compatible with old versions
    if (!opt.qualityCut.enabledFront) {
        opt.qualityCut.enabledFront = cmd.exist("cut_by_quality3");
        if (opt.qualityCut.enabledFront)
            cerr << "WARNING: cut_by_quality3 is deprecated, please use cut_tail instead." << endl;
    }
    opt.qualityCut.enabledRight = cmd.exist("cut_right");
    // back compatible with old versions
    if (!opt.qualityCut.enabledRight) {
        opt.qualityCut.enabledRight = cmd.exist("cut_by_quality_aggressive");
        if (opt.qualityCut.enabledRight)
            cerr << "WARNING: cut_by_quality_aggressive is deprecated, please use cut_right instead." << endl;
    }

    opt.qualityCut.windowSizeShared = cmd.get<int>("cut_window_size");
    opt.qualityCut.qualityShared = cmd.get<int>("cut_mean_quality");

    if (cmd.exist("cut_front_window_size"))
        opt.qualityCut.windowSizeFront = cmd.get<int>("cut_front_window_size");
    else
        opt.qualityCut.windowSizeFront = opt.qualityCut.windowSizeShared;
    if (cmd.exist("cut_front_mean_quality"))
        opt.qualityCut.qualityFront = cmd.get<int>("cut_front_mean_quality");
    else
        opt.qualityCut.qualityFront = opt.qualityCut.qualityShared;

    if (cmd.exist("cut_tail_window_size"))
        opt.qualityCut.windowSizeTail = cmd.get<int>("cut_tail_window_size");
    else
        opt.qualityCut.windowSizeTail = opt.qualityCut.windowSizeShared;
    if (cmd.exist("cut_tail_mean_quality"))
        opt.qualityCut.qualityTail = cmd.get<int>("cut_tail_mean_quality");
    else
        opt.qualityCut.qualityTail = opt.qualityCut.qualityShared;

    if (cmd.exist("cut_right_window_size"))
        opt.qualityCut.windowSizeRight = cmd.get<int>("cut_right_window_size");
    else
        opt.qualityCut.windowSizeRight = opt.qualityCut.windowSizeShared;
    if (cmd.exist("cut_right_mean_quality"))
        opt.qualityCut.qualityRight = cmd.get<int>("cut_right_mean_quality");
    else
        opt.qualityCut.qualityRight = opt.qualityCut.qualityShared;

    // raise a warning if cutting option is not enabled but -W/-M is enabled
    if (!opt.qualityCut.enabledFront && !opt.qualityCut.enabledTail && !opt.qualityCut.enabledRight) {
        if (cmd.exist("cut_window_size") || cmd.exist("cut_mean_quality")
                || cmd.exist("cut_front_window_size") || cmd.exist("cut_front_mean_quality")
                || cmd.exist("cut_tail_window_size") || cmd.exist("cut_tail_mean_quality")
                || cmd.exist("cut_right_window_size") || cmd.exist("cut_right_mean_quality"))
            cerr << "WARNING: you specified the options for cutting by quality, but forogt to enable any of cut_front/cut_tail/cut_right. This will have no effect." << endl;
    }

    // quality filtering
    opt.qualfilter.enabled = !cmd.exist("disable_quality_filtering");
    opt.qualfilter.qualifiedQual = num2qual(cmd.get<int>("qualified_quality_phred"));
    opt.qualfilter.unqualifiedPercentLimit = cmd.get<int>("unqualified_percent_limit");
    opt.qualfilter.avgQualReq = cmd.get<int>("average_qual");
    opt.qualfilter.nBaseLimit = cmd.get<int>("n_base_limit");

    // length filtering
    opt.lengthFilter.enabled = !cmd.exist("disable_length_filtering");
    opt.lengthFilter.requiredLength = cmd.get<int>("length_required");

    if (cmd.get<int>("minlength") == 0) {
        if (opt.transSearch.mode == tGREEDY) {
            opt.transSearch.minAAFragLength = 19;
        } else {
            opt.transSearch.minAAFragLength = 13;
        }
    } else {
        opt.transSearch.minAAFragLength = cmd.get<int>("minlength");
    }
    
    opt.lengthFilter.requiredLength = max(opt.lengthFilter.requiredLength, static_cast<int> (opt.transSearch.minAAFragLength) * 3);
    opt.lengthFilter.maxLength = cmd.get<int>("length_limit");

    // low complexity filter
    opt.complexityFilter.enabled = !cmd.exist("no_low_complexity_filter");
    opt.complexityFilter.threshold = (min(100, max(0, cmd.get<int>("complexity_threshold")))) / 100.0;

    // overlap correction
    opt.correction.enabled = !cmd.exist("disable_correction");
    opt.overlapRequire = cmd.get<int>("overlap_len_require");
    opt.overlapDiffLimit = cmd.get<int>("overlap_diff_limit");
    opt.overlapDiffPercentLimit = cmd.get<int>("overlap_diff_percent_limit");

    // umi
    opt.umi.enabled = cmd.exist("umi");
    opt.umi.length = cmd.get<int>("umi_len");
    opt.umi.prefix = cmd.get<string>("umi_prefix");
    opt.umi.skip = cmd.get<int>("umi_skip");
    if (opt.umi.enabled) {
        string umiLoc = cmd.get<string>("umi_loc");
        str2lower(umiLoc);
        if (umiLoc.empty())
            error_exit("You've enabled UMI by (--umi), you should specify the UMI location by (--umi_loc)");
        if (umiLoc != "index1" && umiLoc != "index2" && umiLoc != "read1" && umiLoc != "read2" && umiLoc != "per_index" && umiLoc != "per_read") {
            error_exit("UMI location can only be index1/index2/read1/read2/per_index/per_read");
        }
        if (!opt.isPaired() && (umiLoc == "index2" || umiLoc == "read2"))
            error_exit("You specified the UMI location as " + umiLoc + ", but the input data is not paired end.");
        if (opt.umi.length == 0 && (umiLoc == "read1" || umiLoc == "read2" || umiLoc == "per_read"))
            error_exit("You specified the UMI location as " + umiLoc + ", but the length is not specified (--umi_len).");
        if (umiLoc == "index1") {
            opt.umi.location = UMI_LOC_INDEX1;
        } else if (umiLoc == "index2") {
            opt.umi.location = UMI_LOC_INDEX2;
        } else if (umiLoc == "read1") {
            opt.umi.location = UMI_LOC_READ1;
        } else if (umiLoc == "read2") {
            opt.umi.location = UMI_LOC_READ2;
        } else if (umiLoc == "per_index") {
            opt.umi.location = UMI_LOC_PER_INDEX;
        } else if (umiLoc == "per_read") {
            opt.umi.location = UMI_LOC_PER_READ;
        }
    }

    // overrepresented sequence analysis
    opt.overRepAnalysis.enabled = cmd.exist("overrepresentation_analysis");
    opt.overRepAnalysis.sampling = cmd.get<int>("overrepresentation_sampling");

    // filtering by index
    string blacklist1 = cmd.get<string>("filter_by_index1");
    string blacklist2 = cmd.get<string>("filter_by_index2");
    int indexFilterThreshold = cmd.get<int>("filter_by_index_threshold");
    opt.initIndexFiltering(blacklist1, blacklist2, indexFilterThreshold);

    //homology search
    opt.mHomoSearchOptions.genemap = cmd.get<string>("genemap");
    opt.mHomoSearchOptions.profiling = cmd.exist("profiling");

    //translated search
    opt.transSearch.tmode = cmd.get<string>("mode");
    if (opt.transSearch.tmode == "tGREEDY") {
        opt.transSearch.mode = tGREEDY;
    } else if (opt.transSearch.tmode == "tMEM") {
        opt.transSearch.mode = tMEM;
    } else {
        error_exit("you must be use either tGREEDY or tMEM mode");
    }

    opt.transSearch.tCodonTable = cmd.get<string>("codontable");
    if (opt.transSearch.tCodonTable == "codontable1") {
        opt.transSearch.codonTable = codontable1;
    } else if (opt.transSearch.tCodonTable == "codontable2") {
        opt.transSearch.codonTable = codontable2;
    } else if (opt.transSearch.tCodonTable == "codontable3") {
        opt.transSearch.codonTable = codontable3;
    } else if (opt.transSearch.tCodonTable == "codontable4") {
        opt.transSearch.codonTable = codontable4;
    } else if (opt.transSearch.tCodonTable == "codontable5") {
        opt.transSearch.codonTable = codontable5;
    } else if (opt.transSearch.tCodonTable == "codontable6") {
        opt.transSearch.codonTable = codontable6;
    } else if (opt.transSearch.tCodonTable == "codontable9") {
        opt.transSearch.codonTable = codontable9;
    } else if (opt.transSearch.tCodonTable == "codontable10") {
        opt.transSearch.codonTable = codontable10;
    } else if (opt.transSearch.tCodonTable == "codontable12") {
        opt.transSearch.codonTable = codontable12;
    } else if (opt.transSearch.tCodonTable == "codontable13") {
        opt.transSearch.codonTable = codontable13;
    } else if (opt.transSearch.tCodonTable == "codontable14") {
        opt.transSearch.codonTable = codontable14;
    } else if (opt.transSearch.tCodonTable == "codontable16") {
        opt.transSearch.codonTable = codontable16;
    } else if (opt.transSearch.tCodonTable == "codontable26") {
        opt.transSearch.codonTable = codontable26;
    } else if (opt.transSearch.tCodonTable == "codontable21") {
        opt.transSearch.codonTable = codontable21;
    } else if (opt.transSearch.tCodonTable == "codontable22") {
        opt.transSearch.codonTable = codontable22;
    } else if (opt.transSearch.tCodonTable == "codontable24") {
        opt.transSearch.codonTable = codontable24;
    } else if (opt.transSearch.tCodonTable == "codontable27") {
        opt.transSearch.codonTable = codontable27;
    } else if (opt.transSearch.tCodonTable == "codontable29") {
        opt.transSearch.codonTable = codontable29;
    } else if (opt.transSearch.tCodonTable == "codontable30") {
        opt.transSearch.codonTable = codontable30;
    } else if (opt.transSearch.tCodonTable == "codontable31") {
        opt.transSearch.codonTable = codontable31;
    } else if (opt.transSearch.tCodonTable == "codontable33") {
        opt.transSearch.codonTable = codontable33;
    } else {
        error_exit("you must select one codon table");
    }

    if (opt.verbose) {
        std::cout << "Codon table of " << opt.transSearch.tCodonTable << " is selected" << std::endl;
    }

    opt.transSearch.misMatches = cmd.get<int>("mismatch");
    opt.transSearch.minScore = cmd.get<int>("minscore");

    opt.transSearch.maxTransLength = cmd.get<int>("maxtranslength");
    opt.transSearch.maxTransLength = max(opt.transSearch.maxTransLength, opt.transSearch.minAAFragLength);
    opt.transSearch.maxTransLength = min((unsigned) 60, opt.transSearch.maxTransLength);
    opt.transSearch.allFragments = cmd.exist("allFragments");
    opt.transSearch.tfmi = cmd.get<string>("tfmi");

    //read all database tables, maps;
    opt.readDB();

    // for selected pathway and customized database
    opt.mHomoSearchOptions.pathway = cmd.get<string>("pathway");
    opt.mHomoSearchOptions.genefa = cmd.get<string>("genefa");

    if (!opt.mHomoSearchOptions.pathway.empty() && !opt.mHomoSearchOptions.genefa.empty()) {
        if (!opt.transSearch.tfmi.empty()) {
            cerr << "you used selected pathway, will ignored the --tfmi file you provided";
        }
        opt.mkSelectedPathwayDB();
    }

    if (opt.transSearch.tfmi.empty()) {
        error_exit("you must provide BWTFMI file using --tfmi file or using --pathway and --genefa");
    } else {
        check_file_valid(opt.transSearch.tfmi);
    }

    BwtFmiDB * tbwtfmiDB = new BwtFmiDB(& opt);

    // I/O
    opt.mHomoSearchOptions.sampleTable = cmd.get<string>("sampletable");
    opt.mHomoSearchOptions.prefix = cmd.get<string>("prefix");
    opt.outputMappedCleanReads = cmd.exist("outputMappedCleanReads");
    opt.outputReadsKOMap = cmd.exist("outputReadsKOMap");

    if (opt.mHomoSearchOptions.prefix.empty() && opt.mHomoSearchOptions.sampleTable.empty()) {
        error_exit("You must specify output file using --prefix or using --sampletable, which contains prefix string");
    } else if (!opt.mHomoSearchOptions.prefix.empty() && !opt.mHomoSearchOptions.sampleTable.empty()) {
        error_exit("You must specify output file using either --prefix or --sampletable");
    }

    //for single file;
    if (!opt.mHomoSearchOptions.prefix.empty()) {
        opt.transSearch.startTime = t_begin;
        opt.in1 = cmd.get<string>("in1");
        opt.in2 = cmd.get<string>("in2");
        std::string outFName;
        opt.htmlFile = opt.mHomoSearchOptions.prefix + "_report.html";
        opt.jsonFile = opt.mHomoSearchOptions.prefix + "_report.json";
        if (opt.outputMappedCleanReads) {
            opt.out1 = opt.mHomoSearchOptions.prefix + "_mapped_R1.fastq.gz";
            if (opt.isPaired()) opt.out2 = opt.mHomoSearchOptions.prefix + "_mapped_R2.fastq.gz";
        }

        if (opt.outputReadsKOMap) {
            opt.outReadsKOMap = opt.mHomoSearchOptions.prefix + "_readsKOMap.txt.gz";
        }

        stringstream ss;
        for (int i = 0; i < argc; i++) {
            ss << argv[i] << " ";
        }
        command = ss.str();
        opt.mHomoSearchOptions.commandStr = command;

        bool supportEvaluation = !opt.inputFromSTDIN && opt.in1 != "/dev/stdin";
        Evaluator eva(&opt);
        if (supportEvaluation) {
            eva.evaluateSeqLen();
            if (opt.overRepAnalysis.enabled)
                eva.evaluateOverRepSeqs();
        }

        long readNum = 0;

        // using evaluator to guess how many reads in total
        if (opt.shallDetectAdapter(false)) {
            if (!supportEvaluation)
                cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
            else {
                cerr << "Detecting adapter sequence for read1..." << endl;
                string adapt = eva.evalAdapterAndReadNum(readNum, false);
                if (adapt.length() > 60)
                    adapt.resize(0, 60);
                if (adapt.length() > 0) {
                    opt.adapter.sequence = adapt;
                    opt.adapter.detectedAdapter1 = adapt;
                } else {
                    cerr << "No adapter detected for read1" << endl;
                    opt.adapter.sequence = "";
                }
                cerr << endl;
            }
        }
        if (opt.shallDetectAdapter(true)) {
            if (!supportEvaluation)
                cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
            else {
                cerr << "Detecting adapter sequence for read2..." << endl;
                string adapt = eva.evalAdapterAndReadNum(readNum, true);
                if (adapt.length() > 60)
                    adapt.resize(0, 60);
                if (adapt.length() > 0) {
                    opt.adapter.sequenceR2 = adapt;
                    opt.adapter.detectedAdapter2 = adapt;
                } else {
                    cerr << "No adapter detected for read2" << endl;
                    opt.adapter.sequenceR2 = "";
                }
                cerr << endl;
            }
        }

        opt.validate();

        // using evaluator to guess how many reads in total
        if (opt.split.needEvaluation && supportEvaluation) {
            // if readNum is not 0, means it is already evaluated by other functions
            if (readNum == 0) {
                eva.evaluateReadNum(readNum);
            }
            opt.split.size = readNum / opt.split.number;
            // one record per file at least
            if (opt.split.size <= 0) {
                opt.split.size = 1;
                cerr << "WARNING: the input file has less reads than the number of files to split" << endl;
            }
        }

        // using evaluator to check if it's two color system
        if (!cmd.exist("trim_poly_g") && !cmd.exist("disable_trim_poly_g") && supportEvaluation) {
            bool twoColorSystem = eva.isTwoColorSystem();
            if (twoColorSystem) {
                opt.polyGTrim.enabled = true;
            }
        }

        Processor p(&opt);
        p.process(tbwtfmiDB);

        cerr << endl << command << endl;
        cerr << endl << "Seq2Fun v" << SEQ2FUNR_VER << ", time used: " << convertSeconds(opt.transSearch.timeLapse) << ", mapping " <<
                opt.transSearch.nTransMappedKOReads << " reads out of " <<
                opt.mHomoSearchOptions.nTotalReads << " (" <<
                getPercentage(long(opt.transSearch.nTransMappedKOReads), opt.mHomoSearchOptions.nTotalReads) << " %); " <<
                "mapped " << opt.transSearch.nTransMappedKOs << " KOs out of " <<
                opt.transSearch.nKODB << " KOs (" <<
                getPercentage(opt.transSearch.nTransMappedKOs, opt.transSearch.nKODB) <<
                " %)" << endl << endl;
        
        opt.transSearch.reset2Default();
        opt.mHomoSearchOptions.reset2Default();
    } else {
        opt.parseSampleTable();
        //std::vector< std::pair<std::string, std::unordered_map<std::string, int> > > sampleKOTableVec;
        int count = 0;
        for (auto & it : opt.samples) {
            time_t t_cycleBegin = time(NULL);
            opt.transSearch.startTime = t_cycleBegin;
            count++;
            std::stringstream msgSS;
            msgSS << "processing sample: " << it.prefix << ", " << count << " out of " << opt.samples.size() << " samples" << "\n";
            loginfo(msgSS.str());
            opt.mHomoSearchOptions.prefix = it.prefix;
            opt.in1 = it.in1;
            opt.in2 = it.in2;
            opt.htmlFile = it.prefix + "_report.html";
            opt.jsonFile = it.prefix + "_report.json";

            if (opt.outputMappedCleanReads) {
                opt.out1 = it.prefix + "_mapped_R1.fastq.gz";
                if (opt.isPaired()) opt.out2 = it.prefix + "_mapped_R2.fastq.gz";
            }

            if (opt.outputReadsKOMap) {
                opt.outReadsKOMap = opt.mHomoSearchOptions.prefix + "_readsKOMap.txt.gz";
            }

            std::stringstream ss;
            for (int i = 0; i < argc; i++) {
                ss << argv[i] << " ";
            }
            command = ss.str();
            opt.mHomoSearchOptions.commandStr = command;

            bool supportEvaluation = !opt.inputFromSTDIN && opt.in1 != "/dev/stdin";
            Evaluator eva(&opt);
            if (supportEvaluation) {
                eva.evaluateSeqLen();
                if (opt.overRepAnalysis.enabled)
                    eva.evaluateOverRepSeqs();
            }

            long readNum = 0;

            // using evaluator to guess how many reads in total
            if (opt.shallDetectAdapter(false)) {
                if (!supportEvaluation)
                    cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
                else {
                    cerr << "Detecting adapter sequence for read1..." << endl;
                    string adapt = eva.evalAdapterAndReadNum(readNum, false);
                    if (adapt.length() > 60)
                        adapt.resize(0, 60);
                    if (adapt.length() > 0) {
                        opt.adapter.sequence = adapt;
                        opt.adapter.detectedAdapter1 = adapt;
                    } else {
                        cerr << "No adapter detected for read1" << endl;
                        opt.adapter.sequence = "";
                    }
                    cerr << endl;
                }
            }
            if (opt.shallDetectAdapter(true)) {
                if (!supportEvaluation)
                    cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
                else {
                    cerr << "Detecting adapter sequence for read2..." << endl;
                    string adapt = eva.evalAdapterAndReadNum(readNum, true);
                    if (adapt.length() > 60)
                        adapt.resize(0, 60);
                    if (adapt.length() > 0) {
                        opt.adapter.sequenceR2 = adapt;
                        opt.adapter.detectedAdapter2 = adapt;
                    } else {
                        cerr << "No adapter detected for read2" << endl;
                        opt.adapter.sequenceR2 = "";
                    }
                    cerr << endl;
                }
            }

            opt.validate();

            // using evaluator to guess how many reads in total
            if (opt.split.needEvaluation && supportEvaluation) {
                // if readNum is not 0, means it is already evaluated by other functions
                if (readNum == 0) {
                    eva.evaluateReadNum(readNum);
                }
                opt.split.size = readNum / opt.split.number;
                // one record per file at least
                if (opt.split.size <= 0) {
                    opt.split.size = 1;
                    cerr << "WARNING: the input file has less reads than the number of files to split" << endl;
                }
            }

            // using evaluator to check if it's two color system
            if (!cmd.exist("trim_poly_g") && !cmd.exist("disable_trim_poly_g") && supportEvaluation) {
                bool twoColorSystem = eva.isTwoColorSystem();
                if (twoColorSystem) {
                    opt.polyGTrim.enabled = true;
                }
            }
            
            Processor p(&opt);
            p.process(tbwtfmiDB);
            cerr << endl << command << endl;
            cerr << endl << "Seq2Fun v" << SEQ2FUNR_VER << ", time used: " << convertSeconds(opt.transSearch.timeLapse) <<
                    ", mapping " << opt.transSearch.nTransMappedKOReads << " reads out of " <<
                    opt.mHomoSearchOptions.nTotalReads << " (" <<
                    getPercentage(long(opt.transSearch.nTransMappedKOReads), opt.mHomoSearchOptions.nTotalReads) << " %); " <<
                    "mapped " << opt.transSearch.nTransMappedKOs << " KOs out of " <<
                    opt.transSearch.nKODB << " KOs (" <<
                    getPercentage(opt.transSearch.nTransMappedKOs, opt.transSearch.nKODB) <<
                    " %)" << endl << endl;
            opt.transSearch.reset2Default();
            opt.mHomoSearchOptions.reset2Default();
        }

        HtmlReporterAll hra(&opt);
        hra.report();
        time_t t_total = time(NULL);
        auto timeUsed = difftime(t_total, t_begin);
        cerr << endl << "Seq2Fun v" << SEQ2FUNR_VER << ", time used: " << convertSeconds(timeUsed) << ", processed " << opt.samples.size() << " samples" << endl << endl;
        opt.samples.clear();
        opt.samples.shrink_to_fit();
    }
    if (tbwtfmiDB) delete tbwtfmiDB;
    return 0;
}
