#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include <time.h>
#include <thread>
#include <valarray>
#include <cmath>

#include "cmdline.h"
#include "options.h"
#include "util.h"
#include "seqtractpeprocessor.h"

std::string command;
mutex logmtx;

int main(int argc, char** argv) {

    if (argc == 1) {
        std::cout << "SeqTract: multi-threaded tool to extract sequences from FASTQ files!";
    }

    cmdline::parser cmd;

    cmd.add<string>("sampleTable", 0, "input table for fastq.gz files with annotation from Seq2Fun output, eg: one column for SE samples or two columns (separated by '\t') of samples for PE sample", true, "");
    cmd.add<string>("geneTable", 0, "input table of target genes for sequences extraction, eg: one columns of s2f ids or KO ids", true, "");
    cmd.add<string>("outputDir", 0, "output folder for extract fastq.gz files, default is current working dir", false, ".");
    //cmd.add<string>("suffix", 0, "the suffix appended to the file name", false, "");
    cmd.add<string>("undetermined", 0, "the file name of undetermined data, default is Undetermined", false, "Undetermined");
    cmd.add<int>("compression", 0, "compression level for gzip output (1 - 9), default is 6", false, 6);
    cmd.add<int>("thread", 'w', "worker thread number, default is 2", false, 2);
    cmd.add("verbose", 'V', "enable verbose");

    cmd.parse_check(argc, argv);

    stringstream ss;
    for (int i = 0; i < argc; i++) {
        ss << argv[i] << " ";
    }
    command = ss.str();
    time_t t1 = time(NULL);

    Options opt;
    opt.verbose = cmd.exist("verbose");
    opt.thread = cmd.get<int>("thread");
    int n_t = std::thread::hardware_concurrency();
    opt.thread = std::min(opt.thread, n_t);

    opt.mSeqExtractions.sampleMappedTableStr = cmd.get<string>("sampleTable");
    opt.mSeqExtractions.targetGeneTableStr = cmd.get<string>("geneTable");
    opt.mSeqExtractions.outputDir = cmd.get<string>("outputDir");
    opt.mSeqExtractions.compression = cmd.get<int>("compression");
    opt.mSeqExtractions.undeterminedFileName = cmd.get<string>("undetermined");
    opt.readSampleExtraction();

    int steps = floor(opt.mSeqExtractions.targetGenesVec.size() / opt.thread);
    int rem = opt.mSeqExtractions.targetGenesVec.size() % opt.thread;

    int i = 0;
    for (i = 0; i < steps; i++) {
        int low = i * opt.thread;
        int up = (i + 1) * opt.thread;
        opt.mSeqExtractions.targetGenesSubVec.clear();
        opt.mSeqExtractions.targetGenesSubVec = sliceVec(opt.mSeqExtractions.targetGenesVec, low, up);

        opt.mSeqExtractions.suffix = "_R1.fastq.gz";
        auto tmp = opt.mSeqExtractions.undeterminedFileName + "_" + std::to_string(i) + opt.mSeqExtractions.suffix;
        opt.mSeqExtractions.undeterminedFileNameOut = joinpath(opt.mSeqExtractions.outputDir, tmp);
        opt.mSeqExtractions.samplesVec.clear();
        if (i == 0) {
            opt.mSeqExtractions.samplesVec = opt.mSeqExtractions.samplesVecF;
        } else {
            auto tmp = opt.mSeqExtractions.undeterminedFileName + "_" + std::to_string(i - 1) + opt.mSeqExtractions.suffix;
            opt.mSeqExtractions.undeterminedFileNameIn = joinpath(opt.mSeqExtractions.outputDir, tmp);
            opt.mSeqExtractions.samplesVec.emplace_back(opt.mSeqExtractions.undeterminedFileNameIn);
        }
        SeqTractPeProcessor sqt(& opt);
        sqt.process();
        if (i > 0) {
            std::string msg = "\ndeleting file: " + opt.mSeqExtractions.undeterminedFileNameIn;
            loginfo(msg);
            remove(opt.mSeqExtractions.undeterminedFileNameIn.c_str());
        }

        if (opt.mSeqExtractions.paired) {
            opt.mSeqExtractions.suffix = "_R2.fastq.gz";
            auto tmp = opt.mSeqExtractions.undeterminedFileName + "_" + std::to_string(i) + opt.mSeqExtractions.suffix;
            opt.mSeqExtractions.undeterminedFileNameOut = joinpath(opt.mSeqExtractions.outputDir, tmp);
            opt.mSeqExtractions.samplesVec.clear();
            if (i == 0) {
                opt.mSeqExtractions.samplesVec = opt.mSeqExtractions.samplesVecR;
            } else {
                auto tmp = opt.mSeqExtractions.undeterminedFileName + "_" + std::to_string(i - 1) + opt.mSeqExtractions.suffix;
                opt.mSeqExtractions.undeterminedFileNameIn = joinpath(opt.mSeqExtractions.outputDir, tmp);
                opt.mSeqExtractions.samplesVec.emplace_back(opt.mSeqExtractions.undeterminedFileNameIn);
            }
            SeqTractPeProcessor sqt(& opt);
            sqt.process();
            if (i > 0) {
                std::string msg = "\ndeleting file: " + opt.mSeqExtractions.undeterminedFileNameIn;
                loginfo(msg);
                remove(opt.mSeqExtractions.undeterminedFileNameIn.c_str());
            }
        }
    }

    if (rem != 0) {
        opt.mSeqExtractions.samplesVec.clear();
        int low = i * opt.thread;
        int up = opt.mSeqExtractions.targetGenesVec.size();
        opt.mSeqExtractions.targetGenesSubVec.clear();
        opt.mSeqExtractions.targetGenesSubVec = sliceVec(opt.mSeqExtractions.targetGenesVec, low, up);

        opt.mSeqExtractions.suffix = "_R1.fastq.gz";

        auto tmp = opt.mSeqExtractions.undeterminedFileName + "_" + std::to_string(i) + opt.mSeqExtractions.suffix;
        opt.mSeqExtractions.undeterminedFileNameOut = joinpath(opt.mSeqExtractions.outputDir, tmp);
        opt.mSeqExtractions.samplesVec.clear();
        if (i == 0) {
            opt.mSeqExtractions.samplesVec = opt.mSeqExtractions.samplesVecF;
        } else {
            auto tmp = opt.mSeqExtractions.undeterminedFileName + "_" + std::to_string(i - 1) + opt.mSeqExtractions.suffix;
            opt.mSeqExtractions.undeterminedFileNameIn = joinpath(opt.mSeqExtractions.outputDir, tmp);
            opt.mSeqExtractions.samplesVec.emplace_back(opt.mSeqExtractions.undeterminedFileNameIn);
        }
        SeqTractPeProcessor sqt(& opt);
        sqt.process();

        std::string msg = "\ndeleting file: " + opt.mSeqExtractions.undeterminedFileNameIn;
        loginfo(msg);
        remove(opt.mSeqExtractions.undeterminedFileNameIn.c_str());

        if (i > 0) {
            std::string msg = "\ndeleting file: " + opt.mSeqExtractions.undeterminedFileNameOut;
            loginfo(msg);
            remove(opt.mSeqExtractions.undeterminedFileNameOut.c_str());
        }

        if (opt.mSeqExtractions.paired) {
            opt.mSeqExtractions.suffix = "_R2.fastq.gz";
            auto tmp = opt.mSeqExtractions.undeterminedFileName + "_" + std::to_string(i) + opt.mSeqExtractions.suffix;
            opt.mSeqExtractions.undeterminedFileNameOut = joinpath(opt.mSeqExtractions.outputDir, tmp);
            opt.mSeqExtractions.samplesVec.clear();

            if (i == 0) {
                opt.mSeqExtractions.samplesVec = opt.mSeqExtractions.samplesVecF;
            } else {
                auto tmp = opt.mSeqExtractions.undeterminedFileName + "_" + std::to_string(i - 1) + opt.mSeqExtractions.suffix;
                opt.mSeqExtractions.undeterminedFileNameIn = joinpath(opt.mSeqExtractions.outputDir, tmp);
                opt.mSeqExtractions.samplesVec.emplace_back(opt.mSeqExtractions.undeterminedFileNameIn);
            }

            SeqTractPeProcessor sqt(& opt);
            sqt.process();

            std::string msg = "\ndeleting file: " + opt.mSeqExtractions.undeterminedFileNameIn;
            loginfo(msg);
            remove(opt.mSeqExtractions.undeterminedFileNameIn.c_str());

            if (i > 0) {
                std::string msg = "\ndeleting file: " + opt.mSeqExtractions.undeterminedFileNameOut;
                loginfo(msg);
                remove(opt.mSeqExtractions.undeterminedFileNameOut.c_str());
            }
        }
    }
    
    time_t t2 = time(NULL);

    cerr << endl << command << endl;
    cerr << endl << "SeqTract v" << SEQ2FUNR_VER << ", time used: " << convertSeconds((t2) - t1) <<
            " seconds processed " << opt.mSeqExtractions.targetGenesVec.size() <<
            " samples with " << opt.mSeqExtractions.numFeatures << " features\n\n";

    return 0;
}

