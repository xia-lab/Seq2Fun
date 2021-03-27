#ifndef HTMLREPORTERALL_H
#define HTMLREPORTERALL_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <chrono>
#include <memory.h>
#include <valarray>
#include <time.h>

#include "util.h"
#include "common.h"
#include "options.h"


using namespace std;

class HtmlReporterAll {
public:
    HtmlReporterAll(Options* opt);
    ~HtmlReporterAll();
    void report();
    static void outputRow(ofstream& ofs, string key, long value);
    static void outputRow(ofstream& ofs, string key, string value);
    static void outputRow(ofstream& ofs, std::vector<Sample> & samplesVec);
    static void outputSummaryTable(ofstream& ofs, std::vector<Sample> & samplesVec);
    
private:
    const string getCurrentSystemTime();
    void printHeader(ofstream& ofs);
    void printCSS(ofstream& ofs);
    void printJS(ofstream& ofs);
    void printFooter(ofstream& ofs);
    
    void printAnnotationResults(ofstream & ofs);
    void reportRarefactionKO(ofstream& ofs);
    void reportRarefactionKO3D(ofstream& ofs);
    void reportKOBarPlot(ofstream& ofs);
    void reportPathwayBarPlot(ofstream& ofs);
    void reportOrgBarPlot(ofstream& ofs);
    void reportReadsQualityPlot3D(ofstream& ofs);
    void reportAllTables();
    static string list2string(std::vector<long> & x_vec, int top);
    static string list2string(std::vector<int> & x_vec, int top);
    static string list2string(std::vector<double> & x_vec, int top);
    static string list2string(std::vector<string> & x_vec, int top);
    static string list2string2(std::vector<string> & x_vec, int top);
    std::vector<std::string> smNmVec;
    std::vector<std::vector<std::string> > koFreqVec;
    std::vector<std::vector<std::string> > pathwayFreqVec;
    std::vector<std::vector<std::string> > orgFreqVec;
    
private:
    Options * mOptions;

};

#endif /* HTMLREPORTERALL_H */

