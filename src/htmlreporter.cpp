#include "htmlreporter.h"
#include <chrono>
#include <memory.h>
#include <valarray>

extern string command;

HtmlReporter::HtmlReporter(Options* & opt){
    mOptions = opt;
    mDupHist = NULL;
    mDupRate = 0.0;
}

HtmlReporter::~HtmlReporter(){
}

void HtmlReporter::setDupHist(int* dupHist, double* dupMeanGC, double dupRate) {
    mDupHist = dupHist;
    mDupMeanGC = dupMeanGC;
    mDupRate = dupRate;
}

void HtmlReporter::setInsertHist(atomic_long* insertHist, int insertSizePeak) {
    mInsertHist = insertHist;
    mInsertSizePeak = insertSizePeak;
}

void HtmlReporter::outputRow(ofstream& ofs, string key, long v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + to_string(v) + "</td></tr>\n";
}

void HtmlReporter::outputRow(ofstream& ofs, string key, string v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + v + "</td></tr>\n";
}

void HtmlReporter::outputLongRow(ofstream& ofs, string key, string v) {
    ofs << "<tr><td class='collarge'>" + key + "</td><td class='collarge'>" + v + "</td></tr>\n";
}

string HtmlReporter::formatNumber(long number) {
    double num = (double)number;
    string unit[6] = {"", "K", "M", "G", "T", "P"};
    int order = 0;
    while (num > 1000.0) {
        order += 1;
        num /= 1000.0;
    }

    if (order == 0)
        return to_string(number);
    else
        return to_string(num) + " " + unit[order];
}

string HtmlReporter::getPercents(long numerator, long denominator) {
    if(denominator == 0)
        return "0.0";
    else
        return to_string((double)numerator * 100.0 / (double)denominator);
}

void HtmlReporter::printSummary(ofstream& ofs, FilterResult* result, Stats* preStats1, Stats* postStats1, Stats* preStats2, Stats* postStats2) {
    long pre_total_reads = preStats1->getReads();
    if(preStats2)
        pre_total_reads += preStats2->getReads();

    long pre_total_bases = preStats1->getBases();
    if(preStats2)
        pre_total_bases += preStats2->getBases();

    long pre_q20_bases = preStats1->getQ20();
    if(preStats2)
        pre_q20_bases += preStats2->getQ20();

    long pre_q30_bases = preStats1->getQ30();
    if(preStats2)
        pre_q30_bases += preStats2->getQ30();

    long pre_total_gc = preStats1->getGCNumber();
    if(preStats2)
        pre_total_gc += preStats2->getGCNumber();

    long post_total_reads = postStats1->getReads();
    if(postStats2)
        post_total_reads += postStats2->getReads();

    long post_total_bases = postStats1->getBases();
    if(postStats2)
        post_total_bases += postStats2->getBases();

    long post_q20_bases = postStats1->getQ20();
    if(postStats2)
        post_q20_bases += postStats2->getQ20();

    long post_q30_bases = postStats1->getQ30();
    if(postStats2)
        post_q30_bases += postStats2->getQ30();

    long post_total_gc = postStats1->getGCNumber();
    if(postStats2)
        post_total_gc += postStats2->getGCNumber();

    string sequencingInfo  = mOptions->isPaired()?"paired end":"single end";
    if(mOptions->isPaired()) {
        sequencingInfo += " (" + to_string(preStats1->getCycles()) + " cycles + " + to_string(preStats2->getCycles()) + " cycles)";
    } else {
        sequencingInfo += " (" + to_string(preStats1->getCycles()) + " cycles)";
    }

    ofs << endl;
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('summary')><a name='summary'>Data QC summary <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='summary'>\n";

    ofs << "<div class='subsection_title' onclick=showOrHide('general')>General</div>\n";
    ofs << "<div id='general'>\n";
    ofs << "<table class='summary_table'>\n";
    outputRow(ofs, "sequencing:", sequencingInfo);

    // report read length change
    if(mOptions->isPaired()) {
        outputRow(ofs, "mean length before filtering:", to_string(preStats1->getMeanLength()) + "bp, " + to_string(preStats2->getMeanLength()) + "bp");
    } else {
        outputRow(ofs, "mean length before filtering:", to_string(preStats1->getMeanLength()) + "bp");
        outputRow(ofs, "mean length after filtering:", to_string(postStats1->getMeanLength()) + "bp");
    }

    if(mOptions->duplicate.enabled) {
        string dupStr = to_string(mDupRate*100) + "%";
        if(!mOptions->isPaired())
            dupStr += " (may be overestimated since this is SE data)";
        outputRow(ofs, "duplication rate:", dupStr);
    }
    if(mOptions->isPaired()) {
        outputRow(ofs, "Insert size peak:", mInsertSizePeak);
    }
    if(mOptions->adapterCuttingEnabled()) {
        if(!mOptions->adapter.detectedAdapter1.empty())
            outputRow(ofs, "Detected read1 adapter:", mOptions->adapter.detectedAdapter1);
        if(!mOptions->adapter.detectedAdapter2.empty())
            outputRow(ofs, "Detected read2 adapter:", mOptions->adapter.detectedAdapter2);
    }
    ofs << "</table>\n";
    ofs << "</div>\n";

    ofs << "<div class='subsection_title' onclick=showOrHide('before_filtering_summary')>Original data</div>\n";
    ofs << "<div id='before_filtering_summary'>\n";
    ofs << "<table class='summary_table'>\n";
    outputRow(ofs, "total reads:", formatNumber(pre_total_reads));
    outputRow(ofs, "total bases:", formatNumber(pre_total_bases));
    outputRow(ofs, "Q20 bases:", formatNumber(pre_q20_bases) + " (" + getPercents(pre_q20_bases,pre_total_bases) + "%)");
    outputRow(ofs, "Q30 bases:", formatNumber(pre_q30_bases) + " (" + getPercents(pre_q30_bases, pre_total_bases) + "%)");
    outputRow(ofs, "GC content:", getPercents(pre_total_gc,pre_total_bases) + "%");
    ofs << "</table>\n";
    ofs << "</div>\n";

    ofs << "<div class='subsection_title' onclick=showOrHide('after_filtering_summary')>Clean data used for detection</div>\n";
    ofs << "<div id='after_filtering_summary'>\n";
    ofs << "<table class='summary_table'>\n";
    outputRow(ofs, "total reads:", formatNumber(post_total_reads));
    outputRow(ofs, "total bases:", formatNumber(post_total_bases));
    outputRow(ofs, "Q20 bases:", formatNumber(post_q20_bases) + " (" + getPercents(post_q20_bases, post_total_bases) + "%)");
    outputRow(ofs, "Q30 bases:", formatNumber(post_q30_bases) + " (" + getPercents(post_q30_bases, post_total_bases) + "%)");
    outputRow(ofs, "GC content:", getPercents(post_total_gc,post_total_bases) + "%");
    ofs << "</table>\n";
    ofs << "</div>\n";

    if(result) {
        ofs << "<div class='subsection_title' onclick=showOrHide('filtering_result')>Filtering result</div>\n";
        ofs << "<div id='filtering_result'>\n";
        result -> reportHtml(ofs, pre_total_reads, pre_total_bases);
        ofs << "</div>\n";
    }

    ofs << "</div>\n";
    ofs << "</div>\n";

    if(result && mOptions->adapterCuttingEnabled()) {
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('adapters')><a name='summary'>Adapters <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
        ofs << "<div id='adapters' style='display:none'>\n";

        result->reportAdapterHtml(ofs, pre_total_bases);

        ofs << "</div>\n";
        ofs << "</div>\n";
    }

    if(mOptions->duplicate.enabled) {
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('duplication')><a name='summary'>Duplication <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
        ofs << "<div id='duplication' style='display:none'>\n";

        reportDuplication(ofs);

        ofs << "</div>\n";
        ofs << "</div>\n";
    }

    if(mOptions->isPaired()) {
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('insert_size')><a name='summary'>Insert size estimation <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
        ofs << "<div id='insert_size' style='display:none'>\n";

        reportInsertSize(ofs, preStats1->getCycles() + preStats2->getCycles() - mOptions->overlapRequire);

        ofs << "</div>\n";
        ofs << "</div>\n";
    }
}

void HtmlReporter::reportInsertSize(ofstream& ofs, int isizeLimit) {
    if(isizeLimit<1)
        isizeLimit = 1;
    int total = min(mOptions->insertSizeMax, isizeLimit);
    long *x = new long[total];
    double allCount = 0;
    for(int i=0; i<total; i++) {
        x[i] = i;
        allCount += mInsertHist[i];
    }
    allCount += mInsertHist[mOptions->insertSizeMax];
    double* percents = new double[total];
    memset(percents, 0, sizeof(double)*total);
    if(allCount > 0) {
        for(int i=0; i<total; i++) {
            percents[i] = (double)mInsertHist[i] * 100.0 / (double)allCount;
        }
    }

    double unknownPercents = (double)mInsertHist[mOptions->insertSizeMax] * 100.0 / (double)allCount;

    ofs << "<div id='insert_size_figure'>\n";
    ofs << "<div class='figure' id='plot_insert_size' style='height:400px;'></div>\n";
    ofs << "</div>\n";

    ofs << "<div class='sub_section_tips'>This estimation is based on paired-end overlap analysis, and there are ";
    ofs << to_string(unknownPercents);
    ofs << "% reads found not overlapped. <br /> The nonoverlapped read pairs may have insert size &lt;" << mOptions->overlapRequire;
    ofs << " or &gt;" << isizeLimit;
    ofs << ", or contain too much sequencing errors to be detected as overlapped.";
    ofs <<"</div>\n";
    
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x, total) + "],";
    json_str += "y:[" + Stats::list2string(percents, total) + "],";
    json_str += "name: 'Percent (%)  ',";
    json_str += "type:'bar',";
    json_str += "line:{color:'rgba(128,0,128,1.0)', width:1}\n";
    json_str += "}";

    json_str += "];\n";

    json_str += "var layout={title:'Insert size distribution (" + to_string(unknownPercents) + "% reads are with unknown length)', xaxis:{title:'Insert size'}, yaxis:{title:'Read percent (%)'}};\n";
    json_str += "Plotly.newPlot('plot_insert_size', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
    delete[] percents;
}

void HtmlReporter::reportDuplication(ofstream& ofs) {

    ofs << "<div id='duplication_figure'>\n";
    ofs << "<div class='figure' id='plot_duplication' style='height:400px;'></div>\n";
    ofs << "</div>\n";
    
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    int total = mOptions->duplicate.histSize - 2;
    long *x = new long[total];
    double allCount = 0;
    for(int i=0; i<total; i++) {
        x[i] = i+1;
        allCount += mDupHist[i+1];
    }
    double* percents = new double[total];
    memset(percents, 0, sizeof(double)*total);
    if(allCount > 0) {
        for(int i=0; i<total; i++) {
            percents[i] = (double)mDupHist[i+1] * 100.0 / (double)allCount;
        }
    }
    int maxGC = total;
    double* gc = new double[total];
    for(int i=0; i<total; i++) {
        gc[i] = (double)mDupMeanGC[i+1] * 100.0;
        // GC ratio will be not accurate if no enough reads to average
        if(percents[i] <= 0.05 && maxGC == total)
            maxGC = i;
    }

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x, total) + "],";
    json_str += "y:[" + Stats::list2string(percents, total) + "],";
    json_str += "name: 'Read percent (%)  ',";
    json_str += "type:'bar',";
    json_str += "line:{color:'rgba(128,0,128,1.0)', width:1}\n";
    json_str += "},";

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x, maxGC) + "],";
    json_str += "y:[" + Stats::list2string(gc, maxGC) + "],";
    json_str += "name: 'Mean GC ratio (%)  ',";
    json_str += "mode:'lines',";
    json_str += "line:{color:'rgba(255,0,128,1.0)', width:2}\n";
    json_str += "}";

    json_str += "];\n";

    json_str += "var layout={title:'duplication rate (" + to_string(mDupRate*100.0) + "%)', xaxis:{title:'duplication level'}, yaxis:{title:'Read percent (%) & GC ratio'}};\n";
    json_str += "Plotly.newPlot('plot_duplication', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
    delete[] percents;
    delete[] gc;
}

void HtmlReporter::reportRarefaction(ofstream& ofs){

    int top = mOptions->transSearch.rarefactionMap.size();
    std::vector<int> x_vec;
    x_vec.reserve(top);
    std::vector<double> y_vec;
    y_vec.reserve(top);

    for (auto & it : mOptions->transSearch.rarefactionMap) {
        x_vec.push_back(it.first);
        y_vec.push_back((double)it.second);
    }
    
    mOptions->transSearch.rarefactionMap.clear();
    
    ofs << "<div id='rare_figure'>\n";
    ofs << "<div class='figure' id='plot_rarefaction_curve' style='height:400px;'></div>\n";
    ofs << "</div>\n";
    
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x_vec, top) + "],";
    json_str += "y:[" + Stats::list2string(y_vec, top) + "],";
    json_str += "name: 'Rarefaction curve  ',";
    json_str += "type:'scatter',";
    json_str += "}";
    json_str += "];\n";

    json_str += "var layout={title:'Rarefaction curve ', xaxis:{title:'Number of reads', automargin: true}, yaxis:{title:'Number of KOs', automargin: true}};\n";
    json_str += "Plotly.newPlot('plot_rarefaction_curve', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    x_vec.clear();
    x_vec.shrink_to_fit();
    y_vec.clear();
    y_vec.shrink_to_fit();
}

void HtmlReporter::reportKOBarPlot(ofstream& ofs){
    std::vector<std::string> x_vec;
    std::vector<double> y_vec;
    std::vector<std::string> y_lable_vec;
    
    int total = mOptions->transSearch.sortedKOFreqTupleVector.size();
    int top = 50;
    top = min(top, total);

    for (int i = 0; i < total; i++) {
        auto it = mOptions->transSearch.sortedKOFreqTupleVector.at(i);
        x_vec.push_back(get<0>(it));
        y_vec.push_back((double) get<1>(it));
        y_lable_vec.push_back(get<2>(it));
    }
    
    mOptions->transSearch.sortedKOFreqTupleVector.clear();
    mOptions->transSearch.sortedKOFreqTupleVector.shrink_to_fit();

    ofs << "<div class='subsection_title' onclick=showOrHide('plot_ko_hits')><a name='summary'>Top-hit KO bar plot<font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div class='figure' id='plot_ko_hits' style='height:400px;'></div>\n";

    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x_vec, top) + "],";
    json_str += "y:[" + Stats::list2string(y_vec, top) + "],";
    json_str += "name: 'Number of hitted kos  ',";
    json_str += "type:'bar',";
    json_str += "line:{color:'rgba(128,0,128,1.0)', width:1}\n";
    json_str += "}";
    json_str += "];\n";
    json_str += "var layout={title:'Top " + to_string(top) + " abundant KOs ', xaxis:{title:'KO', automargin: true}, yaxis:{title:'Number of reads hits', automargin: true}};\n";
    json_str += "Plotly.newPlot('plot_ko_hits', data, layout);\n";
    ofs << json_str;
    ofs << "</script>" << endl;
    
    ofs << "<div class='subsection_title' onclick=showOrHide('top_ko_table')><a name='summary'>KOs Table (full list)<font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='top_ko_table' style='overflow-y:auto; height: 400px;'>\n";
    ofs << "<table class='summary_table'>\n";
    ofs << "<tr><td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << "KO ID" << "</td><td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << "Number" << "</td><td class='collarge' style='font-size:14px;color:#ffffff;background:#008000'>" << "Name" << "</td></tr>\n";
    for (int i = 0; i < total; i++) {
        ofs << "<tr><td class='ko_col'>" << x_vec[i] << "</td><td class='ko_col'>" << y_vec[i] << "</td><td class='collarge'>" << y_lable_vec[i] << "</td></tr>\n";
    }
    
    ofs << "</table>\n";
    ofs << "</div>\n";
    
    x_vec.clear();
    x_vec.shrink_to_fit();
    y_vec.clear();
    y_vec.shrink_to_fit();
    y_lable_vec.clear();
    y_lable_vec.shrink_to_fit();
}

void HtmlReporter::reportRarefactionId(ofstream& ofs){

    int top = mOptions->transSearch.rarefactionIdMap.size();
    std::vector<int> x_vec;
    x_vec.reserve(top);
    std::vector<double> y_vec;
    y_vec.reserve(top);

    for (auto & it : mOptions->transSearch.rarefactionIdMap) {
        x_vec.push_back(it.first);
        y_vec.push_back((double)it.second);
    }
    
    mOptions->transSearch.rarefactionIdMap.clear();
    
    ofs << "<div id='rare_figure'>\n";
    ofs << "<div class='figure' id='plot_rarefaction_curve_id' style='height:400px;'></div>\n";
    ofs << "</div>\n";
    
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x_vec, top) + "],";
    json_str += "y:[" + Stats::list2string(y_vec, top) + "],";
    json_str += "name: 'Rarefaction curve for s2f id  ',";
    json_str += "type:'scatter',";
    json_str += "}";
    json_str += "];\n";

    json_str += "var layout={title:'Rarefaction curve for s2f id', xaxis:{title:'Number of reads', automargin: true}, yaxis:{title:'Number of s2f ids', automargin: true}};\n";
    json_str += "Plotly.newPlot('plot_rarefaction_curve_id', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    x_vec.clear();
    x_vec.shrink_to_fit();
    y_vec.clear();
    y_vec.shrink_to_fit();
}

void HtmlReporter::reportBarPlotId(ofstream& ofs){
//    std::vector<std::string> x_vec;
//    std::vector<double> y_vec;
//    std::vector<std::string> y_lable_vec;
//    
//    int total = mOptions->transSearch.sortedIdFreqTupleVector.size();
//    int top = 50;
//    top = min(top, total);
//
//    for (int i = 0; i < total; i++) {
//        auto it = mOptions->transSearch.sortedIdFreqTupleVector.at(i);
//        x_vec.push_back(get<0>(it));
//        y_vec.push_back((double) get<1>(it));
//        y_lable_vec.push_back(get<2>(it));
//    }
//    
//    mOptions->transSearch.sortedIdFreqTupleVector.clear();
//    mOptions->transSearch.sortedIdFreqTupleVector.shrink_to_fit();
//
//    ofs << "<div class='subsection_title' onclick=showOrHide('plot_id_hits')><a name='summary'>Top-hit s2f id bar plot<font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
//    ofs << "<div class='figure' id='plot_id_hits' style='height:400px;'></div>\n";
//
//    ofs << "\n<script type=\"text/javascript\">" << endl;
//    string json_str = "var data=[";
//
//    json_str += "{";
//    json_str += "x:[" + Stats::list2string(x_vec, top) + "],";
//    json_str += "y:[" + Stats::list2string(y_vec, top) + "],";
//    json_str += "name: 'Number of hitted s2f ids  ',";
//    json_str += "type:'bar',";
//    json_str += "line:{color:'rgba(128,0,128,1.0)', width:1}\n";
//    json_str += "}";
//    json_str += "];\n";
//    json_str += "var layout={title:'Top " + to_string(top) + " abundant s2f ids ', xaxis:{title:'s2f id', automargin: true}, yaxis:{title:'Number of reads hits', automargin: true}};\n";
//    json_str += "Plotly.newPlot('plot_id_hits', data, layout);\n";
//    ofs << json_str;
//    ofs << "</script>" << endl;
    
    ofs << "<div class='subsection_title' onclick=showOrHide('top_id_table')><a name='summary'>S2f ids Table (full list)<font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='top_id_table' style='overflow-y:auto; height: 400px;'>\n";
    ofs << "<table class='summary_table'>\n";
    ofs << "<tr><td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << "s2f_id" 
            << "</td><td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << "Number" 
            << "</td><td class='collarge' style='font-size:14px;color:#ffffff;background:#008000'>" << "KO" 
            << "</td><td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << "GO" 
            << "</td><td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << "Symbol" 
            << "</td><td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << "Gene" 
            << "</td></tr>\n";
    
    std::string ko, go, symbol, gene;
    uint32 count;
    for(const auto & it : mOptions->transSearch.totalIdFreqUMapResults){
        auto itt = mOptions->mHomoSearchOptions.fullDbMap.find(it.first);
        if(itt == mOptions->mHomoSearchOptions.fullDbMap.end()){
            count = 0; ko = "U"; go = "U"; symbol = "U"; gene = "U";
        } else {
            count = it.second;
            ko = itt->second.ko;
            go = itt->second.go;
            symbol = itt->second.symbol;
            gene = itt->second.gene;
        }
        ofs << "<tr><td class='ko_col'>" << "s2f_" << std::setfill('0') << std::setw(10) << *it.first 
                << "</td><td class='ko_col'>" << count 
                << "</td><td class='ko_col'>" << ko 
                << "</td><td class='collarge'>" << go 
                << "</td><td class='ko_col'>" << symbol 
                << "</td><td class='ko_col'>" << gene 
                << "</td></tr>\n";
    }
    ko = ""; go = ""; symbol = ""; gene = "";
    mOptions->transSearch.totalIdFreqUMapResults.clear();
//    for (int i = 0; i < total; i++) {
//        ofs << "<tr><td class='ko_col'>" << x_vec[i] 
//                << "</td><td class='ko_col'>" << y_vec[i] 
//                << "</td><td class='collarge'>" << y_lable_vec[i] 
//                << "</td></tr>\n";
//    }
    
    ofs << "</table>\n";
    ofs << "</div>\n";
    
//    x_vec.clear();
//    x_vec.shrink_to_fit();
//    y_vec.clear();
//    y_vec.shrink_to_fit();
//    y_lable_vec.clear();
//    y_lable_vec.shrink_to_fit();
}

void HtmlReporter::reportPathway(ofstream& ofs){
    
    int total = mOptions->transSearch.sortedPathwayFreqTupleVector.size();
    int top = 30;
    top = min(top, total);
    
    std::vector<std::string> x_vec;
    x_vec.reserve(total);
    std::vector<double> y_vec;
    y_vec.reserve(total);
    std::vector<std::string> y_lable_vec;
    y_lable_vec.reserve(total);
    
    std::vector<int> z_mapped_vec;
    z_mapped_vec.reserve(total);
    std::vector<int> z_total_vec;
    z_total_vec.reserve(total);

    for (int i = 0; i < total; i++) {
        auto it = mOptions->transSearch.sortedPathwayFreqTupleVector.at(i);
        x_vec.push_back(get<0>(it));
        y_vec.push_back((double) get<1>(it));
        y_lable_vec.push_back(get<2>(it));
        z_mapped_vec.push_back(get<3>(it));
        z_total_vec.push_back(get<4>(it));
    }
    
    mOptions->transSearch.sortedPathwayFreqTupleVector.clear();
    mOptions->transSearch.sortedPathwayFreqTupleVector.shrink_to_fit();
    
    ofs << "<div class='subsection_title' onclick=showOrHide('plot_pathway_hits')><a name='summary'>Top-hit pathway bar plot<font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div class='figure' id='plot_pathway_hits' style='height:600px;'></div>\n";
    
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x_vec, top) + "],";
    json_str += "y:[" + Stats::list2string(y_vec, top) + "],";
    json_str += "text: [" + Stats::list2string(y_lable_vec, top) + "],";
    json_str += "name: 'Percentage of KOs in that pathways (%)  ',";
    json_str += "type:'bar',";
    json_str += "line:{color:'rgba(128,0,128,1.0)', width:1}\n";
    json_str += "}";
    json_str += "];\n";

    json_str += "var layout={title:'Top-hit pathways (" + to_string(top) + ")', xaxis:{title:'Pahtways', tickangle: 45}, yaxis:{title:'Percentage of number of hit KOs'}, margin: {b: 200}};\n";
    
    json_str += "Plotly.newPlot('plot_pathway_hits', data, layout);\n";
    ofs << json_str;
    ofs << "</script>" << endl;
    
    ofs << "<div class='subsection_title' onclick=showOrHide('top_pathway_table')><a name='summary'>All hit pathway table<font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='top_pathway_table' style='overflow-y:auto; height: 400px;'>\n";
    ofs << "<table class='summary_table'>\n";
    ofs << "<tr><td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << "Pathway ID" << "</td><td class='collarge' style='font-size:14px;color:#ffffff;background:#008000'>" << "Mapped KOs in pathway (%)" << "</td><td class='collarge' style='font-size:14px;color:#ffffff;background:#008000'>" << "Name" << "</td></tr>\n";

    for (int i = 0; i < total; i++) {
        std::stringstream ss;
        ss << y_vec[i] << " (" << z_mapped_vec[i] << "/" << z_total_vec[i] << ")";
        ofs << "<tr><td class='ko_col'>" <<  x_vec[i] << "</td><td class='ko_col'>" << ss.str() << "</td><td class='collarge'>" << y_lable_vec[i] << "</td></tr>\n";
    }
    
    ofs << "</table>\n";
    ofs << "</div>\n";
    
    x_vec.clear();
    x_vec.shrink_to_fit();
    y_vec.clear();
    y_vec.shrink_to_fit();
    y_lable_vec.clear();
    y_lable_vec.shrink_to_fit();
    z_mapped_vec.clear();
    z_mapped_vec.shrink_to_fit();
    z_total_vec.clear();
    z_total_vec.shrink_to_fit();
}

void HtmlReporter::reportSpecies(ofstream& ofs){

    std::vector<std::string> x_vec;
    std::vector<double> y_vec;
    
    int total = mOptions->transSearch.sortedOrgFreqVec.size();
    int top = 30;
    top = min(top, total);

    for (int i = 0; i < total; i++) {
        auto it = mOptions->transSearch.sortedOrgFreqVec.at(i);
        x_vec.push_back(it.first);
        y_vec.push_back((double) it.second);
    }
    mOptions->transSearch.sortedOrgFreqVec.clear();
    mOptions->transSearch.sortedOrgFreqVec.shrink_to_fit();

    ofs << "<div class='subsection_title' onclick=showOrHide('plot_species_hits')><a name='summary'>Top-hit Species bar plot<font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div class='figure' id='plot_species_hits' style='height:400px;'></div>\n";

    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x_vec, top) + "],";
    json_str += "y:[" + Stats::list2string(y_vec, top) + "],";
    json_str += "name: 'Number of hitted kos  ',";
    json_str += "type:'bar',";
    json_str += "line:{color:'rgba(128,0,128,1.0)', width:1}\n";
    json_str += "}";
    json_str += "];\n";
    json_str += "var layout={title:'Top " + to_string(top) + " hit species ', xaxis:{title:'Species', automargin: true}, yaxis:{title:'Number of KOs', automargin: true}};\n";
    json_str += "Plotly.newPlot('plot_species_hits', data, layout);\n";
    ofs << json_str;
    ofs << "</script>" << endl;
    
    ofs << "<div class='subsection_title' onclick=showOrHide('top_species_table')><a name='summary'>Hit species Table (full list)<font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='top_species_table' style='overflow-y:auto; height: 400px;'>\n";
    ofs << "<table class='summary_table'>\n";
    ofs << "<tr><td class='colmedium' style='font-size:14px;color:#ffffff;background:#008000'>" << "Species" << "</td><td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << "Number of KOs" << "</td></tr>\n";
    for (int i = 0; i < total; i++) {
        ofs << "<tr><td class='colmedium'>" << x_vec[i] << "</td><td class='ko_col'>" << y_vec[i] << "</td></tr>\n";
    }
    
    ofs << "</table>\n";
    ofs << "</div>\n";
    
    x_vec.clear();
    x_vec.shrink_to_fit();
    y_vec.clear();
    y_vec.shrink_to_fit();
}

void HtmlReporter::printAnnotationResults(ofstream& ofs) {
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('result')><a name='result'>Functional quantification results: <I>" << "</I><font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='result'>\n";

    ofs << "<div id='detection_result'>\n";
    ofs << "<table class='summary_table' style='width:1000px'>\n";
    
    outputLongRow(ofs, "Sample", mOptions->mHomoSearchOptions.prefix);

    if (mOptions->samples.size() > 0) {
        int sampleId = mOptions->getWorkingSampleId(mOptions->mHomoSearchOptions.prefix); // get the working sample id;
        outputLongRow(ofs, "Class", mOptions->samples.at(sampleId).feature);
    }
    
    double ratioId = (mOptions->transSearch.nTransMappedIds * 100) / mOptions->transSearch.nIdDB;
    std::string rStrResultId = to_string(mOptions->transSearch.nTransMappedIds) + " / " +  to_string(mOptions->transSearch.nIdDB) + " (" + to_string(ratioId) +  "%)";
    
    outputLongRow(ofs, "Number of s2f id mapped / total s2f ids in database", rStrResultId);
    
//    double ratio = (mOptions->transSearch.nTransMappedKOs * 100) / mOptions->transSearch.nKODB;
//    std::string rStrResult = to_string(mOptions->transSearch.nTransMappedKOs) + " / " +  to_string(mOptions->transSearch.nKODB) + " (" + to_string(ratio) +  "%)";
//    
//    outputLongRow(ofs, "Number of KOs annotated / total KOs in database", rStrResult);
//
//    if (mOptions->mHomoSearchOptions.profiling) {
//        outputLongRow(ofs, "Number of hit species", to_string(mOptions->transSearch.nMappedOrgs));
//        outputLongRow(ofs, "Number of hit pathways", to_string(mOptions->transSearch.nMappedPathways));
//    }
    
//    ratio = (mOptions->transSearch.nTransMappedKOReads * 100) / mOptions->mHomoSearchOptions.nTotalReads;
//    rStrResult = to_string(mOptions->transSearch.nTransMappedKOReads) + " / " +  to_string(mOptions->mHomoSearchOptions.nTotalReads) + " (" + to_string(ratio) +  "%)";
//    outputLongRow(ofs, "Number of annotated reads / total raw reads", rStrResult);
    
    double ratio = (mOptions->mHomoSearchOptions.nCleanReads * 100) / mOptions->mHomoSearchOptions.nTotalReads;
    std::string rStrResult = to_string(mOptions->mHomoSearchOptions.nCleanReads) + " / " +  to_string(mOptions->mHomoSearchOptions.nTotalReads) + " (" + to_string(ratio) +  "%)";
    outputLongRow(ofs, "Number of clean reads / total raw reads", rStrResult);

    auto sTime = ctime(& mOptions->transSearch.startTime);
    auto eTime = ctime(& mOptions->transSearch.endTime);
    rStrResult = convertSeconds(mOptions->transSearch.timeLapse) + "; " + unkown2Str(sTime);
    outputLongRow(ofs, "time used", rStrResult);
    
    outputLongRow(ofs, "command", mOptions->mHomoSearchOptions.commandStr);
    
    ofs << "</table>\n";

    ofs << "</div>\n";
    ofs << "</div>\n";

    if (mOptions->mHomoSearchOptions.profiling) {
        
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('rarefactionId')><a name='summary'>Rarefaction curve for s2f id</a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
        ofs << "<div id='rarefactionid'>\n";
        reportRarefactionId(ofs);
        ofs << "</div>\n";
        ofs << "</div>\n";

        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('top_ids')><a name='summary'>Abundance of s2f ids</a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
        ofs << "<div id='top_ids'>\n";
        reportBarPlotId(ofs);
        ofs << "</div>\n";
        ofs << "</div>\n";
        
//        ofs << "<div class='section_div'>\n";
//        ofs << "<div class='section_title' onclick=showOrHide('rarefaction')><a name='summary'>Rarefaction curve</a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
//        ofs << "<div id='rarefaction'>\n";
//        reportRarefaction(ofs);
//        ofs << "</div>\n";
//        ofs << "</div>\n";
//
//        ofs << "<div class='section_div'>\n";
//        ofs << "<div class='section_title' onclick=showOrHide('top_kos')><a name='summary'>Top abundant KOs</a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
//        ofs << "<div id='top_kos'>\n";
//        reportKOBarPlot(ofs);
//        ofs << "</div>\n";
//        ofs << "</div>\n";
        
//        ofs << "<div class='section_div'>\n";
//        ofs << "<div class='section_title' onclick=showOrHide('top_pathways')><a name='summary'>Top-hit pathways</a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
//        ofs << "<div id='top_pathways'>\n";
//
//        reportPathway(ofs);
//
//        ofs << "</div>\n";
//        ofs << "</div>\n";
//
//        ofs << "<div class='section_div'>\n";
//        ofs << "<div class='section_title' onclick=showOrHide('hitted_species')><a name='summary'>Top-hit species</a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
//        ofs << "<div id='hitted_species'>\n";
//
//        reportSpecies(ofs);
//
//        ofs << "</div>\n";
//        ofs << "</div>\n";
    } else {
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('top_ids')><a name='summary'>Top abundant s2f ids</a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
        ofs << "<div id='top_ids'>\n";
        reportBarPlotId(ofs);
        ofs << "</div>\n";
        ofs << "</div>\n";
    }
}

void HtmlReporter::report(FilterResult* result, Stats* preStats1, Stats* postStats1, Stats* preStats2, Stats* postStats2) {
    ofstream ofs;
    ofs.open(mOptions->htmlFile, ifstream::out);

    printHeader(ofs);

    printAnnotationResults(ofs);
    
    printSummary(ofs, result, preStats1, postStats1, preStats2, postStats2);

    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('before_filtering')><a name='summary'>Original data <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='before_filtering'  style='display:none'>\n";

    if(preStats1) {
        preStats1 -> reportHtml(ofs, "Original data", "read1");
    }

    if(preStats2) {
        preStats2 -> reportHtml(ofs, "Original data", "read2");
    }

    ofs << "</div>\n";
    ofs << "</div>\n";

    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('after_filtering')><a name='summary'>Clean data used for detection <font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='after_filtering'  style='display:none'>\n";

    if(postStats1) {  
        string name = "read1";
        postStats1 -> reportHtml(ofs, "Clean data used for detection", name);
    }

    if(postStats2) {
        postStats2 -> reportHtml(ofs, "Clean data used for detection", "read2");
    }

    ofs << "</div>\n";
    ofs << "</div>\n";

    printFooter(ofs);

}

void HtmlReporter::printHeader(ofstream& ofs){
    ofs << "<html><head><meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />";
    ofs << "<title>Seq2FunR report at " + getCurrentSystemTime() + " </title>";
    printJS(ofs);
    printCSS(ofs);
    ofs << "</head>";
    ofs << "<body><div id='container'>";
}

void HtmlReporter::printCSS(ofstream& ofs){
    ofs << "<style type=\"text/css\">" << endl;
    ofs << "td {border:1px solid #dddddd;padding:5px;font-size:12px;}" << endl;
    ofs << "table {border:1px solid #999999;padding:2x;border-collapse:collapse; width:800px}" << endl;
    ofs << ".col1 {width:400px; font-weight:bold;}" << endl;
    ofs << ".adapter_col {width:800px; font-size:10px;}" << endl;
    ofs << ".ko_col {width:200px; font-weight:bold;}" << endl;
    ofs << ".collarge {width:600px; font-weight:bold;}" << endl;
    ofs << ".colmedium {width:400px; font-weight:bold;}" << endl;
    ofs << "img {padding:30px;}" << endl;
    ofs << "#menu {font-family:Consolas, 'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << "#menu a {color:#0366d6; font-size:18px;font-weight:600;line-height:28px;text-decoration:none;font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol'}" << endl;
    ofs << "a:visited {color: #999999}" << endl;
    ofs << ".alignleft {text-align:left;}" << endl;
    ofs << ".alignright {text-align:right;}" << endl;
    ofs << ".figure {width:1200px;height:600px;}" << endl;
    ofs << ".header {color:#ffffff;padding:1px;height:20px;background:#000000;}" << endl;
    ofs << ".section_title {color:#ffffff;font-size:20px;padding:5px;text-align:left;background:#008000; margin-top:10px;}" << endl;
    ofs << ".subsection_title {font-size:16px;padding:5px;margin-top:10px;text-align:left;color:#008000}" << endl;
    ofs << "#container {text-align:center;padding:3px 3px 3px 10px;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << ".menu_item {text-align:left;padding-top:5px;font-size:18px;}" << endl;
    ofs << ".highlight {text-align:left;padding-top:30px;padding-bottom:30px;font-size:20px;line-height:35px;}" << endl;
    ofs << "#helper {text-align:left;border:1px dotted #fafafa;color:#777777;font-size:12px;}" << endl;
    ofs << "#footer {text-align:left;padding:15px;color:#ffffff;font-size:10px;background:#008000;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << ".kmer_table {text-align:center;font-size:8px;padding:2px;}" << endl;
    ofs << ".kmer_table td{text-align:center;font-size:8px;padding:0px;color:#ffffff}" << endl;
    ofs << ".sub_section_tips {color:#999999;font-size:10px;padding-left:5px;padding-bottom:3px;}" << endl;
    ofs << "</style>" << endl;
}

//void HtmlReporter::printCSS(ofstream& ofs){
//    ofs << "<style type=\"text/css\">" << endl;
//    ofs << "td {border:1px solid #dddddd;padding:5px;font-size:12px;}" << endl;
//    ofs << "table {border:1px solid #999999;padding:2x;border-collapse:collapse; width:800px}" << endl;
//    ofs << ".col1 {width:400px; font-weight:bold;}" << endl;
//    ofs << ".adapter_col {width:800px; font-size:10px;}" << endl;
//    ofs << ".ko_col {width:200px; font-weight:bold;}" << endl;
//    ofs << ".collarge {width:600px; font-weight:bold;}" << endl;
//    ofs << ".colmedium {width:400px; font-weight:bold;}" << endl;
//    ofs << "img {padding:30px;}" << endl;
//    ofs << "#menu {font-family:Consolas, 'Liberation Mono', Menlo, Courier, monospace;}" << endl;
//    ofs << "#menu a {color:#0366d6; font-size:18px;font-weight:600;line-height:28px;text-decoration:none;font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol'}" << endl;
//    ofs << "a:visited {color: #999999}" << endl;
//    ofs << ".alignleft {text-align:left;}" << endl;
//    ofs << ".alignright {text-align:right;}" << endl;
//    ofs << ".figure {width:1200px;height:600px;}" << endl;
//    ofs << ".header {color:#ffffff;padding:1px;height:20px;background:#000000;}" << endl;
//    ofs << ".section_title {color:#ffffff;font-size:20px;padding:5px;text-align:left;background:#008000; margin-top:10px;}" << endl;
//    ofs << ".subsection_title {font-size:16px;padding:5px;margin-top:10px;text-align:left;color:#008000}" << endl;
//    ofs << "#container {text-align:center;padding:3px 3px 3px 10px;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
//    ofs << ".menu_item {text-align:left;padding-top:5px;font-size:18px;}" << endl;
//    ofs << ".highlight {text-align:left;padding-top:30px;padding-bottom:30px;font-size:20px;line-height:35px;}" << endl;
//    ofs << "#helper {text-align:left;border:1px dotted #fafafa;color:#777777;font-size:12px;}" << endl;
//    ofs << "#footer {text-align:left;padding:15px;color:#ffffff;font-size:10px;background:#008000;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
//    ofs << ".kmer_table {text-align:center;font-size:8px;padding:2px;}" << endl;
//    ofs << ".kmer_table td{text-align:center;font-size:8px;padding:0px;color:#ffffff}" << endl;
//    ofs << ".sub_section_tips {color:#999999;font-size:10px;padding-left:5px;padding-bottom:3px;}" << endl;
//    ofs << "</style>" << endl;
//}

const string HtmlReporter::getCurrentSystemTime(){
  auto tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  struct tm* ptm = localtime(&tt);
  char date[60] = {0};
  sprintf(date, "%d-%02d-%02d      %02d:%02d:%02d",
    (int)ptm->tm_year + 1900,(int)ptm->tm_mon + 1,(int)ptm->tm_mday,
    (int)ptm->tm_hour,(int)ptm->tm_min,(int)ptm->tm_sec);
  return std::string(date);
}

void HtmlReporter::printFooter(ofstream& ofs){
    ofs << "\n</div>" << endl;
    ofs << "<div id='footer'> ";
    ofs << "<p>"<<command<<"</p>";
    ofs << "Seq2FunR " << SEQ2FUNR_VER << ", at " << getCurrentSystemTime() << " </div>";
    ofs << "</body></html>";
}

void HtmlReporter::printJS(ofstream& ofs){
    ofs << "<script src='https://www.seq2fun.ca/resources/javascript/plotly-1.2.0.min.js'></script>" << endl;
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << "    function showOrHide(divname) {" << endl;
    ofs << "        div = document.getElementById(divname);" << endl;
    ofs << "        if(div.style.display == 'none')" << endl;
    ofs << "            div.style.display = 'block';" << endl;
    ofs << "        else" << endl;
    ofs << "            div.style.display = 'none';" << endl;
    ofs << "    }" << endl;
    ofs << "</script>" << endl;
}
