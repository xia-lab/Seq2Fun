
#include "htmlreporterall.h"

extern string command;

HtmlReporterAll::HtmlReporterAll(Options * opt) {
    mOptions = opt;
    smNmVec.clear();
}

HtmlReporterAll::~HtmlReporterAll() {
}

void HtmlReporterAll::printHeader(ofstream& ofs){
    ofs << "<html><head><meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />";
    ofs << "<title>Seq2FunM report at " + getCurrentSystemTime() + " </title>";
    printJS(ofs);
    printCSS(ofs);
    ofs << "</head>";
    ofs << "<body><div id='container'>";
}

void HtmlReporterAll::printCSS(ofstream& ofs){
    ofs << "<style type=\"text/css\">" << endl;
    ofs << "td {border:1px solid #dddddd;padding:5px;font-size:12px;}" << endl;
    ofs << "table {border:1px solid #999999;padding:2x;border-collapse:collapse; width:1000px}" << endl;
    ofs << "table thead{position:fixed;}" << endl;
    ofs << ".display"<< endl;
    ofs << ".col1 {width:400px; font-weight:bold;}" << endl;
    ofs << ".adapter_col {width:800px; font-size:10px;}" << endl;
    ofs << ".ko_col {width:200px; font-weight:bold;}" << endl;
    ofs << ".collarge {width:600px; font-weight:bold;}" << endl;
    ofs << ".exlarge {width:1000px; font-weight:bold;}" << endl;
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

const string HtmlReporterAll::getCurrentSystemTime(){
  auto tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  struct tm* ptm = localtime(&tt);
  char date[60] = {0};
  sprintf(date, "%d-%02d-%02d      %02d:%02d:%02d",
    (int)ptm->tm_year + 1900,(int)ptm->tm_mon + 1,(int)ptm->tm_mday,
    (int)ptm->tm_hour,(int)ptm->tm_min,(int)ptm->tm_sec);
  return std::string(date);
}

void HtmlReporterAll::printFooter(ofstream& ofs){
    ofs << "\n</div>" << endl;
    ofs << "<div id='footer'> ";
    ofs << "<p>"<<command<<"</p>";
    ofs << "Seq2FunM " << SEQ2FUNR_VER << ", at " << getCurrentSystemTime() << " </div>";
    ofs << "</body></html>";
}

void HtmlReporterAll::outputRow(ofstream& ofs, string key, long v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + to_string(v) + "</td></tr>\n";
}

void HtmlReporterAll::outputRow(ofstream& ofs, string key, string v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + v + "</td></tr>\n";
}

void HtmlReporterAll::outputRow(ofstream& ofs, std::vector<Sample> & samplesVec) {

    ofs << "<tr>";
    ofs << "<td class='col1'> Id </td>" <<
            "<td class='col1'> sample </td>" <<
            "<td class='col1'> class </td>" <<
            "<td class='col1'> n. S2F ids </td>" <<
            "<td class='col1'> S2F ids rate(%) </td>" <<
            "<td class='col1'> n. S2F ids(db) </td>" <<
            "<td class='col1'> n. reads (S2F id) </td>" <<
            "<td class='col1'> S2F id reads rate(%) </td>" <<
            "<td class='col1'> n. KOs </td>" <<
            "<td class='col1'> KO rate(%) </td>" <<
            "<td class='col1'> n. KOs(db) </td>" <<
            "<td class='col1'> n. reads (KO) </td>" <<
            "<td class='col1'> KO reads rate(%) </td>" <<
            "<td class='col1'> n. clean reads </td>" <<
            "<td class='col1'> clean reads rate(%) </td>" <<
            "<td class='col1'> n. raw reads </td>" <<
            "<td class='col1'> time used</td>" <<
            "<td class='col1'> date </td>";
    ofs << "</tr>\n";   
    
    int id = 0;
    for(Sample & it : samplesVec){
        id++;
        ofs << "<tr>";
        ofs <<  "<td class='col1'>" + to_string(id) + "</td>" << 
                "<td class='col1'>" + basename(it.prefix) + "</td>" << 
                "<td class='col1'>" + it.feature + "</td>" << 
                "<td class='col1'>" + to_string(it.nId) + "</td>" <<
                "<td class='col1'>" + to_string(it.idRate) << "</td>" <<
                "<td class='col1'>" + to_string(it.nIdDb) + "</td>" <<
                "<td class='col1'>" + to_string(it.transSearchMappedIdReads) << "</td>" <<
                "<td class='col1'>" + to_string(it.mappedIdReadsRate) << "</td>" <<
                "<td class='col1'>" + to_string(it.nKO) + "</td>" <<
                "<td class='col1'>" + to_string(it.koRate) << "</td>" <<
                "<td class='col1'>" + to_string(it.nKODb) + "</td>" <<
                "<td class='col1'>" + to_string(it.transSearchMappedKOReads) << "</td>" <<
                "<td class='col1'>" + to_string(it.mappedKOReadsRate) << "</td>" <<
                "<td class='col1'>" + to_string(it.totalCleanReads) << "</td>" <<
                "<td class='col1'>" + to_string(it.cleanReadsRate) << "</td>" <<
                "<td class='col1'>" + to_string(it.totalRawReads) << "</td>" << 
                "<td class='col1'>" + convertSeconds(it.timeLapse) << "</td>" << 
                "<td class='col1'>" + unkown2Str(ctime(&it.startTime)) << "</td>";
        ofs << "</tr>\n";        
    }
}

void HtmlReporterAll::reportAllTables() {
    
    std::set<std::string> koSet, orgSet, pathwaySet, goSet, idSet;
    for (const Sample & sample : mOptions->samples) {
        smNmVec.push_back(basename(sample.prefix));
        for(const auto & it : sample.totalKoFreqUMapResults){
            koSet.insert(it.first);
        }
        for(const auto & it : sample.totalPathwayMap){
            pathwaySet.insert(it.first);
        }
        for(const auto & it : sample.totalOrgKOUMap){
            orgSet.insert(it.first);
        }
        for(const auto & it : sample.totalGoFreqUMapResults){
            goSet.insert(it.first);
        }
        for(const auto & it : sample.totalIdFreqUMapResults){
            idSet.insert(it.first);
        }
    }
    
    //1. for first ko freq;
    std::string fOutNm = joinpath(mOptions->samples.front().path, "All_sample_KO_abundance_table.txt");
    std::ofstream * fOut = NULL;
    fOut = new std::ofstream();
    fOut->open(fOutNm.c_str(), std::ofstream::out);
    if(!fOut->is_open()) error_exit("Can not open all_all_sample_KO_abundance_table.txt");
    if (mOptions->verbose) loginfo("Starting to write all samples KO abundance table");
    
    *fOut << "#NAME\t";
    for(const auto & it : smNmVec){
        *fOut << it << "\t";
    }
    *fOut << "KO_name\n";
    
    *fOut << "#CLASS\t";
    for(const auto & it : mOptions->samples){
        *fOut << it.feature  << "\t";
    }
    *fOut << "class_info\n";
    
    koFreqVec.reserve(koSet.size());
    std::vector<std::string> tmpVec;
    tmpVec.reserve(smNmVec.size() + 2);
    for(const auto & it : koSet){
        *fOut << it << "\t";
        tmpVec.clear();
        tmpVec.push_back(it);
        for (Sample & sample : mOptions->samples) {
            auto itkf = sample.totalKoFreqUMapResults.find(it);
            if(itkf == sample.totalKoFreqUMapResults.end()){
                *fOut << 0 << "\t";
                tmpVec.push_back(to_string(0));
            } else {
                *fOut << itkf->second << "\t";
                tmpVec.push_back(to_string(itkf->second));
            }
        }
        auto itko = mOptions->mHomoSearchOptions.ko_fullname_map.find(it);
        if(itko == mOptions->mHomoSearchOptions.ko_fullname_map.end()){
            *fOut << "UNASSIGNED\n";
            tmpVec.push_back("UNASSIGNED");
        } else {
            *fOut << itko->second << "\n";
            tmpVec.push_back(itko->second);
        }
        koFreqVec.push_back(tmpVec);
        tmpVec.clear();
    }
    fOut->flush();
    fOut->close();
//    if(fOut) delete fOut;
//    fOut = NULL;
    
    fOutNm = joinpath(mOptions->samples.front().path, "All_sample_KO_abundance_table_submit2networkanalyst.txt");
    fOut->open(fOutNm.c_str(), std::ofstream::out);
    if(!fOut->is_open()) error_exit("Can not open All_sample_KO_abundance_table_submit2networkanalyst.txt");
    if (mOptions->verbose) loginfo("Starting to write all samples KO abundance table");
    
    *fOut << "#NAME\t";
    for(const auto & it : smNmVec){
        if(it != smNmVec.back()){
            *fOut << it << "\t";
        } else {
            *fOut << it << "\n";
        }
    }
    
    *fOut << "#CLASS:XX\t";
    for(const auto & it : mOptions->samples){
        if(it.prefix != mOptions->samples.back().prefix){
            *fOut << it.feature  << "\t";
        } else {
            *fOut << it.feature  << "\n";
        }
    }

    for(const auto & it : koSet){
        *fOut << it << "\t";
        for (Sample & sample : mOptions->samples) {
            auto itkf = sample.totalKoFreqUMapResults.find(it);
            if (itkf == sample.totalKoFreqUMapResults.end()) {
                if (sample.prefix != mOptions->samples.back().prefix) {
                    *fOut << 0 << "\t";
                } else {
                    *fOut << 0 << "\n";
                }
            } else {
                if (sample.prefix != mOptions->samples.back().prefix) {
                    *fOut << itkf->second << "\t";
                } else {
                    *fOut << itkf->second << "\n";
                }
            }
        }
    }
    fOut->flush();
    fOut->close();
    
    if (mOptions->verbose) loginfo("Finish to write KO abundance table for all samples");
    
    fOutNm = joinpath(mOptions->samples.front().path, "All_sample_s2fid_abundance_table_submit2networkanalyst.txt");
    fOut->open(fOutNm.c_str(), std::ofstream::out);
    if(!fOut->is_open()) error_exit("Can not open All_sample_s2fid_abundance_table_submit2networkanalyst.txt");
    if (mOptions->verbose) loginfo("Starting to write all samples s2fid abundance table");
    
    *fOut << "#NAME\t";
    for(const auto & it : smNmVec){
        if(it != smNmVec.back()){
            *fOut << it << "\t";
        } else {
            *fOut << it << "\n";
        }
    }
    
    *fOut << "#CLASS:XX\t";
    for(const auto & it : mOptions->samples){
        if(it.prefix != mOptions->samples.back().prefix){
            *fOut << it.feature  << "\t";
        } else {
            *fOut << it.feature  << "\n";
        }
    }

    idFreqVec.reserve(idSet.size());
    for (const auto & it : idSet) {
        *fOut << it << "\t";
        tmpVec.clear();
        tmpVec.push_back(it);
        for (const Sample & sample : mOptions->samples) {
            auto itkf = sample.totalIdFreqUMapResults.find(it);
            if (itkf == sample.totalIdFreqUMapResults.end()) {
                if (sample.prefix != mOptions->samples.back().prefix) {
                    *fOut << 0 << "\t";
                } else {
                    *fOut << 0 << "\n";
                }
                tmpVec.push_back("0");
            } else {
                if (sample.prefix != mOptions->samples.back().prefix) {
                    *fOut << itkf->second << "\t";
                } else {
                    *fOut << itkf->second << "\n";
                }
                tmpVec.push_back(to_string(itkf->second));
            }
        }
        idFreqVec.push_back(tmpVec);
        tmpVec.clear();
    }
    fOut->flush();
    fOut->close();
    
    if (mOptions->verbose) loginfo("Finish to write s2fid abundance table for all samples");
    
    fOutNm = joinpath(mOptions->samples.front().path, "All_sample_GO_abundance_table.txt");
    fOut->open(fOutNm.c_str(), std::ofstream::out);
    if(!fOut->is_open()) error_exit("Can not open All_sample_GO_abundance_table.txt");
    if (mOptions->verbose) loginfo("Starting to write all samples GO abundance table");
    
    *fOut << "#NAME\t";
    for(const auto & it : smNmVec){
        if(it != smNmVec.back()){
            *fOut << it << "\t";
        } else {
            *fOut << it << "\n";
        }
    }
    
    *fOut << "#CLASS:XX\t";
    for(const auto & it : mOptions->samples){
        if(it.prefix != mOptions->samples.back().prefix){
            *fOut << it.feature  << "\t";
        } else {
            *fOut << it.feature  << "\n";
        }
    }

    for(const auto & it : goSet){
        *fOut << it << "\t";
        for (Sample & sample : mOptions->samples) {
            auto itkf = sample.totalGoFreqUMapResults.find(it);
            if (itkf == sample.totalGoFreqUMapResults.end()) {
                if (sample.prefix != mOptions->samples.back().prefix) {
                    *fOut << 0 << "\t";
                } else {
                    *fOut << 0 << "\n";
                }
            } else {
                if (sample.prefix != mOptions->samples.back().prefix) {
                    *fOut << itkf->second << "\t";
                } else {
                    *fOut << itkf->second << "\n";
                }
            }
        }
    }
    fOut->flush();
    fOut->close();
    
    if (mOptions->verbose) loginfo("Finish to write GO abundance table for all samples");
    
    if (mOptions->mHomoSearchOptions.profiling) {
        //2. for pathway
        std::string fOutNm = joinpath(mOptions->samples.front().path, "All_sample_pathway_table.txt");
        //std::ofstream * fOut = new std::ofstream();
        fOut->open(fOutNm, std::ofstream::out);
        if(!fOut->is_open()) error_exit("Can not open All_sample_pathway_table.txt");
        if (mOptions->verbose) loginfo("Starting to write all samples pathway abundance table");

        *fOut << "#Name\t";
        for (auto & it : smNmVec) {
            *fOut << it << "\t";
        }
        *fOut << "n_total_KOs\n";

        *fOut << "#Class\t";
        for (auto & it : mOptions->samples) {
            *fOut << it.feature << "\t";
        }
        *fOut << "class_info\n";
        
        pathwayFreqVec.reserve(pathwaySet.size());
        for (auto & it : pathwaySet) {
            int tmpInt = 0;
            auto itt = mOptions->mHomoSearchOptions.pathway_ko_stats_umap.find(it);
            if (itt == mOptions->mHomoSearchOptions.pathway_ko_stats_umap.end()) {
               tmpInt = 0;
            } else {
               tmpInt = itt->second;
            }
            
            tmpVec.clear();
            *fOut << it << "\t";
            tmpVec.push_back(it);
            
            for (Sample & sample : mOptions->samples) {
                auto itkf = sample.totalPathwayMap.find(it);
                if (itkf == sample.totalPathwayMap.end()) {
                    *fOut << 0 << "\t";
                    tmpVec.push_back("0 (0%)");
                } else {
                    *fOut << itkf->second << "\t";
                    auto str = unkown2Str(itkf->second) + " (" + unkown2Str(double(itkf->second * 100) / double(tmpInt)) + " %)";
                    tmpVec.push_back(str);
                }
            }
            *fOut << tmpInt << "\n";
            tmpVec.push_back(to_string(tmpInt));
//            auto itt = mOptions->mHomoSearchOptions.pathway_ko_stats_umap.find(it);
//            if(itt == mOptions->mHomoSearchOptions.pathway_ko_stats_umap.end()){
//                *fOut << "0\n";
//                tmpVec.push_back("0");
//            } else {
//                *fOut << itt->second << "\n";
//                tmpVec.push_back(to_string(itt->second));
//            }
            pathwayFreqVec.push_back(tmpVec);
            tmpVec.clear();
        }
        fOut->flush();
        fOut->close();
        if (mOptions->verbose) loginfo("Finish to write KO abundance table for all samples");

        //3. for species
        fOutNm.clear();
        fOutNm = joinpath(mOptions->samples.front().path, "All_sample_species_table.txt");
        //std::ofstream * fOut = new std::ofstream();
        fOut->open(fOutNm.c_str(), std::ofstream::out);
        if (!fOut->is_open()) error_exit("Can not open All_sample_species_table.txt");
        if (mOptions->verbose) loginfo("Starting to write all samples species abundance table");

        *fOut << "#Name\t";
        for (auto & it : smNmVec) {
            *fOut << it << "\t";
        }
        *fOut << "\n";

        *fOut << "#Class\t";
        for (auto & it : mOptions->samples) {
            *fOut << it.feature << "\t";
        }
        *fOut << "\n";
        
        orgFreqVec.reserve(orgSet.size());
        for (auto & it : orgSet) {
            tmpVec.clear();
            *fOut << it << "\t";
            tmpVec.push_back(it);
            for (Sample & sample : mOptions->samples) {
                auto itkf = sample.totalOrgKOUMap.find(it);
                if (itkf == sample.totalOrgKOUMap.end()) {
                    *fOut << 0 << "\t";
                    tmpVec.push_back("0");
                } else {
                    *fOut << itkf->second << "\t";
                    tmpVec.push_back(to_string(itkf->second));
                }
            }
            *fOut << "\n";
            orgFreqVec.push_back(tmpVec);
            tmpVec.clear();
        }
        fOut->flush();
        fOut->close();
        if (fOut) delete fOut;
        fOut = NULL;
        if (mOptions->verbose) loginfo("Finish to write spcies table for all samples");
    }
}

void HtmlReporterAll::report(){
    
    reportAllTables();
    
    std::string allHtmlReporter = dirname(mOptions->mHomoSearchOptions.prefix) + "All_samples.html";
    ofstream ofs;
    ofs.open(allHtmlReporter, ifstream::out);
    printHeader(ofs);

    printAnnotationResults(ofs);
    
    printFooter(ofs);
    ofs.close();

}

void HtmlReporterAll::printAnnotationResults(ofstream& ofs) {
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('result')><a name='result'>Functional quantification results: <I>" << "</I><font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
    ofs << "<div id='result'>\n";

    ofs << "<div id='detection_result' style='overflow:auto; height: 400px;'>\n";
    ofs << "<table class='summary_table' style='width:1800px'>\n";
    outputRow(ofs, mOptions->samples);
    ofs << "</table>\n";
    ofs << "</div>\n";
    ofs << "</div>\n";

    if (mOptions->mHomoSearchOptions.profiling) {
        
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('rarefactionS2fid')><a name='summary'>Rarefaction curve (S2F id) <I>" << "</I><font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
        ofs << "<div id='rarefactionS2fid'>\n";
        reportRarefactionS2f(ofs);
        ofs << "</div>\n";
        ofs << "</div>\n";

        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('rarefactionS2fid3d')><a name='summary'>Rarefaction curve (S2F) 3D <I>" << "</I><font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
        ofs << "<div id='rarefactionS2fid3d'>\n";
        reportRarefactionS2f3D(ofs);
        ofs << "</div>\n";
        ofs << "</div>\n";

        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('top_s2fids')><a name='summary'>Top abundant S2F ids</a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
        ofs << "<div id='top_s2fids'>\n";
        reportS2fBarPlot(ofs);
        ofs << "</div>\n";
        ofs << "</div>\n";
        
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('rarefactionKO')><a name='summary'>Rarefaction curve (KO) <I>" << "</I><font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
        ofs << "<div id='rarefactionKO'>\n";
        reportRarefactionKO(ofs);
        ofs << "</div>\n";
        ofs << "</div>\n";

        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('rarefactionKO3d')><a name='summary'>Rarefaction curve (KO) 3D <I>" << "</I><font color='#88CCFF' > (click to show/hide) </font></a></div>\n";
        ofs << "<div id='rarefactionKO3d'>\n";
        reportRarefactionKO3D(ofs);
        ofs << "</div>\n";
        ofs << "</div>\n";

        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('top_kos')><a name='summary'>Top abundant KOs</a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
        ofs << "<div id='top_kos'>\n";
        reportKOBarPlot(ofs);
        ofs << "</div>\n";
        ofs << "</div>\n";

        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('top_pathways')><a name='summary'>Hit pathways</a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
        ofs << "<div id='top_pathways'>\n";
        reportPathwayBarPlot(ofs);
        ofs << "</div>\n";
        ofs << "</div>\n";

        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('top_org')><a name='summary'>Hit Species</a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
        ofs << "<div id='top_org'>\n";
        reportOrgBarPlot(ofs);
        ofs << "</div>\n";
        ofs << "</div>\n";

        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('reads_quality')><a name='summary'>Reads Quality</a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
        ofs << "<div id='reads_quality' style='display:none'>\n";
        reportReadsQualityPlot3D(ofs);
        ofs << "</div>\n";
        ofs << "</div>\n";
        
    } else {

        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('top_s2fids')><a name='summary'>Top abundant S2F ids</a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
        ofs << "<div id='top_s2fids'>\n";
        reportS2fBarPlot(ofs);
        ofs << "</div>\n";
        ofs << "</div>\n";
        
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('top_kos')><a name='summary'>Top abundant KOs</a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
        ofs << "<div id='top_kos'>\n";
        reportKOBarPlot(ofs);
        ofs << "</div>\n";
        ofs << "</div>\n";

        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('reads_quality')><a name='summary'>Reads Quality</a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
        ofs << "<div id='reads_quality' style='display:none'>\n";
        reportReadsQualityPlot3D(ofs);
        ofs << "</div>\n";
        ofs << "</div>\n";
    }
}

void HtmlReporterAll::reportRarefactionKO(ofstream& ofs) {

    ofs << "<div id='ko_hits_figure'>\n";
    ofs << "<div class='figure' id='plot_rarefaction_ko_curve' style='height:400px;'></div>\n";
    ofs << "</div>\n";

    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";
    
    for (Sample & sample : mOptions->samples) {
        
        int top = sample.rarefactionMap.size();
        std::vector<long> x_vec;
        std::vector<double> y_vec;

        for (auto & it : sample.rarefactionMap) {
            x_vec.push_back(it.first);
            y_vec.push_back((double) it.second);
        }
        
        json_str += "{";
        json_str += "x:[" + list2string(x_vec, top) + "],";
        json_str += "y:[" + list2string(y_vec, top) + "],";
        json_str += "name: '" + basename(sample.prefix) + "',";
        json_str += "type:'scatter'";
        json_str += "},";

        x_vec.clear();
        y_vec.clear();
    }
    json_str.pop_back();
    json_str += "];\n";
    json_str += "var layout={title:'Rarefaction curve (KO) ', "
            "xaxis:{title:'Number of reads', automargin: true}, "
            "yaxis:{title:'Number of KOs', automargin: true}};\n";
    json_str += "Plotly.newPlot('plot_rarefaction_ko_curve', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;
}

void HtmlReporterAll::reportRarefactionKO3D(ofstream& ofs) {

    ofs << "<div id='ko3d_hits_figure'>\n";
    ofs << "<div class='figure' id='plot_rarefaction_ko3d_curve' style='height:600px;'></div>\n";
    ofs << "</div>\n";

    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";
    
    for (Sample & sample : mOptions->samples) {
        int top = sample.rarefactionMap.size();
        std::vector<long> x_vec;
        std::vector<double> y_vec;
        std::vector<string> z_vec;

        for (auto & it : sample.rarefactionMap) {
            x_vec.push_back(it.first);
            y_vec.push_back((double) it.second);
            z_vec.push_back(basename(sample.prefix));
        }

        sample.rarefactionMap.clear();
        
        json_str += "{";
        json_str += "x:[" + list2string(z_vec, top) + "],";
        json_str += "y:[" + list2string(x_vec, top) + "],";
        json_str += "z:[" + list2string(y_vec, top) + "],";
        json_str += "name: '" + basename(sample.prefix) + "',";
        json_str += "type: 'scatter3d',";
        json_str += "mode: 'lines',";
        json_str += "line: {width:5},";
        json_str += "showscale: true";
        json_str += "},";

        x_vec.clear();
        y_vec.clear();
        z_vec.clear();
    }
    json_str.pop_back();
    json_str += "];\n";
    json_str += "var layout={title:'Rarefaction curve (KO) 3d', "
            "autosize: 'true',"
            "width: 800, height: 600,"
            "scene: {"
            "xaxis:{title:'Sample', automargin: true}, "
            "yaxis:{title:'Number of reads', automargin: true}, "
            "zaxis:{title:'Number of KOs', automargin: true}}};\n";
    json_str += "Plotly.newPlot('plot_rarefaction_ko3d_curve', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;
}

void HtmlReporterAll::reportKOBarPlot(ofstream& ofs){
    
    ofs << "<div class='subsection_title' onclick=showOrHide('ko_table')><a name='summary'>KOs Table (full list) (click to show/hide) </a></div>\n";
    ofs << "<div id='ko_table' style='overflow:auto; height: 400px;'>\n";
    ofs << "<table class='summary_table'>\n";
    ofs << "<tr><td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << "Sample" << "</td>";
    for(auto & it : smNmVec){
        ofs << "<td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << it << "</td>";
    }
    ofs << "<td class='exlarge' style='font-size:14px;color:#ffffff;background:#008000'>" << "Name" << "</td></tr>\n";

    ofs << "<tr><td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << "Class" << "</td>";
    for (auto & it : mOptions->samples) {
        ofs << "<td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << it.feature << "</td>";
    }
    ofs << "<td class='exlarge' style='font-size:14px;color:#ffffff;background:#008000'>" << "Class_info" << "</td></tr>\n";
    
    for (auto & it : koFreqVec) {
        ofs << "<tr><td class='ko_col'>" << it.at(0) << "</td>";
        if (it.size() > 2) {
            for (int i = 1; i < (it.size() - 1); i++) {
                ofs << "<td class='ko_col'>" << it.at(i) << "</td>";
            }
        }
        ofs << "<td class='exlarge'>" << it.back() << "</td></tr>\n";
    }
    
    ofs << "</table>\n";
    ofs << "</div>\n";
}

void HtmlReporterAll::reportRarefactionS2f(ofstream& ofs) {

    ofs << "<div id='s2f_hits_figure'>\n";
    ofs << "<div class='figure' id='plot_rarefaction_s2f_curve' style='height:400px;'></div>\n";
    ofs << "</div>\n";

    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";
    
    for (Sample & sample : mOptions->samples) {
        
        int top = sample.rarefactionIdMap.size();
        std::vector<long> x_vec;
        std::vector<double> y_vec;

        for (auto & it : sample.rarefactionIdMap) {
            x_vec.push_back(it.first);
            y_vec.push_back((double) it.second);
        }
        
        json_str += "{";
        json_str += "x:[" + list2string(x_vec, top) + "],";
        json_str += "y:[" + list2string(y_vec, top) + "],";
        json_str += "name: '" + basename(sample.prefix) + "',";
        json_str += "type:'scatter'";
        json_str += "},";

        x_vec.clear();
        y_vec.clear();
    }
    json_str.pop_back();
    json_str += "];\n";
    json_str += "var layout={title:'Rarefaction curve (S2F id) ', "
            "xaxis:{title:'Number of reads', automargin: true}, "
            "yaxis:{title:'Number of S2F ids', automargin: true}};\n";
    json_str += "Plotly.newPlot('plot_rarefaction_s2f_curve', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;
}

void HtmlReporterAll::reportRarefactionS2f3D(ofstream& ofs) {

    ofs << "<div id='s2f3d_hits_figure'>\n";
    ofs << "<div class='figure' id='plot_rarefaction_s2f3d_curve' style='height:600px;'></div>\n";
    ofs << "</div>\n";

    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";
    
    for (Sample & sample : mOptions->samples) {
        int top = sample.rarefactionIdMap.size();
        std::vector<long> x_vec;
        std::vector<double> y_vec;
        std::vector<string> z_vec;

        for (auto & it : sample.rarefactionIdMap) {
            x_vec.push_back(it.first);
            y_vec.push_back((double) it.second);
            z_vec.push_back(basename(sample.prefix));
        }

        sample.rarefactionIdMap.clear();
        
        json_str += "{";
        json_str += "x:[" + list2string(z_vec, top) + "],";
        json_str += "y:[" + list2string(x_vec, top) + "],";
        json_str += "z:[" + list2string(y_vec, top) + "],";
        json_str += "name: '" + basename(sample.prefix) + "',";
        json_str += "type: 'scatter3d',";
        json_str += "mode: 'lines',";
        json_str += "line: {width:5},";
        json_str += "showscale: true";
        json_str += "},";

        x_vec.clear();
        y_vec.clear();
        z_vec.clear();
    }
    json_str.pop_back();
    json_str += "];\n";
    json_str += "var layout={title:'Rarefaction curve (S2F ids) 3d', "
            "autosize: 'true',"
            "width: 800, height: 600,"
            "scene: {"
            "xaxis:{title:'Sample', automargin: true}, "
            "yaxis:{title:'Number of reads', automargin: true}, "
            "zaxis:{title:'Number of S2F ids', automargin: true}}};\n";
    json_str += "Plotly.newPlot('plot_rarefaction_s2f3d_curve', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;
}

void HtmlReporterAll::reportS2fBarPlot(ofstream& ofs){
    
    ofs << "<div class='subsection_title' onclick=showOrHide('id_table')><a name='summary'>S2F id Table (full list) </a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
    ofs << "<div id='id_table' style='overflow:auto; height: 400px;'>\n";
    ofs << "<table class='summary_table'>\n";
    ofs << "<tr><td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << "Sample" << "</td>";
    for (auto & it : smNmVec) {
        ofs << "<td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << it << "</td>";
    }

    ofs << "<tr><td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << "Class" << "</td>";
    for (auto & it : mOptions->samples) {
        ofs << "<td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << it.feature << "</td>";
    }

    for (const auto & it : idFreqVec) {
        ofs << "<tr><td class='ko_col'>" << it.at(0) << "</td>";
        if (it.size() > 2) {
            for (int i = 1; i < (it.size() - 1); i++) {
                ofs << "<td class='ko_col'>" << it.at(i) << "</td>";
            }
        }
        ofs << "<td class='ko_col'>" << it.back() << "</td></tr>\n";
    }

    ofs << "</table>\n";
    ofs << "</div>\n";
}

void HtmlReporterAll::reportReadsQualityPlot3D(ofstream& ofs) {
    
    ofs << "<div class='subsection_title' onclick=showOrHide('pre_reads1_quality_3d_curve')><a name='summary'>Prefilter Reads1 (click to show/hide) </a></div>\n";
    ofs << "<div class='figure' id='pre_reads1_quality_3d_curve' style='height:600px;'></div>\n";
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";
    for (Sample & sample : mOptions->samples) {
        auto it = sample.totoalReadsQualityVec.at(0);
        auto total = get<2>(it);
        std::vector<std::string> x_vec(total, basename(sample.prefix));
        json_str += "{";
        json_str += "x:[" + list2string(x_vec, total) + "],";
        json_str += "y:[" + get<0>(it) + "], ";
        json_str += "z:[" + get<1>(it) + "], ";
        json_str += "name: '" + basename(sample.prefix) + "', ";
        json_str += "type: 'scatter3d', ";
        json_str += "mode: 'lines', ";
        json_str += "line: {width:5}, ";
        json_str += "showscale: true";
        json_str += "},";
    }
    
    json_str.pop_back();
    json_str += "];\n";
    json_str += "var layout={title:'Prefilter Reads1 quality', "
            "autosize: 'true',"
            "width: 800, height: 600,"
            "scene: {"
            "xaxis:{title:'Sample', automargin: true}, "
            "yaxis:{title:'Number of cycles', automargin: true}, "
            "zaxis:{title:'Reads quality', automargin: true}}};\n";
    json_str += "Plotly.newPlot('pre_reads1_quality_3d_curve', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;
    
    
    ofs << "<div class='subsection_title' onclick=showOrHide('post_reads1_quality_3d_curve')><a name='summary'>Postfilter Reads1 (click to show/hide) </a></div>\n";
    ofs << "<div class='figure' id='post_reads1_quality_3d_curve' style='height:600px;'></div>\n";
    ofs << "\n<script type=\"text/javascript\">" << endl;
    

    if (mOptions->isPaired()) {
        json_str.clear();
        json_str = "var data=[";
        for (Sample & sample : mOptions->samples) {
            auto it = sample.totoalReadsQualityVec.at(2);
            auto total = get<2>(it);
            std::vector<std::string> x_vec(total, basename(sample.prefix));
            json_str += "{";
            json_str += "x:[" + list2string(x_vec, total) + "],";
            json_str += "y:[" + get<0>(it) + "], ";
            json_str += "z:[" + get<1>(it) + "], ";
            json_str += "name: '" + basename(sample.prefix) + "', ";
            json_str += "type: 'scatter3d', ";
            json_str += "mode: 'lines', ";
            json_str += "line: {width:5}, ";
            json_str += "showscale: true";
            json_str += "},";
        }

        json_str.pop_back();
        json_str += "];\n";
        json_str += "var layout={title:'Postfilter Reads1 quality', "
                "autosize: 'true',"
                "width: 800, height: 600,"
                "scene: {"
                "xaxis:{title:'Sample', automargin: true}, "
                "yaxis:{title:'Number of cycles', automargin: true}, "
                "zaxis:{title:'Reads quality', automargin: true}}};\n";
        json_str += "Plotly.newPlot('post_reads1_quality_3d_curve', data, layout);\n";

        ofs << json_str;
        ofs << "</script>" << endl;
        
    
        ofs << "<div class='subsection_title' onclick=showOrHide('pre_reads2_quality_3d_curve')><a name='summary'>Prefilter Reads2 (click to show/hide) </a></div>\n";
        ofs << "<div class='figure' id='pre_reads2_quality_3d_curve' style='height:600px;'></div>\n";
        ofs << "\n<script type=\"text/javascript\">" << endl;
        
        json_str.clear();
        json_str = "var data=[";
        for (Sample & sample : mOptions->samples) {
            auto it = sample.totoalReadsQualityVec.at(1);
            auto total = get<2>(it);
            std::vector<std::string> x_vec(total, basename(sample.prefix));
            json_str += "{";
            json_str += "x:[" + list2string(x_vec, total) + "],";
            json_str += "y:[" + get<0>(it) + "], ";
            json_str += "z:[" + get<1>(it) + "], ";
            json_str += "name: '" + basename(sample.prefix) + "', ";
            json_str += "type: 'scatter3d', ";
            json_str += "mode: 'lines', ";
            json_str += "line: {width:5}, ";
            json_str += "showscale: true";
            json_str += "},";
        }
        json_str.pop_back();
        json_str += "];\n";
        json_str += "var layout={title:'Prefilter Reads2 quality', "
                "autosize: 'true',"
                "width: 800, height: 600,"
                "scene: {"
                "xaxis:{title:'Sample', automargin: true}, "
                "yaxis:{title:'Number of cycles', automargin: true}, "
                "zaxis:{title:'Reads quality', automargin: true}}};\n";
        json_str += "Plotly.newPlot('pre_reads2_quality_3d_curve', data, layout);\n";

        ofs << json_str;
        ofs << "</script>" << endl;

        ofs << "<div class='subsection_title' onclick=showOrHide('post_reads2_quality_3d_curve')><a name='summary'>Postfilter Reads2 (click to show/hide) </a></div>\n";
        ofs << "<div class='figure' id='post_reads2_quality_3d_curve' style='height:600px;'></div>\n";
        ofs << "\n<script type=\"text/javascript\">" << endl;

        json_str.clear();
        json_str = "var data=[";
        for (Sample & sample : mOptions->samples) {
            auto it = sample.totoalReadsQualityVec.at(3);
            auto total = get<2>(it);
            std::vector<std::string> x_vec(total, basename(sample.prefix));
            json_str += "{";
            json_str += "x:[" + list2string(x_vec, total) + "],";
            json_str += "y:[" + get<0>(it) + "], ";
            json_str += "z:[" + get<1>(it) + "], ";
            json_str += "name: '" + basename(sample.prefix) + "', ";
            json_str += "type: 'scatter3d', ";
            json_str += "mode: 'lines', ";
            json_str += "line: {width:5}, ";
            json_str += "showscale: true";
            json_str += "},";
        }
        json_str.pop_back();
        json_str += "];\n";
        json_str += "var layout={title:'Postfilter Reads2 quality', "
                "autosize: 'true',"
                "width: 800, height: 600,"
                "scene: {"
                "xaxis:{title:'Sample', automargin: true}, "
                "yaxis:{title:'Number of cycles', automargin: true}, "
                "zaxis:{title:'Reads quality', automargin: true}}};\n";
        json_str += "Plotly.newPlot('post_reads2_quality_3d_curve', data, layout);\n";
        ofs << json_str;
        ofs << "</script>" << endl;
    } else {
        json_str.clear();
        json_str = "var data=[";
        for (Sample & sample : mOptions->samples) {
            auto it = sample.totoalReadsQualityVec.at(1);
            auto total = get<2>(it);
            std::vector<std::string> x_vec(total, basename(sample.prefix));
            json_str += "{";
            json_str += "x:[" + list2string(x_vec, total) + "],";
            json_str += "y:[" + get<0>(it) + "], ";
            json_str += "z:[" + get<1>(it) + "], ";
            json_str += "name: '" + basename(sample.prefix) + "', ";
            json_str += "type: 'scatter3d', ";
            json_str += "mode: 'lines', ";
            json_str += "line: {width:5}, ";
            json_str += "showscale: true";
            json_str += "},";
        }

        json_str.pop_back();
        json_str += "];\n";
        json_str += "var layout={title:'Postfilter Reads1 quality', "
                "autosize: 'true',"
                "width: 800, height: 600,"
                "scene: {"
                "xaxis:{title:'Sample', automargin: true}, "
                "yaxis:{title:'Number of cycles', automargin: true}, "
                "zaxis:{title:'Reads quality', automargin: true}}};\n";
        json_str += "Plotly.newPlot('post_reads1_quality_3d_curve', data, layout);\n";

        ofs << json_str;
        ofs << "</script>" << endl;
    }
}

void HtmlReporterAll::reportPathwayBarPlot(ofstream& ofs){
    
    ofs << "<div class='subsection_title' onclick=showOrHide('pathway_table')><a name='summary'> Pathways Table (full list) (click to show/hide) </a></div>\n";
    ofs << "<div id='pathway_table' style='overflow:auto; height: 400px;'>\n";
    ofs << "<table class='summary_table'>\n";
    ofs << "<tr><td class='exlarge' style='font-size:14px;color:#ffffff;background:#008000'>" << "Sample" << "</td>";
    for(auto & it : smNmVec){
        ofs << "<td class='exlarge' style='font-size:14px;color:#ffffff;background:#008000'>" << it << "</td>";
    }
    ofs << "<td class='exlarge' style='font-size:14px;color:#ffffff;background:#008000'>" << "N. of KOs" << "</td></tr>\n";

    ofs << "<tr><td class='exlarge' style='font-size:14px;color:#ffffff;background:#008000'>" << "Class" << "</td>";
    for (auto & it : mOptions->samples) {
        ofs << "<td class='exlarge' style='font-size:14px;color:#ffffff;background:#008000'>" << it.feature << "</td>";
    }
    ofs << "<td class='exlarge' style='font-size:14px;color:#ffffff;background:#008000'>" << "Class_info" << "</td></tr>\n";
    
    for (auto & it : pathwayFreqVec) {
        ofs << "<tr><td class='exlarge'>" << it.at(0) << "</td>";
        if (it.size() > 2) {
            for (int i = 1; i < (it.size() - 1); i++) {
                ofs << "<td class='exlarge'>" << it.at(i) << "</td>";
            }
        }
        ofs << "<td class='exlarge'>" << it.back() << "</td></tr>\n";
    }
    
    ofs << "</table>\n";
    ofs << "</div>\n";
}

void HtmlReporterAll::reportOrgBarPlot(ofstream& ofs){
    
    ofs << "<div class='subsection_title' onclick=showOrHide('org_table')><a name='summary'>Hit Species Table (full list) </a><font color='#88CCFF' > (click to show/hide) </font></div>\n";
    ofs << "<div id='org_table' style='overflow:auto; height: 400px;'>\n";
    ofs << "<table class='summary_table'>\n";
    ofs << "<tr><td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << "Sample" << "</td>";
    for(auto & it : smNmVec){
        ofs << "<td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << it << "</td>";
    }

    ofs << "<tr><td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << "Class" << "</td>";
    for (auto & it : mOptions->samples) {
        ofs << "<td class='ko_col' style='font-size:14px;color:#ffffff;background:#008000'>" << it.feature << "</td>";
    }
   
    for (auto & it : orgFreqVec) {
        ofs << "<tr><td class='ko_col'>" << it.at(0) << "</td>";
        if (it.size() > 2) {
            for (int i = 1; i < (it.size() - 1); i++) {
                ofs << "<td class='ko_col'>" << it.at(i) << "</td>";
            }
        }
        ofs << "<td class='ko_col'>" << it.back() << "</td></tr>\n";
    }
    
    ofs << "</table>\n";
    ofs << "</div>\n";
}

string HtmlReporterAll::list2string(std::vector<long> & x_vec, int top) {
    stringstream ss;
    for(int i = 0; i < top; i++){
        ss << x_vec[i];
        if(i < top - 1){
            ss << ",";
        }
    }
    return ss.str();
}

string HtmlReporterAll::list2string(std::vector<int> & x_vec, int top) {
    stringstream ss;
    for(int i = 0; i < top; i++){
        ss << x_vec[i];
        if(i < top - 1){
            ss << ",";
        }
    }
    return ss.str();
}

string HtmlReporterAll::list2string(std::vector<double> & x_vec, int top) {
    stringstream ss;
    for(int i = 0; i < top; i++){
        ss << x_vec[i];
        if(i < top - 1){
            ss << ",";
        }
    }
    return ss.str();
}

string HtmlReporterAll::list2string(std::vector<string> & x_vec, int top) {
    stringstream ss;
    for(int i = 0; i < top; i++){
        ss << "'";
        ss << x_vec[i];
        if(i < top - 1){
            ss << "',";
        } else {
            ss << "'";
        }
    }
    return ss.str();
}

string HtmlReporterAll::list2string2(std::vector<string> & x_vec, int top) {
    stringstream ss;
    for(int i = 0; i < top; i++){
        ss << x_vec[i];
        if(i < top - 1){
            ss << ",";
        }
    }
    return ss.str();
}

void HtmlReporterAll::printJS(ofstream& ofs){
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
