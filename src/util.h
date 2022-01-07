#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <algorithm>
#include <time.h>
#include <mutex>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <stdio.h>  /* defines FILENAME_MAX */
//#include <limits.h>
//#include <unistd.h>
//#include <boost/algorithm/string.hpp>
#include <tuple>
#include <utility>
#include <numeric>
#include <string>  
#include <sstream>
#include <fstream>
#include <deque>
#include <regex>

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#include <limits.h>
#define GetCurrentDir getcwd
#endif

#ifdef WINDOWS
#include <windows.h>
#define GetExePath _getcwd
#else
#include <unistd.h>
#include <limits.h>
#define GetExePath getcwd
#endif

#include "options.h"


using namespace std;

template <typename T>
void colorCout(const T & str, char color = 'r') {
    if (color == 'r') {
        std::cout << "\033[1;31m" << str << "\033[0m\n";
    } else if (color == 'g') {
        std::cout << "\033[1;32m" << str << "\033[0m\n";
    } else if (color == 'y') {
        std::cout << "\033[1;33m" << str << "\033[0m\n";
    } else if (color == 'b') {
        std::cout << "\033[1;34m" << str << "\033[0m\n";
    } else if (color == 'm') {
        std::cout << "\033[1;35m" << str << "\033[0m\n";
    } else if (color == 'c') {
        std::cout << "\033[1;36m" << str << "\033[0m\n";
    } else if (color == 'w') {
        std::cout << "\033[1;37m" << str << "\033[0m\n";
    } else {
        std::cout << "\033[1;30m" << str << "\033[0m\n";
    }
}

inline char complement(char base) {
    switch (base) {
        case 'A':
        case 'a':
            return 'T';
        case 'T':
        case 't':
            return 'A';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        default:
            return 'N';
    }
}

inline bool starts_with(string const & value, string const & starting) {
    if (starting.size() > value.size()) return false;
    return equal(starting.begin(), starting.end(), value.begin());
}

inline bool ends_with(string const & value, string const & ending) {
    if (ending.size() > value.size()) return false;
    return equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline string trim(const string& str) {
    string::size_type pos = str.find_first_not_of(' ');
    if (pos == string::npos) {
        return string("");
    }
    string::size_type pos2 = str.find_last_not_of(' ');
    if (pos2 != string::npos) {
        return str.substr(pos, pos2 - pos + 1);
    }
    return str.substr(pos);
}

inline string trimStr(const string& str) {
    string::size_type pos = str.find_first_not_of(' ');
    if (pos == string::npos) {
        return string("");
    }
    string::size_type pos2 = str.find_last_not_of(' ');
    if (pos2 != string::npos) {
        return str.substr(pos, pos2 - pos + 1);
    }
    return str.substr(pos);
}

inline string trimEnd(char line[]) {
    int readed = strlen(line);
    if (readed >= 2) {
        if (line[readed - 1] == '\n' || line[readed - 1] == '\r') {
            line[readed - 1] = '\0';
            if (line[readed - 2] == '\r') {
                line[readed - 2] = '\0';
            }
        }
    }
    string str(line);
    return str;
}

inline int splitStr(const string& str, vector<string>& ret_, string sep = "\t") {
    if (str.empty()) {
        return 0;
    }
    string tmp;
    string::size_type pos_begin = str.find_first_not_of(sep);
    string::size_type comma_pos = 0;
    while (pos_begin != string::npos) {
        comma_pos = str.find(sep, pos_begin);
        if (comma_pos != string::npos) {
            tmp = str.substr(pos_begin, comma_pos - pos_begin);
            pos_begin = comma_pos + sep.length();
        } else {
            tmp = str.substr(pos_begin);
            pos_begin = comma_pos;
        }
        ret_.push_back(tmp);
        tmp.clear();
    }
    ret_.shrink_to_fit();
    return 0;
}

inline std::vector<std::string> split2(const std::string & s, const char delim) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delim)) {
        tokens.push_back(token);
    }
    return tokens;
}

inline int countFreq(std::string & s, std::string p){
    int M = s.length();
    int N = p.length();
    int res = 0;
    for(int i = 0; i <= M - N; i ++){
        int j;
        for(j = 0; j < N; j++){
            if(s[i+j] != p[j]) break;
        }
        if(j == N){
            res++;
            j = 0;
        }
    }
    return res;
}

inline string replace(const string& str, const string& src, const string& dest) {
    string ret;

    string::size_type pos_begin = 0;
    string::size_type pos = str.find(src);
    while (pos != string::npos) {
        ret.append(str.data() + pos_begin, pos - pos_begin);
        ret += dest;
        pos_begin = pos + 1;
        pos = str.find(src, pos_begin);
    }
    if (pos_begin < str.length()) {
        ret.append(str.begin() + pos_begin, str.end());
    }
    return ret;
}

inline string reverse(const string& str) {
    string ret(str.length(), 0);
    for (int pos = 0; pos < str.length(); pos++) {
        ret[pos] = str[str.length() - pos - 1];
    }
    return ret;
}

inline string basename(const string& filename) {
    string::size_type pos = filename.find_last_of('/');
    if (pos == string::npos)
        return filename;
    else if (pos == filename.length() - 1)
        return ""; // a bad filename
    else
        return filename.substr(pos + 1, filename.length() - pos - 1);
}

inline string get_current_dir() {
    char buff[FILENAME_MAX]; //create string buffer to hold path
    GetCurrentDir(buff, FILENAME_MAX);
    string current_working_dir(buff);
    return current_working_dir;
}

inline string get_upper_dir() {
    std::string cwd = get_current_dir();
    std::string::size_type bepos = cwd.find_last_of("/");
    cwd.erase(bepos);
    return (cwd);
}

#ifdef WINDOWS

std::string getexepath() {
    char result[ MAX_PATH ];
    return std::string(result, GetModuleFileName(NULL, result, MAX_PATH));
}
#else

inline string GetExePath() {
    char result[ PATH_MAX ];
    ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
    return std::string(result, (count > 0) ? count : 0);
}

inline string get_seq2fun_dir() {
    std::string cwd = GetExePath();
    std::string::size_type bepos = cwd.find("/bin");
    cwd.erase(bepos);
    return (cwd);
}
#endif

inline string dirname(const string& filename) {
    string::size_type pos = filename.find_last_of('/');
    if (pos == string::npos) {
        return "./";
    } else {
        return filename.substr(0, pos + 1);
    }
}

inline string joinpath(const string& dirname, const string& basename) {
    if (dirname[dirname.length() - 1] == '/') {
        return dirname + basename;
    } else {
        return dirname + "/" + basename;
    }
}

inline string checkDirEnd(const string & dirname) {
    return (dirname[dirname.length() - 1] == '/' ? dirname : dirname + '/');
}

//Check if a string is a file or directory

inline bool file_exists(const string& s) {
    bool exists = false;
    if (s.length() > 0) {
        struct stat status;
        int result = stat(s.c_str(), &status);
        if (result == 0) {
            exists = true;
        }
    }
    return exists;
}


// check if a string is a directory

inline bool is_directory(const string& path) {
    bool isdir = false;
    struct stat status;
    // visual studion use _S_IFDIR instead of S_IFDIR
    // http://msdn.microsoft.com/en-us/library/14h5k7ff.aspx
#ifdef _MSC_VER
#define S_IFDIR _S_IFDIR
#endif
    stat(path.c_str(), &status);
    if (status.st_mode & S_IFDIR) {
        isdir = true;
    }
    // #endif
    return isdir;
}

inline void check_file_valid(const string& s) {
    if (!file_exists(s)) {
        cerr << "ERROR: file '" << s << "' doesn't exist, quit now" << endl;
        exit(-1);
    }
    if (is_directory(s)) {
        cerr << "ERROR: '" << s << "' is a folder, not a file, quit now" << endl;
        exit(-1);
    }
}

inline bool check_filename_valid(const string& s) {
    return 0 < trim(s).length() && trim(s).length() <= 255 && regex_match(s, regex("^[A-Za-z0-9_\\.\\-]+$"));
}

inline void check_file_writable(const string& s) {
    string dir = dirname(s);
    if (!file_exists(dir)) {
        cerr << "ERROR: '" << dir << " doesn't exist. Create this folder and run this command again." << endl;
        exit(-1);
    }
    if (is_directory(s)) {
        cerr << "ERROR: '" << s << "' is not a writable file, quit now" << endl;
        exit(-1);
    }
}

// Remove non alphabetic characters from a string

inline string str_keep_alpha(const string& s) {
    string new_str;
    for (size_t it = 0; it < s.size(); it++) {
        if (isalpha(s[it])) {
            new_str += s[it];
        }
    }
    return new_str;
}


// Remove invalid sequence characters from a string

inline void str_keep_valid_sequence(string& s, bool forceUpperCase = false) {
    size_t total = 0;
    const char case_gap = 'a' - 'A';
    for (size_t it = 0; it < s.size(); it++) {
        char c = s[it];
        if (forceUpperCase && c >= 'a' && c <= 'z') {
            c -= case_gap;
        }
        if (isalpha(c) || c == '-' || c == '*') {
            s[total] = c;
            total++;
        }
    }

    s.resize(total);
}

//inline int convert_float2int(int & a, float & b){
//    return((int)round(a * b));
//}

inline int find_with_right_pos(const string& str, const string& pattern, int start = 0) {
    int pos = str.find(pattern, start);
    if (pos < 0)
        return -1;
    else
        return pos + pattern.length();
}

inline void str2upper(string& s) {
    transform(s.begin(), s.end(), s.begin(), (int (*)(int))toupper);
}

inline void str2lower(string& s) {
    transform(s.begin(), s.end(), s.begin(), (int (*)(int))tolower);
}

inline char num2qual(int num) {
    if (num > 127 - 33)
        num = 127 - 33;
    if (num < 0)
        num = 0;

    char c = num + 33;
    return c;
}

inline void error_exit(const string& msg) {
    cerr << "ERROR: " << msg << endl;
    exit(-1);
}

extern mutex logmtx;

inline void loginfo(const string & s, bool next = true) {
    logmtx.lock();
    time_t tt = time(NULL);
    tm* t = localtime(&tt);
    if (next) {
        fprintf(stderr, "[\033[1;35m%02d:%02d:%02d\033[0m] %s\n", t->tm_hour, t->tm_min, t->tm_sec, s.c_str());
    } else {
        fprintf(stderr, "[\033[1;35m%02d:%02d:%02d\033[0m] %s\r", t->tm_hour, t->tm_min, t->tm_sec, s.c_str());
    }
    logmtx.unlock();
}

inline void loginfolong(const string & s){
    logmtx.lock();
    time_t tt = time(NULL);
    tm* t = localtime(&tt);
    fprintf(stderr, "[%02d:%02d:%02d] %s\n", t->tm_hour, t->tm_min, t->tm_sec, s.c_str());
    logmtx.unlock();
}

inline string trimName(string& str) {
    // string strnew;
    str.erase(str.begin());
    string suffixStartCharacters = " /\t\r";
    size_t n = str.find_first_of(suffixStartCharacters);
    if (n != string::npos) {
        return str.erase(n);
    } else {
        return str;
    }
}

template<typename T>
T getMostFreqStrFromVec(std::vector< T > & vectorko) {
    std::map<T, int> freq;

    for (int i = 0; i < vectorko.size(); i++) {
        freq[vectorko[i]]++;
    }

    int max_F = 0;
    T res;
    for (const auto & j : freq) {
        if (max_F < j.second) {
            res = j.first;
            max_F = j.second;
        }
    }
    freq.clear();
    vectorko.clear();
    return res;
}

inline std::string getUniqStrFromVec(std::vector<std::string> & vectorko) {
    sort(vectorko.begin(), vectorko.end());
    vectorko.erase(unique(vectorko.begin(), vectorko.end()), vectorko.end());

    std::string uniq = "";
    if (vectorko.size() == 1) {
        uniq = vectorko[1];
    }
    vectorko.clear();
    return uniq;
}

inline string removeStr(const string &str, const string &src) {
    string ret;
    string::size_type pos_begin = 0;
    string::size_type pos = str.find(src);
    while (pos != string::npos) {
        ret.append(str.data() + pos_begin, pos - pos_begin);
        pos_begin = pos + src.length();
        pos = str.find(src, pos_begin);
    }
    if (pos_begin < str.length()) {
        ret.append(str.begin() + pos_begin, str.end());
    }
    return ret;
}

inline string removeStrs(const string &str) {
    string ret;
    string::size_type pos_begin = 0;
    string::size_type pos = str.find("\tK");
    if (pos != string::npos) {
        ret.append(str.data(), pos);
    } else {
        pos = str.find("\tU");
        if (pos != string::npos) {
            ret.append(str.data(), pos);
        } else {
            return str;
        }
    }
    return ret;
}

template<class T1, class T2>
double getPercentage(T1 v1, T2 v2) {
    double ret = double (v1 * 100) / (double) v2;
    return ret;
}

//inline double getPercetage(int v1, long v2){
//   double ret = double (v1 * 100) / (double) v2;
//   return ret;
//}
//
//inline double getPercetageInt(int v1, int v2){
//   double ret = double (v1 * 100) / (double) v2;
//   return ret;
//}

template<typename T>
std::string unkown2Str(const T & t) {
    std::ostringstream os;
    os << t;
    return os.str();
}


//bool cmp(std::pair<const uint32*, uint32>& l,
//        std::pair<const uint32*, uint32>& r){
//    return l.second > r.second;
//}
//
//std::vector<std::pair<const uint32*, uint32> > sortMapToVector2(std::map<const uint32*, uint32> &unMap) {
//    std::vector<std::pair < const uint32*, uint32 > > sortedVec;
//    sortedVec.reserve(unMap.size());
//    
//    for(const auto & it : unMap){
//        sortedVec.push_back(it);
//    }
//    std::sort(sortedVec.begin(), sortedVec.end(), cmp);
//    
//    for(const auto & it : sortedVec){
//        std::cout << it.first << " : " << it.second << "\n";
//    }
//    return sortedVec;
//}

template <typename T1, typename T2>
std::vector<std::pair<const T1*, T2> > sortMapToVector(std::map<const T1*, T2> &unMap) {
    std::vector<std::pair < const T1*, T2 >> sortedVec;
    sortedVec.reserve(unMap.size());
    std::partial_sort_copy(unMap.begin(),
            unMap.end(),
            sortedVec.begin(),
            sortedVec.end(),
            [](std::pair<const T1*, float> const &l,
            std::pair<const T1*, float> const &r) {
                return l.second > r.second;
            });
    return sortedVec;
}

template <typename T1, typename T2>
std::vector<std::pair<T1, T2> > sortUMapToVector(std::map<T1, T2> &unMap) {
    int size = unMap.size();
    std::vector<std::pair < T1, T2 >> sortedVec(size);
    std::partial_sort_copy(unMap.begin(),
            unMap.end(),
            sortedVec.begin(),
            sortedVec.end(),
            [](std::pair<const T1, float> const &l,
            std::pair<const T1, float> const &r) {
                return l.second > r.second;
            });
    return sortedVec;
}

template <typename T1, typename T2>
std::vector<std::pair<T1, T2> > sortUMapToVector(std::unordered_map<T1, T2> &unMap) {
    int size = unMap.size();
    std::vector<std::pair < T1, T2 >> sortedVec(size);
    std::partial_sort_copy(unMap.begin(),
            unMap.end(),
            sortedVec.begin(),
            sortedVec.end(),
            [](std::pair<const T1, float> const &l,
            std::pair<const T1, float> const &r) {
                return l.second > r.second;
            });
    return sortedVec;
}

//template <const typename T1, typename T2>
//std::vector<std::pair<T1, T2> > sortUMapToVector(std::map<T1, T2> &unMap) {
//    int size = unMap.size();
//    std::vector<std::pair <T1, T2 >> sortedVec(size);
//    std::partial_sort_copy(unMap.begin(),
//            unMap.end(),
//            sortedVec.begin(),
//            sortedVec.end(),
//            [](std::pair<T1, float> const &l,
//            std::pair<T1, float> const &r) {
//                return *(l.second) > *(r.second);
//            });
//    return sortedVec;
//}

inline std::vector<std::tuple<std::string, double, int, int> > sortTupleVector(std::vector<std::tuple<std::string, double, int, int> > & OriVec) {
    int size = OriVec.size();
    std::vector<std::tuple<std::string, double, int, int> > sortedVec(size);
    std::partial_sort_copy(OriVec.begin(),
            OriVec.end(),
            sortedVec.begin(),
            sortedVec.end(),
            [](const std::tuple<std::string, double, int, int > &l,
            const std::tuple<std::string, double, int, int > &r) {
                return get<1>(l) > get<1>(r);
            });
    return sortedVec;
}

inline void getUniqVec(std::vector<std::string> &orgVec) {
    std::unordered_set<std::string> s;
    std::vector<std::string>::iterator itr = orgVec.begin();
    for (auto curr = orgVec.begin(); curr != orgVec.end(); ++curr) {
        if (s.insert(*curr).second)
            *itr++ = *curr;
    }
    orgVec.erase(itr, orgVec.end());
}

inline void stripChar(std::string & s) {
    for (auto it = s.begin(); it != s.end(); ++it) {
        if (!isalpha(*it)) {
            s.erase(it);
            it--;
        }
    }
}

template<typename T>
inline std::string convertSeconds(const T & t) {
    long seconds, minutes, hours;
    seconds = long(t);
    minutes = t / 60;
    hours = minutes / 60;

    std::stringstream ss;
    ss << hours << " hours " << minutes % 60 << " minutes " << seconds % 60 << " seconds";
    return ss.str();
}

template<typename T>
vector<T> sliceVec(vector<T> vec, int x, int y) {
    auto start = vec.begin() + x;
    auto end = vec.begin() + y;
    vector<T> res(y - x);
    copy(start, end, res.begin());
    return res;
}

template<typename T>
int getVecIndex(vector<T> & v, T i) {
    auto it = std::find(v.begin(), v.end(), i);
    return (it != v.end() ? it - v.begin() : -1);
}

#endif /* UTIL_H */
