#include "fragment.h"

Fragment::Fragment(const std::string & s) {
    seq = s;
}

Fragment::Fragment(const std::string & s, bool b) {
 seq = s;
 SEGchecked = true;
}

Fragment:: Fragment(const std::string & s, unsigned int n, unsigned int p, int d){
 seq = s;
 num_mm = n;
 diff = d;
 pos_lastmm = p;
}

Fragment::Fragment(const std::string & s, unsigned int n, unsigned int p, int d, IndexType arg_si0, IndexType arg_si1, int len){
 seq = s;
 num_mm = n;
 diff = d;
 pos_lastmm = p;
 si0 = arg_si0;
 si1 = arg_si1;
 matchlen = len;
 SEGchecked = true;
} // fragments with substitutions have been checked before

Fragment::Fragment(const std::string & s, unsigned int n, unsigned int p, int d, SI * si){
 seq = s;
 num_mm = n;
 diff = d;
 pos_lastmm = p;
 si0 = si->start;
 si1 = si->start + (IndexType) si->len;
 matchlen = si->ql;
}

Fragment::Fragment(const std::string & s, unsigned int n, unsigned int p, SI * si){
 seq = s;
 num_mm = n;
 pos_lastmm = p;
 si0 = si->start;
 si1 = si->start + (IndexType) si->len;
 matchlen = si->ql; 
}

Fragment::Fragment(const std::string & s, unsigned int n, unsigned int p) {
seq = s;
num_mm = n;
pos_lastmm = p; 
}