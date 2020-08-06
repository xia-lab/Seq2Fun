#ifndef FRAGMENT_H
#define FRAGMENT_H

#include <string>

extern "C" {
#include "bwt/bwt.h"
}

using namespace std;

class Fragment {
public:
    std::string seq;
    unsigned int num_mm = 0;
    int diff = 0;
    unsigned int pos_lastmm = 0;
    IndexType si0, si1, arg_si0, arg_si1;
    int matchlen;
    bool SEGchecked = false;

    Fragment(const std::string & s);

    Fragment(const std::string & s, bool b);

    Fragment(const std::string & s, unsigned int n, unsigned int p, int d);

    Fragment(const std::string & s, unsigned int n, unsigned int p, int d, IndexType arg_si0, IndexType arg_si1, int len);

    Fragment(const std::string & s, unsigned int n, unsigned int p, int d, SI * si);

    Fragment(const std::string & s, unsigned int n, unsigned int p, SI * si);

    Fragment(const std::string & s, unsigned int n, unsigned int p);
};

#endif /* FRAGMENT_H */

