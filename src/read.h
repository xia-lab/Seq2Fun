#ifndef READ_H
#define READ_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include "sequence.h"
#include <vector>
#include "util.h"

using namespace std;

class Read{
public:
    Read(string name, string seq, string strand, string quality, bool phred64 = false);
    Read(string name, Sequence seq, string strand, string quality, bool phred64 = false);
    Read(string name, string seq, string strand);
    Read(string name, Sequence seq, string strand);
    Read(Read &r);
    void print();
    void printFile(ofstream& file);
    Read* reverseComplement();
    string firstIndex();
    string lastIndex();
    // default is Q20
    int lowQualCount(int qual = 20);
    int length();
    string toString();
    string toFastaR1();
    string toFastaR2();
    string toStringWithTag(string tag);
    string toStringWithTag(uint32* tag);
    string toStringWithTagRm();
    void resize(int len);
    void convertPhred64To33();
    void trimFront(int len);
    bool fixMGI();

public:
    static bool test();

private:


public:
	string mName;
	Sequence mSeq;
	string mStrand;
	string mQuality;
	bool mHasQuality;
};

class ReadPair{
public:
    ReadPair(Read* left, Read* right);
    ~ReadPair();

    // merge a pair, without consideration of seq error caused false INDEL
    Read* fastMerge();
public:
    Read* mLeft;
    Read* mRight;

public:
    static bool test();
};

class ReadItem{
public:
    std::string name1;
    std::string name2;
    std::string sequence1;
    std::string quality1;
    std::string sequence2;
    std::string quality2;
    bool paired = false;
    ReadItem(const std::string &, const std::string &);
    ReadItem(const std::string &, const std::string &, const std::string &);
//    ReadItem(const std::string & n1, const std::string & s1, const std::string & q1);
//    ReadItem(const std::string & n1, const std::string & s1, const std::string & q1, const std::string & n2, const std::string & s2, const std::string & q2);
//    std::string toStringR1();
    std::string toStringWithTagR1(std::string & tag);
//    std::string toStringR2();
    std::string toStringWithTagR2(std::string & tag);

//    std::string toStringRQ1();
//    std::string toStringWithTagRQ1(std::string & tag);
//    std::string toStringRQ2();
//    std::string toStringWithTagRQ2(std::string & tag);
//    
};

#endif