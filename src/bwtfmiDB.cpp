#include "bwtfmiDB.h"

BwtFmiDB::BwtFmiDB(Options * & opt) {
    mOptions = opt;
    init();
}

BwtFmiDB::~BwtFmiDB() {
    //for trans search
    if (tfmi) {
        delete tfmi;
        tfmi = NULL;
    }

    if (tbwt) {
        delete tbwt;
        tbwt = NULL;
    }

    if (tastruct->trans) free(tastruct->trans);
    if (tastruct->a) free(tastruct->a);
    if (tastruct) free(tastruct);
    if (mOptions->transSearch.SEG) {
        SegParametersFree(tblast_seg_params);
    }
}

void BwtFmiDB::init() {
    if (!mOptions->transSearch.tfmi.empty()) {
        if (mOptions->verbose) {
            std::string msg = "Reading protein (trans search) BWT FMI index from file " + mOptions->transSearch.tfmi;
            loginfo(msg);
        }

        FILE * tfile = fopen(mOptions->transSearch.tfmi.c_str(), "r");
        tbwt = readIndexes(tfile);
        Transsearch = true;
        fclose(tfile);
        tfmi = tbwt->f;
        //if (mOptions->verbose) {
            std::stringstream msgs;
            msgs << "Protein (trans search) BWT of length " << tbwt->len << " has been read with " << tbwt->nseq << " sequences, alphabet = " << tbwt->alphabet;
            loginfo(msgs.str());
        //}

        tdb_length = (double) (tbwt->len - tbwt->nseq);
        if (mOptions->verbose) {
        std::string msg = "Protein (trans search) double length is " + to_string(tdb_length);
            loginfo(msg);
        }

        tastruct = alloc_AlphabetStruct(tbwt->alphabet, 0, 0);

        //need to be conformed.
        if (mOptions->transSearch.SEG) {
            tblast_seg_params = SegParametersNewAa(); //need to be conformed;
            tblast_seg_params->overlaps = TRUE;
        }
    }

    if (mOptions->verbose) loginfo("finish BwtFmiDB initiation");
}

