#include "bwtfmiDB.h"

BwtFmiDB::BwtFmiDB(Options * & opt) {
    mOptions = opt;
    init();
}

BwtFmiDB::~BwtFmiDB() {
    free_BWT();
    if (tastruct->trans) free(tastruct->trans);
    if (tastruct->a) free(tastruct->a);
    if (tastruct) free(tastruct);
    if (mOptions->transSearch.SEG) {
        SegParametersFree(tblast_seg_params);
    }
}


void BwtFmiDB::free_BWT() {
    if (tbwt == NULL)
        return; // Check if bwt is NULL
    if (tbwt->f != NULL) {
        free_FMI(tbwt->f);
        tbwt->f = NULL;
    }
    if (tbwt->s != NULL) {
        free_suffixArray(tbwt->s);
        tbwt->s = NULL;
    }
    // Free dynamically allocated members
    if (tbwt->bwt != NULL) {
        free(tbwt->bwt);
        tbwt->bwt = NULL;
    }
    if (tbwt->alphabet != NULL) {
        free(tbwt->alphabet);
        tbwt->alphabet = NULL;
    }
    // Finally, free the BWT structure itself
    free(tbwt);
}

void BwtFmiDB::free_FMI(FMI*& fmi) {
    if (fmi == NULL) return; // Check if fmi is NULL

    if (fmi->index1 != NULL) {
        for (int i = 0; i < fmi->N1; i++) {
            if (fmi->index1[i] != NULL) {
                free(fmi->index1[i]);
                fmi->index1[i] = NULL;
            }
        }
        free(fmi->index1);
        fmi->index1 = NULL;
    }

    if (fmi->index2 != NULL) {
        for (int i = 0; i < fmi->N2; i++) {
            if (fmi->index2[i] != NULL) {
                free(fmi->index2[i]);
                fmi->index2[i] = NULL;
            }
        }
        free(fmi->index2);
        fmi->index2 = NULL;
    }

    if (fmi->startLcode != NULL) {
        free(fmi->startLcode);
        fmi->startLcode = NULL;
    }
    
        // Free dynamically allocated members
    if (fmi->bwt != NULL) {
        free(fmi->bwt);
        fmi->bwt = NULL;
    }
    // Finally, free the FMI structure itself
    free(fmi);
}

void BwtFmiDB::free_suffixArray(suffixArray*& sa) {
    if (sa == NULL) return; // Check if sa is NULL
    // Free dynamically allocated members
    if (sa->sa != NULL) {
        free(sa->sa);
        sa->sa = NULL;
    }
    if (sa->seqTermOrder != NULL) {
        free(sa->seqTermOrder);
        sa->seqTermOrder = NULL;
    }
    if (sa->seqlengths != NULL) {
        free(sa->seqlengths);
        sa->seqlengths = NULL;
    }
    if (sa->hash != NULL) {
        for (int i = 0; i < sa->nseq; ++i) {
            SEQstruct *cur = sa->hash[i];
            recursive_free_SEQstruct(cur);
        }
        free(sa->hash);
        sa->hash = NULL;
    }
    if (sa->ids != NULL) {
        for (int i = 0; i < sa->nseq; ++i) {
            free(sa->ids[i]);
        }
        free(sa->ids);
        sa->ids = NULL;
    }
    //if (sa->seqstart != NULL) {
      //  free(sa->seqstart);
       // sa->seqstart = NULL;
   // }
    // Finally, free the suffixArray structure itself
    free(sa);
    sa = NULL;
}

void BwtFmiDB::init() {
    if (!mOptions->transSearch.tfmi.empty()) {
        if (mOptions->verbose) {
            std::string msg = "Reading protein (trans search) BWT FMI index from file " + mOptions->transSearch.tfmi;
            mOptions->longlog ? loginfolong(msg) : loginfo(msg);
        }

        FILE * tfile = fopen(mOptions->transSearch.tfmi.c_str(), "r");
        tbwt = readIndexes(tfile);
        Transsearch = true;
        fclose(tfile);
        tfmi = tbwt->f;
        //if (mOptions->verbose) {
            std::stringstream msgs;
            msgs << "Protein (trans search) BWT of length " << tbwt->len << " has been read with " << tbwt->nseq << " sequences, alphabet = " << tbwt->alphabet;
            mOptions->longlog ? loginfolong(msgs.str()) : loginfo(msgs.str());
        //}

        tdb_length = (double) (tbwt->len - tbwt->nseq);
        if (mOptions->verbose) {
        std::string msg = "Protein (trans search) double length is " + to_string(tdb_length);
            mOptions->longlog ? loginfolong(msg) : loginfo(msg);
        }

        tastruct = alloc_AlphabetStruct(tbwt->alphabet, 0, 0);

        //need to be conformed.
        if (mOptions->transSearch.SEG) {
            tblast_seg_params = SegParametersNewAa(); //need to be conformed;
            tblast_seg_params->overlaps = TRUE;
        }
    }

    if (mOptions->verbose) {
        mOptions->longlog ? loginfolong("finish BwtFmiDB initiation") : loginfo("finish BwtFmiDB initiation");
    }
}

