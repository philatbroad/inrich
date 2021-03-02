#ifndef COMMON_H
#define COMMON_H

#include <time.h>
#include <string>

const std::string inrich_ver = "INRICH v.1.0";


enum test_t { TARGET = 0 ,
              INTERVAL = 1 } ;

enum enrich_test_mode {
    INRICH = 1 ,
    ALIGATOR = 2
} ;


// =======================================
// Physical Clustering Test
//
// PL 2010.08.18 ADDED
// =======================================
const static double SIG_STUDY_HIT[] = { 0.001, 0.01, 0.05 };
const static unsigned int N_SIG_STUDY_HIT = 3;

const static int  MAX_PAIR_DIST = 999999999;   // 32bit Maximum Integer: 4,294,967,295
const static bool CAL_DIST_IF_SAME_CHR = false;
const static bool CAL_DIST_IF_ADJACENT = false;

// =========================================
// PL 2010.08.30 ADDED
// =========================================



//void sleep(unsigned int mseconds);




#endif // COMMON_H
