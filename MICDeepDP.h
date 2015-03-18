#ifndef __MIC_DEEP_DP_H__
#define __MIC_DEEP_DP_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>

#include "MIC-SRAAlgnmt.h"
#include "MIC-PEAlgnmt.h"
#include "MIC-DPAlgnmt.h"

#define MIC_DEEP_DP_SEED_PAIR_LIMIT 2

// -----------------------------
// Debug Flag
// -----------------------------
// Define to print output message on stdout
// Not advisable for data set larger than 1
//#define MIC_DP_DEBUG_PRINT_SEED_MATCHING

__attribute__((target(mic)))
int MICDeepDP(
        MICPEArguments * peSeedArgs, 
        unsigned int * peOutputBuffer,
        MICSRAOccMetadata * peMetaBuffer, 
        MICSRAArguments * readArgs,
        MICSRAArguments * mateArgs,
        MICSRAArguments * readSeedSRAArgs,
        MICSRAArguments * mateSeedSRAArgs,
        CPTSRAModel * cpPModels,
        CPTSRAModel * cpNModels,
        HSP * hsp,
        DPWorkMIC * dpWork,
        double seedLength,
        double seedOverlap,
        double orphanTrigger,
        double dpScoreTF,
        int seedOccLimit,
        MICDPOccurrence * dpOcc);


#endif
