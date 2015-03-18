//
//    MIC-PEAlgnmt.h
//
//    SOAP2 / 2BWT
//
//    Copyright (C) 2013, HKU
//
//    This program is free software; you can redistribute it and/or
//    modify it under the terms of the GNU General Public License
//    as published by the Free Software Foundation; either version 2
//    of the License, or (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __MIC_PE_ALIGNMENT_H__
#define __MIC_PE_ALIGNMENT_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "MIC-SRAArguments.h"


#define MIC_PE_MAX_RESULT                       10000
#define MIC_PAIRING_MERGE_LIMIT                 5

#define MIC_PE_OUTPUT_STATUS_OPEN               0
#define MIC_PE_OUTPUT_STATUS_SKIPPED            1
#define MIC_PE_OUTPUT_STATUS_UNHANDLE           2
#define MIC_PE_OUTPUT_STATUS_CLOSED             3
#define MIC_PE_OUTPUT_STATUS_PAIR               4
#define MIC_PE_OUTPUT_STATUS_BAD_PAIR           5
#define MIC_PE_OUTPUT_STATUS_READ_NO_ALIGNMENT  6
#define MIC_PE_OUTPUT_STATUS_MATE_NO_ALIGNMENT  7
#define MIC_PE_OUTPUT_STATUS_BOTH_NO_ALIGNMENT  8
#define MIC_PE_OUTPUT_STATUS_DP_BASE_MATE       9
#define MIC_PE_OUTPUT_STATUS_DP_BASE_READ       10
#define MIC_PE_OUTPUT_STATUS_DP_MIX_BASE        11
#define MIC_PE_OUTPUT_STATUS_DP_SEED            12
//<---- Insert new status here and increase the count
//      (value of MIC_PE_OUTPUT_STATUS_COUNT) below 
//      ; also add the string value to MIC_PE_OUTPUT_STATUS
#define MIC_PE_OUTPUT_STATUS_COUNT              13

__attribute__((target(mic)))
static const char MIC_PE_OUTPUT_STATUS[][128] = {
    "OPEN",
    "SKIPPED",
    "UNHANDLE",
    "CLOSED",
    "PAIR",
    "BAD_PAIR",
    "READ_NO_ALIGNMENT",
    "MATE_NO_ALIGNMENT",
    "BOTH_NO_ALIGNMENT",
    "DP_BASE_MATE",
    "DP_BASE_READ",
    "DP_MIX_BASE",
    "DP_SEED"};

typedef struct MICPEArguments {
    // Input
    MICSRAArguments * readArgs;
    MICSRAArguments * mateArgs;
    unsigned int outputVacancy;
    unsigned int metaVacancy;
    int outputLimit;

    // Upper and lower bound for insert
    unsigned int uBound;
    unsigned int lBound;

    // Merge nearby position, 0 = disable 1 = enable
    uint8_t mergeEnable;

    // Output
    unsigned int * output;
    uint16_t * occCount;
    uint8_t * outputStatus;
    MICSRAOccMetadata * outputMeta;
    
    // DP Arguments
    DPScores * dpScores;

} MICPEArguments;

__attribute__((target(mic)))
void MICPEArgumentsConfig(
        MICPEArguments * peArguments,
        MICSRAArguments * readArgs,
        MICSRAArguments * mateArgs,
        unsigned int outputVacancy,
        unsigned int metaVacancy,
        unsigned int * output,
        MICSRAOccMetadata * outputMeta,
        uint16_t * occCount,
        uint8_t * outputStatus);

__attribute__((target(mic)))
void MICPEMappingInitialise(MICPEArguments * peArguments);

__attribute__((target(mic)))
void MICPEMappingOccurrences(MICPEArguments * peArguments);

__attribute__((target(mic)))
void MICPEMappingComplete(MICPEArguments * peArguments);

__attribute__((target(mic)))
void MICPEArgumentsSetBounds(MICPEArguments * peArguments,
        unsigned int lowerBound,
        unsigned int upperBound);

__attribute__((target(mic)))
void MICPEArgumentsSetDPScores(MICPEArguments * peArguments,
        DPScores * dpScores);


// Set a limit to output, occurrences after the limit will be dropped
__attribute__((target(mic)))
void MICPEArgumentsSetMaxOutput(MICPEArguments * peArguments,
        int outputLimit);

// Set a limit to output, occurrences after the limit will be dropped
__attribute__((target(mic)))
void MICPEArgumentsSetMergeEnable(MICPEArguments * peArguments,
        uint8_t mergeEnable);
#endif

