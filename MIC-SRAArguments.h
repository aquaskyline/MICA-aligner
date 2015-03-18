//
//    MIC-SRAArguments.h
//
//    mica
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

#ifndef __MIC_SRAARGUMENTS_H__
#define __MIC_SRAARGUMENTS_H__

#include "2bwt-flex/SRAArguments.h"

#define MIC_OUTPUT_STATUS_OPEN      0
#define MIC_OUTPUT_STATUS_CLOSE     1
#define MIC_OUTPUT_STATUS_UNHANDLE  2
#define MIC_OUTPUT_STATUS_SKIPPED   3
#define MIC_OUTPUT_STATUS_COMPLETE  4

// The upper bound of the number of alignment per read
#define MIC_SRA_OUTPUT_MAX_ALIGNMENT   10000
#define MIC_SRA_OUTPUT_SIZE_PER_THREAD 100000

#define MIC_SRA_OUTPUT_MAX_META        10000
#define MIC_SRA_OUTPUT_META_PER_THREAD 100000

#define MIC_SRA_MAX_NUM_OF_ERROR    3

typedef struct MICSRAOccMetadata {

    // POS_STRAND:0, NEG_STRAND:1
    uint32_t numPayload:24, strand:1, numOfErr:7;
    SRAError errors[MIC_SRA_MAX_NUM_OF_ERROR];
    uint16_t matchLen;
} MICSRAOccMetadata;

typedef struct MICSRAArguments {
    unsigned char * readCode;
    unsigned char * readCode_Complt;
    uint16_t seedLength;
    uint16_t readLength;
    uint8_t readStrand;
    uint8_t maxNBMismatch;
    
    uint8_t * outputStatus;

    unsigned int * outputBlock;
    uint16_t * occCount;
    MICSRAOccMetadata * metaBlock;
    uint16_t * metaCount;
    
    uint8_t outputType;
    
    BWT * bwt;
    BWT * rev_bwt;
    LT * lt;
    LT * rlt;
    
    SRAWorkingMemory AlgnmtMemory;
    
    // Used by Seed alignment only
    // This paramater needs to be initialised as zero for 
    // non-seed alignment as it has effect on all 
    // occurrences found by MIC-SRA. 
    // 
    // Jeanno: I suggest the Deep DP process should be organised in a struct
    // And then put these seed related variables inside.
    uint16_t seedOffset;
    uint16_t seedOffset_Complt;
} MICSRAArguments;

__attribute__((target(mic)))
MICSRAArguments * MICSRAARGConstruct();
__attribute__((target(mic)))
void MICSRAARGFree(MICSRAArguments * micArgs);
__attribute__((target(mic)))
MICSRAArguments * MICSRAARGMakeCopy(MICSRAArguments * micArgs);
__attribute__((target(mic)))
MICSRAArguments * MICSRAARGCopy(MICSRAArguments * destination, MICSRAArguments * source);

void MICSRAMETADebugPrint(MICSRAOccMetadata * meta);

#endif

