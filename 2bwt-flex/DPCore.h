//
//    DPCore.h
//
//    SOAP2 / 2BWT
//
//    Copyright (C) 2012, HKU
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
/////////////////////////////////////////////////////
/*

   Modification History 
   
   Date   : 3rd February 2012
   Author : Edward MK Wu
   Change : New file.

*/
/////////////////////////////////////////////////////

#ifndef __DP_CORE_H__
#define __DP_CORE_H__

#include "SRACore.h"

#define DP_MAX_REGION_LENGTH        (CHAR_PER_WORD * 4096)
#define DP_MAX_ERROR_COUNT          (SRA_MAX_READ_LENGTH * 2 + 2)
#define DP_2BWT_TRIGGER_THRESHOLD   100

#define DP_2BWT_REG_EXPAND_FACTOR   1.50f
#define DP_2BWT_REG_SPOS_CORRECTION 0.25f

#define DP_ALIGNMENT_COMPLETED    0
#define DP_ALIGNMENT_INITIALISED  1
#define DP_ALIGNMENT_INPUT_ERROR  2
#define DP_ALIGNMENT_MATRIX_READY 3
#define DP_ALIGNMENT_FAILED       4

// This are the 3-bit vector for tracing matrix
// The Dynamics Programming alignment algorithm we employed is the affine
// gap penality 3 matrix models. The algoritm is as below:
//
//               M[i-1,j-1] + MATCH / MISMATCH 
// M[i,j] = max{ I[i-1,j-1] + MATCH / MISMATCH    // match/mismatch
//               D[i-1,j-1] + MATCH / MISMATCH
//
//               M[i,j-1] + GAPOPEN 
// I[i,j] = max{ I[i,j-1] + GAPEXTEND             // insertion on read
//               D[i,j-1] + GAPOPEN 
//
//               M[i-1,j] + GAPOPEN
// D[i,j] = max{ I[i-1,j] + GAPOPEN               // insertion on reference sequence
//               D[i-1,j] + GAPEXTEND
 //
// DP_TRACE_M means the score value can be from Matrix M
// DP_TRACE_I means the score value can be from Matrix I
// DP_TRACE_D means the score value can be from Matrix D
// _____________
// |  M  |  I  |
// |_____|_____| i.e. when DPMatrixCellScore is called for X, 
// |  D  |  X  | we assume the value of A, B, C are ready.
// |_____|_____| 
#define DP_TRACE_BASE             0
#define DP_TRACE_M                1
#define DP_TRACE_I                2
#define DP_TRACE_D                4

#define DP_MATRIX_ID_M            0
#define DP_MATRIX_ID_I            1
#define DP_MATRIX_ID_D            2
#define DP_MATRIX_ID_COUNT        3

// Obsolete ----------------------------
#define DP_TRACE_MISMATCH         1
#define DP_TRACE_REG_INSERT       2
#define DP_TRACE_READ_INSERT      4
// -------------------------------------

#define DP_MM_COORD(regIdx,readIdx)            ((regIdx)*(SRA_MAX_READ_LENGTH+1)+(readIdx))

typedef struct DPMatchElem {
    uint16_t type:4,length:12;
} DPMatchElem;

typedef struct DPOccurrence {
    uint8_t type;
    unsigned long long ambPosition;
    uint8_t strand;
    uint16_t matchElemsCount;
    DPMatchElem matchElems[DP_MAX_ERROR_COUNT];
    uint16_t mismatchCount;
    uint16_t matchLen;
    int16_t matchScore;
    
    // Auxiliary data from DP input
    unsigned int input_regionStartIdx;
    unsigned int input_regionLength;
} DPOccurrence;

typedef struct DPScores {
    int dpMatch;
    int dpMismatch;
    int dpGapOpen;
    int dpGapExtend;
} DPScores;

typedef struct DPMatrixCell {
    unsigned char dpTrace;
    int16_t dpScore;
} DPMatrixCell;

typedef struct DPWork {
    //matrix is allocated to the size of M * N where M is the length of the read
    //and N is the length of the allowed region on reference sequence.
    DPMatrixCell matrix[DP_MATRIX_ID_COUNT][(SRA_MAX_READ_LENGTH+1)*(DP_MAX_REGION_LENGTH+1)];

    unsigned char _regionBuffer[DP_MAX_REGION_LENGTH];
    unsigned char * regionBuffer;
    unsigned int regionStartIdx;
    unsigned int regionLength;
    
    unsigned char * readCode;
    unsigned int readLength;
    
    //leadVerifyLength and tailVerifyLength defines the lengths of the
    //section of region that is only input into Dynamic Programming for
    //verification such that occurrences found at these regions are not
    //considered result. For instance,
    //
    //   A             B            C
    // -----|--------------------|----- region extracted from ref. sequence
    //
    // in this case leadVerifyLength is 5 and the part of region is marked with A;
    // tailVerifyLength is 5 and the part of region is marked with C. Only
    // alignment that starts on region B is output by DPBackTrack.
    // the length of A, B, C should added up equal to regionLength defined above.
    unsigned int leadVerifyLength;
    unsigned int tailVerifyLength;
    
    //leadSoftClipLength and tailSoftClipLength defines the lengths of the
    //section of read that can be clipped to achieve better Dynamic Programming
    //score. 
    unsigned int leadSoftClipLength;
    unsigned int tailSoftClipLength;
    unsigned int totalSoftClipLength;
    
    //leadHardClipLength and tailHardClipLength defines the leading/trailing lengths of the
    //region that must be clipped prior to report. The difference between
    //this and VerifyLength is that the latter attempts to recover the best alignment that
    //fulfils the criteria which no alignment should start/end on that part of the region;
    //while hard-clip is only checked against the best alignment and if hard-clip is needed
    //it will always happen and if the score after clipping falls below the score threshold
    //the occurrence will be suppressed.
    //  NOTE: Verification Length and Hard Clip Length are NOT supposed
    //        to be used together.
    unsigned int leadHardClipLength;
    unsigned int tailHardClipLength;
    
    int flag;
} DPWork;

void DPOccurrenceCopy(DPOccurrence * source, DPOccurrence * destination);


DPWork * DPWorkCreate();
void DPWorkInitialise(DPWork * dpWork);
void DPWorkFree(DPWork * dpWork);
DPScores * DPScoresConstruct();
void DPScoresFree(DPScores * dpScores);

void DEBUGDPPrintWork(DPWork * dpWork);

#endif
