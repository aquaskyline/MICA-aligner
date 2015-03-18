//
//    DPOperations.c
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

#include "DPOperations.h"

// --------------------------------------------------
// DP_DEBUG #1
// Uncomment the below to turn on logging message for
// PE matching of two occurrences.
//#define DP_DEBUG_PRINT_MATCHING_PROGRESS
// --------------------------------------------------
// --------------------------------------------------
// DP_DEBUG #2
// Uncomment the below to turn on visual logging for
// DP match
//#define DP_DEBUG_PRINT_VISUAL_RESULT
// --------------------------------------------------

void DPDebugPrintMatch(DPWork * dpWork, DPOccurrence * dpOcc);

void DPMatrixFetching(DPWork * dpWork, HSP * hsp) {
    unsigned int packedDNAStartBlock =  dpWork->regionStartIdx / CHAR_PER_WORD;
    unsigned int packedDNAEndBlock =  (dpWork->regionStartIdx + dpWork->regionLength + (CHAR_PER_WORD - 1)) / CHAR_PER_WORD - 1;
    unsigned int packedDNABlockIdx = 0;
    int expandRegionIdx = 0;
    
    unsigned char * _regionBuffer = dpWork->_regionBuffer;
    int buffer_i;
    unsigned int word;
    int charMask = (1<<BIT_PER_CHAR)-1;
    
    // Read in the reference sequence into buffer.    
    for (packedDNABlockIdx=packedDNAStartBlock;packedDNABlockIdx<=packedDNAEndBlock;packedDNABlockIdx++) {
        word = hsp->packedDNA[packedDNABlockIdx];
        for (buffer_i=0; buffer_i<CHAR_PER_WORD; buffer_i++) {
            _regionBuffer[expandRegionIdx++] = (unsigned char)((word>>((CHAR_PER_WORD-buffer_i-1)*BIT_PER_CHAR)) & charMask);
        }
    } 
    dpWork->regionBuffer = _regionBuffer + dpWork->regionStartIdx % CHAR_PER_WORD;
}

inline void DPMatrixFillCell(DPWork * dpWork, DPScores * dpScores, unsigned int reg_i, unsigned int read_i) {

    DPMatrixCell * dpMatrixM = dpWork->matrix[DP_MATRIX_ID_M];
    DPMatrixCell * dpMatrixI = dpWork->matrix[DP_MATRIX_ID_I];
    DPMatrixCell * dpMatrixD = dpWork->matrix[DP_MATRIX_ID_D];

    unsigned char reg_c  = dpWork->regionBuffer[reg_i-1];
    unsigned char read_c = dpWork->readCode[read_i-1];
    
    unsigned int mCoordA = DP_MM_COORD(reg_i-1,read_i-1);
    unsigned int mCoordB = DP_MM_COORD(reg_i-1,read_i); // delete on reference sequence
    unsigned int mCoordC = DP_MM_COORD(reg_i,read_i-1); // insertion on reference sequence
    unsigned int mCoord  = DP_MM_COORD(reg_i,read_i);
    
    int dpGapOpen   = dpScores->dpGapOpen;
    int dpGapExtend = dpScores->dpGapExtend;
    
    int A, B, C, max;
    
    int matchScore = (dpScores->dpMismatch * (reg_c!=read_c)) + 
                     (dpScores->dpMatch    * (reg_c==read_c));

    // Matrix M
    // =============================
    
    A = dpMatrixM[mCoordA].dpScore + matchScore;
    B = dpMatrixI[mCoordA].dpScore + matchScore;
    C = dpMatrixD[mCoordA].dpScore + matchScore;
    
    max = A;
    if (B>max) {max = B;}
    if (C>max) {max = C;}
    
    dpMatrixM[mCoord].dpTrace = (DP_TRACE_M    * (max==A)) | 
                                (DP_TRACE_I    * (max==B)) |
                                (DP_TRACE_D    * (max==C));
    dpMatrixM[mCoord].dpScore = max;
    
    #ifdef DP_DEBUG_PRINT_VISUAL_RESULT
        printf("Matrix M = %d/%d\n",dpMatrixM[mCoord].dpScore,dpMatrixM[mCoord].dpTrace);
    #endif
    
    // Matrix I
    // =============================
    
    A = dpMatrixM[mCoordB].dpScore + dpGapOpen;
    B = dpMatrixI[mCoordB].dpScore + dpGapExtend;
    C = dpMatrixD[mCoordB].dpScore + dpGapOpen;
    
    max = A;
    if (B>max) {max = B;}
    if (C>max) {max = C;}
    
    dpMatrixI[mCoord].dpTrace = (DP_TRACE_M    * (max==A)) | 
                                (DP_TRACE_I    * (max==B)) |
                                (DP_TRACE_D    * (max==C));
    dpMatrixI[mCoord].dpScore = max;
    
    #ifdef DP_DEBUG_PRINT_VISUAL_RESULT
        printf("Matrix I = %d/%d\n",dpMatrixI[mCoord].dpScore,dpMatrixI[mCoord].dpTrace);
    #endif
    
    // Matrix D
    // =============================
    
    A = dpMatrixM[mCoordC].dpScore + dpGapOpen;
    B = dpMatrixI[mCoordC].dpScore + dpGapOpen;
    C = dpMatrixD[mCoordC].dpScore + dpGapExtend;
    
    max = A;
    if (B>max) {max = B;}
    if (C>max) {max = C;}
    
    dpMatrixD[mCoord].dpTrace = (DP_TRACE_M    * (max==A)) | 
                                (DP_TRACE_I    * (max==B)) |
                                (DP_TRACE_D    * (max==C));
    dpMatrixD[mCoord].dpScore = max;
    
    #ifdef DP_DEBUG_PRINT_VISUAL_RESULT
        printf("Matrix D = %d/%d\n",dpMatrixD[mCoord].dpScore,dpMatrixD[mCoord].dpTrace);
    #endif
    
}

int DPMatrixFill(DPArguments * dpArgs) {

    DPWork * dpWork = dpArgs->dpWork;
    DPScores * dpScores = dpArgs->dpScores;
    DPMatrixCell * dpMatrixM = dpWork->matrix[DP_MATRIX_ID_M];
    DPMatrixCell * dpMatrixI = dpWork->matrix[DP_MATRIX_ID_I];
    DPMatrixCell * dpMatrixD = dpWork->matrix[DP_MATRIX_ID_D];
    
    unsigned int reg_i, read_i;
    unsigned int mCoord;
    unsigned char * regionBuffer = dpWork->regionBuffer;
    unsigned int regionIdx = 0;
    
    #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
    printf("\nDPMatrixFill -- Invoke\n");
    #endif
    
    #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
    printf("pos = %u\n", dpWork->regionStartIdx);
    printf("offset = %u\n", dpWork->regionStartIdx % CHAR_PER_WORD);
    printf("length = %u\n", dpWork->regionLength);
    #endif
    
    #ifdef DP_DEBUG_PRINT_VISUAL_RESULT
        for (regionIdx=0;regionIdx<dpWork->regionLength;regionIdx++) {
            printf("%c",dnaChar[regionBuffer[regionIdx]]);
        }
        printf("\n");
    #endif
    
    if (dpWork->flag != DP_ALIGNMENT_INITIALISED) {
        printf("Uninitialised DpWork.\n");
        return 0;
    }
    
    //Initialise the matching matrix giving the first column all zeros
    #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
    printf("DPMatrixFill -- Initialising matching matrix\n");
    #endif
    
    #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
        printf("DPMatrixFill -- dpScores->dpMatch = %d\n",dpScores->dpMatch);
        printf("DPMatrixFill -- dpScores->dpMismatch = %d\n",dpScores->dpMismatch);
        printf("DPMatrixFill -- dpScores->dpGapOpen = %d\n",dpScores->dpGapOpen);
        printf("DPMatrixFill -- dpScores->dpGapExtend = %d\n",dpScores->dpGapExtend);
    #endif
    
    if ( dpWork->tailVerifyLength > dpWork->regionLength ) {
        printf("Invalid set up of tailVerifyLength. %u is greater than regionLength %u\n",dpWork->tailVerifyLength,dpWork->regionLength);
        dpWork->flag = DP_ALIGNMENT_FAILED;
        return 0;
    }
    
    int lastScore;
    unsigned int occLength = dpWork->regionLength-dpWork->tailVerifyLength;
    
    // Initialise the matching matrix s.t. it 
    // sets the entire first column as zero.
    for (reg_i=0;reg_i<=occLength;reg_i++) {
        mCoord = DP_MM_COORD(reg_i,0);
        dpMatrixM[mCoord].dpScore = 0;
        dpMatrixI[mCoord].dpScore = dpScores->dpGapOpen;
        dpMatrixD[mCoord].dpScore = dpScores->dpGapOpen;
        dpMatrixM[mCoord].dpTrace = DP_TRACE_BASE;
        dpMatrixI[mCoord].dpTrace = DP_TRACE_BASE;
        dpMatrixD[mCoord].dpTrace = DP_TRACE_BASE;
    }
    
    // Special handling for tailVerifyLength Enhancement
    // ------------------------------------------------
    // Initialise the later part of the first column 
    // to gaps.
    if (dpWork->tailVerifyLength>0) {
    
        mCoord = DP_MM_COORD(reg_i,0);
        dpMatrixM[mCoord].dpScore = dpScores->dpGapOpen;
        dpMatrixI[mCoord].dpScore = dpScores->dpGapOpen;
        dpMatrixD[mCoord].dpScore = dpScores->dpGapOpen;
        dpMatrixM[mCoord].dpTrace = DP_TRACE_I;
        dpMatrixI[mCoord].dpTrace = DP_TRACE_I;
        dpMatrixD[mCoord].dpTrace = DP_TRACE_I;
        lastScore = dpMatrixI[mCoord].dpScore;
        reg_i++;
        
        for (;reg_i<=dpWork->regionLength;reg_i++) {
            mCoord = DP_MM_COORD(reg_i,0);
            dpMatrixM[mCoord].dpScore = lastScore + dpScores->dpGapExtend;
            dpMatrixI[mCoord].dpScore = lastScore + dpScores->dpGapExtend;
            dpMatrixD[mCoord].dpScore = lastScore + dpScores->dpGapExtend;
            dpMatrixM[mCoord].dpTrace = DP_TRACE_I;
            dpMatrixI[mCoord].dpTrace = DP_TRACE_I;
            dpMatrixD[mCoord].dpTrace = DP_TRACE_I;
            lastScore = dpMatrixI[mCoord].dpScore;
        }
    }
    
    #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
        printf("DPMatrixFill\n");
        for (reg_i=0;reg_i<=dpWork->regionLength;reg_i++) {
            mCoord = DP_MM_COORD(reg_i,0);
            printf("%10u : %10d : %10u\n",reg_i,dpMatrixM[mCoord].dpScore,dpMatrixM[mCoord].dpTrace);
        }
    #endif
    
    // End of tailVerifyLength Enhancement
    // ----------------------------------
    
    mCoord = DP_MM_COORD(0,0);
    dpMatrixM[mCoord].dpScore = 0;
    dpMatrixI[mCoord].dpScore = 0;
    dpMatrixD[mCoord].dpScore = 0;
    dpMatrixM[mCoord].dpTrace = DP_TRACE_BASE;
    dpMatrixI[mCoord].dpTrace = DP_TRACE_BASE;
    dpMatrixD[mCoord].dpTrace = DP_TRACE_BASE;
    
    mCoord = DP_MM_COORD(0,1);
    dpMatrixM[mCoord].dpScore = dpScores->dpGapOpen;
    dpMatrixI[mCoord].dpScore = dpScores->dpGapOpen;
    dpMatrixD[mCoord].dpScore = dpScores->dpGapOpen;
    dpMatrixM[mCoord].dpTrace = DP_TRACE_BASE;
    dpMatrixI[mCoord].dpTrace = DP_TRACE_BASE;
    dpMatrixD[mCoord].dpTrace = DP_TRACE_BASE;
    lastScore = dpMatrixD[mCoord].dpScore;
    
    for (read_i=2;read_i<=dpWork->readLength;read_i++) {
        mCoord = DP_MM_COORD(0,read_i);
        dpMatrixM[mCoord].dpScore = lastScore + dpScores->dpGapExtend;
        dpMatrixI[mCoord].dpScore = lastScore + dpScores->dpGapExtend;
        dpMatrixD[mCoord].dpScore = lastScore + dpScores->dpGapExtend;
        dpMatrixM[mCoord].dpTrace = DP_TRACE_BASE;
        dpMatrixI[mCoord].dpTrace = DP_TRACE_BASE;
        dpMatrixD[mCoord].dpTrace = DP_TRACE_BASE;
        lastScore = dpMatrixD[mCoord].dpScore;
    }
    
    #ifdef DP_DEBUG_PRINT_VISUAL_RESULT
        DPDebugPrintMatrix(dpWork);
    #endif
    
    #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
        printf("DPMatrixFill -- Filling matching matrix\n");
    #endif
    
    for (reg_i=1;reg_i<=dpWork->regionLength;reg_i++) {
        for (read_i=1;read_i<=dpWork->readLength;read_i++) {
            DPMatrixFillCell(dpWork,dpScores,reg_i,read_i);
        }
    }
    
    #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
        printf("DPMatrixFill -- DP Matrix completely filled.\n");
    #endif
    
    dpWork->flag = DP_ALIGNMENT_MATRIX_READY;
    return 1;
}


//DPMatrixReadMatchElem retrieves the DP matching elements from the trace.
//It always prefers mismatch over indel.
inline void DPMatrixReadMatchElem(DPWork * dpWork, unsigned int reg_i, unsigned int read_i, unsigned int matrix_i, DPOccurrence * dpOcc) {

    unsigned char reg_c, read_c;
    
    read_c = dpWork->readCode[read_i-1];
    if (reg_i>0) reg_c = dpWork->regionBuffer[reg_i-1];
    
    int trace = dpWork->matrix[matrix_i][DP_MM_COORD(reg_i,read_i)].dpTrace;
    
    int matchType = SRA_CHAR_ERROR_TYPE_MISMATCH;
    
    // Error Type by the Matrix ID
    if ( matrix_i == DP_MATRIX_ID_M ) {
        //Mismatch or Match
        matchType = SRA_CHAR_ERROR_TYPE_MISMATCH;
        dpOcc->mismatchCount += reg_c!=read_c;
        dpOcc->matchLen++;
    } else if ( matrix_i == DP_MATRIX_ID_I ) {
        //Gap Open
        matchType = SRA_CHAR_ERROR_TYPE_INSERT;
        dpOcc->matchLen++;
    } else if ( matrix_i == DP_MATRIX_ID_D ) {
        //Gap Open
        matchType = SRA_CHAR_ERROR_TYPE_DELETE;
    } else {
        printf("[SAFE-GUARD] Unexpected Matrix Idx %d\n",matrix_i);
        printf("Alignment unable to continue.\n");
        exit(1);
    }
    
    if (dpOcc->matchElemsCount == 0) {
        dpOcc->matchElemsCount++;
        dpOcc->matchElems[dpOcc->matchElemsCount-1].length = 1;
        dpOcc->matchElems[dpOcc->matchElemsCount-1].type = matchType;
    } else if (dpOcc->matchElems[dpOcc->matchElemsCount-1].type == matchType) {
        dpOcc->matchElems[dpOcc->matchElemsCount-1].length++;
    } else {
        if (dpOcc->matchElems[dpOcc->matchElemsCount-1].length>0) dpOcc->matchElemsCount++;
        dpOcc->matchElems[dpOcc->matchElemsCount-1].length = 1;
        dpOcc->matchElems[dpOcc->matchElemsCount-1].type = matchType;
    }
}

//DPMatrixBackTrackCell finds the trace of the given cell. This is an atomic operation for DP BackTracing
//Also, this function prefers mismatch over indel greedily for each cell.
inline void DPMatrixBackTrackCell(DPWork * dpWork, unsigned int * reg_i, unsigned int * read_i, unsigned int * matrix_i) {

    int trace = dpWork->matrix[(*matrix_i)][DP_MM_COORD((*reg_i),(*read_i))].dpTrace;
    
    // Position Moving by Matrix Idx
    if ( (*matrix_i) == DP_MATRIX_ID_M ) {
        (*reg_i)--;
        (*read_i)--;
    } else if ( (*matrix_i) == DP_MATRIX_ID_I ) {
        (*reg_i)--;
    } else if ( (*matrix_i) == DP_MATRIX_ID_D ) {
        (*read_i)--;
    } else {
        printf("[SAFE-GUARD] Unexpected Matrix ID %d\n",(*matrix_i));
        printf("Alignment unable to continue.\n");
        exit(1);
    }
    
        
    // Matrix Switch by Tracing Track
    if ( (trace & DP_TRACE_M) > 0 ) {
        //Switching to Matrix M
        (*matrix_i) = DP_MATRIX_ID_M;
    } else if ( (trace & DP_TRACE_I) > 0 ) {
        //Switching to Matrix I
        (*matrix_i) = DP_MATRIX_ID_I;
    } else if ( (trace & DP_TRACE_D) > 0 ) {
        //Switching to Matrix D
        (*matrix_i) = DP_MATRIX_ID_D;
    } else {
        printf("[SAFE-GUARD] Unexpected DP Trace Value %d\n",trace);
        printf("Alignment unable to continue.\n");
        exit(1);
    }
}

int DPBackTrack(DPWork * dpWork, DPScores * dpScores, int minRecoverScore, DPOccurrence * dpOcc) {
    
    unsigned int reg_i, read_i, mCoord, matrix_i;
    
    // Initialise maxScore catcher to smaller than the minRecoverScore
    int maxScore = minRecoverScore-1;
    
    int maxScore_Regi = -1;
    int maxScore_Readi = -1;
    int maxScore_Matrixi = -1;
    
    // Auxiliary Data
    dpOcc->input_regionStartIdx = dpWork->regionStartIdx;
    dpOcc->input_regionLength = dpWork->regionLength;
    
    #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
        printf("dpWork->regionStartIdx = %u\n",dpWork->regionStartIdx);
        printf("dpWork->regionLength = %u\n",dpWork->regionLength);
        printf("dpWork->read = ");
        for (read_i=0;read_i<dpWork->readLength;read_i++) {
            printf("%c",dnaChar[dpWork->readCode[read_i]]);
        }
        printf("\n");
    #endif
    
    #ifdef DP_DEBUG_PRINT_VISUAL_RESULT
        DPDebugPrintMatrix(dpWork);
    #endif
    
    if ( dpWork->leadVerifyLength > dpWork->regionLength ) {
        printf("Invalid set up of leadVerifyLength. %u is greater than regionLength %u\n",dpWork->leadVerifyLength,dpWork->regionLength);
        dpWork->flag = DP_ALIGNMENT_FAILED;
        return 0;
    }
    
    //////////////////////////////////////////////////////////
    //
    // OBTAIN BEST SCORE
    //
    //////////////////////////////////////////////////////////
    
    // A. Default logic to obtain the maximum score from the last column
    maxScore_Regi = 0;
    maxScore_Readi = dpWork->readLength;
    for (reg_i=1+dpWork->leadVerifyLength;reg_i<=dpWork->regionLength;reg_i++) {
        mCoord = DP_MM_COORD(reg_i,dpWork->readLength);
        for (matrix_i=0;matrix_i<DP_MATRIX_ID_COUNT;matrix_i++) {
            if (dpWork->matrix[matrix_i][mCoord].dpScore > maxScore) {
                maxScore = dpWork->matrix[matrix_i][mCoord].dpScore;
                maxScore_Regi = reg_i;
                maxScore_Matrixi = matrix_i;
            }
        }
    }
    
    // B. Soft-Clipping enhancement to obtain the maximum score from the last columns
    // to allow the tail of the read to be clipped ( consider never existed ).
    
    // Set up the softClipRemain variables to denote how many characters can be clipped from tail.
    int softCilpRemain = dpWork->totalSoftClipLength;
    if ( softCilpRemain > dpWork->tailSoftClipLength ) {
        softCilpRemain = dpWork->tailSoftClipLength;
    }
    #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
        printf("\nStarting Maximum Score lookup...\n");
        printf("[HeadSoftClip] softCilpRemain = %d\n",softCilpRemain);
    #endif
    if (softCilpRemain>0) {
        for (reg_i=1+dpWork->leadVerifyLength;reg_i<=dpWork->regionLength;reg_i++) {
            for (read_i=dpWork->readLength-softCilpRemain;read_i<dpWork->readLength;read_i++) {
                mCoord = DP_MM_COORD(reg_i,read_i);
                for (matrix_i=0;matrix_i<DP_MATRIX_ID_COUNT;matrix_i++) {
                    if (dpWork->matrix[matrix_i][mCoord].dpScore > maxScore) {
                        maxScore = dpWork->matrix[matrix_i][mCoord].dpScore;
                        maxScore_Regi = reg_i;
                        maxScore_Readi = read_i;
                        maxScore_Matrixi = matrix_i;
                    }
                }
            }    
        }
    }

    // Set up the softClipRemain variables to denote how many characters can be clipped from head.
    // This is set now but not later because the hard-clip logic might change read_i.
    softCilpRemain = dpWork->totalSoftClipLength - dpWork->readLength - read_i;
    if ( softCilpRemain > dpWork->leadSoftClipLength ) {
        softCilpRemain = dpWork->leadSoftClipLength;
    }

    #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
        printf("DPBackTrack -- Maximum score found is Matrix[%u](%u,%u) = %d\n",maxScore_Matrixi,maxScore_Regi,maxScore_Readi,maxScore);
    #endif

    if (maxScore < minRecoverScore) {
        //Score does not meet minRecoverScore criteria and being discarded
        #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
            printf("Score does not meet minRecoverScore criteria and being discarded. (%d<%d)\n",maxScore,minRecoverScore);
        #endif
        dpOcc->ambPosition = 0;
        dpOcc->matchElemsCount = 0;
        dpOcc->mismatchCount = 0;
        dpOcc->matchLen = 0;
        dpOcc->matchScore = 0;
        return 0;
    }
    
    
    
    
    
    //////////////////////////////////////////////////////////
    //
    // TAIL HARD CLIPPING
    //
    //////////////////////////////////////////////////////////
    
    reg_i = maxScore_Regi;
    read_i = maxScore_Readi;
    matrix_i = maxScore_Matrixi;
    
    dpOcc->mismatchCount = 0;
    dpOcc->matchLen = 0;
    
    #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
        printf("\nStarting Backtracking...\n");
        printf("[HeadSoftClip] softCilpRemain = %d\n",softCilpRemain);
    #endif
    
    // Initial hard-clip backward Trace before the actual one
    // The max score cell found by the above logics have not taken into
    // account the tailHardClipLength.
    if ( dpWork->tailHardClipLength > 0 ) {
    
        int tailLength = dpWork->regionLength - reg_i + 1;
        
        #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
            printf("[HeadSoftClip] tailLength = %d\n",tailLength);
        #endif
        
        mCoord = DP_MM_COORD(reg_i,read_i);
        while ( dpWork->matrix[matrix_i][ mCoord ].dpTrace != DP_TRACE_BASE &&
                tailLength <= dpWork->tailHardClipLength ) {
                
            #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
                printf("DPHardTrack -- Cursor at Matrix[%u](%u,%u) = %d / %d\n",matrix_i,reg_i,read_i,
                    dpWork->matrix[matrix_i][ mCoord ].dpScore,
                    dpWork->matrix[matrix_i][ mCoord ].dpTrace);
            #endif
            
            DPMatrixReadMatchElem(dpWork,reg_i,read_i,matrix_i,dpOcc);
            DPMatrixBackTrackCell(dpWork,&reg_i,&read_i,&matrix_i);
            mCoord = DP_MM_COORD(reg_i,read_i);
            
            #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
                printf("  [Track] New Cursor (%u,%u) matchElemsCount=%d matchElems.type=%d matchElems.length=%d\n",reg_i,read_i,dpOcc->matchElemsCount,dpOcc->matchElems[dpOcc->matchElemsCount-1].type,dpOcc->matchElems[dpOcc->matchElemsCount-1].length);
            #endif
            
            tailLength = dpWork->regionLength - reg_i + 1;
        }

        maxScore = dpWork->matrix[matrix_i][mCoord].dpScore;
        maxScore_Regi = reg_i;
        maxScore_Readi = read_i;
        maxScore_Matrixi = matrix_i;
        
        #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
            printf("DPBackTrack -- Maximum score found is Matrix[%u](%u,%u) = %d\n",maxScore_Matrixi,maxScore_Regi,maxScore_Readi,maxScore);
        #endif

        if (maxScore < minRecoverScore) {
            //Score does not meet minRecoverScore criteria and being discarded
            #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
                printf("Score does not meet minRecoverScore criteria and being discarded. (%d<%d)\n",maxScore,minRecoverScore);
            #endif
            dpOcc->ambPosition = 0;
            dpOcc->matchElemsCount = 0;
            dpOcc->mismatchCount = 0;
            dpOcc->matchLen = 0;
            dpOcc->matchScore = 0;
            return 0;
        }
    }
    
    
    
    
    
    
    //////////////////////////////////////////////////////////
    //
    // BACK-TRACKING
    //
    //////////////////////////////////////////////////////////
    
    // Set up DP occurrences after Soft/Hard-Clipping
    if ( read_i < dpWork->readLength ) {
        dpOcc->matchElemsCount = 1;
        dpOcc->matchElems[0].type = SRA_CHAR_ERROR_TYPE_SOFTCLIP;
        dpOcc->matchElems[0].length = dpWork->readLength - read_i;
    } else {
        dpOcc->matchElemsCount = 0;
    }
    
    // Actual backtracking taking place..
    int minScore_j = 0, matchLen_j, matchElemsCount_j, matchElemsLength_j, read_j = 0, reg_j = 0;
    mCoord = DP_MM_COORD(reg_i,read_i);
    while ( dpWork->matrix[matrix_i][ mCoord ].dpTrace != DP_TRACE_BASE ) {
        #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
            printf("DPBackTrack -- Cursor at Matrix[%u](%u,%u) = %d / %d\n",matrix_i,reg_i,read_i,
                dpWork->matrix[matrix_i][ mCoord ].dpScore,
                dpWork->matrix[matrix_i][ mCoord ].dpTrace);
        #endif
        
        // Head and Total Soft Clipping
        // We are noting down the matchElems information for post-processing
        if ( read_i <= softCilpRemain ) {
            if ( minScore_j > dpWork->matrix[matrix_i][ mCoord ].dpScore ) {
                minScore_j = dpWork->matrix[matrix_i][ mCoord ].dpScore;
                read_j = read_i;
                reg_j = reg_i;
                matchElemsCount_j = dpOcc->matchElemsCount;
                matchElemsLength_j = dpOcc->matchElems[matchElemsCount_j-1].length;
                matchLen_j = dpOcc->matchLen;
                
                #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
                    printf("  [HeadSoftClip] New most-negative score %d\n",minScore_j);
                #endif
            }
        }
        
        // Head Hard-Clipping
        if ( reg_i <= dpWork->leadHardClipLength ) {
            minScore_j = dpWork->matrix[matrix_i][ mCoord ].dpScore;
            read_j = read_i;
            reg_j = reg_i;
            matchElemsCount_j = dpOcc->matchElemsCount;
            matchElemsLength_j = dpOcc->matchElems[matchElemsCount_j-1].length;
            matchLen_j = dpOcc->matchLen;
            
            #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
                printf("  [HeadHardClip] Backtracking reached %u, breached the bound of %u.\n",reg_i,dpWork->leadHardClipLength);
            #endif
            
            break;
        }
    
        DPMatrixReadMatchElem(dpWork,reg_i,read_i,matrix_i,dpOcc);
        
        DPMatrixBackTrackCell(dpWork,&reg_i,&read_i,&matrix_i);
        
        mCoord = DP_MM_COORD(reg_i,read_i);
            
        #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
            printf("  [Track] New Cursor (%u,%u) matchElemsCount=%d matchElems.type=%d matchElems.length=%d\n",reg_i,read_i,dpOcc->matchElemsCount,dpOcc->matchElems[dpOcc->matchElemsCount-1].type,dpOcc->matchElems[dpOcc->matchElemsCount-1].length);
        #endif
    }
    
    // Setup the ambPosition prior to head-clipping correction.
    dpOcc->ambPosition = reg_j + dpWork->regionStartIdx;
    
    #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
        printf("DPBackTrack -- Cursor at Matrix[%u](%u,%u) = %d / %d\n",matrix_i,reg_i,read_i,
            dpWork->matrix[matrix_i][ DP_MM_COORD(reg_i,read_i) ].dpScore,
            dpWork->matrix[matrix_i][ DP_MM_COORD(reg_i,read_i) ].dpTrace);
        printf("Occurrence position at %u (mis=%u matchLen=%d)\n", reg_i + dpWork->regionStartIdx, dpOcc->mismatchCount, dpOcc->matchLen);
        printf("Backtracking completed...\n\n");
    #endif
    
    #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
        DPDebugPrintMatch(dpWork,dpOcc);
    #endif
    
    if ( minScore_j < 0 ) {

        #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
            printf("[HeadSoftClip] Clipping at %d\n",read_j);
        #endif
        
        // In this case the minScore is negative and we should be clipping something...
        maxScore -= minScore_j;
        
        dpOcc->matchElemsCount = matchElemsCount_j;
        
        // Handling special case when the clipping occurs in the 
        // middle of a gap. 
        //                          
        // RefSeq---------------------------------------
        //                 /         / \          \
        // Read         5'/---------/123\----------\3'
        //
        // If soft-clipping happens to occur in the (2) position then when we 
        // clip 5' to 2 we will beed to replace the score of (2) with a
        // gap-open to make sure point (3) remains as a gap-extend.
        if ( dpOcc->matchElems[dpOcc->matchElemsCount-1].type == SRA_CHAR_ERROR_TYPE_DELETE &&
             matchElemsLength_j <= dpOcc->matchElems[dpOcc->matchElemsCount-1].length ) {
             
            #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
                printf("[HeadSoftClip] In the middle of any insertion on reference sequence\n");
            #endif
            
            maxScore += dpScores->dpGapOpen;
            
            dpOcc->matchElems[dpOcc->matchElemsCount-1].length = matchElemsLength_j;
            dpOcc->matchElemsCount++;
            
        } else {
        // Normal Handling
            if ( matchElemsLength_j <=  dpOcc->matchElems[dpOcc->matchElemsCount-1].length ) {
                dpOcc->matchElems[dpOcc->matchElemsCount-1].length = matchElemsLength_j;
                dpOcc->matchElemsCount++;
            }
        }
        
        dpOcc->matchElems[dpOcc->matchElemsCount-1].type = SRA_CHAR_ERROR_TYPE_SOFTCLIP;
        dpOcc->matchElems[dpOcc->matchElemsCount-1].length = read_j;
        dpOcc->matchLen = matchLen_j;
    
        #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
            printf("[HeadSoftClip] Maximum score after clipping is %d\n",maxScore);
        #endif
    }
    
    #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
        DPDebugPrintMatch(dpWork,dpOcc);
    #endif
    
    #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
        printf("DPBackTrack -- Matching Score after head tracking = %d\n",maxScore);
    #endif

    if (maxScore < minRecoverScore) {
        //Score does not meet minRecoverScore criteria and being discarded
        #ifdef DP_DEBUG_PRINT_MATCHING_PROGRESS
            printf("Score does not meet minRecoverScore criteria and being discarded. (%d<%d)\n",maxScore,minRecoverScore);
        #endif
        dpOcc->ambPosition = 0;
        dpOcc->matchElemsCount = 0;
        dpOcc->mismatchCount = 0;
        dpOcc->matchLen = 0;
        dpOcc->matchScore = 0;
        return 0;
    }
    
    // Eventually set up the DP Occurrences Matching Score after all
    // the post-processing.
    dpOcc->matchScore = maxScore;
    
    return 1;
}












//Debug Functions Below This Line ----------------------------------------------------------------------------------

void DPDebugPrintTrace(int dpTrace) {
    char tmp[8];
    tmp[0]='B';
    tmp[1]='M';
    tmp[2]='I';
    tmp[3]='?';
    tmp[4]='D';
    tmp[5]='?';
    tmp[6]='?';
    tmp[7]='?';
    
    char output[4];
    int i = 0;
    if (dpTrace == DP_TRACE_BASE) {
        output[i++] = tmp[dpTrace];
    } else {
        if ((dpTrace & DP_TRACE_M) >0) {
            output[i++] = tmp[DP_TRACE_M];
        }
        if ((dpTrace & DP_TRACE_I) >0) {
            output[i++] = tmp[DP_TRACE_I];
        }
        if ((dpTrace & DP_TRACE_D) >0) {
            output[i++] = tmp[DP_TRACE_D];
        }
    }
    output[i]='\0';
    printf("%-4s",output);
}

void DPDebugPrintMatch(DPWork * dpWork, DPOccurrence * dpOcc) {
    int readi,regi;
    int i,j,k;
    
    unsigned char flag[SRA_MAX_READ_LENGTH*2];
    
    printf("dpWork->regionStartIdx = %u\n",dpWork->regionStartIdx);
    printf("dpWork->regionLength   = %u\n",dpWork->regionLength);
    printf("dpOcc->ambPosition     = %llu\n",dpOcc->ambPosition);
    printf("dpOcc->matchLen        = %u\n",dpOcc->matchLen);
    printf("Occurrences start at w.r.t. ref. seq extract = %llu\n",dpOcc->ambPosition - dpWork->regionStartIdx);
    
    //Reference Sequence
    k=0;
    printf("         >  ");
    for (regi=0;regi<dpOcc->ambPosition - dpWork->regionStartIdx;) {
        printf(".");
        flag[k++]=dnaChar[dpWork->regionBuffer[regi++]];
    }
    for (i=dpOcc->matchElemsCount-1;i>=0;i--) {
        for (j=0;j<dpOcc->matchElems[i].length;j++) {
            if (dpOcc->matchElems[i].type == SRA_CHAR_ERROR_TYPE_SOFTCLIP) {
                printf("S");
                flag[k++]=dnaChar[dpWork->regionBuffer[regi++]];
            } else if (dpOcc->matchElems[i].type == SRA_CHAR_ERROR_TYPE_MISMATCH) {
                printf("M");
                flag[k++]=dnaChar[dpWork->regionBuffer[regi++]];
            } else if (dpOcc->matchElems[i].type == SRA_CHAR_ERROR_TYPE_INSERT) { //deletion on reference
                printf("D");
                flag[k++]=dnaChar[dpWork->regionBuffer[regi++]];
            } else if (dpOcc->matchElems[i].type == SRA_CHAR_ERROR_TYPE_DELETE) { //insertion on reference
                printf("I");
                flag[k++]='-';
            }
        }
    }

    for (;regi<dpWork->regionLength;) {
        printf(".");
        flag[k++]=dnaChar[dpWork->regionBuffer[regi++]];
    }
    printf("\n");
    flag[k]='\0';
    printf("REF.SEQ  =  %s\n",flag);
    
    //Read
    k=0;
    readi=0;
    for (regi=0;regi<dpOcc->ambPosition - dpWork->regionStartIdx;regi++) {
        flag[k++]=' ';
    }
    for (i=dpOcc->matchElemsCount-1;i>=0;i--) {
        for (j=0;j<dpOcc->matchElems[i].length;j++) {
            if (dpOcc->matchElems[i].type == SRA_CHAR_ERROR_TYPE_SOFTCLIP) {
                flag[k++]='*';
                readi++;
            } else if (dpOcc->matchElems[i].type == SRA_CHAR_ERROR_TYPE_MISMATCH) {
                flag[k++]=dnaChar[dpWork->readCode[readi++]];
            } else if (dpOcc->matchElems[i].type == SRA_CHAR_ERROR_TYPE_INSERT) { //deletion on reference
                flag[k++]='-';
            } else if (dpOcc->matchElems[i].type == SRA_CHAR_ERROR_TYPE_DELETE) { //insertion on reference
                flag[k++]=dnaChar[dpWork->readCode[readi++]];
            }
        }
    }
    flag[k]='\0';
    printf("   READ  =  %s\n",flag);
}

void DPDebugPrintMatrix(DPWork * dpWork) {
    unsigned int reg_i, read_i;
    unsigned int regionIdx = 0;
    unsigned char * regionBuffer = dpWork->regionBuffer;
    
    printf("Matrix M - DPScores and Traces\n");
    printf("=================================\n");
    printf("        ");
    printf("        ");
    for (read_i=0;read_i<dpWork->readLength;read_i++) {
        printf("%4c    ",dnaChar[dpWork->readCode[read_i]]);
    }
    printf("\n");
    printf("        ");
    for (reg_i=0;reg_i<=dpWork->regionLength;reg_i++) {
        if (reg_i>0) {
            printf("%8c",dnaChar[regionBuffer[reg_i-1]]);
        }
        for (read_i=0;read_i<=dpWork->readLength;read_i++) {
            int trace = dpWork->matrix[DP_MATRIX_ID_M][DP_MM_COORD(reg_i,read_i)].dpTrace;
            printf("%4d",dpWork->matrix[DP_MATRIX_ID_M][DP_MM_COORD(reg_i,read_i)].dpScore);
            DPDebugPrintTrace(trace);
        }
        printf("\n");
    }
    printf("\n\n");
    
    printf("Matrix I - DPScores and Traces\n");
    printf("=================================\n");
    printf("        ");
    printf("        ");
    for (read_i=0;read_i<dpWork->readLength;read_i++) {
        printf("%4c    ",dnaChar[dpWork->readCode[read_i]]);
    }
    printf("\n");
    printf("        ");
    for (reg_i=0;reg_i<=dpWork->regionLength;reg_i++) {
        if (reg_i>0) {
            printf("%8c",dnaChar[regionBuffer[reg_i-1]]);
        }
        for (read_i=0;read_i<=dpWork->readLength;read_i++) {
            int trace = dpWork->matrix[DP_MATRIX_ID_I][DP_MM_COORD(reg_i,read_i)].dpTrace;
            printf("%4d",dpWork->matrix[DP_MATRIX_ID_I][DP_MM_COORD(reg_i,read_i)].dpScore);
            DPDebugPrintTrace(trace);
        }
        printf("\n");
    }
    printf("\n\n");
    
    printf("Matrix D - DPScores and Traces\n");
    printf("=================================\n");
    printf("        ");
    printf("        ");
    for (read_i=0;read_i<dpWork->readLength;read_i++) {
        printf("%4c    ",dnaChar[dpWork->readCode[read_i]]);
    }
    printf("\n");
    printf("        ");
    for (reg_i=0;reg_i<=dpWork->regionLength;reg_i++) {
        if (reg_i>0) {
            printf("%8c",dnaChar[regionBuffer[reg_i-1]]);
        }
        for (read_i=0;read_i<=dpWork->readLength;read_i++) {
            int trace = dpWork->matrix[DP_MATRIX_ID_D][DP_MM_COORD(reg_i,read_i)].dpTrace;
            printf("%4d",dpWork->matrix[DP_MATRIX_ID_D][DP_MM_COORD(reg_i,read_i)].dpScore);
            DPDebugPrintTrace(trace);
        }
        printf("\n");
    }
    printf("\n\n");
    
    
    printf("Legend\n");
    printf("=================================\n");
    printf(" _____________\n");
    printf(" |  M  |  I  |\n");
    printf(" |_____|_____|\n");
    printf(" |  D  |  X  |\n");
    printf(" |_____|_____|\n");
    printf("\n");
    /*printf("        ");
    for (read_i=0;read_i<dpWork->readLength;read_i++) {
        printf("    %c",dnaChar[dpWork->readCode[read_i]]);
    }
    printf("\n");
    printf("    ");
    for (reg_i=0;reg_i<=dpWork->regionLength;reg_i++) {
        if (reg_i>0) {
            printf("    %c",dnaChar[regionBuffer[reg_i-1]]);
        }
        for (read_i=0;read_i<=dpWork->readLength;read_i++) {
            printf("%4d",dpWork->matrix[DP_MATRIX_ID_M][DP_MM_COORD(reg_i,read_i)].dpTrace);
        }
        printf("\n");
    }*/
    
}
