//
//    DPCore.c
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
   
   Date   : 26th February 2012
   Author : Edward MK Wu
   Change : New file.

*/
/////////////////////////////////////////////////////

#include "DPCore.h"

void DPOccurrenceCopy(DPOccurrence * source, DPOccurrence * destination) {
    int i;
    destination->ambPosition                = source->ambPosition;
    destination->strand                     = source->strand;
    destination->matchElemsCount            = source->matchElemsCount;
    destination->mismatchCount              = source->mismatchCount;
    destination->matchLen                   = source->matchLen;
    destination->matchScore                 = source->matchScore;
    
    for (i=0;i<source->matchElemsCount&&i<DP_MAX_ERROR_COUNT;i++) {
        destination->matchElems[i].type     = source->matchElems[i].type;
        destination->matchElems[i].length   = source->matchElems[i].length;
    }
}

DPWork * DPWorkCreate() {
    int i,k;
    
    DPWork * dpWork = (DPWork *) malloc(sizeof(DPWork));
    dpWork->flag = DP_ALIGNMENT_INITIALISED;
    dpWork->leadVerifyLength = 0;
    dpWork->tailVerifyLength = 0;
    dpWork->leadSoftClipLength = 0;
    dpWork->tailSoftClipLength = 0;
    dpWork->leadHardClipLength = 0;
    dpWork->tailHardClipLength = 0;
    
    return dpWork;
}

void DEBUGDPPrintWork(DPWork * dpWork) {
    printf("[DEBUGDPPrintWork] Printing DP Work structure..\n");
    printf("  FLAG     = %d\n",dpWork->flag);
    int i;
    printf("  REGBUF   = ");
    for (i=0;i<dpWork->regionLength;i++) {
        printf("%c",dnaChar[dpWork->_regionBuffer[i]]);
    } printf(" (%u) \n",dpWork->regionLength);
    printf("  READ     = ");
    for (i=0;i<dpWork->readLength;i++) {
        printf("%c",dnaChar[dpWork->readCode[i]]);
    } printf(" (%u) \n",dpWork->readLength);
    
    printf("  LEAD-VL  = %d\n",dpWork->leadVerifyLength);
    printf("  TAIL-VL  = %d\n",dpWork->tailVerifyLength);
    
    printf("  LEAD-SCL = %d\n",dpWork->leadSoftClipLength);
    printf("  TAIL-SCL = %d\n",dpWork->tailSoftClipLength);
    printf("  TTL-SCL  = %d\n",dpWork->totalSoftClipLength);
    
    printf("  LEAD-HC = %d\n",dpWork->leadHardClipLength);
    printf("  TAIL-HC = %d\n",dpWork->tailHardClipLength);
}

void DPWorkInitialise(DPWork * dpWork) {
    dpWork->flag = DP_ALIGNMENT_INITIALISED;    
}

void DPWorkFree(DPWork * dpWork) {
    free(dpWork);
}

DPScores * DPScoresConstruct() {
    DPScores * dpScores = (DPScores*) malloc (sizeof(DPScores));
    
    dpScores->dpMatch     = 0;
    dpScores->dpMismatch  = 0;
    dpScores->dpGapOpen   = 0;
    dpScores->dpGapExtend = 0;
    
    return dpScores;
}

void DPScoresFree(DPScores * dpScores) {
    free(dpScores);
}
