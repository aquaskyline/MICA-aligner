//
//    MIC-DPCore.h
//
//    MICA 
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

   Date   : 18th September 2013
   Author : shchan 
   Change : Added soft clipping

   Date   : 4th July 2013
   Author : shchan 
   Change : New file.

*/
/////////////////////////////////////////////////////

#include "2bwt-flex/DPCore.h"

#ifndef __MIC_DPCORE_H__
#define __MIC_DPCORE_H__

#define MIC_DP_TRACE_M              1
#define MIC_DP_TRACE_I              2
#define MIC_DP_TRACE_D              3
#define MIC_DP_TRACE_S              0

#define MIC_DP_MAXREGIONLENGTH      2048  // must be multiple of 16
#define MIC_DP_MAXREADLENGTH        256  // must be multiple of 16

#define coord(i,j)                  (((i+j)*MIC_DP_MAXREADLENGTH) + (16*(i+j)) + (j+15))

#define MIC_DP_MATRIXLENGTH         coord(MIC_DP_MAXREGIONLENGTH, MIC_DP_MAXREADLENGTH)
#define MIC_DP_MAXTRACELENGTH       (MIC_DP_MAXREADLENGTH*2)

typedef struct DPParametersMIC {
    int scoreMatch, scoreMismatch, scoreGapOpen, scoreGapExtend;
    int maxLeftClip, maxRightClip, maxTotalClip;
    double scoreThreshold;
} DPParametersMIC;

typedef struct DPArgumentsMIC {
    char *regCode, *readCode;
    unsigned regLength, readLength;
} DPArgumentsMIC;

typedef struct DPMatrixCellMIC {
    unsigned scoreM   : 11;
    unsigned scoreI   : 11;
    unsigned scoreD   : 5;
    unsigned codeRef  : 2;
    unsigned codeRead : 2;
    unsigned flag     : 1;
} DPMatrixCellMIC;

typedef struct DPWorkMIC {

    // working memory
    DPMatrixCellMIC matrix[MIC_DP_MATRIXLENGTH];
    char regBuffer[MIC_DP_MAXREGIONLENGTH];

    // output
    char trace[MIC_DP_MAXTRACELENGTH];
    int traceLength;
    int maxScore, maxScoreStartPos, maxScoreEndPos;
    double scoreThreshold;

    // constant for calculations
    int scoreMatch, scoreMismatch, scoreGapOpen, scoreGapExtend;
    int maxLeftClip, maxRightClip, maxTotalClip;
    __m512i rMatch, rMismatch, rGapOpen, rGapExtend;
    __m512i rMaskTraceM, rMaskTraceI, rMaskTraceD, rMaskCodeRef, rMaskCodeRead, rMaskFlag, rMaskScoreM, rMaskScoreI, rMaskScoreD;

} DPWorkMIC;

__attribute__((target(mic)))
void DPWorkInitMIC(const DPParametersMIC* dpp, DPWorkMIC* dpw);

__attribute__((target(mic)))
int DPMatrixFillMIC(const DPArgumentsMIC* dpa, DPWorkMIC* dpw);

#endif
