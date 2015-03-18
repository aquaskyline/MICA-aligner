//
//    MIC-DPCore.c
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

#include "MIC-DPCore.h"

__attribute__((target(mic)))
void DPWorkInitMIC(const DPParametersMIC * dpp, DPWorkMIC * dpw) {

#ifdef __MIC__
    dpw->maxLeftClip   = dpp->maxLeftClip;
    dpw->maxRightClip  = dpp->maxRightClip;
    dpw->maxTotalClip  = dpp->maxTotalClip;

    dpw->scoreThreshold = dpp->scoreThreshold;

    dpw->scoreMatch    = dpp->scoreMatch;
    dpw->scoreMismatch = dpp->scoreMismatch;
    dpw->scoreGapOpen  = dpp->scoreGapOpen;
    dpw->scoreGapExtend= dpp->scoreGapExtend;

    dpw->rMatch        = _mm512_set1_epi32(dpp->scoreMatch);
    dpw->rMismatch     = _mm512_set1_epi32(dpp->scoreMismatch);
    dpw->rGapOpen      = _mm512_set1_epi32(dpp->scoreGapOpen);
    dpw->rGapExtend    = _mm512_set1_epi32(dpp->scoreGapExtend);

    dpw->rMaskScoreM   = _mm512_set1_epi32(0x000007ff);
    dpw->rMaskScoreI   = _mm512_set1_epi32(0x003ff800);
    dpw->rMaskScoreD   = _mm512_set1_epi32(0x07c00000);
    dpw->rMaskCodeRef  = _mm512_set1_epi32(0x18000000);
    dpw->rMaskCodeRead = _mm512_set1_epi32(0x60000000);
    dpw->rMaskFlag     = _mm512_set1_epi32(0x80000000);
#endif

}


#ifdef __MIC__

__attribute__((target(mic)))
static inline __m512i matrixLoad(void * m, int i){
    uint32_t* p = m;
    //if (i%16==0) return  _mm512_load_epi32(p+i);
    //__m512i r = _mm512_load_epi32(p+i+1);
    //return _mm512_alignr_epi32(r,_mm512_set1_epi32(p[i]),15); 

    //return _mm512_set_epi32(p[i],p[i+1],p[i+2],p[i+3],p[i+4],p[i+5],p[i+6],p[i+7],p[i+8],p[i+9],p[i+10],p[i+11],p[i+12],p[i+13],p[i+14],p[i+15]);
    return _mm512_set_epi32(p[i+15],p[i+14],p[i+13],p[i+12],p[i+11],p[i+10],p[i+9],p[i+8],p[i+7],p[i+6],p[i+5],p[i+4],p[i+3],p[i+2],p[i+1],p[i]);
}

// debug function
__attribute__((target(mic)))
void print512(__m512i var){
   // uint32_t *val = (uint32_t*) &var;
    //printf("Numerical: %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u %u \n", 
   //        val[15], val[14], val[13], val[12], val[11], val[10], val[9], val[8], val[7], val[6], val[5], val[4], val[3], val[2], val[1], val[0]);

    
    DPMatrixCellMIC *c = (DPMatrixCellMIC*) &var;
    printf("%u %u %u\n", c[0].scoreM, c[0].scoreI, c[0].scoreD);
}

__attribute__((target(mic)))
static inline void matrixFill(DPWorkMIC * dpw, int i, int j){

    _mm_prefetch((char*)dpw->matrix+coord(i,j), _MM_HINT_T0);

    // A ^<
    // B <
    // C ^

    // Raw data in A, B, C 
    __m512i rA, rB, rC;
    __m512i rM, rI, rD;

    __m512i rR;
    rR = _mm512_setzero_epi32(); 
   
    rA = matrixLoad(dpw->matrix,(coord(i-1,j-1)));
    rB = matrixLoad(dpw->matrix,coord(i,j-1));
    rC = matrixLoad(dpw->matrix,coord(i-1,j));

    __m512i rcB = _mm512_and_epi32(rB,dpw->rMaskCodeRef);
    __m512i rcC = _mm512_and_epi32(rC,dpw->rMaskCodeRead);

    rR = _mm512_or_epi32(rR, rcB);
    rR = _mm512_or_epi32(rR, rcC);
    
    __m512i r1;

   // set scoreM and traceM
    rM = _mm512_and_epi32(rA, dpw->rMaskScoreM);
    rI = _mm512_srli_epi32(_mm512_and_epi32(rA, dpw->rMaskScoreI),11);
    rD = _mm512_sub_epi32(_mm512_add_epi32(rM, _mm512_srli_epi32(_mm512_and_epi32(rA, dpw->rMaskScoreD),22)), _mm512_set1_epi32(16));

    // r1 = max of rM, rI, rD
    r1 = _mm512_mask_mov_epi32(rI,_mm512_cmpge_epi32_mask(rM,rI),rM);
    r1 = _mm512_mask_mov_epi32(rD,_mm512_cmpge_epi32_mask(r1,rD),r1);

    // set scoreM
    rR = _mm512_or_epi32(rR, _mm512_add_epi32(r1, _mm512_mask_mov_epi32(dpw->rMismatch, _mm512_cmpeq_epi32_mask(_mm512_srli_epi32(rcC,2),rcB), dpw->rMatch)));

    // set scoreI and traceI
    rM = _mm512_add_epi32(dpw->rGapOpen, _mm512_and_epi32(rB, dpw->rMaskScoreM));
    rI = _mm512_add_epi32(dpw->rGapExtend, _mm512_srli_epi32(_mm512_and_epi32(rB, dpw->rMaskScoreI),11));
    rD = _mm512_add_epi32(dpw->rGapOpen, _mm512_sub_epi32(_mm512_add_epi32(rM, _mm512_srli_epi32(_mm512_and_epi32(rB, dpw->rMaskScoreD),22)), _mm512_set1_epi32(16)));
    
    // r1 = max of rM, rI, rD
    r1 = _mm512_mask_mov_epi32(rI,_mm512_cmpge_epi32_mask(rM,rI),rM);
    r1 = _mm512_mask_mov_epi32(rD,_mm512_cmpge_epi32_mask(r1,rD),r1);
  
    // set scoreI
    rR = _mm512_or_epi32(rR, _mm512_slli_epi32(r1, 11));

    // set scoreD and traceD
    rM = _mm512_add_epi32(dpw->rGapOpen, _mm512_and_epi32(rC, dpw->rMaskScoreM));
    rI = _mm512_add_epi32(dpw->rGapOpen, _mm512_srli_epi32(_mm512_and_epi32(rC, dpw->rMaskScoreI),11));
    rD = _mm512_add_epi32(dpw->rGapExtend, _mm512_sub_epi32(_mm512_add_epi32(_mm512_and_epi32(rC, dpw->rMaskScoreM), _mm512_srli_epi32(_mm512_and_epi32(rC, dpw->rMaskScoreD),22)), _mm512_set1_epi32(16)));
    
    // r1 = max of rM, rI, rD
    r1 = _mm512_mask_mov_epi32(rI,_mm512_cmpge_epi32_mask(rM,rI),rM);
    r1 = _mm512_mask_mov_epi32(rD,_mm512_cmpge_epi32_mask(r1,rD),r1);
  
    // set scoreD
    rR = _mm512_or_epi32(rR, _mm512_slli_epi32(_mm512_sub_epi32(_mm512_add_epi32(r1, _mm512_set1_epi32(16)), _mm512_and_epi32(rR,dpw->rMaskScoreM)),22));

    _mm512_store_epi32(dpw->matrix+coord(i,j),rR);
}


// backtrack and store information in the output fields of dpw
__attribute__((target(mic)))
static int backtrackMatrix(const DPArgumentsMIC * dpa, DPWorkMIC * dpw){

    int max = -999999;
    int reg_i, max_reg_i, read_i, max_read_i;
    int i , j, c, t = 0;

    // find the maximum entry in the "potential-to-be-clipped" region
    for (reg_i=1; reg_i<dpa->regLength; reg_i++){
        for (read_i=dpa->readLength-dpw->maxRightClip; read_i<=dpa->readLength; read_i++){
            if (dpw->matrix[coord(reg_i,read_i)].scoreM>max){
                max_reg_i = reg_i;
                max_read_i = read_i;
                t = MIC_DP_TRACE_M;
                max = dpw->matrix[coord(reg_i,read_i)].scoreM;
            }
            if (dpw->matrix[coord(reg_i,read_i)].scoreI>max){
                max_reg_i = reg_i;
                max_read_i = read_i;
                t = MIC_DP_TRACE_I;
                max = dpw->matrix[coord(reg_i,read_i)].scoreI;
            }
        }
    }

    int clipLength = dpa->readLength - max_read_i;
    int hLimit = dpw->maxTotalClip-clipLength<dpw->maxRightClip?dpw->maxTotalClip-clipLength:dpw->maxRightClip;

    // not work for other scores
    // and no visible improvement, so disabled
    //if (max-1000+hLimit < dpw->scoreThreshold * dpa->readLength) return FALSE;

    dpw->traceLength = 0;
    for (i = 0; i<clipLength; i++){
        dpw->trace[dpw->traceLength] = MIC_DP_TRACE_S;
        dpw->traceLength++;
    }

    i = max_reg_i;
    j = max_read_i;
    int score = max;

    int clipScore = 99999;
    int clipI, clipJ, clipMod, clipTrace = 0;

    while(j>0 && i>=0){

        if (t==MIC_DP_TRACE_M){
            i--;
            j--;
            if (score == (dpw->matrix[coord(i,j)].scoreM + (dpa->regCode[i]==dpa->readCode[j]?dpw->scoreMatch:dpw->scoreMismatch))) c = MIC_DP_TRACE_M;
            else if (score == (dpw->matrix[coord(i,j)].scoreI + (dpa->regCode[i]==dpa->readCode[j]?dpw->scoreMatch:dpw->scoreMismatch))) c = MIC_DP_TRACE_I;
            else if (score == (dpw->matrix[coord(i,j)].scoreM + dpw->matrix[coord(i,j)].scoreD - 16 + (dpa->regCode[i]==dpa->readCode[j]?dpw->scoreMatch:dpw->scoreMismatch))) c = MIC_DP_TRACE_D;
            else {
               printf("[DPCore] Unexpected Backtrack Error (M)\n");
               return FALSE;
            }
        }
        else if (t==MIC_DP_TRACE_I){
            j--;
            if (score == (dpw->matrix[coord(i,j)].scoreM+dpw->scoreGapOpen)) c = MIC_DP_TRACE_M;
            else if (score == (dpw->matrix[coord(i,j)].scoreI+dpw->scoreGapExtend)) c = MIC_DP_TRACE_I;
            else if (score == (dpw->matrix[coord(i,j)].scoreM + dpw->matrix[coord(i,j)].scoreD - 16 +dpw->scoreGapOpen)) c = MIC_DP_TRACE_D;
            else {
               printf("[DPCore] Unexpected Backtrack Error (I)\n");
               return FALSE;
            }
        }
        else if (t==MIC_DP_TRACE_D){
            i--;
            if (score == (dpw->matrix[coord(i,j)].scoreM+dpw->scoreGapOpen)) c = MIC_DP_TRACE_M;
            else if (score == (dpw->matrix[coord(i,j)].scoreI+dpw->scoreGapOpen )) c = MIC_DP_TRACE_I;
            else if (score == (dpw->matrix[coord(i,j)].scoreM + dpw->matrix[coord(i,j)].scoreD - 16 +dpw->scoreGapExtend)) c = MIC_DP_TRACE_D;
            else {
               printf("[DPCore] Unexpected Backtrack Error (D)\n");
               return FALSE;
            }
        }
        else {
            printf("[DPCore] Backtrack Trace Error\n");
            return FALSE;
        }

        if (c==MIC_DP_TRACE_M) score = (dpw->matrix[coord(i,j)].scoreM);
        else if (c==MIC_DP_TRACE_I) score = (dpw->matrix[coord(i,j)].scoreI);
        else if (c==MIC_DP_TRACE_D) score = (dpw->matrix[coord(i,j)].scoreM + dpw->matrix[coord(i,j)].scoreD - 16);

        dpw->trace[dpw->traceLength] = t;
        clipMod = (t==2 && c==2)*2; 
        t = c;

        if (t==1){
            if (j<=hLimit && dpw->matrix[coord(i,j)].scoreM-1000<clipScore) {
                clipScore = dpw->matrix[coord(i,j)].scoreM-1000;
                clipI = i;
                clipJ = j;
                clipTrace = dpw->traceLength;
            }
        }
        else if (t==2){
            if (j<=hLimit && dpw->matrix[coord(i,j)].scoreI-1000+clipMod<clipScore) {
                clipScore = dpw->matrix[coord(i,j)].scoreI-1000+clipMod;
                clipI = i;
                clipJ = j;
                clipTrace = dpw->traceLength;
            }
        }
        else if (t==3){
            if (j<=hLimit && dpw->matrix[coord(i,j)].scoreM+dpw->matrix[coord(i,j)].scoreD-1016<clipScore) {
                clipScore = dpw->matrix[coord(i,j)].scoreM+dpw->matrix[coord(i,j)].scoreD-1016;
                clipI = i;
                clipJ = j;
                clipTrace = dpw->traceLength;
            }
        }
        else {
            printf("[DPCore] Backtrack Trace Error\n");
            return FALSE;
        }

        dpw->traceLength++;
    }
        //printf("\n");
     //
    for (j = clipTrace+1; j<dpw->traceLength; j++) dpw->trace[j] = MIC_DP_TRACE_S;
    dpw->traceLength = clipTrace+clipJ+1;
 
    dpw->maxScore = max-1000-clipScore;

    dpw->maxScoreStartPos = clipI;
    dpw->maxScoreEndPos = max_reg_i;

    if (dpw->maxScore < dpw->scoreThreshold * dpa->readLength) return FALSE;
    return TRUE; 
    //printf("tracelength: %d maxscore: %d, pos %d\n", dpw->traceLength, max-1000, max_reg_i);
    //printf("%d %d cliped: %d\n", max, i, clipLength);
}

__attribute__((target(mic)))
int DPMatrixFillMIC(const DPArgumentsMIC * dpa, DPWorkMIC * dpw) {

    int i, j, k;

    dpw->matrix[coord(0,0)].scoreM  = 1000;
    dpw->matrix[coord(0,0)].scoreI  = 300; 
    dpw->matrix[coord(0,0)].scoreD  = 8;

    for (j=1; j<=dpa->readLength; j+=16){
        for (k=0; k<16; k++) {
            dpw->matrix[coord(0,j)+k].scoreM   = 300;
            dpw->matrix[coord(0,j)+k].scoreI   = (k==0)?998+j*dpw->scoreGapExtend:300;
            dpw->matrix[coord(0,j)+k].scoreD   = 16;
            dpw->matrix[coord(0,j)+k].codeRead = (j+k<=dpa->readLength)?(dpa->readCode[j+k-1]):0;
        }
    }

    for (i=1; i<=dpa->regLength+dpa->readLength; i++){
        dpw->matrix[coord(i,0)].scoreM  = (i<=dpa->regLength)?1000:300;
        dpw->matrix[coord(i,0)].scoreI  = 300;
        dpw->matrix[coord(i,0)].scoreD  = 8;
        dpw->matrix[coord(i,0)].codeRef = (i<=dpa->regLength)?dpa->regCode[i-1]:0;
        for (j=1; ((j<i) && (j<=dpa->readLength)); j+=16){
            matrixFill(dpw,i-j,j);
        }
        dpw->matrix[coord(0,i)].scoreM   = 300;
        dpw->matrix[coord(0,i)].scoreI   = 998+i*dpw->scoreGapExtend;
        dpw->matrix[coord(0,i)].scoreD   = 16;
        dpw->matrix[coord(0,i)].codeRead = (i<=dpa->readLength)?(dpa->readCode[i-1]):0;
    }
    return backtrackMatrix(dpa, dpw);
}
#else

// OBSOLETE
// "must-be-correct"  version to verify the solution of the MIC version
int DPMatrixFillMIC(const DPArgumentsMIC * dpa, DPWorkMIC * dpw) {
    return 0;
}

#endif
