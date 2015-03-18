//
//    MIC-SRA2BWTOperations.c
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

#include "MIC-SRA2BWTOperations.h"

// DEBUG #1
// Un-comment below to see message showing the status 
// of an alignment step printing at the beginning and
// end of each function call.
//#define MIC_SRA_2BWT_OPERATIONS_START_END_DEBUG

__attribute__((target(mic)))
unsigned long long MICBWTReportSARange(MICSRAArguments * micArgs, unsigned long long l, unsigned long long r, uint8_t occMismatch) {
    uint8_t * outputStatus = micArgs->outputStatus;
    unsigned int * outputBlock = micArgs->outputBlock;
    MICSRAOccMetadata * metaBlock = micArgs->metaBlock;
    SRAWorkingMemory * workMem  = &(micArgs->AlgnmtMemory);
    int i;
    
    unsigned long long resultCount = 0;
    unsigned int occOutRingIdx = (*micArgs->occCount);
    unsigned int metaOutRingIdx = (*micArgs->metaCount);
    
    // resultCount is independent to the size of the MIC buffer
    // We should be counting the number of occurrences regardless
    // of the size of the MIC SRA buffer.
    resultCount = r - l + 1;
    
    if ((*outputStatus) != MIC_OUTPUT_STATUS_CLOSE) {

        // sh: ATTENTION TODO, the following two lines may overflow
        (*micArgs->occCount) += (r - l + 1);
        (*micArgs->metaCount) += (r - l + 1);
        
        if (resultCount+occOutRingIdx >= MIC_SRA_OUTPUT_MAX_ALIGNMENT) {
            (*micArgs->occCount) = 0;
            (*micArgs->metaCount) = 0;
            (*outputStatus) = MIC_OUTPUT_STATUS_CLOSE;
            workMem->IsClosed = 1;
        } else if (resultCount+metaOutRingIdx >= MIC_SRA_OUTPUT_MAX_META) {
            (*micArgs->occCount) = 0;
            (*micArgs->metaCount) = 0;
            (*outputStatus) = MIC_OUTPUT_STATUS_CLOSE;
            workMem->IsClosed = 1;
        } else {

            ///////////////////////////////////////////////
            // Text Position Per SA Index
            ///////////////////////////////////////////////
            for (;l<=r;l++) {
                
                ///////////////////////////////////////////////
                // Occurrence - Text Position
                ///////////////////////////////////////////////
                outputBlock[occOutRingIdx++] = BWTSaValue(micArgs->bwt,l) - micArgs->seedOffset - workMem->Shared_HeadTrim;

                ///////////////////////////////////////////////
                // Meta Data Per Occurrence
                ///////////////////////////////////////////////
                metaBlock[metaOutRingIdx].numPayload = 1;
                metaBlock[metaOutRingIdx].strand = micArgs->readStrand - 1;
                metaBlock[metaOutRingIdx].numOfErr = occMismatch;
                metaBlock[metaOutRingIdx].matchLen = micArgs->seedLength - workMem->Shared_HeadTrim - workMem->Shared_TailTrim;

                for (i=0;i<occMismatch&&i<MAX_NUM_OF_ERROR;i++) {
                    metaBlock[metaOutRingIdx].errors[i].type = workMem->Shared_AccErrorsPos[i].type;
                    metaBlock[metaOutRingIdx].errors[i].position = workMem->Shared_AccErrorsPos[i].position;
                }
                metaOutRingIdx++;
            }
        }
    }

    return resultCount;
}



// BWTExactModelForward_Lookup lookup your pattern in lookup table, bi-directional and assuming forward
__attribute__((target(mic)))
unsigned long long MICBWTExactModelForward_Lookup(MICSRAArguments * micArgs,
                                    CPTSRACase * alignmentCase, int stepInCase,
                                    unsigned long long * saRanges) {

    SRAWorkingMemory * workMem = &(micArgs->AlgnmtMemory);
    
    unsigned long long l_packedPattern = 0;
    unsigned long long r_packedPattern = 0;

    unsigned long long l = 0;
    unsigned long long r = 0;
    unsigned long long rev_l = 0;
    unsigned long long rev_r = 0;
    unsigned long long j;

    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    
    unsigned char * convertedKey = micArgs->readCode + micArgs->seedOffset;
    
    // Masking out any reference to occ quality interim
    // int * keyQuality = qInfo->ReadQuality;
    LT * lookupTable = micArgs->lt;
    LT * rev_lookupTable = micArgs->rlt;

    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int8_t step = alignmentCase->steps[stepInCase].step;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;
    
    // Temporarily substitute the BWT from model with the 
    // MIC BWT from the Argument.
    BWT * bwt;
    if (step<0) {
        bwt = micArgs->bwt;
    } else {
        bwt = micArgs->rev_bwt;
    }

    int pos = start;
    int k;
    int branchCount = 0;
    char branchChar = 0;
    uint8_t nbMismatch = 0;

    unsigned int lookupLength = (lookupTable->tableSize>len)?len:lookupTable->tableSize;

    #ifdef MIC_SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[MICBWTExactModelForward_Lookup] START STEP=%2d\n",stepInCase);
        printf("    - SARNGR = (%12llu-%12llu)  RANGE=%3d-%3d\n",l,r,start,end);
    #endif
    
    int i;
    for (i = 0; i <lookupLength ; i++) {
        l_packedPattern<<=LOOKUP_BIT_PER_CHAR;
        l_packedPattern |= (convertedKey[pos++]  & LOOKUP_CHAR_MASK );
    }
    r_packedPattern=l_packedPattern;
    for (i = lookupLength; i <lookupTable->tableSize ; i++) {
        l_packedPattern<<=LOOKUP_BIT_PER_CHAR;
        r_packedPattern<<=LOOKUP_BIT_PER_CHAR;
        r_packedPattern |= LOOKUP_CHAR_MASK;
    }
    l = l_packedPattern ? lookupTable->table[l_packedPattern-1]+1 : 1;
    r = lookupTable->table[r_packedPattern];

    l_packedPattern = 0;
    r_packedPattern = 0;
    for (i = 0; i <lookupLength ; i++) {
        l_packedPattern<<=LOOKUP_BIT_PER_CHAR;
        l_packedPattern |= (convertedKey[--pos] & LOOKUP_CHAR_MASK);
    }
    r_packedPattern=l_packedPattern;
    for (i = lookupLength; i <rev_lookupTable->tableSize ; i++) {
        l_packedPattern<<=LOOKUP_BIT_PER_CHAR;
        r_packedPattern<<=LOOKUP_BIT_PER_CHAR;
        r_packedPattern |= LOOKUP_CHAR_MASK;
    }
    rev_l = l_packedPattern ? rev_lookupTable->table[l_packedPattern-1]+1 : 1;
    rev_r = rev_lookupTable->table[r_packedPattern];

    if (r-l<rev_r-rev_l) {
        rev_l=rev_r-(r-l);
    } else if (r-l>rev_r-rev_l) {
        l=r-(rev_r-rev_l);
    }

    i = lookupLength;
    pos = start+lookupLength;

    while (i < len && l<=r) {
        unsigned char c = convertedKey[pos];
        BWTAllOccValue(bwt,rev_l,oL);
        BWTAllOccValue(bwt,rev_r + 1,oR);
        oCount[ALPHABET_SIZE-1]=0;
        branchCount = (oR[ALPHABET_SIZE-1]>oL[ALPHABET_SIZE-1]);
        branchChar = (ALPHABET_SIZE-1) *  branchCount;
        for (k=ALPHABET_SIZE-2;k>=0;k--) {
            oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
            
            branchCount += (oR[k]>oL[k]);
            branchChar |= k * (oR[k]>oL[k]);
        }

        // Non-branching Mismatch Handler
        if ( (errType & SRA_STEP_ERROR_TYPE_NBM) > 0 && 
             oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < micArgs->maxNBMismatch) {
            workMem->Shared_AccErrorsPos[nbMismatch].type = SRA_CHAR_ERROR_TYPE_NBMISMATCH;
            workMem->Shared_AccErrorsPos[nbMismatch].position = pos;
            rev_l = bwt->cumulativeFreq[branchChar] + oL[branchChar] + 1;
            rev_r = bwt->cumulativeFreq[branchChar] + oR[branchChar];
            r = r - oCount[branchChar];
            l = r - (rev_r-rev_l);
            nbMismatch++;
        } else {
            rev_l = bwt->cumulativeFreq[c] + oL[c] + 1;
            rev_r = bwt->cumulativeFreq[c] + oR[c];
            r = r - oCount[c];
            l = r - (rev_r-rev_l);
        }
        
        i++;
        pos++;
    }
    
    if (l<=r) {
        //Next Step
        saRanges[0] = l;
        saRanges[1] = r;
        saRanges[2] = rev_l;
        saRanges[3] = rev_r;
        if (alignmentCase->steps[stepInCase+1].type == SRA_STEP_TYPE_BACKWARD_ONLY_BWT) {
            //printf("[BWTExactModelForward_Lookup] MICBWTModelSwitchBackward\n");
            return MICBWTModelSwitchBackward(micArgs,0,0,
                                            alignmentCase,stepInCase+1,
                                            saRanges, 0, nbMismatch, 0);
        } else {
            //printf("[BWTExactModelForward_Lookup] MICBWTModelSwitchAnyDirection\n");
            return MICBWTModelSwitchAnyDirection(micArgs,0,0,
                                                alignmentCase,stepInCase+1,
                                                saRanges, 0, nbMismatch, 0);
        }

    }

    return 0;
}


// BWTExactModelBackward_Lookup lookup your pattern in lookup table, single direction and assuming backward
__attribute__((target(mic)))
unsigned long long MICBWTExactModelBackward_Lookup(MICSRAArguments * micArgs,
                                    CPTSRACase * alignmentCase, int stepInCase,
                                    unsigned long long * saRanges) {
                                    
    SRAWorkingMemory * workMem = &(micArgs->AlgnmtMemory);
    
    unsigned long long l_packedPattern = 0;
    unsigned long long r_packedPattern = 0;

    unsigned long long l = 0;
    unsigned long long r = 0;
    unsigned long long j = 0;
    
    int k = 0;
    unsigned char * convertedKey = micArgs->readCode + micArgs->seedOffset;
    
    // Masking out any reference to occ quality interim
    // int * keyQuality = qInfo->ReadQuality;
    LT * lookupTable = micArgs->lt;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int8_t step = alignmentCase->steps[stepInCase].step;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;
    
    // Temporarily substitute the BWT from model with the 
    // MIC BWT from the Argument.
    BWT * bwt;
    if (step<0) {
        bwt = micArgs->bwt;
    } else {
        bwt = micArgs->rev_bwt;
    }

    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---

    unsigned int lookupLength = (lookupTable->tableSize>len)?len:lookupTable->tableSize;
    int branchCount = 0;
    unsigned char branchChar = 0;
    uint8_t nbMismatch = 0;

    int i;
    int pos = start-lookupLength+1;

    #ifdef MIC_SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[MICBWTExactModelBackward_Lookup] START STEP=%2d\n",stepInCase);
        printf("    - SARNGR = (%12llu-%12llu)  RANGE=%3d-%3d\n",l,r,start,end);
    #endif
    
    for (i = 0; i <lookupLength ; i++) {
        l_packedPattern<<=LOOKUP_BIT_PER_CHAR;
        l_packedPattern |= (convertedKey[pos++]  & LOOKUP_CHAR_MASK );
    }

    r_packedPattern=l_packedPattern;
    for (i = lookupLength; i <lookupTable->tableSize ; i++) {
        l_packedPattern<<=LOOKUP_BIT_PER_CHAR;
        r_packedPattern<<=LOOKUP_BIT_PER_CHAR;
        r_packedPattern |= LOOKUP_CHAR_MASK;
    }
    l = l_packedPattern ? lookupTable->table[l_packedPattern-1]+1 : 1;
    r = lookupTable->table[r_packedPattern];

    i = lookupLength;
    pos = start-lookupLength;
    while (i<len && l<=r) {

        unsigned char c = convertedKey[pos];
        BWTAllOccValue(bwt,l,oL);
        BWTAllOccValue(bwt,r + 1,oR);
        
        branchCount = 0;
        branchChar = 0;
        for (k=ALPHABET_SIZE-1;k>=0;k--) {
            branchCount += (oR[k]>oL[k]);
            branchChar |= k * (oR[k]>oL[k]);
        }
        // Non-branching Mismatch Handler
        if ( (errType & SRA_STEP_ERROR_TYPE_NBM) > 0 && 
             oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < micArgs->maxNBMismatch) {
            workMem->Shared_AccErrorsPos[nbMismatch].type = SRA_CHAR_ERROR_TYPE_NBMISMATCH;
            workMem->Shared_AccErrorsPos[nbMismatch].position = pos;
            l = bwt->cumulativeFreq[branchChar] + oL[branchChar] + 1;
            r = bwt->cumulativeFreq[branchChar] + oR[branchChar];
            nbMismatch++;
        } else {
            l = bwt->cumulativeFreq[c] + oL[c] + 1;
            r = bwt->cumulativeFreq[c] + oR[c];
        }
        
        i++;
        pos--;
    }

    if (l<=r) {
        //Next Step
        saRanges[0]=l;
        saRanges[1]=r;
        return MICBWTModelSwitchBackward(micArgs,0,0,
                                alignmentCase,stepInCase+1,
                                saRanges, 0, nbMismatch, 0);
    }

    return 0;
}


// BWTExactModelBackward_Lookup lookup your pattern in lookup table, bi-directional and assuming forward
__attribute__((target(mic)))
unsigned long long MICBWTExactModelBackwardAnyDirection_Lookup(MICSRAArguments * micArgs,
                                    CPTSRACase * alignmentCase, int stepInCase,
                                    unsigned long long * saRanges) {
                                    
    SRAWorkingMemory * workMem = &(micArgs->AlgnmtMemory);
    
    unsigned long long l_packedPattern = 0;
    unsigned long long r_packedPattern = 0;

    unsigned long long l = 0;
    unsigned long long r = 0;
    unsigned long long rev_l = 0;
    unsigned long long rev_r = 0;
    unsigned long long j;

    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    
    unsigned char * convertedKey = micArgs->readCode + micArgs->seedOffset;
    
    // Masking out any reference to occ quality interim
    // int * keyQuality = qInfo->ReadQuality;
    LT * lookupTable = micArgs->lt;
    LT * rev_lookupTable = micArgs->rlt;

    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int8_t step = alignmentCase->steps[stepInCase].step;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;
    
    // Temporarily substitute the BWT from model with the 
    // MIC BWT from the Argument.
    BWT * bwt;
    if (step<0) {
        bwt = micArgs->bwt;
    } else {
        bwt = micArgs->rev_bwt;
    }
    
    int k;
    int branchCount = 0;
    char branchChar = 0;
    uint8_t nbMismatch = 0;

    unsigned int lookupLength = (lookupTable->tableSize>len)?len:lookupTable->tableSize;

    int i;
    int pos = start-lookupLength+1;
    
    #ifdef MIC_SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[MICBWTExactModelBackwardAnyDirection_Lookup] START STEP=%2d\n",stepInCase);
        printf("    - SARNGR = (%12llu-%12llu)  RANGE=%3d-%3d\n",l,r,start,end);
    #endif
    
    for (i = 0; i <lookupLength ; i++) {
        l_packedPattern<<=LOOKUP_BIT_PER_CHAR;
        l_packedPattern |= (convertedKey[pos++]  & LOOKUP_CHAR_MASK );
    }
    r_packedPattern=l_packedPattern;
    for (i = lookupLength; i <lookupTable->tableSize ; i++) {
        l_packedPattern<<=LOOKUP_BIT_PER_CHAR;
        r_packedPattern<<=LOOKUP_BIT_PER_CHAR;
        r_packedPattern |= LOOKUP_CHAR_MASK;
    }
    l = l_packedPattern ? lookupTable->table[l_packedPattern-1]+1 : 1;
    r = lookupTable->table[r_packedPattern];

    l_packedPattern = 0;
    r_packedPattern = 0;
    for (i = 0; i <lookupLength ; i++) {
        l_packedPattern<<=LOOKUP_BIT_PER_CHAR;
        l_packedPattern |= (convertedKey[--pos] & LOOKUP_CHAR_MASK);
    }
    r_packedPattern=l_packedPattern;
    for (i = lookupLength; i <rev_lookupTable->tableSize ; i++) {
        l_packedPattern<<=LOOKUP_BIT_PER_CHAR;
        r_packedPattern<<=LOOKUP_BIT_PER_CHAR;
        r_packedPattern |= LOOKUP_CHAR_MASK;
    }
    rev_l = l_packedPattern ? rev_lookupTable->table[l_packedPattern-1]+1 : 1;
    rev_r = rev_lookupTable->table[r_packedPattern];

    if (r-l<rev_r-rev_l) {
        rev_l=rev_r-(r-l);
    } else if (r-l>rev_r-rev_l) {
        l=r-(rev_r-rev_l);
    }

    i = lookupLength;
    pos = start-lookupLength;

    while (i < len && l<=r) {

        unsigned char c = convertedKey[pos];
        BWTAllOccValue(bwt,l,oL);
        BWTAllOccValue(bwt,r + 1,oR);
        oCount[ALPHABET_SIZE-1]=0;
        branchCount = (oR[ALPHABET_SIZE-1]>oL[ALPHABET_SIZE-1]);
        branchChar = (ALPHABET_SIZE-1) *  branchCount;
        for (k=ALPHABET_SIZE-2;k>=0;k--) {
            oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
            
            branchCount += (oR[k]>oL[k]);
            branchChar |= k * (oR[k]>oL[k]);
        }

        // Non-branching Mismatch Handler
        if ( (errType & SRA_STEP_ERROR_TYPE_NBM) > 0 && 
             oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < micArgs->maxNBMismatch) {
            workMem->Shared_AccErrorsPos[nbMismatch].type = SRA_CHAR_ERROR_TYPE_NBMISMATCH;
            workMem->Shared_AccErrorsPos[nbMismatch].position = pos;
            l = bwt->cumulativeFreq[branchChar] + oL[branchChar] + 1;
            r = bwt->cumulativeFreq[branchChar] + oR[branchChar];
            rev_r = rev_r - oCount[branchChar];
            rev_l = rev_r - (r-l);
            nbMismatch++;
        } else {
            l = bwt->cumulativeFreq[c] + oL[c] + 1;
            r = bwt->cumulativeFreq[c] + oR[c];
            rev_r = rev_r - oCount[c];
            rev_l = rev_r - (r-l);
        }
        
        i++;
        pos--;
    }
    
    if (l<=r) {
        //Next Step
        saRanges[0] = l;
        saRanges[1] = r;
        saRanges[2] = rev_l;
        saRanges[3] = rev_r;
        if (alignmentCase->steps[stepInCase+1].type == SRA_STEP_TYPE_BACKWARD_ONLY_BWT) {
            //printf("[BWTExactModelForward_Lookup] MICBWTModelSwitchBackward\n");
            return MICBWTModelSwitchBackward(micArgs,0,0,
                                            alignmentCase,stepInCase+1,
                                            saRanges, 0, nbMismatch, 0);
        } else {
            //printf("[BWTExactModelForward_Lookup] MICBWTModelSwitchAnyDirection\n");
            return MICBWTModelSwitchAnyDirection(micArgs,0,0,
                                                alignmentCase,stepInCase+1,
                                                saRanges, 0, nbMismatch, 0);
        }

    }

    return 0;
}


// MICBWTMismatchModelBackward matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
__attribute__((target(mic)))
unsigned long long MICBWTMismatchModelBackward(MICSRAArguments * micArgs, int i, uint8_t mismatchInserted,
                                    CPTSRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occMismatch, uint8_t nbMismatch, int occQuality) {
    
    SRAWorkingMemory * workMem = &(micArgs->AlgnmtMemory);
    
    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];

    unsigned long long saCount=0;
    int k, pos;
    
    unsigned char * convertedKey = micArgs->readCode + micArgs->seedOffset;
    
    // Masking out any reference to occ quality interim
    //int * keyQuality = micArgs->readQuality;

    uint8_t minMismatch = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxMismatch = alignmentCase->steps[stepInCase].MaxError;
    int8_t step = -1; //alignmentCase->steps[stepInCase].step;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;

    unsigned long long errSaRange[4];
    unsigned long long err_l;
    unsigned long long err_r;
    unsigned char c;

    // Temporarily substitute the BWT from model with the 
    // MIC BWT from the Argument.
    BWT * bwt;
    if (step<0) {
        bwt = micArgs->bwt;
    } else {
        bwt = micArgs->rev_bwt;
    }

    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    
    int branchCount = 0;
    char branchChar = 0;
    
    pos = start + step * i;
    len -= (minMismatch-mismatchInserted-1) * (minMismatch>(mismatchInserted+1));

    #ifdef MIC_SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[MICBWTMismatchModelBackward] START STEP=%2d IDX=%2d\n",stepInCase,i);
        printf("    - SARNGR = (%12llu-%12llu)  RANGE=%3d-%3d   ERR#=(%2d+%2d)\n",l,r,start,end,occMismatch,nbMismatch);
    #endif
    
    while (i < len && l<=r) {
    
        c = convertedKey[pos];
        BWTAllOccValue(bwt,l,oL);
        BWTAllOccValue(bwt,r + 1,oR);

        branchCount = 0;
        branchChar = 0;
        for (k=0;k<ALPHABET_SIZE;k++) {
            branchCount += (oR[k]>oL[k]);
            branchChar |= k * (oR[k]>oL[k]);
        }
        
        if (mismatchInserted<maxMismatch)
        for (k=0;k<ALPHABET_SIZE;k++) {
            if (k!=c) {
                
                errSaRange[0] = bwt->cumulativeFreq[k] + oL[k] + 1;
                errSaRange[1] = bwt->cumulativeFreq[k] + oR[k];

                if (errSaRange[0] <= errSaRange[1]) {
                    workMem->Shared_AccErrorsPos[occMismatch+mismatchInserted].type = SRA_CHAR_ERROR_TYPE_MISMATCH;
                    workMem->Shared_AccErrorsPos[occMismatch+mismatchInserted].position = pos;
                    // Masking out any reference to occ quality interim
                    //saCount+=MICBWTModelSwitchBackward(micArgs,i+1,mismatchInserted+1,
                    //                    alignmentCase,stepInCase,
                    //                    errSaRange, occMismatch, nbMismatch, occQuality+keyQuality[pos]);
                    saCount+=MICBWTModelSwitchBackward(micArgs,i+1,mismatchInserted+1,
                                        alignmentCase,stepInCase,
                                        errSaRange, occMismatch, nbMismatch, 0);
                    if (workMem->IsClosed==1) { return saCount; }
                }
            }
        }
    
        // Non-branching Mismatch Handler
        if ( (errType & SRA_STEP_ERROR_TYPE_NBM) > 0 && 
             oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < micArgs->maxNBMismatch) {
            workMem->Shared_AccErrorsPos[occMismatch+mismatchInserted+nbMismatch].type = SRA_CHAR_ERROR_TYPE_NBMISMATCH;
            workMem->Shared_AccErrorsPos[occMismatch+mismatchInserted+nbMismatch].position = pos;
            l = bwt->cumulativeFreq[branchChar] + oL[branchChar] + 1;
            r = bwt->cumulativeFreq[branchChar] + oR[branchChar];
            nbMismatch++;
        } else {
            l = bwt->cumulativeFreq[c] + oL[c] + 1;
            r = bwt->cumulativeFreq[c] + oR[c];
        }
    
        i += 1;
        pos += step;
    }

    if (mismatchInserted>=minMismatch && l<=r) {
        //Next Step
        saRanges[0] = l;
        saRanges[1] = r;
        saCount+=MICBWTModelSwitchBackward(micArgs,0,0,
                                alignmentCase,stepInCase+1,
                                saRanges, occMismatch+mismatchInserted, nbMismatch, occQuality);
    }

    return saCount;
}


// MICBWTMismatchModelAnyDirection matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
__attribute__((target(mic)))
unsigned long long MICBWTMismatchModelAnyDirection(MICSRAArguments * micArgs, int i, uint8_t mismatchInserted,
                                    CPTSRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occMismatch, uint8_t nbMismatch, int occQuality) {

    SRAWorkingMemory * workMem = &(micArgs->AlgnmtMemory);
    
    int8_t step = alignmentCase->steps[stepInCase].step;
    unsigned long long l = saRanges[(step<0)*2];
    unsigned long long r = saRanges[(step<0)*2+1];
    unsigned long long rev_l = saRanges[(step>0)*2];
    unsigned long long rev_r = saRanges[(step>0)*2+1];

    unsigned long long saCount=0;
    int k, pos;
    
    unsigned char * convertedKey = micArgs->readCode + micArgs->seedOffset;
    
    // Masking out any reference to occ quality interim
    //int * keyQuality = micArgs->readQuality;

    uint8_t minMismatch = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxMismatch = alignmentCase->steps[stepInCase].MaxError;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;

    unsigned long long err_rev_l;
    unsigned long long err_rev_r;
    unsigned long long err_r;
    unsigned long long err_l;
    unsigned long long errSaRange[4];
    unsigned char c;

    // Temporarily substitute the BWT from model with the 
    // MIC BWT from the Argument.
    BWT * bwt;
    if (step<0) {
        bwt = micArgs->bwt;
    } else {
        bwt = micArgs->rev_bwt;
    }

    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    
    int branchCount = 0;
    char branchChar = 0;
    
    pos = start + step * i;
    len -= (minMismatch-mismatchInserted-1) * (minMismatch>(mismatchInserted+1));

    #ifdef MIC_SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[MICBWTMismatchModelAnyDirection] START STEP=%2d IDX=%2d\n",stepInCase,i);
        printf("    - SARNGR = (%12llu-%12llu)  RANGE=%3d-%3d   ERR#=(%2d+%2d)\n",l,r,start,end,occMismatch,nbMismatch);
    #endif
    
    while (i < len && l<=r) {
    
        c = convertedKey[pos];
        BWTAllOccValue(bwt,rev_l,oL);
        BWTAllOccValue(bwt,rev_r + 1,oR);

        oCount[ALPHABET_SIZE-1]=0;
        branchCount = (oR[ALPHABET_SIZE-1]>oL[ALPHABET_SIZE-1]);
        branchChar = (ALPHABET_SIZE-1) *  branchCount;
        for (k=ALPHABET_SIZE-2;k>=0;k--) {
            oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
            branchCount += (oR[k]>oL[k]);
            branchChar |= k * (oR[k]>oL[k]);
        }

        if (mismatchInserted<maxMismatch)
        for (k=0;k<ALPHABET_SIZE;k++) {
            if (k!=c) {
                err_rev_l = bwt->cumulativeFreq[k] + oL[k] + 1;
                err_rev_r = bwt->cumulativeFreq[k] + oR[k];
                err_r = r - oCount[k];
                err_l = err_r - (err_rev_r-err_rev_l);

                if (err_l<=err_r) {
                    errSaRange[(step<0)*2] = err_l;
                    errSaRange[(step<0)*2+1] = err_r;
                    errSaRange[(step>0)*2] = err_rev_l;
                    errSaRange[(step>0)*2+1] = err_rev_r;
                    workMem->Shared_AccErrorsPos[occMismatch+mismatchInserted+nbMismatch].type = SRA_CHAR_ERROR_TYPE_MISMATCH;
                    workMem->Shared_AccErrorsPos[occMismatch+mismatchInserted+nbMismatch].position = pos;
                    //saCount+=BWTModelSwitchAnyDirection(aArgs,i+1,mismatchInserted+1,
                    //                        alignmentCase,stepInCase,
                    //                        errSaRange, occMismatch, nbMismatch, 
                    //                        occQuality+keyQuality[pos]);
                    saCount+=MICBWTModelSwitchAnyDirection(micArgs,i+1,mismatchInserted+1,
                                            alignmentCase,stepInCase,
                                            errSaRange, occMismatch, nbMismatch, 0);
                    if (workMem->IsClosed==1) { return saCount; }
                }
            }
        }
        
        // Non-branching Mismatch Handler
        if ( (errType & SRA_STEP_ERROR_TYPE_NBM) > 0 && 
             oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < micArgs->maxNBMismatch) {
            workMem->Shared_AccErrorsPos[occMismatch+mismatchInserted+nbMismatch].type = SRA_CHAR_ERROR_TYPE_NBMISMATCH;
            workMem->Shared_AccErrorsPos[occMismatch+mismatchInserted+nbMismatch].position = pos;
            rev_l = bwt->cumulativeFreq[branchChar] + oL[branchChar] + 1;
            rev_r = bwt->cumulativeFreq[branchChar] + oR[branchChar];
            r = r - oCount[branchChar];
            l = r - (rev_r-rev_l);
            nbMismatch++;
        } else {
            rev_l = bwt->cumulativeFreq[c] + oL[c] + 1;
            rev_r = bwt->cumulativeFreq[c] + oR[c];
            r = r - oCount[c];
            l = r - (rev_r-rev_l);
        }
        
        i += 1;
        pos += step;
    }

    if (mismatchInserted>=minMismatch && l<=r) {
        //Next Step
        saRanges[(step<0)*2] = l;
        saRanges[(step<0)*2+1] = r;
        saRanges[(step>0)*2] = rev_l;
        saRanges[(step>0)*2+1] = rev_r;
        saCount+=MICBWTModelSwitchAnyDirection(micArgs,0,0,
                                alignmentCase,stepInCase+1,
                                saRanges, occMismatch+mismatchInserted, nbMismatch, occQuality);
    }

    return saCount;
}



// MICBWTOptionalModelAnyDirection matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
// It returns whenever the alignment reaches it ends (before the end of BWT search or end of read)
__attribute__((target(mic)))
unsigned long long MICBWTOptionalModelAnyDirection(MICSRAArguments * micArgs, int i, uint8_t errorInserted,
                                    CPTSRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occError, uint8_t nbMismatch, int occQuality) {
                                    
    SRAWorkingMemory * workMem = &(micArgs->AlgnmtMemory);
    
    int8_t step = alignmentCase->steps[stepInCase].step;
    unsigned long long l = saRanges[(step<0)*2];
    unsigned long long r = saRanges[(step<0)*2+1];
    unsigned long long rev_l = saRanges[(step>0)*2];
    unsigned long long rev_r = saRanges[(step>0)*2+1];

    unsigned long long saCount=0;
    int k, pos;
    
    unsigned char * convertedKey = micArgs->readCode + micArgs->seedOffset;
    
    // Masking out any reference to occ quality interim
    //int * keyQuality = micArgs->readQuality;

    uint8_t minMismatch = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxMismatch = alignmentCase->steps[stepInCase].MaxError;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;

    unsigned char c;

    // Temporarily substitute the BWT from model with the 
    // MIC BWT from the Argument.
    BWT * bwt;
    if (step<0) {
        bwt = micArgs->bwt;
    } else {
        bwt = micArgs->rev_bwt;
    }

    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int oL[ALPHABET_SIZE];
    unsigned int oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    
    int branchCount = 0;
    char branchChar = 0;
    int trimmedLength = 0; 
    
    pos = start + step * i;
    len -= (minMismatch-errorInserted-1) * (minMismatch>(errorInserted+1));

    #ifdef MIC_SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[MICBWTOptionalModelAnyDirection] START STEP=%2d IDX=%2d\n",stepInCase,i);
        printf("    - SARNGR = (%12llu-%12llu)  RANGE=%3d-%3d   ERR#=(%2d+%2d)\n",l,r,start,end,occError,nbMismatch);
    #endif
    
    while (i < len && l<=r) {
    
        c = convertedKey[pos];
        BWTAllOccValue(bwt,rev_l,oL);
        BWTAllOccValue(bwt,rev_r + 1,oR);

        oCount[ALPHABET_SIZE-1]=0;
        branchCount = (oR[ALPHABET_SIZE-1]>oL[ALPHABET_SIZE-1]);
        branchChar = (ALPHABET_SIZE-1) *  branchCount;
        for (k=ALPHABET_SIZE-2;k>=0;k--) {
            oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
            branchCount += (oR[k]>oL[k]);
            branchChar |= k * (oR[k]>oL[k]);
        }

        // Non-branching Mismatch Handler
        if ( (errType & SRA_STEP_ERROR_TYPE_NBM) > 0 && 
             oL[c] + 1 > oR[c] ) {
            if ( branchCount == 1 && nbMismatch < micArgs->maxNBMismatch) {
                workMem->Shared_AccErrorsPos[occError+errorInserted+nbMismatch].type = SRA_CHAR_ERROR_TYPE_NBMISMATCH;
                workMem->Shared_AccErrorsPos[occError+errorInserted+nbMismatch].position = pos;
                rev_l = bwt->cumulativeFreq[branchChar] + oL[branchChar] + 1;
                rev_r = bwt->cumulativeFreq[branchChar] + oR[branchChar];
                r = r - oCount[branchChar];
                l = r - (rev_r-rev_l);
                nbMismatch++;
            } else {
                // In case we are not allowed any more NBM or 
                // NBM is not possible we give up.
                break;
            }
        } else {
            rev_l = bwt->cumulativeFreq[c] + oL[c] + 1;
            rev_r = bwt->cumulativeFreq[c] + oR[c];
            r = r - oCount[c];
            l = r - (rev_r-rev_l);
        }
        
        i += 1;
        pos += step;
    }
    trimmedLength = len - i; 

    if (errorInserted>=minMismatch && l<=r) {
    
        if (step==-1) {
            //Update the Shared_HeadTrim if it's backward search.
            //This is essential for the alignment position calculation.
            //Optional region only have an effect on the alignment position only if it is
            //at the head of the read and backward was performed on that region. In the SRA2BWTMdl we
            //assume optional region only exists at the head/tail and tail alignment have no impact on
            //the alignment position.
            workMem->Shared_HeadTrim = trimmedLength;

            for (i=0;i<occError+errorInserted+nbMismatch&&i<MAX_NUM_OF_ERROR;i++) {
                workMem->Shared_AccErrorsPos[i].position -= trimmedLength;
            }
        } else {
            //Update the Shared_TailTrim if it's forward search.
            //This is essential for the alignment match length calculation. In the SRA2BWTMdl we
            //assume optional region only exists at the head/tail and tail alignment have no impact on
            //the alignment position.
            workMem->Shared_TailTrim = trimmedLength;
        }

        //Next Step
        saRanges[(step<0)*2] = l;
        saRanges[(step<0)*2+1] = r;
        saRanges[(step>0)*2] = rev_l;
        saRanges[(step>0)*2+1] = rev_r;
        saCount+=MICBWTModelSwitchAnyDirection(micArgs,0,0,
                                alignmentCase,stepInCase+1,
                                saRanges, occError+errorInserted, nbMismatch, occQuality);
    }

    return saCount;
}



__attribute__((target(mic)))
unsigned long long MICBWTModelSwitchBackward(MICSRAArguments * micArgs,  int i, uint8_t errorInserted,
                                    CPTSRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occError, uint8_t nbMismatch, int occQuality) {
                                    
    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];
    unsigned long long j;
    
    #ifdef MIC_SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[MICBWTModelSwitchBackward] START STEP=%2d\n",stepInCase);
        printf("    - SARNGR = (%12llu-%12llu)\n",l,r);
    #endif
    
    if (micArgs->AlgnmtMemory.IsClosed==1) {
        #ifdef MIC_SRA_2BWT_OPERATIONS_START_END_DEBUG
            printf("[MICBWTModelSwitchBackward] EXIT NONE\n\n");
        #endif
        return 0;
    }
    
    if (alignmentCase->steps[stepInCase].type == SRA_STEP_TYPE_COMPLETE) {
        #ifdef MIC_SRA_2BWT_OPERATIONS_START_END_DEBUG
            printf("[MICBWTModelSwitchBackward] EXIT = (%12llu-%12llu)  ERR#=(%2d+%2d)\n\n",l,r,occError,nbMismatch);
        #endif
        // ATTENTION, we need a MIC output function to write into output block
        //return OCCReportSARange(micArgs,l,r,occError+errorInserted+nbMismatch,occQuality);
        return MICBWTReportSARange(micArgs,l,r,occError+errorInserted);
    }

    uint8_t errorType = alignmentCase->steps[stepInCase].ErrorType;
    uint8_t maxError = alignmentCase->steps[stepInCase].MaxError;

    if ((errorType & SRA_STEP_ERROR_TYPE_NO_ERROR) > 0) {
        return MICBWTMismatchModelBackward(micArgs,  i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    } else if ((errorType & SRA_STEP_ERROR_TYPE_EDIT_DISTANCE) > 0 || (errorType & SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY) > 0) {
        //Temporarily comment out Edit-Distance alignment function
        //return BWTEditModelBackward(micArgs, i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    } else if (errorInserted>=maxError) {
        return MICBWTMismatchModelBackward(micArgs,  i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    } else if ((errorType & SRA_STEP_ERROR_TYPE_MISMATCH_ONLY) > 0) {
        return MICBWTMismatchModelBackward(micArgs, i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    }

    printf("[BWTModelSwitchBackward] Unknown error type[%u] encountered\n",errorType);
    return 0;
}

__attribute__((target(mic)))
unsigned long long MICBWTModelSwitchAnyDirection(MICSRAArguments * micArgs,  int i, uint8_t errorInserted,
                                    CPTSRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occError, uint8_t nbMismatch, int occQuality) {
                                    
    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];
    
    #ifdef MIC_SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[MICBWTModelSwitchAnyDirection] START STEP=%2d\n",stepInCase);
        printf("    - SARNGR = (%12llu-%12llu)\n",l,r);
    #endif
    
    if (micArgs->AlgnmtMemory.IsClosed==1) {
        #ifdef MIC_SRA_2BWT_OPERATIONS_START_END_DEBUG
            printf("[MICBWTModelSwitchAnyDirection] EXIT NONE\n\n");
        #endif
        return 0;
    }
    
    if (alignmentCase->steps[stepInCase].type == SRA_STEP_TYPE_COMPLETE) {
        #ifdef MIC_SRA_2BWT_OPERATIONS_START_END_DEBUG
            printf("[MICBWTModelSwitchAnyDirection] EXIT = (%12llu-%12llu)  ERR#=(%2d+%2d)\n\n",l,r,occError,nbMismatch);
        #endif
        // ATTENTION, we need a MIC output function to write into output block
        //return OCCReportSARange(micArgs,l,r,occError+errorInserted+nbMismatch,occQuality);
        return MICBWTReportSARange(micArgs,l,r,occError+errorInserted);
    }

    uint8_t errorType = alignmentCase->steps[stepInCase].ErrorType;
    uint8_t maxError = alignmentCase->steps[stepInCase].MaxError;

    if ((errorType & SRA_STEP_ERROR_TYPE_NO_ERROR) > 0) {
        return MICBWTMismatchModelAnyDirection(micArgs,  i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    } else if ((errorType & SRA_STEP_ERROR_TYPE_EDIT_DISTANCE) > 0 || (errorType & SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY) > 0) {
        //Temporarily comment out Edit-Distance alignment function
        //return BWTEditModelBackward(micArgs, i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    } else if ((errorType & SRA_STEP_ERROR_TYPE_OPTIONAL_NBM) > 0) {
        return MICBWTOptionalModelAnyDirection(micArgs,  i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    } else if (errorInserted>=maxError) {
        return MICBWTMismatchModelAnyDirection(micArgs,  i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    } else if ((errorType & SRA_STEP_ERROR_TYPE_MISMATCH_ONLY) > 0) {
        return MICBWTMismatchModelAnyDirection(micArgs, i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    }

    printf("[MICBWTModelSwitchAnyDirection] Unknown error type[%u] encountered\n",errorType);
    return 0;
}
