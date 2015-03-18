//
//    SRA2BWTOperations.c
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
   
   Date   : 8th May 2011
   Author : Edward MK Wu
   Change : New file.

   Date   : 8th February 2011
   Author : Edward MK Wu
   Change : Fixed segmentation fault in HSPLoad and CE func.

*/
/////////////////////////////////////////////////////

#include "SRA2BWTOperations.h"

// DEBUG #1
// Un-comment below to see message showing the status 
// of an alignment step printing at the beginning and
// end of each function call.
//#define SRA_2BWT_OPERATIONS_START_END_DEBUG

// DEBUG #2
// Un-comment below to see message showing the status 
// of an alignment step printing at each step in the case.
//#define SRA_2BWT_OPERATIONS_STEP_DEBUG

void SRAArgumentsPrint(SRAArguments * aArgs) {
    int i;
    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    
    printf("ReadId =                  %llu\n",qInfo->ReadId);
    printf("ReadName =                %s\n",qInfo->ReadName);
    printf("ReadLength =              %llu\n",qInfo->ReadLength);
    printf("ReadCode =                "); for (i=0;i<qInfo->ReadLength;i++) {printf("%c",dnaChar[qInfo->ReadCode[i]]);}printf("\n");
    
    printf("Setting->ReadStrand =     %d\n",aSetting->ReadStrand);
    printf("Setting->MaxError =       %d\n",aSetting->MaxError);
    printf("Setting->ErrorType =      %d\n",aSetting->ErrorType);
    
    printf("Setting->OutputType =     %d\n",aSetting->OutputType);
    printf("Setting->OutputFormat =  %d\n",aSetting->OutputFormat);
}


// BWTxxxModelxxx functions are fully generalised BWT search algorithm that searchs for reads contain any number of edit/mismatch
// However, calling these functions requires user to define themselves a 'searching model'.
// The searching model requires each BWT step to be defined. The search algorithm will then follow the defined model.

// BWTMismatchModelAnyDirection_CE matches steps with check and extend mechanism.
// It allows starting off CE in the middle of a step and recursive CE until SRACase completes.
unsigned long long BWTMismatchModelAnyDirection_CE(SRAArguments * aArgs, int i, uint8_t mismatchInserted, 
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occMismatch, uint8_t nbMismatch, int occQuality) {

    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRAResultCount * aStats  = aArgs->AlgnmtStats;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAIndex * aIndex = aArgs->AlgnmtIndex;
    SRAWorkingMemory * workMem = aArgs->AlgnmtMemory;
    
    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];

    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = (char*)qInfo->ReadQuality;
    uint8_t OutputType = aSetting->OutputType;
    HSP * hsp = aIndex->hsp;
    BWT * ceBwt = aIndex->bwt;

    int * leftMostAligned = alignmentCase->leftMostAligned;
    
    uint8_t * occMismatches = workMem->Shared_AccMismatch;
    unsigned long long * saPositions = workMem->Shared_ReadPositions;
    int * occQualities = workMem->Shared_AccQualities;
    SRAError * occErrorsPos = workMem->Shared_AccErrorsPos;
    
    unsigned long long extensionSeq[MAX_CE_BUFFER_SIZE];
    unsigned long long extensionSeq2[MAX_CE_BUFFER_SIZE];
    SRAError errorsTmp[SRA_MAX_READ_LENGTH];
    SRAError errors[MAX_CE_THRESHOLD][MAX_NUM_OF_ERROR];
    int occNBMismatch[MAX_CE_THRESHOLD];
    unsigned int readPos;
    uint8_t mismatch;
    uint8_t mismatchCe;
    uint8_t minMismatch;
    uint8_t maxMismatch;
    int8_t step;
    int len,start,end,pos,newPosCount,swap,left,right,posCount;
    int bufferCount; 

    int k=0;
    unsigned long long saCount = 0;
    step = alignmentCase->steps[stepInCase].step;
    start = alignmentCase->steps[stepInCase].start;
    pos = start + step * i;
    int leftAlign = leftMostAligned[pos];
    unsigned long long j;
    int param_i = i;
    int param_stepInCase = stepInCase;
    int flag;
    int errorsIdx;
    
    #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[BWTMismatchModelAnyDirection_CE] START CASE=%2d STEP=%2d IDX=%2d\n",alignmentCase->id,stepInCase,i);
        printf("    - SARNGR = (%12llu-%12llu)  ERR#=(%2d+%2d)\n",l,r,occMismatch,nbMismatch);
    #endif

    for (j=l;j<=r;j++) {
        occMismatches[k]=occMismatch+nbMismatch;
        occNBMismatch[k]=nbMismatch;
        occQualities[k]=occQuality;
        saPositions[k]=BWTSaValue(ceBwt,j)-leftAlign;
        for (errorsIdx=0;errorsIdx<occMismatch+mismatchInserted+nbMismatch;errorsIdx++) {
            errors[k][errorsIdx].type = occErrorsPos[errorsIdx].type;
            errors[k][errorsIdx].position = occErrorsPos[errorsIdx].position;
        }
        k++;
    }
    
    posCount = k;
    newPosCount = 0;
    
    while (alignmentCase->steps[stepInCase].type != SRA_STEP_TYPE_COMPLETE) {
        minMismatch = alignmentCase->steps[stepInCase].MinError;
        maxMismatch = alignmentCase->steps[stepInCase].MaxError;
        step = alignmentCase->steps[stepInCase].step;
        len = alignmentCase->steps[stepInCase].len;
        start = alignmentCase->steps[stepInCase].start;
        end = alignmentCase->steps[stepInCase].end;

        pos = start + step * i;
        
        swap = (pos - end) & -(pos < end);
        left = end + swap;
        right = pos - swap;
        len = right - left + 1;
        
        //bufferCount = GetPackedPatternLong(convertedKey,left,len,extensionSeq);
        bufferCount = CEPackPattern(convertedKey,left,len,extensionSeq);

        //printf("> ");
        //for (j=0;j<len;j++) {
        //    printf("%c",dnaChar[convertedKey[left+j]]);
        //    if ((j+1) % 4==0) {printf(" ");}
        //}printf("\n");

        for (j=0; j<posCount; j++) {
            readPos = saPositions[j];
            //mismatch = PackedDifferenceLong(extensionSeq,hsp,readPos+left,len);
            //mismatch = CEPackedMismatchMatching(extensionSeq,len,0,hsp,readPos+left,len);
            if ( readPos+left > ceBwt->textLength ) {
                mismatchCe = MAX_NUM_OF_ERROR + 1;
            } else {
                mismatchCe = CEPackedMismatchMatchingWithQuality(extensionSeq,len,0,hsp,readPos+left,len,errorsTmp);
                        
                /*printf("BWTMismatchModelAnyDirection_CE finds %d mismatch\n",mismatchCe);
                for (errorsIdx=0;errorsIdx<mismatchCe;errorsIdx++) {
                printf("%d/%d/%d ",errorsIdx,errorsTmp[errorsIdx].type,errorsTmp[errorsIdx].position);
                }
                printf("\n");*/
            }
            //int mismatch2 = CEPackedMismatchMatching(extensionSeq,len,0,hsp,readPos+left,len);
            //if (mismatch != mismatch2) {
            //    printf("Return number of mismatch is different : %d %d\n",mismatch,mismatch2);
            //}
            mismatch = mismatchCe + mismatchInserted;
            if (mismatch>=minMismatch && mismatch<=(maxMismatch+aSetting->MaxNBMismatch-occNBMismatch[j])) {
                occNBMismatch[newPosCount]=occNBMismatch[j]+((mismatch-maxMismatch)*(mismatch>maxMismatch));
                occQualities[newPosCount]=occQualities[j];
                for (errorsIdx=0;errorsIdx<occMismatches[j]+mismatchInserted;errorsIdx++) {
                    errors[newPosCount][errorsIdx].type = errors[j][errorsIdx].type;
                    errors[newPosCount][errorsIdx].position = errors[j][errorsIdx].position;
                    occQualities[newPosCount]=occQualities[j];
                }
                for (errorsIdx=0;errorsIdx<mismatchCe;errorsIdx++) {
                    errors[newPosCount][occMismatches[j]+mismatchInserted+errorsIdx].type = errorsTmp[errorsIdx].type;
                    errors[newPosCount][occMismatches[j]+mismatchInserted+errorsIdx].position = errorsTmp[errorsIdx].position;
                    occQualities[newPosCount]+=keyQuality[errorsTmp[errorsIdx].position];
                }
                occMismatches[newPosCount]=occMismatches[j]+mismatch;
                
                saPositions[newPosCount]=readPos;
                newPosCount++;
            }
        }
        posCount = newPosCount;
        newPosCount = 0;
        stepInCase++;
        i=0;mismatchInserted=0;
    }
    
    /*printf("BWTMismatchModelAnyDirection_CE reports %d\n",posCount);
    for (j=0;j<posCount;j++) {
        printf("%d: ",j);
        for (errorsIdx=0;errorsIdx<occMismatches[j];errorsIdx++) {
        printf("%d/%d/%d ",errorsIdx,errors[j][errorsIdx].type,errors[j][errorsIdx].position);
        }
    }
    printf("\n");*/

    return OCCReportTextPositions(aArgs, posCount, errors);
}
/*unsigned long long BWTMismatchModelAnyDirection_CE(SRAArguments * aArgs, int i, uint8_t mismatchInserted, 
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occMismatch, int occQuality) {

    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];

    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    unsigned char * convertedKey = qInfo->QueryInfo->ReadCode;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAWorkingMemory * workMem = aArgs->AlgnmtMemory;
    SRAResultCount * aStats  = aArgs->AlgnmtStats;
    uint8_t OutputType = aSetting->OutputType;
    HSP * hsp = aIndex->hsp;
    BWT * ceBwt = aIndex->bwt;

    int * leftMostAligned = alignmentCase->leftMostAligned;
    uint8_t * occMismatches = workMem->Shared_AccMismatch;
    unsigned long long * saPositions = workMem->Shared_ReadPositions;
    int * occQualities = workMem->Shared_AccQualities;
    unsigned long long extensionSeq[4];
    unsigned int readPos;
    uint8_t mismatch;
    uint8_t minMismatch;
    uint8_t maxMismatch;
    int8_t step;
    int len,start,end,pos,newPosCount,swap,left,right,posCount;
    int bufferCount;

    int k=0;
    unsigned long long saCount = 0;
    step = alignmentCase->steps[stepInCase].step;
    start = alignmentCase->steps[stepInCase].start;
    pos = start + step * i;
    int leftAlign = leftMostAligned[pos];
    unsigned long long j;
    int param_i = i;
    int param_stepInCase = stepInCase;
    int flag;

    for (j=l;j<=r;j++) {
        occMismatches[0]=occMismatch;
        occQualities[0]=occQuality;
        saPositions[0]=BWTSaValue(ceBwt,j)-leftAlign;

        i = param_i;
        stepInCase = param_stepInCase;
        flag = 1;
        mismatch = mismatchInserted;

        while (alignmentCase->steps[stepInCase].type != SRA_STEP_TYPE_COMPLETE) {

            minMismatch = alignmentCase->steps[stepInCase].MinError;
            maxMismatch = alignmentCase->steps[stepInCase].MaxError;
            step = alignmentCase->steps[stepInCase].step;
            len = alignmentCase->steps[stepInCase].len;
            start = alignmentCase->steps[stepInCase].start;
            end = alignmentCase->steps[stepInCase].end;

            pos = start + step * i;
            
            swap = (pos - end) & -(pos < end);
            left = end + swap;
            right = pos - swap;
            len = right - left + 1;
            
            bufferCount = GetPackedPatternLong(convertedKey,left,len,extensionSeq);
            mismatch += PackedDifferenceLong(extensionSeq,hsp,saPositions[0]+left,len);
            if (mismatch>=minMismatch && mismatch<=maxMismatch) {
                flag=1;
                occMismatches[0]+=mismatch;
            } else {
                flag=0;
                break;
            }

            stepInCase++;mismatch=0;
            i=0;
        }

        if (flag==1) {
            saCount+=OCCReportTextPositions(aArgs, 1);

            if (workMem->IsClosed==1) {
                return saCount;
            }
        }
    }
    return saCount;
}//*/



// BWTExactModelForward_Lookup lookup your pattern in lookup table, bi-directional and assuming forward
unsigned long long BWTExactModelForward_Lookup(SRAArguments * aArgs,
                                    SRACase * alignmentCase, int stepInCase,
                                    unsigned long long * saRanges) {
    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRAResultCount * aStats  = aArgs->AlgnmtStats;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAIndex * aIndex = aArgs->AlgnmtIndex;
    SRAWorkingMemory * workMem = aArgs->AlgnmtMemory;
    
    unsigned long long l_packedPattern = 0;
    unsigned long long r_packedPattern = 0;

    unsigned long long l = 0;
    unsigned long long r = 0;
    unsigned long long rev_l = 0;
    unsigned long long rev_r = 0;
    unsigned long long j;

    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int ALIGN_16 oL[ALPHABET_SIZE];
    unsigned int ALIGN_16 oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = (char*)qInfo->ReadQuality;
    LT * lookupTable = aIndex->lookupTable;
    LT * rev_lookupTable = aIndex->rev_lookupTable;

    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;
    int8_t step = alignmentCase->steps[stepInCase].step;
    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;
    int pos = start;
    int k;
    int branchCount = 0;
    char branchChar = 0;
    uint8_t nbMismatch = 0;

    #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[BWTExactModelForward_Lookup] START CASE=%2d STEP=%2d\n",alignmentCase->id,stepInCase);
        printf("    - SARNGR = (%12llu-%12llu)  RANGE=%3d-%3d\n",l,r,start,end);
    #endif

    unsigned int lookupLength = (lookupTable->tableSize>len)?len:lookupTable->tableSize;

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
        if ( r-l < ceThreshold && pos >= ceStart && pos <= ceEnd) {
            saRanges[0]=l;
            saRanges[1]=r;
            return BWTMismatchModelAnyDirection_CE(aArgs,i,0,
                                    alignmentCase,stepInCase,
                                    saRanges, 0, nbMismatch, 0);
        }

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
        if ( ( errType & SRA_STEP_ERROR_TYPE_NBM) > 0 && 
             oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < aSetting->MaxNBMismatch) {
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
            //printf("[BWTExactModelForward_Lookup] BWTModelSwitchBackward\n");
            return BWTModelSwitchBackward(aArgs,0,0,
                                            alignmentCase,stepInCase+1,
                                            saRanges, 0, nbMismatch, 0);
        } else {
            //printf("[BWTExactModelForward_Lookup] BWTModelSwitchAnyDirection\n");
            return BWTModelSwitchAnyDirection(aArgs,0,0,
                                                alignmentCase,stepInCase+1,
                                                saRanges, 0, nbMismatch, 0);
        }

    }

    return 0;
}


// BWTExactModelBackward_Lookup lookup your pattern in lookup table, single direction and assuming backward
unsigned long long BWTExactModelBackward_Lookup(SRAArguments * aArgs,
                                    SRACase * alignmentCase, int stepInCase,
                                    unsigned long long * saRanges) {
                                    
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAResultCount * aStats  = aArgs->AlgnmtStats;
    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRAIndex * aIndex = aArgs->AlgnmtIndex;
    SRAWorkingMemory * workMem = aArgs->AlgnmtMemory;
    
    unsigned long long l_packedPattern = 0;
    unsigned long long r_packedPattern = 0;

    unsigned long long l = 0;
    unsigned long long r = 0;
    unsigned long long j = 0;
    
    int k = 0;
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = (char*)qInfo->ReadQuality;
    LT * lookupTable = aIndex->lookupTable;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;
    int8_t step = alignmentCase->steps[stepInCase].step;
    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;

    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int ALIGN_16 oL[ALPHABET_SIZE];
    unsigned int ALIGN_16 oR[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---

    unsigned int lookupLength = (lookupTable->tableSize>len)?len:lookupTable->tableSize;
    int branchCount = 0;
    unsigned char branchChar = 0;
    uint8_t nbMismatch = 0;

    int i;
    int pos = start-lookupLength+1;

    #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[BWTExactModelBackward_Lookup] START CASE=%2d STEP=%2d\n",alignmentCase->id,stepInCase);
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
        if ( r-l < ceThreshold && pos >= ceStart && pos <= ceEnd) {
            saRanges[0]=l;
            saRanges[1]=r;
            return BWTMismatchModelAnyDirection_CE(aArgs,i,0,
                                    alignmentCase,stepInCase,
                                    saRanges, 0, nbMismatch, 0);
        }

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
        if ( ( errType & SRA_STEP_ERROR_TYPE_NBM) > 0 && 
             oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < aSetting->MaxNBMismatch) {
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
        return BWTModelSwitchBackward(aArgs,0,0,
                                alignmentCase,stepInCase+1,
                                saRanges, 0, nbMismatch, 0);
    }

    return 0;
}


// BWTExactModelBackward_Lookup lookup your pattern in lookup table, bi-directional and assuming forward
unsigned long long BWTExactModelBackwardAnyDirection_Lookup(SRAArguments * aArgs,
                                    SRACase * alignmentCase, int stepInCase,
                                    unsigned long long * saRanges) {
    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRAResultCount * aStats  = aArgs->AlgnmtStats;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAIndex * aIndex = aArgs->AlgnmtIndex;
    SRAWorkingMemory * workMem = aArgs->AlgnmtMemory;
    
    unsigned long long l_packedPattern = 0;
    unsigned long long r_packedPattern = 0;

    unsigned long long l = 0;
    unsigned long long r = 0;
    unsigned long long rev_l = 0;
    unsigned long long rev_r = 0;
    unsigned long long j;

    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int ALIGN_16 oL[ALPHABET_SIZE];
    unsigned int ALIGN_16 oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = (char*)qInfo->ReadQuality;
    LT * lookupTable = aIndex->lookupTable;
    LT * rev_lookupTable = aIndex->rev_lookupTable;

    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;
    int8_t step = alignmentCase->steps[stepInCase].step;
    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;
    
    int k;
    int branchCount = 0;
    char branchChar = 0;
    uint8_t nbMismatch = 0;

    unsigned int lookupLength = (lookupTable->tableSize>len)?len:lookupTable->tableSize;

    int i;
    int pos = start-lookupLength+1;
    
    #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[BWTExactModelBackwardAnyDirection_Lookup] START CASE=%2d STEP=%2d\n",alignmentCase->id,stepInCase);
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
        if ( r-l < ceThreshold && pos >= ceStart && pos <= ceEnd) {
            saRanges[0]=l;
            saRanges[1]=r;
            return BWTMismatchModelAnyDirection_CE(aArgs,i,0,
                                    alignmentCase,stepInCase,
                                    saRanges, 0, nbMismatch, 0);
        }

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
        if ( ( errType & SRA_STEP_ERROR_TYPE_NBM) > 0 && 
             oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < aSetting->MaxNBMismatch) {
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
            //printf("[BWTExactModelForward_Lookup] BWTModelSwitchBackward\n");
            return BWTModelSwitchBackward(aArgs,0,0,
                                            alignmentCase,stepInCase+1,
                                            saRanges, 0, nbMismatch, 0);
        } else {
            //printf("[BWTExactModelForward_Lookup] BWTModelSwitchAnyDirection\n");
            return BWTModelSwitchAnyDirection(aArgs,0,0,
                                                alignmentCase,stepInCase+1,
                                                saRanges, 0, nbMismatch, 0);
        }

    }

    return 0;
}

// BWTExactModelBackward matches pattern on text without using any other aux, e.g. lookup table.
unsigned long long BWTExactModelBackward(SRAArguments * aArgs, int i, uint8_t errorInserted,
                                    SRACase * alignmentCase, int stepInCase,
                                    unsigned long long * saRanges,
                                    uint8_t occError, uint8_t nbMismatch, int occQuality) {

    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRAResultCount * aStats  = aArgs->AlgnmtStats;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAWorkingMemory * workMem = aArgs->AlgnmtMemory;
    
    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];

    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = (char*)qInfo->ReadQuality;

    uint8_t minError = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxError = alignmentCase->steps[stepInCase].MaxError;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;
    int8_t step = -1;
    int pos, k;
    unsigned long long j;
    
    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int ALIGN_16 oL[ALPHABET_SIZE];
    unsigned int ALIGN_16 oR[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    
    int branchCount = 0;
    char branchChar = 0;
    
    pos = start + step*i;
    step = -1;

    #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[BWTExactModelBackward] START CASE=%2d STEP=%2d IDX=%2d\n",alignmentCase->id,stepInCase,i);
        printf("    - SARNGR = (%12llu-%12llu)  RANGE=%3d-%3d   ERR#=(%2d+%2d)\n",l,r,start,end,occError,nbMismatch);
    #endif
    
    while (i<len && l<=r) {
        unsigned char c = convertedKey[pos];

        if ( r-l < ceThreshold  && pos >=ceStart && pos<=ceEnd ) {
            saRanges[0] = l;
            saRanges[1] = r;
            return BWTMismatchModelAnyDirection_CE(aArgs,i,errorInserted,
                                    alignmentCase,stepInCase,
                                    saRanges, occError, nbMismatch, occQuality);
        }
        
        
        BWTAllOccValue(bwt,l,oL);
        BWTAllOccValue(bwt,r + 1,oR);
        branchCount = 0;
        branchChar = 0;
        for (k=ALPHABET_SIZE-1;k>=0;k--) {
            branchCount += (oR[k]>oL[k]);
            branchChar |= k * (oR[k]>oL[k]);
        }

        
        // Non-branching Mismatch Handler
        if ( ( errType & SRA_STEP_ERROR_TYPE_NBM) > 0 && 
             oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < aSetting->MaxNBMismatch) {
            workMem->Shared_AccErrorsPos[occError+errorInserted+nbMismatch].type = SRA_CHAR_ERROR_TYPE_NBMISMATCH;
            workMem->Shared_AccErrorsPos[occError+errorInserted+nbMismatch].position = pos;
            l = bwt->cumulativeFreq[branchChar] + oL[branchChar] + 1;
            r = bwt->cumulativeFreq[branchChar] + oR[branchChar];
            nbMismatch++;
        } else {
            l = bwt->cumulativeFreq[c] + oL[c] + 1;
            r = bwt->cumulativeFreq[c] + oR[c];
        }
        
        #ifdef SRA_2BWT_OPERATIONS_STEP_DEBUG
            int k;
            for (k=0;k<ALPHABET_SIZE;k++) {
                printf("[BWTExactModelBackward] STEP bwt->cumulativeFreq[%d] = %d\n",k,bwt->cumulativeFreq[k]);
            }
            printf("[BWTExactModelBackward] STEP %llu-%llu\n",l,r);
        #endif
        
        i++;
        pos += step;
    }
    
    if (errorInserted>=minError && l<=r) {
        saRanges[0]=l;
        saRanges[1]=r;
        return BWTModelSwitchBackward(aArgs,0,0,
                                alignmentCase,stepInCase+1,
                                saRanges, occError+errorInserted, nbMismatch, occQuality);
    }
    return 0;
}

// BWTMismatchModelAnyDirection matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
unsigned long long BWTExactModelAnyDirection(SRAArguments * aArgs, int i, uint8_t errorInserted,
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occError, uint8_t nbMismatch, int occQuality) {

    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRAResultCount * aStats  = aArgs->AlgnmtStats;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAWorkingMemory * workMem = aArgs->AlgnmtMemory;
    
    int8_t step = alignmentCase->steps[stepInCase].step;
    unsigned long long l = saRanges[(step<0)*2];
    unsigned long long r = saRanges[(step<0)*2+1];
    unsigned long long rev_l = saRanges[(step>0)*2];
    unsigned long long rev_r = saRanges[(step>0)*2+1];

    unsigned long long saCount=0;
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = (char*)qInfo->ReadQuality;
    int k, pos;
    uint8_t minError = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxError = alignmentCase->steps[stepInCase].MaxError;
    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;

    unsigned long long err_rev_l;
    unsigned long long err_rev_r;
    unsigned long long err_r;
    unsigned long long err_l;
    unsigned long long errSaRange[4];
    unsigned char c;
    unsigned long long j;
    int branchCount = 0;
    char branchChar = 0;
    
    //printf("Alignment from %d to %d allowing %u-%u mismatches.\n", start,end, minError, maxError);
    
    //printf("[BWTExactModelAnyDirection] from %d to %d allowing %u<%u<%u mismatches.\n", start,end, minError,errorInserted, maxError);
    //printf("[BWTExactModelAnyDirection] %d/%d\n", i,len);
    
    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int ALIGN_16 oL[ALPHABET_SIZE];
    unsigned int ALIGN_16 oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    
    #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[BWTExactModelAnyDirection] START CASE=%2d STEP=%2d IDX=%2d\n",alignmentCase->id,stepInCase,i);
        printf("    - SARNGR = (%12llu-%12llu)  RANGE=%3d-%3d   ERR#=(%2d+%2d)\n",l,r,start,end,occError,nbMismatch);
    #endif
    
    pos = start + step * i;
    len -= (minError-errorInserted-1) * (minError>errorInserted+1);
    while (i < len && l<=r) {
        if ( r-l < ceThreshold && pos >=ceStart && pos<=ceEnd) {
            saRanges[(step<0)*2] = l;
            saRanges[(step<0)*2+1] = r;
            saRanges[(step>0)*2] = rev_l;
            saRanges[(step>0)*2+1] = rev_r;
            k=0;
            return BWTMismatchModelAnyDirection_CE(aArgs,i,errorInserted,
                                    alignmentCase,stepInCase,
                                    saRanges, occError, nbMismatch, occQuality);
        }


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
        
        /*int testFlag = 0;
        if ( l > r ) {
            for (k=0;k<ALPHABET_SIZE;k++) {
                if ( bwt->cumulativeFreq[c] + oL[c] + 1 <= bwt->cumulativeFreq[c] + oR[c]) {
                    testFlag++;
                }
            }
            if (testFlag==1 && nbMismatch < aSetting->MaxNBMismatch) {
                testFlag = 1;
            } else {
                testFlag = 0;
            }
        }*/

        // Non-branching Mismatch Handler
        if ( ( errType & SRA_STEP_ERROR_TYPE_NBM) > 0 && 
             oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < aSetting->MaxNBMismatch) {
            workMem->Shared_AccErrorsPos[occError+errorInserted+nbMismatch].type = SRA_CHAR_ERROR_TYPE_NBMISMATCH;
            workMem->Shared_AccErrorsPos[occError+errorInserted+nbMismatch].position = pos;
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
        /* else {
            if ( testFlag == 1 ) {
                printf("NBM missed %u-%u\n",l,r);
            }
        }*/
        
        i += 1;
        pos += step;
    }

    if (errorInserted>=minError && l<=r) {
        //Next Step
        saRanges[(step<0)*2] = l;
        saRanges[(step<0)*2+1] = r;
        saRanges[(step>0)*2] = rev_l;
        saRanges[(step>0)*2+1] = rev_r;
        return BWTModelSwitchAnyDirection(aArgs,0,0,
                                alignmentCase,stepInCase+1,
                                saRanges, occError+errorInserted, nbMismatch, occQuality);
    }

    return 0;
}

// BWTMismatchModelAnyDirection matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
unsigned long long BWTMismatchModelAnyDirection(SRAArguments * aArgs, int i, uint8_t mismatchInserted,
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occMismatch, uint8_t nbMismatch, int occQuality) {

    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRAResultCount * aStats  = aArgs->AlgnmtStats;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAWorkingMemory * workMem = aArgs->AlgnmtMemory;
    
    int8_t step = alignmentCase->steps[stepInCase].step;
    unsigned long long l = saRanges[(step<0)*2];
    unsigned long long r = saRanges[(step<0)*2+1];
    unsigned long long rev_l = saRanges[(step>0)*2];
    unsigned long long rev_r = saRanges[(step>0)*2+1];

    unsigned long long saCount=0;
    
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = (char*)qInfo->ReadQuality;

    int k, pos;

    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    uint8_t minMismatch = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxMismatch = alignmentCase->steps[stepInCase].MaxError;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;

    unsigned long long err_rev_l;
    unsigned long long err_rev_r;
    unsigned long long err_r;
    unsigned long long err_l;
    unsigned long long errSaRange[4];
    unsigned long long j;
    unsigned char c;
    int branchCount = 0;
    char branchChar = 0;
    
    //printf("Alignment from %d to %d allowing %u-%u mismatches.\n", start,end, minMismatch, maxMismatch);
    
    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int ALIGN_16 oL[ALPHABET_SIZE];
    unsigned int ALIGN_16 oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    
    #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[BWTMismatchModelAnyDirection] START CASE=%2d STEP=%2d IDX=%2d\n",alignmentCase->id,stepInCase,i);
        printf("    - SARNGR = (%12llu-%12llu)  RANGE=%3d-%3d   ERR#=(%2d+%2d)\n",l,r,start,end,occMismatch,nbMismatch);
    #endif
    
    pos = start + step * i;
    len -= (minMismatch-mismatchInserted-1) * (minMismatch>mismatchInserted+1);
    while (i < len && l<=r) {
        //printf("Current progress: %d(%d) within %d characters. stepping by %d.\n", i, pos, len, step);
        if ( r-l < ceThreshold && pos >=ceStart && pos<=ceEnd ) {
            saRanges[(step<0)*2] = l;
            saRanges[(step<0)*2+1] = r;
            saRanges[(step>0)*2] = rev_l;
            saRanges[(step>0)*2+1] = rev_r;
            saCount+=BWTMismatchModelAnyDirection_CE(aArgs,i,mismatchInserted,
                                    alignmentCase,stepInCase,
                                    saRanges, occMismatch, nbMismatch, occQuality);
            return saCount;
        }


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

        //if (mismatchInserted<maxMismatch) {
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
                    saCount+=BWTModelSwitchAnyDirection(aArgs,i+1,mismatchInserted+1,
                                            alignmentCase,stepInCase,
                                            errSaRange, occMismatch, nbMismatch, 
                                            occQuality+keyQuality[pos]);
                    if (workMem->IsClosed==1) { return saCount; }
                }
            }
        }
        //}

        // Non-branching Mismatch Handler
        if ( ( errType & SRA_STEP_ERROR_TYPE_NBM) > 0 && 
             oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < aSetting->MaxNBMismatch) {
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
        saCount+=BWTModelSwitchAnyDirection(aArgs,0,0,
                                alignmentCase,stepInCase+1,
                                saRanges, occMismatch+mismatchInserted, nbMismatch, occQuality);
    }

    return saCount;
}

// BWTMismatchModelBackward matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
unsigned long long BWTMismatchModelBackward(SRAArguments * aArgs, int i, uint8_t mismatchInserted,
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occMismatch, uint8_t nbMismatch, int occQuality) {

    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRAResultCount * aStats  = aArgs->AlgnmtStats;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAWorkingMemory * workMem = aArgs->AlgnmtMemory;
    
    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];

    unsigned long long saCount=0;
    int k, pos;
    
    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = (char*)qInfo->ReadQuality;

    uint8_t minMismatch = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxMismatch = alignmentCase->steps[stepInCase].MaxError;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;
    int8_t step = -1; //alignmentCase->steps[stepInCase].step;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;

    unsigned long long errSaRange[4];
    unsigned long long err_l;
    unsigned long long err_r;
    unsigned char c;

    //printf("Alignment from %d to %d allowing %u-%u mismatches. i = %d\n", start,end, minMismatch, maxMismatch,i);
    
    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int ALIGN_16 oL[ALPHABET_SIZE];
    unsigned int ALIGN_16 oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    int branchCount = 0;
    char branchChar = 0;
    
    #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[BWTMismatchModelBackward] START CASE=%2d STEP=%2d IDX=%2d\n",alignmentCase->id,stepInCase,i);
        printf("    - SARNGR = (%12llu-%12llu)  RANGE=%3d-%3d   ERR#=(%2d+%2d)\n",l,r,start,end,occMismatch,nbMismatch);
    #endif
    
    pos = start + step * i;
    len -= (minMismatch-mismatchInserted-1) * (minMismatch>(mismatchInserted+1));

    while (i < len && l<=r) {
        if ( r-l < ceThreshold && pos>=ceStart && pos<=ceEnd) {
            saRanges[0] = l;
            saRanges[1] = r;
            saCount += BWTMismatchModelAnyDirection_CE(aArgs,i,mismatchInserted,
                                    alignmentCase,stepInCase,
                                    saRanges, occMismatch, nbMismatch, occQuality);
            return saCount;
        }
        c = convertedKey[pos];
        BWTAllOccValue(bwt,l,oL);
        BWTAllOccValue(bwt,r + 1,oR);

        //if (mismatchInserted<maxMismatch) {
        branchCount = 0;
        branchChar = 0;
        for (k=0;k<ALPHABET_SIZE;k++) {
            branchCount += (oR[k]>oL[k]);
            branchChar |= k * (oR[k]>oL[k]);
        }

        for (k=0;k<ALPHABET_SIZE;k++) {
            if (k!=c) {
                
                errSaRange[0] = bwt->cumulativeFreq[k] + oL[k] + 1;
                errSaRange[1] = bwt->cumulativeFreq[k] + oR[k];

                if (errSaRange[0] <= errSaRange[1]) {
                    workMem->Shared_AccErrorsPos[occMismatch+mismatchInserted].type = SRA_CHAR_ERROR_TYPE_MISMATCH;
                    workMem->Shared_AccErrorsPos[occMismatch+mismatchInserted].position = pos;
                    saCount+=BWTModelSwitchBackward(aArgs,i+1,mismatchInserted+1,
                                        alignmentCase,stepInCase,
                                        errSaRange, occMismatch, nbMismatch, occQuality+keyQuality[pos]);
                    if (workMem->IsClosed==1) { return saCount; }
                }
            }
        }
        //}
        
        // Non-branching Mismatch Handler
        if ( ( errType & SRA_STEP_ERROR_TYPE_NBM) > 0 && 
             oL[c] + 1 > oR[c] && branchCount == 1 && nbMismatch < aSetting->MaxNBMismatch) {
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
        saCount+=BWTModelSwitchBackward(aArgs,0,0,
                                alignmentCase,stepInCase+1,
                                saRanges, occMismatch+mismatchInserted, nbMismatch, occQuality);
    }

    return saCount;
}


// BWTEditModelAnyDirection matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
unsigned long long BWTEditModelAnyDirection(SRAArguments * aArgs, int i, uint8_t editInserted,
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occEdit, uint8_t nbMismatch, int occQuality) {

    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRAResultCount * aStats  = aArgs->AlgnmtStats;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAWorkingMemory * workMem = aArgs->AlgnmtMemory;
    
    int8_t step = alignmentCase->steps[stepInCase].step;
    unsigned long long l = saRanges[(step<0)*2];
    unsigned long long r = saRanges[(step<0)*2+1];
    unsigned long long rev_l = saRanges[(step>0)*2];
    unsigned long long rev_r = saRanges[(step>0)*2+1];

    unsigned long long saCount=0;
    
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = (char*)qInfo->ReadQuality;

    int pos, misLen, delLen;
    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;
    uint8_t minEdit = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxEdit = alignmentCase->steps[stepInCase].MaxError;
    int residueEdit = maxEdit - editInserted;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;

    unsigned long long err_rev_l;
    unsigned long long err_rev_r;
    unsigned long long err_r;
    unsigned long long err_l;
    unsigned long long errSaRange[4];
    unsigned char c,nc;
    int j;
    int k;
    
    //printf("[BWTEditModelAnyDirection] Read %llu Alignment from %d to %d, i = %d\n", aArgs->ReadId,start,end,i);
    //printf("[BWTEditModelAnyDirection] ErrorType=%d  allowing %u-%u edits.\n", errType, minEdit, maxEdit);
    
    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int ALIGN_16 oL[ALPHABET_SIZE];
    unsigned int ALIGN_16 oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    
    nc = ALPHABET_SIZE + 1;
    pos = start + step * i;
    misLen = len - (minEdit-editInserted-1) * (minEdit>(editInserted+1));
    delLen = len - (maxEdit-editInserted-1) * (maxEdit>(editInserted+1));

    #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[BWTEditModelAnyDirection] START CASE=%2d STEP=%2d IDX=%2d\n",alignmentCase->id,stepInCase,i);
        printf("    - SARNGR = (%12llu-%12llu)  RANGE=%3d-%3d   ERR#=(%2d+%2d)\n",l,r,start,end,occEdit,nbMismatch);
    #endif
    
    while (i < len && l<=r) {
        //printf("Current progress: %d(%d) within %d characters. stepping by %d.\n", i, pos, len, step);

        /*if ( r-l < ceThreshold && pos >=ceStart && pos<=ceEnd ) {
            k=0;
            saRanges[(step<0)*2] = l;
            saRanges[(step<0)*2+1] = r;
            saRanges[(step>0)*2] = rev_l;
            saRanges[(step>0)*2+1] = rev_r;
            for (j=saRanges[0];j<=saRanges[1];j++) {
                occEdits[k]=occEdit;
                occQualities[k]=occQuality;
                saPositions[k++]=BWTSaValue(ceBwt,j)-leftMostAligned[pos];

            }
            saCount+=BWTMismatchModelAnyDirection_CE(aArgs,aStats,i,editInserted,
                                    hsp,
                                    alignmentCase,stepInCase,k);
            return saCount;
        }*/


        c = convertedKey[pos];
        if (pos+step>=0 && pos+step<=qInfo->ReadLength-1) {nc=convertedKey[pos+step];} else {nc = ALPHABET_SIZE + 1;}
        BWTAllOccValue(bwt,rev_l,oL);
        BWTAllOccValue(bwt,rev_r + 1,oR);
        oCount[ALPHABET_SIZE-1]=0;
        for (k=ALPHABET_SIZE-2;k>=0;k--) {
            oCount[k]=oCount[k+1]+oR[k+1]-oL[k+1];
        }

        if (editInserted<maxEdit) {

            //Delete
            if (c!=nc) {
                //printf("[BWTEditModelAnyDirection] Delete at position %u (i = %d)\n", pos, i);
                if (i>=0 && i < len - 1) {
                    err_rev_l = bwt->cumulativeFreq[nc] + oL[nc] + 1;
                    err_rev_r = bwt->cumulativeFreq[nc] + oR[nc];
                    err_r = r - oCount[nc];
                    err_l = err_r - (err_rev_r-err_rev_l);
                } else {
                    err_rev_l = rev_l;
                    err_rev_r = rev_r;
                    err_r = r;
                    err_l = l;
                }


                errSaRange[(step<0)*2] = err_l;
                errSaRange[(step<0)*2+1] = err_r;
                errSaRange[(step>0)*2] = err_rev_l;
                errSaRange[(step>0)*2+1] = err_rev_r;
                workMem->Shared_AccErrorsPos[occEdit+editInserted].type = SRA_CHAR_ERROR_TYPE_DELETE;
                workMem->Shared_AccErrorsPos[occEdit+editInserted].position = pos;
                unsigned long long result=BWTModelSwitchAnyDirection(aArgs,i+2,editInserted+1,
                                        alignmentCase,stepInCase,
                                        errSaRange, occEdit, nbMismatch, occQuality);
                //if (result) 
                //    printf("[BWTEditModelAnyDirection:Case%d] Delete at position %u (i = %d)\n", alignmentCase->id, pos, i);
                saCount+=result;
                if (workMem->IsClosed==1) { return saCount; }
            }
            //Long delete
            for (k=residueEdit;k>1;k--) {
                if (len - i >= k + 1 && c!=convertedKey[pos+k*step]) {
                    //printf("[BWTEditModelAnyDirection] Long Delete at position %u (i = %d, gap len = %d)\n", pos, i, k);
                    errSaRange[(step<0)*2] = l;
                    errSaRange[(step<0)*2+1] = r;
                    errSaRange[(step>0)*2] = rev_l;
                    errSaRange[(step>0)*2+1] = rev_r;
                    for (j=0;j<k;j++) {
                        workMem->Shared_AccErrorsPos[occEdit+editInserted+j].type = SRA_CHAR_ERROR_TYPE_DELETE;
                        workMem->Shared_AccErrorsPos[occEdit+editInserted+j].position = pos+(step)*j;
                    }
                    unsigned long long result=BWTModelSwitchAnyDirection(aArgs,i+k,editInserted+k,
                                            alignmentCase,stepInCase,
                                            errSaRange, occEdit, nbMismatch, occQuality);
                    //if (result) 
                    //    printf("[BWTEditModelAnyDirection:Case%d] Long Delete at position %u (i = %d, gap len = %d)\n", alignmentCase->id, pos, i, k);
                    saCount+=result;
                    if (workMem->IsClosed==1) { return saCount; }
                }
            }

            for (k=0;k<ALPHABET_SIZE;k++) {
                err_rev_l = bwt->cumulativeFreq[k] + oL[k] + 1;
                err_rev_r = bwt->cumulativeFreq[k] + oR[k];
                err_r = r - oCount[k];
                err_l = err_r - (err_rev_r-err_rev_l);

                if (err_l<=err_r) {
                    if (k!=c && i < misLen) {
                        //Mismatch

                        //printf("[BWTEditModelAnyDirection] Mismatch at position %u (i = %d)\n", pos, i);
                        errSaRange[(step<0)*2] = err_l;
                        errSaRange[(step<0)*2+1] = err_r;
                        errSaRange[(step>0)*2] = err_rev_l;
                        errSaRange[(step>0)*2+1] = err_rev_r;
                        workMem->Shared_AccErrorsPos[occEdit+editInserted].type = SRA_CHAR_ERROR_TYPE_MISMATCH;
                        workMem->Shared_AccErrorsPos[occEdit+editInserted].position = pos;
                        unsigned long long result=BWTModelSwitchAnyDirection(aArgs,i+1,editInserted+1,
                                                alignmentCase,stepInCase,
                                                errSaRange, occEdit, nbMismatch, occQuality);
                        //if (result) 
                        //    printf("[BWTEditModelAnyDirection:Case%d] Mismatch at position %u (i = %d)\n", alignmentCase->id, pos, i);
                        saCount+=result;
                        if (workMem->IsClosed==1) { return saCount; }
                    }
                    
                    //Insert
                    if ((k!=c && (i!=0 || errType==SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY))) {
                        //printf("[BWTEditModelAnyDirection] Insert at position %u (i = %d)\n", pos, i);
                        errSaRange[(step<0)*2] = err_l;
                        errSaRange[(step<0)*2+1] = err_r;
                        errSaRange[(step>0)*2] = err_rev_l;
                        errSaRange[(step>0)*2+1] = err_rev_r;
                        workMem->Shared_AccErrorsPos[occEdit+editInserted].type = SRA_CHAR_ERROR_TYPE_INSERT;
                        workMem->Shared_AccErrorsPos[occEdit+editInserted].position = pos;
                        unsigned long long result=BWTModelSwitchAnyDirection(aArgs,i,editInserted+1,
                                                alignmentCase,stepInCase,
                                                errSaRange, occEdit, nbMismatch, occQuality);
                        //if (result) 
                        //    printf("[BWTEditModelAnyDirection:Case%d] Insert at position %u (i = %d)\n", alignmentCase->id, pos, i);
                        saCount+=result;
                        if (workMem->IsClosed==1) { return saCount; }
                    }
                }
            }
        }

        rev_l = bwt->cumulativeFreq[c] + oL[c] + 1;
        rev_r = bwt->cumulativeFreq[c] + oR[c];
        r = r - oCount[c];
        l = r - (rev_r-rev_l);

        i += 1;
        pos += step;
    }

    if (editInserted>=minEdit && l<=r) {
        //Next Step
        saRanges[(step<0)*2] = l;
        saRanges[(step<0)*2+1] = r;
        saRanges[(step>0)*2] = rev_l;
        saRanges[(step>0)*2+1] = rev_r;
        saCount+=BWTModelSwitchAnyDirection(aArgs,0,0,
                                alignmentCase,stepInCase+1,
                                saRanges, occEdit+editInserted, nbMismatch, occQuality);
    }

    return saCount;
}

// BWTEditModelBackward matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
unsigned long long BWTEditModelBackward(SRAArguments * aArgs, int i, uint8_t editInserted,
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occEdit, uint8_t nbMismatch, int occQuality) {

    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRAResultCount * aStats  = aArgs->AlgnmtStats;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAWorkingMemory * workMem = aArgs->AlgnmtMemory;
    
    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];

    unsigned long long saCount=0;
    int pos, delLen,misLen;

    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = (char*)qInfo->ReadQuality;

    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;
    uint8_t minEdit = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxEdit = alignmentCase->steps[stepInCase].MaxError;
    int residueEdit = maxEdit - editInserted;
    int8_t step = -1; //alignmentCase->steps[stepInCase].step;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;
    int ceThreshold = alignmentCase->steps[stepInCase].ceThreshold;
    int ceStart = alignmentCase->steps[stepInCase].ceStart;
    int ceEnd = alignmentCase->steps[stepInCase].ceEnd;

    unsigned long long errSaRange[4];
    unsigned long long err_l;
    unsigned long long err_r;
    unsigned char nc,c,k;
    int j;
 
    //printf("[BWTEditModelBackward] Read %llu Alignment from %d to %d, i = %d\n", qInfo->ReadId,start,end,i);
//printf("[BWTEditModelBackward] ErrorType=%d  allowing %u-%u edits.\n", errType, minEdit, maxEdit);
    
    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int ALIGN_16 oL[ALPHABET_SIZE];
    unsigned int ALIGN_16 oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    
    pos = start + step * i;
    misLen = len - (minEdit-editInserted-1) * (minEdit>(editInserted+1));
    delLen = len - (maxEdit-editInserted-1) * (maxEdit>(editInserted+1));

    #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[BWTEditModelBackward] START CASE=%2d STEP=%2d IDX=%2d\n",alignmentCase->id,stepInCase,i);
        printf("    - SARNGR = (%12llu-%12llu)  RANGE=%3d-%3d   ERR#=(%2d+%2d)\n",l,r,start,end,occEdit,nbMismatch);
    #endif
    
    nc = ALPHABET_SIZE + 1;
    while (i < len && l<=r) {
        /*if ( r-l < ceThreshold && pos >=ceStart && pos<=ceEnd) {
            k=0;
            for (err_l=l;err_l<=r;err_l++) {
                occEdits[k]=occEdit;
                occQualities[k]=occQuality;
                saPositions[k++]=BWTSaValue(ceBwt,err_l)-leftMostAligned[pos];
            }
            saCount += BWTMismatchModelAnyDirection_CE(aArgs,aStats,i,editInserted,
                                    hsp,
                                    alignmentCase,stepInCase,k);
            return saCount;
        }*/                
        //printf("Current progress: %d(%d) within %d characters. stepping by %d.\n", i, pos, len, step);
        c = convertedKey[pos];
        if (pos+step>=0 && pos+step<=qInfo->ReadLength-1) {nc=convertedKey[pos+step];} else {nc = ALPHABET_SIZE + 1;}
        BWTAllOccValue(bwt,l,oL);
        BWTAllOccValue(bwt,r + 1,oR);

        if (editInserted<maxEdit) {

            //Delete
            if (c!=nc) {
                if (i>=0 && i < len - 1) {
                    err_l = bwt->cumulativeFreq[nc] + oL[nc] + 1;
                    err_r = bwt->cumulativeFreq[nc] + oR[nc];
                } else {
                    err_r = r;
                    err_l = l;
                }
                errSaRange[0] = err_l;
                errSaRange[1] = err_r;
                workMem->Shared_AccErrorsPos[occEdit+editInserted].type = SRA_CHAR_ERROR_TYPE_DELETE;
                workMem->Shared_AccErrorsPos[occEdit+editInserted].position = pos;
                unsigned long long result = BWTModelSwitchBackward(aArgs,i+2,editInserted+1,
                                    alignmentCase,stepInCase,
                                    errSaRange, occEdit, nbMismatch, occQuality);
                //if (result) 
                //    printf("[BWTEditModelBackward:Case%d] Delete at position %u (i = %d)\n", alignmentCase->id, pos, i);
                saCount+=result;
                if (workMem->IsClosed==1) { return saCount; }
            }
            //Long delete
            for (k=residueEdit;k>1;k--) {
                if (len - i >= k + 1 && c!=convertedKey[pos+k*step]) {
                    errSaRange[0] = l;
                    errSaRange[1] = r;
                    for (j=0;j<k;j++) {
                        workMem->Shared_AccErrorsPos[occEdit+editInserted+j].type = SRA_CHAR_ERROR_TYPE_DELETE;
                        workMem->Shared_AccErrorsPos[occEdit+editInserted+j].position = pos+(step)*j;
                    }
                    unsigned long long result=BWTModelSwitchBackward(aArgs,i+k,editInserted+k,
                                            alignmentCase,stepInCase,
                                            errSaRange, occEdit, nbMismatch, occQuality);
                    //if (result) 
                    //    printf("[BWTEditModelBackward:Case%d] Long Delete at position %u (i = %d, gap len = %d)\n", alignmentCase->id, pos, i, k);
                    saCount+=result;
                    if (workMem->IsClosed==1) { return saCount; }
                }
            }


            for (k=0;k<ALPHABET_SIZE;k++) {
                
                err_l = bwt->cumulativeFreq[k] + oL[k] + 1;
                err_r = bwt->cumulativeFreq[k] + oR[k];
                
                if (err_l <= err_r) {

                    if (k!=c && i < misLen) {
                        //Mismatch
                        //printf("[BWTEditModelBackward] Mismatch at position %u (i = %d)\n", pos, i);
                        errSaRange[0] = err_l;
                        errSaRange[1] = err_r;
                        workMem->Shared_AccErrorsPos[occEdit+editInserted].type = SRA_CHAR_ERROR_TYPE_MISMATCH;
                        workMem->Shared_AccErrorsPos[occEdit+editInserted].position = pos;
                        unsigned long long result=BWTModelSwitchBackward(aArgs,i+1,editInserted+1,
                                            alignmentCase,stepInCase,
                                            errSaRange, occEdit, nbMismatch, occQuality);
                        //if (result) 
                        //    printf("[BWTEditModelBackward:Case%d] Mismatch at position %u (i = %d)\n", alignmentCase->id, pos, i);
                        saCount+=result;
                        if (workMem->IsClosed==1) { return saCount; }
                    }

                    //Insert
                    if ((k!=c && (errType== SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY || i != 0))) {
                        //printf("[BWTEditModelBackward] Insert at position %u (i = %d)\n", pos, i);
                        errSaRange[0] = err_l;
                        errSaRange[1] = err_r;
                        workMem->Shared_AccErrorsPos[occEdit+editInserted].type = SRA_CHAR_ERROR_TYPE_INSERT;
                        workMem->Shared_AccErrorsPos[occEdit+editInserted].position = pos;
                        unsigned long long result=BWTModelSwitchBackward(aArgs,i,editInserted+1,
                                            alignmentCase,stepInCase,
                                            errSaRange, occEdit, nbMismatch, occQuality);
                        
                        //if (result) 
                        //    printf("[BWTEditModelBackward:Case%d] Insert at position %u (i = %d)\n", alignmentCase->id, pos, i);
                        saCount+=result;
                        if (workMem->IsClosed==1) { return saCount; }
                    }
                }
            }
        }

        l = bwt->cumulativeFreq[c] + oL[c] + 1;
        r = bwt->cumulativeFreq[c] + oR[c];

        i += 1;
        pos += step;
    }

    if (editInserted>=minEdit && l<=r) {
        //Next Step
        saRanges[0] = l;
        saRanges[1] = r;
        saCount+=BWTModelSwitchBackward(aArgs,0,0,
                                alignmentCase,stepInCase+1,
                                saRanges, occEdit+editInserted, nbMismatch, occQuality);
    }

    return saCount;
}

// BWTOptionalModelAnyDirection matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
// It returns whenever the alignment reaches it ends (before the end of BWT search or end of read)
unsigned long long BWTOptionalModelAnyDirection(SRAArguments * aArgs, int i, uint8_t errorInserted,
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occError, uint8_t nbMismatch, int occQuality) {

    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRAResultCount * aStats  = aArgs->AlgnmtStats;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAWorkingMemory * workMem = aArgs->AlgnmtMemory;
    
    int8_t step = alignmentCase->steps[stepInCase].step;
    unsigned long long l = saRanges[(step<0)*2];
    unsigned long long r = saRanges[(step<0)*2+1];
    unsigned long long rev_l = saRanges[(step>0)*2];
    unsigned long long rev_r = saRanges[(step>0)*2+1];

    unsigned long long saCount=0;
    unsigned char * convertedKey = qInfo->ReadCode;
    char * keyQuality = (char*)qInfo->ReadQuality;
    int k, pos;
    uint8_t minError = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxError = alignmentCase->steps[stepInCase].MaxError;
    uint8_t errType = alignmentCase->steps[stepInCase].ErrorType;
    BWT * bwt = alignmentCase->steps[stepInCase].bwt;
    int len = alignmentCase->steps[stepInCase].len;
    int start = alignmentCase->steps[stepInCase].start;
    int end = alignmentCase->steps[stepInCase].end;

    unsigned char c;
    unsigned long long j;
    int branchCount = 0;
    char branchChar = 0;

    int trimmedLength = 0;
    
    //MARK_FOR_64_ENHANCEMENT ----
    unsigned int ALIGN_16 oL[ALPHABET_SIZE];
    unsigned int ALIGN_16 oR[ALPHABET_SIZE];
    unsigned int oCount[ALPHABET_SIZE];
    //MARK_FOR_64_ENHANCEMENT ---
    
    pos = start + step * i;
    len -= (minError-errorInserted-1) * (minError>errorInserted+1);
    
    #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[BWTOptionalModelAnyDirection] START CASE=%2d STEP=%2d IDX=%2d\n",alignmentCase->id,stepInCase,i);
        printf("    - SARNGR = (%12llu-%12llu)  RANGE=%3d-%3d   ERR#=(%2d+%2d)\n",l,r,start,end,occError,nbMismatch);
    #endif
    
    while ( i < len && l <= r ) {
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
        if ( ( errType & SRA_STEP_ERROR_TYPE_NBM) > 0 && 
             oL[c] + 1 > oR[c] ) {
            if ( branchCount == 1 && nbMismatch < aSetting->MaxNBMismatch) {
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
        /* else {
            if ( testFlag == 1 ) {
                printf("NBM missed %u-%u\n",l,r);
            }
        }*/
        
        i += 1;
        pos += step;
    }
    trimmedLength = len - i;

    #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[BWTOptionalModelAnyDirection] TERM@%2d/%2d\n",i,len);
    #endif
    
    if (errorInserted>=minError && l <= r) {

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
        return BWTModelSwitchAnyDirection(aArgs,0,0,
                                alignmentCase,stepInCase+1,
                                saRanges, occError+errorInserted, nbMismatch, occQuality);
    }

    return 0;
}

unsigned long long BWTModelSwitchAnyDirection(SRAArguments * aArgs, int i, uint8_t errorInserted,
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occError, uint8_t nbMismatch, int occQuality) {
                                    
    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];

    #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[BWTModelSwitchAnyDirection] START CASE=%2d STEP=%2d\n",alignmentCase->id,stepInCase);
        printf("    - SARNGR = (%12llu-%12llu)\n",l,r);
    #endif
    
    if (aArgs->AlgnmtMemory->IsClosed==1) {
        #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
            printf("[BWTModelSwitchBackward] EXIT NONE\n\n");
        #endif
        return 0;
    }
    
    if (alignmentCase->steps[stepInCase].type == SRA_STEP_TYPE_COMPLETE) {
        #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
            printf("[BWTModelSwitchAnyDirection] EXIT = (%12llu-%12llu)  ERR#=(%2d+%2d)\n\n",l,r,occError,nbMismatch);
        #endif
        return OCCReportSARange(aArgs,l,r,occError+errorInserted+nbMismatch,occQuality);
    }

    uint8_t errorType = alignmentCase->steps[stepInCase].ErrorType;
    uint8_t minError = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxError = alignmentCase->steps[stepInCase].MaxError;
    uint8_t OutputType = aArgs->AlgnmtSetting->OutputType;
    
    if ( ((errorType & SRA_STEP_ERROR_TYPE_NO_ERROR)) > 0) {
        return BWTExactModelAnyDirection(aArgs,i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    } else if ( (errorType & SRA_STEP_ERROR_TYPE_EDIT_DISTANCE) > 0 || (errorType & SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY) > 0) {
        return BWTEditModelAnyDirection(aArgs,i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    } else if ( (errorType & SRA_STEP_ERROR_TYPE_OPTIONAL_NBM) > 0) {
        return BWTOptionalModelAnyDirection(aArgs,i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    } else if (errorInserted>=maxError) {
        return BWTExactModelAnyDirection(aArgs,i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    } else if ((errorType & SRA_STEP_ERROR_TYPE_MISMATCH_ONLY) > 0) {
        return BWTMismatchModelAnyDirection(aArgs,i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    }

    fprintf(stderr,"[BWTModelSwitchAnyDirection] Unknown error type[%u] encountered\n",errorType);
    return 0;
}

unsigned long long BWTModelSwitchBackward(SRAArguments * aArgs,  int i, uint8_t errorInserted,
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occError, uint8_t nbMismatch, int occQuality) {
                                    
    unsigned long long l = saRanges[0];
    unsigned long long r = saRanges[1];
    unsigned long long j;
    
    #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
        printf("[BWTModelSwitchBackward] START CASE=%2d STEP=%2d\n",alignmentCase->id,stepInCase);
        printf("    - SARNGR = (%12llu-%12llu)\n",l,r);
    #endif
    
    if (aArgs->AlgnmtMemory->IsClosed==1) {
        #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
            printf("[BWTModelSwitchBackward] EXIT NONE\n\n");
        #endif
        return 0;
    }
    
    if (alignmentCase->steps[stepInCase].type == SRA_STEP_TYPE_COMPLETE) {
        #ifdef SRA_2BWT_OPERATIONS_START_END_DEBUG
            printf("[BWTModelSwitchBackward] EXIT = (%12llu-%12llu)  ERR#=(%2d+%2d)\n\n",l,r,occError,nbMismatch);
        #endif
        return OCCReportSARange(aArgs,l,r,occError+errorInserted+nbMismatch,occQuality);
    }

    uint8_t errorType = alignmentCase->steps[stepInCase].ErrorType;
    uint8_t minError = alignmentCase->steps[stepInCase].MinError;
    uint8_t maxError = alignmentCase->steps[stepInCase].MaxError;
    uint8_t OutputType = aArgs->AlgnmtSetting->OutputType;

    if ((errorType & SRA_STEP_ERROR_TYPE_NO_ERROR) > 0) {
        return BWTExactModelBackward(aArgs,  i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    } else if ((errorType & SRA_STEP_ERROR_TYPE_EDIT_DISTANCE) > 0 || (errorType & SRA_STEP_ERROR_TYPE_EDIT_DISTANCE_BOUNDARY) > 0) {
        return BWTEditModelBackward(aArgs, i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    } else if (errorInserted>=maxError) {
        return BWTExactModelBackward(aArgs,  i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    } else if ((errorType & SRA_STEP_ERROR_TYPE_MISMATCH_ONLY) > 0) {
        return BWTMismatchModelBackward(aArgs, i,errorInserted,alignmentCase,stepInCase,saRanges, occError, nbMismatch, occQuality);
    } 

    fprintf(stderr,"[BWTModelSwitchBackward] Unknown error type[%u] encountered\n",errorType);
    return 0;
}
