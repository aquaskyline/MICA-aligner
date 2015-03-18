//
//    OutputSAMFormat.c
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
   
   Date   : 18th September 2011
   Author : Edward MK Wu
   Change : New file.
   
*/
/////////////////////////////////////////////////////

#include "OutputSAMFormat.h"

// OutputSAMFormat CTC#1
// -------------
// Set this parameter SAM_READ_NAME_IGNORE_PAIR_FLAG to 1 if you want
// OutputSAMFormat to trim the tailing "/1" and "/2" on the read name
#define SAM_READ_NAME_IGNORE_PAIR_FLAG 1
// -------------

// --------------------------------------------------
// DP_DEBUG #1
// Uncomment the below to turn on logging message for
// SAM output, particular the errors after sorting.
//#define SAM_DEBUG_PRINT_CIGAR_PROGRESS
//#define SAM_DEBUG_PRINT_DP_CIGAR_PROGRESS
// --------------------------------------------------

void SAMCreateCIGAR(int * mis_pos, int num_mis, char * cigar) {
    int i = 0;
    char stack[10];
    int j;
    int k;
    int tmp;
    int l=0;
    
    for (i=0;i<num_mis;i++) {
        tmp = mis_pos[i];
        
        j=0;
        while (tmp>0) {
            stack[j++]='0'+tmp%10;
        }
        
        for (k=j-1;k>=0;k--) {
            cigar[l++]=stack[k];
        }
    }
    
    cigar[l]='\0';
}

void SAMIUint8ConcatUint8(uint8_t * data, int * curSize,
                                uint8_t key) {
    data[(*curSize)++] = key;
}
void SAMIUint8ConcatUint32(uint8_t * data, int * curSize,
                                uint32_t key) {
    int i;
    int len = sizeof(uint32_t) / sizeof(uint8_t);
    uint8_t * key_8 = (uint8_t *) &key;
    for (i=0;i<len;i++) {
        data[(*curSize)++] = key_8[i];
    }
}

__attribute__((target(mic)))
int SAMIUint8ConcatUint32AsString(uint8_t * data, int * curSize,
                                uint32_t key) {
    
    if (key==0) {
        SAMIUint8ConcatString(data,curSize,"0",1);
        return 1;
    }
    
    int i=0;
    int len=0;
    char stack[12];
    char rstack[12];
    while (key>0) {
        stack[len] = '0'+key%10;
        len++;
        key /= 10;
    }
    i=0;
    while (i<len) {
        rstack[i] = stack[len-i-1];
        i++;
    }
    SAMIUint8ConcatString(data,curSize,rstack,len);
    return len;
}

__attribute__((target(mic)))
void SAMIUint8ConcatString(uint8_t * data, int * curSize,
                                char * key, int len) {
    int i;
    for (i=0;i<len;i++) {
        data[(*curSize)++] = key[i];
    }
}

inline int SAMExpandErrorsForCigar(unsigned int readLength, uint8_t readStrand, uint8_t errorCount, SRAError errors[], 
                                    int expandedIndelType[], int expandedIndelLength[], uint8_t errorOrder) {

    //Only insert and delete are selected from the list of errors in the list -- errors
    //These indels are then sorted into the sorted Indel list.
    SRAError sortedIndel[SRA_MAX_READ_LENGTH+1];
    int sortedIndelLen[SRA_MAX_READ_LENGTH];
    int indelCount = 0;
    
    //Above list will then be parsed and merged
    int mergedIndelLen = 0;
    int cigarTotalLength = 0;
        
    int i,j,step;
    int expandedIndelCount = 0;
    
    #ifdef SAM_DEBUG_PRINT_CIGAR_PROGRESS
        printf("Printing mismatch list before expansion,\n%d errors on a %ubp long read.\n",errorCount,readLength);
        for (i=0;i<errorCount;i++) {
            printf("%u:%d ",errors[i].position,errors[i].type);
        }
        printf("\n");
    #endif
    
    if (errorOrder==SAM_ERROR_IN_DESC && readStrand==QUERY_NEG_STRAND) {
        j = errorCount-1;
        step = -1;
    } else {
        j = 0;
        step = 1;
    }
    
    //ATTENTION: THIS LOOP IS MEANINGLESS. EW 20130923
    if (readStrand==QUERY_NEG_STRAND) {
        for (i=0;i<indelCount;i++) {
            sortedIndel[i].position = readLength-sortedIndel[i].position-1;
        }
    }
    
    for (i=0;i<errorCount;i++) {
        if (errors[j].type == SRA_CHAR_ERROR_TYPE_MISMATCH || errors[j].type == SRA_CHAR_ERROR_TYPE_NBMISMATCH) {
            j += step;
            continue;
        }
        sortedIndel[indelCount].type = errors[j].type;
        sortedIndel[indelCount].position = errors[j].position;
        indelCount++;
        
        j += step;
    }

    #ifdef SAM_DEBUG_PRINT_CIGAR_PROGRESS
        printf("Printing mismatch list after eliminating mismatch,\nexpect %d operations.\n",indelCount);
        for (i=0;i<indelCount;i++) {
            printf("%u:%d ",sortedIndel[i].position,sortedIndel[i].type);
        }
        printf("\n");
    #endif
    
    if (indelCount==0 || errorCount==0) {
        return 0;
    }
    
    #ifdef SAM_DEBUG_PRINT_CIGAR_PROGRESS
        printf("Printing mismatch list before sorting\n");
        for (i=0;i<indelCount;i++) {
            printf("%u:%d ",sortedIndel[i].position,sortedIndel[i].type);
        }
        printf("\n");
    #endif
    
    if (errorOrder == SAM_ERROR_RANDOM) {
        int swpTmp, swpPos;
        for (i=0;i<indelCount;i++) {
            for (j=0;j<i;j++) {
                if (sortedIndel[j].position > sortedIndel[j+1].position) {
                    swpTmp = sortedIndel[j].type;
                    swpPos = sortedIndel[j].position;
                    
                    sortedIndel[j].type = sortedIndel[j+1].type;
                    sortedIndel[j].position = sortedIndel[j+1].position;
                    
                    sortedIndel[j+1].type = swpTmp;
                    sortedIndel[j+1].position = swpPos;
                }
            }
        }
    }
    
    #ifdef SAM_DEBUG_PRINT_CIGAR_PROGRESS
        printf("Printing sorted indel only mismatches\n");
        for (i=0;i<indelCount;i++) {
            printf("%u:%d ",sortedIndel[i].position,sortedIndel[i].type);
        }
        printf("\n");
    #endif
    
    sortedIndelLen[mergedIndelLen]=1;
    for (i=1;i<indelCount;i++) {
        if (sortedIndel[i].type == SRA_CHAR_ERROR_TYPE_INSERT &&
            sortedIndel[i-1].type == sortedIndel[i].type &&
            sortedIndel[i-1].position == sortedIndel[i].position) {
            
            sortedIndelLen[mergedIndelLen]++;
        } else if (sortedIndel[i].type == SRA_CHAR_ERROR_TYPE_DELETE &&
            sortedIndel[i-1].type == sortedIndel[i].type &&
            sortedIndel[i-1].position == sortedIndel[i].position - 1) {
        
            sortedIndelLen[mergedIndelLen]++;
        } else if (sortedIndel[i].type == SRA_CHAR_ERROR_TYPE_SOFTCLIP &&
            sortedIndel[i-1].type == sortedIndel[i].type &&
            sortedIndel[i-1].position == sortedIndel[i].position - 1) {
        
            sortedIndelLen[mergedIndelLen]++;
        } else {
            mergedIndelLen++;
            sortedIndel[mergedIndelLen].type = sortedIndel[i].type;
            sortedIndel[mergedIndelLen].position = sortedIndel[i].position;
            sortedIndelLen[mergedIndelLen]=1;
        }
    }
    mergedIndelLen++;
    sortedIndel[mergedIndelLen].position = readLength;
    
    #ifdef SAM_DEBUG_PRINT_CIGAR_PROGRESS
        printf("Printing sorted indel after merging\n");
        for (i=0;i<mergedIndelLen;i++) {
            printf("%u:%d:%d ",sortedIndel[i].position,sortedIndel[i].type,sortedIndelLen[i]);
        }
        printf("\n");
    #endif
    
    if (sortedIndel[0].position>0) {
        expandedIndelType[expandedIndelCount] = SRA_CHAR_ERROR_TYPE_MISMATCH;
        expandedIndelLength[expandedIndelCount] = sortedIndel[0].position;
        cigarTotalLength+=sortedIndel[0].position;
        expandedIndelCount++;
    }
    
    for (i=0;i<mergedIndelLen;i++) {
        expandedIndelType[expandedIndelCount] = sortedIndel[i].type;
        expandedIndelLength[expandedIndelCount] = sortedIndelLen[i];
        cigarTotalLength+=(sortedIndel[i].type!=SRA_CHAR_ERROR_TYPE_INSERT)*expandedIndelLength[expandedIndelCount];
        expandedIndelCount++;
        
        if (sortedIndel[i+1].position - cigarTotalLength>0) {
            expandedIndelType[expandedIndelCount] = SRA_CHAR_ERROR_TYPE_MISMATCH;
            expandedIndelLength[expandedIndelCount] = sortedIndel[i+1].position - cigarTotalLength;
            cigarTotalLength+=expandedIndelLength[expandedIndelCount];
            expandedIndelCount++;
        }
    }
    
	if (sortedIndel[mergedIndelLen-1].type!=SRA_CHAR_ERROR_TYPE_INSERT &&
		sortedIndel[mergedIndelLen-1].position==readLength) {
        expandedIndelCount--;
	}
    
    #ifdef SAM_DEBUG_PRINT_CIGAR_PROGRESS
        printf("Printing sorted indel after expanding\n");
        for (i=0;i<expandedIndelCount;i++) {
            printf("%u:%d ",expandedIndelLength[i],expandedIndelType[i]);
        }
        printf("\n");
    #endif
    
    return expandedIndelCount;
}


int SAMIUint8ConcatCigar(uint8_t * data, int * curSize, 
                        unsigned int readLength, uint8_t readStrand, uint8_t errorCount, SRAError errors[], uint8_t errorOrder) {
    
    int i;
    
    //Then the following list be populated
    int expandedIndelType[SRA_MAX_READ_LENGTH*2+1];
    int expandedIndelLength[SRA_MAX_READ_LENGTH*2+1];
    int expandedIndelCount = SAMExpandErrorsForCigar(readLength,readStrand,errorCount,errors,expandedIndelType,expandedIndelLength,errorOrder);
    
    if (expandedIndelCount==0) {
        SAMIUint8ConcatUint32(data,curSize,readLength<<BAM_CIGAR_SHIFT);
        return 1;
    }
    
	for (i=0;i<expandedIndelCount;i++) {
        if (expandedIndelType[i] == SRA_CHAR_ERROR_TYPE_MISMATCH || expandedIndelType[i] == SRA_CHAR_ERROR_TYPE_NBMISMATCH) {
            SAMIUint8ConcatUint32(data,curSize,expandedIndelLength[i]<<BAM_CIGAR_SHIFT);
        } else if (expandedIndelType[i] == SRA_CHAR_ERROR_TYPE_INSERT) {
            SAMIUint8ConcatUint32(data,curSize,expandedIndelLength[i]<<BAM_CIGAR_SHIFT|BAM_CDEL);
        } else if (expandedIndelType[i] == SRA_CHAR_ERROR_TYPE_DELETE) {
            SAMIUint8ConcatUint32(data,curSize,expandedIndelLength[i]<<BAM_CIGAR_SHIFT|BAM_CINS);
        } else if (expandedIndelType[i] == SRA_CHAR_ERROR_TYPE_SOFTCLIP) {
            SAMIUint8ConcatUint32(data,curSize,expandedIndelLength[i]<<BAM_CIGAR_SHIFT|BAM_CSOFT_CLIP);
        }
    }
	
    return expandedIndelCount;
}


int SAMIUint8ConcatCigarAsString(uint8_t * data, int * curSize, 
                        unsigned int readLength, uint8_t readStrand, uint8_t errorCount, SRAError errors[], uint8_t errorOrder) {
    
    int i;
    //Then the following list be populated
    int expandedIndelType[SRA_MAX_READ_LENGTH*2+1];
    int expandedIndelLength[SRA_MAX_READ_LENGTH*2+1];
    int expandedIndelCount = SAMExpandErrorsForCigar(readLength,readStrand,errorCount,errors,expandedIndelType,expandedIndelLength,errorOrder);
    
    if (expandedIndelCount==0) {
        SAMIUint8ConcatUint32AsString(data,curSize,readLength);
        SAMIUint8ConcatString(data,curSize,"M",1);
        return 1;
    }
	
	for (i=0;i<expandedIndelCount;i++) {
        if (expandedIndelType[i] == SRA_CHAR_ERROR_TYPE_MISMATCH || expandedIndelType[i] == SRA_CHAR_ERROR_TYPE_NBMISMATCH) {
            SAMIUint8ConcatUint32AsString(data,curSize,expandedIndelLength[i]);
            SAMIUint8ConcatString(data,curSize,"M",1);
        } else if (expandedIndelType[i] == SRA_CHAR_ERROR_TYPE_INSERT) {
            SAMIUint8ConcatUint32AsString(data,curSize,expandedIndelLength[i]);
            SAMIUint8ConcatString(data,curSize,"D",1);
        } else if (expandedIndelType[i] == SRA_CHAR_ERROR_TYPE_DELETE) {
            SAMIUint8ConcatUint32AsString(data,curSize,expandedIndelLength[i]);
            SAMIUint8ConcatString(data,curSize,"I",1);
        } else if (expandedIndelType[i] == SRA_CHAR_ERROR_TYPE_SOFTCLIP) {
            SAMIUint8ConcatUint32AsString(data,curSize,expandedIndelLength[i]);
            SAMIUint8ConcatString(data,curSize,"S",1);
        }
    }
	
    return expandedIndelCount;
}


//SAMDPIUint8ConcatCigar expects errors reported by Dynamic Programming.
int SAMDPIUint8ConcatCigar(uint8_t * data, int * curSize, 
                        unsigned int readLength, uint8_t readStrand, uint16_t matchElemsCount, DPMatchElem matchElems[]) {
        
    int i,j,step;
    unsigned int expandedCigarString[DP_MAX_ERROR_COUNT];
    
    #ifdef SAM_DEBUG_PRINT_DP_CIGAR_PROGRESS
        printf("[SAMDPIUint8ConcatCigar] Printing sorted indel before expanding; \nthere are %d matching elements\n",matchElemsCount);
        for (i=0;i<matchElemsCount;i++) {
            printf("%u:%d ",matchElems[i].length, matchElems[i].type);
        }
        printf("\n");
    #endif
    
    //////////////////////////////////////////////////////
    // ATTENTION: Uncertain of how negative strand
    //            affects CIGAR string
    //////////////////////////////////////////////////////
    //if (readStrand==QUERY_NEG_STRAND) {
        j = matchElemsCount-1;
        step = -1;
    //} else {
    //    j = 0;
    //    step = 1;
    //}
    for (i=0;i<matchElemsCount;i++) {
        if (matchElems[j].type == SRA_CHAR_ERROR_TYPE_MISMATCH || matchElems[j].type == SRA_CHAR_ERROR_TYPE_NBMISMATCH) {
            expandedCigarString[i] = matchElems[j].length<<BAM_CIGAR_SHIFT;
        } else if (matchElems[j].type == SRA_CHAR_ERROR_TYPE_INSERT) {
            expandedCigarString[i] = matchElems[j].length<<BAM_CIGAR_SHIFT|BAM_CDEL;
        } else if (matchElems[j].type == SRA_CHAR_ERROR_TYPE_DELETE) {
            expandedCigarString[i] = matchElems[j].length<<BAM_CIGAR_SHIFT|BAM_CINS;
        } else if (matchElems[j].type == SRA_CHAR_ERROR_TYPE_SOFTCLIP) {
            expandedCigarString[i] = matchElems[j].length<<BAM_CIGAR_SHIFT|BAM_CSOFT_CLIP;
        }
        j+=step;
    }
    
    #ifdef SAM_DEBUG_PRINT_DP_CIGAR_PROGRESS
        printf("[SAMDPIUint8ConcatCigar] Printing sorted indel after expanding; \nthere are %d operations\n",matchElemsCount);
        for (i=0;i<matchElemsCount;i++) {
            int lastType = expandedCigarString[i] & ((1<<BAM_CIGAR_SHIFT)-1);
            if (lastType == 0) {
                lastType = SRA_CHAR_ERROR_TYPE_MISMATCH;
            } else if (lastType == BAM_CDEL) {
                lastType = SRA_CHAR_ERROR_TYPE_INSERT;
            } else if (lastType == BAM_CINS) {
                lastType = SRA_CHAR_ERROR_TYPE_DELETE;
            } else if (lastType == BAM_CSOFT_CLIP) {
                lastType = SRA_CHAR_ERROR_TYPE_SOFTCLIP;
            }
            printf("%u:%d ",expandedCigarString[i]>>BAM_CIGAR_SHIFT,lastType);
        }
        printf("\n");
    #endif
    
    int _curSize = (*curSize);
    memcpy(&(data[_curSize]),expandedCigarString,matchElemsCount*sizeof(unsigned int));
    (*curSize) += matchElemsCount*sizeof(unsigned int);
    
    return matchElemsCount;
}


//SAMDPIUint8ConcatCigarAsString expects errors reported by Dynamic Programming.
__attribute__((target(mic)))
int SAMDPIUint8ConcatCigarAsString(uint8_t * data, int * curSize, 
                        unsigned int readLength, uint8_t readStrand, uint16_t matchElemsCount, DPMatchElem matchElems[]) {
        
    int i,j,step;
    unsigned int expandedCigarString[DP_MAX_ERROR_COUNT];
    
    #ifdef SAM_DEBUG_PRINT_DP_CIGAR_PROGRESS
        printf("[SAMDPIUint8ConcatCigarAsString] Printing sorted indel before expanding; \nthere are %d matching elements\n",matchElemsCount);
        for (i=0;i<matchElemsCount;i++) {
            printf("%u:%d ",matchElems[i].length, matchElems[i].type);
        }
        printf("\n");
    #endif
    
    //////////////////////////////////////////////////////
    // ATTENTION: Uncertain of how negative strand
    //            affects CIGAR string
    //////////////////////////////////////////////////////
    //if (readStrand==QUERY_NEG_STRAND) {
        j = matchElemsCount-1;
        step = -1;
    //} else {
    //    j = 0;
    //    step = 1;
    //}
    for (i=0;i<matchElemsCount;i++) {
        if (matchElems[j].type == SRA_CHAR_ERROR_TYPE_MISMATCH || matchElems[j].type == SRA_CHAR_ERROR_TYPE_NBMISMATCH) {
            SAMIUint8ConcatUint32AsString(data,curSize,matchElems[j].length);
            SAMIUint8ConcatString(data,curSize,"M",1);
        } else if (matchElems[j].type == SRA_CHAR_ERROR_TYPE_INSERT) {
            SAMIUint8ConcatUint32AsString(data,curSize,matchElems[j].length);
            SAMIUint8ConcatString(data,curSize,"D",1);
        } else if (matchElems[j].type == SRA_CHAR_ERROR_TYPE_DELETE) {
            SAMIUint8ConcatUint32AsString(data,curSize,matchElems[j].length);
            SAMIUint8ConcatString(data,curSize,"I",1);
        } else if (matchElems[j].type == SRA_CHAR_ERROR_TYPE_SOFTCLIP) {
            SAMIUint8ConcatUint32AsString(data,curSize,matchElems[j].length);
            SAMIUint8ConcatString(data,curSize,"S",1);
        }
        j+=step;
    }
    
    #ifdef SAM_DEBUG_PRINT_DP_CIGAR_PROGRESS
        printf("[SAMDPIUint8ConcatCigarAsString] Printing sorted indel after expanding; \nthere are %d operations\n",matchElemsCount);
        for (i=0;i<matchElemsCount;i++) {
            int lastType = expandedCigarString[i] & ((1<<BAM_CIGAR_SHIFT)-1);
            if (lastType == 0) {
                lastType = SRA_CHAR_ERROR_TYPE_MISMATCH;
            } else if (lastType == BAM_CDEL) {
                lastType = SRA_CHAR_ERROR_TYPE_INSERT;
            } else if (lastType == BAM_CINS) {
                lastType = SRA_CHAR_ERROR_TYPE_DELETE;
            } else if (lastType == BAM_CSOFT_CLIP) {
                lastType = SRA_CHAR_ERROR_TYPE_SOFTCLIP;
            }
            printf("%u:%d ",expandedCigarString[i]>>BAM_CIGAR_SHIFT,lastType);
        }
        printf("\n");
    #endif
    
    return matchElemsCount;
}

// SAMSeqReadNameLength detemine the length of a string from the first character
// to the first separator, which could be a space, tab or any line breaker. 
int SAMSeqReadNameLength(char * str, int maxLen) {
    int i=0;
    while (str[i]!='\0' && i<maxLen) {
    
        if (str[i]==' ' ||
            str[i]=='\t' ||
            str[i]=='\r' ||
            str[i]=='\n') {
            
            break;
        }
        
        i++;
    }
    
    if (SAM_READ_NAME_IGNORE_PAIR_FLAG && i>=2) {
        if ((str[i-1]=='1' || str[i-1]=='2') && str[i-2]=='/') {
            return i-2;
        }
    }
    
    return i;
}

void SAMPackSequenceAndCigar(bam1_t * samAlgnmt, int sharedDataLen, SRAQueryInfo * qInfo, SRAOccurrence * occ) {
    int k;

    samAlgnmt->data_len = sharedDataLen;
    samAlgnmt->core.n_cigar = SAMIUint8ConcatCigar(samAlgnmt->data,&(samAlgnmt->data_len),
                            qInfo->ReadLength,
                            occ->strand,
                            occ->mismatchCount,
                            occ->errors, SAM_ERROR_RANDOM);          //CIGAR and number of CIGAR operations
    if ( occ->strand == QUERY_NEG_STRAND){
        // reverse Compliment for NEG STRAND
        for (k = qInfo->ReadLength - 1; k > 0; k-=2) {                                                                //Read
            unsigned char biChar = 0;
            biChar |= bam_nt16_table[dnaChar[3 - qInfo->ReportingReadCode[k]]] << 4;
            biChar |= bam_nt16_table[dnaChar[3 - qInfo->ReportingReadCode[k-1]]];
            SAMIUint8ConcatUint8(samAlgnmt->data,&(samAlgnmt->data_len),biChar);
        }
        if ( qInfo->ReadLength % 2 == 1 ) {
            uint8_t biChar = bam_nt16_table[dnaChar[3 - qInfo->ReportingReadCode[0]]] << 4;
            SAMIUint8ConcatUint8(samAlgnmt->data,&(samAlgnmt->data_len),biChar);
        }

        for (k =  qInfo->ReadLength - 1; k >= 0; k--) {                                                                  //Quality
            SAMIUint8ConcatUint8(samAlgnmt->data,&(samAlgnmt->data_len), qInfo->ReadQuality[0]==0? 0xff : qInfo->ReadQuality[k] - 33);
       // SAMIUint8ConcatUint8(samAlgnmt->data,&(samAlgnmt->data_len), 0xff);
        }
    }else{
        for (k = 0; k < qInfo->ReadLength/2; k++) {                                                                //Read
            unsigned char biChar = 0;
            biChar |= bam_nt16_table[dnaChar[qInfo->ReportingReadCode[k*2]]] << 4;
            biChar |= bam_nt16_table[dnaChar[qInfo->ReportingReadCode[k*2+1]]];
            SAMIUint8ConcatUint8(samAlgnmt->data,&(samAlgnmt->data_len),biChar);
        }
        if ( qInfo->ReadLength % 2 == 1 ) {
            uint8_t biChar = bam_nt16_table[dnaChar[qInfo->ReportingReadCode[qInfo->ReadLength-1]]] << 4;
            SAMIUint8ConcatUint8(samAlgnmt->data,&(samAlgnmt->data_len),biChar);
        }

        for (k = 0; k < qInfo->ReadLength; k++) {                                                                  //Quality
            SAMIUint8ConcatUint8(samAlgnmt->data,&(samAlgnmt->data_len), qInfo->ReadQuality[0]==0 ? 0xff : qInfo->ReadQuality[k] - 33 );
       // SAMIUint8ConcatUint8(samAlgnmt->data,&(samAlgnmt->data_len), 0xff);
        }

    }
}

void SAMDPPackSequenceAndCigar(bam1_t * samAlgnmt, int sharedDataLen, SRAQueryInfo * qInfo, DPOccurrence * dpOcc) {
    int k;
    
    samAlgnmt->data_len = sharedDataLen;
    samAlgnmt->core.n_cigar = SAMDPIUint8ConcatCigar(samAlgnmt->data,&(samAlgnmt->data_len),
                            qInfo->ReadLength,
                            dpOcc->strand,
                            dpOcc->matchElemsCount,
                            dpOcc->matchElems);          //CIGAR and number of CIGAR operations
                            
    unsigned char readBuffer[SRA_MAX_READ_LENGTH+2];
    unsigned char qualBuffer[SRA_MAX_READ_LENGTH+2];
    int readLength = qInfo->ReadLength;
    unsigned char * readCode = qInfo->ReportingReadCode;
    int byteCount = (readLength+1)/2;
    if ( dpOcc->strand == QUERY_NEG_STRAND){
        // NEG STRAND
        int Position_counter = 0;
        for (k = readLength - 1; k > 0; k-=2) {
            unsigned char biChar = 0;
            biChar |= bam_nt16_table[dnaChar[3 - readCode[k]]] << 4 | bam_nt16_table[dnaChar[3 - readCode[k-1]]];
            readBuffer[Position_counter++] = biChar;
        }

        for (k = readLength - 1; k >= 0; k--) {
            qualBuffer[readLength - 1 - k] = qInfo->ReadQuality[0]==0? 0xff : qInfo->ReadQuality[k] - 33 ;
            // qualBuffer[k] = 0xff;
        }
    }else{
        // POS STRAND
        for (k = 0; k < byteCount; k++) {
            unsigned char biChar = 0;
            biChar |= bam_nt16_table[dnaChar[readCode[k*2]]] << 4 | bam_nt16_table[dnaChar[readCode[k*2+1]]];
            readBuffer[k] = biChar;
        }

        for (k = 0; k < readLength; k++) {
            qualBuffer[k] = qInfo->ReadQuality[0]==0? 0xff : qInfo->ReadQuality[k] - 33 ;
            // qualBuffer[k] = 0xff;
        }
    }
    memcpy(&(samAlgnmt->data[samAlgnmt->data_len]),readBuffer,byteCount);
    samAlgnmt->data_len += byteCount;
    memcpy(&(samAlgnmt->data[samAlgnmt->data_len]),qualBuffer,readLength);
    samAlgnmt->data_len += readLength;
    /*for (k = 0; k < qInfo->ReadLength/2; k++) {                                                                //Read
        unsigned char biChar = 0;
        biChar |= bam_nt16_table[dnaChar[qInfo->ReportingReadCode[k*2]]] << 4;
        biChar |= bam_nt16_table[dnaChar[qInfo->ReportingReadCode[k*2+1]]];
        SAMIUint8ConcatUint8(samAlgnmt.data,&(samAlgnmt.data_len),biChar);
    }
    if ( qInfo->ReadLength % 2 == 1 ) {
        uint8_t biChar = bam_nt16_table[dnaChar[qInfo->ReportingReadCode[qInfo->ReadLength-1]]] << 4;
        SAMIUint8ConcatUint8(samAlgnmt.data,&(samAlgnmt.data_len),biChar);
    }
    for (k = 0; k < qInfo->ReadLength; k++) {                                                                  //Quality                
        SAMIUint8ConcatUint8(samAlgnmt.data,&(samAlgnmt.data_len),0xff);
    }*/
}

// old SAMOutputHeaderConstruct is kept for MICA-SE usage
void SAMOutputHeaderConstruct(bam_header_t * sheader, HSP * hsp) {
    unsigned short * ambiguityMap = hsp->ambiguityMap;
    Translate * translate = hsp->translate;
    unsigned int tp, approxIndex, approxValue;
    int i,j;
    
    //Number of sequences
    sheader->n_targets = hsp->numOfSeq;
    sheader->target_name = (char**) malloc(sizeof(char*)*hsp->numOfSeq);
    sheader->target_len = (uint32_t*) malloc(sizeof(uint32_t)*hsp->numOfSeq);
    sheader->text = (char*) malloc(sizeof(char)*1024);

    
    for (i=0;i<hsp->numOfSeq;i++) {
        j = SAMSeqReadNameLength(hsp->annotation[i].text,SAM_MAX_CHROMOSOME_NAME_LENGTH);
        hsp->annotation[i].text[j]='\0';

        //Fill in the Ref.Seq names and lengths.
        sheader->target_name[i] = hsp->annotation[i].text;
        sheader->target_len[i] = hsp->seqActualOffset[i].endPos - hsp->seqActualOffset[i].startPos;
    }
    sprintf(sheader->text,"@PG\tID:%s\tVN:v%d.%d.%d (%s)\n",PROJECT_NAME,PROJECT_MAJOR,PROJECT_MINOR,PROJECT_REV,PROJECT_SPECIAL);
    sheader->text = sheader->text;
    sheader->l_text=strlen(sheader->text);

    //Given up unknown parameters.
    //If someone ever found out what the hell is this pls update.
    sheader->hash = NULL;
    sheader->rg2lib = NULL;
}

// new SAMOutputHeaderConstruct containing listConsumer information like readGroup, sampleName and readGrpOption
// used by mica-pe for new batch list
void SAMOutputHeaderConstruct_ListReader(bam_header_t * sheader, HSP * hsp, ListReader * listreader) {
    unsigned short * ambiguityMap = hsp->ambiguityMap;
    Translate * translate = hsp->translate;
    unsigned int tp, approxIndex, approxValue;
    int i,j;
    
    //Number of sequences
    sheader->n_targets = hsp->numOfSeq;
    sheader->target_name = (char**) malloc(sizeof(char*)*hsp->numOfSeq);
    sheader->target_len = (uint32_t*) malloc(sizeof(uint32_t)*hsp->numOfSeq);
    sheader->text = (char*) malloc(sizeof(char)*1024);

    //char * snTags = ( char * ) malloc ( MAX_CHROMOSOME_NUM * MAX_SEQ_NAME_LENGTH * 2 ); // maximum chromosomes * (number of chars max each*2)
    //memset ( snTags, '\0', MAX_CHROMOSOME_NUM * MAX_SEQ_NAME_LENGTH * 2 );
    //char * snTagsTmp = ( char * ) malloc ( MAX_SEQ_NAME_LENGTH * 2 );
    //memset ( snTagsTmp, '\0', MAX_SEQ_NAME_LENGTH * 2 );
    
    for (i=0;i<hsp->numOfSeq;i++) {
        j = SAMSeqReadNameLength(hsp->annotation[i].text,SAM_MAX_CHROMOSOME_NAME_LENGTH);
        hsp->annotation[i].text[j]='\0';

        //Fill in the Ref.Seq names and lengths.
        sheader->target_name[i] = hsp->annotation[i].text;
        sheader->target_len[i] = hsp->seqActualOffset[i].endPos - hsp->seqActualOffset[i].startPos;
       // sprintf ( snTagsTmp, "@SQ\tSN:%s\tLN:%u\n", hsp->annotation[i].text, hsp->seqActualOffset[i].endPos - hsp->seqActualOffset[i].startPos );
       // strcat ( snTags, snTagsTmp );
    }

    char optional_field_buffer[1024];
    if ( ( listreader->sampleName != NULL && strlen ( listreader->sampleName ) != 0 ) || ( listreader->readGroup != NULL && strlen ( listreader->readGroup ) )
         || ( listreader->readGrpOption != NULL && strlen ( listreader->readGrpOption ) ) ){
        strcpy( optional_field_buffer , "@RG");
        if ( listreader->readGroup != NULL && strlen ( listreader->readGroup ) != 0 ){
            sprintf( optional_field_buffer, "%s\tID:%s", optional_field_buffer , listreader->readGroup );
        }
        if ( listreader->sampleName != NULL && strlen ( listreader->sampleName ) != 0 ){
            sprintf( optional_field_buffer, "%s\tSM:%s", optional_field_buffer , listreader->sampleName );
        }
        if ( listreader->readGrpOption != NULL && strlen ( listreader->readGrpOption ) != 0 ){
            sprintf( optional_field_buffer, "%s\t%s", optional_field_buffer , listreader->readGrpOption );
        }
        strcat( optional_field_buffer , "\n");
    } else {
        optional_field_buffer[0]='\0';
    }
    char programInfo[1024];
    sprintf(programInfo,"@PG\tID:%s\tVN:v%d.%d.%d (%s)\n",PROJECT_NAME,PROJECT_MAJOR,PROJECT_MINOR,PROJECT_REV,PROJECT_SPECIAL);
    int textLen = strlen ( optional_field_buffer ) + strlen ( programInfo )  + 15 + 23;
    sheader->text = ( char * ) malloc ( textLen * sizeof ( char ) );
    sprintf ( sheader->text, "@HD\tVN:1.3\tSO:unsorted\n%s%s", optional_field_buffer, programInfo );
    sheader->l_text = strlen ( sheader->text );

    //Given up unknown parameters.
    //If someone ever found out what the hell is this pls update.
    sheader->hash = NULL;
    sheader->rg2lib = NULL;
}

void SAMOutputHeaderDestruct(bam_header_t * sheader) {
    if ( sheader->text != NULL){
        free(sheader->text);
        sheader->text = NULL;
    }
    if ( sheader->target_name != NULL){
        free(sheader->target_name);
        sheader->target_name = NULL;
    }
    if ( sheader->target_len != NULL){
        free(sheader->target_len);
        sheader->target_len = NULL;
    }
}

