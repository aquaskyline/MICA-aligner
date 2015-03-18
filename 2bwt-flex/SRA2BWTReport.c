//
//    SRA2BWTReport.c
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
   
   Date   : 3rd July 2011
   Author : Edward MK Wu
   Change : Add support to pair-ending alignment mode.
            In this mode, the alignment output is not translated to
            (chromosome id/offset) tuple and output as 
            (mismatch count+1/ambiguous position). occMismatch is 
            incremented to avoid colliding with the delimitors.

   Date   : 24th July 2011
   Author : Edward MK Wu
   Change : Add Output Format support.

   Date   : 20th August 2011
   Author : Edward MK Wu
   Change : Created Output Format 20110820
            - Added a header length count (32-bit integer) to 
              indicate the length of the header. This might 
              possibly introduce backward compatibility of future
              output file format.
            - Added support to multi-threading output dump.

   Date   : 4th September 2011
   Author : Edward MK Wu
   Change : Added support to DEFAULT output format for debug.
   
*/
/////////////////////////////////////////////////////

#include "SRA2BWTReport.h"

// 2BWT-PE CTC#1
// -------------
// Uncomment this parameter SRA_OCC_PLAIN_DUMP_MISMATCH_POSITION if you want
// 2BWT-Combine to print the exact position of each mismatch
// found in individual alignment in Plain output mode.
//#define SRA_OCC_PLAIN_DUMP_MISMATCH_POSITION
// -------------

unsigned int OCCWriteOutputHeader(SRAOutput * aOutput, SRASetting * aSetting, SRAIndex * aIndex) {
    
    FILE * outFilePtr = aOutput->OutFilePtr;
    HSP * hsp = aIndex->hsp;
    unsigned short * ambiguityMap = hsp->ambiguityMap;
    Translate * translate = hsp->translate;
    unsigned long long tp;
    unsigned int tp_32;
    uint8_t chromId;
    int i;

    time_t mytime = time(NULL);
        
    fprintf(outFilePtr,"%s v%d.%d.%d (%s) OUTPUT\n",PROJECT_NAME,PROJECT_MAJOR,PROJECT_MINOR,PROJECT_REV,PROJECT_SPECIAL);
    fprintf(outFilePtr,"FORMAT_STR=%u ",OCC_OUTPUT_FORMAT);
    fprintf(outFilePtr,"CALLTIME=%s",ctime(&mytime));
    
    fprintf(outFilePtr,"%u\n",aSetting->OutputFormat);
    
    fprintf(outFilePtr,"%u\n",hsp->numOfSeq);
    for (i=0;i<hsp->numOfSeq;i++) {
        OCCTranslateOccurrence(translate, ambiguityMap, (unsigned long long) hsp->seqOffset[i].endPos, &chromId, &tp);
        tp_32 = tp;
        fprintf(outFilePtr,"%u %u %s\n",hsp->annotation[i].gi,tp_32,hsp->annotation[i].text);
    }
                
    return 0;
}

void OCCFlushCachePlain(SRAArguments * aArgs) {
    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRAOutput * aOutput = aArgs->AlgnmtOutput;
    SRAIndex * aIndex = aArgs->AlgnmtIndex;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    
    SRAOCCCollector * occCollector = aOutput->occCollector;
    FILE * outFilePtr = aOutput->OutFilePtr;
    HSP * hsp = aIndex->hsp;
    Translate * occTranslate = hsp->translate;
    unsigned short * occAmbiguityMap = hsp->ambiguityMap;
    
    DPArguments * dpArgs = aArgs->dpArguments;
    DPOCCCollector * dpOccCollector = dpArgs->dpOccCollector;
    DPOCCCollectorPointer * dpOccIterator = DPOCCPTCreate(dpOccCollector);
    
    DPOccurrence * dpOcc_1;
    
    char strandStr[4] = "?+-";
    SRAOccurrence lastDelimEntry;
    int lastDelimValid;
    uint8_t chromId;
    unsigned long long tp;
    
    SRAOccurrence_ll * myNode = occCollector->occRootNode;
    if (myNode==NULL || myNode->size==0) return;
    unsigned int i = 0;
    
    #ifdef SRA_OCC_PLAIN_DUMP_MISMATCH_POSITION
        int j;
    #endif
    
    while (myNode != NULL) {
        for (i=0;i<myNode->size;i++) {
            if (myNode->bucket[i].type==SRAOCC_TYPE_DELIMITOR_READ) {
            
                lastDelimEntry.type             = myNode->bucket[i].type;
                lastDelimEntry.ambPosition      = myNode->bucket[i].ambPosition;
                lastDelimEntry.strand           = myNode->bucket[i].strand;
                lastDelimEntry.mismatchCount    = myNode->bucket[i].mismatchCount;
                lastDelimValid = 1;
                
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_NO_ALIGNMENT) {
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_AWAIT_TRANSLATE) {
            
                OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                        myNode->bucket[i].ambPosition,
                                        &chromId,
                                        &tp);
                                        
                #ifndef DEBUG_2BWT_NO_WRITE_FILE
                fprintf(outFilePtr,"%llu %u %llu %c %d",  lastDelimEntry.ambPosition,
                                                            chromId,
                                                            tp,
                                                            strandStr[myNode->bucket[i].strand],
                                                            myNode->bucket[i].mismatchCount);
                                                            
                #ifdef SRA_OCC_PLAIN_DUMP_MISMATCH_POSITION
                    fprintf(outFilePtr," %d SRA",myNode->bucket[i].matchLen);
                    for (j=0;j<myNode->bucket[i].mismatchCount;j++) {
                        if (myNode->bucket[i].errors[j].type==SRA_CHAR_ERROR_TYPE_MISMATCH) {
                            fprintf(outFilePtr," M%d",myNode->bucket[i].errors[j].position);
                        } else if (myNode->bucket[i].errors[j].type==SRA_CHAR_ERROR_TYPE_INSERT) {
                            fprintf(outFilePtr," I%d",myNode->bucket[i].errors[j].position);
                        } else if (myNode->bucket[i].errors[j].type==SRA_CHAR_ERROR_TYPE_DELETE) {
                            fprintf(outFilePtr," D%d",myNode->bucket[i].errors[j].position);
                        } else if (myNode->bucket[i].errors[j].type==SRA_CHAR_ERROR_TYPE_NBMISMATCH) {
                            fprintf(outFilePtr," N%d",myNode->bucket[i].errors[j].position);
                        } else {
                            fprintf(outFilePtr," U%d",myNode->bucket[i].errors[j].position);
                        }
                    }
                #endif
                
                fprintf(outFilePtr,"\n");
                #endif
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_DP_SEED_OCC) {

                dpOcc_1 = DPOCCPTRead(dpOccIterator);DPOCCPTNext(dpOccIterator);
                
                if (dpOcc_1==NULL) {
                    #ifdef PE_OCC_PRINT_OVERFLOW_MESSAGE
                        printf("[PEOccFlushCache/Read#%llu] DP Occurrence Collector does not have this occurrence.\n",qInfo->ReadId);
                    #endif
                } else {
                    
                    OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                            dpOcc_1->ambPosition,
                                            &chromId,
                                            &tp);
                                            
                    #ifndef DEBUG_2BWT_NO_WRITE_FILE
                    fprintf(outFilePtr,"%llu %u %llu %c %d",  lastDelimEntry.ambPosition,
                                                                chromId,
                                                                tp,
                                                                strandStr[dpOcc_1->strand],
                                                                dpOcc_1->mismatchCount);

                    #ifdef SRA_OCC_PLAIN_DUMP_MISMATCH_POSITION
                        fprintf(outFilePtr," %d DP",dpOcc_1->matchLen);
                        for (j=dpOcc_1->matchElemsCount-1;j>=0;j--) {
                            if (dpOcc_1->matchElems[j].type==SRA_CHAR_ERROR_TYPE_MISMATCH) {
                                fprintf(outFilePtr," %dM",dpOcc_1->matchElems[j].length);
                            } else if (dpOcc_1->matchElems[j].type==SRA_CHAR_ERROR_TYPE_INSERT) {
                                fprintf(outFilePtr," %dI",dpOcc_1->matchElems[j].length);
                            } else if (dpOcc_1->matchElems[j].type==SRA_CHAR_ERROR_TYPE_DELETE) {
                                fprintf(outFilePtr," %dD",dpOcc_1->matchElems[j].length);
                            } else if (dpOcc_1->matchElems[j].type==SRA_CHAR_ERROR_TYPE_NBMISMATCH) {
                                fprintf(outFilePtr," %dN",dpOcc_1->matchElems[j].length);
                            } else {
                                fprintf(outFilePtr," %dU",dpOcc_1->matchElems[j].length);
                            }
                        }
                    #endif
                    
                    fprintf(outFilePtr,"\n");
                    #endif
                    
                }
            }
        }
        myNode = myNode->next;
    }
    
    DPOCCPTFree(dpOccIterator);
    SRAOCCInitialise(occCollector);
    DPOCCInitialise(dpOccCollector);
    
    if (lastDelimValid) {
        SRAOCCAddOccurrenceToBuckets(occCollector, &lastDelimEntry);
    }
}

void OCCFlushCacheSAM(SRAArguments * aArgs) {
    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRAOutput * aOutput = aArgs->AlgnmtOutput;
    SRAIndex * aIndex = aArgs->AlgnmtIndex;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    
    SRAOCCCollector * occCollector = aOutput->occCollector;
    FILE * outFilePtr = aOutput->OutFilePtr;
    HSP * hsp = aIndex->hsp;
    Translate * occTranslate = hsp->translate;
    unsigned short * occAmbiguityMap = hsp->ambiguityMap;

    DPArguments * dpArgs = aArgs->dpArguments;
    DPOCCCollector * dpOccCollector = dpArgs->dpOccCollector;
    DPOCCCollectorPointer * dpOccIterator = DPOCCPTCreate(dpOccCollector);
    
    DPOccurrence * dpOcc_1;
    
    samfile_t * samFilePtr = aOutput->SAMOutFilePtr;
    
    bam1_t samAlgnmt;
    
    unsigned long long i=0, k=0;
    uint8_t chromId;
    unsigned long long tp;
    
    if (qInfo->ReportingReadCode==NULL) { qInfo->ReportingReadCode = qInfo->ReadCode; }
    
    SRAOccurrence_ll * myNode = occCollector->occRootNode;
    if (myNode==NULL || myNode->size==0) return;
    
    unsigned int readNumOfOcc = SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_AWAIT_TRANSLATE);
    unsigned int readNumOfDpOcc = SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_PE_DP_SEED_OCC);
    
    unsigned long long samAuxMultiOccBlockSize = SAM_MDATA_SIZE_PER_READ + 
                                                (readNumOfOcc) * (SAM_MDATA_SIZE_PER_OCC) + 
                                                (readNumOfDpOcc) * (SAM_MDATA_SIZE_PER_DPOCC);
    
    int samPrimaryOccAvailable = 0;
    int samXAZtagCreated = 0;
    
    if (samAuxMultiOccBlockSize>aOutput->SAMAuxDataBlockSize) {
    
        aOutput->SAMAuxDataBlock = (uint8_t *)realloc(aOutput->SAMAuxDataBlock,
                                                        sizeof(uint8_t)*samAuxMultiOccBlockSize);
        aOutput->SAMAuxDataBlockSize = samAuxMultiOccBlockSize;
    }
    
    //Assuming the answer being handled by OCCFlushCacheSAMAPI are all answers of qInfo
    // Each occurrences reported by SAM_API will share the following aux data.
    //----------------------------------------------------------------------------------------------------------------------------
        samAlgnmt.core.bin = bam_reg2bin(0,0);                     //SAM: bin calculated by bam_reg2bin()
        samAlgnmt.core.l_qseq = qInfo->ReadLength;                 //SAM: length of the query sequence (read)
        samAlgnmt.core.l_qname = SAMSeqReadNameLength(qInfo->ReadName,MAX_SEQ_NAME_LENGTH)+1;        //SAM: length of the query name
        
        samAlgnmt.core.mtid = -1;
        samAlgnmt.core.mpos = -1;
        samAlgnmt.core.isize = 0;
        
        samAlgnmt.l_aux = 0;
        samAlgnmt.data_len = 0;
        samAlgnmt.data = aOutput->SAMAuxDataBlock;
        samAlgnmt.m_data = samAuxMultiOccBlockSize;
        SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),qInfo->ReadName,samAlgnmt.core.l_qname);    //Name
        int sharedDataLen = samAlgnmt.data_len;
    //----------------------------------------------------------------------------------------------------------------------------                      
    
    while (myNode != NULL) {
        for (i=0;i<myNode->size;i++) {
            if (myNode->bucket[i].type==SRAOCC_TYPE_DELIMITOR_READ) {
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_DELIMITOR_STEP) {
            } else if (myNode->bucket[k].type==SRAOCC_TYPE_NO_ALIGNMENT) {
            
                samAlgnmt.data_len = sharedDataLen;
                SAMPackSequenceAndCigar(&samAlgnmt,sharedDataLen,qInfo,&(myNode->bucket[i]));
                samAlgnmt.l_aux = 0;
                
                unsigned int samFlag = 0;
                samFlag |= SAM_FLAG_READ_UNMAPPED;
            
                // ---------------------->  |
                                            samAlgnmt.core.tid = -1;                                   //SAM: chromosome ID, defined by bam_header_t
                                            samAlgnmt.core.pos = -1;                                   //SAM: 0-based leftmost coordinate4
                                            samAlgnmt.core.qual = SAM_MAPQ_UNALIGNED;                  //SAM: mapping quality
                                            samAlgnmt.core.flag = samFlag;                             //SAM: bitwise flag
                                            
                                            #ifndef DEBUG_2BWT_NO_WRITE_FILE
                                            samwrite(samFilePtr,&samAlgnmt);
                                            #endif
                // ---------------------->  |
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_AWAIT_TRANSLATE) {
        
                OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                        myNode->bucket[i].ambPosition,
                                        &chromId,
                                        &tp);

                if (samPrimaryOccAvailable) {
                    unsigned int auxStart = samAlgnmt.data_len;
                    if (!samXAZtagCreated) {
                        SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),"XAZ",3);
                        samXAZtagCreated=1;
                    }
                    SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),
                                            hsp->annotation[chromId-1].text,
                                            SAMSeqReadNameLength(hsp->annotation[chromId-1].text,MAX_SEQ_NAME_LENGTH));
                    if (myNode->bucket[i].strand==QUERY_NEG_STRAND) {
                        SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",-",2);
                    } else if (myNode->bucket[i].strand==QUERY_POS_STRAND) {
                        SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",+",2);
                    } else {
                        SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",?",2);
                    }
                    SAMIUint8ConcatUint32AsString(samAlgnmt.data,&(samAlgnmt.data_len),tp);
                    SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",",1);
                    SAMIUint8ConcatCigarAsString(samAlgnmt.data,&(samAlgnmt.data_len),
                                                qInfo->ReadLength,
                                                myNode->bucket[i].strand,
                                                myNode->bucket[i].mismatchCount,
                                                myNode->bucket[i].errors, SAM_ERROR_RANDOM);
                    SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",",1);
                    SAMIUint8ConcatUint32AsString(samAlgnmt.data,&(samAlgnmt.data_len),myNode->bucket[i].mismatchCount);
                    SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),";",1);
                        
                    samAlgnmt.l_aux += samAlgnmt.data_len - auxStart;
                    
                } else {

                    SAMPackSequenceAndCigar(&samAlgnmt,sharedDataLen,qInfo,&myNode->bucket[i]);
                    samAlgnmt.l_aux = 0;
                    
                    unsigned int samFlag = 0;
                    if (myNode->bucket[i].strand==QUERY_NEG_STRAND) {
                        samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
                    }
                    
                    // ---------------------->  |
                                                samAlgnmt.core.tid = chromId-1;                           //SAM: chromosome ID, defined by bam_header_t
                                                samAlgnmt.core.pos = tp-1;                                 //SAM: 0-based leftmost coordinate
                                                samAlgnmt.core.qual = SAM_MAPQ_UNAVAILABLE;                //SAM: mapping quality
                                                samAlgnmt.core.flag = samFlag;                             //SAM: bitwise flag

                                                samPrimaryOccAvailable = 1;
                    // ---------------------->  |
                }
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_DP_SEED_OCC) {
            
                // 
                // The occurrences are found by PE-DP and the occurrences were stored in
                // dpOcc_1 and dpOcc_2.
                // The occurrences will then be stored into samAlgnmt and samAlgnmt_mate.
                //
                
                dpOcc_1 = DPOCCPTRead(dpOccIterator);DPOCCPTNext(dpOccIterator);
                
                if (dpOcc_1==NULL) {
                    #ifdef SRA_OCC_PRINT_OVERFLOW_MESSAGE
                        printf("[PEOccFlushCache/Read#%llu] DP Occurrence Collector does not have this occurrence.\n",qInfo->ReadId);
                    #endif
                } else {
                
                    OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                            dpOcc_1->ambPosition,
                                            &chromId,
                                            &tp);
                                            
                    if (samPrimaryOccAvailable) {

                        //--------------------------
                            unsigned int auxStart = samAlgnmt.data_len;
                            if (!samXAZtagCreated) {
                                SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),"XAZ",3);
                                samXAZtagCreated=1;
                            }
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),
                                                    hsp->annotation[chromId-1].text,
                                                    SAMSeqReadNameLength(hsp->annotation[chromId-1].text,MAX_SEQ_NAME_LENGTH));
                            if (dpOcc_1->strand==QUERY_NEG_STRAND) {
                                SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",-",2);
                            } else if (dpOcc_1->strand==QUERY_POS_STRAND) {
                                SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",+",2);
                            } else {
                                SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",?",2);
                            }
                            SAMIUint8ConcatUint32AsString(samAlgnmt.data,&(samAlgnmt.data_len),tp);
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",",1);
                            SAMDPIUint8ConcatCigarAsString(samAlgnmt.data,&(samAlgnmt.data_len),
                                                        qInfo->ReadLength,
                                                        dpOcc_1->strand,
                                                        dpOcc_1->matchElemsCount,
                                                        dpOcc_1->matchElems);
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",",1);
                            SAMIUint8ConcatUint32AsString(samAlgnmt.data,&(samAlgnmt.data_len),dpOcc_1->mismatchCount);
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),";",1);
                                
                            samAlgnmt.l_aux += samAlgnmt.data_len - auxStart;
                        //--------------------------
                    } else {

                        SAMDPPackSequenceAndCigar(&samAlgnmt,sharedDataLen,qInfo,dpOcc_1);
                        samAlgnmt.l_aux = 0;
                                                
                        unsigned int samFlag = 0;
                        if (dpOcc_1->strand==QUERY_NEG_STRAND) {
                            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
                        }
                        
                        // ---------------------->  |
                                                    samAlgnmt.core.tid = chromId-1;                           //SAM: chromosome ID, defined by bam_header_t
                                                    samAlgnmt.core.pos = tp-1;                                 //SAM: 0-based leftmost coordinate
                                                    samAlgnmt.core.qual = SAM_MAPQ_UNAVAILABLE;                //SAM: mapping quality
                                                    samAlgnmt.core.flag = samFlag;                             //SAM: bitwise flag

                                                    samPrimaryOccAvailable = 1;
                        // ---------------------->  |
                    }
                }
            
            }
        }
        myNode = myNode->next;
    }

    if (samPrimaryOccAvailable) {
        unsigned int auxStart = samAlgnmt.data_len;
        SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),"\0",1);
        samAlgnmt.l_aux += samAlgnmt.data_len - auxStart;
        #ifndef DEBUG_2BWT_NO_WRITE_FILE
        samwrite(samFilePtr,&samAlgnmt);
        #endif
    }
    DPOCCPTFree(dpOccIterator);
    SRAOCCInitialise(occCollector);
    DPOCCInitialise(dpOccCollector);
}

void OCCFlushCache(SRAArguments * aArgs) {
    SRAOutput * aOutput = aArgs->AlgnmtOutput;
    SRAOCCCollector * occCollector = aOutput->occCollector;

    switch (aOutput->OutFileFormat) {
        case SRA_OUTPUT_FORMAT_PLAIN:
            OCCFlushCachePlain(aArgs);
            break;
        case SRA_OUTPUT_FORMAT_SAM:
        case SRA_OUTPUT_FORMAT_BAM:
            OCCFlushCacheSAM(aArgs);
            break;
        case SRA_OUTPUT_FORMAT_SAM_STORE:
            break;
        default:
            printf("Occurrences are dropped as the output format is invalid\n");
            SRAOCCInitialise(occCollector);
    }
}

uint16_t OCCComputeMatchLength(uint16_t readLength,
                                uint8_t occMismatch,
                                const SRAError errorPos[],
                                SRAWorkingMemory * workMem) {

    int i;
    uint16_t matchLen = readLength;
    
    for (i=0;i<occMismatch;i++) {
        if (errorPos[i].type == SRA_CHAR_ERROR_TYPE_INSERT) {
            matchLen++;
        } else if (errorPos[i].type == SRA_CHAR_ERROR_TYPE_DELETE) {
            matchLen--;
        }
    }
    
    matchLen -= workMem->Shared_HeadTrim;
    matchLen -= workMem->Shared_TailTrim;
    
    //printf("[OCCComputeMatchLength] Trimmed = %d / %d\n",
    //    workMem->Shared_HeadTrim,workMem->Shared_TailTrim);
    
    return matchLen;
}

void OCCFlushBestIntoCache(SRAArguments * aArgs) {
    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAWorkingMemory * workMem  = aArgs->AlgnmtMemory;
    SRAResultCount * aStats = aArgs->AlgnmtStats;
    
    uint16_t matchLen = qInfo->ReadLength;
    
    matchLen = OCCComputeMatchLength(matchLen,workMem->BestAlignmentNumOfError,workMem->BestAlignmentErrorPos,workMem);
    
    OCCAddTextPositionToCache(aArgs,SRAOCC_TYPE_DELIMITOR_STEP,0,0,matchLen,
                                workMem->BestAlignmentNumOfError,0,
                                workMem->BestAlignmentErrorPos);
    
    OCCAddTextPositionToCache(aArgs,SRAOCC_TYPE_AWAIT_TRANSLATE,
                                workMem->BestAlignmentPosition,
                                workMem->BestAlignmentStrand,
                                matchLen,
                                workMem->BestAlignmentNumOfError,
                                workMem->BestAlignmentQuality,
                                workMem->BestAlignmentErrorPos);
                                
    aStats->WithError[workMem->BestAlignmentNumOfError]++;
}

//OCCReportSARange : asumming l<=r
unsigned long long OCCReportSARange(SRAArguments * aArgs, 
                                    unsigned long long l, unsigned long long r, 
                                    int occMismatch, int occQuality) {
    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAOutput * aOutput = aArgs->AlgnmtOutput;
    SRAIndex * aIndex = aArgs->AlgnmtIndex;
    SRAWorkingMemory * workMem  = aArgs->AlgnmtMemory;
    SRAResultCount * aStats = aArgs->AlgnmtStats;
    SRAOCCCollector * occCollector = aOutput->occCollector;
    int i;
    
#ifdef DEBUG_2BWT_NO_OUTPUT
    if (l>r) return 0;
    aStats->WithError[occMismatch] += r-l+1;
    return r-l+1;
#endif

    BWT * bwt = aIndex->bwt;
    HOCC * highOcc = aIndex->highOcc;
    uint8_t OutputType = aSetting->OutputType;
    unsigned long long j,k;
    unsigned long long saCount = r-l+1;
    uint8_t readStrand = qInfo->ReadStrand;

    if (OutputType == SRA_REPORT_UNIQUE_BEST) {
        if (l==r && workMem->IsOpened==0) {
            workMem->IsUnique = 1;
            workMem->IsOpened = 1;
            workMem->BestAlignmentPosition = BWTSaValue(bwt,l) - workMem->Shared_HeadTrim;
            workMem->BestAlignmentNumOfError = occMismatch;
            workMem->BestAlignmentStrand = readStrand;
            workMem->BestAlignmentQuality = occQuality;
            for (i=0;i<occMismatch;i++) {
                workMem->BestAlignmentErrorPos[i].type = workMem->Shared_AccErrorsPos[i].type;
                workMem->BestAlignmentErrorPos[i].position = workMem->Shared_AccErrorsPos[i].position;
            }
            
            aStats->RetrievedBySa++;
        } else {
            workMem->IsOpened = 1;
            workMem->IsUnique = 0;
            workMem->IsClosed = 1;
        }
        return 0;
    } else if (OutputType == SRA_REPORT_RANDOM_BEST) {
        if (workMem->IsClosed==1) return 0;
        workMem->BestAlignmentPosition = BWTSaValue(bwt,l) - workMem->Shared_HeadTrim;
        workMem->BestAlignmentNumOfError = occMismatch;
        workMem->BestAlignmentStrand = readStrand;
        workMem->BestAlignmentQuality = occQuality;
        for (i=0;i<occMismatch;i++) {
            workMem->BestAlignmentErrorPos[i].type = workMem->Shared_AccErrorsPos[i].type;
            workMem->BestAlignmentErrorPos[i].position = workMem->Shared_AccErrorsPos[i].position;
        }
        workMem->IsOpened = 1;
        workMem->IsClosed=1;
        
        aStats->RetrievedBySa++;
        return 1;
    } else if (OutputType == SRA_REPORT_BEST_QUALITY) {
        if (workMem->IsOpened==0 || 
           (workMem->IsOpened==1 && occQuality < workMem->BestAlignmentQuality)) {
            workMem->BestAlignmentPosition = BWTSaValue(bwt,l) - workMem->Shared_HeadTrim;
            workMem->BestAlignmentNumOfError = occMismatch;
            workMem->BestAlignmentStrand = readStrand;
            workMem->BestAlignmentQuality = occQuality;
            for (i=0;i<occMismatch;i++) {
                workMem->BestAlignmentErrorPos[i].type = workMem->Shared_AccErrorsPos[i].type;
                workMem->BestAlignmentErrorPos[i].position = workMem->Shared_AccErrorsPos[i].position;
            }
            workMem->IsOpened = 1;
            
            aStats->RetrievedBySa++;
        }
        return 0;
    } else if (OutputType == SRA_REPORT_ALL_BEST || OutputType == SRA_REPORT_ALL || OutputType == SRA_REPORT_ALL_SORTED) {
    
        /////////////////////////////////////////
        // Checking MaxResult setting
        unsigned int occCount = SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_AWAIT_TRANSLATE);
        unsigned long long ll = l;
        unsigned long long lr = r;
        unsigned long long lsaCount = saCount;

        if (aSetting->MaxResult!=-1) {
            if (occCount>aSetting->MaxResult) {
                workMem->IsClosed=1;
                return 0;
            } else if (occCount+saCount>aSetting->MaxResult) {
                lr = ll + (aSetting->MaxResult - occCount) - 1;
                lsaCount = lr - ll + 1;
                workMem->IsClosed=1;
            }
        }
        //
        /////////////////////////////////////////
    
        uint16_t matchLen = qInfo->ReadLength;
        
        matchLen = OCCComputeMatchLength(matchLen,occMismatch,workMem->Shared_AccErrorsPos,workMem);
            
        OCCAddTextPositionToCache(aArgs,SRAOCC_TYPE_DELIMITOR_STEP,0,0,0,occMismatch,0,workMem->Shared_AccErrorsPos);
        
        //Retrieve SA values from either SA or HOCC. Depending on the size of the SA Range.
        if (highOcc!=NULL && saCount>=4) {
            //By HOCC data structure.
            unsigned long long hOccCachedL,hOccCachedR,tp;
            unsigned long long approxL = HOCCApproximateRank(highOcc,ll+1);
            
            HOCCGetCachedSARange(highOcc, approxL,&hOccCachedL,&hOccCachedR);
            while (hOccCachedL>ll) {
                approxL--;
                HOCCGetCachedSARange(highOcc, approxL,&hOccCachedL,&hOccCachedR);
            }

            //Case 1: When the entire SA range falls within a cached range of HOOC
            if (hOccCachedL<=ll && hOccCachedR>=lr) {
                //Get the first index on HOCC
                unsigned long long hOccAnsIndex = highOcc->table[approxL*2+1];
                hOccAnsIndex+= (ll-hOccCachedL);

                j=0; k=hOccAnsIndex;
                while (j<lsaCount) {
                    OCCAddTextPositionToCache(aArgs,SRAOCC_TYPE_AWAIT_TRANSLATE,
                                                highOcc->occ[k] - workMem->Shared_HeadTrim,
                                                readStrand,matchLen,
                                                occMismatch,occQuality,
                                                workMem->Shared_AccErrorsPos);
                    k++;
                    j++;
                }
                aStats->RetrievedByHocc+=lsaCount;

            //Case 2: When the reported SA range enclosed some cached ranges of HOOC
            } else if (hOccCachedR<lr) {
                unsigned long long hOccAnsIndex = highOcc->table[approxL*2+1];
                j=ll; k=hOccAnsIndex;
                while (j<=lr) {
                    if (j>hOccCachedR && approxL<highOcc->ttlPattern) {
                        approxL++;
                        HOCCGetCachedSARange(highOcc, approxL,&hOccCachedL,&hOccCachedR);
                        hOccAnsIndex = highOcc->table[approxL*2+1];
                        k=hOccAnsIndex;
                    }

                    if (hOccCachedL<=j && j<=hOccCachedR) {
                        OCCAddTextPositionToCache(aArgs,SRAOCC_TYPE_AWAIT_TRANSLATE,
                                                    highOcc->occ[k] - workMem->Shared_HeadTrim,
                                                    readStrand,matchLen,
                                                    occMismatch,occQuality,
                                                    workMem->Shared_AccErrorsPos);
                        aStats->RetrievedByHocc++;
                        k++;
                    } else {
                        OCCAddTextPositionToCache(aArgs,SRAOCC_TYPE_AWAIT_TRANSLATE,
                                                    BWTSaValue(bwt,j) - workMem->Shared_HeadTrim,
                                                    readStrand,matchLen,
                                                    occMismatch,occQuality,
                                                    workMem->Shared_AccErrorsPos);
                        aStats->RetrievedBySa++;
                    }
                    j++;
                }
            //Case 3: Cannot find the records in HOCC (Should not happen)
            //In this case, we compute the occ by Suffix array
            } else {
                j=ll;
                while (j<=lr) {
                    OCCAddTextPositionToCache(aArgs,SRAOCC_TYPE_AWAIT_TRANSLATE,
                                                BWTSaValue(bwt,j) - workMem->Shared_HeadTrim,
                                                readStrand,matchLen,
                                                occMismatch,occQuality,
                                                workMem->Shared_AccErrorsPos);
                    j++;
                    aStats->RetrievedBySa++;
                }
            }
            //*/

        } else {
            //By Suffix array.
            j=ll;
            while (j<=lr) {
                OCCAddTextPositionToCache(aArgs,SRAOCC_TYPE_AWAIT_TRANSLATE,
                                        BWTSaValue(bwt,j) - workMem->Shared_HeadTrim,
                                        readStrand,matchLen,
                                        occMismatch,occQuality,
                                        workMem->Shared_AccErrorsPos);
                j++;
            }
            aStats->RetrievedBySa+=lsaCount;
        }
        aStats->WithError[occMismatch] += lsaCount;
        return lsaCount;
    
    } else if (OutputType == SRA_REPORT_NONE) {
        aStats->WithError[occMismatch] += saCount;
        return saCount;
    }

    return 0;
}

void OCCReportNoAlignment(SRAArguments * aArgs) {
    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRAOutput * aOutput = aArgs->AlgnmtOutput;
    
    switch (aOutput->OutFileFormat) {
        case SRA_OUTPUT_FORMAT_PLAIN:
        case SRA_OUTPUT_FORMAT_SAM:
        case SRA_OUTPUT_FORMAT_BAM:
        case SRA_OUTPUT_FORMAT_SAM_STORE:
            OCCAddTextPositionToCache(aArgs,SRAOCC_TYPE_NO_ALIGNMENT,
                                        qInfo->ReadId,
                                        0,0,
                                        0,0,
                                        NULL);
            break;
    }
    
}

void OCCReportDelimitor(SRAArguments * aArgs) {
    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRAOutput * aOutput = aArgs->AlgnmtOutput;
    
    switch (aOutput->OutFileFormat) {
        case SRA_OUTPUT_FORMAT_PLAIN: 
        case SRA_OUTPUT_FORMAT_SAM_STORE:
            OCCAddTextPositionToCache(aArgs,SRAOCC_TYPE_DELIMITOR_READ,
                                        qInfo->ReadId,
                                        0,0,
                                        0,0,
                                        NULL);
            break;
        case SRA_OUTPUT_FORMAT_SAM:
        case SRA_OUTPUT_FORMAT_BAM:
            break;
    }
}

unsigned long long OCCReportTextPositions(SRAArguments * aArgs, int posCount, SRAError errors[][MAX_NUM_OF_ERROR]) {
    
    SRAQueryInfo * qInfo = aArgs->QueryInfo;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAOutput * aOutput = aArgs->AlgnmtOutput;
    SRAIndex * aIndex = aArgs->AlgnmtIndex;
    SRAWorkingMemory * workMem  = aArgs->AlgnmtMemory;
    SRAResultCount * aStats = aArgs->AlgnmtStats;
    SRAOCCCollector * occCollector = aOutput->occCollector;

    int i;
    int k;
    
#ifdef DEBUG_2BWT_NO_OUTPUT
    for (i=0;i<posCount;i++) {
        aStats->WithError[workMem->Shared_AccMismatch[i]]++;
    }
    return posCount;
#endif

    if (posCount<=0) return 0;
   
    unsigned long long * saPositions = workMem->Shared_ReadPositions;
    uint8_t * occMismatches = workMem->Shared_AccMismatch;
    int * occQualities = workMem->Shared_AccQualities;
    
    uint8_t OutputType = aSetting->OutputType;
    FILE * outFilePtr = aOutput->OutFilePtr;
        
    if (OutputType == SRA_REPORT_UNIQUE_BEST) {
        if (posCount==1 && workMem->IsOpened == 0) {
            workMem->IsUnique = 1;
            workMem->BestAlignmentPosition = saPositions[0] - workMem->Shared_HeadTrim;
            workMem->BestAlignmentNumOfError = occMismatches[0];
            workMem->BestAlignmentQuality = occQualities[0];
            workMem->BestAlignmentStrand = qInfo->ReadStrand;
            for (i=0;i<occMismatches[0];i++) {
                workMem->BestAlignmentErrorPos[i].type = errors[0][i].type;
                workMem->BestAlignmentErrorPos[i].position = errors[0][i].position;
            }
            
            aStats->RetrievedByCe++;
        } else {
            workMem->IsUnique = 0;
            workMem->IsClosed = 1;
        }
        workMem->IsOpened = 1;
        return 0;
    } else if (OutputType == SRA_REPORT_RANDOM_BEST) {
        if (workMem->IsClosed==1) return 0;
        workMem->BestAlignmentPosition = saPositions[0] - workMem->Shared_HeadTrim;
        workMem->BestAlignmentNumOfError = occMismatches[0];
        workMem->BestAlignmentQuality = occQualities[0];
        workMem->BestAlignmentStrand = qInfo->ReadStrand;
        for (i=0;i<occMismatches[0];i++) {
            workMem->BestAlignmentErrorPos[i].type = errors[0][i].type;
            workMem->BestAlignmentErrorPos[i].position = errors[0][i].position;
        }
        workMem->IsOpened=1;
        workMem->IsClosed=1;
        
        aStats->RetrievedByCe++;
        return 1;
    } else if (OutputType == SRA_REPORT_BEST_QUALITY) {
        for (k=0;k<posCount;k++) {
            if (workMem->IsOpened==0 || 
                (workMem->IsOpened==1 && occQualities[k] < workMem->BestAlignmentQuality)) {
                workMem->BestAlignmentPosition = saPositions[k] - workMem->Shared_HeadTrim;
                workMem->BestAlignmentNumOfError = occMismatches[k];
                workMem->BestAlignmentQuality = occQualities[k];
                workMem->BestAlignmentStrand = qInfo->ReadStrand;
                for (i=0;i<occMismatches[k];i++) {
                    workMem->BestAlignmentErrorPos[i].type = errors[k][i].type;
                    workMem->BestAlignmentErrorPos[i].position = errors[k][i].position;
                }
                workMem->IsOpened = 1;
                
                aStats->RetrievedByCe++;
            }
        }
        return 0;
    } else if (OutputType == SRA_REPORT_ALL_BEST || OutputType == SRA_REPORT_ALL || OutputType == SRA_REPORT_ALL_SORTED) {
    
        /////////////////////////////////////////
        // Checking MaxResult setting
        unsigned int occCount = SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_AWAIT_TRANSLATE);
        unsigned long long lposCount = posCount;
        if (aSetting->MaxResult!=-1) {
            if (occCount>aSetting->MaxResult) {
                workMem->IsClosed=1;
                return 0;
            } else if (occCount+posCount>aSetting->MaxResult) {
                lposCount = aSetting->MaxResult - occCount;
                workMem->IsClosed=1;
            }
        }
        //
        /////////////////////////////////////////
    
        for (k=0;k<lposCount;k++) {
        
            uint16_t matchLen = qInfo->ReadLength;
            
            matchLen = OCCComputeMatchLength(matchLen,occMismatches[k],errors[k],workMem);
                
            OCCAddTextPositionToCache(aArgs,SRAOCC_TYPE_DELIMITOR_STEP,
                                        saPositions[k] - workMem->Shared_HeadTrim,
                                        qInfo->ReadStrand,matchLen,
                                        occMismatches[k],occQualities[k],
                                        errors[k]);
        
            OCCAddTextPositionToCache(aArgs,SRAOCC_TYPE_AWAIT_TRANSLATE,
                                        saPositions[k] - workMem->Shared_HeadTrim,
                                        qInfo->ReadStrand,matchLen,
                                        occMismatches[k],occQualities[k],
                                        errors[k]);
                                        
            aStats->WithError[occMismatches[k]]++;
        }
        aStats->RetrievedByCe+=lposCount;
        return lposCount;
    } else if (OutputType == SRA_REPORT_NONE) {
        for (k=0;k<posCount;k++) {
            aStats->WithError[occMismatches[k]]++;
        }
        return posCount;
    }
    return 0;
}


void OCCTranslateOccurrence(Translate * occTranslate, unsigned short * occAmbiguityMap,
                            unsigned long long ambPos,
                            uint8_t * ChromId, unsigned long long * unambPos) {
                            
    unsigned int approxIndex,approxValue;
    
    unsigned long long tp = ambPos;
    unsigned int correctPosition = ambPos;
    
    #ifndef DEBUG_2BWT_NO_TRANSLATION
    approxIndex = ambPos>>GRID_SAMPLING_FACTOR_2_POWER;
    approxValue = occAmbiguityMap[approxIndex];
    while (occTranslate[approxValue].startPos>ambPos) {
        approxValue--;
    }
    correctPosition-=occTranslate[approxValue].correction;
    tp=correctPosition;
    #endif
    
    (*ChromId) = occTranslate[approxValue].chrID;
    (*unambPos) = tp;
}

void OCCAddTextPositionToCache(SRAArguments * aArgs, 
                                    uint8_t type,
                                    unsigned long long ambPosition, 
                                    uint8_t strand, uint16_t matchLen,
                                    uint8_t occMismatch, int occQuality,
                                    const SRAError errorPos[]) { 
                                    
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAOutput * aOutput = aArgs->AlgnmtOutput;
    SRAOCCCollector * occCollector = aOutput->occCollector;
    SRAOccurrence occ;
    int i;
    
    occ.type             = type;
    occ.ambPosition      = ambPosition;
    occ.strand           = strand;
    occ.matchLen         = matchLen;
    occ.mismatchCount    = occMismatch;

    for (i=0;i<occMismatch&&i<MAX_NUM_OF_ERROR;i++) {
        occ.errors[i].type = errorPos[i].type;
        occ.errors[i].position = errorPos[i].position;
    }
        
    if (!SRAOCCAddOccurrenceToBuckets(occCollector,&occ)) {
        OCCFlushCache(aArgs);
        if (!SRAOCCAddOccurrenceToBuckets(occCollector,&occ)) {
            printf("[OCCAddTextPositionToCache] Unexpected error. OccCollector not clean. \n");
            exit(1);
        }
    }
}
