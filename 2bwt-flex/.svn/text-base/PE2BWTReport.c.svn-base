//
//    PE2BWTReport.c
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
   
   Date   : 19th October 2011
   Author : Edward MK Wu
   Change : New file created based on 2BWT-AlgnmtRpt.c
   
*/
/////////////////////////////////////////////////////

#include "PE2BWTReport.h"

// 2BWT-PE CTC#1
// -------------
// Uncomment this parameter PE_OCC_PLAIN_DUMP_MISMATCH_POSITION if you want
// 2BWT-Combine to print the exact position of each mismatch
// found in individual alignment in Plain output mode.
//#define PE_OCC_PLAIN_DUMP_MISMATCH_POSITION
// -------------
// 2BWT-PE CTC#2
// Uncomment the below to turn on logging message for
// overflowing dynamic programming result
#define PE_OCC_PRINT_OVERFLOW_MESSAGE
// --------------------------------------------------

unsigned int PEOCCWriteOutputHeader(PEArguments * peArgs) {
    SRAOutput * peOutput_OccHandler = peArgs->PEAlgnmtOutput_OccHandle;
    SRAIndex * aIndex = peArgs->sraArgsPos->AlgnmtIndex;
    PEInput * peInput = peArgs->PEAlgnmtInput;
    
    FILE * outFilePtr = peOutput_OccHandler->OutFilePtr;
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
    
    fprintf(outFilePtr,"%u\n",peOutput_OccHandler->OutFileFormat);
    
    fprintf(outFilePtr,"%u\n",hsp->numOfSeq);
    for (i=0;i<hsp->numOfSeq;i++) {
        OCCTranslateOccurrence(translate, ambiguityMap, (unsigned long long) hsp->seqOffset[i].endPos, &chromId, &tp);
        tp_32 = tp;
        fprintf(outFilePtr,"%u %u %s\n",hsp->annotation[i].gi,tp_32,hsp->annotation[i].text);
    }
    return 0;
}

void PEOCCFlushCachePlain(PEArguments * peArgs) {
    SRAQueryInfo * qInfo            = peArgs->sraArgsPos->QueryInfo;
    SRAOutput * peOutput_OccHandler = peArgs->PEAlgnmtOutput_OccHandle;
    SRAIndex * aIndex = peArgs->sraArgsPos->AlgnmtIndex;
    SRAOCCCollector * occCollector = peOutput_OccHandler->occCollector;
    FILE * outFilePtr = peOutput_OccHandler->OutFilePtr;
    HSP * hsp = aIndex->hsp;
    Translate * occTranslate = hsp->translate;
    unsigned short * occAmbiguityMap = hsp->ambiguityMap;
        
    DPArguments * dpArgs = peArgs->dpArguments;
    DPOCCCollector * dpOccCollector = dpArgs->dpOccCollector;
    DPOCCCollectorPointer * dpOccIterator = DPOCCPTCreate(dpOccCollector);
    
    DPOccurrence * dpOcc_1, * dpOcc_2;
    
    char strandStr[4] = "?+-";
    SRAOccurrence lastDelimEntry;
    int lastDelimValid;
    uint8_t chromId;
    unsigned long long tp;
    
    SRAOccurrence_ll * myNode = occCollector->occRootNode;
    if (myNode==NULL || myNode->size==0) return;
    unsigned int i = 0;

    #ifdef PE_OCC_PLAIN_DUMP_MISMATCH_POSITION
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
                
            /*} else if (myNode->bucket[i].type==SRAOCC_TYPE_PE_MATE_READ_NO_ALIGNMENT ||
                       myNode->bucket[i].type==SRAOCC_TYPE_PE_READ_MATE_NO_ALIGNMENT) {
                       
                OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                        myNode->bucket[i].ambPosition,
                                        &chromId,
                                        &tp);
                                        
                fprintf(outFilePtr,"%llu %u %llu %c %d",  lastDelimEntry.ambPosition,
                                                            chromId,
                                                            tp,
                                                            strandStr[myNode->bucket[i].strand],
                                                            myNode->bucket[i].mismatchCount);

                #ifdef PE_OCC_PLAIN_DUMP_MISMATCH_POSITION
                    int j;
                    for (j=0;j<myNode->bucket[i].mismatchCount;j++) {
                        if (myNode->bucket[i].errors[j].type==SRA_CHAR_ERROR_TYPE_MISMATCH) {
                            fprintf(outFilePtr," M %d",myNode->bucket[i].errors[j].position);
                        } else if (myNode->bucket[i].errors[j].type==SRA_CHAR_ERROR_TYPE_INSERT) {
                            fprintf(outFilePtr," I %d",myNode->bucket[i].errors[j].position);
                        } else if (myNode->bucket[i].errors[j].type==SRA_CHAR_ERROR_TYPE_DELETE) {
                            fprintf(outFilePtr," D %d",myNode->bucket[i].errors[j].position);
                        } else {
                            fprintf(outFilePtr," U %d",myNode->bucket[i].errors[j].position);
                        }
                    }
                #endif
                
                fprintf(outFilePtr,"\n");
            */} else if (myNode->bucket[i].type==SRAOCC_TYPE_PE_PAIR) {
            
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
                                                            
                #ifdef PE_OCC_PLAIN_DUMP_MISMATCH_POSITION
                    fprintf(outFilePtr," %d SRA",myNode->bucket[i].matchLen);
                    for (j=0;j<myNode->bucket[i].mismatchCount;j++) {
                        if (myNode->bucket[i].errors[j].type==SRA_CHAR_ERROR_TYPE_MISMATCH) {
                            fprintf(outFilePtr," M%d",myNode->bucket[i].errors[j].position);
                        } else if (myNode->bucket[i].errors[j].type==SRA_CHAR_ERROR_TYPE_INSERT) {
                            fprintf(outFilePtr," I%d",myNode->bucket[i].errors[j].position);
                        } else if (myNode->bucket[i].errors[j].type==SRA_CHAR_ERROR_TYPE_DELETE) {
                            fprintf(outFilePtr," D%d",myNode->bucket[i].errors[j].position);
                        } else {
                            fprintf(outFilePtr," U%d",myNode->bucket[i].errors[j].position);
                        }
                    }
                #endif
                
                fprintf(outFilePtr,"\n");
                #endif
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_PE_DP_BASE_READ) {
            
                dpOcc_1 = DPOCCPTRead(dpOccIterator);DPOCCPTNext(dpOccIterator);
                
                if (dpOcc_1==NULL) {
                    #ifdef PE_OCC_PRINT_OVERFLOW_MESSAGE
                        printf("[PEOccFlushCache/Read#%llu] DP Occurrence Collector does not have this occurrence.\n",qInfo->ReadId);
                    #endif
                } else {
                    
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
                                                                                                                                
                    #ifdef PE_OCC_PLAIN_DUMP_MISMATCH_POSITION
                        fprintf(outFilePtr," %d SRA",myNode->bucket[i].matchLen);
                        for (j=0;j<myNode->bucket[i].mismatchCount;j++) {
                            if (myNode->bucket[i].errors[j].type==SRA_CHAR_ERROR_TYPE_MISMATCH) {
                                fprintf(outFilePtr," M%d",myNode->bucket[i].errors[j].position);
                            } else if (myNode->bucket[i].errors[j].type==SRA_CHAR_ERROR_TYPE_INSERT) {
                                fprintf(outFilePtr," I%d",myNode->bucket[i].errors[j].position);
                            } else if (myNode->bucket[i].errors[j].type==SRA_CHAR_ERROR_TYPE_DELETE) {
                                fprintf(outFilePtr," D%d",myNode->bucket[i].errors[j].position);
                            } else {
                                fprintf(outFilePtr," U%d",myNode->bucket[i].errors[j].position);
                            }
                        }
                    #endif
                    
                    fprintf(outFilePtr,"\n");
                    #endif

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
                    
                    #ifdef PE_OCC_PLAIN_DUMP_MISMATCH_POSITION
                        fprintf(outFilePtr," %d DP",dpOcc_1->matchLen);
                        for (j=dpOcc_1->matchElemsCount-1;j>=0;j--) {
                            if (dpOcc_1->matchElems[j].type==SRA_CHAR_ERROR_TYPE_MISMATCH) {
                                fprintf(outFilePtr," %dM",dpOcc_1->matchElems[j].length);
                            } else if (dpOcc_1->matchElems[j].type==SRA_CHAR_ERROR_TYPE_INSERT) {
                                fprintf(outFilePtr," %dI",dpOcc_1->matchElems[j].length);
                            } else if (dpOcc_1->matchElems[j].type==SRA_CHAR_ERROR_TYPE_DELETE) {
                                fprintf(outFilePtr," %dD",dpOcc_1->matchElems[j].length);
                            } else {
                                fprintf(outFilePtr," %dU",dpOcc_1->matchElems[j].length);
                            }
                        }
                    #endif
                    
                    fprintf(outFilePtr,"\n");
                    #endif
                }
                
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_PE_DP_BASE_MATE) {
            
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
                                                                
                                                                
                    #ifdef PE_OCC_PLAIN_DUMP_MISMATCH_POSITION
                        fprintf(outFilePtr," %d DP",dpOcc_1->matchLen);
                        for (j=dpOcc_1->matchElemsCount-1;j>=0;j--) {
                            if (dpOcc_1->matchElems[j].type==SRA_CHAR_ERROR_TYPE_MISMATCH) {
                                fprintf(outFilePtr," %dM",dpOcc_1->matchElems[j].length);
                            } else if (dpOcc_1->matchElems[j].type==SRA_CHAR_ERROR_TYPE_INSERT) {
                                fprintf(outFilePtr," %dI",dpOcc_1->matchElems[j].length);
                            } else if (dpOcc_1->matchElems[j].type==SRA_CHAR_ERROR_TYPE_DELETE) {
                                fprintf(outFilePtr," %dD",dpOcc_1->matchElems[j].length);
                            } else {
                                fprintf(outFilePtr," %dU",dpOcc_1->matchElems[j].length);
                            }
                        }
                    #endif
                
                    fprintf(outFilePtr,"\n");
                    #endif
                    
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
                                                                                                                                
                    #ifdef PE_OCC_PLAIN_DUMP_MISMATCH_POSITION
                        fprintf(outFilePtr," %d SRA",myNode->bucket[i].matchLen);
                        for (j=0;j<myNode->bucket[i].mismatchCount;j++) {
                            if (myNode->bucket[i].errors[j].type==SRA_CHAR_ERROR_TYPE_MISMATCH) {
                                fprintf(outFilePtr," M%d",myNode->bucket[i].errors[j].position);
                            } else if (myNode->bucket[i].errors[j].type==SRA_CHAR_ERROR_TYPE_INSERT) {
                                fprintf(outFilePtr," I%d",myNode->bucket[i].errors[j].position);
                            } else if (myNode->bucket[i].errors[j].type==SRA_CHAR_ERROR_TYPE_DELETE) {
                                fprintf(outFilePtr," D%d",myNode->bucket[i].errors[j].position);
                            } else {
                                fprintf(outFilePtr," U%d",myNode->bucket[i].errors[j].position);
                            }
                        }
                    #endif
                    
                    fprintf(outFilePtr,"\n");
                    #endif
                }
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_PE_DP_SEED_OCC) {
            
                dpOcc_1 = DPOCCPTRead(dpOccIterator);DPOCCPTNext(dpOccIterator);
                dpOcc_2 = DPOCCPTRead(dpOccIterator);DPOCCPTNext(dpOccIterator);
                
                if (dpOcc_1==NULL || dpOcc_2==NULL) {
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

                                                                
                    #ifdef PE_OCC_PLAIN_DUMP_MISMATCH_POSITION
                        fprintf(outFilePtr," %d DP",dpOcc_1->matchLen);
                        for (j=dpOcc_1->matchElemsCount-1;j>=0;j--) {
                            if (dpOcc_1->matchElems[j].type==SRA_CHAR_ERROR_TYPE_MISMATCH) {
                                fprintf(outFilePtr," %dM",dpOcc_1->matchElems[j].length);
                            } else if (dpOcc_1->matchElems[j].type==SRA_CHAR_ERROR_TYPE_INSERT) {
                                fprintf(outFilePtr," %dI",dpOcc_1->matchElems[j].length);
                            } else if (dpOcc_1->matchElems[j].type==SRA_CHAR_ERROR_TYPE_DELETE) {
                                fprintf(outFilePtr," %dD",dpOcc_1->matchElems[j].length);
                            } else {
                                fprintf(outFilePtr," %dU",dpOcc_1->matchElems[j].length);
                            }
                        }
                    #endif
                    
                    fprintf(outFilePtr,"\n");
                    #endif
                                        
                    OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                            dpOcc_2->ambPosition,
                                            &chromId,
                                            &tp);
                                            
                    #ifndef DEBUG_2BWT_NO_WRITE_FILE
                    fprintf(outFilePtr,"%llu %u %llu %c %d",  lastDelimEntry.ambPosition,
                                                                chromId,
                                                                tp,
                                                                strandStr[dpOcc_2->strand],
                                                                dpOcc_2->mismatchCount);

                                                                
                    #ifdef PE_OCC_PLAIN_DUMP_MISMATCH_POSITION
                        fprintf(outFilePtr," %d DP",dpOcc_2->matchLen);
                        for (j=dpOcc_2->matchElemsCount-1;j>=0;j--) {
                            if (dpOcc_2->matchElems[j].type==SRA_CHAR_ERROR_TYPE_MISMATCH) {
                                fprintf(outFilePtr," %dM",dpOcc_2->matchElems[j].length);
                            } else if (dpOcc_2->matchElems[j].type==SRA_CHAR_ERROR_TYPE_INSERT) {
                                fprintf(outFilePtr," %dI",dpOcc_2->matchElems[j].length);
                            } else if (dpOcc_2->matchElems[j].type==SRA_CHAR_ERROR_TYPE_DELETE) {
                                fprintf(outFilePtr," %dD",dpOcc_2->matchElems[j].length);
                            } else {
                                fprintf(outFilePtr," %dU",dpOcc_2->matchElems[j].length);
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

void PEOCCFlushCacheSAM(PEArguments * peArgs) {
    SRAQueryInfo * qInfo            = peArgs->sraArgsPos->QueryInfo;
    SRAQueryInfo * qInfo_mate       = peArgs->sraArgsPos_mate->QueryInfo;
    SRAIndex * aIndex               = peArgs->sraArgsPos->AlgnmtIndex;
    
    SRAOutput * peOutput_OccHandler = peArgs->PEAlgnmtOutput_OccHandle;
    SRAOCCCollector * occCollector = peOutput_OccHandler->occCollector;
    samfile_t * samFilePtr = peOutput_OccHandler->SAMOutFilePtr;
    
    DPArguments * dpArgs = peArgs->dpArguments;
    DPOCCCollector * dpOccCollector = dpArgs->dpOccCollector;
    DPOCCCollectorPointer * dpOccIterator = DPOCCPTCreate(dpOccCollector);
    
    DPOccurrence * dpOcc_1, * dpOcc_2;
    
    HSP * hsp = aIndex->hsp;
    Translate * occTranslate = hsp->translate;
    unsigned short * occAmbiguityMap = hsp->ambiguityMap;

    bam1_t samAlgnmt;
    bam1_t samAlgnmt_mate;
    
    unsigned long long i=0, k=0;
    uint8_t chromId_1,chromId_2;
    unsigned long long tp_1,tp_2;
    
    SRAOccurrence_ll * myNode = occCollector->occRootNode;
    if (myNode==NULL || myNode->size==0) return;
    
    if (qInfo->ReportingReadCode==NULL) { qInfo->ReportingReadCode = qInfo->ReadCode; }
    
    unsigned int readNumOfOcc = SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_PE_PAIR)/2 + 
                                SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_PE_READ_READ_NO_ALIGNMENT) + 
                                SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_PE_READ_MATE_NO_ALIGNMENT);
    unsigned int readNumOfDpOcc = SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_PE_DP_BASE_MATE) + 
                                  SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_PE_DP_SEED_OCC);
                                  
    unsigned int readNumOfOcc_mate = SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_PE_PAIR)/2 + 
                                SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_PE_MATE_READ_NO_ALIGNMENT) + 
                                SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_PE_MATE_MATE_NO_ALIGNMENT);
    unsigned int readNumOfDpOcc_mate = SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_PE_DP_BASE_READ) + 
                                       SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_PE_DP_SEED_OCC);
                                
    unsigned long long samAuxMultiOccBlockSize = SAM_MDATA_SIZE_PER_READ + 
                                                (readNumOfOcc) * (SAM_MDATA_SIZE_PER_OCC) + 
                                                (readNumOfDpOcc) * (SAM_MDATA_SIZE_PER_DPOCC);
    unsigned long long samAuxMultiOccBlockSize_mate = SAM_MDATA_SIZE_PER_READ + 
                                                (readNumOfOcc_mate) * (SAM_MDATA_SIZE_PER_OCC) + 
                                                (readNumOfDpOcc_mate) * (SAM_MDATA_SIZE_PER_DPOCC);
    
    int samPrimaryOccAvailable = 0;
    int samPrimaryOccAvailable_mate = 0;
    int samXAZtagCreated = 0;
    int samXAZtagCreated_mate = 0;
    
    if (samAuxMultiOccBlockSize+samAuxMultiOccBlockSize_mate+MAX_FIELD_LEN >peOutput_OccHandler->SAMAuxDataBlockSize) {
    
        peOutput_OccHandler->SAMAuxDataBlock = (uint8_t *)realloc(peOutput_OccHandler->SAMAuxDataBlock,
                                               sizeof(uint8_t)*(samAuxMultiOccBlockSize+samAuxMultiOccBlockSize_mate + MAX_FIELD_LEN ));
        peOutput_OccHandler->SAMAuxDataBlockSize = samAuxMultiOccBlockSize+samAuxMultiOccBlockSize_mate+ MAX_FIELD_LEN;
        
        printf("[PEOCCFlushCacheSAM] Re-sizing the SAM buffer to %llu bytes.\n",
                (unsigned long long) sizeof(uint8_t)*(samAuxMultiOccBlockSize+samAuxMultiOccBlockSize_mate+ MAX_FIELD_LEN));
    }
    
    samAlgnmt.data = peOutput_OccHandler->SAMAuxDataBlock;
    samAlgnmt.m_data = samAuxMultiOccBlockSize;
    samAlgnmt_mate.data = (peOutput_OccHandler->SAMAuxDataBlock)+samAuxMultiOccBlockSize;
    samAlgnmt_mate.m_data = samAuxMultiOccBlockSize_mate;
    
    //Assuming the answer being handled by OCCFlushCacheSAMAPI are all answers of qInfo
    // Each occurrences reported by SAM_API will share the following aux data.
    //----------------------------------------------------------------------------------------------------------------------------
        samAlgnmt.core.bin = bam_reg2bin(0,0);                     //SAM: bin calculated by bam_reg2bin()
        samAlgnmt.core.l_qseq = qInfo->ReadLength;                 //SAM: length of the query sequence (read)
        samAlgnmt.core.l_qname = SAMSeqReadNameLength(qInfo->ReadName,MAX_SEQ_NAME_LENGTH)+1;        //SAM: length of the query name
        samAlgnmt.core.qual = SAM_MAPQ_UNAVAILABLE;                //MAPQ: Unavailable by default
        
        samAlgnmt.l_aux = 0;
        samAlgnmt.data_len = 0;
        SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),qInfo->ReadName,samAlgnmt.core.l_qname);    //Name
        
        int sharedDataLen = samAlgnmt.data_len;
        
    //----------------------------------------------------------------------------------------------------------------------------                      
        samAlgnmt_mate.core.bin = bam_reg2bin(0,0);                     //SAM: bin calculated by bam_reg2bin()
        samAlgnmt_mate.core.l_qseq = qInfo_mate->ReadLength;            //SAM: length of the query sequence (read)
        samAlgnmt_mate.core.l_qname = SAMSeqReadNameLength(qInfo_mate->ReadName,MAX_SEQ_NAME_LENGTH)+1;        //SAM: length of the query name
        samAlgnmt_mate.core.qual = SAM_MAPQ_UNAVAILABLE;                //MAPQ: Unavailable by default
        
        samAlgnmt_mate.l_aux = 0;
        samAlgnmt_mate.data_len = 0;
        SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),qInfo_mate->ReadName,samAlgnmt_mate.core.l_qname);    //Name
        int sharedDataLen_mate = samAlgnmt_mate.data_len;

    //----------------------------------------------------------------------------------------------------------------------------                      

    int validPair = 0;
    int invalidPair = 0;
    SRAOccurrence * occIdx[2];
    long long insertion;
    while (myNode != NULL) {
        for (i=0;i<myNode->size;i++) {
            if (myNode->bucket[i].type==SRAOCC_TYPE_DELIMITOR_READ) {
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_PE_READ_READ_NO_ALIGNMENT) {
                            
                OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                        myNode->bucket[i].ambPosition,
                                        &chromId_2,
                                        &tp_2);
                                        
                if (samPrimaryOccAvailable) {
                    /*if (!samXAZtagCreated) {
                        SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),"XAZ",3);
                        samXAZtagCreated=1;
                    }
                    //--------------------------
                        unsigned int auxStart = samAlgnmt.data_len;
                        SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),
                                                hsp->annotation[chromId_2-1].text,
                                                SAMSeqReadNameLength(hsp->annotation[chromId_2-1].text,MAX_SEQ_NAME_LENGTH));
                        if (myNode->bucket[i].strand==QUERY_NEG_STRAND) {
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",-",2);
                        } else if (myNode->bucket[i].strand==QUERY_POS_STRAND) {
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",+",2);
                        } else {
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",?",2);
                        }
                        SAMIUint8ConcatUint32AsString(samAlgnmt.data,&(samAlgnmt.data_len),tp_2);
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
                    //--------------------------*/
                } else {
                    
                    unsigned int samFlag = 0;
                    samFlag |= SAM_FLAG_IS_PAIR_READ;
                    samFlag |= SAM_FLAG_READ_UNMAPPED;
                    samFlag |= SAM_FLAG_FIRST_IN_PAIR;
                
                    if (myNode->bucket[i].strand==QUERY_NEG_STRAND) {
                        samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
                    }
                    
                    SAMPackSequenceAndCigar(&samAlgnmt,sharedDataLen,qInfo,&(myNode->bucket[i]));
                    samAlgnmt.l_aux = 0;
                    // Read Group at end of each line
                    if (peOutput_OccHandler->ReadGroup!= NULL && strlen( peOutput_OccHandler->ReadGroup ) >0){
                        unsigned int auxStart = samAlgnmt.data_len; 
                        bam_aux_append ( &samAlgnmt, "RG", 'Z', strlen (peOutput_OccHandler->ReadGroup  ) + 1, ( uint8_t * ) peOutput_OccHandler->ReadGroup  );
                        samAlgnmt.l_aux += samAlgnmt.data_len - auxStart;
                    }
                    // ---------------------->  |
                                                samAlgnmt.core.tid = -1;                                   //SAM: chromosome ID, defined by bam_header_t
                                                samAlgnmt.core.pos = -1;                                   //SAM: 0-based leftmost coordinate4
                                                samAlgnmt.core.qual = SAM_MAPQ_UNALIGNED;                  //SAM: mapping quality
                                                samAlgnmt.core.flag = samFlag;                             //SAM: bitwise flag
                                                
                                                samAlgnmt.core.mtid = chromId_2-1;
                                                samAlgnmt.core.mpos = tp_2-1;
                                                samAlgnmt.core.isize = 0;
                                                
                                                samPrimaryOccAvailable = 1;
                    // ---------------------->  |
                }
                
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_PE_MATE_READ_NO_ALIGNMENT) {
            
                OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                        myNode->bucket[i].ambPosition,
                                        &chromId_2,
                                        &tp_2);
                                        
                if (samPrimaryOccAvailable_mate) {
                    /*if (!samXAZtagCreated_mate) {
                        SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),"XAZ",3);
                        samXAZtagCreated_mate=1;
                    }
                    //--------------------------
                        unsigned int auxStart_mate = samAlgnmt_mate.data_len;
                        SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),
                                                hsp->annotation[chromId_2-1].text,
                                                SAMSeqReadNameLength(hsp->annotation[chromId_2-1].text,MAX_SEQ_NAME_LENGTH));
                        if (myNode->bucket[i].strand==QUERY_NEG_STRAND) {
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",-",2);
                        } else if (myNode->bucket[i].strand==QUERY_POS_STRAND) {
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",+",2);
                        } else {
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",?",2);
                        }
                        SAMIUint8ConcatUint32AsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),tp_2);
                        SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",",1);
                        SAMIUint8ConcatCigarAsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),
                                                    qInfo->ReadLength,
                                                    myNode->bucket[i].strand,
                                                    myNode->bucket[i].mismatchCount,
                                                    myNode->bucket[i].errors, SAM_ERROR_RANDOM);
                        SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",",1);
                        SAMIUint8ConcatUint32AsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),myNode->bucket[i].mismatchCount);
                        SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),";",1);
                            
                        samAlgnmt_mate.l_aux += samAlgnmt_mate.data_len - auxStart_mate;
                    //--------------------------*/
                } else {
                
                    unsigned int samFlag = 0;
                    samFlag |= SAM_FLAG_IS_PAIR_READ;
                    samFlag |= SAM_FLAG_READ_UNMAPPED;
                    samFlag |= SAM_FLAG_SECOND_IN_PAIR;
                
                    if (myNode->bucket[i].strand==QUERY_NEG_STRAND) {
                        samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
                    }
                    
                    SAMPackSequenceAndCigar(&samAlgnmt_mate,sharedDataLen_mate,qInfo_mate,&(myNode->bucket[i]));
                    samAlgnmt_mate.l_aux = 0;
                    // Read Group at end of each line
                    if (peOutput_OccHandler->ReadGroup!= NULL && strlen( peOutput_OccHandler->ReadGroup ) >0){
                        unsigned int auxStart = samAlgnmt_mate.data_len; 
                        bam_aux_append ( &samAlgnmt_mate, "RG", 'Z', strlen (peOutput_OccHandler->ReadGroup  ) + 1, ( uint8_t * ) peOutput_OccHandler->ReadGroup  );
                        samAlgnmt_mate.l_aux += samAlgnmt_mate.data_len - auxStart;
                    }
                    // ---------------------->  |
                                                samAlgnmt_mate.core.tid = -1;                                   //SAM: chromosome ID, defined by bam_header_t
                                                samAlgnmt_mate.core.pos = -1;                                   //SAM: 0-based leftmost coordinate4
                                                samAlgnmt_mate.core.qual = SAM_MAPQ_UNALIGNED;                  //SAM: mapping quality
                                                samAlgnmt_mate.core.flag = samFlag;                             //SAM: bitwise flag
                                                
                                                samAlgnmt_mate.core.mtid = chromId_2-1;
                                                samAlgnmt_mate.core.mpos = tp_2-1;
                                                samAlgnmt_mate.core.isize = 0;
                                                
                                                samPrimaryOccAvailable_mate = 1;
                    // ---------------------->  |
                }
                
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_PE_READ_MATE_NO_ALIGNMENT) {
            
                OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                        myNode->bucket[i].ambPosition,
                                        &chromId_1,
                                        &tp_1);
                if (samPrimaryOccAvailable) {
                    if (!samXAZtagCreated) {
                        SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),"XAZ",3);
                        samXAZtagCreated=1;
                    }
                    //--------------------------
                        unsigned int auxStart = samAlgnmt.data_len;
                        SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),
                                                hsp->annotation[chromId_1-1].text,
                                                SAMSeqReadNameLength(hsp->annotation[chromId_1-1].text,MAX_SEQ_NAME_LENGTH));
                        if (myNode->bucket[i].strand==QUERY_NEG_STRAND) {
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",-",2);
                        } else if (myNode->bucket[i].strand==QUERY_POS_STRAND) {
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",+",2);
                        } else {
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",?",2);
                        }
                        SAMIUint8ConcatUint32AsString(samAlgnmt.data,&(samAlgnmt.data_len),tp_1);
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
                    //--------------------------
                } else {
                
                    unsigned int samFlag = 0;
                    samFlag |= SAM_FLAG_IS_PAIR_READ;
                    samFlag |= SAM_FLAG_MATE_UNMAPPED;
                    samFlag |= SAM_FLAG_FIRST_IN_PAIR;
                    
                    if (myNode->bucket[i].strand==QUERY_NEG_STRAND) {
                        samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
                    }
                
                    SAMPackSequenceAndCigar(&samAlgnmt,sharedDataLen,qInfo,&(myNode->bucket[i]));
                    samAlgnmt.l_aux = 0;
                    // Read Group at end of each line
                    if (peOutput_OccHandler->ReadGroup!= NULL && strlen( peOutput_OccHandler->ReadGroup ) >0){
                        unsigned int auxStart = samAlgnmt.data_len; 
                        bam_aux_append ( &samAlgnmt, "RG", 'Z', strlen (peOutput_OccHandler->ReadGroup  ) + 1, ( uint8_t * ) peOutput_OccHandler->ReadGroup  );
                        samAlgnmt.l_aux += samAlgnmt.data_len - auxStart;
                    }
                    // ---------------------->  |
                                                samAlgnmt.core.tid = chromId_1-1;                          //SAM: chromosome ID, defined by bam_header_t
                                                samAlgnmt.core.pos = tp_1-1;                               //SAM: 0-based leftmost coordinate
                                                samAlgnmt.core.qual = SAM_MAPQ_UNALIGNED;                  //SAM: mapping quality
                                                samAlgnmt.core.flag = samFlag;                             //SAM: bitwise flag
                                                
                                                samAlgnmt.core.mtid = -1;
                                                samAlgnmt.core.mpos = -1;
                                                samAlgnmt.core.isize = 0;
                                                
                                                samPrimaryOccAvailable = 1;
                    // ---------------------->  |
                }
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_PE_MATE_MATE_NO_ALIGNMENT) {
            
                OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                        myNode->bucket[i].ambPosition,
                                        &chromId_1,
                                        &tp_1);
                                        
                if (samPrimaryOccAvailable_mate) {
                    //--------------------------
                        unsigned int auxStart_mate = samAlgnmt_mate.data_len;
                        
                        if (!samXAZtagCreated_mate) {
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),"XAZ",3);
                            samXAZtagCreated_mate=1;
                        }
                        SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),
                                                hsp->annotation[chromId_1-1].text,
                                                SAMSeqReadNameLength(hsp->annotation[chromId_1-1].text,MAX_SEQ_NAME_LENGTH));
                        if (myNode->bucket[i].strand==QUERY_NEG_STRAND) {
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",-",2);
                        } else if (myNode->bucket[i].strand==QUERY_POS_STRAND) {
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",+",2);
                        } else {
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",?",2);
                        }
                        SAMIUint8ConcatUint32AsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),tp_1);
                        SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",",1);
                        SAMIUint8ConcatCigarAsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),
                                                    qInfo->ReadLength,
                                                    myNode->bucket[i].strand,
                                                    myNode->bucket[i].mismatchCount,
                                                    myNode->bucket[i].errors, SAM_ERROR_RANDOM);
                        SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",",1);
                        SAMIUint8ConcatUint32AsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),myNode->bucket[i].mismatchCount);
                        SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),";",1);
                            
                        samAlgnmt_mate.l_aux += samAlgnmt_mate.data_len - auxStart_mate;
                    //--------------------------
                } else {
                    
                    unsigned int samFlag = 0;
                    samFlag |= SAM_FLAG_IS_PAIR_READ;
                    samFlag |= SAM_FLAG_MATE_UNMAPPED;
                    samFlag |= SAM_FLAG_SECOND_IN_PAIR;
                    
                    if (myNode->bucket[i].strand==QUERY_NEG_STRAND) {
                        samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
                    }
                    
                    SAMPackSequenceAndCigar(&samAlgnmt_mate,sharedDataLen_mate,qInfo_mate,&(myNode->bucket[i]));
                    samAlgnmt_mate.l_aux = 0;
                    // Read Group at end of each line
                    if (peOutput_OccHandler->ReadGroup!= NULL && strlen( peOutput_OccHandler->ReadGroup ) >0){
                        unsigned int auxStart = samAlgnmt_mate.data_len; 
                        bam_aux_append ( &samAlgnmt_mate, "RG", 'Z', strlen (peOutput_OccHandler->ReadGroup  ) + 1, ( uint8_t * ) peOutput_OccHandler->ReadGroup  );
                        samAlgnmt_mate.l_aux += samAlgnmt_mate.data_len - auxStart;
                    }
                    // ---------------------->  |
                                                samAlgnmt_mate.core.tid = chromId_1-1;                          //SAM: chromosome ID, defined by bam_header_t
                                                samAlgnmt_mate.core.pos = tp_1-1;                               //SAM: 0-based leftmost coordinate
                                                samAlgnmt_mate.core.qual = SAM_MAPQ_UNALIGNED;                  //SAM: mapping quality
                                                samAlgnmt_mate.core.flag = samFlag;                             //SAM: bitwise flag
                                                
                                                samAlgnmt_mate.core.mtid = -1;
                                                samAlgnmt_mate.core.mpos = -1;
                                                samAlgnmt_mate.core.isize = 0;
                                                
                                                samPrimaryOccAvailable_mate = 1;
                    // ---------------------->  |
                }
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_PE_READ_BOTH_NO_ALIGNMENT) {
                unsigned int samFlag = 0;
                samFlag |= SAM_FLAG_IS_PAIR_READ;
                samFlag |= SAM_FLAG_READ_UNMAPPED;
                samFlag |= SAM_FLAG_MATE_UNMAPPED;
                samFlag |= SAM_FLAG_FIRST_IN_PAIR;
            
                SAMPackSequenceAndCigar(&samAlgnmt,sharedDataLen,qInfo,&(myNode->bucket[i]));
                samAlgnmt.l_aux = 0;
                // Read Group at end of each line
                if (peOutput_OccHandler->ReadGroup!= NULL && strlen( peOutput_OccHandler->ReadGroup ) >0){
                    unsigned int auxStart = samAlgnmt.data_len; 
                    bam_aux_append ( &samAlgnmt, "RG", 'Z', strlen (peOutput_OccHandler->ReadGroup  ) + 1, ( uint8_t * ) peOutput_OccHandler->ReadGroup  );
                    samAlgnmt.l_aux += samAlgnmt.data_len - auxStart;
                }
                // ---------------------->  |
                                            samAlgnmt.core.tid = -1;                                   //SAM: chromosome ID, defined by bam_header_t
                                            samAlgnmt.core.pos = -1;                                   //SAM: 0-based leftmost coordinate4
                                            samAlgnmt.core.qual = SAM_MAPQ_UNALIGNED;                  //SAM: mapping quality
                                            samAlgnmt.core.flag = samFlag;                             //SAM: bitwise flag
                                            
                                            samAlgnmt.core.mtid = -1;
                                            samAlgnmt.core.mpos = -1;
                                            samAlgnmt.core.isize = 0;
                                            
                                            #ifndef DEBUG_2BWT_NO_WRITE_FILE
                                            samwrite(samFilePtr,&samAlgnmt);
                                            #endif
                // ---------------------->  |
                
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_PE_MATE_BOTH_NO_ALIGNMENT) {
                unsigned int samFlag = 0;
                samFlag |= SAM_FLAG_IS_PAIR_READ;
                samFlag |= SAM_FLAG_READ_UNMAPPED;
                samFlag |= SAM_FLAG_MATE_UNMAPPED;
                samFlag |= SAM_FLAG_SECOND_IN_PAIR;
                
                SAMPackSequenceAndCigar(&samAlgnmt_mate,sharedDataLen_mate,qInfo_mate,&(myNode->bucket[i]));
                samAlgnmt_mate.l_aux = 0;
                // Read Group at end of each line
                if (peOutput_OccHandler->ReadGroup!= NULL && strlen( peOutput_OccHandler->ReadGroup ) >0){
                    unsigned int auxStart = samAlgnmt_mate.data_len; 
                    bam_aux_append ( &samAlgnmt_mate, "RG", 'Z', strlen (peOutput_OccHandler->ReadGroup  ) + 1, ( uint8_t * ) peOutput_OccHandler->ReadGroup  );
                    samAlgnmt_mate.l_aux += samAlgnmt_mate.data_len - auxStart;
                }
                // ---------------------->  |
                                            samAlgnmt_mate.core.tid = -1;                                   //SAM: chromosome ID, defined by bam_header_t
                                            samAlgnmt_mate.core.pos = -1;                                   //SAM: 0-based leftmost coordinate4
                                            samAlgnmt_mate.core.qual = SAM_MAPQ_UNALIGNED;                  //SAM: mapping quality
                                            samAlgnmt_mate.core.flag = samFlag;                             //SAM: bitwise flag
                                            
                                            samAlgnmt_mate.core.mtid = -1;
                                            samAlgnmt_mate.core.mpos = -1;
                                            samAlgnmt_mate.core.isize = 0;
                                            
                                            #ifndef DEBUG_2BWT_NO_WRITE_FILE
                                            samwrite(samFilePtr,&samAlgnmt_mate);
                                            #endif
                // ---------------------->  |
                
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_PE_PAIR) {
                occIdx[validPair++] = &(myNode->bucket[i]);
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_PE_BAD_PAIR) {
                occIdx[invalidPair++] = &(myNode->bucket[i]);
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_PE_DP_BASE_READ) {
            
                // 
                // The occurrence from the READ is found by PE-SRA and the occurrence was stored
                // in myNode->bucket[i] as SRAOccurrence. The occurrence will then be stored into
                // samAlgnmt.
                //
                // The occurrence from the MATE is found by PE-DP and the occurrence was stored in
                // dpOcc_1. The occurrence will then be stored into samAlgnmt_mate.
                //
                
                dpOcc_1 = DPOCCPTRead(dpOccIterator);DPOCCPTNext(dpOccIterator);
                
                if (dpOcc_1==NULL) {
                    #ifdef PE_OCC_PRINT_OVERFLOW_MESSAGE
                        printf("[PEOccFlushCache/Read#%llu] DP Occurrence Collector does not have this occurrence.\n",qInfo->ReadId);
                    #endif
                } else {
                    
                    OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                            myNode->bucket[i].ambPosition,
                                            &chromId_1,
                                            &tp_1);
                                            
                    OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                            dpOcc_1->ambPosition,
                                            &chromId_2,
                                            &tp_2);
                                            
                    if (samPrimaryOccAvailable || samPrimaryOccAvailable_mate) {
                        
                        //--------------------------
                            unsigned int auxStart = samAlgnmt.data_len;
                            if (!samXAZtagCreated) {
                                SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),"XAZ",3);
                                samXAZtagCreated=1;
                            }
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),
                                                    hsp->annotation[chromId_1-1].text,
                                                    SAMSeqReadNameLength(hsp->annotation[chromId_1-1].text,MAX_SEQ_NAME_LENGTH));
                            if (myNode->bucket[i].strand==QUERY_NEG_STRAND) {
                                SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",-",2);
                            } else if (myNode->bucket[i].strand==QUERY_POS_STRAND) {
                                SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",+",2);
                            } else {
                                SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",?",2);
                            }
                            SAMIUint8ConcatUint32AsString(samAlgnmt.data,&(samAlgnmt.data_len),tp_1);
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
                        //--------------------------
                            unsigned int auxStart_mate = samAlgnmt_mate.data_len;
                            if (!samXAZtagCreated_mate) {
                                SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),"XAZ",3);
                                samXAZtagCreated_mate=1;
                            }
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),
                                                    hsp->annotation[chromId_2-1].text,
                                                    SAMSeqReadNameLength(hsp->annotation[chromId_2-1].text,MAX_SEQ_NAME_LENGTH));
                            if (dpOcc_1->strand==QUERY_NEG_STRAND) {
                                SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",-",2);
                            } else if (dpOcc_1->strand==QUERY_POS_STRAND) {
                                SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",+",2);
                            } else {
                                SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",?",2);
                            }
                            SAMIUint8ConcatUint32AsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),tp_2);
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",",1);
                            SAMDPIUint8ConcatCigarAsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),
                                                        qInfo->ReadLength,
                                                        dpOcc_1->strand,
                                                        dpOcc_1->matchElemsCount,
                                                        dpOcc_1->matchElems);
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",",1);
                            SAMIUint8ConcatUint32AsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),dpOcc_1->mismatchCount);
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),";",1);
                                
                            samAlgnmt_mate.l_aux += samAlgnmt_mate.data_len - auxStart_mate;
                        //--------------------------
                    } else {
                        if (chromId_1==chromId_2) {
                            if (myNode->bucket[i].ambPosition > dpOcc_1->ambPosition) {
                                insertion = - (myNode->bucket[i].ambPosition + myNode->bucket[i].matchLen - dpOcc_1->ambPosition);
                            } else {
                                insertion = dpOcc_1->ambPosition + dpOcc_1->matchLen - myNode->bucket[i].ambPosition;
                            }
                        } else {
                            insertion = 0;
                        }
                                                
                        unsigned int samFlag = 0;
                        
                        samFlag |= SAM_FLAG_IS_PAIR_READ;
                        samFlag |= SAM_FLAG_PROPER_PAIR;
                        samFlag |= SAM_FLAG_FIRST_IN_PAIR;
                        if (myNode->bucket[i].strand==QUERY_NEG_STRAND) {
                            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
                        }
                        if (dpOcc_1->strand==QUERY_NEG_STRAND) {
                            samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
                        }
                        
                        SAMPackSequenceAndCigar(&samAlgnmt,sharedDataLen,qInfo,&(myNode->bucket[i]));
                        samAlgnmt.l_aux = 0;
                        // Read Group at end of each line
                        if (peOutput_OccHandler->ReadGroup!= NULL && strlen( peOutput_OccHandler->ReadGroup ) >0){
                            unsigned int auxStart = samAlgnmt.data_len; 
                            bam_aux_append ( &samAlgnmt, "RG", 'Z', strlen (peOutput_OccHandler->ReadGroup  ) + 1, ( uint8_t * ) peOutput_OccHandler->ReadGroup  );
                            samAlgnmt.l_aux += samAlgnmt.data_len - auxStart;
                        }
                        // ---------------------->  |
                                                    samAlgnmt.core.tid = chromId_1-1;                           //SAM: chromosome ID, defined by bam_header_t
                                                    samAlgnmt.core.pos = tp_1-1;                                 //SAM: 0-based leftmost coordinate
                                                    samAlgnmt.core.flag = samFlag;                             //SAM: bitwise flag
                                                    
                                                    samAlgnmt.core.mtid = chromId_2-1;
                                                    samAlgnmt.core.mpos = tp_2-1;
                                                    samAlgnmt.core.isize = insertion;
                        // ---------------------->  |
                        
                        samFlag = 0;
                        
                        samFlag |= SAM_FLAG_IS_PAIR_READ;
                        samFlag |= SAM_FLAG_PROPER_PAIR;
                        samFlag |= SAM_FLAG_SECOND_IN_PAIR;
                        if (myNode->bucket[i].strand==QUERY_NEG_STRAND) {
                            samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
                        }
                        if (dpOcc_1->strand==QUERY_NEG_STRAND) {
                            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
                        }
                        
                        SAMDPPackSequenceAndCigar(&samAlgnmt_mate,sharedDataLen_mate,qInfo_mate,dpOcc_1);
                        samAlgnmt_mate.l_aux = 0;
                        //SAMPackSequenceAndCigar(&samAlgnmt_mate,sharedDataLen_mate,qInfo_mate,&(myNode->bucket[i]));
                        // Read Group at end of each line
                        if (peOutput_OccHandler->ReadGroup!= NULL && strlen( peOutput_OccHandler->ReadGroup ) >0){
                            unsigned int auxStart = samAlgnmt_mate.data_len; 
                            bam_aux_append ( &samAlgnmt_mate, "RG", 'Z', strlen (peOutput_OccHandler->ReadGroup  ) + 1, ( uint8_t * ) peOutput_OccHandler->ReadGroup  );
                            samAlgnmt_mate.l_aux += samAlgnmt_mate.data_len - auxStart;
                        }
                        // ---------------------->  |
                                                    samAlgnmt_mate.core.tid = chromId_2-1;                           //SAM: chromosome ID, defined by bam_header_t
                                                    samAlgnmt_mate.core.pos = tp_2-1;                                 //SAM: 0-based leftmost coordinate
                                                    samAlgnmt_mate.core.flag = samFlag;                             //SAM: bitwise flag
                                                    
                                                    samAlgnmt_mate.core.mtid = chromId_1-1;
                                                    samAlgnmt_mate.core.mpos = tp_1-1;
                                                    samAlgnmt_mate.core.isize = - insertion;
                        // ---------------------->  |
                        
                        samPrimaryOccAvailable = 1;
                        samPrimaryOccAvailable_mate = 1;
                    }
                    
                }

            } else if (myNode->bucket[i].type==SRAOCC_TYPE_PE_DP_BASE_MATE) {
            
                // 
                // The occurrence from the READ is found by PE-SRA and the occurrence was stored
                // in myNode->bucket[i] as SRAOccurrence. The occurrence will then be stored into
                // samAlgnmt.
                //
                // The occurrence from the MATE is found by PE-DP and the occurrence was stored in
                // dpOcc_1. The occurrence will then be stored into samAlgnmt_mate.
                //
                
                dpOcc_1 = DPOCCPTRead(dpOccIterator);DPOCCPTNext(dpOccIterator);
                
                if (dpOcc_1==NULL) {
                    #ifdef PE_OCC_PRINT_OVERFLOW_MESSAGE
                        printf("[PEOccFlushCache/Read#%llu] DP Occurrence Collector does not have this occurrence.\n",qInfo->ReadId);
                    #endif
                } else {
                
                    OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                            myNode->bucket[i].ambPosition,
                                            &chromId_2,
                                            &tp_2);
                                            
                    OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                            dpOcc_1->ambPosition,
                                            &chromId_1,
                                            &tp_1);

                    if (samPrimaryOccAvailable || samPrimaryOccAvailable_mate) {
                        //--------------------------
                            unsigned int auxStart = samAlgnmt.data_len;
                            if (!samXAZtagCreated) {
                                SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),"XAZ",3);
                                samXAZtagCreated=1;
                            }
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),
                                                    hsp->annotation[chromId_1-1].text,
                                                    SAMSeqReadNameLength(hsp->annotation[chromId_1-1].text,MAX_SEQ_NAME_LENGTH));
                            if (dpOcc_1->strand==QUERY_NEG_STRAND) {
                                SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",-",2);
                            } else if (dpOcc_1->strand==QUERY_POS_STRAND) {
                                SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",+",2);
                            } else {
                                SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",?",2);
                            }
                            SAMIUint8ConcatUint32AsString(samAlgnmt.data,&(samAlgnmt.data_len),tp_1);
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
                            unsigned int auxStart_mate = samAlgnmt_mate.data_len;
                            if (!samXAZtagCreated_mate) {
                                SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),"XAZ",3);
                                samXAZtagCreated_mate=1;
                            }
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),
                                                    hsp->annotation[chromId_2-1].text,
                                                    SAMSeqReadNameLength(hsp->annotation[chromId_2-1].text,MAX_SEQ_NAME_LENGTH));
                            if (myNode->bucket[i].strand==QUERY_NEG_STRAND) {
                                SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",-",2);
                            } else if (myNode->bucket[i].strand==QUERY_POS_STRAND) {
                                SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",+",2);
                            } else {
                                SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",?",2);
                            }
                            SAMIUint8ConcatUint32AsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),tp_2);
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",",1);
                            SAMIUint8ConcatCigarAsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),
                                                        qInfo->ReadLength,
                                                        myNode->bucket[i].strand,
                                                        myNode->bucket[i].mismatchCount,
                                                        myNode->bucket[i].errors,SAM_ERROR_RANDOM);
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",",1);
                            SAMIUint8ConcatUint32AsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),myNode->bucket[i].mismatchCount);
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),";",1);
                                
                            samAlgnmt_mate.l_aux += samAlgnmt_mate.data_len - auxStart_mate;
                        //--------------------------
                    } else {
                        if (chromId_1==chromId_2) {
                            if (dpOcc_1->ambPosition > myNode->bucket[i].ambPosition) {
                                insertion = - (dpOcc_1->ambPosition + dpOcc_1->matchLen - myNode->bucket[i].ambPosition);
                            } else {
                                insertion = myNode->bucket[i].ambPosition + myNode->bucket[i].matchLen - dpOcc_1->ambPosition;
                            }
                        } else {
                            insertion = 0;
                        }
                                                
                        unsigned int samFlag = 0;
                        
                        samFlag |= SAM_FLAG_IS_PAIR_READ;
                        samFlag |= SAM_FLAG_PROPER_PAIR;
                        samFlag |= SAM_FLAG_FIRST_IN_PAIR;
                        if (dpOcc_1->strand==QUERY_NEG_STRAND) {
                            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
                        }
                        if (myNode->bucket[i].strand==QUERY_NEG_STRAND) {
                            samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
                        }
                        
                        SAMDPPackSequenceAndCigar(&samAlgnmt,sharedDataLen,qInfo,dpOcc_1);
                        samAlgnmt.l_aux = 0;
                        // Read Group at end of each line
                        if (peOutput_OccHandler->ReadGroup!= NULL && strlen( peOutput_OccHandler->ReadGroup ) >0){
                            unsigned int auxStart = samAlgnmt.data_len; 
                            bam_aux_append ( &samAlgnmt, "RG", 'Z', strlen (peOutput_OccHandler->ReadGroup  ) + 1, ( uint8_t * ) peOutput_OccHandler->ReadGroup  );
                            samAlgnmt.l_aux += samAlgnmt.data_len - auxStart;
                        }
                        // ---------------------->  |
                                                    samAlgnmt.core.tid = chromId_1-1;                           //SAM: chromosome ID, defined by bam_header_t
                                                    samAlgnmt.core.pos = tp_1-1;                                 //SAM: 0-based leftmost coordinate
                                                    samAlgnmt.core.flag = samFlag;                             //SAM: bitwise flag
                                                    
                                                    samAlgnmt.core.mtid = chromId_2-1;
                                                    samAlgnmt.core.mpos = tp_2-1;
                                                    samAlgnmt.core.isize = insertion;
                        // ---------------------->  |
                        
                        samFlag = 0;
                        
                        samFlag |= SAM_FLAG_IS_PAIR_READ;
                        samFlag |= SAM_FLAG_PROPER_PAIR;
                        samFlag |= SAM_FLAG_SECOND_IN_PAIR;
                        if (dpOcc_1->strand==QUERY_NEG_STRAND) {
                            samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
                        }
                        if (myNode->bucket[i].strand==QUERY_NEG_STRAND) {
                            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
                        }
                        
                        SAMPackSequenceAndCigar(&samAlgnmt_mate,sharedDataLen_mate,qInfo_mate,&(myNode->bucket[i]));
                        samAlgnmt_mate.l_aux = 0;
                        // Read Group at end of each line
                        if (peOutput_OccHandler->ReadGroup!= NULL && strlen( peOutput_OccHandler->ReadGroup ) >0){
                            unsigned int auxStart = samAlgnmt_mate.data_len; 
                            bam_aux_append ( &samAlgnmt_mate, "RG", 'Z', strlen (peOutput_OccHandler->ReadGroup  ) + 1, ( uint8_t * ) peOutput_OccHandler->ReadGroup  );
                            samAlgnmt_mate.l_aux += samAlgnmt_mate.data_len - auxStart;
                        }
                        // ---------------------->  |
                                                    samAlgnmt_mate.core.tid = chromId_2-1;                           //SAM: chromosome ID, defined by bam_header_t
                                                    samAlgnmt_mate.core.pos = tp_2-1;                                 //SAM: 0-based leftmost coordinate
                                                    samAlgnmt_mate.core.flag = samFlag;                             //SAM: bitwise flag
                                                    
                                                    samAlgnmt_mate.core.mtid = chromId_1-1;
                                                    samAlgnmt_mate.core.mpos = tp_1-1;
                                                    samAlgnmt_mate.core.isize = - insertion;
                        // ---------------------->  |
                        
                        samPrimaryOccAvailable = 1;
                        samPrimaryOccAvailable_mate = 1;
                    }
                }
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_PE_DP_SEED_OCC) {
            
                // 
                // The occurrences are found by PE-DP and the occurrences were stored in
                // dpOcc_1 and dpOcc_2.
                // The occurrences will then be stored into samAlgnmt and samAlgnmt_mate.
                //
                
                dpOcc_1 = DPOCCPTRead(dpOccIterator);DPOCCPTNext(dpOccIterator);
                dpOcc_2 = DPOCCPTRead(dpOccIterator);DPOCCPTNext(dpOccIterator);
                
                if (dpOcc_1==NULL || dpOcc_2==NULL) {
                    #ifdef PE_OCC_PRINT_OVERFLOW_MESSAGE
                        printf("[PEOccFlushCache/Read#%llu] DP Occurrence Collector does not have this occurrence.\n",qInfo->ReadId);
                    #endif
                } else {
                
                                            
                    OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                            dpOcc_1->ambPosition,
                                            &chromId_1,
                                            &tp_1);

                    OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                            dpOcc_2->ambPosition,
                                            &chromId_2,
                                            &tp_2);
                                            
                    if (samPrimaryOccAvailable || samPrimaryOccAvailable_mate) {
                    
                        //--------------------------
                            unsigned int auxStart = samAlgnmt.data_len;
                            if (!samXAZtagCreated) {
                                SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),"XAZ",3);
                                samXAZtagCreated=1;
                            }
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),
                                                    hsp->annotation[chromId_1-1].text,
                                                    SAMSeqReadNameLength(hsp->annotation[chromId_1-1].text,MAX_SEQ_NAME_LENGTH));
                            if (dpOcc_1->strand==QUERY_NEG_STRAND) {
                                SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",-",2);
                            } else if (dpOcc_1->strand==QUERY_POS_STRAND) {
                                SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",+",2);
                            } else {
                                SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",?",2);
                            }
                            SAMIUint8ConcatUint32AsString(samAlgnmt.data,&(samAlgnmt.data_len),tp_1);
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
                            unsigned int auxStart_mate = samAlgnmt_mate.data_len;
                            if (!samXAZtagCreated_mate) {
                                SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),"XAZ",3);
                                samXAZtagCreated_mate=1;
                            }
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),
                                                    hsp->annotation[chromId_2-1].text,
                                                    SAMSeqReadNameLength(hsp->annotation[chromId_2-1].text,MAX_SEQ_NAME_LENGTH));
                            if (dpOcc_2->strand==QUERY_NEG_STRAND) {
                                SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",-",2);
                            } else if (dpOcc_2->strand==QUERY_POS_STRAND) {
                                SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",+",2);
                            } else {
                                SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",?",2);
                            }
                            SAMIUint8ConcatUint32AsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),tp_2);
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",",1);
                            SAMDPIUint8ConcatCigarAsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),
                                                        qInfo->ReadLength,
                                                        dpOcc_2->strand,
                                                        dpOcc_2->matchElemsCount,
                                                        dpOcc_2->matchElems);
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",",1);
                            SAMIUint8ConcatUint32AsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),dpOcc_2->mismatchCount);
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),";",1);
                                
                            samAlgnmt_mate.l_aux += samAlgnmt_mate.data_len - auxStart_mate;
                        //--------------------------
                    } else {
                        if (chromId_1==chromId_2) {
                            if (dpOcc_1->ambPosition > dpOcc_2->ambPosition) {
                                insertion = - (dpOcc_1->ambPosition + dpOcc_1->matchLen - dpOcc_2->ambPosition);
                            } else {
                                insertion = dpOcc_2->ambPosition + dpOcc_2->matchLen - dpOcc_1->ambPosition;
                            }
                        } else {
                            insertion = 0;
                        }
                                                
                        unsigned int samFlag = 0;
                        
                        samFlag |= SAM_FLAG_IS_PAIR_READ;
                        samFlag |= SAM_FLAG_PROPER_PAIR;
                        samFlag |= SAM_FLAG_FIRST_IN_PAIR;
                        if (dpOcc_1->strand==QUERY_NEG_STRAND) {
                            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
                        }
                        if (dpOcc_2->strand==QUERY_NEG_STRAND) {
                            samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
                        }
                        
                        SAMDPPackSequenceAndCigar(&samAlgnmt,sharedDataLen,qInfo,dpOcc_1);
                        samAlgnmt.l_aux = 0;
                        // Read Group at end of each line
                        if (peOutput_OccHandler->ReadGroup!= NULL && strlen( peOutput_OccHandler->ReadGroup ) >0){
                            unsigned int auxStart = samAlgnmt.data_len; 
                            bam_aux_append ( &samAlgnmt, "RG", 'Z', strlen (peOutput_OccHandler->ReadGroup  ) + 1, ( uint8_t * ) peOutput_OccHandler->ReadGroup  );
                            samAlgnmt.l_aux += samAlgnmt.data_len - auxStart;
                        }
                        // ---------------------->  |
                                                    samAlgnmt.core.tid = chromId_1-1;                           //SAM: chromosome ID, defined by bam_header_t
                                                    samAlgnmt.core.pos = tp_1-1;                                 //SAM: 0-based leftmost coordinate
                                                    samAlgnmt.core.flag = samFlag;                             //SAM: bitwise flag
                                                    
                                                    samAlgnmt.core.mtid = chromId_2-1;
                                                    samAlgnmt.core.mpos = tp_2-1;
                                                    samAlgnmt.core.isize = insertion;
                        // ---------------------->  |
                        
                        samFlag = 0;
                        
                        samFlag |= SAM_FLAG_IS_PAIR_READ;
                        samFlag |= SAM_FLAG_PROPER_PAIR;
                        samFlag |= SAM_FLAG_SECOND_IN_PAIR;
                        if (dpOcc_1->strand==QUERY_NEG_STRAND) {
                            samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
                        }
                        if (dpOcc_2->strand==QUERY_NEG_STRAND) {
                            samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
                        }
                        
                        SAMDPPackSequenceAndCigar(&samAlgnmt_mate,sharedDataLen_mate,qInfo_mate,dpOcc_2);
                        samAlgnmt_mate.l_aux = 0;
                        // Read Group at end of each line
                        if (peOutput_OccHandler->ReadGroup!= NULL && strlen( peOutput_OccHandler->ReadGroup ) >0){
                            unsigned int auxStart = samAlgnmt_mate.data_len; 
                            bam_aux_append ( &samAlgnmt_mate, "RG", 'Z', strlen (peOutput_OccHandler->ReadGroup  ) + 1, ( uint8_t * ) peOutput_OccHandler->ReadGroup  );
                            samAlgnmt_mate.l_aux += samAlgnmt_mate.data_len - auxStart;
                        }
                        // ---------------------->  |
                                                    samAlgnmt_mate.core.tid = chromId_2-1;                           //SAM: chromosome ID, defined by bam_header_t
                                                    samAlgnmt_mate.core.pos = tp_2-1;                                 //SAM: 0-based leftmost coordinate
                                                    samAlgnmt_mate.core.flag = samFlag;                             //SAM: bitwise flag
                                                    
                                                    samAlgnmt_mate.core.mtid = chromId_1-1;
                                                    samAlgnmt_mate.core.mpos = tp_1-1;
                                                    samAlgnmt_mate.core.isize = - insertion;
                        // ---------------------->  |
                        
                        samPrimaryOccAvailable = 1;
                        samPrimaryOccAvailable_mate = 1;
                    }
                }
            } else if (myNode->bucket[i].type==SRAOCC_TYPE_MAPQ_VALUE) {
                unsigned long long packedMapqValue = myNode->bucket[i].ambPosition;
                int value0 = packedMapqValue >> 32;
                int value1 = packedMapqValue & UINT32_MAX;
                samAlgnmt.core.qual = value0;                       //SAM: mapping quality
                samAlgnmt_mate.core.qual = value1;                  //SAM: mapping quality
            }
            
            if (validPair==2 || invalidPair==2) {
                OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                        occIdx[0]->ambPosition,
                                        &chromId_1,
                                        &tp_1);

                OCCTranslateOccurrence(occTranslate, occAmbiguityMap, 
                                        occIdx[1]->ambPosition,
                                        &chromId_2,
                                        &tp_2);
                                        
                if (samPrimaryOccAvailable || samPrimaryOccAvailable_mate) {
                    //--------------------------
                        unsigned int auxStart = samAlgnmt.data_len;
                        if (!samXAZtagCreated) {
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),"XAZ",3);
                            samXAZtagCreated=1;
                        }
                        SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),
                                                hsp->annotation[chromId_1-1].text,
                                                SAMSeqReadNameLength(hsp->annotation[chromId_1-1].text,MAX_SEQ_NAME_LENGTH));
                        if (occIdx[0]->strand==QUERY_NEG_STRAND) {
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",-",2);
                        } else if (occIdx[0]->strand==QUERY_POS_STRAND) {
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",+",2);
                        } else {
                            SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",?",2);
                        }
                        SAMIUint8ConcatUint32AsString(samAlgnmt.data,&(samAlgnmt.data_len),tp_1);
                        SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",",1);
                        SAMIUint8ConcatCigarAsString(samAlgnmt.data,&(samAlgnmt.data_len),
                                                    qInfo->ReadLength,
                                                    occIdx[0]->strand,
                                                    occIdx[0]->mismatchCount,
                                                    occIdx[0]->errors,SAM_ERROR_RANDOM);
                        SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),",",1);
                        SAMIUint8ConcatUint32AsString(samAlgnmt.data,&(samAlgnmt.data_len),occIdx[0]->mismatchCount);
                        SAMIUint8ConcatString(samAlgnmt.data,&(samAlgnmt.data_len),";",1);
                            
                        samAlgnmt.l_aux += samAlgnmt.data_len - auxStart;
                    //--------------------------
                        unsigned int auxStart_mate = samAlgnmt_mate.data_len;
                        if (!samXAZtagCreated_mate) {
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),"XAZ",3);
                            samXAZtagCreated_mate=1;
                        }
                        SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),
                                                hsp->annotation[chromId_2-1].text,
                                                SAMSeqReadNameLength(hsp->annotation[chromId_2-1].text,MAX_SEQ_NAME_LENGTH));
                        if (occIdx[1]->strand==QUERY_NEG_STRAND) {
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",-",2);
                        } else if (occIdx[1]->strand==QUERY_POS_STRAND) {
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",+",2);
                        } else {
                            SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",?",2);
                        }
                        SAMIUint8ConcatUint32AsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),tp_2);
                        SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",",1);
                        SAMIUint8ConcatCigarAsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),
                                                    qInfo_mate->ReadLength,
                                                    occIdx[1]->strand,
                                                    occIdx[1]->mismatchCount,
                                                    occIdx[1]->errors,SAM_ERROR_RANDOM);
                        SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),",",1);
                        SAMIUint8ConcatUint32AsString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),occIdx[1]->mismatchCount);
                        SAMIUint8ConcatString(samAlgnmt_mate.data,&(samAlgnmt_mate.data_len),";",1);
                            
                        samAlgnmt_mate.l_aux += samAlgnmt_mate.data_len - auxStart_mate;
                    //--------------------------
                } else {
                    
                    if (chromId_1==chromId_2) {
                        if (occIdx[0]->ambPosition > occIdx[1]->ambPosition) {
                            insertion = - (occIdx[0]->ambPosition + occIdx[0]->matchLen - occIdx[1]->ambPosition);
                        } else {
                            insertion = occIdx[1]->ambPosition + occIdx[1]->matchLen - occIdx[0]->ambPosition;
                        }
                    } else {
                        insertion = 0;
                    }
                                            
                    unsigned int samFlag = 0;
                    
                    samFlag |= SAM_FLAG_IS_PAIR_READ;
                    samFlag |= SAM_FLAG_PROPER_PAIR * (validPair==2);
                    samFlag |= SAM_FLAG_FIRST_IN_PAIR;
                    if (occIdx[0]->strand==QUERY_NEG_STRAND) {
                        samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
                    }
                    if (occIdx[1]->strand==QUERY_NEG_STRAND) {
                        samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
                    }
                    
                    SAMPackSequenceAndCigar(&samAlgnmt,sharedDataLen,qInfo,occIdx[0]);
                    samAlgnmt.l_aux = 0;
                    // Read Group at end of each line
                    if (peOutput_OccHandler->ReadGroup!= NULL && strlen( peOutput_OccHandler->ReadGroup ) >0){
                        unsigned int auxStart = samAlgnmt.data_len; 
                        bam_aux_append ( &samAlgnmt, "RG", 'Z', strlen (peOutput_OccHandler->ReadGroup  ) + 1, ( uint8_t * ) peOutput_OccHandler->ReadGroup  );
                        samAlgnmt.l_aux += samAlgnmt.data_len - auxStart;
                    }
                    // ---------------------->  |
                                                samAlgnmt.core.tid = chromId_1-1;                           //SAM: chromosome ID, defined by bam_header_t
                                                samAlgnmt.core.pos = tp_1-1;                                 //SAM: 0-based leftmost coordinate
                                                samAlgnmt.core.flag = samFlag;                             //SAM: bitwise flag
                                                
                                                samAlgnmt.core.mtid = chromId_2-1;
                                                samAlgnmt.core.mpos = tp_2-1;
                                                samAlgnmt.core.isize = insertion;
                    // ---------------------->  |
                    
                    samFlag = 0;
                    
                    samFlag |= SAM_FLAG_IS_PAIR_READ;
                    samFlag |= SAM_FLAG_PROPER_PAIR * (validPair==2);
                    samFlag |= SAM_FLAG_SECOND_IN_PAIR;
                    if (occIdx[0]->strand==QUERY_NEG_STRAND) {
                        samFlag |= SAM_FLAG_MATE_ALGNMT_STRAND;
                    }
                    if (occIdx[1]->strand==QUERY_NEG_STRAND) {
                        samFlag |= SAM_FLAG_READ_ALGNMT_STRAND;
                    }
                    
                    SAMPackSequenceAndCigar(&samAlgnmt_mate,sharedDataLen_mate,qInfo_mate,occIdx[1]);
                    samAlgnmt_mate.l_aux = 0;
                    // Read Group at end of each line
                    if (peOutput_OccHandler->ReadGroup!= NULL && strlen( peOutput_OccHandler->ReadGroup ) >0){
                        unsigned int auxStart = samAlgnmt_mate.data_len; 
                        bam_aux_append ( &samAlgnmt_mate, "RG", 'Z', strlen (peOutput_OccHandler->ReadGroup  ) + 1, ( uint8_t * ) peOutput_OccHandler->ReadGroup  );
                        samAlgnmt_mate.l_aux += samAlgnmt_mate.data_len - auxStart;
                    }
                    // ---------------------->  |
                                                samAlgnmt_mate.core.tid = chromId_2-1;                           //SAM: chromosome ID, defined by bam_header_t
                                                samAlgnmt_mate.core.pos = tp_2-1;                                 //SAM: 0-based leftmost coordinate
                                                samAlgnmt_mate.core.flag = samFlag;                             //SAM: bitwise flag
                                                
                                                samAlgnmt_mate.core.mtid = chromId_1-1;
                                                samAlgnmt_mate.core.mpos = tp_1-1;
                                                samAlgnmt_mate.core.isize = - insertion;
                    // ---------------------->  |
                    
                    samPrimaryOccAvailable = 1;
                    samPrimaryOccAvailable_mate = 1;
                }

                validPair = 0;
                invalidPair = 0;
            }
        }
        myNode = myNode->next;
    }
    
    if (samPrimaryOccAvailable) {
        unsigned int auxStart = samAlgnmt.data_len;
        int tmpDataLen = samAlgnmt.data_len;
        SAMIUint8ConcatString(samAlgnmt.data,&(tmpDataLen),"\0",1);
        samAlgnmt.l_aux += samAlgnmt.data_len - auxStart;
        #ifndef DEBUG_2BWT_NO_WRITE_FILE
        samwrite(samFilePtr,&samAlgnmt);
        #endif
    }
    if (samPrimaryOccAvailable_mate) {
        unsigned int auxStart_mate = samAlgnmt_mate.data_len;
        int tmpDataLen_mate = samAlgnmt_mate.data_len;
        SAMIUint8ConcatString(samAlgnmt_mate.data,&(tmpDataLen_mate),"\0",1);
        samAlgnmt_mate.l_aux += samAlgnmt_mate.data_len - auxStart_mate;
        #ifndef DEBUG_2BWT_NO_WRITE_FILE
        samwrite(samFilePtr,&samAlgnmt_mate);
        #endif
    }
    
    DPOCCPTFree(dpOccIterator);
    SRAOCCInitialise(occCollector);
    DPOCCInitialise(dpOccCollector);
}

void PEOCCFlushCache(PEArguments * peArgs) {
    SRAOutput * peOutput_OccHandler = peArgs->PEAlgnmtOutput_OccHandle;
    SRAOCCCollector * occCollector = peOutput_OccHandler->occCollector;
    
    DPArguments * dpArgs = peArgs->dpArguments;
    DPOCCCollector * dpOccCollector = dpArgs->dpOccCollector;
    
    switch (peOutput_OccHandler->OutFileFormat) {
        case SRA_OUTPUT_FORMAT_PLAIN:
            PEOCCFlushCachePlain(peArgs);
            break;
        case SRA_OUTPUT_FORMAT_SAM:
        case SRA_OUTPUT_FORMAT_BAM:
            PEOCCFlushCacheSAM(peArgs);
            break;
        case SRA_OUTPUT_FORMAT_SAM_STORE:
            PEOCCFlushCacheSAM(peArgs);
            break;
        default:
            printf("Occurrences are dropped as the output format is invalid\n");
            SRAOCCInitialise(occCollector);
            DPOCCInitialise(dpOccCollector);
    }
    
}

void PEOCCReportNoAlignment(PEArguments * peArgs,uint8_t type) {

    SRAQueryInfo * qInfo = peArgs->sraArgsPos->QueryInfo;
    SRAQueryInfo * qInfo_mate = peArgs->sraArgsPos_mate->QueryInfo;
    SRAOutput * peOutput_OccHandler = peArgs->PEAlgnmtOutput_OccHandle;
    
    switch (peOutput_OccHandler->OutFileFormat) {
        case SRA_OUTPUT_FORMAT_PLAIN:
        case SRA_OUTPUT_FORMAT_SAM:
        case SRA_OUTPUT_FORMAT_BAM:
        case SRA_OUTPUT_FORMAT_SAM_STORE:
            PEOCCAddTextPositionToCache(peArgs,type,qInfo->ReadId,
                                        0,0,0,0,NULL);
            break;
    }
    
}

unsigned int PEOCCCountAllPEAlignment(PEArguments * peArgs) {

    SRAOutput * peOutput_OccHandler = peArgs->PEAlgnmtOutput_OccHandle;
    SRAOCCCollector * occCollector = peOutput_OccHandler->occCollector;
    
    DPArguments * dpArgs = peArgs->dpArguments;
    DPOCCCollector * dpOccCollector = dpArgs->dpOccCollector;
    
    unsigned int numOfOcc = SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_PE_PAIR)/2 + 
                            SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_PE_BAD_PAIR)/2;
                                
    unsigned int numOfDpOcc = SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_PE_DP_BASE_READ) + 
                              SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_PE_DP_BASE_MATE) + 
                              SRAOCCCountOccurrencesByType(occCollector,SRAOCC_TYPE_PE_DP_SEED_OCC);
                              
    return numOfOcc+numOfDpOcc;
}

void PEOCCReportSubValidAlignment(PEArguments * peArgs, SRAOccurrence * occ_1, SRAOccurrence * occ_2) {
    SRAQueryInfo * qInfo = peArgs->sraArgsPos->QueryInfo;
    SRAQueryInfo * qInfo_mate = peArgs->sraArgsPos_mate->QueryInfo;
    SRAOutput * peOutput_OccHandler = peArgs->PEAlgnmtOutput_OccHandle;
    
    switch (peOutput_OccHandler->OutFileFormat) {
        case SRA_OUTPUT_FORMAT_PLAIN:
        case SRA_OUTPUT_FORMAT_SAM:
        case SRA_OUTPUT_FORMAT_BAM:
        case SRA_OUTPUT_FORMAT_SAM_STORE:
            PEOCCAddTextPositionToCache(peArgs,SRAOCC_TYPE_PE_BAD_PAIR,occ_1->ambPosition,
                                        occ_1->strand,occ_1->matchLen,occ_1->mismatchCount,0,occ_1->errors);
            PEOCCAddTextPositionToCache(peArgs,SRAOCC_TYPE_PE_BAD_PAIR,occ_2->ambPosition,
                                        occ_2->strand,occ_2->matchLen,occ_2->mismatchCount,0,occ_2->errors);
            break;
            break;
    }
    
}

void PEOCCReportDelimitor(PEArguments * peArgs) {

    SRAQueryInfo * qInfo = peArgs->sraArgsPos->QueryInfo;
    SRAOutput * peOutput_OccHandler = peArgs->PEAlgnmtOutput_OccHandle;
    
    switch (peOutput_OccHandler->OutFileFormat) {
        case SRA_OUTPUT_FORMAT_PLAIN: 
        case SRA_OUTPUT_FORMAT_SAM_STORE:
            PEOCCAddTextPositionToCache(peArgs,SRAOCC_TYPE_DELIMITOR_READ,qInfo->ReadId,
                                        0,0,0,0,NULL);
            break;
        case SRA_OUTPUT_FORMAT_SAM:
        case SRA_OUTPUT_FORMAT_BAM:
            break;
    }
}

unsigned long long PEOCCDumpOneAlignment(PEArguments * peArgs, SRAOccurrence * occ_1, SRAOccurrence * occ_2) {
    PEOCCAddTextPositionToCache(peArgs,SRAOCC_TYPE_PE_PAIR,occ_1->ambPosition,
                                occ_1->strand,occ_1->matchLen,occ_1->mismatchCount,0,occ_1->errors);
    PEOCCAddTextPositionToCache(peArgs,SRAOCC_TYPE_PE_PAIR,occ_2->ambPosition,
                                occ_2->strand,occ_2->matchLen,occ_2->mismatchCount,0,occ_2->errors);
    return 1;
}

unsigned long long PEOCCDumpAlignments(PEArguments * peArgs) {
    SRAOutput * peOutput_OccHandler = peArgs->PEAlgnmtOutput_OccHandle;
    PEOutput * peOutput = peArgs->PEAlgnmtOutput;
    PEPairList * pairList = peOutput->root;
    PEInput * peInput = peArgs->PEAlgnmtInput;
    //printf("[PEOCCDumpAlignments] %u\n",pairList->pairsCount);
    
    unsigned long long occCount = 0;
    unsigned int i;
    
    while (pairList!=NULL && pairList->pairsCount>0) {
        for (i=0;i<pairList->pairsCount;i++) {
            PEPairs * pePair = &(pairList->pairs[i]);
            SRAOccurrence * occ_1 = pePair->occ_1;
            SRAOccurrence * occ_2 = pePair->occ_2;
            switch (peOutput_OccHandler->OutFileFormat) {
                case SRA_OUTPUT_FORMAT_PLAIN: 
                case SRA_OUTPUT_FORMAT_SAM:
                case SRA_OUTPUT_FORMAT_BAM:
                case SRA_OUTPUT_FORMAT_SAM_STORE:
                    PEOCCAddTextPositionToCache(peArgs,SRAOCC_TYPE_PE_PAIR,occ_1->ambPosition,
                                                occ_1->strand,occ_1->matchLen,occ_1->mismatchCount,0,occ_1->errors);
                    PEOCCAddTextPositionToCache(peArgs,SRAOCC_TYPE_PE_PAIR,occ_2->ambPosition,
                                                occ_2->strand,occ_2->matchLen,occ_2->mismatchCount,0,occ_2->errors);
                    if (peInput->OutputType == PE_REPORT_RANDOM_BEST) return 1;
                
                    break;
            }
        }
        occCount += pairList->pairsCount;
        pairList = pairList->next;
    }

    return occCount;
}

unsigned long long PEOCCDumpSingleAlignments(PEArguments * peArgs, SRAOccurrence * occ, unsigned int count, uint8_t type) {
    SRAOutput * peOutput_OccHandler = peArgs->PEAlgnmtOutput_OccHandle;
    PEInput * peInput = peArgs->PEAlgnmtInput;
    PEOutput * peOutput = peArgs->PEAlgnmtOutput;
    
    unsigned int i;
    for (i=0;i<count;i++) {
        switch (peOutput_OccHandler->OutFileFormat) {
            case SRA_OUTPUT_FORMAT_PLAIN: 
            case SRA_OUTPUT_FORMAT_SAM:
            case SRA_OUTPUT_FORMAT_BAM:
            case SRA_OUTPUT_FORMAT_SAM_STORE:
                PEOCCAddTextPositionToCache(peArgs,type,occ[i].ambPosition,
                                            occ[i].strand,occ[i].matchLen,
                                            occ[i].mismatchCount,0,occ[i].errors);
                if (peInput->OutputType == PE_REPORT_RANDOM_BEST) return 1;
            
                break;
        }
    }
    return count;
}


unsigned long long PEOCCMAPQValue(PEArguments * peArgs, unsigned long long value0, unsigned long long value1) {
    unsigned long long packedMapqValue = value0<<32 | value1;
    PEOCCAddTextPositionToCache(peArgs,SRAOCC_TYPE_MAPQ_VALUE,packedMapqValue,
                                0,0,0,0,NULL);
    return 1;
}


void PEOCCAddTextPositionToCache(PEArguments * peArgs, 
                                    uint8_t type,
                                    unsigned long long ambPosition, 
                                    uint8_t strand, uint16_t matchLen,
                                    uint8_t occMismatch, int occQuality,
                                    const SRAError errorPos[]) { 
                                    
    SRAOutput * peOutput = peArgs->PEAlgnmtOutput_OccHandle;
    SRAOCCCollector * occCollector = peOutput->occCollector;
    SRAOccurrence occ;
    int i;

    occ.type             = type;
    occ.ambPosition      = ambPosition;
    occ.strand           = strand;
    occ.mismatchCount    = occMismatch;
    occ.matchLen         = matchLen;

    for (i=0;i<occMismatch&&i<MAX_NUM_OF_ERROR;i++) {
        occ.errors[i].type = errorPos[i].type;
        occ.errors[i].position = errorPos[i].position;
    }
        
    if (!SRAOCCAddOccurrenceToBuckets(occCollector,&occ)) {
        PEOCCFlushCache(peArgs);
        if (!SRAOCCAddOccurrenceToBuckets(occCollector,&occ)) {
            printf("[PEOCCAddTextPositionToCache] Unexpected error. OccCollector not clean. \n");
            exit(1);
        }
    }
}


void MAPQPEGetScore (PEArguments * peArg ,int * mapScore0, int *mapScore1) {
    SRAArguments * sraArgsPos      = peArg->sraArgsPos;
    SRAArguments * sraArgsPos_mate = peArg->sraArgsPos_mate;
    
    MAPQCalculator * readMapqCalc = sraArgsPos->MapqCalc;
    MAPQCalculator * mateMapqCalc = sraArgsPos_mate->MapqCalc;
    MAPQCalculator * peMapqCalc = peArg->MapqCalc;
    
    // x0_0 = sraArgsPos->MapqCalc->rankedCount[0]
    // x1_0 = sraArgsPos_mate->MapqCalc->rankedCount[0]
    // x0_1 = sraArgsPos->MapqCalc->rankedCount[1]
    // x1_1 = sraArgsPos_mate->MapqCalc->rankedCount[1]
    
    // g_log_n = peArg->MapqCalc->g_log_n
    
    // op_score         = peArg->MapqCalc->rankedScore[0]
    // op_num           = peArg->MapqCalc->rankedCount[0]
    // subop_score      = peArg->MapqCalc->rankedScore[1]
    // subop_num        = peArg->MapqCalc->rankedCount[1]
    
    // readlen_0        = sraArgsPos->QueryInfo->ReadLength
    // readlen_1        = sraArgsPos_mate->QueryInfo->ReadLength
    
    // map_score0 map_score1..
    
    bwaLikePairQualScore ( 
        readMapqCalc->rankedCount[0],
        mateMapqCalc->rankedCount[0],
        readMapqCalc->rankedCount[1],
        mateMapqCalc->rankedCount[1],
        
        peMapqCalc->g_log_n,
        
        peMapqCalc->rankedScore[0],
        peMapqCalc->rankedCount[0],
        peMapqCalc->rankedScore[1],
        peMapqCalc->rankedCount[1],
        
        sraArgsPos->QueryInfo->ReadLength,
        sraArgsPos_mate->QueryInfo->ReadLength,
        
        mapScore0,mapScore1
    );

}
