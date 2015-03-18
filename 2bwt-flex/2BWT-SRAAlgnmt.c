//
//    2BWT-SRAAlgnmt.c
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
   
   Date   : 11th April 2012
   Author : Edward MK Wu
   Change : New file.

*/
/////////////////////////////////////////////////////

#include "2BWT-SRAAlgnmt.h"

// --------------------------------------------------
// SRA_DEBUG #1
// Uncomment the below to turn on logging message for
// MAPQ Calculation process.
//#define SRA_DEBUG_PRINT_MAPQ_CALCULATION
// --------------------------------------------------
// SRA_DEBUG #2
// Uncomment the below to turn on logging message for
// SRA matching for seed enhancement.
//#define SRA_DEBUG_PRINT_SEED_MATCHING
// --------------------------------------------------


unsigned long long ProcessStage(SRAArguments * aArgs, SRAModel * sraModel, int * caseId) {
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAIndex * aIndex = aArgs->AlgnmtIndex;
    SRAWorkingMemory * workMem = aArgs->AlgnmtMemory;
    BWT * bwt = aIndex->bwt;
    int OutputType = aSetting->OutputType;
    
    unsigned long long saRanges[4];
    unsigned long long saCount = 0;
    SRACase * thisCase = &(sraModel->cases[(*caseId)]);
    
    while (thisCase->type != SRA_CASE_TYPE_NOT_INITALISED && 
           thisCase->type != SRA_CASE_TYPE_NEXT_STAGE) {
           
        saRanges[0] = 0;
        saRanges[1] = bwt->textLength;
        saRanges[2] = 0;
        saRanges[3] = bwt->textLength;

        if (thisCase->steps[0].type == SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP) {
            saCount += BWTExactModelBackward_Lookup(aArgs,thisCase,0,saRanges);
        } else if (thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP) {
            saCount += BWTExactModelForward_Lookup(aArgs,thisCase,0,saRanges);
        } else if (thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_BACKWARD_LOOKUP) {
            saCount += BWTExactModelBackwardAnyDirection_Lookup(aArgs,thisCase,0,saRanges);
        } else if (thisCase->steps[0].type == SRA_STEP_TYPE_BACKWARD_ONLY_BWT) {
            saCount += BWTModelSwitchBackward(aArgs,0,0,
                                                thisCase,0,
                                                saRanges, 0, 0, 0);
        } else if (thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_BWT) {
            saCount += BWTModelSwitchAnyDirection(aArgs,0,0,
                                                    thisCase,0,
                                                    saRanges, 0, 0, 0);
        }
        if (workMem->IsClosed==1) {
            if (OutputType==SRA_REPORT_RANDOM_BEST) {
                OCCFlushBestIntoCache(aArgs);
            }
            return saCount;
        }
        (*caseId)++;
        thisCase = &(sraModel->cases[(*caseId)]);
    }
    
    return saCount;
}

unsigned long long ProcessReadDoubleStrand(SRAArguments * aArgs, SRAArguments * aArgs_neg, SRAModel * sraModel, SRAModel * sraModel_neg) {
    unsigned long long saCount = 0;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAWorkingMemory * workMem = aArgs->AlgnmtMemory;
    int OutputType = aSetting->OutputType;
    int j,k;
    
    SRAWorkingMemoryInitialise(workMem);

    j=0;
    while (1) {

        //Report unique best alignment
        //Report all best alignment
        k=j;
        saCount += ProcessStage(aArgs,sraModel,&k);
        
        if (OutputType==SRA_REPORT_RANDOM_BEST &&
            workMem->IsClosed==1) {
            
            return saCount;
        }

        k=j;
        saCount += ProcessStage(aArgs_neg,sraModel_neg,&k);
        
        if (OutputType==SRA_REPORT_RANDOM_BEST &&
            workMem->IsClosed==1) {
            
            return saCount;
        }

        if (OutputType==SRA_REPORT_UNIQUE_BEST) {
            if (workMem->IsUnique==1) {
                OCCFlushBestIntoCache(aArgs);
                workMem->IsClosed = 1;
                return 1;
            } else if (workMem->IsOpened==1) {
                workMem->IsClosed = 1;
                return 0;
            }
        } 

        if (OutputType==SRA_REPORT_ALL_BEST && saCount>0) {
            return saCount;
        }
        
        if (sraModel->cases[k].type==SRA_CASE_TYPE_NOT_INITALISED ||
            sraModel_neg->cases[k].type==SRA_CASE_TYPE_NOT_INITALISED) {
            if (sraModel->cases[k].type!=SRA_CASE_TYPE_NOT_INITALISED ||
                sraModel_neg->cases[k].type!=SRA_CASE_TYPE_NOT_INITALISED) {
                fprintf(stderr,"[ERROR] Imbalance SRA Alignment model!\n");
                exit(1);
            } else {
                break;
            }
        }

        j=k+1;
    }

    if (OutputType==SRA_REPORT_BEST_QUALITY && workMem->IsOpened==1) {
        OCCFlushBestIntoCache(aArgs);
        return 1;
    }

    return saCount;

}


unsigned long long ProcessReadSingleStrand(SRAArguments * aArgs, SRAModel * sraModel) {
    unsigned long long saCount = 0;
    SRASetting * aSetting = aArgs->AlgnmtSetting;
    SRAWorkingMemory * workMem = aArgs->AlgnmtMemory;
    int OutputType = aSetting->OutputType;
    int j;
    

    SRAWorkingMemoryInitialise(workMem);

    j=0;
    while (1) {

        saCount += ProcessStage(aArgs,sraModel,&j);
        
        if (OutputType==SRA_REPORT_RANDOM_BEST &&
            workMem->IsClosed==1) {
            
            return saCount;
        }
        
        if (OutputType==SRA_REPORT_UNIQUE_BEST) {
            if (workMem->IsUnique==1) {
                OCCFlushBestIntoCache(aArgs);
                workMem->IsClosed = 1;
                return 1;
            } else if (workMem->IsOpened==1) {
                workMem->IsClosed = 1;
                return 0;
            }
        }
        
        if (OutputType==SRA_REPORT_ALL_BEST && saCount>0) {
            return saCount;
        }
        
        if (sraModel->cases[j].type==SRA_CASE_TYPE_NOT_INITALISED) {
            break;
        }
        
        j++;
    }

    if (OutputType==SRA_REPORT_BEST_QUALITY && workMem->IsOpened==1) {
        OCCFlushBestIntoCache(aArgs);
        return 1;
    }

    return saCount;

}


unsigned long long SRASeedProcessReadDbl(SRAArguments * aArgsPos, SRAArguments * aArgsNeg, 
                                        SRAModel * sraModel, SRAModel * sraModel_neg,
                                        int seedLength, int seedUniqeLength, SRAOCCCollector * sraOccCollector ) {
                                        
    SRAQueryInfo * qInfoPos = aArgsPos->QueryInfo;
    SRAQueryInfo * qInfoNeg = aArgsNeg->QueryInfo;
    
    int properlySeeded = 1;
    int seedCorrection;
    int seedOffset = 0;
    SRAOCCCollector * seedOccCollector = aArgsPos->AlgnmtOutput->occCollector;
    SRAOCCCollectorPointer * seedReadOccPtr = SRAOCCPTCreate(seedOccCollector);
    SRAOccurrence * seedOcc;
    
    int readLength = qInfoPos->ReadLength;
    
    // Initial state
    qInfoPos->ReadLength         = seedLength;
    // Initial state
    qInfoNeg->ReadLength         = seedLength;
    qInfoNeg->ReadCode          += readLength - seedLength;
    
    while (properlySeeded && seedOffset + seedLength <= readLength) {
    
        #ifdef SRA_DEBUG_PRINT_SEED_MATCHING
            printf("[SRASeedProcessReadDbl] seedOffset=%d seedLength=%d seedUniqeLength=%d\n",seedOffset,seedLength,seedUniqeLength);
            DEBUGSRAQueryInfoPrint(qInfoPos,"SEED-POS");
            DEBUGSRAQueryInfoPrint(qInfoNeg,"SEED-NEG");
        #endif

        ProcessReadDoubleStrand(aArgsPos,aArgsNeg,sraModel,sraModel_neg);
        
        #ifdef SRA_DEBUG_PRINT_SEED_MATCHING
            printf("[SRASeedProcessReadDbl] %u alignments were found.\n",SRAOCCCountOccurrencesByType(seedOccCollector,SRAOCC_TYPE_AWAIT_TRANSLATE));
        #endif
        
        if (aArgsPos->AlgnmtMemory->IsClosed==1) {
            properlySeeded = 0;
            SRAOCCInitialise(seedOccCollector);
            break;
        }
        
        SRAOCCPTRoot(seedOccCollector,seedReadOccPtr);
        seedOcc = SRAOCCPTRead(seedReadOccPtr);
        while (seedOcc!=NULL) {
        
            if (seedOcc->type == SRAOCC_TYPE_AWAIT_TRANSLATE) {
                #ifdef SRA_DEBUG_PRINT_SEED_MATCHING
                    printf("[sraQueryInputPos] Found seed at this position = %llu\n",seedOcc->ambPosition);
                #endif
                
                if (seedOcc->strand == QUERY_POS_STRAND) {
                    seedCorrection = seedOffset;
                } else {
                    seedCorrection = (readLength - seedOffset - seedLength);
                }
                
                // Ignoring the seed alignment in case of boundary
                if (seedOcc->ambPosition >= seedCorrection) {
                
                    seedOcc->ambPosition -= seedCorrection;
                        
                    #ifdef SRA_DEBUG_PRINT_SEED_MATCHING
                        uint8_t seedChrId;
                        SRAIndex * sraIndex = aArgsPos->AlgnmtIndex;
                        unsigned long long seedUnambPos;
                        OCCTranslateOccurrence(sraIndex->hsp->translate,sraIndex->hsp->ambiguityMap,seedOcc->ambPosition,&seedChrId,&seedUnambPos);
                        printf("[SRASeedProcessReadDbl] Found corrected at this position = %llu (c%do%llu)\n",seedOcc->ambPosition,seedChrId,seedUnambPos);
                    #endif
                    
                    SRAOCCAddOccurrenceToBuckets(sraOccCollector,seedOcc);
                }
            }
        
            SRAOCCPTNext(seedReadOccPtr);
            seedOcc = SRAOCCPTRead(seedReadOccPtr);
        }
        
        SRAOCCInitialise(seedOccCollector);
        seedOffset += seedUniqeLength;
        
        qInfoPos->ReadCode               += seedUniqeLength;
        qInfoNeg->ReadCode               -= seedUniqeLength;
        
        qInfoPos->ReportingReadCode      += seedUniqeLength;
        qInfoNeg->ReportingReadCode      += seedUniqeLength;
    }

    SRAOCCPTFree(seedReadOccPtr);
    
    return properlySeeded;
}

void SRAPopulateMAPQCalculator(SRAArguments * aArgsPos, SRAArguments * aArgsNeg, DPScores * dpScores, SRAModel * sraModels_Extend, SRAModel * sraModelsNeg_Extend) {

    SRASetting * sraAlgnmtSetting = aArgsPos->AlgnmtSetting;
    SRAResultCount * sraAlgnmtStats = aArgsPos->AlgnmtStats;
    MAPQCalculator * mapqCalc = aArgsPos->MapqCalc;

    
    // Get SRA Optimal and Sub-Optimal
    int resultRank = 0;
    int firstHit = 0;
    int i;
    for (i=0;i<=MAX_NUM_OF_ERROR&&i<=sraAlgnmtSetting->MaxError;i++) {
        if (sraAlgnmtStats->WithError[i]>0 || firstHit) {
            mapqCalc->rankedMismatch[resultRank] = i;
            mapqCalc->rankedCount[resultRank] = sraAlgnmtStats->WithError[i];
            mapqCalc->rankedScore[resultRank] = MAPQNormaliseSRAScore(aArgsPos->QueryInfo->ReadLength,i,dpScores);
            resultRank++;
            firstHit=1;
        }
        
        if (resultRank>=2) {
            break;
        }
    }
    
    if (resultRank == 1) {
    
        #ifdef SRA_DEBUG_PRINT_MAPQ_CALCULATION
            printf("[SRAPopulateMAPQCalculator] Only found optimal. Enriching results.\n");
            SRAPrintResultCount(aArgsPos->AlgnmtStats);
            MAPQPrintParameters(mapqCalc);
        #endif

        // If there isn't an sub-optimal alignment results
        // we will need to enrich the quality data with a
        // one-more-mismatch SRA alignment.
        
        // SRASetting for "One-More-Mismatch" 
        SRASetting * extSetting = SRASettingMakeClone(sraAlgnmtSetting);
        extSetting->OutputType = SRA_REPORT_NONE;
        
        SRAResultCount * extStats       = SRAResultCountConstruct();
        SRAArguments * extPosArg        = SRAARGMakeClone(aArgsPos);
        SRAArguments * extNegArg        = SRAARGMakeClone(aArgsNeg);
        
        extPosArg->AlgnmtSetting = extSetting;
        extNegArg->AlgnmtSetting = extSetting;
        
        extPosArg->AlgnmtStats = extStats;
        extNegArg->AlgnmtStats = extStats;
        
        // Get SRA Count
        ProcessReadDoubleStrand(extPosArg,extNegArg,sraModels_Extend,sraModels_Extend);
        
        // Steal figures from SRA AlgnmtStats
        i = mapqCalc->rankedMismatch[resultRank-1] + 1;
        mapqCalc->rankedMismatch[resultRank] = i;
        mapqCalc->rankedCount[resultRank] = extStats->WithError[i];
        mapqCalc->rankedScore[resultRank] = MAPQNormaliseSRAScore(aArgsPos->QueryInfo->ReadLength,i,dpScores);
        resultRank++;
            
        #ifdef SRA_DEBUG_PRINT_MAPQ_CALCULATION
            printf("[SRAPopulateMAPQCalculator] After enrichment:\n");
            SRAPrintResultCount(extPosArg->AlgnmtStats);
            MAPQPrintParameters(mapqCalc);
        #endif
        
        // Memory Clean-up.
        SRAResultCountFree(extStats);
        SRAARGCloneFree(extPosArg);
        SRAARGCloneFree(extNegArg);
        SRASettingFree(extSetting);
    } else {
        
        #ifdef SRA_DEBUG_PRINT_MAPQ_CALCULATION
            printf("[SRAPopulateMAPQCalculator] MAPQ parameters found. No enrichment needed.\n");
            SRAPrintResultCount(aArgsPos->AlgnmtStats);
            MAPQPrintParameters(mapqCalc);
        #endif
    }
    
    mapqCalc->status = MAPQ_CALC_FIGURE_READY;
    
}
