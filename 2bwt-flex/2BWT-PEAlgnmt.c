//
//    2BWT-PEAlgnmt.c
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
   
   Date   : 11th March 2012
   Author : Edward MK Wu
   Change : New file.

*/
/////////////////////////////////////////////////////

#include "2BWT-PEAlgnmt.h"
#include "MappingQuality.h"
void SRAPopulateMAPQCalculator(SRAArguments * aArgsPos, SRAArguments * aArgsNeg, DPScores * dpScores, SRAModel * sraModels_Extend, SRAModel * sraModelsNeg_Extend);
unsigned long long SRASeedProcessReadDbl(SRAArguments * aArgsPos, SRAArguments * aArgsNeg,
                                        SRAModel * sraModel, SRAModel * sraModel_neg,
                                        int seedLength, int seedUniqeLength, SRAOCCCollector * sraOccCollector ) ;
// --------------------------------------------------
// PE_DEBUG #1
// Uncomment the below to turn on logging message for
// MAPQ Calculation process.
//#define PE_DEBUG_PRINT_MAPQ_CALCULATION
// --------------------------------------------------

// --------------------------------------------------
// DP_DEBUG #1
// Uncomment the below to turn on logging message for
// PE matching for seed enhancement.
//#define DP_DEBUG_PRINT_SEED_MATCHING
// --------------------------------------------------

// 2BWT-PE CTC#1
// -------------
// Set this parameter PE_DUMP_SRA_OUTPUT to 1 if you want
// 2BWT-PE to also dump out all SRA alignment along with
// the pair-end alignment results. Dumping SRA output only
// works with All-Valid Alignment. Nothing will be output
// for the other type of alignment.
#define PE_DUMP_SRA_OUTPUT 0
// -------------
// 2BWT-PE CTC#2
// -------------
// Define PE_DP_MERGE_OCCURRENCES to allow the Seed-
// Enhancement module to remove seemingly duplicated or very close
// occurrences returned by seed-alignment to speed up
// Dynamic Programming and reduce redundant output.
#define PE_DP_MERGE_OCCURRENCES
// -------------

void _PEPopulateMAPQCalculator(PEArguments * peArguments);

unsigned long long PEProcessReadDoubleStrand(PEArguments * peArguments, 
                            SRAModelSet * sraModelSet, SRAModelSet * sraModelSet_Seed, SRAModelSet * sraModelSet_Extend) {
    unsigned long long peCount = 0;
    unsigned long long runSaCount = 0;
    unsigned long long runSaCount_mate = 0;
    unsigned long long saCount = 0;
    unsigned long long saCount_mate = 0;
    unsigned long long peRunOccCount=0;
    int isAligned = 0;
    int isAligned_mate = 0;
    int isAligned_pair = 0;
    
    PEInput * peInput = peArguments->PEAlgnmtInput;
    int OutputType = peInput->OutputType;
    int i,j,k;
    uint8_t stageIdx        = 0;
    uint8_t stageIdx_mate   = 0;
    
    SRAArguments * sraQueryInputPos = peArguments->sraArgsPos;
    SRAArguments * sraQueryInputNeg = peArguments->sraArgsNeg;
    SRAArguments * sraQueryInputPos_mate = peArguments->sraArgsPos_mate;
    SRAArguments * sraQueryInputNeg_mate = peArguments->sraArgsNeg_mate;
    
    SRASetting * sraAlgnmtSetting = sraQueryInputPos->AlgnmtSetting;
    
    SRAModel * sraModel = SRAModelSetGetModel(sraModelSet,sraQueryInputPos->QueryInfo->ReadLength,QUERY_POS_STRAND);
    SRAModel * sraModel_mate = SRAModelSetGetModel(sraModelSet,sraQueryInputPos_mate->QueryInfo->ReadLength,QUERY_POS_STRAND);
    
    SRAModel * sraModel_neg = SRAModelSetGetModel(sraModelSet,sraQueryInputPos->QueryInfo->ReadLength,QUERY_NEG_STRAND);
    SRAModel * sraModel_mate_neg = SRAModelSetGetModel(sraModelSet,sraQueryInputPos_mate->QueryInfo->ReadLength,QUERY_NEG_STRAND);
    
    DPArguments * dpArguments = peArguments->dpArguments;
    DPWork * dpWork = dpArguments->dpWork;
    SRAIndex * sraIndex = peArguments->sraArgsPos->AlgnmtIndex;
    
    PEOutput * peOutput = peArguments->PEAlgnmtOutput;
    PEStats * peStats = peArguments->PEAlgnmtStats;
    SRAOutput * peOutput_OccHandler = peArguments->PEAlgnmtOutput_OccHandle;
    
    PEDPSetting * pedpSetting = peArguments->pedpSetting;
    
    PEOutputInitialise(peOutput);
    PESTATInitialise(peStats);
    MAPQCalculatorInitialise(peArguments->MapqCalc);
    
    OCCReportDelimitor(sraQueryInputPos);
    OCCReportDelimitor(sraQueryInputPos_mate);
    PEOCCReportDelimitor(peArguments);
    
    SRAWorkingMemoryInitialise(sraQueryInputPos->AlgnmtMemory);
    SRAWorkingMemoryInitialise(sraQueryInputPos_mate->AlgnmtMemory);
    
    unsigned int occCounts_1[MAX_NUM_OF_MISMATCH+1];
    unsigned int occCounts_2[MAX_NUM_OF_MISMATCH+1];
    SRAOccurrence * occLists_1[MAX_NUM_OF_MISMATCH+1];
    SRAOccurrence * occLists_2[MAX_NUM_OF_MISMATCH+1];
    
    for (i=0;i<MAX_NUM_OF_MISMATCH+1;i++) {
        occCounts_1[i]=0;
        occCounts_2[i]=0;
        occLists_1[i]=NULL;
        occLists_2[i]=NULL;
    }
    
    int InputSGAOrphanTriggerThreshold = (int) (pedpSetting->SGAOrphanTriggerTF*(double)sraQueryInputPos->QueryInfo->ReadLength);
    int InputSGAOrphanTriggerThreshold_mate = (int) (pedpSetting->SGAOrphanTriggerTF*(double)sraQueryInputPos_mate->QueryInfo->ReadLength);
    
    int InputSGAOrphanExtendTriggerThreshold = (int) (pedpSetting->SGAOrphanExtendTriggerTF*(double)sraQueryInputPos->QueryInfo->ReadLength);
    int InputSGAOrphanExtendTriggerThreshold_mate = (int) (pedpSetting->SGAOrphanExtendTriggerTF*(double)sraQueryInputPos_mate->QueryInfo->ReadLength);
    
    j=0;
    while (1) {

        //////////////////////////////////////////////////////////////////////////
        // Read from first input file
        //////////////////////////////////////////////////////////////////////////
            k=j;                                                                //
            runSaCount  = ProcessStage(sraQueryInputPos,sraModel,&k);           //
            k=j;                                                                //
            runSaCount += ProcessStage(sraQueryInputNeg,sraModel_neg,&k);       //
            saCount += runSaCount;                                              //
            isAligned |= runSaCount>0;                                          //
        //////////////////////////////////////////////////////////////////////////
        
        
        //////////////////////////////////////////////////////////////////////////
        // Read from second input file
        //////////////////////////////////////////////////////////////////////////
            k=j;                                                                //
            runSaCount_mate  = ProcessStage(sraQueryInputPos_mate,sraModel_mate,&k); //
            k=j;                                                                //
            runSaCount_mate += ProcessStage(sraQueryInputNeg_mate,sraModel_mate_neg,&k); //
            saCount_mate += runSaCount_mate;                                    //
            isAligned_mate |= runSaCount_mate>0;                                //
        //////////////////////////////////////////////////////////////////////////
        
        
        if (OutputType==PE_REPORT_RANDOM_BEST) {

            if (runSaCount) {
                occCounts_1[stageIdx] = SRAOCCCountOccurrences(sraQueryInputPos->AlgnmtOutput->occCollector);
                occLists_1[stageIdx] = (SRAOccurrence *) malloc( sizeof(SRAOccurrence) * occCounts_1[stageIdx] );
                occCounts_1[stageIdx] = SRAOCCPopulateSRAOccList(sraQueryInputPos->AlgnmtOutput->occCollector, occLists_1[stageIdx]);
            }
            SRAOCCInitialise(sraQueryInputPos->AlgnmtOutput->occCollector);
            stageIdx++;
            
            if (runSaCount_mate) {
                occCounts_2[stageIdx_mate] = SRAOCCCountOccurrences(sraQueryInputPos_mate->AlgnmtOutput->occCollector);
                occLists_2[stageIdx_mate] = (SRAOccurrence *) malloc( sizeof(SRAOccurrence) * occCounts_2[stageIdx_mate] );
                occCounts_2[stageIdx_mate] = SRAOCCPopulateSRAOccList(sraQueryInputPos_mate->AlgnmtOutput->occCollector, occLists_2[stageIdx_mate]);
            }
            SRAOCCInitialise(sraQueryInputPos_mate->AlgnmtOutput->occCollector);
            stageIdx_mate++;
            

            i = 0;
            while (1) {
            
                if (i<stageIdx_mate && occCounts_1[stageIdx-1]>0 && occCounts_2[i]>0) {
                    PEMappingOccurrences(peInput,peOutput,peStats,occLists_1[stageIdx-1],occCounts_1[stageIdx-1],occLists_2[i],occCounts_2[i],PE_INPUTOCC_UNSORTED);
                    if (peOutput->flag!=PE_ALIGNMENT_COMPLETED) {
                        printf("[PairEnd] Module terminated abnormally. errCode = %u\n",peOutput->flag);
                    }
                }
                
                ////////////////////////////////////////////////////////////
                // Reporting Pairs
                //////////////////////////////////////////////////////////////////////////////////////
                if (PEAligned(peOutput)) {                                                          //
                    peRunOccCount = PEOCCDumpAlignments(peArguments);                               //
                    isAligned_pair=1;                                                               //
                    peCount += peRunOccCount;                                                       //
                    break;                                                                          //
                }                                                                                   //
                //////////////////////////////////////////////////////////////////////////////////////
                
                if (i<stageIdx && occCounts_2[stageIdx_mate-1]>0 && occCounts_1[i]>0) {
                    PEMappingOccurrences(peInput,peOutput,peStats,occLists_1[i],occCounts_1[i],occLists_2[stageIdx_mate-1],occCounts_2[stageIdx_mate-1],PE_INPUTOCC_UNSORTED);
                    if (peOutput->flag!=PE_ALIGNMENT_COMPLETED) {
                        printf("[PairEnd] Module terminated abnormally. errCode = %u\n",peOutput->flag);
                    }
                }
                
                ////////////////////////////////////////////////////////////
                // Reporting Pairs
                //////////////////////////////////////////////////////////////////////////////////////
                if (PEAligned(peOutput)) {                                                          //
                    peRunOccCount = PEOCCDumpAlignments(peArguments);                               //
                    isAligned_pair=1;                                                               //
                    peCount += peRunOccCount;                                                       //
                    break;                                                                          //
                }                                                                                   //
                //////////////////////////////////////////////////////////////////////////////////////
                
                i++;
                if (i>=stageIdx_mate && i>=stageIdx) break;
            }
            
            if (isAligned_pair) break;
            
        }
        
        if (OutputType==PE_REPORT_ALL_BEST) {

            if (runSaCount) {
                occCounts_1[stageIdx] = SRAOCCCountOccurrences(sraQueryInputPos->AlgnmtOutput->occCollector);
                occLists_1[stageIdx] = (SRAOccurrence *) malloc( sizeof(SRAOccurrence) * occCounts_1[stageIdx] );
                occCounts_1[stageIdx] = SRAOCCPopulateSRAOccList(sraQueryInputPos->AlgnmtOutput->occCollector, occLists_1[stageIdx]);
            }
            SRAOCCInitialise(sraQueryInputPos->AlgnmtOutput->occCollector);
            stageIdx++;
            
            if (runSaCount_mate) {
                occCounts_2[stageIdx_mate] = SRAOCCCountOccurrences(sraQueryInputPos_mate->AlgnmtOutput->occCollector);
                occLists_2[stageIdx_mate] = (SRAOccurrence *) malloc( sizeof(SRAOccurrence) * occCounts_2[stageIdx_mate] );
                occCounts_2[stageIdx_mate] = SRAOCCPopulateSRAOccList(sraQueryInputPos_mate->AlgnmtOutput->occCollector, occLists_2[stageIdx_mate]);
            }
            SRAOCCInitialise(sraQueryInputPos_mate->AlgnmtOutput->occCollector);
            stageIdx_mate++;
            

            i = 0;
            while (1) {
            
                if (i<stageIdx_mate && occCounts_1[stageIdx-1]>0 && occCounts_2[i]>0) {
                    PEMappingOccurrences(peInput,peOutput,peStats,occLists_1[stageIdx-1],occCounts_1[stageIdx-1],occLists_2[i],occCounts_2[i],PE_INPUTOCC_UNSORTED);
                    if (peOutput->flag!=PE_ALIGNMENT_COMPLETED) {
                        printf("[PairEnd] Module terminated abnormally. errCode = %u\n",peOutput->flag);
                    }
                }
                
                ////////////////////////////////////////////////////////////
                // Reporting Pairs
                //////////////////////////////////////////////////////////////////////////////////////
                if (PEAligned(peOutput)) {                                                          //
                    peRunOccCount = PEOCCDumpAlignments(peArguments);                               //
                    isAligned_pair=1;                                                               //
                    peCount += peRunOccCount;                                                       //
                }                                                                                   //
                //////////////////////////////////////////////////////////////////////////////////////
                
                // In order to avoid duplicated PE.
                // 0-mismatch vs 0-mismatch in case 0; 1-m vs 1-m in case 1; etc..
                // The mate processing can be skipped if it's the last round.
                if (!(i==stageIdx-1 && i==stageIdx_mate-1)) {
                
                    if (i<stageIdx && occCounts_2[stageIdx_mate-1]>0 && occCounts_1[i]>0) {
                        PEMappingOccurrences(peInput,peOutput,peStats,occLists_1[i],occCounts_1[i],occLists_2[stageIdx_mate-1],occCounts_2[stageIdx_mate-1],PE_INPUTOCC_UNSORTED);
                        if (peOutput->flag!=PE_ALIGNMENT_COMPLETED) {
                            printf("[PairEnd] Module terminated abnormally. errCode = %u\n",peOutput->flag);
                        }
                    }
                    
                    ////////////////////////////////////////////////////////////
                    // Reporting Pairs
                    //////////////////////////////////////////////////////////////////////////////////////
                    if (PEAligned(peOutput)) {                                                          //
                        peRunOccCount = PEOCCDumpAlignments(peArguments);                               //
                        isAligned_pair=1;                                                               //
                        peCount += peRunOccCount;                                                       //
                    }                                                                                   //
                    //////////////////////////////////////////////////////////////////////////////////////
                }
                
                i++;
                if (isAligned_pair) break;
                if (i>=stageIdx_mate && i>=stageIdx) break;
            }
            
            if (isAligned_pair) {

                ///////////////////////////////////////////////////////////////////
                // SRA-MAPQ Computation
                ///////////////////////////////////////////////////////////////////////////////////////////
                SRAModel * sraModelExtend = SRAModelSetGetModel(sraModelSet_Extend,sraQueryInputPos->QueryInfo->ReadLength,QUERY_POS_STRAND);
                SRAModel * sraModelExtend_mate = SRAModelSetGetModel(sraModelSet_Extend,sraQueryInputPos_mate->QueryInfo->ReadLength,QUERY_POS_STRAND);
                SRAModel * sraModelExtendNeg = SRAModelSetGetModel(sraModelSet_Extend,sraQueryInputPos->QueryInfo->ReadLength,QUERY_NEG_STRAND);
                SRAModel * sraModelExtendNeg_mate = SRAModelSetGetModel(sraModelSet_Extend,sraQueryInputPos_mate->QueryInfo->ReadLength,QUERY_NEG_STRAND);
                
                SRAPopulateMAPQCalculator(sraQueryInputPos,sraQueryInputNeg,dpArguments->dpScores,sraModelExtend,sraModelExtendNeg);
                SRAPopulateMAPQCalculator(sraQueryInputPos_mate,sraQueryInputNeg_mate,dpArguments->dpScores,sraModelExtend_mate,sraModelExtendNeg_mate);

                break;
            }
            
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
    
    //Dynamic programming enhancement to random best alignment
    //Dynamic programming enhancement to all best alignment
    if (OutputType==PE_REPORT_RANDOM_BEST || 
        OutputType==PE_REPORT_ALL_BEST) {
    
        /////////////////////////////////////////////////
        // Determine whether we need to run round Alpha and Beta
        int roundAlphaNeeded = 0;
        int roundBetaNeeded = 0;
        
        if (!isAligned_pair) {
            if (pedpSetting->SGAOrphanEnhancement) {
                if (isAligned && !isAligned_mate && saCount<=InputSGAOrphanTriggerThreshold) {
                    roundAlphaNeeded = 1;
                } else if (!isAligned && isAligned_mate && saCount_mate<=InputSGAOrphanTriggerThreshold_mate) {
                    roundBetaNeeded = 1;
                }
            }
            if (pedpSetting->SGAOrphanExtendEnhancement) {
                if (isAligned && isAligned_mate && 
                    saCount<=InputSGAOrphanExtendTriggerThreshold &&
                    saCount_mate<=InputSGAOrphanExtendTriggerThreshold_mate ) {
                    
                    roundAlphaNeeded = 1;
                    roundBetaNeeded = 1;
                }
            }
        }
        
        /////////////////////////////////////////////////
        // Round Alpha
        // Matching the alignments from the first read file (result from all rounds) by dynamic programming
        if (roundAlphaNeeded && !isAligned_pair) {
            for (i=0;i<stageIdx;i++) {
                if (occCounts_1[i]>0) {
                    peRunOccCount = DPPEOrphanMappingOccurrences(peArguments,
                                                            occLists_1[i],occCounts_1[i],
                                                            sraQueryInputPos->QueryInfo->ReadLength,
                                                            sraQueryInputPos_mate->QueryInfo->ReadCode,
                                                            sraQueryInputNeg_mate->QueryInfo->ReadCode,
                                                            sraQueryInputPos_mate->QueryInfo->ReadLength,
                                                            sraQueryInputPos_mate->MapqCalc,
                                                            SRAOCC_TYPE_PE_DP_BASE_READ);
                    
                    ////////////////////////////////////////////////////////////
                    // Reporting Pairs
                    //////////////////////////////////////////////////////////////
                    if (peRunOccCount) {                                        //
                        isAligned_pair=1;                                       //
                        peCount += peRunOccCount;                               //
                        break;                                                  //
                    }                                                           //
                    //////////////////////////////////////////////////////////////
                }
            }
        }
        /////////////////////////////////////////////////
        // Round Beta
        // Matching the alignments from the second read file (result from all rounds) by dynamic programming
        if (roundBetaNeeded && !isAligned_pair) {
            for (i=0;i<stageIdx_mate;i++) {
                if (occCounts_2[i]>0) {
                    peRunOccCount = DPPEOrphanMappingOccurrences(peArguments,
                                                            occLists_2[i],occCounts_2[i],
                                                            sraQueryInputPos_mate->QueryInfo->ReadLength,
                                                            sraQueryInputPos->QueryInfo->ReadCode,
                                                            sraQueryInputNeg->QueryInfo->ReadCode,
                                                            sraQueryInputPos->QueryInfo->ReadLength,
                                                            sraQueryInputPos->MapqCalc,
                                                            SRAOCC_TYPE_PE_DP_BASE_MATE);
                    
                    ////////////////////////////////////////////////////////////
                    // Reporting Pairs
                    //////////////////////////////////////////////////////////////
                    if (peRunOccCount) {                                        //
                        isAligned_pair=1;                                       //
                        peCount += peRunOccCount;                               //
                        break;                                                  //
                    }                                                           //
                    //////////////////////////////////////////////////////////////
                }
            }
        }
    }
    
    if (OutputType==PE_REPORT_ALL) {
        occCounts_1[stageIdx] = SRAOCCCountOccurrences(sraQueryInputPos->AlgnmtOutput->occCollector);
        occLists_1[stageIdx] = (SRAOccurrence *) malloc( sizeof(SRAOccurrence) * occCounts_1[stageIdx] );
        occCounts_1[stageIdx] = SRAOCCPopulateSRAOccList(sraQueryInputPos->AlgnmtOutput->occCollector, occLists_1[stageIdx]);
        stageIdx++;
    
        occCounts_2[stageIdx_mate] = SRAOCCCountOccurrences(sraQueryInputPos_mate->AlgnmtOutput->occCollector);
        occLists_2[stageIdx_mate] = (SRAOccurrence *) malloc( sizeof(SRAOccurrence) * occCounts_2[stageIdx_mate] );
        occCounts_2[stageIdx_mate] = SRAOCCPopulateSRAOccList(sraQueryInputPos_mate->AlgnmtOutput->occCollector, occLists_2[stageIdx_mate]);
        stageIdx_mate++;
        
        if (occCounts_1[stageIdx-1] > 0 && occCounts_2[stageIdx_mate-1] > 0) {
            //Calling Normal Core
            PEMappingOccurrences(peInput,peOutput,peStats,occLists_1[stageIdx-1],occCounts_1[stageIdx-1],occLists_2[stageIdx_mate-1],occCounts_2[stageIdx_mate-1],PE_INPUTOCC_UNSORTED);
            peRunOccCount = PEOCCDumpAlignments(peArguments);
            
            if (peOutput->flag!=PE_ALIGNMENT_COMPLETED) {
                printf("[PairEnd] Module terminated abnormally. errCode = %u\n",peOutput->flag);
            }
            
            if (peRunOccCount) {
                isAligned_pair=1;
                peCount += peRunOccCount;
            }
        }
        
        ///////////////////////////////////////////////////////////////////
        // SRA-MAPQ Computation
        ///////////////////////////////////////////////////////////////////////////////////////////
        if (isAligned_pair) {
            SRAModel * sraModelExtend = SRAModelSetGetModel(sraModelSet_Extend,sraQueryInputPos->QueryInfo->ReadLength,QUERY_POS_STRAND);
            SRAModel * sraModelExtend_mate = SRAModelSetGetModel(sraModelSet_Extend,sraQueryInputPos_mate->QueryInfo->ReadLength,QUERY_POS_STRAND);
            SRAModel * sraModelExtendNeg = SRAModelSetGetModel(sraModelSet_Extend,sraQueryInputPos->QueryInfo->ReadLength,QUERY_NEG_STRAND);
            SRAModel * sraModelExtendNeg_mate = SRAModelSetGetModel(sraModelSet_Extend,sraQueryInputPos_mate->QueryInfo->ReadLength,QUERY_NEG_STRAND);
            
            SRAPopulateMAPQCalculator(sraQueryInputPos,sraQueryInputNeg,dpArguments->dpScores,sraModelExtend,sraModelExtendNeg);
            SRAPopulateMAPQCalculator(sraQueryInputPos_mate,sraQueryInputNeg_mate,dpArguments->dpScores,sraModelExtend_mate,sraModelExtendNeg_mate);
        }

        /////////////////////////////////////////////////
        // Determine whether we need to run round Alpha and Beta
        int roundAlphaNeeded = 0;
        int roundBetaNeeded = 0;
        
        if (!isAligned_pair) {
            if (pedpSetting->SGAOrphanEnhancement && 
                occCounts_1[stageIdx-1] == 0 && 
                occCounts_2[stageIdx_mate-1] <= InputSGAOrphanTriggerThreshold_mate) {
                
                roundAlphaNeeded = 1;
            } else if (pedpSetting->SGAOrphanEnhancement && 
                occCounts_2[stageIdx_mate-1] == 0 && 
                occCounts_1[stageIdx-1] <= InputSGAOrphanTriggerThreshold) {
                
                roundBetaNeeded = 1;
            }

            if (pedpSetting->SGAOrphanExtendEnhancement && 
                occCounts_1[stageIdx-1] > 0 && 
                occCounts_2[stageIdx_mate-1] > 0 && 
                occCounts_1[stageIdx-1] <= InputSGAOrphanExtendTriggerThreshold &&
                occCounts_2[stageIdx_mate-1] <= InputSGAOrphanExtendTriggerThreshold_mate) {
                
                roundAlphaNeeded = 1;
                roundBetaNeeded = 1;
            }
        }
        
        /////////////////////////////////////////////////
        // Round Alpha
        // Matching the alignments from the first read file (result from all rounds) by dynamic programming
        if (roundAlphaNeeded) {
            //Calling SGA Core Case A
            peRunOccCount = DPPEOrphanMappingOccurrences(peArguments,occLists_2[stageIdx_mate-1],occCounts_2[stageIdx_mate-1],sraQueryInputPos_mate->QueryInfo->ReadLength,
                                    sraQueryInputPos->QueryInfo->ReadCode,sraQueryInputNeg->QueryInfo->ReadCode,sraQueryInputPos->QueryInfo->ReadLength,
                                    sraQueryInputPos->MapqCalc,
                                    SRAOCC_TYPE_PE_DP_BASE_MATE);
            
            if (peRunOccCount) {
                isAligned_pair=1;
                peCount += peRunOccCount;
            }
        }

        /////////////////////////////////////////////////
        // Round Beta
        // Matching the alignments from the second read file (result from all rounds) by dynamic programming
        if (roundBetaNeeded) {
            //Calling SGA Core Case B
            peRunOccCount = DPPEOrphanMappingOccurrences(peArguments,occLists_1[stageIdx-1],occCounts_1[stageIdx-1],sraQueryInputPos->QueryInfo->ReadLength,
                                    sraQueryInputPos_mate->QueryInfo->ReadCode,sraQueryInputNeg_mate->QueryInfo->ReadCode,sraQueryInputPos_mate->QueryInfo->ReadLength,
                                    sraQueryInputPos_mate->MapqCalc,
                                    SRAOCC_TYPE_PE_DP_BASE_READ);
            
            if (peRunOccCount) {
                isAligned_pair=1;
                peCount += peRunOccCount;
            }
        }
        
    } 

    /////////////////////////////////////////////////
    // Round Gamma-1
    // Attempt to produce Pair-end alignment result by seeding the reads for candidate
    if (pedpSetting->SGASeedEnhancement && !isAligned_pair) {
    
        peRunOccCount = DPPESeedRecoverOccurrences(peArguments, sraModelSet_Seed,
                                        pedpSetting->SGASeedLength, pedpSetting->SGASeedLengthOverlap,
                                        pedpSetting->SGASeedLooseCriteria, pedpSetting->SGASeedOccLimitation);

        if (peRunOccCount) {
            isAligned_pair=1;
            peCount += peRunOccCount;
        }
    }
    
    /////////////////////////////////////////////////
    // Round Gamma-2
    // Attempt to produce Pair-end alignment result by seeding the reads for candidate
    if (pedpSetting->SGASeedEnhancement_1 && !isAligned_pair) {
    
        peRunOccCount = DPPESeedRecoverOccurrences(peArguments, sraModelSet_Seed,
                                        pedpSetting->SGASeedLength_1, pedpSetting->SGASeedLengthOverlap_1,
                                        pedpSetting->SGASeedLooseCriteria_1, pedpSetting->SGASeedOccLimitation_1);

        if (peRunOccCount) {
            isAligned_pair=1;
            peCount += peRunOccCount;
        }
    }
    
    ///////////////////////////////////////////////////////////////////
    // Reporting output line that does not involve a valid pair
    // this part is completely based on the assumption
    // feel free to edit it to correct the output.
    // For any pair of read that does not produce a valid pair-end alignment..
    ///////////////////////////////////////////////////////////////////////////////////////////
    if (!isAligned_pair) {
        /////////////////////////////////////////////////
        //Report a sub-valid pair-end if both read aligned to something
        if (isAligned && isAligned_mate) {
            unsigned int nonEmptyList_1, nonEmptyList_2;
            for (i=0;i<stageIdx;i++) {
                if (occLists_1[i]!=NULL) nonEmptyList_1 = i;
            }
            for (i=0;i<stageIdx_mate;i++) {
                if (occLists_2[i]!=NULL) nonEmptyList_2 = i;
            }
            PEOCCReportSubValidAlignment(peArguments, occLists_1[nonEmptyList_1], occLists_2[nonEmptyList_2]);
            
        /////////////////////////////////////////////////
        //Report that both read does not produce any short read alignment
        } else if (!isAligned && !isAligned_mate) {
            PEOCCReportNoAlignment(peArguments,SRAOCC_TYPE_PE_READ_BOTH_NO_ALIGNMENT);
            PEOCCReportNoAlignment(peArguments,SRAOCC_TYPE_PE_MATE_BOTH_NO_ALIGNMENT);
            
        /////////////////////////////////////////////////
        //For the cases that one of the read produce some short read alignment
        // while the other does not.
        } else {
            /////////////////////////////////////////////////
            //The first of pair
            if (!isAligned) {
                for (i=0;i<stageIdx_mate;i++) {
                    // Output only the first maxResult of them if defined
                    unsigned int occCount = occCounts_2[i];
                    if (peInput->maxResult != -1 && 
                        occCount > peInput->maxResult) {
                        occCount = peInput->maxResult;
                    }
                    if (occLists_2[i]!=NULL) PEOCCDumpSingleAlignments(peArguments,occLists_2[i],occCount,SRAOCC_TYPE_PE_READ_READ_NO_ALIGNMENT);
                }
            } else {
                for (i=0;i<stageIdx;i++) {
                    // Output only the first maxResult of them if defined
                    unsigned int occCount = occCounts_1[i];
                    if (peInput->maxResult != -1 && 
                        occCount > peInput->maxResult) {
                        occCount = peInput->maxResult;
                    }
                    if (occLists_1[i]!=NULL) PEOCCDumpSingleAlignments(peArguments,occLists_1[i],occCount,SRAOCC_TYPE_PE_READ_MATE_NO_ALIGNMENT);
                }
            }
            
            /////////////////////////////////////////////////
            //The second of pair
            if (!isAligned_mate) {
                for (i=0;i<stageIdx;i++) {
                    // Output only the first maxResult of them if defined
                    unsigned int occCount = occCounts_1[i];
                    if (peInput->maxResult != -1 && 
                        occCount > peInput->maxResult) {
                        occCount = peInput->maxResult;
                    }
                    if (occLists_1[i]!=NULL) PEOCCDumpSingleAlignments(peArguments,occLists_1[i],occCount,SRAOCC_TYPE_PE_MATE_READ_NO_ALIGNMENT);
                }
            } else {
                for (i=0;i<stageIdx_mate;i++) {
                    // Output only the first maxResult of them if defined
                    unsigned int occCount = occCounts_2[i];
                    if (peInput->maxResult != -1 && 
                        occCount > peInput->maxResult) {
                        occCount = peInput->maxResult;
                    }
                    if (occLists_2[i]!=NULL) PEOCCDumpSingleAlignments(peArguments,occLists_2[i],occCount,SRAOCC_TYPE_PE_MATE_MATE_NO_ALIGNMENT);
                }
            }
        }
    }
    ///////////////////////////////////////////////////////////////////
    // PE-MAPQ Computation
    ///////////////////////////////////////////////////////////////////////////////////////////
    if (isAligned_pair) {
        if (peArguments->MapqCalc->status != MAPQ_CALC_FIGURE_READY) {
            _PEPopulateMAPQCalculator(peArguments);
        }
        
        int value0,value1;
        
        #ifdef PE_DEBUG_PRINT_MAPQ_CALCULATION
            MAPQPrintParameters(sraQueryInputPos->MapqCalc);
            MAPQPrintParameters(sraQueryInputPos_mate->MapqCalc);
            MAPQPrintParameters(peArguments->MapqCalc);
        #endif
        
        MAPQPEGetScore(peArguments,&value0,&value1);
        PEOCCMAPQValue(peArguments,value0,value1);
        
        #ifdef PE_DEBUG_PRINT_MAPQ_CALCULATION
            printf("Value0 = %u             Value1 = %u\n",value0,value1);
        
            if (peArguments->MapqCalc->status != MAPQ_CALC_FIGURE_READY ||
                sraQueryInputPos->MapqCalc->status != MAPQ_CALC_FIGURE_READY ||
                sraQueryInputPos_mate->MapqCalc->status != MAPQ_CALC_FIGURE_READY) {
                printf("DEBUG FAILED %s\n",sraQueryInputPos->QueryInfo->ReadName);
            }
        #endif
    }
    

    if (PE_DUMP_SRA_OUTPUT==1 && OutputType==PE_REPORT_ALL) {
        OCCFlushCache(sraQueryInputPos);
        OCCFlushCache(sraQueryInputPos_mate);
    } else {
        SRAOCCInitialise(sraQueryInputPos->AlgnmtOutput->occCollector);
        SRAOCCInitialise(sraQueryInputPos_mate->AlgnmtOutput->occCollector);
    }

    PEOCCFlushCache(peArguments);

    for (i=0;i<stageIdx;i++) {
        if (occLists_1[i]!=NULL) free(occLists_1[i]);
    }
    for (i=0;i<stageIdx_mate;i++) {
        if (occLists_2[i]!=NULL) free(occLists_2[i]);
    }

    return peCount;
}


unsigned int DPPESeedRecoverOccurrences(PEArguments * peArguments,
                                        SRAModelSet * sraModelSet_Seed, 
                                        double SGASeedLength, double SGASeedLengthOverlap,
                                        int SGASeedLooseCriteria, int SGASeedOccLimitation) {
    unsigned long long peCount = 0;
    int i;

    SRAArguments * sraQueryInputPos = peArguments->sraArgsPos;
    SRAArguments * sraQueryInputNeg = peArguments->sraArgsNeg;
    SRAArguments * sraQueryInputPos_mate = peArguments->sraArgsPos_mate;
    SRAArguments * sraQueryInputNeg_mate = peArguments->sraArgsNeg_mate;
    
    DPArguments * dpArguments = peArguments->dpArguments;
    DPWork * dpWork = dpArguments->dpWork;
    SRAIndex * sraIndex = peArguments->sraArgsPos->AlgnmtIndex;
    
    PEInput * peInput = peArguments->PEAlgnmtInput;
    PEOutput * peOutput = peArguments->PEAlgnmtOutput;
    
    PEDPSetting * pedpSetting = peArguments->pedpSetting;
    
    //1. For the read, perform SRA on the seeds. Store the result in sraQueryInputPos->AlgnmtOutput->occCollector
    //   tweak the occurrences to store the anticipated occurrence position.
    //2. For the mate, perform SRA on the seeds. Store the result in sraQueryInputPos_mate->AlgnmtOutput->occCollector
    //   tweak the occurrences to store the anticipated occurrence position.
    //3. Check the number of occ reported by either read/mate should not exceed SGASeedOccLimitation.
    //4. Do a PE on the occurrences to find the candidate. The result will be stored in peArguments->PEAlgnmtOutput_OccHandle.
    //4a.We merge the results that are too close to each other. Too close is defined as distance within 25% of the read length.
    //5. For each pair of result it indicates that it could be a possible occurrence. We then verify the candidate by dynamic programming.
    //6. Post-verification validation to make sure the occurrence falls between the allowed insertion size.
    //7. Report the occurrences

    #ifdef DP_DEBUG_PRINT_SEED_MATCHING
        printf("\n\n[DPPESeedRecoverOccurrences] Seed Recover Module Initiated.\n");
    #endif

    // Setting up the seed-lengths
    int seedLength,      seedLength_mate;
    int seedUniqeLength, seedUniqeLength_mate;
    
    // If seedLength > 1, treat input as constant seedLength or seedOverlap
    // If seedLength <= 1, treat input as a ratio to readLen set Length to zero when readLength or mateLength = 0 as a ratio to readLen; + 0.5 to do rounding
    seedLength           =    ( SGASeedLength > 1.0 ? SGASeedLength : (double)sraQueryInputPos->QueryInfo->ReadLength * SGASeedLength ) + 0.5;
    seedLength_mate      =    ( SGASeedLength > 1.0 ? SGASeedLength : (double)sraQueryInputPos_mate->QueryInfo->ReadLength * SGASeedLength ) + 0.5;
    if (seedLength < SRA_MIN_READ_LENGTH) { seedLength = SRA_MIN_READ_LENGTH; }
    else if (seedLength > SRA_MAX_READ_LENGTH) { seedLength = SRA_MAX_READ_LENGTH; }
    
    if (seedLength_mate < SRA_MIN_READ_LENGTH) { seedLength_mate = SRA_MIN_READ_LENGTH; }
    else if (seedLength_mate > SRA_MAX_READ_LENGTH) { seedLength_mate = SRA_MAX_READ_LENGTH; }
                
    // If seedOverlap > 1, treat input as constant seedLength or seedOverlap
    // If seedOverlap <= 1, treat input as a ratio to readLen; + 0.5 to do rounding
    seedUniqeLength      =   ( seedLength      - ( SGASeedLengthOverlap > 1.0 ? SGASeedLengthOverlap : (double)sraQueryInputPos->QueryInfo->ReadLength      * SGASeedLengthOverlap ) ) + 0.5;
    seedUniqeLength_mate =   ( seedLength_mate - ( SGASeedLengthOverlap > 1.0 ? SGASeedLengthOverlap : (double)sraQueryInputPos_mate->QueryInfo->ReadLength * SGASeedLengthOverlap ) ) + 0.5;
    
    //Step 1
    //=========================================================================
    
    SRAArguments * seedSraArgPos = SRAARGMakeMate(sraQueryInputPos);
    SRAArguments * seedSraArgNeg = SRAARGMakeSlave(seedSraArgPos);
    SRAArguments * seedSraArgPos_mate = SRAARGMakeMate(sraQueryInputPos_mate);
    SRAArguments * seedSraArgNeg_mate = SRAARGMakeSlave(seedSraArgPos_mate);
    //Setting the QueryInfo structures
    SRAQueryInfoCopy(sraQueryInputPos->QueryInfo, seedSraArgPos->QueryInfo);
    SRAQueryInfoCopy(sraQueryInputNeg->QueryInfo, seedSraArgNeg->QueryInfo);
    SRAQueryInfoCopy(sraQueryInputPos_mate->QueryInfo, seedSraArgPos_mate->QueryInfo);
    SRAQueryInfoCopy(sraQueryInputNeg_mate->QueryInfo, seedSraArgNeg_mate->QueryInfo);

    #ifdef DP_DEBUG_PRINT_SEED_MATCHING
        DEBUGSRAQueryInfoPrint(seedSraArgPos->QueryInfo,"POS");
        DEBUGSRAQueryInfoPrint(seedSraArgNeg->QueryInfo,"NEG");
        DEBUGSRAQueryInfoPrint(seedSraArgPos_mate->QueryInfo,"MPOS");
        DEBUGSRAQueryInfoPrint(seedSraArgNeg_mate->QueryInfo,"MNEG");
    #endif
    
    //Setting the AlgnmtOutput structures
    SRAAlgnmtOutputInitNone(seedSraArgPos->AlgnmtOutput);
    SRAAlgnmtOutputInitNone(seedSraArgPos_mate->AlgnmtOutput);
    
    // Setting the outputType be report all-best
    seedSraArgPos->AlgnmtSetting = SRASettingMakeClone(seedSraArgPos->AlgnmtSetting);
    seedSraArgPos->AlgnmtSetting->OutputType = SRA_REPORT_ALL_BEST;
    seedSraArgPos->AlgnmtSetting->MaxResult  = SGASeedOccLimitation;
    
    seedSraArgNeg->AlgnmtSetting = seedSraArgPos->AlgnmtSetting;
    seedSraArgPos_mate->AlgnmtSetting = seedSraArgPos->AlgnmtSetting;
    seedSraArgNeg_mate->AlgnmtSetting = seedSraArgPos->AlgnmtSetting;
    
    //Step 1
    //=========================================================================
    
    int properlySeeded = 1;
    
    // Re-using the SRA Occ Collector for seed alignment corrected ambPosition collection
    SRAOCCCollector * sraOccCollector = sraQueryInputPos->AlgnmtOutput->occCollector;
    SRAOCCCollector * sraOccCollector_mate = sraQueryInputPos_mate->AlgnmtOutput->occCollector;

    unsigned long long seedOccCounts_1 = SRAOCCCountOccurrencesByType(sraOccCollector     ,SRAOCC_TYPE_AWAIT_TRANSLATE);
    unsigned long long seedOccCounts_2 = SRAOCCCountOccurrencesByType(sraOccCollector_mate,SRAOCC_TYPE_AWAIT_TRANSLATE);
    
    /////////////////////////////////////////////////////////////////////
    // Process the seeds from READ
    /////////////////////////////////////////////////////////////////////
    
    if (properlySeeded && seedOccCounts_1==0) {
        #ifdef DP_DEBUG_PRINT_SEED_MATCHING
            printf("[SGASeedEnhancement] SRA Seed Alignment for READ.\n");
        #endif
        SRAOCCInitialise(sraOccCollector);
        SRAModel * sraModels_Seeding     = SRAModelSetGetModel(sraModelSet_Seed,seedLength,QUERY_POS_STRAND);
        SRAModel * sraModels_Seeding_neg = SRAModelSetGetModel(sraModelSet_Seed,seedLength,QUERY_POS_STRAND);
        properlySeeded = SRASeedProcessReadDbl(seedSraArgPos,seedSraArgNeg,
                          sraModels_Seeding,sraModels_Seeding_neg,
                          seedLength,seedUniqeLength,sraOccCollector);
        seedOccCounts_1 = SRAOCCCountOccurrencesByType(sraOccCollector     ,SRAOCC_TYPE_AWAIT_TRANSLATE);
    } else {
        //SAFE-GUARD
        //printf("[SGASeedEnhancement] Skipping the Seeding SRA alignment -- %llu alignments available.\n",seedOccCounts_1);
    }
    
    /////////////////////////////////////////////////////////////////////
    // Process the seeds from MATE
    /////////////////////////////////////////////////////////////////////
    
    if (properlySeeded && seedOccCounts_2==0) {
        #ifdef DP_DEBUG_PRINT_SEED_MATCHING
            printf("[SGASeedEnhancement] SRA Seed Alignment for MATE.\n");
        #endif
        SRAOCCInitialise(sraOccCollector_mate);
        SRAModel * sraModels_Seeding     = SRAModelSetGetModel(sraModelSet_Seed,seedLength_mate,QUERY_POS_STRAND);
        SRAModel * sraModels_Seeding_neg = SRAModelSetGetModel(sraModelSet_Seed,seedLength_mate,QUERY_POS_STRAND);
        properlySeeded = SRASeedProcessReadDbl(seedSraArgPos_mate,seedSraArgNeg_mate,
                          sraModels_Seeding,sraModels_Seeding_neg,
                          seedLength_mate,seedUniqeLength_mate,sraOccCollector_mate);
        seedOccCounts_2 = SRAOCCCountOccurrencesByType(sraOccCollector_mate,SRAOCC_TYPE_AWAIT_TRANSLATE);
    } else {
        //SAFE-GUARD
        //printf("[SGASeedEnhancement] Skipping the Seeding SRA alignment -- %llu alignments available.\n",seedOccCounts_2);
    }
    
    /////////////////////////////////////////////////////////////////////
    // Check the seeding completed successfully
    /////////////////////////////////////////////////////////////////////
    if (!properlySeeded) {
    
        #ifdef DP_DEBUG_PRINT_SEED_MATCHING
            printf("[SGASeedEnhancement] Read skipped due to limitation exceeded. 1=%u 2=%u b=%d\n",
                SRAOCCCountOccurrencesByType(sraOccCollector     ,SRAOCC_TYPE_AWAIT_TRANSLATE),
                SRAOCCCountOccurrencesByType(sraOccCollector_mate,SRAOCC_TYPE_AWAIT_TRANSLATE),
                SGASeedOccLimitation);
        #endif
        
    } else {
    
        /////////////////////////////////////////////////////////////////////
        // Pair-End Alignment On sraOccCollector and sraOccCollector_mate
        /////////////////////////////////////////////////////////////////////
        
        SRAOccurrence * seedOccLists_1 = (SRAOccurrence *) malloc( sizeof(SRAOccurrence) * seedOccCounts_1 );
        SRAOccurrence * seedOccLists_2 = (SRAOccurrence *) malloc( sizeof(SRAOccurrence) * seedOccCounts_2 );
        
        SRAOCCPopulateSRAOccList(sraOccCollector     ,seedOccLists_1);
        SRAOCCPopulateSRAOccList(sraOccCollector_mate,seedOccLists_2);
        
        #ifdef DP_DEBUG_PRINT_SEED_MATCHING
            printf("[SGASeedEnhancement] -- READ Seed Alignment gives the following result (%llu)\n",seedOccCounts_1);
            SRAOCCDebugPrint(sraOccCollector);
            printf("[SGASeedEnhancement] -- MATE Seed Alignment gives the following result (%llu)\n",seedOccCounts_2);
            SRAOCCDebugPrint(sraOccCollector_mate);
        #endif
        
        ////////////////////////////////////////////////////
        // Sort and Unique
        ////////////////////////////////////////////////////
        SRAOccurrencesSort(seedOccLists_1,&seedOccCounts_1);
        SRAOccurrencesSort(seedOccLists_2,&seedOccCounts_2);
        
        ////////////////////////////////////////////////////
        // Merging the occurrences that are too close to each other
        // This is a heuristic that help removing duplicated 
        // resulting pair-end alignment
        ////////////////////////////////////////////////////
        #ifdef PE_DP_MERGE_OCCURRENCES
            if (seedOccCounts_1>0) {
                unsigned long long mergedLastOccIdx = 0;
                unsigned long long mergedListIdx = 0;
                int mergeDist = (double)sraQueryInputPos->QueryInfo->ReadLength * 0.25f;
                
                SRAOccurrenceCopy(&(seedOccLists_1[0]),&(seedOccLists_1[mergedListIdx++]));
                for (i=1;i<seedOccCounts_1;i++) {
                    if (seedOccLists_1[i].ambPosition > seedOccLists_1[mergedLastOccIdx].ambPosition + mergeDist) {
                        SRAOccurrenceCopy(&(seedOccLists_1[i]),&(seedOccLists_1[mergedListIdx++]));
                        mergedLastOccIdx = i;
                    }
                }
                seedOccCounts_1 = mergedListIdx;
            }
            if (seedOccCounts_2>0) {
                unsigned long long mergedLastOccIdx = 0;
                unsigned long long mergedListIdx = 0;
                int mergeDist = (double)sraQueryInputPos_mate->QueryInfo->ReadLength * 0.25f;
                
                SRAOccurrenceCopy(&(seedOccLists_2[0]),&(seedOccLists_2[mergedListIdx++]));
                for (i=1;i<seedOccCounts_2;i++) {
                    if (seedOccLists_2[i].ambPosition > seedOccLists_2[mergedLastOccIdx].ambPosition + mergeDist) {
                        SRAOccurrenceCopy(&(seedOccLists_2[i]),&(seedOccLists_2[mergedListIdx++]));
                        mergedLastOccIdx = i;
                    }
                }
                seedOccCounts_2 = mergedListIdx;
            }
        #endif
        
        #ifdef DP_DEBUG_PRINT_SEED_MATCHING
            printf("[SGASeedEnhancement] -- READ Seed Alignment gives the following unique results (%llu)\n",seedOccCounts_1);
            SRAOccurrencesPrint(seedOccLists_1,seedOccCounts_1);
            printf("[SGASeedEnhancement] -- MATE Seed Alignment gives the following unique results (%llu)\n",seedOccCounts_2);
            SRAOccurrencesPrint(seedOccLists_2,seedOccCounts_2);
        #endif
            
        /////////////////////////////////////////////////////////////////////
        // Check the number of occurrences does not exceed the boundary
        /////////////////////////////////////////////////////////////////////
        if (SGASeedOccLimitation != -1 &&
           (seedOccCounts_1 > SGASeedOccLimitation ||
            seedOccCounts_2 > SGASeedOccLimitation)) {
            
            #ifdef DP_DEBUG_PRINT_SEED_MATCHING
                printf("[SGASeedEnhancement] Read skipped due to limitation exceeded. 1=%llu 2=%llu b=%d\n",
                    seedOccCounts_1,seedOccCounts_2,SGASeedOccLimitation);
            #endif
            
        } else {
            
            PEInput * seedPEInput = PEInputConstruct(sraIndex->bwt,sraIndex->hsp);
            PEStats * seedPEStats = PESTATConstruct();
            
            memcpy(seedPEInput,peInput,sizeof(PEInput));
            seedPEInput->OutputType = PE_REPORT_ALL;
            seedPEInput->maxResult = -1;
            
            PEMappingOccurrences(seedPEInput,peOutput,seedPEStats,seedOccLists_1,seedOccCounts_1,seedOccLists_2,seedOccCounts_2,PE_INPUTOCC_SORTED);
            if (peOutput->flag!=PE_ALIGNMENT_COMPLETED) {
                printf("[PairEnd] Module terminated abnormally. errCode = %u\n",peOutput->flag);
            }

            PEInputFree(seedPEInput);
            PESTATFree(seedPEStats);
            
            /////////////////////////////////////////////////////////////////////
            // SGA on candidate Pair-end alignment
            /////////////////////////////////////////////////////////////////////
            DPOccurrence dpOcc, dpOcc_mate;
            int dpAlgnmtOK = TRUE;
                            
            int seedScore, seedScore_mate;
            PEPairList * myNode = peOutput->root;
            
            seedScore = pedpSetting->SGAScoreTF * sraQueryInputPos->QueryInfo->ReadLength;
            seedScore_mate = pedpSetting->SGAScoreTF * sraQueryInputPos_mate->QueryInfo->ReadLength;

            MAPQCalculator mapqCalcOverlay_read;
            MAPQCalculator mapqCalcOverlay_mate;
            MAPQCalculatorInitialise(&mapqCalcOverlay_read);
            MAPQCalculatorInitialise(&mapqCalcOverlay_mate);

            OCCReportDelimitor(sraQueryInputPos);
            OCCReportDelimitor(sraQueryInputPos_mate);
            while (myNode!=NULL && dpAlgnmtOK==TRUE) {
            
                for (i=0;i<myNode->pairsCount;i++) {
                
                            
                    //////////////////////////////
                    // occ_1
                    //////////////////////////////
                    DPWorkInitialise(dpWork);
                    dpWork->leadVerifyLength = 0;
                    dpWork->tailVerifyLength = 0;
                    dpWork->totalSoftClipLength = pedpSetting->SGASoftTotalClipLength * (double)sraQueryInputPos->QueryInfo->ReadLength;
                    dpWork->leadSoftClipLength = pedpSetting->SGASoftHeadClipLength * (double)sraQueryInputPos->QueryInfo->ReadLength;
                    dpWork->tailSoftClipLength = pedpSetting->SGASoftTailClipLength * (double)sraQueryInputPos->QueryInfo->ReadLength;
                    dpWork->leadHardClipLength = 0;
                    dpWork->tailHardClipLength = 0;
                    dpOcc.type = DPOCC_TYPE_AWAIT_TRANSLATE;
                    dpOcc.strand = myNode->pairs[i].occ_1->strand;
                    
                    dpWork->readLength = sraQueryInputPos->QueryInfo->ReadLength;
                    dpWork->regionLength = (double)sraQueryInputPos->QueryInfo->ReadLength * DP_2BWT_REG_EXPAND_FACTOR;
                    dpWork->regionStartIdx = myNode->pairs[i].occ_1->ambPosition - ((double)sraQueryInputPos->QueryInfo->ReadLength * DP_2BWT_REG_SPOS_CORRECTION);
                    
                    //Handle boundary case
                    if ((double)sraQueryInputPos->QueryInfo->ReadLength * DP_2BWT_REG_SPOS_CORRECTION > myNode->pairs[i].occ_1->ambPosition) {
                        dpWork->regionStartIdx = 0;
                    } else if (dpWork->regionStartIdx + dpWork->regionLength > sraIndex->hsp->dnaLength) {
                        dpWork->regionLength = sraIndex->hsp->dnaLength - dpWork->regionStartIdx;
                    }

                    if (myNode->pairs[i].occ_1->strand == QUERY_POS_STRAND) {
                        dpWork->readCode = sraQueryInputPos->QueryInfo->ReadCode;
                    } else if (myNode->pairs[i].occ_1->strand == QUERY_NEG_STRAND) {
                        dpWork->readCode = sraQueryInputNeg->QueryInfo->ReadCode;
                    }
                    
                    #ifdef DP_DEBUG_PRINT_SEED_MATCHING
                        printf("\n[SGASeedEnhancement] DP-verifying occ_1 by %llu (region from s%u-l%u)\n",myNode->pairs[i].occ_1->ambPosition,dpWork->regionStartIdx,dpWork->regionLength);
                    #endif
                    
                    DPMatrixFetching(dpWork,sraIndex->hsp);
                    
                    #ifdef DP_DEBUG_PRINT_SEED_MATCHING
                        DEBUGDPPrintWork(dpWork);
                    #endif
                    
                    if (!DPMatrixFill(dpArguments)) {
                        printf("[SGASeedEnhancement] Unknown error occurred during SeedEnhancement DP-Verification Stage\n");
                        continue;
                    }
                    
                    // DPBackTrack only returns success if the there is an occurrence with
                    // sufficiently high score is found.
                    if (!DPBackTrack(dpWork,dpArguments->dpScores,seedScore,&dpOcc)) {
                        #ifdef DP_DEBUG_PRINT_SEED_MATCHING
                            printf("[SGASeedEnhancement] Occurrence ignored due to insufficient score (%d<%d).\n",dpOcc.matchScore,seedScore);
                        #endif
                        continue;
                    } else {
                        // Log this occurrence if this is a optimal
                        if (dpOcc.matchScore > mapqCalcOverlay_read.rankedScore[0]) {
                            mapqCalcOverlay_read.rankedScore[1] = mapqCalcOverlay_read.rankedScore[0];
                            mapqCalcOverlay_read.rankedScore[0] = dpOcc.matchScore;
                            mapqCalcOverlay_read.rankedCount[0] = 1;
                        } else if (dpOcc.matchScore == mapqCalcOverlay_read.rankedScore[0]) {
                            mapqCalcOverlay_read.rankedCount[0]++;
                        } else if (dpOcc.matchScore > mapqCalcOverlay_read.rankedScore[1]) {
                            mapqCalcOverlay_read.rankedScore[1] = dpOcc.matchScore;
                        }
                    }
                    
                    
                    //////////////////////////////
                    // occ_2
                    //////////////////////////////
                    DPWorkInitialise(dpWork);
                    dpWork->leadVerifyLength = 0;
                    dpWork->tailVerifyLength = 0;
                    dpWork->totalSoftClipLength = pedpSetting->SGASoftTotalClipLength * (double)sraQueryInputPos_mate->QueryInfo->ReadLength;
                    dpWork->leadSoftClipLength = pedpSetting->SGASoftHeadClipLength * (double)sraQueryInputPos_mate->QueryInfo->ReadLength;
                    dpWork->tailSoftClipLength = pedpSetting->SGASoftTailClipLength * (double)sraQueryInputPos_mate->QueryInfo->ReadLength;
                    dpWork->leadHardClipLength = 0;
                    dpWork->tailHardClipLength = 0;
                    dpOcc_mate.type = DPOCC_TYPE_AWAIT_TRANSLATE;
                    dpOcc_mate.strand = myNode->pairs[i].occ_2->strand;
                    
                    dpWork->readLength = sraQueryInputPos_mate->QueryInfo->ReadLength;
                    dpWork->regionLength = (double)sraQueryInputPos_mate->QueryInfo->ReadLength * DP_2BWT_REG_EXPAND_FACTOR;
                    dpWork->regionStartIdx = myNode->pairs[i].occ_2->ambPosition - ((double)sraQueryInputPos_mate->QueryInfo->ReadLength * DP_2BWT_REG_SPOS_CORRECTION);
                    
                    // Hard Clipping is only performed on occ_2 because when occ_1 is being found there is no hard bound..
                    if ( myNode->pairs[i].occ_1->ambPosition < myNode->pairs[i].occ_2->ambPosition ) {
                        // occ_1 being left leg, occ_2 being right leg
                        // right leg cannot extend to the left of the left leg
                        if ( dpWork->regionStartIdx < dpOcc.ambPosition ) {
                            dpWork->leadHardClipLength = dpOcc.ambPosition - dpWork->regionStartIdx;
                        }
                    } else {
                        // occ_1 being right leg, occ_2 being left leg
                        // left leg cannot extend to the right of the right leg
                        if ( dpWork->regionStartIdx + dpWork->regionLength - 1 > dpOcc.ambPosition + dpOcc.matchLen - 1 ) {
                            dpWork->tailHardClipLength = dpWork->regionStartIdx + dpWork->regionLength - dpOcc.ambPosition - dpOcc.matchLen;
                        }
                    }

                    //Handle boundary case
                    if ((double)sraQueryInputPos_mate->QueryInfo->ReadLength * DP_2BWT_REG_SPOS_CORRECTION > myNode->pairs[i].occ_2->ambPosition) {
                        dpWork->regionStartIdx = 0;
                    } else if (dpWork->regionStartIdx + dpWork->regionLength > sraIndex->hsp->dnaLength) {
                        dpWork->regionLength = sraIndex->hsp->dnaLength - dpWork->regionStartIdx;
                    }
                    
                    if (myNode->pairs[i].occ_2->strand == QUERY_POS_STRAND) {
                        dpWork->readCode = sraQueryInputPos_mate->QueryInfo->ReadCode;
                    } else if (myNode->pairs[i].occ_2->strand == QUERY_NEG_STRAND) {
                        dpWork->readCode = sraQueryInputNeg_mate->QueryInfo->ReadCode;
                    }
                    
                    #ifdef DP_DEBUG_PRINT_SEED_MATCHING
                        printf("\n[SGASeedEnhancement] DP-verifying occ_2 by %llu (region from s%u-l%u)\n",myNode->pairs[i].occ_2->ambPosition,dpWork->regionStartIdx,dpWork->regionLength);
                    #endif
                    
                    DPMatrixFetching(dpWork,sraIndex->hsp);
                    
                    #ifdef DP_DEBUG_PRINT_SEED_MATCHING
                        DEBUGDPPrintWork(dpWork);
                    #endif
                    
                    if (!DPMatrixFill(dpArguments)) {
                        printf("[SGASeedEnhancement_mate] Unknown error occurred during SeedEnhancement DP-Verification Stage\n");
                        continue;
                    }
                    
                    // DPBackTrack only returns success if the there is an occurrence with
                    // sufficiently high score is found.
                    if (!DPBackTrack(dpWork,dpArguments->dpScores,seedScore_mate,&dpOcc_mate)) {
                        #ifdef DP_DEBUG_PRINT_SEED_MATCHING
                            printf("[SGASeedEnhancement] Occurrence ignored due to insufficient score.\n");
                        #endif
                        continue;
                    } else {
                        // Log this occurrence if this is a optimal
                        if (dpOcc_mate.matchScore > mapqCalcOverlay_mate.rankedScore[0]) {
                            mapqCalcOverlay_mate.rankedScore[1] = mapqCalcOverlay_mate.rankedScore[0];
                            mapqCalcOverlay_mate.rankedScore[0] = dpOcc_mate.matchScore;
                            mapqCalcOverlay_mate.rankedCount[0] = 1;
                        } else if (dpOcc_mate.matchScore == mapqCalcOverlay_mate.rankedScore[0]) {
                            mapqCalcOverlay_mate.rankedCount[0]++;
                        } else if (dpOcc_mate.matchScore > mapqCalcOverlay_mate.rankedScore[1]) {
                            mapqCalcOverlay_mate.rankedScore[1] = dpOcc_mate.matchScore;
                        }
                    }
                    
                    //////////////////////////////
                    // Reporting Pairs
                    //////////////////////////////
                    
                    #ifdef PE_DP_FORCE_INSERTION_SIZE
                        unsigned int insertion;
                        if (PEDPIsPairEndMatch(peInput,&dpOcc,&dpOcc_mate,&insertion) || 
                            PEDPIsPairEndMatch(peInput,&dpOcc_mate,&dpOcc,&insertion)) {
                            continue;
                        }
                    #endif
                
                    #ifdef DP_DEBUG_PRINT_SEED_MATCHING
                        printf("[SGASeedEnhancement] ReadID#%llu:%llu Found DP Occ %llu and %llu\n",
                            sraQueryInputPos->QueryInfo->ReadId, sraQueryInputPos_mate->QueryInfo->ReadId, 
                            dpOcc.ambPosition,dpOcc_mate.ambPosition);
                    #endif
                            
                    PEDPOCCDumpOneSeedAlignment(peArguments,&dpOcc,&dpOcc_mate);
                    peCount++;
                    
                    // Early terminate if maxResult is reached
                    if (peInput->maxResult != -1 && 
                        peCount >= peInput->maxResult) {
                        dpAlgnmtOK = FALSE;
                        break;
                    }
                    
                    // Early terminate if one alignment is already found for random best alignment
                    if (peInput->OutputType == PE_REPORT_RANDOM_BEST) {
                        dpAlgnmtOK = FALSE;
                        break;
                    }
                }
                
                myNode = myNode->next;
            }
            
            if (sraQueryInputPos->MapqCalc->status == MAPQ_CALC_FIGURE_INITIALISED) {
                MAPQCalculatorOverlay(&mapqCalcOverlay_read, sraQueryInputPos->MapqCalc);
                sraQueryInputPos->MapqCalc->status == MAPQ_CALC_FIGURE_READY;
            }
            if (sraQueryInputPos_mate->MapqCalc->status == MAPQ_CALC_FIGURE_INITIALISED) {
                MAPQCalculatorOverlay(&mapqCalcOverlay_mate, sraQueryInputPos_mate->MapqCalc);
                sraQueryInputPos->MapqCalc->status == MAPQ_CALC_FIGURE_READY;
            }
        }

        free(seedOccLists_1);
        free(seedOccLists_2);
    }
        
    /////////////////////////////////////////////////////////////////////
    // Step 8 Memory Clean Up
    /////////////////////////////////////////////////////////////////////
    SRAAlgnmtOutputFreeNone(seedSraArgPos->AlgnmtOutput);
    SRAAlgnmtOutputFreeNone(seedSraArgPos_mate->AlgnmtOutput);
    SRAARGSlaveFree(seedSraArgNeg_mate);
    SRAARGSlaveFree(seedSraArgNeg);
    SRAARGMateFree(seedSraArgPos_mate);
    SRAARGMateFree(seedSraArgPos);
    
    return peCount;
//}
}

unsigned int DPPEOrphanMappingOccurrences(PEArguments * peArguments, SRAOccurrence * occList_base, unsigned long long occCount_base, unsigned int readLength_base,
                            unsigned char * readCode, unsigned char * readCode_neg, unsigned int readLength, MAPQCalculator * mapqCalc,
                            int dpBase) {

    PEInput * peInput = peArguments->PEAlgnmtInput;
        
    int peCount = 0;
    SRAOutput * peOutput_OccHandler = peArguments->PEAlgnmtOutput_OccHandle;
    DPArguments * dpArguments = peArguments->dpArguments;
    DPWork * dpWork = dpArguments->dpWork;
    SRAIndex * sraIndex = peArguments->sraArgsPos->AlgnmtIndex;
    int strandLeftLeg = peInput->strandLeftLeg;
    int strandRightLeg = peInput->strandRightLeg;
    
    PEDPSetting * pedpSetting = peArguments->pedpSetting;
    
    unsigned int i;
    unsigned int expectedPos = 0;
    unsigned int expectedLength = peInput->insertUbound - peInput->insertLbound + readLength - 1;
    
    int peAlgnmtOk = TRUE;
    
    unsigned int reportedPeCount = PEOCCCountAllPEAlignment(peArguments);
    if ( peInput->maxResult != -1 ) {
        if ( reportedPeCount >= peInput->maxResult ) {
            peAlgnmtOk = FALSE;
        }
    }

    DPOccurrence dpOcc_found;
    
    MAPQCalculator mapqCalcOverlay;
    MAPQCalculatorInitialise(&mapqCalcOverlay);

    dpWork->regionLength = expectedLength;
    dpWork->readLength = readLength;
    
    int InputSGAOrphanRecoverScore = pedpSetting->SGAScoreTF * readLength;
    
    // -----------------------------------------------------------------------------
    //    =====================                      [                    ]
    //        Aligned (SRA)                              Un-aligned (DP)
    //
    if (peAlgnmtOk) {
        dpOcc_found.type = DPOCC_TYPE_AWAIT_TRANSLATE;
        dpOcc_found.strand = strandRightLeg;
        if (strandRightLeg == QUERY_POS_STRAND) {
            dpWork->readCode = readCode;
        } else if (strandRightLeg == QUERY_NEG_STRAND) {
            dpWork->readCode = readCode_neg;
        }

        dpWork->leadVerifyLength = 0;
        dpWork->tailVerifyLength = 0;

        dpWork->totalSoftClipLength = pedpSetting->SGASoftTotalClipLength * (double)readLength;
        dpWork->leadSoftClipLength = pedpSetting->SGASoftHeadClipLength * (double)readLength;
        dpWork->tailSoftClipLength = pedpSetting->SGASoftTailClipLength * (double)readLength;

        dpWork->leadHardClipLength = 0;
        dpWork->tailHardClipLength = 0;
        
        for (i=0;i<occCount_base;i++) {
            if (occList_base[i].strand == strandLeftLeg) {
                if (occList_base[i].ambPosition + peInput->insertLbound + expectedLength < sraIndex->hsp->dnaLength) {
                    DPWorkInitialise(dpWork);
                    
                    expectedPos = occList_base[i].ambPosition + peInput->insertLbound - readLength;
                    dpWork->regionStartIdx = expectedPos;
                    DPMatrixFetching(dpWork,sraIndex->hsp);
                    
                    // right leg cannot extend to the left of the left leg
                    if ( dpWork->regionStartIdx < occList_base[i].ambPosition ) {
                        dpWork->leadHardClipLength = occList_base[i].ambPosition - dpWork->regionStartIdx;
                    }
                    
                    if (!DPMatrixFill(dpArguments)) {
                        printf("[DPPEOrphanMappingOccurrences] Unknown error occurred during OrphanRecover\n");
                        continue;
                    }
                    
                    // DPBackTrack only returns success if the there is an occurrence with
                    // sufficiently high score is found.
                    if (DPBackTrack(dpWork,dpArguments->dpScores,InputSGAOrphanRecoverScore,&dpOcc_found)) {
                    
                        #ifdef PE_DP_FORCE_INSERTION_SIZE
                            unsigned int insertion;
                            if (!PEIsPairEndMatch(peInput,
                                            occList_base[i].ambPosition,occList_base[i].strand,occList_base[i].matchLen,
                                            dpOcc_found.ambPosition,dpOcc_found.strand,dpOcc_found.matchLen,
                                            &insertion)) {
                                continue;
                            }
                        #endif
                        
                        PEDPOCCDumpOneOrphanAlignment(peArguments,&(occList_base[i]),&dpOcc_found,dpBase);
                        reportedPeCount++;
                        peCount++;
                        
                        // Log this occurrence if this is a optimal
                        if (dpOcc_found.matchScore > mapqCalcOverlay.rankedScore[0]) {
                            mapqCalcOverlay.rankedScore[1] = mapqCalcOverlay.rankedScore[0];
                            mapqCalcOverlay.rankedScore[0] = dpOcc_found.matchScore;
                            mapqCalcOverlay.rankedCount[0] = 1;
                        } else if (dpOcc_found.matchScore == mapqCalcOverlay.rankedScore[0]) {
                            mapqCalcOverlay.rankedCount[0]++;
                        } else if (dpOcc_found.matchScore > mapqCalcOverlay.rankedScore[1]) {
                            mapqCalcOverlay.rankedScore[1] = dpOcc_found.matchScore;
                        }
                        
                        // Early terminate if maxResult is reached
                        if (peInput->maxResult != -1 && 
                            reportedPeCount >= peInput->maxResult) {
                            peAlgnmtOk = FALSE;
                            break;
                        }
                        
                        // Early terminate if one alignment is already found for random best alignment
                        if (peInput->OutputType == PE_REPORT_RANDOM_BEST) {
                            peAlgnmtOk = FALSE;
                            break;
                        }
                    }
                }
            }
        }
    }
    
    // -----------------------------------------------------------------------------
    //    [                    ]                       ======================
    //        Un-aligned (DP)                               Aligned (SRA)
    //
    if (peAlgnmtOk) {
        dpOcc_found.type = DPOCC_TYPE_AWAIT_TRANSLATE;
        dpOcc_found.strand = strandLeftLeg;
        if (strandLeftLeg == QUERY_POS_STRAND) {
            dpWork->readCode = readCode;
        } else if (strandLeftLeg == QUERY_NEG_STRAND) {
            dpWork->readCode = readCode_neg;
        }
        
        dpWork->leadVerifyLength = 0;
        dpWork->tailVerifyLength = 0;
        
        dpWork->totalSoftClipLength = pedpSetting->SGASoftTotalClipLength * (double)readLength;
        dpWork->leadSoftClipLength = pedpSetting->SGASoftHeadClipLength * (double)readLength;
        dpWork->tailSoftClipLength = pedpSetting->SGASoftTailClipLength * (double)readLength;
        
        dpWork->leadHardClipLength = 0;
        dpWork->tailHardClipLength = 0;
        
        for (i=0;i<occCount_base;i++) {
            if (occList_base[i].strand == strandRightLeg) {
                if (occList_base[i].ambPosition >= peInput->insertUbound + readLength) {
                
                    DPWorkInitialise(dpWork);
                    
                    expectedPos = occList_base[i].ambPosition - peInput->insertUbound + readLength;
                    dpWork->regionStartIdx = expectedPos;
                    DPMatrixFetching(dpWork,sraIndex->hsp);
                    
                    // left leg cannot extend to the right of the right leg
                    if ( dpWork->regionStartIdx + dpWork->regionLength - 1 > occList_base[i].ambPosition + readLength - 1 ) {
                        dpWork->tailHardClipLength = dpWork->regionStartIdx + dpWork->regionLength - occList_base[i].ambPosition - readLength;
                    }
                    
                    if (!DPMatrixFill(dpArguments)) {
                        printf("[DPPEOrphanMappingOccurrences] Unknown error occurred during OrphanRecover\n");
                        continue;
                    }
                    
                    // DPBackTrack only returns success if the there is an occurrence with
                    // sufficiently high score is found.
                    if (DPBackTrack(dpWork,dpArguments->dpScores,InputSGAOrphanRecoverScore,&dpOcc_found)) {
                    
                        #ifdef PE_DP_FORCE_INSERTION_SIZE
                            unsigned int insertion;
                            if (!PEIsPairEndMatch(peInput,
                                            dpOcc_found.ambPosition,dpOcc_found.strand,dpOcc_found.matchLen,
                                            occList_base[i].ambPosition,occList_base[i].strand,occList_base[i].matchLen,
                                            &insertion)) {
                                continue;
                            }
                        #endif
                        
                        PEDPOCCDumpOneOrphanAlignment(peArguments,&(occList_base[i]),&dpOcc_found,dpBase);
                        reportedPeCount++;
                        peCount++;
                        
                        // Log this occurrence if this is a optimal
                        if (dpOcc_found.matchScore > mapqCalcOverlay.rankedScore[0]) {
                            mapqCalcOverlay.rankedScore[1] = mapqCalcOverlay.rankedScore[0];
                            mapqCalcOverlay.rankedScore[0] = dpOcc_found.matchScore;
                            mapqCalcOverlay.rankedCount[0] = 1;
                        } else if (dpOcc_found.matchScore == mapqCalcOverlay.rankedScore[0]) {
                            mapqCalcOverlay.rankedCount[0]++;
                        } else if (dpOcc_found.matchScore > mapqCalcOverlay.rankedScore[1]) {
                            mapqCalcOverlay.rankedScore[1] = dpOcc_found.matchScore;
                        }
                    
                        // Early terminate if maxResult is reached
                        if (peInput->maxResult != -1 && 
                            reportedPeCount >= peInput->maxResult) {
                            peAlgnmtOk = FALSE;
                            break;
                        }
                        
                        // Early terminate if one alignment is already found for random best alignment
                        if (peInput->OutputType == PE_REPORT_RANDOM_BEST) {
                            peAlgnmtOk = FALSE;
                            break;
                        }
                    }
                }
            }
        }
    }
    
    // Wrapping up alignments
    if (mapqCalc->status == MAPQ_CALC_FIGURE_INITIALISED) {
        MAPQCalculatorOverlay(&mapqCalcOverlay, mapqCalc);
        mapqCalc->status == MAPQ_CALC_FIGURE_READY;
    }
    
    return peCount;
}


void _PEPopulateMAPQCalculator(PEArguments * peArguments) {

    MAPQCalculator * mapqCalc = peArguments->MapqCalc;
    
    int readLength = peArguments->sraArgsPos->QueryInfo->ReadLength;
    int mateLength = peArguments->sraArgsPos_mate->QueryInfo->ReadLength;
    DPScores * dpScores = peArguments->dpArguments->dpScores;
    
    // Setting up the PE SRA Occ Collector
    SRAOutput * peOutput_OccHandler = peArguments->PEAlgnmtOutput_OccHandle;
    SRAOCCCollector * occCollector = peOutput_OccHandler->occCollector;
    SRAOCCCollectorPointer * sraOccIterator = SRAOCCPTCreate(occCollector);
    SRAOccurrence * peOcc_1 = SRAOCCPTRead(sraOccIterator);
    SRAOccurrence * peOcc_2;
    
    // Setting up the PE DP Occ Collector
    DPOCCCollector * dpOccCollector = peArguments->dpArguments->dpOccCollector;
    DPOCCCollectorPointer * dpOccIterator = DPOCCPTCreate(dpOccCollector);
    DPOccurrence * dpOcc_1 = DPOCCPTRead(dpOccIterator);
    DPOccurrence * dpOcc_2;
    
    int16_t optimalScore = 0;
    int16_t suboptimalScore = 0;
    unsigned int optimalCount = 0;
    int16_t peOccScore = 0;
    int16_t peTtlOccCount = 0;
    
    while (peOcc_1!=NULL) {
    
        if (peOcc_1->type == SRAOCC_TYPE_DELIMITOR_READ) {
        } else if (peOcc_1->type == SRAOCC_TYPE_PE_PAIR) {
            // Proper SRA-SRA PE Alignment Results
            SRAOCCPTNext(sraOccIterator);
            peOcc_2 = SRAOCCPTRead(sraOccIterator);
            
            peOccScore = MAPQNormaliseSRAScore(readLength,peOcc_1->mismatchCount,dpScores) + 
                            MAPQNormaliseSRAScore(readLength,peOcc_2->mismatchCount,dpScores);
                            
            peTtlOccCount++;
            
        } else if (peOcc_1->type==SRAOCC_TYPE_PE_DP_BASE_READ ||
                    peOcc_1->type==SRAOCC_TYPE_PE_DP_BASE_MATE) {
                    
            // Orphan-DP PE Alignment Results
            DPOCCPTNext(dpOccIterator);
            dpOcc_2 = DPOCCPTRead(dpOccIterator);
            
            peOccScore = MAPQNormaliseSRAScore(readLength,peOcc_1->mismatchCount,dpScores) + 
                            dpOcc_1->matchScore;
                            
            peTtlOccCount++;
            
        } else if (peOcc_1->type==SRAOCC_TYPE_PE_DP_SEED_OCC) {
            // Seed-DP PE Alignment Results
            DPOCCPTNext(dpOccIterator);
            dpOcc_2 = DPOCCPTRead(dpOccIterator);
            
            peOccScore = dpOcc_1->matchScore + 
                            dpOcc_2->matchScore;
            
            peTtlOccCount++;
            
            DPOCCPTNext(dpOccIterator);
            dpOcc_1 = DPOCCPTRead(dpOccIterator);
        } else {
            // ATTENTION
            printf("[_PEPopulateMAPQCalculator] Unknown Pair-End Alignment Type %d\n",peOcc_1->type);
            peOccScore = 0;
        }
        
        // Loop-through all PE alignment results to find
        // the optimal and sub-optimal alignments
        if (peOccScore > optimalScore) {
            suboptimalScore = optimalScore;
            optimalScore = peOccScore;
            optimalCount = 1;
        } else if (peOccScore == optimalScore) {
            optimalCount++;
        } else if (peOccScore > suboptimalScore) {
            suboptimalScore = peOccScore;
        }
        
        SRAOCCPTNext(sraOccIterator);
        peOcc_1 = SRAOCCPTRead(sraOccIterator);
    }
    
    
    // Get SRA Optimal and Sub-Optimal
    int resultRank = 0;
    mapqCalc->rankedMismatch[resultRank] = 0;
    mapqCalc->rankedCount[resultRank] = optimalCount;
    mapqCalc->rankedScore[resultRank] = optimalScore;
    
    resultRank++;
    mapqCalc->rankedMismatch[resultRank] = 0;
    mapqCalc->rankedCount[resultRank] = peTtlOccCount - optimalCount;
    mapqCalc->rankedScore[resultRank] = suboptimalScore;
    
    #ifdef PE_DEBUG_PRINT_MAPQ_CALCULATION
        printf("[_PEPopulateMAPQCalculator] SRAOcc = %u / DPOcc = %u\n",SRAOCCCountOccurrences(occCollector),DPOCCCountOccurrences(dpOccCollector));
        printf("[_PEPopulateMAPQCalculator] detects %u alignment results\n",peTtlOccCount);
        MAPQPrintParameters(mapqCalc);
    #endif
    
    mapqCalc->status = MAPQ_CALC_FIGURE_READY;
}
