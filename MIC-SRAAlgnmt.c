//
//    MIC-SRAAlgnmt.c
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

#include "MIC-SRAAlgnmt.h"

// --------------------------------------------------
// DP_DEBUG #1
// Uncomment the below to turn on logging message for
// PE matching for seed enhancement.
//#define DP_DEBUG_PRINT_SEED_MATCHING
// --------------------------------------------------

__attribute__((target(mic)))
void _MICDebugPrintArgument (MICSRAArguments * micArgs);
__attribute__((target(mic)))
void _MICDebugPrintOccurrences (MICSRAArguments * micArgs);

__attribute__((target(mic)))
unsigned long long MICProcessStage(MICSRAArguments * micArgs, CPTSRAModel * cpModel, int * caseId) {
    SRAWorkingMemory * workMem = &(micArgs->AlgnmtMemory);
    BWT * bwt = micArgs->bwt;
    int OutputType = micArgs->outputType;
    
    unsigned long long saRanges[4];
    unsigned long long saCount = 0;
    CPTSRACase * thisCase = &(cpModel->cases[(*caseId)]);
    
    while (thisCase->type != SRA_CASE_TYPE_NOT_INITALISED && 
           thisCase->type != SRA_CASE_TYPE_NEXT_STAGE) {
           
        saRanges[0] = 0;
        saRanges[1] = bwt->textLength;
        saRanges[2] = 0;
        saRanges[3] = bwt->textLength;

        if (thisCase->steps[0].type == SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP) {
            saCount += MICBWTExactModelBackward_Lookup(micArgs,
                                                thisCase,0,
                                                saRanges);
        } else if (thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP) {
            saCount += MICBWTExactModelForward_Lookup(micArgs,
                                                thisCase,0,
                                                saRanges);
        } else if (thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_BACKWARD_LOOKUP) {
            saCount += MICBWTExactModelBackwardAnyDirection_Lookup(micArgs,
                                                thisCase,0,
                                                saRanges);
        } else if (thisCase->steps[0].type == SRA_STEP_TYPE_BACKWARD_ONLY_BWT) {
            saCount += MICBWTModelSwitchBackward(micArgs,0,0,
                                                thisCase,0,
                                                saRanges, 0, 0, 0);
        } else if (thisCase->steps[0].type == SRA_STEP_TYPE_BI_DIRECTIONAL_BWT) {
            saCount += MICBWTModelSwitchAnyDirection(micArgs,0,0,
                                                thisCase,0,
                                                saRanges, 0, 0, 0);
        }
        
        // Comment out random-best alignment interim.
        //if (workMem->IsClosed==1) {
        //    if (OutputType==SRA_REPORT_RANDOM_BEST) {
        //        OCCFlushBestIntoCache(micArgs);
        //    }
        //    return saCount;
        //}
        (*caseId)++;
        thisCase = &(cpModel->cases[(*caseId)]);
    }
    
    return saCount;
}


// Function supports all-valid alignment and all-best alignment
// It will terminate when the all work required to obtain all-best/all-valid alignment is done.
// It will also early terminate if the remaining of the alignment reported will be discarded anyway. (all-best scenario)
__attribute__((target(mic)))
unsigned long long MICProcessStageDoubleStrand(MICSRAArguments * micArgs, CPTSRAModel * cpPModels, CPTSRAModel * cpNModels, int * caseId) {

    unsigned long long saCount = 0;
    int k;
    unsigned char * readCode = micArgs->readCode;
    uint16_t seedOffset = micArgs->seedOffset;
    
    k = (*caseId);
    micArgs->readStrand = QUERY_POS_STRAND;
    micArgs->readCode = readCode;
    micArgs->seedOffset = seedOffset;
    saCount += MICProcessStage(micArgs,cpPModels,&k);

    k = (*caseId);
    micArgs->readStrand = QUERY_NEG_STRAND;
    micArgs->readCode = micArgs->readCode_Complt;
    micArgs->seedOffset = micArgs->seedOffset_Complt;
    saCount += MICProcessStage(micArgs,cpNModels,&k);
    
    (*caseId) = k;

    // Restoring the readStrand and readCode for any downstream process
    micArgs->readStrand = QUERY_POS_STRAND;
    micArgs->readCode = readCode;
    micArgs->seedOffset = seedOffset;
    
    return saCount;
}


// Function supports all-valid alignment and all-best alignment
// It will terminate when the all work required to obtain all-best/all-valid alignment is done.
// It will also early terminate if the remaining of the alignment reported will be discarded anyway. (all-best scenario)
__attribute__((target(mic)))
unsigned long long MICProcessReadDoubleStrand(MICSRAArguments * micArgs, CPTSRAModel * cpPModels, CPTSRAModel * cpNModels) {

    unsigned long long saCount = 0;
    SRAWorkingMemory * workMem = &(micArgs->AlgnmtMemory);
    int OutputType = micArgs->outputType;
    int j,k;
    
    SRAWorkingMemoryInitialise(workMem);
    
    unsigned char * readCode = micArgs->readCode;
    uint16_t seedOffset = micArgs->seedOffset;
    
    j=0;
    while (1) {

        //Report unique best alignment
        //Report all best alignment
        k=j;
        micArgs->readStrand = QUERY_POS_STRAND;
        micArgs->readCode = readCode;
        micArgs->seedOffset = seedOffset;
        saCount += MICProcessStage(micArgs,cpPModels,&k);
        
        // Comment out random-best alignment interim.
        //if (OutputType==SRA_REPORT_RANDOM_BEST &&
        //    workMem->IsClosed==1) {
        //    
        //    return saCount;
        //}

        k=j;
        micArgs->readStrand = QUERY_NEG_STRAND;
        micArgs->readCode = micArgs->readCode_Complt;
        micArgs->seedOffset = micArgs->seedOffset_Complt;
        saCount += MICProcessStage(micArgs,cpNModels,&k);
        
        // Comment out random-best alignment interim.
        //if (OutputType==SRA_REPORT_RANDOM_BEST &&
        //    workMem->IsClosed==1) {
        //    
        //    return saCount;
        //}
 
        if (OutputType==SRA_REPORT_ALL_BEST && saCount>0) {
            break;
        }
        
        if (cpPModels->cases[k].type==SRA_CASE_TYPE_NOT_INITALISED ||
            cpNModels->cases[k].type==SRA_CASE_TYPE_NOT_INITALISED) {
            if (cpPModels->cases[k].type!=SRA_CASE_TYPE_NOT_INITALISED ||
                cpNModels->cases[k].type!=SRA_CASE_TYPE_NOT_INITALISED) {
                printf("[ERROR] Imbalance SRA Alignment model!\n");
                exit(1);
            } else {
                break;
            }
        }

        j=k+1;
    }

    // Comment out best-quality alignment interim.
    //if (OutputType==SRA_REPORT_BEST_QUALITY && workMem->IsOpened==1) {
    //    OCCFlushBestIntoCache(micArgs);
    //    return 1;
    //}

    // Restoring the readStrand and readCode for any downstream process
    micArgs->readStrand = QUERY_POS_STRAND;
    micArgs->readCode = readCode;
    micArgs->seedOffset = seedOffset;
    
    return saCount;
}

__attribute__((target(mic)))
void _MICDebugPrintArgument (MICSRAArguments * micArgs) {
    int i;
    
    printf("Align + ");
    for (i=0;i<micArgs->seedLength;i++) {
        if (micArgs->readCode[i + micArgs->seedOffset]==0) printf("A");
        if (micArgs->readCode[i + micArgs->seedOffset]==1) printf("C");
        if (micArgs->readCode[i + micArgs->seedOffset]==2) printf("G");
        if (micArgs->readCode[i + micArgs->seedOffset]==3) printf("T");
    }
    printf(" (%d)\n",micArgs->seedLength);
    
    printf("Align - ");
    for (i=0;i<micArgs->seedLength;i++) {
        if (micArgs->readCode_Complt[i + micArgs->seedOffset_Complt]==0) printf("A");
        if (micArgs->readCode_Complt[i + micArgs->seedOffset_Complt]==1) printf("C");
        if (micArgs->readCode_Complt[i + micArgs->seedOffset_Complt]==2) printf("G");
        if (micArgs->readCode_Complt[i + micArgs->seedOffset_Complt]==3) printf("T");
    }
    printf(" (%d)\n",micArgs->seedLength);
    printf("maxNBMismatch = %u\n",micArgs->maxNBMismatch);
    printf("seedOffset = %u\n",micArgs->seedOffset);
    printf("seedOffset_Complt = %u\n",micArgs->seedOffset_Complt);
}

__attribute__((target(mic)))
void _MICDebugPrintOccurrences (MICSRAArguments * micArgs) {

    int i;
    uint8_t * outputStatus = micArgs->outputStatus;
    unsigned int * outputBlock = micArgs->outputBlock;
    MICSRAOccMetadata * metaBlock = micArgs->metaBlock;
    
    for (i=0;i<(*micArgs->occCount);i++) {
        printf("Occurrences : pos=%15u strand=%d mis=%2d\n",outputBlock[i],metaBlock[i].strand,metaBlock[i].numOfErr);
    }
}

__attribute__((target(mic)))
unsigned long long MICSeedAlgnmtDoubleStrand(MICSRAArguments * micArgs, CPTSRAModel * cpSeedModel, CPTSRAModel * cpSeedModel_neg,
                                                    uint16_t seedLength, uint16_t seedUniqeLength, int occLimit) {
                                             
    unsigned long long saCount = 0;
    
    MICSRAArguments micSeedArgs;
    MICSRAARGCopy(&micSeedArgs,micArgs);
    
    // Debug printing of the SRA model and the SRA argument
    //_MICDebugPrintArgument(&micSeedArgs);
    //MICSRADebugPrintModel(cpSeedModel,stdout);
    
    uint16_t seedOffset = 0;
    uint16_t readLength = micArgs->seedLength;
    
    // Position strand starts from the left-most position
    // Setup the starting position of the negative strand as it starts from the right-most position minus the seedLength
    micSeedArgs.seedLength      = seedLength;
    micSeedArgs.readCode_Complt = micArgs->readCode_Complt;
    		//+ micArgs->readLength - seedLength;
    
    while (seedOffset + seedLength <= readLength) {

        ///////////////////////////////////////////////
        // Prepare Seed Alignment Correction
        ///////////////////////////////////////////////
        micSeedArgs.seedOffset         = seedOffset;
        micSeedArgs.seedOffset_Complt  = readLength - seedOffset - seedLength;
        
        saCount += MICProcessReadDoubleStrand(&micSeedArgs, cpSeedModel, cpSeedModel_neg);
    
        if (occLimit != -1 && saCount>occLimit) {
            // We don't need to proceed further as the number of seed alignment has
            // exceeded the limitation defined in configuration.
            return saCount;
        }
    
//        micSeedArgs.readCode        += seedUniqeLength;
//        micSeedArgs.readCode_Complt -= seedUniqeLength;
        seedOffset += seedUniqeLength;
    }
    
    return saCount;
}
