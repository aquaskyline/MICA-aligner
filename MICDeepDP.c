#include "MICDeepDP.h" 

    ///////////////////////////////////////
    // Deep-DP Flow (sh: obsolete)
    ///////////////////////////////////////
    //1. For the read, perform SRA on the seeds. Store the result.
    //   tweak the occurrences to store the anticipated occurrence position.
    //2. For the mate, perform SRA on the seeds. Store the result.
    //   tweak the occurrences to store the anticipated occurrence position.
    //3. Check the number of occ reported by either read/mate should not exceed SGASeedOccLimitation.
    //4. Do a PE on the occurrences to find the candidate. The result will be stored in peArguments->PEAlgnmtOutput_OccHandle.
    //4a.We merge the results that are too close to each other. Too close is defined as distance within 25% of the read length.
    //5. For each pair of result it indicates that it could be a possible occurrence. We then verify the candidate by dynamic programming.
    //6. Post-verification validation to make sure the occurrence falls between the allowed insertion size.
    //7. Report the occurrences
    
    //1. Performed in MICController by MIC-SRAAlgnmt
    
    //2. Performed in MICController by MIC-SRAAlgnmt
    
    //3. Performed in MICController
    
    //4. Performed in MICController
    
    //4a.Not yet implemented ATTENTION
    
    //5.
    
__attribute__((target(mic)))
static inline int MICDPDeepAlignment(MICDPOccurrence * dpOutputPtr, DPWorkMIC * dpWork, const MICPEArguments * peArgs, const HSP * hsp) {

    DPArgumentsMIC dpArgs;
    int dpOccCount = 0;
    int peOccIdx = 0;
    
    MICSRAArguments * readArgs = peArgs->readArgs;
    MICSRAArguments * mateArgs = peArgs->mateArgs;
    
    unsigned int readLength = readArgs->seedLength;
    unsigned int mateLength = mateArgs->seedLength;
    unsigned int regionStartPos;
    
    MICSRAOccMetadata * peMeta     = peArgs->outputMeta;
    unsigned int *      peOutput   = peArgs->output;
    int j;

    for (peOccIdx=0; peOccIdx<(*peArgs->occCount); peOccIdx+=2) {

        ///////////////////////////////////////
        // Setting up parameters re. READ
        ///////////////////////////////////////
        dpArgs.readLength = readLength;
        if (peMeta[peOccIdx].strand == 0) {
            dpArgs.readCode = readArgs->readCode + readArgs->seedOffset ;
        } else {
            dpArgs.readCode = readArgs->readCode_Complt + readArgs->seedOffset_Complt;
        }
        ///////////////////////////////////////
        // Setting up parameters for the anticipated region on Reference Sequence
        ///////////////////////////////////////
        dpArgs.regLength = (double)readLength * 1.50f;
        regionStartPos = peOutput[peOccIdx] - ((double)readArgs->seedLength * 0.25f);
        
        //Handle boundary case
        if (peOutput[peOccIdx] < ((double)readArgs->seedLength * 0.25f)) regionStartPos = 0;
        if (regionStartPos + dpArgs.regLength > hsp->dnaLength) dpArgs.regLength = hsp->dnaLength - regionStartPos;

        // Extract the reference sequence from HSP with the estimated location from the BASE
        dpArgs.regCode = MICDPExtractCodeFromHSP(
            dpWork->regBuffer,
            regionStartPos,
            dpArgs.regLength,
            hsp
        );

        // Invocation of the MIC-DP Module
        // Report DP occurrence of the score is greater than the threshold
        if (DPMatrixFillMIC(&dpArgs, dpWork)){
            if (dpOccCount >= MIC_DP_OUTPUT_MAX_ALIGNMENT) {
                return MIC_DP_STATUS_BREACH_LIMIT;
            }
            if (!_MICDPPopulateOccurrence(dpOutputPtr+dpOccCount, dpWork, regionStartPos, peMeta[peOccIdx].strand)) {
                return MIC_DP_STATUS_ERROR; 
            }
        }
        else continue;
        #ifdef MIC_DP_DEBUG_PRINT_SEED_MATCHING
        printf("[SeedDp] DP was carried out [%u-%u] by seed position %u\n",
            regionStartPos,regionStartPos+dpArgs.regLength-1,peOutput[peOccIdx+1]);
        printf("DP Alignment found at %u (SCORE %d > BOUND %.2f)\n",
            regionStartPos + dpWork->maxScoreStartPos,dpWork->maxScore,dpWork->scoreThreshold * dpArgs.readLength);
        #endif

        ///////////////////////////////////////
        // Setting up parameters re. MATE
        ///////////////////////////////////////
        dpArgs.readLength = mateLength;
        if (peMeta[peOccIdx+1].strand == 0) {
            dpArgs.readCode = mateArgs->readCode + mateArgs->seedOffset ;
        } else {
            dpArgs.readCode = mateArgs->readCode_Complt + mateArgs->seedOffset_Complt ;
        }
        ///////////////////////////////////////
        // Setting up parameters for the anticipated region on Reference Sequence
        ///////////////////////////////////////
        dpArgs.regLength = (double)mateLength * 1.50f;
        regionStartPos = peOutput[peOccIdx+1] - ((double)mateArgs->seedLength * 0.25f);
        
        //Handle boundary case
        if (peOutput[peOccIdx+1] < ((double)mateArgs->seedLength * 0.25f)) regionStartPos = 0;
        if (regionStartPos + dpArgs.regLength > hsp->dnaLength) dpArgs.regLength = hsp->dnaLength - regionStartPos;

        // Extract the reference sequence from HSP with the
        // estimated location from the BASE
        dpArgs.regCode = MICDPExtractCodeFromHSP(
            dpWork->regBuffer,
            regionStartPos,
            dpArgs.regLength,
            hsp
        );

        // Invocation of the MIC-DP Module
        // Report DP occurrence of the score is greater than the threshold
        if (DPMatrixFillMIC(&dpArgs, dpWork)){
            if (dpOccCount+1 >= MIC_DP_OUTPUT_MAX_ALIGNMENT) {
                return MIC_DP_STATUS_BREACH_LIMIT;
            }
            if (!_MICDPPopulateOccurrence(dpOutputPtr+dpOccCount+1, dpWork, regionStartPos, peMeta[peOccIdx+1].strand)) {
                return MIC_DP_STATUS_ERROR;
            }
            dpOccCount+=2;
        }
        #ifdef MIC_DP_DEBUG_PRINT_SEED_MATCHING
        printf("[SeedDp] DP was carried out [%u-%u] by seed position %u\n",
            regionStartPos,regionStartPos+dpArgs.regLength-1,peOutput[peOccIdx+1]);
        printf("DP Alignment found at %u (SCORE %d > BOUND %.2f)\n",
            regionStartPos + dpWork->maxScoreStartPos,dpWork->maxScore,dpWork->scoreThreshold * dpArgs.readLength);
        #endif
        
      /*  #ifdef MIC_DP_FORCE_INSERTION_SIZE
            if (dpAlgnmtOK) {
                MICDPOccurrence * dpOcc1 = dpOutputPtr+dpOccCount;
                unsigned int dpAmbPos1 = dpOcc1->ambPosition;
                unsigned int dpAmbPos2 = regionStartPos + dpWork->maxScoreStartPos;
                
                unsigned int insertion_1 = dpAmbPos2 + dpWork->maxScoreEndPos - dpWork->maxScoreStartPos + 1 - dpAmbPos1;
                unsigned int insertion_2 = dpAmbPos1 + dpOcc1->matchLen - dpAmbPos2;
                
                // ATTENTION : This checking logic is a short-term solution before
                // the StrandLeftLeg, StrandRightLeg option being implemented.
                if ((insertion_1 <= peArgs->lBound || insertion_1 >= peArgs->uBound) &&
                    (insertion_2 <= peArgs->lBound || insertion_2 >= peArgs->uBound)) {
                    dpAlgnmtOK = FALSE;
                }
            }
       #endif
        */
    }
    return dpOccCount;
}


__attribute__((target(mic)))
int MICDeepDP(
        MICPEArguments * peSeedArgs, 
        unsigned int * peOutputBuffer,
        MICSRAOccMetadata * peMetaBuffer, 
        MICSRAArguments * readArgs,
        MICSRAArguments * mateArgs,
        MICSRAArguments * readSeedSRAArgs,
        MICSRAArguments * mateSeedSRAArgs,
        CPTSRAModel * cpPModels,
        CPTSRAModel * cpNModels,
        HSP * hsp,
        DPWorkMIC * dpWork,
        double seedLength,
        double seedOverlap,
        double orphanTrigger,
        double dpScoreTF,
        int seedOccLimit,
        MICDPOccurrence * dpOcc
    ) {
    // Seed DP - Only kick in when the DP is enabled
    // aka Deep DP
    
    #ifdef MIC_DP_DEBUG_PRINT_SEED_MATCHING
    printf("[SeedDp] MIC Deep DP initiated.\n");
    #endif

    // SeedDp - Initialisation of the round based result
    uint16_t peSeedOccCount = 0;
    uint8_t  peSeedStatus = MIC_OUTPUT_STATUS_OPEN;

    //char seedOverLimit = 0;
    unsigned long long saCount;

    //The following MUST be initialised as the down stream process might
    //make decision based on them.
    *readSeedSRAArgs->outputStatus   = MIC_OUTPUT_STATUS_OPEN; 
    *readSeedSRAArgs->occCount       = 0;
    *readSeedSRAArgs->metaCount      = 0;
    *mateSeedSRAArgs->outputStatus      = MIC_OUTPUT_STATUS_OPEN;
    *mateSeedSRAArgs->occCount          = 0;
    *mateSeedSRAArgs->metaCount         = 0;

    if (*readArgs->occCount > 0 ){
        readSeedSRAArgs = readArgs;
        //TODO
    } else {
        ////////////////////////////////////////////
        // Seed-Round READ
        ////////////////////////////////////////////
        readSeedSRAArgs->seedLength      = readArgs->seedLength;
        readSeedSRAArgs->readLength      = readArgs->readLength;
        readSeedSRAArgs->readCode        = readArgs->readCode;
        readSeedSRAArgs->readCode_Complt = readArgs->readCode_Complt;
        // Initialisation of the round based result
        uint16_t readLen                = readSeedSRAArgs->seedLength;
        // If seedLength or seedOverlap > 1, treat input as constant seedLength or seedOverlap
        // If seedLength or seedOverlap <= 1, treat input as a ratio to readLen + 0.5 to do rounding
        uint16_t readSeedLength      =  ( seedLength > 1.0 ? seedLength : (double)readLen * seedLength ) + 0.5 ;
        uint16_t readSeedUniqeLength =  ( readSeedLength - ( seedOverlap > 1.0 ? seedOverlap : (double)readLen * seedOverlap ) ) + 0.5 ;
        // Calling the model aligner
        saCount = MICSeedAlgnmtDoubleStrand(readSeedSRAArgs,&(cpPModels[readSeedLength]),&(cpNModels[readSeedLength]),
                readSeedLength,readSeedUniqeLength,seedOccLimit);

        if (*readSeedSRAArgs->outputStatus == MIC_OUTPUT_STATUS_OPEN) {
            *readSeedSRAArgs->outputStatus = MIC_OUTPUT_STATUS_COMPLETE;
        }

        #ifdef MIC_DP_DEBUG_PRINT_SEED_MATCHING
        printf("[SeedDp] READ SRA Found Results = %llu (FLAG %u)\n",saCount,*readSeedSRAArgs->outputStatus);
        #endif

        ////////////////////////////////////////////
        // Checkpoint : SRA Results Checking
        ////////////////////////////////////////////
        // If the SRA output buffer of the MIC has flooded
        // MIC is deemed to be incapable of handling the
        // seed alignment.
        if (*readSeedSRAArgs->outputStatus == MIC_OUTPUT_STATUS_CLOSE) {
            return MIC_DP_STATUS_TOO_MANY_RESULT;
        }

        // If either side of the SRA results in no
        // alignment result, cannot proceed
        if (saCount == 0) {
            return 0;
        }

        // sh : i believe it is already checked, so no need to check again here
        // If either side of the SRA results exceeds
        // the alignment result boundard, cannot proceed
        //if (seedOccLimit != -1 && *readSeedSRAArgs->occCount > seedOccLimit) {
        //    seedOverLimit = 1;
        //    *readSeedSRAArgs->occCount = seedOccLimit - 1;
        //}
    }

    if (*mateArgs->occCount > 0 ){
        mateSeedSRAArgs = mateArgs;
        //TODO
    } else {
        ////////////////////////////////////////////
        // Seed-Round MATE
        ////////////////////////////////////////////
        mateSeedSRAArgs->seedLength      = mateArgs->seedLength;
        mateSeedSRAArgs->readLength      = mateArgs->readLength;
        mateSeedSRAArgs->readCode        = mateArgs->readCode;
        mateSeedSRAArgs->readCode_Complt = mateArgs->readCode_Complt;
        // Initialisation of the round based result
        uint16_t mateLen                    = mateSeedSRAArgs->seedLength;
        // If seedLength or seedOverlap > 1, treat input as constant seedLength or seedOverlap
        // If seedLength or seedOverlap <= 1, treat input as a ratio to mateLen; + 0.5 to do rounding
        uint16_t mateSeedLength             =  ( seedLength > 1.0 ? seedLength : (double)mateLen * seedLength ) + 0.5 ;
        uint16_t mateSeedUniqeLength        =  ( mateSeedLength - ( seedOverlap > 1.0 ? seedOverlap :  (double)mateLen * seedOverlap ) ) + 0.5 ;

        // Calling the model aligner
        saCount = MICSeedAlgnmtDoubleStrand(mateSeedSRAArgs,&(cpPModels[mateSeedLength]),&(cpNModels[mateSeedLength]),
                mateSeedLength,mateSeedUniqeLength,seedOccLimit);

        if (*mateSeedSRAArgs->outputStatus == MIC_OUTPUT_STATUS_OPEN) {
            *mateSeedSRAArgs->outputStatus = MIC_OUTPUT_STATUS_COMPLETE;
        }

        #ifdef MIC_DP_DEBUG_PRINT_SEED_MATCHING
        printf("[SeedDp] MATE SRA Found Results = %llu (FLAG %u)\n",saCount,*mateSeedSRAArgs->outputStatus);
        #endif

        ////////////////////////////////////////////
        // Checkpoint : SRA Results Checking
        //////////////////////////////////////////// 
        // If the SRA output buffer of the MIC has flooded
        // MIC is deemed to be incapable of handling the
        // seed alignment.
        if (*mateSeedSRAArgs->outputStatus == MIC_OUTPUT_STATUS_CLOSE) {
            return MIC_DP_STATUS_TOO_MANY_RESULT;
        }

        // If either side of the SRA results in no
        // alignment result, cannot proceed
        if (saCount == 0) {
            return 0;
        }

        // If either side of the SRA results exceeds
        // the alignment result boundard, cannot proceed
        //if (seedOccLimit != -1 && *mateSeedSRAArgs->occCount > seedOccLimit) {
        //    seedOverLimit = 1;
        //    *mateSeedSRAArgs->occCount = seedOccLimit - 1;
        //}
    }

    ////////////////////////////////////////////
    // Seed-Round PE
    ////////////////////////////////////////////

    MICPEArgumentsConfig(peSeedArgs, readSeedSRAArgs, mateSeedSRAArgs, MIC_PE_MAX_RESULT,
            MIC_PE_MAX_RESULT, peOutputBuffer, peMetaBuffer, &peSeedOccCount, &peSeedStatus);

    MICPEArgumentsSetMaxOutput(peSeedArgs, MIC_DEEP_DP_SEED_PAIR_LIMIT);
    MICPEMappingInitialise(peSeedArgs);
    MICPEMappingOccurrences(peSeedArgs);
    MICPEMappingComplete(peSeedArgs);

    #ifdef MIC_DP_DEBUG_PRINT_SEED_MATCHING
    printf("[SeedDp] Pair-ing Engine Found Result Flag = %u\n",*peSeedArgs->outputStatus);
    #endif
    
    // proceeed only if paring success
    if (*(peSeedArgs->outputStatus) != MIC_PE_OUTPUT_STATUS_PAIR) {
       return 0;
    }

    ////////////////////////////////////////////
    // Seed-Round DP
    ////////////////////////////////////////////

    return MICDPDeepAlignment(dpOcc, dpWork, peSeedArgs, hsp);

}
