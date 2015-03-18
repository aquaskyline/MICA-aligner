/*
 *
 *    
 *      MicMappingQuality.c
 *    Soap3(gpu)
 *
 *    Copyright (C) 2011, HKU
 *
 *    This program is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License
 *    as published by the Free Software Foundation; either version 2
 *    of the License, or (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */

#include "MicMappingQuality.h"

__attribute__((target(mic)))
static void MicCalculateSraResult(
        MICSRAArguments * sraArgs,
        CPTSRAModel * cpPModels, CPTSRAModel * cpNModels,
        int * numBestOutput, int * num2ndBestOutput) {
    int NO_ENTRY = 999;

    int size = *sraArgs->occCount;
    int numBest = 0;
    int bestMismatch = NO_ENTRY;
    int numSecondBest = 0;
    int secondBestMismatch = NO_ENTRY;

    int i;
    for (i = 0; i < size; i++) {
        int currentMismatch = sraArgs->metaBlock[i].numOfErr;
        if (currentMismatch < bestMismatch) {
            // New best is found, push the ranking down
            // Original best becomes the second best
            secondBestMismatch = bestMismatch;
            numSecondBest = numBest;

            // New best
            numBest = 1;
            bestMismatch = currentMismatch;
        } else if (currentMismatch == bestMismatch) {
            numBest++;
        } else if (currentMismatch < secondBestMismatch) {
            numSecondBest = 1;
            secondBestMismatch = currentMismatch;
        } else if (currentMismatch == secondBestMismatch) {
            numSecondBest++;
        }
    }

    // No second best is found, we need to perform oneMoreMismatch
    if (secondBestMismatch == NO_ENTRY) {
        MICProcessReadDoubleStrand(sraArgs, cpPModels,
                cpNModels);

        numSecondBest = *sraArgs->occCount - size;
        if (numSecondBest < 0) {
            numSecondBest = 0;
        }
        // TODO: Handle the case where oneMoreMismatch may give too many results
    }

    (* numBestOutput) = numBest;
    (* num2ndBestOutput) = numSecondBest;
}


__attribute__((target(mic)))
static void MicCalculatePeOptimal(
        MICPEArguments * peArgs,
        DPScores * dpScores,
        int * optimalScoreOutput,
        int * numOptimalOutput,
        int * subOptimalScoreOutput,
        int * numSubOptimalOutput) {

    int NO_ENTRY = -99999;

    int size = *peArgs->occCount;
    int numBest = 0;
    int bestScore = NO_ENTRY;
    int secondBestScore = NO_ENTRY;

    int readMateTotalLegnth = peArgs->readArgs->seedLength +
            peArgs->mateArgs->seedLength;

    int i;
    for (i = 0; i < size; i += 2) {
        int currentError = peArgs->outputMeta[i].numOfErr +
                peArgs->outputMeta[i + 1].numOfErr;
        int currentScore = readMateTotalLegnth -
                currentError * (dpScores->dpMismatch - dpScores->dpMatch);
        if (currentScore > bestScore) {
            // New best is found, push the ranking down
            // Original best becomes the second best
            secondBestScore = bestScore;

            // New best
            numBest = 1;
            bestScore = currentScore;
        } else if (currentScore == bestScore) {
            numBest++;
        } else if (currentScore > secondBestScore) {
            secondBestScore = currentScore;
        }
    }
    (* numOptimalOutput) = numBest;
    (* optimalScoreOutput) = bestScore;
    (* numSubOptimalOutput) = size / 2 - numBest;
    (* subOptimalScoreOutput) = secondBestScore;


}
// A more generic function for getting the num of best
// and second best alignment from DP, which can be used by
// Default DP as well as Deep DP
__attribute__((target(mic)))
static void _MicCalculateDpResult(
        MICDPOccurrence * dpOccurrences, uint32_t count,        
        int * numBestOutput, int * num2ndBestOutput,
        uint32_t start, uint32_t step) {

    int NO_ENTRY = -99999;

    int numBest = 0;
    int bestScore = NO_ENTRY;

    int i;
    for (i = start; i < count; i += step) {
        int currentScore = dpOccurrences[i].score;
        if (currentScore > bestScore) {
            // New best
            numBest = 1;
            bestScore = currentScore;
        } else if (currentScore == bestScore) {
            numBest++;
        }
    }

    (* numBestOutput) = numBest;
    (* num2ndBestOutput) = count - numBest;
}

// This function gives number of best alignment and 
// number of second best alignment from DP results,
// without considering the mixing of read and mate DP results
__attribute__((target(mic)))
static void MicCalculateDpResult(
        MICDPOccurrence * dpOccurrences, uint32_t count,        
        int * numBestOutput, int * num2ndBestOutput) {

    _MicCalculateDpResult(dpOccurrences, count,        
            numBestOutput, num2ndBestOutput, 0, 1);

}


__attribute__((target(mic)))
static void MicCalculateDeepDpReadResult(
        MICDPOccurrence * dpOccurrences, uint32_t count,        
        int * numBestOutput, int * num2ndBestOutput) {

    _MicCalculateDpResult(dpOccurrences, count,        
            numBestOutput, num2ndBestOutput, 0 /* start */, 2 /* step */);

}

__attribute__((target(mic)))
static void MicCalculateDeepDpMateResult(
        MICDPOccurrence * dpOccurrences, uint32_t count,        
        int * numBestOutput, int * num2ndBestOutput) {

    _MicCalculateDpResult(dpOccurrences, count,        
            numBestOutput, num2ndBestOutput, 1 /* start */, 2 /* step */);
}

// Jeanno: Currently the result given by (New) Default DP are stored
// in MICPEArguments (Sra Side result) and in MICDPOccurrence array
// (DP side result). Looking forward to redesigning it in the future.
__attribute__((target(mic)))
static void MicCalculateDpOptimal(
        MICPEArguments * peArgs,
        MICSRAArguments * sraArgs,            // The base side sraArg
        MICDPOccurrence * dpOccurrences,
        DPScores * dpScores,
        int * optimalScoreOutput,
        int * numOptimalOutput,
        int * subOptimalScoreOutput,
        int * numSubOptimalOutput) {

    int NO_ENTRY = -99999;
    int count = *peArgs->occCount;
    int numBest = 0;
    int bestScore = NO_ENTRY;
    int secondBestScore = NO_ENTRY;

    int i;
    for (i = 0; i < count; i++) {
        int currentError = peArgs->outputMeta[i].numOfErr;
        int currentScore = sraArgs->seedLength +
                currentError * (dpScores->dpMismatch - dpScores->dpMatch) +
                dpOccurrences[i].score;
        if (currentScore > bestScore) {
            // New best is found, push the ranking down
            // Original best becomes the second best
            secondBestScore = bestScore;

            // New best
            numBest = 1;
            bestScore = currentScore;
        } else if (currentScore == bestScore) {
            numBest++;
        } else if (currentScore > secondBestScore) {
            secondBestScore = currentScore;
        }
    }
    (* numOptimalOutput) = numBest;
    (* optimalScoreOutput) = bestScore;
    (* numSubOptimalOutput) = count - numBest;
    (* subOptimalScoreOutput) = secondBestScore;
}

__attribute__((target(mic)))
static void MicCalculateDeepDpOptimal(
        MICDPOccurrence * dpOccurrences,
        uint32_t count,
        int * optimalScoreOutput,
        int * numOptimalOutput,
        int * subOptimalScoreOutput,
        int * numSubOptimalOutput) {

    int NO_ENTRY = -99999;
    int numBest = 0;
    int bestScore = NO_ENTRY;
    int secondBestScore = NO_ENTRY;

    int i;
    for (i = 0; i < count; i += 2) {
        int currentScore = dpOccurrences[i].score + dpOccurrences[i + 1].score;
        if (currentScore > bestScore) {
            // New best is found, push the ranking down
            // Original best becomes the second best
            secondBestScore = bestScore;

            // New best
            numBest = 1;
            bestScore = currentScore;
        } else if (currentScore == bestScore) {
            numBest++;
        } else if (currentScore > secondBestScore) {
            secondBestScore = currentScore;
        }
    }
    (* numOptimalOutput) = numBest;
    (* optimalScoreOutput) = bestScore;
    (* numSubOptimalOutput) = count - numBest;
    (* subOptimalScoreOutput) = secondBestScore;
}


// This is used to get optimal scores after performing mix base DP
__attribute__((target(mic)))
static void MicCalculateMixDpOptimal(
        MICPEArguments * peArgs,
        MICDPOccurrence * dpOccurrences,
        DPScores * dpScores,
        int * optimalScoreOutput,
        int * numOptimalOutput,
        int * subOptimalScoreOutput,
        int * numSubOptimalOutput) {

    int size = *peArgs->occCount;
    // Find the starting point of Mate base
    int mateStart;
    for (mateStart = 0; mateStart < size; mateStart++) {
        MICDPOccurrence * dpOcc = dpOccurrences + mateStart;
        if (dpOccurrences->ambPosition == 0 &&
                dpOccurrences->matchLen == 0) {
            break;
        }
    }

    if (mateStart == size) {
        // In this case readBase result doesn't exist
        // We can then treat it as single-based mate based DP result
        MicCalculateDpOptimal(peArgs, peArgs->mateArgs, dpOccurrences, dpScores,
                optimalScoreOutput, numOptimalOutput,
                subOptimalScoreOutput, numSubOptimalOutput);
        return;
    }
    mateStart++;
    // dpOccurrences[mateStart] should be at the start of mate base results here

    // Prepare peArgs for read and mate separately
    // to mimic the argument for one side base DP optimal calculation
    MICPEArguments readBasePeArgs = *peArgs;
    uint16_t readBaseOcc = mateStart - 1;
    readBasePeArgs.occCount = &readBaseOcc; 
    MICPEArguments mateBasePeArgs = *peArgs;
    mateBasePeArgs.output += mateStart;
    mateBasePeArgs.outputMeta += mateStart;
    uint16_t mateBaseOcc = *peArgs->occCount - readBaseOcc;
    mateBasePeArgs.occCount = &mateBaseOcc; 


    // Perform a read base Dp Optimal calculation
    int readBasePeOptimalScore, readBasePeNumOptimal,
            readBasePeSubOptScore, readBasePeNumSubOpt;
    MicCalculateDpOptimal(&readBasePeArgs, readBasePeArgs.readArgs,
            dpOccurrences, dpScores,
            &readBasePeOptimalScore, &readBasePeNumOptimal,
            &readBasePeSubOptScore, &readBasePeNumSubOpt);

    // Perform a mate base Dp Optimal calculation
    int mateBasePeOptimalScore, mateBasePeNumOptimal,
            mateBasePeSubOptScore, mateBasePeNumSubOpt;
    MicCalculateDpOptimal(&mateBasePeArgs, mateBasePeArgs.mateArgs,
            dpOccurrences, dpScores,
            &mateBasePeOptimalScore, &mateBasePeNumOptimal,
            &mateBasePeSubOptScore, &mateBasePeNumSubOpt);

    // Calculate the final result
    // Combining the separated result from readBase and mateBase DP
    // into the final opt and subopt result
    (* optimalScoreOutput) = readBasePeOptimalScore;
    (* numOptimalOutput) = readBasePeNumOptimal;
    (* subOptimalScoreOutput) = readBasePeSubOptScore;
    (* numSubOptimalOutput) = readBasePeNumSubOpt;

    // Merge mateBaseOptimal result into final output
    if (mateBasePeOptimalScore > readBasePeOptimalScore) {
        (* subOptimalScoreOutput) = readBasePeOptimalScore;
        (* numSubOptimalOutput) = readBasePeNumOptimal;

        (* optimalScoreOutput) = mateBasePeOptimalScore;
        (* numOptimalOutput) = mateBasePeNumOptimal;
    } else {
        if (mateBasePeOptimalScore == readBasePeOptimalScore) {
            (* numOptimalOutput) += mateBasePeNumOptimal;
        } else {
            if (mateBasePeOptimalScore > readBasePeSubOptScore) {
                (* subOptimalScoreOutput) = mateBasePeOptimalScore;
                (* numSubOptimalOutput) = mateBasePeNumOptimal;
            } else if (mateBasePeOptimalScore == readBasePeSubOptScore) {
                (* numSubOptimalOutput) += mateBasePeNumOptimal;
            }
        }
    }

    // Merge mateBaseSubOptimal result into final output
    // Since mateBaseSubOpt must not be opt in the final result,
    // we can safely skip comparing it to the current opt score
    if (mateBasePeSubOptScore > (* subOptimalScoreOutput)) {
        (* subOptimalScoreOutput) = mateBasePeSubOptScore;
        (* numSubOptimalOutput) = mateBasePeNumSubOpt;
    } else if (mateBasePeSubOptScore == (* subOptimalScoreOutput)) {
        (* numSubOptimalOutput) += mateBasePeNumSubOpt;
    }

}

// g_log_n: this should be initiate by bwase_initialize
//      the size of the array is 256
__attribute__((target(mic)))
PEMappingQuality MicCalculatePEMappingQuality(
        MICSRAArguments * readSraArgs, 
        MICSRAArguments * mateSraArgs,
        MICPEArguments * peArgs,
        MICDPOccurrence * dpOccurrences,
        CPTSRAModel * cpPModels,
        CPTSRAModel * cpNModels,
        int * g_log_n) {
    DPScores * dpScores = peArgs->dpScores;
    PEMappingQuality ret;
    ret.readQuality = 0;
    ret.mateQuality = 0;

    if (*peArgs->outputStatus == MIC_PE_OUTPUT_STATUS_PAIR) {

        // Find readNumBest & readNum2ndBest
        int readNumBest, readNum2ndBest;
        MicCalculateSraResult(readSraArgs, cpPModels, cpNModels,
            &readNumBest, &readNum2ndBest);
        // Find mateNumBest & mateNum2ndBest
        int mateNumBest, mateNum2ndBest;
        MicCalculateSraResult(mateSraArgs, cpPModels, cpNModels,
            &mateNumBest, &mateNum2ndBest);

        // Find peOptimalScore & peNumOptimal
        // Find peSubOptScore & peNumSubOpt
        int peOptimalScore, peNumOptimal, peSubOptScore, peNumSubOpt;
        MicCalculatePeOptimal(peArgs, dpScores, &peOptimalScore, &peNumOptimal,
                &peSubOptScore, &peNumSubOpt);

        // Plug the parameter needed for the scores
        bwaLikePairQualScore(readNumBest, readNum2ndBest, mateNumBest, mateNum2ndBest,
                g_log_n, peOptimalScore, peNumOptimal, peSubOptScore, peNumSubOpt,
                readSraArgs->seedLength, mateSraArgs->seedLength,
                &ret.readQuality, &ret.mateQuality);
    } else if (*peArgs->outputStatus == MIC_PE_OUTPUT_STATUS_DP_BASE_READ) {
        // Read alignment is from SRA
        // Mate alignment is from DP

        // Find readNumBest & readNum2ndBest
        int readNumBest, readNum2ndBest;
        MicCalculateSraResult(readSraArgs, cpPModels, cpNModels,
            &readNumBest, &readNum2ndBest);

        // Find mateNumBest & mateNum2ndBest
        int mateNumBest, mateNum2ndBest;
        MicCalculateDpResult(dpOccurrences, *peArgs->occCount,
                &mateNumBest, &mateNum2ndBest);


        // Find peOptimalScore & peNumOptimal
        // Find peSubOptScore & peNumSubOpt
        int peOptimalScore, peNumOptimal, peSubOptScore, peNumSubOpt;
        MicCalculateDpOptimal(peArgs, readSraArgs, dpOccurrences, dpScores,
                &peOptimalScore, &peNumOptimal,
                &peSubOptScore, &peNumSubOpt);

        // Plug the parameter needed for the scores
        bwaLikePairQualScore(readNumBest, readNum2ndBest, mateNumBest, mateNum2ndBest,
                g_log_n, peOptimalScore, peNumOptimal, peSubOptScore, peNumSubOpt,
                readSraArgs->seedLength, mateSraArgs->seedLength,
                &ret.readQuality, &ret.mateQuality);
    } else if (*peArgs->outputStatus == MIC_PE_OUTPUT_STATUS_DP_BASE_MATE) {
        // Read alignment is from DP
        // Mate alignment is from SRA

        // Find readNumBest & readNum2ndBest
        int readNumBest, readNum2ndBest;
        MicCalculateDpResult(dpOccurrences, *peArgs->occCount,
                &readNumBest, &readNum2ndBest);

        // Find mateNumBest & mateNum2ndBest
        int mateNumBest, mateNum2ndBest;
        MicCalculateSraResult(mateSraArgs, cpPModels, cpNModels,
            &mateNumBest, &mateNum2ndBest);

        // Find peOptimalScore & peNumOptimal
        // Find peSubOptScore & peNumSubOpt
        int peOptimalScore, peNumOptimal, peSubOptScore, peNumSubOpt;
        MicCalculateDpOptimal(peArgs, mateSraArgs, dpOccurrences, dpScores,
                &peOptimalScore, &peNumOptimal,
                &peSubOptScore, &peNumSubOpt);

        // Plug the parameter needed for the scores
        bwaLikePairQualScore(readNumBest, readNum2ndBest, mateNumBest, mateNum2ndBest,
                g_log_n, peOptimalScore, peNumOptimal, peSubOptScore, peNumSubOpt,
                readSraArgs->seedLength, mateSraArgs->seedLength,
                &ret.readQuality, &ret.mateQuality);
    } else if (*peArgs->outputStatus == MIC_PE_OUTPUT_STATUS_DP_MIX_BASE) {
        // Find readNumBest & readNum2ndBest
        int readNumBest, readNum2ndBest;
        MicCalculateSraResult(readSraArgs, cpPModels, cpNModels,
            &readNumBest, &readNum2ndBest);
        // Find mateNumBest & mateNum2ndBest
        int mateNumBest, mateNum2ndBest;
        MicCalculateSraResult(mateSraArgs, cpPModels, cpNModels,
            &mateNumBest, &mateNum2ndBest);

        // Find peOptimalScore & peNumOptimal
        // Find peSubOptScore & peNumSubOpt
        int peOptimalScore, peNumOptimal, peSubOptScore, peNumSubOpt;
        MicCalculateMixDpOptimal(peArgs, dpOccurrences, dpScores,
                &peOptimalScore, &peNumOptimal,
                 &peSubOptScore, &peNumSubOpt);
    } else if (*peArgs->outputStatus == MIC_PE_OUTPUT_STATUS_DP_SEED) {
        // // Case of Deep DP
        // // Find readNumBest & readNum2ndBest
        int readNumBest, readNum2ndBest;
        MicCalculateDeepDpReadResult(dpOccurrences, *peArgs->output,
                &readNumBest, &readNum2ndBest);

        // Find mateNumBest & mateNum2ndBest
        int mateNumBest, mateNum2ndBest;
        MicCalculateDeepDpMateResult(dpOccurrences, *peArgs->output,
                &mateNumBest, &mateNum2ndBest);
        
        // Find peOptimalScore & peNumOptimal
        // Find peSubOptScore & peNumSubOpt
        int peOptimalScore, peNumOptimal, peSubOptScore, peNumSubOpt;
        MicCalculateDeepDpOptimal(dpOccurrences, *peArgs->occCount,
                &peOptimalScore, &peNumOptimal,
                &peSubOptScore, &peNumSubOpt);

        // Plug the parameter needed for the scores
        bwaLikePairQualScore(readNumBest, readNum2ndBest, mateNumBest, mateNum2ndBest,
                g_log_n, peOptimalScore, peNumOptimal, peSubOptScore, peNumSubOpt,
                readSraArgs->seedLength, mateSraArgs->seedLength,
                &ret.readQuality, &ret.mateQuality);
    }

    return ret;
}
