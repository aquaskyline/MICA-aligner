//
//    MICCT-ShortHandler.c
//
//    mica
//
//    Copyright (C) 2014, HKU
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

#include "MICCT-ShortHandler.h"

__attribute__((target(mic)))
void MICSHProcessRead(HSP * hsp,
                     MICSRAArguments * readMicArgs, MICSRAArguments * mateMicArgs,
                     CPTSRAModel * cpPModels, CPTSRAModel * cpNModels,
                     PEInput * peAlgnmtInput, MICPEArguments * peArgs,
                     
                     MICSRAArguments * readSeedMicArgs, MICSRAArguments * mateSeedMicArgs,
                     CPTSRAModel * cpPModels_seed, CPTSRAModel * cpNModels_seed,
                     MICPEArguments * peSeedArgs,
                     unsigned int * peSeedOutputPtr, MICSRAOccMetadata * peSeedMetaPtr,
                     
                     PEDPSetting * pedpSetting, DPWorkMIC * dpWork,
                     
                     unsigned int * outputPtr, uint8_t * outputBufferStatus, uint16_t * occCount,
                     MICDPOccurrence * dpOcc, unsigned dpOccCount, uint32_t * dpOccCountInc) {
    
    int caseIdx = 0;
    int stageCaseIdx = 0;
    while (1) {
    
        ////////////////////////////////////////////
        // READ
        ////////////////////////////////////////////
        stageCaseIdx = caseIdx;
        // Calling the model aligner
        MICProcessStageDoubleStrand(readMicArgs,
                                    &(cpPModels[readMicArgs->readLength]),
                                    &(cpNModels[readMicArgs->readLength]),
                                    &stageCaseIdx);
        
        ////////////////////////////////////////////
        // MATE
        ////////////////////////////////////////////
        stageCaseIdx = caseIdx;
        // Calling the model aligner
        MICProcessStageDoubleStrand(mateMicArgs,
                                    &(cpPModels[mateMicArgs->readLength]),
                                    &(cpNModels[mateMicArgs->readLength]),
                                    &stageCaseIdx);
        
        // Case that either of the leg alignment has been flooded
        if ((*readMicArgs->outputStatus) != MIC_OUTPUT_STATUS_OPEN || 
            (*mateMicArgs->outputStatus) != MIC_OUTPUT_STATUS_OPEN) {
            break;
        }

        // All-best alignment
        // This placed after the last step check because we want to carry out the last PE
        // outside of this loop. The next call to the PE engine do not check the SRA occ count
        // hence it will carry out singleEndAlignmentCopy.
        if (peAlgnmtInput->OutputType == PE_REPORT_ALL_BEST) {
            MICPEMappingInitialise(peArgs);
            MICPEMappingOccurrences(peArgs);
            
            if ( (*outputBufferStatus) == MIC_PE_OUTPUT_STATUS_SKIPPED ||
                 (*outputBufferStatus) == MIC_PE_OUTPUT_STATUS_UNHANDLE ||
                 (*outputBufferStatus) == MIC_PE_OUTPUT_STATUS_CLOSED ||
                 (*outputBufferStatus) == MIC_PE_OUTPUT_STATUS_PAIR ) {
                break;
            }
            
        }
        
        // Case that all the cases have been completed.
        if (cpPModels[readMicArgs->readLength].cases[stageCaseIdx].type==SRA_CASE_TYPE_NOT_INITALISED) {
            break;
        }
        
        caseIdx = stageCaseIdx + 1;
    }
    
    if ((*readMicArgs->outputStatus) == MIC_OUTPUT_STATUS_OPEN) {
        (*readMicArgs->outputStatus) = MIC_OUTPUT_STATUS_COMPLETE;
    }
    if ((*mateMicArgs->outputStatus) == MIC_OUTPUT_STATUS_OPEN) {
        (*mateMicArgs->outputStatus) = MIC_OUTPUT_STATUS_COMPLETE;
    }
    
    ////////////////////////////////////////////
    // Pairing 
    ////////////////////////////////////////////
    // This is where we will call the PE module.
    // We need to carry out the following here
    // 1. Check that the READ and MATE alignment
    //    count does not exceed a defined boundary
    // 2. Pipe READ and MATE alignments into the
    //    PEMappingCore
    // 3. Restrict the PEMappingCore to report only
    //    if the number of pair-end alignment is
    //    less than some defined boundary.
    // 4. If orphan => Default DP. Symmetric
    // 5. If still no PE result, seed and deep DP
    
    // All-valid alignment
    if (peAlgnmtInput->OutputType == PE_REPORT_ALL) {
        MICPEMappingInitialise(peArgs);
        MICPEMappingOccurrences(peArgs);
    }
    
    MICPEMappingComplete(peArgs);
    
    ////////////////////////////////////////////
    // DP Engine
    ////////////////////////////////////////////

    // skip if no vacancy
    if ( dpOccCount+ MIC_DP_OUTPUT_MAX_ALIGNMENT > MIC_DP_OUTPUT_SIZE_PER_THREAD ) {
         (*outputBufferStatus) = MIC_PE_OUTPUT_STATUS_CLOSED;
         (*occCount) = 0;
    }

    (*dpOccCountInc) = 0;

    // Recovery DP - Only kick in when the DP is enabled
    // aka Default DP
    if ( (pedpSetting->SGAOrphanEnhancement) &&
            ((*outputBufferStatus)==MIC_PE_OUTPUT_STATUS_READ_NO_ALIGNMENT ||
            (*outputBufferStatus)==MIC_PE_OUTPUT_STATUS_MATE_NO_ALIGNMENT ||
            (*outputBufferStatus)==MIC_PE_OUTPUT_STATUS_BAD_PAIR )
    ) {
        
        int defaultDPCount = MICDPOrphanAlignment(dpOcc + dpOccCount, dpWork, peArgs, hsp, pedpSetting->SGAOrphanTriggerTF);

        #ifdef MIC_DP_DEBUG_PRINT_ALIGNMENT_FLOW
            printf("[MICCT] Recovery DP is invoked. Returned %d\n",defaultDPCount);
        #endif
        
        if (defaultDPCount > 0){
            if ((*outputBufferStatus) == MIC_PE_OUTPUT_STATUS_READ_NO_ALIGNMENT){
                (*outputBufferStatus) = MIC_PE_OUTPUT_STATUS_DP_BASE_MATE;
            } else if ((*outputBufferStatus) == MIC_PE_OUTPUT_STATUS_MATE_NO_ALIGNMENT){
                (*outputBufferStatus) = MIC_PE_OUTPUT_STATUS_DP_BASE_READ;
            } else if ((*outputBufferStatus) == MIC_PE_OUTPUT_STATUS_BAD_PAIR){
                (*outputBufferStatus) = MIC_PE_OUTPUT_STATUS_DP_MIX_BASE;
            }
            (*dpOccCountInc) = defaultDPCount;
            (*occCount) = defaultDPCount;
        } else if (defaultDPCount == MIC_DP_STATUS_BREACH_LIMIT) {
            // Too many DP results reported, requires CPU assistance
            (*outputBufferStatus) = MIC_PE_OUTPUT_STATUS_CLOSED;
            (*occCount) = 0;
        } else if (defaultDPCount == 0 || defaultDPCount == MIC_DP_STATUS_TOO_MANY_RESULT){
            // processed but no result or too many result ... do nothing
            //(*outputBufferStatus) = MIC_PE_OUTPUT_STATUS_CLOSED;
            //(*occCount) = 0;
        } else {
            // unprocessed or error
            (*outputBufferStatus) = MIC_PE_OUTPUT_STATUS_CLOSED;
            (*occCount) = 0;
        }

    }

    // Seed DP - Only kick in when the DP is enabled
    // aka Deep DP
    if ( (pedpSetting->SGASeedEnhancement) &&
            ((*outputBufferStatus) == MIC_PE_OUTPUT_STATUS_READ_NO_ALIGNMENT ||
            (*outputBufferStatus) == MIC_PE_OUTPUT_STATUS_MATE_NO_ALIGNMENT ||
            (*outputBufferStatus) == MIC_PE_OUTPUT_STATUS_BOTH_NO_ALIGNMENT)
    ) {
        int deepDPCount = MICDeepDP(peSeedArgs, 
                peSeedOutputPtr, peSeedMetaPtr, 
                readMicArgs, mateMicArgs,
                readSeedMicArgs, mateSeedMicArgs,
                cpPModels_seed, cpNModels_seed,
                hsp,
                dpWork,
                pedpSetting->SGASeedLength,
                pedpSetting->SGASeedLengthOverlap,
                pedpSetting->SGAOrphanTriggerTF,
                pedpSetting->SGAScoreTF,
                pedpSetting->SGASeedOccLimitation,
                dpOcc + dpOccCount);

        #ifdef MIC_DP_DEBUG_PRINT_ALIGNMENT_FLOW
            printf("[MICCT] Deep DP is invoked. Returned %d\n",deepDPCount);
        #endif
        
        if (deepDPCount > 0) {
            (*outputBufferStatus) = MIC_PE_OUTPUT_STATUS_DP_SEED;
            (*dpOccCountInc) = deepDPCount;
            (*occCount) = 1;
            *outputPtr = deepDPCount;
        } else if (deepDPCount == MIC_DP_STATUS_BREACH_LIMIT) {
            // Too many DP results reported, requires CPU assistance
            (*outputBufferStatus) = MIC_PE_OUTPUT_STATUS_CLOSED;
            (*occCount) = 0;
        } else if (pedpSetting->SGASeedEnhancement_1) {
            
            deepDPCount = MICDeepDP(peSeedArgs, 
                    peSeedOutputPtr, peSeedMetaPtr, 
                    readMicArgs, mateMicArgs,
                    readSeedMicArgs, mateSeedMicArgs,
                    cpPModels_seed, cpNModels_seed,
                    hsp,
                    dpWork,
                    pedpSetting->SGASeedLength_1,
                    pedpSetting->SGASeedLengthOverlap_1,
                    pedpSetting->SGAOrphanTriggerTF,
                    pedpSetting->SGAScoreTF,
                    pedpSetting->SGASeedOccLimitation_1,
                    dpOcc + dpOccCount);
                    
            #ifdef MIC_DP_DEBUG_PRINT_ALIGNMENT_FLOW
                printf("[MICCT] Deep DP (Level-2) is invoked. Returned %d\n",deepDPCount);
            #endif
            
            if (deepDPCount > 0) {
                (*outputBufferStatus) = MIC_PE_OUTPUT_STATUS_DP_SEED;
                (*dpOccCountInc) = deepDPCount;
                (*occCount) = 1;
                *outputPtr = deepDPCount;
            } else if (deepDPCount == MIC_DP_STATUS_BREACH_LIMIT) {
                // Too many DP results reported, requires CPU assistance
                (*outputBufferStatus) = MIC_PE_OUTPUT_STATUS_CLOSED;
                (*occCount) = 0;
            } else if (deepDPCount == 0 || deepDPCount == MIC_DP_STATUS_TOO_MANY_RESULT){
                // processed but no result or too many result ... do nothing
                //(*outputBufferStatus) = MIC_PE_OUTPUT_STATUS_SKIPPED;
                //(*occCount) = 0;
            } else {
                // unprocessed or error
                (*outputBufferStatus) = MIC_PE_OUTPUT_STATUS_CLOSED;
                (*occCount) = 0;
            }
        } else if (deepDPCount == 0 || deepDPCount == MIC_DP_STATUS_TOO_MANY_RESULT){
            // processed but no result or too many result ... do nothing
            //(*outputBufferStatus) = MIC_PE_OUTPUT_STATUS_SKIPPED;
            //(*occCount) = 0;
        } else {
            // unprocessed or error
            (*outputBufferStatus) = MIC_PE_OUTPUT_STATUS_CLOSED;
            (*occCount) = 0;
        }
    }
}
