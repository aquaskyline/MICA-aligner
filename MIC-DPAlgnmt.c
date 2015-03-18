//
//    MIC-DPAlgnmt.c
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

#include "MIC-DPAlgnmt.h"

__attribute__((target(mic)))
char * MICDPExtractCodeFromHSP(char * buffer, const unsigned startIndex, const unsigned length, const HSP * hsp){
    unsigned int * packedDNA = hsp->packedDNA;
    unsigned int packedDNAStartBlock =  startIndex / CHAR_PER_WORD;
    unsigned int packedDNAEndBlock =  (startIndex + length + (CHAR_PER_WORD - 1)) / CHAR_PER_WORD - 1;
    unsigned int packedDNABlockIdx = 0;
    int expandRegionIdx = 0;
    
    int buffer_i;
    unsigned int word;
    int charMask = (1<<BIT_PER_CHAR)-1;
    
    // Read in the reference sequence into buffer.    
    for (packedDNABlockIdx=packedDNAStartBlock;packedDNABlockIdx<=packedDNAEndBlock;packedDNABlockIdx++) {
        word = packedDNA[packedDNABlockIdx];
        for (buffer_i=0; buffer_i<CHAR_PER_WORD; buffer_i++) {
            buffer[expandRegionIdx++] = (unsigned char)((word>>((CHAR_PER_WORD-buffer_i-1)*BIT_PER_CHAR)) & charMask);
        }
    } 
    return  buffer + startIndex % CHAR_PER_WORD;
}


__attribute__((target(mic)))
int _MICDPPopulateOccurrence(MICDPOccurrence * dpoToWrite, DPWorkMIC * dpWork, unsigned int startPos, int dpStrand) {

    int j;
    dpoToWrite -> ambPosition = startPos + dpWork->maxScoreStartPos;
    dpoToWrite -> strand = dpStrand;
    dpoToWrite -> matchElemsCount = 0;
   

    //NOTE: MIC_TRACE_I/D refers to REF insert/delete, SRA_CHAR_ERROR_TYPE is the opposite
 
    int currentMatchType = -232;
    int matchType;

    //for (j = dpWork->traceLength-1; j>=0; j--){
    for (j =0; j<= dpWork->traceLength-1;  j++){

        if (dpWork->trace[j]==MIC_DP_TRACE_M) {
            matchType = SRA_CHAR_ERROR_TYPE_MISMATCH;
        } else if (dpWork->trace[j]==MIC_DP_TRACE_I) {
            matchType = SRA_CHAR_ERROR_TYPE_DELETE;
        } else if (dpWork->trace[j]==MIC_DP_TRACE_D) {
            matchType = SRA_CHAR_ERROR_TYPE_INSERT;
        } else if (dpWork->trace[j]==MIC_DP_TRACE_S) {
            matchType = SRA_CHAR_ERROR_TYPE_SOFTCLIP;
        } else {
            //ATTENTION
            printf("[SAFE-GUARD] Unable Matching Type %u\n",dpWork->trace[j]);
            return 0;
        }
        
        if (matchType != currentMatchType){
        
            //if (dpoToWrite->matchElemsCount >= MIC_DP_MAX_ERROR_COUNT) {
            //    return 0;
            //}
            currentMatchType = matchType; 
            dpoToWrite->matchElems[dpoToWrite->matchElemsCount].type = matchType;
            dpoToWrite->matchElems[dpoToWrite->matchElemsCount].length = 0;
            dpoToWrite->matchElemsCount++;
        }
        
        dpoToWrite->matchElems[dpoToWrite->matchElemsCount-1].length++;
    }
    dpoToWrite -> matchLen = dpWork->maxScoreEndPos - dpWork->maxScoreStartPos + 1;
    dpoToWrite->score = dpWork->maxScore;
    return 1;
}


__attribute__((target(mic)))
static int _MICDPOrphanRecovery(MICDPOccurrence * dpOutputPtr, int dpOccCount, DPWorkMIC * dpWork, const MICPEArguments * peArgs, 
          const HSP * hsp,  MICSRAArguments * mappedRead, MICSRAArguments * unmappedRead) {

    DPArgumentsMIC dpArgs;
    unsigned int * mappedOutputBlock      = mappedRead->outputBlock;
    MICSRAOccMetadata * mappedMetaBlock   = mappedRead->metaBlock;
    
    // Setting up values that is not occurrence specific
    dpArgs.regLength = peArgs->uBound - peArgs->lBound + unmappedRead->seedLength - 1;
    dpArgs.readLength = unmappedRead->seedLength;

    int i, j;
    int startPos, dpStrand;
    for (i=0; i<mappedRead->occCount[0] ; i++){

        // Estimating location from the BASE
        startPos = 0;
        if (mappedMetaBlock[i].strand==(QUERY_POS_STRAND - 1)){
            if ( mappedOutputBlock[i] + peArgs->lBound > unmappedRead->seedLength ) {
                startPos = mappedOutputBlock[i] + peArgs->lBound - unmappedRead->seedLength;
            }
            dpArgs.readCode = unmappedRead->readCode_Complt + unmappedRead->seedOffset_Complt;
            dpStrand = QUERY_NEG_STRAND - 1;
        } else {
            if ( mappedOutputBlock[i] + mappedRead->seedLength > peArgs->uBound ) {
                startPos = mappedOutputBlock[i] + mappedRead->seedLength - peArgs->uBound;
            }
            dpArgs.readCode = unmappedRead->readCode + unmappedRead->seedOffset;
            dpStrand = QUERY_POS_STRAND - 1;
        }
        
        #ifdef MIC_DP_DEBUG_PRINT_ORPHAN_MATCHING
            printf("Base on occurrence %u, estimating %u-%u\n",mappedOutputBlock[i],startPos,dpArgs.regLength);
        #endif
        
        // Extract the reference sequence from HSP with the
        // estimated location from the BASE
        if ( startPos + dpArgs.regLength > hsp->dnaLength ) {
            dpArgs.regLength = hsp->dnaLength - startPos;
        }
        dpArgs.regCode = MICDPExtractCodeFromHSP(
            dpWork->regBuffer,
            startPos,
            dpArgs.regLength,
            hsp
        );
        
        // Invocation of the MIC-DP Module
        // Report DP occurrence of the score is greater than the threshold
        if (DPMatrixFillMIC(&dpArgs, dpWork)){
            #ifdef MIC_DP_DEBUG_PRINT_ORPHAN_MATCHING
                printf("Reporting position %u\n",startPos+dpWork->maxScoreStartPos);
            #endif
            
            // If the number of occurrences found by DP has exceed the defined parameter MIC_DP_OUTPUT_MAX_ALIGNMENT.
            if (dpOccCount >= MIC_DP_OUTPUT_MAX_ALIGNMENT) {
                return MIC_DP_STATUS_TOO_MANY_RESULT;
            }
            if (!_MICDPPopulateOccurrence(dpOutputPtr+dpOccCount, dpWork, startPos, dpStrand)) {
                return MIC_DP_STATUS_ERROR;
            }
            // Write down the PE-SRA occurrence at the position dpOccCount
            peArgs->output[dpOccCount] = mappedOutputBlock[i];
            peArgs->outputMeta[dpOccCount] = mappedMetaBlock[i];
            
            dpOccCount++;
        }

        /*#ifdef MIC_DP_FORCE_INSERTION_SIZE
            if (dpAlgnmtOK) {
                unsigned int dpAmbPos = startPos + dpWork->maxScoreStartPos;
                unsigned int sraAmbPos = mappedRead->outputBlock[i];
                
                unsigned int insertion_1 = sraAmbPos + mappedRead->readLength - dpAmbPos;
                unsigned int insertion_2 = dpAmbPos + dpWork->maxScoreEndPos - dpWork->maxScoreStartPos - sraAmbPos + 1;
                
                // ATTENTION : This checking logic is a short-term solution before
                // the StrandLeftLeg, StrandRightLeg option being implemented.
                if ((insertion_1 <= peArgs->lBound || insertion_1 >= peArgs->uBound) &&
                    (insertion_2 <= peArgs->lBound || insertion_2 >= peArgs->uBound)) {
                    dpAlgnmtOK = FALSE;
                }
            }
        #endif*/
    }
    return dpOccCount;
}


__attribute__((target(mic)))
int MICDPOrphanAlignment(MICDPOccurrence * dpOutputPtr, DPWorkMIC * dpWork, const MICPEArguments * peArgs, 
        const HSP * hsp, double InputSGAOrphanTriggerTF){

    MICSRAArguments * mappedRead, * unmappedRead;
    int dpOccCount = 0;

    ///////////////////////////////////////
    // Default-DP Flow
    ///////////////////////////////////////

    // Perform DP by estimating the location of MATE based on occurrences of READ
    #ifdef MIC_DP_DEBUG_PRINT_ORPHAN_MATCHING
    printf("[ORPHAN] New-Default-DP Flow Initiated. READ base sequence.\n");
    #endif

    if ((peArgs->readArgs->occCount[0] > 0 ) && (peArgs->readArgs->occCount[0] <= (int)(InputSGAOrphanTriggerTF*(double)peArgs->readArgs->seedLength))){
        dpOccCount = _MICDPOrphanRecovery(dpOutputPtr,0,dpWork,peArgs,hsp, peArgs->readArgs,peArgs->mateArgs);
    }

    // check if error occured
    if (dpOccCount<0) return dpOccCount;

    //TODO omg?
    if (*peArgs->outputStatus==MIC_PE_OUTPUT_STATUS_BAD_PAIR && dpOccCount > 0) {
        // Insert a zero'ed cell to separate the output from base READ
        // from  base MATE. The separator only exists of base READ is aligned.
        // if only base MATE is aligned there will not be any separator. 
        // If base READ is aligned then there will be always be a separator.
        if (dpOccCount >= MIC_DP_OUTPUT_MAX_ALIGNMENT) return MIC_DP_STATUS_TOO_MANY_RESULT;

        // Write down the zero'd occurrence at the position dpOccCount
        MICDPOccurrence * dpoToWrite = dpOutputPtr+dpOccCount;
        dpoToWrite -> ambPosition = 0;
        dpoToWrite -> strand = 0;
        dpoToWrite -> matchElemsCount = 0;
        dpoToWrite -> matchLen = 0;
        // Write down the zero'd occurrence at the position dpOccCount
        peArgs->output[dpOccCount] = 0;
        dpOccCount++;
    }

    // Perform DP by estimating the location of READ based on occurrences of MATE
    #ifdef MIC_DP_DEBUG_PRINT_ORPHAN_MATCHING
    printf("[ORPHAN] New-Default-DP Flow Initiated. MATE base sequence.\n");
    #endif

    if ((peArgs->mateArgs->occCount[0] > 0 ) && (peArgs->mateArgs->occCount[0] <= (int)(InputSGAOrphanTriggerTF*(double)peArgs->mateArgs->seedLength))){
        dpOccCount = _MICDPOrphanRecovery(dpOutputPtr,dpOccCount,dpWork,peArgs,hsp, peArgs->mateArgs,peArgs->readArgs);
    }

    return dpOccCount;
}

void MICDPOccurrenceConvert(MICDPOccurrence  * source, DPOccurrence * destination) {
    int i;
    destination->ambPosition                = source->ambPosition;
    destination->strand                     = source->strand + 1;
    destination->matchElemsCount            = source->matchElemsCount;
    destination->mismatchCount              = 0;
    destination->matchLen                   = source->matchLen;
    
    for (i=0;i<source->matchElemsCount&&i<DP_MAX_ERROR_COUNT;i++) {
        destination->matchElems[i].type     = source->matchElems[i].type;
        destination->matchElems[i].length   = source->matchElems[i].length;
        if (source->matchElems[i].type == SRA_CHAR_ERROR_TYPE_MISMATCH) {
            destination->mismatchCount++;
        }
    }
}
