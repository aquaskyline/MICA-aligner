//
//    MIC-PEAlgnmt.c
//
//    SOAP2 / 2BWT
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


#include "MIC-PEAlgnmt.h"



__attribute__((target(mic)))
void MICPEArgumentsConfig(
    MICPEArguments * peArguments,
    MICSRAArguments * readArgs,
    MICSRAArguments * mateArgs,
    unsigned int outputVacancy,
    unsigned int metaVacancy,
    unsigned int * output,
    MICSRAOccMetadata * outputMeta,
    uint16_t * occCount,
    uint8_t * outputStatus) {

    peArguments->readArgs = readArgs;
    peArguments->mateArgs = mateArgs;

    peArguments->output = output;
    peArguments->outputMeta = outputMeta;
    peArguments->outputVacancy = outputVacancy;
    peArguments->metaVacancy = metaVacancy;
    peArguments->outputStatus = outputStatus;
    peArguments->occCount = occCount;

    peArguments->mergeEnable = 0;
}

__attribute__((target(mic)))
static void MICQuickSort(unsigned int * nums, MICSRAOccMetadata * metadata, uint16_t size) {
    if (size <= 1) return;

    unsigned int pivot = nums[size/2];
    unsigned int temp;
    MICSRAOccMetadata tempMeta;
    int left = 0;
    int right = size - 1;
    
    while (left <= right) {
        while (nums[left] < pivot) {
            left++;
        }
        while (nums[right] > pivot) {
            right--;
        }
        if (left <= right) {
            // Swap the elements
            temp = nums[left];
            nums[left] = nums[right];
            nums[right] = temp;
            
            tempMeta = metadata[left];
            metadata[left] = metadata[right];
            metadata[right] = tempMeta;

            left++;
            right--;
        } 
    }

    MICQuickSort(nums, metadata, right + 1);
    MICQuickSort(nums + left, metadata + left, size - left);
}

// Input parameter pos, should be sorted before passing into this routine
__attribute__((target(mic)))
static uint16_t mergePositions(unsigned int * pos, MICSRAOccMetadata * metadata, uint16_t size) {
    if (size == 0) {
        return 0;
    }

    int i;
    uint16_t newSize = 1;
    
    for (i = 1; i < size; i++) {
        if (pos[i] - pos[newSize - 1] < MIC_PAIRING_MERGE_LIMIT) {
            pos[newSize] = pos[i];
            metadata[newSize] = metadata[i];
            newSize++;
        }
    }
    
    return newSize;
}

__attribute__((target(mic)))
static void copySingleEndResults(MICPEArguments * peArgs, MICSRAArguments * sraArgs) {
    if (*sraArgs->occCount > peArgs->outputVacancy ||
            *sraArgs->occCount > peArgs->metaVacancy) {
       (*peArgs->outputStatus) = MIC_PE_OUTPUT_STATUS_UNHANDLE;
       return;
    }

    (*peArgs->occCount) = *sraArgs->occCount;

    // Limit the output
    if (*peArgs->occCount > peArgs->outputLimit &&
            peArgs->outputLimit > -1) {
        *peArgs->occCount = peArgs->outputLimit;
    }

    unsigned int i;
    for (i=0; i<(*peArgs->occCount); i++) {
        peArgs->output[i] = sraArgs->outputBlock[i];
        peArgs->outputMeta[i] = sraArgs->metaBlock[i];
    }
}


__attribute__((target(mic)))
static void pairOneStrand(unsigned int * readOcc, unsigned int * mateOcc,
    MICSRAOccMetadata * readMeta, MICSRAOccMetadata * mateMeta,
    uint32_t readOccSize, uint32_t mateOccSize,
    MICPEArguments * peArguments, char readStrand) {

    // Perform merging
    unsigned int ub = peArguments->uBound;
    unsigned int lb = peArguments->lBound;
    unsigned int outputIdx = *peArguments->occCount;
    uint16_t leftIdx = 0;
    uint16_t rightStart = 0;
    uint16_t rightEnd = 0;
    unsigned int vacancy = peArguments->outputVacancy;
    char shouldEnd = 0;

    uint32_t leftOccSize;
    unsigned int * leftOcc;
    MICSRAOccMetadata * leftMeta;
    uint32_t rightOccSize;
    unsigned int * rightOcc;
    MICSRAOccMetadata * rightMeta;
    uint16_t rightLength;

    // Correctly assign left and right legs
    if (readStrand == 0) {
        // case of read should be leftLeg
        leftOccSize = readOccSize;
        leftOcc = readOcc;
        leftMeta = readMeta;
        rightOccSize = mateOccSize;
        rightOcc = mateOcc;
        rightMeta = mateMeta;
        rightLength = peArguments->mateArgs->seedLength;
    } else {
        // case of read should be rightLeg
        leftOccSize = mateOccSize;
        leftOcc = mateOcc;
        leftMeta = mateMeta;
        rightOccSize = readOccSize;
        rightOcc = readOcc;
        rightMeta = readMeta;
        rightLength = peArguments->readArgs->seedLength;
    }

    // The range is [rightStart, rightEnd)
    // rightEnd is not included.

    // Pair each leftOcc with rightOcc
    // Strands are checked
    for (leftIdx=0; leftIdx<leftOccSize; leftIdx++) {
    
        // This check aims at ignoring any duplicated occurrences on the list
        if (leftIdx > 0 && leftOcc[leftIdx] == leftOcc[leftIdx-1]) {
            continue;
        }
    
        if (leftMeta[leftIdx].strand != 0) {
            continue;
        }

        // Move rightStart into correct index
        while (rightStart < rightOccSize &&
                (leftOcc[leftIdx] + lb > rightOcc[rightStart] + rightLength ||
                rightMeta[rightStart].strand != 1)) {
            rightStart++;
        }
        if (rightStart >= rightOccSize) {
            break;
        }
        // Move rightEnd into correct index
        rightEnd = rightStart;
        while (rightEnd < rightOccSize &&
                (leftOcc[leftIdx] + ub >= rightOcc[rightEnd] + rightLength ||
                rightMeta[rightStart].strand != 1)) {
            rightEnd++;
        }
        //assert(rightEnd <= rightOccSize);
        // Write to output
        int i;
        for (i=rightStart; i<rightEnd &&
                (*peArguments->outputStatus) != MIC_PE_OUTPUT_STATUS_UNHANDLE; i++) {
            if (rightMeta[i].strand != 1) {
                continue;
            }
            if (vacancy - outputIdx >= 2) {
                // Put result into output

                // Check output limit per PE alignment
                if (outputIdx >= MIC_PE_MAX_RESULT) {
                    shouldEnd = 1;
                    (*peArguments->outputStatus) = MIC_PE_OUTPUT_STATUS_CLOSED; 
                    break;
                } else if (outputIdx >= peArguments->outputLimit &&
                        peArguments->outputLimit > -1) {
                    shouldEnd = 1;
                    break;

                }
                //printf("Pair: %d, %d; Gap: %d\n", leftOcc[leftIdx], rightOcc[i], rightOcc[i] - leftOcc[leftIdx]);
                
                // Handle output order
                if (readStrand == 0) {
                    peArguments->outputMeta[outputIdx] = leftMeta[leftIdx];
                    peArguments->output[outputIdx++] = leftOcc[leftIdx];
                }
                peArguments->outputMeta[outputIdx] = rightMeta[i];
                peArguments->output[outputIdx++] = rightOcc[i];

                if (readStrand == 1) {
                    peArguments->outputMeta[outputIdx] = leftMeta[leftIdx];
                    peArguments->output[outputIdx++] = leftOcc[leftIdx];
                }
            } else {
                (*peArguments->outputStatus) = MIC_PE_OUTPUT_STATUS_UNHANDLE;
                shouldEnd = 1;
                break;
            }
        }

        if (shouldEnd) {
            break;
        }
    }

    if ((*peArguments->outputStatus) == MIC_PE_OUTPUT_STATUS_OPEN) {
        (*peArguments->occCount) = outputIdx;
    } else {
        (*peArguments->occCount) = 0;
    }

}


__attribute__((target(mic)))
static void outputLeastMismatchResult(MICPEArguments * peArgs, MICSRAArguments * sraArgs) {
    uint8_t leastError = 100;
    uint16_t selectedOccIdx = 0;
    uint16_t selectedMetaIdx = 0;
    uint16_t occIdx = 0;
    uint16_t metaIdx = 0;
    for (metaIdx=0; metaIdx<(*sraArgs->metaCount); metaIdx++) {
        uint8_t currentErrorNum = sraArgs->metaBlock[metaIdx].numOfErr;
        if (currentErrorNum < leastError) {
            selectedOccIdx = occIdx;
            selectedMetaIdx = metaIdx;
            leastError = currentErrorNum;
        }

        occIdx++;
    }

    // Now selectedOccIdx and selectedMetaIdx should have the least mismatches.
    peArgs->output[(*peArgs->occCount)] = sraArgs->outputBlock[selectedOccIdx];
    peArgs->outputMeta[(*peArgs->occCount)] = sraArgs->metaBlock[selectedMetaIdx];
    
    (*peArgs->occCount)++;
}

__attribute__((target(mic)))
static void pairBadPair(MICPEArguments * peArgs) {
    outputLeastMismatchResult(peArgs, peArgs->readArgs);
    outputLeastMismatchResult(peArgs, peArgs->mateArgs);
    (*peArgs->outputStatus) = MIC_PE_OUTPUT_STATUS_BAD_PAIR;
}


/**
 * Returns number of output occurrences
 **/
__attribute__((target(mic)))
void MICPEMappingInitialise(MICPEArguments * peArguments) {
    (*peArguments->occCount) = 0;
}

/**
 * Returns number of output occurrences
 **/
__attribute__((target(mic)))
void MICPEMappingOccurrences(MICPEArguments * peArguments) {
    MICSRAArguments * readArgs = peArguments->readArgs;
    MICSRAArguments * mateArgs = peArguments->mateArgs;

    if ((*readArgs->outputStatus) == MIC_OUTPUT_STATUS_SKIPPED ||
            (*mateArgs->outputStatus) == MIC_OUTPUT_STATUS_SKIPPED) {
        (*peArguments->outputStatus) = MIC_PE_OUTPUT_STATUS_SKIPPED;
        return;
    }

    if ((*readArgs->outputStatus) == MIC_OUTPUT_STATUS_UNHANDLE ||
            (*mateArgs->outputStatus) == MIC_OUTPUT_STATUS_UNHANDLE ||
            (*readArgs->outputStatus) == MIC_OUTPUT_STATUS_CLOSE ||
            (*mateArgs->outputStatus) == MIC_OUTPUT_STATUS_CLOSE) {
        (*peArguments->outputStatus) = MIC_PE_OUTPUT_STATUS_CLOSED;
        return;
    }

    if ((*readArgs->occCount == 0) || (*mateArgs->occCount == 0)) {
        return;
    }

    if (peArguments->outputVacancy < 2 || peArguments->metaVacancy < 2) {
        (*peArguments->outputStatus) = MIC_PE_OUTPUT_STATUS_UNHANDLE;
        return;
    }

    // Finished handling error cases
    // Remaining possible cases of status
    // PAIR, BAD_PAIR, UNHANDLE, CLOSED

    // Expand the read and mate metaBlock into the length of occCount
    uint32_t readOccSize = *readArgs->occCount;
    uint32_t mateOccSize = *mateArgs->occCount;
    MICSRAOccMetadata * readMeta = readArgs->metaBlock;
    MICSRAOccMetadata * mateMeta = mateArgs->metaBlock;
    // Copy out read and mate occ
    unsigned int * readOcc = readArgs->outputBlock;
    unsigned int * mateOcc = mateArgs->outputBlock;

    // Enhanced quicksort to move extended metaBlock
    // Sort Read Occ
    MICQuickSort(readOcc, readMeta, readOccSize);
    // Sort Mate Occ
    MICQuickSort(mateOcc, mateMeta, mateOccSize);

    if (peArguments->mergeEnable) {
        readOccSize = mergePositions(readOcc, readMeta, readOccSize);
        mateOccSize = mergePositions(mateOcc, mateMeta, mateOccSize);
    }

    pairOneStrand(readOcc, mateOcc, readMeta, mateMeta,
            readOccSize, mateOccSize, peArguments, 0);
    if ((*peArguments->outputStatus) == MIC_PE_OUTPUT_STATUS_OPEN ||
        (*peArguments->outputStatus) == MIC_PE_OUTPUT_STATUS_PAIR) {
        pairOneStrand(readOcc, mateOcc, readMeta, mateMeta,
                readOccSize, mateOccSize, peArguments, 1);
    }
    

    if ((*peArguments->outputStatus) == MIC_PE_OUTPUT_STATUS_OPEN &&
        (*peArguments->occCount) > 0) {
        (*peArguments->outputStatus) = MIC_PE_OUTPUT_STATUS_PAIR;
    }
}


/**
 * Returns number of output occurrences
 **/
__attribute__((target(mic)))
void MICPEMappingComplete(MICPEArguments * peArguments) {
    MICSRAArguments * readArgs = peArguments->readArgs;
    MICSRAArguments * mateArgs = peArguments->mateArgs;

    if ((*peArguments->outputStatus) == MIC_PE_OUTPUT_STATUS_SKIPPED || 
        (*peArguments->outputStatus) == MIC_PE_OUTPUT_STATUS_CLOSED ||
        (*peArguments->outputStatus) == MIC_PE_OUTPUT_STATUS_UNHANDLE) {
        return;
    }

    if ((*readArgs->occCount == 0) && (*mateArgs->occCount == 0)) {
        (*peArguments->outputStatus) = MIC_PE_OUTPUT_STATUS_BOTH_NO_ALIGNMENT;
        return;
    }

    if (*readArgs->occCount == 0) {
        (*peArguments->outputStatus) = MIC_PE_OUTPUT_STATUS_READ_NO_ALIGNMENT;
        copySingleEndResults(peArguments, mateArgs);
        return;
    }

    if (*mateArgs->occCount == 0) {
        (*peArguments->outputStatus) = MIC_PE_OUTPUT_STATUS_MATE_NO_ALIGNMENT;
        copySingleEndResults(peArguments, readArgs);
        return;
    }

    if ((*peArguments->outputStatus) == MIC_PE_OUTPUT_STATUS_OPEN) {
        if ((*peArguments->occCount) == 0) {
            // Proceed to BAD_PAIR
            pairBadPair(peArguments);
        } else {
            (*peArguments->outputStatus) = MIC_PE_OUTPUT_STATUS_PAIR;
        }
    }

}


__attribute__((target(mic)))
void MICPEArgumentsSetBounds(MICPEArguments * peArguments,
    unsigned int lowerBound,
    unsigned int upperBound) {

    peArguments->lBound = lowerBound;
    peArguments->uBound = upperBound;
}

// Setter function for DP related parameters
__attribute__((target(mic)))
void MICPEArgumentsSetDPScores(MICPEArguments * peArguments,
    DPScores * dpScores) {

    peArguments->dpScores = dpScores;
}


// Set a limit to output, occurrences after the limit will be dropped
__attribute__((target(mic)))
void MICPEArgumentsSetMaxOutput(MICPEArguments * peArguments,
        int outputLimit) {

    peArguments->outputLimit = outputLimit;
}



// Set a limit to output, occurrences after the limit will be dropped
__attribute__((target(mic)))
void MICPEArgumentsSetMergeEnable(MICPEArguments * peArguments,
        uint8_t mergeEnable) {

    peArguments->mergeEnable = mergeEnable;
}
