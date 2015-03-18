//
//    DPOccCollector.c
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
   
   Date   : 3rd October 2011
   Author : Edward MK Wu
   Change : New file.

*/
/////////////////////////////////////////////////////

#include "DPOccCollector.h"

inline void DPOCCBucketsInitialise(DPOccurrence_ll * node) {
    if (node==NULL) return;
    DPOCCBucketsInitialise(node->next);
    node->size = 0;
}


inline DPOccurrence_ll * DPOCCBucketCreate() {
    DPOccurrence_ll * newNode = (DPOccurrence_ll*) malloc(sizeof(DPOccurrence_ll));
    newNode->size = 0;
    newNode->next = NULL;
    return newNode;
}

inline DPOccurrence_ll * DPOCCBucketsAddRawOccurrence(DPOccurrence_ll * node, 
                                                        uint8_t type, unsigned long long ambPosition, 
                                                        uint8_t strand, uint16_t matchLen,
                                                        uint16_t mismatchCount, int16_t matchScore,
                                                        uint16_t matchElemsCount, const DPMatchElem matchElems[],
                                                        unsigned int input_regionStartIdx, unsigned int input_regionLength) {
    if (node==NULL) return node;
    
    DPOccurrence_ll * addingIntoNode = node;
    unsigned int size = node->size;
    int i;
    
    //If the current bucket is full.
    if (size>=DPOCC_BUCKET_SIZE) {
        if (node->next==NULL) {
            //If the next bucket has not been created/initialised.
            addingIntoNode = DPOCCBucketCreate();
            node->next = addingIntoNode;
        } else {
            //If the next bucket has already been created/initialised.
            addingIntoNode = node->next;
        }
        size = 0;
    }
    
    addingIntoNode->bucket[size].type             = type;
    addingIntoNode->bucket[size].ambPosition      = ambPosition;
    addingIntoNode->bucket[size].strand           = strand;
    addingIntoNode->bucket[size].mismatchCount    = mismatchCount;
    addingIntoNode->bucket[size].matchLen         = matchLen;
    addingIntoNode->bucket[size].matchElemsCount  = matchElemsCount;
    addingIntoNode->bucket[size].matchScore       = matchScore;
    for (i=0;i<matchElemsCount&&i<DP_MAX_ERROR_COUNT;i++) {
        addingIntoNode->bucket[size].matchElems[i].type     = matchElems[i].type;
        addingIntoNode->bucket[size].matchElems[i].length   = matchElems[i].length;
    }    
    addingIntoNode->bucket[size].input_regionStartIdx = input_regionStartIdx;
    addingIntoNode->bucket[size].input_regionLength = input_regionLength;
    
    addingIntoNode->size++;
    
    return addingIntoNode;
}

inline DPOccurrence_ll * DPOCCBucketsAddOccurrence(DPOccurrence_ll * node, DPOccurrence * occ1) {
    return DPOCCBucketsAddRawOccurrence(
                    node, occ1->type,
                    occ1->ambPosition, 
                    occ1->strand, occ1->matchLen, occ1->mismatchCount,
                    occ1->matchScore, occ1->matchElemsCount, occ1->matchElems,
                    occ1->input_regionStartIdx, occ1->input_regionLength);
}

inline unsigned int DPOCCBucketsCountOccurrences(DPOccurrence_ll * node) {
    if (node==NULL) return 0;
    return node->size + DPOCCBucketsCountOccurrences(node->next);
}


inline unsigned int DPOCCBucketsPopulateDPOccList(DPOccurrence_ll * root, DPOccurrence * sraOcc) {
    unsigned int i;
    unsigned int sraOccIdx;
    DPOccurrence_ll * myNode = root;
    
    sraOccIdx = 0;
    while (myNode != NULL) {
        for (i=0;i<myNode->size;i++) {
            if (myNode->bucket[i].type!=DPOCC_TYPE_DELIMITOR_READ &&
                myNode->bucket[i].type!=DPOCC_TYPE_DELIMITOR_STEP &&
                myNode->bucket[i].type!=DPOCC_TYPE_NO_ALIGNMENT) {
                
                DPOccurrenceCopy(&(myNode->bucket[i]),&(sraOcc[sraOccIdx]));
                sraOccIdx++;
            }
        }
        myNode = myNode->next;
    }
    return sraOccIdx;
}


inline void DPOCCBucketsFree(DPOccurrence_ll * node) {
    if (node==NULL) return;
    DPOCCBucketsFree(node->next);
    free(node);
}


DPOCCCollector * DPOCCCreate(uint8_t floodHandling) {
    int i;
    DPOCCCollector * newOccCollector = (DPOCCCollector*) malloc(sizeof(DPOCCCollector));
    newOccCollector->floodHandling = floodHandling;
    for (i=0;i<256;i++) { newOccCollector->occStats[i] = 0; }
    newOccCollector->occRootNode = DPOCCBucketCreate();
    newOccCollector->occCurrNode = newOccCollector->occRootNode;
    return newOccCollector;
}

void DPOCCInitialise(DPOCCCollector * occCollector) {
    int i;
    DPOCCBucketsInitialise(occCollector->occRootNode);
    for (i=0;i<256;i++) { occCollector->occStats[i] = 0; }
    occCollector->occCurrNode = occCollector->occRootNode;
}

int DPOCCAddOccurrenceToBuckets(DPOCCCollector * occCollector, DPOccurrence * occ1) {
    if (occCollector->floodHandling==DPOCC_FLOOD_TYPE_MAINTAIN_2_SLOT &&
        occCollector->occCurrNode->size >= DPOCC_BUCKET_SIZE-2) {
        return 0;
    } else if (occCollector->floodHandling==DPOCC_FLOOD_TYPE_FLUSH && 
        occCollector->occCurrNode->size == DPOCC_BUCKET_SIZE) {
        return 0;
    }
    occCollector->occStats[occ1->type]++;
    occCollector->occCurrNode = DPOCCBucketsAddOccurrence(occCollector->occCurrNode,occ1);
    return 1;
}

unsigned int DPOCCPopulateDPOccList(DPOCCCollector * occCollector, DPOccurrence * sraOcc) {
    return DPOCCBucketsPopulateDPOccList(occCollector->occRootNode,sraOcc);
}

unsigned int DPOCCCountOccurrences(DPOCCCollector * occCollector) {
    return DPOCCBucketsCountOccurrences(occCollector->occRootNode);
}

unsigned int DPOCCCountOccurrencesByType(DPOCCCollector * occCollector, uint8_t type) {
    return occCollector->occStats[type];
}

void DPOCCFree(DPOCCCollector * occCollector) {
    DPOCCBucketsFree(occCollector->occRootNode);
    free(occCollector);
}





DPOCCCollectorPointer * DPOCCPTCreate(DPOCCCollector * occCollector) {
    DPOCCCollectorPointer * occPtr = (DPOCCCollectorPointer*) malloc(sizeof(DPOCCCollectorPointer));
    DPOCCPTRoot(occCollector,occPtr);
    return occPtr;
}

void DPOCCPTRoot(DPOCCCollector * occCollector, DPOCCCollectorPointer * occPtr) {
    occPtr->bucketPtr = occCollector->occRootNode;
    occPtr->bucketIdx   = 0;
    occPtr->bucketCount = 0;
    if (occPtr->bucketPtr!=NULL) occPtr->bucketCount = occPtr->bucketPtr->size;
}

int DPOCCPTNext(DPOCCCollectorPointer * occPtr) {
    occPtr->bucketIdx++;
    if (occPtr->bucketIdx >= occPtr->bucketCount) {
    
        occPtr->bucketPtr = occPtr->bucketPtr->next;
        occPtr->bucketIdx = 0;
        occPtr->bucketCount = 0;
        if (occPtr->bucketPtr==NULL) return 0;
        if (occPtr->bucketPtr->size==0) return 0;
        
        occPtr->bucketCount = occPtr->bucketPtr->size;
    }
    return 1;
}

DPOccurrence * DPOCCPTRead(DPOCCCollectorPointer * occPtr) {
    if (occPtr->bucketIdx<occPtr->bucketCount) {
        return &(occPtr->bucketPtr->bucket[occPtr->bucketIdx]);
    } else {
        return NULL;
    }
}

void DPOCCPTFree(DPOCCCollectorPointer * occPtr) {
    free(occPtr);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//      DEBUG FUNCTIONS BELOW
//////////////////////////////////////////////////////////////////////////////////////////////////////////

void DPOCCBucketsDebugPrint(DPOccurrence_ll * node) {
    if (node==NULL) return;
    
    unsigned int i;
    char strandChar[5] = ".+-?";
    
    for (i=0;i<node->size;i++) {
        printf("%d\t%c(%d)\t%llu %u\n",node->bucket[i].mismatchCount,strandChar[node->bucket[i].strand],node->bucket[i].strand,node->bucket[i].ambPosition,node->bucket[i].matchLen);
    }
    
    DPOCCBucketsDebugPrint(node->next);
}

void DPOCCDebugPrint(DPOCCCollector * occCollector) {
    if (occCollector==NULL) return;
    DPOCCBucketsDebugPrint(occCollector->occRootNode);
}

unsigned int DPOCCBucketsCountBucket(DPOccurrence_ll * node) {
    if (node==NULL) return 0;
    return 1 + DPOCCBucketsCountBucket(node->next);
}

