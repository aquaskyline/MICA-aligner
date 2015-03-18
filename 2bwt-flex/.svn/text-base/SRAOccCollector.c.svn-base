//
//    SRAOccCollector.c
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

#include "SRAOccCollector.h"

inline void SRAOCCBucketsInitialise(SRAOccurrence_ll * node) {
    if (node==NULL) return;
    SRAOCCBucketsInitialise(node->next);
    node->size = 0;
}


inline SRAOccurrence_ll * SRAOCCBucketCreate() {
    SRAOccurrence_ll * newNode = (SRAOccurrence_ll*) malloc(sizeof(SRAOccurrence_ll));
    newNode->size = 0;
    newNode->next = NULL;
    return newNode;
}

inline SRAOccurrence_ll * SRAOCCBucketsAddRawOccurrence(SRAOccurrence_ll * node, uint8_t type, unsigned long long ambPosition, uint8_t strand, uint16_t matchLen, uint8_t mismatchCount, const SRAError errorPos[]) {
    if (node==NULL) return node;
    
    SRAOccurrence_ll * addingIntoNode = node;
    unsigned int size = node->size;
    int i;
    
    //If the current bucket is full.
    if (size>=SRAOCC_BUCKET_SIZE) {
        if (node->next==NULL) {
            //If the next bucket has not been created/initialised.
            addingIntoNode = SRAOCCBucketCreate();
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
    for (i=0;i<mismatchCount&&i<MAX_NUM_OF_ERROR;i++) {
        addingIntoNode->bucket[size].errors[i].type = errorPos[i].type;
        addingIntoNode->bucket[size].errors[i].position = errorPos[i].position;
    }    
    addingIntoNode->size++;
    
    return addingIntoNode;
}

inline SRAOccurrence_ll * SRAOCCBucketsAddOccurrence(SRAOccurrence_ll * node, SRAOccurrence * occ1) {
    return SRAOCCBucketsAddRawOccurrence(node,occ1->type,occ1->ambPosition,occ1->strand,occ1->matchLen,occ1->mismatchCount,occ1->errors);
}

inline unsigned int SRAOCCBucketsCountOccurrences(SRAOccurrence_ll * node) {
    if (node==NULL) return 0;
    return node->size + SRAOCCBucketsCountOccurrences(node->next);
}


inline unsigned int SRAOCCBucketsPopulateSRAOccList(SRAOccurrence_ll * root, SRAOccurrence * sraOcc) {
    unsigned int i;
    unsigned int sraOccIdx;
    SRAOccurrence_ll * myNode = root;
    
    sraOccIdx = 0;
    while (myNode != NULL) {
        for (i=0;i<myNode->size;i++) {
            if (myNode->bucket[i].type!=SRAOCC_TYPE_DELIMITOR_READ &&
                myNode->bucket[i].type!=SRAOCC_TYPE_DELIMITOR_STEP &&
                myNode->bucket[i].type!=SRAOCC_TYPE_NO_ALIGNMENT) {
                
                SRAOccurrenceCopy(&(myNode->bucket[i]),&(sraOcc[sraOccIdx]));
                sraOccIdx++;
            }
        }
        myNode = myNode->next;
    }
    return sraOccIdx;
}


inline void SRAOCCBucketsFree(SRAOccurrence_ll * node) {
    if (node==NULL) return;
    SRAOCCBucketsFree(node->next);
    free(node);
}


SRAOCCCollector * SRAOCCCreate(uint8_t floodHandling) {
    int i;
    SRAOCCCollector * newOccCollector = (SRAOCCCollector*) malloc(sizeof(SRAOCCCollector));
    newOccCollector->floodHandling = floodHandling;
    for (i=0;i<256;i++) { newOccCollector->occStats[i] = 0; }
    newOccCollector->occRootNode = SRAOCCBucketCreate();
    newOccCollector->occCurrNode = newOccCollector->occRootNode;
    return newOccCollector;
}

void SRAOCCInitialise(SRAOCCCollector * occCollector) {
    int i;
    SRAOCCBucketsInitialise(occCollector->occRootNode);
    for (i=0;i<256;i++) { occCollector->occStats[i] = 0; }
    occCollector->occCurrNode = occCollector->occRootNode;
}

int SRAOCCAddOccurrenceToBuckets(SRAOCCCollector * occCollector, SRAOccurrence * occ1) {
    if (occCollector->floodHandling==SRAOCC_FLOOD_TYPE_MAINTAIN_2_SLOT &&
        occCollector->occCurrNode->size >= SRAOCC_BUCKET_SIZE-2) {
        return 0;
    } else if (occCollector->floodHandling==SRAOCC_FLOOD_TYPE_FLUSH && 
        occCollector->occCurrNode->size == SRAOCC_BUCKET_SIZE) {
        return 0;
    }
    occCollector->occStats[occ1->type]++;
    occCollector->occCurrNode = SRAOCCBucketsAddOccurrence(occCollector->occCurrNode,occ1);
    return 1;
}

unsigned int SRAOCCPopulateSRAOccList(SRAOCCCollector * occCollector, SRAOccurrence * sraOcc) {
    return SRAOCCBucketsPopulateSRAOccList(occCollector->occRootNode,sraOcc);
}

unsigned int SRAOCCCountOccurrences(SRAOCCCollector * occCollector) {
    return SRAOCCBucketsCountOccurrences(occCollector->occRootNode);
}

unsigned int SRAOCCCountOccurrencesByType(SRAOCCCollector * occCollector, uint8_t type) {
    return occCollector->occStats[type];
}

void SRAOCCFree(SRAOCCCollector * occCollector) {
    SRAOCCBucketsFree(occCollector->occRootNode);
    free(occCollector);
}





SRAOCCCollectorPointer * SRAOCCPTCreate(SRAOCCCollector * occCollector) {
    SRAOCCCollectorPointer * occPtr = (SRAOCCCollectorPointer*) malloc(sizeof(SRAOCCCollectorPointer));
    SRAOCCPTRoot(occCollector,occPtr);
    return occPtr;
}

void SRAOCCPTRoot(SRAOCCCollector * occCollector, SRAOCCCollectorPointer * occPtr) {
    occPtr->bucketPtr = occCollector->occRootNode;
    occPtr->bucketIdx   = 0;
    occPtr->bucketCount = 0;
    if (occPtr->bucketPtr!=NULL) occPtr->bucketCount = occPtr->bucketPtr->size;
}

int SRAOCCPTNext(SRAOCCCollectorPointer * occPtr) {
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

SRAOccurrence * SRAOCCPTRead(SRAOCCCollectorPointer * occPtr) {
    if (occPtr->bucketIdx<occPtr->bucketCount) {
        return &(occPtr->bucketPtr->bucket[occPtr->bucketIdx]);
    } else {
        return NULL;
    }
}

void SRAOCCPTFree(SRAOCCCollectorPointer * occPtr) {
    free(occPtr);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//      DEBUG FUNCTIONS BELOW
//////////////////////////////////////////////////////////////////////////////////////////////////////////

void SRAOCCBucketsDebugPrint(SRAOccurrence_ll * node) {
    if (node==NULL) return;
    
    unsigned int i;
    char strandChar[5] = ".+-?";
    
    for (i=0;i<node->size;i++) {
        printf("%d\t%d\t%c(%d)\t%llu %u\n",node->bucket[i].type,node->bucket[i].mismatchCount,strandChar[node->bucket[i].strand],node->bucket[i].strand,node->bucket[i].ambPosition,node->bucket[i].matchLen);
    }
    
    SRAOCCBucketsDebugPrint(node->next);
}

void SRAOCCDebugPrint(SRAOCCCollector * occCollector) {
    if (occCollector==NULL) return;
    printf("[SRAOCCDebugPrint] type mCnt str() ambPos mLen\n");
    SRAOCCBucketsDebugPrint(occCollector->occRootNode);
    printf("[SRAOCCDebugPrint] Finished\n");
}

unsigned int SRAOCCBucketsCountBucket(SRAOccurrence_ll * node) {
    if (node==NULL) return 0;
    return 1 + SRAOCCBucketsCountBucket(node->next);
}

