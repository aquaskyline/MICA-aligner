//
//    DPOccCollector.h
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


#ifndef __DP_OCC_COLLECTOR_H__
#define __DP_OCC_COLLECTOR_H__

#include "DPCore.h"

#define DPOCC_BUCKET_SIZE      8192

#define DPOCC_FLOOD_TYPE_FLUSH                0
#define DPOCC_FLOOD_TYPE_MAINTAIN_2_SLOT      1
#define DPOCC_FLOOD_TYPE_EXPAND               2

#define DPOCC_TYPE_DELIMITOR_READ                  0
#define DPOCC_TYPE_DELIMITOR_STEP                  1
#define DPOCC_TYPE_NO_ALIGNMENT                    5
#define DPOCC_TYPE_AWAIT_TRANSLATE                 6

#define DPOCC_TYPE_PE_PAIR                         10
#define DPOCC_TYPE_PE_BAD_PAIR                     11
#define DPOCC_TYPE_PE_READ_READ_NO_ALIGNMENT       12
#define DPOCC_TYPE_PE_READ_MATE_NO_ALIGNMENT       13
#define DPOCC_TYPE_PE_READ_BOTH_NO_ALIGNMENT       14
#define DPOCC_TYPE_PE_MATE_READ_NO_ALIGNMENT       15
#define DPOCC_TYPE_PE_MATE_MATE_NO_ALIGNMENT       16
#define DPOCC_TYPE_PE_MATE_BOTH_NO_ALIGNMENT       17

#define DPOCC_TYPE_PE_DP_BASE_READ                 20
#define DPOCC_TYPE_PE_DP_BASE_MATE                 21
#define DPOCC_TYPE_PE_DP_SEED_OCC                  22

typedef struct DPOccurrence_ll {
    DPOccurrence bucket[DPOCC_BUCKET_SIZE];
    unsigned int size;
    struct DPOccurrence_ll * next;
} DPOccurrence_ll;


typedef struct DPOCCCollector {
    uint8_t floodHandling;
    unsigned int occStats[256];
    
    DPOccurrence_ll * occRootNode;
    DPOccurrence_ll * occCurrNode;
} DPOCCCollector;

typedef struct DPOCCCollectorPointer {

    DPOccurrence_ll * bucketPtr;
    
    //The current and the count of the bucket being
    //pointed at by bucketPtr.
    unsigned int bucketIdx;
    unsigned int bucketCount;
} DPOCCCollectorPointer;



//////////////////////////////////////////////////////////////////////////////////////////////////////////
//      SRA OCC COLLECTOR FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// DPOCCCreate creates and initialise a occ collector.
DPOCCCollector * DPOCCCreate(uint8_t floodHandling);

// DPOCCInitialise reinitialises the occ collector so it
// could be re-used in new alignment.
void DPOCCInitialise(DPOCCCollector * root);

// DPOCCAddOccurrenceToBuckets tries adding one occurrence into the collector
// It returns 1 if the addition was successful. 0 otherwise;
int DPOCCAddOccurrenceToBuckets(DPOCCCollector * node, DPOccurrence * occ1);

// DPOCCCountOccurrences gives a count of the total number of occurrences being stored 
// inside the collector.
unsigned int DPOCCCountOccurrences(DPOCCCollector * occCollector);

// DPOCCCountOccurrences gives a count of the total number of occurrences being stored 
// inside the collector with specific type.
unsigned int DPOCCCountOccurrencesByType(DPOCCCollector * occCollector, uint8_t type);

// DPOCCPopulateDPOccList populates the occurrences inside the collector
// into a DPOccurrence list. The memory required should be initialised and
// handled by the caller function.
unsigned int DPOCCPopulateDPOccList(DPOCCCollector * occCollector, DPOccurrence * sraOcc);

// DPOCCFree frees all the memory used by the collector.
void DPOCCFree(DPOCCCollector * root);







//////////////////////////////////////////////////////////////////////////////////////////////////////////
//      SRA OCC COLLECTOR POINTER FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////

//DPOCCPTCreate construct an instance of the DPOCCCollectorPointer
DPOCCCollectorPointer * DPOCCPTCreate(DPOCCCollector * occCollector);

// DPOCCPTRoot move the pointer to the first record in the occCollector
void DPOCCPTRoot(DPOCCCollector * occCollector, DPOCCCollectorPointer * occPtr);

// DPOCCPTNext move the pointer to the next record in the occCollector
int DPOCCPTNext(DPOCCCollectorPointer * occPtr);

// DPOCCPTRead returns the pointer to the record in DPOCCCollector
DPOccurrence * DPOCCPTRead(DPOCCCollectorPointer * occPtr);

// DPOCCPTFree frees all the memory used by the pointer
void DPOCCPTFree(DPOCCCollectorPointer * occPtr);






//////////////////////////////////////////////////////////////////////////////////////////////////////////
//      DEBUG FUNCTIONS BELOW
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// DPOCCBucketsDebugPrint prints all the items in DPOCCBuckets.
// This should only be used for debug purpose.
void DPOCCBucketsDebugPrint(DPOccurrence_ll * root);

// DPOCCDebugPrint prints all the items in DPOCCCollector.
// This should only be used for debug purpose.
void DPOCCDebugPrint(DPOCCCollector * occCollector);

// DPOCCBucketsCountBucket counts the number of buckets
// used by the ll.
unsigned int DPOCCBucketsCountBucket(DPOccurrence_ll * root);


#endif