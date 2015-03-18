//
//    SRAOccCollector.h
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


#ifndef __SRA_OCC_COLLECTOR_H__
#define __SRA_OCC_COLLECTOR_H__

#include "SRACore.h"

#define SRAOCC_BUCKET_SIZE      8192

#define SRAOCC_FLOOD_TYPE_FLUSH                0
#define SRAOCC_FLOOD_TYPE_MAINTAIN_2_SLOT      1
#define SRAOCC_FLOOD_TYPE_EXPAND               2

#define SRAOCC_TYPE_DELIMITOR_READ                  0
#define SRAOCC_TYPE_DELIMITOR_STEP                  1
#define SRAOCC_TYPE_NO_ALIGNMENT                    5
#define SRAOCC_TYPE_AWAIT_TRANSLATE                 6
#define SRAOCC_TYPE_DP_SEED_OCC                     7
#define SRAOCC_TYPE_MAPQ_VALUE                      8

#define SRAOCC_TYPE_PE_PAIR                         10
#define SRAOCC_TYPE_PE_BAD_PAIR                     11
#define SRAOCC_TYPE_PE_READ_READ_NO_ALIGNMENT       12
#define SRAOCC_TYPE_PE_READ_MATE_NO_ALIGNMENT       13
#define SRAOCC_TYPE_PE_READ_BOTH_NO_ALIGNMENT       14
#define SRAOCC_TYPE_PE_MATE_READ_NO_ALIGNMENT       15
#define SRAOCC_TYPE_PE_MATE_MATE_NO_ALIGNMENT       16
#define SRAOCC_TYPE_PE_MATE_BOTH_NO_ALIGNMENT       17

#define SRAOCC_TYPE_PE_DP_BASE_READ                 20
#define SRAOCC_TYPE_PE_DP_BASE_MATE                 21
#define SRAOCC_TYPE_PE_DP_SEED_OCC                  22
#define SRAOCC_TYPE_PE_DP_SINGLE_END_OCC            23


typedef struct SRAOccurrence_ll {
    SRAOccurrence bucket[SRAOCC_BUCKET_SIZE];
    unsigned int size;
    struct SRAOccurrence_ll * next;
} SRAOccurrence_ll;


typedef struct SRAOCCCollector {
    uint8_t floodHandling;
    unsigned int occStats[256];
    
    SRAOccurrence_ll * occRootNode;
    SRAOccurrence_ll * occCurrNode;
} SRAOCCCollector;

typedef struct SRAOCCCollectorPointer {

    SRAOccurrence_ll * bucketPtr;
    
    //The current and the count of the bucket being
    //pointed at by bucketPtr.
    unsigned int bucketIdx;
    unsigned int bucketCount;
} SRAOCCCollectorPointer;



//////////////////////////////////////////////////////////////////////////////////////////////////////////
//      SRA OCC COLLECTOR FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// SRAOCCCreate creates and initialise a occ collector.
SRAOCCCollector * SRAOCCCreate(uint8_t floodHandling);

// SRAOCCInitialise reinitialises the occ collector so it
// could be re-used in new alignment.
void SRAOCCInitialise(SRAOCCCollector * root);

// SRAOCCAddOccurrenceToBuckets tries adding one occurrence into the collector
// It returns 1 if the addition was successful. 0 otherwise;
int SRAOCCAddOccurrenceToBuckets(SRAOCCCollector * node, SRAOccurrence * occ1);

// SRAOCCCountOccurrences gives a count of the total number of occurrences being stored 
// inside the collector.
unsigned int SRAOCCCountOccurrences(SRAOCCCollector * occCollector);

// SRAOCCCountOccurrences gives a count of the total number of occurrences being stored 
// inside the collector with specific type.
unsigned int SRAOCCCountOccurrencesByType(SRAOCCCollector * occCollector, uint8_t type);

// SRAOCCPopulateSRAOccList populates the occurrences inside the collector
// into a SRAOccurrence list. The memory required should be initialised and
// handled by the caller function.
unsigned int SRAOCCPopulateSRAOccList(SRAOCCCollector * occCollector, SRAOccurrence * sraOcc);

// SRAOCCFree frees all the memory used by the collector.
void SRAOCCFree(SRAOCCCollector * root);







//////////////////////////////////////////////////////////////////////////////////////////////////////////
//      SRA OCC COLLECTOR POINTER FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////

//SRAOCCPTCreate construct an instance of the SRAOCCCollectorPointer
SRAOCCCollectorPointer * SRAOCCPTCreate(SRAOCCCollector * occCollector);

// SRAOCCPTRoot move the pointer to the first record in the occCollector
void SRAOCCPTRoot(SRAOCCCollector * occCollector, SRAOCCCollectorPointer * occPtr);

// SRAOCCPTNext move the pointer to the next record in the occCollector
int SRAOCCPTNext(SRAOCCCollectorPointer * occPtr);

// SRAOCCPTRead returns the pointer to the record in SRAOCCCollector
SRAOccurrence * SRAOCCPTRead(SRAOCCCollectorPointer * occPtr);

// SRAOCCPTFree frees all the memory used by the pointer
void SRAOCCPTFree(SRAOCCCollectorPointer * occPtr);






//////////////////////////////////////////////////////////////////////////////////////////////////////////
//      DEBUG FUNCTIONS BELOW
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// SRAOCCBucketsDebugPrint prints all the items in SRAOCCBuckets.
// This should only be used for debug purpose.
void SRAOCCBucketsDebugPrint(SRAOccurrence_ll * root);

// SRAOCCDebugPrint prints all the items in SRAOCCCollector.
// This should only be used for debug purpose.
void SRAOCCDebugPrint(SRAOCCCollector * occCollector);

// SRAOCCBucketsCountBucket counts the number of buckets
// used by the ll.
unsigned int SRAOCCBucketsCountBucket(SRAOccurrence_ll * root);


#endif
