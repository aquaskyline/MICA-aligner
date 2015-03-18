//
//    OutputSAMFormat.h
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
   
   Date   : 25th July 2011
   Author : Edward MK Wu
   Change : New file.
   
*/
/////////////////////////////////////////////////////

#ifndef __OUTPUT_SAM_FORMAT_H__
#define __OUTPUT_SAM_FORMAT_H__

#include "SRACore.h"
#include "DPCore.h"

#include "samtools-0.1.19/sam.h"
#include "ListConsumer.h"

#define SAM_FLAG_IS_PAIR_READ           1
#define SAM_FLAG_PROPER_PAIR            2
#define SAM_FLAG_READ_UNMAPPED          4
#define SAM_FLAG_MATE_UNMAPPED          8

#define SAM_FLAG_READ_ALGNMT_STRAND     16
#define SAM_FLAG_MATE_ALGNMT_STRAND     32
#define SAM_FLAG_FIRST_IN_PAIR          64
#define SAM_FLAG_SECOND_IN_PAIR         128

#define SAM_FLAG_SECONDARY_ALGNMT       256
#define SAM_FLAG_QC_CHECK_FAILED        512
#define SAM_FLAG_DUPLICATE              1024

#define SAM_MDATA_INITIAL_ALLOWED_OCC   1024
#define SAM_MDATA_INITIAL_ALLOWED_DPOCC 100
#define SAM_MDATA_SIZE_PER_READ         2048
#define SAM_MDATA_SIZE_PER_OCC          (MAX_SEQ_NAME_LENGTH + 4*(MAX_NUM_OF_MISMATCH*2+2)) //Size in byte
#define SAM_MDATA_SIZE_PER_DPOCC        (MAX_SEQ_NAME_LENGTH + 4*(DP_MAX_ERROR_COUNT)) //Size in byte

#define SAM_ERROR_IN_ASC                 0
#define SAM_ERROR_IN_DESC                1
#define SAM_ERROR_RANDOM                 2

#define SAM_MAX_CHROMOSOME_NAME_LENGTH  256

// SAM_MAPQ_UNAVAILABLE should be set to 255 
// according to SAM v1.4 specification.
// However we are making an exception here
// as some of the downstream programs are not
// capable of handling 255.
#define SAM_MAPQ_UNAVAILABLE            255
#define SAM_MAPQ_UNALIGNED              0

int SAMSeqReadNameLength(char * str, int maxLen);

void SAMIUint8ConcatUint8(uint8_t * data, int * curSize,
                                uint8_t key);
void SAMIUint8ConcatUint32(uint8_t * data, int * curSize,
                                uint32_t key);

__attribute__((target(mic)))
int SAMIUint8ConcatUint32AsString(uint8_t * data, int * curSize,
                                uint32_t key);

__attribute__((target(mic)))
void SAMIUint8ConcatString(uint8_t * data, int * curSize,
                                char * key, int len);
                                
int SAMIUint8ConcatCigar(uint8_t * data, int * curSize, 
                                unsigned int readLength, uint8_t readStrand, 
                                uint8_t errorCount, SRAError errors[], uint8_t errorOrder);
int SAMIUint8ConcatCigarAsString(uint8_t * data, int * curSize, 
                        unsigned int readLength, uint8_t readStrand, 
                        uint8_t errorCount, SRAError errors[], uint8_t errorOrder);
void SAMPackSequenceAndCigar(bam1_t * samAlgnmt, int sharedDataLen, SRAQueryInfo * qInfo, SRAOccurrence * occ);

//Dynamic Programming
void SAMDPPackSequenceAndCigar(bam1_t * samAlgnmt, int sharedDataLen, SRAQueryInfo * qInfo, DPOccurrence * dpOcc);

__attribute__((target(mic)))
int SAMDPIUint8ConcatCigarAsString(uint8_t * data, int * curSize, 
                        unsigned int readLength, uint8_t readStrand, uint16_t matchElemsCount, DPMatchElem matchElems[]);
                        
void SAMOutputHeaderConstruct(bam_header_t * sheader, HSP * hsp);
void SAMOutputHeaderConstruct_ListReader(bam_header_t * sheader, HSP * hsp, ListReader * listreader);
void SAMOutputHeaderDestruct(bam_header_t * sheader);

#endif
