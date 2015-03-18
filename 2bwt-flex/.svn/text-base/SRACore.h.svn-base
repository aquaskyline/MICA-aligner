//
//    SRACore.h
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
   
   Date   : 30th October 2011
   Author : Edward MK Wu
   Change : New file.

*/
/////////////////////////////////////////////////////

#ifndef __SRA_CORE_H__
#define __SRA_CORE_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "../Release.h"
#include "2bwt-lib/BWT.h"
#include "2bwt-lib/HSP.h"
#include "HOCC.h"
#include "LT.h"

#define SRA_MIN_READ_LENGTH                         25
#define SRA_MAX_READ_LENGTH                         200
#define MAX_NUM_OF_MISMATCH                         5
#define MAX_NUM_OF_INDEL                            5
#define MAX_NUM_OF_NBM_ERROR                        10

#define SRA_TYPE_MISMATCH_ONLY                      0
#define SRA_TYPE_EDIT_DISTANCE                      1

__attribute__((target(mic)))
static const char SRA_TYPE[][128] = {
    "SRA_TYPE_MISMATCH_ONLY",
    "SRA_TYPE_EDIT_DISTANCE" };
    
// MAX_NUM_OF_ERROR must be greater than both
// MAX(MAX_NUM_OF_MISMATCH and MAX_NUM_OF_INDEL) + MAX_NUM_OF_NBM_ERROR
#define MAX_NUM_OF_ERROR                            15

#define SRA_CHAR_ERROR_TYPE_NONE                    0
#define SRA_CHAR_ERROR_TYPE_MISMATCH                1
#define SRA_CHAR_ERROR_TYPE_INSERT                  2
#define SRA_CHAR_ERROR_TYPE_DELETE                  3
#define SRA_CHAR_ERROR_TYPE_NBMISMATCH              4
#define SRA_CHAR_ERROR_TYPE_SOFTCLIP                5

// SRA_REPORT_NONE - Report nothing to the output buffer
//                   so only counting how many occurrences there are.
#define SRA_REPORT_UNIQUE_BEST                      3
#define SRA_REPORT_RANDOM_BEST                      4
#define SRA_REPORT_ALL                              1
#define SRA_REPORT_ALL_BEST                         2
#define SRA_REPORT_BEST_QUALITY                     5
#define SRA_REPORT_DETAIL                           6
#define SRA_REPORT_EXACT_NUM_ERROR                  7
#define SRA_REPORT_NONE                             8
#define SRA_REPORT_ALL_SORTED                       9

//#define SRA_OUTPUT_FORMAT_DEFAULT                   0
// ATTENTION: Constant is obsolete and should be removed.
#define SRA_OUTPUT_FORMAT_PLAIN                     1
#define SRA_OUTPUT_FORMAT_SAM                       2
#define SRA_OUTPUT_FORMAT_BAM                       3
#define SRA_OUTPUT_FORMAT_UNDEF                     10
#define SRA_OUTPUT_FORMAT_NONE                      11
#define SRA_OUTPUT_FORMAT_DISCARD                   12
#define SRA_OUTPUT_FORMAT_SAM_STORE                 13
    
typedef struct SRAError {
    uint16_t type:4, position:12;
} SRAError;

typedef struct SRAOccurrence {
    uint8_t type;
    unsigned long long ambPosition;
    uint8_t strand;
    uint8_t mismatchCount;
    SRAError errors[MAX_NUM_OF_ERROR];
    uint16_t matchLen;
} SRAOccurrence;

typedef struct SRAQueryInfo {
    unsigned long long ReadId;
    char * ReadName;
    unsigned long long ReadLength;
    unsigned char * ReportingReadCode; //This is identical to ReadCode in QueryInputPos; But different in QueryInputNeg
    unsigned char * ReadCode;
    unsigned char ReadStrand;
    unsigned char * ReadQuality;
} SRAQueryInfo;

typedef struct SRASetting {
    uint8_t ReadStrand;
    uint8_t MaxNBMismatch;
    uint8_t MaxError;
    uint8_t ErrorType;
    uint8_t OutputType;
    uint8_t OutputFormat;
    int16_t MaxResult;
    double MaxOptionalRegionPct;
} SRASetting;

typedef struct SRAIndex {
    BWT * bwt;
    BWT * rev_bwt;
    HSP * hsp;
    HOCC * highOcc;
    LT * lookupTable;
    LT * rev_lookupTable;
} SRAIndex;

/////////////////////////////////////////////////////
/*  Function SRAOccurrenceCopy --
    Copy the SRAOccurrence fully into destination
*////////////////////////////////////////////////////
void SRAOccurrenceCopy(SRAOccurrence * source, SRAOccurrence * destination);


/////////////////////////////////////////////////////
/*  Function SRARadixSort --
    Sort the list of SRAOccurrences by Radix Sort
    assuming the buffer1 and buffer2 has sufficient
    length and allocated memory.
*////////////////////////////////////////////////////
void SRARadixSort(SRAOccurrence * list, unsigned long long * _count,
                 SRAOccurrence * buffer1, SRAOccurrence * buffer2,
                 SRAOccurrence * sortedList);
                 
/////////////////////////////////////////////////////
/*  Function SRAOccurrencesSort --
    Sort the list of SRAOccurrences by Radix Sort.
    A wrapper for above so user doesn't need to
    worry about memory allocation.
*////////////////////////////////////////////////////
void SRAOccurrencesSort(SRAOccurrence * occList_1, unsigned long long * _occCount_1);


/////////////////////////////////////////////////////
/*  Function SRAOccurrencesPrint --
    Debug function to output a list of SRAOccurrences
    to standard output.
*////////////////////////////////////////////////////
void SRAOccurrencesPrint(SRAOccurrence * occList_1, unsigned long long occCount_1);


/////////////////////////////////////////////////////
/*  Function SRAOccurrencesPrint --
    Construction and Destruction
*////////////////////////////////////////////////////
SRAQueryInfo * SRAQueryInfoConstruct();
void SRAQueryInfoPopulate(SRAQueryInfo * qInfo, 
                            unsigned long long ReadId, char * ReadName,
                            unsigned long long ReadLength,
                            unsigned char * ReadCode, unsigned char * ReportingReadCode,
                            unsigned char ReadStrand,
                            char * ReadQuality);
SRAQueryInfo * SRAQueryInfoMakeClone(SRAQueryInfo * qInfo);
SRAQueryInfo * SRAQueryInfoCopy(SRAQueryInfo * source, SRAQueryInfo * destination);
void SRAQueryInfoFree(SRAQueryInfo * qInfo);
SRAIndex * SRAIndexConstruct();
void SRAIndexPopulate(SRAIndex * aIndex, BWT * bwt, BWT * rev_bwt, HSP * hsp, 
                        HOCC * highOcc, LT * lookupTable, LT * rev_lookupTable);
void SRAIndexFree(SRAIndex * aIndex);
SRASetting * SRASettingConstruct();
SRASetting * SRASettingMakeClone(SRASetting * source);
void SRASettingFree(SRASetting * aSettings);

void SRAFlipPattern(unsigned char * charMap, unsigned char * ReadCode, unsigned long long ReadLength, unsigned char * FlipReadCode);
void SRAFlipQuality(unsigned char * charMap,unsigned char * ReadQuality, unsigned long long ReadLength, unsigned char * FlipReadQuality);
void DEBUGSRAQueryInfoPrint(SRAQueryInfo * qInfo, char * tag);
void DEBUGSRASequencePrint(const char * tag, unsigned char * convertedPattern, int length);

#endif
