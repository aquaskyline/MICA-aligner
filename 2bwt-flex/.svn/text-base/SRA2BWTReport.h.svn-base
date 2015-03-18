//
//    SRA2BWTReport.h
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
   
   Date   : 8th May 2011
   Author : Edward MK Wu
   Change : New file.
   
*/
/////////////////////////////////////////////////////

#ifndef __SRA_REPORT_H__
#define __SRA_REPORT_H__

#include "SRAArguments.h"
#include <time.h>

//----------------------------
// DEBUG FLAGS
//----------------------------
//Define the below parameter to output the alignment result(text position)
// on screen instead of writing into output file
//#define DEBUG_2BWT_OUTPUT_TO_SCREEN

//Define the below parameter to stop cache reported SA range and text position.
// This is turned on when we are trying to obtain alignment timing only.
//#define DEBUG_2BWT_NO_OUTPUT

//Define the below parameter to skip writing the alignment result(text position)
// into disk. This is turned on when we are trying to obtain alignment timing only.
//#define DEBUG_2BWT_NO_WRITE_FILE

//Define the below parameter to skip translating the alignment result(text position).
// This is turned on when we are trying to obtain alignment timing only.
//#define DEBUG_2BWT_NO_TRANSLATION
//----------------------------

unsigned int OCCWriteOutputHeader(SRAOutput * aOutput, SRASetting * aSetting, SRAIndex * aIndex);

//OCCReportSARange : asumming l<=r
unsigned long long OCCReportSARange(SRAArguments * aArgs,
                                    unsigned long long l, unsigned long long r, 
                                    int occMismatch, int occQuality);

void OCCReportDelimitor(SRAArguments * aArgs);

void OCCReportNoAlignment(SRAArguments * aArgs);

unsigned long long OCCReportTextPositions(SRAArguments * aArgs, int posCount, SRAError errors[][MAX_NUM_OF_ERROR]);

void OCCFlushCache(SRAArguments * aArgs);

void OCCFlushBestIntoCache(SRAArguments * aArgs);

// OCCTranslateOccurrence translates a position w.r.t. to the ambiguous-processed position to a
// chromosomeId,Offset pair. The Translate structure however only supports 32-bit operations;
// therefore unsigned int correctPosition is needed to correct the pos/neg values of 32-bit
// unsigned integer.
void OCCTranslateOccurrence(Translate * occTranslate, unsigned short * occAmbiguityMap,
                            unsigned long long ambPos,
                            unsigned char * ChromId, unsigned long long * unambPos);

// OCCAddTextPositionToCache is called by OCCReport* functions to add 1 record to the position cache.
void OCCAddTextPositionToCache(SRAArguments * aArgs, 
                                    uint8_t type, 
                                    unsigned long long ambPosition, 
                                    uint8_t strand, uint16_t matchLen,
                                    uint8_t occMismatch, int occQuality,
                                    const SRAError errorPos[]);











#endif
