//
//    PE2BWTReport.h
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
   
   Date   : 19th October 2011
   Author : Edward MK Wu
   Change : New file created based on 2BWT-AlgnmtRpt.h
   
*/
/////////////////////////////////////////////////////

#ifndef __PE_REPORT_H__
#define __PE_REPORT_H__

#include "PEArguments.h"
#include "SRA2BWTReport.h"
#include "PEDPReport.h"

#define PE_OCC_OUTPUT_FORMAT            20111019

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

unsigned int PEOCCWriteOutputHeader(PEArguments * peArgs);

void PEOCCFlushCache(PEArguments * peArgs);

//PEOCCCountAllPEAlignment counts how many PE alignments are there
//stored in the occCollector, including SRA PE and DP PE results.
unsigned int PEOCCCountAllPEAlignment(PEArguments * peArgs);

//PEOCCDumpAlignments dumps alignment from PEOuput to OccCollector.
unsigned long long PEOCCDumpAlignments(PEArguments * peArgs);
//PEOCCDumpAlignments dumps one alignment to OccCollector.
unsigned long long PEOCCDumpOneAlignment(PEArguments * peArgs, SRAOccurrence * occ_1, SRAOccurrence * occ_2);

//PEOCCDumpSingleAlignments dumps alignment from SRAOccurrence list to OccCollector.
unsigned long long PEOCCDumpSingleAlignments(PEArguments * peArgs, SRAOccurrence * occ, unsigned int count, uint8_t type);

//PEOCCReportSubValidAlignment for reporting a invalid pair-end alignment for random-best alignment
// if a valid pair-end alignment is not available then a sub-valid pair-end alignment is reported.
void PEOCCReportSubValidAlignment(PEArguments * peArgs, SRAOccurrence * occ_1, SRAOccurrence * occ_2);

//PEOCCReportNoAlignment for reporting short read alignment results without a valid pair-end alignment found
// when one or both of the reads does not produce any short read alignment.
void PEOCCReportNoAlignment(PEArguments * peArgs,uint8_t type);

//PEOCCReportDelimitor helps plain/default format to separate results.
void PEOCCReportDelimitor(PEArguments * peArgs);

void PEOCCAddTextPositionToCache(PEArguments * peArgs, 
                                    uint8_t type,
                                    unsigned long long ambPosition, 
                                    uint8_t strand, uint16_t matchLen, uint8_t occMismatch, int occQuality,
                                    const SRAError errorPos[]);

void MAPQPEGetScore (PEArguments * peArg ,int * mapScore0, int *mapScore1);
unsigned long long PEOCCMAPQValue(PEArguments * peArgs, unsigned long long value0, unsigned long long value1);

#endif
