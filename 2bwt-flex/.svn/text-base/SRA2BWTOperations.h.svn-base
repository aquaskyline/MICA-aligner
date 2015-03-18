//
//    SRA2BWTOperations.h
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

   Date   : 8th February 2011
   Author : Edward MK Wu
   Change : Fixed segmentation fault in HSPLoad and CE func.

   Date   : 19th June 2011
   Author : Edward MK Wu
   Change : Packaging 2BWT library as a separate product.
             Thus, changing all references to 2bwt lib to subdirectory.

*/
/////////////////////////////////////////////////////

#ifndef __SRA_2BWT_OPERATIONS_H__
#define __SRA_2BWT_OPERATIONS_H__

#include "SRA2BWTMdl.h"
#include "SRA2BWTReport.h"
#include "SRA2BWTCheckAndExtend.h"

void SRAAlgnmtInputPrint(SRAArguments * aArgs);

// BWTxxxModelxxx functions are fully generalised BWT search algorithm that searchs for reads contain any number of edit/mismatch
// However, calling these functions requires user to define themselves a 'searching model'.
// The searching model requires each BWT step to be defined. The search algorithm will then follow the defined model.

// BWTMismatchModelAnyDirection_CE matches steps with check and extend mechanism.
// It allows starting off CE in the middle of a step and recursive CE until SRACase completes.
unsigned long long BWTMismatchModelAnyDirection_CE(SRAArguments * aArgs, int i, uint8_t mismatchInserted, 
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occMismatch, uint8_t nbMismatch, int occQuality);


// BWTExactModelForward_Lookup lookup your pattern in lookup table, bi-directional and assuming forward
unsigned long long BWTExactModelForward_Lookup(SRAArguments * aArgs,
                                    SRACase * alignmentCase, int stepInCase,
                                    unsigned long long * saRanges);

// BWTExactModelBackward_Lookup lookup your pattern in lookup table, single direction and assuming backward
unsigned long long BWTExactModelBackward_Lookup(SRAArguments * aArgs,
                                    SRACase * alignmentCase, int stepInCase,
                                    unsigned long long * saRanges);
                                    
// BWTExactModelBackward_Lookup lookup your pattern in lookup table, bi-directional and assuming forward
unsigned long long BWTExactModelBackwardAnyDirection_Lookup(SRAArguments * aArgs,
                                    SRACase * alignmentCase, int stepInCase,
                                    unsigned long long * saRanges);
// BWTExactModelBackward matches pattern on text without using any other aux, e.g. lookup table.
unsigned long long BWTExactModelBackward(SRAArguments * aArgs, int i, uint8_t errorInserted,
                                    SRACase * alignmentCase, int stepInCase,
                                    unsigned long long * saRanges,
                                    uint8_t occError, uint8_t nbMismatch, int occQuality);
                                    
// BWTMismatchModelAnyDirection matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
unsigned long long BWTExactModelAnyDirection(SRAArguments * aArgs, int i, uint8_t errorInserted,
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occError, uint8_t nbMismatch, int occQuality);
                                    
// BWTMismatchModelAnyDirection matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
unsigned long long BWTMismatchModelAnyDirection(SRAArguments * aArgs, int i, uint8_t mismatchInserted,
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occMismatch, uint8_t nbMismatch, int occQuality);

// BWTMismatchModelBackward matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
unsigned long long BWTMismatchModelBackward(SRAArguments * aArgs, int i, uint8_t mismatchInserted,
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occMismatch, uint8_t nbMismatch, int occQuality);


// BWTEditModelAnyDirection matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
unsigned long long BWTEditModelAnyDirection(SRAArguments * aArgs, int i, uint8_t editInserted,
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occEdit, uint8_t nbMismatch, int occQuality);
                                    
// BWTEditModelBackward matches pattern on text without using any other aux, e.g. lookup table, with mismatches.
unsigned long long BWTEditModelBackward(SRAArguments * aArgs, int i, uint8_t editInserted,
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occEdit, uint8_t nbMismatch, int occQuality);

unsigned long long BWTModelSwitchAnyDirection(SRAArguments * aArgs, int i, uint8_t errorInserted,
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occError, uint8_t nbMismatch, int occQuality);

unsigned long long BWTModelSwitchBackward(SRAArguments * aArgs,  int i, uint8_t errorInserted,
                                    SRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occError, uint8_t nbMismatch, int occQuality);

#endif
