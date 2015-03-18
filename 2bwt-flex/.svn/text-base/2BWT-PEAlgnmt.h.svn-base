//
//    2BWT-PEAlgnmt.h
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
   
   Date   : 11th March 2012
   Author : Edward MK Wu
   Change : New file.

*/
/////////////////////////////////////////////////////

#ifndef __2BWT_PE_ALIGNMENT__
#define __2BWT_PE_ALIGNMENT__

#include "2BWT-SRAAlgnmt.h"

#include "PEArguments.h"
#include "PE2BWTReport.h"
#include "PEDPReport.h"

// -----------------------------
// Compile-time Configurable #1
// Ticket#78 [MICA] Default-DP not as sensitive as SOAP3, even SOAP2
// -----------------------------
// Enable the following flag to force the DP results to fall within the
// insertion size by imposing extra checking prior to report
//#define PE_DP_FORCE_INSERTION_SIZE

unsigned long long PEProcessReadDoubleStrand(PEArguments * peArguments, 
                            SRAModelSet * sraModelSet, SRAModelSet * sraModelSet_Seed, SRAModelSet * sraModelSet_Extend);

unsigned int DPPEOrphanMappingOccurrences(PEArguments * peArguments, SRAOccurrence * occList_base, unsigned long long occCount_base, unsigned int readLength_base,
                            unsigned char * readCode, unsigned char * readCode_neg, unsigned int readLength, MAPQCalculator * mapqCalc,
                            int dpBase);
unsigned int DPPESeedRecoverOccurrences(PEArguments * peArguments,
                                        SRAModelSet * sraModelSet_Seed, 
                                        double SGASeedLength, double SGASeedLengthOverlap,
                                        int SGASeedLooseCriteria, int SGASeedOccLimitation);
#endif

