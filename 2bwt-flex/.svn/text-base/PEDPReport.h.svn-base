//
//    PEDPReport.h
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
   
   Date   : 26th February 2012
   Author : Edward MK Wu
   Change : New file.

*/
/////////////////////////////////////////////////////

#ifndef __PE_DP_REPORT_H__
#define __PE_DP_REPORT_H__

#include "PE2BWTReport.h"
#include "PEArguments.h"
#include "DPArguments.h"

//PEDPOCCDumpOneOrphanAlignment reports a pair of occurrences where one of them
//being result from PE-SRA and the other from PE-DP.
unsigned long long PEDPOCCDumpOneOrphanAlignment(PEArguments * peArguments, SRAOccurrence * occ_found, DPOccurrence * dpOcc_found, int dpBase);

//PEDPOCCDumpOneSeedAlignment reports a pair of dynamic programming occurrences from PE-DP.
unsigned long long PEDPOCCDumpOneSeedAlignment(PEArguments * peArguments, DPOccurrence * dpOcc_1, DPOccurrence * dpOcc_2);

//PEDPOCCDumpOneSingleEndAlignment reports a single dynamic programming occurrences from PE-DP.
unsigned long long PEDPOCCDumpOneSingleEndAlignment(PEArguments * peArguments, DPOccurrence * dpOcc_1);

#endif

