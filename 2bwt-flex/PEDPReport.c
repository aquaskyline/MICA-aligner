//
//    PEDPReport.c
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

#include "PEDPReport.h"

unsigned long long PEDPOCCDumpOneOrphanAlignment(PEArguments * peArguments, SRAOccurrence * occ_found, DPOccurrence * dpOcc_found, int dpBase) {

    DPArguments * dpArguments = peArguments->dpArguments;
    DPOCCCollector * dpOccCollector = dpArguments->dpOccCollector;
    
    // Report the PE-DP Occurrence
    DPOCCAddOccurrenceToBuckets(dpOccCollector,dpOcc_found);
    
    // Report the PE-SRA Occurrence. The type being the dpBase.
    PEOCCAddTextPositionToCache(peArguments,dpBase,
                                occ_found->ambPosition,
                                occ_found->strand,
                                occ_found->matchLen,
                                occ_found->mismatchCount,0,
                                occ_found->errors);
    return 1;
}

unsigned long long PEDPOCCDumpOneSeedAlignment(PEArguments * peArguments, DPOccurrence * dpOcc_1, DPOccurrence * dpOcc_2) {

    DPArguments * dpArguments = peArguments->dpArguments;
    DPOCCCollector * dpOccCollector = dpArguments->dpOccCollector;
    
    // Report the PE-DP Occurrence
    DPOCCAddOccurrenceToBuckets(dpOccCollector,dpOcc_1);
    DPOCCAddOccurrenceToBuckets(dpOccCollector,dpOcc_2);
    
    // Report the PE-SRA Occurrence. The type being the dpBase.
    PEOCCAddTextPositionToCache(peArguments,SRAOCC_TYPE_PE_DP_SEED_OCC,0,0,0,0,0,NULL);
    return 1;
}

unsigned long long PEDPOCCDumpOneSingleEndAlignment(PEArguments * peArguments, DPOccurrence * dpOcc_1) {

    DPArguments * dpArguments = peArguments->dpArguments;
    DPOCCCollector * dpOccCollector = dpArguments->dpOccCollector;

    // Report the PE-DP Occurrence
    DPOCCAddOccurrenceToBuckets(dpOccCollector,dpOcc_1);

    // Report the PE-SRA Occurrence. The type being the dpBase.
    PEOCCAddTextPositionToCache(peArguments,SRAOCC_TYPE_PE_DP_SINGLE_END_OCC,0,0,0,0,0,NULL);
    return 1;
}

