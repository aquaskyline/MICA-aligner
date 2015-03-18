//
//    SRADPReport.c
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
   
   Date   : 12th April 2012
   Author : Edward MK Wu
   Change : New file.

*/
/////////////////////////////////////////////////////

#include "SRADPReport.h"

unsigned long long DPOCCDumpOneSeedAlignment(SRAArguments * sraArguments, DPOccurrence * dpOcc_1) {

    DPArguments * dpArguments = sraArguments->dpArguments;
    DPOCCCollector * dpOccCollector = dpArguments->dpOccCollector;
    
    // Report the PE-DP Occurrence
    DPOCCAddOccurrenceToBuckets(dpOccCollector,dpOcc_1);
    
    // Report the PE-SRA Occurrence. The type being the dpBase.
    OCCAddTextPositionToCache(sraArguments,SRAOCC_TYPE_DP_SEED_OCC,0,0,0,0,0,NULL);
    
    return 1;
}
