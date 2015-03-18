//
//    DPArguments.c
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
   
   Date   : 16th June 2013
   Author : Edward MK Wu
   Change : New file.

*/
/////////////////////////////////////////////////////

#include <stdlib.h>
#include <string.h>

#include "DPArguments.h"

DPArguments * DPARGConstruct() {
    DPArguments * dpArgs = (DPArguments*) malloc ( sizeof(DPArguments) );
    
    dpArgs->dpWork = DPWorkCreate();
    dpArgs->dpScores = DPScoresConstruct();
    dpArgs->dpOccCollector = DPOCCCreate(DPOCC_FLOOD_TYPE_FLUSH);
    
    return dpArgs;
}

// DPARGMakeClone function creates a copy of DPArguments based on a populated one
// the slave will share only the static data like dpScores.
DPArguments * DPARGMakeMate(DPArguments * source) {
    DPArguments * destination = (DPArguments*) malloc ( sizeof(DPArguments) );
    memcpy(destination,source,sizeof(DPArguments));
    
    destination->dpWork = DPWorkCreate();
    destination->dpOccCollector = DPOCCCreate(DPOCC_FLOOD_TYPE_FLUSH);
    
    return destination;
}

void DPARGMateFree(DPArguments * dpArgs) {
    if (dpArgs->dpWork!=NULL) {
        DPWorkFree(dpArgs->dpWork);
        dpArgs->dpWork=NULL;
    }
    if (dpArgs->dpOccCollector!=NULL) {
        DPOCCFree(dpArgs->dpOccCollector);
        dpArgs->dpOccCollector=NULL;
    }
    free(dpArgs);
}

void DPARGFree(DPArguments * dpArgs) {
    if (dpArgs->dpWork!=NULL) {
        DPWorkFree(dpArgs->dpWork);
        dpArgs->dpWork=NULL;
    }
    if (dpArgs->dpScores!=NULL) {
        DPScoresFree(dpArgs->dpScores);
        dpArgs->dpScores=NULL;
    }
    if (dpArgs->dpOccCollector!=NULL) {
        DPOCCFree(dpArgs->dpOccCollector);
        dpArgs->dpOccCollector=NULL;
    }
    free(dpArgs);
}

