//
//    PEArguments.h
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
   
   Date   : 4th July 2013
   Author : Edward MK Wu
   Change : New file.
   
*/
/////////////////////////////////////////////////////

#include "PEArguments.h"

void PEARGPopulateSRAARG(PEArguments * peArgs, SRAArguments * readArg, SRAArguments * readArg_neg, SRAArguments * mateArg, SRAArguments * mateArg_neg) {
    peArgs->sraArgsPos   = readArg;
    peArgs->sraArgsNeg   = readArg_neg;
    
    peArgs->sraArgsPos_mate   = mateArg;
    peArgs->sraArgsNeg_mate   = mateArg_neg;
}

PEArguments * PEARGConstruct(BWT * bwt, HSP * hsp) {

    PEArguments * peArgs = (PEArguments*) malloc( sizeof(PEArguments) );
    
    // For each elements in PEArguments
    peArgs->sraArgsPos   = NULL;
    peArgs->sraArgsNeg   = NULL;
    
    peArgs->sraArgsPos_mate   = NULL;
    peArgs->sraArgsNeg_mate   = NULL;
    
    peArgs->PEAlgnmtInput            = PEInputConstruct(bwt,hsp);
    peArgs->PEAlgnmtOutput           = PEOutputConstruct();
    peArgs->PEAlgnmtStats            = PESTATConstruct();
    peArgs->PEAlgnmtOutput_OccHandle = SRAAlgnmtOutputConstruct();
    
    peArgs->dpArguments   = DPARGConstruct();
    peArgs->pedpSetting   = PEDPSettingConstruct();
    
    peArgs->MapqCalc      = MAPQCalculatorCreate();
    
    return peArgs;
}

// PEARGMakeMate function creates a copy of PEArguments based on a populated one
// the slave will share only the static data like PEAlgnmtInput, pedpSetting.
PEArguments * PEARGMakeMate(PEArguments * source) {
    PEArguments * destination = (PEArguments*) malloc( sizeof(PEArguments) );
    memcpy(destination,source,sizeof(PEArguments));
    
    // Specific structure that is not supposed to be shared.
    destination->PEAlgnmtOutput           = PEOutputConstruct();
    destination->PEAlgnmtStats            = PESTATConstruct();
    destination->PEAlgnmtOutput_OccHandle = SRAAlgnmtOutputConstruct();
    
    destination->dpArguments   = DPARGMakeMate(source->dpArguments);
    destination->MapqCalc      = MAPQCalculatorCreate();
    
    return destination;
}

// PEARGMakeClone function creates an exact copy of PEArguments based on a populated one
PEArguments * PEARGMakeClone(PEArguments * source) {
    PEArguments * destination = (PEArguments*) malloc( sizeof(PEArguments) );
    memcpy(destination,source,sizeof(PEArguments));    
    return destination;
}

void PEARGCloneFree(PEArguments * peArgs) {
    free(peArgs);
}

void PEARGMateFree(PEArguments * peArgs) {
    if (peArgs->PEAlgnmtOutput!=NULL) {
        PEOutputFree(peArgs->PEAlgnmtOutput);
    }
    if (peArgs->PEAlgnmtStats!=NULL) {
        PESTATFree(peArgs->PEAlgnmtStats);
    }
    if (peArgs->PEAlgnmtOutput_OccHandle!=NULL) {
        SRAAlgnmtOutputFree(peArgs->PEAlgnmtOutput_OccHandle);
    }
    if (peArgs->dpArguments!=NULL) {
        DPARGMateFree(peArgs->dpArguments);
    }
    if (peArgs->MapqCalc!=NULL) {
        MAPQCalculatorFree(peArgs->MapqCalc);
    }
    
    free(peArgs);
}


void PEARGFree(PEArguments * peArgs) {
    if (peArgs->PEAlgnmtInput!=NULL) {
        PEInputFree(peArgs->PEAlgnmtInput);
    }
    if (peArgs->PEAlgnmtOutput!=NULL) {
        PEOutputFree(peArgs->PEAlgnmtOutput);
    }
    if (peArgs->PEAlgnmtStats!=NULL) {
        PESTATFree(peArgs->PEAlgnmtStats);
    }
    if (peArgs->PEAlgnmtOutput_OccHandle!=NULL) {
        SRAAlgnmtOutputFree(peArgs->PEAlgnmtOutput_OccHandle);
    }
    if (peArgs->dpArguments!=NULL) {
        DPARGFree(peArgs->dpArguments);
    }
    if (peArgs->pedpSetting!=NULL) {
        PEDPSettingFree(peArgs->pedpSetting);
    }
    if (peArgs->MapqCalc!=NULL) {
        MAPQCalculatorFree(peArgs->MapqCalc);
    }
    
    free(peArgs);
}
