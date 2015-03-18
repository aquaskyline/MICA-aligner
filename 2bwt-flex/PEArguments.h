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
   
   Date   : 8th May 2011
   Author : Edward MK Wu
   Change : New file.
   
*/
/////////////////////////////////////////////////////

#ifndef __PE_ARGUMENTS_H__
#define __PE_ARGUMENTS_H__

#include "PECore.h"
#include "SRAArguments.h"
#include "PEDPArguments.h"

typedef struct PEArguments{

    SRAArguments * sraArgsPos;
    SRAArguments * sraArgsNeg;
    
    SRAArguments * sraArgsPos_mate;
    SRAArguments * sraArgsNeg_mate;
    
    PEInput * PEAlgnmtInput;
    PEOutput * PEAlgnmtOutput;
    // Accumulative Statistics. User responsibility to reset if needed.
    PEStats * PEAlgnmtStats;
    //Handle output from the pair-end alignment
    SRAOutput * PEAlgnmtOutput_OccHandle;
    
    DPArguments * dpArguments;
    PEDPSetting * pedpSetting;
    
    MAPQCalculator * MapqCalc;
} PEArguments;

PEArguments * PEARGConstruct(BWT * bwt, HSP * hsp);
void PEARGPopulateSRAARG(PEArguments * peArgs, SRAArguments * readArg, SRAArguments * readArg_neg, SRAArguments * mateArg, SRAArguments * mateArg_neg);
PEArguments * PEARGMakeMate(PEArguments * source);
void PEARGMateFree(PEArguments * peArgs);
PEArguments * PEARGMakeClone(PEArguments * source);
void PEARGCloneFree(PEArguments * peArgs);
void PEARGFree(PEArguments * peArgs);

#endif
