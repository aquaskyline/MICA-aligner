//
//    MIC-SRAAlgnmt.h
//
//    mica
//
//    Copyright (C) 2013, HKU
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

#ifndef __MIC_SRAALGNMT_H__
#define __MIC_SRAALGNMT_H__

#include "MIC-SRA2BWTOperations.h"
#include "MIC-SRA2BWTMdl.h"
#include "MIC-SRAArguments.h"

__attribute__((target(mic)))
unsigned long long MICProcessStage(MICSRAArguments * micArgs, CPTSRAModel * cpModel, int * caseId);

__attribute__((target(mic)))
unsigned long long MICProcessStageDoubleStrand(MICSRAArguments * micArgs, CPTSRAModel * cpModel, CPTSRAModel * cpModel_neg, int * caseId);

__attribute__((target(mic)))
unsigned long long MICProcessReadDoubleStrand(MICSRAArguments * micArgs, CPTSRAModel * cpModel, CPTSRAModel * cpModel_neg);

__attribute__((target(mic)))
unsigned long long MICSeedAlgnmtDoubleStrand(MICSRAArguments * micArgs, CPTSRAModel * cpSeedModel, CPTSRAModel * cpSeedModel_neg,
                                                    uint16_t seedLength, uint16_t seedUniqeLength, int occLimit);
#endif

