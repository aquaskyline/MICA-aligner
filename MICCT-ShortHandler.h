//
//    MICCT-ShortHandler.h
//
//    mica
//
//    Copyright (C) 2014, HKU
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

#ifndef __MICCT_SHORTHANDLER_H__
#define __MICCT_SHORTHANDLER_H__

///////////////////////////////////////
// Including Standard Libraries
///////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>

///////////////////////////////////////
// Including MIC Libraries
///////////////////////////////////////
#include "MIC-SRA2BWTMdl.h"
#include "MIC-SRAAlgnmt.h"
#include "MIC-PEAlgnmt.h"
#include "MIC-DPAlgnmt.h"
#include "MICA-PE-ReadThread.h"
#include "MICA-PE-WriteThread.h"
#include "MicMappingQuality.h"
#include "MICDeepDP.h"
#include "MemMan.h"

__attribute__((target(mic)))
void MICSHProcessRead(HSP * hsp,
                     MICSRAArguments * readMicArgs, MICSRAArguments * mateMicArgs,
                     CPTSRAModel * cpPModels, CPTSRAModel * cpNModels,
                     PEInput * peAlgnmtInput, MICPEArguments * peArgs,
                     
                     MICSRAArguments * readSeedMicArgs, MICSRAArguments * mateSeedMicArgs,
                     CPTSRAModel * cpPModels_seed, CPTSRAModel * cpNModels_seed,
                     MICPEArguments * peSeedArgs,
                     unsigned int * peSeedOutputPtr, MICSRAOccMetadata * peSeedMetaPtr,
                     
                     PEDPSetting * pedpSetting, DPWorkMIC * dpWork,
                     
                     unsigned int * outputPtr, uint8_t * outputBufferStatus, uint16_t * occCount,
                     MICDPOccurrence * dpOcc, unsigned dpOccCount, uint32_t * dpOccCountInc);

#endif

