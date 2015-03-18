/*
 *
 *    
 *      MicMappingQuality.h
 *    Soap3(gpu)
 *
 *    Copyright (C) 2011, HKU
 *
 *    This program is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License
 *    as published by the Free Software Foundation; either version 2
 *    of the License, or (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */

#ifndef __MIC_MAPPING_QUALITY_H__
#define __MIC_MAPPING_QUALITY_H__

#include "2bwt-flex/MappingQuality.h"
#include "2bwt-flex/DPCore.h"

#include "MIC-SRAArguments.h"
#include "MIC-SRA2BWTMdl.h"
#include "MIC-SRAAlgnmt.h"
#include "MIC-PEAlgnmt.h"
#include "MIC-DPAlgnmt.h"


typedef int MappingQuality;

typedef struct PEMappingQuality {
    MappingQuality readQuality;
    MappingQuality mateQuality;
} PEMappingQuality;

// g_log_n: this should be initiate by bwase_initialize
//          in MappingQuality.h
//          the size of the array is 256
//          output buffer inside sraArgs will be used if need
__attribute__((target(mic)))
PEMappingQuality MicCalculatePEMappingQuality(
        MICSRAArguments * readSraArgs, 
        MICSRAArguments * mateSraArgs,
        MICPEArguments * peArgs,
        MICDPOccurrence * dpOccurrences,
        CPTSRAModel * cpPModels,
        CPTSRAModel * cpNModels,
        int * g_log_n);

#endif
