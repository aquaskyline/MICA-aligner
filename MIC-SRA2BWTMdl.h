//
//    MIC-SRA2BWTMdl.h
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

#ifndef __MIC_SRA2BWTMDL_H__
#define __MIC_SRA2BWTMDL_H__

#include "2bwt-flex/SRA2BWTMdl.h"

#define MAX_NUM_OF_CPT_SRA_STEPS                         8
#define MAX_NUM_OF_CPT_SRA_CASES                         18

// -----------------------------
// Compile-time Configurable #1
// -----------------------------
// Define this token to disable the use of lookup table
// in MIC-SRA alignment. Please note that any related
// memory operations are still carried out.
//#define MIC_2BWT_MODEL_DISABLE_LOOKUP

typedef struct CPTSRAStep {
    uint8_t type;
    int start;
    int end;

    uint8_t ErrorType;
    uint8_t MinError;
    uint8_t MaxError;

    char step;         //computed
    int len;           //computed
} CPTSRAStep;

typedef struct CPTSRACase {
    int id;
    uint8_t type;
    CPTSRAStep steps[MAX_NUM_OF_CPT_SRA_STEPS];
} CPTSRACase;

typedef struct CPTSRAModel {
    CPTSRACase     cases[MAX_NUM_OF_CPT_SRA_CASES];
} CPTSRAModel;

void MICSRA2BWTModelPopulate(CPTSRAModel * cpModel, SRAModel * sraModel);

__attribute__((target(mic)))
void MICSRADebugPrintModel(CPTSRAModel * cpModel, FILE * stream);

#endif

