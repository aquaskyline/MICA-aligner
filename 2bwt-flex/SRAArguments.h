//
//    SRAArguments.h
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

#ifndef __SRA_ARG_H__
#define __SRA_ARG_H__

#include "DPArguments.h"

#include "SRA2BWTMdl.h"
#include "SRAOutputFile.h"

#include "MappingQuality.h"

#define MAPQ_CALC_FIGURE_INITIALISED 0
#define MAPQ_CALC_FIGURE_READY       1

// MAPQCalculator is strictly implemented for the calculation of 
// MAPQ values from SRA / PE / DP results. This structures will be
// read by the output handler in *Report library.
typedef struct MAPQCalculator {
    uint8_t rankedMismatch[2];
    int16_t rankedScore[2];
    unsigned int rankedCount[2];
    int * g_log_n;
    int status;
} MAPQCalculator;

typedef struct SRAResultCount{
    //WithError stores the number of occurrences found with the corresponding number of mismatch/edits
    unsigned int WithError[MAX_NUM_OF_ERROR+1];

    //Position retrieved by SA, CE or HOCC
    unsigned int RetrievedBySa;
    unsigned int RetrievedByCe;
    unsigned int RetrievedByHocc;
} SRAResultCount;

typedef struct SRAWorkingMemory {
    // Working memory item that is initialised as they
    // are used.
    unsigned long long Shared_ReadPositions[MAX_CE_THRESHOLD];
    uint8_t Shared_AccMismatch[MAX_CE_THRESHOLD];
    int Shared_AccQualities[MAX_CE_THRESHOLD];
    SRAError Shared_AccErrorsPos[MAX_NUM_OF_ERROR];
    unsigned long long BestAlignmentPosition;
    uint8_t BestAlignmentStrand;
    uint8_t BestAlignmentNumOfError;
    int BestAlignmentQuality;
    SRAError BestAlignmentErrorPos[MAX_NUM_OF_ERROR];


    // Working memory item that MUST be initialised independently
    // before invoking the 2BWT alignment model search engine.
    int Shared_HeadTrim;
    int Shared_TailTrim;
    
    //IsOpened      1 : Yes
    //              0 : No
    int8_t IsOpened;

    //IsClosed      1 : Yes
    //              0 : No
    int8_t IsClosed;

    //IsUnique      1 : Yes
    //              0 : No
    int8_t IsUnique;

    //Best alignment position
    //IsBest        0 : Unknown
    //              1 : Exists
    int8_t IsBest;

} SRAWorkingMemory;

// For any arguments added to this struct, constructor functions
// must be updated to reflect the changes
typedef struct SRAArguments {
    SRAQueryInfo * QueryInfo;
    SRASetting * AlgnmtSetting;
    SRAIndex * AlgnmtIndex;
    SRAWorkingMemory * AlgnmtMemory;
    SRAResultCount * AlgnmtStats;
    SRAOutput * AlgnmtOutput;
    DPArguments * dpArguments;
    MAPQCalculator * MapqCalc;
} SRAArguments;



SRAArguments * SRAARGConstruct();
SRAArguments * SRAARGMakeClone(SRAArguments * source);
SRAArguments * SRAARGMakeSlave(SRAArguments * source);
SRAArguments * SRAARGMakeMate(SRAArguments * source);
void SRAARGFree ( SRAArguments * aArgs );
void SRAARGSlaveFree ( SRAArguments * aArgs );
void SRAARGMateFree ( SRAArguments * aArgs );
void SRAARGCloneFree ( SRAArguments * aArgs );

SRAResultCount * SRAResultCountConstruct();
void SRAResultCountInitialise(SRAResultCount * aStats);
void SRAResultCountFree(SRAResultCount * aStats);
SRAWorkingMemory * SRAWorkingMemoryConstruct();
__attribute__((target(mic)))
void SRAWorkingMemoryInitialise(SRAWorkingMemory * aMem);
void SRAWorkingMemoryFree(SRAWorkingMemory * aMem);

MAPQCalculator * MAPQCalculatorCreate();
void MAPQCalculatorInitialise(MAPQCalculator * mapqCalc);
void MAPQCalculatorOverlay(MAPQCalculator * source, MAPQCalculator * destination);
void MAPQCalculatorFree(MAPQCalculator * mapqCalc);

////////////////////////////////////////////////////////
// DEBUG FUNCTIONS
////////////////////////////////////////////////////////

void SRAPrintResultCount(SRAResultCount * aStats);
void MAPQPrintParameters(MAPQCalculator * mapqCalc);

#endif
