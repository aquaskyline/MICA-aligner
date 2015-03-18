//
//    SRAArguments.c
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

#include "SRAArguments.h"

SRAArguments * SRAARGConstruct() {
    SRAArguments * aArgs = (SRAArguments*) malloc( sizeof(SRAArguments) );
    
    // For each elements in SRAArguments
    aArgs->QueryInfo     = SRAQueryInfoConstruct();
    aArgs->AlgnmtSetting = SRASettingConstruct();
    aArgs->AlgnmtIndex   = SRAIndexConstruct();
    aArgs->AlgnmtMemory  = SRAWorkingMemoryConstruct();
    aArgs->AlgnmtStats   = SRAResultCountConstruct();
    aArgs->AlgnmtOutput  = SRAAlgnmtOutputConstruct();
    aArgs->dpArguments   = DPARGConstruct();
    aArgs->MapqCalc      = MAPQCalculatorCreate();
    
    return aArgs;
}

// SRAARGMakeClone function creates a copy of SRAArgument based on a populated one
// the clone will share everything with the source
SRAArguments * SRAARGMakeClone(SRAArguments * source) {

    SRAArguments * destination = (SRAArguments*) malloc( sizeof(SRAArguments) );
    memcpy(destination,source,sizeof(SRAArguments));
    
    return destination;
}

// SRAARGMakeSlave function creates a copy of SRAArgument based on a populated one
// the slave will share pretty much everything with the source except QueryInfo
SRAArguments * SRAARGMakeSlave(SRAArguments * source) {

    SRAArguments * destination = (SRAArguments*) malloc( sizeof(SRAArguments) );
    memcpy(destination,source,sizeof(SRAArguments));
    
    // Specific structure that is not supposed to be shared.
    destination->QueryInfo     = SRAQueryInfoConstruct();
    
    return destination;
}


// SRAARGMakeMate function creates a copy of SRAArgument based on a populated one
// the slave will share only the static data like AlgnmtSetting, AlgnmtIndex.
SRAArguments * SRAARGMakeMate(SRAArguments * source) {

    SRAArguments * destination = (SRAArguments*) malloc( sizeof(SRAArguments) );
    memcpy(destination,source,sizeof(SRAArguments));
    
    // Specific structure that is not supposed to be shared.
    destination->QueryInfo     = SRAQueryInfoConstruct();
    destination->AlgnmtMemory  = SRAWorkingMemoryConstruct();
    destination->AlgnmtStats   = SRAResultCountConstruct();
    destination->AlgnmtOutput  = SRAAlgnmtOutputConstruct();
    destination->dpArguments   = DPARGMakeMate(source->dpArguments);
    destination->MapqCalc      = MAPQCalculatorCreate();
    
    return destination;
}


void SRAARGFree ( SRAArguments * aArgs ) {

    // For each elements in SRAArguments
    SRAQueryInfoFree(aArgs->QueryInfo);
    SRASettingFree(aArgs->AlgnmtSetting);
    SRAIndexFree(aArgs->AlgnmtIndex);
    SRAWorkingMemoryFree(aArgs->AlgnmtMemory);
    SRAResultCountFree(aArgs->AlgnmtStats);
    SRAAlgnmtOutputFree(aArgs->AlgnmtOutput);
    DPARGFree(aArgs->dpArguments);
    MAPQCalculatorFree(aArgs->MapqCalc);
    
    free(aArgs);
}


void SRAARGCloneFree ( SRAArguments * aArgs ) {
    free(aArgs);
}

void SRAARGSlaveFree ( SRAArguments * aArgs ) {

    // For each elements in SRAArguments
    SRAQueryInfoFree(aArgs->QueryInfo);
    
    free(aArgs);
}

void SRAARGMateFree ( SRAArguments * aArgs ) {

    // For each elements in SRAArguments
    SRAQueryInfoFree(aArgs->QueryInfo);
    SRAWorkingMemoryFree(aArgs->AlgnmtMemory);
    SRAResultCountFree(aArgs->AlgnmtStats);
    SRAAlgnmtOutputFree(aArgs->AlgnmtOutput);
    DPARGMateFree(aArgs->dpArguments);
    MAPQCalculatorFree(aArgs->MapqCalc);
    
    free(aArgs);
}


SRAResultCount * SRAResultCountConstruct() {
    SRAResultCount * aStats = (SRAResultCount*) malloc ( sizeof(SRAResultCount) );
    
    SRAResultCountInitialise(aStats);
    
    return aStats;
}

void SRAResultCountInitialise(SRAResultCount * aStats) {
    memset(aStats,0,sizeof(SRAResultCount));
}

void SRAResultCountFree(SRAResultCount * aStats) {
    free(aStats);
}

MAPQCalculator * MAPQCalculatorCreate() {
    MAPQCalculator * mapqCalc = (MAPQCalculator*) malloc ( sizeof(MAPQCalculator) );
    mapqCalc->g_log_n = (int *) malloc(sizeof(int) * 256) ;
    bwase_initialize(mapqCalc->g_log_n);
    MAPQCalculatorInitialise(mapqCalc);
    return mapqCalc;
}

void MAPQCalculatorInitialise(MAPQCalculator * mapqCalc) {
    int i;
    for (i=0;i<2;i++) {
        mapqCalc->rankedMismatch[i] = 0;
        mapqCalc->rankedScore[i] = 0;
        mapqCalc->rankedCount[i] = 0;
    }
    mapqCalc->status = MAPQ_CALC_FIGURE_INITIALISED;
}

void MAPQCalculatorOverlay(MAPQCalculator * source, MAPQCalculator * destination) {
    int i;
    for (i=0;i<2;i++) {
        destination->rankedMismatch[i] = source->rankedMismatch[i];
        destination->rankedScore[i] = source->rankedScore[i];
        destination->rankedCount[i] = source->rankedCount[i];
    }
}

void MAPQCalculatorFree(MAPQCalculator * mapqCalc) {
    free(mapqCalc->g_log_n);
	free(mapqCalc);
}

SRAWorkingMemory * SRAWorkingMemoryConstruct() {
    SRAWorkingMemory * aMem = (SRAWorkingMemory*) malloc ( sizeof(SRAWorkingMemory) );
    SRAWorkingMemoryInitialise(aMem);
    return aMem;
}

__attribute__((target(mic)))
void SRAWorkingMemoryInitialise(SRAWorkingMemory * aMem) {
    aMem->IsOpened = 0;
    aMem->IsClosed = 0;
    aMem->IsUnique = 0;
    aMem->IsBest = 0;
    aMem->Shared_HeadTrim = 0;
    aMem->Shared_TailTrim = 0;
}

void SRAWorkingMemoryFree(SRAWorkingMemory * aMem) {
    free(aMem);
}    

////////////////////////////////////////////////////////
// DEBUG FUNCTIONS
////////////////////////////////////////////////////////

void SRAPrintResultCount(SRAResultCount * aStats) {
    unsigned int count = 0;
    int i;
    for (i=0;i<=MAX_NUM_OF_ERROR;i++) {
        if (aStats->WithError[i]>0)
            printf("  - %2u mismatch/edit alignments  = %u\n", i, aStats->WithError[i]);
        count+=aStats->WithError[i];
    }
    printf("Total                            = %u\n", count);
}

void MAPQPrintParameters(MAPQCalculator * mapqCalc) {
    printf("Optimal Occurrences     Count = %u      Score = %u (%u-mm)\n",
        mapqCalc->rankedCount[0],
        mapqCalc->rankedScore[0],
        mapqCalc->rankedMismatch[0]);
    printf("Sub-Optimal Occurrences Count = %u      Score = %u (%u-mm)\n",
        mapqCalc->rankedCount[1],
        mapqCalc->rankedScore[1],
        mapqCalc->rankedMismatch[1]);
    printf("\n");
}
