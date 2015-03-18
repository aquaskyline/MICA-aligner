//
//    MIC-SRAArguments.c
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

#include "MIC-SRAArguments.h"

__attribute__((target(mic)))
MICSRAArguments * MICSRAARGConstruct() {
    MICSRAArguments * micArgs = (MICSRAArguments*) malloc( sizeof(MICSRAArguments) );
    
    micArgs->readCode = NULL;
    micArgs->readCode_Complt = NULL;
    micArgs->seedLength = 0;
    micArgs->readLength = 0;
    micArgs->readStrand = 0;
    micArgs->maxNBMismatch = 0;
    micArgs->seedOffset = 0;
    micArgs->seedOffset_Complt = 0;
    micArgs->outputBlock = NULL;
    micArgs->outputStatus = NULL;
    micArgs->occCount = NULL;
    micArgs->metaBlock = NULL;
    micArgs->metaCount = NULL;
    
    micArgs->outputType = SRA_REPORT_ALL;
    
    micArgs->bwt = NULL;
    micArgs->rev_bwt = NULL;
    
    micArgs->lt = NULL;
    micArgs->rlt = NULL;
    
    return micArgs;
}

__attribute__((target(mic)))
void MICSRAARGFree(MICSRAArguments * micArgs) {
    free(micArgs);
}

__attribute__((target(mic)))
MICSRAArguments * MICSRAARGMakeCopy(MICSRAArguments * micArgs) {
    MICSRAArguments * copyArgs = (MICSRAArguments*) malloc( sizeof(MICSRAArguments) );
    memcpy(copyArgs,micArgs,sizeof(MICSRAArguments));
    return copyArgs;
}

__attribute__((target(mic)))
MICSRAArguments * MICSRAARGCopy(MICSRAArguments * destination, MICSRAArguments * source) {
    memcpy(destination,source,sizeof(MICSRAArguments));
    return destination;
}

void MICSRAMETADebugPrint(MICSRAOccMetadata * meta) {
    int i;
    printf("Meta (%u) %c #M=%u |M|=%u\n",meta->numPayload,(meta->strand==0)?'+':'-',meta->numOfErr,meta->matchLen);
    for (i=0;i<meta->numOfErr&&i<MAX_NUM_OF_ERROR;i++) {
        printf("  - @%u TYPE %d\n",meta->errors[i].position,meta->errors[i].type);
    }
}