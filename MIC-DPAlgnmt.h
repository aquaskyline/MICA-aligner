//
//    MIC-DPAlgnmt.h
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

#ifndef __MIC_DPALGNMT_H__
#define __MIC_DPALGNMT_H__

#include "2bwt-flex/SRACore.h"
#include "MIC-PEAlgnmt.h"
#include "MIC-DPCore.h"

#define MIC_DP_MAX_ERROR_COUNT          (SRA_MAX_READ_LENGTH/2)
#define MIC_DP_OUTPUT_MAX_ALIGNMENT     1000
#define MIC_DP_OUTPUT_SIZE_PER_THREAD   10000

#define MIC_DP_STATUS_BREACH_LIMIT      -1
#define MIC_DP_STATUS_ERROR             -2
#define MIC_DP_STATUS_TOO_MANY_RESULT   -3

// -----------------------------
// Debug Flag
// -----------------------------
// Define to print output message on stdout
// Not advisable for data set larger than 1
//#define MIC_DP_DEBUG_PRINT_ORPHAN_MATCHING

// -----------------------------
// Compile-time Configurable #1
// Ticket#78 [MICA] Default-DP not as sensitive as SOAP3, even SOAP2
// -----------------------------
// Enable the following flag to force the DP results to fall within the
// insertion size by imposing extra checking prior to report
//
//sh : TODO  dont enable it, not tested, looks like have bug
//#define MIC_DP_FORCE_INSERTION_SIZE

typedef struct MICDPOccurrence {
    unsigned int ambPosition;
    uint8_t strand;
    uint16_t matchElemsCount;
    DPMatchElem matchElems[MIC_DP_MAX_ERROR_COUNT];
    uint16_t matchLen;
    uint16_t score;
} MICDPOccurrence; 

__attribute__((target(mic)))
char * MICDPExtractCodeFromHSP(char * buffer, const unsigned startIndex, const unsigned length, const HSP * hsp);

__attribute__((target(mic)))
int _MICDPPopulateOccurrence(MICDPOccurrence * dpoToWrite, DPWorkMIC * dpWork, unsigned int startPos, int dpStrand);

__attribute__((target(mic)))
int MICDPOrphanAlignment(MICDPOccurrence * dpOutputPtr, DPWorkMIC * dpWork, const MICPEArguments * peArgs, 
        const HSP * hsp, double InputSGAOrphanTriggerTF);

__attribute__((target(mic)))
void MICDPOccurrenceConvert(MICDPOccurrence  * source, DPOccurrence * destination);

#endif

