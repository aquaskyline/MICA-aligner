//
//    MIC-SRA2BWTOperations.h
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

#ifndef __MIC_SRA2BWTOPERATIONS_H__
#define __MIC_SRA2BWTOPERATIONS_H__

#include "MIC-SRAArguments.h"
#include "MIC-SRA2BWTMdl.h"



// BWTExactModelForward_Lookup lookup your pattern in lookup table, bi-directional and assuming forward
__attribute__((target(mic)))
unsigned long long MICBWTExactModelForward_Lookup(MICSRAArguments * micArgs,
                                    CPTSRACase * alignmentCase, int stepInCase,
                                    unsigned long long * saRanges);
__attribute__((target(mic)))
unsigned long long MICBWTExactModelBackward_Lookup(MICSRAArguments * micArgs,
                                    CPTSRACase * alignmentCase, int stepInCase,
                                    unsigned long long * saRanges);
__attribute__((target(mic)))
unsigned long long MICBWTExactModelBackwardAnyDirection_Lookup(MICSRAArguments * micArgs,
                                    CPTSRACase * alignmentCase, int stepInCase,
                                    unsigned long long * saRanges);
__attribute__((target(mic)))
unsigned long long MICBWTModelSwitchBackward(MICSRAArguments * micArgs,  int i, uint8_t errorInserted,
                                    CPTSRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occError, uint8_t nbMismatch, int occQuality);

__attribute__((target(mic)))
unsigned long long MICBWTModelSwitchAnyDirection(MICSRAArguments * micArgs,  int i, uint8_t errorInserted,
                                    CPTSRACase * alignmentCase, int stepInCase, 
                                    unsigned long long * saRanges,
                                    uint8_t occError, uint8_t nbMismatch, int occQuality);

#endif

