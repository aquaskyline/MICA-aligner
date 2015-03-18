//
//    MICControllerThread.h
//
//    MICA / 2BWT
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



#ifndef __MIC_CONTROLLER_THREAD__
#define __MIC_CONTROLLER_THREAD__

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
// Including 2BWT Libraries
///////////////////////////////////////
#include "2bwt-flex/2bwt-lib/iniparser.h"
#include "2bwt-flex/2bwt-lib/2BWT-Interface.h"
#include "2bwt-flex/2bwt-lib/Timing.h"
#include "2bwt-flex/utilities/Logging.h"

///////////////////////////////////////
// Including SOAP2-DP Libraries
///////////////////////////////////////
#include "2bwt-flex/2BWT-PEAlgnmt.h"
#include "2bwt-flex/2BWT-SRAAlgnmt.h"
#include "2bwt-flex/SRAQueryParser.h"
#include "2bwt-flex/MappingQuality.h"

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
#include "MIC-MEMControl.h"

///////////////////////////////////////
// Including MICCT Read Handler Libraries
///////////////////////////////////////
#include "MICCT-ShortHandler.h"

#define MIC_CONTROL_WRITE_THREAD_BUFFER_COUNT 3
#define MIC_CONTROL_HEALTH_STRING_MAX_LENGTH 20

// -----------------------------
// Compile-time Configurable #1
// Ticket#70 [MICA] An infrastructure change for the health checking of multi-threading processes
// -----------------------------
// Enable the core health checker to poll
// multiple component and get health status on stdout.
//#define MICA_PE_ENABLE_HEALTH_CHECK

// -----------------------------
// Debug Flag
// -----------------------------
// Define to print alignment flow messages on stdout
// Not advisable for data set larger than 1
//#define MIC_DP_DEBUG_PRINT_ALIGNMENT_FLOW

// -----------------------------
// Debug Flag
// -----------------------------
// Define to print input parameter per batch.
//#define MIC_PRINT_INPUT_PARAM

// -----------------------------
// Debug Flag
// -----------------------------
// Define to disable MICController
// by causing it to terminate once initialised.
//#define MICCT_DISABLE_PROCESSING

// -----------------------------
// Debug Flag
// -----------------------------
// Define to enable confirmation on MICController
// flow and pause at significant milestone
//#define MICCT_CONFIRMATION

typedef struct MICControllerThread {
    uint8_t gThreadId;
    uint8_t threadId;

    Idx2BWT * micIdx2BWT;
    char * charMap;

    // SRA Settings
    int alignmentModel;
    SRAArguments * cpuSraArgTemplate;

    // PE Settings
    PEArguments * cpuPeArgTemplate;
    
    // MIC Thread Settings
    unsigned int micMaxBatchSize;
    int numMicThreads;

    // Read Thread
    PEReadThread * peReadThread;
    
    // Write Thread
    PEWriteThread * writeThread;
    
    // Thread
    pthread_t * t;
    
    // Health string
    char             healthString[MIC_CONTROL_HEALTH_STRING_MAX_LENGTH];
    char             consolidHealth[MIC_CONTROL_HEALTH_STRING_MAX_LENGTH+MICA_PE_WRITE_HEALTH_STRING_MAX_LENGTH];
    
    // Statistics
    double           statsTotalAlignmentTime;
    double           statsTotalReadLoadTime;
    double           statsTotalModelBuildTime;
    double           statsTotalWaitTime;
    unsigned long long statsQueryHandled;
    
} MICControllerThread;

MICControllerThread * MICCTCreate(uint8_t gThreadId, Idx2BWT * micIdx2BWT, char * charMap,
    unsigned int batchSize, uint8_t threadId, int numMicThreads);
void MICCTSetSRA(MICControllerThread * thread, int alignmentModel, SRAArguments * cpuSraArgTemplate);
void MICCTSetPE(MICControllerThread * thread, PEArguments * cpuPeArgTemplate);
void MICCTSetReadThread(MICControllerThread * thread, PEReadThread * peReadThread);
void MICCTSetWriteThread(MICControllerThread * thread, PEWriteThread * peWriteThread);
void MICCTStart(MICControllerThread * thread);
void MICCTJoin(MICControllerThread * thread);
void MICCTFree(MICControllerThread * thread);

void MICCTPollHealthString(MICControllerThread * thread);
char * MICCTGetHealthText(MICControllerThread * thread);

#endif
