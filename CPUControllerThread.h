//
//    CPUControllerThread.h
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

#ifndef __CPUCONTROLLERTHREAD_H__
#define __CPUCONTROLLERTHREAD_H__

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

///////////////////////////////////////
// Including MIC Libraries
///////////////////////////////////////
#include "MIC-SRA2BWTMdl.h"
#include "MIC-SRAAlgnmt.h"
#include "MIC-PEAlgnmt.h"
#include "MIC-DPAlgnmt.h"
#include "MICA-PE-ReadThread.h"
#include "MICA-PE-WriteThread.h"
#include "MemMan.h"

#define CPU_CONTROL_WRITE_THREAD_BUFFER_COUNT 1
#define CPU_CONTROL_HEALTH_STRING_MAX_LENGTH 20

// -----------------------------
// Compile-time Configurable #1
// Ticket#70 [MICA] An infrastructure change for the health checking of multi-threading processes
// -----------------------------
// Enable the core health checker to poll
// multiple component and get health status on stdout.
//#define MICA_PE_ENABLE_HEALTH_CHECK
// -----------------------------
// Compile-time Configurable #2
// -----------------------------
// Limit CPU Controller to process only a certain number of batches
#define CPU_CONTROL_MAX_BATCH        1

typedef struct CPUControllerThread {
    uint8_t threadId;

    char * charMap;

    // SRA Settings
    int alignmentModel;
    SRAArguments * cpuSraArgTemplate;

    // PE Settings
    PEArguments * cpuPeArgTemplate;
    
    // CPU Thread Settings
    unsigned int cpuMaxBatchSize;
    int numOfMICThreads;
    
    // Read Thread
    PEReadThread * peReadThread;
    
    // Write Thread
    PEWriteThread * writeThread;
    
    // Thread
    pthread_t * t;
    
    // Health string
    char             healthString[CPU_CONTROL_HEALTH_STRING_MAX_LENGTH];
    char             consolidHealth[CPU_CONTROL_HEALTH_STRING_MAX_LENGTH+MICA_PE_WRITE_HEALTH_STRING_MAX_LENGTH];
    
    // Statistics
    double           statsTotalAlignmentTime;
    double           statsTotalReadLoadTime;
    double           statsTotalModelBuildTime;
    double           statsTotalWaitTime;
    unsigned long long statsQueryHandled;
    
} CPUControllerThread;

CPUControllerThread * CPUCTCreate(uint8_t threadId, char * charMap,
    unsigned int batchSize);
void CPUCTSetSRA(CPUControllerThread * thread, int alignmentModel, SRAArguments * cpuSraArgTemplate);
void CPUCTSetPE(CPUControllerThread * thread, PEArguments * cpuPeArgTemplate);
void CPUCTSetReadThread(CPUControllerThread * thread, PEReadThread * peReadThread);
void CPUCTSetWriteThread(CPUControllerThread * thread, PEWriteThread * peWriteThread);
void CPUCTStart(CPUControllerThread * thread);
void CPUCTJoin(CPUControllerThread * thread);
void CPUCTFree(CPUControllerThread * thread);

void CPUCTPollHealthString(CPUControllerThread * thread);
char * CPUCTGetHealthText(CPUControllerThread * thread);

#endif

