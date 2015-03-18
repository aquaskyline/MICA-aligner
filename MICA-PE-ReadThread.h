//
//    MICA-PE-ReadThread.h
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


/**
 * This module separate a thread of reading input from file.
 * Currently only support one consumer.
 */

#ifndef __MICA_PE_READ_THREAD_H__
#define __MICA_PE_READ_THREAD_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#include <pthread.h>
#include <unistd.h>
#include <errno.h>

#include "2bwt-flex/ListConsumer.h"
#include "2bwt-flex/SRAQueryParser.h"
#include "2bwt-flex/PEArguments.h"
#include "MemMan.h"

#define PERT_MAX_NUM_OF_BATCH         3
#define PERT_HEALTH_STRING_MAX_LENGTH 20

#define PERT_OUTPUT_STATUS_INITIALISED    0
#define PERT_OUTPUT_STATUS_RUNNING        1
#define PERT_OUTPUT_STATUS_DEPLETED       2
#define PERT_OUTPUT_STATUS_FREE           3

#define MIC_ALGNMT_FOR_UNDEFINED          0
#define MIC_ALGNMT_FOR_OUTPUT             1

//-----------------------------------
// DEBUG FLAGS #1
//-----------------------------------
// Define to print output for flow
//#define MICA_PERT_DEBUG_FLOW
//-----------------------------------
// DEBUG FLAGS #2
//-----------------------------------
// Define to print output statistics
//#define MICA_PERT_DEBUG_INPUT_FILE

// This Input Parameters structures are initialised for each set of
// input file, and maintained by ReadThread.
// ReadThread initialise all the member variables at the beginning.
typedef struct PEInputParam {
    int InputPEInsertionUpperBound;
    int InputPEInsertionLowerBound;
} PEInputParam;

// This Read Output structures are initialised for each set of
// input file, and maintained by ReadThread.
// ReadThread initialise all the member variables at the beginning,
// and continuously update the status and numBatchProduced as 
// it reads the input file. Whenever a write-thread completed a
// batch it will call a callback function PEReadThread_ClosingOutput
// which will then update the numBatchProcessed. If it matches up with
// numBatchProduced and status is DEPLETED, the the output Handler will be
// closed and freed by the callback function.
typedef struct PEReadOutput {
    uint8_t status;
    unsigned int numBatchProduced;
    unsigned int numBatchProcessed;
    SRAOutput * outputHandler;
} PEReadOutput;


// A batch of output from PEReadThread
typedef struct PEReadBatch {
    char * readNames;
    char * mateNames;
    unsigned char * readPatterns;
    unsigned char * matePatterns;
    uint16_t * readLengths;
    uint16_t * mateLengths;
    unsigned readSize;
    unsigned mateSize;
    
    // Store the inputParam here too for WriteThread to use
    PEInputParam * batchParam;
    // Store the outputHandler here too for WriteThread to use
    PEReadOutput * batchOutput;

    char * readQuality;
    char * mateQuality;

    int * readUncertaintyNumber;
    int * mateUncertaintyNumber;

    // Store the flag that indicates
    // the reason for the alignment to be done on this batch
    int algnFuncReasonFlag;
} PEReadBatch;


typedef struct PEReadThread PEReadThread;

PEReadThread * PEReadThread_create(
    unsigned char * charMap,
    unsigned int numOfQuery,
    
    unsigned int maxNameLength,
    unsigned int maxPatternLength,
    
    HSP * hsp,
    
    ListConsumer * pairListConsumer,

    int outputFormat,
    int numSAMThreads
);

void PEReadThread_start(PEReadThread * readThread);

// Wait for the thread to give the caller a new batch
PEReadBatch * PEReadThread_wait(PEReadThread * readThread);

// Signal the thread that the caller has finishing using a batch.
void PEReadThread_signal(PEReadThread * readThread);
void PEReadThread_quit(PEReadThread * readThread);
void PEReadThread_free(PEReadThread * readThread);
void PEReadThread_join(PEReadThread * readThread);

// Health information collector
void PEReadThread_PollHealthString(PEReadThread * readThread);
char * PEReadThread_GetHealthText(PEReadThread * readThread);

// Output Handler
void PEReadThread_ClosingOutput(PEReadOutput * batchOutput);

#endif

