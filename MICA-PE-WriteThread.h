//
//    MICA-PE-WriteThread.h
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

#ifndef __MICA_PE_WRITETHREAD_H__
#define __MICA_PE_WRITETHREAD_H__

///////////////////////////////////////
// Including Standard Libraries
///////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#include <pthread.h>
#include <unistd.h>
#include <errno.h>

///////////////////////////////////////
// Including SOAP2-DP Libraries
///////////////////////////////////////
#include "2bwt-flex/ListConsumer.h"
#include "2bwt-flex/2BWT-PEAlgnmt.h"
#include "2bwt-flex/2bwt-lib/Timing.h"
#include "2bwt-flex/utilities/Logging.h"
#include "2bwt-flex/utilities/SimpleQueue.h"

///////////////////////////////////////
// Including MIC Libraries
///////////////////////////////////////
#include "MIC-SRA2BWTMdl.h"
#include "MIC-SRAAlgnmt.h"
#include "MIC-PEAlgnmt.h"
#include "MIC-DPAlgnmt.h"
#include "MemMan.h"
#include "MicMappingQuality.h"
#include "MICA-PE-ReadThread.h"

///////////////////////////////////////
// DEBUG FLAGS
///////////////////////////////////////
//Define the below parameter to output the write thread
//audit events into standard output
//#define DEBUG_MICA_PEWT_LOGGING

//Define the below parameter to output the write thread
//read handling event into standard output
//#define DEBUG_MICA_PEWT_READ_HANDLING

#define MICA_PE_WRITE_THREAD_MAX_BATCH          3

#define MICA_PE_WRITE_THREAD_HEALTH_OK          0
#define MICA_PE_WRITE_THREAD_HEALTH_DEPLETED    1

#define MICA_PE_WRITE_THREAD_PAYLOAD_OUTPUT_TYPE_UNKNOWN     0
#define MICA_PE_WRITE_THREAD_PAYLOAD_OUTPUT_TYPE_MIC_MIX     1
#define MICA_PE_WRITE_THREAD_PAYLOAD_OUTPUT_TYPE_SKIP_MIC    2

#define MICA_PE_WRITE_HEALTH_STRING_MAX_LENGTH               20

typedef struct PEWTArgs {
    SRAArguments  * sraArgTemplate;
    PEArguments   * peArgTemplate;
    unsigned char * charMap;
    unsigned int maxNumReadPerCore;
    int maxNumOfMICThreads;
    int InputAlignmentModel;
} PEWTArgs;

typedef struct PEWTPerfStats {
    unsigned long long readNumHandled;
    unsigned long long readNumByHandlers[MIC_PE_OUTPUT_STATUS_COUNT];
    double             processTime;
} PEWTPerfStats;

typedef struct PEATPerfStats {
    unsigned long long readNumHandled;
    unsigned long long readNumAligned;
    double             processTime;
} PEATPerfStats;

typedef struct PEAlgnmtThreadArg {
    void *        writeThread;
    int           threadIdx;
    PEATPerfStats algnmtThreadStats;
} PEAlgnmtThreadArg;

typedef struct PEWriteThread {

    // Global MICA-PE Argument Configuration
    PEWTArgs pewtArguments;
    
    // Statistics
    unsigned long long consPairCount;
    unsigned long long consPairRead;
    unsigned long long cpuConsSaCount;
    unsigned long long cpuConsReadCount;
    unsigned long long consReadCount;
    double totalOutputTime;
    
    // POSIX Thread
    int numAlignmentThread;
    pthread_t       outputThreadBody;
    pthread_t     * alignmentThreadBody;
    PEAlgnmtThreadArg * algnmtThreadArg;
    
    // POSIX Mutex and Condition Variables
    pthread_mutex_t mutexThreadHealth;
    pthread_mutex_t mutexFreeQueue;
    pthread_cond_t  condFreeQueue;
    pthread_mutex_t mutexAlignmentQueue;
    pthread_cond_t  condAlignmentQueue;
    pthread_mutex_t mutexOutputQueue;
    pthread_cond_t  condOutputQueue;
    
    ////////////////////////////////////////
    // Thread Payload - Configuration
    ////////////////////////////////////////
    unsigned int          queryInBatch[MICA_PE_WRITE_THREAD_MAX_BATCH];
    unsigned int          batchSize[MICA_PE_WRITE_THREAD_MAX_BATCH];
    unsigned int          numThreads[MICA_PE_WRITE_THREAD_MAX_BATCH];
    unsigned int          firstReadIdx[MICA_PE_WRITE_THREAD_MAX_BATCH];
    int                   algnFuncReasonFlag[MICA_PE_WRITE_THREAD_MAX_BATCH];
    
    ////////////////////////////////////////
    // Thread Payload - Input
    ////////////////////////////////////////
    char                * readNameBufferFrame[MICA_PE_WRITE_THREAD_MAX_BATCH];
    unsigned char       * readBodyBufferFrame[MICA_PE_WRITE_THREAD_MAX_BATCH];
    uint16_t            * readLengthBufferFrame[MICA_PE_WRITE_THREAD_MAX_BATCH];
    char                * readQualityBufferFrame[MICA_PE_WRITE_THREAD_MAX_BATCH];
    char                * mateNameBufferFrame[MICA_PE_WRITE_THREAD_MAX_BATCH];
    unsigned char       * mateBodyBufferFrame[MICA_PE_WRITE_THREAD_MAX_BATCH];
    uint16_t            * mateLengthBufferFrame[MICA_PE_WRITE_THREAD_MAX_BATCH];
    char                * mateQualityBufferFrame[MICA_PE_WRITE_THREAD_MAX_BATCH];
    PEReadOutput        * peBatchOutput[MICA_PE_WRITE_THREAD_MAX_BATCH];
    PEInputParam        * peInputParam[MICA_PE_WRITE_THREAD_MAX_BATCH];

    
    ////////////////////////////////////////
    // Thread Payload - Output
    ////////////////////////////////////////
    uint8_t             outputType[MICA_PE_WRITE_THREAD_MAX_BATCH];
    uint16_t            * occCount[MICA_PE_WRITE_THREAD_MAX_BATCH];
    uint8_t             * outputBufferStatus[MICA_PE_WRITE_THREAD_MAX_BATCH];
    unsigned int        * outputBufferFrame[MICA_PE_WRITE_THREAD_MAX_BATCH];
    MICSRAOccMetadata   * outputBufferMeta[MICA_PE_WRITE_THREAD_MAX_BATCH];
    MICDPOccurrence     * dpOutputBufferFrame[MICA_PE_WRITE_THREAD_MAX_BATCH];
    PEMappingQuality    * mappingQualities[MICA_PE_WRITE_THREAD_MAX_BATCH];
    SRAOCCCollector     * cpuSraOccCollector[MICA_PE_WRITE_THREAD_MAX_BATCH];
    DPOCCCollector      * cpuDpOccCollector[MICA_PE_WRITE_THREAD_MAX_BATCH];
    
    // Ring Buffer Information & Thread Health
    int bufferSize;
    SimpleQueue     * freeQueue;
    SimpleQueue     * alignmentQueue;
    SimpleQueue     * outputQueue;
    uint8_t         threadHealth;
    int             parentControllerIdx;
    
    // Health string and stats
    char             healthString[MICA_PE_WRITE_HEALTH_STRING_MAX_LENGTH];
    PEWTPerfStats    performanceStats;
    
} PEWriteThread;

PEWriteThread * PEWTCreate(unsigned int maxNumReadPerCore, int maxNumOfMICThreads, int numWriteThreadBuffer, int numAlignmentThread, int parentControllerIdx);

void PEWTRegisterArguments(PEWriteThread * writeThread, 
                        SRAArguments * sraArgs, PEArguments * peArgs, 
                        int InputAlignmentModel,
                        unsigned char * charMap);


void PEWTCloneBuffer(PEWriteThread * writeThread, 
                        unsigned int numThreads, unsigned int batchSize, unsigned int firstReadIdx, unsigned int queryInBatch,
                        char * readNameBufferFrame, unsigned char * readBodyBufferFrame, uint16_t * readLengthBufferFrame,
                        char * mateNameBufferFrame, unsigned char * mateBodyBufferFrame, uint16_t * mateLengthBufferFrame,
                        char * readQualityBufferFrame, char * mateQualityBufferFrame,
                        PEInputParam * peInputParam, PEReadOutput * peBatchOutput,
                        uint8_t payloadOutputType,
                        uint16_t * occCount, uint8_t * outputBufferStatus, unsigned int * outputBufferFrame,
                        MICSRAOccMetadata * outputBufferMeta, MICDPOccurrence * dpOutputBufferFrame, PEMappingQuality * mappingQualities,
                        int algnFuncReasonFlag);

void PEWTRoll(PEWriteThread * writeThread);
void PEWTSignal(PEWriteThread * writeThread, int threadHealth);
void PEWTThreadJoin(PEWriteThread * writeThread);
void PEWTPollHealthString(PEWriteThread * writeThread);
char * PEWTGetHealthText(PEWriteThread * writeThread);
void PEWTFree(PEWriteThread * writeThread);

#endif

