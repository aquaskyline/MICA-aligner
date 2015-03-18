//
//    CPUControllerThread.c
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

#include "CPUControllerThread.h"

void CPUCTWriteStatus(CPUControllerThread * thread, char *status) {
    strcpy(thread->healthString,status);
}

void CPUCTSetSRA(CPUControllerThread * thread, int alignmentModel, SRAArguments * cpuSraArgTemplate) {
    thread->alignmentModel = alignmentModel;
    thread->cpuSraArgTemplate = cpuSraArgTemplate;
}

void CPUCTSetPE(CPUControllerThread * thread, PEArguments * cpuPeArgTemplate) {
    thread->cpuPeArgTemplate = cpuPeArgTemplate;
}

void CPUCTSetReadThread(CPUControllerThread * thread, PEReadThread * peReadThread) {
    thread->peReadThread = peReadThread;
}

void CPUCTSetWriteThread(CPUControllerThread * thread, PEWriteThread * peWriteThread) {
    thread->writeThread = peWriteThread;
}

CPUControllerThread * CPUCTCreate(uint8_t threadId, char * charMap,
    unsigned int cpuMaxBatchSize) {

    CPUControllerThread * ret = (CPUControllerThread *) MEMManMalloc(sizeof(CPUControllerThread), MEMORY_TYPE_CPU);
    
    // Initialize some fields in CPUControllerThread
    ret->threadId = threadId;
    ret->charMap = charMap;
    ret->cpuMaxBatchSize = cpuMaxBatchSize;
    ret->t = MEMManMalloc(sizeof(pthread_t), MEMORY_TYPE_CPU);
    ret->writeThread = NULL;
    
    ret->statsTotalAlignmentTime = 0;
    ret->statsTotalReadLoadTime = 0;
    ret->statsTotalModelBuildTime = 0;
    ret->statsTotalWaitTime = 0;
    ret->statsQueryHandled = 0;
    
    #ifdef MICA_PE_ENABLE_HEALTH_CHECK
        strcpy(ret->healthString,"INIT");
    #else
        strcpy(ret->healthString,"DISABLED");
    #endif
    MICCTPollHealthString(ret);
    
    return ret;
}

static void * start(void * cpuControllerThread) {
    CPUControllerThread * thread = (CPUControllerThread *) cpuControllerThread;
    
    uint8_t threadId = thread->threadId;
    PEReadThread * peReadThread = thread->peReadThread;

    // Timestamp'ed -----------
    double startTime = setStartTime();
    double lastEventTime = 0;
    // double timestamp = getElapsedTime(startTime);
    // printf("CPUC#%d - Elapsed time (Index Loading) : %9.4f seconds\n\n", timestamp - lastEventTime);
    // lastEventTime = timestamp;
    // ------------------------
    
    ////////////////////////////////////////////////////////////
    // CPU Allocation of Memory Frame
    // Fake pointer for memory straight copy.
    ////////////////////////////////////////////////////////////
    unsigned int numReadPerCore = thread->cpuMaxBatchSize;
    
    char * readNameBufferFrame;
    unsigned char * readBodyBufferFrame;
    uint16_t * readLengthBufferFrame;
    char * mateNameBufferFrame;
    unsigned char * mateBodyBufferFrame;
    uint16_t * mateLengthBufferFrame;
    
    PEInputParam * peInputParamBufferFrame;
    PEReadOutput * peReadOutputBufferFrame;
    char * readQualityBufferFrame;
    char * mateQualityBufferFrame;
    int algnFuncReasonFlag;

    // Timestamp'ed -----------
    double timestamp = getElapsedTime(startTime);
    printf("CPUC#%d - Elapsed time (CPU Memory Allocation) : %9.4f seconds\n\n", thread->threadId, timestamp - lastEventTime);
    lastEventTime = timestamp;
    // ------------------------
    
    // Write Thread
    PEWriteThread * writeThread = thread->writeThread;

    double processingTime = 0;
    double alignmentTime = 0;
    double readLoadTime = 0;
    double modelTime = 0;
    double waitTime = 0;

    unsigned int queryInBatch = 1;
    unsigned int readInBatch = 1;
    unsigned int mateInBatch = 1;
    unsigned long long queryHandled = 0;

    #ifdef MIC_DEBUG_PRINT_OUTPUT_COUNTS
    unsigned int numOfClosed = 0;
    unsigned int numOfSaReq = 0;
    unsigned int numOfSaContent = 0;
    unsigned int numOfOccContent = 0;
    #endif    

    char batchHandled = 0;
    while ( queryInBatch != 0 && ( batchHandled < CPU_CONTROL_MAX_BATCH || CPU_CONTROL_MAX_BATCH == -1 ) ) {

        #ifdef MICA_PE_ENABLE_HEALTH_CHECK
            // Update status to WAIT
            CPUCTWriteStatus(thread,"WAIT");
        #endif
        
        // Read in the query for the buckets for all threads
        // The current stamp is being retrieved
        PEReadBatch * peReadBatch = PEReadThread_wait(peReadThread);

        #ifdef MICA_PE_ENABLE_HEALTH_CHECK
            // Update status to WAIT
            CPUCTWriteStatus(thread,"START");
        #endif
        
        // Timestamp'ed -----------
        timestamp = getElapsedTime(startTime);
        waitTime += timestamp - lastEventTime;
        processingTime += timestamp - lastEventTime;
        lastEventTime = timestamp;
        // ------------------------
        
        readInBatch = peReadBatch->readSize;
        mateInBatch = peReadBatch->mateSize;
        
        readNameBufferFrame = peReadBatch->readNames;
        mateNameBufferFrame = peReadBatch->mateNames;
        readBodyBufferFrame = peReadBatch->readPatterns;
        mateBodyBufferFrame = peReadBatch->matePatterns;
        readLengthBufferFrame = peReadBatch->readLengths;
        mateLengthBufferFrame = peReadBatch->mateLengths;
        peInputParamBufferFrame = peReadBatch->batchParam;
        peReadOutputBufferFrame = peReadBatch->batchOutput;
        readQualityBufferFrame = peReadBatch->readQuality;
        mateQualityBufferFrame = peReadBatch->mateQuality;
        algnFuncReasonFlag = peReadBatch->algnFuncReasonFlag;
        queryInBatch = readInBatch > mateInBatch ? mateInBatch : readInBatch;

        if (readInBatch != mateInBatch) {
            printf("CPUC#%d - [ERROR] Number of reads are different between READ(%u) and MATE(%u). Only processing %u..\n", thread->threadId,readInBatch,mateInBatch,queryInBatch);
        }

        if (queryInBatch == 0) {
        
            #ifdef MICA_PE_ENABLE_HEALTH_CHECK
                // Update status to WAIT
                CPUCTWriteStatus(thread,"EXITING");
            #endif
            
            PEReadThread_signal(peReadThread);
            break;
        }
        
        // Timestamp'ed -----------
        timestamp = getElapsedTime(startTime);
        //printf("CPUC#%d - Elapsed time (Read Loading) : %9.4f seconds\n\n", thread->threadId, timestamp - lastEventTime);
        readLoadTime += timestamp - lastEventTime;
        processingTime += timestamp - lastEventTime;
        lastEventTime = timestamp;
        // ------------------------
        
        ////////////////////////////////////////////////////////////
        // Handling Output
        // This following pick results from MIC and output with CPU
        ////////////////////////////////////////////////////////////
    
        #ifdef MICA_PE_ENABLE_HEALTH_CHECK
            // Update status to WAIT
            CPUCTWriteStatus(thread,"ALGNMT");
        #endif
            
        #ifndef MICA_PE_DISCARD_OUTPUT_FROM_MIC
            PEWTCloneBuffer(writeThread,
                            1,queryInBatch,queryHandled,queryInBatch,
                            readNameBufferFrame,readBodyBufferFrame,readLengthBufferFrame,
                            mateNameBufferFrame,mateBodyBufferFrame,mateLengthBufferFrame,
                            readQualityBufferFrame, mateQualityBufferFrame,
                            peInputParamBufferFrame, peReadOutputBufferFrame,
                            MICA_PE_WRITE_THREAD_PAYLOAD_OUTPUT_TYPE_SKIP_MIC,
                            NULL,NULL,NULL,
                            NULL,NULL,NULL,
                            algnFuncReasonFlag);
        #endif

        // We do a direct copy from Read Thread into Write Thread
        // Therefore the read completed signal cannot be fired too early
        PEReadThread_signal(peReadThread);
        
        queryHandled+=queryInBatch;

        // Timestamp'ed -----------
        timestamp = getElapsedTime(startTime);
        //printf("CPUC#%d - Elapsed time (Output) : %9.4f seconds\n\n", thread->threadId, timestamp - lastEventTime);
        processingTime += timestamp - lastEventTime;
        lastEventTime = timestamp;
        // ------------------------

        batchHandled++;
        printf("CPUC#%d - Number of Query Handled = %llu\n", thread->threadId,queryHandled);
        printf("\n");
    }
    
    ////////////////////////////////////////////////////////////
    // Statistics
    ////////////////////////////////////////////////////////////
    thread->statsTotalAlignmentTime = alignmentTime;
    thread->statsTotalReadLoadTime = readLoadTime;
    thread->statsTotalModelBuildTime = modelTime;
    thread->statsTotalWaitTime = waitTime;
    thread->statsQueryHandled = queryHandled;
    
    thread->writeThread = NULL;
    
    return NULL;

}

void CPUCTPrintStats (CPUControllerThread * thread, Logging * logger) {
    
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "CPUC#%d - Number of Query Handled = %llu\n", thread->threadId,thread->statsQueryHandled);
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "CPUC#%d - Read Loading Time       = %9.4f seconds\n", thread->threadId,thread->statsTotalReadLoadTime);
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "CPUC#%d - Model Building Time     = %9.4f seconds\n", thread->threadId,thread->statsTotalModelBuildTime);
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "CPUC#%d - Wait Time               = %9.4f seconds\n", thread->threadId,thread->statsTotalWaitTime);
//    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "CPUC#%d - Alignment Time          = %9.4f seconds\n", thread->threadId,alignmentTime);
//    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "CPUC#%d - Processing (Total) Time = %9.4f seconds\n", thread->threadId,processingTime);
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "\n");
}

void CPUCTStart(CPUControllerThread * thread) {
    pthread_create(thread->t, NULL, start, thread);
}

void CPUCTJoin(CPUControllerThread * thread) {
    pthread_join(*thread->t, NULL);
}

void CPUCTFree(CPUControllerThread * thread) {
    free(thread->t);
    free(thread);
}

void CPUCTPollHealthString(CPUControllerThread * thread) {
    sprintf(thread->consolidHealth,"%s",thread->healthString);
}

char * CPUCTGetHealthText(CPUControllerThread * thread) {
    return thread->consolidHealth;
}
