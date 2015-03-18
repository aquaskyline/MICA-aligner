//
//    MICA-PE-WriteThread.c
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

#include "MICA-PE-WriteThread.h"

PEWriteThread * PEWTCreate(unsigned int maxNumReadPerCore, int maxNumOfMICThreads,
        int numWriteThreadBuffer, int numAlignmentThread, int parentControllerIdx) {

    int i,j;
    
    PEWriteThread * writeThread = (PEWriteThread*) MEMManMalloc(sizeof(PEWriteThread),MEMORY_TYPE_CPU);

    unsigned int memsize_readNameBufferFrame    = MAX_SEQ_NAME_LENGTH*maxNumReadPerCore*maxNumOfMICThreads;
    unsigned int memsize_readBodyBufferFrame    = SRA_MAX_READ_LENGTH*maxNumReadPerCore*maxNumOfMICThreads;
    unsigned int memsize_outputBufferFrame      = MIC_SRA_OUTPUT_SIZE_PER_THREAD*maxNumOfMICThreads;
    unsigned int memsize_outputBufferMeta       = MIC_SRA_OUTPUT_META_PER_THREAD*maxNumOfMICThreads;
    unsigned int memsize_bufferPerRead          = maxNumReadPerCore*maxNumOfMICThreads;
    unsigned int memsize_dpOutputBufferFrame    = MIC_DP_OUTPUT_SIZE_PER_THREAD*maxNumOfMICThreads; 
    
    PEWTArgs * pewtArguments = &(writeThread->pewtArguments);
    pewtArguments->sraArgTemplate = NULL;
    pewtArguments->peArgTemplate = NULL;
    pewtArguments->charMap = NULL;
    pewtArguments->maxNumReadPerCore = maxNumReadPerCore;
    pewtArguments->maxNumOfMICThreads = maxNumOfMICThreads;
    
    writeThread->consPairCount = 0;
    writeThread->consPairRead = 0;
    writeThread->cpuConsSaCount = 0;
    writeThread->cpuConsReadCount = 0;
    writeThread->consReadCount = 0;
    
    for (i=0;i<numWriteThreadBuffer;i++) {
        
        writeThread->readNameBufferFrame[i]         = (char*)          MEMManMalloc(sizeof(char)*memsize_readNameBufferFrame,MEMORY_TYPE_CPU);
        writeThread->readBodyBufferFrame[i]         = (unsigned char*) MEMManMalloc(sizeof(unsigned char)*memsize_readBodyBufferFrame,MEMORY_TYPE_CPU);
        writeThread->readQualityBufferFrame[i]      = (char*)          MEMManMalloc(sizeof(char)*memsize_readBodyBufferFrame,MEMORY_TYPE_CPU);
        writeThread->readLengthBufferFrame[i]       = (uint16_t*)      MEMManMalloc(sizeof(uint16_t)*memsize_bufferPerRead,MEMORY_TYPE_CPU);
        writeThread->mateNameBufferFrame[i]         = (char*)          MEMManMalloc(sizeof(char)*memsize_readNameBufferFrame,MEMORY_TYPE_CPU);
        writeThread->mateBodyBufferFrame[i]         = (unsigned char*) MEMManMalloc(sizeof(unsigned char)*memsize_readBodyBufferFrame,MEMORY_TYPE_CPU);
        writeThread->mateQualityBufferFrame[i]      = (char*)          MEMManMalloc(sizeof(char)*memsize_readBodyBufferFrame,MEMORY_TYPE_CPU);
        writeThread->mateLengthBufferFrame[i]       = (uint16_t*)      MEMManMalloc(sizeof(uint16_t)*memsize_bufferPerRead,MEMORY_TYPE_CPU);
        
        writeThread->outputType[i]                  = MICA_PE_WRITE_THREAD_PAYLOAD_OUTPUT_TYPE_UNKNOWN;
        writeThread->occCount[i]                    = (uint16_t*)       MEMManMalloc(sizeof(uint16_t)*memsize_bufferPerRead,MEMORY_TYPE_CPU);
        writeThread->outputBufferStatus[i]          = (uint8_t*)       MEMManMalloc(sizeof(uint8_t)*memsize_bufferPerRead,MEMORY_TYPE_CPU);
        writeThread->outputBufferFrame[i]           = (unsigned int*)  MEMManMalloc(sizeof(unsigned int)*memsize_outputBufferFrame,MEMORY_TYPE_CPU);
        writeThread->outputBufferMeta[i]            = (MICSRAOccMetadata*) MEMManMalloc(sizeof(MICSRAOccMetadata)*memsize_outputBufferMeta,MEMORY_TYPE_CPU);
        writeThread->dpOutputBufferFrame[i]         = (MICDPOccurrence*) MEMManMalloc(sizeof(MICDPOccurrence) * memsize_dpOutputBufferFrame,MEMORY_TYPE_CPU);
        writeThread->mappingQualities[i]            = (PEMappingQuality*) MEMManMalloc(sizeof(PEMappingQuality) * memsize_bufferPerRead,MEMORY_TYPE_CPU);
        
        writeThread->cpuSraOccCollector[i]          = SRAOCCCreate(SRAOCC_FLOOD_TYPE_EXPAND);
        writeThread->cpuDpOccCollector[i]           = DPOCCCreate(SRAOCC_FLOOD_TYPE_EXPAND);
    }
    
    writeThread->algnmtThreadArg = (PEAlgnmtThreadArg*) MEMManMalloc(sizeof(PEAlgnmtThreadArg) * numAlignmentThread,MEMORY_TYPE_CPU);
    writeThread->alignmentThreadBody = (pthread_t*) MEMManMalloc(sizeof(pthread_t) * numAlignmentThread,MEMORY_TYPE_CPU);
    
    writeThread->outputThreadBody = 0;
    for (i=0;i<numAlignmentThread;i++) {
        writeThread->alignmentThreadBody[i] = 0;
        
        writeThread->algnmtThreadArg[i].writeThread = (void*) writeThread;
        writeThread->algnmtThreadArg[i].threadIdx = i;
        
        writeThread->algnmtThreadArg[i].algnmtThreadStats.readNumHandled = 0;
        writeThread->algnmtThreadArg[i].algnmtThreadStats.readNumAligned = 0;
        writeThread->algnmtThreadArg[i].algnmtThreadStats.processTime = 0;
    }
    writeThread->performanceStats.readNumHandled = 0;
    for (j=0;j<MIC_PE_OUTPUT_STATUS_COUNT;j++) {
        writeThread->performanceStats.readNumByHandlers[j] = 0;
    }
    writeThread->performanceStats.processTime = 0;
    writeThread->numAlignmentThread = numAlignmentThread;
    
    // Set up the bufferSize according to the input
    writeThread->bufferSize = numWriteThreadBuffer;
    
    // Initialise the queue for buffer cloning
    writeThread->freeQueue = SQCreate(numWriteThreadBuffer);
    SQInitialiseFullQueue(writeThread->freeQueue);
    
    // Initialise the queue for alignment thread
    writeThread->alignmentQueue = SQCreate(numWriteThreadBuffer);
    SQInitialiseEmptyQueue(writeThread->alignmentQueue);
    
    // Initialise the queue for output thread
    writeThread->outputQueue = SQCreate(numWriteThreadBuffer);
    SQInitialiseEmptyQueue(writeThread->outputQueue);
    
    writeThread->parentControllerIdx = parentControllerIdx;
    
    writeThread->threadHealth = MICA_PE_WRITE_THREAD_HEALTH_OK;
    
    // Initialise the pthread conditional variables
    pthread_mutex_init(&(writeThread->mutexThreadHealth), NULL);
    pthread_mutex_init(&(writeThread->mutexFreeQueue), NULL);
    pthread_cond_init(&(writeThread->condFreeQueue), NULL);
    pthread_mutex_init(&(writeThread->mutexAlignmentQueue), NULL);
    pthread_cond_init(&(writeThread->condAlignmentQueue), NULL);
    pthread_mutex_init(&(writeThread->mutexOutputQueue), NULL);
    pthread_cond_init(&(writeThread->condOutputQueue), NULL);
    
    strcpy(writeThread->healthString,"INIT");
    
    
    return writeThread;
}

void PEWTRegisterArguments(PEWriteThread * writeThread, 
                        SRAArguments * sraArgTemplate, PEArguments * peArgTemplate, 
                        int InputAlignmentModel,
                        unsigned char * charMap) {

    PEWTArgs * pewtArguments = &(writeThread->pewtArguments);

    if ( pewtArguments->sraArgTemplate != NULL || pewtArguments->peArgTemplate != NULL ) {
        printf("[SAFE-GUARD] WriteThread can only register arguments once.\n");

    } else {
        pewtArguments->sraArgTemplate = sraArgTemplate;
        pewtArguments->peArgTemplate  = peArgTemplate;
        
        pewtArguments->charMap = charMap;
        pewtArguments->InputAlignmentModel = InputAlignmentModel;
    }
}

void PEWTCloneBuffer(PEWriteThread * writeThread, 
                        unsigned int numThreads, unsigned int batchSize, unsigned int firstReadIdx, unsigned int queryInBatch,
                        char * readNameBufferFrame, unsigned char * readBodyBufferFrame, uint16_t * readLengthBufferFrame,
                        char * mateNameBufferFrame, unsigned char * mateBodyBufferFrame, uint16_t * mateLengthBufferFrame,
                        char * readQualityBufferFrame, char * mateQualityBufferFrame,
                        PEInputParam * peInputParam, PEReadOutput * peBatchOutput,
                        uint8_t payloadOutputType,
                        uint16_t * occCount, uint8_t * outputBufferStatus, unsigned int * outputBufferFrame,
                        MICSRAOccMetadata * outputBufferMeta, MICDPOccurrence * dpOutputBufferFrame, PEMappingQuality * mappingQualities,
                        int algnFuncReasonFlag) {

    #ifdef DEBUG_MICA_PEWT_LOGGING
        printf ( "PEWT  -  PEWTCloneBuffer   -  Invoked.\n");
    #endif
    
    pthread_mutex_t * mutexFreeQueue = &(writeThread->mutexFreeQueue);
    pthread_cond_t * condFreeQueue = &(writeThread->condFreeQueue);
    
    pthread_mutex_t * mutexAlignmentQueue = &(writeThread->mutexAlignmentQueue);
    pthread_cond_t * condAlignmentQueue = &(writeThread->condAlignmentQueue);
    
    int bufferIdx = 0;
    
    PEWTArgs * pewtArguments = &(writeThread->pewtArguments);
    
    if ( numThreads * batchSize > pewtArguments->maxNumOfMICThreads * pewtArguments->maxNumReadPerCore ) {
        printf("[SAFE-GUARD] Supplied clone subject is larger than the write-thread allocated buffer size.\n");
        printf("%u x %u is larger than %u x %u\n",numThreads,batchSize,pewtArguments->maxNumOfMICThreads,pewtArguments->maxNumReadPerCore);
    }
    
    unsigned int memsize_outputBufferFrame      = MIC_SRA_OUTPUT_SIZE_PER_THREAD*numThreads;
    unsigned int memsize_outputBufferMeta       = MIC_SRA_OUTPUT_META_PER_THREAD*numThreads;
    unsigned int memsize_bufferPerRead          = batchSize*numThreads;
    unsigned int memsize_dpOutputBufferFrame    = MIC_DP_OUTPUT_SIZE_PER_THREAD*numThreads; 
    
    // Block until there is an empty buffer
    pthread_mutex_lock(mutexFreeQueue);
        
        #ifdef DEBUG_MICA_PEWT_LOGGING
            printf ( "PEWT  -  PEWTCloneBuffer   -  Mutex lock (mutexFreeQueue) obtained. Waiting for condition (condFreeQueue)\n");
        #endif
        
        while (SQIsEmpty(writeThread->freeQueue)) {
            printf("[WARNING] Blocking has occurred -- Please revisit pipelining parameters to achieve better performance.\n");
            pthread_cond_wait(condFreeQueue,mutexFreeQueue);
        }

        if (!SQDequeue(writeThread->freeQueue,&bufferIdx)) {
            printf("[SAFE-GUARD] Unexpected mutex error.\n");
            printf("The freeQueue is supposed to be non-empty.\n");
        }
        
        #ifdef DEBUG_MICA_PEWT_LOGGING
            printf ( "PEWT  -  PEWTCloneBuffer   -  Using buffer %d.\n",bufferIdx);
            printf ( "PEWT  -  PEWTCloneBuffer   -  Unlock Mutex(mutexFreeQueue)\n");
        #endif
    
    pthread_mutex_unlock(mutexFreeQueue);

    #ifdef DEBUG_MICA_PEWT_LOGGING
        printf ( "PEWT  -  PEWTCloneBuffer   -  Free Queue Status\n" );
        SQPrintQueue(writeThread->freeQueue);
        printf ( "PEWT  -  PEWTCloneBuffer   -  Copying Payload - Configuration\n");
    #endif
    
    // Clone the payload
    
    ////////////////////////////////////////
    // Thread Payload - Configuration
    ////////////////////////////////////////
    writeThread->queryInBatch[bufferIdx]        = queryInBatch;
    writeThread->batchSize[bufferIdx]           = batchSize;
    writeThread->numThreads[bufferIdx]          = numThreads;
    writeThread->firstReadIdx[bufferIdx]        = firstReadIdx;
    writeThread->peInputParam[bufferIdx]        = peInputParam;
    writeThread->peBatchOutput[bufferIdx]       = peBatchOutput;
    writeThread->algnFuncReasonFlag[bufferIdx]  = algnFuncReasonFlag;
    
    #ifdef DEBUG_MICA_PEWT_LOGGING
        printf ( "PEWT  -  PEWTCloneBuffer   -  Copying Payload - Input\n");
    #endif
    
    ////////////////////////////////////////
    // Thread Payload - Input
    ////////////////////////////////////////
    uint32_t readNameCopySize = MAX_SEQ_NAME_LENGTH * queryInBatch;
    uint32_t readBodyCopySize = SRA_MAX_READ_LENGTH * queryInBatch; 
    uint32_t readLengthCopySize = sizeof(uint16_t) * queryInBatch;
    memcpy(writeThread->readNameBufferFrame[bufferIdx],readNameBufferFrame, readNameCopySize);
    memcpy(writeThread->readBodyBufferFrame[bufferIdx],readBodyBufferFrame, readBodyCopySize);
    memcpy(writeThread->readQualityBufferFrame[bufferIdx], readQualityBufferFrame, readBodyCopySize);
    memcpy(writeThread->readLengthBufferFrame[bufferIdx], readLengthBufferFrame, readLengthCopySize);
    memset(writeThread->readLengthBufferFrame[bufferIdx] + queryInBatch,
        0, sizeof(uint16_t) * memsize_bufferPerRead - readLengthCopySize);
    memcpy(writeThread->mateNameBufferFrame[bufferIdx],mateNameBufferFrame, readNameCopySize);
    memcpy(writeThread->mateBodyBufferFrame[bufferIdx],mateBodyBufferFrame, readBodyCopySize);
    memcpy(writeThread->mateQualityBufferFrame[bufferIdx], mateQualityBufferFrame, readBodyCopySize);
    memcpy(writeThread->mateLengthBufferFrame[bufferIdx], mateLengthBufferFrame, readLengthCopySize);
    memset(writeThread->mateLengthBufferFrame[bufferIdx] + queryInBatch,
        0, sizeof(uint16_t) * memsize_bufferPerRead - readLengthCopySize);

    #ifdef DEBUG_MICA_PEWT_LOGGING
        printf ( "PEWT  -  PEWTCloneBuffer   -  Copying Payload - Output\n");
    #endif
    
    ////////////////////////////////////////
    // Thread Payload - Output
    ////////////////////////////////////////
    writeThread->outputType[bufferIdx]=payloadOutputType;
    if ( payloadOutputType == MICA_PE_WRITE_THREAD_PAYLOAD_OUTPUT_TYPE_MIC_MIX) {
        memcpy(writeThread->occCount[bufferIdx],occCount,sizeof(uint16_t)*memsize_bufferPerRead);
        memcpy(writeThread->outputBufferStatus[bufferIdx],outputBufferStatus,sizeof(uint8_t)*memsize_bufferPerRead);
        memcpy(writeThread->outputBufferFrame[bufferIdx],outputBufferFrame,sizeof(unsigned int)*memsize_outputBufferFrame);
        memcpy(writeThread->outputBufferMeta[bufferIdx],outputBufferMeta,sizeof(MICSRAOccMetadata)*memsize_outputBufferMeta);
        memcpy(writeThread->dpOutputBufferFrame[bufferIdx],dpOutputBufferFrame,sizeof(MICDPOccurrence) * memsize_dpOutputBufferFrame);
        memcpy(writeThread->mappingQualities[bufferIdx],mappingQualities,sizeof(PEMappingQuality) * memsize_bufferPerRead);
    }

    // Block until there is an empty buffer
    pthread_mutex_lock(mutexAlignmentQueue);
        
        #ifdef DEBUG_MICA_PEWT_LOGGING
            printf ( "PEWT  -  PEWTCloneBuffer   -  Mutex lock (mutexAlignmentQueue) obtained. Waiting for condition (condAlignmentQueue)\n");
        #endif
        
        // EW: I am thinking the condition wait is a bit pointless as
        // the free queue was not empty meaning the alignment must have vacancy.
        //while (SQIsFull(writeThread->alignmentQueue)) {
        //    pthread_cond_wait(condAlignmentQueue,mutexAlignmentQueue);
        //}

        if (!SQEnqueue(writeThread->alignmentQueue,bufferIdx)) {
            printf("[SAFE-GUARD] Unexpected mutex error.\n");
            printf("The alignmentQueue is supposed to have vacancy.\n");
        }
        
        #ifdef DEBUG_MICA_PEWT_LOGGING
            printf ( "PEWT  -  PEWTCloneBuffer   -  Signal Listener of condition(condAlignmentQueue)\n");
        #endif
        
        pthread_cond_signal(condAlignmentQueue);
    
    #ifdef DEBUG_MICA_PEWT_LOGGING
        printf ( "PEWT  -  PEWTCloneBuffer   -  Unlock Mutex(mutexAlignmentQueue)\n");
    #endif
    
    pthread_mutex_unlock(mutexAlignmentQueue);
    
    #ifdef DEBUG_MICA_PEWT_LOGGING
        printf ( "PEWT  -  PEWTCloneBuffer   -  AlignmentQueue Status\n" );
        SQPrintQueue(writeThread->alignmentQueue);
    #endif
}

static PEArguments * _PEWTCreatePEArgument(PEWriteThread * writeThread) {
    
    PEWTArgs * pewtArguments = &(writeThread->pewtArguments);

    // Make a mate of the SRA arguments
    SRAArguments * pewtSraReadArgs      = SRAARGMakeMate(pewtArguments->sraArgTemplate);
    SRAArguments * pewtSraReadArgs_neg  = SRAARGMakeSlave(pewtSraReadArgs);
    SRAArguments * pewtSraMateArgs      = SRAARGMakeMate(pewtSraReadArgs);
    SRAArguments * pewtSraMateArgs_neg  = SRAARGMakeSlave(pewtSraMateArgs);
    
    SRAAlgnmtOutputInitNone(pewtSraReadArgs->AlgnmtOutput);
    SRAAlgnmtOutputInitNone(pewtSraMateArgs->AlgnmtOutput);
    
    // Make a mate of the PE arguments
    PEArguments * pewtPeArg = PEARGMakeMate(pewtArguments->peArgTemplate);
    pewtPeArg->PEAlgnmtInput = PEInputMakeClone(pewtPeArg->PEAlgnmtInput);
    //Free-ing the output handler as it will be replaced when a write buffer is being processed
    // by the output handler embedded the read buffer, and passed along.
    SRAAlgnmtOutputFree(pewtPeArg->PEAlgnmtOutput_OccHandle);
    pewtPeArg->PEAlgnmtOutput_OccHandle = NULL;

    // Popluate the PE Argument
    PEARGPopulateSRAARG(pewtPeArg,
                        pewtSraReadArgs,pewtSraReadArgs_neg,
                        pewtSraMateArgs,pewtSraMateArgs_neg);
    
    return pewtPeArg;
}

static void _PEWTFreePEArgument(PEArguments * peArgs) {

    SRAAlgnmtOutputFreeNone(peArgs->sraArgsPos->AlgnmtOutput);
    SRAAlgnmtOutputFreeNone(peArgs->sraArgsPos_mate->AlgnmtOutput);

    SRAARGSlaveFree(peArgs->sraArgsNeg);
    SRAARGSlaveFree(peArgs->sraArgsNeg_mate);
    SRAARGMateFree(peArgs->sraArgsPos);
    SRAARGMateFree(peArgs->sraArgsPos_mate);

    PEInputCloneFree(peArgs->PEAlgnmtInput);
    PEARGMateFree(peArgs);    
}
    
    
void PEWTFree(PEWriteThread * writeThread) {
    int i;
    
    #ifdef DEBUG_MICA_PEWT_LOGGING
        printf ( "PEWT  -  PEWTFree   -  Free-ing the Write Thread\n");
    #endif

    PEWTSignal(writeThread, MICA_PE_WRITE_THREAD_HEALTH_DEPLETED);
    PEWTThreadJoin(writeThread);

    for (i=0;i<writeThread->bufferSize;i++) {
        free(writeThread->readNameBufferFrame[i]);
        free(writeThread->readBodyBufferFrame[i]);
        free(writeThread->readLengthBufferFrame[i]);
        free(writeThread->mateNameBufferFrame[i]);
        free(writeThread->mateBodyBufferFrame[i]);
        free(writeThread->mateLengthBufferFrame[i]);
        free(writeThread->readQualityBufferFrame[i]);
        free(writeThread->mateQualityBufferFrame[i]);
                
        free(writeThread->occCount[i]);
        free(writeThread->outputBufferStatus[i]);
        free(writeThread->outputBufferFrame[i]);
        free(writeThread->outputBufferMeta[i]);
        free(writeThread->dpOutputBufferFrame[i]);
        free(writeThread->mappingQualities[i]);
        
        SRAOCCFree(writeThread->cpuSraOccCollector[i]);
        DPOCCFree(writeThread->cpuDpOccCollector[i]);
    }

    free(writeThread->alignmentThreadBody);
    free(writeThread->algnmtThreadArg);
    
    SQFree(writeThread->freeQueue);
    SQFree(writeThread->alignmentQueue);
    SQFree(writeThread->outputQueue);
    
    free(writeThread);
}

/////////////////////////////////////////////////////////////////////////
//
// ALIGNMENT THREAD BODY
//
/////////////////////////////////////////////////////////////////////////

void * _PEWTAlignmentRollBody(void * _threadArg) {

    PEAlgnmtThreadArg * threadArg = (PEAlgnmtThreadArg*) _threadArg;

    PEWriteThread * writeThread = (PEWriteThread*) threadArg->writeThread;
    int threadIdx = threadArg->threadIdx;
    
    pthread_mutex_t * mutexFreeQueue = &(writeThread->mutexFreeQueue);
    pthread_cond_t * condFreeQueue = &(writeThread->condFreeQueue);
    
    pthread_mutex_t * mutexAlignmentQueue = &(writeThread->mutexAlignmentQueue);
    pthread_cond_t * condAlignmentQueue = &(writeThread->condAlignmentQueue);
    
    pthread_mutex_t * mutexOutputQueue = &(writeThread->mutexOutputQueue);
    pthread_cond_t * condOutputQueue = &(writeThread->condOutputQueue);
    
    PEWTArgs * pewtArguments = &(writeThread->pewtArguments);
    
    PEArguments * peArg = _PEWTCreatePEArgument(writeThread);
    
    // Set up the temporary output collector for the intermediate results
    peArg->PEAlgnmtOutput_OccHandle = SRAAlgnmtOutputConstruct();
    SRAAlgnmtOutputInitSAMStore(peArg->PEAlgnmtOutput_OccHandle);
    
    unsigned char * charMap = pewtArguments->charMap;
    
    PEDPSetting * pedpSetting = peArg->pedpSetting;
    
    // Loosen the matching criteria
    SRASetting * seedSetting = SRASettingMakeClone(peArg->sraArgsPos->AlgnmtSetting);
    seedSetting->OutputType = SRA_REPORT_ALL_BEST;
    if (pedpSetting->SGASeedLooseCriteria>=0) {
        seedSetting->MaxError = pedpSetting->SGASeedLooseCriteria;
    } else {
        seedSetting->MaxError += pedpSetting->SGASeedLooseCriteria;
    }
    
    // oneMoreMismatchSetting increase mismatch by 1
    SRASetting * oneMoreMismatchSetting = SRASettingMakeClone(peArg->sraArgsPos->AlgnmtSetting);
    oneMoreMismatchSetting->OutputType == SRA_REPORT_ALL;
    oneMoreMismatchSetting->MaxError += 1;
    
    // Construct SRA Model from SRA Settings
    SRAModelSet * sraModelSet = SRAModelSetConstruct(
                                    peArg->sraArgsPos->AlgnmtSetting,
                                    peArg->sraArgsPos->AlgnmtIndex,
                                    SRA_MODEL_16G,
                                    SRA_MIN_READ_LENGTH,
                                    SRA_MAX_READ_LENGTH);
    SRAModelSet * sraModelSet_seed = SRAModelSetConstruct(
                                    seedSetting,
                                    peArg->sraArgsPos->AlgnmtIndex,
                                    SRA_MODEL_16G,
                                    SRA_MIN_READ_LENGTH,
                                    SRA_MAX_READ_LENGTH);
    SRAModelSet * sraModelSet_extend = SRAModelSetConstruct(
                                    oneMoreMismatchSetting,
                                    peArg->sraArgsPos->AlgnmtIndex,
                                    SRA_MODEL_16G,
                                    SRA_MIN_READ_LENGTH,
                                    SRA_MAX_READ_LENGTH);
    
    unsigned long long i, j, k;
    int invalidStatus = FALSE;
    //int dummyQuality[SRA_MAX_READ_LENGTH];
    //memset(dummyQuality,0,sizeof(int)*SRA_MAX_READ_LENGTH);
    
    #ifdef DEBUG_MICA_PEWT_LOGGING
        printf ( "PEWT  -  _PEWTAlignmentRollBody  -  Invoked and Entering Main-Loop\n");
    #endif
    
    // Timestamp'ed -----------
    double startTime = setStartTime();
    double lastEventTime = 0;
    double totalElipsedTime = 0;
    // ------------------------

    while (1) {

        int bufferIdx;
        
        // Block until there is an non-empty buffer
        pthread_mutex_lock(mutexAlignmentQueue);
        
            // Loop wait as long as
            // - alignment queue is empty
            // - AND free queue is not full OR thread is healthy
            while (
                SQIsEmpty(writeThread->alignmentQueue) &&
                (
                    !SQIsFull(writeThread->freeQueue) ||
                    writeThread->threadHealth==MICA_PE_WRITE_THREAD_HEALTH_OK
                )
            ) {
                pthread_cond_wait(condAlignmentQueue,mutexAlignmentQueue);
            }
            // Break in case
            // - free queue is full
            // - AND alignment queue is empty
            // - AND thread is not healthy
            if (
                SQIsFull(writeThread->freeQueue) &&
                SQIsEmpty(writeThread->alignmentQueue) &&
                writeThread->threadHealth!=MICA_PE_WRITE_THREAD_HEALTH_OK
            ) {

                #ifdef DEBUG_MICA_PEWT_LOGGING
                    printf ( "PEWT  -  _PEWTAlignmentRollBody  -  Thread is not healthy and Buffer is empty. Exiting.\n");
                #endif
                
                pthread_cond_signal(condAlignmentQueue);
                pthread_mutex_unlock(mutexAlignmentQueue);
                break;
            }
            
            if (!SQDequeue(writeThread->alignmentQueue,&bufferIdx)) {
                printf("[SAFE-GUARD] Unexpected mutex error.\n");
                printf("The alignmentQueue is supposed to be non-empty.\n");
            }
        pthread_mutex_unlock(mutexAlignmentQueue);

        // Timestamp'ed -----------
        double timestamp = getElapsedTime(startTime);
        lastEventTime = timestamp;
        // ------------------------
        
        #ifdef DEBUG_MICA_PEWT_LOGGING
            printf ( "PEWT  -  _PEWTAlignmentRollBody  -  Received Signal on condition(condAlignmentQueue)\n");
        #endif
        
        
        //Process Payload
        //
        #ifdef DEBUG_MICA_PEWT_LOGGING
            printf ( "PEWT  -  _PEWTAlignmentRollBody  -  Processing Payload from buffer %d.\n", bufferIdx);
        #endif
        
        unsigned int          queryInBatch          = writeThread->queryInBatch[bufferIdx];
        unsigned int          batchSize             = writeThread->batchSize[bufferIdx];
        unsigned int          numThreads            = writeThread->numThreads[bufferIdx];
        unsigned int          firstReadIdx          = writeThread->firstReadIdx[bufferIdx];
        int                   algnFuncReasonFlag    = writeThread->algnFuncReasonFlag[bufferIdx];
        
        #ifdef DEBUG_MICA_PEWT_LOGGING
            printf ( "PEWT  -  _PEWTAlignmentRollBody  -  #Thread %u / BatchSize %u / #Query %u.\n",numThreads,batchSize,queryInBatch);
        #endif
        
        char                * readNameBufferFrame   = writeThread->readNameBufferFrame[bufferIdx];
        unsigned char       * readBodyBufferFrame   = writeThread->readBodyBufferFrame[bufferIdx];
        char                * readQualityBufferFrame= writeThread->readQualityBufferFrame[bufferIdx];
        uint16_t            * readLengthBufferFrame = writeThread->readLengthBufferFrame[bufferIdx];
        char                * mateNameBufferFrame   = writeThread->mateNameBufferFrame[bufferIdx];
        unsigned char       * mateBodyBufferFrame   = writeThread->mateBodyBufferFrame[bufferIdx];
        char                * mateQualityBufferFrame= writeThread->mateQualityBufferFrame[bufferIdx];
        uint16_t            * mateLengthBufferFrame = writeThread->mateLengthBufferFrame[bufferIdx];
        
        PEInputParam        * peInputParam          = writeThread->peInputParam[bufferIdx];
        PEReadOutput        * peBatchOutput         = writeThread->peBatchOutput[bufferIdx];
            
        uint8_t               outputType            = writeThread->outputType[bufferIdx];
        uint16_t            * occCount              = writeThread->occCount[bufferIdx];
        uint8_t             * outputBufferStatus    = writeThread->outputBufferStatus[bufferIdx];
        
        // Set up the temporary output collector for the intermediate results
        SRAOCCCollector     * cpuSraOccCollector    = writeThread->cpuSraOccCollector[bufferIdx];
        DPOCCCollector      * cpuDpOccCollector     = writeThread->cpuDpOccCollector[bufferIdx];
        
        SRAArguments * sraReadArgs = peArg->sraArgsPos;
        SRAArguments * sraReadArgs_neg = peArg->sraArgsNeg;
        SRAArguments * sraMateArgs = peArg->sraArgsPos_mate;
        SRAArguments * sraMateArgs_neg = peArg->sraArgsNeg_mate;
        
        unsigned char readNegPattern[SRA_MAX_READ_LENGTH], mateNegPattern[SRA_MAX_READ_LENGTH];
        unsigned char readNegQuality[SRA_MAX_READ_LENGTH], mateNegQuality[SRA_MAX_READ_LENGTH];
        // Populate input parameter and output handler for batch
        peArg->PEAlgnmtInput->insertUbound = peInputParam->InputPEInsertionUpperBound;
        peArg->PEAlgnmtInput->insertLbound = peInputParam->InputPEInsertionLowerBound;
        
        for (i=0;i<numThreads;i++) {
            unsigned int offset = i * batchSize;
                
            // Build up the SRA Model in case there are new read lengths
            for (j=0;j<batchSize;j++) {
                unsigned int readIdxInBatch = offset + j;
                // buildSRAModels with each read and mate into
                // sraModel, seedModel and oneMoreMismatchModel
                int readLength = readLengthBufferFrame[readIdxInBatch];
                int mateLength = mateLengthBufferFrame[readIdxInBatch];
                // If seedLength > 1, treat input as constant seedLength or seedOverlap
                // If seedLength <= 1, treat input as a ratio to readLen
                // + 0.5 to do rounding, set Length to zero when readLength or mateLength = 0
                int readSeedLength = ( pedpSetting->SGASeedLength > 1.0 && readLength != 0 ? pedpSetting->SGASeedLength
                    : (double) readLength * pedpSetting->SGASeedLength ) + 0.5;
                int mateSeedLength = ( pedpSetting->SGASeedLength > 1.0 && mateLength != 0 ? pedpSetting->SGASeedLength
                    : (double) mateLength * pedpSetting->SGASeedLength ) + 0.5;
                int readSeedLength2 = ( pedpSetting->SGASeedLength_1 > 1.0 && readLength != 0 ? pedpSetting->SGASeedLength_1
                    : (double) readLength * pedpSetting->SGASeedLength_1 ) + 0.5;
                int mateSeedLength2 = ( pedpSetting->SGASeedLength_1 > 1.0 && mateLength != 0 ? pedpSetting->SGASeedLength_1
                    :(double) mateLength * pedpSetting->SGASeedLength_1 ) + 0.5;
                SRAModelConstruct(sraModelSet, readLength);
                SRAModelConstruct(sraModelSet, mateLength);
                SRAModelConstruct(sraModelSet_seed, readSeedLength);
                SRAModelConstruct(sraModelSet_seed, mateSeedLength);
                SRAModelConstruct(sraModelSet_seed, readSeedLength2);
                SRAModelConstruct(sraModelSet_seed, mateSeedLength2);
                SRAModelConstruct(sraModelSet_extend, readLength);
                SRAModelConstruct(sraModelSet_extend, mateLength);
            }
            
            for (j=0;j<batchSize;j++) {
                unsigned int readIdxInBatch = offset + j;

                // Buffer could be partially filled.
                // Early termination w.r.t to the size of the particular buffer
                if (readIdxInBatch >= queryInBatch) {
                    break;
                }
                unsigned int readIdx = firstReadIdx + readIdxInBatch;
            
                #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                    printf("[INFO] Handling Read %u\n",readIdx);
                #endif
                
                unsigned char * readPatternPtr = &readBodyBufferFrame[readIdxInBatch*SRA_MAX_READ_LENGTH];
                char * readPatternName         = &readNameBufferFrame[readIdxInBatch*MAX_SEQ_NAME_LENGTH];
                char * readQualityPtr          = &readQualityBufferFrame[readIdxInBatch*SRA_MAX_READ_LENGTH];
                uint16_t readPatternLen        = readLengthBufferFrame[readIdxInBatch];
                SRAFlipPattern(charMap,readPatternPtr,readPatternLen,readNegPattern);
                SRAFlipQuality(charMap,readQualityPtr,readPatternLen,readNegQuality);

                unsigned char * matePatternPtr = &mateBodyBufferFrame[readIdxInBatch*SRA_MAX_READ_LENGTH];
                char * matePatternName         = &mateNameBufferFrame[readIdxInBatch*MAX_SEQ_NAME_LENGTH];
                char * mateQualityPtr          = &mateQualityBufferFrame[readIdxInBatch*SRA_MAX_READ_LENGTH]; ;
                uint16_t matePatternLen        = mateLengthBufferFrame[readIdxInBatch];
                SRAFlipPattern(charMap,matePatternPtr,matePatternLen,mateNegPattern);
                SRAFlipQuality(charMap,mateQualityPtr,readPatternLen,mateNegQuality);
                
                // Populate the SRA information
                SRAQueryInfoPopulate(sraReadArgs->QueryInfo,readIdx,readPatternName,readPatternLen,readPatternPtr,readPatternPtr,QUERY_POS_STRAND,readQualityPtr);
                SRAQueryInfoPopulate(sraReadArgs_neg->QueryInfo,readIdx,readPatternName,readPatternLen,readNegPattern,readPatternPtr,QUERY_NEG_STRAND,readNegQuality);

                SRAQueryInfoPopulate(sraMateArgs->QueryInfo,readIdx,matePatternName,matePatternLen,matePatternPtr,matePatternPtr,QUERY_POS_STRAND,mateQualityPtr);
                SRAQueryInfoPopulate(sraMateArgs_neg->QueryInfo,readIdx,matePatternName,matePatternLen,mateNegPattern,matePatternPtr,QUERY_NEG_STRAND,mateNegQuality);
                
                // Initialise the result counter and MAPQ buffers
                SRAResultCountInitialise(sraReadArgs->AlgnmtStats);
                SRAResultCountInitialise(sraMateArgs->AlgnmtStats);
                MAPQCalculatorInitialise(sraReadArgs->MapqCalc);
                MAPQCalculatorInitialise(sraMateArgs->MapqCalc);
                MAPQCalculatorInitialise(peArg->MapqCalc);
                
                if ( outputType == MICA_PE_WRITE_THREAD_PAYLOAD_OUTPUT_TYPE_UNKNOWN) {
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Payload Error - Unknown - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    
                    printf("[INFO] Payload Error - Unknown - %s / %s!\n",readPatternName,matePatternName);
                    // In this case the read has never been initialised.
                    // therefore the read always requires CPU processing.

                    unsigned int runPairOcc = PEProcessReadDoubleStrand(peArg,
                                                                        sraModelSet,sraModelSet_seed,sraModelSet_extend);
                    // As the output buffer are invalid
                    occCount[readIdxInBatch] = 0;
                    threadArg->algnmtThreadStats.readNumAligned++;
                } else if ( outputType == MICA_PE_WRITE_THREAD_PAYLOAD_OUTPUT_TYPE_SKIP_MIC) {
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Payload unaligned by MIC - MIC unaligned - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    // In this case the read has never been passed into MIC
                    // therefore the read always requires CPU processing.
                    
                    unsigned int runPairOcc = PEProcessReadDoubleStrand(peArg,
                                                                        sraModelSet,sraModelSet_seed,sraModelSet_extend);
                    // As the output buffer are invalid
                    occCount[readIdxInBatch] = 0;
                    threadArg->algnmtThreadStats.readNumAligned++;
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_CLOSED) {
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Unaligned closed read - MIC unaligned - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    // In this case the output buffer is flooded as the number of output
                    // generated is larger than the MIC_SRA_OUTPUT_MAX_ALIGNMENT; or the number of header
                    // generated is larger than MIC_SRA_OUTPUT_MAX_META. MIC was not
                    // able to handle all the alignment result hence declared the read 
                    // unaligned on MIC.
                    //cpuConsSaCount += ProcessReadDoubleStrand(sraReadArgs,sraReadArgs_neg,sraModelSet,sraModelSet);
                    unsigned int runPairOcc = PEProcessReadDoubleStrand(peArg,
                                                                        sraModelSet,sraModelSet_seed,sraModelSet_extend);

                    threadArg->algnmtThreadStats.readNumAligned++;
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_UNHANDLE) {
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Unhandled read - MIC unaligned - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    // In this case the entire core shared output buffer is flooded as the number of output
                    // generated is larger than the MIC_SRA_OUTPUT_SIZE_PER_THREAD. MIC was not
                    // able to handle all the reads in the input read block hence declared all consecutive read 
                    // unaligned on MIC.
                    //cpuConsSaCount += ProcessReadDoubleStrand(sraReadArgs,sraReadArgs_neg,sraModelSet,sraModelSet);
                    unsigned int runPairOcc = PEProcessReadDoubleStrand(peArg,
                                                                        sraModelSet,sraModelSet_seed,sraModelSet_extend);

                    threadArg->algnmtThreadStats.readNumAligned++;
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_SKIPPED) {
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_OUTPUT_STATUS_OPEN) {
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_PAIR) {
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_BAD_PAIR) {
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_READ_NO_ALIGNMENT) {
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_MATE_NO_ALIGNMENT) {
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_BOTH_NO_ALIGNMENT) {
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_DP_BASE_MATE) {
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_DP_BASE_READ) {
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_DP_MIX_BASE) {
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_DP_SEED) {
                } else {
                }
                
                //Copy result into intermediate buffer from read buffer
                SRAOCCCollectorPointer * occptr = SRAOCCPTCreate(peArg->PEAlgnmtOutput_OccHandle->occCollector);
                SRAOccurrence * sraOcc = SRAOCCPTRead(occptr);
                while ( sraOcc != NULL ) {
                    SRAOCCAddOccurrenceToBuckets(cpuSraOccCollector,sraOcc);
                    
                    SRAOCCPTNext(occptr);
                    sraOcc = SRAOCCPTRead(occptr);
                }
                SRAOCCPTFree(occptr);

                DPOCCCollectorPointer * dpoccptr = DPOCCPTCreate(peArg->dpArguments->dpOccCollector);
                DPOccurrence * dpOcc = DPOCCPTRead(dpoccptr);
                while ( dpOcc != NULL ) {
                    DPOCCAddOccurrenceToBuckets(cpuDpOccCollector,dpOcc);
                    
                    DPOCCPTNext(dpoccptr);
                    dpOcc = DPOCCPTRead(dpoccptr);
                }
                DPOCCPTFree(dpoccptr);
                
                // Clear result from read buffer
                SRAOCCInitialise(peArg->PEAlgnmtOutput_OccHandle->occCollector);
                DPOCCInitialise(peArg->dpArguments->dpOccCollector);
                
                threadArg->algnmtThreadStats.readNumHandled++;
            }
        }
        
        #ifdef DEBUG_MICA_PEWT_LOGGING
            printf ( "PEWT  -  _PEWTAlignmentRollBody  -  Finished Processing Payload.\n");
            printf ( "PEWT  -  _PEWTAlignmentRollBody  -  Handler Statistics. Total Number of %llu/%llu Reads Aligned.\n", 
                threadArg->algnmtThreadStats.readNumAligned,threadArg->algnmtThreadStats.readNumHandled);
        #endif

        // Timestamp'ed -----------
        timestamp = getElapsedTime(startTime);
        totalElipsedTime += timestamp - lastEventTime;
        threadArg->algnmtThreadStats.processTime += timestamp - lastEventTime;
        lastEventTime = timestamp;
        // ------------------------
        
        pthread_mutex_lock(mutexOutputQueue);

            #ifdef DEBUG_MICA_PEWT_LOGGING
                printf ( "PEWT  -  _PEWTAlignmentRollBody  -  Mutex Lock obtained(mutexOutputQueue). Updating Buffer Pointer.\n");
            #endif
            
            if (!SQEnqueue(writeThread->outputQueue,bufferIdx)) {
                printf("[SAFE-GUARD] Unexpected mutex error.\n");
                printf("The outputQueue is supposed to be non-full.\n");
            }
            
            pthread_cond_signal(condOutputQueue);

            #ifdef DEBUG_MICA_PEWT_LOGGING
                printf ( "PEWT  -  _PEWTAlignmentRollBody  -  Unlock mutex(mutexOutputQueue).\n");
                printf ( "PEWT  -  _PEWTAlignmentRollBody  -  Output Queue Status\n" );
                SQPrintQueue(writeThread->outputQueue);
            #endif
            
        pthread_mutex_unlock(mutexOutputQueue);
        fflush(0);
    }
    
    SRASettingFree(seedSetting);
    SRASettingFree(oneMoreMismatchSetting);
    
    SRAModelSetFree(sraModelSet);
    SRAModelSetFree(sraModelSet_seed);
    SRAModelSetFree(sraModelSet_extend);
    
    SRAAlgnmtOutputFreeSAMStore(peArg->PEAlgnmtOutput_OccHandle);
    SRAAlgnmtOutputFree(peArg->PEAlgnmtOutput_OccHandle);
    peArg->PEAlgnmtOutput_OccHandle = NULL;
    
    _PEWTFreePEArgument(peArg);
    
    writeThread->alignmentThreadBody[threadIdx] = 0;
    pthread_exit(0);
    return 0;
}


/////////////////////////////////////////////////////////////////////////
//
// OUTPUT THREAD BODY
//
/////////////////////////////////////////////////////////////////////////

unsigned int _PEWTPumpOutputFromBuffer(PEArguments * peArg, 
                                SRAOCCCollectorPointer * occptr, DPOCCCollectorPointer * dpoccptr) {

    unsigned int runPairOcc = 0;
                    
    SRAOCCCollector * resultOccCollector = peArg->PEAlgnmtOutput_OccHandle->occCollector;
    DPOCCCollector * resultDpOccCollector = peArg->dpArguments->dpOccCollector;

    SRAOccurrence * sraOcc = SRAOCCPTRead(occptr);
    DPOccurrence * dpOcc = DPOCCPTRead(dpoccptr);
    
    if ( sraOcc != NULL && sraOcc->type == SRAOCC_TYPE_DELIMITOR_READ) {
        SRAOCCPTNext(occptr);
        sraOcc = SRAOCCPTRead(occptr);
    }
    
    while ( sraOcc != NULL && sraOcc->type != SRAOCC_TYPE_DELIMITOR_READ ) {
        
        SRAOCCAddOccurrenceToBuckets(resultOccCollector,sraOcc);

        switch (sraOcc->type) {
            case SRAOCC_TYPE_PE_PAIR:
                runPairOcc++;
                break;
            case SRAOCC_TYPE_PE_DP_BASE_READ:
            case SRAOCC_TYPE_PE_DP_BASE_MATE:
                runPairOcc+=2;
                if ( dpOcc == NULL ) {
                    printf("[SAFE-GUARD] Insufficient number of DP occurrences for SRA output.\n");
                } else {
                    DPOCCAddOccurrenceToBuckets(resultDpOccCollector,dpOcc);
                    DPOCCPTNext(dpoccptr);
                    dpOcc = DPOCCPTRead(dpoccptr);
                }
                break;
            case SRAOCC_TYPE_PE_DP_SEED_OCC:
                runPairOcc+=2;
                if ( dpOcc == NULL ) {
                    printf("[SAFE-GUARD] Insufficient number of DP occurrences for SRA output.\n");
                } else {
                    DPOCCAddOccurrenceToBuckets(resultDpOccCollector,dpOcc);
                    DPOCCPTNext(dpoccptr);
                    dpOcc = DPOCCPTRead(dpoccptr);
                }
                
                if ( dpOcc == NULL ) {
                    printf("[SAFE-GUARD] Insufficient number of DP occurrences for SRA output.\n");
                } else {
                    DPOCCAddOccurrenceToBuckets(resultDpOccCollector,dpOcc);
                    DPOCCPTNext(dpoccptr);
                    dpOcc = DPOCCPTRead(dpoccptr);
                }
                break;
        }
        
        SRAOCCPTNext(occptr);
        sraOcc = SRAOCCPTRead(occptr);
    }
    
    PEOCCFlushCache(peArg);
    runPairOcc /= 2;
    
    return runPairOcc;
}

void * _PEWTOutputRollBody(void * _writeThread) {
    PEWriteThread * writeThread = (PEWriteThread*) _writeThread;
    
    pthread_mutex_t * mutexFreeQueue = &(writeThread->mutexFreeQueue);
    pthread_cond_t * condFreeQueue = &(writeThread->condFreeQueue);
    
    pthread_mutex_t * mutexAlignmentQueue = &(writeThread->mutexAlignmentQueue);
    pthread_cond_t * condAlignmentQueue = &(writeThread->condAlignmentQueue);
    
    pthread_mutex_t * mutexOutputQueue = &(writeThread->mutexOutputQueue);
    pthread_cond_t * condOutputQueue = &(writeThread->condOutputQueue);
    
    PEWTArgs * pewtArguments = &(writeThread->pewtArguments);
    
    unsigned char * charMap = pewtArguments->charMap;
    
    PEArguments * peArg = _PEWTCreatePEArgument(writeThread);
    
    unsigned long long i, j, k;
    int invalidStatus = FALSE;
    //int dummyQuality[SRA_MAX_READ_LENGTH];
    //memset(dummyQuality,0,sizeof(int)*SRA_MAX_READ_LENGTH);
    
    #ifdef DEBUG_MICA_PEWT_LOGGING
        printf ( "PEWT  -  _PEWTOutputRollBody  -  Invoked and Entering Main-Loop\n");
    #endif

    //ATTENTION
    writeThread->consPairCount = 0;
    writeThread->consPairRead = 0;
    writeThread->cpuConsSaCount = 0;
    writeThread->cpuConsReadCount = 0;
    writeThread->consReadCount = 0;
    
    // Timestamp'ed -----------
    double startTime = setStartTime();
    double lastEventTime = 0;
    double totalElipsedTime = 0;
    // ------------------------

    while (1) {

        int bufferIdx;
        
        // Block until there is an non-empty buffer
        pthread_mutex_lock(mutexOutputQueue);
        
            // Loop wait as long as
            // - alignment queue is empty
            // - AND free queue is not full OR thread is healthy
            while (
                SQIsEmpty(writeThread->outputQueue) &&
                (
                    !SQIsFull(writeThread->freeQueue) ||
                    !SQIsEmpty(writeThread->alignmentQueue) ||
                    writeThread->threadHealth==MICA_PE_WRITE_THREAD_HEALTH_OK
                )
            ) {
                pthread_cond_wait(condOutputQueue,mutexOutputQueue);
            }
            // Break in case
            // - free queue is full
            // - AND alignment queue is empty
            // - AND thread is not healthy
            if (
                SQIsFull(writeThread->freeQueue) &&
                SQIsEmpty(writeThread->alignmentQueue) &&
                SQIsEmpty(writeThread->outputQueue) &&
                writeThread->threadHealth!=MICA_PE_WRITE_THREAD_HEALTH_OK
            ) {

                #ifdef DEBUG_MICA_PEWT_LOGGING
                    printf ( "PEWT  -  _PEWTOutputRollBody  -  Thread is not healthy and Buffer is empty. Exiting.\n");
                #endif
                
                pthread_mutex_unlock(mutexOutputQueue);
                break;
            }
            
            if (!SQDequeue(writeThread->outputQueue,&bufferIdx)) {
                printf("[SAFE-GUARD] Unexpected mutex error.\n");
                printf("The outputQueue is supposed to be non-empty.\n");
            }
        pthread_mutex_unlock(mutexOutputQueue);

        // Timestamp'ed -----------
        double timestamp = getElapsedTime(startTime);
        lastEventTime = timestamp;
        // ------------------------
        
        #ifdef DEBUG_MICA_PEWT_LOGGING
            printf ( "PEWT  -  _PEWTOutputRollBody  -  Received Signal on condition(condOutputQueue)\n");
        #endif
        
        
        //Process Payload
        //
        #ifdef DEBUG_MICA_PEWT_LOGGING
            printf ( "PEWT  -  _PEWTOutputRollBody  -  Processing Payload from buffer %d.\n", bufferIdx);
        #endif
        
        unsigned int          queryInBatch          = writeThread->queryInBatch[bufferIdx];
        unsigned int          batchSize             = writeThread->batchSize[bufferIdx];
        unsigned int          numThreads            = writeThread->numThreads[bufferIdx];
        unsigned int          firstReadIdx          = writeThread->firstReadIdx[bufferIdx];
        int                   algnFuncReasonFlag    = writeThread->algnFuncReasonFlag[bufferIdx];
        
        if ( algnFuncReasonFlag == MIC_ALGNMT_FOR_UNDEFINED ) {
            printf("[SAFE-GUARD] Alignment function flag indicates undefined reason for batch!\n");
        }

        #ifdef DEBUG_MICA_PEWT_LOGGING
            printf ( "PEWT  -  _PEWTOutputRollBody  -  #Thread %u / BatchSize %u / #Query %u.\n",numThreads,batchSize,queryInBatch);
        #endif
        
        char                * readNameBufferFrame   = writeThread->readNameBufferFrame[bufferIdx];
        unsigned char       * readBodyBufferFrame   = writeThread->readBodyBufferFrame[bufferIdx];
        char                * readQualityBufferFrame= writeThread->readQualityBufferFrame[bufferIdx];
        uint16_t            * readLengthBufferFrame = writeThread->readLengthBufferFrame[bufferIdx];
        char                * mateNameBufferFrame   = writeThread->mateNameBufferFrame[bufferIdx];
        unsigned char       * mateBodyBufferFrame   = writeThread->mateBodyBufferFrame[bufferIdx];
        char                * mateQualityBufferFrame= writeThread->mateQualityBufferFrame[bufferIdx];
        uint16_t            * mateLengthBufferFrame = writeThread->mateLengthBufferFrame[bufferIdx];
        
        PEReadOutput        * peBatchOutput         = writeThread->peBatchOutput[bufferIdx];
            
        uint8_t               outputType            = writeThread->outputType[bufferIdx];
        uint16_t            * occCount              = writeThread->occCount[bufferIdx];
        uint8_t             * outputBufferStatus    = writeThread->outputBufferStatus[bufferIdx];
        unsigned int        * outputBufferFrame     = writeThread->outputBufferFrame[bufferIdx];
        MICSRAOccMetadata   * outputBufferMeta      = writeThread->outputBufferMeta[bufferIdx];
        MICDPOccurrence     * dpOutputBufferFrame   = writeThread->dpOutputBufferFrame[bufferIdx];
        PEMappingQuality    * mappingQualities      = writeThread->mappingQualities[bufferIdx];
        
        // Set up the temporary output collector for the intermediate results
        SRAOCCCollector     * cpuSraOccCollector    = writeThread->cpuSraOccCollector[bufferIdx];
        DPOCCCollector      * cpuDpOccCollector     = writeThread->cpuDpOccCollector[bufferIdx];

        SRAOCCCollectorPointer * occptr = SRAOCCPTCreate(cpuSraOccCollector);
        DPOCCCollectorPointer * dpoccptr = DPOCCPTCreate(cpuDpOccCollector);
        
        SRAArguments * sraReadArgs = peArg->sraArgsPos;
        SRAArguments * sraReadArgs_neg = peArg->sraArgsNeg;
        SRAArguments * sraMateArgs = peArg->sraArgsPos_mate;
        SRAArguments * sraMateArgs_neg = peArg->sraArgsNeg_mate;
        
        unsigned char readNegPattern[SRA_MAX_READ_LENGTH], mateNegPattern[SRA_MAX_READ_LENGTH];
        unsigned char readNegQuality[SRA_MAX_READ_LENGTH], mateNegQuality[SRA_MAX_READ_LENGTH];
        DPOccurrence dpOcc_1, dpOcc_2;
        DPOCCCollector * dpOccCollector = peArg->dpArguments->dpOccCollector;
        
        // Populate input parameter and output handler for batch
        peArg->PEAlgnmtOutput_OccHandle = peBatchOutput->outputHandler;
        
        for (i=0;i<numThreads;i++) {
            unsigned int offset = i * batchSize;
            unsigned int * outputPtr     = &outputBufferFrame[i*MIC_SRA_OUTPUT_SIZE_PER_THREAD];
            MICSRAOccMetadata * metaPtr  = &outputBufferMeta[i*MIC_SRA_OUTPUT_META_PER_THREAD];
            MICDPOccurrence * dpOutPtr   = &dpOutputBufferFrame[i*MIC_DP_OUTPUT_SIZE_PER_THREAD];
            
            for (j=0;j<batchSize;j++) {
                unsigned int readIdxInBatch = offset + j;

                // Buffer could be partially filled.
                // Early termination w.r.t to the size of the particular buffer
                if (readIdxInBatch >= queryInBatch) {
                    break;
                }
                unsigned int readIdx = firstReadIdx + readIdxInBatch;
            
                #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                    printf("[INFO] Handling Read %u\n",readIdx);
                #endif
                
                unsigned char * readPatternPtr = &readBodyBufferFrame[readIdxInBatch*SRA_MAX_READ_LENGTH];
                char * readPatternName         = &readNameBufferFrame[readIdxInBatch*MAX_SEQ_NAME_LENGTH];
                uint16_t readPatternLen        = readLengthBufferFrame[readIdxInBatch];
                char * readQualityPtr         = &readQualityBufferFrame[readIdxInBatch*SRA_MAX_READ_LENGTH];
                SRAFlipPattern(charMap,readPatternPtr,readPatternLen,readNegPattern);
                SRAFlipQuality(charMap,readQualityPtr,readPatternLen,readNegQuality);
                
                unsigned char * matePatternPtr = &mateBodyBufferFrame[readIdxInBatch*SRA_MAX_READ_LENGTH];
                char * matePatternName         = &mateNameBufferFrame[readIdxInBatch*MAX_SEQ_NAME_LENGTH];
                uint16_t matePatternLen        = mateLengthBufferFrame[readIdxInBatch];
                char * mateQualityPtr         = &mateQualityBufferFrame[readIdxInBatch*SRA_MAX_READ_LENGTH];
                SRAFlipPattern(charMap,matePatternPtr,matePatternLen,mateNegPattern);
                SRAFlipQuality(charMap,mateQualityPtr,matePatternLen,mateNegQuality);
                
                // Populate the SRA information
                SRAQueryInfoPopulate(sraReadArgs->QueryInfo,readIdx,readPatternName,readPatternLen,readPatternPtr,readPatternPtr,QUERY_POS_STRAND,readQualityPtr);
                SRAQueryInfoPopulate(sraReadArgs_neg->QueryInfo,readIdx,readPatternName,readPatternLen,readNegPattern,readPatternPtr,QUERY_NEG_STRAND,readNegQuality);
                
                SRAQueryInfoPopulate(sraMateArgs->QueryInfo,readIdx,matePatternName,matePatternLen,matePatternPtr,matePatternPtr,QUERY_POS_STRAND,mateQualityPtr);
                SRAQueryInfoPopulate(sraMateArgs_neg->QueryInfo,readIdx,matePatternName,matePatternLen,mateNegPattern,matePatternPtr,QUERY_NEG_STRAND,mateNegQuality);
                
                // Initialise the result counter and MAPQ buffers
                SRAResultCountInitialise(sraReadArgs->AlgnmtStats);
                SRAResultCountInitialise(sraMateArgs->AlgnmtStats);

                //Assuming status is valid prior to parsing
                invalidStatus = FALSE;
                
                if ( outputType == MICA_PE_WRITE_THREAD_PAYLOAD_OUTPUT_TYPE_UNKNOWN) {
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Payload Error - Unknown - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    
                    printf("[INFO] Payload Error - Unknown - %s / %s!\n",readPatternName,matePatternName);
                    // In this case the read has never been initialised.
                    // therefore the read always requires CPU processing.

                    unsigned int runPairOcc = _PEWTPumpOutputFromBuffer(peArg,occptr,dpoccptr);
                    
                    if (runPairOcc) {
                        writeThread->consPairCount += runPairOcc;
                        writeThread->cpuConsSaCount += runPairOcc;
                        writeThread->consPairRead++;
                        // Mark this handled by CPU
                        writeThread->cpuConsReadCount++;
                    }
                    
                    // Mark this read as invalid status
                    invalidStatus = TRUE;
                    
                    // As the output buffer are invalid
                    occCount[readIdxInBatch] = 0;
                } else if ( outputType == MICA_PE_WRITE_THREAD_PAYLOAD_OUTPUT_TYPE_SKIP_MIC) {
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Payload unaligned by MIC - MIC unaligned - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    // In this case the read has never been passed into MIC
                    // therefore the read always requires CPU processing.

                    unsigned int runPairOcc = _PEWTPumpOutputFromBuffer(peArg,occptr,dpoccptr);
                    
                    if (runPairOcc) {
                        writeThread->consPairCount += runPairOcc;
                        writeThread->cpuConsSaCount += runPairOcc;
                        writeThread->consPairRead++;
                        // Mark this handled by CPU
                        writeThread->cpuConsReadCount++;
                    }
                    
                    // As the output buffer are invalid
                    occCount[readIdxInBatch] = 0;
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_CLOSED) {
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Unaligned closed read - MIC unaligned - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    // In this case the output buffer is flooded as the number of output
                    // generated is larger than the MIC_SRA_OUTPUT_MAX_ALIGNMENT; or the number of header
                    // generated is larger than MIC_SRA_OUTPUT_MAX_META. MIC was not
                    // able to handle all the alignment result hence declared the read 
                    // unaligned on MIC.
                    //cpuConsSaCount += ProcessReadDoubleStrand(sraReadArgs,sraReadArgs_neg,sraModelSet,sraModelSet);

                    unsigned int runPairOcc = _PEWTPumpOutputFromBuffer(peArg,occptr,dpoccptr);
                    
                    if (runPairOcc) {
                        writeThread->consPairCount += runPairOcc;
                        writeThread->cpuConsSaCount += runPairOcc;
                        writeThread->consPairRead++;
                        // Mark this handled by CPU
                        writeThread->cpuConsReadCount++;
                    }
                    
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_UNHANDLE) {
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Unhandled read - MIC unaligned - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    // In this case the entire core shared output buffer is flooded as the number of output
                    // generated is larger than the MIC_SRA_OUTPUT_SIZE_PER_THREAD. MIC was not
                    // able to handle all the reads in the input read block hence declared all consecutive read 
                    // unaligned on MIC.
                    //cpuConsSaCount += ProcessReadDoubleStrand(sraReadArgs,sraReadArgs_neg,sraModelSet,sraModelSet);

                    unsigned int runPairOcc = _PEWTPumpOutputFromBuffer(peArg,occptr,dpoccptr);
                    
                    if (runPairOcc) {
                        writeThread->consPairCount += runPairOcc;
                        writeThread->cpuConsSaCount += runPairOcc;
                        writeThread->consPairRead++;
                        // Mark this handled by CPU
                        writeThread->cpuConsReadCount++;
                    }
                    
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_SKIPPED) {
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Skipped read - MIC unaligned - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    if ( readIdxInBatch < queryInBatch ) {
                        printf("[INFO] Skipped read - MIC unaligned - %s / %s!\n",readPatternName,matePatternName);
                    }
                    // In this case the MIC actively skipped the alignment of this read. Reason being the 
                    // read contains invalid character or unknown reason. CPU should also skip it.
                    // These sample code gives output when these reads are found from the output block and attempt to process it.
                    // cpuConsSaCount += ProcessReadDoubleStrand(sraReadArgs,sraReadArgs_neg,sraModelSet,sraModelSet);
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_OUTPUT_STATUS_OPEN) {
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Skipped read - MIC still open - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    printf("[ERROR] Skipped read - MIC still open - %s / %s!\n",readPatternName,matePatternName);
                    // In this case the MIC actively skipped the alignment of this read. For unknown reason.
                    //cpuConsSaCount += ProcessReadDoubleStrand(sraReadArgs,sraReadArgs_neg,sraModelSet,sraModelSet);
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_PAIR) {
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Proper pair - MIC aligned - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    // Ticket#58
                    // This indicates that the pair-end matching engine returns some valid pair-end alignment. 
                    // The output buffer should be populated with the alignment results.

                    if (occCount[readIdxInBatch]>0) {
                    
                        // ATTENTION SAFE GUARD
                        if ( occCount[readIdxInBatch] % 2 != 0 ) {
                            printf("[SAFEGUARD] MIC_PE_OUTPUT_STATUS_PAIR Read acts abnormal - %s / %s!\n",readPatternName,matePatternName);
                            printf("Invalid number of output in buffer occCount[readIdxInBatch] = %u\n",occCount[readIdxInBatch]);
                        }
                        
                        for (k=0;k<occCount[readIdxInBatch];) {
                            PEOCCAddTextPositionToCache(peArg,SRAOCC_TYPE_PE_PAIR,
                                                    outputPtr[k], metaPtr[k].strand + 1,
                                                    metaPtr[k].matchLen, metaPtr[k].numOfErr, 0,
                                                    metaPtr[k].errors);
                            k++;
                            
                            PEOCCAddTextPositionToCache(peArg,SRAOCC_TYPE_PE_PAIR,
                                                    outputPtr[k], metaPtr[k].strand + 1,
                                                    metaPtr[k].matchLen, metaPtr[k].numOfErr, 0,
                                                    metaPtr[k].errors);
                            k++;
                        }
                        PEOCCMAPQValue(peArg,mappingQualities[readIdxInBatch].readQuality,mappingQualities[readIdxInBatch].mateQuality);
                        PEOCCFlushCache(peArg);
                        
                        writeThread->consPairCount += occCount[readIdxInBatch]/2;
                        writeThread->consPairRead++;
                    } else {
                        // ATTENTION SAFE GUARD
                        printf("[SAFEGUARD] MIC_PE_OUTPUT_STATUS_PAIR Read acts abnormal - %s / %s!\n",readPatternName,matePatternName);
                        printf("Invalid number of output in buffer occCount[readIdxInBatch] = %u\n",occCount[readIdxInBatch]);
                    }
                    
                    // Mark this handled by MIC
                    writeThread->consReadCount++;
                    
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_BAD_PAIR) {
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Both SRA aligned but unable to PE pair - MIC unaligned - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    // Ticket#58
                    // This indicates that the pair-end matching engine does NOT give any valid pair-end alignment; 
                    // yet both of the read and mate are aligned to the reference sequence. In this case, 
                    // we randomly output one of the single-end alignment from each read and mate with the least 
                    // number of mismatch. The output buffer should be populated such pair of alignment.
                
                    k = 0;
                    PEOCCAddTextPositionToCache(peArg,SRAOCC_TYPE_PE_BAD_PAIR,
                                            outputPtr[k], metaPtr[k].strand + 1,
                                            metaPtr[k].matchLen, metaPtr[k].numOfErr, 0,
                                            metaPtr[k].errors);
                    k++;
                    PEOCCAddTextPositionToCache(peArg,SRAOCC_TYPE_PE_BAD_PAIR,
                                            outputPtr[k], metaPtr[k].strand + 1,
                                            metaPtr[k].matchLen, metaPtr[k].numOfErr, 0,
                                            metaPtr[k].errors);
                    PEOCCMAPQValue(peArg,mappingQualities[readIdxInBatch].readQuality,mappingQualities[readIdxInBatch].mateQuality);
                    PEOCCFlushCache(peArg);
                                   
                    // Mark this handled by MIC
                    writeThread->consReadCount++;
                    
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_READ_NO_ALIGNMENT) {
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Left SRA unaligned hence unable to PE pair - MIC unaligned - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    
                    // Ticket#58
                    // This indicates that the READ is unaligned. Hence implies PE is unperformed. 
                    // The output buffer should be populated with the SRA alignment results from MATE.
                    if (occCount[readIdxInBatch]>0) {
                        for (k=0;k<occCount[readIdxInBatch];k++) {
                            PEOCCAddTextPositionToCache(peArg,SRAOCC_TYPE_PE_READ_READ_NO_ALIGNMENT,
                                                    outputPtr[k], metaPtr[k].strand + 1,
                                                    metaPtr[k].matchLen, metaPtr[k].numOfErr, 0,
                                                    metaPtr[k].errors);
                            PEOCCAddTextPositionToCache(peArg,SRAOCC_TYPE_PE_MATE_MATE_NO_ALIGNMENT,
                                                    outputPtr[k], metaPtr[k].strand + 1,
                                                    metaPtr[k].matchLen, metaPtr[k].numOfErr, 0,
                                                    metaPtr[k].errors);
                        }
                        PEOCCMAPQValue(peArg,mappingQualities[readIdxInBatch].readQuality,mappingQualities[readIdxInBatch].mateQuality);
                        PEOCCFlushCache(peArg);
                    } else {
                        // ATTENTION SAFE GUARD
                        printf("[SAFEGUARD] MIC_PE_OUTPUT_STATUS_READ_NO_ALIGNMENT Read acts abnormal - %s / %s!\n",readPatternName,matePatternName);
                        printf("Invalid number of output in buffer occCount[readIdxInBatch] = %u\n",occCount[readIdxInBatch]);
                    }
                    
                    // Mark this handled by MIC
                    writeThread->consReadCount++;
                    
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_MATE_NO_ALIGNMENT) {
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Right SRA unaligned hence unable to PE pair - MIC unaligned - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    
                    // Ticket#58
                    // This indicates that the MATE is unaligned. Hence implies PE is unperformed. 
                    // The output buffer should be populated with the SRA alignment results from READ.

                    if (occCount[readIdxInBatch]>0) {
                        for (k=0;k<occCount[readIdxInBatch];k++) {
                            PEOCCAddTextPositionToCache(peArg,SRAOCC_TYPE_PE_READ_MATE_NO_ALIGNMENT,
                                                    outputPtr[k], metaPtr[k].strand + 1,
                                                    metaPtr[k].matchLen, metaPtr[k].numOfErr, 0,
                                                    metaPtr[k].errors);
                            PEOCCAddTextPositionToCache(peArg,SRAOCC_TYPE_PE_MATE_READ_NO_ALIGNMENT,
                                                    outputPtr[k], metaPtr[k].strand + 1,
                                                    metaPtr[k].matchLen, metaPtr[k].numOfErr, 0,
                                                    metaPtr[k].errors);
                        }
                        PEOCCMAPQValue(peArg,mappingQualities[readIdxInBatch].readQuality,mappingQualities[readIdxInBatch].mateQuality);
                        PEOCCFlushCache(peArg);
                    } else {
                        // ATTENTION SAFE GUARD
                        printf("[SAFEGUARD] MIC_PE_OUTPUT_STATUS_MATE_NO_ALIGNMENT Read acts abnormal - %s / %s!\n",readPatternName,matePatternName);
                        printf("Invalid number of output in buffer occCount[readIdxInBatch] = %u\n",occCount[readIdxInBatch]);
                    }
                    
                    // Mark this handled by MIC
                    writeThread->consReadCount++;
                    
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_BOTH_NO_ALIGNMENT) {
                
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Both SRA unaligned hence unable to PE pair - MIC unaligned - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    
                    // Ticket#58
                    // This indicates none of READ or MATE is aligned. No output is written to the buffer.
                    if (occCount[readIdxInBatch]>0) {
                        // ATTENTION SAFE GUARD
                        printf("[SAFEGUARD] MIC_PE_OUTPUT_STATUS_MATE_NO_ALIGNMENT Read acts abnormal - %s / %s!\n",readPatternName,matePatternName);
                        printf("Invalid number of output in buffer occCount[readIdxInBatch] = %u\n",occCount[readIdxInBatch]);
                    }
                    
                    PEOCCReportNoAlignment(peArg,SRAOCC_TYPE_PE_READ_BOTH_NO_ALIGNMENT);
                    PEOCCReportNoAlignment(peArg,SRAOCC_TYPE_PE_MATE_BOTH_NO_ALIGNMENT);
                    PEOCCMAPQValue(peArg,mappingQualities[readIdxInBatch].readQuality,mappingQualities[readIdxInBatch].mateQuality);
                    PEOCCFlushCache(peArg);
                    
                    // Mark this handled by MIC
                    writeThread->consReadCount++;
                    
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_DP_BASE_MATE) {
                
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Left SRA unaligned PE found by DefaultDP - MIC aligned - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    
                    // Ticket#56
                    // This indicates the READ is a result of DP based on estimation made with MATE
                    
                    if (occCount[readIdxInBatch]>0) {
                        for (k=0;k<occCount[readIdxInBatch];k++) {
                            MICDPOccurrenceConvert(&dpOutPtr[k],&dpOcc_1);
                            DPOCCAddOccurrenceToBuckets(dpOccCollector,&dpOcc_1);
                            PEOCCAddTextPositionToCache(peArg,SRAOCC_TYPE_PE_DP_BASE_MATE,
                                                    outputPtr[k], metaPtr[k].strand + 1,
                                                    metaPtr[k].matchLen, metaPtr[k].numOfErr, 0,
                                                    metaPtr[k].errors);
                        }
                        PEOCCMAPQValue(peArg,mappingQualities[readIdxInBatch].readQuality,mappingQualities[readIdxInBatch].mateQuality);
                        PEOCCFlushCache(peArg);
                        
                        writeThread->consPairCount += occCount[readIdxInBatch];
                        writeThread->consPairRead++;
                    } else {
                        // ATTENTION SAFE GUARD
                        printf("[SAFEGUARD] MIC_PE_OUTPUT_STATUS_DP_BASE_MATE Read acts abnormal - %s / %s!\n",readPatternName,matePatternName);
                        printf("Invalid number of output in buffer occCount[readIdxInBatch] = %u\n",occCount[readIdxInBatch]);
                    }
                    
                    dpOutPtr += occCount[readIdxInBatch];
                    
                    // Mark this handled by MIC
                    writeThread->consReadCount++;
                    
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_DP_BASE_READ) {
                
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] Right SRA unaligned PE found by DefaultDP - MIC aligned - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    
                    // Ticket#56
                    // This indicates the MATE is a result of DP based on estimation made with READ

                    if (occCount[readIdxInBatch]>0) {
                        for (k=0;k<occCount[readIdxInBatch];k++) {
                            MICDPOccurrenceConvert(&dpOutPtr[k],&dpOcc_1);
                            DPOCCAddOccurrenceToBuckets(dpOccCollector,&dpOcc_1);
                            PEOCCAddTextPositionToCache(peArg,SRAOCC_TYPE_PE_DP_BASE_READ,
                                                    outputPtr[k], metaPtr[k].strand + 1,
                                                    metaPtr[k].matchLen, metaPtr[k].numOfErr, 0,
                                                    metaPtr[k].errors);
                        }
                        PEOCCMAPQValue(peArg,mappingQualities[readIdxInBatch].readQuality,mappingQualities[readIdxInBatch].mateQuality);
                        PEOCCFlushCache(peArg);
                        
                        writeThread->consPairCount += occCount[readIdxInBatch];
                        writeThread->consPairRead++;
                    } else {
                        // ATTENTION SAFE GUARD
                        printf("[SAFEGUARD] MIC_PE_OUTPUT_STATUS_DP_BASE_READ Read acts abnormal - %s / %s!\n",readPatternName,matePatternName);
                        printf("Invalid number of output in buffer occCount[readIdxInBatch] = %u\n",occCount[readIdxInBatch]);
                    }
                    
                    dpOutPtr += occCount[readIdxInBatch];
                    
                    // Mark this handled by MIC
                    writeThread->consReadCount++;
                    
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_DP_MIX_BASE) {
                
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] PE found by New-Default - MIC aligned - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    
                    // Ticket#56
                    // This indicates the combination of the above two cases. The output contains the following:
                    // - DP found occurrences based on the estimated position by READ; and
                    // - Separater which is all zero
                    // - DP found occurrences based on the estimated position by MATE
                    
                    if (occCount[readIdxInBatch]==0) {
                        // ATTENTION SAFE GUARD
                        printf("[SAFEGUARD] MIC_PE_OUTPUT_STATUS_DP_MIX_BASE Read acts abnormal - %s / %s!\n",readPatternName,matePatternName);
                        printf("Invalid number of output in buffer occCount[readIdxInBatch] = %u\n",occCount[readIdxInBatch]);
                    }
                    
                    int dpBase;
                    int readCount = -1;
                    for (k=0;k<occCount[readIdxInBatch];k++) {
                        if (dpOutPtr[k].ambPosition == 0 &&
                            dpOutPtr[k].strand == 0 &&
                            dpOutPtr[k].matchElemsCount == 0 &&
                            dpOutPtr[k].matchLen == 0 &&
                            outputPtr[k] == 0) {
                            
                            readCount = k;
                            break;
                            
                        }
                    }
                    int mateCount = occCount[readIdxInBatch] - readCount - 1;
                    
                    if (readCount == 0 && mateCount == 0) {
                        // ATTENTION SAFE GUARD
                        printf("[SAFEGUARD] MIC_PE_OUTPUT_STATUS_DP_MIX_BASE Read acts abnormal - %s / %s!\n",readPatternName,matePatternName);
                        printf("[SAFEGUARD] MIC_PE_OUTPUT_STATUS_DP_MIX_BASE Read acts abnormal - %s / %s!\n",readPatternName,matePatternName);
                        printf("Invalid number of output in buffer occCount[readIdxInBatch] = %u\n",occCount[readIdxInBatch]);
                        printf("Yet found the count of occurrences of both sides be %u and %u\n",readCount,mateCount);
                    }
                    
                    // Base READ DP occurrences
                    k=0;
                    if (readCount>0) {
                        dpBase = SRAOCC_TYPE_PE_DP_BASE_READ;
                        for (;k<readCount;k++) {
                            MICDPOccurrenceConvert(&dpOutPtr[k],&dpOcc_1);
                            DPOCCAddOccurrenceToBuckets(dpOccCollector,&dpOcc_1);
                            PEOCCAddTextPositionToCache(peArg,SRAOCC_TYPE_PE_DP_BASE_READ,
                                                    outputPtr[k], metaPtr[k].strand + 1,
                                                    metaPtr[k].matchLen, metaPtr[k].numOfErr, 0,
                                                    metaPtr[k].errors);
                        }
                        k++;
                    }
                    
                    // Base MATE DP occurrences
                    dpBase = SRAOCC_TYPE_PE_DP_BASE_MATE;
                    for (;k<mateCount;k++) {
                        MICDPOccurrenceConvert(&dpOutPtr[k],&dpOcc_1);
                        DPOCCAddOccurrenceToBuckets(dpOccCollector,&dpOcc_1);
                        PEOCCAddTextPositionToCache(peArg,SRAOCC_TYPE_PE_DP_BASE_MATE,
                                                outputPtr[k], metaPtr[k].strand + 1,
                                                metaPtr[k].matchLen, metaPtr[k].numOfErr, 0,
                                                metaPtr[k].errors);
                    }
                    PEOCCMAPQValue(peArg,mappingQualities[readIdxInBatch].readQuality,mappingQualities[readIdxInBatch].mateQuality);
                    PEOCCFlushCache(peArg);
                    
                    writeThread->consPairCount += occCount[readIdxInBatch] - 1;
                    writeThread->consPairRead++;
                    
                    dpOutPtr += occCount[readIdxInBatch];
                    
                    // Mark this handled by MIC
                    writeThread->consReadCount++;
                    
                } else if ( outputBufferStatus[readIdxInBatch] == MIC_PE_OUTPUT_STATUS_DP_SEED) {
                    
                    #ifdef DEBUG_MICA_PEWT_READ_HANDLING
                        printf("[DEBUG] PE found by DeepDp - MIC aligned - %s / %s!\n",readPatternName,matePatternName);
                    #endif
                    
                    // Ticket#56
                    // This case is the SeedDP (DeepDP) output. The output contains the following:
                    // - Exactly one integer on output buffer and it stores the number of alignment (2 for each pair-end alignment)
                    // All PE alignment can be retrieved from the DP buffer.
                    
                    unsigned int peCount = outputPtr[0];
                    
                    if (occCount[readIdxInBatch] != 1) {
                        // ATTENTION SAFE GUARD
                        printf("[SAFEGUARD] MIC_PE_OUTPUT_STATUS_DP_SEED Read acts abnormal - %s / %s!\n",readPatternName,matePatternName);
                        printf("Invalid number of output in buffer occCount[readIdxInBatch] = %u. Expect 1.\n",occCount[readIdxInBatch]);
                        break;
                    }
                    
                    if (peCount>MIC_DP_OUTPUT_MAX_ALIGNMENT) {
                        // ATTENTION SAFE GUARD
                        printf("[SAFEGUARD] MIC_PE_OUTPUT_STATUS_DP_SEED Read acts abnormal - %s / %s!\n",readPatternName,matePatternName);
                        printf("Invalid number of output in buffer peCount = %u. Expect <=%u.\n",peCount,MIC_DP_OUTPUT_MAX_ALIGNMENT);
                        break;
                    }
                    
                    
                    for (k=0;k<peCount;){
                        MICDPOccurrenceConvert(&dpOutPtr[k],&dpOcc_1);
                        MICDPOccurrenceConvert(&dpOutPtr[k+1],&dpOcc_2);
                        PEDPOCCDumpOneSeedAlignment(peArg,&dpOcc_1,&dpOcc_2);
                        k+=2;
                    }
                    PEOCCMAPQValue(peArg,mappingQualities[readIdxInBatch].readQuality,mappingQualities[readIdxInBatch].mateQuality);
                    PEOCCFlushCache(peArg);
                    writeThread->consPairCount += peCount/2;
                    writeThread->consPairRead++;
                    
                    dpOutPtr += peCount;
                    
                    // Mark this handled by MIC
                    writeThread->consReadCount++;
                    
                } else {
                    printf("[SAFEGUARD] Read return unexpected output status - %s / %s - %u!\n",readPatternName,matePatternName,outputBufferStatus[readIdxInBatch]);
                    invalidStatus = TRUE;
                }
                
                if (invalidStatus==FALSE) {
                    // Increase status counter
                    writeThread->performanceStats.readNumHandled++;
                    writeThread->performanceStats.readNumByHandlers[outputBufferStatus[readIdxInBatch]]++;
                }
                
                outputPtr += occCount[readIdxInBatch];
                metaPtr += occCount[readIdxInBatch];
            }

        }
        
        SRAOCCInitialise(cpuSraOccCollector);
        DPOCCInitialise(cpuDpOccCollector);
        
        SRAOCCPTFree(occptr);
        DPOCCPTFree(dpoccptr);
        
        PEReadThread_ClosingOutput(peBatchOutput);
        
        #ifdef DEBUG_MICA_PEWT_LOGGING
            printf ( "PEWT  -  _PEWTOutputRollBody  -  Finished Processing Payload.\n");
            printf ( "PEWT  -  _PEWTOutputRollBody  -  Handler Statistics. Total Number of %llu Reads Handled.\n", writeThread->performanceStats.readNumHandled);
            printf ( "PEWT  -  _PEWTOutputRollBody  -  Breakdown by handler.\n");
            
            for (k=0;k<MIC_PE_OUTPUT_STATUS_COUNT;k++) {
                printf("%20s %llu\n",MIC_PE_OUTPUT_STATUS[k],writeThread->performanceStats.readNumByHandlers[k]);
            }
        #endif

        // Timestamp'ed -----------
        timestamp = getElapsedTime(startTime);
        totalElipsedTime += timestamp - lastEventTime;
        writeThread->performanceStats.processTime += timestamp - lastEventTime;
        lastEventTime = timestamp;
        // ------------------------
        
        pthread_mutex_lock(mutexFreeQueue);

            #ifdef DEBUG_MICA_PEWT_LOGGING
                printf ( "PEWT  -  _PEWTOutputRollBody  -  Mutex Lock obtained(mutexFreeQueue). Updating Buffer Pointer.\n");
            #endif
            
            if (!SQEnqueue(writeThread->freeQueue,bufferIdx)) {
                printf("[SAFE-GUARD] Unexpected mutex error.\n");
                printf("The freeQueue is supposed to be non-full.\n");
            }
            
            pthread_cond_signal(condAlignmentQueue);
            pthread_cond_signal(condFreeQueue);

            #ifdef DEBUG_MICA_PEWT_LOGGING
                printf ( "PEWT  -  _PEWTOutputRollBody  -  Unlock mutex(mutexFreeQueue).\n");
                printf ( "PEWT  -  _PEWTOutputRollBody  -  Free Queue Status\n" );
                SQPrintQueue(writeThread->freeQueue);
            #endif
            
        pthread_mutex_unlock(mutexFreeQueue);
        fflush(0);
    }

    peArg->PEAlgnmtOutput_OccHandle = NULL;
    
    _PEWTFreePEArgument(peArg);
    
    writeThread->outputThreadBody = 0;
    pthread_exit(0);
    return 0;
}

void PEWTPrintStats (PEWriteThread * writeThread, Logging * logger) {
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "PEWT - Consld-PE Pair Count    = %15llu / %15llu\n", writeThread->consPairCount, writeThread->consPairRead);
    // Found by MIC figures are the occurrences/reads that are solely handled by MIC(in MICCT) without CPU assistant
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "     -   Found by MIC          = %15llu / %15llu\n", writeThread->consPairCount - writeThread->cpuConsSaCount, writeThread->consPairRead - writeThread->cpuConsReadCount);
    // Found by CPU figures are the occurrences/reads that are either unhandled by MICCT/or allocated by scheduler into CPUCT directly(no MICCT involvement).
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "     -   Found by CPU          = %15llu / %15llu\n", writeThread->cpuConsSaCount,writeThread->cpuConsReadCount);
    
    int i;
    for (i=0;i<writeThread->numAlignmentThread;i++) {
        //Simplified view
        LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "PEWT - Alignment Thread Effort = %15llu ( %9.4f rps )\n",
            writeThread->algnmtThreadArg[i].algnmtThreadStats.readNumHandled,
            writeThread->algnmtThreadArg[i].algnmtThreadStats.readNumHandled / writeThread->algnmtThreadArg[i].algnmtThreadStats.processTime
            );
        
        //Detailed view
        /*LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "PEWT - Alignment Thread Effort = %15llu\n", writeThread->algnmtThreadArg[i].algnmtThreadStats.readNumHandled);
        LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "PEWT - CPU Aligner Enrich Time = %9.4f seconds ( %9.4f rps )\n", 
            writeThread->algnmtThreadArg[i].algnmtThreadStats.processTime,
            writeThread->algnmtThreadArg[i].algnmtThreadStats.readNumHandled / writeThread->algnmtThreadArg[i].algnmtThreadStats.processTime);*/
    }
    //Simplified view
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "PEWT - Output Thread Effort    = %15llu ( %9.4f rps )\n", 
        writeThread->performanceStats.readNumHandled,
        writeThread->performanceStats.readNumHandled / writeThread->performanceStats.processTime
    );
    
    //Detailed view
    /*LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "PEWT - Output Thread Effort    = %15llu\n", writeThread->performanceStats.readNumHandled);
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "PEWT - Output Time             = %9.4f seconds ( %9.4f rps )\n\n", 
        writeThread->performanceStats.processTime, 
        writeThread->performanceStats.readNumHandled / writeThread->performanceStats.processTime);*/
    
}

void PEWTRoll(PEWriteThread * writeThread) {
    
    pthread_t * alignmentThreadBody = writeThread->alignmentThreadBody;
    pthread_t * outputThreadBody    = &(writeThread->outputThreadBody);
    int i;
    
    for (i=0;i<writeThread->numAlignmentThread;i++) {
        if (alignmentThreadBody[i]!=0) {
            printf("[WriteThread] Unable to re-roll a PE Write Thread when an alignment thread running.\n");
            return;
        }
    }
    
    if ((*outputThreadBody)!=0) {
        printf("[WriteThread] Unable to re-roll a PE Write Thread when an output thread is running.\n");
        return;
    }

    PEWTSignal(writeThread,MICA_PE_WRITE_THREAD_HEALTH_OK);
    
    for (i=0;i<writeThread->numAlignmentThread;i++) {
        if (pthread_create(&alignmentThreadBody[i], NULL, _PEWTAlignmentRollBody, (void*) &(writeThread->algnmtThreadArg[i]))) {
            printf("[WriteThread] Unable to create thread!\n");
            exit(1);
        }
    }
    
    if (pthread_create(outputThreadBody, NULL, _PEWTOutputRollBody, (void*) writeThread)) {
        printf("[WriteThread] Unable to create thread!\n");
        exit(1);
    }
}

void PEWTSignal(PEWriteThread * writeThread, int threadHealth) {
    
    pthread_mutex_t * mutexThreadHealth = &(writeThread->mutexThreadHealth);
    pthread_cond_t * condAlignmentQueue = &(writeThread->condAlignmentQueue);
    pthread_cond_t * condOutputQueue = &(writeThread->condOutputQueue);
    
    pthread_mutex_lock(mutexThreadHealth);
    
        #ifdef DEBUG_MICA_PEWT_LOGGING
            printf ( "PEWT  -  PEWTSignal  -  Mutex Lock obtained(mutexThreadHealth). Updating Thread Health to %d.\n",threadHealth);
        #endif
        
        writeThread->threadHealth = threadHealth;
        pthread_cond_signal(condAlignmentQueue);
        pthread_cond_signal(condOutputQueue);
        
        #ifdef DEBUG_MICA_PEWT_LOGGING
            printf ( "PEWT  -  PEWTSignal  -  Unlock mutex(mutexThreadHealth).\n");
        #endif
        
    pthread_mutex_unlock(mutexThreadHealth);

    #ifdef DEBUG_MICA_PEWT_LOGGING
        printf ( "PEWT  -  PEWTSignal  -  Roll Body has returned.\n");
        printf ( "PEWT  -  PEWTSignal  -  Free Queue Status\n" );
        SQPrintQueue(writeThread->freeQueue);
        printf ( "PEWT  -  PEWTSignal  -  Alignment Queue Status\n" );
        SQPrintQueue(writeThread->alignmentQueue);
        printf ( "PEWT  -  PEWTSignal  -  Output Queue Status\n" );
        SQPrintQueue(writeThread->outputQueue);
    #endif
}

void PEWTThreadJoin(PEWriteThread * writeThread) {

    #ifdef DEBUG_MICA_PEWT_LOGGING
        printf ( "PEWT  -  PEWTThreadJoin  -  Attempting to join ALL threads.\n");
    #endif
    
    pthread_t * alignmentThreadBody = writeThread->alignmentThreadBody;
    pthread_t * outputThreadBody    = &(writeThread->outputThreadBody);
    int i;
    
    for (i=0;i<writeThread->numAlignmentThread;i++) {
        if (alignmentThreadBody[i]!=0) {
            #ifdef DEBUG_MICA_PEWT_LOGGING
                printf ( "PEWT  -  PEWTThreadJoin  -  Attempting to join the alignment thread.\n");
            #endif
            pthread_join(alignmentThreadBody[i],NULL);
        }
        alignmentThreadBody[i] = 0;
    }
    
    if ((*outputThreadBody)!=0) {
        #ifdef DEBUG_MICA_PEWT_LOGGING
            printf ( "PEWT  -  PEWTThreadJoin  -  Attempting to join the output thread.\n");
        #endif
        pthread_join((*outputThreadBody),NULL);
    }
    
    writeThread->outputThreadBody = 0;
    
    #ifdef DEBUG_MICA_PEWT_LOGGING
        printf ( "PEWT  -  PEWTThreadJoin  -  ALL threads joined.\n");
    #endif
}

void PEWTPollHealthString(PEWriteThread * writeThread) {
    // healthString is expected to be at least MICA_PE_WRITE_HEALTH_STRING_MAX_LENGTH character long.
    sprintf(writeThread->healthString,"%d/%d FREI",writeThread->freeQueue->count,writeThread->freeQueue->size);
}

char * PEWTGetHealthText(PEWriteThread * writeThread) {
    return writeThread->healthString;
}
