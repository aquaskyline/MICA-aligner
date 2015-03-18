//
//    MICA-PE-ReadThread.c
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

#include "MICA-PE-ReadThread.h"

typedef struct PEReadThread {
    // Input
    unsigned char * charMap;
    unsigned int numOfQuery;
    unsigned int maxNameLength;
    unsigned int maxPatternLength;

    // Buffer is a part of ReadThread
    // such that ReadThread can provide input file switching
    // when user supplied a list of input files
    UTBFRBuffer * readQueryBuffer;
    UTBFRBuffer * mateQueryBuffer;
    ListConsumer * listConsumer;
    ListReader * listReader;
    PEInputParam inputParam[LIST_CONS_MAX_NUM_OF_ITEM];
    PEReadOutput readOutput[LIST_CONS_MAX_NUM_OF_ITEM];

    // Output
    PEReadBatch batches[PERT_MAX_NUM_OF_BATCH];
    char nextWrite;  // Which slot to write next
    char nextRead;   // Which slot to read next

    int outputFormat, outputNumSamThreads;

    // Concurrent related
    char batchReady;
    pthread_mutex_t * batchReadyLock;
    pthread_cond_t * batchProduced;
    pthread_cond_t * batchConsumed;
    pthread_mutex_t * quit;
    pthread_t * t;
    
    // Health string
    char healthString[PERT_HEALTH_STRING_MAX_LENGTH];
    
    // All OutputHandlers
    bam_header_t samOutputHeader;
    HSP * PEReadThread_hsp;
} PEReadThread;

static int shouldQuit(pthread_mutex_t * mtx) {

    switch(pthread_mutex_trylock(mtx)) {
    case 0: /* if we got the lock, unlock and return 1 (true) */
      pthread_mutex_unlock(mtx);
      return 1;
    case EBUSY: /* return 0 (false) if the mutex was locked */
      return 0;
    }
    return 1;
}


static void mallocBatch(PEReadThread * t, PEReadBatch * batch) {
    batch->readNames = MEMManMalloc(sizeof(char) * t->numOfQuery * t->maxNameLength, MEMORY_TYPE_CPU);
    batch->mateNames = MEMManMalloc(sizeof(char) * t->numOfQuery * t->maxNameLength, MEMORY_TYPE_CPU);
    batch->readPatterns = MEMManMalloc(sizeof(unsigned char) * t->numOfQuery * t->maxPatternLength, MEMORY_TYPE_CPU);
    batch->matePatterns = MEMManMalloc(sizeof(unsigned char) * t->numOfQuery * t->maxPatternLength, MEMORY_TYPE_CPU);
    batch->readQuality = MEMManMalloc(sizeof(char) * t->numOfQuery * t->maxPatternLength, MEMORY_TYPE_CPU);
    batch->mateQuality = MEMManMalloc(sizeof(char) * t->numOfQuery * t->maxPatternLength, MEMORY_TYPE_CPU);
    batch->readLengths = MEMManMalloc(sizeof(uint16_t) * t->numOfQuery, MEMORY_TYPE_CPU);
    batch->mateLengths = MEMManMalloc(sizeof(uint16_t) * t->numOfQuery, MEMORY_TYPE_CPU);
    batch->readSize = 0;
    batch->mateSize = 0;
    batch->algnFuncReasonFlag = MIC_ALGNMT_FOR_UNDEFINED;

    batch->readUncertaintyNumber = MEMManMalloc(sizeof(int) * t->numOfQuery, MEMORY_TYPE_CPU);
    batch->mateUncertaintyNumber = MEMManMalloc(sizeof(int) * t->numOfQuery, MEMORY_TYPE_CPU);
}

static void freeBatch(PEReadBatch * batch) {
    free(batch->readNames);
    free(batch->mateNames);
    free(batch->readPatterns);
    free(batch->matePatterns);
    free(batch->readLengths);
    free(batch->mateLengths);
    free(batch->readQuality);
    free(batch->mateQuality);

    free(batch->readUncertaintyNumber);
    free(batch->mateUncertaintyNumber);
}

static int fetchNextListItem(PEReadThread * t) {
    ListReader * listReader = t->listReader;
    ListConsumer * listConsumer = t->listConsumer;
    
    int retVal = 0;
    
    char outputFilename[LIST_CONS_MAX_INPUT_IDEN_LEN+15]; // Max file name length + payload enough to store file extension.
    
    //Free previous buffer if loaded
    if (t->readQueryBuffer!=NULL) {
        UTBFRFree(t->readQueryBuffer);
        t->readQueryBuffer = NULL;
    }
    if (t->mateQueryBuffer!=NULL) {
        UTBFRFree(t->mateQueryBuffer);
        t->mateQueryBuffer = NULL;
    }
    
    //If there is an existing item, closing it
    if ( listReader->listIndex>=0 ) {
        t->readOutput[listReader->listIndex].status = PERT_OUTPUT_STATUS_DEPLETED;
    }
    
    //Get next item on list
    while ( !retVal ) {
        retVal = LRGetNextPairEndRead(listConsumer, listReader);

        // List depleted
        if (!retVal) {
            break;
        }

        // If file names are returned. Check them.
        if (!UTBFRIsFile(listReader->readFilename)) {
            printf("[ReadThread] %s cannot be opened!\n",listReader->readFilename);
            retVal = 0;
        } else if (!UTBFRIsFile(listReader->mateFilename)) {
            printf("[ReadThread] %s cannot be opened!\n",listReader->mateFilename);
            retVal = 0;
        } else if (listReader->ubound < listReader->lbound) {
            printf("[ReadThread] %s / %s upper bound smaller than lower bound",
                listReader->readFilename,listReader->mateFilename);
            retVal = 0;
        }
    }
    
    if (retVal) {
        #ifdef MICA_PERT_DEBUG_INPUT_FILE
            printf("\n[MICA_PERT_DEBUG_INPUT_FILE] fetchNextListItem is invoked.\n");
            printf("  %s / %s files will be processed.\n",listReader->readFilename,listReader->mateFilename);
            printf("  %5u / %5u insertion size.\n",listReader->lbound,listReader->ubound);
            LCDebugPrint(t->listConsumer);
        #endif
        
        //Load up the file buffers
        t->readQueryBuffer = UTBFRLoad(listReader->readFilename);
        t->mateQueryBuffer = UTBFRLoad(listReader->mateFilename);
        
        // Create new output for this input set
        if ( listReader->listIndex >= LIST_CONS_MAX_NUM_OF_ITEM ) {
            printf("[SAFE-GUARD] This item on the list has index %u which is greater than/equal to the max number of allowed item on list(%u).\n",
                listReader->listIndex,LIST_CONS_MAX_NUM_OF_ITEM);
            exit(1);
        }
        // reConstruct outputHeader using listReader, which may contain optional arguments for file header
        SAMOutputHeaderConstruct_ListReader(&(t->samOutputHeader),t->PEReadThread_hsp,listReader);
        SRAOutput * outputHandler = SRAAlgnmtOutputConstruct();
        if (t->outputFormat == SRA_OUTPUT_FORMAT_SAM){
            if ( listReader->outputPrefix!=NULL && strlen(listReader->outputPrefix) > 0 ){
                sprintf(outputFilename,"%s.sam",listReader->outputPrefix );
            } else {
                sprintf(outputFilename,"%s_pe.sam",listReader->inputIdentifier);
            }
            // SRAAlgnmtOutputInitSAM(outputHandler,outputFilename,&(t->samOutputHeader));
            SRAAlgnmtOutputInitSAM(outputHandler,outputFilename,&(t->samOutputHeader),listReader->readGroup==NULL? "" : listReader->readGroup );
        }
        else if (t->outputFormat == SRA_OUTPUT_FORMAT_BAM){
            if ( listReader->outputPrefix!=NULL && strlen(listReader->outputPrefix) > 0 ){
                sprintf(outputFilename,"%s.bam",listReader->outputPrefix );
            } else {
                sprintf(outputFilename,"%s_pe.bam",listReader->inputIdentifier);
            }
            //SRAAlgnmtOutputInitBAM(outputHandler,outputFilename,&(t->samOutputHeader), t->outputNumSamThreads );
            SRAAlgnmtOutputInitBAM(outputHandler,outputFilename,&(t->samOutputHeader), t->outputNumSamThreads, listReader->readGroup==NULL? "" : listReader->readGroup );
        }
        else {
            printf("[SAFE-GUARD] [ReadThread] Output Format must be SAM or BAM");
            exit(1);
        }
        // Destroy the samOutputHeader
        SAMOutputHeaderDestruct(&(t->samOutputHeader));
        t->inputParam[listReader->listIndex].InputPEInsertionUpperBound = listReader->ubound;
        t->inputParam[listReader->listIndex].InputPEInsertionLowerBound = listReader->lbound;
        
        t->readOutput[listReader->listIndex].status = PERT_OUTPUT_STATUS_RUNNING;
        t->readOutput[listReader->listIndex].numBatchProduced = 0;
        t->readOutput[listReader->listIndex].numBatchProcessed = 0;
        t->readOutput[listReader->listIndex].outputHandler = outputHandler;
        
        
    } else {
        #ifdef MICA_PERT_DEBUG_INPUT_FILE
            printf("\n[MICA_PERT_DEBUG_INPUT_FILE] fetchNextListItem is invoked.\n");
            printf("  Failed to obtained any item from the list.\n");
            LCDebugPrint(t->listConsumer);
        #endif
    }
    
    return retVal;
}

PEReadThread * PEReadThread_create(
    unsigned char * charMap,
    unsigned int numOfQuery,
    
    unsigned int maxNameLength,
    unsigned int maxPatternLength,
    
    HSP * hsp,

    ListConsumer * pairListConsumer,

    int outputFormat,
    int numSAMThreads) {

    int i;
    
    PEReadThread * ret = MEMManMalloc(sizeof(PEReadThread), MEMORY_TYPE_CPU);

    ret->charMap = charMap;
    ret->numOfQuery = numOfQuery;
    ret->maxNameLength = maxNameLength;
    ret->maxPatternLength = maxPatternLength;
    ret->readQueryBuffer = NULL;
    ret->mateQueryBuffer = NULL;
    ret->listConsumer = pairListConsumer;
    ret->listReader = LRCreate(pairListConsumer);
    ret->nextWrite = 0;
    ret->nextRead = 0;
    ret->outputFormat = outputFormat;
    ret->outputNumSamThreads = numSAMThreads;
    ret->PEReadThread_hsp = hsp;
    for (i=0; i<LIST_CONS_MAX_NUM_OF_ITEM; i++) {
        ret->inputParam[i].InputPEInsertionUpperBound = 0;
        ret->inputParam[i].InputPEInsertionLowerBound = 0;
        
        ret->readOutput[i].status = PERT_OUTPUT_STATUS_INITIALISED;
        ret->readOutput[i].numBatchProduced = 0;
        ret->readOutput[i].numBatchProcessed = 0;
        ret->readOutput[i].outputHandler = NULL;
    }
    
    // Setting up SAM header necessary,
    // SAM header will be generated in fetchNextListItem() for different readfiles
    // SAMOutputHeaderConstruct(&(ret->samOutputHeader),hsp);
    
    // Ready the buffer for the first file
    if (!fetchNextListItem(ret)) {
        printf("[Safe-Guard] Supplied input cannot be loaded.\n");
        exit(1);
    }

    for (i=0; i<PERT_MAX_NUM_OF_BATCH; i++) {
        mallocBatch(ret, ret->batches + i);
    }

    ret->t = MEMManMalloc(sizeof(pthread_t), MEMORY_TYPE_CPU);
    ret->batchReadyLock = MEMManMalloc(sizeof(pthread_mutex_t), MEMORY_TYPE_CPU);
    ret->quit = MEMManMalloc(sizeof(pthread_mutex_t), MEMORY_TYPE_CPU);
    ret->batchProduced = MEMManMalloc(sizeof(pthread_cond_t), MEMORY_TYPE_CPU);
    ret->batchConsumed = MEMManMalloc(sizeof(pthread_cond_t), MEMORY_TYPE_CPU);
    pthread_mutex_init(ret->quit, NULL);
    pthread_mutex_init(ret->batchReadyLock, NULL);
    pthread_cond_init(ret->batchProduced, NULL);
    pthread_cond_init(ret->batchConsumed, NULL);
    pthread_mutex_lock(ret->quit);

    ret->batchReady = 0;

    strcpy(ret->healthString,"INIT");
    
    return ret;
}


static void writeBuffer(PEReadThread * t, char idx) {


#ifdef MICA_PERT_DEBUG_FLOW
    printf("[READ THREAD] Write batch %d\n", idx);
#endif
    PEReadBatch * batch = t->batches + idx;
    ListConsumer * listConsumer = t->listConsumer;
    ListReader * listReader = t->listReader;

    // Clear size information
    batch->readSize = 0;
    batch->mateSize = 0;
    
    //First Attempt is to try the current input file
    if (t->readQueryBuffer != NULL && t->mateQueryBuffer != NULL) {
        batch->readSize = SRAQueryAndNameGetBatchFromFASTAQLength(t->readQueryBuffer,
            t->charMap, t->numOfQuery, t->maxNameLength, batch->readNames,
            t->maxPatternLength, batch->readPatterns, batch->readQuality, batch->readLengths, batch->readUncertaintyNumber);

        batch->mateSize = SRAQueryAndNameGetBatchFromFASTAQLength(t->mateQueryBuffer,
            t->charMap, t->numOfQuery, t->maxNameLength, batch->mateNames,
            t->maxPatternLength, batch->matePatterns, batch->mateQuality, batch->mateLengths, batch->mateUncertaintyNumber);
            
        batch->batchParam = &(t->inputParam[listReader->listIndex]);
        batch->batchOutput = &(t->readOutput[listReader->listIndex]);
        batch->algnFuncReasonFlag = MIC_ALGNMT_FOR_OUTPUT;
    }
    
    // Second Attempt is the next file on the input list
    if (batch->readSize==0 && batch->mateSize == 0) {
    
        if (!LCEndOfList(listConsumer) && fetchNextListItem(t)) {
                
            batch->readSize = SRAQueryAndNameGetBatchFromFASTAQLength(t->readQueryBuffer,
                t->charMap, t->numOfQuery, t->maxNameLength, batch->readNames,
                t->maxPatternLength, batch->readPatterns, batch->readQuality, batch->readLengths, batch->readUncertaintyNumber);

            batch->mateSize = SRAQueryAndNameGetBatchFromFASTAQLength(t->mateQueryBuffer,
                t->charMap, t->numOfQuery, t->maxNameLength, batch->mateNames,
                t->maxPatternLength, batch->matePatterns, batch->mateQuality, batch->mateLengths, batch->mateUncertaintyNumber);
                
            batch->batchParam = &(t->inputParam[listReader->listIndex]);
            batch->batchOutput = &(t->readOutput[listReader->listIndex]);
            batch->algnFuncReasonFlag = MIC_ALGNMT_FOR_OUTPUT;
        }
    }
    
    // If the buffer is valid add 1 to the stats
    if (batch->readSize!=0 || batch->mateSize != 0) {
        t->readOutput[listReader->listIndex].numBatchProduced++;
    }
}

static void * start(void * readThread) {
    PEReadThread * thread = (PEReadThread *) readThread;

    // thread->batchReady should be 0 now
    while (!shouldQuit(thread->quit)) {

        // Read Here
#ifdef MICA_PERT_DEBUG_FLOW
        printf("[READ THREAD] Start Reading\n");
#endif

        writeBuffer(thread, thread->nextWrite);

#ifdef MICA_PERT_DEBUG_FLOW
        printf("[READ THREAD] End Reading\n"); fflush(0);
#endif
        pthread_mutex_lock(thread->batchReadyLock);
        thread->nextWrite = (thread->nextWrite + 1) % PERT_MAX_NUM_OF_BATCH;
        thread->batchReady++;

        pthread_cond_signal(thread->batchProduced);
        while (thread->batchReady >= PERT_MAX_NUM_OF_BATCH && !shouldQuit(thread->quit)) {
            pthread_cond_wait(thread->batchConsumed, thread->batchReadyLock);
        }
        pthread_mutex_unlock(thread->batchReadyLock);
    }

    return NULL;
}


// Wait until a batch is ready for reading
PEReadBatch * PEReadThread_wait(PEReadThread * readThread) {

#ifdef MICA_PERT_DEBUG_FLOW
    printf("[READ THREAD] Start Waiting\n");
#endif

    pthread_mutex_lock(readThread->batchReadyLock);
    while (readThread->batchReady == 0) {
        pthread_cond_wait(readThread->batchProduced, readThread->batchReadyLock);
    }

#ifdef MICA_PERT_DEBUG_FLOW
    printf("[READ THREAD] End Waiting\n");
#endif

#ifdef MICA_PERT_DEBUG_FLOW
    printf("[READ THREAD] Get batch %d\n", readThread->nextRead);
#endif
    PEReadBatch * ret = readThread->batches + readThread->nextRead;
    return ret;
}

// Notify the thread that the batch has been finished reading
// by caller thread. Signal the thread to continue reading the next batch
// which will overwrite the current batch.
void PEReadThread_signal(PEReadThread * readThread) {

#ifdef MICA_PERT_DEBUG_FLOW
        printf("[READ THREAD] Signal\n");
#endif
    readThread->nextRead = (readThread->nextRead + 1) % PERT_MAX_NUM_OF_BATCH;
    readThread->batchReady--;
    pthread_cond_signal(readThread->batchConsumed);
    pthread_mutex_unlock(readThread->batchReadyLock);
}


void PEReadThread_start(PEReadThread * readThread) {
    pthread_create(readThread->t, NULL, start, readThread);
}

// Should be called before the last PEReadThread_signal
void PEReadThread_quit(PEReadThread * readThread) {
    pthread_mutex_unlock(readThread->quit);
    pthread_cond_signal(readThread->batchConsumed);
}

void PEReadThread_ClosingOutput(PEReadOutput * batchOutput) {
    
    batchOutput->numBatchProcessed++;
    
#ifdef MICA_PERT_DEBUG_FLOW
    printf("[READ THREAD] %u/%u batch processed on file with status %d\n",
        batchOutput->numBatchProcessed,batchOutput->numBatchProduced,
        batchOutput->status);
#endif
        
    if ( batchOutput->status == PERT_OUTPUT_STATUS_DEPLETED &&
         batchOutput->numBatchProcessed >= batchOutput->numBatchProduced) {
         
#ifdef MICA_PERT_DEBUG_FLOW
        printf("[READ THREAD] Closing output\n");
#endif

        if ( batchOutput->outputHandler != NULL ) {
            if ( batchOutput->outputHandler->OutFileFormat == SRA_OUTPUT_FORMAT_SAM ) {
                SRAAlgnmtOutputFreeSAM(batchOutput->outputHandler);
            } else if ( batchOutput->outputHandler->OutFileFormat == SRA_OUTPUT_FORMAT_BAM ) {
                SRAAlgnmtOutputFreeBAM(batchOutput->outputHandler);
            }
            SRAAlgnmtOutputFree(batchOutput->outputHandler);
            batchOutput->outputHandler = NULL;
            batchOutput->status = PERT_OUTPUT_STATUS_FREE;
        }
    }

}

void PEReadThread_free(PEReadThread * readThread) {
    int i;
    for (i=0; i<PERT_MAX_NUM_OF_BATCH; i++) {
        freeBatch(readThread->batches + i);
    }
    for (i=0; i<LIST_CONS_MAX_NUM_OF_ITEM; i++) {
        if ( readThread->readOutput[i].outputHandler != NULL ) {
            if (readThread->readOutput[i].outputHandler->OutFileFormat == SRA_OUTPUT_FORMAT_SAM) {
                SRAAlgnmtOutputFreeSAM(readThread->readOutput[i].outputHandler);
            } else if (readThread->readOutput[i].outputHandler->OutFileFormat == SRA_OUTPUT_FORMAT_BAM) {
                SRAAlgnmtOutputFreeBAM(readThread->readOutput[i].outputHandler);
            }
            SRAAlgnmtOutputFree(readThread->readOutput[i].outputHandler);
            readThread->readOutput[i].outputHandler = NULL;
            readThread->readOutput[i].status = PERT_OUTPUT_STATUS_FREE;
        }
    }
    SAMOutputHeaderDestruct(&(readThread->samOutputHeader));
    pthread_mutex_destroy(readThread->quit);
    pthread_mutex_destroy(readThread->batchReadyLock);
    pthread_cond_destroy(readThread->batchConsumed);
    free(readThread->batchConsumed);
    pthread_cond_destroy(readThread->batchProduced);
    free(readThread->batchProduced);
    free(readThread->quit);
    free(readThread->t);
    free(readThread->batchReadyLock);
    if (readThread->readQueryBuffer!=NULL) {
        UTBFRFree(readThread->readQueryBuffer);
        readThread->readQueryBuffer = NULL;
    }
    if (readThread->mateQueryBuffer!=NULL) {
        UTBFRFree(readThread->mateQueryBuffer);
        readThread->mateQueryBuffer = NULL;
    }
    LRFree(readThread->listReader);
    free(readThread);
}

void PEReadThread_join(PEReadThread * readThread) {
    pthread_join(*(readThread->t), NULL);
}

void PEReadThread_PollHealthString(PEReadThread * readThread) {
    // healthString is expected to be at least PERT_HEALTH_STRING_MAX_LENGTH character long.
    sprintf(readThread->healthString,"%d/%d REIF",readThread->batchReady,PERT_MAX_NUM_OF_BATCH);
}

char * PEReadThread_GetHealthText(PEReadThread * readThread) {
    return readThread->healthString;
}
