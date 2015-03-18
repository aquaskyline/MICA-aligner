//
//    MIC-MEMControl.c
//
//    mica
//
//    Copyright (C) 2014, HKU
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

#include "MIC-MEMControl.h"

void * MICMemLoadBlock(MICMemBlock * block, int threadId, int threadCount) {

    void * micPtr = NULL;
    int i;
    
    if (block->status == MIC_MEM_STATUS_READY) {
    
        uint8_t * hbytes = (uint8_t*) block->hostPtr;
        uint64_t  hsize = block->atomSize * block->atomCount;
        
        #ifdef MIC_MEM_DEBUG_PRINT_EVENT
            printf("[MIC-MEM] Offloading buffer#%s..\n",block->tag);
            printf("          MIC Pointer Data Size = %llu\n",hsize);
            fflush(0);
        #endif
        
        #pragma offload target(mic : threadId)
        {
            // Allocate memory block of the same size
            micPtr = _mm_malloc(hsize, 64);
            
            #ifdef MIC_MEM_DEBUG_PRINT_EVENT
                printf("          MIC Pointer Address = %llu\n",(unsigned long long) micPtr);
                fflush(0);
            #endif
        }
        
        int cellSize[MIC_MEMCPY_MAX_THREAD];
        if ( threadCount > MIC_MEMCPY_MAX_THREAD ) {
            threadCount = MIC_MEMCPY_MAX_THREAD;
        }
        int evenBlockSize = hsize/threadCount;
        for (i=0;i<threadCount;i++) {
            cellSize[i]=evenBlockSize;
        }
        cellSize[threadCount-1] += hsize - evenBlockSize * threadCount;
        
        #pragma offload target(mic : threadId) \
            in(hbytes:length(hsize) alloc_if(1) free_if(1))
        #pragma omp parallel for private(i) num_threads(threadCount)
        for (i = 0; i < threadCount; i++) {
            int skip = i * evenBlockSize;
            void * t = micPtr + skip;
            uint8_t * h = hbytes + skip;
            memcpy(t,h,cellSize[i]);
        }
    
        #ifdef MIC_MEM_DEBUG_PRINT_EVENT
            printf("DONE\n");
        #endif
        
        block->micPtr = micPtr;
        block->status == MIC_MEM_STATUS_LOADED;
    }
    
    return micPtr;
}

void * MICMemAddBlock(MICMemBlockArray * micMem,
                    void * hostPtr, uint64_t atomSize, uint64_t atomCount, 
                    char * tag) {
    if (micMem->blockUsed >=  micMem->blockMaxSize) {
        printf("[SAFE-GUARD] There isn't enough MIC memory block for your new block %s\n",tag);
        printf("Consider top-up the value MIC_MEM_MAX_BLOCK_SIZE or review all allocated block to reduce count.\n");
        return NULL;
    } else {
    
        int idx = micMem->blockUsed;
        
        MICMemBlock * block = &(micMem->blocks[idx]);
        
        block->atomSize = atomSize;
        block->atomCount = atomCount;
        block->hostPtr = hostPtr;
        block->micPtr = NULL;
        
        strcpy(block->tag,tag);
        block->status = MIC_MEM_STATUS_READY;
        
        void * micPtr = MICMemLoadBlock(block,micMem->threadId,micMem->threadCount);
        block->micPtr = micPtr;
    
        micMem->blockUsed++;
        return block->micPtr;
    }
}

MICMemBlockArray * MICMemCreate(int threadId, int threadCount) {
    MICMemBlockArray * micMem = (MICMemBlockArray*) malloc(sizeof(MICMemBlockArray));
    
    micMem->blocks = (MICMemBlock*) malloc(sizeof(MICMemBlock) * MIC_MEM_MAX_BLOCK_SIZE);
    micMem->blockMaxSize = MIC_MEM_MAX_BLOCK_SIZE;
    micMem->blockUsed = 0;
    micMem->threadId = threadId;
    micMem->threadCount = threadCount;
    
    int i;
    for (i=0;i<micMem->blockMaxSize;i++) {
        micMem->blocks[i].status = MIC_MEM_STATUS_EMPTY;
    }
    
    return micMem;
}

/*
void MICMemLoad(MICMemBlockArray * micMem, int threadId) {

    int i;
    for (i=0;i<micMem->blockMaxSize;i++) {
        
        MICMemBlock * block = &(micMem->blocks[i]);
        MICMemLoadBlock(block,threadId);
        void * micPtr = (*(block->micPtr));
printf("ASD\n");
        
////// TEMP HACK BEGIN
if (i==0) {
#pragma offload target(mic : threadId) \
                inout(block:length(1) alloc_if(1) free_if(1))
{
    printf("nbytes : %llu\n",micPtr);
    int j;
    uint64_t * m = (uint64_t *) micPtr;
    for (j=0;j<100;j++) {
        printf("%u ",m[j]);
    }printf("\n");
    fflush(0);
}
}
////// TEMP HACK END

    }
}

void MICMemUnload(MICMemBlockArray * micMem) {
}
*/

void MICMemFreeBlock(MICMemBlock * block, int threadId) {

    if (block->status == MIC_MEM_STATUS_READY) {
        void * micPtr = block->micPtr;
        #pragma offload target(mic : threadId)
        {
            _mm_free(micPtr);
        }
    }
}

void MICMemFree(MICMemBlockArray * micMem) {

    int i;
    for (i=0;i<micMem->blockMaxSize;i++) {
        MICMemBlock * block = &(micMem->blocks[i]);
        MICMemFreeBlock(block,micMem->threadId);
    }
    
    free(micMem->blocks);
    micMem->blocks = NULL;
    free(micMem);
}

void MICMemDebugPrint(MICMemBlockArray * micMem) {
    int i;
    for (i=0;i<micMem->blockMaxSize;i++) {
        MICMemBlock * block = &(micMem->blocks[i]);
        if (block->status != MIC_MEM_STATUS_EMPTY) {
            printf("BUFFER#%-4d - %s - %d\n",i,block->tag,block->status);
            printf("    SIZE = %llu x %llu = %llu\n",block->atomSize,block->atomCount,block->atomSize*block->atomCount);
        }
    }
}
