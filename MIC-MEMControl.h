//
//    MIC-MEMControl.h
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

#ifndef __MIC_MEMCONTROL_H__
#define __MIC_MEMCONTROL_H__

///////////////////////////////////////
// DEBUG FLAGS
///////////////////////////////////////
// Define the below debug flag to output
// debug messages when there is operations on
// MIC Memory Controller.
//#define MIC_MEM_DEBUG_PRINT_EVENT

///////////////////////////////////////
// Including Standard Libraries
///////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#define MIC_MEM_MAX_BLOCK_SIZE  256

#define MIC_MEM_MAX_TAG_LENGTH  256

#define MIC_MEM_STATUS_EMPTY    0
#define MIC_MEM_STATUS_READY    1
#define MIC_MEM_STATUS_LOADED   2

#define MIC_MEMCPY_MAX_THREAD   224

typedef struct MICMemBlock {
    uint64_t    atomSize; // in bytes
    uint64_t    atomCount;
    void *      hostPtr;
    void *      micPtr;
    char        tag[MIC_MEM_MAX_TAG_LENGTH];
    uint8_t     status;
} MICMemBlock;

typedef struct MICMemBlockArray {
    MICMemBlock * blocks;
    int blockUsed;
    int blockMaxSize;
    int threadId;
    int threadCount;
} MICMemBlockArray;

MICMemBlockArray * MICMemCreate(int threadId,int threadCount);
void * MICMemLoadBlock(MICMemBlock * block, int threadId, int threadCount);
void * MICMemAddBlock(MICMemBlockArray * micMem,
                    void * hostPtr, uint64_t atomSize, uint64_t atomCount, 
                    char * tag);
void MICMemFree(MICMemBlockArray * micMem);
void MICMemFreeBlock(MICMemBlock * block, int threadId);
void MICMemFree(MICMemBlockArray * micMem) ;

#endif

