//
//    ListConsumer.h
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

#ifndef __LISTCONSUMER_H__
#define __LISTCONSUMER_H__

#include <stdint.h>
#include "utilities/BufferedFileReader.h"

#define LIST_CONS_STATUS_INITIALISED   0
#define LIST_CONS_STATUS_LOADED        1
#define LIST_CONS_STATUS_DEPLETED      2

#define LIST_CONS_TYPE_UNKNOWN         0
#define LIST_CONS_TYPE_LIST_FILE       0
#define LIST_CONS_TYPE_STATIC          1

#define LIST_CONS_MAX_INPUT_IDEN_LEN   1024

#define LIST_CONS_MAX_NUM_OF_ITEM      1024

typedef struct ListConsumer {
    uint8_t status;
    uint8_t type;
    unsigned int listIndex;
    
    // List file
    UTBFRBuffer * listBufferReader;
    int maxFilenameLength;
    int maxFieldLength;
    // Static Definition
    char * readFilenames;
    char * mateFilenames;
    int * lbounds;
    int * ubounds;
    unsigned int listCount;

    // additional information, same as soap3-dp input_options
    char * outputPrefix;
    char * readGroup;
    char * sampleName;
    char * readGrpOption;
} ListConsumer;

typedef struct ListReader {
    char * readFilename;
    char * mateFilename;
    int lbound;
    int ubound;
    
    unsigned int listIndex;
    char * inputIdentifier;
    
    // additional information, same as soap3-dp input_options
    char * outputPrefix;
    char * readGroup;
    char * sampleName;
    char * readGrpOption;
} ListReader;

ListConsumer * LCCreate();
void LCLoadList(ListConsumer * listConsumer, char * listFilename, int maxFilenameLength, int maxFieldLength);
void LCCreateList(ListConsumer * listConsumer, 
                char * readFilenames, char * mateFilenames,
                int * lbounds, int * ubounds, 
                int listSize, int maxFilenameLength, int maxFieldLength);
void LCFree(ListConsumer * listConsumer);
int LCEndOfList(ListConsumer * listConsumer);
int LCGetNextPairEndRead(ListConsumer * listConsumer, char * readFilename, char * mateFilename, int * lbound, int * ubound,
                                        char * outputPrefix, char * readGroup, char * sampleName, char * readGrpOption );

ListReader * LRCreate(ListConsumer * listConsumer);
int LRGetNextPairEndRead(ListConsumer * listConsumer, ListReader * listReader);
void LRFree(ListReader * listReader);

void LCDebugPrint(ListConsumer * listConsumer);
#endif

