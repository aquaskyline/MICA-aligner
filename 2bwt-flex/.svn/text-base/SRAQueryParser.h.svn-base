//
//    SRAQueryParser.h
//
//    SOAP2 / 2BWT
//
//    Copyright (C) 2012, HKU
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
/////////////////////////////////////////////////////
/*
    File identifier: file:///usr/local/svn/project/2bwt-gpu/latest/SRAQuery.h
    GPU-BWT Library for 2BWT-Flex
    SRAQuery is a FASTA/FASTQ parser that buffer the read file loading.
    Last updated: 2011/03/05

    Modification History

    Date   : 5th March 2011
    Author : Edward MK Wu
    Change : New file.
    
    Date   : 7th March 2011
    Author : Edward MK Wu
    Change : New function SRAFetchReadFromPacked to take
             SRA read from a packed buffer.

    Date   : 19th June 2011
    Author : Edward MK Wu
    Change : Packaging 2BWT library as a separate product.
             Thus, changing all references to 2bwt lib to subdirectory.

*/
/////////////////////////////////////////////////////

#ifndef __SRAQUERY_H__
#define __SRAQUERY_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "2bwt-lib/HSP.h"
#include "utilities/BufferedFileReader.h"

#define SANGER_QUALITY_SCALE 33
#define SOLEXA_QUALITY_SCALE 64

// ------------------------------------------------------
// SRAQUERY CTC#1 SRAQUERY_READ_NAME_STOP_AT_SPACE
// Configurable flag to control whether the read name
// parser should stop storing characters after the first
// space on the read name. Define to stop.
#define SRAQUERY_READ_NAME_STOP_AT_SPACE
// ------------------------------------------------------

char _SRAFetchReadName(UTBFRBuffer * queryBuffer, char * queryName, unsigned int maxReadNameLength);
char _SRAFetchReadSequenceNoMap(UTBFRBuffer * queryBuffer,
                                char * queryPattern,
                                unsigned int maxReadLength,
                                int * readLength);

unsigned int SRAQueryGetBatchFromFASTA(UTBFRBuffer * queryBuffer,
                                unsigned char * charMap,
                                unsigned int numOfQuery,
                                unsigned int maxReadLength,
                                unsigned char * queryPattern);
                                
unsigned int SRAQueryAndNameGetBatchFromFASTA(UTBFRBuffer * queryBuffer,
                                            unsigned char * charMap,
                                            unsigned int numOfQuery,
                                            unsigned int maxReadNameLength,
                                            char * queryName,
                                            unsigned int maxReadLength,
                                            unsigned char * queryPattern,
                                            char terminalChar);
                                            
unsigned int SRAQueryAndNameGetBatchFromFASTAQ(UTBFRBuffer * queryBuffer,
                                            unsigned char * charMap,
                                            unsigned int numOfQuery,
                                            unsigned int maxReadNameLength,
                                            char * queryName,
                                            unsigned int maxReadLength,
                                            unsigned char * queryPattern,
                                            char terminalChar);
                                            
unsigned int SRAQueryAndNameGetBatchFromFASTAQLength(UTBFRBuffer * queryBuffer,
                                            unsigned char * charMap,
                                            unsigned int numOfQuery,
                                            unsigned int maxReadNameLength,
                                            char * queryName,
                                            unsigned int maxReadLength,
                                            unsigned char * queryPattern,
                                            char * queryQuality,
                                            uint16_t * queryLength,
                                            int * uncertaintyNumber);

unsigned int SRAQueryAndNameGetBatchFromFASTAQNoMap(UTBFRBuffer * queryBuffer,
                                            unsigned int numOfQuery,
                                            unsigned int maxReadNameLength,
                                            char * queryName,
                                            unsigned int maxReadLength,
                                            char * queryPattern);
                                            
/*void SRAFetchReadFromPacked(unsigned int * packedQueries,
                            unsigned int readId, unsigned int readLength,
                            char * output);
*/
                            
#endif
