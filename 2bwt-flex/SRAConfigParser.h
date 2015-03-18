//
//    SRAConfigParser.h
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
    File identifier: file:///usr/local/svn/project/2bwt-flex/latest/SRAConfigParser.h
    2BWT-Flex Library
    SRAConfigParser is a CFG parser that reads in a presumed format CFG

    Modification History

    Date   : 16th December 2012
    Author : Edward MK Wu
    Change : New file.
    
*/
/////////////////////////////////////////////////////

#ifndef __SRA_CONFIG_PARSE_H__
#define __SRA_CONFIG_PARSE_H__

#define SRA_CFGPARSE_MAX_LEVEL        6
#define SRA_CFGPARSE_MAX_TOKEN_LENGTH 256
#define SRA_CFGPARSE_MAX_TOKEN_COUNT  20480
#define SRA_CFGPARSE_MAX_TOKEN_PER_LEVEL SRA_CFGPARSE_MAX_TOKEN_COUNT / SRA_CFGPARSE_MAX_LEVEL

#define SRA_CFGPARSE_TYPE_NONE        0
#define SRA_CFGPARSE_TYPE_TOKEN       1
#define SRA_CFGPARSE_TYPE_VALUE       2
#define SRA_CFGPARSE_TYPE_BLOCK_ENTRY 3
#define SRA_CFGPARSE_TYPE_BLOCK_EXIT  4

#define SRA_CFGPARSE_READER_UNREAD    0
#define SRA_CFGPARSE_READER_READ      1

#include "SRAArguments.h"
#include "utilities/BufferedFileReader.h"

typedef struct SRACFGNode {
    int8_t type;
    int16_t index;
} SRACFGNode;

typedef struct SRACFGLevel {
    int16_t nodeCount;
    SRACFGNode nodes[SRA_CFGPARSE_MAX_TOKEN_PER_LEVEL];
} SRACFGLevel;

typedef struct SRACFG {
    SRACFGLevel levels[SRA_CFGPARSE_MAX_LEVEL];
    char tokens[SRA_CFGPARSE_MAX_TOKEN_COUNT][SRA_CFGPARSE_MAX_TOKEN_LENGTH];
} SRACFG;


typedef struct SRACFGReaderNode {
    int8_t readerFlag;
} SRACFGReaderNode;

typedef struct SRACFGReaderLevel {
    int16_t readerCount;
    SRACFGReaderNode nodes[SRA_CFGPARSE_MAX_TOKEN_PER_LEVEL];
} SRACFGReaderLevel;

typedef struct SRACFGReader {
    int16_t readerLevelIdx;
    SRACFGReaderLevel levels[SRA_CFGPARSE_MAX_LEVEL];
} SRACFGReader;




SRACFG * SRACPLoad ( UTBFRBuffer * cfgBuffer);
void SRACPFree ( SRACFG * cfg );
void SRACPLoadModel ( SRACFG * cfg, SRAModel * SRAMismatchModel,
                    unsigned int ReadLength,
                    SRASetting * aSettings, SRAIndex * aIndex, int modelId );

void SRACPDebugPrint ( SRACFG * cfg );
int8_t _SRACPReaderEnter ( SRACFG * cfg, SRACFGReader * cfgReader, char * blockName );
void _SRACPReaderExit ( SRACFG * cfg, SRACFGReader * cfgReader );
int8_t _SRACPReaderGetValue ( SRACFG * cfg, SRACFGReader * cfgReader, char * tokenName, char * tokenValue, char * defaultValue );

#endif
