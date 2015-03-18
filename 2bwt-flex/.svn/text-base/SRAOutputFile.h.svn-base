//
//    SRAOutputFile.h
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

#ifndef __SRAOUTPUTFILE_H__
#define __SRAOUTPUTFILE_H__

#include "SRAOccCollector.h"

#include "OutputAnyFormat.h"
#include "OutputBinaryFormat.h"
#include "OutputSAMFormat.h"

typedef struct SRAOutput {
    char OutFileName[MAX_FILENAME_LEN*2+1];
    uint8_t OutFileFormat;
    FILE * OutFilePtr;
    char ReadGroup[MAX_FIELD_LEN+1];
    //Parameter for SAM
    samfile_t * SAMOutFilePtr;
    uint8_t * SAMAuxDataBlock;
    unsigned long long SAMAuxDataBlockSize;
    
    //Parameter for new OCC
    SRAOCCCollector * occCollector;
} SRAOutput;


SRAOutput * SRAAlgnmtOutputConstruct();
void SRAAlgnmtOutputFree(SRAOutput * sraOutput);

// INITIALISER 
// Initialiser prepare the structure for different type of output format
// It does not allow overriding the same structure already occupied by 
// another format
void SRAAlgnmtOutputInitNone(SRAOutput * sraOutput);
void SRAAlgnmtOutputInitDiscard(SRAOutput * sraOutput);
void SRAAlgnmtOutputInitPlain(SRAOutput * sraOutput, char * OutFileName);
void SRAAlgnmtOutputInitSAMStore(SRAOutput * sraOutput);
void SRAAlgnmtOutputInitSAM(SRAOutput * sraOutput, char * OutFileName, bam_header_t * samOutputHeader, char * readGroup);
void SRAAlgnmtOutputInitBAM(SRAOutput * sraOutput, char * OutFileName, bam_header_t * samOutputHeader, int numthreads,  char * readGroup);
// FREE
// Free-er clean up the structure and re-initialise it after cleaning up
// format-specific data. It however does not flush the output left in
// the occ-collector. It's the aligner responsibility.
void SRAAlgnmtOutputFreeNone(SRAOutput * sraOutput);
void SRAAlgnmtOutputFreeDiscard(SRAOutput * sraOutput);
void SRAAlgnmtOutputFreePlain(SRAOutput * sraOutput);
void SRAAlgnmtOutputFreeSAMStore(SRAOutput * sraOutput);
void SRAAlgnmtOutputFreeSAM(SRAOutput * sraOutput);
void SRAAlgnmtOutputFreeBAM(SRAOutput * sraOutput);

#endif

