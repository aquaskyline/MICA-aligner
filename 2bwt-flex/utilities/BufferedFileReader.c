//
//    BufferedFileReader.c
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

#include "BufferedFileReader.h"

int UTBFRIsFile(const char * filename) {
    int retVal = 1;
    FILE * fp = fopen(filename, "r");
    if (fp == NULL) {
        retVal = 0;
    } else {
        fclose(fp);
    }
    return retVal;
}

UTBFRBuffer * UTBFRLoad(const char * filename) {
    UTBFRBuffer * result = (UTBFRBuffer*) malloc (sizeof(UTBFRBuffer));
    result->queryFile = gzopen(filename, "r");
    if (result->queryFile==NULL) { 
        fprintf(stderr,"[UTBFRLoad] Cannot open file %s\n",filename); 
        exit(1);
    }
    UTBFRGetBufferFromFile(result);
    return result;
}

void UTBFRResetBuffer(UTBFRBuffer * queryBuffer) {
    gzseek(queryBuffer->queryFile, 0, SEEK_SET);
    UTBFRGetBufferFromFile(queryBuffer);
}

void UTBFRGetBufferFromFile(UTBFRBuffer * queryBuffer) {
    queryBuffer->bufferSize = gzread(queryBuffer->queryFile,&(queryBuffer->queryFileBuffer),INPUT_BUFFER_SIZE);
    queryBuffer->bufferPtr = 0;
}

char UTBFRGetCharacter(UTBFRBuffer * queryBuffer) {
    if (queryBuffer->bufferPtr >= queryBuffer->bufferSize) {
        UTBFRGetBufferFromFile(queryBuffer);
    }
    char character = queryBuffer->queryFileBuffer[queryBuffer->bufferPtr];
    queryBuffer->bufferPtr++;
    return character;
}

int UTBFREOF(UTBFRBuffer * queryBuffer) {
    if (queryBuffer->bufferPtr>=queryBuffer->bufferSize && gzeof(queryBuffer->queryFile)) {
        return 1;
    }
    return 0;
}

void UTBFRFree(UTBFRBuffer * queryBuffer) {
    if (queryBuffer==NULL) return;
    if (queryBuffer->queryFile != NULL) {gzclose(queryBuffer->queryFile);}
    free(queryBuffer);
}
