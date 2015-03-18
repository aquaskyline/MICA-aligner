//
//    ListConsumer.c
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

#include "ListConsumer.h"
static int _LCGetFileName(UTBFRBuffer * listBufferReader, char * filename, int maxLength);
static int _LCGetInteger(UTBFRBuffer * listBufferReader, int * num);
static int _LCGetOptionalField(UTBFRBuffer * listBufferReader, char * optionalfield, int maxLength);
static int _LCSkipRest(UTBFRBuffer * listBufferReader);

void _LRBuildInputIdentifier(ListReader * listReader);

ListConsumer * LCCreate() {
    ListConsumer * listConsumer = (ListConsumer *) malloc ( sizeof(ListConsumer) );
    listConsumer->listIndex = 0;
    listConsumer->listBufferReader = NULL;
    listConsumer->maxFilenameLength = 0;
    listConsumer->maxFieldLength = 0;
    listConsumer->status = LIST_CONS_STATUS_INITIALISED;
    listConsumer->type = LIST_CONS_TYPE_UNKNOWN;
    return listConsumer;
}

void LCLoadList(ListConsumer * listConsumer, char * listFilename, int maxFilenameLength, int maxFieldLength) {

    if ( listConsumer->status != LIST_CONS_STATUS_INITIALISED ) {
        printf("[SAFE-GUARD] LCLoadList found unexpected status from list consumer slot.\n");
    } else {
        listConsumer->listCount = 0;
    
        listConsumer->readFilenames = NULL;
        listConsumer->mateFilenames = NULL;
        listConsumer->lbounds = NULL;
        listConsumer->ubounds = NULL;
        
        listConsumer->listBufferReader = UTBFRLoad(listFilename);
        listConsumer->maxFilenameLength = maxFilenameLength;
        listConsumer->maxFieldLength = maxFieldLength;
        
        listConsumer->status = LIST_CONS_STATUS_LOADED;
        listConsumer->type = LIST_CONS_TYPE_LIST_FILE;
        
        // additional information, same as soap3-dp input_options
        listConsumer->outputPrefix= NULL;
        listConsumer->readGroup = NULL;
        listConsumer->sampleName = NULL;
        listConsumer->readGrpOption = NULL;
    }
}

// readFilenames is a 2-D matrix [listSize x maxFilenameLength] ; so as mateFilenames
void LCCreateList(ListConsumer * listConsumer, 
                char * readFilenames, char * mateFilenames,
                int * lbounds, int * ubounds, 
                int listSize, int maxFilenameLength, int maxFieldLength) {
    if ( listConsumer->status != LIST_CONS_STATUS_INITIALISED ) {
        printf("[SAFE-GUARD] LCCreateList found unexpected status from list consumer slot.\n");
    } else {
        listConsumer->listCount = listSize;
        
        listConsumer->readFilenames = (char*) malloc( sizeof(char) * maxFilenameLength * listSize );
        listConsumer->mateFilenames = (char*) malloc( sizeof(char) * maxFilenameLength * listSize );
        listConsumer->lbounds = (int*) malloc( sizeof(int) * listSize );
        listConsumer->ubounds = (int*) malloc( sizeof(int) * listSize );
        
        memcpy(listConsumer->readFilenames,readFilenames,sizeof(char) * maxFilenameLength * listSize);
        memcpy(listConsumer->mateFilenames,mateFilenames,sizeof(char) * maxFilenameLength * listSize);
        memcpy(listConsumer->lbounds,lbounds,sizeof(int) * listSize);
        memcpy(listConsumer->ubounds,ubounds,sizeof(int) * listSize);
        
        listConsumer->listBufferReader = NULL;
        listConsumer->maxFilenameLength = maxFilenameLength;
        listConsumer->maxFieldLength = maxFieldLength;
        listConsumer->status = LIST_CONS_STATUS_LOADED;
        listConsumer->type = LIST_CONS_TYPE_STATIC;
        
        // additional information, same as soap3-dp input_options
        listConsumer->outputPrefix = (char*) malloc( sizeof(char) * maxFieldLength * listSize );
        listConsumer->readGroup = (char*) malloc( sizeof(char) * maxFieldLength * listSize );
        listConsumer->sampleName = (char*) malloc( sizeof(char) * maxFieldLength * listSize );
        listConsumer->readGrpOption = (char*) malloc( sizeof(char) * maxFieldLength * listSize );
    }
}

int LCEndOfList(ListConsumer * listConsumer) {
    return (listConsumer->status == LIST_CONS_STATUS_DEPLETED);
}

int LCGetNextPairEndRead(ListConsumer * listConsumer, char * readFilename, char * mateFilename, int * lbound, int * ubound, 
                         char * outputPrefix, char * readGroup, char * sampleName, char * readGrpOption ) {

    int lineOk = 1;
    
    if (listConsumer->status == LIST_CONS_STATUS_DEPLETED) {
    
        printf("[SAFE-GUARD] List consumer does not have a valid status. This implies\n");
        printf("List consumer could have been corrupted or invalid list file.\n");
        lineOk = 0;
        
    } else if ( listConsumer->listIndex >= LIST_CONS_MAX_NUM_OF_ITEM ) {
        
        printf("[SAFE-GUARD] List exceeds maximum number of items (%u).\n",LIST_CONS_MAX_NUM_OF_ITEM);
        printf("Review compile time definition of LIST_CONS_MAX_NUM_OF_ITEM.\n");
        lineOk = 0;
        
    } else {
        if ( listConsumer->type == LIST_CONS_TYPE_LIST_FILE ) {
            // Assuming the format is as below:
            // <Read file> <Mate file> <insertion Lower Bound> <Upper Bound> ( optional: <Output Prefix> <Read Group ID> <Sample Name> <Read Group Option> )
            lineOk = lineOk && _LCGetFileName(listConsumer->listBufferReader,readFilename,listConsumer->maxFilenameLength);
            lineOk = lineOk && _LCGetFileName(listConsumer->listBufferReader,mateFilename,listConsumer->maxFilenameLength);
            lineOk = lineOk && _LCGetInteger(listConsumer->listBufferReader,lbound);
            lineOk = lineOk && _LCGetInteger(listConsumer->listBufferReader,ubound);
            
            
            //int isEndOfLine = 0;
            //_LCGetOptionalField(listConsumer->listBufferReader, char * optionalfield, MAX_FIELD_LENGTH, int * isEndOfLine) 
            
            _LCGetOptionalField(listConsumer->listBufferReader,  outputPrefix, listConsumer->maxFieldLength );
            _LCGetOptionalField(listConsumer->listBufferReader,  readGroup, listConsumer->maxFieldLength );
            _LCGetOptionalField(listConsumer->listBufferReader,  sampleName, listConsumer->maxFieldLength );
            _LCGetOptionalField(listConsumer->listBufferReader,  readGrpOption, listConsumer->maxFieldLength );
            lineOk = lineOk && _LCSkipRest(listConsumer->listBufferReader);
        } else if ( listConsumer->type == LIST_CONS_TYPE_STATIC ) {
            lineOk = ( listConsumer->listIndex < listConsumer->listCount );
            if (lineOk) {
                strcpy( readFilename , listConsumer->readFilenames + listConsumer->listIndex * listConsumer->maxFilenameLength );
                strcpy( mateFilename , listConsumer->mateFilenames + listConsumer->listIndex * listConsumer->maxFilenameLength );
                (*lbound) = listConsumer->lbounds[listConsumer->listIndex];
                (*ubound) = listConsumer->ubounds[listConsumer->listIndex];
                if (outputPrefix!=NULL){
                    outputPrefix[0]='\0';
                }
                if (readGroup!=NULL){
                    readGroup[0]='\0';
                }
                if (sampleName!=NULL){
                    sampleName[0]='\0';
                }
                if (readGrpOption!=NULL){
                    readGrpOption[0]='\0';
                }
            }
        }
        
        listConsumer->listIndex++;

        if (!lineOk) {
            listConsumer->status = LIST_CONS_STATUS_DEPLETED;
        } else {
            if (readFilename[0]=='\0'||mateFilename[0]=='\0') {
                printf("[SAFE-GUARD] LCGetNextPairEndRead returns empty file names ");
                printf("yet the function returns normally.\n");
            }
        }
    }
    
    return lineOk;
}

void LCFree(ListConsumer * listConsumer) {
    if ( listConsumer->readFilenames != NULL ) {
        free(listConsumer->readFilenames);
    }
    if ( listConsumer->mateFilenames != NULL ) {
        free(listConsumer->mateFilenames);
    }
    if ( listConsumer->lbounds != NULL ) {
        free(listConsumer->lbounds);
    }
    if ( listConsumer->ubounds != NULL ) {
        free(listConsumer->ubounds);
    }
    if ( listConsumer->listBufferReader != NULL ) {
        UTBFRFree(listConsumer->listBufferReader);
    }
    if ( listConsumer->readGroup!= NULL ) {
        free(listConsumer->readGroup);
    }
    if ( listConsumer->sampleName!= NULL ) {
        free(listConsumer->sampleName);
    }
    if ( listConsumer->outputPrefix != NULL ) {
        free(listConsumer->outputPrefix);
    }
    if ( listConsumer->readGrpOption!= NULL ) {
        free(listConsumer->readGrpOption);
    }
    free(listConsumer);
}

ListReader * LRCreate(ListConsumer * listConsumer) {
    if ( listConsumer->status != LIST_CONS_STATUS_LOADED ) {
        printf("[SAFE-GUARD] LRCreate found unexpected status from list consumer slot.\n");
        return NULL;
    } else {
        ListReader * listReader = (ListReader*) malloc (sizeof(ListReader));
        listReader->readFilename = (char*) malloc( sizeof(char) * listConsumer->maxFilenameLength );
        listReader->mateFilename = (char*) malloc( sizeof(char) * listConsumer->maxFilenameLength );
        listReader->inputIdentifier = (char*) malloc( sizeof(char) * (LIST_CONS_MAX_INPUT_IDEN_LEN + 1) );
        listReader->lbound = 0;
        listReader->ubound = 0;
        listReader->listIndex = 0;
        listReader ->outputPrefix = (char*) malloc( sizeof(char) * listConsumer->maxFieldLength );
        listReader ->readGroup= (char*) malloc( sizeof(char) * listConsumer->maxFieldLength );
        listReader ->sampleName= (char*) malloc( sizeof(char) * listConsumer->maxFieldLength );
        listReader ->readGrpOption= (char*) malloc( sizeof(char) * listConsumer->maxFieldLength );
        return listReader;
    }
}

int LRGetNextPairEndRead(ListConsumer * listConsumer, ListReader * listReader) {
    listReader->listIndex++;
    
    int retVal = LCGetNextPairEndRead(listConsumer,
                                listReader->readFilename,
                                listReader->mateFilename,
                                &(listReader->lbound),
                                &(listReader->ubound),
                                listReader->outputPrefix,
                                listReader->readGroup,
                                listReader->sampleName,
                                listReader->readGrpOption );
                                
    if (retVal) {
        _LRBuildInputIdentifier(listReader);
    } else {
        sprintf(listReader->inputIdentifier,"SafeGuard.Error");
    }
    
    return retVal;
}

void LRFree(ListReader * listReader) {
    if (listReader->readFilename != NULL) {
        free(listReader->readFilename);
        listReader->readFilename = NULL;
    }
    if (listReader->mateFilename != NULL) {
        free(listReader->mateFilename);
        listReader->mateFilename = NULL;
    }
    if (listReader->inputIdentifier != NULL) {
        free(listReader->inputIdentifier);
        listReader->inputIdentifier = NULL;
    }
    if ( listReader->readGroup!= NULL ) {
        free(listReader->readGroup);
        listReader->readGroup = NULL;
    }
    if ( listReader->sampleName!= NULL ) {
        free(listReader->sampleName);
        listReader->sampleName = NULL;
    }
    if ( listReader->outputPrefix!= NULL ) {
        free(listReader->outputPrefix);
        listReader->outputPrefix = NULL;
    }
    if ( listReader->readGrpOption!= NULL ) {
        free(listReader->readGrpOption);
        listReader->readGrpOption = NULL;
    }
    free(listReader);
}

/////////////////////////////////////////////////////////////////////////////
// Internal Functions
/////////////////////////////////////////////////////////////////////////////

void _LRBuildInputIdentifier(ListReader * listReader) {
    char idBuffer[LIST_CONS_MAX_INPUT_IDEN_LEN+1];
    int charLeft = LIST_CONS_MAX_INPUT_IDEN_LEN;
    int charIdx = 0;
    
    // remove the filePath before the input filenames to get the output filename
    char * readFilenameShort = strrchr ( listReader->readFilename , '/' ); 
    char * mateFilenameShort = strrchr ( listReader->mateFilename , '/' ); 
    if ( readFilenameShort == NULL ) {
        readFilenameShort = listReader->readFilename;
    } else {
        readFilenameShort++;
    }
    if ( mateFilenameShort == NULL ) {
        mateFilenameShort = listReader->mateFilename;
    } else {
        mateFilenameShort++;
    }
    
    int len1 = strlen(readFilenameShort);
    int len2 = strlen(mateFilenameShort);
    
    int strcpyLen = 0;
    
    // Part one of the string
    if (charLeft < len1) {
        strcpyLen = charLeft;
    } else {
        strcpyLen = len1;
    }
    
    if (strcpyLen > 0) {
        memcpy(idBuffer+charIdx,readFilenameShort,strcpyLen);
        charIdx  += strcpyLen;
        charLeft -= strcpyLen;
    }
    
    // Underscore
    if (charLeft < 1) {
        strcpyLen = charLeft;
    } else {
        strcpyLen = 1;
    }
    
    if (strcpyLen > 0) {
        idBuffer[charIdx] = '_';
        charIdx += strcpyLen;
        charLeft -= strcpyLen;
    }

    // Part two of the string
    if (charLeft < len2) {
        strcpyLen = charLeft;
    } else {
        strcpyLen = len2;
    }
    
    if (strcpyLen > 0) {
        memcpy(idBuffer+charIdx,mateFilenameShort,strcpyLen);
        charIdx  += strcpyLen;
        charLeft -= strcpyLen;
    }
    
    idBuffer[charIdx++] = '\0';
    
    memcpy(listReader->inputIdentifier,idBuffer,charIdx);
}

// Read 1 filename from the buffered reader. Caller responsibility to 
// initialise the filename memory
static int _LCGetFileName(UTBFRBuffer * listBufferReader, char * filename, int maxLength) {
    char c = UTBFRGetCharacter(listBufferReader);

    // Skip any padding
    while (!UTBFREOF(listBufferReader) && ( c == '\n' || c == '\t' || c == ' ' ) ) {
        c = UTBFRGetCharacter(listBufferReader);
    }
    
    int i = 0;
    int retVal = 1;

    while (!UTBFREOF(listBufferReader) && c != '\n' && c != '\t' && c != ' ') {
        if (i<maxLength) {
            filename[i++] = c;
        } else {
            retVal = 0;
        }
        c = UTBFRGetCharacter(listBufferReader);
    }
    filename[i] = '\0';
    return (retVal && i>0);
}

// Read 1 integer from the buffered reader.
static int _LCGetInteger(UTBFRBuffer * listBufferReader, int * num) {
    char c = UTBFRGetCharacter(listBufferReader);
    int sum = 0;
    int retVal = 1;
    int i = 0;
    
    // Skip any padding
    while (!UTBFREOF(listBufferReader) && ( c == '\n' || c == '\t' || c == ' ' ) ) {
        c = UTBFRGetCharacter(listBufferReader);
    }
    
    if ( c < '0' || c > '9' ) {
        retVal = 0;
    } else {
        while (!UTBFREOF(listBufferReader) && c >= '0' && c <= '9') {
            sum = ( sum * 10 ) + ( c - '0' );
            i++;
            c = UTBFRGetCharacter(listBufferReader);
        }
        
        (*num) = sum;
    }
    
    return (retVal && i>0);
}

// Read 1 optional field from the buffered reader. Caller responsibility to 
// initialise the field memory. with isEndOfLine indicator
static int _LCGetOptionalField(UTBFRBuffer * listBufferReader, char * optionalfield, int maxLength) {
    
    if (optionalfield == NULL || UTBFREOF(listBufferReader) )
        return 0;
    listBufferReader->bufferPtr--;
    char c = UTBFRGetCharacter(listBufferReader);
    
    // if c == '\n' then the optional field does not exist
    if ( UTBFREOF(listBufferReader) || c == '\n' ){
        optionalfield[0]='\0';
        //isEndOfLine = 1 ;
        return 0;
    }
    
    // Skip any padding
    while (!UTBFREOF(listBufferReader) && ( c == '\n' || c == '\t' || c == ' ' ) ) {
        c = UTBFRGetCharacter(listBufferReader);
    }
    
    int i = 0;
    int retVal = 1;
    int not_preskipchar = 1;
    while (!UTBFREOF(listBufferReader) && c != '\n' && c != '\t' && c != ' ') {
        if (i<maxLength) {
            if (not_preskipchar){
                optionalfield[i++] = c;
                if (c=='\\') { not_preskipchar = 0;}
            }else{
                if (c=='t') { optionalfield[i-1] = '\t';}
                else { optionalfield[i++] = c; }
                not_preskipchar = 1;
            }
            
        } else {
            retVal = 0;
        }
        c = UTBFRGetCharacter(listBufferReader);
    }
    optionalfield[i] = '\0';
    
    //isEndOfLine = ( c == '\n' ? 1 : isEndOfLine );
    return (retVal && i>0);
}


// Skip everything until new-line
static int _LCSkipRest(UTBFRBuffer * listBufferReader) {
    int retVal = 1;

    listBufferReader->bufferPtr--;
    char c = UTBFRGetCharacter(listBufferReader);
    
    while (!UTBFREOF(listBufferReader) && c != '\n') {
        c = UTBFRGetCharacter(listBufferReader);
    }
    
    return retVal;
}

void LCDebugPrint(ListConsumer * listConsumer) {
    if ( listConsumer->status != LIST_CONS_STATUS_LOADED ) {
        printf("[SAFE-GUARD] LCDebugPrint found unexpected status from list consumer slot.\n");
    } else {
        int i;
        unsigned int nameOffset = 0;
        for (i=0;i<listConsumer->listCount;i++) {
            printf("%s %s %u %u",listConsumer->readFilenames+nameOffset,
                                   listConsumer->readFilenames+nameOffset,
                                   listConsumer->lbounds[i],
                                   listConsumer->ubounds[i]);
            if (i==listConsumer->listIndex){ 
                printf("  <");
            }
            printf("\n");
            nameOffset += listConsumer->maxFilenameLength;
        }
    }
}
