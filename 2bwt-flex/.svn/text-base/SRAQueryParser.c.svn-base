//
//    SRAQueryParser.c
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
    File identifier: file:///usr/local/svn/project/2bwt-gpu/latest/SRAQuery.c
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
             
    Date   : 2nd October 2011
    Author : Edward MK Wu
    Change : Add support to FASTQ file.
             
*/
/////////////////////////////////////////////////////

#include "SRAQueryParser.h"

char _SRAFetchReadName(UTBFRBuffer * queryBuffer, char * queryName, unsigned int maxReadNameLength) {
    char queryChar = UTBFRGetCharacter(queryBuffer);
    int headerCharId=0;
    #ifdef SRAQUERY_READ_NAME_STOP_AT_SPACE
    int reachedSpace=0;
    #endif
    while (!UTBFREOF(queryBuffer) && queryChar != '\n') {
        #ifdef SRAQUERY_READ_NAME_STOP_AT_SPACE
            if (queryChar==' ') {
                reachedSpace = 1;
            }
            if (reachedSpace==0)
        #endif
        if (headerCharId<maxReadNameLength)
            queryName[headerCharId++]=queryChar;
        queryChar = UTBFRGetCharacter(queryBuffer);
    }
    queryName[headerCharId++]='\0';
    return queryChar;
}

inline char _SRAFetchReadSequence(UTBFRBuffer * queryBuffer, 
                                unsigned char * charMap,
                                unsigned char * queryPattern, 
                                unsigned int maxReadLength, 
                                char terminalChar,
                                int * readLength) {

    char queryChar = UTBFRGetCharacter(queryBuffer);
    int nucleoId=0;
    int invalidChar=0;
    while (!UTBFREOF(queryBuffer) && queryChar != '>' && queryChar != '@' && queryChar != '+') {
        if (queryChar != '\n') {
            if (queryChar >= 'a' && queryChar <= 'z') {
                queryChar += (unsigned char)('A' - 'a');
            }
            if (nucleoId<maxReadLength) {
                if (charMap[queryChar] >= ALPHABET_SIZE) invalidChar=1;
                queryPattern[nucleoId++]=charMap[queryChar];
            }
        }
        queryChar = UTBFRGetCharacter(queryBuffer);
    }
    if (invalidChar) {
        queryPattern[0] = terminalChar;
        nucleoId = 0;
    }
    (*readLength) = nucleoId;
    return queryChar;
}

inline char _SRAFetchReadSequenceUncertaintyNumber(UTBFRBuffer * queryBuffer,
                                unsigned char * charMap,
                                unsigned char * queryPattern,
                                unsigned int maxReadLength,
                                char terminalChar,
                                int * readLength,
                                int * uncertaintyNumber) {

    char queryChar = UTBFRGetCharacter(queryBuffer);
    int nucleoId=0;
    int invalidChar=0;
    int numberOfN=0;
    while (!UTBFREOF(queryBuffer) && queryChar != '>' && queryChar != '@' && queryChar != '+') {
        if (queryChar != '\n') {
            if (queryChar >= 'a' && queryChar <= 'z') {
                queryChar += (unsigned char)('A' - 'a');
            }
            if (queryChar == 'N') numberOfN++;
            if (nucleoId<maxReadLength) {
                if (charMap[queryChar] >= ALPHABET_SIZE) invalidChar=1;
                queryPattern[nucleoId++]=charMap[queryChar];
            }
        }
        queryChar = UTBFRGetCharacter(queryBuffer);
    }
    if (invalidChar) {
        queryPattern[0] = terminalChar;
        nucleoId = 0;
    }
    (*readLength) = nucleoId;
    *uncertaintyNumber = numberOfN;
    return queryChar;
}

char _SRAFetchReadSequenceNoMap(UTBFRBuffer * queryBuffer,
                                char * queryPattern, 
                                unsigned int maxReadLength,
                                int * readLength) {

    char queryChar = UTBFRGetCharacter(queryBuffer);
    int nucleoId=0;
    int invalidChar=0;
    while (!UTBFREOF(queryBuffer) && queryChar != '>' && queryChar != '@' && queryChar != '+') {
        if (queryChar != '\n') {
            if (queryChar >= 'a' && queryChar <= 'z') {
                queryChar += (unsigned char)('A' - 'a');
            }
            if (nucleoId<maxReadLength) {
                queryPattern[nucleoId++]=queryChar;
            }
        }
        queryChar = UTBFRGetCharacter(queryBuffer);
    }
    if (invalidChar) {
        queryPattern[0] = '\0';
    }
    (*readLength) = nucleoId;
    return queryChar;
}

unsigned int SRAQueryGetBatchFromFASTA(UTBFRBuffer * queryBuffer,
                                unsigned char * charMap,
                                unsigned int numOfQuery,
                                unsigned int maxReadLength,
                                unsigned char * queryPattern) {
                                
    memset(queryPattern,ALPHABET_SIZE+1,maxReadLength*numOfQuery);

    int queryLoaded = 0;
    char queryChar = 0;
    int readLength;
    
    while (queryLoaded<numOfQuery && !UTBFREOF(queryBuffer)) {
    
        //Read everything before the entry point of a read, character ">"
        while (!UTBFREOF(queryBuffer) && queryChar != '>') {
            queryChar = UTBFRGetCharacter(queryBuffer);
        }

        //Read the header of a read
        if (UTBFREOF(queryBuffer)) break;
        _SRAFetchReadName(queryBuffer,NULL,-1);

        //Read the pattern body of a read
        if (UTBFREOF(queryBuffer)) break;
        queryChar = _SRAFetchReadSequence(queryBuffer,charMap, queryPattern+(queryLoaded*maxReadLength), maxReadLength, ALPHABET_SIZE+1, &readLength);
        
        queryLoaded++;
        
    }
    queryBuffer->bufferPtr--;
    
    return queryLoaded;
}


unsigned int SRAQueryAndNameGetBatchFromFASTA(UTBFRBuffer * queryBuffer,
                                            unsigned char * charMap,
                                            unsigned int numOfQuery,
                                            unsigned int maxReadNameLength,
                                            char * queryName,
                                            unsigned int maxReadLength,
                                            unsigned char * queryPattern,
                                            char terminalChar) {
                                
    memset(queryName,'\0',maxReadNameLength*numOfQuery);
    memset(queryPattern,terminalChar,maxReadLength*numOfQuery);

    int queryLoaded = 0;
    char queryChar = 0;
    int readLength;
    
    while (queryLoaded<numOfQuery && !UTBFREOF(queryBuffer)) {
    
        //Read everything before the entry point of a read, character ">"
        while (!UTBFREOF(queryBuffer) && queryChar != '>') {
            queryChar = UTBFRGetCharacter(queryBuffer);
        }

        //Read the header of a read
        if (UTBFREOF(queryBuffer)) break;
        _SRAFetchReadName(queryBuffer,queryName+(queryLoaded*maxReadNameLength),maxReadNameLength);

        //Read the pattern body of a read
        if (UTBFREOF(queryBuffer)) break;
        queryChar = _SRAFetchReadSequence(queryBuffer,charMap, queryPattern+(queryLoaded*maxReadLength), maxReadLength, terminalChar, &readLength);
        
        queryLoaded++;
    }
    queryBuffer->bufferPtr--;
    
    return queryLoaded;
}


unsigned int SRAQueryAndNameGetBatchFromFASTAQ(UTBFRBuffer * queryBuffer,
                                            unsigned char * charMap,
                                            unsigned int numOfQuery,
                                            unsigned int maxReadNameLength,
                                            char * queryName,
                                            unsigned int maxReadLength,
                                            unsigned char * queryPattern,
                                            char terminalChar) {
                                
    memset(queryName,'\0',maxReadNameLength*numOfQuery);
    memset(queryPattern,terminalChar,maxReadLength*numOfQuery);

    int queryLoaded = 0;
    char queryChar = 0;
    int isFastq = 0;
    int i;
    int readLength;
    
    while (queryLoaded<numOfQuery && !UTBFREOF(queryBuffer)) {
    
        //Read everything before the entry point of a read, character ">"
        while (!UTBFREOF(queryBuffer) && queryChar != '>' && queryChar != '@') {
            queryChar = UTBFRGetCharacter(queryBuffer);
        }

        isFastq=0;
        if (queryChar=='@') {isFastq=1;}

        //Read the header of a read
        if (UTBFREOF(queryBuffer)) break;
        _SRAFetchReadName(queryBuffer,queryName+(queryLoaded*maxReadNameLength),maxReadNameLength);

        //Read the pattern body of a read
        if (UTBFREOF(queryBuffer)) break;
        queryChar = _SRAFetchReadSequence(queryBuffer, charMap, queryPattern+(queryLoaded*maxReadLength), maxReadLength, terminalChar, &readLength);

        if (isFastq) {
            //Read everything before the entry point of a read, character ">"
            while (!UTBFREOF(queryBuffer) && queryChar != '+') {
                queryChar = UTBFRGetCharacter(queryBuffer);
            }
            //Read the header of a read
            if (UTBFREOF(queryBuffer)) break;
            queryChar = UTBFRGetCharacter(queryBuffer);
            while (!UTBFREOF(queryBuffer) && queryChar != '\n') {
                queryChar = UTBFRGetCharacter(queryBuffer);
            }
            //Read the quality of a read
            if (UTBFREOF(queryBuffer)) break;
            queryChar = UTBFRGetCharacter(queryBuffer);
            for (i=0;i<readLength;i++) {
                queryChar = UTBFRGetCharacter(queryBuffer);
                //In SOAP2 we do not have handling for read quality.
                //therefore the read quality is CURRENTLY read from the input file
                //and discarded immediately.
            }
        }

        queryLoaded++;
    }
    queryBuffer->bufferPtr--;
    
    return queryLoaded;
}


unsigned int SRAQueryAndNameGetBatchFromFASTAQLength(UTBFRBuffer * queryBuffer,
                                            unsigned char * charMap,
                                            unsigned int numOfQuery,
                                            unsigned int maxReadNameLength,
                                            char * queryName,
                                            unsigned int maxReadLength,
                                            unsigned char * queryPattern,
                                            char * queryQuality,
                                            uint16_t * queryLength,
                                            int * uncertaintyNumber) {
                                
    memset(queryName,'\0',maxReadNameLength*numOfQuery);
    memset(queryLength,0,sizeof(uint16_t)*numOfQuery);

    int queryLoaded = 0;
    char queryChar = 0;
    int isFastq = 0;
    int i;
    int readLength;
    
    while (queryLoaded<numOfQuery && !UTBFREOF(queryBuffer)) {
    
        //Read everything before the entry point of a read, character ">"
        while (!UTBFREOF(queryBuffer) && queryChar != '>' && queryChar != '@') {
            queryChar = UTBFRGetCharacter(queryBuffer);
        }

        isFastq=0;
        if (queryChar=='@') {isFastq=1;}

        //Read the header of a read
        if (UTBFREOF(queryBuffer)) break;
        _SRAFetchReadName(queryBuffer,queryName+(queryLoaded*maxReadNameLength),maxReadNameLength);

        //Read the pattern body of a read
        if (UTBFREOF(queryBuffer)) break;
        queryChar = _SRAFetchReadSequenceUncertaintyNumber(queryBuffer, charMap, queryPattern+(queryLoaded*maxReadLength), maxReadLength, '\0', &readLength, uncertaintyNumber+queryLoaded);
        
        queryLength[queryLoaded] = readLength;

        if (isFastq) {
            //Read everything before the entry point of a read, character ">"
            while (!UTBFREOF(queryBuffer) && queryChar != '+') {
                queryChar = UTBFRGetCharacter(queryBuffer);
            }
            //Read the header of a read
            if (UTBFREOF(queryBuffer)) break;
            queryChar = UTBFRGetCharacter(queryBuffer);
            while (!UTBFREOF(queryBuffer) && queryChar != '\n') {
                queryChar = UTBFRGetCharacter(queryBuffer);
            }
            //Read the quality of a read
            if (UTBFREOF(queryBuffer)) break;
            _SRAFetchReadName(queryBuffer,queryQuality+(queryLoaded*maxReadLength),maxReadLength);

            /*queryChar = UTBFRGetCharacter(queryBuffer);
            for (i=0;i<readLength;i++) {
                queryChar = UTBFRGetCharacter(queryBuffer);
                In SOAP2 we do not have handling for read quality.
                therefore the read quality is CURRENTLY read from the input file
                and discarded immediately.
            }*/
        }else{
            queryQuality[queryLoaded*maxReadLength] = 0;
        }

        queryLoaded++;
    }
    queryBuffer->bufferPtr--;
    
    return queryLoaded;
}


unsigned int SRAQueryAndNameGetBatchFromFASTAQNoMap(UTBFRBuffer * queryBuffer,
                                            unsigned int numOfQuery,
                                            unsigned int maxReadNameLength,
                                            char * queryName,
                                            unsigned int maxReadLength,
                                            char * queryPattern) {
                                            
    memset(queryName,'\0',maxReadNameLength*numOfQuery);
    memset(queryPattern,'\0',maxReadLength*numOfQuery);

    int queryLoaded = 0;
    char queryChar = 0;
    int isFastq = 0;
    int i;
    int readLength;
    
    while (queryLoaded<numOfQuery && !UTBFREOF(queryBuffer)) {
    
        //Read everything before the entry point of a read, character ">"
        while (!UTBFREOF(queryBuffer) && queryChar != '>' && queryChar != '@') {
            queryChar = UTBFRGetCharacter(queryBuffer);
        }

        isFastq=0;
        if (queryChar=='@') {isFastq=1;}

        //Read the header of a read
        if (UTBFREOF(queryBuffer)) break;
        _SRAFetchReadName(queryBuffer,queryName+(queryLoaded*maxReadNameLength),maxReadNameLength);

        //Read the pattern body of a read
        if (UTBFREOF(queryBuffer)) break;
        queryChar = _SRAFetchReadSequenceNoMap(queryBuffer, queryPattern+(queryLoaded*maxReadLength), maxReadLength, &readLength);

        if (isFastq) {
            //Read everything before the entry point of a read, character ">"
            while (!UTBFREOF(queryBuffer) && queryChar != '+') {
                queryChar = UTBFRGetCharacter(queryBuffer);
            }
            //Read the header of a read
            if (UTBFREOF(queryBuffer)) break;
            queryChar = UTBFRGetCharacter(queryBuffer);
            while (!UTBFREOF(queryBuffer) && queryChar != '\n') {
                queryChar = UTBFRGetCharacter(queryBuffer);
            }
            //Read the quality of a read
            if (UTBFREOF(queryBuffer)) break;
            queryChar = UTBFRGetCharacter(queryBuffer);
            for (i=0;i<readLength;i++) {
                queryChar = UTBFRGetCharacter(queryBuffer);
                //In SOAP2 we do not have handling for read quality.
                //therefore the read quality is CURRENTLY read from the input file
                //and discarded immediately.
            }
        }
        
        queryLoaded++;

    }
    queryBuffer->bufferPtr--;
    
    return queryLoaded;
}

/*
void SRAFetchReadFromPacked(unsigned int * packedQueries,
                            unsigned int readId, unsigned int readLength,
                            char * output) {

    unsigned int idx = readId % 32;
    unsigned int offset = readId / 32 * 32 * GPU_WORD_PER_QUERY + idx;
    unsigned int *addr = packedQueries + offset;

       unsigned int chars16;

       unsigned int i = 0;
       for (i = 0; i < readLength; ++i) { // backward search
               unsigned int m = i & 15;
               if (m == 0) {
                 chars16 = *(addr + (i << 1));
               }
               unsigned int c = ((chars16 >> (m*2)) & 3);
               output[i] = dnaChar[c];
       }
       output[i]='\0';
       return;
}
*/

