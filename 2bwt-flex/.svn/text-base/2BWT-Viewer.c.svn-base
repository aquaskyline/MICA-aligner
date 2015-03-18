//
//    2BWT-Viewer.c
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
/* 

    SOAP2 Viewer for SOAP2 Aligner
    Date:    2010-11-13
    Author:  Edward Wu
    Website: http://www.cs.hku.hk/~mkewu/2bwt

    Designed for output format 20110820

    Modification History 

    Date   : 8th May 2011
    Author : Edward MK Wu
    Change : New file.

    Date   : 16th May 2011
    Author : Edward MK Wu
    Change : Re-design of 2BWT-Viewer. No output format change.
             The read file is now entirely loaded in the memory.
             With the assumption that the number of reads does not
             exceed SRA_QUERY_BATCH_SIZE * SRA_MAX_NUM_BATCH.
    
    Date   : 19th June 2011
    Author : Edward MK Wu
    Change : Packaging 2BWT library as a separate product.
             Thus, changing all references to 2bwt lib to subdirectory.

*/

#include "2BWT-Viewer.h"

//Define the debug flag to turn on extra output.
//#define DEBUG_2BWT_PRINT_MESSAGE
//#define DEBUG_2BWT_PRINT_CELL

unsigned char charMap[256];
char AlignmentTypeString[2][5] = {"MISM","EDIT"};

void SRAQueryBatchInitialise(SRAQueryBatch * sraQueryBatch) {
    sraQueryBatch->startIdx = 0;
    sraQueryBatch->endIdx = 0;
}

SRAQueryBatch * SRAQueryBatchConstruct(unsigned int batchSize) {
    SRAQueryBatch * sraQueryBatch = (SRAQueryBatch*) malloc(sizeof(SRAQueryBatch));
    sraQueryBatch->names = (char*) malloc(sizeof(char) * batchSize * (MAX_SEQ_NAME_LENGTH));
    sraQueryBatch->patterns = (unsigned char*) malloc(sizeof(unsigned char) * batchSize * (SRA_MAX_READ_LENGTH));
    sraQueryBatch->batchSize = batchSize;
    SRAQueryBatchInitialise(sraQueryBatch);
    return sraQueryBatch;
}

void SRAQueryBatchFree(SRAQueryBatch * sraQueryBatch) {
    free(sraQueryBatch->names);
    free(sraQueryBatch->patterns);
    free(sraQueryBatch);
}

void SRAQueryBatchFill(SRAQueryBatch * sraQueryBatch, UTBFRBuffer * queryBuffer) {

    unsigned int queryInBatch = SRAQueryAndNameGetBatchFromFASTAQ(queryBuffer,
                                charMap,
                                sraQueryBatch->batchSize,
                                MAX_SEQ_NAME_LENGTH,
                                sraQueryBatch->names,
                                SRA_MAX_READ_LENGTH,
                                sraQueryBatch->patterns,
                                ALPHABET_SIZE+1);
                                
    sraQueryBatch->startIdx = sraQueryBatch->endIdx + 1;
    sraQueryBatch->endIdx += queryInBatch;
    #ifdef DEBUG_2BWT_PRINT_MESSAGE
        fprintf(stdout, "[SRAQueryBatchFill] invoked serving %llu to %llu\n",sraQueryBatch->startIdx,sraQueryBatch->endIdx);
    #endif
}

inline int GetReadPatternLength(unsigned char * queryPattern) {
    int i;
    int length = 0;
    for (i=0;i<SRA_MAX_READ_LENGTH;i++) {
        if (queryPattern[i]==ALPHABET_SIZE+1) {
            break;
        }
        length++;
    }
    return length;
}

void SRAFillCharMap(unsigned char * charMap);

int main(int argc, char ** argv) {
    #ifdef DEBUG_2BWT_PRINT_MESSAGE
        fprintf(stdout, "SOAP2 Viewer initiated.\n");
    #endif

    char * readFileName;
    char * readFileName_mate;
    char * outputFileName;
    OCCPositionCacheToDisk OUTResultBuffer;

    unsigned long long OUTNumOfOcc=0;
    
    char *          OUTCurrentQueryName[2];
    unsigned char * OUTCurrentQuerySequence[2];
    int             OUTCurrentQueryLength[2];
        
    unsigned char       cacheFlag1,cacheFlag2;
    unsigned long long  cacheFlag3;

    unsigned char sraReadStrand;
    unsigned long long sraChromId, sraReadId, sraOffset;
        
    unsigned int loadedBatchCount = 0;
    unsigned int loadedBatchCount_mate = 0;
    UTBFRBuffer * queryBuffer, * queryBuffer_mate;
    
    unsigned long long i,j,k;
    uint8_t argImpliedQueryMode = 0;
    
    if (argc==3) {
        fprintf(stderr, "Output Result Type = Single-end Alignment\n");
        argImpliedQueryMode = SHORT_READ_ALIGNMENT;
        readFileName = (char*) argv[1];
        outputFileName = (char*) argv[2];
    } else if (argc==4) {
        fprintf(stderr, "Output Result Type = Pair-end Alignment\n");
        argImpliedQueryMode = PAIR_END_ALIGNMENT;
        readFileName = (char*) argv[1];
        readFileName_mate = (char*) argv[2];
        outputFileName = (char*) argv[3];
    } else {
        printf("\n[Main] %s Viewer v%d.%d.%d (%s) - Usage guide:\n",PROJECT_NAME,PROJECT_MAJOR,PROJECT_MINOR,PROJECT_REV,PROJECT_SPECIAL);
        printf("\n");
        printf("For single-end alignment result:\n");
        printf("  Windows Platform:\n\t %s <query_file> <output_file>\n",PROJECT_VIEWER_BINARY);
        printf("  Linux Platform:\n\t ./%s <query_file> <output_file>\n",PROJECT_VIEWER_BINARY);
        printf("For pair-end alignment result:\n");
        printf("  Windows Platform:\n\t %s <query_file_1> <query_file_2> <output_file>\n",PROJECT_VIEWER_BINARY);
        printf("  Linux Platform:\n\t ./%s <query_file_1> <query_file_2> <output_file>\n",PROJECT_VIEWER_BINARY);
        printf("Note: <output_file> must be the corresponding output file of running %s with <query_file>\n\n",PROJECT_ALIGNER_BINARY);
        return 1;
    }
    
    SRAFillCharMap(charMap);
    
    
    
    
    // =======================================================
    // | Open various input and output files
    // =======================================================
    //Open read file
    fprintf(stderr, "Opening %s .. \n",readFileName);
    queryBuffer = UTBFRLoad(readFileName);
    
    //Open read file
    if (argImpliedQueryMode == PAIR_END_ALIGNMENT) {
        fprintf(stderr, "Opening %s .. \n",readFileName_mate);
        queryBuffer_mate = UTBFRLoad(readFileName_mate);
    }
    
    //Open output file
    SRAOutputFileInfo * sraOutputFileInfo = SRAOutputOpen(outputFileName);
    
    
    
    // =======================================================
    // | Initialise the SRA query batch
    // =======================================================
    unsigned long long sraQueryBatchStartIdx = 0;
    unsigned long long sraQueryBatchEndIdx = 0;
    SRAQueryBatch * sraQueryBatch      = SRAQueryBatchConstruct(sraOutputFileInfo->queryBatch);
    SRAQueryBatch * sraQueryBatch_mate = SRAQueryBatchConstruct(sraOutputFileInfo->queryBatch);
     
    
    
    // =======================================================
    // | Checking the output file
    // =======================================================
    if (argImpliedQueryMode != sraOutputFileInfo->queryMode) {
        fprintf(stderr, "Viewer call style does not match output file.\n");
        fprintf(stderr, "Output file contains result of ");
        if (sraOutputFileInfo->queryMode == SHORT_READ_ALIGNMENT) {
            fprintf(stderr, "single-end alignment");
        } else if (sraOutputFileInfo->queryMode == PAIR_END_ALIGNMENT) {
            fprintf(stderr, "pair-end alignment");
        } else {
            fprintf(stderr, "unknown alignment");
        }
        fprintf(stderr, "\n");
        return 1;
    }
    
    if (sraOutputFileInfo->outputFormat == SRA_OUTPUT_FORMAT_PLAIN ||
                sraOutputFileInfo->outputFormat == SRA_OUTPUT_FORMAT_SAM) {
        fprintf(stderr, "Output Type = Plain / SAM\n");
        char * dumpFileNamePtr = sraOutputFileInfo->dumpFileNames;
        unsigned int fileId = 0;
        for (fileId=0;fileId<sraOutputFileInfo->dumpNumOfFiles;fileId++) {
            fprintf(stderr, "[SRAOutputOpen] Opening %s ..",dumpFileNamePtr);
            FILE * dumpFile = (FILE*)fopen64(dumpFileNamePtr,"r");
            if (dumpFile==NULL || dumpFile ==0) {
                fprintf(stderr, "Failed!\n");
                exit(1);
            } else {
                fprintf(stderr, "Done!\n");
            }
            if (sraOutputFileInfo->outputFormat == SRA_OUTPUT_FORMAT_SAM && fileId>0) {
                size_t lineBufferSize = OCC_OUTPUT_MAX_LINE_LENGTH;
                char * lineBuffer = (char*) malloc(sizeof(char) * lineBufferSize);
                unsigned int charInLine;
                //Ignore the header lines
                for (i=0;i<1+sraOutputFileInfo->dbNumOfSeqs;i++) {
                    charInLine = getline(&lineBuffer, &lineBufferSize, dumpFile);
                }
                free(lineBuffer);
            }
            char c = fgetc(dumpFile);
            while (!feof(dumpFile)) {
                fputc(c,stdout);
                c = fgetc(dumpFile);
            }
            fclose(dumpFile);
            dumpFileNamePtr += (OCC_OUTPUT_MAX_LINE_LENGTH+1);
        }
        return 0;
    } else {
        fprintf(stderr, "Output Type = Unknown\n");
        return 1;
    }
    
    // =======================================================
    // | Parsing the payload
    // =======================================================
    
    char * dumpFileNamePtr = sraOutputFileInfo->dumpFileNames;
    
    printf("#Read_Pos\tRead_Name\tRead\tChromosome\tOffset\tRead_Length\tStrand\tType");
        
    unsigned int fileId = 0;
    for (fileId=0;fileId<sraOutputFileInfo->dumpNumOfFiles;fileId++) {
    
        SRAQueryBatchInitialise(sraQueryBatch);
        SRAQueryBatchFill(sraQueryBatch,queryBuffer);
        if (argImpliedQueryMode == PAIR_END_ALIGNMENT) SRAQueryBatchFill(sraQueryBatch_mate,queryBuffer_mate);
    
        fprintf(stderr, "[SRAOutputOpen] Opening %s ..",dumpFileNamePtr);
        FILE * dumpFile = (FILE*)fopen64(dumpFileNamePtr,"rb");
        if (dumpFile==NULL || dumpFile ==0) {
            fprintf(stderr, "Failed!\n");
            exit(1);
        } else {
            fprintf(stderr, "Done!\n");
        }
        
        ////////////////////////////////////////////////////
        //Check file body size
        ////////////////////////////////////////////////////
        fseek(dumpFile, 0, SEEK_END);
        unsigned long long dumpFileSize = ftell(dumpFile);
        if (dumpFileSize % sizeof(OCCPositionCacheToDisk)>0) {
            fprintf(stderr,"Bad file format in %s.\n[Content:%llu LeftOver:%llu]!\n",
                        dumpFileNamePtr,
                        dumpFileSize,
                        dumpFileSize % sizeof(OCCPositionCacheToDisk));
            
            return 0;
        }
        fseek(dumpFile, 0, SEEK_SET);
        unsigned long long recordInDumpFile = dumpFileSize / sizeof(OCCPositionCacheToDisk);
        
        ////////////////////////////////////////////////////
        //Parse file body
        ////////////////////////////////////////////////////
        int occEven = 0;
        for (i=0;i<recordInDumpFile;i++) {
            if (SRAOutputGetNextRow(dumpFile,&OUTResultBuffer)==0) break;

            cacheFlag1 = OUTResultBuffer.cell[0];
            cacheFlag2 = OUTResultBuffer.cell[1];
            cacheFlag3=(*(unsigned long long*)&OUTResultBuffer.cell[2]);
            #ifdef DEBUG_2BWT_PRINT_CELL
                fprintf(stdout,"\nRecord is read : %u %u %llu\n",cacheFlag1, cacheFlag2, cacheFlag3);
            #endif

            if (cacheFlag2==0) {
                //Delimitor
                sraReadId=cacheFlag3;
                
                while ( sraReadId>sraQueryBatch->endIdx ) {
                        
                    SRAQueryBatchFill(sraQueryBatch,queryBuffer);
                    if (argImpliedQueryMode == PAIR_END_ALIGNMENT) SRAQueryBatchFill(sraQueryBatch_mate,queryBuffer_mate);
                    //if (sraQueryBatch->startIdx<=sraQueryBatch->endIdx) {
                    //    fprintf(stdout,"\nIncorrect file body. Unexpected read Id %llu.\n",sraReadId);
                    //}
                }
                
                OUTCurrentQueryName[0] = sraQueryBatch->names + (sraReadId - sraQueryBatch->startIdx) * (MAX_SEQ_NAME_LENGTH);
                OUTCurrentQuerySequence[0] = sraQueryBatch->patterns + (sraReadId-sraQueryBatch->startIdx) * (SRA_MAX_READ_LENGTH);
                OUTCurrentQueryLength[0] = GetReadPatternLength(OUTCurrentQuerySequence[0]);
                
                if (argImpliedQueryMode == PAIR_END_ALIGNMENT) {
                    OUTCurrentQueryName[1] = sraQueryBatch_mate->names + (sraReadId - sraQueryBatch->startIdx) * (MAX_SEQ_NAME_LENGTH);
                    OUTCurrentQuerySequence[1] = sraQueryBatch_mate->patterns + (sraReadId-sraQueryBatch->startIdx) * (SRA_MAX_READ_LENGTH);
                    OUTCurrentQueryLength[1] = GetReadPatternLength(OUTCurrentQuerySequence[1]);
                }

                #ifdef DEBUG_2BWT_PRINT_CELL
                    fprintf(stdout,"Query Name = %s\n",OUTCurrentQueryName[0]);
                    fprintf(stdout,"Query Pattern = ");
                    for (j=0;j<OUTCurrentQueryLength[0];j++) {
                        printf("%c",dnaChar[OUTCurrentQuerySequence[0][j]]);
                    }
                    fprintf(stdout,"\nQuery Length = %d\n",OUTCurrentQueryLength[0]);
                    
                    if (argImpliedQueryMode == PAIR_END_ALIGNMENT) {
                        fprintf(stdout,"Query Name = %s\n",OUTCurrentQueryName[1]);
                        fprintf(stdout,"Query Pattern = ");
                        for (j=0;j<OUTCurrentQueryLength[1];j++) {
                            printf("%c",dnaChar[OUTCurrentQuerySequence[1][j]]);
                        }
                        fprintf(stdout,"\nQuery Length = %d\n",OUTCurrentQueryLength[1]);
                    }
                #endif
            } else {
                //Occurrence position
                sraReadStrand = cacheFlag1;
                sraChromId = cacheFlag2-1;
                sraOffset  = cacheFlag3;
                
                printf("\n%llu\t",sraReadId);
                for (j=0;j<MAX_SEQ_NAME_LENGTH;j++) {
                    if (OUTCurrentQueryName[occEven][j]=='\t' || 
                        OUTCurrentQueryName[occEven][j]==' ' || 
                        OUTCurrentQueryName[occEven][j]=='\n' || 
                        OUTCurrentQueryName[occEven][j]=='\0' || 
                        OUTCurrentQueryName[occEven][j]=='\r') {break;}
                    printf("%c",OUTCurrentQueryName[occEven][j]);
                }
                printf("\t");
                for (j=0;j<OUTCurrentQueryLength[occEven];j++) {
                    printf("%c",dnaChar[OUTCurrentQuerySequence[occEven][j]]);
                }
                printf("\t");
                if (sraOutputFileInfo->dbAnnotation[sraChromId].gi>0) printf("%u|",sraOutputFileInfo->dbAnnotation[sraChromId].gi);
                printf("%s\t%llu\t%u\t",sraOutputFileInfo->dbAnnotation[sraChromId].text,sraOffset,OUTCurrentQueryLength[occEven]);
                if (sraReadStrand==QUERY_POS_STRAND) {printf("+\t");} 
                else if (sraReadStrand==QUERY_NEG_STRAND) {printf("-\t");}
                else {printf("?\t");}
                OUTNumOfOcc++;
                occEven = (++occEven) % 2;
                occEven &= (argImpliedQueryMode == PAIR_END_ALIGNMENT);
            }
        }
        fprintf(stderr,"\nThere were %llu alignments.\n",OUTNumOfOcc);
        fclose(dumpFile);
        
        dumpFileNamePtr += (OCC_OUTPUT_MAX_LINE_LENGTH+1);
        
        UTBFRResetBuffer(queryBuffer);
        if (argImpliedQueryMode == PAIR_END_ALIGNMENT) UTBFRResetBuffer(queryBuffer_mate);
    }
    printf("\n");
    
    // =======================================================
    // | Clean up
    // =======================================================
    SRAOutputFree(sraOutputFileInfo);
    UTBFRFree(queryBuffer);
    if (argImpliedQueryMode == PAIR_END_ALIGNMENT)  UTBFRFree(queryBuffer_mate);
    free(sraQueryBatch);
    free(sraQueryBatch_mate);
    return 0;
}

void SRAFillCharMap(unsigned char * charMap) {
    int i;
    for (i=0;i<256;i++) { charMap[i] = 2; };
    charMap['A'] = 0;
    charMap['C'] = 1;
    charMap['G'] = 2;
    charMap['T'] = 3;
}

