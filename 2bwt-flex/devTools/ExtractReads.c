//
//    devTools/ExtractReads.c
//
//    soap2-dp
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

#include "ExtractReads.h"

int main ( int argc, char ** argv ) {

    int i;
    if ( argc != 2 ) {
        printf("Usage\n% ./extractReads <FASTA File>\n");
        exit(1);
    }
    printf("Read file %s.\n",argv[1]);
    
    SRAOCCCollector * occCollector = SRAOCCCreate(SRAOCC_FLOOD_TYPE_EXPAND);
    unsigned int readId;
    unsigned long long readCount;
    SRAOccurrence sraOcc;
    memset(&sraOcc,0,sizeof(SRAOccurrence));
    
    printf("Reading Read Index (1-based) (0 to escape)\n");
    while (scanf("%u", &readId)) {
        if (readId==0) break;
        sraOcc.type = SRAOCC_TYPE_AWAIT_TRANSLATE;
        sraOcc.ambPosition = readId;
        SRAOCCAddOccurrenceToBuckets(occCollector,&sraOcc);
    }
    
    readCount = SRAOCCCountOccurrences(occCollector);
    printf("%llu read items to be extracted.\n",readCount);
    SRAOccurrence * allReadIds = (SRAOccurrence*) malloc(sizeof(SRAOccurrence)*SRAOCCCountOccurrences(occCollector));
    SRAOCCPopulateSRAOccList(occCollector,allReadIds);
    //for (i=0;i<readCount;i++) {
    //    printf("%llu\n",allReadIds[i].ambPosition);
    //}
    printf("Sorting read indexes..\n");
    SRAOccurrencesSort(allReadIds,&readCount);

    printf("Opening read file %s.\n",argv[1]);
    UTBFRBuffer * queryBuffer = UTBFRLoad(argv[1]);
    char tmpPatternNames[MAX_SEQ_NAME_LENGTH+1];
    char tmpPatternReads[SRA_MAX_READ_LENGTH+1];
    unsigned int queryInBatch = 1;
    
    unsigned int traversingReadIdx = 0;
    unsigned int printingReadIdx = 0;
    while (queryInBatch>0 && printingReadIdx<readCount) {
        queryInBatch = SRAQueryAndNameGetBatchFromFASTAQNoMap(queryBuffer,
                                                                1,
                                                                MAX_SEQ_NAME_LENGTH,
                                                                tmpPatternNames,
                                                                SRA_MAX_READ_LENGTH,
                                                                tmpPatternReads);
        traversingReadIdx++;
        if (allReadIds[printingReadIdx].ambPosition==traversingReadIdx) {
            printf(">%u/%llu-%s\n%s\n",printingReadIdx,readCount,tmpPatternNames,tmpPatternReads);
            printingReadIdx++;
        } else if (allReadIds[printingReadIdx].ambPosition<traversingReadIdx) {
            printingReadIdx++;
        }
    }
    free(allReadIds);
    SRAOCCFree(occCollector);
    UTBFRFree(queryBuffer);
    return 0;
}
