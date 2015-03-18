/*

   Date   : 5th September 2013
   Author : Edward MK Wu
   Change : 2BWT Index Utilities

*/

///////////////////////////////////////
// Including Standard Libraries
///////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>

///////////////////////////////////////
// Including 2BWT Libraries
///////////////////////////////////////
#include "../2bwt-lib/iniparser.h"
#include "../2bwt-lib/2BWT-Interface.h"
#include "../2bwt-lib/Timing.h"

///////////////////////////////////////
// Including SOAP2-DP Libraries
///////////////////////////////////////
#include "../2BWT-PEAlgnmt.h"
#include "../2BWT-SRAAlgnmt.h"
#include "../SRAQueryParser.h"

#define BWT_OPS_CONVERT_PATT_TO_SARANGE   1
#define SA_OPS_CONVERT_SAI_TO_SAV         10
#define HSP_OPS_CONVERT_AMB_TO_CHROFFSET  20
#define FASTA_OPS_CONVERT_COORD_TO_OUTPUT 30
#define SEQ_OPS_BRUTE_FORCE_COMPARE       40

unsigned char charMap[256];
unsigned char complementMap[256];
char DatabaseName[MAX_FILENAME_LEN+1] = "";
char FASTAName[MAX_FILENAME_LEN+1] = "";
char InputSaValueFileName[MAX_FILENAME_LEN+1] = ".sa8";

void ParseIniFile(char *iniFileName);
dictionary *ParseInput(int argc, char** argv);
void PrintHelp();
void SRAFillCharMap(unsigned char * charMap);

typedef struct FastaSequence {
    int length;
    char * body;
} FastaSequence;

typedef struct FastaFile {
    FastaSequence * seq;
} FastaFile;

void BwtConvertPatternToSaRange(Idx2BWT * idx2BWT) {
    unsigned int ui;
    char c = '\0';
    int readLength = 0;
    
    unsigned char pattern[SRA_MAX_READ_LENGTH];
    
    printf("Pattern (End at EOL) = "); fflush(stdout);
    c = getchar();
    while ( c != '\n' ) {
        pattern[readLength++] = c;
        c = getchar();
    }
    pattern[readLength] = '\0';
    
    BWT * bwt = idx2BWT->bwt;
    
    
    printf("Read Length = %d\n",readLength);
}

void SaConvertSaiToSaValue(Idx2BWT * idx2BWT) {
    unsigned int ui;
    
    printf("SA Index = "); fflush(stdout);
    scanf("%u",&ui);
    
    BWT * bwt = idx2BWT->bwt;
    
    printf("SA Value = %u\n",BWTSaValue(bwt,ui));
}

void HspConvertAmbToChrOffset(Idx2BWT * idx2BWT) {
    unsigned int ui;
    unsigned char chrId;
    unsigned long long unambPos;
    
    printf("Ambitious Location = "); fflush(stdout);
    scanf("%u",&ui);
    
    Translate * translate = idx2BWT->hsp->translate;
    unsigned short * ambiguityMap = idx2BWT->hsp->ambiguityMap;
    
    OCCTranslateOccurrence(translate,ambiguityMap,ui,&chrId,&unambPos);
    printf("Chromosome %u Offset %llu\n",chrId,unambPos);
}

void FastaExtractSeqBody(Idx2BWT * idx2BWT, FastaFile * fastaFile) {
    unsigned char chrId;
    unsigned long long unambPos;
    unsigned long long length;
    printf("Chromosome ID (Based 1) = "); fflush(stdout);
    scanf("%u",&chrId);chrId--;
    printf("Offset (Based 1) = "); fflush(stdout);
    scanf("%llu",&unambPos);unambPos--;
    printf("Length = "); fflush(stdout);
    scanf("%llu",&length);
    
    unsigned long long i;
    printf("> Extracted %s (+ %llu-%llu)\n",idx2BWT->hsp->annotation[chrId].decoratedText,unambPos,unambPos+length-1);
    for (i=0;i<length;i++) {
        printf("%c",fastaFile->seq[chrId].body[unambPos+i]);
    }
    printf("\n");
    printf("> Extracted %s (- %llu-%llu)\n",idx2BWT->hsp->annotation[chrId].decoratedText,unambPos,unambPos+length-1);
    for (i=0;i<length;i++) {
        printf("%c",complementMap[fastaFile->seq[chrId].body[unambPos+length-i-1]]);
    }
    printf("\n");
}

void SeqCompareStrings() {

    unsigned int maxLength = 1024;
    
    int i, lenA, lenB;
    char * stringA = (char*) malloc (sizeof(char)*1024);
    char * stringB = (char*) malloc (sizeof(char)*1024);
    
    printf("String A = "); fflush(stdout);
    scanf("%s",stringA);
    lenA = strlen(stringA);
    
    printf("String B = "); fflush(stdout);
    scanf("%s",stringB);
    lenB = strlen(stringB);
    
    printf("\n%s\n",stringA);
    for (i=0;i<lenA;i++) {
        if (i>=lenB) {
            printf("?");
        } else if (stringA[i]==stringB[i]) {
            printf(" ");
        } else {
            printf("X");
        }
    }printf("\n");
    printf("%s\n",stringB);
    
    free(stringA);
    free(stringB);
}

int main(int argc, char** argv) {

    ////////////////////////////////////////////////////////////
    // Declaration of Variables
    ////////////////////////////////////////////////////////////
    // Program input
    dictionary *programInput;
    unsigned long long i,j,k;
    unsigned int u,v,w;

    HSPFillCharMap(charMap);
    HSPFillComplementMap(complementMap);
    
    ////////////////////////////////////////////////////////////
    // Ini Configuration
    ////////////////////////////////////////////////////////////
    char iniFilename[MAX_FILENAME_LEN];
    sprintf(iniFilename, "%s.ini", argv[0]);
    ParseIniFile(iniFilename);

    // Command Argument - Override
    programInput = ParseInput(argc, argv);
    ////////////////////////////////////////////////////////////
    // Index Handling
    ////////////////////////////////////////////////////////////
    printf("Loading index %s ... ", DatabaseName); fflush(stdout);
    Idx2BWT * idx2BWT = BWTLoad2BWT(DatabaseName,InputSaValueFileName);
    printf("DONE\n\n");
    
    UTBFRBuffer * fastaFilePtr = UTBFRLoad(FASTAName);
    FastaFile fastaFile;
    fastaFile.seq = (FastaSequence*) malloc(sizeof(FastaSequence) * idx2BWT->hsp->numOfSeq);
    char RefSeqName[1024];
    for (i=0;i<idx2BWT->hsp->numOfSeq;i++) {
        fastaFile.seq[i].length = idx2BWT->hsp->seqActualOffset[i].endPos - idx2BWT->hsp->seqActualOffset[i].startPos + 1;
        fastaFile.seq[i].body = (char*)malloc(sizeof(char) * (fastaFile.seq[i].length+1));
        
        _SRAFetchReadName(fastaFilePtr,RefSeqName,1024);
        printf("Reading %s expecting length %u .. ",RefSeqName,fastaFile.seq[i].length); fflush(stdout);
        
        _SRAFetchReadSequenceNoMap(fastaFilePtr,fastaFile.seq[i].body,fastaFile.seq[i].length,&(fastaFile.seq[i].length));
        printf("DONE(%u)\n",fastaFile.seq[i].length); fflush(stdout);
    }
    UTBFRFree(fastaFilePtr);
    
    unsigned int ui;
    while (1) {
    
        printf("\n=================================\n");
        printf("  DEV INDEX UTILITIES\n");
        printf("=================================\n");
        printf("Please select the operations:\n");
        printf("0. Exit\n");
        //printf("- BWT Operations\n");
        //printf("%d. Convert Pattern to SA Range (Exact-Matching)\n",BWT_OPS_CONVERT_PATT_TO_SARANGE);
        printf("- SA Retrieve\n");
        printf("%d. Convert SA Index to Value.\n",SA_OPS_CONVERT_SAI_TO_SAV);
        printf("- HSP Data Mine\n");
        printf("%d. Convert ambPosition to ChrID-Offset.\n",HSP_OPS_CONVERT_AMB_TO_CHROFFSET);
        printf("- REF-SEQ Data Mine\n");
        printf("%d. Extract ChrID-Offset from sequence Body.\n",FASTA_OPS_CONVERT_COORD_TO_OUTPUT);
        printf("- SEQ Manipulation\n");
        printf("%d. Brute-force comparison of two sequence bodies.\n",SEQ_OPS_BRUTE_FORCE_COMPARE);
        
        printf("\nChoice = ");fflush(stdout);
        scanf("%u",&ui);
        printf("\n");
        
        if (ui==0) {
            break;
        } else if (ui==BWT_OPS_CONVERT_PATT_TO_SARANGE) {
            BwtConvertPatternToSaRange(idx2BWT);
        } else if (ui==SA_OPS_CONVERT_SAI_TO_SAV) {
            SaConvertSaiToSaValue(idx2BWT);
        } else if (ui==HSP_OPS_CONVERT_AMB_TO_CHROFFSET) {
            HspConvertAmbToChrOffset(idx2BWT);
        } else if (ui==FASTA_OPS_CONVERT_COORD_TO_OUTPUT) {
            FastaExtractSeqBody(idx2BWT,&fastaFile);
        } else if (ui==SEQ_OPS_BRUTE_FORCE_COMPARE) {
            SeqCompareStrings();
        } else {
            printf("[ERROR] Do not understand input!\n");
        }
    }
    
    ////////////////////////////////////////////////////////////
    // Free Memory Allocated
    ////////////////////////////////////////////////////////////
    BWTFree2BWT(idx2BWT);
    
    return 0;
}


void ParseIniFile(char *iniFileName) {

    dictionary *ini;

    printf("Loading %s ..", iniFileName);
    ini = iniparser_load(iniFileName, FALSE);
    if (ini == NULL) {
        printf("not found.\n");
        return;
    }
    printf("done.\n");

    iniparser_freedict(ini);

}

dictionary *ParseInput(int argc, char** argv) {
    dictionary *programInput;
    char t1[3] = "-c";    // specify that this is a boolean type parameter; no following argument
    char t2[3] = "-U";    // specify that this is a boolean type parameter; no following argument
    char t3[3] = "-A";    // specify that this is a boolean type parameter; no following argument
    char *d[3];

    char *tempString;
    int len;
    int i;

    d[0] = t1;
    d[1] = t2;
    d[2] = t3;

    programInput = paraparser_load(argc, argv, 3, d);    // 4 parameters are boolean type

    if (!iniparser_find_entry(programInput, "argument:1")) {
        PrintHelp();
        exit(1);
    }
    iniparser_copystring(programInput, "argument:1", FASTAName, FASTAName, MAX_FILENAME_LEN);

    
    if (!iniparser_find_entry(programInput, "argument:2")) {
        PrintHelp();
        exit(1);
    }
    iniparser_copystring(programInput, "argument:2", DatabaseName, DatabaseName, MAX_FILENAME_LEN);

    return programInput;

}

void PrintHelp() {
    printf("Usage Syntax:\n");
    printf("    ./devIndexUtils <ref seq> <2bwt index>\n");
    printf("\n");
}

void SRAFillCharMap(unsigned char * charMap) {
    int i;
    for (i=0;i<256;i++) {
        charMap[i] = i;
    }
}
