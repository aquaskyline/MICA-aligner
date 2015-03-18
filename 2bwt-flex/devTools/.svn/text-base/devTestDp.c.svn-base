#include "../2bwt-lib/iniparser.h"
#include "../2bwt-lib/MiscUtilities.h"
#include "../2bwt-lib/MemManager.h"
#include "../2bwt-lib/Timing.h"
#include "../2bwt-lib/TextConverter.h"
#include "../2bwt-lib/BWT.h"
#include "../2bwt-lib/HSP.h"
#include "../2bwt-lib/HSPstatistic.h"
#include "../2bwt-lib/blast_dust.h"
#include "../2bwt-lib/Socket.h"

#include "../DPOperations.h"
#include "../PEArguments.h"
#include "../PE2BWTReport.h"
#include "../DPArguments.h"
#include "../PEDPReport.h"

/*
#define PackedDNAFileName  "db2/ncbi.genome1-5.fa.index.pac"
#define AnnotationFileName "db2/ncbi.genome1-5.fa.index.ann"
#define AmbiguityFileName  "db2/ncbi.genome1-5.fa.index.amb"
#define TranslateFileName  "db2/ncbi.genome1-5.fa.index.tra"
*/

#define PackedDNAFileName  "db6/genome.37.1.fa.index.pac"
#define AnnotationFileName "db6/genome.37.1.fa.index.ann"
#define AmbiguityFileName  "db6/genome.37.1.fa.index.amb"
#define TranslateFileName  "db6/genome.37.1.fa.index.tra"

// Debug Flag 1: Enable the testing of tailVerifyLength
//#define DEV_TEST_DP_TEST_TAIL_VERIFY_LENGTH

// Debug Flag 2: Enable the testing of TailSoftClip
//#define DEV_TEST_DP_TEST_TAIL_SOFT_CLIP

// Debug Flag 2: Enable the testing of HeadSoftClip
//#define DEV_TEST_DP_TEST_HEAD_SOFT_CLIP


void DevTestDpPrintDpOcc(DPOccurrence * dpOcc) {
    printf("DPOccurrence is found at %llu\n",dpOcc->ambPosition);
    int i; for (i=0;i<dpOcc->matchElemsCount;i++) {
        printf("%d-%d\t\t",dpOcc->matchElems[i].length ,dpOcc->matchElems[i].type);
    }
}

int main() {
    printf("testDp -- invoked\n");
    
    MMPool *mmPool;
    int PoolSize = 2097152;                // 2M  - fixed; not configurable through ini
    MMMasterInitialize(1, 0, FALSE, NULL);
    mmPool = MMPoolCreate(PoolSize);
    HSP * hsp = HSPLoad(mmPool, PackedDNAFileName, AnnotationFileName, AmbiguityFileName,TranslateFileName, 1);
    //=======================
    DPOccurrence dpOcc_found;
    dpOcc_found.matchElemsCount=0;
    
    printf("testDp -- begining to initialise the matrix\n");
    
    DPArguments * dpArgument = DPARGConstruct();
    dpArgument->dpScores->dpMatch = 1;
    dpArgument->dpScores->dpMismatch = -2;
    dpArgument->dpScores->dpGapOpen = -3;
    dpArgument->dpScores->dpGapExtend = -1;
    
    DPWork * dpWork = dpArgument->dpWork;
    dpWork->regionStartIdx = 0;
    /*
    dpWork->regionLength = 100;
    char * pattern = "AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTGGCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC";
    //*/
    
    /*
    dpWork->regionLength = 150;
    char * pattern = "ATGGAATGGAATGGAATGGAATTAAACCGAATAGAATGGAATGGAATGGAATGGAACGGAACAGAACGGAACGGAACGGAATGGAGTGGAATGGAATGGA";
    //*/
    
    // Testing DP-Tweak
    #ifdef DEV_TEST_DP_TEST_TAIL_SOFT_CLIP
        dpWork->totalSoftClipLength = 100;
        dpWork->tailSoftClipLength = 2;
    #endif
    
    // Testing DP-Tweak
    #ifdef DEV_TEST_DP_TEST_HEAD_SOFT_CLIP
        dpWork->totalSoftClipLength = 100;
        dpWork->leadSoftClipLength = 2;
    #endif
    
    // Soft-Clip Head
    /*
    dpWork->regionStartIdx = 0;
    dpWork->regionLength = 11;
    char * pattern = "CGGCCCTAAC";
    //*/
    /*
    // Soft-Clip Tail
    dpWork->regionStartIdx = 0;
    dpWork->regionLength = 11;
    char * pattern = "TAACCCTTTT";
    //*/
    /*
    // Soft-Clip Total
    dpWork->regionStartIdx = 0;
    dpWork->regionLength = 11;
    char * pattern = "CGACCCTTTT";
    //*/
    /*
    // Soft-Clip Head - Boundary Case
    dpWork->regionStartIdx = 0;
    dpWork->regionLength = 11;
    char * pattern = "CCCCCCTA";
    // NOTE: Set dpWork->tailSoftClipLength = 0; dpWork->totalSoftClipLength = 2;
    
    dpWork->leadSoftClipLength = 2;
    dpWork->tailSoftClipLength = 0;
    dpWork->totalSoftClipLength = 2;
    //*/
    /*
    // Soft-Clip Tail - Boundary Case
    dpWork->regionStartIdx = 0;
    dpWork->regionLength = 11;
    char * pattern = "TAACCCGGGG";
    
    dpWork->leadSoftClipLength = 0;
    dpWork->tailSoftClipLength = 2;
    dpWork->totalSoftClipLength = 100;
    //*/
    /*
    // Soft-Clip Both - Boundary Case
    dpWork->regionStartIdx = 0;
    dpWork->regionLength = 11;
    char * pattern = "CCACCCGGGG";
    
    dpWork->leadSoftClipLength = 2;
    dpWork->tailSoftClipLength = 2;
    dpWork->totalSoftClipLength = 100;
    //*/
    /*
    // Hard-Clip Head
    dpWork->regionStartIdx = 0;
    dpWork->regionLength = 11;
    char * pattern = "AAACCCTTAA";
    
    dpWork->leadSoftClipLength = 0;
    dpWork->tailSoftClipLength = 0;
    dpWork->totalSoftClipLength = 0;
    
    dpWork->leadHardClipLength = 2;
    dpWork->tailHardClipLength = 0;
    //*/
    
    // Tail-Clip Tail
    dpWork->regionStartIdx = 0;
    dpWork->regionLength = 11;
    char * pattern = "AAACCCTTAA";
    
    dpWork->leadSoftClipLength = 0;
    dpWork->tailSoftClipLength = 0;
    dpWork->totalSoftClipLength = 0;
    
    dpWork->leadHardClipLength = 0;
    dpWork->tailHardClipLength = 5;
    //*/
    
    dpWork->readCode = (unsigned char*) malloc(sizeof(unsigned char)*(strlen(pattern)+1));
    int i;
    for (i=0;i<strlen(pattern);i++) {
        if (pattern[i]=='A') dpWork->readCode[i] = 0;
        if (pattern[i]=='C') dpWork->readCode[i] = 1;
        if (pattern[i]=='G') dpWork->readCode[i] = 2;
        if (pattern[i]=='T') dpWork->readCode[i] = 3;
    } pattern[i]='\0';
    dpWork->readLength = i;
    
    // Testing DP-Tweak
    #ifdef DEV_TEST_DP_TEST_TAIL_VERIFY_LENGTH
        dpWork->tailVerifyLength = 4;
    #endif
    
    DPMatrixFetching(dpWork,hsp);
    printf("DP Work Status = %d\n",dpWork->flag);
    
    printf("testDp -- begining to fill the matrix - DPMatrixFill\n");
    DPMatrixFill(dpArgument);
    printf("DP Work Status = %d\n",dpWork->flag);
    
    printf("testDp -- Debug Print - DPDebugPrintMatrix\n");
    DPDebugPrintMatrix(dpWork);
    printf("DP Work Status = %d\n",dpWork->flag);
    
    printf("testDp -- BackTrack - DPBackTrack\n");
    DPBackTrack(dpWork,dpArgument->dpScores,-5,&dpOcc_found);
    printf("DP Work Status = %d\n",dpWork->flag);

    printf("testDp -- output occurrence\n");
    DevTestDpPrintDpOcc(&dpOcc_found);
    
    printf("DP Work Status = %d\n",dpWork->flag);
    free(dpWork->readCode);
    
    DPWorkFree(dpWork);
    
    printf("testDp -- terminating\n");
    
    //=======================
    HSPFree(mmPool, hsp, 1);
    MMPoolFree(mmPool);
}
