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

#include "../SRAQueryParser.h"
#include "../2BWT-SRAAlgnmt.h"

#include <string.h>

//#define PRINTBWT_OPTION_PRINT_HEADER
//#define PRINTBWT_OPTION_PRINT_BWT
//#define PRINTBWT_OPTION_PRINT_OCC_DEBUG
//#define PRINTBWT_OPTION_PRINT_OCC_VALUE
//#define PRINTBWT_OPTION_PRINT_ALL_BWT_DECODE
//#define PRINTBWT_OPTION_PRINT_OCC_VALUE_MAJOR

int _BWTBackwardSearch(BWT * bwt, unsigned long long rawPattern, int len, unsigned long long * saL, unsigned long long * saR, char * pattern) {

    unsigned long long l = 0;
    unsigned long long r = bwt->textLength;
    int charMask = ((1<<LOOKUP_BIT_PER_CHAR)-1);
    
    int i = 0;
    pattern[len]='\0';
    while (i<len && l<=r) {
        unsigned char c = rawPattern & charMask;
        pattern[len-i-1]=dnaChar[c];

        l = bwt->cumulativeFreq[c] + BWTOccValue(bwt, l, c) + 1;
        r = bwt->cumulativeFreq[c] + BWTOccValue(bwt, r + 1, c);
        
        i++;rawPattern>>=LOOKUP_BIT_PER_CHAR;
    }
    
    (*saL) = l;
    (*saR) = r;
    
    return (l <= r);
}

int main(int argc, char** argv) {

    BWT *bwt;
    MMPool *mmPool;
    LT * lookup;

    int PoolSize = 2097152;                // 2M  - fixed; not configurable through ini
    MMMasterInitialize(3, 0, FALSE, NULL);
    mmPool = MMPoolCreate(PoolSize);

    char * indexPrefix = argv[1];
    char bwtFile[1024];
    char fmvFile[1024];
    char ltFile[1024];
    long long i,j;
    unsigned char c;

    printf("Index prefix                 = %s\n",indexPrefix);
    strcpy(bwtFile,indexPrefix);
    strcat(bwtFile,".bwt");
    strcpy(fmvFile,indexPrefix);
    strcat(fmvFile,".fmv");
    strcpy(ltFile,indexPrefix);
    strcat(ltFile,".lkt");
    printf("BWT file                     = %s\n",bwtFile);
    printf("FMV file                     = %s\n",fmvFile);
    printf("LT file                      = %s\n",ltFile);

    bwt = BWTLoad(mmPool, bwtFile, fmvFile, NULL, NULL, NULL, NULL);
    lookup = LTLoad(ltFile);
    
    printf("lookup->tableSize            = %d\n",lookup->tableSize);
    printf("lookup->ltSizeInWord         = %d\n",lookup->ltSizeInWord);
    
    #ifdef PRINTBWT_OPTION_PRINT_HEADER
    printf("\n");
    printf("bwt->inverseSa0              = %20llu\n",(unsigned long long) bwt->inverseSa0);
    printf("bwt->textLength              = %20llu\n",(unsigned long long) bwt->textLength);
    printf("bwt->saInterval              = %20llu\n",(unsigned long long) bwt->saInterval);
    printf("bwt->inverseSaInterval       = %20llu\n",(unsigned long long) bwt->inverseSaInterval);
    printf("bwt->decodeTableGenerated    = %20d\n",(int) bwt->decodeTableGenerated);
    printf("bwt->bwtSizeInWord           = %20llu\n",(unsigned long long) bwt->bwtSizeInWord);
    printf("bwt->occSizeInWord           = %20llu\n",(unsigned long long) bwt->occSizeInWord);
    printf("bwt->occMajorSizeInWord      = %20llu\n",(unsigned long long) bwt->occMajorSizeInWord);
    printf("bwt->saValueSizeInWord       = %20llu\n",(unsigned long long) bwt->saValueSizeInWord);
    printf("bwt->inverseSaSizeInWord     = %20llu\n",(unsigned long long) bwt->inverseSaSizeInWord);
    printf("bwt->cachedSaIndexSizeInWord = %20llu\n",(unsigned long long) bwt->cachedSaIndexSizeInWord);
    for (i=0;i<ALPHABET_SIZE;i++) {
        printf("bwt->cumulativeFreq[%lld]      = %20llu\n",i,(unsigned long long)bwt->cumulativeFreq[i]);
    }
    #endif
    
    unsigned long long rawPattern;
    unsigned long long l,r;
    
    char pattern[1024];
    char verdict[3] = " !";
    for (rawPattern=0;rawPattern<lookup->ltSizeInWord;rawPattern++) {
        if (_BWTBackwardSearch(bwt,rawPattern,lookup->tableSize,&l,&r,pattern)) {
            printf("Patt. = %s      ",pattern);
            printf("LT = %llu            BWT = %llu (%llu)    %c\n",(unsigned long long)lookup->table[rawPattern],l,r,verdict[r!=lookup->table[rawPattern]]);
        }
    }
    
    #ifdef PRINTBWT_OPTION_PRINT_OCC_DEBUG
    #define POS_ARRAY_SIZE 16
    printf("\n");
    unsigned long long oL[ALPHABET_SIZE];
    static const long long posArray[POS_ARRAY_SIZE] = { 0, 1, 128, 256, 259, 300,
                                                        512, 1024, 2048, 3000, 30000,
                                                        150000, 159918, 159919, 159920, 159921 };
    for (i=0;i<POS_ARRAY_SIZE;i++) {
        unsigned long long pos = posArray[i];
        BWTAllOccValue(bwt,pos,oL);
        for (c=0; c<ALPHABET_SIZE; c++) {
            printf("BWTOccValue(%15llu,%d)    = %15lld / %15lld (ALL)\n",pos,c,(unsigned long long)BWTOccValue(bwt,pos,c),(unsigned long long)oL[c]);
        }
    }
    #endif
    
    #ifdef PRINTBWT_OPTION_PRINT_ALL_BWT_DECODE
    printf("\n");
    for (i=0;i<10240;i++) {
        unsigned long long pos = i;
        BWTAllOccValue(bwt,pos,oL);
        for (c=0; c<ALPHABET_SIZE; c++) {
            printf("BWTOccValue(%15llu,%d)    = %15lld / %15lld (ALL)\n",pos,c,(unsigned long long)BWTOccValue(bwt,pos,c),(unsigned long long)oL[c]);
        }
    }
    #endif

    #ifdef PRINTBWT_OPTION_PRINT_BWT
    printf("\n");
    printf("BWT Code Value Debug =========================== \n");
    for (i=0;i<bwt->textLength/CHAR_PER_WORD;i++) {
        printf("%c",dnaChar[(bwt->bwtCode[i]>>28)&15]);
        printf("%c",dnaChar[(bwt->bwtCode[i]>>24)&15]);
        printf("%c",dnaChar[(bwt->bwtCode[i]>>20)&15]);
        printf("%c",dnaChar[(bwt->bwtCode[i]>>16)&15]);
        
        printf("%c",dnaChar[(bwt->bwtCode[i]>>12)&15]);
        printf("%c",dnaChar[(bwt->bwtCode[i]>>8)&15]);
        printf("%c",dnaChar[(bwt->bwtCode[i]>>4)&15]);
        printf("%c",dnaChar[(bwt->bwtCode[i])&15]);
        
        if (i%10==9) {printf("\n");}
    }
    #endif

    #ifdef PRINTBWT_OPTION_PRINT_OCC_VALUE
    printf("\n");
    printf("BWT Occ Value Debug =========================== \n");
    printf("NOTE: Alignment of this values are expected to be different between 32-bit and 64-bit BWT index\n");
    for (i=0;i<bwt->occSizeInWord/(ALPHABET_SIZE*2);i++) {
        for (j=0;j<ALPHABET_SIZE;j++) {
            printf("%15llu",(unsigned long long)bwt->occValue[i*(ALPHABET_SIZE*2)+j*2]);
            if (j%4==3) {printf("\n");}
        }
        for (j=0;j<ALPHABET_SIZE;j++) {
            printf("%15llu",(unsigned long long)bwt->occValue[i*(ALPHABET_SIZE*2)+j*2+1]);
            if (j%4==3) {printf("\n");}
        }
    }
    printf("\n");
    #endif

    #ifdef PRINTBWT_OPTION_PRINT_OCC_VALUE_MAJOR
    printf("\n");
    printf("BWT Occ Major Value Debug =========================== \n");
    for (i=0;i<bwt->occMajorSizeInWord;i++) {
        printf("%15llu",(unsigned long long)bwt->occValueMajor[i]);
        if (i%5==4) {printf("\n");}
    }
    printf("\n");
    #endif

    BWTFree(mmPool, bwt);
    MMPoolFree(mmPool);
    LTFree(lookup);

    return 0;
}
