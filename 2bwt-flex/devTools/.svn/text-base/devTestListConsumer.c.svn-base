

#include "../ListConsumer.h"
#include "../2bwt-lib/TypeNLimit.h"

int main () {
    char readFilename[1025];
    char mateFilename[1025];
    int ub, lb;
    char outputPrefix[1025];
    char readGroup[1025];
    char sampleName[1025];
    char readGrpOption[1025];
    
    ListConsumer * listConsumer = LCCreate();
    LCLoadList(listConsumer,"devTestListConsumer.list",MAX_FILENAME_LEN,MAX_FIELD_LEN);
    
    int i = 0;
    while (LCGetNextPairEndRead(listConsumer,readFilename,mateFilename,&lb,&ub,outputPrefix,readGroup,sampleName,readGrpOption)) {
        printf("%d-th line = %s / %s / %d / %d\n",i,readFilename,mateFilename,lb,ub);
        i++;
    }
    printf("Reading of the %d-th line... FAILED\n",i);
    printf("JUST PLAYING..\n");
    LCGetNextPairEndRead(listConsumer,readFilename,mateFilename,&lb,&ub,outputPrefix,readGroup,sampleName,readGrpOption);
    
    LCFree(listConsumer);
}
