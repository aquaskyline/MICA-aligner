
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>

#include "../MIC-MEMControl.h"

#define SIZE 1000000000

void SetupBuffer(char * buffer, int factor, int size) {
    int i;
    for (i=0;i<size;i++) {
        buffer[i] = i*factor;
    }
}

__attribute__ (( target (mic)))
void PrintBuffer(char * buffer, int size) {
    int i;
    for (i=0;i<size;i++) {
        printf("%5d ",buffer[i]);
        if ( i % 5 == 4 ) { printf("\n"); }
    }
}

void BruteForceAllocation() {

    char c;
    printf("Allocating CPU host memory..");fflush(stdout);
    char * buffer1 = malloc(SIZE);
    char * micBuffer1;
    char * micBuffer2 = malloc(1);
    SetupBuffer(buffer1,1,5);

    printf("DONE\nHit return to continue.");
    scanf("%c",&c);

    printf("Allocating MIC device memory..");fflush(stdout);
    #pragma offload target(mic : 0) mandatory \
        nocopy(micBuffer1) \
        in(micBuffer2:length(1) alloc_if(1) free_if(0))
    { 
        micBuffer1 = (char*) malloc(SIZE);
        micBuffer2=micBuffer1;
        printf("micBuffer1 = %llu\n",(unsigned long long) micBuffer1);fflush(0);
        printf("micBuffer2 = %llu\n",(unsigned long long) micBuffer2);fflush(0);
    }
    printf("DONE\nHit return to continue.");
    scanf("%c",&c);
    
    printf("Offloading MIC device memory #1..");fflush(stdout);
    #pragma offload target(mic : 0) mandatory \
        nocopy(micBuffer1) \
        nocopy(micBuffer2:length(1) alloc_if(0) free_if(0))
    { 
        printf("micBuffer2 = %llu\n",(unsigned long long) micBuffer2);fflush(0);
    }
    printf("DONE\nHit return to continue.");
    scanf("%c",&c);
    
    printf("DONE\nHit return to continue.");
    scanf("%c",&c);
    printf("Offloading MIC device memory #3..");fflush(stdout);
    #pragma offload target(mic : 0) mandatory \
        in(buffer1:length(SIZE) into(micBuffer2) alloc_if(0))
    { PrintBuffer(buffer1,5); }
    printf("DONE\nHit return to continue.");
    scanf("%c",&c);

    printf("Free-ing CPU host memory..");fflush(stdout);
    free(buffer1);
    printf("DONE\nHit return to continue.");
    scanf("%c",&c);

}


void MICMemAllocation() {
    MICMemBlockArray * micMem = MICMemCreate(0);
    int i;char c;
    
    int size = 100000000;
    
    // Create Host Memory Block
    uint64_t * testU64 = malloc(sizeof(uint64_t) * size);
    for (i=0;i<1024;i++) {
        testU64[i] = i;
    }
    
    // Create Host Memory Block
    uint32_t * testU32 = malloc(sizeof(uint32_t) * size);
    for (i=0;i<1024;i++) {
        testU32[i] = i*10;
    }
    
    // MIC Memory Pointer
    void * micPtr1 = NULL;
    void * micPtr2 = NULL;
    printf("CPU micPtr1 : %llu\n", micPtr1);fflush(0);
    printf("CPU micPtr2 : %llu\n", micPtr2);fflush(0);
    
    micPtr1 = MICMemAddBlock(micMem,(void*)testU64,sizeof(uint64_t),size,"TEST_U64");
    printf("Successfully added memory block #1 into MIC Memory Cache\n");
    MICMemDebugPrint(micMem);

        printf("CPU %llu\n",micPtr1);
    micPtr2 = MICMemAddBlock(micMem,(void*)testU32,sizeof(uint32_t),size,"TEST_U32");
    printf("Successfully added memory block #2 into MIC Memory Cache\n");
    MICMemDebugPrint(micMem);
    
    printf("Offloaded host memory..DONE\n");fflush(stdout);
    printf("Hit return to continue.");fflush(stdout);
    scanf("%c",&c);
    
    printf("CPU micPtr1 : %llu\n", micPtr1);fflush(0);
    printf("CPU micPtr2 : %llu\n", micPtr2);fflush(0);
    //printf("Offloading host memory..");fflush(stdout);
    //MICMemLoad(micMem,0);
    //printf("DONE\nHit return to continue.");fflush(stdout);
    //scanf("%c",&c);
    

////// TEMP HACK BEGIN
/*MICMemBlock * block = &(micMem->blocks[0]);
#pragma offload target(mic : 0) \
                inout(block:length(1) alloc_if(1) free_if(1))
{
    printf("nbytes : %llu\n",block->micPtr);
    int j;
    uint64_t * m = (uint64_t * ) block->micPtr;
    for (j=0;j<100;j++) {
        printf("%llu ",m[j]);
    }printf("\n");
    fflush(0);
}*/
////// TEMP HACK END

////// TEMP HACK BEGIN
#pragma offload target(mic : 0)
{
    printf("nbytes : %llu\n",micPtr1);
    int j;
    uint64_t * m = (uint64_t *) micPtr1;
    for (j=0;j<1024;j++) {
        printf("%llu ",m[j]);
    }printf("\n");
    fflush(0);
}
////// TEMP HACK END

////// TEMP HACK BEGIN
#pragma offload target(mic : 0)
{
    printf("nbytes : %llu\n",micPtr2);
    int j;
    uint32_t * m = (uint32_t *) micPtr2;
    for (j=0;j<1024;j++) {
        printf("%u ",m[j]);
    }printf("\n");
    fflush(0);
}
////// TEMP HACK END

    printf("Free MIC memory..");fflush(stdout);
    MICMemFree(micMem);
    printf("DONE\nHit return to continue.");fflush(stdout);
    scanf("%c",&c);
}

int main() {

    if (0) BruteForceAllocation();
    
    MICMemAllocation();

    return 0;

}
