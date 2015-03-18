
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>

#define SIZE 1000000000


__attribute__ (( target (mic)))
void SetupBuffer(char * buffer, int factor, unsigned long long size) {
    unsigned long long i;
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
__attribute__((target(mic)))
char * testPassing;

void Offload1(char * buffer1) {
    #pragma offload target(mic : 0) mandatory \
        in(buffer1:length(SIZE)                alloc_if(1) free_if(1) align(64)) \
        nocopy(testPassing:length(1))
    {
        testPassing = malloc(sizeof(char) * SIZE);
        memcpy(testPassing,buffer1,SIZE);
        
        printf("%u\n",testPassing[0]);
        printf("%u\n",testPassing[1]);
        printf("%u\n",testPassing[2]);
        printf("%u\n",testPassing[3]);
        printf("%u\n",testPassing[4]);fflush(0);
    }
}

void Offload2(char * buffer1) {
    #pragma offload target(mic : 0) mandatory \
        nocopy(testPassing:length(1))
    {
        printf("%u\n",testPassing[0]);
        printf("%u\n",testPassing[1]);
        printf("%u\n",testPassing[2]);
        printf("%u\n",testPassing[3]);
        printf("%u\n",testPassing[4]);fflush(0);
    }
}


int main() {

    char c;

    printf("Allocating CPU host memory..");fflush(stdout);
    char * buffer1 = malloc(SIZE);
    char * buffer2 = malloc(SIZE*2);
    char * tmp;
    SetupBuffer(buffer1,1,5);

    printf("DONE\nHit return to continue.");
    scanf("%c",&c);
    
    //Offload1(buffer1);
    Offload2(buffer1);
    

    printf("Transfer MIC device memory.. #1");fflush(stdout);
    #pragma offload target(mic : 0) mandatory \
        in(buffer1:length(SIZE)                alloc_if(1) free_if(1) align(64))
    { }
    printf("DONE\nHit return to continue.");
    scanf("%c",&c);

    printf("Transfer MIC device memory.. #2");fflush(stdout);
    #pragma offload target(mic : 0) mandatory \
        in(buffer2:length(SIZE*2)                alloc_if(1) free_if(1))
    { }
    printf("DONE\nHit return to continue.");
    scanf("%c",&c);

    printf("Transfer MIC device memory.. #1");fflush(stdout);
    #pragma offload target(mic : 0) mandatory \
        in(buffer1:length(SIZE)                alloc_if(1) free_if(1))
    { }
    printf("DONE\nHit return to continue.");
    scanf("%c",&c);

    printf("Free-ing CPU host memory..");fflush(stdout);
    free(buffer1);
    free(buffer2);
    printf("DONE\nHit return to continue.");
    scanf("%c",&c);

    return 0;

}
