//
//    devExtractSeq.c
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

#include "devExtractSeq.h"


char getRealChar(unsigned char c) {
    if (c==0) return 'A';
    if (c==1) return 'C';
    if (c==2) return 'G';
    if (c==3) return 'T';
    return '?';
}

unsigned char getCharInWord(unsigned int word, int offset) {
    word>>=(15-offset)*2;
    word &= 3;
    return (unsigned char)word;    
}

int main (int argc, char ** argv) {
	FILE *inputFile;
	unsigned char tempChar[4];
	unsigned int writeBuffer=0;
	unsigned int writeBufferIndex=0;
	unsigned int *packedText;
	unsigned int packedFileLen;
	unsigned char lastByteLength;
	unsigned int wordToProcess;
    unsigned int mk2,textLength;
    char seq[32];
	int i,j,k;
	int queryNum=1,queryLen=35;
    
    if (argc!=3) {
        printf("%s <PackedDNAFileName> <FASTASequenceFileName>\n",argv[0]);
        exit(1);
    }

	inputFile = (FILE*)(FILE*)fopen64(argv[1], "rb");

	if (inputFile == NULL) {
		fprintf(stderr, "DNALoadPacked() : Cannot open inputFileName!\n");
		exit(1);
	}

	fseek(inputFile, -1, SEEK_END);
	packedFileLen = ftell(inputFile);
	if ((int)packedFileLen < 0) {
		fprintf(stderr, "DNALoadPacked(): Cannot determine file length!\n");
		exit(1);
	}
	fread(&lastByteLength, sizeof(unsigned char), 1, inputFile);

	textLength = (packedFileLen - 1) * 4 + lastByteLength;

	wordToProcess = (textLength + 16 - 1) / 16;
	fprintf(stderr,"Allocating memory..\n");
	packedText = (unsigned int *) malloc((wordToProcess + 1) * sizeof(unsigned int));	// allocate 1 more word at end
	fprintf(stderr,"Allocated memory..\n");
	packedText[wordToProcess - 1] = 0;
	packedText[wordToProcess] = 0;

	fseek(inputFile, 0, SEEK_SET);
	fread(packedText, 1, packedFileLen, inputFile);
	
	for (i=0; i<wordToProcess; i++) {
		*(unsigned int*)tempChar = packedText[i];
		packedText[i] = (tempChar[0] << 24) | (tempChar[1] << 16) | (tempChar[2] << 8) | tempChar[3];
    }
    
    
	FILE * fin = fopen(argv[2], "r");
	fseek(fin, -1, SEEK_END);
	unsigned int fileLen = ftell(fin);
	fprintf(stderr,"Allocating memory..%u\n",fileLen);
	char * fullText = (char *) malloc(fileLen * sizeof(char));
	fprintf(stderr,"Allocated memory..\n");
	unsigned int seqSP[25];
	unsigned int seqCount=0;
	unsigned int count = 0;
	char c;
	seqSP[0]=0;
	
	fseek(fin, 0, SEEK_SET);
	fprintf(stderr,"Reading sequence file..\n");
	while (!feof(fin)) {
		c = getc(fin);
		if (c=='>') {
			while (!feof(fin) && c != '\n') {
				c = getc(fin);
			}
			seqSP[seqCount]=count;
			seqCount++;
		}
		if (c!='\n' && c!='\t') {
            fullText[count]=c;
			count++;
			/*if (count>tp) {
				printf("%c",c);
			}
			if (count>tp+len-1) {
				
				printf("\nChromosome = %u\n",seqCount);
			}	*/
		}
	}
	seqSP[seqCount]=count;
	seqCount++;
    
    while (1) {
        unsigned int word,len,location;
        printf("By Chromosome? [1/0] ");int c;scanf("%d",&c);
        if (c==0) {
            printf("[1-based] Location on Packed text?[0 to disable] ");scanf("%u",&word);
            printf("[1-based] Location on (unambiguous)Full text?[0 to disable] ");scanf("%u",&location);
            printf("Length of substring? ");scanf("%u",&len);
	        printf("Positive[1]/Negative[0] strand? ");int d;scanf("%d",&d);
            if (word!=0) {
				word--;
				unsigned int tp = word+1;
                unsigned int offset = word % 16;
                word/=16;
                j=0;
            	printf("%u\t",tp);
                while (j<len) {
					char c = getRealChar(getCharInWord(packedText[word],offset));
					
					if (d==1) {
		                printf("%c",c);
					} else {
						if (c=='A') {
							printf("T");
						} else if (c=='C') {
							printf("G");
						} else if (c=='G') {
							printf("C");
						} else if (c=='T') {
							printf("A");
						}
		                
					}
                    j++;
                    offset++;
                    if (offset>=16) {offset=0;word++;}
                }
            	printf("\t%u",tp+len);
                printf("\n");
            }
            
            if (location!=0) {
				location--;
                unsigned int i;
            	printf("%u\t",location+1);
                for (i=location;i<location+len;i++) {
					if (d==1) {
		                printf("%c",fullText[i]);
					} else {
						if (fullText[i]=='A') {
							printf("T");
						} else if (fullText[i]=='C') {
							printf("G");
						} else if (fullText[i]=='G') {
							printf("C");
						} else if (fullText[i]=='T') {
							printf("A");
						}
		                
					}
                    //printf("%c",fullText[i]);
                }
            	printf("\t%u",location+len);
                printf("\n");
            }
        } else {
            printf("Input chromosome#? ");scanf("%u",&word);
            printf("[1-based] Substring location? ");scanf("%u",&location);location--;
            printf("Length of substring? ");scanf("%u",&len);
	        printf("Positive[1]/Negative[0] strand? ");int d;scanf("%d",&d);
            unsigned int i;
            location+=seqSP[word-1];
            printf("%u(%u)\t",location+1,seqSP[word-1]);
            for (i=location;i<location+len;i++) {
				if (d==1) {
	                printf("%c",fullText[i]);
				} else {
					if (fullText[i]=='A') {
						printf("T");
					} else if (fullText[i]=='C') {
						printf("G");
					} else if (fullText[i]=='G') {
						printf("C");
					} else if (fullText[i]=='T') {
						printf("A");
					}
				}
            }
            printf("\t%u",location+len);
            printf("\n");
        }
    }
    /*fprintf(stderr,"Generating Queries\n");
	for (i=0;i<queryNum;i++) {
        printf("> Test Queries Number = %u\n",i);
        unsigned int word = genRand(textLength/16-1);
        unsigned int offset = rand() % 16;
        j=0;
        while (j<queryLen) {
            printf("%c", getRealChar(getCharInWord(packedText[word],offset)));
            j++;
            offset++;
            if (offset>=16) {offset=0;word++;}
            
        }
        printf("\n");
    }*/
	
	fclose(fin);
	free(fullText);
	
	free(packedText);
	fclose(inputFile);
	
    fprintf(stderr,"textLength = %u\n",textLength);
}
