#include "../2bwt-flex/OutputSAMFormat.h"
#include "stdio.h"
#include "stdlib.h"

char samFlagTranslate[12][1024] = { 
    "0x1 template having multiple segments in sequencing",
    "0x2 each segment properly aligned according to the aligner",
    "0x4 segment unmapped",
    "0x8 next segment in the template unmapped",
    "0x10 SEQ being reverse complemented",
    "0x20 SEQ of the next segment in the template being reversed",
    "0x40 the first segment in the template",
    "0x80 the last segment in the template",
    "0x100 secondary alignment",
    "0x200 not passing quality controls",
    "0x400 PCR or optical duplicate",
    "0x800 supplementary alignment" };

int main() {
    int flag;
    
    while (1) {
        printf("Please enter a SAM flag = ");fflush(stdout);
        scanf("%d",&flag);
        
        int msgIdx = 0;
        int flagIdx = 1;
        while (flagIdx <= flag) {
            if ( (flag & flagIdx) > 0 ) {
                printf("%9d %s\n",flagIdx,samFlagTranslate[msgIdx]);
            }
            flagIdx *= 2;
            msgIdx++;
        }
    }
    return 0;
}