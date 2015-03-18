
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>

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
#include "../SRA2BWTOperations.h"

#include "../PEArguments.h"
#include "../PE2BWTReport.h"
#include "../DPOperations.h"

void printSize(char *name, int sizeInByte) {
    printf(" %30s",name);
    printf("       = %d Byte / %d KByte / %d MByte",
            sizeInByte, 
            (int)((double)(sizeInByte + 1023.0)/1024.0), 
            (int)((double)(sizeInByte + 1023.0)/1024.0/1024.0));
    printf("\n");
}


int main() {
    printSize("char",(int)sizeof(char));
    printf("\n");
    printSize("SRAError",(int)sizeof(SRAError));
    printSize("SRAOccurrence",(int)sizeof(SRAOccurrence));
    printf("\n");
    
    printSize("SRAArguments",(int)sizeof(SRAArguments));
    printSize("SRAQueryInfo",(int)sizeof(SRAQueryInfo));
    printf("\n");
    
    printSize("SRAOCCCollector",(int)sizeof(SRAOCCCollector));
    printSize("SRAOccurrence_ll",(int)sizeof(SRAOccurrence_ll));
    printf("\n");
    
    printSize("DPWork",(int)sizeof(DPWork));
    printSize("DPOccurrence",(int)sizeof(DPOccurrence));
    printf("\n");
    
    printf("Constants\n");
    printf("SRA_MAX_READ_LENGTH = %u\n",SRA_MAX_READ_LENGTH);
    
    printf("DP_MAX_REGION_LENGTH = %u\n",DP_MAX_REGION_LENGTH);
    printf("DP_MAX_ERROR_COUNT = %u\n",DP_MAX_ERROR_COUNT);
    
    return 0;
}
