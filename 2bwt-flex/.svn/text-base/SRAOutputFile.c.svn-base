//
//    SRAOutputFile.c
//
//    mica
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

#include "SRAOutputFile.h"

SRAOutput * SRAAlgnmtOutputConstruct() {

    SRAOutput * sraOutput = (SRAOutput*) malloc (sizeof(SRAOutput));
    
    // Common:
    sraOutput->OutFileFormat = SRA_OUTPUT_FORMAT_UNDEF;
    sraOutput->OutFileName[0] = '\0';
    sraOutput->occCollector = NULL;

    // For Plain output usage:
    sraOutput->OutFilePtr = NULL;
    
    // For SAM output usage:
    sraOutput->SAMOutFilePtr = NULL;
    sraOutput->SAMAuxDataBlockSize = 0;
    sraOutput->SAMAuxDataBlock = NULL;
    sraOutput->ReadGroup[0]='\0';
    return sraOutput;
}

void SRAAlgnmtOutputFree(SRAOutput * sraOutput) {

    if (sraOutput->OutFileFormat != SRA_OUTPUT_FORMAT_UNDEF) {
        printf("[SAFEGUARD] SRAAlgnmtOutputFree is destroying a SRAOutput block that is occupied by output type %u.\n", sraOutput->OutFileFormat);
    }
    
    if (sraOutput->SAMOutFilePtr != NULL) {
        samclose(sraOutput->SAMOutFilePtr);
    }
    
    if (sraOutput->OutFilePtr != NULL) {
        fclose(sraOutput->OutFilePtr);
    }

    if (sraOutput->occCollector != NULL) {
        SRAOCCFree(sraOutput->occCollector);
    }

    if (sraOutput->SAMAuxDataBlock != NULL) {
        free(sraOutput->SAMAuxDataBlock);
    }
    
    free(sraOutput);
}

//-------------------------------------------------
// Output Format dictated procedures - INITIALISATION
//-------------------------------------------------
void SRAAlgnmtOutputInitNone(SRAOutput * sraOutput) {

    if (sraOutput->OutFileFormat == SRA_OUTPUT_FORMAT_UNDEF) {
    
        // Setup the status
        sraOutput->OutFileFormat = SRA_OUTPUT_FORMAT_NONE;
        
        // Setup the necessary auxiliary data structure
        sraOutput->occCollector = SRAOCCCreate(SRAOCC_FLOOD_TYPE_EXPAND);
        
    } else {
        printf("[SAFEGUARD] SRAAlgnmtOutputInitNone is initialising SRAOutput block that is already occupied by output type %u.\n", sraOutput->OutFileFormat);
        printf("Function is invoked for output type:%u\n",SRA_OUTPUT_FORMAT_NONE);
    }
}

void SRAAlgnmtOutputInitDiscard(SRAOutput * sraOutput) {

    if (sraOutput->OutFileFormat == SRA_OUTPUT_FORMAT_UNDEF) {
    
        // Setup the status
        sraOutput->OutFileFormat = SRA_OUTPUT_FORMAT_DISCARD;
        
        // Setup the necessary auxiliary data structure
        sraOutput->occCollector = SRAOCCCreate(SRAOCC_FLOOD_TYPE_FLUSH);
        
    } else {
        printf("[SAFEGUARD] SRAAlgnmtOutputInitDiscard is initialising SRAOutput block that is already occupied by output type %u.\n", sraOutput->OutFileFormat);
        printf("Function is invoked for output type:%u\n",SRA_OUTPUT_FORMAT_DISCARD);
    }
}

void SRAAlgnmtOutputInitPlain(SRAOutput * sraOutput, char * OutFileName) {

    if (sraOutput->OutFileFormat == SRA_OUTPUT_FORMAT_UNDEF) {
    
        // Setup the status
        sraOutput->OutFileFormat = SRA_OUTPUT_FORMAT_PLAIN;
        strcpy(sraOutput->OutFileName,OutFileName);
        
        // Setup the necessary auxiliary data structure
        sraOutput->occCollector = SRAOCCCreate(SRAOCC_FLOOD_TYPE_EXPAND);
        sraOutput->OutFilePtr = (FILE*)fopen64(OutFileName, "w");
        
        printf("Preparing Output File = %s\n", OutFileName);
    } else {
        printf("[SAFEGUARD] SRAAlgnmtOutputInitPlain is initialising SRAOutput block that is already occupied by output type %u.\n", sraOutput->OutFileFormat);
        printf("Function is invoked for output type:%u filename:%s\n",SRA_OUTPUT_FORMAT_PLAIN,OutFileName);
    }
}

void SRAAlgnmtOutputInitSAMStore(SRAOutput * sraOutput) {

    if (sraOutput->OutFileFormat == SRA_OUTPUT_FORMAT_UNDEF) {
    
        // Setup the status
        sraOutput->OutFileFormat = SRA_OUTPUT_FORMAT_SAM_STORE;
        
        // Setup the necessary auxiliary data structure
        sraOutput->occCollector = SRAOCCCreate(SRAOCC_FLOOD_TYPE_EXPAND);
        
    } else {
        printf("[SAFEGUARD] SRAAlgnmtOutputInitSAMStore is initialising SRAOutput block that is already occupied by output type %u.\n", sraOutput->OutFileFormat);
        printf("Function is invoked for output type:%u\n",SRA_OUTPUT_FORMAT_SAM_STORE);
    }
}

void SRAAlgnmtOutputInitSAM(SRAOutput * sraOutput, char * OutFileName, bam_header_t * samOutputHeader, char * readGroup) {

    if (sraOutput->OutFileFormat == SRA_OUTPUT_FORMAT_UNDEF) {
    
        // Setup the status
        sraOutput->OutFileFormat = SRA_OUTPUT_FORMAT_SAM;
        strcpy(sraOutput->OutFileName,OutFileName);
        
        // Setup the necessary auxiliary data structure
        sraOutput->occCollector = SRAOCCCreate(SRAOCC_FLOOD_TYPE_EXPAND);
        sraOutput->SAMOutFilePtr = samopen(OutFileName,"wh",samOutputHeader);
        sraOutput->SAMAuxDataBlockSize = SAM_MDATA_SIZE_PER_READ + SAM_MDATA_INITIAL_ALLOWED_OCC * SAM_MDATA_SIZE_PER_OCC;
        sraOutput->SAMAuxDataBlock = (uint8_t *)malloc(sizeof(uint8_t)*sraOutput->SAMAuxDataBlockSize);

        strcpy(sraOutput->ReadGroup,readGroup);
        printf("Preparing Output File = %s\n", OutFileName);
    } else {
        printf("[SAFEGUARD] SRAAlgnmtOutputInitSAM is initialising SRAOutput block that is already occupied by output type %u.\n", sraOutput->OutFileFormat);
        printf("Function is invoked for output type:%u filename:%s\n",SRA_OUTPUT_FORMAT_SAM,OutFileName);
    }
}

void SRAAlgnmtOutputInitBAM(SRAOutput * sraOutput, char * OutFileName, bam_header_t * samOutputHeader, int numthreads, char * readGroup) {

    if (sraOutput->OutFileFormat == SRA_OUTPUT_FORMAT_UNDEF) {
    
        // Setup the status
        sraOutput->OutFileFormat = SRA_OUTPUT_FORMAT_BAM;
        strcpy(sraOutput->OutFileName,OutFileName);
        sraOutput->ReadGroup[0]='\0';
        // Setup the necessary auxiliary data structure
        sraOutput->occCollector = SRAOCCCreate(SRAOCC_FLOOD_TYPE_EXPAND);
        sraOutput->SAMOutFilePtr = samopen(OutFileName,"wb",samOutputHeader);
        samthreads(sraOutput->SAMOutFilePtr, numthreads, 256);
        sraOutput->SAMAuxDataBlockSize = SAM_MDATA_SIZE_PER_READ + SAM_MDATA_INITIAL_ALLOWED_OCC * SAM_MDATA_SIZE_PER_OCC;
        sraOutput->SAMAuxDataBlock = (uint8_t *)malloc(sizeof(uint8_t)*sraOutput->SAMAuxDataBlockSize);
        strcpy(sraOutput->ReadGroup,readGroup);
        printf("Preparing Output File = %s\n", OutFileName);
    } else {
        printf("[SAFEGUARD] SRAAlgnmtOutputInitBAM is initialising SRAOutput block that is already occupied by output type %u.\n", sraOutput->OutFileFormat);
        printf("Function is invoked for output type:%u filename:%s\n",SRA_OUTPUT_FORMAT_BAM,OutFileName);
    }
}

//-------------------------------------------------
// Output Format dictated procedures - FREE
//-------------------------------------------------
void SRAAlgnmtOutputFreeNone(SRAOutput * sraOutput) {

    if (sraOutput->OutFileFormat == SRA_OUTPUT_FORMAT_NONE) {
    
        // Setup the status
        sraOutput->OutFileFormat = SRA_OUTPUT_FORMAT_UNDEF;

        if (sraOutput->occCollector != NULL) {
            SRAOCCFree(sraOutput->occCollector);
            sraOutput->occCollector = NULL;
        } else {
            printf("[SAFEGUARD] SRAAlgnmtOutputFreeNone is free-ing occCollector that is already NULL!\n");
        }
        
    } else {
        printf("[SAFEGUARD] SRAAlgnmtOutputFreeNone is free-ing SRAOutput block that is already occupied by output type %u.\n", sraOutput->OutFileFormat);
        printf("Function is invoked for output type:%u\n",SRA_OUTPUT_FORMAT_NONE);
    }
}

void SRAAlgnmtOutputFreeDiscard(SRAOutput * sraOutput) {

    if (sraOutput->OutFileFormat == SRA_OUTPUT_FORMAT_NONE) {
    
        // Setup the status
        sraOutput->OutFileFormat = SRA_OUTPUT_FORMAT_UNDEF;

        if (sraOutput->occCollector != NULL) {
            SRAOCCFree(sraOutput->occCollector);
            sraOutput->occCollector = NULL;
        } else {
            printf("[SAFEGUARD] SRAAlgnmtOutputFreeDiscard is free-ing occCollector that is already NULL!\n");
        }
        
    } else {
        printf("[SAFEGUARD] SRAAlgnmtOutputFreeDiscard is free-ing SRAOutput block that is already occupied by output type %u.\n", sraOutput->OutFileFormat);
        printf("Function is invoked for output type:%u\n",SRA_OUTPUT_FORMAT_NONE);
    }
}

void SRAAlgnmtOutputFreePlain(SRAOutput * sraOutput) {

    if (sraOutput->OutFileFormat == SRA_OUTPUT_FORMAT_PLAIN) {
    
        // Setup the status
        sraOutput->OutFileFormat = SRA_OUTPUT_FORMAT_UNDEF;
        sraOutput->OutFileName[0] = '\0';
        sraOutput->ReadGroup[0] = '\0';
        // Setup the necessary auxiliary data structure
        if (sraOutput->OutFilePtr != NULL) {
            fclose(sraOutput->OutFilePtr);
            sraOutput->OutFilePtr = NULL;
        } else {
            printf("[SAFEGUARD] SRAAlgnmtOutputFreePlain is free-ing OutFilePtr that is already NULL!\n");
        }

        if (sraOutput->occCollector != NULL) {
            SRAOCCFree(sraOutput->occCollector);
            sraOutput->occCollector = NULL;
        } else {
            printf("[SAFEGUARD] SRAAlgnmtOutputFreePlain is free-ing occCollector that is already NULL!\n");
        }
        
    } else {
        printf("[SAFEGUARD] SRAAlgnmtOutputFreePlain is free-ing SRAOutput block that is already occupied by output type %u.\n", sraOutput->OutFileFormat);
        printf("Function is invoked for output type:%u\n",SRA_OUTPUT_FORMAT_PLAIN);
    }
}

void SRAAlgnmtOutputFreeSAMStore(SRAOutput * sraOutput) {

    if (sraOutput->OutFileFormat == SRA_OUTPUT_FORMAT_SAM_STORE) {
    
        // Setup the status
        sraOutput->OutFileFormat = SRA_OUTPUT_FORMAT_UNDEF;

        if (sraOutput->occCollector != NULL) {
            SRAOCCFree(sraOutput->occCollector);
            sraOutput->occCollector = NULL;
        } else {
            printf("[SAFEGUARD] SRAAlgnmtOutputFreeSAMStore is free-ing occCollector that is already NULL!\n");
        }
        
    } else {
        printf("[SAFEGUARD] SRAAlgnmtOutputFreeSAMStore is free-ing SRAOutput block that is already occupied by output type %u.\n", sraOutput->OutFileFormat);
        printf("Function is invoked for output type:%u\n",SRA_OUTPUT_FORMAT_SAM_STORE);
    }
}

void SRAAlgnmtOutputFreeSAM(SRAOutput * sraOutput) {

    if (sraOutput->OutFileFormat == SRA_OUTPUT_FORMAT_SAM) {
    
        // Setup the status
        sraOutput->OutFileFormat = SRA_OUTPUT_FORMAT_UNDEF;
        sraOutput->OutFileName[0] = '\0';
        sraOutput->ReadGroup[0] = '\0';
        // Setup the necessary auxiliary data structure
        if (sraOutput->SAMOutFilePtr != NULL) {
            samclose(sraOutput->SAMOutFilePtr);
            sraOutput->SAMOutFilePtr = NULL;
        } else {
            printf("[SAFEGUARD] SRAAlgnmtOutputFreeSAM is free-ing SAMOutFilePtr that is already NULL!\n");
        }
        
        sraOutput->SAMAuxDataBlockSize = 0;
        if ( sraOutput->SAMAuxDataBlock != NULL ) {
            free(sraOutput->SAMAuxDataBlock);
            sraOutput->SAMAuxDataBlock = NULL;
        } else {
            printf("[SAFEGUARD] SRAAlgnmtOutputFreeSAM is free-ing SAMAuxDataBlock that is already NULL!\n");
        }

        if (sraOutput->occCollector != NULL) {
            SRAOCCFree(sraOutput->occCollector);
            sraOutput->occCollector = NULL;
        } else {
            printf("[SAFEGUARD] SRAAlgnmtOutputFreeSAM is free-ing occCollector that is already NULL!\n");
        }
        
    } else {
        printf("[SAFEGUARD] SRAAlgnmtOutputFreeSAM is free-ing SRAOutput block that is already occupied by output type %u.\n", sraOutput->OutFileFormat);
        printf("Function is invoked for output type:%u\n",SRA_OUTPUT_FORMAT_SAM);
    }
}

void SRAAlgnmtOutputFreeBAM(SRAOutput * sraOutput) {

    if (sraOutput->OutFileFormat == SRA_OUTPUT_FORMAT_BAM) {
    
        // Setup the status
        sraOutput->OutFileFormat = SRA_OUTPUT_FORMAT_UNDEF;
        sraOutput->OutFileName[0] = '\0';
        sraOutput->ReadGroup[0] = '\0';
        // Setup the necessary auxiliary data structure
        if (sraOutput->SAMOutFilePtr != NULL) {
            samclose(sraOutput->SAMOutFilePtr);
            sraOutput->SAMOutFilePtr = NULL;
        } else {
            printf("[SAFEGUARD] SRAAlgnmtOutputFreeBAM is free-ing SAMOutFilePtr that is already NULL!\n");
        }
        
        sraOutput->SAMAuxDataBlockSize = 0;
        if ( sraOutput->SAMAuxDataBlock != NULL ) {
            free(sraOutput->SAMAuxDataBlock);
            sraOutput->SAMAuxDataBlock = NULL;
        } else {
            printf("[SAFEGUARD] SRAAlgnmtOutputFreeBAM is free-ing SAMAuxDataBlock that is already NULL!\n");
        }

        if (sraOutput->occCollector != NULL) {
            SRAOCCFree(sraOutput->occCollector);
            sraOutput->occCollector = NULL;
        } else {
            printf("[SAFEGUARD] SRAAlgnmtOutputFreeBAM is free-ing occCollector that is already NULL!\n");
        }
        
    } else {
        printf("[SAFEGUARD] SRAAlgnmtOutputFreeBAM is free-ing SRAOutput block that is already occupied by output type %u.\n", sraOutput->OutFileFormat);
        printf("Function is invoked for output type:%u\n",SRA_OUTPUT_FORMAT_BAM);
    }
}
//-------------------------------------------------
// End of output Format dictated procedures
//-------------------------------------------------
