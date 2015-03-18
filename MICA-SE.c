//
//    MICA.c
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
//    Modification History
//
//    Date    : 4th June 2013
//    Changes : Correctly distribute the bucket. Currently the numReadPerCore for each core is
//              hardcoded causing inefficient load balancing while the number of reads is not
//              a multiple of the numReadPerCore.
//
//    Date    : 6th June 2013
//    Changes : Add Occurrence Collector to gather output but discard them to avoid IO
//
//    To-do   : 1. Construct the BWT structure outside the pragma offload.
//              2. Don't need inout saCount everytime
//
///////////////////////////////////////////////////////////////////////////////////////////////

#include "MICA-SE.h"
#include "Release.h"

#define MIC_MAIN_NON_INPUT_VALUE    99999999
#define MIC_BATCH_SIZE              20000

#define SRA_READ_SKIP_INVALID      0
#define SRA_READ_REPLACE_INVALID   1

//-----------------------------------
// DEBUG FLAGS
//-----------------------------------
// Define to print output statistics
//#define MIC_DEBUG_PRINT_OUTPUT_COUNTS


// Declaration of Global Variables
// Global Variables should NEVER be referenced outside of
// this file though they can be.
__attribute__((target(mic)))
int InputNumOfMICThreads = 240;
int InputNumOfCPUThreads = 1;

char ReadFileName[MAX_FILENAME_LEN+1] = "";
char OutputFileName[MAX_FILENAME_LEN+1] = "*.out";
char DatabaseName[MAX_FILENAME_LEN+1] = "";
char InputSaValueFileName[MAX_FILENAME_LEN+1] = ".sa8";
char InputCPUSaValueFileName[MAX_FILENAME_LEN+1] = ".sa";

int InputSRAInvalidReadHandling = SRA_READ_SKIP_INVALID;
int InputQueryStrand = QUERY_BOTH_STRAND;
int InputAlignmentModel = SRA_MODEL_8G;
int InputMaxNumOfAlignment = -1;

int InputMaxError = 0;
int InputMaxNBMismatch = MIC_MAIN_NON_INPUT_VALUE;
int InputErrorType = MIC_MAIN_NON_INPUT_VALUE;
int InputOutputFileName = MIC_MAIN_NON_INPUT_VALUE;
int OutputFormat = MIC_MAIN_NON_INPUT_VALUE;

unsigned int InputMemoryCPUMaxUsageMBytes = -1;
unsigned int InputMemoryMICMaxUsageMBytes = 3600;
    
void PrintHelp();
dictionary *ParseInput(int argc, char** argv);
void ParseIniFile(char *iniFileName);
void SRAFillCharMap(unsigned char * charMap);
void PrintAlignmentSettings(SRAArguments * sArgs);
void SRAGetFileSuffix(int threadId, char * suffixStr);

// GetBatchSizing function takes the number of reads and
// calculates required number of cores; and the number of reads assigned to each core.
// The output number (batchSize * numThreads) can be greater than (totalNumRead).
void GetBatchSizing(unsigned int totalNumRead,
                    unsigned int maxNumThreads, unsigned int minBatchSize,
                    unsigned int * batchSize, unsigned int * numThreads) {

    // If totalNumRead is empty, the output should all be zero
    // i.e. 0 number of CPU and 0 read in each core
    if (totalNumRead==0) {
        (*batchSize) = 0;
        (*numThreads) = 0;
        return;
    }

    // If totalNumRead is smaller than the number of threads available,
    // each core gets at most 1 read.
    if (totalNumRead<=maxNumThreads) {
        (*batchSize) = 1;
        (*numThreads) = totalNumRead;
        return;
    }

    // Otherwise, distribute the big bucket of reads to the cores
    unsigned int tmp = (totalNumRead + (maxNumThreads - 1)) / maxNumThreads;
    (*batchSize) = tmp;
    (*numThreads) = maxNumThreads;
    return;
}

int main(int argc, char** argv) {

    printf("\n%s v%d.%d.%d (%s)\n",PROJECT_NAME,PROJECT_MAJOR,PROJECT_MINOR,PROJECT_REV,PROJECT_SPECIAL);
    
    ////////////////////////////////////////////////////////////
    // Declaration of Variables
    ////////////////////////////////////////////////////////////
    // Program input
    dictionary *programInput;
    double startTime, lastEventTime, timestamp;
    startTime = setStartTime();
    lastEventTime = 0;
    unsigned long long i,j,k;
    unsigned int u,v,w;

    unsigned char charMap[256];
    unsigned char complementMap[256];
    char dummyQuality[SRA_MAX_READ_LENGTH];
    memset(dummyQuality,0,sizeof(int)*SRA_MAX_READ_LENGTH);
    
    //Initial Environmental Parameters for Intel MIC
    MICSetParameters();

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
    printf("Loading index %s ... ",DatabaseName); fflush(stdout);
    Idx2BWT * idx2BWT = BWTLoad2BWT(DatabaseName,InputSaValueFileName);
    
    char LookupTableFileName[MAX_FILENAME_LEN+1];
    char RevLookupTableFileName[MAX_FILENAME_LEN+1];
    strcpy(LookupTableFileName,DatabaseName);
    strcpy(RevLookupTableFileName,DatabaseName);
    strcat(LookupTableFileName,".lkt");
    strcat(RevLookupTableFileName,".rev.lkt");
    LT * lookup = LTLoad(LookupTableFileName);
    LT * rev_lookup = LTLoad(RevLookupTableFileName);
    
    BWT * bwt = idx2BWT->bwt;
    BWT * rev_bwt = idx2BWT->rev_bwt;
    HSP * hsp = idx2BWT->hsp;
    printf("DONE\n\n");

    // Timestamp'ed -----------
    timestamp = getElapsedTime(startTime);
    printf("(Elapsed time (Index Loading) : %9.4f seconds)\n\n", timestamp - lastEventTime);
    lastEventTime = timestamp;
    // ------------------------
    
    
    
    
    
    
    






    ////////////////////////////////////////////////////////////
    // CPU Allocation of Memory Frame
    ////////////////////////////////////////////////////////////
    unsigned int numReadPerCore = MIC_BATCH_SIZE;
    unsigned int memsize_readNameBufferFrame    = MAX_SEQ_NAME_LENGTH*numReadPerCore*InputNumOfMICThreads;
    unsigned int memsize_readBodyBufferFrame    = SRA_MAX_READ_LENGTH*numReadPerCore*InputNumOfMICThreads;
    unsigned int memsize_outputBufferFrame      = MIC_SRA_OUTPUT_SIZE_PER_THREAD*InputNumOfMICThreads;
    unsigned int memsize_outputBufferMeta       = MIC_SRA_OUTPUT_META_PER_THREAD*InputNumOfMICThreads;
    unsigned int memsize_bufferPerRead          = numReadPerCore*InputNumOfMICThreads;

    char * readNameBufferFrame                  = (char*)               MEMManMalloc(sizeof(char)*memsize_readNameBufferFrame,MEMORY_TYPE_CPU);
    unsigned char * readBodyBufferFrame         = (unsigned char*)      MEMManMalloc(sizeof(unsigned char)*memsize_readBodyBufferFrame,MEMORY_TYPE_SHARED);
    char * readQualityBufferFrame               = (char*)      MEMManMalloc(sizeof(char)*memsize_readBodyBufferFrame,MEMORY_TYPE_SHARED);;
    uint16_t * readLengthBufferFrame            = (uint16_t*)           MEMManMalloc(sizeof(uint16_t)*memsize_bufferPerRead,MEMORY_TYPE_SHARED);
    int * uncertaintyNumberBufferFrame          = (int*)                MEMManMalloc(sizeof(int)*memsize_bufferPerRead,MEMORY_TYPE_SHARED);

    unsigned int * threadAlgmt                  = (unsigned int*)       MEMManMalloc(sizeof(unsigned int)*InputNumOfMICThreads,MEMORY_TYPE_SHARED);
    uint16_t * occCount                         = (uint16_t*)           MEMManMalloc(sizeof(uint16_t)*memsize_bufferPerRead,MEMORY_TYPE_SHARED);
    unsigned int * outputBufferFrame            = (unsigned int*)       MEMManMalloc(sizeof(unsigned int)*memsize_outputBufferFrame,MEMORY_TYPE_SHARED);
    uint8_t * outputBufferStatus                = (uint8_t*)            MEMManMalloc(sizeof(uint8_t)*memsize_bufferPerRead,MEMORY_TYPE_SHARED);
    MICSRAOccMetadata * outputBufferMeta        = (MICSRAOccMetadata*)  MEMManMalloc(sizeof(MICSRAOccMetadata)*memsize_outputBufferMeta,MEMORY_TYPE_SHARED);

    CPTSRAModel * cpModels                      = (CPTSRAModel*)        MEMManMalloc(sizeof(CPTSRAModel) * SRA_MAX_READ_LENGTH,MEMORY_TYPE_SHARED);

    printf("[INFO] Total Shared Working Memory allocated = %9.2f Mbytes\n",MEMManGetUsage(MEMORY_TYPE_SHARED)/1024.0/1024.0);
    if (InputMemoryMICMaxUsageMBytes!=-1 && MEMManGetUsage(MEMORY_TYPE_SHARED)>(unsigned long long) InputMemoryMICMaxUsageMBytes*1024*1024) {
        printf("[ERROR] Exceeded defined MIC limitation of %u MBytes\n",InputMemoryMICMaxUsageMBytes);
    }
    
    printf("[INFO] Total CPU Working Memory allocated = %9.2f Mbytes\n",MEMManGetUsage(MEMORY_TYPE_CPU)/1024.0/1024.0);
    if (InputMemoryCPUMaxUsageMBytes!=-1 && MEMManGetUsage(MEMORY_TYPE_CPU)>(unsigned long long) InputMemoryCPUMaxUsageMBytes*1024*1024) {
        printf("[ERROR] Exceeded defined CPU limitation of %u MBytes\n",InputMemoryMICMaxUsageMBytes);
    }
    
    printf("\n");

    // Timestamp'ed -----------
    timestamp = getElapsedTime(startTime);
    printf("(Elapsed time (CPU Memory Allocation) : %9.4f seconds)\n\n", timestamp - lastEventTime);
    lastEventTime = timestamp;
    // ------------------------











    ////////////////////////////////////////////////////////////
    // MIC Memory Allocation
    ////////////////////////////////////////////////////////////
    
    // Forward BWT
    unsigned int bwtCodeLengthInFile = BWTFileSizeInWord(bwt->textLength);
    unsigned int bwtTextLength = bwt->textLength;
    unsigned int bwtOccSizeInWord = bwt->occSizeInWord;
    unsigned int bwtInverseSa0 = bwt->inverseSa0;

    unsigned int * bwtBwtCode = bwt->bwtCode;
    unsigned int * bwtCumulativeFreq = bwt->cumulativeFreq;
    unsigned int * bwtOccValue = bwt->occValue;
    unsigned int * bwtOccValueMajor = bwt->occValueMajor;
    unsigned int * micOccValue = bwt->micOccValue;
    
    unsigned int micBwtSaValueSizeInWord = bwt->saValueSizeInWord;
    unsigned int * micBwtSaValue = bwt->saValue;
    
    #pragma offload target(mic) \
        nocopy(bwt:length(1)                                                alloc_if(1) free_if(0) align(2*1024*1024)) \
        nocopy(bwtBwtCode:length(bwtCodeLengthInFile)                       alloc_if(1) free_if(0) align(2*1024*1024)) \
        nocopy(bwtCumulativeFreq:length(ALPHABET_SIZE+1)                    alloc_if(1) free_if(0) align(2*1024*1024)) \
        nocopy(bwtOccValue:length(bwt->occSizeInWord)                       alloc_if(1) free_if(0) align(2*1024*1024)) \
        nocopy(bwtOccValueMajor:length(bwt->occMajorSizeInWord)             alloc_if(1) free_if(0) align(2*1024*1024)) \
        nocopy(micBwtSaValue:length(micBwtSaValueSizeInWord)                alloc_if(1) free_if(0) align(2*1024*1024)) \
        nocopy(micOccValue:length(bwt->micOccSizeInWord)                    alloc_if(1) free_if(0) align(2*1024*1024))
        
    {}
    
    // Reverse BWT
    unsigned int rbwtCodeLengthInFile = BWTFileSizeInWord(rev_bwt->textLength);
    unsigned int rbwtTextLength = rev_bwt->textLength;
    unsigned int rbwtOccSizeInWord = rev_bwt->occSizeInWord;
    unsigned int rbwtInverseSa0 = rev_bwt->inverseSa0;

    unsigned int * rbwtBwtCode = rev_bwt->bwtCode;
    unsigned int * rbwtCumulativeFreq = rev_bwt->cumulativeFreq;
    unsigned int * rbwtOccValue = rev_bwt->occValue;
    unsigned int * rbwtOccValueMajor = rev_bwt->occValueMajor;
    unsigned int * rmicOccValue = rev_bwt->micOccValue;

    #pragma offload target(mic) \
        nocopy(rev_bwt:length(1)                                            alloc_if(1) free_if(0) align(2*1024*1024)) \
        nocopy(rbwtBwtCode:length(rbwtCodeLengthInFile)                     alloc_if(1) free_if(0) align(2*1024*1024)) \
        nocopy(rbwtCumulativeFreq:length(ALPHABET_SIZE+1)                   alloc_if(1) free_if(0) align(2*1024*1024)) \
        nocopy(rbwtOccValue:length(rev_bwt->occSizeInWord)                  alloc_if(1) free_if(0) align(2*1024*1024)) \
        nocopy(rbwtOccValueMajor:length(rev_bwt->occMajorSizeInWord)        alloc_if(1) free_if(0) align(2*1024*1024)) \
        nocopy(rmicOccValue:length(rev_bwt->micOccSizeInWord)               alloc_if(1) free_if(0) align(2*1024*1024))
        
    {}
    
    // Input and Output Data
    #pragma offload target(mic) \
        nocopy(readBodyBufferFrame:length(memsize_readBodyBufferFrame)      alloc_if(1) free_if(0) align(2*1024*1024)) \
        nocopy(readLengthBufferFrame:length(memsize_bufferPerRead)          alloc_if(1) free_if(0) align(2*1024*1024)) \
        \
        nocopy(cpModels:length(SRA_MAX_READ_LENGTH)                         alloc_if(1) free_if(0) align(2*1024*1024)) \
        nocopy(occCount:length(memsize_bufferPerRead)                       alloc_if(1) free_if(0) align(2*1024*1024)) \
        nocopy(threadAlgmt:length(InputNumOfMICThreads)                     alloc_if(1) free_if(0) align(2*1024*1024)) \
        \
        nocopy(outputBufferFrame:length(memsize_outputBufferFrame)          alloc_if(1) free_if(0) align(2*1024*1024)) \
        nocopy(outputBufferMeta:length(memsize_outputBufferMeta)            alloc_if(1) free_if(0) align(2*1024*1024)) \
        nocopy(outputBufferStatus:length(memsize_bufferPerRead)             alloc_if(1) free_if(0) align(2*1024*1024))
    {}

    // Timestamp'ed -----------
    timestamp = getElapsedTime(startTime);
    printf("(Elapsed time (MIC Memory Allocation) : %9.4f seconds)\n\n", timestamp - lastEventTime);
    lastEventTime = timestamp;
    // ------------------------

    
    
    
    
    
    
    
    
    
    ////////////////////////////////////////////////////////////
    // MIC Memory Copy
    ////////////////////////////////////////////////////////////

    // Forward BWT
    #pragma offload target(mic) \
        in(bwt:length(1)                                                    alloc_if(0) free_if(0)) \
        in(bwtBwtCode:length(bwtCodeLengthInFile)                           alloc_if(0) free_if(0)) \
        in(bwtCumulativeFreq:length(ALPHABET_SIZE+1)                        alloc_if(0) free_if(0)) \
        in(bwtOccValue:length(bwt->occSizeInWord)                           alloc_if(0) free_if(0)) \
        in(bwtOccValueMajor:length(bwt->occMajorSizeInWord)                 alloc_if(0) free_if(0)) \
        in(micBwtSaValue:length(micBwtSaValueSizeInWord)                    alloc_if(0) free_if(0)) \
        in(micOccValue:length(bwt->micOccSizeInWord)                        alloc_if(0) free_if(0))
    {
        bwt->cumulativeFreq = bwtCumulativeFreq;
        bwt->bwtCode = bwtBwtCode;
        bwt->occValue = bwtOccValue;
        bwt->occValueMajor = bwtOccValueMajor;
        bwt->saValue = micBwtSaValue;
        bwt->micOccValue = micOccValue;
    }

    // Reverse BWT
    #pragma offload target(mic) \
        in(rev_bwt:length(1)                                                alloc_if(0) free_if(0)) \
        in(rbwtBwtCode:length(rbwtCodeLengthInFile)                         alloc_if(0) free_if(0)) \
        in(rbwtCumulativeFreq:length(ALPHABET_SIZE+1)                       alloc_if(0) free_if(0)) \
        in(rbwtOccValue:length(rev_bwt->occSizeInWord)                      alloc_if(0) free_if(0)) \
        in(rbwtOccValueMajor:length(rev_bwt->occMajorSizeInWord)            alloc_if(0) free_if(0)) \
        in(rmicOccValue:length(rev_bwt->micOccSizeInWord)                   alloc_if(0) free_if(0))
    {
        rev_bwt->cumulativeFreq = rbwtCumulativeFreq;
        rev_bwt->bwtCode = rbwtBwtCode;
        rev_bwt->occValue = rbwtOccValue;
        rev_bwt->occValueMajor = rbwtOccValueMajor;
        rev_bwt->micOccValue = rmicOccValue;
    }

    // Timestamp'ed -----------
    timestamp = getElapsedTime(startTime);
    printf("(Elapsed time (MIC Memory Copy) : %9.4f seconds)\n\n", timestamp - lastEventTime);
    lastEventTime = timestamp;
    // ------------------------

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    ////////////////////////////////////////////////////////////
    // CPU Index Re-Loading
    ////////////////////////////////////////////////////////////
    
    char saFilename[MAX_INDEX_FILENAME_LENGTH]; 
    strcpy(saFilename,DatabaseName);
    strcat(saFilename,InputCPUSaValueFileName);
    printf("Loading Full Suffix-Array for CPU %s ... ",saFilename); fflush(stdout);
    BWTSAFree(idx2BWT->mmPool,idx2BWT->bwt);
    BWTSALoad(idx2BWT->mmPool,idx2BWT->bwt,saFilename,NULL);
    printf("DONE\n\n");

    // Timestamp'ed -----------
    timestamp = getElapsedTime(startTime);
    printf("(Elapsed time (Index Loading) : %9.4f seconds)\n\n", timestamp - lastEventTime);
    lastEventTime = timestamp;
    // ------------------------
    
    
    
    
    
    
    
    
    
    
    
    


    
    SRAArguments * sraArgs = SRAARGConstruct();
    
    // Setting up SRA Indexes
    SRAIndexPopulate(sraArgs->AlgnmtIndex,idx2BWT->bwt,idx2BWT->rev_bwt,idx2BWT->hsp,NULL,lookup,rev_lookup);
    
    // Setting up SRA Settings
    sraArgs->AlgnmtSetting->MaxError = InputMaxError;
    sraArgs->AlgnmtSetting->MaxResult= InputMaxNumOfAlignment;
    
    // Setting up output file
    bam_header_t samOutputHeader;
    SAMOutputHeaderConstruct(&samOutputHeader,hsp);
    sraArgs->AlgnmtOutput->occCollector  = SRAOCCCreate(SRAOCC_FLOOD_TYPE_EXPAND);
    sraArgs->AlgnmtOutput->OutFileFormat = SRA_OUTPUT_FORMAT_SAM;
    sraArgs->AlgnmtSetting->OutputFormat = SRA_OUTPUT_FORMAT_SAM;
    sraArgs->AlgnmtOutput->SAMOutFilePtr = samopen("out","wh",&samOutputHeader);
    
    SRAArguments * sraArgs_neg = SRAARGMakeSlave(sraArgs);
    
    SRAModelSet * sraModelSet = SRAModelSetConstruct(
                                    sraArgs->AlgnmtSetting,
                                    sraArgs->AlgnmtIndex,
                                    SRA_MODEL_16G,
                                    SRA_MIN_READ_LENGTH,
                                    SRA_MAX_READ_LENGTH);
    
    
    


    //Printing Settings
    PrintAlignmentSettings(sraArgs);
    
    
    
    
    
    

    UTBFRBuffer * readQueryBuffer = UTBFRLoad(ReadFileName);
    HSPFillCharMap(charMap);
    if (InputSRAInvalidReadHandling == SRA_READ_REPLACE_INVALID) SRAFillCharMap(charMap);
    HSPFillComplementMap(complementMap);

    SRAOCCCollector * occCollector = sraArgs->AlgnmtOutput->occCollector;
    
    double processingTime = 0;
    double alignmentTime = 0;
    double outputTime = 0;
    double readLoadTime = 0;
    double modelTime = 0;

    unsigned int queryInBatch = 1;
    unsigned long long consSaCount = 0;
    unsigned long long cpuConsSaCount = 0;
    unsigned long long queryHandled = 0;

    #ifdef MIC_DEBUG_PRINT_OUTPUT_COUNTS
    unsigned int numOfClosed = 0;
    unsigned int numOfSaReq = 0;
    unsigned int numOfSaContent = 0;
    unsigned int numOfOccContent = 0;
    #endif
    
    while (queryInBatch!=0) {

        // Read in the query for the buckets for all threads
        // The current stamp is being retrieved
        queryInBatch = SRAQueryAndNameGetBatchFromFASTAQLength(readQueryBuffer,
                                                        charMap,
                                                        numReadPerCore * InputNumOfMICThreads,
                                                        MAX_SEQ_NAME_LENGTH,
                                                        readNameBufferFrame,
                                                        SRA_MAX_READ_LENGTH,
                                                        readBodyBufferFrame,
                                                        readQualityBufferFrame,
                                                        readLengthBufferFrame,
                                                        uncertaintyNumberBufferFrame);

        if (queryInBatch == 0) break;
        
        // Timestamp'ed -----------
        timestamp = getElapsedTime(startTime);
        //printf("(Elapsed time (Read Loading) : %9.4f seconds)\n\n", timestamp - lastEventTime);
        readLoadTime += timestamp - lastEventTime;
        processingTime += timestamp - lastEventTime;
        lastEventTime = timestamp;
        // ------------------------
        
        // Build up the SRA Model in case there are new read lengths
        for (i=0;i<queryInBatch;i++) {
            unsigned int readLength = readLengthBufferFrame[i];
            
            SRAModel * sraModels = SRAModelSetGetModel(sraModelSet,readLength,QUERY_POS_STRAND);
    
            if (sraModels==NULL) {
                SRAModelConstruct(sraModelSet,readLength);
                sraModels = SRAModelSetGetModel(sraModelSet,readLength,QUERY_POS_STRAND);
                MICSRA2BWTModelPopulate(&(cpModels[readLength]),sraModels);
                
                //Print the SRAModel for Debug Purposes
                // DebugPrintModel(sraModels,stdout);
                // MICSRADebugPrintModel(&cpModels[readLength],stdout);
            }
        }
        
        // Timestamp'ed -----------
        timestamp = getElapsedTime(startTime);
        //printf("(Elapsed time (Read Loading) : %9.4f seconds)\n\n", timestamp - lastEventTime);
        modelTime += timestamp - lastEventTime;
        processingTime += timestamp - lastEventTime;
        lastEventTime = timestamp;
        // ------------------------

        /*#pragma offload target(mic) \
            in(readBodyBufferFrame:length(memsize_readBodyBufferFrame)          alloc_if(0) free_if(0)) \
            in(readLengthBufferFrame:length(memsize_bufferPerRead)      alloc_if(0) free_if(0))
        {}

        // Timestamp'ed -----------
        timestamp = getElapsedTime(startTime);
        printf("(Elapsed time (MIC Read Copy) : %9.4f seconds)\n\n", timestamp - lastEventTime);
        lastEventTime = timestamp;
        // ------------------------*/





        unsigned int batchSize, numThreads;
        GetBatchSizing(queryInBatch,InputNumOfMICThreads,1,&batchSize,&numThreads);
        memset(threadAlgmt,0,sizeof(unsigned int)*InputNumOfMICThreads);
        memset(occCount,0,sizeof(uint16_t)*memsize_bufferPerRead);
        memset(outputBufferStatus,MIC_OUTPUT_STATUS_OPEN,sizeof(unsigned char)*memsize_bufferPerRead);




        printf("MIC Core initiated (%u) ... ",numThreads); fflush(stdout);
        ////////////////////////////////////////////////////////////
        // Alignment Body
        ////////////////////////////////////////////////////////////

        #pragma offload target(mic) \
            nocopy(bwt:length(1)                                                alloc_if(0) free_if(0)) \
            nocopy(rev_bwt:length(1)                                            alloc_if(0) free_if(0)) \
            \
            in(readBodyBufferFrame:length(memsize_readBodyBufferFrame)          alloc_if(0) free_if(0)) \
            in(readLengthBufferFrame:length(memsize_bufferPerRead)              alloc_if(0) free_if(0)) \
            in(cpModels:length(SRA_MAX_READ_LENGTH)                             alloc_if(0) free_if(0)) \
            \
            inout(threadAlgmt:length(InputNumOfMICThreads)                      alloc_if(0) free_if(0)) \
            inout(occCount:length(memsize_bufferPerRead)                        alloc_if(0) free_if(0)) \
            \
            out(outputBufferFrame:length(memsize_outputBufferFrame)             alloc_if(0) free_if(0)) \
            out(outputBufferMeta:length(memsize_outputBufferMeta)               alloc_if(0) free_if(0)) \
            inout(outputBufferStatus:length(memsize_bufferPerRead)              alloc_if(0) free_if(0))

        #pragma omp parallel for private(i,j) num_threads(InputNumOfMICThreads)
        for (i=0;i<numThreads;i++) {

            unsigned int offset = i * batchSize;
            int charIdx;
            unsigned char c;
            uint16_t metaCount;
            
            // Shared in core
            unsigned int * outputPtr     = &outputBufferFrame[i*MIC_SRA_OUTPUT_SIZE_PER_THREAD];
            MICSRAOccMetadata * metaPtr  = &outputBufferMeta[i*MIC_SRA_OUTPUT_META_PER_THREAD];
            unsigned int outputVacancy   = MIC_SRA_OUTPUT_SIZE_PER_THREAD;
            unsigned int metaVacancy     = MIC_SRA_OUTPUT_META_PER_THREAD;

            // Multiple cell per read
            unsigned char * patternPtr   = &readBodyBufferFrame[offset*SRA_MAX_READ_LENGTH];
            
            MICSRAArguments * micArgs = MICSRAARGConstruct();
            // ATTENTION. BWT is hacked into MICSRAArgs before
            // proper implementation of MICSRAMdl.
            micArgs->outputType      = SRA_REPORT_ALL;
            micArgs->bwt             = bwt;
            micArgs->rev_bwt         = rev_bwt;
            
            #pragma ivdep
            for (j=0;j<batchSize;j++) {

                // One per read
                uint16_t patternLen          = readLengthBufferFrame[offset + j];
                uint8_t * outputStatus       = &outputBufferStatus[offset + j];
                uint16_t * _occCount         = &occCount[offset + j];
                metaCount                    = 0;
                
                ////////////////////////////////////////////
                // READ
                ////////////////////////////////////////////
                if (outputVacancy<MIC_SRA_OUTPUT_MAX_ALIGNMENT) {
                    //Skip if the global output buffer is flooded
                    (*outputStatus) = MIC_OUTPUT_STATUS_UNHANDLE;
                } else if (metaVacancy<MIC_SRA_OUTPUT_MAX_META) {
                    //Skip if the global output buffer is flooded
                    (*outputStatus) = MIC_OUTPUT_STATUS_UNHANDLE;
                } else if (patternLen==0) {
                    //Skip if the read is an invalid read
                    (*outputStatus) = MIC_OUTPUT_STATUS_SKIPPED;
                } else {
                    
                    micArgs->seedLength     = patternLen;
                    micArgs->readLength     = patternLen;
                    micArgs->readCode       = patternPtr;
                    micArgs->outputStatus   = outputStatus;
                    micArgs->outputBlock    = outputPtr;
                    micArgs->occCount       = _occCount;
                    micArgs->metaBlock      = metaPtr;
                    micArgs->metaCount      = &metaCount;
                    
                    // Calling the model aligner
                    threadAlgmt[i] += MICProcessReadDoubleStrand(micArgs,&(cpModels[patternLen]),&(cpModels[patternLen]));
                    
                    // This is equivelance to if (*outputStatus == MIC_OUTPUT_STATUS_OPEN) then (*outputStatus) = MIC_OUTPUT_STATUS_COMPLETE
                    (*outputStatus) |= MIC_OUTPUT_STATUS_COMPLETE * (*outputStatus == MIC_OUTPUT_STATUS_OPEN);
                }
                
                patternPtr += SRA_MAX_READ_LENGTH;
                outputPtr += (*_occCount);
                outputVacancy -= (*_occCount);
                metaPtr += (metaCount);
                metaVacancy -= (metaCount);
            }
            
            MICSRAARGFree(micArgs);
        }
        printf("DONE\n"); 

        // Timestamp'ed -----------
        timestamp = getElapsedTime(startTime);
        alignmentTime += timestamp - lastEventTime;
        processingTime += timestamp - lastEventTime;
        lastEventTime = timestamp;
        // ------------------------

        for (i=0;i<numThreads;i++) {
            consSaCount += threadAlgmt[i];
        }

        ////////////////////////////////////////////////////////////
        // Produce Statistics of Output
        ////////////////////////////////////////////////////////////
        #ifdef MIC_DEBUG_PRINT_OUTPUT_COUNTS
        for (i=0;i<queryInBatch;i++) {
            if ( outputBufferStatus[i] == MIC_OUTPUT_STATUS_CLOSE) {
                numOfClosed++;
                printf("Closed Item size = %u\n",occCount[i]);
            }  
            unsigned int * outputPtr = &outputBufferFrame[i*MIC_SRA_OUTPUT_SIZE_PER_THREAD];
            numOfOccContent+=occCount[i];
        }
        printf("[STATISTICS] Number of Occurrences found = %u\n",numOfOccContent);
        printf("[STATISTICS] Number of Closed Output Collector = %u\n",numOfClosed);
        #endif
        
        ////////////////////////////////////////////////////////////
        // Handling Output
        // This following pick results from MIC and output with CPU
        ////////////////////////////////////////////////////////////
        unsigned char negPattern[SRA_MAX_READ_LENGTH];
        for (i=0;i<numThreads;i++) {
        
            unsigned int offset = i * batchSize;
            unsigned int * outputPtr     = &outputBufferFrame[i*MIC_SRA_OUTPUT_SIZE_PER_THREAD];
            MICSRAOccMetadata * metaPtr  = &outputBufferMeta[i*MIC_SRA_OUTPUT_META_PER_THREAD];
            unsigned int metaCount       = metaPtr->numPayload;
            
            for (j=0;j<batchSize;j++) {     

                unsigned int bufferIdx = offset + j;

                unsigned int readIdx = queryHandled + bufferIdx;
            
                unsigned char * patternPtr = &readBodyBufferFrame[bufferIdx*SRA_MAX_READ_LENGTH];
                char * patternName         = &readNameBufferFrame[bufferIdx*MAX_SEQ_NAME_LENGTH];
                uint16_t patternLen        = readLengthBufferFrame[bufferIdx];
                SRAFlipPattern(charMap,patternPtr,patternLen,negPattern);
                
                SRAQueryInfoPopulate(sraArgs->QueryInfo,readIdx,patternName,patternLen,patternPtr,patternPtr,QUERY_POS_STRAND,dummyQuality);
                SRAQueryInfoPopulate(sraArgs_neg->QueryInfo,readIdx,patternName,patternLen,negPattern,patternPtr,QUERY_NEG_STRAND,dummyQuality);
                
                SRAModel * sraModels = SRAModelSetGetModel(sraModelSet,patternLen,QUERY_POS_STRAND);
                
                if ( outputBufferStatus[bufferIdx] == MIC_OUTPUT_STATUS_CLOSE) {
                    printf("[INFO] Unaligned closed read - MIC unaligned - %s!\n",patternName);
                    // In this case the output buffer is flooded as the number of output
                    // generated is larger than the MIC_SRA_OUTPUT_MAX_ALIGNMENT; or the number of header
                    // generated is larger than MIC_SRA_OUTPUT_MAX_META. MIC was not
                    // able to handle all the alignment result hence declared the read 
                    // unaligned on MIC.
                    cpuConsSaCount += ProcessReadDoubleStrand(sraArgs,sraArgs_neg,sraModels,sraModels);
                } else if ( outputBufferStatus[bufferIdx] == MIC_OUTPUT_STATUS_UNHANDLE) {
                    printf("[INFO] Unhandled read - MIC unaligned - %s!\n",patternName);
                    // In this case the entire core shared output buffer is flooded as the number of output
                    // generated is larger than the MIC_SRA_OUTPUT_SIZE_PER_THREAD. MIC was not
                    // able to handle all the reads in the input read block hence declared all consecutive read 
                    // unaligned on MIC.
                    cpuConsSaCount += ProcessReadDoubleStrand(sraArgs,sraArgs_neg,sraModels,sraModels);
                } else if ( outputBufferStatus[bufferIdx] == MIC_OUTPUT_STATUS_SKIPPED) {
                    if ( bufferIdx < queryInBatch ) {
                        printf("[INFO] Skipped read - MIC unaligned - %s!\n",patternName);
                    }
                    // In this case the MIC actively skipped the alignment of this read. Reason being the 
                    // read contains invalid character or unknown reason. CPU should also skip it.
                    // These sample code gives output when these reads are found from the output block and attempt to process it.
                    // cpuConsSaCount += ProcessReadDoubleStrand(sraArgs,sraArgs_neg,sraModels,sraModels);
                } else if ( outputBufferStatus[bufferIdx] == MIC_OUTPUT_STATUS_OPEN) {
                    printf("[ERROR] Skipped read - MIC still open - %s!\n",patternName);
                    // In this case the MIC actively skipped the alignment of this read. For unknown reason.
                    cpuConsSaCount += ProcessReadDoubleStrand(sraArgs,sraArgs_neg,sraModels,sraModels);
                } else if (occCount[bufferIdx]==0) {
                    OCCReportNoAlignment(sraArgs);
                } else {
                    if (occCount[bufferIdx]>0) {
                        
                        for (k=0;k<occCount[bufferIdx];k++) {
                            if (metaCount==0) {
                                metaPtr++;
                                metaCount = metaPtr->numPayload;
                            }
                            
                            OCCAddTextPositionToCache(sraArgs,SRAOCC_TYPE_AWAIT_TRANSLATE,
                                                    outputPtr[k], (metaPtr->strand + 1),
                                                    patternLen, metaPtr->numOfErr, 0,
                                                    metaPtr->errors);
                            metaCount--;
                        }
                    }
                }
                outputPtr += occCount[bufferIdx];;
                OCCFlushCache(sraArgs);
            }
        }
        
        queryHandled+=queryInBatch;

        // Timestamp'ed -----------
        timestamp = getElapsedTime(startTime);
        //printf("(Elapsed time (Output) : %9.4f seconds)\n\n", timestamp - lastEventTime);
        outputTime += timestamp - lastEventTime;
        processingTime += timestamp - lastEventTime;
        lastEventTime = timestamp;
        // ------------------------

        printf("Number of Query Handled = %llu\n",queryHandled);
        printf("Consolidated SA Count (MIC-Accm) = %llu\n",consSaCount);
        printf("Consolidated SA Count (CPU-Accm) = %llu\n",cpuConsSaCount);
        printf("Read Loading Time       = %9.4f seconds\n",readLoadTime);
        printf("Model Building Time     = %9.4f seconds\n",modelTime);
        printf("Alignment Time          = %9.4f seconds\n",alignmentTime);
        printf("Output Time             = %9.4f seconds\n",outputTime);
        printf("Processing (Total) Time = %9.4f seconds\n",processingTime);
        printf("\n");
        printf("\n");

    }


    SRAModelSetFree(sraModelSet);
    UTBFRFree(readQueryBuffer);
    SRAARGFree(sraArgs);
    SRAARGSlaveFree(sraArgs_neg);
    SAMOutputHeaderDestruct(&samOutputHeader);

    
    ////////////////////////////////////////////////////////////
    // Free MIC Memory Allocated
    ////////////////////////////////////////////////////////////
    #pragma offload target(mic) \
        nocopy(bwt:length(1)                                                alloc_if(0) free_if(1)) \
        nocopy(bwtBwtCode:length(bwtCodeLengthInFile)                       alloc_if(0) free_if(1)) \
        nocopy(bwtCumulativeFreq:length(ALPHABET_SIZE+1)                    alloc_if(0) free_if(1)) \
        nocopy(bwtOccValue:length(bwt->occSizeInWord)                       alloc_if(0) free_if(1)) \
        nocopy(bwtOccValueMajor:length(bwt->occMajorSizeInWord)             alloc_if(0) free_if(1)) \
        nocopy(micBwtSaValue:length(micBwtSaValueSizeInWord)                alloc_if(0) free_if(1)) \
        nocopy(micOccValue:length(bwt->micOccSizeInWord)                    alloc_if(0) free_if(1))
    {}
    
    #pragma offload target(mic) \
        nocopy(rev_bwt:length(1)                                            alloc_if(0) free_if(1)) \
        nocopy(rbwtBwtCode:length(rbwtCodeLengthInFile)                     alloc_if(0) free_if(1)) \
        nocopy(rbwtCumulativeFreq:length(ALPHABET_SIZE+1)                   alloc_if(0) free_if(1)) \
        nocopy(rbwtOccValue:length(rev_bwt->occSizeInWord)                  alloc_if(0) free_if(1)) \
        nocopy(rbwtOccValueMajor:length(rev_bwt->occMajorSizeInWord)        alloc_if(0) free_if(1)) \
        nocopy(rmicOccValue:length(rev_bwt->micOccSizeInWord)               alloc_if(0) free_if(1))
    {}
    
    // Input and Output Data
    #pragma offload target(mic) \
        nocopy(readBodyBufferFrame:length(memsize_readBodyBufferFrame)      alloc_if(0) free_if(1)) \
        nocopy(readLengthBufferFrame:length(memsize_bufferPerRead)          alloc_if(0) free_if(1)) \
        \
        nocopy(cpModels:length(SRA_MAX_READ_LENGTH)                         alloc_if(0) free_if(1)) \
        nocopy(occCount:length(memsize_bufferPerRead)                       alloc_if(0) free_if(1)) \
        nocopy(threadAlgmt:length(InputNumOfMICThreads)                     alloc_if(0) free_if(1)) \
        \
        nocopy(outputBufferFrame:length(memsize_outputBufferFrame)          alloc_if(0) free_if(1)) \
        nocopy(outputBufferMeta:length(memsize_outputBufferMeta)            alloc_if(0) free_if(1)) \
        nocopy(outputBufferStatus:length(memsize_bufferPerRead)             alloc_if(0) free_if(1))
    {}
    
    ////////////////////////////////////////////////////////////
    // Free Memory Allocated
    ////////////////////////////////////////////////////////////
    printf("\nFree index ... ");
    fflush(stdout);
    LTFree(lookup);
    LTFree(rev_lookup);
    BWTFree2BWT(idx2BWT);
    printf("DONE\n");
    free(readNameBufferFrame);
    free(readBodyBufferFrame);
    free(readLengthBufferFrame);
    free(threadAlgmt);
    free(occCount);
    free(outputBufferFrame);
    free(outputBufferStatus);
    free(outputBufferMeta);
    free(cpModels);
    
    return 0;
}


void PrintAlignmentSettings(SRAArguments * sArgs) {

    SRASetting * sraSettings = sArgs->AlgnmtSetting;
    SRAIndex * sraIndex = sArgs->AlgnmtIndex;
    
    switch (sraSettings->OutputType) {
        case SRA_REPORT_UNIQUE_BEST:
            printf("Report         = Unique Best Alignment.\n");
            break;
        case SRA_REPORT_RANDOM_BEST:
            printf("Report         = Random Best Alignment.\n");
            break;
        case SRA_REPORT_BEST_QUALITY:
            printf("Report         = Best Quality Alignment.\n");
            break;
        case SRA_REPORT_ALL_BEST:
            if (sraSettings->MaxResult==-1) {
                printf("Report         = All Best Alignment (All).\n");
            } else {
                printf("Report         = All Best Alignment (#<=%d).\n",sraSettings->MaxResult);
            }
            break;
        case SRA_REPORT_ALL:
            if (sraSettings->MaxResult==-1) {
                printf("Report         = All Alignment (All).\n");
            } else {
                printf("Report         = All Alignment (#<=%d).\n",sraSettings->MaxResult);
            }
            break;
        case SRA_REPORT_DETAIL:
            printf("Report         = Unique Best Alignment /w Detail (For SNP Analysis).\n");
            break;
    }
    
    if (sraSettings->ErrorType==SRA_TYPE_MISMATCH_ONLY) {
        printf("Max #Mismatch  = %d\n",sraSettings->MaxError);
    } else if  (sraSettings->ErrorType==SRA_TYPE_EDIT_DISTANCE) {
        printf("Max #Edit      = %d\n",sraSettings->MaxError);
    }
    printf("Max #NBMism    = %d\n",sraSettings->MaxNBMismatch);
    printf("Ref.Seq Length = %u\n",sraIndex->bwt->textLength);
    printf("Strand         = ");
    switch (sraSettings->ReadStrand) {
        case QUERY_POS_STRAND:
            printf("Positive\n");
            break;
        case QUERY_NEG_STRAND:
            printf("Negative\n");
            break;
        case QUERY_BOTH_STRAND:
            printf("Both\n");
            break;
    }
    printf("Threads        = %d / %d\n",InputNumOfCPUThreads,InputNumOfMICThreads);
    
    /*
    printf("Pre-Align Trim = ");
    if (TrimPrealignRead>0) {
        printf("Enabled (trimming=%u%c)\n",TrimPrealignRead ,'%');
    } else {
        printf("Disabled\n");
    }
    printf("2-Pass         = ");
    if (SRASecondPassEnabled) {
        printf("Enabled (trimming=%u%c loosen=%d forceIndel=%d)\n",TrimUnalignRead,'%',LooseCriteriaUnalignRead,ForceIndelUnalignedRead);
    } else {
        printf("Disabled\n");
    }
    printf("SGA Seed Rvy   = ");
    if (InputSGASeedEnhancement) {
        printf("%.2f / %.2f\n",InputSGASeedLength, InputSGASeedLengthOverlap);
        printf("               = loosen=%d\n",InputSGASeedLooseCriteria);
        if (InputSGASeedOccLimitation!=-1)
            printf("               = limit'd=%d\n",InputSGASeedOccLimitation);
    } else {
        printf("Disabled\n");
    }*/
    
    switch (sraSettings->OutputFormat) {
        case SRA_OUTPUT_FORMAT_PLAIN:
            printf("Output Format  = Succinct (Plain).\n");
            break;
        case SRA_OUTPUT_FORMAT_SAM:
            printf("Output Format  = SAM v1.4 (Plain).\n");
            break;
        case SRA_OUTPUT_FORMAT_BAM:
            printf("Output Format  = BAM v1.4 (Binary).\n");
            break;
    }
    printf("\n");
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

    // Memory parameters
    InputMemoryMICMaxUsageMBytes = iniparser_getuint(ini, "Memory:MemoryMICMaxUsageMBytes", InputMemoryMICMaxUsageMBytes);
    InputMemoryCPUMaxUsageMBytes = iniparser_getuint(ini, "Memory:MemoryCPUMaxUsageMBytes", InputMemoryCPUMaxUsageMBytes);

    // Database parameters
    iniparser_copystring(ini, "Database:SaValueFileName", InputSaValueFileName, InputSaValueFileName, MAX_FILENAME_LEN);
    iniparser_copystring(ini, "Database:CPUSaValueFileName", InputCPUSaValueFileName, InputCPUSaValueFileName, MAX_FILENAME_LEN);

    //Query Parameters
    InputNumOfMICThreads = iniparser_getuint(ini, "MultipleThreading:NumOfMICThreads", InputNumOfMICThreads );
    InputNumOfCPUThreads = iniparser_getuint(ini, "MultipleThreading:NumOfCPUThreads", InputNumOfCPUThreads );
    
    // Short Read Alignement Parameters
    char tmp[10];
    iniparser_copystring(ini, "AlignmentModel:AlignmentModel", tmp, tmp, 4);
    if (strcmp(tmp,"16G")==0) {
        InputAlignmentModel = SRA_MODEL_16G;
    } else {
        InputAlignmentModel = SRA_MODEL_8G;
    }

    // Read input handling
    iniparser_copystring(ini, "AlignmentModel:InvalidReadHandling", tmp, tmp, 8);
    if (strcmp(tmp,"REPLACE")==0) {
        InputSRAInvalidReadHandling = SRA_READ_REPLACE_INVALID;
    } else {
        InputSRAInvalidReadHandling = SRA_READ_SKIP_INVALID;
    } 
    
    InputMaxNumOfAlignment = iniparser_getdouble(ini, "SingleEnd:MaxNumOfAlignment", InputMaxNumOfAlignment );
    
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

    if (argc<4) {
        PrintHelp();
        exit(1);
    }

    if (strcmp(argv[1],"single")!=0) {
        PrintHelp();
        exit(1);
    }

    // Get database, query name and output file name
    if (!iniparser_find_entry(programInput, "argument:2")) {
        PrintHelp();
        exit(1);
    }
    iniparser_copystring(programInput, "argument:2", DatabaseName, DatabaseName, MAX_FILENAME_LEN);

    if (!iniparser_find_entry(programInput, "argument:3")) {
        PrintHelp();
        exit(1);
    }
    iniparser_copystring(programInput, "argument:3", ReadFileName, ReadFileName, MAX_FILENAME_LEN);

    InputQueryStrand = iniparser_getint(programInput, "parameter:-S", InputQueryStrand);

    if (iniparser_find_entry(programInput, "parameter:-m")) {
        tempString = iniparser_getstr(programInput, "parameter:-m");
        len = (int)strlen(tempString);
        InputErrorType = SRA_TYPE_MISMATCH_ONLY;
        if ( len > 0 && strncmp(tempString + len - 1, "e", 1) == 0) {
            len--;
            tempString[len] = '\0';
            InputErrorType = SRA_TYPE_EDIT_DISTANCE;
        }
        if ( len > 0 ) {
            for (i=0;i<len;i++) {
                if (tempString[i]<'0' || tempString[i]>'9') {
                    PrintHelp();
                    exit(1);
                }
            }
            InputMaxError = atoi(tempString);
        }
    }

    InputMaxNBMismatch = iniparser_getint(programInput,"parameter:-n", InputMaxNBMismatch);

    InputOutputFileName = iniparser_find_entry(programInput, "parameter:-o");
    if (InputOutputFileName) {
        iniparser_copystring(programInput, "parameter:-o", OutputFileName, OutputFileName, MAX_FILENAME_LEN);
    } else {
        InputOutputFileName=MIC_MAIN_NON_INPUT_VALUE;
    }
    OutputFormat = SRA_OUTPUT_FORMAT_PLAIN;
    OutputFormat = iniparser_getint(programInput, "parameter:-b", OutputFormat);

    return programInput;

}


void PrintHelp() {
    printf("Syntax:\n");
    printf("    %s load <2bwt index>\n",PROJECT_ALIGNER_BINARY);
    printf("    %s unload\n",PROJECT_ALIGNER_BINARY);
    printf("\n");
    printf("Short Read Alignment Syntax:\n");
    printf("    %s single <2bwt index> <read file> [Options]\n",PROJECT_ALIGNER_BINARY);
    printf("\n");
    printf("    [Options]\n");
    printf("        -m: Maximum #errors allowed. [Def=Auto]\n");
    printf("            Expect value in this format \"-m [e]<intNum>\".\n");
    printf("            [e] is an optional flag; If supplied, error is defined as edit (mis+indel).\n");
    printf("            Otherwise, as mismatch. <intNum> is the maximum number of errors\n");
    printf("            allowed. Valid values are 0,..,%d or e0,..,e%d\n",MAX_NUM_OF_MISMATCH,MAX_NUM_OF_INDEL);
    printf("            If above is not supplied, SOAP2 will auto-default\n");
    printf("            the alignment criteria by inspecting the read files.\n");
    printf("            See Readme for more information.\n");
    printf("        -S: %d - Upper; %d - Lower; %d - Both [Def=%d]\n",
                        QUERY_POS_STRAND,
                        QUERY_NEG_STRAND,
                        QUERY_BOTH_STRAND,
                        QUERY_BOTH_STRAND);
    printf("        -h: Alignment type. [Def=%d]\n",SRA_REPORT_RANDOM_BEST);
    printf("            %d : All Valid Alignment\n", SRA_REPORT_ALL);
    printf("            %d : All Best Alignment\n", SRA_REPORT_ALL_BEST);
    printf("            %d : Unique Best Alignment\n", SRA_REPORT_UNIQUE_BEST);
    printf("            %d : Random Best Alignment\n", SRA_REPORT_RANDOM_BEST);
    printf("            %d : Best Quality Alignment\n", SRA_REPORT_BEST_QUALITY);
    printf("        -b: Output format. [Def=%d]\n",SRA_OUTPUT_FORMAT_PLAIN);
    printf("            %d : Succinct (Plain)\n", SRA_OUTPUT_FORMAT_PLAIN);
    printf("            %d : SAM 1.4 (Plain)\n", SRA_OUTPUT_FORMAT_SAM);
    printf("\n");
}

void SRAFillCharMap(unsigned char * charMap) {
    int i;
    memset(charMap,2,256);
    for (i=0;i<ALPHABET_SIZE;i++) {
        charMap[dnaChar[i]] = i;
    }
}

void SRAGetFileSuffix(int threadId, char * suffixStr) {
    if (threadId<0 || threadId>=1000) {
        threadId=0;
    }
    suffixStr[0]='.';
    suffixStr[1]='0'+threadId / 100;
    suffixStr[2]='0'+(threadId % 100) / 10;
    suffixStr[3]='0'+(threadId % 10);
    suffixStr[4]='\0';
}
