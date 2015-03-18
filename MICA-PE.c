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

#include "MICA-PE.h"
#include "Release.h"

#define MIC_MAIN_NON_INPUT_VALUE    99999999

// Input Arguments for MICA
#define MICA_INPUT_FILE_ARG_FILE    0
#define MICA_INPUT_FILE_ARG_LIST    1

// -----------------------------
// Compile-time Configurable #1
// -----------------------------
// Maximum number of reads to be processed by
// each MIC controller batch run.
#define MIC_MAX_BATCH_SIZE              1000000

// -----------------------------
// Compile-time Configurable #2
// -----------------------------
// Number of MIC card equipped in the machine.
// User may set the numOfMICThread to 0 in mica-pe.ini to disable
// a particular MIC card. 
#define NUM_MIC_EQUIPPED            3

// -----------------------------
// Compile-time Configurable #3
// -----------------------------
// Enable the core health checker to poll
// multiple component and get health status on stdout.
//#define MICA_PE_ENABLE_HEALTH_CHECK

//-----------------------------------
// DEBUG FLAGS
//-----------------------------------
// Define to print output statistics
//#define MIC_DEBUG_PRINT_OUTPUT_COUNTS


// Declaration of Global Variables
// Global Variables should NEVER be referenced outside of
// this file though they can be.

// Number of CPU thread allocated to aid MIC alignment
// It's common for nowadays computer to equip with a multi-
// core CPU. MICA supports to definition of extra aligner
// utilising CPU processing power to perform alignment
// alongside the MIC processing unit.
int InputNumOfCPUThreads = 1;

int InputFileArgument = MICA_INPUT_FILE_ARG_FILE;
char ReadFileName[MAX_FILENAME_LEN+1] = "";
char MateFileName[MAX_FILENAME_LEN+1] = "";
char ListFileName[MAX_FILENAME_LEN+1] = "";
char OutputFileName[MAX_FILENAME_LEN+1] = "*.out";
char DatabaseName[MAX_FILENAME_LEN+1] = "";
char InputSaValueFileName[MAX_FILENAME_LEN+1] = ".sa8";
char InputCPUSaValueFileName[MAX_FILENAME_LEN+1] = ".sa";

int InputCPULoadSAValue = 1;

int InputSRAInvalidReadHandling = SRA_READ_SKIP_INVALID;
int InputQueryStrand = QUERY_BOTH_STRAND;
int InputAlignmentModel = SRA_MODEL_8G;
int InputMaxNumOfAlignment = -1;

int InputReportType = PE_REPORT_ALL_BEST;
int InputErrorType = SRA_TYPE_MISMATCH_ONLY;
int InputMaxError = 2;
int InputMaxNBMismatch = 0;
int InputOutputFileName = MIC_MAIN_NON_INPUT_VALUE;
int InputOutputFormat = SRA_OUTPUT_FORMAT_SAM;
int InputOutputNumSamThreads = 1;

int InputPEInsertionUpperBound = MIC_MAIN_NON_INPUT_VALUE;
int InputPEInsertionLowerBound = MIC_MAIN_NON_INPUT_VALUE;
int InputPEStrandLeftLeg = QUERY_POS_STRAND;
int InputPEStrandRightLeg = QUERY_NEG_STRAND;

int    InputSGAOrphanEnhancement = 0;
double InputSGAOrphanTriggerTF = 0.30f;
double InputSGAScoreTF = 0.30f;

int    InputSGAOrphanExtendEnhancement = 0;
double InputSGAOrphanExtendTriggerTF = 0.30f;

int    InputSGAScoreMatch     = 1;
int    InputSGAScoreMismatch  = -2;
int    InputSGAScoreGapOpen   = -3;
int    InputSGAScoreGapExtend = -1;

double InputSGASoftHeadClipLength = 0.0;
double InputSGASoftTailClipLength = 0.40;
double InputSGASoftTotalClipLength = 0.40;

int    InputSGASeedEnhancement = 0;
double InputSGASeedLength = 0.30f;
double InputSGASeedLengthOverlap = 0.15f;
int    InputSGASeedLooseCriteria = -1;
int    InputSGASeedOccLimitation = 4096;

int    InputSGASeedEnhancement_1 = 0;
double InputSGASeedLength_1 = 0.30f;
double InputSGASeedLengthOverlap_1 = 0.15f;
int    InputSGASeedLooseCriteria_1 = -1;
int    InputSGASeedOccLimitation_1 = 4096;

unsigned int InputMemoryCPUMaxUsageMBytes = -1;
unsigned int InputMemoryMICMaxUsageMBytes = 3600;

unsigned int InputMICThreads[NUM_MIC_EQUIPPED];
    
void PrintHelp();
dictionary *ParseInput(int argc, char** argv);
void ParseIniFile(char *iniFileName);
void SRAFillCharMap(unsigned char * charMap);
void PrintAlignmentSettings(Logging * logger, SRAArguments * sArgs, PEArguments * peArgs);
void SRAGetFileSuffix(int threadId, char * suffixStr);
void ParameterVerification();

int main(int argc, char** argv) {

    // Start Audit Logging
    Logging * logger = LOGCreate("mica-pe.log",LOGGING_MODE_RUNNING);

    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL,
                "\n%s v%d.%d.%d (%s)\n",PROJECT_NAME,PROJECT_MAJOR,PROJECT_MINOR,PROJECT_REV,PROJECT_SPECIAL);
    
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

    // Parameter Control
    ParameterVerification();

    ////////////////////////////////////////////////////////////
    // Initialise the threads information
    ////////////////////////////////////////////////////////////
    int gThreadId = 0;
    int micThreadId = 0;
    int cpuThreadId = 0;
    int maxNumOfThreads = 0;
    int numMICActivated = 0;
    for (micThreadId = 0; micThreadId < NUM_MIC_EQUIPPED; micThreadId++) {
        numMICActivated += InputMICThreads[micThreadId] != 0;
        if ( InputMICThreads[micThreadId] > maxNumOfThreads ) {
            maxNumOfThreads = InputMICThreads[micThreadId];
        }
    }

    ////////////////////////////////////////////////////////////
    // Index Handling
    ////////////////////////////////////////////////////////////
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL,
                "Loading index %s ... ", DatabaseName); fflush(stdout);
    Idx2BWT * micIdx2BWT = BWTLoad2BWT(DatabaseName,InputSaValueFileName);
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL,
                "DONE\n\n");

    // Jeanno: The CPUIdx2BWT Index creation here is a bit clumsy.
    // Maybe we can work on a module for this later.
    Idx2BWT * cpuIdx2BWT = NULL;
    if (InputCPULoadSAValue) {
        // Only load CPU index if the SA being used is different to the
        // one being used by MIC.
        cpuIdx2BWT = (Idx2BWT *) MEMManMalloc(sizeof(Idx2BWT), MEMORY_TYPE_CPU);
        (*cpuIdx2BWT) = (*micIdx2BWT);
        cpuIdx2BWT->bwt = (BWT *) MEMManMalloc(sizeof(BWT), MEMORY_TYPE_CPU);
        (*cpuIdx2BWT->bwt) = (*micIdx2BWT->bwt);

        // Full Suffix Array for CPU
        char saFilename[MAX_INDEX_FILENAME_LENGTH]; 
        strcpy(saFilename, DatabaseName);
        strcat(saFilename, InputCPUSaValueFileName);
        printf("Loading Full Suffix-Array for CPU %s ... ",saFilename); fflush(stdout);
        BWTSALoad(cpuIdx2BWT->mmPool, cpuIdx2BWT->bwt, saFilename,NULL);
        printf("DONE\n\n");
    } else {
        cpuIdx2BWT = micIdx2BWT;
    }

    // Lookup Table for CPU
    char LookupTableFileName[MAX_FILENAME_LEN+1];
    char RevLookupTableFileName[MAX_FILENAME_LEN+1];
    strcpy(LookupTableFileName,DatabaseName);
    strcpy(RevLookupTableFileName,DatabaseName);
    strcat(LookupTableFileName,".lkt");
    strcat(RevLookupTableFileName,".rev.lkt");
    LT * lookup = LTLoad(LookupTableFileName);
    LT * revLookup = LTLoad(RevLookupTableFileName);    

    // Timestamp'ed -----------
    timestamp = getElapsedTime(startTime);
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL,
                "(Elapsed time (Index Loading) : %9.4f seconds)\n\n", timestamp - lastEventTime);
    lastEventTime = timestamp;
    // ------------------------

    HSPFillCharMap(charMap);
    if (InputSRAInvalidReadHandling == SRA_READ_REPLACE_INVALID) SRAFillCharMap(charMap);
    HSPFillComplementMap(complementMap);

    #ifdef MICA_PE_ENABLE_HEALTH_CHECK
        // Construct the health checker procedure
        MICAHealthParm * micaHealth = MICAHealthParmCreate();
    #endif
    
    printf("Initialising Read Thread ... "); fflush(stdout);
    ////////////////////////////////////////////////////////////
    // Read Thread
    ////////////////////////////////////////////////////////////
    ListConsumer * listConsumer = LCCreate();
    if ( InputFileArgument == MICA_INPUT_FILE_ARG_FILE ) {
        LCCreateList(listConsumer,
            ReadFileName, MateFileName,
            &InputPEInsertionLowerBound, &InputPEInsertionUpperBound,
            1, MAX_FILENAME_LEN,MAX_FIELD_LEN);
    } else if ( InputFileArgument == MICA_INPUT_FILE_ARG_LIST ) {
        LCLoadList(listConsumer,ListFileName,MAX_FILENAME_LEN,MAX_FIELD_LEN);
    } else {
        printf("[SAFE-GUARD] List consumer not initialised due to unexpected input type.\n");
    }
    
    PEReadThread * readThread = PEReadThread_create(charMap, 
        MIC_MAX_BATCH_SIZE, MAX_SEQ_NAME_LENGTH,
        SRA_MAX_READ_LENGTH, 
        cpuIdx2BWT->hsp, listConsumer, InputOutputFormat, InputOutputNumSamThreads);
    PEReadThread_start(readThread);
    printf("DONE\n");

    #ifdef MICA_PE_ENABLE_HEALTH_CHECK
        MICAAddMonitor_ReadThread(micaHealth,readThread);
    #endif
    
    MICControllerThread * micThreads[NUM_MIC_EQUIPPED];
    CPUControllerThread ** cpuThreads = (CPUControllerThread**) malloc(sizeof(CPUControllerThread*)*InputNumOfCPUThreads);
    
    ////////////////////////////////////////////////////////////
    // Set up shared arguments and parameters for SRA and PE - CPU
    ////////////////////////////////////////////////////////////
    SRAArguments * cpuSraArgTemplate = SRAARGConstruct();
    PEArguments * cpuPeArgTemplate = PEARGConstruct(cpuIdx2BWT->bwt,cpuIdx2BWT->hsp);

    // Setting up SRA Indexes
    SRAIndexPopulate(cpuSraArgTemplate->AlgnmtIndex, cpuIdx2BWT->bwt, cpuIdx2BWT->rev_bwt,
        cpuIdx2BWT->hsp, NULL, lookup, revLookup);
        
    // Setting up SRA Settings
    // The Pair-end aligner requires SRA to behave differently in different cases.
    // For example, when we are performing all-valid PE we don't need SRA to be staged
    // by different number of mismatch; Thus in that case we merely need to perform a
    // SRA_REPORT_ALL. In cases like all-best alignment, it's a requirement to have it
    // set toe SRA_REPORT_ALL_SORTED to allow the pair-end alignment to early terminate.
    if ( InputReportType == PE_REPORT_ALL_BEST ) {
        cpuSraArgTemplate->AlgnmtSetting->OutputType = SRA_REPORT_ALL_SORTED;
    } else {
        cpuSraArgTemplate->AlgnmtSetting->OutputType = SRA_REPORT_ALL;
    }
    cpuSraArgTemplate->AlgnmtSetting->ErrorType = InputErrorType;
    cpuSraArgTemplate->AlgnmtSetting->MaxError = InputMaxError;
    cpuSraArgTemplate->AlgnmtSetting->MaxNBMismatch = InputMaxNBMismatch;
    
    // Setting up SRA Output Settings
    cpuSraArgTemplate->AlgnmtSetting->OutputFormat = InputOutputFormat;

    // Setting up PE Settings
    cpuPeArgTemplate->PEAlgnmtInput->OutputType        = InputReportType;
    cpuPeArgTemplate->PEAlgnmtInput->maxResult         = InputMaxNumOfAlignment;
    cpuPeArgTemplate->PEAlgnmtInput->insertLbound      = InputPEInsertionLowerBound;
    cpuPeArgTemplate->PEAlgnmtInput->insertUbound      = InputPEInsertionUpperBound;
    cpuPeArgTemplate->PEAlgnmtInput->strandLeftLeg     = InputPEStrandLeftLeg;
    cpuPeArgTemplate->PEAlgnmtInput->strandRightLeg    = InputPEStrandRightLeg;
    
    // Setting up PE-DP Settings
    cpuPeArgTemplate->pedpSetting->SGAScoreTF              = InputSGAScoreTF;
    cpuPeArgTemplate->pedpSetting->SGASoftHeadClipLength   = InputSGASoftHeadClipLength;
    cpuPeArgTemplate->pedpSetting->SGASoftTailClipLength   = InputSGASoftTailClipLength;
    cpuPeArgTemplate->pedpSetting->SGASoftTotalClipLength  = InputSGASoftTotalClipLength;
    
    cpuPeArgTemplate->pedpSetting->SGAOrphanEnhancement    = InputSGAOrphanEnhancement;
    cpuPeArgTemplate->pedpSetting->SGAOrphanTriggerTF      = InputSGAOrphanTriggerTF;
    
    cpuPeArgTemplate->pedpSetting->SGAOrphanExtendEnhancement = InputSGAOrphanExtendEnhancement;
    cpuPeArgTemplate->pedpSetting->SGAOrphanExtendTriggerTF   = InputSGAOrphanExtendTriggerTF;
    
    cpuPeArgTemplate->dpArguments->dpScores->dpMatch       = InputSGAScoreMatch;
    cpuPeArgTemplate->dpArguments->dpScores->dpMismatch    = InputSGAScoreMismatch;
    cpuPeArgTemplate->dpArguments->dpScores->dpGapOpen     = InputSGAScoreGapOpen;
    cpuPeArgTemplate->dpArguments->dpScores->dpGapExtend   = InputSGAScoreGapExtend;
    
    cpuPeArgTemplate->pedpSetting->SGASeedEnhancement      = InputSGASeedEnhancement;
    cpuPeArgTemplate->pedpSetting->SGASeedLength           = InputSGASeedLength;
    cpuPeArgTemplate->pedpSetting->SGASeedLengthOverlap    = InputSGASeedLengthOverlap;
    cpuPeArgTemplate->pedpSetting->SGASeedLooseCriteria    = InputSGASeedLooseCriteria;
    cpuPeArgTemplate->pedpSetting->SGASeedOccLimitation    = InputSGASeedOccLimitation;
    
    cpuPeArgTemplate->pedpSetting->SGASeedEnhancement_1    = InputSGASeedEnhancement_1;
    cpuPeArgTemplate->pedpSetting->SGASeedLength_1         = InputSGASeedLength_1;
    cpuPeArgTemplate->pedpSetting->SGASeedLengthOverlap_1  = InputSGASeedLengthOverlap_1;
    cpuPeArgTemplate->pedpSetting->SGASeedLooseCriteria_1  = InputSGASeedLooseCriteria_1;
    cpuPeArgTemplate->pedpSetting->SGASeedOccLimitation_1  = InputSGASeedOccLimitation_1;
    
    //Printing Settings
    if ( InputFileArgument == MICA_INPUT_FILE_ARG_FILE ) {
        LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL,
            "Input          = %s / %s\n",ReadFileName,MateFileName);
    } else if ( InputFileArgument == MICA_INPUT_FILE_ARG_LIST ) {
        LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL,
            "Input          = %s\n",ListFileName);
    } else if ( InputFileArgument == MICA_INPUT_FILE_ARG_LIST ) {
        LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL,
            "Input          = Unknown\n");
    }
    PrintAlignmentSettings(logger,cpuSraArgTemplate,cpuPeArgTemplate);

    ////////////////////////////////////////////////////////////
    // Write Thread
    ////////////////////////////////////////////////////////////
    printf("Initialising Write Thread ... "); fflush(stdout);
    unsigned int maxNumReadPerCore = (MIC_MAX_BATCH_SIZE + maxNumOfThreads - 1) / maxNumOfThreads;
    PEWriteThread * gWriteThread = PEWTCreate(maxNumReadPerCore, maxNumOfThreads, 
                                              MIC_CONTROL_WRITE_THREAD_BUFFER_COUNT, 
                                              InputNumOfCPUThreads + NUM_MIC_EQUIPPED,
                                              0);

    PEWTRegisterArguments(gWriteThread, cpuSraArgTemplate, cpuPeArgTemplate,
                            InputAlignmentModel,
                            charMap);
    
    PEWTRoll(gWriteThread);
    
    printf("DONE\n");
    
    #ifdef MICA_PE_ENABLE_HEALTH_CHECK
        MICAAddMonitor_WriteThread(micaHealth,gWriteThread);
    #endif
    
    ////////////////////////////////////////////////////////////
    // Create multiple MIC Threads
    ////////////////////////////////////////////////////////////
    printf("Initialising MIC Controller (%u) ... ", numMICActivated); fflush(stdout);
    for (micThreadId = 0; micThreadId < NUM_MIC_EQUIPPED; micThreadId++) {
        if (InputMICThreads[micThreadId] == 0) {
            continue;
        }
        micThreads[micThreadId] = MICCTCreate(gThreadId, micIdx2BWT, charMap,
            MIC_MAX_BATCH_SIZE, micThreadId, InputMICThreads[micThreadId]);

        MICControllerThread * micThread = micThreads[micThreadId];
        MICCTSetSRA(micThread, InputAlignmentModel, cpuSraArgTemplate);
        MICCTSetPE(micThread, cpuPeArgTemplate);
        MICCTSetReadThread(micThread, readThread);
        MICCTSetWriteThread(micThread, gWriteThread);
        MICCTStart(micThread);
        #ifdef MICA_PE_ENABLE_HEALTH_CHECK
            MICAAddMonitor_MICControl(micaHealth,micThread);
        #endif
        
        gThreadId++;
    }
    printf("DONE\n");

    ////////////////////////////////////////////////////////////
    // Create multiple CPU Threads
    ////////////////////////////////////////////////////////////
    printf("Initialising CPU Controller (%u) ... ",InputNumOfCPUThreads); fflush(stdout);
    for (cpuThreadId = 0; cpuThreadId < InputNumOfCPUThreads; cpuThreadId++) {
        cpuThreads[cpuThreadId] = CPUCTCreate(gThreadId, charMap, MIC_MAX_BATCH_SIZE);

        CPUControllerThread * cpuThread = cpuThreads[cpuThreadId];
        CPUCTSetSRA(cpuThread, InputAlignmentModel, cpuSraArgTemplate);
        CPUCTSetPE(cpuThread, cpuPeArgTemplate);
        CPUCTSetReadThread(cpuThread, readThread);
        CPUCTSetWriteThread(cpuThread, gWriteThread);
        CPUCTStart(cpuThread);
        #ifdef MICA_PE_ENABLE_HEALTH_CHECK
            MICAAddMonitor_CPUControl(micaHealth,cpuThread);
        #endif
        
        gThreadId++;
    }
    printf("DONE\n");
    
    printf("\n");
    
    #ifdef MICA_PE_ENABLE_HEALTH_CHECK
        MICAHealthParmRoll(micaHealth);
    #endif
    
    // Join and Free micThreads
    for (micThreadId = 0; micThreadId < NUM_MIC_EQUIPPED; micThreadId++) {
        if (InputMICThreads[micThreadId] == 0) {
            continue;
        }

        MICCTJoin(micThreads[micThreadId]);
    }

    // Join and Free cpuThreads
    for (cpuThreadId = 0; cpuThreadId < InputNumOfCPUThreads; cpuThreadId++) {
        CPUCTJoin(cpuThreads[cpuThreadId]);
    }

    // Join and Free micThreads
    for (micThreadId = 0; micThreadId < NUM_MIC_EQUIPPED; micThreadId++) {
        if (InputMICThreads[micThreadId] == 0) {
            continue;
        }

        MICCTPrintStats(micThreads[micThreadId],logger);
        MICCTFree(micThreads[micThreadId]);
    }

    // Join and Free cpuThreads
    for (cpuThreadId = 0; cpuThreadId < InputNumOfCPUThreads; cpuThreadId++) {
        CPUCTPrintStats(cpuThreads[cpuThreadId],logger);
        CPUCTFree(cpuThreads[cpuThreadId]);
    }

    free(cpuThreads);
    
    #ifdef MICA_PE_ENABLE_HEALTH_CHECK
        // Free the health checker
        MICAHealthParmJoin(micaHealth);
        MICAHealthParmFree(micaHealth);
    #endif

    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL,
                "[INFO] Total Shared Working Memory allocated = %9.2f Mbytes\n",MEMManGetUsage(MEMORY_TYPE_SHARED)/1024.0/1024.0);
    if (InputMemoryMICMaxUsageMBytes!=-1 && MEMManGetUsage(MEMORY_TYPE_SHARED)>(unsigned long long) InputMemoryMICMaxUsageMBytes*1024*1024) {
        LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL,
                    "[ERROR] Exceeded defined MIC limitation of %u MBytes\n",InputMemoryMICMaxUsageMBytes);
    }
    
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL,
                "[INFO] Total CPU Working Memory allocated = %9.2f Mbytes\n",MEMManGetUsage(MEMORY_TYPE_CPU)/1024.0/1024.0);
    if (InputMemoryCPUMaxUsageMBytes!=-1 && MEMManGetUsage(MEMORY_TYPE_CPU)>(unsigned long long) InputMemoryCPUMaxUsageMBytes*1024*1024) {
        LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL,
                    "[ERROR] Exceeded defined CPU limitation of %u MBytes\n",InputMemoryMICMaxUsageMBytes);
    }
    
    printf("\n");

    PEWTSignal(gWriteThread, MICA_PE_WRITE_THREAD_HEALTH_DEPLETED);
    PEWTThreadJoin(gWriteThread);
    PEWTPrintStats(gWriteThread,logger);
    
    PEReadThread_quit(readThread);
    PEReadThread_join(readThread);

    // Timestamp'ed -----------
    timestamp = getElapsedTime(startTime);
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL,
                "Elapsed time (Total MIC and CPU Processing) : %9.4f seconds\n\n", timestamp - lastEventTime);
    lastEventTime = timestamp;
    // ------------------------


    ////////////////////////////////////////////////////////////
    // Free Memory Allocated
    ////////////////////////////////////////////////////////////
    printf("Free index ... ");
    fflush(stdout);
    
    PEARGFree(cpuPeArgTemplate);
    SRAARGFree(cpuSraArgTemplate);
    
    LTFree(lookup);
    LTFree(revLookup);
    if (InputCPULoadSAValue) {
        // Free up CPU memory only if it's branched off from MIC Index.
        BWTSAFree(cpuIdx2BWT->mmPool, cpuIdx2BWT->bwt);
        free(cpuIdx2BWT->bwt);
        free(cpuIdx2BWT);
    }
    BWTFree2BWT(micIdx2BWT);
    printf("DONE\n");

    printf("Free Read/Write Thead ... ");
    fflush(stdout);
    PEReadThread_free(readThread);
    PEWTFree(gWriteThread);
    LCFree(listConsumer);
    
    printf("DONE\n");

    LOGFree(logger);
    iniparser_freedict(programInput); 
    return 0;
}


void PrintAlignmentSettings(Logging * logger, SRAArguments * sArgs, PEArguments * peArgs) {

    SRASetting * sraSettings = sArgs->AlgnmtSetting;
    SRAIndex * sraIndex = sArgs->AlgnmtIndex;
    
    switch (peArgs->PEAlgnmtInput->OutputType) {
        case PE_REPORT_ALL_BEST:
            if (peArgs->PEAlgnmtInput->maxResult==-1) {
                LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Report         = All Best Alignment (All).\n");
            } else {
                LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Report         = All Best Alignment (#<=%d).\n",peArgs->PEAlgnmtInput->maxResult);
            }
            break;
        case PE_REPORT_ALL:
            if (peArgs->PEAlgnmtInput->maxResult==-1) {
                LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Report         = All Alignment (All).\n");
            } else {
                LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Report         = All Alignment (#<=%d).\n",peArgs->PEAlgnmtInput->maxResult);
            }
            break;
    }
    
    if (sraSettings->ErrorType==SRA_TYPE_MISMATCH_ONLY) {
        LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Max #Mismatch  = %d\n",sraSettings->MaxError);
    } else if  (sraSettings->ErrorType==SRA_TYPE_EDIT_DISTANCE) {
        LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Max #Edit      = %d\n",sraSettings->MaxError);
    }
    LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Max #NBMism    = %d\n",sraSettings->MaxNBMismatch);
    LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Ref.Seq Length = %u\n",sraIndex->bwt->textLength);
    LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Strand         = ");
    switch (sraSettings->ReadStrand) {
        case QUERY_POS_STRAND:
            LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Positive\n");
            break;
        case QUERY_NEG_STRAND:
            LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Negative\n");
            break;
        case QUERY_BOTH_STRAND:
            LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Both\n");
            break;
    }

    char tmp[4] = "?+-";
    LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Pair strand    = %c/%c\n",tmp[InputPEStrandLeftLeg],tmp[InputPEStrandRightLeg]);
    
    LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"CPU Threads    = %d\n",InputNumOfCPUThreads);
    
    if (InputSGAOrphanEnhancement || InputSGAOrphanExtendEnhancement || InputSGASeedEnhancement) {
        LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"SGA Soft-clip  = %.2f / %.2f / %.2f \n", InputSGASoftHeadClipLength, InputSGASoftTailClipLength, InputSGASoftTotalClipLength);
    }
    if (InputSGAOrphanEnhancement) LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"SGA Orphan Rvy = %.2f / %.2f \n", InputSGAOrphanTriggerTF, InputSGAScoreTF);
    if (InputSGAOrphanExtendEnhancement) LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"SGA OrpExt Rvy = %.2f / %.2f \n", InputSGAOrphanExtendTriggerTF, InputSGAScoreTF);
    //if (InputSGAOrphanEnhancement) LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"SGA Orphan Rvy = %.2f / %.2f\n",InputSGAOrphanTriggerTF, InputSGAScoreTF);
    /*
    LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Pre-Align Trim = ");
    if (TrimPrealignRead>0) {
        LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Enabled (trimming=%u%c)\n",TrimPrealignRead ,'%');
    } else {
        LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Disabled\n");
    }
    LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"2-Pass         = ");
    if (SRASecondPassEnabled) {
        LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Enabled (trimming=%u%c loosen=%d forceIndel=%d)\n",TrimUnalignRead,'%',LooseCriteriaUnalignRead,ForceIndelUnalignedRead);
    } else {
        LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Disabled\n");
    }
    */
    if (InputSGASeedEnhancement) {
        LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"SGA Seed Rvy   = %.2f / %.2f\n",InputSGASeedLength, InputSGASeedLengthOverlap);
        LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"               = loosen=%d\n",InputSGASeedLooseCriteria);
        if (InputSGASeedOccLimitation!=-1)
            LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"               = limit'd=%d\n",InputSGASeedOccLimitation);
    }
    if (InputSGASeedEnhancement_1) {
        LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"SGA Seed Rvy-2 = %.2f / %.2f\n",InputSGASeedLength_1, InputSGASeedLengthOverlap_1);
        LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"               = loosen=%d\n",InputSGASeedLooseCriteria_1);
        if (InputSGASeedOccLimitation_1!=-1)
            LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"               = limit'd=%d\n",InputSGASeedOccLimitation_1);
    }
    
    switch (sraSettings->OutputFormat) {
        case SRA_OUTPUT_FORMAT_PLAIN:
            LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Output Format  = Succinct (Plain).\n");
            break;
        case SRA_OUTPUT_FORMAT_SAM:
            LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Output Format  = SAM v1.4 (Plain).\n");
            break;
        case SRA_OUTPUT_FORMAT_BAM:
            LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"Output Format  = BAM v1.4 (Binary).\n");
            break;
    }
    LOGWriteLine(logger,LOGGING_MSG_TYPE_ALL,"\n");
}


void ParseIniFile(char *iniFileName) {

    dictionary *ini;
    char tmp[10];

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
    
    if (strcmp(InputSaValueFileName,InputCPUSaValueFileName)==0) {
        InputCPULoadSAValue = 0;
    } else {
        InputCPULoadSAValue = 1;
    }

    //Query Parameters
    InputNumOfCPUThreads = iniparser_getuint(ini, "MultipleThreading:NumOfCPUThreads", InputNumOfCPUThreads );

    //Output Format Parameters 
    InputOutputNumSamThreads = iniparser_getuint(ini, "BAMSpecific:NumOfCompressionThreads", InputOutputNumSamThreads );

    // Short Read Alignement Parameters
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
    
    InputMaxNumOfAlignment = iniparser_getdouble(ini, "PairEnd:MaxNumOfAlignment", InputMaxNumOfAlignment );
    
    iniparser_copystring(ini, "PairEnd:StrandArrangement", tmp, tmp, 10);
    if (strcmp(tmp,"+/+")==0) {
        InputPEStrandLeftLeg = QUERY_POS_STRAND;
        InputPEStrandRightLeg = QUERY_POS_STRAND;
    } else if (strcmp(tmp,"-/+")==0) {
        InputPEStrandLeftLeg = QUERY_NEG_STRAND;
        InputPEStrandRightLeg = QUERY_POS_STRAND;
    } else if (strcmp(tmp,"-/-")==0) {
        InputPEStrandLeftLeg = QUERY_NEG_STRAND;
        InputPEStrandRightLeg = QUERY_NEG_STRAND;
    } else {
        InputPEStrandLeftLeg = QUERY_POS_STRAND;
        InputPEStrandRightLeg = QUERY_NEG_STRAND;
    }
    
    // DP Common
    InputSGAScoreTF = iniparser_getdouble(ini, "PairEndSGAEnhancement:SGAScoreTF", InputSGAScoreTF );
    InputSGAScoreMatch     = iniparser_getint(ini, "PairEndSGAEnhancement:SGAScoreMatch", InputSGAScoreMatch );
    InputSGAScoreMismatch  = iniparser_getint(ini, "PairEndSGAEnhancement:SGAScoreMismatch", InputSGAScoreMismatch );
    InputSGAScoreGapOpen   = iniparser_getint(ini, "PairEndSGAEnhancement:SGAScoreGapOpen", InputSGAScoreGapOpen );
    InputSGAScoreGapExtend = iniparser_getint(ini, "PairEndSGAEnhancement:SGAScoreGapExtend", InputSGAScoreGapExtend );
    
    InputSGASoftHeadClipLength = iniparser_getdouble(ini, "PairEndSGAEnhancement:SGASoftHeadClipLength", InputSGASoftHeadClipLength );
    InputSGASoftTailClipLength = iniparser_getdouble(ini, "PairEndSGAEnhancement:SGASoftTailClipLength", InputSGASoftTailClipLength );
    InputSGASoftTotalClipLength = iniparser_getdouble(ini, "PairEndSGAEnhancement:SGASoftTotalClipLength", InputSGASoftTotalClipLength );
    
    // Orhpan Recovery
    iniparser_copystring(ini, "PairEndSGAEnhancement:SGAOrphanEnhancement", tmp, tmp, 6);
    if (strcmp(tmp,"TRUE")==0) {
        InputSGAOrphanEnhancement = 1;
    } else {
        InputSGAOrphanEnhancement = 0;
    } 
    InputSGAOrphanTriggerTF = iniparser_getdouble(ini, "PairEndSGAEnhancement:SGAOrphanTriggerTF", InputSGAOrphanTriggerTF );
    
    iniparser_copystring(ini, "PairEndSGAEnhancement:SGAOrphanExtendEnhancement", tmp, tmp, 6);
    if (strcmp(tmp,"TRUE")==0) {
        InputSGAOrphanExtendEnhancement = 1;
    } else {
        InputSGAOrphanExtendEnhancement = 0;
    } 
    InputSGAOrphanExtendTriggerTF = iniparser_getdouble(ini, "PairEndSGAEnhancement:SGAOrphanExtendTriggerTF", InputSGAOrphanExtendTriggerTF );
    
    // Seed Enhancement
    iniparser_copystring(ini, "PairEndSGAEnhancement:SGASeedEnhancement", tmp, tmp, 6);
    if (strcmp(tmp,"TRUE")==0) {
        InputSGASeedEnhancement = 1;
    } else {
        InputSGASeedEnhancement = 0;
    }
    InputSGASeedLength = iniparser_getdouble(ini, "PairEndSGAEnhancement:SGASeedLength", InputSGASeedLength );
    InputSGASeedLengthOverlap = iniparser_getdouble(ini, "PairEndSGAEnhancement:SGASeedLengthOverlap", InputSGASeedLengthOverlap );
    InputSGASeedLooseCriteria = iniparser_getint(ini, "PairEndSGAEnhancement:SGASeedLooseCriteria", InputSGASeedLooseCriteria );
    InputSGASeedOccLimitation = iniparser_getint(ini, "PairEndSGAEnhancement:SGASeedOccLimitation", InputSGASeedOccLimitation );

    // Seed Enhancement (2nd-level)
    iniparser_copystring(ini, "PairEndSGAEnhancement:SGASeedEnhancement_1", tmp, tmp, 6);
    if (strcmp(tmp,"TRUE")==0) {
        InputSGASeedEnhancement_1 = 1;
    } else {
        InputSGASeedEnhancement_1 = 0;
    }
    InputSGASeedLength_1 = iniparser_getdouble(ini, "PairEndSGAEnhancement:SGASeedLength_1", InputSGASeedLength_1 );
    InputSGASeedLengthOverlap_1 = iniparser_getdouble(ini, "PairEndSGAEnhancement:SGASeedLengthOverlap_1", InputSGASeedLengthOverlap_1 );
    InputSGASeedLooseCriteria_1 = iniparser_getint(ini, "PairEndSGAEnhancement:SGASeedLooseCriteria_1", InputSGASeedLooseCriteria_1 );
    InputSGASeedOccLimitation_1 = iniparser_getint(ini, "PairEndSGAEnhancement:SGASeedOccLimitation_1", InputSGASeedOccLimitation_1 );

    // MIC threads parameters
    InputMICThreads[0] = iniparser_getuint(ini, "MultipleThreading:NumOfMICThreads", 0);
    InputMICThreads[1] = iniparser_getuint(ini, "MultipleThreading:NumOfMICThreads_1", 0);
    InputMICThreads[2] = iniparser_getuint(ini, "MultipleThreading:NumOfMICThreads_2", 0);

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

    if (argc<5) {
        PrintHelp();
        exit(1);
    }

    if (strcmp(argv[1],"pair")!=0) {
        PrintHelp();
        exit(1);
    }
    
    // Get database, query name and output file name
    if (!iniparser_find_entry(programInput, "argument:2")) {
        PrintHelp();
        exit(1);
    }
    iniparser_copystring(programInput, "argument:2", DatabaseName, DatabaseName, MAX_FILENAME_LEN);
    
    if (iniparser_find_entry(programInput, "parameter:-i")) {
        // User input is a list of reads
        InputFileArgument = MICA_INPUT_FILE_ARG_LIST;
        iniparser_copystring(programInput, "parameter:-i", ListFileName, ListFileName, MAX_FILENAME_LEN);
    } else {
        // User input is a pair of file
        InputFileArgument = MICA_INPUT_FILE_ARG_FILE;
        if (!iniparser_find_entry(programInput, "argument:3")) {
            PrintHelp();
            exit(1);
        }
        iniparser_copystring(programInput, "argument:3", ReadFileName, ReadFileName, MAX_FILENAME_LEN);
        iniparser_copystring(programInput, "argument:4", MateFileName, MateFileName, MAX_FILENAME_LEN);
        InputPEInsertionLowerBound = iniparser_getint(programInput,"parameter:-v", InputPEInsertionLowerBound);
        InputPEInsertionUpperBound = iniparser_getint(programInput,"parameter:-u", InputPEInsertionUpperBound);
    }

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
    InputReportType = iniparser_getint(programInput,"parameter:-h", InputReportType);

    InputOutputFileName = iniparser_find_entry(programInput, "parameter:-o");
    if (InputOutputFileName) {
        iniparser_copystring(programInput, "parameter:-o", OutputFileName, OutputFileName, MAX_FILENAME_LEN);
    } else {
        InputOutputFileName=MIC_MAIN_NON_INPUT_VALUE;
    }
    InputOutputFormat = iniparser_getint(programInput, "parameter:-b", InputOutputFormat);

    return programInput;

}


void PrintHelp() {
    printf("Syntax:\n");
    printf("Pair End Alignment Syntax:\n");
    printf("    %s pair <2bwt index> <read file 1> <read file 2> -v <l-insert> -u <u-insert> [Options]\n",PROJECT_ALIGNER_BINARY);
    printf("    %s pair <2bwt index> -i <input list file> [Options]\n", PROJECT_ALIGNER_BINARY);
    printf("\n");
    printf("    [Mandatory]\n");
    printf("        -v: Lower bound;\n");
    printf("        -u: Upper bound of the insertion size between a pair.\n");
    printf("        A pair will be reported if the insertion size falls\n");
    printf("        between [v,u] and their strands match.\n");
    printf("\n");
    printf("    [Options]\n");
    printf("        -m: Maximum #errors allowed. [Def=2]\n");
    printf("            Expect value in this format \"-m <intNum>\".\n");
    printf("            <intNum> is the maximum number of errors\n");
    printf("            allowed. Valid values are 0,..,%d\n", MAX_NUM_OF_MISMATCH);
    printf("        -h: Alignment type. [Def=%d]\n",SRA_REPORT_ALL_BEST);
    printf("            %d : All Valid Alignment\n", SRA_REPORT_ALL);
    printf("            %d : All Best Alignment\n", SRA_REPORT_ALL_BEST);
    //printf("            %d : Unique Best Alignment\n", SRA_REPORT_UNIQUE_BEST);
    //printf("            %d : Random Best Alignment\n", SRA_REPORT_RANDOM_BEST);
    //printf("            %d : Best Quality Alignment\n", SRA_REPORT_BEST_QUALITY);
    printf("        -b: Output format. [Def=%d]\n",SRA_OUTPUT_FORMAT_SAM);
    //printf("            %d : Succinct (Plain)\n", SRA_OUTPUT_FORMAT_PLAIN);
    printf("            %d : SAM 1.4 (Plain)\n", SRA_OUTPUT_FORMAT_SAM);
    printf("            %d : BAM 1.4 (Compress)\n", SRA_OUTPUT_FORMAT_BAM);
    printf("\n");
}

void ParameterVerification() {
    // Internal Checking
    if (InputSGASeedOccLimitation==-1 || InputSGASeedOccLimitation>MIC_SRA_OUTPUT_MAX_ALIGNMENT) {
        printf("Potential Configuration Error:\n");
        printf("Seed alignment expects boundary %u which cannot be satisfied \n",InputSGASeedOccLimitation);
        printf("by MIC-SRA alignment model (capped at %u). \n",MIC_SRA_OUTPUT_MAX_ALIGNMENT);
        printf("Seed alignment with more than %u results will be passed to CPU\n",InputSGASeedOccLimitation);
        printf("for Seed alignment which might degrade performance\n");
    }
    
    // Internal Checking
    if (MIC_DP_OUTPUT_MAX_ALIGNMENT>MIC_PE_MAX_RESULT) {
        printf("Potential Configuration Error:\n");
        printf("DP Enhancement could reporting occurrences through PE output buffer.\n");
        printf("The number of output per read is configured as:\n");
        printf("SRA-PE <= %u and DP <= %u, which could result in overflow.\n",MIC_PE_MAX_RESULT,MIC_DP_OUTPUT_MAX_ALIGNMENT);
    }
    
    if ( InputMaxNumOfAlignment < -1 || InputMaxNumOfAlignment == 0 ) {
        printf("MaxNumOfAlignment can only be set to -1; or any positive integer.\n");
        exit(1);
    }
    
    if ( InputMaxNumOfAlignment > MIC_PE_MAX_RESULT ) {
        printf("Potential Configuration Error:\n");
        printf("MaxNumOfAlignment(%d) exceeded the capacity of MIC allocated buffer(%d)\n",InputMaxNumOfAlignment,MIC_PE_MAX_RESULT);
        printf("Any read pairs have number of output exceeding the MIC capacity are handled by CPU.\n");
        printf("That could potentially affect alignment performance.\n");
    }

    if (InputReportType!=MIC_MAIN_NON_INPUT_VALUE && 
        InputReportType!=PE_REPORT_ALL &&
        InputReportType!=PE_REPORT_ALL_BEST ) {
        fprintf(stderr, "Report mode must be %d=All-Valid %d=All-Best.\n", PE_REPORT_ALL, PE_REPORT_ALL_BEST);
        exit(1);
    }
    
    /* ATTENTION : AWAITING IMPLEMENTATION OF ALL-BEST ALIGNMENT
    if (InputReportType!=MIC_MAIN_NON_INPUT_VALUE && 
        InputReportType!=PE_REPORT_ALL &&
        InputReportType!=REPORT_ALL_BEST ) {
        
        fprintf(stderr, "Report mode must be either %d=All-Valid, %d=All-Best.\n", PE_REPORT_ALL, PE_REPORT_ALL_BEST);
        exit(1);
    }
    */

    if (InputErrorType == SRA_TYPE_MISMATCH_ONLY) {
        if (InputMaxError!=MIC_MAIN_NON_INPUT_VALUE && (InputMaxError<0 || InputMaxError>MAX_NUM_OF_MISMATCH)) {
            fprintf(stderr, "Maximum mismatches must be within the range [0-%d].\n",MAX_NUM_OF_MISMATCH);
            exit(1);
        }
    } else if (InputErrorType == SRA_TYPE_EDIT_DISTANCE) {
        if (InputMaxError!=MIC_MAIN_NON_INPUT_VALUE && (InputMaxError<0 || InputMaxError>MAX_NUM_OF_INDEL)) {
            fprintf(stderr, "Maximum edit distance must be within the range [0-%d].\n",MAX_NUM_OF_INDEL);
            exit(1);
        }
    }

    /* ATTENTION : AWAITING IMPLEMENTATION OF NBM
    if (InputMaxNBMismatch!=MIC_MAIN_NON_INPUT_VALUE && (InputMaxNBMismatch < 0 || InputMaxNBMismatch > MAX_NUM_OF_NBM_ERROR)) {
        fprintf(stderr, "Maximum Non-Branching error must be within the range [0-%d].\n",MAX_NUM_OF_NBM_ERROR);
        exit(1);
    }
    */
    
    if (InputOutputFormat != SRA_OUTPUT_FORMAT_SAM && InputOutputFormat != SRA_OUTPUT_FORMAT_BAM) {
        fprintf(stderr, "Output format must be %d=SAM or %d=BAM.\n", SRA_OUTPUT_FORMAT_SAM, SRA_OUTPUT_FORMAT_BAM);
        exit(1);
    }

    if (InputOutputFormat==SRA_OUTPUT_FORMAT_BAM && InputOutputNumSamThreads < 1) {
        fprintf(stderr, "Number of output compression threads must be at least 1\n");
        exit(1);
    }
    
    if ( InputFileArgument == MICA_INPUT_FILE_ARG_FILE) {
        if (InputPEInsertionUpperBound==MIC_MAIN_NON_INPUT_VALUE || InputPEInsertionLowerBound==MIC_MAIN_NON_INPUT_VALUE) {
            fprintf(stderr, "Lower and upper bound of the insertion must be supplied!\n");
            exit(1);
        } else if (InputPEInsertionUpperBound < InputPEInsertionLowerBound) {
            fprintf(stderr, "Supplied upper bound is smaller than lower bound!\n");
            exit(1);
        }
    }
    
    // Dynamics Program Parameters
    if ( !InputSGASeedEnhancement && InputSGASeedEnhancement_1 ) {
        fprintf(stderr, "Second level SRA Seed Enhancement is enabled but the first level is disabled.\n");
        exit(1);
    }
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

