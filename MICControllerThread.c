//
//    MICControllerThread.c
//
//    MICA / 2BWT
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

#include "MICControllerThread.h"

void MICCTWriteStatus(MICControllerThread * thread, char * status) {
    strcpy(thread->healthString,status);
}

void MICCTSetSRA(MICControllerThread * thread, int alignmentModel, SRAArguments * cpuSraArgTemplate) {
    thread->alignmentModel = alignmentModel;
    thread->cpuSraArgTemplate = cpuSraArgTemplate;
}


void MICCTSetPE(MICControllerThread * thread, PEArguments * cpuPeArgTemplate) {
    thread->cpuPeArgTemplate = cpuPeArgTemplate;
}

void MICCTSetReadThread(MICControllerThread * thread, PEReadThread * peReadThread) {
    thread->peReadThread = peReadThread;
}

void MICCTSetWriteThread(MICControllerThread * thread, PEWriteThread * peWriteThread) {
    thread->writeThread = peWriteThread;
}

MICControllerThread * MICCTCreate(uint8_t gThreadId, Idx2BWT * micIdx2BWT, char * charMap,
    unsigned int micMaxBatchSize, uint8_t threadId, int numMicThreads) {

    MICControllerThread * ret = (MICControllerThread *) MEMManMalloc(sizeof(MICControllerThread), MEMORY_TYPE_CPU);
    
    // Initialize some fields in MICControllerThread
    ret->threadId = threadId;
    ret->gThreadId = gThreadId;
    ret->micIdx2BWT = micIdx2BWT;
    ret->charMap = charMap;
    ret->micMaxBatchSize = micMaxBatchSize;
    ret->numMicThreads = numMicThreads;
    ret->t = MEMManMalloc(sizeof(pthread_t), MEMORY_TYPE_CPU);
    ret->writeThread = NULL;
    
    ret->statsTotalAlignmentTime = 0;
    ret->statsTotalReadLoadTime = 0;
    ret->statsTotalModelBuildTime = 0;
    ret->statsTotalWaitTime = 0;
    ret->statsQueryHandled = 0;
    
    #ifdef MICA_PE_ENABLE_HEALTH_CHECK
        strcpy(ret->healthString,"INIT");
    #else
        strcpy(ret->healthString,"DISABLED");
    #endif
    MICCTPollHealthString(ret);
    
    return ret;
}

// GetBatchSizing function takes the number of reads and
// calculates required number of cores; and the number of reads assigned to each core.
// The output number (batchSize * numThreads) can be greater than (totalNumRead).
static void getBatchSizing(unsigned int totalNumRead,
        unsigned int maxNumThreads, unsigned int minBatchSize,
        unsigned int * batchSize, unsigned int * numThreads) {

    // If totalNumRead is empty, the output should all be zero
    // i.e. 0 number of CPU and 0 read in each core
    if (totalNumRead == 0) {
        (*batchSize) = 0;
        (*numThreads) = 0;
        return;
    }

    // If totalNumRead is smaller than the number of threads available,
    // each core gets at most 1 read.
    if (totalNumRead <= maxNumThreads) {
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

static void buildSRAModel (SRAModelSet * sraModelSet, int readLength, 
        CPTSRAModel * cpPModels, CPTSRAModel * cpNModels) {

    #ifndef MICCT_DISABLE_PROCESSING
    
    // Determine if the SRA model exists
    SRAModel * sraPModels = SRAModelSetGetModel(sraModelSet,readLength,QUERY_POS_STRAND);
    
    if ( sraPModels == NULL ) {
    
        SRAModelConstruct(sraModelSet,readLength);
        
        SRAModel * sraPModels = SRAModelSetGetModel(sraModelSet,readLength,QUERY_POS_STRAND);
        MICSRA2BWTModelPopulate(&(cpPModels[readLength]), sraPModels);
        
        SRAModel * sraNModels = SRAModelSetGetModel(sraModelSet,readLength,QUERY_NEG_STRAND);
        MICSRA2BWTModelPopulate(&(cpNModels[readLength]), sraNModels);
    }
    #endif
}

static int getLocalConfigNumThreads (int threadId) {
    int numThreads = 0;
    #pragma offload target(mic : threadId)
    #pragma omp parallel
    #pragma omp master
        numThreads = omp_get_num_threads();

    return numThreads;
}

static unsigned copyReadBatch(PEReadBatch * peReadBatch, char * readNames, char * mateNames,
        unsigned char * readBodies, unsigned char * mateBodies,
        char * readQuality, char * mateQuality,
        uint16_t * readLengths, uint16_t * mateLengths, 
        int * readUncertaintyNumber, int * mateUncertaintyNumber,
        PEInputParam ** peInputParamBufferFrame,
        PEReadOutput ** peReadOutputBufferFrame, int * algnFuncReasonFlag, unsigned int readMateLengthsSize) {

    unsigned readInBatch = peReadBatch->readSize;
    unsigned mateInBatch = peReadBatch->mateSize;
    uint32_t readNameCopySize = MAX_SEQ_NAME_LENGTH * readInBatch;
    uint32_t mateNameCopySize = MAX_SEQ_NAME_LENGTH * mateInBatch;
    memcpy(readNames, peReadBatch->readNames, readNameCopySize);
    memcpy(mateNames, peReadBatch->mateNames, mateNameCopySize);

    uint32_t readBodyCopySize = SRA_MAX_READ_LENGTH * readInBatch; 
    uint32_t mateBodyCopySize = SRA_MAX_READ_LENGTH * mateInBatch; 
    memcpy(readBodies, peReadBatch->readPatterns, readBodyCopySize);
    memcpy(mateBodies, peReadBatch->matePatterns, mateBodyCopySize);
    memcpy(readQuality, peReadBatch->readQuality, readBodyCopySize);
    memcpy(mateQuality, peReadBatch->mateQuality, mateBodyCopySize);

    uint32_t readLengthCopySize = sizeof(uint16_t) * readInBatch;
    uint32_t mateLengthCopySize = sizeof(uint16_t) * mateInBatch;
    memcpy(readLengths, peReadBatch->readLengths, readLengthCopySize);
    memset(readLengths + readInBatch,
        0, sizeof(uint16_t) * readMateLengthsSize - readLengthCopySize);
    memcpy(mateLengths, peReadBatch->mateLengths, mateLengthCopySize);
    memset(mateLengths + mateInBatch,
        0, sizeof(uint16_t) * readMateLengthsSize - mateLengthCopySize);

    uint32_t readUncertaintyNumberCopySize = sizeof(int) * readInBatch;
    uint32_t mateUncertaintyNumberCopySize = sizeof(int) * mateInBatch;
    memcpy(readUncertaintyNumber, peReadBatch->readUncertaintyNumber, readUncertaintyNumberCopySize);
    memcpy(mateUncertaintyNumber, peReadBatch->mateUncertaintyNumber, mateUncertaintyNumberCopySize);

    (*peInputParamBufferFrame) = peReadBatch->batchParam;
    (*peReadOutputBufferFrame) = peReadBatch->batchOutput;
    (*algnFuncReasonFlag) = peReadBatch->algnFuncReasonFlag;
    
    unsigned queryInBatch = readInBatch > mateInBatch ? mateInBatch : readInBatch;
    if (readInBatch != mateInBatch) {
        printf("MICC - [ERROR] Number of reads are different between READ(%u) and MATE(%u). Only processing %u..\n", readInBatch, mateInBatch, queryInBatch);
    }
    return queryInBatch;
}

static void * start(void * micControllerThread) {
    MICControllerThread * thread = (MICControllerThread *) micControllerThread;

    uint8_t threadId = thread->threadId;
    Idx2BWT * micIdx2BWT = thread->micIdx2BWT;
    BWT * bwt = micIdx2BWT->bwt;
    BWT * rev_bwt = micIdx2BWT->rev_bwt;
    HSP * hsp = micIdx2BWT->hsp;
    LT * lt = thread->cpuSraArgTemplate->AlgnmtIndex->lookupTable;
    LT * rlt = thread->cpuSraArgTemplate->AlgnmtIndex->rev_lookupTable;
    PEReadThread * peReadThread = thread->peReadThread;

    // Timestamp'ed -----------
    double startTime = setStartTime();
    double lastEventTime = 0;
    // ------------------------

    int numOfMICThreads;
    if (thread->numMicThreads != -1) {
        numOfMICThreads = thread->numMicThreads;
    } else {
        numOfMICThreads = getLocalConfigNumThreads(threadId);
    }

    ////////////////////////////////////////////////////////////
    // CPU Allocation of Memory Frame
    ////////////////////////////////////////////////////////////
    unsigned int numReadPerCore = (thread->micMaxBatchSize + numOfMICThreads - 1) / numOfMICThreads;
    unsigned int memsize_readNameBufferFrame    = MAX_SEQ_NAME_LENGTH*numReadPerCore*numOfMICThreads;
    unsigned int memsize_readBodyBufferFrame    = SRA_MAX_READ_LENGTH*numReadPerCore*numOfMICThreads;
    unsigned int memsize_outputBufferFrame      = MIC_SRA_OUTPUT_SIZE_PER_THREAD*numOfMICThreads;
    unsigned int memsize_outputBufferMeta       = MIC_SRA_OUTPUT_META_PER_THREAD*numOfMICThreads;
    unsigned int memsize_sraOutputBufferFrame   = MIC_SRA_OUTPUT_MAX_ALIGNMENT*numOfMICThreads;
    unsigned int memsize_sraOutputBufferMeta    = MIC_SRA_OUTPUT_MAX_META*numOfMICThreads;
    unsigned int memsize_bufferPerRead          = numReadPerCore*numOfMICThreads;
    unsigned int memsize_dpOutputBufferFrame    = MIC_DP_OUTPUT_SIZE_PER_THREAD*numOfMICThreads; 
    unsigned int memsize_dpWorkingBuffer        = numOfMICThreads;
    
    unsigned int memsize_seedPeOutputBuffer     = MIC_PE_MAX_RESULT*numOfMICThreads;
    unsigned int memsize_seedSraOutputBuffer    = MIC_SRA_OUTPUT_MAX_ALIGNMENT*numOfMICThreads;
    
    char * readNameBufferFrame                  = (char*)               MEMManMalloc(sizeof(char)*memsize_readNameBufferFrame,MEMORY_TYPE_CPU);
    unsigned char * readBodyBufferFrame         = (unsigned char*)      MEMManMalloc(sizeof(unsigned char)*memsize_readBodyBufferFrame,MEMORY_TYPE_SHARED);
    uint16_t * readLengthBufferFrame            = (uint16_t*)           MEMManMalloc(sizeof(uint16_t)*memsize_bufferPerRead,MEMORY_TYPE_SHARED);
    char * readQualityBufferFrame               = (char*)               MEMManMalloc(sizeof(char)*memsize_readBodyBufferFrame,MEMORY_TYPE_CPU);
    int * readUncertaintyNumberBufferFrame      = (int*)                MEMManMalloc(sizeof(int)*memsize_bufferPerRead,MEMORY_TYPE_SHARED);
    char * mateNameBufferFrame                  = (char*)               MEMManMalloc(sizeof(char)*memsize_readNameBufferFrame,MEMORY_TYPE_CPU);
    unsigned char * mateBodyBufferFrame         = (unsigned char*)      MEMManMalloc(sizeof(unsigned char)*memsize_readBodyBufferFrame,MEMORY_TYPE_SHARED);
    uint16_t * mateLengthBufferFrame            = (uint16_t*)           MEMManMalloc(sizeof(uint16_t)*memsize_bufferPerRead,MEMORY_TYPE_SHARED);
    char * mateQualityBufferFrame               = (char*)               MEMManMalloc(sizeof(char)*memsize_readBodyBufferFrame,MEMORY_TYPE_CPU);
    int * mateUncertaintyNumberBufferFrame      = (int*)                MEMManMalloc(sizeof(int)*memsize_bufferPerRead,MEMORY_TYPE_SHARED);
    PEInputParam * peInputParamBufferFrame;
    PEReadOutput * peReadOutputBufferFrame;
    int algnFuncReasonFlag;

    uint16_t * occCount                         = (uint16_t*)           MEMManMalloc(sizeof(uint16_t)*memsize_bufferPerRead,MEMORY_TYPE_SHARED);
    unsigned int * outputBufferFrame            = (unsigned int*)       MEMManMalloc(sizeof(unsigned int)*memsize_outputBufferFrame,MEMORY_TYPE_SHARED);
    uint8_t * outputBufferStatus                = (uint8_t*)            MEMManMalloc(sizeof(uint8_t)*memsize_bufferPerRead,MEMORY_TYPE_SHARED);
    MICSRAOccMetadata * outputBufferMeta        = (MICSRAOccMetadata*)  MEMManMalloc(sizeof(MICSRAOccMetadata)*memsize_outputBufferMeta,MEMORY_TYPE_SHARED);
    
    unsigned int * readOutputBufferFrame        = (unsigned int*)       MEMManMalloc(sizeof(unsigned int)*memsize_sraOutputBufferFrame,MEMORY_TYPE_SHARED);
    MICSRAOccMetadata * readOutputBufferMeta    = (MICSRAOccMetadata*)  MEMManMalloc(sizeof(MICSRAOccMetadata)*memsize_sraOutputBufferMeta,MEMORY_TYPE_SHARED);
    
    unsigned int * mateOutputBufferFrame        = (unsigned int*)       MEMManMalloc(sizeof(unsigned int)*memsize_sraOutputBufferFrame,MEMORY_TYPE_SHARED);
    MICSRAOccMetadata * mateOutputBufferMeta    = (MICSRAOccMetadata*)  MEMManMalloc(sizeof(MICSRAOccMetadata)*memsize_sraOutputBufferMeta,MEMORY_TYPE_SHARED);

    CPTSRAModel * cpPModels                     = (CPTSRAModel*)        MEMManMalloc(sizeof(CPTSRAModel) * SRA_MAX_READ_LENGTH, MEMORY_TYPE_SHARED);
    CPTSRAModel * cpPModels_seed                = (CPTSRAModel*)        MEMManMalloc(sizeof(CPTSRAModel) * SRA_MAX_READ_LENGTH, MEMORY_TYPE_SHARED);
    CPTSRAModel * cpPModels_ext                 = (CPTSRAModel*)        MEMManMalloc(sizeof(CPTSRAModel) * SRA_MAX_READ_LENGTH, MEMORY_TYPE_SHARED);

    CPTSRAModel * cpNModels                     = (CPTSRAModel*)        MEMManMalloc(sizeof(CPTSRAModel) * SRA_MAX_READ_LENGTH, MEMORY_TYPE_SHARED);
    CPTSRAModel * cpNModels_seed                = (CPTSRAModel*)        MEMManMalloc(sizeof(CPTSRAModel) * SRA_MAX_READ_LENGTH, MEMORY_TYPE_SHARED);
    CPTSRAModel * cpNModels_ext                 = (CPTSRAModel*)        MEMManMalloc(sizeof(CPTSRAModel) * SRA_MAX_READ_LENGTH, MEMORY_TYPE_SHARED);

    MICDPOccurrence * dpOutputBufferFrame       = (MICDPOccurrence*)    MEMManMalloc(sizeof(MICDPOccurrence) * memsize_dpOutputBufferFrame, MEMORY_TYPE_SHARED);
    DPWorkMIC * dpWorkingBuffer                 = (DPWorkMIC*)          MEMManMalloc(sizeof(DPWorkMIC) * memsize_dpWorkingBuffer, MEMORY_TYPE_SHARED);

    unsigned int * readSeedOutputBufferFrame    = (unsigned int*)       MEMManMalloc(sizeof(unsigned int)*memsize_seedSraOutputBuffer,MEMORY_TYPE_SHARED);
    unsigned int * mateSeedOutputBufferFrame    = (unsigned int*)       MEMManMalloc(sizeof(unsigned int)*memsize_seedSraOutputBuffer,MEMORY_TYPE_SHARED);
    
    MICSRAOccMetadata * readSeedOutputBufferMeta = (MICSRAOccMetadata*) MEMManMalloc(sizeof(MICSRAOccMetadata)*memsize_seedSraOutputBuffer,MEMORY_TYPE_SHARED);
    MICSRAOccMetadata * mateSeedOutputBufferMeta = (MICSRAOccMetadata*) MEMManMalloc(sizeof(MICSRAOccMetadata)*memsize_seedSraOutputBuffer,MEMORY_TYPE_SHARED);
    
    unsigned int * seedPeOutputBufferFrame      = (unsigned int*)       MEMManMalloc(sizeof(unsigned int)*memsize_seedPeOutputBuffer,MEMORY_TYPE_SHARED);
    MICSRAOccMetadata * seedPeOutputBufferMeta  = (MICSRAOccMetadata*)  MEMManMalloc(sizeof(MICSRAOccMetadata)*memsize_seedPeOutputBuffer,MEMORY_TYPE_SHARED);

    PEMappingQuality * mappingQualities         = (PEMappingQuality *)  MEMManMalloc(sizeof(PEMappingQuality)*memsize_bufferPerRead, MEMORY_TYPE_SHARED);
    int * g_log_n                               = (int *)               MEMManMalloc(sizeof(int) * 256, MEMORY_TYPE_SHARED) ;
    
    bwase_initialize(g_log_n);
    
    SRAArguments * sraArg = thread->cpuSraArgTemplate;
    PEArguments * peArg = PEARGMakeMate(thread->cpuPeArgTemplate);
    peArg->PEAlgnmtInput = PEInputMakeClone(peArg->PEAlgnmtInput);
   
    // Loosen the matching criteria
    PEDPSetting * pedpSetting = peArg->pedpSetting;
    DPScores * dpScores = peArg->dpArguments->dpScores;
    SRASetting * sraSetting = sraArg->AlgnmtSetting;
    SRASetting * seedSetting = SRASettingMakeClone(sraArg->AlgnmtSetting);
    seedSetting->OutputType = SRA_REPORT_ALL_BEST;
    if (pedpSetting->SGASeedLooseCriteria>=0) {
        seedSetting->MaxError = pedpSetting->SGASeedLooseCriteria;
    } else {
        seedSetting->MaxError += pedpSetting->SGASeedLooseCriteria;
    }
   
    // oneMoreMismatchSetting increase mismatch by 1
    SRASetting * oneMoreMismatchSetting = SRASettingMakeClone(sraArg->AlgnmtSetting);
    oneMoreMismatchSetting->OutputType = SRA_REPORT_ALL;
    oneMoreMismatchSetting->MaxError += 1;

    // Timestamp'ed -----------
    double timestamp = getElapsedTime(startTime);
    printf("MICC#%d - Elapsed time (CPU Memory Allocation) : %9.4f seconds\n\n", thread->threadId, timestamp - lastEventTime);
    lastEventTime = timestamp;
    // ------------------------
    
    #ifdef MICCT_CONFIRMATION
        char c;
        printf("CONFIRMATION - CPU Memory Allocation Completed - Hit return to continue.");fflush(stdout);
        scanf("%c",&c);
    #endif

    // ////////////////////////////////////////////////////////////
    // // MIC Memory Allocation
    ////////////////////////////////////////////////////////////
    
    // MIC Memory Control for large memory block usage.
    MICMemBlockArray * micMem = MICMemCreate(threadId,numOfMICThreads);
    
    // Forward BWT
    unsigned int micBwtSaValueSizeInWord = bwt->saValueSizeInWord;

    void * _bwtBwtCode = MICMemAddBlock(micMem,(void*)bwt->bwtCode,sizeof(unsigned int),bwt->bwtSizeInWord,"BWT");
    void * _revBwtBwtCode = MICMemAddBlock(micMem,(void*)rev_bwt->bwtCode,sizeof(unsigned int),rev_bwt->bwtSizeInWord,"R-BWT");
    void * _saValue = MICMemAddBlock(micMem,(void*)bwt->saValue,sizeof(unsigned int),micBwtSaValueSizeInWord,"SA");
    
    unsigned int bwtTextLength = bwt->textLength;
    unsigned int bwtOccSizeInWord = bwt->occSizeInWord;
    unsigned int bwtInverseSa0 = bwt->inverseSa0;

    unsigned int * bwtCumulativeFreq = bwt->cumulativeFreq;
    unsigned int * bwtOccValue = bwt->occValue;
    unsigned int * bwtOccValueMajor = bwt->occValueMajor;
    unsigned int * micOccValue = bwt->micOccValue;

    #ifndef MICCT_DISABLE_PROCESSING
    #pragma offload target(mic : threadId) \
        in(bwt:length(1)                                                alloc_if(1) free_if(0)) \
        in(bwtCumulativeFreq:length(ALPHABET_SIZE+1)                    alloc_if(1) free_if(0)) \
        in(bwtOccValue:length(bwt->occSizeInWord)                       alloc_if(1) free_if(0) align(64)) \
        in(bwtOccValueMajor:length(bwt->occMajorSizeInWord)             alloc_if(1) free_if(0) align(64)) \
        in(micOccValue:length(bwt->micOccSizeInWord)                    alloc_if(1) free_if(0) align(64))
    {
        bwt->cumulativeFreq = bwtCumulativeFreq;
        bwt->bwtCode = _bwtBwtCode;
        bwt->occValue = bwtOccValue;
        bwt->occValueMajor = bwtOccValueMajor;
        bwt->saValue = _saValue;
        bwt->micOccValue = micOccValue;
    }
    #endif
    
    // Reverse BWT
    unsigned int rbwtTextLength = rev_bwt->textLength;
    unsigned int rbwtOccSizeInWord = rev_bwt->occSizeInWord;
    unsigned int rbwtInverseSa0 = rev_bwt->inverseSa0;

    unsigned int * rbwtCumulativeFreq = rev_bwt->cumulativeFreq;
    unsigned int * rbwtOccValue = rev_bwt->occValue;
    unsigned int * rbwtOccValueMajor = rev_bwt->occValueMajor;
    unsigned int * rmicOccValue = rev_bwt->micOccValue;

    #ifndef MICCT_DISABLE_PROCESSING
    #pragma offload target(mic : threadId) \
        in(rev_bwt:length(1)                                            alloc_if(1) free_if(0)) \
        in(rbwtCumulativeFreq:length(ALPHABET_SIZE+1)                   alloc_if(1) free_if(0)) \
        in(rbwtOccValue:length(rev_bwt->occSizeInWord)                  alloc_if(1) free_if(0) align(64)) \
        in(rbwtOccValueMajor:length(rev_bwt->occMajorSizeInWord)        alloc_if(1) free_if(0) align(64)) \
        in(rmicOccValue:length(rev_bwt->micOccSizeInWord)               alloc_if(1) free_if(0) align(64))
    {
        rev_bwt->cumulativeFreq = rbwtCumulativeFreq;
        rev_bwt->bwtCode = _revBwtBwtCode;
        rev_bwt->occValue = rbwtOccValue;
        rev_bwt->occValueMajor = rbwtOccValueMajor;
        rev_bwt->micOccValue = rmicOccValue;
    
    }
    #endif
    
    // HSP
    unsigned int hspLengthInFile = (hsp->dnaLength + CHAR_PER_WORD - 1) / 16;
    void * _packedDNA = MICMemAddBlock(micMem,(void*)hsp->packedDNA,sizeof(unsigned int),hspLengthInFile,"PACKDNA");
    
    unsigned int hspDnaLength = hsp->dnaLength;
    
    #ifndef MICCT_DISABLE_PROCESSING
    #pragma offload target(mic : threadId) \
        in(hsp:length(1)                                                alloc_if(1) free_if(0))
    {
        hsp->packedDNA = _packedDNA;
        hsp->dnaLength = hspDnaLength;
    
    }
    #endif

    // LT
    unsigned int ltSizeInWord = lt->ltSizeInWord;
    void * _ltTable  = MICMemAddBlock(micMem,(void*)lt->table,sizeof(unsigned int),ltSizeInWord,"LKTB");
    void * _rltTable = MICMemAddBlock(micMem,(void*)rlt->table,sizeof(unsigned int),ltSizeInWord,"R-LKTB");
    
    unsigned int ltSize = lt->tableSize;
    #ifndef MICCT_DISABLE_PROCESSING
    #pragma offload target(mic : threadId) \
        in(lt:length(1)                                                 alloc_if(1) free_if(0)) \
        in(rlt:length(1)                                                alloc_if(1) free_if(0))
    {
        lt->tableSize = ltSize;
        lt->ltSizeInWord = ltSizeInWord;
        lt->table = _ltTable;
        
        rlt->tableSize = ltSize;
        rlt->ltSizeInWord = ltSizeInWord;
        rlt->table = _rltTable;
    }
    #endif
    
    // Miscellaneous
    PEInput * peAlgnmtInput = peArg->PEAlgnmtInput;
    #ifndef MICCT_DISABLE_PROCESSING
    #pragma offload_transfer target(mic : threadId) \
        in(sraSetting:length(1)                         alloc_if(1) free_if(0)) \
        in(peAlgnmtInput:length(1)                      alloc_if(1) free_if(0)) \
        in(pedpSetting:length(1)                        alloc_if(1) free_if(0)) \
        in(dpScores:length(1)                           alloc_if(1) free_if(0)) \
        in(g_log_n: length(256)                         alloc_if(1) free_if(0))
    #endif


    // Timestamp'ed -----------
    timestamp = getElapsedTime(startTime);
    printf("MICC#%d - Elapsed time (MIC Index Copy) : %9.4f seconds\n\n", thread->threadId, timestamp - lastEventTime);
    lastEventTime = timestamp;

    #ifdef MICCT_CONFIRMATION
        printf("CONFIRMATION - 2BWT Index MIC Memory Allocation Completed - Hit return to continue.");fflush(stdout);
        scanf("%c",&c);
    #endif
    
    ////////////////////////////////////////////////////////////
    // Input and Output data Allocation
    ////////////////////////////////////////////////////////////
    #ifndef MICCT_DISABLE_PROCESSING
    #pragma offload_transfer target(mic : threadId) \
        nocopy(readBodyBufferFrame:length(memsize_readBodyBufferFrame)          alloc_if(1) free_if(0)) \
        nocopy(readLengthBufferFrame:length(memsize_bufferPerRead)              alloc_if(1) free_if(0)) \
        nocopy(readUncertaintyNumberBufferFrame:length(memsize_bufferPerRead)   alloc_if(1) free_if(0)) \
        nocopy(mateBodyBufferFrame:length(memsize_readBodyBufferFrame)          alloc_if(1) free_if(0)) \
        nocopy(mateLengthBufferFrame:length(memsize_bufferPerRead)              alloc_if(1) free_if(0)) \
        nocopy(mateUncertaintyNumberBufferFrame:length(memsize_bufferPerRead)   alloc_if(1) free_if(0)) \
        \
        nocopy(cpPModels:length(SRA_MAX_READ_LENGTH)                            alloc_if(1) free_if(0)) \
        nocopy(cpPModels_seed:length(SRA_MAX_READ_LENGTH)                       alloc_if(1) free_if(0)) \
        nocopy(cpPModels_ext:length(SRA_MAX_READ_LENGTH)                        alloc_if(1) free_if(0)) \
        nocopy(cpNModels:length(SRA_MAX_READ_LENGTH)                            alloc_if(1) free_if(0)) \
        nocopy(cpNModels_seed:length(SRA_MAX_READ_LENGTH)                       alloc_if(1) free_if(0)) \
        nocopy(cpNModels_ext:length(SRA_MAX_READ_LENGTH)                        alloc_if(1) free_if(0)) \
        \
        nocopy(occCount:length(memsize_bufferPerRead)                           alloc_if(1) free_if(0)) \
        \
        nocopy(outputBufferFrame:length(memsize_outputBufferFrame)              alloc_if(1) free_if(0)) \
        nocopy(outputBufferMeta:length(memsize_outputBufferMeta)                alloc_if(1) free_if(0)) \
        nocopy(outputBufferStatus:length(memsize_bufferPerRead)                 alloc_if(1) free_if(0)) \
        \
        nocopy(readOutputBufferFrame:length(memsize_sraOutputBufferFrame)       alloc_if(1) free_if(0)) \
        nocopy(readOutputBufferMeta:length(memsize_sraOutputBufferMeta)         alloc_if(1) free_if(0)) \
        nocopy(mateOutputBufferFrame:length(memsize_sraOutputBufferFrame)       alloc_if(1) free_if(0)) \
        nocopy(mateOutputBufferMeta:length(memsize_sraOutputBufferMeta)         alloc_if(1) free_if(0)) \
        \
        nocopy(dpOutputBufferFrame:length(memsize_dpOutputBufferFrame)          alloc_if(1) free_if(0)) \
        nocopy(dpWorkingBuffer:length(memsize_dpWorkingBuffer)                  alloc_if(1) free_if(0)) \
        \
        nocopy(readSeedOutputBufferFrame:length(memsize_seedSraOutputBuffer)    alloc_if(1) free_if(0)) \
        nocopy(mateSeedOutputBufferFrame:length(memsize_seedSraOutputBuffer)    alloc_if(1) free_if(0)) \
        nocopy(readSeedOutputBufferMeta:length(memsize_seedSraOutputBuffer)     alloc_if(1) free_if(0)) \
        nocopy(mateSeedOutputBufferMeta:length(memsize_seedSraOutputBuffer)     alloc_if(1) free_if(0)) \
        \
        nocopy(seedPeOutputBufferFrame:length(memsize_seedPeOutputBuffer)       alloc_if(1) free_if(0)) \
        nocopy(seedPeOutputBufferMeta:length(memsize_seedPeOutputBuffer)        alloc_if(1) free_if(0)) \
        \
        nocopy(mappingQualities: length(memsize_bufferPerRead)                  alloc_if(1) free_if(0))
    #endif

    #ifdef MICCT_CONFIRMATION
        printf("CONFIRMATION - In/Output Buffer Allocation Completed - Hit return to continue.");fflush(stdout);
        scanf("%c",&c);
    #endif

    // Timestamp'ed -----------
    timestamp = getElapsedTime(startTime);
    printf("MICC#%d - Elapsed time (MIC Input Output Allocation) : %9.4f seconds\n\n", thread->threadId, timestamp - lastEventTime);
    lastEventTime = timestamp;

    // Construct SRA Model from SRA Settings
    SRAModelSet * sraModelSet = SRAModelSetConstruct(
                                    sraArg->AlgnmtSetting,
                                    sraArg->AlgnmtIndex,
                                    SRA_MODEL_16G,
                                    SRA_MIN_READ_LENGTH,
                                    SRA_MAX_READ_LENGTH);
    SRAModelSet * sraModelSet_Seed = SRAModelSetConstruct(
                                    seedSetting,
                                    sraArg->AlgnmtIndex,
                                    SRA_MODEL_16G,
                                    SRA_MIN_READ_LENGTH,
                                    SRA_MAX_READ_LENGTH);
    SRAModelSet * sraModelSet_Extend = SRAModelSetConstruct(
                                    oneMoreMismatchSetting,
                                    sraArg->AlgnmtIndex,
                                    SRA_MODEL_16G,
                                    SRA_MIN_READ_LENGTH,
                                    SRA_MAX_READ_LENGTH);
    
    
    // Write Thread
    PEWriteThread * writeThread = thread->writeThread;
    
    double processingTime = 0;
    double alignmentTime = 0;
    double readLoadTime = 0;
    double modelTime = 0;
    double waitTime = 0;

    unsigned int queryInBatch = 1;
    unsigned long long queryHandled = 0;

    #ifndef MICCT_DISABLE_PROCESSING
    while (queryInBatch!=0) {

        #ifdef MICA_PE_ENABLE_HEALTH_CHECK
            // Update status to WAIT
            MICCTWriteStatus(thread,"WAIT");
        #endif

        // Wait for the read batch
        PEReadBatch * peReadBatch = PEReadThread_wait(peReadThread);

        #ifdef MICA_PE_ENABLE_HEALTH_CHECK
            // Update status to WAIT
            MICCTWriteStatus(thread,"START");
        #endif
        
        // Timestamp'ed -----------
        timestamp = getElapsedTime(startTime);
        //printf("MICC#%d - Elapsed time (Read Loading) : %9.4f seconds\n\n",thread->threadId , timestamp - lastEventTime);
        waitTime += timestamp - lastEventTime;
        processingTime += timestamp - lastEventTime;
        lastEventTime = timestamp;
        // ------------------------
        
        //////////////////////////////////////
        // Copy the read batch
        //////////////////////////////////////
        queryInBatch = copyReadBatch(peReadBatch, readNameBufferFrame, mateNameBufferFrame,
                readBodyBufferFrame, mateBodyBufferFrame, readQualityBufferFrame, mateQualityBufferFrame,
                readLengthBufferFrame, mateLengthBufferFrame, readUncertaintyNumberBufferFrame, mateUncertaintyNumberBufferFrame,
                &peInputParamBufferFrame, &peReadOutputBufferFrame, &algnFuncReasonFlag,
                memsize_bufferPerRead);

        PEReadThread_signal(peReadThread);
        if (queryInBatch == 0) {
        
            #ifdef MICA_PE_ENABLE_HEALTH_CHECK
                // Update status to WAIT
                MICCTWriteStatus(thread,"EXITING");
            #endif
            
            break;
        }
        
        // Timestamp'ed -----------
        timestamp = getElapsedTime(startTime);
        //printf("MICC#%d - Elapsed time (Read Loading) : %9.4f seconds\n\n",thread->threadId , timestamp - lastEventTime);
        readLoadTime += timestamp - lastEventTime;
        processingTime += timestamp - lastEventTime;
        lastEventTime = timestamp;
        // ------------------------
                
        // Build up the SRA Model in case there are new read lengths
        unsigned int i;
        for (i=0;i<queryInBatch;i++) {
            // buildSRAModels with each read and mate into
            // sraModel, seedModel and oneMoreMismatchModel
            int readLength = readLengthBufferFrame[i];
            int mateLength = mateLengthBufferFrame[i];
            
            // If seedLength > 1, treat input as constant seedLength or seedOverlap
            // If seedLength <= 1, treat input as a ratio to readLen
            // + 0.5 to do rounding, set Length to zero when readLength or mateLength = 0
            int readSeedLength = ( pedpSetting->SGASeedLength > 1.0 && readLength != 0 ? pedpSetting->SGASeedLength
            : (double) readLength * pedpSetting->SGASeedLength ) + 0.5;
            int mateSeedLength = ( pedpSetting->SGASeedLength > 1.0 && mateLength != 0 ? pedpSetting->SGASeedLength
            : (double) mateLength * pedpSetting->SGASeedLength ) + 0.5;
            int readSeedLength2 = ( pedpSetting->SGASeedLength_1 > 1.0 && readLength != 0 ? pedpSetting->SGASeedLength_1
            : (double) readLength * pedpSetting->SGASeedLength_1 ) + 0.5;
            int mateSeedLength2 = ( pedpSetting->SGASeedLength_1 > 1.0 && mateLength != 0 ? pedpSetting->SGASeedLength_1
            : (double) mateLength * pedpSetting->SGASeedLength_1 ) + 0.5;
            buildSRAModel(sraModelSet, readLength, cpPModels, cpNModels);
            buildSRAModel(sraModelSet, mateLength, cpPModels, cpNModels);
            buildSRAModel(sraModelSet_Seed, readSeedLength, cpPModels_seed, cpNModels_seed);
            buildSRAModel(sraModelSet_Seed, mateSeedLength, cpPModels_seed, cpNModels_seed);
            buildSRAModel(sraModelSet_Seed, readSeedLength2, cpPModels_seed, cpNModels_seed);
            buildSRAModel(sraModelSet_Seed, mateSeedLength2, cpPModels_seed, cpNModels_seed);
            buildSRAModel(sraModelSet_Extend, readLength, cpPModels_ext, cpNModels_ext);
            buildSRAModel(sraModelSet_Extend, mateLength, cpPModels_ext, cpNModels_ext);

        }

        // Timestamp'ed -----------
        timestamp = getElapsedTime(startTime);
        //printf("MICC#%d - Elapsed time (Read Loading) : %9.4f seconds\n\n", timestamp - lastEventTime);
        modelTime += timestamp - lastEventTime;
        processingTime += timestamp - lastEventTime;
        lastEventTime = timestamp;
        // ------------------------

        unsigned int batchSize, numThreads;
        getBatchSizing(queryInBatch,numOfMICThreads,1,&batchSize,&numThreads);
        
        memset(occCount,0,sizeof(uint16_t)*memsize_bufferPerRead);
        
        memset(outputBufferStatus,MIC_OUTPUT_STATUS_OPEN,sizeof(unsigned char)*memsize_bufferPerRead);
        
        // Populate input parameter for batch
        peAlgnmtInput->insertUbound = peInputParamBufferFrame->InputPEInsertionUpperBound;
        peAlgnmtInput->insertLbound = peInputParamBufferFrame->InputPEInsertionLowerBound;

        #ifdef MICA_PE_ENABLE_HEALTH_CHECK
            // Update status to WAIT
            MICCTWriteStatus(thread,"MIC-ALGNMT");
        #endif
        
        ////////////////////////////////////////////////////////////
        // Alignment Body
        ////////////////////////////////////////////////////////////
        #ifdef MIC_PRINT_INPUT_PARAM
            printf("MICC#%d - MIC Core Input Parameters - Insertion %u-%u\n", 
                thread->threadId,
                peAlgnmtInput->insertLbound,
                peAlgnmtInput->insertUbound);
            fflush(stdout);
        #endif
        printf("MICC#%d - MIC Core initiated (%u)\n", thread->threadId,numThreads); fflush(stdout);
        unsigned long long j,k;
        #pragma offload target(mic : threadId) \
            nocopy(bwt:length(1)                                                    alloc_if(0) free_if(0)) \
            nocopy(rev_bwt:length(1)                                                alloc_if(0) free_if(0)) \
            nocopy(hsp:length(1)                                                    alloc_if(0) free_if(0)) \
            nocopy(lt:length(1)                                                     alloc_if(0) free_if(0)) \
            nocopy(rlt:length(1)                                                    alloc_if(0) free_if(0)) \
            \
            in(readBodyBufferFrame: length(memsize_readBodyBufferFrame)             alloc_if(0) free_if(0)) \
            in(readLengthBufferFrame: length(memsize_bufferPerRead)                 alloc_if(0) free_if(0)) \
            in(readUncertaintyNumberBufferFrame: length(memsize_bufferPerRead)      alloc_if(0) free_if(0)) \
            in(mateBodyBufferFrame: length(memsize_readBodyBufferFrame)             alloc_if(0) free_if(0)) \
            in(mateLengthBufferFrame: length(memsize_bufferPerRead)                 alloc_if(0) free_if(0)) \
            in(mateUncertaintyNumberBufferFrame: length(memsize_bufferPerRead)      alloc_if(0) free_if(0)) \
            \
            in(cpPModels:length(SRA_MAX_READ_LENGTH)                                alloc_if(0) free_if(0)) \
            in(cpPModels_seed:length(SRA_MAX_READ_LENGTH)                           alloc_if(0) free_if(0)) \
            in(cpPModels_ext:length(SRA_MAX_READ_LENGTH)                            alloc_if(0) free_if(0)) \
            in(cpNModels:length(SRA_MAX_READ_LENGTH)                                alloc_if(0) free_if(0)) \
            in(cpNModels_seed:length(SRA_MAX_READ_LENGTH)                           alloc_if(0) free_if(0)) \
            in(cpNModels_ext:length(SRA_MAX_READ_LENGTH)                            alloc_if(0) free_if(0)) \
            \
            inout(occCount:length(memsize_bufferPerRead)                            alloc_if(0) free_if(0)) \
            \
            out(outputBufferFrame:length(memsize_outputBufferFrame)                 alloc_if(0) free_if(0)) \
            out(outputBufferMeta:length(memsize_outputBufferMeta)                   alloc_if(0) free_if(0)) \
            inout(outputBufferStatus:length(memsize_bufferPerRead)                  alloc_if(0) free_if(0)) \
            \
            nocopy(readOutputBufferFrame:length(memsize_sraOutputBufferFrame)       alloc_if(0) free_if(0)) \
            nocopy(readOutputBufferMeta:length(memsize_sraOutputBufferMeta)         alloc_if(0) free_if(0)) \
            nocopy(mateOutputBufferFrame:length(memsize_sraOutputBufferFrame)       alloc_if(0) free_if(0)) \
            nocopy(mateOutputBufferMeta:length(memsize_sraOutputBufferMeta)         alloc_if(0) free_if(0)) \
            \
            out(dpOutputBufferFrame:length(memsize_dpOutputBufferFrame)             alloc_if(0) free_if(0)) \
            nocopy(dpWorkingBuffer:length(memsize_dpWorkingBuffer)                  alloc_if(0) free_if(0)) \
            \
            in(sraSetting:length(1)                                                 alloc_if(0) free_if(0)) \
            in(peAlgnmtInput:length(1)                                              alloc_if(0) free_if(0)) \
            nocopy(pedpSetting:length(1)                                            alloc_if(0) free_if(0)) \
            nocopy(dpScores:length(1)                                               alloc_if(0) free_if(0)) \
            \
            nocopy(readSeedOutputBufferFrame:length(memsize_seedSraOutputBuffer)    alloc_if(0) free_if(0)) \
            nocopy(mateSeedOutputBufferFrame:length(memsize_seedSraOutputBuffer)    alloc_if(0) free_if(0)) \
            nocopy(readSeedOutputBufferMeta:length(memsize_seedSraOutputBuffer)     alloc_if(0) free_if(0)) \
            nocopy(mateSeedOutputBufferMeta:length(memsize_seedSraOutputBuffer)     alloc_if(0) free_if(0)) \
            \
            nocopy(seedPeOutputBufferFrame:length(memsize_seedPeOutputBuffer)       alloc_if(0) free_if(0)) \
            nocopy(seedPeOutputBufferMeta:length(memsize_seedPeOutputBuffer)        alloc_if(0) free_if(0)) \
            \
            out(mappingQualities: length(memsize_bufferPerRead)                     alloc_if(0) free_if(0)) \
            nocopy(g_log_n: length(256)                                             alloc_if(0) free_if(0))
        
        #pragma omp parallel for private(i,j,k) num_threads(numThreads)
        for (i = 0; i < numThreads; i++) {
            unsigned int offset = i * batchSize;
            int charIdx;
            unsigned char c;
            uint16_t readMetaCount,mateMetaCount;
            unsigned char readComplement[SRA_MAX_READ_LENGTH];
            unsigned char mateComplement[SRA_MAX_READ_LENGTH];
            
            // Shared in core
            unsigned int * outputPtr     = &outputBufferFrame[i*MIC_SRA_OUTPUT_SIZE_PER_THREAD];
            unsigned int * readOutputPtr = &readOutputBufferFrame[i*MIC_SRA_OUTPUT_MAX_ALIGNMENT];
            unsigned int * mateOutputPtr = &mateOutputBufferFrame[i*MIC_SRA_OUTPUT_MAX_ALIGNMENT];

            MICSRAOccMetadata * metaPtr      = &outputBufferMeta[i*MIC_SRA_OUTPUT_META_PER_THREAD];
            MICSRAOccMetadata * readMetaPtr  = &readOutputBufferMeta[i*MIC_SRA_OUTPUT_MAX_META];
            MICSRAOccMetadata * mateMetaPtr  = &mateOutputBufferMeta[i*MIC_SRA_OUTPUT_MAX_META];
            
            unsigned int * readSeedOutputPtr = &readSeedOutputBufferFrame[i*MIC_SRA_OUTPUT_MAX_ALIGNMENT];
            unsigned int * mateSeedOutputPtr = &mateSeedOutputBufferFrame[i*MIC_SRA_OUTPUT_MAX_ALIGNMENT];
            
            MICSRAOccMetadata * readSeedMetaPtr = &readSeedOutputBufferMeta[i*MIC_SRA_OUTPUT_MAX_ALIGNMENT];
            MICSRAOccMetadata * mateSeedMetaPtr = &mateSeedOutputBufferMeta[i*MIC_SRA_OUTPUT_MAX_ALIGNMENT];
            
            unsigned int * peSeedOutputPtr      = &seedPeOutputBufferFrame[i*MIC_PE_MAX_RESULT];
            MICSRAOccMetadata * peSeedMetaPtr   = &seedPeOutputBufferMeta[i*MIC_PE_MAX_RESULT];
            
            unsigned int peOutputVacancy     = MIC_SRA_OUTPUT_SIZE_PER_THREAD;
            unsigned int peMetaVacancy         = MIC_SRA_OUTPUT_META_PER_THREAD;

            uint16_t readOccCount            = 0;
            uint16_t mateOccCount            = 0;
            
            uint8_t readOutputStatus         = MIC_OUTPUT_STATUS_OPEN;
            uint8_t mateOutputStatus         = MIC_OUTPUT_STATUS_OPEN;
            
            // Multiple cell per read
            unsigned char * readPatternPtr   = &readBodyBufferFrame[offset*SRA_MAX_READ_LENGTH];
            unsigned char * matePatternPtr   = &mateBodyBufferFrame[offset*SRA_MAX_READ_LENGTH];

            MICSRAArguments * readMicArgs = MICSRAARGConstruct();
            MICSRAArguments * mateMicArgs = MICSRAARGConstruct();
            MICPEArguments peArgs;
            MICPEArguments peSeedArgs;
            MICPEArgumentsSetBounds(&peArgs, peAlgnmtInput->insertLbound, peAlgnmtInput->insertUbound);
            MICPEArgumentsSetBounds(&peSeedArgs, peAlgnmtInput->insertLbound, peAlgnmtInput->insertUbound);
            MICPEArgumentsSetDPScores(&peArgs, dpScores);
            MICPEArgumentsSetDPScores(&peSeedArgs, dpScores);
            MICPEArgumentsSetMaxOutput(&peArgs, peAlgnmtInput->maxResult);
            MICPEArgumentsSetMergeEnable(&peArgs, 1);
            MICPEArgumentsSetMergeEnable(&peSeedArgs, 1);

            // ATTENTION. BWT is hacked into MICSRAArgs before
            // proper implementation of MICSRAMdl.
            if (peAlgnmtInput->OutputType == PE_REPORT_ALL) {
                readMicArgs->outputType     = SRA_REPORT_ALL;
            } else {
                readMicArgs->outputType     = SRA_REPORT_ALL_SORTED;
            }
            readMicArgs->bwt                = bwt;
            readMicArgs->rev_bwt            = rev_bwt;
            readMicArgs->lt                 = lt;
            readMicArgs->rlt                = rlt;
            readMicArgs->outputStatus       = &readOutputStatus;
            readMicArgs->outputBlock        = readOutputPtr;
            readMicArgs->occCount           = &readOccCount;
            readMicArgs->metaBlock          = readMetaPtr;
            readMicArgs->metaCount          = &readMetaCount;
            readMicArgs->readCode_Complt    = readComplement;
            readMicArgs->seedOffset         = 0;
            readMicArgs->seedOffset_Complt  = 0;
            readMicArgs->maxNBMismatch      = sraSetting->MaxNBMismatch;

            if (peAlgnmtInput->OutputType == PE_REPORT_ALL) {
                mateMicArgs->outputType     = SRA_REPORT_ALL;
            } else {
                mateMicArgs->outputType     = SRA_REPORT_ALL_SORTED;
            }
            mateMicArgs->bwt                = bwt;
            mateMicArgs->rev_bwt            = rev_bwt;
            mateMicArgs->lt                 = lt;
            mateMicArgs->rlt                = rlt;
            mateMicArgs->outputStatus       = &mateOutputStatus;
            mateMicArgs->outputBlock        = mateOutputPtr;
            mateMicArgs->occCount           = &mateOccCount;
            mateMicArgs->metaBlock          = mateMetaPtr;
            mateMicArgs->metaCount          = &mateMetaCount;
            mateMicArgs->readCode_Complt    = mateComplement;
            mateMicArgs->seedOffset         = 0;
            mateMicArgs->seedOffset_Complt  = 0;
            mateMicArgs->maxNBMismatch      = sraSetting->MaxNBMismatch;
            //////////////////////////////////////////
            // Create Copy for Seed Alignment
            //////////////////////////////////////////
            MICSRAArguments * readSeedMicArgs = MICSRAARGMakeCopy(readMicArgs);
            MICSRAArguments * mateSeedMicArgs = MICSRAARGMakeCopy(mateMicArgs);
            
            uint16_t readSeedOccCount         = 0;
            uint16_t mateSeedOccCount         = 0;
            
            uint8_t readSeedOutputStatus      = MIC_OUTPUT_STATUS_OPEN;
            uint8_t mateSeedOutputStatus      = MIC_OUTPUT_STATUS_OPEN;
            
            uint16_t readSeedMetaCount        = 0;
            uint16_t mateSeedMetaCount        = 0;
            
            readSeedMicArgs->outputType       = SRA_REPORT_ALL_BEST;
            readSeedMicArgs->outputStatus     = &readSeedOutputStatus;
            readSeedMicArgs->outputBlock      = readSeedOutputPtr;
            readSeedMicArgs->occCount         = &readSeedOccCount;
            readSeedMicArgs->metaBlock        = readSeedMetaPtr;
            readSeedMicArgs->metaCount        = &readSeedMetaCount;
           
            mateSeedMicArgs->outputType       = SRA_REPORT_ALL_BEST;
            mateSeedMicArgs->outputStatus     = &mateSeedOutputStatus;
            mateSeedMicArgs->outputBlock      = mateSeedOutputPtr;
            mateSeedMicArgs->occCount         = &mateSeedOccCount;
            mateSeedMicArgs->metaBlock        = mateSeedMetaPtr;
            mateSeedMicArgs->metaCount        = &mateSeedMetaCount;
            
            // DP related stuff start 
            // TODO this is put here temporarily, should be moved to somewhere else 
            // ATTENTION
            MICDPOccurrence * dpOcc = &dpOutputBufferFrame[i*MIC_DP_OUTPUT_SIZE_PER_THREAD];
            DPWorkMIC * dpWork = &dpWorkingBuffer[i];
            DPParametersMIC dpp;
            dpp.scoreMatch = dpScores->dpMatch;
            dpp.scoreMismatch = dpScores->dpMismatch;
            dpp.scoreGapOpen = dpScores->dpGapOpen;
            dpp.scoreGapExtend = dpScores->dpGapExtend;
            dpp.maxLeftClip = pedpSetting->SGASoftHeadClipLength * 100;
            dpp.maxRightClip = pedpSetting->SGASoftTailClipLength * 100;
            dpp.maxTotalClip = pedpSetting->SGASoftTotalClipLength * 100;
            dpp.scoreThreshold = pedpSetting->SGAScoreTF;
            DPWorkInitMIC(&dpp, dpWork);
            unsigned dpOccCount = 0;
            
            // end of DP related stuff
            #pragma ivdep
            for (j=0;j<batchSize;j++) {
            
                // Only proceed if the empty slot is absolutely enough 
                // for another round of PE alignment
                if (peOutputVacancy<MIC_PE_MAX_RESULT) {
                    outputBufferStatus[offset + j] = MIC_PE_OUTPUT_STATUS_UNHANDLE;
                    continue;
                } else if (peMetaVacancy<MIC_PE_MAX_RESULT) {
                    outputBufferStatus[offset + j] = MIC_PE_OUTPUT_STATUS_UNHANDLE;
                    continue;
                }
            
                // Initialisation of the round based result
                uint16_t readLen             = readLengthBufferFrame[offset + j];
                readOutputStatus             = MIC_OUTPUT_STATUS_OPEN; 
                readOccCount                 = 0;
                readMetaCount                = 0;
                // Initialisation of the round based result
                uint16_t mateLen             = mateLengthBufferFrame[offset + j];
                mateOutputStatus             = MIC_OUTPUT_STATUS_OPEN;
                mateOccCount                 = 0;
                mateMetaCount                = 0;

                // SeedDp - Initialisation of the round based result
                uint16_t peSeedOccCount = 0;
                uint8_t  peSeedStatus = MIC_OUTPUT_STATUS_OPEN;
                
                if (readLen==0 || mateLen==0) {
                    //Skip if the read is an invalid read
                    readOutputStatus = MIC_OUTPUT_STATUS_SKIPPED;
                    mateOutputStatus = MIC_OUTPUT_STATUS_SKIPPED;
                    outputBufferStatus[offset + j] = MIC_PE_OUTPUT_STATUS_SKIPPED;
                    occCount[offset + j] = 0;
                } else {
                    readMicArgs->seedLength     = readLen;
                    readMicArgs->readLength     = readLen;
                    readMicArgs->readCode       = readPatternPtr;
                    
                    SRAWorkingMemoryInitialise(&(readMicArgs->AlgnmtMemory));

                    mateMicArgs->seedLength     = mateLen;
                    mateMicArgs->readLength     = mateLen;
                    mateMicArgs->readCode       = matePatternPtr;
                    
                    SRAWorkingMemoryInitialise(&(mateMicArgs->AlgnmtMemory));
                    
                    /////////////////////////////////////////////////
                    // ATTENTION. The really bad read flip.
                    /////////////////////////////////////////////////
                    for (k=0;k<readLen;k++) {
                        readMicArgs->readCode_Complt[readLen-k-1] = 3 - readPatternPtr[k];
                    }
                    for (k=0;k<mateLen;k++) {
                        mateMicArgs->readCode_Complt[mateLen-k-1] = 3 - matePatternPtr[k];
                    }
            
                    MICPEArgumentsConfig(&peArgs, readMicArgs, mateMicArgs, peOutputVacancy,
                        peMetaVacancy, outputPtr, metaPtr, &occCount[offset + j],
                        &outputBufferStatus[offset + j]);
                    
                    uint32_t dpOccCountInc = 0;
                    
                    MICSHProcessRead(hsp, readMicArgs, mateMicArgs,
                                     cpPModels, cpNModels, peAlgnmtInput, &peArgs,
                                     readSeedMicArgs,mateSeedMicArgs,
                                     cpPModels_seed, cpNModels_seed, &peSeedArgs,
                                     peSeedOutputPtr, peSeedMetaPtr,
                                     pedpSetting, dpWork,
                                     outputPtr, &(outputBufferStatus[offset + j]), &(occCount[offset + j]),
                                     dpOcc, dpOccCount, &dpOccCountInc);
                    
                    mappingQualities[offset + j] = MicCalculatePEMappingQuality(
                            readMicArgs, mateMicArgs, &peArgs,
                            dpOcc + dpOccCount,
                            cpPModels_ext, cpNModels_ext, g_log_n);
                            
                    dpOccCount += dpOccCountInc;

                    outputPtr += *peArgs.occCount;
                    metaPtr += *peArgs.occCount;
                    peMetaVacancy -= *peArgs.occCount;
                    peOutputVacancy -= *peArgs.occCount;

                }

                readPatternPtr += SRA_MAX_READ_LENGTH;
                matePatternPtr += SRA_MAX_READ_LENGTH;
            }

            MICSRAARGFree(readMicArgs);
            MICSRAARGFree(mateMicArgs);
            MICSRAARGFree(readSeedMicArgs);
            MICSRAARGFree(mateSeedMicArgs);
        } // End offload

        // Timestamp'ed -----------
        timestamp = getElapsedTime(startTime);
        alignmentTime += timestamp - lastEventTime;
        processingTime += timestamp - lastEventTime;
        printf("MICC#%d - Elapsed time (Offload) : %9.4f seconds ( %9.4f rps ) \n", thread->threadId, timestamp - lastEventTime, (queryHandled+queryInBatch) / alignmentTime);
        fflush(0);
        lastEventTime = timestamp;
        // ------------------------

        ////////////////////////////////////////////////////////////
        // Produce Statistics of Output
        ////////////////////////////////////////////////////////////
        #ifdef MIC_DEBUG_PRINT_OUTPUT_COUNTS
        unsigned numOfClosed = 0;
        unsigned numOfOccContent = 0;
        for (i=0;i<queryInBatch;i++) {
            if (outputBufferStatus[i] == MIC_PE_OUTPUT_STATUS_CLOSED) {
                numOfClosed++;
                printf("MICC#%d - Closed Item size = %u\n", thread->threadId, occCount[i]);
            }  
            unsigned int * outputPtr = &outputBufferFrame[i*MIC_SRA_OUTPUT_SIZE_PER_THREAD];
            numOfOccContent+=occCount[i];
        }
        printf("MICC#%d - [STATISTICS] Number of Occurrences found = %u\n", thread->threadId, numOfOccContent);
        printf("MICC#%d - [STATISTICS] Number of Closed Output Collector = %u\n", thread->threadId, numOfClosed);
        #endif
        
        ////////////////////////////////////////////////////////////
        // Handling Output
        // This following pick results from MIC and output with CPU
        ////////////////////////////////////////////////////////////
        
        #ifdef MICA_PE_ENABLE_HEALTH_CHECK
            // Update status to WAIT
            MICCTWriteStatus(thread,"CPU-ALGNMT");
        #endif
        
        #ifndef MICA_PE_DISCARD_OUTPUT_FROM_MIC
            PEWTCloneBuffer(writeThread,
                            numThreads,batchSize,queryHandled,queryInBatch,
                            readNameBufferFrame,readBodyBufferFrame,readLengthBufferFrame,
                            mateNameBufferFrame,mateBodyBufferFrame,mateLengthBufferFrame,
                            readQualityBufferFrame, mateQualityBufferFrame,
                            peInputParamBufferFrame, peReadOutputBufferFrame,
                            MICA_PE_WRITE_THREAD_PAYLOAD_OUTPUT_TYPE_MIC_MIX,
                            occCount,outputBufferStatus,outputBufferFrame,
                            outputBufferMeta,dpOutputBufferFrame,mappingQualities,
                            algnFuncReasonFlag);
        #endif

        queryHandled+=queryInBatch;

        // Timestamp'ed -----------
        timestamp = getElapsedTime(startTime);
        //printf("MICC#%d - Elapsed time (Output) : %9.4f seconds\n\n", timestamp - lastEventTime);
        waitTime += timestamp - lastEventTime;
        processingTime += timestamp - lastEventTime;
        lastEventTime = timestamp;
        // ------------------------

        printf("MICC#%d - Number of Query Handled = %llu\n\n", thread->threadId,queryHandled);
        fflush(0);
    }
    #endif
    
    #ifdef MICCT_CONFIRMATION
        printf("CONFIRMATION - Read Processing Completed - Hit return to continue.");fflush(stdout);
        scanf("%c",&c);
    #endif
    
    ////////////////////////////////////////////////////////////
    // Statistics
    ////////////////////////////////////////////////////////////
    thread->statsTotalAlignmentTime = alignmentTime;
    thread->statsTotalReadLoadTime = readLoadTime;
    thread->statsTotalModelBuildTime = modelTime;
    thread->statsTotalWaitTime = waitTime;
    thread->statsQueryHandled = queryHandled;

    thread->writeThread = NULL;
    
    SRAModelSetFree(sraModelSet);
    SRAModelSetFree(sraModelSet_Seed);
    SRAModelSetFree(sraModelSet_Extend);
    
    SRASettingFree(seedSetting);
    SRASettingFree(oneMoreMismatchSetting);
    
    PEInputCloneFree(peArg->PEAlgnmtInput);
    PEARGMateFree(peArg);
    
    ////////////////////////////////////////////////////////////
    // Free MIC Memory Allocated
    ////////////////////////////////////////////////////////////
    printf("\nMICC#%d - Free index from MIC Controller ... ", thread->threadId);
    fflush(stdout);
    #ifndef MICCT_DISABLE_PROCESSING
    #pragma offload_transfer target(mic : threadId) \
        nocopy(bwt:length(1)                                                alloc_if(0) free_if(1)) \
        nocopy(bwtCumulativeFreq:length(ALPHABET_SIZE+1)                    alloc_if(0) free_if(1)) \
        nocopy(bwtOccValue:length(bwt->occSizeInWord)                       alloc_if(0) free_if(1)) \
        nocopy(bwtOccValueMajor:length(bwt->occMajorSizeInWord)             alloc_if(0) free_if(1)) \
        nocopy(micOccValue:length(bwt->micOccSizeInWord)                    alloc_if(0) free_if(1))
    
    #pragma offload_transfer target(mic : threadId) \
        nocopy(rev_bwt:length(1)                                            alloc_if(0) free_if(1)) \
        nocopy(rbwtCumulativeFreq:length(ALPHABET_SIZE+1)                   alloc_if(0) free_if(1)) \
        nocopy(rbwtOccValue:length(rev_bwt->occSizeInWord)                  alloc_if(0) free_if(1)) \
        nocopy(rbwtOccValueMajor:length(rev_bwt->occMajorSizeInWord)        alloc_if(0) free_if(1)) \
        nocopy(rmicOccValue:length(rev_bwt->micOccSizeInWord)               alloc_if(0) free_if(1))

    // HSP
    #pragma offload_transfer target(mic : threadId) \
        nocopy(hsp:length(1)                                                alloc_if(0) free_if(1))
    
    // LT
    #pragma offload_transfer target(mic : threadId) \
        nocopy(lt:length(1)                                                 alloc_if(0) free_if(1)) \
        nocopy(rlt:length(1)                                                alloc_if(0) free_if(1))
    
    // Miscellaneous
    #pragma offload_transfer target(mic : threadId) \
        nocopy(sraSetting: length(1)                                        alloc_if(0) free_if(1)) \
        nocopy(peAlgnmtInput: length(1)                                     alloc_if(0) free_if(1)) \
        nocopy(dpScores: length(1)                                          alloc_if(0) free_if(1)) \
        nocopy(g_log_n: length(256)                                         alloc_if(0) free_if(1))

    // Input and Output Data
    #pragma offload_transfer target(mic : threadId) \
        nocopy(readBodyBufferFrame:length(memsize_readBodyBufferFrame)          alloc_if(0) free_if(1)) \
        nocopy(readLengthBufferFrame:length(memsize_bufferPerRead)              alloc_if(0) free_if(1)) \
        nocopy(mateBodyBufferFrame:length(memsize_readBodyBufferFrame)          alloc_if(0) free_if(1)) \
        nocopy(mateLengthBufferFrame:length(memsize_bufferPerRead)              alloc_if(0) free_if(1)) \
        \
        nocopy(cpPModels:length(SRA_MAX_READ_LENGTH)                            alloc_if(0) free_if(1)) \
        nocopy(cpPModels_seed:length(SRA_MAX_READ_LENGTH)                       alloc_if(0) free_if(1)) \
        nocopy(cpPModels_ext:length(SRA_MAX_READ_LENGTH)                        alloc_if(0) free_if(1)) \
        nocopy(cpNModels:length(SRA_MAX_READ_LENGTH)                            alloc_if(0) free_if(1)) \
        nocopy(cpNModels_seed:length(SRA_MAX_READ_LENGTH)                       alloc_if(0) free_if(1)) \
        nocopy(cpNModels_ext:length(SRA_MAX_READ_LENGTH)                        alloc_if(0) free_if(1)) \
        \
        nocopy(occCount:length(memsize_bufferPerRead)                           alloc_if(0) free_if(1)) \
        \
        nocopy(outputBufferFrame:length(memsize_outputBufferFrame)              alloc_if(0) free_if(1)) \
        nocopy(outputBufferStatus:length(memsize_bufferPerRead)                 alloc_if(0) free_if(1)) \
        nocopy(outputBufferMeta:length(memsize_outputBufferMeta)                alloc_if(0) free_if(1)) \
        \
        nocopy(readOutputBufferFrame:length(memsize_sraOutputBufferFrame)       alloc_if(0) free_if(1)) \
        nocopy(readOutputBufferMeta:length(memsize_sraOutputBufferMeta)         alloc_if(0) free_if(1)) \
        nocopy(mateOutputBufferFrame:length(memsize_sraOutputBufferFrame)       alloc_if(0) free_if(1)) \
        nocopy(mateOutputBufferMeta:length(memsize_sraOutputBufferMeta)         alloc_if(0) free_if(1)) \
        \
        nocopy(dpOutputBufferFrame:length(memsize_dpOutputBufferFrame)          alloc_if(0) free_if(1)) \
        nocopy(dpWorkingBuffer:length(memsize_dpWorkingBuffer)                  alloc_if(0) free_if(1)) \
        \
        nocopy(readSeedOutputBufferFrame:length(memsize_seedSraOutputBuffer)    alloc_if(0) free_if(1)) \
        nocopy(mateSeedOutputBufferFrame:length(memsize_seedSraOutputBuffer)    alloc_if(0) free_if(1)) \
        nocopy(readSeedOutputBufferMeta:length(memsize_seedSraOutputBuffer)     alloc_if(0) free_if(1)) \
        nocopy(mateSeedOutputBufferMeta:length(memsize_seedSraOutputBuffer)     alloc_if(0) free_if(1)) \
        \
        nocopy(seedPeOutputBufferFrame:length(memsize_seedPeOutputBuffer)       alloc_if(0) free_if(1)) \
        nocopy(seedPeOutputBufferMeta:length(memsize_seedPeOutputBuffer)        alloc_if(0) free_if(1)) \
        \
        nocopy(mappingQualities: length(memsize_bufferPerRead)                  alloc_if(0) free_if(1))
    #endif
    
    MICMemFree(micMem);
    
    ////////////////////////////////////////////////////////////
    // Free Memory Allocated
    ////////////////////////////////////////////////////////////
    free(readNameBufferFrame);
    free(readBodyBufferFrame);
    free(readLengthBufferFrame);
    free(readQualityBufferFrame);
    free(mateNameBufferFrame);
    free(mateBodyBufferFrame);
    free(mateLengthBufferFrame);
    free(mateQualityBufferFrame);
    free(occCount);
    free(outputBufferFrame);
    free(outputBufferStatus);
    free(outputBufferMeta);
    free(readOutputBufferFrame);
    free(readOutputBufferMeta);
    free(mateOutputBufferFrame);
    free(mateOutputBufferMeta);
    free(cpPModels);
    free(cpPModels_seed);
    free(cpPModels_ext);
    free(cpNModels);
    free(cpNModels_seed);
    free(cpNModels_ext);
    free(dpOutputBufferFrame);
    free(dpWorkingBuffer);
    free(readSeedOutputBufferFrame);
    free(mateSeedOutputBufferFrame);
    free(readSeedOutputBufferMeta);
    free(mateSeedOutputBufferMeta);
    free(seedPeOutputBufferFrame);
    free(seedPeOutputBufferMeta);
    free(mappingQualities);
    free(g_log_n);
    printf("DONE\n");

    return NULL;

}

void MICCTPrintStats (MICControllerThread * thread, Logging * logger) {
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "MICC#%d - Number of Query Handled = %llu\n", thread->threadId,thread->statsQueryHandled);
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "MICC#%d - Read Loading Time       = %9.4f seconds\n", thread->threadId,thread->statsTotalReadLoadTime);
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "MICC#%d - Model Building Time     = %9.4f seconds\n", thread->threadId,thread->statsTotalModelBuildTime);
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "MICC#%d - Wait Time               = %9.4f seconds\n", thread->threadId,thread->statsTotalWaitTime);
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "MICC#%d - Alignment Time          = %9.4f seconds\n", thread->threadId,thread->statsTotalAlignmentTime);
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "MICC#%d -                         ( %9.4f rps ) \n", thread->threadId, thread->statsQueryHandled / thread->statsTotalAlignmentTime);
    LOGWriteLine(logger, LOGGING_MSG_TYPE_ALL, "\n");
}

void MICCTStart(MICControllerThread * thread) {
    pthread_create(thread->t, NULL, start, thread);
}

void MICCTJoin(MICControllerThread * thread) {
    pthread_join(*thread->t, NULL);
}

void MICCTFree(MICControllerThread * thread) {
    free(thread->t);
    free(thread);
}

void MICCTPollHealthString(MICControllerThread * thread) {
    sprintf(thread->consolidHealth,"%s",thread->healthString);
}

char * MICCTGetHealthText(MICControllerThread * thread) {
    return thread->consolidHealth;
}
