//
//    SRAConfigParser.c
//
//    SOAP2 / 2BWT
//
//    Copyright (C) 2012, HKU
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
/////////////////////////////////////////////////////
/*
    File identifier: file:///usr/local/svn/project/2bwt-flex/latest/SRAConfigParser.c
    2BWT-Flex Library
    SRAConfigParser is a CFG parser that reads in a presumed format CFG

    Modification History

    Date   : 16th December 2012
    Author : Edward MK Wu
    Change : New file.
    
*/
/////////////////////////////////////////////////////

#include "SRAConfigParser.h"

inline uint8_t _SRACPGetNextToken ( UTBFRBuffer * cfgBuffer, char * lastCharacter, char * token ) {
    char c = *lastCharacter;
    
    int tokenLength = 0;
    int isInsideQuote = 0;
    
    while (!UTBFREOF(cfgBuffer)) {
    
        token[tokenLength]=c;
    
        if ( c == '{' ) {
            if (tokenLength==0) {
                (*lastCharacter) = UTBFRGetCharacter(cfgBuffer);
                token[tokenLength] = '\0';
                return SRA_CFGPARSE_TYPE_BLOCK_ENTRY;
            } else {
                break;
            }
        } else if ( c == '}' ) {
            if (tokenLength==0) {
                (*lastCharacter) = UTBFRGetCharacter(cfgBuffer);
                token[tokenLength] = '\0';
                return SRA_CFGPARSE_TYPE_BLOCK_EXIT;
            } else {
                break;
            }
        } else if ( c == ' ' || c == '\t' ||  c == '\n'  ||  c == '\r' ) {
            if (isInsideQuote) {
                tokenLength++;
            } else {
                if (tokenLength!=0) {
                    break;
                }
            }
        } else {
            tokenLength++;
        }
    
        c = UTBFRGetCharacter(cfgBuffer);
    }
    
    (*lastCharacter) = c;
    token[tokenLength] = '\0';
    
    if (tokenLength==0) {
        return SRA_CFGPARSE_TYPE_NONE;
    }
    
    return 1;
}

inline void _SRACPEnqueueNode ( SRACFG * cfg, int16_t levelIdx, int16_t tokenType, int16_t tokenIdx )
{
    SRACFGLevel * level = &(cfg->levels[levelIdx]);
    SRACFGNode * node = &(level->nodes[level->nodeCount]);
    node->type = tokenType;
    node->index = tokenIdx;
    level->nodeCount++;
}

SRACFG * SRACPLoad ( UTBFRBuffer * cfgBuffer) {

    int i,j;
    char c = UTBFRGetCharacter(cfgBuffer);
    char token[SRA_CFGPARSE_MAX_TOKEN_LENGTH];
    
    SRACFG * newCfg = ( SRACFG * ) malloc(sizeof(SRACFG));
    for (i=0;i<SRA_CFGPARSE_MAX_LEVEL;i++) {
        newCfg->levels[i].nodeCount=0;
        
        for (j=0;j<SRA_CFGPARSE_MAX_TOKEN_PER_LEVEL;j++) {
            newCfg->levels[i].nodes[j].type = SRA_CFGPARSE_TYPE_NONE;
            newCfg->levels[i].nodes[j].index = -1;
        }
    }
    
    int16_t tokenIdx = 0;
    int16_t levelIdx = 0;
    
    uint8_t tokenType = _SRACPGetNextToken(cfgBuffer,&c,token);
    while (tokenType!=SRA_CFGPARSE_TYPE_NONE) {
        //printf("[SRACFG] Token (%u/%u) obtained %d %s\n",newCfg->levels[levelIdx].nodeCount,SRA_CFGPARSE_MAX_TOKEN_PER_LEVEL,tokenType,token);
        if (tokenType==SRA_CFGPARSE_TYPE_TOKEN) {
        
            _SRACPEnqueueNode(newCfg,levelIdx,tokenType,tokenIdx);
            strcpy(newCfg->tokens[tokenIdx],token);
            tokenIdx++;
            
            if (newCfg->levels[levelIdx].nodeCount>SRA_CFGPARSE_MAX_TOKEN_PER_LEVEL) {
                fprintf(stderr,"[SRACFG] Error while reading %s.\n",token);
                fprintf(stderr,"[SRACFG] Number of token exceeds SRA_CFGPARSE_MAX_TOKEN_PER_LEVEL.\n");
                exit(1);
            }
            
        } else if (tokenType==SRA_CFGPARSE_TYPE_BLOCK_ENTRY) {
        
            if (levelIdx+1>=SRA_CFGPARSE_MAX_LEVEL) {
                fprintf(stderr,"[SRACFG] Number of level exceeds SRA_CFGPARSE_MAX_LEVEL.\n");
                exit(1);
            }
        
            _SRACPEnqueueNode(newCfg,levelIdx,tokenType,newCfg->levels[levelIdx+1].nodeCount);
            levelIdx++;
            
        } else if (tokenType==SRA_CFGPARSE_TYPE_BLOCK_EXIT) {
        
            if (levelIdx-1<0) {
                fprintf(stderr,"[SRACFG] Number of level reduced beyond zero.\n");
                exit(1);
            }
        
            _SRACPEnqueueNode(newCfg,levelIdx,tokenType,0);
            levelIdx--;
        } else {
            fprintf(stderr,"[SRACFG] Unrecognised token type.\n");
            exit(1);
        }
        tokenType = _SRACPGetNextToken(cfgBuffer,&c,token);
    }
    
    return newCfg;
}

void SRACPFree ( SRACFG * cfg ) {
    free(cfg);
}

void SRACPDebugPrint ( SRACFG * cfg ) {
    printf("[SRACPDebugPrint] Debug printing of SRACFG structure\n");
    printf("Size of the structure = %u Bytes\n",(unsigned int)sizeof(SRACFG));
    printf("SRA_CFGPARSE_MAX_LEVEL = %u\n",SRA_CFGPARSE_MAX_LEVEL);
    printf("SRA_CFGPARSE_MAX_TOKEN_COUNT = %u\n",SRA_CFGPARSE_MAX_TOKEN_COUNT);
    printf("SRA_CFGPARSE_MAX_TOKEN_LENGTH = %u\n",SRA_CFGPARSE_MAX_TOKEN_LENGTH);
    
    int i;
    for (i=0;i<SRA_CFGPARSE_MAX_LEVEL;i++) {
        printf("Printing Level #%d\n",i);
        unsigned nodeCount = cfg->levels[i].nodeCount;
        printf("    nodeCount = %u\n",nodeCount);
        int j = 0;
        for (j=0;j<nodeCount;j++) {
            int type = cfg->levels[i].nodes[j].type;
            int index = cfg->levels[i].nodes[j].index;
            if (type==SRA_CFGPARSE_TYPE_TOKEN) {
                printf("    Node[%d] = (TOKEN / %d) %s\n",j,index,cfg->tokens[index]);
            } else if (type==SRA_CFGPARSE_TYPE_BLOCK_ENTRY) {
                printf("    Node[%d] = (ENTRY / %d)\n",j,index);
            } else if (type==SRA_CFGPARSE_TYPE_BLOCK_EXIT) {
                printf("    Node[%d] = (EXIT  / %d)\n",j,index);
            } else {
                printf("    Node[%d] = (UNKNOWN %d / %d)\n",j,type,index);
            }
        }
    }
    
    for (i=0;i<SRA_CFGPARSE_MAX_TOKEN_COUNT;i++) {
        printf("%d : %s\n",i,cfg->tokens[i]);
    }
}

int16_t _SRACPReaderFindToken ( SRACFG * cfg, SRACFGReader * cfgReader, char * tokenName ) {
    int i = 0;
    int16_t levelIdx = cfgReader->readerLevelIdx;
    int16_t nodeIdx = cfgReader->levels[levelIdx].readerCount;
    SRACFGNode * nodes = cfg->levels[levelIdx].nodes;
    SRACFGReaderNode * readerNodes = cfgReader->levels[levelIdx].nodes;
    while ( nodes[nodeIdx].type != SRA_CFGPARSE_TYPE_BLOCK_EXIT &&
            nodes[nodeIdx].type != SRA_CFGPARSE_TYPE_NONE) {
            
        if (nodes[nodeIdx].type == SRA_CFGPARSE_TYPE_TOKEN && 
            readerNodes[nodeIdx].readerFlag != SRA_CFGPARSE_READER_READ &&
            strcmp(cfg->tokens[nodes[nodeIdx].index],tokenName)==0) {
            return nodeIdx;
        }

        nodeIdx++;
    }
    return -1;
}

int8_t _SRACPReaderEnter ( SRACFG * cfg, SRACFGReader * cfgReader, char * blockName ) {
    int16_t nodeIdx = _SRACPReaderFindToken(cfg,cfgReader,blockName);
    if (nodeIdx==-1) {
        //fprintf(stderr,"[SRACFG] Failed to find token %s.\n",blockName);
        return 0;
    } else {
        if (cfg->levels[cfgReader->readerLevelIdx].nodes[nodeIdx+1].type != SRA_CFGPARSE_TYPE_BLOCK_ENTRY) {
            fprintf(stderr,"[SRACFG] token %s is not a block.\n",blockName);
            exit(1);
        } else {
            cfgReader->levels[cfgReader->readerLevelIdx].nodes[nodeIdx].readerFlag = SRA_CFGPARSE_READER_READ;
            int16_t blockNodeIdx = cfg->levels[cfgReader->readerLevelIdx].nodes[nodeIdx+1].index;
            cfgReader->readerLevelIdx++;
            cfgReader->levels[cfgReader->readerLevelIdx].readerCount = blockNodeIdx;
        }
    }
    return 1;
}

void _SRACPReaderExit ( SRACFG * cfg, SRACFGReader * cfgReader ) {
    if (cfgReader->readerLevelIdx<=0) {
        fprintf(stderr,"[SRACFG] Failed to exit block.\n");
        exit(1);
    } else {
        cfgReader->readerLevelIdx--;
    }
}

int8_t _SRACPReaderGetValue ( SRACFG * cfg, SRACFGReader * cfgReader, char * tokenName, char * tokenValue, char * defaultValue ) {
    int16_t nodeIdx = _SRACPReaderFindToken(cfg,cfgReader,tokenName);
    if (nodeIdx==-1 && defaultValue==NULL) {
        fprintf(stderr,"[SRACFG] Failed to find token %s.\n",tokenName);
        exit(1);
    } else if (defaultValue!=NULL) {
        strcpy(tokenValue,defaultValue);
        return 1;
    } else {
        if (cfg->levels[cfgReader->readerLevelIdx].nodes[nodeIdx+1].type != SRA_CFGPARSE_TYPE_TOKEN) {
            fprintf(stderr,"[SRACFG] token %s does not carry a value.\n",tokenName);
            exit(1);
        } else {
            int16_t tokenIdx = cfg->levels[cfgReader->readerLevelIdx].nodes[nodeIdx+1].index;
            strcpy(tokenValue,cfg->tokens[tokenIdx]);
        }
    }
    return 1;
}


int8_t SRACPConstEnumModelType ( char * constantKey ) {
    int8_t index = 0;
    while ( index<2 ) {
        if (strcmp(SRA_TYPE[index],constantKey)==0) {
            return index;
        }
        index++;
    }
    fprintf(stderr,"[SRACPConstEnumModelType] Cannot map key %s to enum.\n",constantKey);
    exit (1);
    return -1;
}

int8_t SRACPConstEnumErrorType ( char * constantKey ) {
    int8_t index = 0;
    while ( index<4 ) {
        if (strcmp(SRA_STEP_ERROR_TYPE[index],constantKey)==0 ) {
            return index;
        }
        index++;
    }
    fprintf(stderr,"[SRACPConstEnumErrorType] Cannot map key %s to enum.\n",constantKey);
    exit (1);
    return -1;
}

int8_t SRACPConstEnumCaseType ( char * constantKey ) {
    int8_t index = 0;
    while ( index<3 ) {
        if (strcmp(SRA_CASE_TYPE[index],constantKey)==0) {
            return index;
        }
        index++;
    }
    fprintf(stderr,"[SRACPConstEnumCaseType] Cannot map key %s to enum.\n",constantKey);
    exit (1);
    return -1;
}

int8_t SRACPConstEnumStepType ( char * constantKey ) {
    int8_t index = 0;
    while ( index<6 ) {
        if (strcmp(SRA_STEP_TYPE[index],constantKey)==0) {
            return index;
        }
        index++;
    }
    fprintf(stderr,"[SRACPConstEnumStepType] Cannot map key %s to enum.\n",constantKey);
    exit (1);
    return -1;
}

int8_t _SRACPReadRegion ( char * region, int * start, int * end ) {
    int i = 0;
    int j = 0;
    char buffer[1024];
    while ( i < strlen(region) && region[i]!='-' ) {
        buffer[j++] = region[i];
        i++;
    }
    if ( i >= strlen(region) || j == 0) {
        fprintf(stderr,"[SRACPReadRegion] Cannot parse region string %s.\n",region);
        exit (1);
        return 0;
    }
    buffer[j] = '\0';
    (*start) = atoi(buffer);
    i++;
    j=0;
    while ( i < strlen(region) ) {
        buffer[j++] = region[i];
        i++;
    }
    if ( j == 0) {
        fprintf(stderr,"[SRACPReadRegion] Cannot parse region string %s.\n",region);
        exit (1);
        return 0;
    }
    buffer[j] = '\0';
    (*end) = atoi(buffer);
    return 1;
}


void SRACPLoadModel ( SRACFG * cfg, SRAModel * SRAMismatchModel,
                    unsigned int ReadLength,
                    SRASetting * aSettings, SRAIndex * aIndex, int modelId ) {
                    
    BWT * bwt = aIndex->bwt;
    BWT * rev_bwt = aIndex->rev_bwt;
    unsigned int readRegion[100];
    uint8_t ErrorType = aSettings->ErrorType;
    uint8_t MaxError = aSettings->MaxError;
    uint8_t smallErrorFirst = 0;
    if ( aSettings->OutputType == SRA_REPORT_UNIQUE_BEST ||
         aSettings->OutputType == SRA_REPORT_RANDOM_BEST ||
         aSettings->OutputType == SRA_REPORT_ALL_BEST ||
         aSettings->OutputType == SRA_REPORT_ALL_SORTED) {
        smallErrorFirst = 1;
    }
    
    /////////////////////////////////////////////////////////
    // Model Initialisation
    /////////////////////////////////////////////////////////
    
    int i,j,k,l;
    SRAMismatchModel->ModelReadLength = ReadLength;
    memset(readRegion,0,sizeof(unsigned int)*100);
    if (ReadLength==0) return;
    int caseIdx = 0;
    int stepIdx = 0;
    
    for (caseIdx=0;caseIdx<MAX_NUM_OF_SRA_CASES;caseIdx++) {
        SRACase * thisCase = &(SRAMismatchModel->cases[caseIdx]);
        thisCase->id=0;
        thisCase->type=SRA_CASE_TYPE_NOT_INITALISED;
        for (stepIdx=0;stepIdx<MAX_NUM_OF_SRA_STEPS;stepIdx++) {
            SRAStep * thisStep = &(thisCase->steps[stepIdx]);
            thisStep->type = SRA_STEP_TYPE_COMPLETE;
            thisStep->start = 0;
            thisStep->end = 0;
            thisStep->ceThreshold = 0;
            thisStep->ceStart = 0;
            thisStep->ceEnd = ReadLength-1;
            
            thisStep->ErrorType = SRA_STEP_ERROR_TYPE_NO_ERROR | SRA_STEP_ERROR_TYPE_NBM;
            thisStep->MinError = 0;
            thisStep->MaxError = 0;
        }
    }
    
    /////////////////////////////////////////////////////////
    // Loading Model from File
    /////////////////////////////////////////////////////////
    SRACFGReader cfgReader;
    memset(&cfgReader,0,sizeof(SRACFGReader));
    char tokenValue[SRA_CFGPARSE_MAX_TOKEN_LENGTH];
    if (!_SRACPReaderEnter(cfg,&cfgReader,"alignmentModel")) {
        fprintf(stderr,"[SRACPLoadModel] Failed to load SRA Alignment Model. Missing mandatory block alignmentModel\n");
        exit(1);
    }
    caseIdx = 0;
    while (_SRACPReaderEnter(cfg,&cfgReader,"alignment")) {
        _SRACPReaderGetValue(cfg,&cfgReader,"errorType",tokenValue,NULL);
        
        int8_t blockErrorType = SRACPConstEnumModelType(tokenValue);
        
        if (blockErrorType==ErrorType) {
            
            _SRACPReaderGetValue(cfg,&cfgReader,"maxNumOfError",tokenValue,NULL);
            int8_t blockNumOfError = atoi(tokenValue);
            
            int8_t readModel = (smallErrorFirst && MaxError >= blockNumOfError) ||
                                (!smallErrorFirst && MaxError == blockNumOfError);
                                
            if (readModel) {
                while (_SRACPReaderEnter(cfg,&cfgReader,"case")) {
                    SRACase * thisCase = &(SRAMismatchModel->cases[caseIdx]);
                
                    _SRACPReaderGetValue(cfg,&cfgReader,"type",tokenValue,NULL);
                    thisCase->type = SRACPConstEnumCaseType(tokenValue);
                
                    int stepIdx = 0;
                    while (_SRACPReaderEnter(cfg,&cfgReader,"step")) {
                        SRAStep * thisStep = &(thisCase->steps[stepIdx]);
                        
                        _SRACPReaderGetValue(cfg,&cfgReader,"stepType",tokenValue,NULL);
                        thisStep->type = SRACPConstEnumStepType(tokenValue);
                        
                        _SRACPReaderGetValue(cfg,&cfgReader,"searchRegion",tokenValue,NULL);
                        _SRACPReadRegion(tokenValue,&(thisStep->start),&(thisStep->end));
                        thisStep->start--;
                        thisStep->end--;
                        
                        _SRACPReaderGetValue(cfg,&cfgReader,"maxError",tokenValue,NULL);
                        thisStep->MaxError = atoi(tokenValue);
                        
                        _SRACPReaderGetValue(cfg,&cfgReader,"minError",tokenValue,NULL);
                        thisStep->MinError = atoi(tokenValue);
                        
                        _SRACPReaderGetValue(cfg,&cfgReader,"errorType",tokenValue,NULL);
                        thisStep->ErrorType = SRACPConstEnumErrorType(tokenValue);
                        
                        if (_SRACPReaderEnter(cfg,&cfgReader,"checkExtend")) {
                                
                            _SRACPReaderGetValue(cfg,&cfgReader,"threshold",tokenValue,NULL);
                            thisStep->ceThreshold = atoi(tokenValue);
                            
                            _SRACPReaderGetValue(cfg,&cfgReader,"region",tokenValue,NULL);
                            _SRACPReadRegion(tokenValue,&(thisStep->ceStart),&(thisStep->ceEnd));
                            thisStep->ceStart--;
                            thisStep->ceEnd--;
                            
                            _SRACPReaderExit(cfg,&cfgReader);
                        }
                        _SRACPReaderExit(cfg,&cfgReader);
                        stepIdx++;
                    }
                    _SRACPReaderExit(cfg,&cfgReader);
                    caseIdx++;
                }
            }
        }
        _SRACPReaderExit(cfg,&cfgReader);
    }
    
    /////////////////////////////////////////////////////////
    // Model Validation
    /////////////////////////////////////////////////////////
    
    for (caseIdx=0;caseIdx<MAX_NUM_OF_SRA_CASES;caseIdx++) {
    
        SRACase * thisCase = &(SRAMismatchModel->cases[caseIdx]);
        
        if (thisCase->type == SRA_CASE_TYPE_NOT_INITALISED) {
            break;
        }
    
        for (stepIdx=0;stepIdx<MAX_NUM_OF_SRA_STEPS;stepIdx++) {
            SRAStep * thisStep = &(thisCase->steps[stepIdx]);
            
            if (thisStep->type == SRA_STEP_TYPE_COMPLETE) {
                break;
            }
            
            // Marking the right-most end of each region
            // Compute the compressed/expanded range index w.r.t to ReadLength
            int max = thisStep->start;
            if (thisStep->end>max) {
                max = thisStep->end;
            }
            readRegion[max] = ReadLength * ((float)max/100.0);
            readRegion[max+1] = readRegion[max] + 1;
            
        }
    }
    readRegion[0] = 0;
    readRegion[ReadLength-1] = ReadLength-1;
    

    /////////////////////////////////////////////////////////
    // Model Post Processing
    /////////////////////////////////////////////////////////
    
    for (caseIdx=0;caseIdx<MAX_NUM_OF_SRA_CASES;caseIdx++) {
    
        SRACase * thisCase = &(SRAMismatchModel->cases[caseIdx]);
        if (thisCase->type != SRA_CASE_TYPE_NOT_INITALISED) {
        
            int leftMost = ReadLength;
            for (stepIdx=0;stepIdx<MAX_NUM_OF_SRA_STEPS;stepIdx++) {
                SRAStep * thisStep = &(thisCase->steps[stepIdx]);
            
                if (thisStep->type != SRA_STEP_TYPE_COMPLETE) {
                    thisStep->start = readRegion[thisStep->start];
                    thisStep->end = readRegion[thisStep->end];
                    
                    if ( thisStep->start == 0 || thisStep->end == 0 ) {
                        printf("[SRACPLoadModel] Error computing the range %u-%u w.r.t to readLength %u\n",
                            thisStep->start,thisStep->end,ReadLength);
                    }
                
                    int start = thisStep->start;
                    int end = thisStep->end;
                    int step = -2 * (start>end) + 1;
                    thisStep->step = step;
                    thisStep->len = (end - start) * step + 1;
                    if (step > 0) {
                        //Forward alignment
                        thisStep->bwt = rev_bwt;
                        if (start<leftMost) {leftMost = start;}
                        for (l=start;l<=end;l+=step) {
                            thisCase->leftMostAligned[l]=leftMost;
                        }
                    } else {
                        //Backward alignment
                        thisStep->bwt = bwt;
                        for (l=start;l>=end;l+=step) {
                            thisCase->leftMostAligned[l]=leftMost--;
                        }
                    }
                    if (smallErrorFirst) {
                        thisStep->MinError = thisStep->MaxError;
                    }
                    if (thisStep->ceThreshold>MAX_CE_THRESHOLD) {
                        thisStep->ceThreshold=MAX_CE_THRESHOLD;
                    }
                    
                    #ifdef DEBUG_2BWT_DISABLE_CHECK_AND_EXTEND
                        thisStep->ceThreshold=0;
                    #endif
                    
                    if (thisStep->ceThreshold>0) {
                        thisStep->ceStart = ReadLength * ((float)thisStep->ceStart/100.0);
                        thisStep->ceEnd = ReadLength * ((float)thisStep->ceEnd/100.0);
                    
                        if (thisStep->ceStart>thisStep->ceEnd) {
                            int swap = thisStep->ceEnd;
                            thisStep->ceEnd=thisStep->ceStart;
                            thisStep->ceStart=swap;
                        }
                    }
                    
                    #ifdef DEBUG_2BWT_DISABLE_LOOKUP_TABLE
                        if (thisStep->type == SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP ||
                            thisStep->type == SRA_STEP_TYPE_BI_DIRECTIONAL_BACKWARD_LOOKUP ) {
                            thisStep->type = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
                        } else if (thisStep->type == SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP) {
                            thisStep->type = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
                        }
                    #endif
                }

            }
        }
    }
}


