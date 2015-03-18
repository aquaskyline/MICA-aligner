#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "../SRAConfigParser.h"

int main(int argc, char ** argv) {
    printf("[MAIN] Configuration Reader is invoked.\n");
    printf("Size of SRACFGReader           %llu bytes\n", (unsigned long long)sizeof(SRACFGReader));
    printf("Size of SRACFG                 %llu bytes\n", (unsigned long long)sizeof(SRACFG));
    
    UTBFRBuffer * buf = UTBFRLoad("../sraModel.8G.cfg");
    SRACFG * cfg = SRACPLoad(buf);
    //SRACPDebugPrint(cfg);
    
    /*char tokenValue[SRA_CFGPARSE_MAX_TOKEN_LENGTH];
    _SRACPReaderEnter(cfg,"alignmentModel");
    while (_SRACPReaderEnter(cfg,"alignment")) {
        _SRACPReaderGetValue(cfg,"errorType",tokenValue);
        printf("token value = %s\n",tokenValue);
        _SRACPReaderGetValue(cfg,"maxNumOfError",tokenValue);
        printf("token value = %s\n",tokenValue);
        _SRACPReaderEnter(cfg,"case");
        _SRACPReaderGetValue(cfg,"type",tokenValue);
        printf("token value = %s\n",tokenValue);
        _SRACPReaderExit(cfg);
        _SRACPReaderExit(cfg);
    }*/
    
    SRASetting sraQuerySettingSocket;
    sraQuerySettingSocket.ReadStrand=QUERY_BOTH_STRAND;
    sraQuerySettingSocket.ErrorType=SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
    sraQuerySettingSocket.OutputType=SRA_REPORT_ALL;
    sraQuerySettingSocket.OutputFormat=SRA_OUTPUT_FORMAT_SAM;
    sraQuerySettingSocket.MaxResult=1000;
    
    SRAIndex sraIndex;
    
    SRAModel * SRAMismatchModel;
    SRAModel * SRAMismatchModel_cfg = (SRAModel*) malloc(sizeof(SRAModel));
    
    //////////////////////////////////////
    // OutputType
    //////////////////////////////////////
    //sraQuerySettingSocket.OutputType=SRA_REPORT_ALL_BEST;
    sraQuerySettingSocket.OutputType=SRA_REPORT_ALL;
    sraQuerySettingSocket.MaxError=4;
    
    //////////////////////////////////////
    // ReadLength
    //////////////////////////////////////
    int ReadLength = 55;
    
    //////////////////////////////////////
    // ErrorType
    //////////////////////////////////////
    printf("Creating SRA Model\n");
    sraQuerySettingSocket.ErrorType=SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
    SRAModelSet * modelSet = SRAModelSetConstruct(&sraQuerySettingSocket,&sraIndex,SRA_MODEL_8G,SRA_MIN_READ_LENGTH,SRA_MAX_READ_LENGTH); 
    SRAModelConstruct(modelSet,ReadLength);
    SRAMismatchModel = SRAModelSetGetModel(modelSet,ReadLength,QUERY_POS_STRAND);
    sraQuerySettingSocket.ErrorType=SRA_TYPE_MISMATCH_ONLY;
    SRACPLoadModel(cfg,SRAMismatchModel_cfg,ReadLength,&sraQuerySettingSocket,&sraIndex,SRA_MODEL_8G); 
    
    FILE * ofp = fopen("cfgOutput.old","w");
    DebugPrintModel(SRAMismatchModel,ofp);
    
    FILE * nfp = fopen("cfgOutput.new","w");
    DebugPrintModel(SRAMismatchModel_cfg,nfp);
    
    //////////////////////////////////////
    // Exhaustion
    //////////////////////////////////////
    #ifndef DEBUG_DISABLE_EXHAUSTION
    int outputTypes[5] = {SRA_REPORT_UNIQUE_BEST, SRA_REPORT_RANDOM_BEST, SRA_REPORT_ALL,
                            SRA_REPORT_ALL_BEST, SRA_REPORT_BEST_QUALITY };
    int outputIdx = 0;
    int errorIdx;
    //for (outputIdx=0;outputIdx<5;outputIdx++) {
    //    for (ReadLength=35;ReadLength<200;ReadLength++) {
            for (errorIdx=0;errorIdx<=5;errorIdx++) {
                sraQuerySettingSocket.OutputType=outputIdx;
                sraQuerySettingSocket.MaxError = errorIdx;
                ReadLength=100;
                
                printf("Creating SRA Model for Out:%d Length:%d Err:%d\n",outputIdx,ReadLength,errorIdx);
                sraQuerySettingSocket.ErrorType=SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
                //SRAMismatchModel = SRAModelConstruct(ReadLength,QUERY_POS_STRAND,&sraQuerySettingSocket,&sraIndex,SRA_MODEL_8G); 
                SRAModelConstruct(modelSet,ReadLength);
                SRAMismatchModel = SRAModelSetGetModel(modelSet,ReadLength,QUERY_POS_STRAND);
                sraQuerySettingSocket.ErrorType=SRA_TYPE_MISMATCH_ONLY;
                SRACPLoadModel(cfg,SRAMismatchModel_cfg,ReadLength,&sraQuerySettingSocket,&sraIndex,SRA_MODEL_8G); 
                
                fprintf(ofp,"Creating SRA Model for Out:%d Length:%d Err:%d\n",outputIdx,ReadLength,errorIdx);
                DebugPrintModel(SRAMismatchModel,ofp);
                fprintf(nfp,"Creating SRA Model for Out:%d Length:%d Err:%d\n",outputIdx,ReadLength,errorIdx);
                DebugPrintModel(SRAMismatchModel_cfg,nfp);
            }
    //    }
    //}
    #endif
    
    fclose(ofp);
    fclose(nfp);
    SRACPFree(cfg);
    UTBFRFree(buf);
    
    return 0;
    // Statements below generates sraModel.cfg from SRA2BWTMdl.c hardcodes.
    // --------------------------------------------------------------------

    int errorCount,errorType;
    
    printf("alignmentModel\n");
    printf("{\n");
    for (errorType=0;errorType<2;errorType++) {
    for (errorCount=0;errorCount<=MAX_NUM_OF_MISMATCH;errorCount++) {
        
        printf("    alignment\n");
        printf("    {\n");
        if (errorType==0) {
            printf("        errorType     SRA_TYPE_MISMATCH_ONLY\n");
        } else {
            printf("        errorType     SRA_TYPE_EDIT_DISTANCE\n");
        }
        printf("        maxNumOfError %u\n",errorCount);

        if (errorType==0) {
            sraQuerySettingSocket.ErrorType=SRA_STEP_ERROR_TYPE_MISMATCH_ONLY;
        } else {
            sraQuerySettingSocket.ErrorType=SRA_STEP_ERROR_TYPE_EDIT_DISTANCE;
        }
        
        sraQuerySettingSocket.MaxError=errorCount;
        //SRAMismatchModel = SRAModelConstruct(100,QUERY_POS_STRAND,&sraQuerySettingSocket,&sraIndex,SRA_MODEL_8G); 
        SRAModelConstruct(modelSet,100);
        SRAMismatchModel = SRAModelSetGetModel(modelSet,100,QUERY_POS_STRAND);
        
        int caseIdx;
        int stepIdx;
        for (caseIdx=0;caseIdx<MAX_NUM_OF_SRA_CASES;caseIdx++) {
            SRACase * thisCase = &(SRAMismatchModel->cases[caseIdx]);
            if (thisCase->type == SRA_CASE_TYPE_NOT_INITALISED) {
                caseIdx = MAX_NUM_OF_SRA_CASES;
                continue;
            }
            
            printf("        case\n");
            printf("        {\n");
            printf("            type %s\n",SRA_CASE_TYPE[thisCase->type]);
            for (stepIdx=0;stepIdx<MAX_NUM_OF_SRA_STEPS;stepIdx++) {
                SRAStep * thisStep = &(thisCase->steps[stepIdx]);
                if (thisStep->type == SRA_STEP_TYPE_COMPLETE) {
                    stepIdx = MAX_NUM_OF_SRA_STEPS;
                    continue;
                }
                
                printf("            step\n");
                printf("            {\n");
                printf("                stepType %s\n",SRA_STEP_TYPE[thisStep->type]);
                printf("                searchRegion %u-%u\n",thisStep->start+1,thisStep->end+1);
                printf("                minError %u\n",thisStep->MinError);
                printf("                maxError %u\n",thisStep->MaxError);
                printf("                errorType %s\n",SRA_STEP_ERROR_TYPE[thisStep->ErrorType]);
                
                if (thisStep->ceThreshold!=0) {
                    printf("                checkExtend\n");
                    printf("                {\n");
                    printf("                    threshold %u\n",thisStep->ceThreshold);
                    int lowRegion, highRegion;
                    lowRegion = thisStep->ceStart;
                    highRegion = thisStep->ceStart;
                    if (thisStep->ceEnd<lowRegion) {
                        lowRegion = thisStep->ceEnd;
                    } else {
                        highRegion = thisStep->ceEnd;
                    }
                    printf("                    region %u-%u\n",lowRegion+1,highRegion+1);
                    printf("                }\n");
                }
                
                printf("            }\n");
            }
            printf("        }\n");
            
        }
        printf("    }\n");
    }}
    
    printf("}\n");

    SRAModelSetFree(modelSet);
}
