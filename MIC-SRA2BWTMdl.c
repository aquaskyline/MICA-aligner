//
//    MIC-SRA2BWTMdl.c
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

#include "MIC-SRA2BWTMdl.h"


void MICSRA2BWTStepPopulate(CPTSRAStep * cpStep, SRAStep * sraStep) {

    // Step Level
    cpStep->type = sraStep->type;
    
#ifdef MIC_2BWT_MODEL_DISABLE_LOOKUP
    // If the compiler configuration item is defined,
    // disable Lookup Table Temporarily on MIC card.
    if (cpStep->type == SRA_STEP_TYPE_BI_DIRECTIONAL_FORWARD_LOOKUP ||
        cpStep->type == SRA_STEP_TYPE_BI_DIRECTIONAL_BACKWARD_LOOKUP ) {
        cpStep->type = SRA_STEP_TYPE_BI_DIRECTIONAL_BWT;
    } else if (cpStep->type == SRA_STEP_TYPE_BACKWARD_ONLY_LOOKUP) {
        cpStep->type = SRA_STEP_TYPE_BACKWARD_ONLY_BWT;
    }
#endif

    cpStep->start = sraStep->start;
    cpStep->end = sraStep->end;
    
    cpStep->ErrorType = sraStep->ErrorType;
    cpStep->MinError = sraStep->MinError;
    cpStep->MaxError = sraStep->MaxError;
    
    cpStep->step = sraStep->step;
    cpStep->len = sraStep->len;
}

void MICSRA2BWTCasePopulate(CPTSRACase * cpCase, SRACase * sraCase) {

    int stepIdx;
    int reachStop = 0;
    
    // Case Level
    cpCase->id = sraCase->id;
    cpCase->type = sraCase->type;
    
    for (stepIdx=0;stepIdx<MAX_NUM_OF_CPT_SRA_STEPS;stepIdx++) {
    
        CPTSRAStep * cpStep = &cpCase->steps[stepIdx];
        SRAStep * sraStep = &sraCase->steps[stepIdx];
        
        MICSRA2BWTStepPopulate(cpStep,sraStep);
        
        if (sraStep->type == SRA_STEP_TYPE_COMPLETE) {
            reachStop = 1;
            break;
        }
    }
    
    if ( !reachStop ) {
        printf("[ERROR] CPTSRAModel is not large enough to take in the entire SRAModel.\n");
        printf("Please review parameter MAX_NUM_OF_CPT_SRA_STEPS(%d)\n",MAX_NUM_OF_CPT_SRA_STEPS);
        for (;stepIdx<MAX_NUM_OF_SRA_STEPS;stepIdx++) {
            SRAStep * sraStep = &sraCase->steps[stepIdx];
            if ( sraStep->type == SRA_STEP_TYPE_COMPLETE ) {
                break;
            }
        }
        printf("The SRA full case has %d steps.\n",stepIdx);
    }

}

void MICSRA2BWTModelPopulate(CPTSRAModel * cpModel, SRAModel * sraModel) {

    int caseIdx;
    int reachStop = 0;
    
    // Model Level
    // -- Empty currently
    
    for (caseIdx=0;caseIdx<MAX_NUM_OF_CPT_SRA_CASES;caseIdx++) {
    
        CPTSRACase * cpCase = &cpModel->cases[caseIdx];
        SRACase * sraCase   = &sraModel->cases[caseIdx];
        
        MICSRA2BWTCasePopulate(cpCase,sraCase);
        
        if ( sraCase->type == SRA_CASE_TYPE_NOT_INITALISED ) {
            reachStop = 1;
            break;
        }
    }
    
    if ( !reachStop ) {
        printf("[ERROR] CPTSRAModel is not large enough to take in the entire SRAModel.\n");
        printf("Please review parameter MAX_NUM_OF_CPT_SRA_CASES(%d)\n",MAX_NUM_OF_CPT_SRA_CASES);
        for (;caseIdx<MAX_NUM_OF_SRA_CASES;caseIdx++) {
            SRACase * sraCase   = &sraModel->cases[caseIdx];
            if ( sraCase->type == SRA_CASE_TYPE_NOT_INITALISED ) {
                break;
            }
        }
        printf("The SRA full model has %d cases.\n",caseIdx);
    }
}

__attribute__((target(mic)))
void MICSRADebugPrintModel(CPTSRAModel * cpModel, FILE * stream) {
    int i,j,k;
    fprintf(stream,"Debug printing alignment model:\n");
    for (j=0;j<MAX_NUM_OF_CPT_SRA_CASES;j++) {
        if (cpModel->cases[j].type==SRA_CASE_TYPE_NOT_INITALISED) break;
        fprintf(stream,"    Case %u\n", j);
        fprintf(stream,"    Case Type %s\n", SRA_CASE_TYPE[cpModel->cases[j].type]);
        for (k=0;k<MAX_NUM_OF_CPT_SRA_STEPS;k++) {
        
            fprintf(stream,"        Step type : %s\n",SRA_STEP_TYPE[cpModel->cases[j].steps[k].type]);
            fprintf(stream,"        Step ErrorType : %s\n",SRA_STEP_ERROR_TYPE[cpModel->cases[j].steps[k].ErrorType]);
            
            
            if (cpModel->cases[j].steps[k].type == SRA_STEP_TYPE_COMPLETE) 
            break;
            fprintf(stream,"        From %d to %d allowing mismatch %d-%d\n", 
            cpModel->cases[j].steps[k].start,
            cpModel->cases[j].steps[k].end,
            cpModel->cases[j].steps[k].MinError,
            cpModel->cases[j].steps[k].MaxError);
                
            fprintf(stream,"        Computed value :\n");
            fprintf(stream,"        - step %u\n",cpModel->cases[j].steps[k].step);
            fprintf(stream,"        - len %u\n",cpModel->cases[j].steps[k].len);
        }
        
        // Print left aligned for CE
        //for (k=0;k<cpModel->ModelReadLength;k++) {
        //    fprintf(stream,"%d ",cpModel->cases[j].leftMostAligned[k]);
        //}
        //fprintf(stream,"\n");
                    
        
    }
}
