//
//    MICA-Health.c
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

#include "MICA-Health.h"

MICAHealthParm * MICAHealthParmCreate() {
    MICAHealthParm * micaHealthParm = (MICAHealthParm*) MEMManMalloc(sizeof(MICAHealthParm), MEMORY_TYPE_CPU);
    
    micaHealthParm->readThreadCount = 0;
    micaHealthParm->micControllerCount = 0;
    micaHealthParm->cpuControllerCount = 0;
    
    micaHealthParm->threadBody = 0;
    micaHealthParm->exitSignal = 0;
    
    return micaHealthParm;
}

void MICAHealthParmFree(MICAHealthParm * micaHealthParm) {
    free(micaHealthParm);
}

void MICAAddMonitor_ReadThread(MICAHealthParm * micaHealthParm, PEReadThread * readThread) {
    micaHealthParm->readThreads[micaHealthParm->readThreadCount++] = readThread;
}

void MICAAddMonitor_WriteThread(MICAHealthParm * micaHealthParm, PEWriteThread * writeThread) {
    micaHealthParm->writeThreads[micaHealthParm->writeThreadCount++] = writeThread;
}

void MICAAddMonitor_MICControl(MICAHealthParm * micaHealthParm, MICControllerThread * micController) {
    micaHealthParm->micControllers[micaHealthParm->micControllerCount++] = micController;
}

void MICAAddMonitor_CPUControl(MICAHealthParm * micaHealthParm, CPUControllerThread * cpuController) {
    micaHealthParm->cpuControllers[micaHealthParm->cpuControllerCount++] = cpuController;
}


void * _MICAHealthParmRollBody (void * _micaHealthParm) {
    MICAHealthParm * micaHealthParm = (MICAHealthParm*) _micaHealthParm;
    
    int i;
    while (!micaHealthParm->exitSignal) {
    
        printf("\n[MICA Health Report]\n");
        printf("  ReadThd(%d) = ",micaHealthParm->readThreadCount); fflush(stdout);
        for (i=0;i<micaHealthParm->readThreadCount;i++) {
            if (i>0) printf(" / ");
            PEReadThread_PollHealthString(micaHealthParm->readThreads[i]);
            printf("%s",PEReadThread_GetHealthText(micaHealthParm->readThreads[i]));
        }
        printf("\n");
        
        printf("  MICCtl(%d)  = ",micaHealthParm->micControllerCount); fflush(stdout);
        for (i=0;i<micaHealthParm->micControllerCount;i++) {
            if (i>0) printf(" / ");
            MICCTPollHealthString(micaHealthParm->micControllers[i]);
            printf("%s",MICCTGetHealthText(micaHealthParm->micControllers[i]));
        }
        printf("\n");
        
        printf("  CPUCtl(%d)  = ",micaHealthParm->cpuControllerCount); fflush(stdout);
        for (i=0;i<micaHealthParm->cpuControllerCount;i++) {
            if (i>0) printf(" / ");
            CPUCTPollHealthString(micaHealthParm->cpuControllers[i]);
            printf("%s",CPUCTGetHealthText(micaHealthParm->cpuControllers[i]));
        }
        printf("\n");

        printf("  WteThd(%d) = ",micaHealthParm->writeThreadCount); fflush(stdout);
        for (i=0;i<micaHealthParm->writeThreadCount;i++) {
            if (i>0) printf(" / ");
            PEWTPollHealthString(micaHealthParm->writeThreads[i]);
            printf("%s",PEWTGetHealthText(micaHealthParm->writeThreads[i]));
        }
        
        printf("\n");
        
        // Sleep
        sleep(MICA_HEALTH_POLL_SECOND);
    }
    return 0;
}

void MICAHealthParmRoll (MICAHealthParm * micaHealthParm) {

    pthread_t * threadBody = &(micaHealthParm->threadBody);
    
    if ((*threadBody)!=0) {
        printf("[micaHealthParm] Unable to re-roll a MICA Health when it is running.\n");
        return;
    }
    
    if (pthread_create(threadBody, NULL, _MICAHealthParmRollBody, (void*) micaHealthParm)) {
        printf("[micaHealthParm] Unable to create thread!\n");
        exit(1);
    }
    
}

void MICAHealthParmJoin(MICAHealthParm * micaHealthParm) {

    pthread_t * threadBody = &(micaHealthParm->threadBody);
    
    if ((*threadBody)!=0) {
        micaHealthParm->exitSignal = 1;
        pthread_join((*threadBody),NULL);
    }
    
    micaHealthParm->threadBody = 0;
    micaHealthParm->exitSignal = 0;
}
