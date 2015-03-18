//
//    MICA-Health.h
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

#ifndef __MICA_HEALTH_H__
#define __MICA_HEALTH_H__

///////////////////////////////////////
// Including Standard Libraries
///////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>
#include <pthread.h>
#include <unistd.h>
#include <errno.h>

///////////////////////////////////////
// Including MIC Libraries
///////////////////////////////////////
#include "MICA-PE-ReadThread.h"
#include "MICControllerThread.h"
#include "CPUControllerThread.h"

#define MICA_HEALTH_MAX_READ_THREAD    1
#define MICA_HEALTH_MAX_WRITE_THREAD   1
#define MICA_HEALTH_MAX_MIC_CONTROL    10
#define MICA_HEALTH_MAX_CPU_CONTROL    10

#define MICA_HEALTH_POLL_SECOND        2

typedef struct MICAHealthParm {
    PEReadThread * readThreads[MICA_HEALTH_MAX_READ_THREAD];
    PEWriteThread * writeThreads[MICA_HEALTH_MAX_WRITE_THREAD];
    MICControllerThread * micControllers[MICA_HEALTH_MAX_MIC_CONTROL];
    CPUControllerThread * cpuControllers[MICA_HEALTH_MAX_CPU_CONTROL];
    
    int readThreadCount;
    int writeThreadCount;
    int micControllerCount;
    int cpuControllerCount;
    
    // POSIX Thread
    pthread_t       threadBody;
    int             exitSignal;
} MICAHealthParm;


MICAHealthParm * MICAHealthParmCreate();
void MICAHealthParmFree(MICAHealthParm * micaHealthParm);

void MICAAddMonitor_ReadThread(MICAHealthParm * micaHealthParm, PEReadThread * readThread);
void MICAAddMonitor_WriteThread(MICAHealthParm * micaHealthParm, PEWriteThread * writeThread);
void MICAAddMonitor_MICControl(MICAHealthParm * micaHealthParm, MICControllerThread * micController);
void MICAAddMonitor_CPUControl(MICAHealthParm * micaHealthParm, CPUControllerThread * cpuController);

void MICAHealthParmRoll (MICAHealthParm * micaHealthParm);
void MICAHealthParmJoin(MICAHealthParm * micaHealthParm);

#endif

