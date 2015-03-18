//
//    testWriteThread.c
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

#include "testWriteThread.h"


int main(int argc, char** argv) {

    PEWriteThread * myThread = PEWTCreate(1000,240,1,1,1);
    
    printf("[MAIN] Roll the Write Thread.\n");
    PEWTRoll(myThread);
    
    int sleepSec = 5;
    printf("[MAIN] Sleeping for %d seconds.\n",sleepSec);
    sleep(sleepSec);
    printf("[MAIN] Awake!\n");
    
    printf("\n");
    PEWTCloneBuffer(myThread,1000,240,0,0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0,NULL,NULL,NULL,NULL,NULL,NULL);
    PEWTCloneBuffer(myThread,1000,240,0,0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0,NULL,NULL,NULL,NULL,NULL,NULL);
    PEWTCloneBuffer(myThread,1000,240,0,0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0,NULL,NULL,NULL,NULL,NULL,NULL);
    PEWTCloneBuffer(myThread,1000,240,0,0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0,NULL,NULL,NULL,NULL,NULL,NULL);
    PEWTCloneBuffer(myThread,1000,240,0,0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0,NULL,NULL,NULL,NULL,NULL,NULL);
    printf("\n");
    
    printf("[MAIN] Sleeping for %d seconds.\n",sleepSec);
    sleep(sleepSec);
    printf("[MAIN] Awake!\n");
    
    printf("\n");
    PEWTFree(myThread);

    return 0;
}
