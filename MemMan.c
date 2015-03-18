//
//    MemMan.c
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

#include "MemMan.h"

// Singleton
// The freeing of this instance relies on OS.
// It is not done explicitly in the code.
static MEMMan * instance = NULL;

static MEMMan * getInstance() {
    if (instance == NULL) {
        instance = MEMManCreate();
    }

    return instance;
}

// A wrapper to record the total memory allocated for
// the aligner use.
MEMMan * MEMManCreate() {
    MEMMan * memMan = (MEMMan *) malloc (sizeof(MEMMan));
    int i;
    for (i=0;i<MEMORY_TYPE_COUNT;i++) {
        memMan->allocMemoryByte[i] = 0;
    }
    return memMan;
}
void * MEMManMalloc(size_t size, int type) {
    MEMMan * memMan = getInstance();
    memMan->allocMemoryByte[type]+=size;
    return malloc(size);
}

unsigned long long MEMManGetUsage(int type) {
    MEMMan * memMan = getInstance();
    return memMan->allocMemoryByte[type];
}

void MEMManFree(MEMMan * memMan) {
    free(memMan);
}
