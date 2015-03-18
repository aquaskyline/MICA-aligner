//
//    MemMan.h
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

#ifndef __MEMMAN_H__
#define __MEMMAN_H__

#include <stdio.h>
#include <stdlib.h>

#define MEMORY_TYPE_CPU         0
#define MEMORY_TYPE_SHARED      1
#define MEMORY_TYPE_COUNT       2

typedef struct MEMMan {
    unsigned long long allocMemoryByte[MEMORY_TYPE_COUNT];
} MEMMan;

MEMMan * MEMManCreate();
unsigned long long MEMManGetUsage(int type);
void * MEMManMalloc(size_t size, int type);

#endif

