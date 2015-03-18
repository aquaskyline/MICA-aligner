//
//    2BWT-Viewer.h
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
/* 

    SOAP2 Viewer for SOAP2 Aligner
    Date:    2010-11-13
    Author:  Edward Wu
    Website: http://www.cs.hku.hk/~mkewu/2bwt

    Designed for output format 20101113

    Modification History 

    Date   : 3rh July 2011
    Author : Edward MK Wu
    Change : New file.

*/

#ifndef __2BWT_VIEWER_H__
#define __2BWT_VIEWER_H__

#include "SRACore.h"
#include "SRAQueryParser.h"
#include "OutputAnyFormat.h"
#include "OutputBinaryFormat.h"

#define SOCKET_OPERATION      0
#define LOCAL_ALIGNMENT       1
#define SHORT_READ_ALIGNMENT  2
#define PAIR_END_ALIGNMENT    3

typedef struct SRAQueryBatch{
    char * names;
    unsigned char * patterns;
    unsigned int batchSize;
    unsigned long long startIdx;
    unsigned long long endIdx;
} SRAQueryBatch;

#endif
