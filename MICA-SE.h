//
//    MICA.h
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

#ifndef __MICA_SE_H__
#define __MICA_SE_H__

///////////////////////////////////////
// Including Standard Libraries
///////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <time.h>
#include <omp.h>

///////////////////////////////////////
// Including 2BWT Libraries
///////////////////////////////////////
#include "2bwt-flex/2bwt-lib/iniparser.h"
#include "2bwt-flex/2bwt-lib/2BWT-Interface.h"
#include "2bwt-flex/2bwt-lib/Timing.h"

///////////////////////////////////////
// Including SOAP2-DP Libraries
///////////////////////////////////////
#include "2bwt-flex/2BWT-SRAAlgnmt.h"
#include "2bwt-flex/SRAQueryParser.h"

///////////////////////////////////////
// Including MIC Libraries
///////////////////////////////////////
#include "MIC-SRA2BWTMdl.h"
#include "MIC-SRAAlgnmt.h"
#include "MemMan.h"

#endif

