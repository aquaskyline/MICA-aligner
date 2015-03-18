//
//    2BWT-SRAAlgnmt.h
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

   Modification History 
   
   Date   : 11th April 2012
   Author : Edward MK Wu
   Change : New file.

*/
/////////////////////////////////////////////////////

#ifndef __2BWT_SRA_ALIGNMENT__
#define __2BWT_SRA_ALIGNMENT__

#include "SRA2BWTOperations.h"
#include "DPOperations.h"
#include "SRADPReport.h"
#include "MappingQuality.h"

unsigned long long ProcessStage(SRAArguments * aArgs, SRAModel * sraModel, int * caseId);

unsigned long long ProcessReadDoubleStrand(SRAArguments * aArgs, SRAArguments * aArgs_neg, SRAModel * sraModel, SRAModel * sraModel_neg);

unsigned long long ProcessReadSingleStrand(SRAArguments * aArgs, SRAModel * sraModel);

void SRAPopulateMAPQParameters(SRAArguments * aArgsPos, SRAArguments * aArgsNeg, DPScores * dpScores, SRAModel * sraModels_Extend, SRAModel * sraModelsNeg_Extend);

#endif

