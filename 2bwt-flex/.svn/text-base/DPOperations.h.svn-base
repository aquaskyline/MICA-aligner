//
//    DPOperations.h
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
   
   Date   : 3rd February 2012
   Author : Edward MK Wu
   Change : New file.

*/
/////////////////////////////////////////////////////

#ifndef __PE_DP_OPERATIONS_H__
#define __PE_DP_OPERATIONS_H__

#include "DPArguments.h"

void DPMatrixFetching(DPWork * dpWork, HSP * hsp);

// Assumption 1 -- Coordinates on Matching Matrix
// The operations assumes the first operand to coordinate parameter be
// reference on region on reference sequence and the second be on read.
//
// Assumption 2 -- DPMatrixCellScore called by DPMatrixFill
// DPMatrixCellScore computes the value of one cell 
// assuming the "base" cells are correctly computed.
// _____________
// |  A  |  B  |
// |_____|_____| i.e. when DPMatrixCellScore is called for X, 
// |  C  |  X  | we assume the value of A, B, C are ready.
// |_____|_____| 
//
// DPMatrixFill calls DPMatrixCellScore to fill up the entire matching matrix
int DPMatrixFill(DPArguments * dpArgs);

// DPBackTrack output the alignments from a completed matching matrix
int DPBackTrack(DPWork * dpWork, DPScores * dpScores, int minRecoverScore, DPOccurrence * dpOcc);

void DPDebugPrintMatrix(DPWork * dpWork);

#endif
