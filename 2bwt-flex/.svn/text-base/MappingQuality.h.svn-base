/*
 *
 *    MappingQuality.h
 *    Soap3(gpu)
 *
 *    Copyright (C) 2011, HKU
 *
 *    This program is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License
 *    as published by the Free Software Foundation; either version 2
 *    of the License, or (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 */

#ifndef __MAPPING_QUALITY_H__
#define __MAPPING_QUALITY_H__

#include "DPCore.h"

__attribute__((target(mic)))
void bwaLikePairQualScore (
	int x0_0, int x1_0, int x0_1, int x1_1, int * g_log_n,
	int op_score, int op_num, int subop_score, int subop_num,
	int readlen_0, int readlen_1, int * map_score0, int * map_score1);

__attribute__((target(mic)))
void bwase_initialize ( int * g_log_n );
int MAPQNormaliseSRAScore ( int readLength, int numMismatch, DPScores * dpScores );

#endif
