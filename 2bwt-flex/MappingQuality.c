/*
 *
 *    MappingQuality.c
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
 #include <math.h>

 #include "MappingQuality.h"

__attribute__((target(mic)))
int bwaLikeSingleQualScore ( int x0, int x1, int * g_log_n )
{
    // fprintf(stderr, "x0 : %i; x1 : %i\n", x0, x1);
    int score;

    if ( x0 > 1 )
    { score = 0; }
    else if ( x1 == 0 )
    { score = 37; }
    else
    {
        int cappedX1 = ( x1 > 255 ) ? 255 : x1;
        int n = g_log_n[cappedX1];
        score = ( 23 < n ) ? 0 : 23 - n;
    }

    return score;
}


__attribute__((target(mic)))
void bwaLikePairQualScore ( int x0_0, int x1_0, int x0_1, int x1_1, int * g_log_n, int op_score, int op_num, int subop_score, int subop_num, int readlen_0, int readlen_1, int * map_score0, int * map_score1 )
{
    // fprintf(stderr, "x0_0:%i x1_0:%i x0_1:%i x1_1:%i op_score:%i op_num:%i subop_score:%i subop_num:%i readlen_0:%i readlen_1:%i \n", x0_0, x1_0, x0_1, x1_1, op_score, op_num, subop_score, subop_num, readlen_0, readlen_1);
    int mapq0 = bwaLikeSingleQualScore ( x0_0, x1_0, g_log_n );
    int mapq1 = bwaLikeSingleQualScore ( x0_1, x1_1, g_log_n );
    // fprintf(stderr, "mapq0 : %i ; mapq1 : %i\n", mapq0, mapq1);
    op_score = op_score * 10;
    subop_score = subop_score * 10;
    int mapq_p = 0;

    if ( mapq0 > 0 && mapq1 > 0 )
    {
        mapq_p = mapq0 + mapq1;

        if ( mapq_p > 60 ) { mapq_p = 60; }

        mapq0 = mapq_p;
        mapq1 = mapq_p;
    }
    else
    {
        if ( op_num == 1 )
        {
            if ( subop_num == 0 )
            { mapq_p = 29; }
            else if ( op_score - subop_score > ( 0.3 * ( ( readlen_0 + readlen_1 ) / 2 ) ) ) { mapq_p = 23; }
            else
            {
                subop_num = ( subop_num > 255 ) ? 255 : subop_num;
                mapq_p = ( op_score - subop_score ) / 2 - g_log_n[subop_num];

                if ( mapq_p < 0 ) { mapq_p = 0; }
            }
        }

        if ( mapq0 == 0 )
        {
            mapq0 = ( mapq_p + 7 < mapq1 ) ? mapq_p + 7 : mapq1;
        }

        if ( mapq1 == 0 )
        {
            mapq1 = ( mapq_p + 7 < mapq0 ) ? mapq_p + 7 : mapq0;
        }
    }

    ( *map_score0 ) = mapq0;
    ( *map_score1 ) = mapq1;
}


__attribute__((target(mic)))
void bwase_initialize ( int * g_log_n )
{
    int i;

    for ( i = 1; i != 256; ++i ) { g_log_n[i] = ( int ) ( 4.343 * log ( i ) + 0.5 ); }
}

int MAPQNormaliseSRAScore ( int readLength, int numMismatch, DPScores * dpScores ) {
    return dpScores->dpMismatch * numMismatch + dpScores->dpMatch * ( readLength - numMismatch );
}
