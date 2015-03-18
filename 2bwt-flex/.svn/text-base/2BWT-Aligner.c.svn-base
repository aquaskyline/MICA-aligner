//
//    2BWT-Aligner.c
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

    2BWT-Aligner.c
    
    This application is only a wrapper class that cakls 
    the corresponding 2BWT/SOAP2 module.
    
    Modification History 

    Date   : 28th December 2011
    Author : Edward MK Wu
    Change : New file.

*/
/////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "Release.h"

#define WRAPPER_MAX_ARGUMENT 256

#define WRAPPER_SINGLE_END_ALIGNMENT_APP    "soap2-cpu-se"
#define WRAPPER_PAIR_END_ALIGNMENT_APP      "soap2-cpu-pe"

void PrintUsage(const char * arg_0) {
    printf("\n");
    printf("Usage Guide:\n");
    printf("    This software performs single-end and pair-end alignment\n");
    printf("    with pre-processed FASTA reference sequence and FASTA/FASTQ short queries.\n");
    printf("    For a full how-to reference please refer to the readme document.\n");
    printf("\n");
    printf("For usage guide of Single-End alignment run below:\n");
    printf("    %s single\n",arg_0);
    printf("\n");
    printf("For usage guide of Pair-End alignment run below:\n");
    printf("    %s pair\n",arg_0);
    printf("\n");
}


int main (int argc, char ** argv) {

    printf("\n[Main] %s v%d.%d.%d (%s)\n",PROJECT_NAME,PROJECT_MAJOR,PROJECT_MINOR,PROJECT_REV,PROJECT_SPECIAL);

    if (argc<=1) {
        PrintUsage(argv[0]);
        return 1;
    }
    
    time_t now;
    time(&now);
    printf("[Main] invoked at %s", ctime(&now));
    
    char * tmpArgv[WRAPPER_MAX_ARGUMENT];
    int i;
    printf("[Main] with arguments -- %s ",argv[0]);
    for (i=1;i<argc;i++) {
        tmpArgv[i] = argv[i];
        printf("%s ",argv[i]);
    }
    printf("\n\n");
    tmpArgv[i] = NULL;

    ////////////////////////////////////////////////////////////
    // Wrapper program will now replace itself with the 
    // corresponding child program -- by execv.
    ////////////////////////////////////////////////////////////
    if (strcmp(argv[1],"single")==0) {
        tmpArgv[0] = argv[0];
        execv(WRAPPER_SINGLE_END_ALIGNMENT_APP, tmpArgv);
    } else if (strcmp(argv[1],"pair")==0) {
        tmpArgv[0] = argv[0];
        execv(WRAPPER_PAIR_END_ALIGNMENT_APP, tmpArgv);
    } else {
        PrintUsage(argv[0]);
        return 1;
    }
    ////////////////////////////////////////////////////////////
    // Anything below this line is unreachable in any situation
    ////////////////////////////////////////////////////////////
    return 0;
}
