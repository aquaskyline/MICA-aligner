//
//    Logging.c
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

#include "Logging.h"

Logging * LOGCreate( char * logFilename, int mode ) {
    Logging * logging = (Logging*) malloc ( sizeof(Logging) );
    
    if ( mode == LOGGING_MODE_RUNNING ) {
        logging->logFp = (FILE*)(FILE*)fopen64(logFilename,"a");
    } else if ( mode == LOGGING_MODE_STANDALONE ) {
        logging->logFp = (FILE*)(FILE*)fopen64(logFilename,"w");
    } else {
        printf("[SAFE-GUARD] Unknown logging mode %d\n",mode);
        logging->logFp = NULL;
    }
    
    if ( logging->logFp == NULL ) {
        printf("[ERROR] Unable to open file %s for log writing.\n",logFilename);
        printf("Log messages will be printed to screen instead.\n");
    }
    
    return logging;
}

void LOGFree( Logging * logging ) {
    fclose(logging->logFp);
    logging->logFp = NULL;
    free(logging);
}

////////////////////////////////////////////////////
// Log Writer
////////////////////////////////////////////////////

int LOGWriteLine( Logging * logging, int msgType, const char *format, ...) {
    va_list args;
    
    // Print to screen if applicable
    va_start(args, format);
    if( msgType == LOGGING_MSG_TYPE_ALL || 
        msgType == LOGGING_MSG_TYPE_SCREEN_ONLY) {
        vprintf(format, args);
        fflush(stdout);
    }
    va_end(args);
    
    // Print to file if applicable
    va_start(args, format);
    if( msgType == LOGGING_MSG_TYPE_ALL || 
        msgType == LOGGING_MSG_TYPE_LOG_ONLY) {
        if (logging->logFp != NULL) {
            vfprintf(logging->logFp, format, args);
            fflush(logging->logFp);
        }
    }
    va_end(args);
    
    return 0;
}

