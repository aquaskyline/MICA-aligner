//
//    devNegSeq.c
//
//    soap2-dp
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

#include "stdio.h"
#include "stdlib.h"
#include "../2bwt-lib/HSP.h"


void SRAFillCharMap(unsigned char * charMap) {
    int i;
    for (i=0;i<256;i++) { charMap[i] = 2; };
    for (i=0;i<ALPHABET_SIZE;i++) {
        charMap[dnaChar[i]] = i;
    }
}


int main (int argc, char ** argv) {

    unsigned char tempChar[1024];
    unsigned char charMap[256];
    char c;
    
    SRAFillCharMap(charMap);
    
    printf("Read Body = ");
    fflush(stdout);
    c = '0';
    int i =0;
    while(c!='\n') {
        c = fgetc(stdin);
        tempChar[i++] = charMap[c];
    }
    tempChar[--i] = '\0';
    
    int j;
    printf("Pos Body = ");
    for (j=0;j<i;j++) {
        printf("%c",dnaChar[tempChar[j]]);
    }
    printf("\n");
    printf("Neg Body = ");
    for (j=i-1;j>=0;j--) {
        int nc = 3 - tempChar[j];
        printf("%c",dnaChar[nc]);
    }
    printf("\n");
}
