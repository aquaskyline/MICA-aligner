#!/bin/bash
#
# This script strip off all auxiliary data on a SAM output format and retain
# only the chromosome name and offset. The output will be printed to standard
# out and will be in the format of <SIDE>: <READ-NAME> <CHROME-NAME> <+/-><OFFSET>
# The usage of the scripts is ./stripSAM.sh <SAM file>
#

cat $1 | awk '{if (substr($1,0,1)!="@") { \
if (and($2,16)>0) {temp="-"} else {temp="+"}  \
if (NR%2==0) {side="READ"} else {side="MATE"} \
printf("%s: %s %s %s%s\n",side,$1,$3,temp,$4); \
secOccStr=substr($12,6,length($12)); \
secOccCount=split(secOccStr,secOcc,";"); \
for (i=1;i<=secOccCount;i++) { \
  if (secOcc[i]=="") continue; \
  split(secOcc[i],occElem,","); \
  printf("%s: %s %s %s\n",side,$1,occElem[1],occElem[2]); \
}\
}\
}'
