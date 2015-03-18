#
#    Makefile
#
#    SOAP2 / 2BWT
#
#    Copyright (C) 2012, HKU
#
#    This program is free software; you can redistribute it and/or
#    modify it under the terms of the GNU General Public License
#    as published by the Free Software Foundation; either version 2
#    of the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
####################################################################
#
#   Modification History
#
#   Date   : 17th May 2013
#   Author : Edward MK Wu
#   Change : MICA New Makefile
#
####################################################################

CC = icc
DEFINE = 
CFLAGS = -wd2568 -openmp -O3 -funroll-loops -I../thirdparty/zlib-1.2.8/ -vec-report=3
LIBS = -L/miclib64 -lz
THREADLIBS = -lpthread

UTILLIB = 2bwt-flex/utilities
UTILOBJLIBS = $(UTILLIB)/BufferedFileReader.o $(UTILLIB)/Logging.o $(UTILLIB)/SimpleQueue.o

BWTLIB = 2bwt-flex/2bwt-lib
BWTOBJLIBS = $(BWTLIB)/BWT.o $(BWTLIB)/dictionary.o $(BWTLIB)/blast_dust.o $(BWTLIB)/DNACount.o \
			$(BWTLIB)/HSP.o $(BWTLIB)/HSPstatistic.o $(BWTLIB)/iniparser.o $(BWTLIB)/inistrlib.o \
			$(BWTLIB)/karlin.o $(BWTLIB)/MemManager.o $(BWTLIB)/MiscUtilities.o $(BWTLIB)/QSufSort.o \
			$(BWTLIB)/r250.o $(BWTLIB)/TextConverter.o $(BWTLIB)/Timing.o $(BWTLIB)/Socket.o $(BWTLIB)/2BWT-Interface.o

PELIB = 2bwt-flex/pair-end
PEOBJLIBS = $(PELIB)/PEAlgnmt.o

SAMLIB = 2bwt-flex/samtools-0.1.19
SAMOBJLIBS = $(SAMLIB)/sam.o $(SAMLIB)/bam.o $(SAMLIB)/bgzf.o $(SAMLIB)/kstring.o $(SAMLIB)/bam_import.o \
				$(SAMLIB)/faidx.o $(SAMLIB)/bam_pileup.o $(SAMLIB)/bam_aux.o $(SAMLIB)/sam_header.o $(SAMLIB)/razf.o

SRALIB = 2bwt-flex/
SRAOBJLIBS = $(SRALIB)/LT.o $(SRALIB)/HOCC.o $(SRALIB)/OutputAnyFormat.o $(SRALIB)/OutputBinaryFormat.o \
			$(SRALIB)/OutputSAMFormat.o $(SRALIB)/SRACore.o $(SRALIB)/SRAOccCollector.o $(SRALIB)/SRA2BWTReport.o \
			$(SRALIB)/SRA2BWTMdl.o $(SRALIB)/SRAConfigParser.o $(SRALIB)/SRA2BWTCheckAndExtend.o $(SRALIB)/SRA2BWTOperations.o \
			$(SRALIB)/SRAQueryParser.o $(SRALIB)/2BWT-SRAAlgnmt.o $(SRALIB)/DPCore.o $(SRALIB)/DPOperations.o $(SRALIB)/DPOccCollector.o \
			$(SRALIB)/SRADPReport.o $(SRALIB)/MICCore.o $(SRALIB)/SRAArguments.o $(SRALIB)/SRAOutputFile.o $(SRALIB)/PECore.o $(SRALIB)/PEArguments.o \
			$(SRALIB)/DPArguments.o $(SRALIB)/PEDPArguments.o $(SRALIB)/MappingQuality.o
PEAOBJLIBS = $(SRALIB)/2BWT-PEAlgnmt.o $(SRALIB)/PE2BWTReport.o $(SRALIB)/PEDPReport.o \
            $(SRALIB)/ListConsumer.o

all:	MICA-SE MICA-PE SOAP2-Builder index-builder-package

clean:
	rm $(UTILLIB)/*.o $(BWTLIB)/*.o $(PELIB)/*.o $(SRALIB)/*.o *.o

fastclean:
	rm *.o

$(UTILLIB) : force_look
	cd $(UTILLIB); $(MAKE)
$(SAMLIB) : force_look
	cd $(SAMLIB); $(MAKE)
$(BWTLIB) : force_look
	cd $(BWTLIB); $(MAKE)
$(PELIB) : force_look
	cd $(PELIB); $(MAKE)
$(SRALIB) : force_look
	cd $(SRALIB); $(MAKE)

# MIC Library

# MICA-SE
MICA-SE.o:				MICA-SE.c MICA-SE.h Makefile
MIC-SRA2BWTMdl.o:			MIC-SRA2BWTMdl.c MIC-SRA2BWTMdl.h Makefile
MIC-SRAArguments.o:			MIC-SRAArguments.c MIC-SRAArguments.h Makefile
MIC-SRA2BWTOperations.o:		MIC-SRA2BWTOperations.c MIC-SRA2BWTOperations.h Makefile
MIC-SRAAlgnmt.o:		MIC-SRAAlgnmt.c MIC-SRAAlgnmt.h Makefile
MIC-PEAlgnmt.o:	MIC-PEAlgnmt.c MIC-PEAlgnmt.h
MIC-DPCore.o: MIC-DPCore.h MIC-DPCore.c
MIC-DPAlgnmt.o: MIC-DPAlgnmt.h MIC-DPAlgnmt.c
MicMappingQuality.o: MicMappingQuality.h MicMappingQuality.c
MICDeepDP.o:		MICDeepDP.c MICDeepDP.h
MemMan.o: MemMan.c MemMan.h

# MICA-PE
MICA-PE-ReadThread.o: MICA-PE-ReadThread.c MICA-PE-ReadThread.h
MICA-PE-WriteThread.o: MICA-PE-WriteThread.c MICA-PE-WriteThread.h
MICControllerThread.o: MICControllerThread.c MICControllerThread.h 
CPUControllerThread.o: CPUControllerThread.c CPUControllerThread.h 
MICCT-ShortHandler.o: MICCT-ShortHandler.c MICCT-ShortHandler.h 
MIC-MEMControl.o: MIC-MEMControl.c MIC-MEMControl.h
MICA-Health.o: MICA-Health.c MICA-Health.h 

MICA_SE_OBJECTS = MIC-SRA2BWTMdl.o MIC-SRAArguments.o MIC-SRA2BWTOperations.o \
			MIC-SRAAlgnmt.o MemMan.o MICA-PE-ReadThread.o \
			MicMappingQuality.o \
			$(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(PEAOBJLIBS) $(UTILOBJLIBS)
MICA-SE: MICA-SE.o $(MICA_SE_OBJECTS) Makefile
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) MICA-SE.o $(MICA_SE_OBJECTS) -o mica-se


MICA_PE_OBJECTS = MIC-PEAlgnmt.o MIC-DPCore.o MIC-DPAlgnmt.o MICDeepDP.o \
	MICA-Health.o \
	MICA-PE-WriteThread.o \
	MICControllerThread.o CPUControllerThread.o \
	MICCT-ShortHandler.o \
	MIC-MEMControl.o \
	$(MICA_SE_OBJECTS)
MICA-PE: MICA-PE.o $(MICA_PE_OBJECTS) Makefile
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) MICA-PE.o $(MICA_PE_OBJECTS) -o mica-pe

SOAP2-Builder:	$(SRALIB)/2BWT-Builder.o $(BWTOBJLIBS) $(BWTLIB)/BWTConstruct.o $(SRALIB)/LTConstruct.o $(SRALIB)/HOCCConstruct.o Makefile
	$(CC) $(CFLAGS) $(LIBS) $(SRALIB)/2BWT-Builder.o $(BWTOBJLIBS) $(BWTLIB)/BWTConstruct.o $(SRALIB)/LTConstruct.o $(SRALIB)/HOCCConstruct.o -o soap2-dp-builder
index-builder-package: SOAP2-Builder
	cp soap2-dp-builder index-builder-step1;
	cp soap2-dp-builder index-builder-step2;
	
