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
CFLAGS = -wd2568 -openmp -O3 -funroll-loops -I/home/rbluo/test/zlib/zlib-1.2.8/ -vec-report=3
LIBS = -L/miclib64 -lz
THREADLIBS = -lpthread

UTILLIB = ../2bwt-flex/utilities
UTILOBJLIBS = $(UTILLIB)/BufferedFileReader.o $(UTILLIB)/Logging.o $(UTILLIB)/SimpleQueue.o

BWTLIB = ../2bwt-flex/2bwt-lib
BWTOBJLIBS = $(BWTLIB)/BWT.o $(BWTLIB)/dictionary.o $(BWTLIB)/blast_dust.o $(BWTLIB)/DNACount.o \
			$(BWTLIB)/HSP.o $(BWTLIB)/HSPstatistic.o $(BWTLIB)/iniparser.o $(BWTLIB)/inistrlib.o \
			$(BWTLIB)/karlin.o $(BWTLIB)/MemManager.o $(BWTLIB)/MiscUtilities.o $(BWTLIB)/QSufSort.o \
			$(BWTLIB)/r250.o $(BWTLIB)/TextConverter.o $(BWTLIB)/Timing.o $(BWTLIB)/Socket.o $(BWTLIB)/2BWT-Interface.o

PELIB = ../2bwt-flex/pair-end
PEOBJLIBS = $(PELIB)/PEAlgnmt.o

SAMLIB = ../2bwt-flex/samtools-0.1.19
SAMOBJLIBS = $(SAMLIB)/sam.o $(SAMLIB)/bam.o $(SAMLIB)/bgzf.o $(SAMLIB)/kstring.o $(SAMLIB)/bam_import.o \
				$(SAMLIB)/faidx.o $(SAMLIB)/bam_pileup.o $(SAMLIB)/bam_aux.o $(SAMLIB)/sam_header.o $(SAMLIB)/razf.o

SRALIB = ../2bwt-flex/
SRAOBJLIBS = $(SRALIB)/LT.o $(SRALIB)/HOCC.o $(SRALIB)/OutputAnyFormat.o $(SRALIB)/OutputBinaryFormat.o \
			$(SRALIB)/OutputSAMFormat.o $(SRALIB)/SRACore.o $(SRALIB)/SRAOccCollector.o $(SRALIB)/SRA2BWTReport.o \
			$(SRALIB)/SRA2BWTMdl.o $(SRALIB)/SRAConfigParser.o $(SRALIB)/SRA2BWTCheckAndExtend.o $(SRALIB)/SRA2BWTOperations.o \
			$(SRALIB)/SRAQueryParser.o $(SRALIB)/2BWT-SRAAlgnmt.o $(SRALIB)/DPCore.o $(SRALIB)/DPOperations.o $(SRALIB)/DPOccCollector.o \
			$(SRALIB)/SRADPReport.o $(SRALIB)/MICCore.o $(SRALIB)/SRAArguments.o $(SRALIB)/SRAOutputFile.o $(SRALIB)/PECore.o $(SRALIB)/PEArguments.o \
			$(SRALIB)/DPArguments.o $(SRALIB)/PEDPArguments.o $(SRALIB)/MappingQuality.o

PEAOBJLIBS = $(SRALIB)/2BWT-PEAlgnmt.o $(SRALIB)/PE2BWTReport.o $(SRALIB)/PEDPReport.o \
		$(SRALIB)/ListConsumer.o

MICLIB = ../
MICOBJLIBS = $(MICLIB)/MIC-SRA2BWTMdl.o $(MICLIB)/MIC-SRAArguments.o $(MICLIB)/MIC-SRA2BWTOperations.o \
			$(MICLIB)/MIC-SRAAlgnmt.o $(MICLIB)/MIC-PEAlgnmt.o $(MICLIB)/MIC-DPCore.o \
			$(MICLIB)/MIC-DPAlgnmt.o $(MICLIB)/MICDeepDP.o $(MICLIB)/MicMappingQuality.o \
			$(MICLIB)/MICA-PE-ReadThread.o $(MICLIB)/MICA-PE-WriteThread.o $(MICLIB)/MemMan.o \
			$(MICLIB)/MICControllerThread.o $(MICLIB)/CPUControllerThread.o \
			$(MICLIB)/MICCT-ShortHandler.o \
			$(MICLIB)/MIC-MEMControl.o

all:	testWriteThread samFlagTranslator testOffload testOffloadMalloc

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
$(MICLIB) : force_look
	cd $(MICLIB); $(MAKE)

MICA_PE_OBJECTS = $(MICOBJLIBS) $(SAMOBJLIBS) $(PEOBJLIBS) \
	$(BWTOBJLIBS) $(SRAOBJLIBS) $(PEAOBJLIBS) $(UTILOBJLIBS)
testWriteThread: testWriteThread.o $(MICA_PE_OBJECTS)
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) testWriteThread.o $(MICA_PE_OBJECTS) -o testWriteThread
samFlagTranslator: samFlagTranslator.o $(MICA_PE_OBJECTS)
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) samFlagTranslator.o $(MICA_PE_OBJECTS) -o samFlagTranslator
	
testOffload: $(MICA_PE_OBJECTS) testOffload.o
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) $(MICA_PE_OBJECTS) testOffload.o -o testOffload

testOffloadMalloc: $(MICA_PE_OBJECTS) testOffloadMalloc.o
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) $(MICA_PE_OBJECTS) testOffloadMalloc.o -o testOffloadMalloc
