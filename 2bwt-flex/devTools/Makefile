####################################################################
#
#   Modification History
#
#   Date   : 19th February 2012
#   Author : Edward MK Wu
#   Change : New file.
#
####################################################################

CC = icpc
DEFINE = 
CFLAGS = -wd2568 -openmp -O3 -funroll-loops -I/home/rbluo/test/zlib/zlib-1.2.8/ -vec-report=3
LIBS = -lm -lz
THREADLIBS = -lpthread

UTILLIB = ../utilities
UTILOBJLIBS = $(UTILLIB)/BufferedFileReader.o

BWTLIB = ../2bwt-lib
BWTOBJLIBS = $(BWTLIB)/BWT.o $(BWTLIB)/dictionary.o $(BWTLIB)/blast_dust.o $(BWTLIB)/DNACount.o \
			$(BWTLIB)/HSP.o $(BWTLIB)/HSPstatistic.o $(BWTLIB)/iniparser.o $(BWTLIB)/inistrlib.o \
			$(BWTLIB)/karlin.o $(BWTLIB)/MemManager.o $(BWTLIB)/MiscUtilities.o $(BWTLIB)/QSufSort.o \
			$(BWTLIB)/r250.o $(BWTLIB)/TextConverter.o $(BWTLIB)/Timing.o $(BWTLIB)/Socket.o $(BWTLIB)/2BWT-Interface.o

PELIB = ../pair-end
PEOBJLIBS = $(PELIB)/PEAlgnmt.o

SAMLIB = ../samtools-0.1.19
SAMOBJLIBS = $(SAMLIB)/sam.o $(SAMLIB)/bam.o $(SAMLIB)/bgzf.o $(SAMLIB)/kstring.o $(SAMLIB)/bam_import.o \
				$(SAMLIB)/faidx.o $(SAMLIB)/bam_pileup.o $(SAMLIB)/bam_aux.o $(SAMLIB)/sam_header.o $(SAMLIB)/razf.o

SRALIB = ../
SRAOBJLIBS = $(SRALIB)/LT.o $(SRALIB)/HOCC.o $(SRALIB)/OutputAnyFormat.o $(SRALIB)/OutputBinaryFormat.o \
			$(SRALIB)/OutputSAMFormat.o $(SRALIB)/SRACore.o $(SRALIB)/SRAOccCollector.o $(SRALIB)/SRA2BWTReport.o \
			$(SRALIB)/SRA2BWTMdl.o $(SRALIB)/SRAConfigParser.o $(SRALIB)/SRA2BWTCheckAndExtend.o $(SRALIB)/SRA2BWTOperations.o \
			$(SRALIB)/SRAQueryParser.o $(SRALIB)/2BWT-SRAAlgnmt.o $(SRALIB)/DPCore.o $(SRALIB)/DPOperations.o $(SRALIB)/DPOccCollector.o \
			$(SRALIB)/SRADPReport.o $(SRALIB)/MICCore.o $(SRALIB)/SRAArguments.o $(SRALIB)/SRAOutputFile.o $(SRALIB)/PECore.o $(SRALIB)/PEArguments.o \
			$(SRALIB)/DPArguments.o $(SRALIB)/PEDPArguments.o $(SRALIB)/MappingQuality.o
PEAOBJLIBS = $(SRALIB)/2BWT-PEAlgnmt.o $(SRALIB)/PE2BWTReport.o $(SRALIB)/PEDPReport.o $(SRALIB)/ListConsumer.o

all:	DEV_STRUCT_SIZE DEV_TEST_DP DEV_PRINT_LT DEV_CFG_READER EXTRACT_READ EXTRACT_SEQ DEV_INDEX_UTILS DEV_LIST_CONS_UTILS DEV_NEG_BODY

$(BWTLIB) : force_look
	cd $(BWTLIB); $(MAKE)

$(PELIB) : force_look
	cd $(PELIB); $(MAKE)

$(SAMLIB) : force_look
	cd $(SAMLIB); $(MAKE)

# SOAP3-CPU-PE
PE2BWTReport.o:	            ../PE2BWTReport.c ../PE2BWTReport.h Makefile
DPOperations.o:				../DPOperations.c ../DPOperations.h Makefile
DPOccCollector.o:	        ../DPOccCollector.c ../DPOccCollector.h Makefile
DPCore.o:					../DPCore.c ../DPCore.h Makefile
PEDPReport.o:	            ../PEDPReport.c ../PEDPReport.h Makefile
devStructSize.o:	        devStructSize.c Makefile
DEV_STRUCT_SIZE: devStructSize.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) Makefile
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) devStructSize.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) -o devStructSize
	
devTestDp.o:				devTestDp.c Makefile
DEV_TEST_DP: 				devTestDp.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) Makefile
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) devTestDp.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) -o devTestDp
	
PrintBWT.o:					PrintBWT.c Makefile
DEV_PRINT_BWT:				PrintBWT.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) Makefile
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) PrintBWT.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) -o printBwt

PrintSA.o:					PrintSA.c Makefile
DEV_PRINT_SA:				PrintSA.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) Makefile
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) PrintSA.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) -o printSa

PrintLT.o:					PrintLT.c Makefile
DEV_PRINT_LT:				PrintLT.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) Makefile
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) PrintLT.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) -o printLt

cfgReader.o:				cfgReader.c Makefile
DEV_CFG_READER:				cfgReader.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) Makefile
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) cfgReader.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) -o cfgReader

ExtractReads.o:				ExtractReads.c Makefile
EXTRACT_READ:				ExtractReads.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) Makefile
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) ExtractReads.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) -o extractReads
	
devExtractSeq.o:				devExtractSeq.c Makefile
EXTRACT_SEQ:				devExtractSeq.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) Makefile
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) devExtractSeq.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) -o devExtractSeq
	
devIndexUtils.o:				devIndexUtils.c Makefile
DEV_INDEX_UTILS:			devIndexUtils.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) Makefile
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) devIndexUtils.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(UTILOBJLIBS) -o devIndexUtils

devTestListConsumer.o:	    devTestListConsumer.c Makefile
DEV_LIST_CONS_UTILS:			devTestListConsumer.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(PEAOBJLIBS) $(UTILOBJLIBS) Makefile
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) devTestListConsumer.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(PEAOBJLIBS) $(UTILOBJLIBS) -o devTestListConsumer

devNegSeq.o:	    devNegSeq.c Makefile
DEV_NEG_BODY:			devNegSeq.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(PEAOBJLIBS) $(UTILOBJLIBS) Makefile
	$(CC) $(CFLAGS) $(LIBS) $(THREADLIBS) devNegSeq.o $(SAMOBJLIBS) $(PEOBJLIBS) $(BWTOBJLIBS) $(SRAOBJLIBS) $(PEAOBJLIBS) $(UTILOBJLIBS) -o devNegSeq
