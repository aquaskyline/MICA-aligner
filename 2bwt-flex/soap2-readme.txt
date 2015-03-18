soap2-dp - README
Version 2.5.265
==================

Introduction
------------
SOAP2 - is a short read alignment software that is based on the 2BWT index.
The soap2-dp extends the functionality of SOAP2 to support alignment
with up to 5 mismatches or 5 insert/delete of single character. Furthermore, 
multi-threading is also supported in soap2-dp to make full use of the 
multi-core CPU to speed up the alignment process.

soap2-dp supports both single-end and pair-end alignment. 

The alignment software in this package handles reference sequences 
in the commonly known FASTA format and short reads in FASTA/FASTQ format.
The alignment output is in humanly readable manner.

soap2-dp normally operates with 16 GB main memory. 

The 2BWT index and search algorithms used by soap2-dp were developed by the 
algorithms reserach group of the University of Hong Kong (T.W. Lam, C.M. Liu
Alan Tam, Simon Wong, Thomas Wong, Edward Wu and S.M. Yiu).


Hardware Requirement
--------------------
Recommended system requirement
- A quad core CPU
- 16GB main memory
- 64-bit environment


Installation and Configuration
------------------------------
soap2-dp can be downloaded from this website : http://www.cs.hku.hk/2bwt-tools/

Follow the following steps to install the software after you have 
obtained the software archive.

1. Move the archive to the directory you wish to install it.

2. Unzip the archive with gunzip and tar
    % gunzip soap2-dp-v2.5.<rev>-x86-64bit.tar.gz
    % tar -xvvf soap2-dp-v2.5.<rev>-x86-64bit.tar
    % cd soap2-dp-v2.5.<rev>\

3. Inside the archive you should find 11 files extracted.
   They are:
   Index building software:                      soap2-dp-builder soap2-dp-builder.ini
   Alignment software:                           soap2-dp soap2-dp.ini
   Core alignment modules (advanced users only): soap2-cpu-se soap2-cpu-pe
   Bespoke Binary Output viewing software:       soap2-dp-viewer

4. For the usages of these software please refer to the next section.

Upgrade Instruction
-------------------
Upgrading into this version of software requires the following actions:
 - 05/09/2012 Upgrading from revision earlier than 218,
    There is an update with the translation mechanism such that pre-processed chromosome
    lengths can be outputted correctly. Previously, the SAM output contains chromosome 
    lengths that are only good as guess. To take advantage of the correct chromosome
    lengths in SAM output, the translation index needs to be re-built. To do that,
    one may set all construction step to N, except ParseFASTA remains as Y, in
    soap2-dp-builder.ini and re-run soap2-dp-builder against the reference sequence.

Software Usage Guide - Index Builder
----------------------------------------------
Index builder is the essential first step to preprocess a reference
sequence (FASTA format) for the alignment software. It
generates a 2BWT Index.

    > Limitations
    -------------
    Please take note of the following restriction on the index builder
    1.  Builder assumes the input is a well formatted FASTA file.
    2.  There can be multiple sequences within the FASTA file; however,
        the total length of the reference sequences cannot exceed 2^32 - 1.
    3.  No more than 254 sequences in the FASTA file.
    4.  Characters other than A, C, G and T will be considered invalid characters.
        Any segment of more than 10 consecutive invalid characters will be removed.
        All the rest of the invalid characters will be replaced by "G".

To invoke the builder, follow the step:

    % ./soap2-dp-builder <FASTA sequence file>

E.g.,
    % ./soap2-dp-builder ncbi.genome.fa

In this example,soap2-dp-builder will build up a 2BWT index for the sequence 
ncbi.genome.fa. The output index will be named ncbi.genome.fa.index.*. 
There should be 16 files with the same prefix created.

Software Usage Guide - Short Read Alignment
-------------------------------------------------------
Assume that the 2BWT index of a reference sequence has been built.
The alignment software takes in the 2BWT index and a FASTA/FASTQ
file contains any number of short reads for alignment.
There are several options to control the alignment process.

    > Limitations
    -------------
    1.  It finds alignment with up to 5 mismatches or 5 insert/delete of single character.
    2.  It supports reads of any length between 20-200.
    3.  It supports reads contains only A, C, G, T.

Please invoke the alignment software as follows:

    % ./soap2-dp single <index>.index <FASTA/FASTQ read file> [option] [option] [option] ...

<index> is the filename of the FASTA reference sequence file you built the index on
<FASTA/FASTQ read file> is the filename of the short reads in FASTA/FASTQ format

    [Options]
        -m: Maximum #mismatch allowed. [1-5;Def=Auto]
        -g: Maximum #edit allowed.  [1-5;Def=Auto]
            If none of above is supplied, aligner will auto-default
            the alignment criteria by inspecting the read files.
            See Readme for more information.
        -S: 1 - Upper; 2 - Lower; 3 - Both (default)
        -h: Alignment type. [Def=4]
            1 : All Valid Alignment
            2 : All Best Alignment
            3 : Unique Best Alignment
            4 : Random Best Alignment
            5 : Best Quality Alignment
        -b: Output format. [Def=1]
            0 : Bespoke (Binary)
            1 : Default (Plain)
            2 : SAM 1.4 (Plain)

Examples:

    % ./soap2-dp single ncbi.genome.fa.index reads.fa -m 2 -h 1

In this example, it will report all-valid alignment of reads.fa 
with at most 2-mismatches. See Output Format section for description of
the output.

    % ./soap2-dp single ncbi.genome.fa.index reads.fa -g 4 -h 1
In this example, it will report all-valid alignment of reads.fa
with at most 4-indels.

    > Default Parameters for Alignment
    ----------------------------
    If neither -m or -g is supplied, the alignment parameters will be
    defaulted. Aligner will inspect the first 10 reads and
    determine the parameter with their max lengths. If the max lengths
    is smaller than 50bp then the following default will apply:
        -m 2
    Otherwise, the following default will apply:
        -m 3

    > Multi-threading Support
    -------------------------
    Aligner supports spawning alignment threads to speed up the alignment
    when using a multi-core CPU. Inside soap2-dp.ini, 
    the parameter NumOfCPUThreads controls the number CPU threads to be used.

    [MultipleThreading]
    NumOfCPUThreads = 1;

    By setting this value to k, aligner will spawns k threads from the main program to
    perform alignment. When this value is set to k, the alignment results will be written
    onto k separate files.
    reads.fa.out.0, reads.fa.out.1, reads.fa.out.2, ... , reads.fa.out.k
    For simple text format or SAM v1.4 output format, each of these result files are read
    individually. For bespoked binary format, these results file are not humanly readable and
    required additional processing. Please refer to the section -
    "SOAP2 Bespoke Binary Output Format".

Pre-Alignment Trimming
----------------------
The aligner supports trimming the read by a percentage before feeding the reads into
the alignment program. Please read soap2-dp.ini for guide on how to configure
it. This is by-default disabled.

Post-Alignment Trimming and 2nd-Pass Alignment
----------------------------------------------
The aligner supports trimming each read by a percentage after the alignment programm cannot
find any possible alignment result for the same read in the first round of alignment.
The read will then be trimmed and passed into the second round of alignment. The alignment program also
supports modifiying the alignment parameters of the second round alignment, e.g. loosen the #of mismatch.
Please read soap2-dp.ini for guide on how to configure it. This is by-default disabled.


Software Usage Guide - Pair-End Alignment
----------------------------------------------------
The alignment software takes in the 2BWT index and two FASTA/FASTQ
file contains equal number of short reads for pair-end alignment.
There are several options to control the alignment process.

    > Limitations
    -------------
    1.  It finds alignment with up to 5 mismatches or 5 insert/delete of single character.
    2.  It supports reads of any length between 20-200.
    3.  It supports reads contains only A, C, G, T.

Please invoke the alignment software as follows:

    % ./soap2-dp_ pair <index>.index <FASTA/FASTQ read file 1> <FASTA/FASTQ read file 2>  -v <l-insert> -u <u-insert> [option] ...

<index> is the filename of the FASTA reference sequence file you built the index on
<FASTA/FASTQ read file 1/2> is the filenames of the short reads in FASTA/FASTQ format

    [Mandatory]
        -v: Lower bound;
        -u: Upper bound of the insertion size between a pair.
        A pair will be reported if the insertion size falls
        between [v,u] and their strands match.

    [Options]
        -m: Maximum #mismatch allowed. [1-5;Def=Auto]
        -g: Maximum #edit allowed.  [1-5;Def=Auto]
        If none of above is supplied, aligner will auto-default
        the alignment criteria by inspecting the read files.
        See Readme for more information.
        -h: Alignment type. [Def=4]
            1 : All Valid Pair-end Alignment
            4 : Arbitrary Best Pair-end Alignment
        -b: Output format. [Def=1]
            0 : Bespoke (Binary)
            1 : Default (Plain)
            2 : SAM 1.4 (Plain)


Examples:

    % ./soap2-dp pair ncbi.genome.fa.index reads_A.fa reads_B.fa -g 4 -h 1 -v 0 -u 3000
In this example, it will report all-valid pair-end alignment of reads_A.fa 
nd reads_B.fa with at most 4-indels; tolerating insertion size of [0,3000].

    > Default Parameters for Alignment
    ----------------------------
    Similar to Single-End Alignment, yet max read length is defined by the max
    among the first 10 reads of both read files.

    > Multi-threading Support
    -------------------------
    Similar to Single-End Alignment


Output Format
----------------------------------------------
soap2-dp supports up to three output formats.
(0): Bespoke format for high performance (Binary)
(1): Simple text format (Plain)
(2): SAM v1.4

Bespoke Binary 
----------------------
In order to maintain reasonable performance upon huge dataset, we have prepared a binary formatted 
output to minimize possible system I/O. Aligner can be configured to report alignment positions 
into its bespoke binary format. 
(In Pair-End Alignment, two consecutive entries form one pair-end alignment in this format.)

In order to convert the binary file into human readable format. We need to parse the output
with the output viewer. The human readable output will be sent to standard output by the
viewer. If you need it into a file, then you will have to pipe it to a file.

Please invoke the output viewer by following the below guide.
    % ./soap2-dp-viewer <FASTA/FASTQ read file> <output file>

e.g.
    % ./soap2-dp-viewer reads.fa reads.fa.out >  reads.fa.out.txt
In this example, the output files of a 4-threads alignment will be output to a file.

Simple Text Format
----------------------
This output format target for a concise output of each alignment result. Except for the header lines,
each line in the output file corresponds to one alignment output. Read Multi-threading section for 
detail of the output files.
(In Pair-End Alignment, two consecutive entries form one pair-end alignment in this format.)

SAM v1.4 by SAM-tools
----------------------
SAM-tools v0.1.18 is included in soap2-dp package to faciliate outputting alignment result 
into SAM output format. We have slightly modified the original code of SAM-tools to make it 
compilable under g++. Please see http://samtools.sourceforge.net/ for more detail of this package.
Also read Multi-threading section for detail of the output files.

Reference
----------------------------------------------
T.W. Lam, R. Li, Alan Tam, Simon Wong, Edward Wu, S.M. Yiu: High Throughput Short Read Alignment via Bi-directional BWT. BIBM 2009: 31-36


