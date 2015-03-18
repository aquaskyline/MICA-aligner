//
//    PEAlgnmt.h
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
/*

    Modification History 

    Date  : 22nd June 2011
    Author: Edward MK Wu
    Change: New file.

            Each pair-end result is considered to have two legs
            which each of them is a position on the reference sequence.
            
            The left leg is the position of the left-most position of
            the match (which must be a %strandLeftLeg% strand alignment) ;
            the right leg is the position of its mate 
            (which must be a %strandRightLeg% strand alignment);
            
                | left-leg                  | right-leg
                v                           v
            ___|XXXXXXXXXXXXX|_____________|XXXXXXXXXXXXX|_

                |<----------- insertion size ---------->|
*/
/////////////////////////////////////////////////////

#ifndef __PE_ALIGNMENT_H__
#define __PE_ALIGNMENT_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "../2bwt-lib/BWT.h"
#include "../2bwt-lib/HSP.h"
#include "../SRACore.h"
#include "../DPCore.h"

#define PE_MAX_BUCKET_SIZE 1024

#define PE_REPORT_RANDOM_BEST     4
#define PE_REPORT_ONE             7
#define PE_REPORT_ALL             1
#define PE_REPORT_ALL_BEST        2

#define PE_BEST_NOT_AVAILABLE     0
#define PE_BEST_AVAILABLE         1

#define PE_ALIGNMENT_COMPLETED    0
#define PE_ALIGNMENT_INITIALISED  1
#define PE_ALIGNMENT_INPUT_ERROR  2
#define PE_ALIGNMENT_FAILED       3

#define PE_INPUTOCC_SORTED        0
#define PE_INPUTOCC_UNSORTED      1

#define PE_NON_INPUT_VALUE 99999999

#define PE_MAX_NUM_OF_ERROR       2*MAX_NUM_OF_ERROR+1

/////////////////////////////////////////////////////
/*
    Structure SRAOccurrence / PESRAAlignmentResult (Input)
    They are the structures to hold the input to Pair-End alignment.
    
    SRAOccurrence is supposed to hold a list of occurrences.
    The function PEMappingOccurrences() takes two lists of SRAOccurrence as input.
    
    PESRAAlignmentResult is supposed to hold a list of SA ranges.
    The function PEMapping() takes two lists of PESRAAlignmentResult as input.
*/
/////////////////////////////////////////////////////    

typedef struct PESRAAlignmentResult {
    unsigned long long saIndexLeft;
    unsigned long long saIndexRight;
    uint8_t strand;
    uint8_t mismatchCount;
} PESRAAlignmentResult;

/////////////////////////////////////////////////////
/*
    Structure PEInput (Parameters)
    It is the structure to hold the query parameters to Pair-End alignment.
    Some of the parameters are compulsory; some of them are optional.
    PEInput SHOULD ALWAYS BE constructed by PEInputConstruct() function.
*/
/////////////////////////////////////////////////////    
typedef struct PEInput {
    BWT * bwt;                  //Compulsory
    HSP * hsp;                  //Compulsory
    uint8_t OutputType;         //Compulsory. Valid value being PE_REPORT_*.
        
    uint8_t strandLeftLeg;      //Compulsory
    uint8_t strandRightLeg;     //Compulsory
    
    int16_t maxResult;          //Compulsory. Only report the first maxResult
                                //result if set to >=0. Not apply if -1.

    // Compulsory
    ///////////////////////////////////////////////////
    // Either option A or B should be defined and it is compulsory
    // Option B takes precedence over option A.
    //
    // Option A: Insertion size defined by Mean and StdDev.
    int insertSizeMean;         //Optional.
    int insertSizeStdDev;       //Optional.
    // Option B: Insertion size defined by Lower bound and Upper bound.
    int insertLbound;           //Optional.
    int insertUbound;           //Optional.
    ///////////////////////////////////////////////////
} PEInput;


/////////////////////////////////////////////////////
/*
    Structure PEPairs (Output)
    It is the structure to hold 1 Pair-End alignment result.
    The result being held is w.r.t to the read set source.
    The alignment/strand from the first read set is stored in
    occ_1 while those from the second read set is 
    stored in occ_2.
*/
/////////////////////////////////////////////////////    
typedef struct PEPairs {

    SRAOccurrence * occ_1;
    SRAOccurrence * occ_2;
    
    int insertion;
    
} PEPairs;


/////////////////////////////////////////////////////
/*
    Structure PEPairList (Output)
    It is the structure to hold a bucket of alignment result(PEPairs).
*/
/////////////////////////////////////////////////////    
typedef struct PEPairList {
    // The bucket of pair-end results
    PEPairs pairs[PE_MAX_BUCKET_SIZE];
    unsigned long long pairsCount;
    
    // Link to the next bucket
    struct PEPairList * next;
} PEPairList;


/////////////////////////////////////////////////////
/*
    Structure PEOutput (Output)
    It is the structure to enclose a linked-list of bucket(PEPairList).
*/
/////////////////////////////////////////////////////    
typedef struct PEOutput {
    int flag;
    struct PEPairList * root;
    struct PEPairList * tail;
} PEOutput;

/////////////////////////////////////////////////////
/*
    Structure PEStats (Output Counting and Statistics)
    It is the structure to store the statistics for output counting
    and the counts are stored by number of total mismatches on the 
    pair-end results
*/
/////////////////////////////////////////////////////    
typedef struct PEStats {
    unsigned long long NumPEAlignments;
    unsigned long long WithError[PE_MAX_NUM_OF_ERROR];
} PEStats;

/////////////////////////////////////////////////////
/*
    Function PEMapping
    This function takes in 2 lists of SA ranges, which
    are assumed to be enriched by SRAEnrichSARanges.
    
    This function will perform the following,
    1. Initialise output collector
    2. Retrieve the occurrences of all SA ranges
    3. Allocate enough memory for all occurrences
    4. Sort the 2 lists of occurrences
    5. Merge the occurrences to generate a list of PE 
       alignment
       
    Parameters:
        peInput - The input parameter for PE
        peOutput - The output of the PE alignment, which is
                   expected to be constructed by the caller.
                   It is initialised again by PEMapping to allow
                   reusing of constructed PEOutput.
        peStats - The statistics of the PE alignment. 
                  User responsibility to initialise if needed.
        resultListA, 
        resultCountA - The first list of SA ranges being 
                       passed into PEMapping.
        resultListB, 
        resultCountB - The second list of SA ranges being 
                       passed into PEMapping.
                
*/
/////////////////////////////////////////////////////
void PEMapping(PEInput * peInput, PEOutput * peOutput, PEStats * peStats,
                PESRAAlignmentResult * resultListA, unsigned long long resultCountA,
                PESRAAlignmentResult * resultListB, unsigned long long resultCountB,
                int resultListsQualifier);

                
/////////////////////////////////////////////////////
/*
    Function PEMappingOccurrences
    This function basically behaves the same as the above function,
    however it takes in 2 lists of Occurrences, instead of two lists of 
    SA Ranges.
    
    This function will perform the following,
    1. Allocate enough memory for sorting all occurrences
    2. Sort the 2 lists of occurrences
    3. Merge the occurrences to generate a list of PE 
       alignment
       
    Parameters:
        peInput  - The input parameter for PE
        peOutput - The output of the PE alignment, which is
                   expected to be constructed by the caller.
                   It is initialised again by PEMapping to allow
                   reusing of constructed PEOutput.
        peStats - The statistics of the PE alignment. 
                  User responsibility to initialise if needed.
        occList_1, 
        occCount_1 - The first list of occurrences
                   being passed into PEMapping.
        occList_2, 
        occCount_2 - The second list of occurrences
                   being passed into PEMapping.
                
*/
/////////////////////////////////////////////////////
void PEMappingOccurrences(PEInput * peInput, PEOutput * peOutput, PEStats * peStats,
                        SRAOccurrence * occList_1, unsigned long long occCount_1,
                        SRAOccurrence * occList_2, unsigned long long occCount_2,
                        int occListsQualifier);

/////////////////////////////////////////////////////
/*
    Function PEOutputInitialise
    This function cleans up the Pair-End alignment output buffer.
       
    Parameters:
        peOutput - The output of the PE alignment, which is
                   expected to be constructed by the caller.
                   It will be initialised by this function.
                
*/
/////////////////////////////////////////////////////
void PEOutputInitialise(PEOutput * peOutput);

/////////////////////////////////////////////////////
/*
    Function PEReportPairResult
    This function stores a found pair-end alignment result to
    PEOutput buffer.
       
    Parameters:
        peInput - The input of the PE alignment, which is
                   expected to be constructed by the caller.
        peOutput - The output of the PE alignment, which is
                   expected to be constructed by the caller.
        peStats - The statistics of the PE alignment. 
                  It will be updated accordingly.
        occ_1 - The occurrence from the alignment of the read in the
                first read file.
        occ_2 - The occurrence from the alignment of the read in the
                second read file.
        insertion - The gap between the pair
                
*/
/////////////////////////////////////////////////////
void PEReportPairResult(PEInput * peInput, PEOutput * peOutput, PEStats * peStats, SRAOccurrence * occ_1, SRAOccurrence * occ_2, unsigned int insertion);



/////////////////////////////////////////////////////
/*
    Function PEIsPairEndMatch
    This function takes in 2 occurrence and determind
    if the two occurrences form a pair end mapping.
    Assume left <= right.
    
    Input:
        peInput - The input parameter for PE
        ambPos_left - The ambiguous position on left leg
        strand_left - The strand on left leg
        matchLen_left - The match length on left leg
        ambPos_right - The ambiguous position on left leg
        strand_right - The strand on left leg
        matchLen_right - The match length on left leg
        insertion - The insertion size output
        
    Output:
        1 if the two occurrences form a pair end mapping
        0 otherwise
        
*/
/////////////////////////////////////////////////////
int PEIsPairEndMatch(PEInput * peInput, 
                    unsigned long long ambPos_left, int strand_left, unsigned int matchLen_left, 
                    unsigned long long ambPos_right, int strand_right, unsigned int matchLen_right, 
                    unsigned int * insertion);
int PESRAIsPairEndMatch(PEInput * peInput, SRAOccurrence * occ_1, SRAOccurrence * occ_2, unsigned int * insertion);
int PEDPIsPairEndMatch(PEInput * peInput, DPOccurrence * occ_1, DPOccurrence * occ_2, unsigned int * insertion);

/////////////////////////////////////////////////////
// Constructor and Destructor
/////////////////////////////////////////////////////
PEInput * PEInputConstruct(BWT * bwt, HSP * hsp);
PEInput * PEInputMakeClone(PEInput * source);
void PEInputFree(PEInput * peInput);
void PEInputCloneFree(PEInput * peInput);
PEOutput * PEOutputConstruct();
void PEOutputFree(PEOutput * peOutput);
PEStats * PESTATConstruct();
void PESTATInitialise(PEStats * peStats);
void PESTATFree(PEStats * peStats);

/////////////////////////////////////////////////////
// Utilities
/////////////////////////////////////////////////////
unsigned long long PECountPEOutput ( PEOutput * peOutput );
unsigned long long PEStatsPEOutput ( PEOutput * peOutput, PEStats * peStats, PEPairs ** optimal, unsigned int * mismatchStats );

unsigned long long PECountPEPairList(PEPairList * pairList);
unsigned long long PECountPEOutput(PEOutput * peOutput) ;
int PEAligned(PEOutput * peOutput);
#endif

