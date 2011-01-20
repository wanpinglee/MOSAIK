// ***************************************************************************
// CPairwiseUtilities - a centralized location for all parameters related to
//                      pairwise alignment routines.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "PairwiseUtilities.h"

// sets the Smith-Waterman match score
float CPairwiseUtilities::MatchScore              = 10.0f;

// sets the Smith-Waterman mismatch score
float CPairwiseUtilities::MismatchScore           = -9.0f;

// sets the Smith-Waterman gap open penalty
float CPairwiseUtilities::GapOpenPenalty          = 15.0f;

// sets the Smith-Waterman gap extend penalty
float CPairwiseUtilities::GapExtendPenalty        = 6.66f;

// toggles if we should filter sequences that match too little of the original sequence
unsigned int CPairwiseUtilities::MinAlignment  = 0;
bool CPairwiseUtilities::UseMinAlignmentFilter = false;

// toggles if we should filter sequences that match too little of the original sequence
double CPairwiseUtilities::MinPercentAlignment        = 0.0;
bool CPairwiseUtilities::UseMinAlignmentPercentFilter = false;

// toggles if we should filter sequences that have low alignment qualities
unsigned char CPairwiseUtilities::MinAlignmentQuality = 0;
bool CPairwiseUtilities::UseMinAlignmentQualityFilter = false;

// toggles if we should filter sequences with too many mismatches (count)
unsigned int CPairwiseUtilities::MaxNumMismatches = -1;
bool CPairwiseUtilities::UseMismatchFilter        = false;

// toggles if we should filter sequences with too many mismatches (percent)
double CPairwiseUtilities::MaxMismatchPercent     = 1.0;
bool CPairwiseUtilities::UseMismatchPercentFilter = false;

// toggles if we should compare mismatches against the aligned read length or the original read length
bool CheckMismatchesAgainstReadLength             = false;

// sets the Smith-Waterman homo-polymer gap open penalty
bool CPairwiseUtilities::UseHomoPolymerGapOpenPenalty = false;
float CPairwiseUtilities::HomoPolymerGapOpenPenalty   = 4.0f;
