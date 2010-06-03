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

#pragma once

class CPairwiseUtilities {
public:
	// toggles if we should filter sequences with too many mismatches (count)
	static bool UseMismatchFilter;
	// toggles if we should filter sequences with too many mismatches (percent)
	static bool UseMismatchPercentFilter;
	// sets the maximum number of mismatches to allow if the mismatch filter is enabled
	static unsigned int MaxNumMismatches;
	// sets the maximum mismatch percent to allow if the mismatch percent filter is enabled
	static double MaxMismatchPercent;
	// toggles if we should filter sequences that match too little of the original sequence
	static bool UseMinAlignmentFilter;
	// sets the minimum number of aligned nucleotides to allow if the minimum alignment filter is enabled
	static unsigned int MinAlignment;
	// toggles if we should filter sequences that match too little of the original sequence
	static bool UseMinAlignmentPercentFilter;
	// sets the minimum percentage of alignment to allow if the minimum alignment filter is enabled
	static double MinPercentAlignment;
	// toggles if we should filter sequences that have low alignment qualities
	static bool UseMinAlignmentQualityFilter;
	// sets the minimum aligment quality to allow if the minimum alignment quality filter is enabled
	static unsigned char MinAlignmentQuality;
	// toggles if we should compare mismatches against the aligned read length or the original read length
	static bool CheckMismatchesAgainstReadLength;
	// sets the Smith-Waterman match score
	static float MatchScore;
	// sets the Smith-Waterman mismatch score
	static float MismatchScore;
	// sets the Smith-Waterman gap open penalty
	static float GapOpenPenalty;
	// sets the Smith-Waterman gap extend penalty
	static float GapExtendPenalty;
	// toggles the use of the Smith-Waterman homo-polymer gap open penalty
	static bool UseHomoPolymerGapOpenPenalty;
	// sets the Smith-Waterman homo-polymer gap open penalty
	static float HomoPolymerGapOpenPenalty;
};
