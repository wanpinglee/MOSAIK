// ***************************************************************************
// CSmithWatermanGotoh - aligns reads using a Smith-Waterman algorithm.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iostream>
#include <algorithm>
#include <memory>
#include "Alignment.h"
#include "Mosaik.h"

using namespace std;

#define MOSAIK_NUM_NUCLEOTIDES 26
#define GAP '-'

class CSmithWatermanGotoh {
public:
	// constructor
	CSmithWatermanGotoh(float matchScore, float mismatchScore, float gapOpenPenalty, float gapExtendPenalty);
	// destructor
	~CSmithWatermanGotoh(void);
	// aligns the query sequence to the reference using the Smith Waterman Gotoh algorithm
	void Align(Alignment& alignment, const char* s1, const unsigned int s1Length, const char* s2, const unsigned int s2Length);
	// enables homo-polymer scoring
	void EnableHomoPolymerGapPenalty(float hpGapOpenPenalty);

#ifndef WINUNIT
private:
#endif
	// creates a simple scoring matrix to align the nucleotides and the ambiguity code N
	void CreateScoringMatrix(void);
	// corrects the homopolymer gap order for forward alignments
	static void CorrectHomopolymerGapOrder(Alignment& al);
	// returns the maximum floating point number
	static inline float MaxFloats(const float& a, const float& b, const float& c);
	// our simple scoring matrix
	float mScoringMatrix[MOSAIK_NUM_NUCLEOTIDES][MOSAIK_NUM_NUCLEOTIDES];
	// keep track of maximum initialized sizes
	unsigned int mCurrentMatrixSize;
	unsigned int mCurrentAnchorSize;
	unsigned int mCurrentQuerySize;
	unsigned int mCurrentAQSumSize;
	// define our traceback directions
	// N.B. This used to be defined as an enum, but gcc doesn't like being told
	// which storage class to use
	const static char Directions_STOP;
	const static char Directions_LEFT;
	const static char Directions_DIAGONAL;
	const static char Directions_UP;
	// define scoring constants
	const float mMatchScore;
	const float mMismatchScore;
	const float mGapOpenPenalty;
	const float mGapExtendPenalty;
	// store the backtrace pointers
	char* mPointers;
	// store the vertical gap sizes - assuming gaps are not longer than 32768 bases long
	short* mSizesOfVerticalGaps;
	// store the horizontal gap sizes - assuming gaps are not longer than 32768 bases long
	short* mSizesOfHorizontalGaps;	
	// score if xi aligns to a gap after yi
	float* mQueryGapScores;
	// best score of alignment x1...xi to y1...yi
	float* mBestScores;
	// our reversed alignment
	char* mReversedAnchor;
	char* mReversedQuery;
	// define static constants
	static const float FLOAT_NEGATIVE_INFINITY;
	// toggles the use of the homo-polymer gap open penalty
	bool mUseHomoPolymerGapOpenPenalty;
	// specifies the homo-polymer gap open penalty
	float mHomoPolymerGapOpenPenalty;
};

// returns the maximum floating point number
inline float CSmithWatermanGotoh::MaxFloats(const float& a, const float& b, const float& c) {
	float max = 0.0f;
	if(a > max) max = a;
	if(b > max) max = b;
	if(c > max) max = c;
	return max;
}
