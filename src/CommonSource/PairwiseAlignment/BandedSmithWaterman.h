// ***************************************************************************
// CBandedSmithWaterman - aligns reads using a banded Smith-Waterman algorithm.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Wan-Ping Lee & Michael Strömberg
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
#include "HashRegion.h"

using namespace std;

#define MOSAIK_NUM_NUCLEOTIDES 26
#define GAP '-'

typedef unsigned char DirectionType;
typedef unsigned char PositionType;

struct ElementInfo {
	unsigned int Direction             : 2;
	unsigned int mSizeOfVerticalGaps   : 15;
	unsigned int mSizeOfHorizontalGaps : 15;
};

class CBandedSmithWaterman {
public:
	// constructor
	CBandedSmithWaterman(float matchScore, float mismatchScore, float gapOpenPenalty, float gapExtendPenalty, unsigned int bandWidth);
	// destructor
	~CBandedSmithWaterman(void);
	// aligns the query sequence to the anchor using the Smith Waterman Gotoh algorithm
	void Align(Alignment& alignment, const char* s1, const unsigned int s1Length, const char* s2, const unsigned int s2Length, HashRegion& hr);
	// enables homo-polymer scoring
	void EnableHomoPolymerGapPenalty(float hpGapOpenPenalty);
private:
	// calculates the score during the forward algorithm
	float CalculateScore(const char* s1, const char* s2, const unsigned int rowNum, const unsigned int columnNum, float& currentQueryGapScore, const unsigned int rowOffset, const unsigned int columnOffset);
	// creates a simple scoring matrix to align the nucleotides and the ambiguity code N
	void CreateScoringMatrix(void);
	// corrects the homopolymer gap order for forward alignments
	static void CorrectHomopolymerGapOrder(Alignment& al);
	// returns the maximum floating point number
	static inline float MaxFloats(const float& a, const float& b, const float& c);
	// reinitializes the matrices
	void ReinitializeMatrices(const PositionType& positionType, const unsigned int& s1Length, const unsigned int& s2Length, const HashRegion& hr);
	// performs the backtrace algorithm
	void Traceback(Alignment& alignment, const char* s1, const char* s2, const unsigned int s2Length, unsigned int bestRow, unsigned int bestColumn, const unsigned int rowOffset, const unsigned int columnOffset);
	// updates the best score during the forward algorithm
	inline void UpdateBestScore(unsigned int& bestRow, unsigned int& bestColumn, float& bestScore, const unsigned int rowNum, const unsigned int columnNum, const float score);
	// our simple scoring matrix
	float mScoringMatrix[MOSAIK_NUM_NUCLEOTIDES][MOSAIK_NUM_NUCLEOTIDES];
	// keep track of maximum initialized sizes
	unsigned int mCurrentMatrixSize;
	unsigned int mCurrentAnchorSize;
	unsigned int mCurrentAQSumSize;
	unsigned int mBandwidth;
	// define our backtrace directions
	const static DirectionType Directions_STOP;
	const static DirectionType Directions_LEFT;
	const static DirectionType Directions_DIAGONAL;
	const static DirectionType Directions_UP;
	// store the backtrace pointers
	ElementInfo* mPointers;
	// define our position types
	const static PositionType Position_REF_AND_QUERY_ZERO;
	const static PositionType Position_REF_ZERO;
	const static PositionType Position_QUERY_ZERO;
	const static PositionType Position_REF_AND_QUERO_NONZERO;
	// define scoring constants
	const float mMatchScore;
	const float mMismatchScore;
	const float mGapOpenPenalty;
	const float mGapExtendPenalty;
	// score if xi aligns to a gap after yi
	float* mAnchorGapScores;
	// best score of alignment x1...xi to y1...yi
	float* mBestScores;
	// our reversed alignment
	char* mReversedAnchor;
	char* mReversedQuery;
	// define static constants
	static const float FLOAT_NEGATIVE_INFINITY;
	// toggles the use of the homo-polymer gap open penalty
	bool mUseHomoPolymerGapOpenPenalty;
	float mHomoPolymerGapOpenPenalty;
};

// returns the maximum floating point number
inline float CBandedSmithWaterman::MaxFloats(const float& a, const float& b, const float& c) {
	float max = 0.0f;
	if(a > max) max = a;
	if(b > max) max = b;
	if(c > max) max = c;
	return max;
}

// updates the best score during the forward algorithm
inline void CBandedSmithWaterman::UpdateBestScore(unsigned int& bestRow, unsigned int& bestColumn, float& bestScore, const unsigned int rowNum, const unsigned int columnNum, const float score) {

	//const unsigned int row    = rowNum + rowOffset;
	//const unsigned int column = columnOffset - rowNum + columnNum;

	if(score >= bestScore) {
		bestRow    = rowNum;
		bestColumn = columnNum;
		bestScore  = score;
	}
}
