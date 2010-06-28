// ***************************************************************************
// CSmithWatermanGotoh - aligns reads using a Smith-Waterman algorithm.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "SmithWatermanGotoh.h"

const float CSmithWatermanGotoh::FLOAT_NEGATIVE_INFINITY = (float)-1e+30;

const char CSmithWatermanGotoh::Directions_STOP     = 0;
const char CSmithWatermanGotoh::Directions_LEFT     = 1;
const char CSmithWatermanGotoh::Directions_DIAGONAL = 2;
const char CSmithWatermanGotoh::Directions_UP       = 3;

CSmithWatermanGotoh::CSmithWatermanGotoh(float matchScore, float mismatchScore, float gapOpenPenalty, float gapExtendPenalty) 
: mCurrentMatrixSize(0)
, mCurrentAnchorSize(0)
, mCurrentQuerySize(0)
, mCurrentAQSumSize(0)
, mMatchScore(matchScore)
, mMismatchScore(mismatchScore)
, mGapOpenPenalty(gapOpenPenalty)
, mGapExtendPenalty(gapExtendPenalty)
, mPointers(NULL)
, mSizesOfVerticalGaps(NULL)
, mSizesOfHorizontalGaps(NULL)
, mQueryGapScores(NULL)
, mBestScores(NULL)
, mReversedAnchor(NULL)
, mReversedQuery(NULL)
, mUseHomoPolymerGapOpenPenalty(false)
{
	CreateScoringMatrix();
}

CSmithWatermanGotoh::~CSmithWatermanGotoh(void) {
	if(mPointers)              delete [] mPointers;
	if(mSizesOfVerticalGaps)   delete [] mSizesOfVerticalGaps;
	if(mSizesOfHorizontalGaps) delete [] mSizesOfHorizontalGaps;
	if(mQueryGapScores)        delete [] mQueryGapScores;
	if(mBestScores)            delete [] mBestScores;
	if(mReversedAnchor)        delete [] mReversedAnchor;
	if(mReversedQuery)         delete [] mReversedQuery;
}

// aligns the query sequence to the reference using the Smith Waterman Gotoh algorithm
void CSmithWatermanGotoh::Align(Alignment& alignment, const char* s1, const unsigned int s1Length, const char* s2, const unsigned int s2Length) {

	if((s1Length == 0) || (s2Length == 0)) {
		cout << "ERROR: Found a read with a zero length." << endl;
		exit(1);
	}

	unsigned int referenceLen      = s1Length + 1;
	unsigned int queryLen          = s2Length + 1;
	unsigned int sequenceSumLength = s1Length + s2Length;

	// reinitialize our matrices

	if((referenceLen * queryLen) > mCurrentMatrixSize) {

		// calculate the new matrix size
		mCurrentMatrixSize = referenceLen * queryLen;

		// delete the old arrays
		if(mPointers)              delete [] mPointers;
		if(mSizesOfVerticalGaps)   delete [] mSizesOfVerticalGaps;
		if(mSizesOfHorizontalGaps) delete [] mSizesOfHorizontalGaps;

		try {

			// initialize the arrays
			mPointers              = new char[mCurrentMatrixSize];
			mSizesOfVerticalGaps   = new short[mCurrentMatrixSize];
			mSizesOfHorizontalGaps = new short[mCurrentMatrixSize];

		} catch(bad_alloc) {
			cout << "ERROR: Unable to allocate enough memory for the Smith-Waterman algorithm." << endl;
			exit(1);
		}
	}

	// initialize the traceback matrix to STOP
	memset((char*)mPointers, 0, SIZEOF_CHAR * queryLen);
	for(unsigned int i = 1; i < referenceLen; i++) mPointers[i * queryLen] = 0;

	// initialize the gap matrices to 1
	uninitialized_fill(mSizesOfVerticalGaps, mSizesOfVerticalGaps + mCurrentMatrixSize, 1);
	uninitialized_fill(mSizesOfHorizontalGaps, mSizesOfHorizontalGaps + mCurrentMatrixSize, 1);

	//
	// construct
	//

	// reinitialize our query-dependent arrays
	if(s2Length > mCurrentQuerySize) {

		// calculate the new query array size
		mCurrentQuerySize = s2Length;

		// delete the old arrays
		if(mQueryGapScores) delete [] mQueryGapScores;
		if(mBestScores)     delete [] mBestScores;

		// initialize the arrays
		try {

			mQueryGapScores = new float[mCurrentQuerySize + 1];
			mBestScores     = new float[mCurrentQuerySize + 1];

		} catch(bad_alloc) {
			cout << "ERROR: Unable to allocate enough memory for the Smith-Waterman algorithm." << endl;
			exit(1);
		}
	}

	// reinitialize our reference+query-dependent arrays
	if(sequenceSumLength > mCurrentAQSumSize) {

		// calculate the new reference array size
		mCurrentAQSumSize = sequenceSumLength;

		// delete the old arrays
		if(mReversedAnchor) delete [] mReversedAnchor;
		if(mReversedQuery)  delete [] mReversedQuery;

		// initialize the arrays
		try {

			mReversedAnchor = new char[mCurrentAQSumSize + 1];	// reversed sequence #1
			mReversedQuery  = new char[mCurrentAQSumSize + 1];	// reversed sequence #2

		} catch(bad_alloc) {
			cout << "ERROR: Unable to allocate enough memory for the Smith-Waterman algorithm." << endl;
			exit(1);
		}
	}

	// initialize the gap score and score vectors
	uninitialized_fill(mQueryGapScores, mQueryGapScores + queryLen, FLOAT_NEGATIVE_INFINITY);
	memset((char*)mBestScores, 0, SIZEOF_FLOAT * queryLen);

	float similarityScore, totalSimilarityScore, bestScoreDiagonal;
	float queryGapExtendScore, queryGapOpenScore;
	float referenceGapExtendScore, referenceGapOpenScore, currentAnchorGapScore;

	unsigned int BestColumn = 0;
	unsigned int BestRow    = 0;
	float BestScore         = FLOAT_NEGATIVE_INFINITY;

	for(unsigned int i = 1, k = queryLen; i < referenceLen; i++, k += queryLen) {

		currentAnchorGapScore = FLOAT_NEGATIVE_INFINITY;
		bestScoreDiagonal = mBestScores[0];

		for(unsigned int j = 1, l = k + 1; j < queryLen; j++, l++) {

			// calculate our similarity score
			similarityScore = mScoringMatrix[s1[i - 1] - 'A'][s2[j - 1] - 'A'];

			// fill the matrices
			totalSimilarityScore = bestScoreDiagonal + similarityScore;

			//cout << "i: " << i << ", j: " << j << ", totalSimilarityScore: " << totalSimilarityScore << endl;

			queryGapExtendScore = mQueryGapScores[j] - mGapExtendPenalty;
			queryGapOpenScore   = mBestScores[j] - mGapOpenPenalty;

			// compute the homo-polymer gap score if enabled
			if(mUseHomoPolymerGapOpenPenalty)
				if((j > 1) && (s2[j - 1] == s2[j - 2]))
					queryGapOpenScore = mBestScores[j] - mHomoPolymerGapOpenPenalty;

			if(queryGapExtendScore > queryGapOpenScore) {
				mQueryGapScores[j] = queryGapExtendScore;
				mSizesOfVerticalGaps[l] = (short)(mSizesOfVerticalGaps[l - queryLen] + 1);
			} else mQueryGapScores[j] = queryGapOpenScore;

			referenceGapExtendScore = currentAnchorGapScore - mGapExtendPenalty;
			referenceGapOpenScore   = mBestScores[j - 1] - mGapOpenPenalty;

			// compute the homo-polymer gap score if enabled
			if(mUseHomoPolymerGapOpenPenalty)
				if((i > 1) && (s1[i - 1] == s1[i - 2]))
					referenceGapOpenScore = mBestScores[j - 1] - mHomoPolymerGapOpenPenalty;

			if(referenceGapExtendScore > referenceGapOpenScore) {
				currentAnchorGapScore = referenceGapExtendScore;
				mSizesOfHorizontalGaps[l] = (short)(mSizesOfHorizontalGaps[l - 1] + 1);
			} else currentAnchorGapScore = referenceGapOpenScore;

			bestScoreDiagonal = mBestScores[j];
			mBestScores[j] = MaxFloats(totalSimilarityScore, mQueryGapScores[j], currentAnchorGapScore);

			// determine the traceback direction
			// diagonal (445364713) > stop (238960195) > up (214378647) > left (166504495)
			if(mBestScores[j] == 0)                         mPointers[l] = Directions_STOP;
			else if(mBestScores[j] == totalSimilarityScore) mPointers[l] = Directions_DIAGONAL;
			else if(mBestScores[j] == mQueryGapScores[j])   mPointers[l] = Directions_UP;
			else                                            mPointers[l] = Directions_LEFT;

			// set the traceback start at the current cell i, j and score
			if(mBestScores[j] > BestScore) {
				BestRow    = i;
				BestColumn = j;
				BestScore  = mBestScores[j];
			}
		}
	}

	//
	// traceback
	//

	// aligned sequences
	int gappedAnchorLen  = 0;   // length of sequence #1 after alignment
	int gappedQueryLen   = 0;   // length of sequence #2 after alignment
	int numMismatches    = 0;   // the mismatched nucleotide count

	char c1, c2;

	int ci = BestRow;
	int cj = BestColumn;
	int ck = ci * queryLen;

	// traceback flag
	bool keepProcessing = true;
	bool hasGap = false;

	while(keepProcessing) {

		// diagonal (445364713) > stop (238960195) > up (214378647) > left (166504495)
		switch(mPointers[ck + cj]) {

			case Directions_DIAGONAL:
				c1 = s1[--ci];
				c2 = s2[--cj];
				ck -= queryLen;

				mReversedAnchor[gappedAnchorLen++] = c1;
				mReversedQuery[gappedQueryLen++]   = c2;

				// increment our mismatch counter
				if(mScoringMatrix[c1 - 'A'][c2 - 'A'] == mMismatchScore) numMismatches++;	
				break;

			case Directions_STOP:
				keepProcessing = false;
				break;

			case Directions_UP:
				for(unsigned int l = 0, len = mSizesOfVerticalGaps[ck + cj]; l < len; l++) {
					mReversedAnchor[gappedAnchorLen++] = s1[--ci];
					mReversedQuery[gappedQueryLen++]   = GAP;
					ck -= queryLen;
					numMismatches++;
				}
				hasGap = true;
				break;

			case Directions_LEFT:
				for(unsigned int l = 0, len = mSizesOfHorizontalGaps[ck + cj]; l < len; l++) {
					mReversedAnchor[gappedAnchorLen++] = GAP;
					mReversedQuery[gappedQueryLen++]   = s2[--cj];
					numMismatches++;
				}
				hasGap = true;
				break;
		}
	}

	// define the reference and query sequences
	mReversedAnchor[gappedAnchorLen] = 0;
	mReversedQuery[gappedQueryLen]   = 0;

	// catch sequences with different lengths
	if(gappedAnchorLen != gappedQueryLen) {
		cout << "ERROR: The aligned sequences have different lengths after Smith-Waterman-Gotoh algorithm." << endl;
		exit(1);
	}

	// reverse the strings and assign them to our alignment structure
	reverse(mReversedAnchor, mReversedAnchor + gappedAnchorLen);
	reverse(mReversedQuery,  mReversedQuery  + gappedQueryLen);

	alignment.Reference = mReversedAnchor;
	alignment.Query     = mReversedQuery;

	// set the reference endpoints
	alignment.ReferenceBegin = ci;
	alignment.ReferenceEnd   = BestRow - 1;

	// set the query endpoints
	if(alignment.IsReverseStrand) {
		alignment.QueryBegin = s2Length - BestColumn;
		alignment.QueryEnd   = s2Length - cj - 1;
	} else {
		alignment.QueryBegin = cj;
		alignment.QueryEnd   = BestColumn - 1;
	}

	// set the query length and number of mismatches
	alignment.QueryLength = alignment.QueryEnd - alignment.QueryBegin + 1;
	alignment.NumMismatches  = numMismatches;

	// fix the gap order
	if(hasGap) CorrectHomopolymerGapOrder(alignment);
}

// creates a simple scoring matrix to align the nucleotides and the ambiguity code N
void CSmithWatermanGotoh::CreateScoringMatrix(void) {

	unsigned int nIndex = 13;
	unsigned int xIndex = 23;

	// define the N score to be 1/4 of the span between mismatch and match
	//const short nScore = mMismatchScore + (short)(((mMatchScore - mMismatchScore) / 4.0) + 0.5);

	// calculate the scoring matrix
	for(unsigned char i = 0; i < MOSAIK_NUM_NUCLEOTIDES; i++) {
		for(unsigned char j = 0; j < MOSAIK_NUM_NUCLEOTIDES; j++) {

			// N.B. matching N to everything (while conceptually correct) leads to some
			// bad alignments, lets make N be a mismatch instead.

			// add the matches or mismatches to the hashtable (N is a mismatch)
			if((i == nIndex) || (j == nIndex)) mScoringMatrix[i][j] = mMismatchScore;
			else if((i == xIndex) || (j == xIndex)) mScoringMatrix[i][j] = mMismatchScore;
			else if(i == j) mScoringMatrix[i][j] = mMatchScore;
			else mScoringMatrix[i][j] = mMismatchScore;
		}
	}

	// add ambiguity codes
	mScoringMatrix['M' - 'A']['A' - 'A'] = mMatchScore;	// M - A
	mScoringMatrix['A' - 'A']['M' - 'A'] = mMatchScore;
	mScoringMatrix['M' - 'A']['C' - 'A'] = mMatchScore; // M - C
	mScoringMatrix['C' - 'A']['M' - 'A'] = mMatchScore;

	mScoringMatrix['R' - 'A']['A' - 'A'] = mMatchScore;	// R - A
	mScoringMatrix['A' - 'A']['R' - 'A'] = mMatchScore;
	mScoringMatrix['R' - 'A']['G' - 'A'] = mMatchScore; // R - G
	mScoringMatrix['G' - 'A']['R' - 'A'] = mMatchScore;

	mScoringMatrix['W' - 'A']['A' - 'A'] = mMatchScore;	// W - A
	mScoringMatrix['A' - 'A']['W' - 'A'] = mMatchScore;
	mScoringMatrix['W' - 'A']['T' - 'A'] = mMatchScore; // W - T
	mScoringMatrix['T' - 'A']['W' - 'A'] = mMatchScore;

	mScoringMatrix['S' - 'A']['C' - 'A'] = mMatchScore;	// S - C
	mScoringMatrix['C' - 'A']['S' - 'A'] = mMatchScore;
	mScoringMatrix['S' - 'A']['G' - 'A'] = mMatchScore; // S - G
	mScoringMatrix['G' - 'A']['S' - 'A'] = mMatchScore;

	mScoringMatrix['Y' - 'A']['C' - 'A'] = mMatchScore;	// Y - C
	mScoringMatrix['C' - 'A']['Y' - 'A'] = mMatchScore;
	mScoringMatrix['Y' - 'A']['T' - 'A'] = mMatchScore; // Y - T
	mScoringMatrix['T' - 'A']['Y' - 'A'] = mMatchScore;

	mScoringMatrix['K' - 'A']['G' - 'A'] = mMatchScore;	// K - G
	mScoringMatrix['G' - 'A']['K' - 'A'] = mMatchScore;
	mScoringMatrix['K' - 'A']['T' - 'A'] = mMatchScore; // K - T
	mScoringMatrix['T' - 'A']['K' - 'A'] = mMatchScore;

	mScoringMatrix['V' - 'A']['A' - 'A'] = mMatchScore;	// V - A
	mScoringMatrix['A' - 'A']['V' - 'A'] = mMatchScore;
	mScoringMatrix['V' - 'A']['C' - 'A'] = mMatchScore; // V - C
	mScoringMatrix['C' - 'A']['V' - 'A'] = mMatchScore;
	mScoringMatrix['V' - 'A']['G' - 'A'] = mMatchScore; // V - G
	mScoringMatrix['G' - 'A']['V' - 'A'] = mMatchScore;

	mScoringMatrix['H' - 'A']['A' - 'A'] = mMatchScore;	// H - A
	mScoringMatrix['A' - 'A']['H' - 'A'] = mMatchScore;
	mScoringMatrix['H' - 'A']['C' - 'A'] = mMatchScore; // H - C
	mScoringMatrix['C' - 'A']['H' - 'A'] = mMatchScore;
	mScoringMatrix['H' - 'A']['T' - 'A'] = mMatchScore; // H - T
	mScoringMatrix['T' - 'A']['H' - 'A'] = mMatchScore;

	mScoringMatrix['D' - 'A']['A' - 'A'] = mMatchScore;	// D - A
	mScoringMatrix['A' - 'A']['D' - 'A'] = mMatchScore;
	mScoringMatrix['D' - 'A']['G' - 'A'] = mMatchScore; // D - G
	mScoringMatrix['G' - 'A']['D' - 'A'] = mMatchScore;
	mScoringMatrix['D' - 'A']['T' - 'A'] = mMatchScore; // D - T
	mScoringMatrix['T' - 'A']['D' - 'A'] = mMatchScore;

	mScoringMatrix['B' - 'A']['C' - 'A'] = mMatchScore;	// B - C
	mScoringMatrix['C' - 'A']['B' - 'A'] = mMatchScore;
	mScoringMatrix['B' - 'A']['G' - 'A'] = mMatchScore; // B - G
	mScoringMatrix['G' - 'A']['B' - 'A'] = mMatchScore;
	mScoringMatrix['B' - 'A']['T' - 'A'] = mMatchScore; // B - T
	mScoringMatrix['T' - 'A']['B' - 'A'] = mMatchScore;
}

// enables homo-polymer scoring
void CSmithWatermanGotoh::EnableHomoPolymerGapPenalty(float hpGapOpenPenalty) {
	mUseHomoPolymerGapOpenPenalty = true;
	mHomoPolymerGapOpenPenalty    = hpGapOpenPenalty;
}

// corrects the homopolymer gap order for forward alignments
void CSmithWatermanGotoh::CorrectHomopolymerGapOrder(Alignment& al) {

	// this is only required for alignments with mismatches
	if(al.NumMismatches == 0) return;

	// initialize
	char* pReference = al.Reference.Data();
	char* pQuery     = al.Query.Data();

	char* pGapString  = NULL;
	char* pBaseString = NULL;

	const unsigned int pairwiseLen = al.Reference.Length();
	const unsigned int lastIndex   = pairwiseLen - 1;

	// identify the gapped regions
	for(unsigned int i = 0; i < pairwiseLen; ++i) {

		// find the gaps in the reference
		bool foundGap = false;

		if(pReference[i] == '-') {
			foundGap    = true;
			pGapString  = pReference;
			pBaseString = pQuery;
		}

		// find gaps in the query
		if(pQuery[i] == '-') {
			foundGap    = true;
			pBaseString = pReference;
			pGapString  = pQuery;
		}

		// investigate the gap region
		if(foundGap) {

			unsigned int index  = i;
			while((pGapString[index] == '-') && (index < lastIndex)) index++;
			const unsigned int length = index - i;

			// check if this stretch is a homopolymer
			bool isHomopolymer = true;
			const char firstBase = pBaseString[i];
			for(unsigned int j = 1; j < length; ++j) {
				if(pBaseString[i + j] != firstBase) {
					isHomopolymer = false;
					break;
				}
			}

			// if this is a homopolymer, extend the coordinates
			if(isHomopolymer) {

				// find the left endpoint
				unsigned int numExtraBases = 0;
				unsigned int lowIndex = i;
				while((pGapString[lowIndex] == firstBase) && (pBaseString[lowIndex] == firstBase) && (lowIndex != 0)) {
					--lowIndex;
					++numExtraBases;
				}

				// find the right endpoint
				unsigned int highIndex = index;
				while((pGapString[highIndex] == firstBase) && (pBaseString[highIndex] == firstBase) && (highIndex < lastIndex)) {
					++highIndex;
					++numExtraBases;
				}

				// write the reordered gap
				if(numExtraBases > 0) {
					memset(pGapString + lowIndex, firstBase, numExtraBases);
					memset(pGapString + lowIndex + numExtraBases, GAP, highIndex - lowIndex - numExtraBases);
				}
			} 

			i = index;
		}
	}
}
