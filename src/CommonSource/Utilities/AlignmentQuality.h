// ***************************************************************************
// CAlignmentQuality - calculates the alignment qualities from a technology-
//                     specific logistic regression model.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <cmath>
#include <vector>
#include "Alignment.h"

using namespace std;

// preprocess pattern inputs
#define MIN_BITS_ILLUMINA            0.0f
#define MIN_MM_BQ_PERCENT_ILLUMINA   0.0f
#define MIN_LN_READ_LEN_ILLUMINA     3.583518938456110f
#define MIN_LN_REF_LEN_ILLUMINA      19.231609711279024f

#define MIN_BITS_454                 0.0f
#define MIN_MM_BQ_PERCENT_454        0.0f
#define MIN_LN_READ_LEN_454          3.583518938456110f
#define MIN_LN_REF_LEN_454           19.231609711279024f

#define RANGE_BITS_ILLUMINA          20.0f
#define RANGE_MM_BQ_PERCENT_ILLUMINA 20.0f
#define RANGE_LN_READ_LEN_ILLUMINA   0.747214401830221f
#define RANGE_LN_REF_LEN_ILLUMINA    2.532772211451253f

#define RANGE_BITS_454               20.0f
#define RANGE_MM_BQ_PERCENT_454      20.0f
#define RANGE_LN_READ_LEN_454        2.631089159966082f
#define RANGE_LN_REF_LEN_454         2.532772211451253f

#define RANGE_AQ_ILLUMINA            64.599999999999994f
#define RANGE_AQ_454                 66.930000000000007f

#define LAYER_BIAS_ILLUMINA          -0.564142895586588f
#define LAYER_BIAS_454               -1.058189629447166f

#define NUM_MMBQ_BINS                21

class CAlignmentQuality {
public:
	// constructor
	CAlignmentQuality(bool usingIllumina, unsigned int refLen);
	// destructor
	~CAlignmentQuality(void);
	// calculates the alignment quality
	void CalculateQuality(vector<Alignment>::iterator& alIter) const;

private:
	// denormalizes the output pattern
	inline unsigned char Denormalize(float alignmentQuality) const;
	// normalizes the input pattern
	inline void Normalize(float& bits, float& mmPercent, float& readLen, float& refLen) const;
	// our reference length
	float mLnReferenceLength;
	// denotes whether or not we're using the 454 or Illumina model
	bool mUsingIllumina; 
	// our bias vectors
	static const float mB_Illumina[30];
	static const float mB_454[40];
	// our input weight matrices
	static const float mIW_Illumina[30][4];
	static const float mIW_454[40][4];
	// our layer weight matrices
	static const float mLW_Illumina[30];
	static const float mLW_454[40];
	// our mismatch bin boundaries
	static const double mMismatchBQBins[21];
};

// denormalizes the output pattern
inline unsigned char CAlignmentQuality::Denormalize(float alignmentQuality) const {
	float aq = (mUsingIllumina ? RANGE_AQ_ILLUMINA : RANGE_AQ_454) * (alignmentQuality + 1.0f) * 0.5f;
	if(aq < 0.0f)  aq = 0.0f;
	if(aq > 99.0f) aq = 99.0f;
	return (unsigned char)(aq + 0.5);
}

// normalizes the input pattern
inline void CAlignmentQuality::Normalize(float& bits, float& mmPercent, float& readLen, float& refLen) const {
	if(mUsingIllumina) {
		bits      = 2.0f * (bits      - MIN_BITS_ILLUMINA)          / RANGE_BITS_ILLUMINA          - 1.0f;
		mmPercent = 2.0f * (mmPercent - MIN_MM_BQ_PERCENT_ILLUMINA) / RANGE_MM_BQ_PERCENT_ILLUMINA - 1.0f;
		readLen   = 2.0f * (readLen   - MIN_LN_READ_LEN_ILLUMINA)   / RANGE_LN_READ_LEN_ILLUMINA   - 1.0f;
		refLen    = 2.0f * (refLen    - MIN_LN_REF_LEN_ILLUMINA)    / RANGE_LN_REF_LEN_ILLUMINA    - 1.0f;
	} else {
		bits      = 2.0f * (bits      - MIN_BITS_454)          / RANGE_BITS_454          - 1.0f;
		mmPercent = 2.0f * (mmPercent - MIN_MM_BQ_PERCENT_454) / RANGE_MM_BQ_PERCENT_454 - 1.0f;
		readLen   = 2.0f * (readLen   - MIN_LN_READ_LEN_454)   / RANGE_LN_READ_LEN_454   - 1.0f;
		refLen    = 2.0f * (refLen    - MIN_LN_REF_LEN_454)    / RANGE_LN_REF_LEN_454    - 1.0f;
	}
}
