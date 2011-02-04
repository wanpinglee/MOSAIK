// ***************************************************************************
// CColorspaceUtilities - conversion to colorspace from basespace & vice versa
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "Alignment.h"
#include "MosaikString.h"

using namespace std;

#define PACK_SHORT(a,b) (((a) << 8) | (b))
#define ARRAY_EXTENSION 10

// define our region data structure
struct RegionT {
	unsigned short Begin;
	unsigned short Length;

	// constructor
	RegionT(unsigned short beg)
		: Begin(beg)
		, Length(0)
	{}

	RegionT(void)
	        : Begin(0)
		, Length(0)
	{}
};

struct csAlignment {

	unsigned int csAlignmentLength;
	
	char* csReference;
	char* csQuery;
	char* bsReference;
	char* bsQuery;

	// 0: mismatch
	// 1: the position after a dash region in reference
	// 2: the position after a dash region in query
	// 3: the preceding position of a dash region in reference
	// 4: the preceding position of a dash region in query
	// 5-8: the position has been confirmed a SNP, whose original type is 1-4, respectively
	// 9: the position has been confirmed an error
	unsigned short* type;

	unsigned int nDashReference;
	unsigned int nDashQuery;
	RegionT* dashReference;
	RegionT* dashQuery;

	unsigned int  nMismatch;
	unsigned int* mismatch;

	unsigned int nIdentical;
	RegionT* identical;


	csAlignment(void)
		: csAlignmentLength(0)
		, csReference(NULL)
		, csQuery(NULL)
		, bsReference(NULL)
		, bsQuery(NULL)
		, type(NULL)
	
		, nDashReference(0)
		, nDashQuery(0)
		, dashReference(NULL)
		, dashQuery(NULL)

		, nMismatch(0)
		, mismatch(NULL)

		, nIdentical(0)
		, identical(NULL)
	{}

};

typedef map<short, char> CS_MAP_t;
typedef map<short, char> BS_MAP_t;
typedef vector<RegionT> RegionVector;

class CColorspaceUtilities {
public:
	// constructor
	CColorspaceUtilities(unsigned int nAllowedMismatch);
	CColorspaceUtilities(void);
	// destructor
	~CColorspaceUtilities(void);
	// converts the supplied alignment from colorspace to basespace
	void ConvertAlignmentToBasespace(Alignment& al);
	// converts the supplied read from basespace to pseudo-colorspace
	void ConvertReadBasespaceToPseudoColorspace(CMosaikString& s);
	// converts the supplied read from colorspace to pseudo-colorspace
	void ConvertReadColorspaceToPseudoColorspace(CMosaikString& s);
	// converts the supplied read from pseudo-colorspace to colorspace
	void ConvertReadPseudoColorspaceToColorspace(CMosaikString& s);
	// sets the reference sequences
	void SetReferenceSequences(char** pBsRefSeqs);
	// sets mNAllowedMismatch used for sequence errors
	void SetNumAllowedMismatch(unsigned int allowedMismatch);

private:
	// copy constructor 
	CColorspaceUtilities( const CColorspaceUtilities& copy );
	// assing operator
	CColorspaceUtilities& operator=( const CColorspaceUtilities& copy );
	// converts a colorspace sequence with provided seed base into basespace
	void ConvertColorspaceToBasespace(char seed, const string& colorspaceSeq, string& basespaceSeq);
	// records regions of contiguous identity in the alignment
	void FindIndenticalRegions(char* pReference, char* pQuery, const unsigned short pairwiseLen, RegionVector& rv);
	// returns the simplified version of the IUPAC ambiguity code - based on base frequencies in the human genome 36.3
	static inline char GetSimplifiedBase(char c);
	// adds the colorspace to basespace conversions
	void InitializeBasespaceMap(void);
	// adds the basespace to colorspace conversions
	void InitializeColorspaceMap(void);
	// replaces the gaps in the pairwise alignment with a dibase transition code '4'
	void PatchColorspaceGaps(char* pReference, char* pQuery, const unsigned short pairwiseLen);
	// updates the observed gaps array and returns true if gaps are found in the alignment
	bool UpdateObservedGaps(const char* pReference, const char* pQuery, const unsigned short pairwiseLen);
	// detect sequencing errors
	void FindSequencingError(unsigned int pairwiseLen);
	// adjust positions of insertions or deletion
	void AdjustDash(const char* csSequence, const char* csSequenceOpp, const RegionT* dashRegion, const unsigned int nDashRegion, char* bsSequence);
	// convert cs sequence to bs sequence
	void ConvertCs2Bs (const char* csSequence, char* bsSequence, const unsigned int start, const unsigned int end, const char startBase);
	// get the complement base
	inline char GetComplementBase(char c) ;
	// our colorspace to basespace conversion map
	BS_MAP_t mBSMap;
	// our basespace to colorspace conversion map
	CS_MAP_t mCSMap;
	// our basespace reference sequence
	char** mpBsRefSeqs;
	
	csAlignment mCsAl;
	unsigned int mNAllowedMismatch;
};

// returns the simplified version of the IUPAC ambiguity code - based on base frequencies in the human genome 36.3
inline char CColorspaceUtilities::GetSimplifiedBase(char c) {

	switch(c) {
		case 'M':
		case 'R':
		case 'V':
			c = 'A';
			break;

		case 'S':
			c = 'G';
			break;

		case 'B':
		case 'D':
		case 'H':
		case 'K':
		case 'W':
		case 'Y':
			c = 'T';
			break;

		// Marked by Lee on 1/19/2009
		//case 'N':
		//	c = 'X';
		//	break;
	}

	return c;
}

// get the complement base
inline char CColorspaceUtilities::GetComplementBase(char c) {
	
	switch(c) {
		case 'A':
			c = 'T';
			break;
		case 'C':
			c = 'G';
			break;
		case 'G':
			c = 'C';
			break;
		case 'T':
			c = 'A';
		default:
			break;
	}

	return c;
}
