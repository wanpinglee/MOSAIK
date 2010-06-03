// ***************************************************************************
// CColorspaceUtilities - conversion to colorspace from basespace & vice versa
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "ColorspaceUtilities.h"

// constructor
CColorspaceUtilities::CColorspaceUtilities(unsigned int nAllowedMismatch) 
: mpBsRefSeqs(NULL)
, mNAllowedMismatch(nAllowedMismatch)
{
	// initialize our basespace and colorspace conversion maps
	InitializeBasespaceMap();
	InitializeColorspaceMap();
}

// constructor
CColorspaceUtilities::CColorspaceUtilities(void) 
: mpBsRefSeqs(NULL)
, mNAllowedMismatch(3)
{
	// initialize our basespace and colorspace conversion maps
	InitializeBasespaceMap();
	InitializeColorspaceMap();
}

// destructor
CColorspaceUtilities::~CColorspaceUtilities(void) {
	//mNumGapsObserved = 0;
	//delete [] mNumGapsObserved;

	mCsAl.csAlignmentLength = 0;
	if (mCsAl.csReference) delete [] mCsAl.csReference;
	if (mCsAl.csQuery)     delete [] mCsAl.csQuery;
	if (mCsAl.bsReference) delete [] mCsAl.bsReference;
	if (mCsAl.bsQuery)     delete [] mCsAl.bsQuery;
	if (mCsAl.type)        delete [] mCsAl.type;
	
	mCsAl.nDashReference = 0;
	mCsAl.nDashQuery     = 0;
	if (mCsAl.dashReference) delete [] mCsAl.dashReference;
	if (mCsAl.dashQuery)     delete [] mCsAl.dashQuery;

	mCsAl.nMismatch = 0;
	if (mCsAl.mismatch)    delete [] mCsAl.mismatch;

	mCsAl.nIdentical = 0;
	if (mCsAl.identical)   delete [] mCsAl.identical;
}

// converts the supplied alignment from colorspace to basespace
void CColorspaceUtilities::ConvertAlignmentToBasespace(Alignment& al) {

	
	// convert the alignment to character arrays
	const unsigned int pairwiseLen = al.Reference.Length();
	//char* pReference = al.Reference.Data();
	//char* pQuery     = al.Query.Data();

	// re-allocate mBsRef & mBsQuery if the reversed space is insufficient
	if(  pairwiseLen > mCsAl.csAlignmentLength ) {
		
		if ( mCsAl.csReference ) delete [] mCsAl.csReference;
		if ( mCsAl.csQuery )     delete [] mCsAl.csQuery;
		if ( mCsAl.bsReference ) delete [] mCsAl.bsReference;
		if ( mCsAl.bsQuery )     delete [] mCsAl.bsQuery;
		if ( mCsAl.type )        delete [] mCsAl.type;
	
		if ( mCsAl.dashReference ) delete [] mCsAl.dashReference;
		if ( mCsAl.dashQuery )     delete [] mCsAl.dashQuery;

		if ( mCsAl.mismatch )      delete [] mCsAl.mismatch;

		if ( mCsAl.identical )     delete [] mCsAl.identical;
		
		try {
			mCsAl.csAlignmentLength = pairwiseLen;
			mCsAl.csReference = new char [ pairwiseLen ];
			mCsAl.csQuery     = new char [ pairwiseLen ];
			mCsAl.bsReference = new char [ pairwiseLen + 2 ];
		        mCsAl.bsQuery     = new char [ pairwiseLen + 2 ];
			mCsAl.type        = new unsigned short [pairwiseLen];

			mCsAl.dashReference = new RegionT [ pairwiseLen ];
			mCsAl.dashQuery     = new RegionT [ pairwiseLen ];

			mCsAl.mismatch      = new unsigned int [ pairwiseLen ];

			mCsAl.identical     = new RegionT [ pairwiseLen ];

		}
		catch( bad_alloc ){
		        cout << "ERROR: Unable to allocate enough memory for converting colorspace ." << endl;
			exit(1);
		}

	}

	// initialize the counters
	mCsAl.nDashReference = 0;
	mCsAl.nDashQuery     = 0;
	mCsAl.nMismatch      = 0;
	mCsAl.nIdentical     = 0;
	

	// convert cs to bs
	char bsBase = mpBsRefSeqs[al.ReferenceIndex][al.ReferenceBegin];
	memcpy ( mCsAl.csReference, al.Reference.Data(), pairwiseLen );
	memcpy ( mCsAl.csQuery,     al.Query.Data(),     pairwiseLen );

	
	mCsAl.bsReference[0] = bsBase;
	mCsAl.bsQuery[0]     = bsBase;
	ConvertCs2Bs(mCsAl.csReference, mCsAl.bsReference, 0, pairwiseLen-1, bsBase);
	ConvertCs2Bs(mCsAl.csQuery, mCsAl.bsQuery, 0, pairwiseLen-1, bsBase);


	// search the dash regions & mismatches
	bool continuedDReference = false;
	bool continuedDQuery     = false;
	unsigned int nMatch      = 0;
	BS_MAP_t::const_iterator bsIter;
	for ( unsigned int i = 0; i < pairwiseLen; i++ ) {
		
		// determine identical region
		const bool isEndIdentity = ( mCsAl.csQuery[i] != mCsAl.csReference[i] ) && ( nMatch >= mNAllowedMismatch );
		if ( isEndIdentity ) {
			mCsAl.identical[ mCsAl.nIdentical ].Begin  = i - nMatch;
			mCsAl.identical[ mCsAl.nIdentical ].Length = nMatch;
			mCsAl.nIdentical++;
		}

		if ( mCsAl.csQuery[i] == mCsAl.csReference[i])
			nMatch++;
		else 
			nMatch = 0;
		
		// determine mismatches
		bool isN = (mCsAl.csReference[i] != 'A') && (mCsAl.csReference[i] != 'C') && (mCsAl.csReference[i] != 'G') && (mCsAl.csReference[i] != 'T');
		const bool isMismatch = (mCsAl.csReference[i] != '-') && (mCsAl.csQuery[i] != '-') && !isN && (mCsAl.csReference[i] != mCsAl.csQuery[i]);
		if ( isMismatch ) {
			// set the position
			mCsAl.mismatch[mCsAl.nMismatch] = i;
			mCsAl.nMismatch++;
			// set the mismatch flag
			mCsAl.type[i] = 0;
		}

		
		// for reference convertion
		if ( mCsAl.csReference[i] != '-' ) {
			// end the current dash region
			if ( continuedDReference ) {
				mCsAl.nDashReference++;
				// the current position could be a mismatch
				mCsAl.mismatch[mCsAl.nMismatch] = i;
				mCsAl.nMismatch++;
				mCsAl.type[i] = 1;
			}
			continuedDReference = false;
		}
		else {
			// start a dash region
			if ( !continuedDReference ) {
				mCsAl.dashReference[mCsAl.nDashReference].Begin  = i;
				mCsAl.dashReference[mCsAl.nDashReference].Length = 0;
				// the preceding position could be a mismatch
				mCsAl.mismatch[mCsAl.nMismatch] = i;
				mCsAl.nMismatch++;
				mCsAl.type[i] = 3;
			}
			mCsAl.dashReference[mCsAl.nDashReference].Length++;
			continuedDReference = true;
		}

		// for query convertion
		if ( mCsAl.csQuery[i] != '-' ) {
			// end the current dash region
			if ( continuedDQuery ) {
				mCsAl.nDashQuery++;
				// the current position could be a mismatch
				mCsAl.mismatch[mCsAl.nMismatch] = i;
				mCsAl.nMismatch++;
				mCsAl.type[i] = 2;
			}
			continuedDQuery = false;
		}
		else {
			// start a dash region
			if ( !continuedDQuery ) {
				mCsAl.dashQuery[mCsAl.nDashQuery].Begin  = i;
				mCsAl.dashQuery[mCsAl.nDashQuery].Length = 0;
				// the preceding position could be a mismatch
				mCsAl.mismatch[mCsAl.nMismatch] = i;
				mCsAl.nMismatch++;
				mCsAl.type[i] = 4;
			}
			mCsAl.dashQuery[mCsAl.nDashQuery].Length++;
			continuedDQuery = true;
		}


	}

	if ( nMatch > 0 ) {
		mCsAl.identical[ mCsAl.nIdentical ].Begin  = pairwiseLen - nMatch;
		mCsAl.identical[ mCsAl.nIdentical ].Length = nMatch;
		mCsAl.nIdentical++;
	}
	
	
	if ( mCsAl.identical[mCsAl.nIdentical - 1].Begin != 0 ) {
		// find sequencing errors
		if ( mCsAl.nMismatch > 0 )
			FindSequencingError(pairwiseLen);

        	if ( mCsAl.nDashReference > 0 ) {
			for ( unsigned int i = 0; i < mCsAl.nDashReference; i++ ) {
				unsigned int curPosition = mCsAl.dashReference[i].Begin + mCsAl.dashReference[i].Length;
				if ( mCsAl.bsReference[ curPosition ] != mCsAl.bsQuery[ curPosition ] ) {
					curPosition = mCsAl.dashReference[i].Begin;
					for ( unsigned int j = 0; j < mCsAl.dashReference[i].Length; j++ )
						mCsAl.bsQuery[ curPosition + j + 1 ] = 'N';
				}
			}

			AdjustDash(mCsAl.csReference, mCsAl.csQuery, mCsAl.dashReference, mCsAl.nDashReference, mCsAl.bsReference);
		}
		
	        if ( mCsAl.nDashQuery > 0 )
		        AdjustDash(mCsAl.csQuery, mCsAl.csReference, mCsAl.dashQuery, mCsAl.nDashQuery, mCsAl.bsQuery);

		
		// deal with the indentical region
		for ( unsigned int i = 0; i < mCsAl.nIdentical; i++ ) {
			unsigned int csEnd = mCsAl.identical[i].Begin + mCsAl.identical[i].Length - 1;
			unsigned int curPosition = mCsAl.identical[i].Begin;
			bsBase = mCsAl.bsReference[ curPosition ];
			while ( ( bsBase == '-' ) || ( bsBase == 'N' ) ) {
				curPosition--;
				bsBase = mCsAl.bsReference[ curPosition ];
				if ( curPosition == 0 )
					break;
			}

			bool isGoodBsBase = false;
			if ( ( bsBase != '-' ) && ( bsBase != 'N' ) )
				isGoodBsBase = true;
			
			if ( isGoodBsBase )
				ConvertCs2Bs(mCsAl.csQuery, mCsAl.bsQuery, mCsAl.identical[i].Begin, csEnd, bsBase);
		}
		

	}

	// end up the sequences
	mCsAl.bsReference[ pairwiseLen + 1 ] = 0;
	mCsAl.bsQuery[ pairwiseLen + 1 ]     = 0;

        al.Reference = mCsAl.bsReference;
	al.Query     = mCsAl.bsQuery;
	
	++al.ReferenceEnd;
	++al.QueryEnd;
	al.QueryLength = al.QueryEnd - al.QueryBegin + 1;

        
	// ------------------------------------------------------------------------------------------
	// convert the colorspace transition qualities to base qualities
	// NOTE: this algorithm will simply take the minimum of the two qualities that overlap a base
	// ------------------------------------------------------------------------------------------
	const unsigned short numColorspaceQualities = al.BaseQualities.Length();
	const unsigned short lastCSQIndex = numColorspaceQualities - 1;

	CMosaikString csQualities = al.BaseQualities;
	al.BaseQualities.Reserve(numColorspaceQualities + 1);
	al.BaseQualities.SetLength(numColorspaceQualities + 1);

	const char* pCSQual = csQualities.CData();
	char* pBSQual       = al.BaseQualities.Data();

	// handle the first base quality
	*pBSQual = *pCSQual;
	++pBSQual;

	// handle the internal base qualities
	for(unsigned short i = 1; i < numColorspaceQualities; ++i, ++pBSQual)
		*pBSQual = min(pCSQual[i - 1], pCSQual[i]);

	// handle the final base quality
           *pBSQual = pCSQual[lastCSQIndex];



	// update the number of mismatches
	// TODO: This should be augmented to support IUPAC ambiguity codes
	const unsigned int bsPairwiseLen = al.Reference.Length();

	al.NumMismatches = 0;
	for(unsigned short i = 0; i < bsPairwiseLen; ++i) {
		if(mCsAl.bsReference[i] != mCsAl.bsQuery[i]) al.NumMismatches++;
	}

}

// convert cs sequence to bs sequence
// a '-' converter would produce a '-'
void CColorspaceUtilities::ConvertCs2Bs (const char* csSequence, char* bsSequence, const unsigned int start, const unsigned int end, const char startBase) {

	char lastQueryBase = startBase;
	//bsSequence[ start ] = startBase;
	
	BS_MAP_t::const_iterator bsIter;
	// Note: the bs sequence has one more char than cd sequence
	for ( unsigned int i = start; i < end + 1; i++ ) {
		if ( csSequence[i] != '-' ) {
			bsIter = mBSMap.find(PACK_SHORT(lastQueryBase, csSequence[i]));
			if(bsIter == mBSMap.end()) {
		        	printf("ERROR: Unknown combination found when converting to basespace: [%c] & [%c]\n", lastQueryBase, csSequence[i]);
				exit(1);
			}
			bsSequence[ i + 1 ] = bsIter->second;
			lastQueryBase   = bsIter->second;
		}
		else
			bsSequence[ i + 1 ] = '-';
	}
}



// adjust positions of insertions or deletion
void CColorspaceUtilities::AdjustDash(const char* csSequence, const char* csSequenceOpp, const RegionT* dashRegion, const unsigned int nDashRegion, char* bsSequence) {

	// a dash appears after a '-' translation
	// sometimes, it should be in front of a '-' translation
	for ( unsigned int i = 0; i < nDashRegion; i++ ) {

		// if a, only one, cs match between two dash regions, the succeeding dash region would be moved one prevous position.
		bool hasPreviousDashRegion = false;
		if ( i > 0 )
		        hasPreviousDashRegion = ( dashRegion[i - 1].Begin + dashRegion[i - 1].Length ) == ( dashRegion[i].Begin - 1 );
			
		//unsigned short position = dashRegion[i].Begin;
		//bool isFrontDash = hasNextDashRegion || snpDRegion[i] || ( csSequence[position] != csSequenceOpp[position] );
		
		bool isFrontDash = hasPreviousDashRegion;
		if ( isFrontDash ) {
			
			unsigned short position = dashRegion[ i ].Begin;
			char ch  = bsSequence[ position ];
			
			bsSequence[ position ] = '-';

			position = dashRegion[ i ].Begin + dashRegion[ i ].Length;
			bsSequence[ position ] = ch;
		}
	}
}


// detect sequencing errors
void CColorspaceUtilities::FindSequencingError(const unsigned int pairwiseLen) {

	
	for (unsigned int i = 0; i < mCsAl.nMismatch; i++) {
	
		unsigned short curPosition = mCsAl.mismatch[ i ];
		// Assumption: the mismatch before or after a dash region couldn't be a sequencing error
		if ( mCsAl.type[ curPosition ] != 0 )
		        continue;

		// the mismatch doesn't have the same basespace characters in reference and query
		if ( mCsAl.bsReference[ curPosition ] != mCsAl.bsQuery[ curPosition ] )
		        continue;

		// try to find the end of SNPs
		bool isSnp = false;
		unsigned int nSnp = 0;
		
		for (unsigned int j = i + 1; j < mCsAl.nMismatch; j++) {
		        unsigned short nextPosition = mCsAl.mismatch[ j ];
			char nextBsReference = mCsAl.bsReference[ nextPosition + 1 ];
			char nextBsQuery     = mCsAl.bsQuery[ nextPosition + 1 ];

			// the # of SNPS is larger than the given # of mismatchs
			// in this case, the mismatch is determined as a sequence error
			nSnp = mCsAl.mismatch[ j ] - mCsAl.mismatch[ i ];
			if ( nSnp >  mNAllowedMismatch ) {
				i = j - 1;
			
				isSnp = false;
			
				break;
			}
			
			// find SNPs
			if ( nextBsReference == nextBsQuery ) {
				nSnp = mCsAl.mismatch[ j ] - mCsAl.mismatch[ i ];
				i = j - 1;
				
				isSnp = true;
				
				//if ( mCsAl.type[ nextPosition ]   == 1 )
				//	mCsAl.type[ nextPosition ] = 5;
				
				//if ( mCsAl.type[ nextPosition ]   == 2 )
				//	mCsAl.type[ nextPosition ] = 6;

				if ( mCsAl.type[ nextPosition ]   == 3 )
					mCsAl.type[ nextPosition ] = 7;
				
				if ( mCsAl.type[ nextPosition ]   == 4 )
					mCsAl.type[ nextPosition ] = 8;
				break;
			}

			// when meeting a beginning of dash region
			// we get a chance to look at one more
			if ( (mCsAl.type[ j ] == 3) || (mCsAl.type[ j ] == 4) ) {
				nextPosition = mCsAl.mismatch[ j + 1 ];
				nextBsReference = mCsAl.bsReference[ nextPosition + 1 ];
				nextBsQuery     = mCsAl.bsQuery[ nextPosition + 1 ];

				if ( nextBsReference == nextBsQuery ) {
					// the length of the dash region shouldn't be counted
					nSnp = mCsAl.mismatch[ j ] - curPosition;
					i = j;

					isSnp = true;
					
					if ( mCsAl.type[ nextPosition ]   == 1 )
						mCsAl.type[ nextPosition ] = 5;

					if ( mCsAl.type[ nextPosition ]   == 2 )
						mCsAl.type[ nextPosition ] = 6;
				}

				// SNPs should end before the dash region
				break;
			}

		}


		// the current position is a sequencing error
		bool isN = (mCsAl.csReference[ curPosition ] != 'A') && (mCsAl.csReference[ curPosition ] != 'C') && (mCsAl.csReference[ curPosition ] != 'G') && (mCsAl.csReference[ curPosition ] != 'T');
		if ( !isSnp && !isN) {
		        
			mCsAl.type[ curPosition ] = 9;
			
			// correct the error
			mCsAl.csQuery[ curPosition ] = mCsAl.csReference[ curPosition ];
			
			unsigned int curPosition2 = curPosition;
			char lastQueryBase = mCsAl.bsReference[ curPosition2 ];
			while ( ( lastQueryBase == '-' ) || ( lastQueryBase == 'N' ) ) {
				curPosition2--;
				lastQueryBase = mCsAl.bsReference[ curPosition2 ];
				if ( curPosition2 == 0 )
					break;
			}

			bool isGoodBsBase = false;
			if ( ( lastQueryBase != '-' ) && ( lastQueryBase != 'N' ) )
				isGoodBsBase = true;
			
			if ( isGoodBsBase )
				ConvertCs2Bs(mCsAl.csQuery, mCsAl.bsQuery, curPosition, pairwiseLen - 1, lastQueryBase);
			
		} // end of if ( !isSnp )

		// if the number of snps is larger than 2
		// we use Ns to present the SNPs
		if ( isSnp && (nSnp > 2) ) {
			for ( unsigned int i = 0; i < nSnp; i++ ) {
				mCsAl.bsQuery[ curPosition + i + 1 ] = 'N';
			}
		}
		

	}

}


// converts a colorspace sequence with provided seed base into basespace
void CColorspaceUtilities::ConvertColorspaceToBasespace(char seed, const string& colorspaceSeq, string& basespaceSeq) {

	// make the basespace sequence just as long as the colorspace sequence
	const unsigned short csLen = colorspaceSeq.size();
	basespaceSeq.resize(csLen);

	// create traversal pointers to our strings
	const char* pCS = colorspaceSeq.data();
	char* pBS       = (char*)basespaceSeq.data();

	// convert each colorspace/seed combo into a basespace nucleotide
	BS_MAP_t::const_iterator bsIter;
	for(unsigned int i = 0; i < csLen; ++i, ++pCS, ++pBS) {

		// find the appropriate seed/colorspace transition combination
		bsIter = mBSMap.find(PACK_SHORT(seed, *pCS));
		if(bsIter == mBSMap.end()) {
			printf("ERROR: Unknown combination found when converting to basespace: [%c] & [%c]\n", seed, *pCS);
			exit(1);
		}

		seed = *pBS = bsIter->second;
	}
}

// converts the supplied read from basespace to pseudo-colorspace
void CColorspaceUtilities::ConvertReadBasespaceToPseudoColorspace(CMosaikString& s) {

	char* pPrev   = s.Data();
	char* pString = pPrev + 1;

	// simplify various ambiguity codes
	*pPrev = GetSimplifiedBase(*pPrev);

	CS_MAP_t::const_iterator csIter;
	for(unsigned int i = 1; i < s.Length(); ++i, ++pString, ++pPrev) {

		// simplify various ambiguity codes
		*pString = GetSimplifiedBase(*pString);

		csIter = mCSMap.find(PACK_SHORT(*pPrev, *pString));
		if(csIter == mCSMap.end()) {
			printf("ERROR: Unknown combination found when converting to colorspace: [%c] & [%c]\n", *pPrev, *pString);
			exit(1);
		}

		*pPrev = csIter->second;
	}

	// adjust the read
	s.TrimEnd(1);
}

// converts the supplied read from colorspace to pseudo-colorspace
void CColorspaceUtilities::ConvertReadColorspaceToPseudoColorspace(CMosaikString& s) {
	char* pBases = s.Data();
	for(unsigned int i = 0; i < s.Length(); ++i, ++pBases) {
		switch(*pBases) {
			case '0':
				*pBases = 'A';
				break;
			case '1':
				*pBases = 'C';
				break;
			case '2':
				*pBases = 'G';
				break;
			case '3':
				*pBases = 'T';
				break;
			case 'X':
				break;
			case '-':
				*pBases = 'N';
				break;
			case '.':
				// here we pick an arbitrary colorspace transition, this will have at
				// least 25 % of being correct as opposed to specifying an 'N'.
				*pBases = 'A';
				break;
			default:
				printf("ERROR: Unrecognized nucleotide (%c) when converting read to pseudo-colorspace.\n", pBases[i]);
				exit(1);
				break;
		}
	}
}

// converts the supplied read from pseudo-colorspace to colorspace
// The function is used by the unaligned-read writer.
void CColorspaceUtilities::ConvertReadPseudoColorspaceToColorspace(CMosaikString& s) {
	char* pBases = s.Data();
	for(unsigned int i = 0; i < s.Length(); ++i, ++pBases) {
		switch(*pBases) {
			case 'A':
				*pBases = '0';
				break;
			case 'C':
				*pBases = '1';
				break;
			case 'G':
				*pBases = '2';
				break;
			case 'T':
				*pBases = '3';
				break;
			case 'X':
			case 'N':
				break;
			default:
				printf("ERROR: Unrecognized nucleotide (%c) when converting read to colorspace.\n", pBases[i]);
				exit(1);
				break;
		}
	}
}

// records regions of contiguous identity in the alignment
void CColorspaceUtilities::FindIndenticalRegions(char* pReference, char* pQuery, const unsigned short pairwiseLen, RegionVector& rv) {

	for(unsigned short i = 0; i < pairwiseLen; i++) {
		if(pReference[i] == pQuery[i]) {
			RegionT r(i);
			unsigned short end = i;
			while((end < pairwiseLen) && (pReference[end] == pQuery[end])) {
				++r.Length;
				++end;
			}
			rv.push_back(r);
			i = end;
		}
	}
}

// adds the pseudo-colorspace to basespace conversions
void CColorspaceUtilities::InitializeBasespaceMap(void) {

	// modified by Lee on 1/19/2009
	
	mBSMap[PACK_SHORT('A','A')] = 'A';
	mBSMap[PACK_SHORT('A','C')] = 'C';
	mBSMap[PACK_SHORT('A','G')] = 'G';
	mBSMap[PACK_SHORT('A','T')] = 'T';
	mBSMap[PACK_SHORT('A','E')] = 'N';
	//mBSMap[PACK_SHORT('A','M')] = '-';

	mBSMap[PACK_SHORT('C','A')] = 'C';
	mBSMap[PACK_SHORT('C','C')] = 'A';
	mBSMap[PACK_SHORT('C','G')] = 'T';
	mBSMap[PACK_SHORT('C','T')] = 'G';
	mBSMap[PACK_SHORT('C','F')] = 'N';
	//mBSMap[PACK_SHORT('C','U')] = '-';

	mBSMap[PACK_SHORT('G','A')] = 'G';
	mBSMap[PACK_SHORT('G','C')] = 'T';
	mBSMap[PACK_SHORT('G','G')] = 'A';
	mBSMap[PACK_SHORT('G','T')] = 'C';
	mBSMap[PACK_SHORT('G','I')] = 'N';
	//mBSMap[PACK_SHORT('G','D')] = '-';

	mBSMap[PACK_SHORT('T','A')] = 'T';
	mBSMap[PACK_SHORT('T','C')] = 'G';
	mBSMap[PACK_SHORT('T','G')] = 'C';
	mBSMap[PACK_SHORT('T','T')] = 'A';
	mBSMap[PACK_SHORT('T','L')] = 'N';
	//mBSMap[PACK_SHORT('T','B')] = '-';

	//mBSMap[PACK_SHORT('-','M')] = 'A';
	//mBSMap[PACK_SHORT('-','U')] = 'C';
	//mBSMap[PACK_SHORT('-','D')] = 'G';
	//mBSMap[PACK_SHORT('-','B')] = 'T';
	//mBSMap[PACK_SHORT('-','R')] = 'N';
	//mBSMap[PACK_SHORT('-','K')] = '-';
	
	mBSMap[PACK_SHORT('N','E')] = 'A';
	mBSMap[PACK_SHORT('N','F')] = 'C';
	mBSMap[PACK_SHORT('N','I')] = 'G';
	mBSMap[PACK_SHORT('N','L')] = 'T';
	mBSMap[PACK_SHORT('N','O')] = 'N';
	//mBSMap[PACK_SHORT('N','R')] = '-';
}

// adds the basespace to pseudo-colorspace conversions
void CColorspaceUtilities::InitializeColorspaceMap(void) {

	// modified by Lee on 1/19/2009
	
	mCSMap[PACK_SHORT('A','A')] = 'A';
	mCSMap[PACK_SHORT('A','C')] = 'C';
	mCSMap[PACK_SHORT('A','G')] = 'G';
	mCSMap[PACK_SHORT('A','T')] = 'T';
	mCSMap[PACK_SHORT('A','N')] = 'E';
	//mCSMap[PACK_SHORT('A','-')] = 'M';

	mCSMap[PACK_SHORT('C','A')] = 'C';
	mCSMap[PACK_SHORT('C','C')] = 'A';
	mCSMap[PACK_SHORT('C','G')] = 'T';
	mCSMap[PACK_SHORT('C','T')] = 'G';
	mCSMap[PACK_SHORT('C','N')] = 'F';
	//mCSMap[PACK_SHORT('C','-')] = 'U';

	mCSMap[PACK_SHORT('G','A')] = 'G';
	mCSMap[PACK_SHORT('G','C')] = 'T';
	mCSMap[PACK_SHORT('G','G')] = 'A';
	mCSMap[PACK_SHORT('G','T')] = 'C';
	mCSMap[PACK_SHORT('G','N')] = 'I';
	//mCSMap[PACK_SHORT('G','-')] = 'D';

	mCSMap[PACK_SHORT('T','A')] = 'T';
	mCSMap[PACK_SHORT('T','C')] = 'G';
	mCSMap[PACK_SHORT('T','G')] = 'C';
	mCSMap[PACK_SHORT('T','T')] = 'A';
	mCSMap[PACK_SHORT('T','N')] = 'L';
	//mCSMap[PACK_SHORT('T','-')] = 'B';

	//mCSMap[PACK_SHORT('-','A')] = 'M';
	//mCSMap[PACK_SHORT('-','C')] = 'U';
	//mCSMap[PACK_SHORT('-','G')] = 'D';
	//mCSMap[PACK_SHORT('-','T')] = 'B';
	//mCSMap[PACK_SHORT('-','N')] = 'R';
	//mCSMap[PACK_SHORT('-','-')] = 'K';

	mCSMap[PACK_SHORT('N','A')] = 'E';
	mCSMap[PACK_SHORT('N','C')] = 'F';
	mCSMap[PACK_SHORT('N','G')] = 'I';
	mCSMap[PACK_SHORT('N','T')] = 'L';
	mCSMap[PACK_SHORT('N','N')] = 'O';
	//mCSMap[PACK_SHORT('N','X')] = 'R';
}

// replaces the gaps in the pairwise alignment with a dibase transition code '4'
// NOTE: we arbitrarily add an extra gap transition - this should be modified so
//       that we create a true transition to the next base.
void CColorspaceUtilities::PatchColorspaceGaps(char* pReference, char* pQuery, const unsigned short pairwiseLen) {

	const unsigned int lastIndex = pairwiseLen - 1;

	// here we take advantage of the fact that gaps should NEVER occur in the
	// reference and query sequence at the same position
	for(unsigned short i = 0; i < pairwiseLen; ++i) {

		// patch the gaps in the reference
		if(pReference[i] == '-') {
			unsigned int index = i;
			while((pReference[index] == '-') && (index < lastIndex)) index++;
			const unsigned int length = index - i + 1;
			memset(pReference + i, 'N', length);
			i = index;
		}

		// patch the gaps in the query
		if(pQuery[i] == '-') {
			unsigned int index = i;
			while((pQuery[index] == '-') && (index < lastIndex)) index++;
			const unsigned int length = index - i + 1;
			memset(pQuery + i, 'N', length);
			i = index;
		}
	}
}

// sets the reference sequences
void CColorspaceUtilities::SetReferenceSequences(char** pBsRefSeqs) {
	mpBsRefSeqs = pBsRefSeqs;
}

void CColorspaceUtilities::SetNumAllowedMismatch(unsigned int allowedMismatch) {
	mNAllowedMismatch = allowedMismatch;
}

