// ***************************************************************************
// CReferenceSequenceReader - loads reference sequences from the MOSAIK 
//                            reference sequence archive.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "FastLZIO.h"
#include "Mosaik.h"
#include "LargeFileSupport.h"
#include "ReferenceSequence.h"
#include "ReferenceSequenceStatus.h"
#include "UnorderedMap.h"

using namespace std;

namespace MosaikReadFormat {
	class CReferenceSequenceReader {
	public:
		// constructor
		CReferenceSequenceReader(void);
		// destructor
		~CReferenceSequenceReader(void);
		// checks to see if this is truly a MOSAIK reference sequence archive
		static bool CheckFile(const string& filename, const bool showError);
		// closes the reference sequence archive
		void Close(void);
		// copies the reference sequences from this archive into the supplied character array
		// NOTE: caller frees the memory
		void CopyReferenceSequences(char** &pSeqs);
		// returns the number of reference sequences in this archive
		unsigned int GetNumReferenceSequences(void) const;
		// returns the reference sequence length
		unsigned int GetReferenceSequenceLength(void) const;
		// retrieves the desired reference sequence and places it in the specified string
		void GetReferenceSequence(const string& name, string& bases);
		// adds the reference sequences to the supplied vector
		void GetReferenceSequences(vector<ReferenceSequence>& referenceSequences);
		// returns the reference sequence status
		ReferenceSequenceStatus GetStatus(void) const;
		// returns true if the reference sequences in this archive match those from the supplied vector
		bool HasSameReferenceSequences(vector<ReferenceSequence>& otherSeqs);
		// initializes the supplied pointer with the concatenated reference sequence
		void LoadConcatenatedSequence(char* &referenceSequence);
		//initializes the supplied pointer with the 2-bit concatenated reference sequence
		void Load2BitConcatenatedSequence(char* &referenceSequence, char* &maskSequence, unsigned int& numMaskedPositions);
		// opens the reference sequence archive
		void Open(const string& filename);

	private:
		// define a comparison function for sorting our alignment positions (ascending)
		struct SortReferenceSequencesByBeginAsc {
			bool operator()(const ReferenceSequence& ar1, const ReferenceSequence& ar2) {
				return ar1.Begin < ar2.Begin;
			}
		};
		// stores the file state
		bool mIsOpen;
		// our input file stream
		FILE* mInStream;
		// our offsets
		off_type mConcatenatedOffset;
		off_type mConcatenated2bOffset;
		off_type mIndexOffset;
		off_type mReferenceBasesOffset;
		off_type mMaskedRegionsOffset;
		// the number of reference sequences contained in the reference archive
		unsigned int mNumReferenceSequences;
		// the concatenated sequence length
		unsigned int mConcatenatedLen;
		// the concatenated 2-bit sequence length
		unsigned int mConcatenated2bLen;
		// our file status
		ReferenceSequenceStatus mStatus;
		// stores the index for our reference sequences
		unordered_map<string, ReferenceSequence> mIndex;
		// our FastLZ IO reader and writer
		CFastLZIO mFIO;
	};
}
