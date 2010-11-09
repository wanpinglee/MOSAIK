// ***************************************************************************
// CAlignmentReader - loads alignments from the MOSAIK alignment archive.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <list>
#include <map>
#include <vector>
#include "fastlz.h"
#include "AlignedRead.h"
#include "AlignmentStatus.h"
#include "FastLZIO.h"
#include "GapInfo.h"
#include "LargeFileSupport.h"
#include "MemoryUtilities.h"
#include "Mosaik.h"
#include "ReadGroup.h"
#include "ReferenceSequence.h"
#include "SequenceUtilities.h"
#include "SequencingTechnologies.h"

using namespace std;

namespace MosaikReadFormat {
	class CAlignmentReader {
	public:
		// constructor
		CAlignmentReader(void);
		// destructor
		~CAlignmentReader(void);
		// checks to see if this is truly a MOSAIK alignment archive
		static bool CheckFile(const string& filename, SequencingTechnologies& st, AlignmentStatus& as, const bool showError);
		// closes the alignment archive
		void Close(void);
		// returns the a pointer to the header tags map
		map<unsigned char, Tag>* GetHeaderTags(void);
		// returns the number of bases in the archive
		uint64_t GetNumBases(void) const;
		// returns the number of reads in the archive
		uint64_t GetNumReads(void) const;
		// retrieves the read group data given a read group code
		ReadGroup GetReadGroupFromCode(const unsigned int code);
		// retrieves the read groups vector
		void GetReadGroups(vector<ReadGroup>& readGroups) const;
		// retrieves the reference sequence gaps
		vector<vector<GapInfo> >* GetReferenceSequenceGaps(void);
		// if the reference name is found, the referenceIndex will be set and the function will return true
		bool GetReferenceSequenceIndex(const string& referenceName, unsigned int& referenceIndex) const;
		// retrieves the reference sequence data
		vector<ReferenceSequence>* GetReferenceSequences(void);
		void GetReferenceSequences( vector<ReferenceSequence>& refVec);
		// gets the alignment archive sequencing technology
		SequencingTechnologies GetSequencingTechnology(void) const;
		// retrieves the file status
		AlignmentStatus GetStatus(void) const;
		// retrieves the signature
		void GetSignature ( char*& signature );
		// jumps to the block containing the specified reference index and position
		void Jump(const unsigned int referenceIndex, const unsigned int referencePosition);
		// loads the next alignment from the alignment archive
		bool LoadNextAlignment(Alignment& al);
		// loads the next read from the alignment archive
		bool LoadNextRead(Mosaik::AlignedRead& ar);
		// opens the alignment archive
		void Open(const string& filename);
		// sets the file pointer to the beginning of the read data
		void Rewind(void);

	private:
		// load the read header from disk
		void LoadReadHeader(CMosaikString& readName, unsigned int& readGroupCode, unsigned char& readStatus, unsigned int& numMate1Alignments, unsigned int& numMate2Alignments);
		// deserializes each alignment and stores them in the supplied vector
		void ReadAlignments(vector<Alignment>& alignments, const bool isLongRead, const bool isPairedInSequencing, const bool isResolvedAsPair, const unsigned int readGroupCode);
		// deserialize the alignment
		void ReadAlignment(Alignment& al, const bool isLongRead, const bool isPairedInSequencing, const bool isResolvedAsPair);
		// reads a new compressed partition (returns false if EOF occurs)
		bool ReadPartition(void);
		// reads the tag from disk
		void ReadTag(Tag& tag);
		// denotes the status of the output stream
		bool mIsOpen;
		// our compressed output stream
		FILE* mInStream;
		// stores the archive read count
		uint64_t mNumReads;
		uint64_t mNumBases;
		// stores the current read number
		uint64_t mCurrentRead;
		// stores the file offsets
		off_type mReadsOffset;
		off_type mReferenceGapOffset;
		off_type mIndexOffset;
		// our input buffer
		char* mBuffer;
		char* mBufferPtr;
		unsigned int mBufferLen;
		// our input compression buffer
		unsigned char* mCompressionBuffer;
		unsigned int mCompressionBufferLen;
		// our input filename
		string mInputFilename;
		// our partitioning setup
		unsigned short mPartitionSize;
		unsigned short mPartitionMembers;
		// our reference sequence LUT
		char** mRefSeqLUT;
		unsigned int mNumRefSeqs;
		// our reference sequences
		vector<ReferenceSequence> mReferenceSequences;
		// our reference sequence gap vector
		vector<vector<GapInfo> > mRefSeqGaps;
		// our read groups
		vector<ReadGroup> mReadGroups;
		// our file status
		AlignmentStatus mStatus;
		// our sequencing technology
		SequencingTechnologies mSeqTech;
		// our header tags
		map<unsigned char, Tag> mHeaderTags;
		// our read group LUT
		map<unsigned int, ReadGroup> mReadGroupLUT;
		// our file signature
		//static const char* MOSAIK_SIGNATURE;
		//static const unsigned char SIGNATURE_LENGTH;
		char* MosaikSignature;
	};
}
