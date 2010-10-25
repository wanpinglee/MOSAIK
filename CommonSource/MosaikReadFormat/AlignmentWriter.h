// ***************************************************************************
// CAlignmentWriter - stores reads in a MOSAIK alignment archive.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#ifndef _AlignmentWriter_H_
#define _AlignmentWriter_H_

#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include <cstdio>
#include "fastlz.h"
#include "AlignedRead.h"
#include "AlignmentStatus.h"
#include "FastLZIO.h"
#include "FileUtilities.h"
#include "GapInfo.h"
#include "MemoryUtilities.h"
#include "Mosaik.h"
#include "NaiveAlignmentSet.h"
#include "Read.h"
#include "ReadGroup.h"
#include "ReferenceSequence.h"
#include "SequenceUtilities.h"
#include "SequencingTechnologies.h"
#include "TimeSupport.h"
#include "UnorderedMap.h"

#define FASTLZ_BETTER_SPEED       1
#define FASTLZ_BETTER_COMPRESSION 2

#define ALS_IS_REVERSE_STRAND      1
#define ALS_IS_MATE_REVERSE_STRAND 2

// here we assume that we'll need space for Sanger length reads (reference bases,
// query bases, query base qualities)
#define MEMORY_BUFFER_SIZE        3000

#define NUM_READS_OFFSET 25

namespace MosaikReadFormat {
	class CAlignmentWriter {
	public:
		// constructor
		CAlignmentWriter(void);
		// destructor
		~CAlignmentWriter(void);
		// adds a header tag
		void AddHeaderTag(const Tag& tag);
		// closes the alignment archive
		void Close(void);
		// retrieves the number of bases written
		uint64_t GetNumBases(void) const;
		// retrieves the number of reads written
		uint64_t GetNumReads(void) const;
		// opens the alignment archive
		void Open(const string& filename, const vector<ReferenceSequence>& referenceSequences, const vector<ReadGroup>& readGroups, const AlignmentStatus as);
		// saves the read to the alignment archive
		void SaveAlignedRead(const Mosaik::AlignedRead& ar);
		// saves the alignment to the alignment archive
		void SaveAlignment(Alignment* pAl);
		// saves the read to the alignment archive
		void SaveRead(const Mosaik::Read& mr, CNaiveAlignmentSet& mate1Alignments, CNaiveAlignmentSet& mate2Alignments);
		// adds a header tag (only works before opening the file)
		void AddHeaderTag(const unsigned char tagID, const TagType& tagType);
		// set reference gaps vector
		void SetReferenceGaps(vector<unordered_map<unsigned int, unsigned short> >* pRefGapVector);
		// the lengths of references should + 1 
		// since there would be one more base after converting references from solid (colorspace) to basespace
		void AdjustSolidReferenceBases(void);
		// adjust the size of partition; the default is 20000
		void AdjustPartitionSize(unsigned short size);

	private:
		// specifies our index entry
		struct IndexEntry {
			off_type Offset;
			unsigned int Position;
			unsigned int ReferenceIndex;
		};
		// adjusts the buffer
		void AdjustBuffer(void);
		// serializes the specified alignment
		void WriteAlignment(const Alignment* pAl, const bool isLongRead, const bool isPairedEnd, const bool isFirstMate, const bool isResolvedAsPair);
		// write partition to disk
		void WritePartition(void);
		// write the read header to disk
		void WriteReadHeader(const CMosaikString& readName, const unsigned int readGroupCode, const unsigned char readStatus, const unsigned int numMate1Alignments, const unsigned int numMate2Alignments);
		// writes the tag to disk
		void WriteTag(const map<unsigned char, Tag>::const_iterator& htIter);
		// denotes the status of the output stream
		bool mIsOpen;
		// our compressed output stream
		FILE* mOutStream;
		// stores the number of sequences that have been written
		uint64_t mNumReads;
		uint64_t mNumBases;
		// our output buffer
		unsigned char* mBuffer;
		unsigned int mBufferLen;
		unsigned int mBufferPosition;
		unsigned int mBufferThreshold;
		// our output compression buffer
		unsigned char* mCompressionBuffer;
		unsigned int mCompressionBufferLen;
		// our partitioning setup
		unsigned short mPartitionSize;
		unsigned short mPartitionMembers;
		// our output filename
		string mOutputFilename;
		// our reference sequence vector
		vector<ReferenceSequence> mReferenceSequences;
		unsigned int mNumRefSeqs;
		// our reference sequence gap vector
		vector<unordered_map<unsigned int, unsigned short> >* mpRefGapVector;
		// our alignment status
		AlignmentStatus mStatus;
		// denotes that this alignment archive is paired-end (used in SaveRead)
		bool mIsPairedEndArchive;
		// our block index
		vector<IndexEntry> mIndex;
		unsigned int mLastReferenceIndex;
		unsigned int mLastReferencePosition;
		bool mStoreIndex;
		// our header tags
		map<unsigned char, Tag> mHeaderTags;
	};
}

#endif
