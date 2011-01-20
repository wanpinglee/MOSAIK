// ***************************************************************************
// CReadWriter - stores reads in a MOSAIK read archive.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iostream>
#include <cstdio>
#include "fastlz.h"
#include "Mosaik.h"
#include "FileUtilities.h"
#include "ReadGroup.h"
#include "MemoryUtilities.h"
#include "Read.h"
#include "ReadStatus.h"
#include "SequencingTechnologies.h"
#include "SequenceUtilities.h"
#include "TimeSupport.h"

#define FASTLZ_BETTER_SPEED       1
#define FASTLZ_BETTER_COMPRESSION 2

// here we assume that we'll need space for Sanger length reads (reference bases,
// query bases, query base qualities)
#define MEMORY_BUFFER_SIZE        3000

#define UPDATE_HEADER_OFFSET 16

namespace MosaikReadFormat {
	class CReadWriter {
	public:
		// constructor
		CReadWriter(void);
		// destructor
		~CReadWriter(void);
		// closes the read archive
		void Close(void);
		// retrieves the number of bases written
		uint64_t GetNumBases(void) const;
		// retrieves the number of reads written
		uint64_t GetNumReads(void) const;
		// opens the read archive
		void Open(const string& filename, const ReadStatus rs, const ReadGroup& md);
		// saves the read to the read archive
		void SaveRead(const Mosaik::Read& mr);
	private:
		// adjusts the buffer
		void AdjustBuffer(void);
		// write partition to disk
		void WritePartition(void);
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
		// our AB SOLiD flag
		bool mIsSOLiD;
		// our reads offset
		off_type mReadsOffset;
	};
}
