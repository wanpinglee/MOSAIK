// ***************************************************************************
// CReadReader - loads reads from the MOSAIK read archive.
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
#include "MemoryUtilities.h"
#include "ReadGroup.h"
#include "Read.h"
#include "ReadStatus.h"
#include "SequencingTechnologies.h"
#include "SequenceUtilities.h"
#include "TimeSupport.h"

using namespace std;

namespace MosaikReadFormat {
	class CReadReader {
	public:
		// constructor
		CReadReader(void);
		// destructor
		~CReadReader(void);
		// validates the supplied read archive file
		static bool CheckFile(const string& filename, SequencingTechnologies& st, ReadStatus& rs, const bool showError);
		// closes the read archive
		void Close(void);
		// gets the metadata object
		ReadGroup GetReadGroup(void) const;
		// gets the archive read count
		uint64_t GetNumReads(void) const;
		// gets the read archive sequencing technology
		SequencingTechnologies GetSequencingTechnology(void) const;
		// gets the read archive status
		ReadStatus GetStatus(void) const;
		// loads the next read from the read archive
		bool LoadNextRead(Mosaik::Read& mr);
		// opens the read archive
		void Open(const string& filename);
		// sets the file pointer to the beginning of the read data
		void Rewind(void);

	private:
		// denotes the status of the output stream
		bool mIsOpen;
		// our compressed output stream
		FILE* mInStream;
		// stores the archive read and base count
		uint64_t mNumReads;
		uint64_t mNumBases;
		// stores the current read number
		uint64_t mCurrentRead;
		// our input buffer
		unsigned char* mBuffer;
		unsigned char* mBufferPtr;
		unsigned int mBufferLen;
		// our input compression buffer
		unsigned char* mCompressionBuffer;
		unsigned int mCompressionBufferLen;
		// our output filename
		string mInputFilename;
		// our partitioning setup
		unsigned short mPartitionSize;
		unsigned short mPartitionMembers;
		// read archive status
		ReadStatus mStatus;
		// our AB SOLiD flag
		bool mIsSOLiD;
		// our reads offset
		off_type mReadsOffset;
		// our metadata
		ReadGroup mReadGroup;
	};
}
