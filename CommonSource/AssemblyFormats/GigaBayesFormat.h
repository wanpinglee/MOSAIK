// ***************************************************************************
// CGigaBayesFormat - exports MOSAIK assemblies into the GigaBayes (GIG)
//                    format.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <cstdio>
#include "AbstractAssemblyFormat.h"
#include "FileUtilities.h"
#include "LargeFileSupport.h"
#include "MemoryUtilities.h"

class CGigaBayesFormat : public CAbstractAssemblyFormat {
public:
	// constructor
	CGigaBayesFormat(unsigned char referenceBaseQuality);
	// destructor
	~CGigaBayesFormat(void);
	// checks if the base quality buffer is large enough to accomodate the requested size
	void CheckBQBufferSize(const unsigned int requestedLength);
	// opens the GIG file
	void Open(const CMosaikString& filename);
	// closes the GIG file
	void Close(void);
	// saves the specified reference sequence to the header
	void SaveHeader(CMosaikString& reference, const string& referenceName, const unsigned int ungappedRefLength, const uint64_t& numSequences);
	// saves the specified read to the index and reads files
	void SaveRead(const Alignment& al, CMosaikString& gappedRead);

private:
	// our file output streams
	FILE* mOutputStream;
	FILE* mReadIndexStream;
	// our output buffer
	unsigned char* mBuffer;
	unsigned int mBufferLen;
	// our base qualities buffer
	unsigned short* mBQBuffer;
	unsigned int mBQBufferLen;
	// our temporary file
	string mTempReadIndexFilename;
	// the contig name and offset
	string mContigName;
	// the locations where store the file offsets
	off_type mContigIndexOffsetLocus;
	off_type mReadIndexOffsetLocus;
	// our file offsets
	off_type mContigIndexOffset;
	off_type mContigFileOffset;
	off_type mReadIndexOffset;
	// our reference sequence base quality
	unsigned char mReferenceSequenceBaseQuality;
};
