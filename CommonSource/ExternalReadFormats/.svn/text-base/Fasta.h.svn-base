// ***************************************************************************
// CFasta - imports reads from the FASTA file format.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iostream>
#include <string>
#include <sstream>
#include "zlib.h"
#include "LargeFileSupport.h"
#include "Mosaik.h"
#include "MosaikString.h"
#include "MemoryUtilities.h"
#include "Read.h"
#include "SequenceUtilities.h"
#include "RegexUtilities.h"

using namespace std;

struct FastaTags {
	CMosaikString Name;
	CMosaikString Species;
	CMosaikString GenomeAssemblyID;
	CMosaikString URI;
};

class CFasta {
public:
	// constructor
	CFasta(void);
	// destructor
	~CFasta(void);
	// Checks if a file is truly a FASTA file
	static bool CheckFile(const string& filename, const bool showError);
	// closes the FASTA file
	void Close(void);
	// enables parsing of the base quality file
	void EnableBaseQualityFile(const string& filename);
	// loads the next read from the FASTA file
	bool LoadNextMate(FastaTags& ft, Mosaik::Mate& m);
	// opens the FASTA file
	void Open(const string& filename);
	// sets the file pointer to the beginning of the read data
	void Rewind(void);
	// sets the assigned base quality if a base quality file is not specified
	void SetAssignedBaseQuality(const unsigned char baseQuality);

private:
	// denotes the status of the output stream
	bool mIsOpen;
	// denotes if the file is compressed
	bool mAreBasesCompressed;
	bool mAreBaseQualitiesCompressed;
	// denotes if we have a base quality file
	bool mHasBaseQualityFile;
	// our compressed output stream
	FILE* mInStream;
	FILE* mInQualityStream;
	gzFile mInZStream;
	gzFile mInQualityZStream;
	// our input buffers
	char* mBaseBuffer;
	unsigned int mBaseBufferLen;
	char* mQualityBuffer;
	unsigned int mQualityBufferLen;
	// stores the start of the read data (handles csfasta case)
	off_type mReadDataBaseOffset;
	off_type mReadDataQualityOffset;
	// our base quality FASTA filename
	string mBaseQualityFilename;
	// our assigned base quality
	unsigned char mAssignedBaseQuality;
};
