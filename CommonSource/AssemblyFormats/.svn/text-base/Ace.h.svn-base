// ***************************************************************************
// CAce - exports MOSAIK assemblies into the ACE format.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// This code is dual licenced under the GNU General Public License 2.0+ or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include "AbstractAssemblyFormat.h"
#include "FileUtilities.h"
#include "MemoryUtilities.h"

using namespace std;

#define ACE_LINE_LENGTH 50

class CAce : public CAbstractAssemblyFormat {
public:
	// constructor
	CAce(unsigned char referenceBaseQuality);
	// destructor
	~CAce(void);
	// opens the ace file
	void Open(const CMosaikString& filename);
	// closes the ace file
	void Close(void);
	// saves the specified reference sequence to the header
	void SaveHeader(CMosaikString& reference, const string& referenceName, const unsigned int ungappedRefLength, const uint64_t& numSequences);
	// saves the specified read to the index and reads files
	void SaveRead(const Alignment& al, CMosaikString& gappedRead);

private:
	// write a sequence to the specified file stream in rows of ACE_LINE_LENGTH bases
	static void ChopSequence(FILE* out, const CMosaikString& reference);
	// retrieves the current time in the format used by consed (slightly modified asctime)
	static void GetTime(string& time);
	// our file output streams
	FILE* mHeaderStream;
	FILE* mReadStream;
	// our temporary files
	string mTempHeaderFilename;
	string mTempReadFilename;
	// our time string
	string mTimeString;
	// our reference sequence base quality
	unsigned char mReferenceSequenceBaseQuality;
};
