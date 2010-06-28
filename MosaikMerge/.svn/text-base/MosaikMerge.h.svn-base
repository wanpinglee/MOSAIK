// ***************************************************************************
// CMosaikMerge - merges two or more sorted alignment archives.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <list>
#include <string>
#include <vector>
#include <cstdio>
#include "Alignment.h"
#include "AlignmentReader.h"
#include "AlignmentWriter.h"
#include "Mosaik.h"
#include "ConsoleUtilities.h"
#include "ProgressBar.h"
#include "ReadGroup.h"
#include "UnorderedMap.h"

using namespace std;

// define our serialization status flags
#define MERGE_IS_LONG_READ               1
#define MERGE_IS_REVERSE_COMPLEMENT      2
#define MERGE_IS_MATE_REVERSE_COMPLEMENT 4

class CMosaikMerge {
public:
	// constructor
	CMosaikMerge(const unsigned int numCachedAlignments);
	// destructor
	~CMosaikMerge(void);
	// merges the files contained in the file vector and stores them in a specified output file
	void MergeFiles(vector<string>& fileVector, string& outputFilename);
private:
	// retrieves an alignment from the specified temporary file and adds it to the specified list
	void AddAlignment(FILE* tempFile, const unsigned int owner, list<Alignment>& alignments);
	// retrieves an alignment from the specified temporary file
	bool GetAlignment(FILE* tempFile, const unsigned int owner, Alignment& al);
	// records the observed gaps in the specified reference 
	void RecordReferenceGaps(const unsigned short refIndex, const unsigned int refBegin, const CMosaikString& refSeq);
	// serializes the specified list to a temporary file
	uint64_t Serialize(list<Alignment>& alignmentCache, const unsigned int numEntries);
	// denotes the alignment cache size
	unsigned int mNumCachedAlignments;
	// our temporary file vector
	vector<string> mTempFiles;
	// our output buffer
	unsigned char* mBuffer;
	unsigned int mBufferLen;
	// our reference gap hash map vector and associated iterator
	vector<unordered_map<unsigned int, unsigned short> > mRefGapVector;
	unordered_map<unsigned int, unsigned short>::iterator mRefGapIter;
};
