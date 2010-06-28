// ***************************************************************************
// CSingleEndSort - sorts single-end reads by reference sequence position.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iostream>
#include <set>
#include <string>
#include <vector>
#include "sqlite3.h"
#include "Alignment.h"
#include "AlignmentReader.h"
#include "AlignmentWriter.h"
#include "ReadGroup.h"
#include "Mosaik.h"
#include "ConsoleUtilities.h"
#include "ProgressBar.h"
#include "UnorderedMap.h"

using namespace std;

#define READ_NAME_BUFFER_SIZE 16
#define GAP '-'

// define our serialization status flags
#define SE_IS_LONG_READ      1
#define SE_IS_REVERSE_STRAND 2

#define SQL_BUFFER_SIZE 10240

class CSingleEndSort {
public:
	// constructor
	CSingleEndSort(const unsigned int numCachedAlignments);
	// destructor
	~CSingleEndSort(void);
	// enables consed renaming
	void EnableConsedRenaming(void);
	// enables duplicate read filtering
	void EnableDuplicateFiltering(const string& duplicateDirectory);
	// processes multiply aligned reads
	void EnableNonUniqueMode(void);
	// sorts the input alignments and saves them to the output file
	void SaveAlignmentsOrderedByPosition(const string& inputFilename, const string& outputFilename);
private:
	// retrieves an alignment from the specified temporary file and adds it to the specified list
	void AddAlignment(FILE* tempFile, const unsigned int owner, list<Alignment>& alignments);
	// corrects the homopolymer gap order for reverse alignments
	//void CorrectHomopolymerGapOrder(Alignment& al);
	// retrieves an alignment from the specified temporary file
	bool GetAlignment(FILE* tempFile, const unsigned int owner, Alignment& al);
	// records the observed gaps in the specified reference 
	void RecordReferenceGaps(Alignment& al);
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
	// toggles if we want to use non-unique reads
	bool mSortNonUniqueMates;
	// our duplicate removal variables
	bool mRemoveDuplicates;
	string mDuplicateDirectory;
	// toggles if we want to append the alignment count to the read name
	bool mRenameReads;
};
