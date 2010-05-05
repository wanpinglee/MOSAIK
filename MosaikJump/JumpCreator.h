// ***************************************************************************
// CJumpCreator - creates a jump database for use with MosaikAligner.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iostream>
#include <memory>
#include "ConsoleUtilities.h"
#include "FileUtilities.h"
#include "MemoryUtilities.h"
#include "ProgressBar.h"
#include "ProgressCounter.h"
#include "ReferenceSequenceReader.h"

using namespace std;

#define KEY_LENGTH 5

class CJumpCreator {
public:
	// constructor
	CJumpCreator(const unsigned char hashSize, const string& filenameStub, const unsigned char sortingMemoryGB, const bool keepKeysInMemory, const unsigned int hashPositionThreshold);
	// destructor
	~CJumpCreator(void);
	// builds the jump database
	void BuildJumpDatabase(void);
	// enables hash position logging
	//void EnableHashPositionsLogging(const string& filename);
	// hashes the reference and stores the results in sorted temporary files
	void HashReference(const string& referenceFilename);
	// saves the metadata to the jump database
	void WriteMetadata(void);
private:
	struct HashPosition {
		uint64_t Hash;
		unsigned int Position;
		unsigned char Owner;

		// deserialize this object from the supplied file stream
		bool Deserialize(FILE* temp) {
			fread((char*)&Hash,     SIZEOF_UINT64, 1, temp);
			fread((char*)&Position, SIZEOF_INT,       1, temp);

			if(feof(temp)) return false;
			return true;
		}

		// serialize this object to the supplied file stream
		void Serialize(FILE* temp) {
			fwrite((char*)&Hash,     SIZEOF_UINT64, 1, temp);
			fwrite((char*)&Position, SIZEOF_INT,       1, temp);
		}
	};
	// define a comparison function for sorting our hash positions (ascending)
	struct SortHashPositionAsc {
		bool operator()(const HashPosition& hp1, const HashPosition& hp2) {
			if(hp1.Hash == hp2.Hash) return hp1.Position < hp2.Position;
			return hp1.Hash < hp2.Hash;
		}
	};
	// define a comparison function for sorting our hash positions (descending)
	struct SortHashPositionDesc {
		bool operator()(const HashPosition& hp1, const HashPosition& hp2) {
			if(hp1.Hash == hp2.Hash) return hp2.Position < hp1.Position;
			return hp2.Hash < hp1.Hash;
		}
	};
	// creates the hash for a supplied fragment
	static void CreateHash(const char* fragment, const unsigned char fragmentLen, uint64_t& key);
	// serializes the sorting vector to temporary files
	void SerializeSortingVector(vector<HashPosition>& hashPositions);
	// stores the supplied hash positions in the jump database
	void StoreHash(vector<HashPosition>& hashPositions);
	// stores the all of the serialized filenames used
	vector<string> mSerializedPositionsFilenames;
	// our hash size
	unsigned char mHashSize;
	// the number of GB RAM allocated for sorting
	unsigned char mSortingMemoryGB;
	// our jump database file handles
	FILE* mKeys;
	FILE* mPositions;
	//gzFile mHashPositionLog;
	// our output buffer
	unsigned char* mBuffer;
	// the output buffer size
	unsigned int mBufferLen;
	// the total number of hash positions sorted
	unsigned int mNumHashPositions;
	// sets the limit for how many hash positions should be retrieved
	unsigned int mMaxHashPositions;
	// toggles if keys should be kept in memory until processing is finished
	bool mKeepKeysInMemory;
	// toggles whether or not we return all genome positions or just a subset
	bool mLimitPositions;
	// toggles whether or not hash positions should be logged
	bool mLogHashPositions;
	// our key buffer when keys are kept in memory
	uint64_t* mKeyBuffer;
	// the key buffer size
	uint64_t mKeyBufferLen;
};
