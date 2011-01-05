// ***************************************************************************
// CJumpDnaHash - a disk/memory agnostic genome hash map used in the all 
//                algorithm. Much more memory efficient than a standard hash
//                map when storing mammalian genomes. (unlimited pos / hash)
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "AbstractDnaHash.h"
#include "FileUtilities.h"
#include "LargeFileSupport.h"
#include "MemoryUtilities.h"
#include "MruCache.h"

using namespace std;

#define KEY_LENGTH 5

class CJumpDnaHash : public CAbstractDnaHash {
public:
	// constructor
	CJumpDnaHash(const unsigned char hashSize, const string& filenameStub, const unsigned short numPositions, const bool keepKeysInMemory, const bool keepPositionsInMemory, const unsigned int numCachedElements, const unsigned int begin, const unsigned int end, const unsigned int offset, const unsigned int expectedMemory, const bool useLowMemory, const bool bubbleSpecialHashes, const uint64_t specialBegin, const unsigned int specialPercent);
	// destructor
	~CJumpDnaHash(void);
	// dummy function
	void Add(const uint64_t& key, const unsigned int genomePosition);
	// dummy function
	void Clear(void);
	// retrieves the genome location of the fragment
	void Get(const uint64_t& key, const unsigned int& queryPosition, CHashRegionTree& hrt, double& mhpOccupancy);
	// returns the numbers of jump database cache hits and misses
	void GetCacheStatistics(uint64_t& cacheHits, uint64_t& cacheMisses);
	// dumps the contents of the hash table to standard output
	void Dump();
	// close the jump database
	void FreeMemory(void);
	// dummy function
	void RandomizeAndTrimHashPositions(unsigned short numHashPositions);
	// load hash keys and positions from the file to memory
	void LoadKeysNPositions(void);
	void GetHashStatistics(const vector<pair<unsigned int, unsigned int> >& referenceSequences, vector<unsigned int>& nHashs, vector<unsigned int>& expectedMemories);

private:
	// loads the keys database into memory
	void LoadKeys(void);
	// loads the positions database into memory
	void LoadPositions(void);
	// dummy function
	void Resize(void);
	// specifies how many hash positions should be retrieved
	unsigned short mNumPositions;
	// toggles whether or not we return all hash positions or just a subset
	bool mLimitPositions;
	// toggles if the keys should be stored in memory
	bool mKeepKeysInMemory;
	// toggles if the positions should be stored in memory
	bool mKeepPositionsInMemory;
	// toggles if the hash table cache should be used
	bool mUseCache;
	// our jump database file handles
	FILE* mKeys;
	FILE* mMeta;
	FILE* mPositions;
	// our input/output buffer
	unsigned char* mBuffer;
	// the current buffer size
	unsigned int mBufferLen;
	// sets the limit for how many hash positions should be retrieved
	unsigned int mMaxHashPositions;
	// our input/output key buffer
	char* mKeyBuffer;
	uint64_t mKeyBufferLen;
	uintptr_t mKeyBufferPtr;
	// our input/output position buffer
	char* mPositionBuffer;
	uint64_t mPositionBufferLen;
	uintptr_t mPositionBufferPtr;
	// caches the most recently used hashes
	CMruCache<uint64_t, vector<unsigned int> > mMruCache;
	// load a block of hash positions
	inline void LoadBlockPositions( char* blockPosition, uint64_t& bytesLeft, const unsigned int& fillBufferSize );
	// Store hash positions
	inline void StorePositions ( off_type& curFilePosition, off_type& left, vector<unsigned int>& positions, const off_type keyOffest);
	// determine the chromosome which positions locating in
	void SetPositionDistribution(const vector<pair<unsigned int, unsigned int> >& referenceSequences, vector<unsigned int>& nHashs, vector<unsigned int>& expectedMemories, const vector<unsigned int>& positions);
	// the begining of current chromosome
	unsigned int _begin;
	// the end of current chromosome
	unsigned int _end;
	unsigned int _offset;
	unsigned int _expectedMemory;
	bool hasKeysNPositions;
	bool _useLowMemory;
	bool _bubbleSpecialHashes;
	uint64_t _specialBegin;
	unsigned int  _nSpecialHash;
};
