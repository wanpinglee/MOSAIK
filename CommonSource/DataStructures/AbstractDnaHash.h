// ***************************************************************************
// CAbstractDnaHash - superclass to the genome hash maps.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iostream>
#include "Mosaik.h"
#include "PosixThreads.h"
#include "HashRegionTree.h"

using namespace std;
using namespace AVLTree;

class CAbstractDnaHash {
public:
	CAbstractDnaHash(void);
	virtual ~CAbstractDnaHash(void) = 0;
	// adds a fragment to the hash table
	virtual void Add(const uint64_t& key, const unsigned int genomePosition) = 0;
	// resets the counter and hash positions values
	virtual void Clear(void) = 0;
	// retrieves the genome location of the fragment
	virtual void Get(const uint64_t& key, const unsigned int& queryPosition, CHashRegionTree& hrt, double& mhpOccupancy) = 0;
	// dumps the contents of the hash table to standard output
	virtual void Dump(void) = 0;
	// redimension the hash table to the specified size
	virtual void FreeMemory(void) = 0;
	// randomize and trim hash positions
	virtual void RandomizeAndTrimHashPositions(unsigned short numHashPositions) = 0;
	// register our thread mutexes
	static pthread_mutex_t mJumpCacheMutex;
	static pthread_mutex_t mJumpKeyMutex;
	static pthread_mutex_t mJumpPositionMutex;
	
protected:
	// translates the supplied hash to a position in the hash table
	inline unsigned int IndexFor(uint64_t index) const;
	// runs when we need to resize the hash table
	virtual void Resize(void) = 0;
	// stores the hashes
	uint64_t* mHashes;
	// registers how many elements our hash can handle
	unsigned int mCapacity;
	// specifies the mask used when calculating index
	unsigned int mMask;
	// stores the maximum hash table load
	float mLoad;
	// keeps a threshold for when we should increase the size of the hash table
	unsigned int mThreshold;
	// registers how many elements are actually present in our hash table
	unsigned int mCount;
	// stores the status of our allocated memory
	bool mMemoryAllocated;
	// defines the code for the empty dna hash code
	const static uint64_t DNA_HASH_EMPTY_KEY;
	// specifies the largest resizeable hash table size
	const static unsigned int LargestResizeableSize;
	// stores the current hash size
	unsigned char mHashSize;
};

// translates the supplied hash to a position in the hash table
inline unsigned int CAbstractDnaHash::IndexFor(uint64_t index) const {
	index = (~index) + (index << 21);
	index = index ^ (index >> 24);
	index = (index + (index << 3)) + (index << 8);
	index = index ^ (index >> 14);
	index = (index + (index << 2)) + (index << 4);
	index = index ^ (index >> 28);
	index = index + (index << 31);
	return index & mMask;
}

