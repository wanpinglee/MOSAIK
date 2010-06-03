// ***************************************************************************
// CMultiDnaHash - genome hash map used in the multi algorithm. (9 pos / hash)
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
#include <climits>
#include "AbstractDnaHash.h"
#include "MemoryUtilities.h"

using namespace std;

// indicate the number of positions stored per hash position
#define DNA_HASH_NUM_STORED		9

class CMultiDnaHash : public CAbstractDnaHash {
public:
	CMultiDnaHash(const unsigned char bitCapacity, const unsigned char hashSize);
	~CMultiDnaHash(void);
	// adds a fragment to the hash table
	void Add(const uint64_t& key, const unsigned int genomePosition);
	// resets the counter and hash positions values
	void Clear(void);
	// retrieves the genome location of the fragment
	void Get(const uint64_t& key, const unsigned int& queryPosition, CHashRegionTree& hrt, double& mhpOccupancy);
	// dumps the contents of the hash table to standard output
	void Dump();
	// frees all memory used by the hash table
	void FreeMemory(void);
	// randomize and trim hash positions
	void RandomizeAndTrimHashPositions(unsigned short numHashPositions);

private:	
	// runs when we need to resize the hash table
	void Resize(void);
	// stores track of hash positions
	unsigned int* mHashPositions;
	// defines the code for an empty hash position
	const static unsigned int DNA_EMPTY_HASH_POSITION;
};
