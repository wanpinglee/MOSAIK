// ***************************************************************************
// CDnaHash - genome hash map used in the fast algorithm. (1 position / hash)
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
#include <cmath>
#include <cstdlib>
#include "AbstractDnaHash.h"
#include "MemoryUtilities.h"

using namespace std;

class CDnaHash : public CAbstractDnaHash {
public:
	CDnaHash(const unsigned char bitCapacity, const unsigned char hashSize);
	~CDnaHash(void);
	// adds a fragment to the hash table
	void Add(const uint64_t& key, const unsigned int genomePosition);
	// increments a counter every time the hash is seen
	void AddCount(const uint64_t& key);
	// resets the counter and hash positions values
	void Clear(void);
	// retrieves the genome location of the fragment
	void Get(const uint64_t& key, const unsigned int& queryPosition, CHashRegionTree& hrt, double& mhpOccupancy);
	// returns statistics about the hash table
	void GetStatistics(unsigned int& numUsedHashes, unsigned int& numUniqueHashes, unsigned int& numNonUniqueHashes, unsigned int& numUsedHashesCount, unsigned int& numUniqueHashesCount, unsigned int& numNonUniqueHashesCount, double& mean, double& stddev);
	// dumps the contents of the hash table to standard output
	void Dump();
	// frees all memory used by the hash table
	void FreeMemory(void);
	// randomize and trim hash positions
	void RandomizeAndTrimHashPositions(unsigned short numHashPositions);

private:
	// runs when we need to resize the hash table
	void Resize(void);
	// stores our hash positions
	unsigned int* mHashPositions;
	// defines the code for a non-unique key
	const static unsigned int DNA_HASH_NON_UNIQUE_KEY;
};
