// ***************************************************************************
// CUbiqDnaHash - genome hash map used in the all algorithm. 
//                (unlimited genome positions / hash)
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <climits>
#include <cmath>
#include "AbstractDnaHash.h"
#include "MemoryUtilities.h"
#include "ProgressBar.h"

using namespace std;

class CUbiqDnaHash : public CAbstractDnaHash {
public:
	// constructor
	CUbiqDnaHash(const unsigned char bitCapacity, const unsigned char hashSize);
	// deconstructor
	~CUbiqDnaHash(void);
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
	// keeps track of hash positions
	vector<unsigned int>* mHashPositions;
	// stores the number of positions in the hash table
	unsigned int mPositions;
};
