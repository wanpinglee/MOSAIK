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

#include "UbiqDnaHash.h"

CUbiqDnaHash::CUbiqDnaHash(const unsigned char bitCapacity, const unsigned char hashSize)
: mHashPositions(NULL)
, mPositions(0)
{
	if(bitCapacity == 32) {
		mCapacity  = UINT_MAX;
		mMask      = UINT_MAX;
		mLoad      = 1.0f;
		mThreshold = UINT_MAX;
	} else {
		mCapacity = 1 << bitCapacity;
		mMask     = mCapacity - 1;
		mLoad     = 0.8f;
		mThreshold = (unsigned int)(mCapacity * mLoad);
	}

	mCount     = 0;
	mHashSize  = hashSize;

	// create our hash table
	try {

		mHashes        = new uint64_t[mCapacity];
		mHashPositions = new vector<unsigned int>[mCapacity];

	} catch(bad_alloc) {
		cout << "ERROR: Unable to allocate enough memory for the DNA hash map." << endl;
		exit(1);
	}

	// set the default settings for each element
	mMemoryAllocated = true;
	Clear();
}

CUbiqDnaHash::~CUbiqDnaHash(void) {
	if(mMemoryAllocated) FreeMemory();
}

// frees all memory used by the hash table
void CUbiqDnaHash::FreeMemory(void) {
	mMemoryAllocated = false;
	delete [] mHashes;
	delete [] mHashPositions;
}

// adds a fragment to the hash table
void CUbiqDnaHash::Add(const uint64_t& key, const unsigned int genomePosition) {

	// check to see if we need to resize the hash table
	if((mCount + 1) > mThreshold) Resize();

	// retrieve the array position for this hash
	unsigned int position = IndexFor(key);
	if(position >= mCapacity) position = 0;

	// set our found key variable to false
	bool alreadyExists = false;

	// find an unused element
	while(mHashes[position] != DNA_HASH_EMPTY_KEY) {

		// check to see if it already exists
		if(mHashes[position] == key) {
			alreadyExists = true;
			break;
		}

		// get the next position
		position++;

		// wrap around if needed
		if(position >= mCapacity) position = 0;
	}

	if(!alreadyExists) {

		// assign the key
		mHashes[position] = key;

		// increase the counter
		mCount++;
	}

	mPositions++;

	// add the genome position to our vector
	mHashPositions[position].push_back(genomePosition);
}

// resets the counter and hash positions values
void CUbiqDnaHash::Clear(void) {

	// set all of the elements to their default values
	uninitialized_fill(mHashes, mHashes + mCapacity, DNA_HASH_EMPTY_KEY);
	for(unsigned int i = 0; i < mCapacity; i++) mHashPositions[i].clear();

	// reset the counter and collisions variables
	mCount = 0;
}

// retrieves the genome location of the fragment
void CUbiqDnaHash::Get(const uint64_t& key, const unsigned int& queryPosition, CHashRegionTree& hrt, double& mhpOccupancy) {

	// use a fixed mhp occupancy
	mhpOccupancy = 1.0;

	// retrieve the array position for this hash
	unsigned int position = IndexFor(key);
	if(position >= mCapacity) position = 0;

	// set our found key variable to false
	bool foundKey = false;

	// find an unused element
	while(mHashes[position] != DNA_HASH_EMPTY_KEY) {

		// check to see if it already exists
		if(mHashes[position] == key) {
			foundKey = true;
			break;
		}

		// get the next position
		position++;

		// wrap around if needed
		if(position >= mCapacity) position = 0;
	}

	// create new hash regions and add them to the tree
	if(foundKey) {

		// check if we have any empty hash positions		
		for(unsigned int i = 0; i < (unsigned int)mHashPositions[position].size(); i++) {

			// create a new island and add it to the island list
			HashRegion island;
			island.Begin         = mHashPositions[position][i];
			island.End           = mHashPositions[position][i] + mHashSize - 1;
			island.QueryBegin    = queryPosition;
			island.QueryEnd      = queryPosition + mHashSize - 1;
			hrt.Insert(island);
		}
	}
}

// runs when we need to resize the hash table
void CUbiqDnaHash::Resize(void) {

	// check to see if we're already at maximum capacity
	if(mCapacity == UINT_MAX) {
		cout << "ERROR: Cannot resize hash table. Already at maximum capacity." << endl;
		exit(1);
	}

	try {

		//
		// copy the elements to a new temporary hash table
		//

		// create temporary arrays
		uint64_t* tHashes                  = new uint64_t[mCapacity];
		vector<unsigned int>* tHashPositions = new vector<unsigned int>[mCapacity];

		// copy the hash keys and delete the old hash keys
		// N.B. copy integrity checked
		memcpy(tHashes, mHashes, SIZEOF_UINT64 * mCapacity);
		delete [] mHashes;

		// copy the hash positions and delete the old hash values
		for(unsigned int i = 0; i < mCapacity; i++) {

			// reserve the right amount of space
			tHashPositions[i].reserve(mHashPositions[i].size());

			// copy the vector elements
			copy(mHashPositions[i].begin(), mHashPositions[i].end(), back_inserter(tHashPositions[i]));
		}

		// delete the old hash positions
		delete [] mHashPositions;

		//
		// calculate the new hash table size
		//

		// save the old capacity
		unsigned int oldCapacity = mCapacity;

		// increase the capacity by a factor of 2
		if(mCapacity < LargestResizeableSize) {
			mCapacity = mCapacity << 1;
			mMask     = mCapacity - 1;
		} else {
			mCapacity = UINT_MAX;
			mMask     = UINT_MAX;
			mLoad     = 1.0;
		}

		// increase the threshold
		mThreshold = (unsigned int)(mCapacity * mLoad);

		//
		// populate the new hash table
		//

		mHashes        = new uint64_t[mCapacity];
		mHashPositions = new vector<unsigned int>[mCapacity];

		// set the default values for the new table
		// N.B. erase integrity checked
		uninitialized_fill(mHashes, mHashes + mCapacity, DNA_HASH_EMPTY_KEY);

		// copy the old values
		for(unsigned int i = 0; i < oldCapacity; i++) {

			// if it was an active element, add it to the new hash
			if(tHashes[i] != DNA_HASH_EMPTY_KEY) {

				// retrieve the array position for this hash
				unsigned int position = IndexFor(tHashes[i]);
				if(position >= mCapacity) position = 0;

				// find an unused element
				while(mHashes[position] != DNA_HASH_EMPTY_KEY) {

					// get the next position
					position++;

					// wrap around if needed
					if(position >= mCapacity) position = 0;
				}

				// copy the information from the old table to the new table
				mHashes[position] = tHashes[i];

				mHashPositions[position].reserve(tHashPositions[i].size());
				copy(tHashPositions[i].begin(), tHashPositions[i].end(), back_inserter(mHashPositions[position]));
			}
		}

		// delete the old tables
		delete [] tHashes;
		delete [] tHashPositions;

	} catch(bad_alloc &ba) {

		cout << "ERROR: Could not allocate enough memory to resize the hash table: " << ba.what() << endl;
		exit(1);
	}
}

// dumps the contents of the hash table to standard output
void CUbiqDnaHash::Dump() {

	cout << "Ubiq DNA hash table contents:" << endl;
	cout << "=============================" << endl;

	unsigned int numDisplayedKeys = 0, numPositions = 0;
	for(unsigned int i = 0; i < mCapacity; i++) {
		if(mHashes[i] != DNA_HASH_EMPTY_KEY) {

			cout << "key: " << mHashes[i] << ", positions:";

			// check if we have any empty hash positions
			for(unsigned int hashPos = 0; hashPos < (unsigned int)mHashPositions[i].size(); hashPos++, numPositions++)
				cout << " " << mHashPositions[i][hashPos];

			cout << endl;

			numDisplayedKeys++;
		}
	}

	cout << endl;
	cout << "keys found in hash table: " << numDisplayedKeys << ", positions found: " << numPositions << ", hash table count: " << mCount << endl;
}

// randomize and trim hash positions
void CUbiqDnaHash::RandomizeAndTrimHashPositions(unsigned short numHashPositions) {

	vector<unsigned int>* pHashPositions = NULL;

	// calculate the threshold if a zero parameter is given
	if(numHashPositions == 0) {

		cout << endl << "- calculating genome position threshold... ";
		cout.flush();

		double sum = 0.0;
		unsigned int numHashes = mCount;

		// calculate the sum
		for(unsigned int i = 0; i < mCapacity; i++)
			if(mHashes[i] != DNA_HASH_EMPTY_KEY)
				sum += mHashPositions[i].size();

		// calculate the mean number of hash positions
		double mean = sum / (double)numHashes;

		// calculate the standard deviation
		double diffSumSquare = 0.0;

		for(unsigned int i = 0; i < mCapacity; i++) {
			if(mHashes[i] != DNA_HASH_EMPTY_KEY) {
				double diffSum = mHashPositions[i].size() - mean;
				diffSumSquare += diffSum * diffSum;
			}
		}

		double variance = diffSumSquare / (numHashes - 1.0);
		double stddev   = sqrt(variance);

		numHashPositions = (unsigned short)(mean + 4.0 * stddev);

		cout << "finished." << endl;
		cout << "- setting the max number of hash positions per hash (" << numHashPositions << ")" << endl;
	}

	// randomize and trim
	for(unsigned int i = 0; i < mCapacity; i++) {
		if(mHashes[i] != DNA_HASH_EMPTY_KEY) {
			if(mHashPositions[i].size() > numHashPositions) {
				pHashPositions = &mHashPositions[i];
				random_shuffle(pHashPositions->begin(), pHashPositions->end());
				pHashPositions->erase(pHashPositions->begin() + numHashPositions, pHashPositions->end());
			}
		}
	}
}
