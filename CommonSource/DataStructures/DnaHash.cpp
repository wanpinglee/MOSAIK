// ***************************************************************************
// CDnaHash - genome hash map used in the fast algorithm. (1 position / hash)
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "DnaHash.h"

// defines the code for a non-unique key
const unsigned int CDnaHash::DNA_HASH_NON_UNIQUE_KEY = -1;

CDnaHash::CDnaHash(const unsigned char bitCapacity, const unsigned char hashSize)
: mHashPositions(NULL)
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
		mHashPositions = new unsigned int[mCapacity];

	} catch(bad_alloc) {
		cout << "ERROR: Unable to allocate enough memory for the DNA hash map." << endl;
		exit(1);
	}

	// set the default settings for each element
	mMemoryAllocated = true;
	Clear();
}

CDnaHash::~CDnaHash(void) {
	if(mMemoryAllocated) FreeMemory();
}

// redimension the hash table to the specified size
void CDnaHash::FreeMemory(void) {
	mMemoryAllocated = false;
	delete [] mHashes;
	delete [] mHashPositions;
}

// adds a fragment to the hash table
void CDnaHash::Add(const uint64_t& key, const unsigned int genomePosition) {

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

	if(alreadyExists) {

		// set the hash position element to NON_UNIQUE_KEY
		mHashPositions[position] = DNA_HASH_NON_UNIQUE_KEY;

	} else {

		// add the keys to the hash table
		mHashes[position]        = key;
		mHashPositions[position] = genomePosition;

		// increase the counter
		mCount++;
	}
}

// increments a counter every time the hash is seen
void CDnaHash::AddCount(const uint64_t& key) {

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

	if(alreadyExists) {

		// increment the counter
		mHashPositions[position]++;

	} else {

		// add the keys to the hash table
		mHashes[position]        = key;
		mHashPositions[position] = 1;

		// increase the counter
		mCount++;
	}
}

// resets the counter and hash positions values
void CDnaHash::Clear(void) {

	// set all of the elements to their default values
	uninitialized_fill(mHashes, mHashes + mCapacity, DNA_HASH_EMPTY_KEY);
	uninitialized_fill(mHashPositions, mHashPositions + mCapacity, 0);

	// reset the counter and collisions variables
	mCount = 0;
}

// retrieves the genome location of the fragment
void CDnaHash::Get(const uint64_t& key, const unsigned int& queryPosition, CHashRegionTree& hrt, double& mhpOccupancy) {

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

	// create a new hash region and add it to the tree
	if(foundKey && (mHashPositions[position] != DNA_HASH_NON_UNIQUE_KEY)) {
		HashRegion island;
		island.Begin         = mHashPositions[position];
		island.End           = mHashPositions[position] + mHashSize - 1;
		island.QueryBegin    = queryPosition;
		island.QueryEnd      = queryPosition + mHashSize - 1;
		hrt.Insert(island);
	}
}

// runs when we need to resize the hash table
void CDnaHash::Resize(void) {

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
		uint64_t* tHashes = new uint64_t[mCapacity];
		unsigned int* tHashPositions  = new unsigned int[mCapacity];

		// copy the hash keys and delete the old hash keys
		// N.B. copy integrity checked
		memcpy(tHashes, mHashes, SIZEOF_UINT64 * mCapacity);
		delete [] mHashes;

		// copy the hash values and delete the old hash values
		memcpy(tHashPositions, mHashPositions, SIZEOF_INT * mCapacity);
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
		mHashPositions = new unsigned int[mCapacity];

		// set the default values for the new table
		uninitialized_fill(mHashes, mHashes + mCapacity, DNA_HASH_EMPTY_KEY);
		uninitialized_fill(mHashPositions, mHashPositions + mCapacity, 0);

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
				mHashes[position]        = tHashes[i];
				mHashPositions[position] = tHashPositions[i];
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
void CDnaHash::Dump() {

	cout << "DNA hash table contents:" << endl;
	cout << "========================" << endl;

	unsigned int numDisplayedKeys = 0, numPositions = 0;
	for(unsigned int i = 0; i < mCapacity; i++)
		if(mHashes[i] != DNA_HASH_EMPTY_KEY) {
			cout << "key: " << mHashes[i] << ", position: ";

			if(mHashPositions[i] == DNA_HASH_NON_UNIQUE_KEY) {
				cout << "non-unique" << endl;
			} else {
				cout << mHashPositions[i] << endl;
				numPositions++;
			}

			numDisplayedKeys++;
		}

		cout << endl;
		cout << "keys found in hash table: " << numDisplayedKeys << ", positions found: " << numPositions << ", hash table count: " << mCount << endl;
}

// returns statistics about the hash table
void CDnaHash::GetStatistics(unsigned int& numUsedHashes, unsigned int& numUniqueHashes, unsigned int& numNonUniqueHashes, unsigned int& numUsedHashesCount, unsigned int& numUniqueHashesCount, unsigned int& numNonUniqueHashesCount, double& mean, double& stddev) {

	// initialization
	numUsedHashes           = 0;
	numUniqueHashes         = 0;
	numNonUniqueHashes      = 0;
	numUsedHashesCount      = 0;
	numUniqueHashesCount    = 0;
	numNonUniqueHashesCount = 0;

	for(unsigned int i = 0; i < mCapacity; i++)
		if(mHashes[i] != DNA_HASH_EMPTY_KEY) {
			numUsedHashes++;
			numUsedHashesCount += mHashPositions[i];

			if(mHashPositions[i] != 1) {
				numNonUniqueHashes++;
				numNonUniqueHashesCount += mHashPositions[i];
			} else {
				numUniqueHashes++;
				numUniqueHashesCount++;
			}
		}

		// calculate the mean number of hash positions
		mean = (double)numUsedHashesCount / (double)numUsedHashes;

		// calculate the standard deviation
		double diffSumSquare = 0.0;

		for(unsigned int i = 0; i < mCapacity; i++)
			if(mHashes[i] != DNA_HASH_EMPTY_KEY) {
				double diffSum = mHashPositions[i] - mean;
				diffSumSquare += diffSum * diffSum;
			}

			double variance = diffSumSquare / (numUsedHashes - 1.0);
			stddev          = sqrt(variance);
}

// dummy function
void CDnaHash::RandomizeAndTrimHashPositions(unsigned short numHashPositions) {}
