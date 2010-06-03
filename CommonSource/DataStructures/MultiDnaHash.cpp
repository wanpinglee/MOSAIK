// ***************************************************************************
// CMultiDnaHash - genome hash map used in the multi algorithm. (9 pos / hash)
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "MultiDnaHash.h"

// defines the code for an empty hash position
const unsigned int CMultiDnaHash::DNA_EMPTY_HASH_POSITION = -1;

CMultiDnaHash::CMultiDnaHash(const unsigned char bitCapacity, const unsigned char hashSize)
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
		mHashPositions = new unsigned int[mCapacity * DNA_HASH_NUM_STORED];

	} catch(bad_alloc) {
		cout << "ERROR: Unable to allocate enough memory for the DNA hash map." << endl;
		exit(1);
	}

	// set the default settings for each element
	mMemoryAllocated = true;
	Clear();
}

CMultiDnaHash::~CMultiDnaHash(void) {
	if(mMemoryAllocated) FreeMemory();
}

// redimension the hash table to the specified size
void CMultiDnaHash::FreeMemory(void) {
	mMemoryAllocated = false;
	delete [] mHashes;
	delete [] mHashPositions;
}

// adds a fragment to the hash table
void CMultiDnaHash::Add(const uint64_t& key, const unsigned int genomePosition) {

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

		// check if we have any empty hash positions
		unsigned int startPos = position * DNA_HASH_NUM_STORED;
		unsigned int endPos   = startPos + DNA_HASH_NUM_STORED;

		for(unsigned int hashPos = startPos; hashPos < endPos; hashPos++) {

			// if there is a position available, add the current position
			if(mHashPositions[hashPos] == DNA_EMPTY_HASH_POSITION) {
				mHashPositions[hashPos] = genomePosition;
				break;
			}
		}

	} else {

		// add the keys to the hash table
		mHashes[position]                              = key;
		mHashPositions[position * DNA_HASH_NUM_STORED] = genomePosition;

		// increase the counter
		mCount++;
	}
}

// resets the counter and hash positions values
void CMultiDnaHash::Clear(void) {

	// set all of the elements to their default values
	unsigned int numHashPositions = mCapacity * DNA_HASH_NUM_STORED;
	uninitialized_fill(mHashes, mHashes + mCapacity, DNA_HASH_EMPTY_KEY);
	uninitialized_fill(mHashPositions, mHashPositions + numHashPositions, DNA_EMPTY_HASH_POSITION);

	// reset the counter and collisions variables
	mCount = 0;
}

// retrieves the genome location of the fragment
void CMultiDnaHash::Get(const uint64_t& key, const unsigned int& queryPosition, CHashRegionTree& hrt, double& mhpOccupancy) {

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
		unsigned int startPos = position * DNA_HASH_NUM_STORED;
		unsigned int endPos   = startPos + DNA_HASH_NUM_STORED;

		for(unsigned int hashPos = startPos; hashPos < endPos; hashPos++) {

			// if there is a position available, add the current position
			if(mHashPositions[hashPos] == DNA_EMPTY_HASH_POSITION) break;

			HashRegion island;
			island.Begin         = mHashPositions[hashPos];
			island.End           = mHashPositions[hashPos] + mHashSize - 1;
			island.QueryBegin    = queryPosition;
			island.QueryEnd      = queryPosition + mHashSize - 1;
			hrt.Insert(island);
		}
	}
}

// runs when we need to resize the hash table
void CMultiDnaHash::Resize(void) {

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
		uint64_t* tHashes        = new uint64_t[mCapacity];
		unsigned int* tHashPositions = new unsigned int[mCapacity * DNA_HASH_NUM_STORED];

		// copy the hash keys and delete the old hash keys
		memcpy(tHashes, mHashes, SIZEOF_UINT64 * mCapacity);
		delete [] mHashes;

		// copy the hash values and delete the old hash values
		memcpy(tHashPositions, mHashPositions, SIZEOF_INT * mCapacity * DNA_HASH_NUM_STORED);
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

		unsigned int numHashPositions = mCapacity * DNA_HASH_NUM_STORED;
		mHashes        = new uint64_t[mCapacity];
		mHashPositions = new unsigned int[numHashPositions];

		// set the default values for the new table
		uninitialized_fill(mHashes, mHashes + mCapacity, DNA_HASH_EMPTY_KEY);
		uninitialized_fill(mHashPositions, mHashPositions + numHashPositions, DNA_EMPTY_HASH_POSITION);

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

				unsigned char numPos = 0;
				unsigned int startPos = position * DNA_HASH_NUM_STORED;
				unsigned int endPos   = startPos + DNA_HASH_NUM_STORED;

				for(unsigned int hashPos = startPos; hashPos < endPos; hashPos++, numPos++) {
					if(hashPos < numHashPositions) mHashPositions[hashPos] = tHashPositions[i * DNA_HASH_NUM_STORED + numPos];
					else {
						cout << "ERROR: An invalid hash position was found when resizing the hash table." << endl;
						exit(1);
					}
				}
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
void CMultiDnaHash::Dump() {

	cout << "Multi DNA hash table contents:" << endl;
	cout << "==============================" << endl;

	unsigned int numDisplayedKeys = 0, numPositions = 0;
	for(unsigned int i = 0; i < mCapacity; i++)
		if(mHashes[i] != DNA_HASH_EMPTY_KEY) {

			cout << "key: " << mHashes[i] << ", positions:";

			// check if we have any empty hash positions
			unsigned int startPos = i * DNA_HASH_NUM_STORED;
			unsigned int endPos   = startPos + DNA_HASH_NUM_STORED;

			for(unsigned int hashPos = startPos; hashPos < endPos; hashPos++) {
				if(mHashPositions[hashPos] == DNA_EMPTY_HASH_POSITION) break;
				cout << " " << mHashPositions[hashPos];
				numPositions++;
			}

			cout << endl;

			numDisplayedKeys++;
		}

		cout << endl;
		cout << "keys found in hash table: " << numDisplayedKeys << ", positions found: " << numPositions << ", hash table count: " << mCount << endl;
}

// dummy function
void CMultiDnaHash::RandomizeAndTrimHashPositions(unsigned short numHashPositions) {}
