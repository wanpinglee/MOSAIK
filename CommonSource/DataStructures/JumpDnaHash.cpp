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

#include "JumpDnaHash.h"

// constructor
CJumpDnaHash::CJumpDnaHash(const unsigned char hashSize, const string& filenameStub, const unsigned short numPositions, const bool keepKeysInMemory, const bool keepPositionsInMemory, const unsigned int numCachedElements, const unsigned int begin, const unsigned int end, const unsigned int offset, const unsigned int expectedMemory, const bool useLowMemory, const bool bubbleSpecialHashes, const uint64_t specialBegin, const unsigned int nSpecialHash)
: mNumPositions(numPositions)
, mLimitPositions(false)
, mKeepKeysInMemory(keepKeysInMemory)
, mKeepPositionsInMemory(keepPositionsInMemory)
, mUseCache(false)
, mKeys(NULL)
, mPositions(NULL)
, mBuffer(NULL)
, mBufferLen(4096)
, mKeyBuffer(NULL)
, mKeyBufferLen(0)
, mKeyBufferPtr(0)
, mPositionBuffer(NULL)
, mPositionBufferLen(0)
, mPositionBufferPtr(0)
, mMruCache(numCachedElements)
, _begin(begin)
, _end(end)
, _offset(offset)
, _expectedMemory(expectedMemory)
, _useLowMemory(useLowMemory)
, _bubbleSpecialHashes(bubbleSpecialHashes)
, _specialBegin(specialBegin)
, _nSpecialHash(nSpecialHash)
{
	mHashSize = hashSize;

	// generate our filenames
	string keysFilename      = filenameStub + "_keys.jmp";
	string metaFilename      = filenameStub + "_meta.jmp";
	string positionsFilename = filenameStub + "_positions.jmp";

	// get the file sizes
	CFileUtilities::GetFileSize(keysFilename, mKeyBufferLen);
	CFileUtilities::GetFileSize(positionsFilename, mPositionBufferLen);

	// open our files
	fopen_s(&mKeys, keysFilename.c_str(), "rb");

	if(!mKeys) {
		cout << "ERROR: Unable to open the keys file (" << keysFilename << ") for reading." << endl;
		exit(1);
	}

	fopen_s(&mMeta, metaFilename.c_str(), "rb");

	if(!mMeta) {
		cout << "ERROR: Unable to open the metadata file (" << metaFilename << ") for reading." << endl;
		exit(1);
	}

	fopen_s(&mPositions, positionsFilename.c_str(), "rb");

	if(!mPositions) {
		cout << "ERROR: Unable to open the positions file (" << positionsFilename << ") for reading." << endl;
		exit(1);
	}

	// initialize the file buffer
	try {
		mBuffer = new unsigned char[mBufferLen];
	} catch(bad_alloc) {
		cout << "ERROR: Unable to allocate enough memory for the jump database buffer." << endl;
		exit(1);
	}

	// check the hash size
	unsigned char jumpHashSize = fgetc(mMeta);
	if(jumpHashSize != hashSize) {
		cout << "ERROR: The supplied hash size (" << (short)hashSize << ") is different from the hash size of the jump database (" << (short)jumpHashSize << "). Please create another jump database or select a different hash size." << endl;
		exit(1);
	}

	// close the metadata file
	fclose(mMeta);

	// place the keys and positions in memory
	//if(keepKeysInMemory)      LoadKeys();
	//if(keepPositionsInMemory) LoadPositions();
	hasKeysNPositions = false;


	// activate the MRU cache
	if(numCachedElements > 0) mUseCache = true;
	if(keepKeysInMemory && keepPositionsInMemory) mUseCache = false;

	// limit the number of hash positions
	if(numPositions > 0) RandomizeAndTrimHashPositions(numPositions);

	// set the default settings for each element
	mMemoryAllocated = true;
}

// deconstructor
CJumpDnaHash::~CJumpDnaHash(void) {
	if(mMemoryAllocated) FreeMemory();
}

// dummy function
void CJumpDnaHash::Add(const uint64_t& key, const unsigned int genomePosition) {
	cout << "ERROR: This function has not been implemented. Please use MosaikJump to create jump databases." << endl;
	exit(1);
}

// dummy function
void CJumpDnaHash::Clear(void) {}

// dummy function
void CJumpDnaHash::Dump() {
	cout << "ERROR: This function has not been implemented." << endl;
	exit(1);
}

// close the jump database
void CJumpDnaHash::FreeMemory(void) {
	if(mBuffer) {
		delete [] mBuffer;
		mBuffer = NULL;
	}

	if(mKeyBuffer) {
		delete [] mKeyBuffer;
		mKeyBuffer = NULL;
	}
	
	if(mPositionBuffer) {
		delete [] mPositionBuffer;
		mPositionBuffer = NULL;
	}

	hasKeysNPositions = false;
}

// load hash keys and positions form file to memory
void CJumpDnaHash::LoadKeysNPositions() {
	
	if(mKeepKeysInMemory) {
		cout << "- loading jump key database into memory... ";
		cout.flush();
		LoadKeys();
	}
	cout << "finished." << endl;
	
	if(mKeepPositionsInMemory) {
		cout << "- loading jump positions database into memory... ";
		cout.flush();
		LoadPositions();
	}
	cout << "finished." << endl;

	hasKeysNPositions = true;
}

// retrieves the genome location of the fragment
void CJumpDnaHash::Get(const uint64_t& key, const unsigned int& queryPosition, CHashRegionTree& hrt, double& mhpOccupancy) {

	if ( !hasKeysNPositions ) {
		cout << "Have not loaded hash keys and positions before using them" << endl;
		exit(1);
	}
	
	// find the correct position in the keys database
	const off_type offset = key * KEY_LENGTH;
	off_type position = 0;

	// initialize the mhp occupancy
	mhpOccupancy = 1.0;

	// ===================
	// check the MRU cache
	// ===================

	if(mUseCache) {
		vector<unsigned int> positionVector;

		pthread_mutex_lock(&mJumpCacheMutex);
		bool isCached = mMruCache.Get(key, positionVector);
		pthread_mutex_unlock(&mJumpCacheMutex);

		// TODO: handle the mhp occupancy. How do we get the mhp occupancy when using the cache?

		if(isCached) {
			for(unsigned int i = 0; i < positionVector.size(); i++) {
				HashRegion island;
				island.Begin      = positionVector[i];
				island.End        = positionVector[i] + mHashSize - 1;
				island.QueryBegin = queryPosition;
				island.QueryEnd   = queryPosition + mHashSize - 1;
				hrt.Insert(island);
			}
			return;
		}
	}

	// ==========================
	// retrieve the file position
	// ==========================

	if(mKeepKeysInMemory) {
		memcpy((char*)&position, (char*)(mKeyBufferPtr + offset), KEY_LENGTH);
	} else {
		pthread_mutex_lock(&mJumpKeyMutex);
		fseek64(mKeys, offset, SEEK_SET);
		fread((char*)&position, KEY_LENGTH, 1, mKeys);
		pthread_mutex_unlock(&mJumpKeyMutex);
	}

	// return if the key is undefined
	if(position == 0xffffffffffULL) return;

	if((uint64_t)position > mPositionBufferLen) {
		cout << "ERROR: A position (" << position << ") was specified that is larger than the jump positions database (" << mPositionBufferLen << ")." << endl;
		exit(1);
	}

	// ===========================
	// retrieve the hash positions
	// ===========================

	if(mKeepPositionsInMemory) {

		char* pPositions = (char*)(mPositionBufferPtr + position);

		unsigned int numPositions = 0;
		memcpy((char*)&numPositions, pPositions, SIZEOF_INT);
		unsigned int bufferOffset = SIZEOF_INT;

		// load positions
		//vector<unsigned int> hashPositions;
		//hashPositions.resize(numPositions);
		//unsigned int hashPosition = 0;
		//for(unsigned int i = 0; i < numPositions; ) {
		//	memcpy((char*)&hashPosition, pPositions + bufferOffset, SIZEOF_INT);
		//	bufferOffset += SIZEOF_INT;
		//	if ( hashPosition == 0xff ) {
				// skip next int which indicates how many continuous chromosomes not having hash hits
		//		bufferOffset += SIZEOF_INT;
		//	} 
		//	else {
		//		hashPositions[i] = hashPosition;
		//		i++;
		//	}
		//}

		// random_shuffle
		//random_shuffle(hashPositions.begin(), hashPositions.end());

		
		// set the mhp occupancy
		bool found = false;
		if(mLimitPositions && (numPositions > mMaxHashPositions)) {
			mhpOccupancy = (double)mMaxHashPositions / (double)numPositions;
			numPositions = mMaxHashPositions;
			found = true;
		}

		unsigned int hashPosition = 0;
		for(unsigned int i = 0; i < numPositions; i++) {
			memcpy((char*)&hashPosition, pPositions + bufferOffset, SIZEOF_INT);
			bufferOffset += SIZEOF_INT;

			HashRegion island;
			island.Begin         = hashPosition;
			island.End           = hashPosition  + mHashSize - 1;
			island.QueryBegin    = queryPosition;
			island.QueryEnd      = queryPosition + mHashSize - 1;
			hrt.Insert(island);
			//if ( hrt.Insert(island) && found )
			//	i--;

		}

	} else {

		pthread_mutex_lock(&mJumpPositionMutex);
		fseek64(mPositions, position, SEEK_SET);

		unsigned int numPositions = 0;
		fread((char*)&numPositions, SIZEOF_INT, 1, mPositions);

		unsigned int entrySize = numPositions * SIZEOF_INT;
		CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, entrySize);

		fread(mBuffer, entrySize, 1, mPositions);
		pthread_mutex_unlock(&mJumpPositionMutex);

		// set the mhp occupancy
		if(mLimitPositions && (numPositions > mMaxHashPositions)) {
			mhpOccupancy = (double)mMaxHashPositions / (double)numPositions;
			numPositions = mMaxHashPositions;
		}

		unsigned int bufferOffset = 0;
		unsigned int hashPosition = 0;

		vector<unsigned int> positionVector;
		if(mUseCache) positionVector.resize(numPositions);

		for(unsigned int i = 0; i < numPositions; i++) {
			memcpy((char*)&hashPosition, mBuffer + bufferOffset, SIZEOF_INT);
			bufferOffset += SIZEOF_INT;

			if(mUseCache) positionVector[i] = hashPosition;

			HashRegion island;
			island.Begin         = hashPosition;
			island.End           = hashPosition  + mHashSize - 1;
			island.QueryBegin    = queryPosition;
			island.QueryEnd      = queryPosition + mHashSize - 1;
			hrt.Insert(island);
		}

		if(mUseCache) {
			pthread_mutex_lock(&mJumpCacheMutex);
			mMruCache.Insert(key, positionVector);
			pthread_mutex_unlock(&mJumpCacheMutex);
		}
	}
}

// returns the numbers of jump database cache hits and misses
void CJumpDnaHash::GetCacheStatistics(uint64_t& cacheHits, uint64_t& cacheMisses) {
	mMruCache.GetStatistics(cacheHits, cacheMisses);
}

// get the distribution of # hashs aginst the chromosomes
void CJumpDnaHash::GetHashStatistics(
	const vector<pair<unsigned int, unsigned int> >& referenceSequences, 
	vector<unsigned int>& nHashs, 
	vector<unsigned int>& expectedMemories,
	const bool&           hasSpecial,
	const unsigned int&   specialBegin,
	const unsigned int&   specialMaxHashPositions) {
	
	LoadKeys();

	const unsigned int fillBufferSize  = 536870912; // 500 MB
	char*        blockPosition;
	uintptr_t    blockPositionPtr;
	unsigned int nBlock = 0; // indicates how many blocks have been handled
	blockPosition = new char[ (size_t) fillBufferSize ];

	if( !blockPosition ) {
		cout << "ERROR: Memory allocation for the temporary jump positions failed." << endl;
		exit(1);
	}

	uint64_t bytesLeft           = mPositionBufferLen;
	//uint64_t mPositionBufferLen1 = _end - _begin + 1;
	LoadBlockPositions( blockPosition, bytesLeft, fillBufferSize );

	off_type offset          = 0;                   // for mKeyBufferPtr
	//off_type curFilePosition = 0;                   // for mPositionBufferPtr
	//off_type left            = mPositionBufferLen1; // for mPositionBuffer full detection

	//cout << "- loading jump positions database into memory... ";
	//cout.flush();

	//off_type noHash = 0;
	while ( (uint64_t)offset < mKeyBufferLen ) {

		// get the position of the position file for the current key
		// pointer of key
		off_type filePosition = 0;
		memcpy((char*)&filePosition, (char*)(mKeyBufferPtr + offset), KEY_LENGTH);

		// no hash hits
		if ( filePosition ==  0xffffffffffULL ) {
			offset += KEY_LENGTH;
			//noHash++;
			continue;
		}

		
		// if the required filePosition isn't within the current block,
		// then load the next block of positions
		if ( filePosition >= (off_type)( nBlock + 1 ) * fillBufferSize ) {
			LoadBlockPositions( blockPosition, bytesLeft, fillBufferSize );
			nBlock++;
		}

		// load number of hash hits
		blockPositionPtr   = (uintptr_t)&blockPosition[0];
		off_type posOffset = filePosition - (off_type) nBlock * fillBufferSize;
		unsigned int numPositions;
		memcpy((char*)&numPositions, (char*)(blockPositionPtr + posOffset), SIZEOF_INT);
		

		vector <unsigned int> positions;
		positions.reserve(numPositions);
		// load all positions of hash hits
		for ( unsigned int i = 0; i < numPositions; i++ ) {
			
			filePosition += SIZEOF_INT;
			// if the required filePosition isn't within the current block,
			// then load the next block of positions
			if ( filePosition >= (off_type)( nBlock + 1 ) * fillBufferSize ) {
				LoadBlockPositions( blockPosition, bytesLeft, fillBufferSize );
				nBlock++;
			}

			posOffset = filePosition - (off_type) nBlock * fillBufferSize;
			unsigned int hashPosition;
			memcpy((char*)&hashPosition, (char*)(blockPositionPtr + posOffset), SIZEOF_INT);
			
			// the hash position is not within the current chromosome
			if ( hashPosition > _end )
				break;
			if ( hashPosition < _begin )
				continue;
			
			// the hash position is within the current chromosome, and keep it in the vector
			
			if ( hashPosition < _offset ) {
				cout << "ERROR: The hash position is smaller than offset." << endl;
				exit(1);
			}
			
			hashPosition -= _offset;
			
			positions.push_back(hashPosition);
		}

		
		// has hash hits
		if ( positions.size() != 0 ) {

			// handle sepcail references
			if ( hasSpecial ) {
				bool found = false;
				unsigned int nSpecial = 0;
				for ( vector<unsigned int>::reverse_iterator rite = positions.rbegin(); rite != positions.rend(); ++rite ) {
					
					if ( *rite < specialBegin )
						break;
					
					if ( nSpecial <= specialMaxHashPositions ) {
						found = true;
						// the last slot is for special reference
						nHashs[ nHashs.size() - 1 ]++;
					}

					nSpecial++;
				}
				if ( found )
					// the last slot is for special reference
					expectedMemories[ expectedMemories.size() - 1 ]++;

				// remove special positions
				if ( nSpecial > 0 ) {
					unsigned int eraseBegin = positions.size() - nSpecial;
					positions.erase( positions.begin() + eraseBegin, positions.end() );
				}

			}
			
			// handle regular references
			if ( positions.size() != 0 ) {
				random_shuffle( positions.begin(), positions.end(), randomGenerator );
				SetPositionDistribution(referenceSequences, nHashs, expectedMemories, positions);
			}
			
		}

		positions.clear();

		offset += KEY_LENGTH;
	}

	delete [] blockPosition;

	//cout << "finished." << endl;
	//
	//cout << endl << noHash << endl;

	fclose(mPositions);


}

// determine the chromosome which positions locating in
void CJumpDnaHash::SetPositionDistribution(
	const vector<pair<unsigned int, unsigned int> >& referenceSequences, 
	vector<unsigned int>&       nHashs, 
	vector<unsigned int>&       expectedMemories, 
	const vector<unsigned int>& positions
	//const bool&                 hasSpecial,
	//const unsigned int&         specialBegin,
	//const unsigned int&         specialMaxHashPositions
) {

	vector <bool> hasPositions;
	hasPositions.resize(nHashs.size(), false);

	unsigned int nPositions = 0;

	if ( mLimitPositions && ( positions.size() > mMaxHashPositions ) )
		nPositions = mMaxHashPositions;
	else
		nPositions = positions.size();
	
	for ( unsigned int i = 0; i < nPositions; i++ ) {
		unsigned int refNo = 0;
		// search the position belonging to which group
		while( positions[i] > referenceSequences[refNo].second ) refNo++;
		
		nHashs[refNo]++;

		if ( !hasPositions[refNo] ) {
			hasPositions[refNo] = true;
			expectedMemories[refNo]++;
		}
		
		//expectedMemories[refNo]++;

	}
}


// loads the keys database into memory
void CJumpDnaHash::LoadKeys(void) {

	// check if we can allocate enough memory
	// TODO: find a platform-independent fix for this check
	//if(mKeyBufferLen > SIZE_MAX) {
	//	cout << "ERROR: Cannot allocate enough memory to store the jump keys database (" << mKeyBufferLen 
	//		<< " bytes). The largest allocation size is " << SIZE_MAX 
	//		<< " bytes. Try using a smaller hash size or leaving the jump keys database on disk." << endl;
	//	exit(1);
	//}

	mKeyBuffer = new char[(size_t)mKeyBufferLen];
	mKeyBufferPtr = (uintptr_t)&mKeyBuffer[0];

	if(!mKeyBuffer) {
		cout << "ERROR: Memory allocation for the jump keys database failed." << endl;
		exit(1);
	}

	uint64_t bytesLeft = mKeyBufferLen;
	const unsigned int fillBufferSize = 2147483648ULL; // 2 GB

	//cout << "- loading jump keys database into memory... ";
	//cout.flush();

	char* pKeys = (char*)mKeyBuffer;
	while(bytesLeft > fillBufferSize) {
		fread(pKeys, fillBufferSize, 1, mKeys);
		pKeys     += fillBufferSize;
		bytesLeft -= fillBufferSize;
	}

	fread(pKeys, (size_t)bytesLeft, 1, mKeys);
	//cout << "finished." << endl;

	fclose(mKeys);
}

// load partial positions into memory for low-memory algorithm
inline void CJumpDnaHash::LoadBlockPositions( char* blockPosition, uint64_t& bytesLeft, const unsigned int& fillBufferSize ) {
	
	// clear the buffer
	memset(blockPosition, 0, (size_t)fillBufferSize);
	
	if ( bytesLeft < fillBufferSize ) {
		fread(blockPosition, (size_t)bytesLeft, 1, mPositions);
		bytesLeft -= bytesLeft;
	}
	else {
		fread(blockPosition, (size_t)fillBufferSize, 1, mPositions);
		bytesLeft -= fillBufferSize;
	}
}


// loads the positions database into memory
void CJumpDnaHash::LoadPositions(void) {

	// check if we can allocate enough memory
	// TODO: find a platform-independent fix for this check
	//if(mPositionBufferLen > SIZE_MAX) {
	//	cout << "ERROR: Cannot allocate enough memory to store the jump positions database (" << mPositionBufferLen 
	//		<< " bytes). The largest allocation size is " << SIZE_MAX 
	//		<< " bytes. Try using a smaller reference sequence or leaving the jump positions database on disk." << endl;
	//	exit(1);
	//}
	

	// prepare the buffer for loading positions from the file
	// for low-memory usage, each tile we load fillBufferSize positions
	//const unsigned int fillBufferSize  = 1073741824ULL; // 1 GB
	const unsigned int fillBufferSize  = 536870912; // 500 MB
	char*        blockPosition;
	uintptr_t    blockPositionPtr;
	unsigned int nBlock = 0; // indicates how many blocks have been handled
	blockPosition = new char[(size_t)fillBufferSize];

	if( !blockPosition ) {
		cout << "ERROR: Memory allocation for the temporary jump positions failed." << endl;
		exit(1);
	}


	// initialize positions memory
	uint64_t mPositionBufferLen1 = _useLowMemory ? (_expectedMemory * SIZEOF_INT) : mPositionBufferLen;

	mPositionBuffer    = new char[mPositionBufferLen1];
	mPositionBufferPtr = (uintptr_t)&mPositionBuffer[0];

	if ( !mPositionBuffer ) {
		cout << "ERROR: Memory allocation for the jump positions failed." << endl;
		exit(1);
	}

	// load the first block of positions from the file
	// bytesLeft indicates how many positions are left in the file
	uint64_t bytesLeft = mPositionBufferLen;
	LoadBlockPositions( blockPosition, bytesLeft, fillBufferSize );

	
	off_type offset          = 0;                   // for mKeyBufferPtr
	off_type curFilePosition = 0;                   // for mPositionBufferPtr
	off_type left            = mPositionBufferLen1; // for mPositionBuffer full detection

	while ( (uint64_t)offset < mKeyBufferLen ) {

		// get the position of the position file for the current key
		// pointer of key
		off_type filePosition = 0;
		memcpy((char*)&filePosition, (char*)(mKeyBufferPtr + offset), KEY_LENGTH);

		// no hash hits
		if ( filePosition ==  0xffffffffffULL ) {
			offset += KEY_LENGTH;
			continue;
		}

		
		// if the required filePosition isn't within the current block,
		// then load the next block of positions
		if ( filePosition >= (off_type)(nBlock+1)*fillBufferSize ) {
			LoadBlockPositions( blockPosition, bytesLeft, fillBufferSize );
			nBlock++;
		}

		// load number of hash hits
		blockPositionPtr   = (uintptr_t)&blockPosition[0];
		off_type posOffset = filePosition - (off_type)nBlock*fillBufferSize;
		unsigned int numPositions;
		memcpy((char*)&numPositions, (char*)(blockPositionPtr + posOffset), SIZEOF_INT);
		

		vector <unsigned int> positions;
		positions.reserve(numPositions);
		// load all positions of hash hits
		for ( unsigned int i = 0; i < numPositions; i++ ) {
			
			filePosition += SIZEOF_INT;
			if ( filePosition >= (off_type)(nBlock+1)*fillBufferSize ) {
				LoadBlockPositions( blockPosition, bytesLeft, fillBufferSize );
				nBlock++;
			}

			posOffset = filePosition - (off_type)nBlock*fillBufferSize;
			unsigned int hashPosition;
			memcpy((char*)&hashPosition, (char*)(blockPositionPtr + posOffset), SIZEOF_INT);
			
			// the hash position is not within the current chromosome
			if ( hashPosition > _end )
				break;
			if ( hashPosition < _begin )
				continue;
			
			// the hash position is within the current chromosome, and keep it in the vector
			
			if ( hashPosition < _offset ) {
				cout << "ERROR: The hash position is smaller than offset." << endl;
				exit(1);
			}
			
			hashPosition -= _offset;
			
			positions.push_back(hashPosition);
		}

		
		// no hash hit
		if ( positions.size() == 0 ) {
			// revise the key pointer
			memset((char*)(mKeyBufferPtr + offset), 0xff, KEY_LENGTH);
		}
		else {
			// NOTE: that low-memory won't go into here
			if ( _bubbleSpecialHashes ) {
				unsigned int totalPos     = positions.size();
				unsigned int totalSpecial = 0;
				for ( vector <unsigned int>::reverse_iterator rit = positions.rbegin(); rit != positions.rend(); ++rit ) {
					if ( *rit > _specialBegin )
						totalSpecial++;
				}
				
				//if ( totalSpecial > 0 ) {
				//	unsigned short i = 0;
				//	uint64_t keyt = offset / KEY_LENGTH;
				//	string keys;
				//	keys.resize(15);
				//	while ( i < 15 ) {
				//		uint64_t keyc = keyt & 0x0000000000000003ULL;
				//		if ( keyc == 0 ) keys[14-i] = 'A';
				//		else if ( keyc == 1 ) keys[14-i] = 'C';
				//		else if ( keyc == 2 ) keys[14-i] = 'G';
				//		else if ( keyc == 3 ) keys[14-i] = 'T';
				//		keyt = keyt >> 2;
				//		++i;
				//	}
				//	cerr << keys << "\t" << totalSpecial << "\t" << totalPos << endl;
				//}
				
				unsigned int movedSpecial = ( totalSpecial > _nSpecialHash ) ? _nSpecialHash : totalSpecial;
				if ( ( totalSpecial != totalPos ) && ( movedSpecial > 0 ) ) {
					unsigned int specialBegin = totalPos - totalSpecial;
					unsigned int normalEnd    = totalPos - totalSpecial - 1;
					random_shuffle( positions.begin(), positions.begin() + normalEnd, randomGenerator );
					random_shuffle( positions.begin() + specialBegin, positions.end(), randomGenerator );
					unsigned int temp;
					for ( unsigned int i = 0; i < totalSpecial; ++i ) {
						temp = positions[i];
						positions[i] = positions[ specialBegin + i ];
						positions[ specialBegin + i ] = temp;
					}

				} else
					random_shuffle( positions.begin(), positions.end(), randomGenerator );
			} else
				random_shuffle( positions.begin(), positions.end(), randomGenerator );
			StorePositions(curFilePosition, left, positions, offset);
		}

		positions.clear();

		offset += KEY_LENGTH;
	}

	delete [] blockPosition;

	//for ( uint64_t i = curFilePosition; i < mPositionBufferLen1; i++ ) {
	//	delete &mPositionBuffer[i];
	//}

	//cout << "finished." << endl;

	fclose(mPositions);


}

// store the positions in the memory
inline void CJumpDnaHash::StorePositions ( off_type& curFilePosition, off_type& left, vector<unsigned int>& positions, const off_type keyOffset) {
	
	//if ( mLimitPositions && (positions.size() > mMaxHashPositions) ) {
	//	cout << "ERROR: The amount of hash positions is incorrect." << endl;
	//	exit(1);
	//}
	
	// revise the key pointer
	memcpy((char*)(mKeyBufferPtr + keyOffset), (char*)&curFilePosition, KEY_LENGTH);
	

	// mPositionBuffer is full
	if ( left < SIZEOF_INT ) {
		cout << "ERROR: Run out the allocated position memory." << endl;
		exit(1);
	}

	// store number of hash hits
	unsigned int nPositions = 0;
	if ( mLimitPositions && (positions.size() > mMaxHashPositions) )
		nPositions = mMaxHashPositions;
	else
		nPositions = positions.size();
	
	memcpy((char*)(mPositionBufferPtr + curFilePosition), (char*)&nPositions, SIZEOF_INT);
	curFilePosition += SIZEOF_INT;
	left            -= SIZEOF_INT;

	for ( unsigned int i = 0; i < nPositions; i++ ) {
		// mPositionBuffer is full
		if ( left < SIZEOF_INT ) {
			cout << "ERROR: Run out the allocated position memory." << endl;
			exit(1);
		}

		unsigned int position = positions[i];
		memcpy((char*)(mPositionBufferPtr + curFilePosition), (char*)&position, SIZEOF_INT);
		curFilePosition += SIZEOF_INT;
		left            -= SIZEOF_INT;
	}


}

// randomize and trim hash positions
void CJumpDnaHash::RandomizeAndTrimHashPositions(unsigned short numHashPositions) {
	mLimitPositions   = true;
	mMaxHashPositions = numHashPositions;
}

// dummy function
void CJumpDnaHash::Resize(void) {}
