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
CJumpDnaHash::CJumpDnaHash(const unsigned char hashSize, const string& filenameStub, const unsigned short numPositions, const bool keepKeysInMemory, const bool keepPositionsInMemory, const unsigned int numCachedElements, const unsigned int begin, const unsigned int end, const unsigned int offset)
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
	if(mBuffer)         delete [] mBuffer;
	if(mKeyBuffer)      delete [] mKeyBuffer;
	if(mPositionBuffer) delete [] mPositionBuffer;

	hasKeysNPositions = false;
}

// load hash keys and positions form file to memory
void CJumpDnaHash::LoadKeysNPositions() {
	if(mKeepKeysInMemory)      LoadKeys();
	if(mKeepPositionsInMemory) LoadPositions();

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
		if(mLimitPositions && (numPositions > mMaxHashPositions)) {
			mhpOccupancy = (double)mMaxHashPositions / (double)numPositions;
			numPositions = mMaxHashPositions;
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
void CJumpDnaHash::GetHashStatistics(const vector<pair<unsigned int, unsigned int> > referenceSequences, vector<unsigned int>& nHashs) {
	
	LoadKeys();

	const unsigned int fillBufferSize  = 536870912; // 500 MB
	char*        blockPosition;
	uintptr_t    blockPositionPtr;
	unsigned int nBlock = 0; // indicates how many blocks have been handled
	blockPosition = new char[(size_t)fillBufferSize];

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

	cout << "- loading jump positions database into memory... ";
	cout.flush();

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
		if ( positions.size() != 0 ) {
			random_shuffle( positions.begin(), positions.end() );
			//StorePositions(curFilePosition, left, positions, offset);
			SetPositionDistribution(referenceSequences, nHashs, positions);
		}

		positions.clear();

		offset += KEY_LENGTH;
	}

	delete [] blockPosition;

	cout << "finished." << endl;

	fclose(mPositions);


}

// determine the chromosome which positions locating in
void CJumpDnaHash::SetPositionDistribution(const vector<pair<unsigned int, unsigned int> > referenceSequences, vector<unsigned int>& nHashs, const vector<unsigned int> positions) {

	for ( unsigned int i = 0; i < positions.size(); i++ ) {
		unsigned int refNo = 0;
		while( positions[i] > referenceSequences[refNo].second ) refNo++;
		nHashs[refNo]++;
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

	cout << "- loading jump keys database into memory... ";
	cout.flush();

	char* pKeys = (char*)mKeyBuffer;
	while(bytesLeft > fillBufferSize) {
		fread(pKeys, fillBufferSize, 1, mKeys);
		pKeys     += fillBufferSize;
		bytesLeft -= fillBufferSize;
	}

	fread(pKeys, (size_t)bytesLeft, 1, mKeys);
	cout << "finished." << endl;

	fclose(mKeys);
}

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
	// TODO: we may ask too large memory space
	uint64_t mPositionBufferLen1 = _end - _begin + 1;
	if ( ( mPositionBufferLen / ( 2 * SIZEOF_INT) ) < mPositionBufferLen1 )
		mPositionBufferLen1 = mPositionBufferLen;
	else
		mPositionBufferLen1 = mPositionBufferLen1 * 2 * SIZEOF_INT;

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

	cout << "- loading jump positions database into memory... ";
	cout.flush();

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
			random_shuffle( positions.begin(), positions.end() );
			StorePositions(curFilePosition, left, positions, offset);
		}

		positions.clear();

		offset += KEY_LENGTH;
	}

	delete [] blockPosition;

	cout << "finished." << endl;

	fclose(mPositions);


}

// store the positions in the memory
inline void CJumpDnaHash::StorePositions ( off_type& curFilePosition, off_type& left, vector<unsigned int>& positions, const off_type keyOffset) {
	
	// revise the key pointer
	memcpy((char*)(mKeyBufferPtr + keyOffset), (char*)&curFilePosition, KEY_LENGTH);
	

	// mPositionBuffer is full
	if ( left < SIZEOF_INT ) {
		cout << "ERROR: Out of the allocated position memory." << endl;
		exit(1);
	}

	// store number of hash hits
	unsigned int nPositions = positions.size();
	memcpy((char*)(mPositionBufferPtr + curFilePosition), (char*)&nPositions, SIZEOF_INT);
	curFilePosition += SIZEOF_INT;
	left            -= SIZEOF_INT;
	
	for ( vector<unsigned int>::iterator ptr = positions.begin(); ptr != positions.end(); ptr++ ) {
		// mPositionBuffer is full
		if ( left < SIZEOF_INT ) {
			cout << "ERROR: Out of the allocated position memory." << endl;
			exit(1);
		}

		unsigned int position = *ptr;
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
