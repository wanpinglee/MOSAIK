// ***************************************************************************
// CJumpCreator - creates a jump database for use with MosaikAligner.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "JumpCreator.h"

// constructor
CJumpCreator::CJumpCreator(const unsigned char hashSize, const string& filenameStub, const unsigned char sortingMemoryGB, const bool keepKeysInMemory, const unsigned int hashPositionThreshold)
: mHashSize(hashSize)
, mSortingMemoryGB(sortingMemoryGB)
, mKeys(NULL)
, mPositions(NULL)
, mBuffer(NULL)
, mBufferLen(4096)
, mNumHashPositions(0)
, mMaxHashPositions(0)
, mKeepKeysInMemory(keepKeysInMemory)
, mLimitPositions(false)
, mLogHashPositions(false)
, mKeyBuffer(NULL)
, mKeyBufferLen(0)
{
	// initialize the file buffer
	try {
		mBuffer = new unsigned char[mBufferLen];
	} catch(bad_alloc) {
		cout << "ERROR: Unable to allocate enough memory for the jump database buffer." << endl;
		exit(1);
	}

	// initialize the key buffer
	mKeyBufferLen = (uint64_t)(pow(4.0, (double)mHashSize) * KEY_LENGTH);

	if(keepKeysInMemory) {
		try {
			unsigned int num64uint = (unsigned int)(mKeyBufferLen / (double)SIZEOF_UINT64);
			mKeyBuffer = new uint64_t[num64uint];
			const uint64_t EMPTY = 0xffffffffffffffffULL;
			uninitialized_fill(mKeyBuffer, mKeyBuffer + num64uint, EMPTY);
		} catch(bad_alloc) {
			cout << "ERROR: Unable to allocate enough memory for the jump database key buffer." << endl;
			exit(1);
		}
	}

	string keyFilename       = filenameStub + "_keys.jmp";
	string metaFilename      = filenameStub + "_meta.jmp";
	string positionsFilename = filenameStub + "_positions.jmp";

	fopen_s(&mKeys, keyFilename.c_str(), "wb");

	if(!mKeys) {
		cout << "ERROR: Unable to open the keys file (" << keyFilename << ") for writing." << endl;
		exit(1);
	}

	FILE* meta = NULL;
	fopen_s(&meta, metaFilename.c_str(), "wb");

	if(!meta) {
		cout << "ERROR: Unable to open the metadata file (" << metaFilename << ") for writing." << endl;
		exit(1);
	}

	putc(hashSize, meta);
	fclose(meta);

	fopen_s(&mPositions, positionsFilename.c_str(), "wb");

	if(!mPositions) {
		cout << "ERROR: Unable to open the positions file (" << positionsFilename << ") for writing." << endl;
		exit(1);
	}

	if(hashPositionThreshold > 0) {
		cout << "- setting hash position threshold to " << hashPositionThreshold << endl;
		mLimitPositions   = true;
		mMaxHashPositions = hashPositionThreshold;
	}
}

// destructor
CJumpCreator::~CJumpCreator(void) {

	if(mBuffer)                         delete [] mBuffer;
	if(mKeepKeysInMemory && mKeyBuffer) delete [] mKeyBuffer;

	// close the jump database files
	fclose(mKeys);
	fclose(mPositions);
	//if(mLogHashPositions) gzclose(mHashPositionLog);

	// delete our temporary files
	for(unsigned int i = 0; i < mSerializedPositionsFilenames.size(); i++)
		rm(mSerializedPositionsFilenames[i].c_str());
}

// builds the jump database
void CJumpCreator::BuildJumpDatabase(void) {

	// ------------------------------
	// create the jump database files
	// ------------------------------

	if(!mKeepKeysInMemory) {

		uint64_t bytesLeft = mKeyBufferLen;

		unsigned int fillBufferSize = 314572800; // 300 MB
		unsigned int numBuffers = (unsigned int)(mKeyBufferLen / (double)fillBufferSize);
		unsigned int currentBuffer = 0;

		CConsole::Heading(); 
		cout << endl << "- creating jump keys database (" << numBuffers << " blocks):" << endl;
		CConsole::Reset(); 

		char* fillBuffer = new char[fillBufferSize];
		uninitialized_fill(fillBuffer, fillBuffer + fillBufferSize, 0xff);

		bool isRunning = true;
		CProgressCounter<unsigned int>::StartThread(&currentBuffer, &isRunning, "blocks");

		for(; currentBuffer < numBuffers; currentBuffer++) {
			fwrite(fillBuffer, fillBufferSize, 1, mKeys);
			bytesLeft -= fillBufferSize;
		}

		fwrite(fillBuffer, (size_t)bytesLeft, 1, mKeys);

		// stop the progress counter
		isRunning = false;
		CProgressCounter<unsigned int>::WaitThread();

		// clean up
		delete [] fillBuffer;
	}

	// ---------------------------------------
	// open all of the temporary sorting files
	// ---------------------------------------

	unsigned int numSortingFiles = mSerializedPositionsFilenames.size();
	FILE** sortHandles = new FILE*[numSortingFiles];

	for(unsigned int i = 0; i < numSortingFiles; i++) {
		sortHandles[i] = NULL;
		fopen_s(&sortHandles[i], mSerializedPositionsFilenames[i].c_str(), "rb");

		if(!sortHandles[i]) {
			cout << "ERROR: Unable to open temporary file (" << mSerializedPositionsFilenames[i] << ") for reading." << endl;
			exit(1);
		}
	}

	// -------------------------
	// initialize sorting vector
	// -------------------------

	double memoryAllocated = mSortingMemoryGB * 1073741824.0;
	unsigned int maxSortingElements = (int)(memoryAllocated / (double)sizeof(HashPosition));

	vector<HashPosition> hashCounts;
	hashCounts.reserve(maxSortingElements);

	// ---------------
	// get the top row
	// ---------------

	vector<HashPosition> sameHash;
	vector<HashPosition> topRow;
	topRow.reserve(numSortingFiles);

	for(unsigned int i = 0; i < numSortingFiles; i++) {
		HashPosition hp;
		hp.Owner = i;
		if(hp.Deserialize(sortHandles[i])) topRow.push_back(hp);
	}

	sort(topRow.begin(), topRow.end(), SortHashPositionDesc());

	// -----------------
	// process the files
	// -----------------

	CConsole::Heading(); 
	cout << endl << "- writing jump positions database:" << endl;
	CConsole::Reset(); 

	unsigned int numProcessed = 0;
	CProgressBar<unsigned int>::StartThread(&numProcessed, 0, mNumHashPositions, "hash positions");

	while(true) {

		// no more hash positions
		if(topRow.empty()) break;
		HashPosition bestPosition = topRow.back();
		topRow.pop_back();

		if(sameHash.empty()) {

			sameHash.push_back(bestPosition);

		} else {

			if(bestPosition.Hash == sameHash[0].Hash) {
				sameHash.push_back(bestPosition);
			} else {
				numProcessed += sameHash.size();
				StoreHash(sameHash);

				sameHash.clear();
				sameHash.push_back(bestPosition);
			}
		}

		// get the next hash position from the appropriate temp file
		{
			HashPosition hp;
			hp.Owner = bestPosition.Owner;
			if(hp.Deserialize(sortHandles[bestPosition.Owner])) topRow.push_back(hp);
		}

		sort(topRow.begin(), topRow.end(), SortHashPositionDesc());
	}

	// store the last hash
	numProcessed += sameHash.size();
	StoreHash(sameHash);

	// stop the progress bar
	CProgressBar<unsigned int>::WaitThread();

	// write the keys to file
	if(mKeepKeysInMemory) {

		uint64_t bytesLeft = mKeyBufferLen;

		unsigned int fillBufferSize = 314572800 ; // 300 MB
		unsigned int numBuffers = (unsigned int)(mKeyBufferLen / (double)fillBufferSize);
		unsigned int currentBuffer = 0;
		char* pKeys = (char*)mKeyBuffer;

		CConsole::Heading(); 
		cout << endl << "- serializing jump keys database (" << numBuffers << " blocks):" << endl;
		CConsole::Reset(); 

		bool isRunning = true;
		CProgressCounter<unsigned int>::StartThread(&currentBuffer, &isRunning, "blocks");

		for(; currentBuffer < numBuffers; currentBuffer++) {
			fwrite(pKeys, fillBufferSize, 1, mKeys);
			pKeys     += fillBufferSize;
			bytesLeft -= fillBufferSize;
		}

		fwrite(pKeys, (size_t)bytesLeft, 1, mKeys);

		// stop the progress counter
		isRunning = false;
		CProgressCounter<unsigned int>::WaitThread();
	}

	// close our files
	for(unsigned int i = 0; i < numSortingFiles; i++) fclose(sortHandles[i]);

	// clean up
	delete [] sortHandles;
}

// creates the hash for a supplied fragment
void CJumpCreator::CreateHash(const char* fragment, const unsigned char fragmentLen, uint64_t& key) {

	// set the key to zero
	key = 0;
	const char translation[26] = { 0, 3, 1, 3, -1, -1, 2, 3, -1, -1, 3, -1, 0, 3, -1, -1, -1, 0, 2, 3, -1, 0, 3, 1, 3, -1 };

	if(fragmentLen > 32) {
		cout << "ERROR: This hash table can only handle fragments smaller or equal to 32 bases." << endl;
		exit(1);
	}	

	// convert each nucleotide to its 2-bit representation
	for(unsigned char i = 0; i < fragmentLen; i++) {

		// convert [A,C,G,T] to [0,1,2,3]
		char tValue = translation[fragment[i] - 'A'];

		// catch any unrecognized nucleotides
		if(tValue < 0) {
			cout << "ERROR: Unrecognized nucleotide in hash table: " << fragment[i] << endl;
			cout << "- fragment: ";
			for(unsigned j = 0; j < fragmentLen; j++) cout << fragment[j];
			cout << endl;
			exit(1);
		}

		// shift the key and add the new value
		key = key << 2 | tValue;
	}
}

// enables hash position logging
//void CJumpCreator::EnableHashPositionsLogging(const string& filename) {
//
//	mLogHashPositions = true;
//
//	mHashPositionLog = gzopen(filename.c_str(), "wb3");
//
//	if(!mHashPositionLog) {
//		cout << "ERROR: Unable to open the hash positions log (" << filename << ") for writing." << endl;
//		exit(1);
//	}
//}

// hashes the reference and stores the results in sorted temporary files
void CJumpCreator::HashReference(const string& referenceFilename) {

	// --------------------------------------------
	// retrieve the concatenated reference sequence
	// --------------------------------------------

	cout << endl << "- retrieving reference sequence... ";
	cout.flush();

	MosaikReadFormat::CReferenceSequenceReader refseq;
	refseq.Open(referenceFilename);
	unsigned int referenceLength = refseq.GetReferenceSequenceLength();

	char* pReference = NULL;
	refseq.LoadConcatenatedSequence(pReference);

	cout << "finished." << endl << endl;
	refseq.Close();

	// -------------------------
	// initialize sorting vector
	// -------------------------

	double memoryAllocated = mSortingMemoryGB * 1073741824.0;
	unsigned int maxSortingElements = (int)(memoryAllocated / (double)sizeof(HashPosition));

	vector<HashPosition> hashPositions;
	hashPositions.reserve(maxSortingElements);

	// ----------------------------------------
	// hash the concatenated reference sequence
	// ----------------------------------------

	CConsole::Heading(); 
	cout << "- hashing reference sequence:" << endl;
	CConsole::Reset(); 

	char* pAnchor = pReference;

	unsigned int maxPositions = referenceLength - mHashSize + 1;
	unsigned int i = 0;
	mNumHashPositions = 0;

	CProgressBar<unsigned int>::StartThread(&i, 0, maxPositions, "hashes");

	for(; i < maxPositions; i++, pAnchor++) {

		bool skipHash = false;

		for(unsigned int j = 0; j < mHashSize; j++) {
			char anchorChar = *(pAnchor + j);
			// E, F, I, L, O are used for SOLiD reference
			skipHash = (anchorChar == 'J') || (anchorChar == 'X') || (anchorChar == 'N') || (anchorChar == 'E') || (anchorChar == 'F') || (anchorChar == 'I') || (anchorChar == 'L') || (anchorChar == 'O');
			if( skipHash ) {
				break;
			}
		}

		if( skipHash ) continue;

		HashPosition hp;
		hp.Position = i;
		CreateHash(pAnchor, mHashSize, hp.Hash);
		hashPositions.push_back(hp);

		// dump our sorting vector
		if(hashPositions.size() >= maxSortingElements) {
			mNumHashPositions += hashPositions.size();
			SerializeSortingVector(hashPositions);
		}
	}

	CProgressBar<unsigned int>::WaitThread();

	cout << endl << "- serializing final sorting vector... ";
	cout.flush();

	mNumHashPositions += hashPositions.size();
	SerializeSortingVector(hashPositions);

	cout << "finished." << endl;

	// clean up
	delete [] pReference;
}

// serializes the sorting vector to temporary files
void CJumpCreator::SerializeSortingVector(vector<HashPosition>& hashPositions) {

	// retrieve a temporary filename
	string tempFilename;
	CFileUtilities::GetTempFilename(tempFilename);

	mSerializedPositionsFilenames.push_back(tempFilename);

	// open the temporary file
	FILE* temp = NULL;
	fopen_s(&temp, tempFilename.c_str(), "wb");

	if(!temp) {
		cout << "ERROR: Unable to open temporary file (" << tempFilename << ") for writing." << endl;
		exit(1);
	}

	// sort the hash positions
	sort(hashPositions.begin(), hashPositions.end(), SortHashPositionAsc());

	// serialize
	for(unsigned int i = 0; i < hashPositions.size(); i++) 
		hashPositions[i].Serialize(temp);

	// close the temporary file
	fclose(temp);

	// clear the vector
	hashPositions.clear();
}

// stores the supplied hash positions in the jump database
void CJumpCreator::StoreHash(vector<HashPosition>& hashPositions) {

	unsigned int numHashes = hashPositions.size();

	// store the hash positions
	//if(mLogHashPositions) gzwrite(mHashPositionLog, (char*)&numHashes, SIZEOF_INT);

	if(numHashes == 0) {
		cout << "ERROR: Tried to store an empty hash." << endl;
		exit(1);
	}

	// limit the number of hashes that will be written
	if(mLimitPositions && (numHashes > mMaxHashPositions)) numHashes = mMaxHashPositions;

	// shuffle the vector
	random_shuffle(hashPositions.begin(), hashPositions.end());

	// localize the hash
	uint64_t hash = hashPositions[0].Hash;

	// write the position file offset in the keys file
	off_type offset = hash * KEY_LENGTH;
	off_type positionStart = ftell64(mPositions);

	if(mKeepKeysInMemory) {

		unsigned int num64uint = (unsigned int)(offset / (double)SIZEOF_UINT64);
		uint64_t* p64uint = &mKeyBuffer[num64uint];

		unsigned int numBytes = (unsigned int)(offset - num64uint * SIZEOF_UINT64);
		char* pKeys = (char*)p64uint;

		memcpy(pKeys + numBytes, (char*)&positionStart, KEY_LENGTH);

	} else {

		fseek64(mKeys, offset, SEEK_SET);
		fwrite((char*)&positionStart, KEY_LENGTH, 1, mKeys);
	}

	// write the hash positions
	unsigned int entrySize = (numHashes + 1) * SIZEOF_INT;
	CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, entrySize);

	unsigned int bufferOffset = 0;
	memcpy(mBuffer + bufferOffset, (char*)&numHashes, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	for(unsigned int i = 0; i < numHashes; i++) {
		memcpy(mBuffer + bufferOffset, (char*)&hashPositions[i].Position, SIZEOF_INT);
		bufferOffset += SIZEOF_INT;
	}

	fwrite(mBuffer, bufferOffset, 1, mPositions);
}
