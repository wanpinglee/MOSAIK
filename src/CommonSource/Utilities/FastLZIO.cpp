// ***************************************************************************
// CFastLZIO - embeds compressed data blocks into open file streams.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "FastLZIO.h"

#define FASTLZ_BETTER_SPEED       1
#define FASTLZ_BETTER_COMPRESSION 2

// constructor
CFastLZIO::CFastLZIO(void)
: mBuffer(NULL)
, mBufferLen(0)
{
	try {
		mBufferLen = FASTLZ_IO_BUFFER_LEN;
		mBuffer = new char[mBufferLen];
	} catch(bad_alloc) {
		printf("ERROR: Unable to initialize the FastLZ buffer.\n");
		exit(1);
	}
}

// destructor
CFastLZIO::~CFastLZIO(void) {
	Clear();
}

// clears the buffer
void CFastLZIO::Clear(void) {
	mBufferLen = 0;
	if(mBuffer) {
		delete [] mBuffer;
		mBuffer = NULL;
	}
}

// our input method
void CFastLZIO::Read(char* &buffer, unsigned int& bufferLen, FILE* stm) {

	// read the buffer length
	unsigned int newBufferLen;
	fread((char*)&newBufferLen, SIZEOF_INT, 1, stm);

	// allocate memory if necessary
	if(newBufferLen > bufferLen) {
		try {
			bufferLen = newBufferLen;
			if(buffer) delete [] buffer;
			buffer = new char[bufferLen + 1];
		} catch(bad_alloc) {
			printf("ERROR: Unable to initialize the uncompressed FastLZ buffer (Read).\n");
			exit(1);
		}
	}

	// calculate the number of blocks
	const unsigned short numBlocksRead = (unsigned short)ceil(bufferLen / (double)FASTLZ_IO_OUTPUT_BUFFER_LEN);

	// read each block
	char* pBuffer = buffer;
	int numCompressedBytes;

	for(unsigned int i = 0; i < numBlocksRead; ++i) {

		// read the block length
		fread((char*)&numCompressedBytes, SIZEOF_INT, 1, stm);

		// read the compressed block
		fread(mBuffer, numCompressedBytes, 1, stm);

		// uncompress the block
		int numUncompressedBytes = fastlz_decompress(mBuffer, numCompressedBytes, (void*)pBuffer, FASTLZ_IO_BUFFER_LEN);
		pBuffer += numUncompressedBytes;
	}

	// add the null termination
	*pBuffer = 0;
}

// our input method (STL string)
void CFastLZIO::Read(string& s, FILE* stm) {

	// read the buffer length
	unsigned int bufferLen;
	fread((char*)&bufferLen, SIZEOF_INT, 1, stm);

	// handle an empty packet
	if(bufferLen == 0) {
		s.clear();
		return;
	}

	// calculate the number of blocks
	const unsigned short numBlocksRead = (unsigned short)ceil(bufferLen / (double)FASTLZ_IO_OUTPUT_BUFFER_LEN);

	// resize the string
	s.resize(bufferLen);

	// read each block
	char* pBuffer = (char*)s.data();
	int numCompressedBytes;

	for(unsigned int i = 0; i < numBlocksRead; ++i) {
		fread((char*)&numCompressedBytes, SIZEOF_INT, 1, stm);
		fread(mBuffer, numCompressedBytes, 1, stm);
		int numUncompressedBytes = fastlz_decompress(mBuffer, numCompressedBytes, (void*)pBuffer, FASTLZ_IO_BUFFER_LEN);
		pBuffer += numUncompressedBytes;
	}
}

// our output method
void CFastLZIO::Write(const char* buffer, const unsigned int bufferLen, FILE* stm) {

	// write the buffer length
	fwrite((char*)&bufferLen, SIZEOF_INT, 1, stm);

	// calculate the number of blocks
	const unsigned short numBlocksWritten = (unsigned short)ceil(bufferLen / (double)FASTLZ_IO_OUTPUT_BUFFER_LEN);

	// write each block
	const char* pBuffer    = buffer;
	unsigned int bytesLeft = bufferLen;

	for(unsigned int i = 0; i < numBlocksWritten; ++i) {

		// compress the block
		unsigned int numUncompressedBytes = (bytesLeft > FASTLZ_IO_OUTPUT_BUFFER_LEN ? FASTLZ_IO_OUTPUT_BUFFER_LEN : bytesLeft);
		int numCompressedBytes = fastlz_compress_level(FASTLZ_BETTER_COMPRESSION, pBuffer, numUncompressedBytes, mBuffer);
		pBuffer   += numUncompressedBytes;
		bytesLeft -= numUncompressedBytes;

		// write the block length
		fwrite((char*)&numCompressedBytes, SIZEOF_INT, 1, stm);
		//printf("uncompressed bytes: %u, compressed bytes: %u\n", numUncompressedBytes, numCompressedBytes);

		// write the compressed block
		fwrite(mBuffer, numCompressedBytes, 1, stm);
	}
}
