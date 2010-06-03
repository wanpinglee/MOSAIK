// ***************************************************************************
// CFastLZIO - embeds compressed data blocks into open file streams.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <new>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "fastlz.h"
#include "LargeFileSupport.h"
#include "Mosaik.h"

using namespace std;

// the buffer is currently set to 10 MB
#define FASTLZ_IO_BUFFER_LEN        10485760

// the buffer must be at least 5 % larger
#define FASTLZ_IO_OUTPUT_BUFFER_LEN  9986438

class CFastLZIO {
public:
	// constructor
	CFastLZIO(void);
	// destructor
	~CFastLZIO(void);
	// clears the buffer
	void Clear(void);
	// our input method
	void Read(char* &buffer, unsigned int& bufferLen, FILE* stm);
	// our input method (STL string)
	void Read(string& s, FILE* stm);
	// our output method
	void Write(const char* buffer, const unsigned int bufferLen, FILE* stm);
private:
	// our buffer
	char* mBuffer;
	unsigned int mBufferLen;
};
