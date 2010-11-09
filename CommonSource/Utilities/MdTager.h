
#ifndef _MdTager_H_
#define _MdTager_H_

#include <stdlib.h>
#include <string.h>
#include "FileUtilities.h"

class CMdTager {
	public:
		CMdTager();
		char* GetMdTag( const char* reference, const char* query, const unsigned int& referenceLen );
	private:
		// our MD tag buffer
		unsigned int bufferLen;
		char* buffer;
		char* tempBases;

		// extend the buffer
		void ExtendBuffer( const unsigned int& length );
		// initialize used variable
		inline void InitializeVar( const char* reference, const char* query, const unsigned int& referenceLen );

		// our used variables
		char* pMd;
		char* pReference;
		char* pQuery;
		unsigned int numBases;
		unsigned int currentPos;
		unsigned int numBufferBytes;
		char zeroChar;
};

#endif
