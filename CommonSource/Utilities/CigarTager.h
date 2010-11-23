
#ifndef _CigarTager_H_
#define _CigarTager_H_

#include <stdlib.h>
#include <string.h>
#include "FileUtilities.h"

class CCigarTager {
	public:
		CCigarTager();
		~CCigarTager();
		const char* GetCigarTag( const char* reference, const char* query, const unsigned int& referenceLen );
	private:
		// our cigar buffer
		unsigned int bufferLen;
		char* buffer;

		// extend the buffer
		void ExtendBuffer( const unsigned int& length );
		// initialize used variable
		inline void InitializeVar( const char* reference, const char* query, const unsigned int& referenceLen );

		// our used variables
		char* pCigar;
		char* pReference;
		char* pQuery;
		unsigned int numBases;
		unsigned int currentPos;
		unsigned int numBufferBytes;
};

#endif
