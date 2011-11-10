
#ifndef CIGARTAGER_H_
#define CIGARTAGER_H_

#include <stdlib.h>
#include <string.h>
#include "Alignment.h"
#include "BamHeader.h"
#include "FileUtilities.h"

class CCigarTager {
	public:
		CCigarTager();
		~CCigarTager();
		const char* GetCigarTag( const char* reference, 
			const char* query, 
			const unsigned int& referenceLen,
			const unsigned int& leadingClip = 0,
			const unsigned int& laggingClip = 0,
			const bool packCigar = false);
		// copy constructor
		CCigarTager ( const CCigarTager & copy )
		    : bufferLen(0)
		    , buffer(NULL)
		    , pCigar(NULL)
		    , pPackCigar(NULL)
		    , pReference(NULL)
		    , pQuery(NULL)
		    , numBases(0)
		    , currentPos(0)
		    , numBufferBytes(0)
		{
			bufferLen = copy.bufferLen;
			buffer    = new char [ bufferLen ];
			memcpy( buffer, copy.buffer, bufferLen );

			numBases       = copy.numBases;
			currentPos     = copy.currentPos;
			numBufferBytes = copy.numBufferBytes;

		};
		// assign operator
		CCigarTager& operator=( const CCigarTager & copy ) {
			char * temp = new char [ copy.bufferLen ];
			delete [] buffer;
			bufferLen = copy.bufferLen;
			buffer    = temp;
			memcpy( buffer, copy.buffer, bufferLen );

			numBases       = copy.numBases;
			currentPos     = copy.currentPos;
			numBufferBytes = copy.numBufferBytes;
			return *this;
		};


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
		unsigned int* pPackCigar;
		char* pReference;
		char* pQuery;
		unsigned int numBases;
		unsigned int currentPos;
		unsigned int numBufferBytes;
};

#endif // CIGARTAGER_H_
