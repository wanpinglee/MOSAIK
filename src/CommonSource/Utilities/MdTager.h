
#ifndef _MdTager_H_
#define _MdTager_H_

#include <stdlib.h>
#include <string.h>
#include "FileUtilities.h"

class CMdTager {
	public:
		CMdTager();
		~CMdTager();
		const char* GetMdTag( const char* reference, const char* query, const unsigned int& referenceLen );

		// copy constructor
		/*
		CMdTager ( const CMdTager & copy ) {
			bufferLen = copy.bufferLen;
			buffer    = new char [ bufferLen ];
			memcpy( buffer, copy.buffer, bufferLen );

			numBases       = copy.numBases;
			currentPos     = copy.currentPos;
			numBufferBytes = copy.numBufferBytes;
			zeroChar       = '0';
		};
		// assign operator
		CMdTager& operator=( CMdTager const & copy )
		    : bufferLen(0)
		    , buffer(NULL)
		    , tempBases(NULL)
		    , pMd(NULL)
		    , pReference(NULL)
		    , pQuery(NULL)
		    , numBases(0)
		    , currentPos(0)
		    , numBufferBytes(0)
		    , zeroChar('0')
		    {
			char * temp = new char [ copy.bufferLen ];
			delete [] buffer;
			bufferLen = copy.bufferLen;
			buffer    = temp;
			memcpy( buffer, copy.buffer, bufferLen );

			numBases       = copy.numBases;
			currentPos     = copy.currentPos;
			numBufferBytes = copy.numBufferBytes;
			zeroChar       = '0';
			return *this;
		};
		*/
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

		CMdTager (const CMdTager&);
		CMdTager& operator= (const CMdTager&);
};

#endif
