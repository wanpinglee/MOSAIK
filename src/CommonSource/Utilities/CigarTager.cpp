#include "CigarTager.h"

CCigarTager::CCigarTager( void )
		    : bufferLen(512)
		    , buffer(NULL)
		    , pCigar(NULL)
		    , pPackCigar(NULL)
		    , pReference(NULL)
		    , pQuery(NULL)
		    , numBases(0)
		    , currentPos(0)
		    , numBufferBytes(0)
{
	buffer    = new char [ bufferLen ];
	//packBuffer= new unsigned int [ bufferLen ];
	memset(buffer, 0, bufferLen);
}

CCigarTager::~CCigarTager( void ) {
	if ( buffer )     delete [] buffer;
	buffer = NULL;
	//if ( packBuffer ) delete [] packBuffer;
}

void CCigarTager::ExtendBuffer( const unsigned int& length ) {
	if ( buffer )     delete [] buffer;
	//if ( packBuffer ) delete [] packBuffer;

	bufferLen = length + 10;
	
	buffer = new char [ bufferLen ];
	//packBuffer = new unsigned int [ bufferLen ];
	memset(buffer, 0, bufferLen);
	//memset(packBuffer, 0, bufferLen);
}

inline void CCigarTager::InitializeVar( const char* reference, const char* query, const unsigned int& referenceLen ) {
	memset(buffer, 0, bufferLen);
	pCigar  = buffer;
	pPackCigar = (unsigned int*) buffer;
	pReference = (char*) reference;
	pQuery     = (char*) query;

	numBases       = referenceLen;
	currentPos     = 0;
	numBufferBytes = 0;
}

const char* CCigarTager::GetCigarTag( 
	const char* reference, 
	const char* query, 
	const unsigned int& referenceLen,
	const unsigned int& leadingClip,
	const unsigned int& laggingClip,
	const bool getPackCigar) {
	
	// check the buffer size
	if ( ( referenceLen * 2 ) > bufferLen )
		ExtendBuffer( referenceLen * 2 );

	InitializeVar( reference, query, referenceLen );

	int numWritten = 0;

	if ( leadingClip > 0 ) {
		if ( getPackCigar ) {
			*pPackCigar = leadingClip << BAM_CIGAR_SHIFT | BAM_CSOFT_CLIP;
			pPackCigar++;
		} else {
			numWritten = 0;
			numWritten = sprintf_s(pCigar, bufferLen, "%uS", leadingClip);
			pCigar += numWritten;
		}
	}
	
	while(currentPos < numBases) {

		unsigned short testPos = currentPos;
		unsigned short operationLength = 0;
		numWritten = 0;

		if( (pReference[currentPos] != '-') && (pQuery[currentPos] != '-') && (pReference[currentPos] != 'Z') ) {

			while((pReference[testPos] != '-') && (pQuery[testPos] != '-') && (pReference[testPos] != 'Z') && (testPos < numBases)) {
				++testPos;					
				++operationLength;
			}

			if ( getPackCigar )
				*pPackCigar = operationLength << BAM_CIGAR_SHIFT | BAM_CMATCH;
				
			else
				numWritten = sprintf_s(pCigar, bufferLen, "%uM", operationLength);

		} else if ( pReference[currentPos] == 'Z' ) {

			while( ( pReference[testPos] == 'Z' ) && ( testPos < numBases ) ){
				++testPos;
				++operationLength;
			}

			if ( getPackCigar )
				*pPackCigar = operationLength << BAM_CIGAR_SHIFT | BAM_CSOFT_CLIP;
			else
				numWritten = sprintf_s(pCigar, bufferLen, "%uS", operationLength);

		
		} else if(pReference[currentPos] == '-') {

			while((pReference[testPos] == '-') && (testPos < numBases)) {
				++testPos;					
				++operationLength;
			}

			if ( getPackCigar )
				*pPackCigar = operationLength << BAM_CIGAR_SHIFT | BAM_CINS;
			else
				numWritten = sprintf_s(pCigar, bufferLen, "%uI", operationLength);

		} else if(pQuery[currentPos] == '-') {

			while((pQuery[testPos] == '-') && (testPos < numBases)) {
				++testPos;					
				++operationLength;
			}

			if ( getPackCigar )
				*pPackCigar = operationLength << BAM_CIGAR_SHIFT | BAM_CDEL;
			else
				numWritten = sprintf_s(pCigar, bufferLen, "%uD", operationLength);

		} else {
			cout << "ERROR: CIGAR string generation failed." << endl;
			exit(1);
		}

		// increment our position
		if ( getPackCigar ) {
			pPackCigar++;
			numBufferBytes++;
			currentPos     += operationLength;
		} else {
			pCigar         += numWritten;
			numBufferBytes += numWritten;
			currentPos     += operationLength;
		}

		// make sure aren't creating a buffer overflow
		if(numBufferBytes >= bufferLen) {
			printf("ERROR: buffer overflow detected when creating the cigar string.\n");
			exit(1);
		}
	}


	if ( laggingClip > 0 ) {
		if ( getPackCigar ) {
			*pPackCigar = laggingClip << BAM_CIGAR_SHIFT | BAM_CSOFT_CLIP;
			pPackCigar++;
		} else {
			numWritten = 0;
			numWritten = sprintf_s(pCigar, bufferLen, "%uS", laggingClip);
			pCigar += numWritten;
		}
	}

	return buffer;
}
