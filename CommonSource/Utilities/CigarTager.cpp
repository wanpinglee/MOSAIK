#include "CigarTager.h"

CCigarTager::CCigarTager( void )
	: bufferLen(1024)
	, buffer(NULL)
{
	buffer    = new char [ bufferLen ];
	memset(buffer, 0, bufferLen);
}

CCigarTager::~CCigarTager( void ) {
	if ( buffer )    delete [] buffer;
}

void CCigarTager::ExtendBuffer( const unsigned int& length ) {
	if ( buffer ) delete [] buffer;

	buffer = new char [ length + 10 ];
	memset(buffer, 0, bufferLen);
}

inline void CCigarTager::InitializeVar( const char* reference, const char* query, const unsigned int& referenceLen ) {
	pCigar  = buffer;
	pReference = (char*) reference;
	pQuery     = (char*) query;

	numBases       = referenceLen;
	currentPos     = 0;
	numBufferBytes = 0;
}

const char* CCigarTager::GetCigarTag( const char* reference, const char* query, const unsigned int& referenceLen ) {
	
	// check the buffer size
	if ( ( referenceLen * 2 ) > bufferLen )
		ExtendBuffer( referenceLen * 2 );

	InitializeVar( reference, query, referenceLen );

	while(currentPos < numBases) {

		unsigned short testPos = currentPos;
		unsigned short operationLength = 0;
		int numWritten = 0;

		if( (pReference[currentPos] != '-') && (pQuery[currentPos] != '-') && (pReference[currentPos] != 'Z') ) {

			while((pReference[testPos] != '-') && (pQuery[testPos] != '-') && (pReference[testPos] != 'Z') && (testPos < numBases)) {
				++testPos;					
				++operationLength;
			}

			numWritten = sprintf_s(pCigar, bufferLen, "%uM", operationLength);

		} else if ( pReference[currentPos] == 'Z' ) {

			while( ( pReference[testPos] == 'Z' ) && ( testPos < numBases ) ){
				++testPos;
				++operationLength;
			}

			numWritten = sprintf_s(pCigar, bufferLen, "%uS", operationLength);

		
		} else if(pReference[currentPos] == '-') {

			while((pReference[testPos] == '-') && (testPos < numBases)) {
				++testPos;					
				++operationLength;
			}

			numWritten = sprintf_s(pCigar, bufferLen, "%uI", operationLength);

		} else if(pQuery[currentPos] == '-') {

			while((pQuery[testPos] == '-') && (testPos < numBases)) {
				++testPos;					
				++operationLength;
			}

			numWritten = sprintf_s(pCigar, bufferLen, "%uD", operationLength);

		} else {
			cout << "ERROR: CIGAR string generation failed." << endl;
			exit(1);
		}

		// increment our position
		pCigar         += numWritten;
		numBufferBytes += numWritten;
		currentPos     += operationLength;

		// make sure aren't creating a buffer overflow
		if(numBufferBytes >= bufferLen) {
			printf("ERROR: buffer overflow detected when creating the cigar string.\n");
			exit(1);
		}
	}

	*pCigar = 0;

	return buffer;
}
