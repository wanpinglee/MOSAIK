#include "ZaTager.h"

CZaTager::CZaTager( void )
	: bufferLen(512)
	, buffer(NULL)
	, cigarTager()
	, mdTager()
{
	buffer = new char [ bufferLen ];
	memset(buffer, 0, bufferLen);
	//buffer.resize(512);
}

CZaTager::~CZaTager( void ) {
	if ( buffer ) delete [] buffer;
	buffer = NULL;
}

void CZaTager::ExtendBuffer( const unsigned int& length ) {

	char* newBuffer = new char [ length ];
	memcpy( newBuffer, buffer, bufferLen );
	delete [] buffer;
	buffer = newBuffer;
	
	bufferLen = length;
}

const char* CZaTager::GetZaTag( const Alignment& query, const Alignment& mate, const bool& isFirstMate, const bool& isSingleton, const bool& isMateUnmapped ) {
	
	char* zaPtr = buffer;
	unsigned int len = 0;

	Alignment al1, al2;
	if ( isFirstMate ) {
		al1 = query;
		al2 = mate;
	} else {
		al1 = mate;
		al2 = query;
	}

	if ( !isSingleton || ( isSingleton && isFirstMate ) ) {
	// read 1
	len = sprintf( zaPtr, "<");
	zaPtr += len;
	if ( isFirstMate )
		len = sprintf( zaPtr, "@;");
	else
		len = sprintf( zaPtr, "&;");
	zaPtr += len;
	len = sprintf( zaPtr, "%u;%u;", al1.Quality, al1.NextBestQuality );
	zaPtr += len;
	if ( !al1.SpecialCode.empty() )
		len = sprintf( zaPtr, "%s;%u;", al1.SpecialCode.c_str(), al1.NumMapped );
	else
		len = sprintf( zaPtr, ";%u;", al1.NumMapped );
	zaPtr += len;
	if ( isFirstMate || isMateUnmapped ) {
		len = sprintf( zaPtr, ";>" );
		zaPtr += len;
	}
	else { 
		len = sprintf( zaPtr, "%s;", cigarTager.GetCigarTag( al1.Reference.CData(), al1.Query.CData(), al1.Reference.Length() ) );
		zaPtr += len;
		len = sprintf( zaPtr, "%s>", mdTager.GetMdTag( al1.Reference.CData(), al1.Query.CData(), al1.Reference.Length() ) );
		zaPtr += len;
	}
	}
	

	// read 2
	if ( !isSingleton || ( isSingleton && !isFirstMate ) ) {
	len = sprintf( zaPtr, "<");
	zaPtr += len;
	if ( !isFirstMate )
		len = sprintf( zaPtr, "@;");
	else
		len = sprintf( zaPtr, "&;");
	zaPtr += len;
	len = sprintf( zaPtr, "%u;%u;", al2.Quality, al2.NextBestQuality );
	zaPtr += len;
	if ( !al2.SpecialCode.empty() )
		len = sprintf( zaPtr, "%s;%u;", al2.SpecialCode.c_str(), al2.NumMapped );
	else
		len = sprintf( zaPtr, ";%u;", al2.NumMapped );
	zaPtr += len;
	if ( !isFirstMate || isMateUnmapped ) {
		len = sprintf( zaPtr, ";>" );
		zaPtr += len;
	}
	else { 
		len = sprintf( zaPtr, "%s;", cigarTager.GetCigarTag( al2.Reference.CData(), al2.Query.CData(), al2.Reference.Length() ) );
		zaPtr += len;
		len = sprintf( zaPtr, "%s>", mdTager.GetMdTag( al2.Reference.CData(), al2.Query.CData(), al2.Reference.Length() ) );
		zaPtr += len;
	}
	}
	
	*zaPtr = 0;

	return buffer;
}
