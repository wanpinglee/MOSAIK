#include "ZaTager.h"

CZaTager::CZaTager( void )
	: bufferLen(1024)
	, buffer(NULL)
{
	buffer = new char [ bufferLen ];
	memset(buffer, 0, bufferLen);
}

CZaTager::~CZaTager( void ) {
	if ( buffer ) delete [] buffer;
}

void CZaTager::ExtendBuffer( const unsigned int& length ) {

	if ( buffer ) delete [] buffer;

	buffer = new char [ length + 10 ];
	memset(buffer, 0, bufferLen);
}

const char* CZaTager::GetZaTag( vector<Alignment>& ar1, vector<Alignment>& ar2 ) {
	
	char* zaPtr = buffer;
	unsigned int len = 0;
	len = sprintf( zaPtr, "1");
	zaPtr += len;
	
	for ( vector<Alignment>::iterator ite = ar1.begin(); ite != ar1.end(); ++ite ) {
		len = sprintf( zaPtr, "<");
		zaPtr += len;
		len = sprintf( zaPtr, "%s;", ite->ReferenceName );
		zaPtr += len;
		len = sprintf( zaPtr, "%u;", ite->ReferenceBegin + 1 );
		zaPtr += len;
		len = sprintf( zaPtr, "%u;", ite->Quality);
		zaPtr += len;
		char strand = ( ite->IsReverseStrand ) ? '-' : '+';
		len = sprintf( zaPtr, "%c;", strand );
		zaPtr += len;
		const char* pCigar = cigarTager.GetCigarTag( ite->Reference.CData(), ite->Query.CData(), ite->Reference.Length() );
		len = strlen( pCigar );
		memcpy( zaPtr, pCigar, len );
		zaPtr += len;
		len = sprintf( zaPtr, ">");
		zaPtr += len;
	}


	len = sprintf( zaPtr, "2");
	zaPtr += len;

	for ( vector<Alignment>::iterator ite = ar2.begin(); ite != ar2.end(); ++ite ) {
		len = sprintf( zaPtr, "<");
		zaPtr += len;
		len = sprintf( zaPtr, "%s;", ite->ReferenceName );
		zaPtr += len;
		len = sprintf( zaPtr, "%u;", ite->ReferenceBegin + 1 );
		zaPtr += len;
		len = sprintf( zaPtr, "%u;", ite->Quality);
		zaPtr += len;
		char strand = ( ite->IsReverseStrand ) ? '-' : '+';
		len = sprintf( zaPtr, "%c;", strand );
		zaPtr += len;
		const char* pCigar = cigarTager.GetCigarTag( ite->Reference.CData(), ite->Query.CData(), ite->Reference.Length() );
		len = strlen( pCigar );
		memcpy( zaPtr, pCigar, len );
		zaPtr += len;
		len = sprintf( zaPtr, ">");
		zaPtr += len;
	}
	
	*zaPtr = 0;

	return buffer;
}
