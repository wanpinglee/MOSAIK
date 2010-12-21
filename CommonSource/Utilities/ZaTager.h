#ifndef _ZATAGER_H_
#define _ZATAGER_H_

#include "AlignedRead.h"
#include "Alignment.h"
#include "CigarTager.h"
#include "vector"

using namespace std;

class CZaTager {
	public:
		CZaTager();
		~CZaTager();
		const char* GetZaTag( vector<Alignment>& ar1, vector<Alignment>& ar2 );
		const char* GetZaTag( const Alignment& query, const Alignment& mate, const bool& isFirstMate );

	private:
		// our ZA tag buffer
		unsigned int bufferLen;
		char* buffer;

		// extend the buffer
		void ExtendBuffer( const unsigned int& length );

		// cigar tag
		CCigarTager cigarTager;
};

#endif
