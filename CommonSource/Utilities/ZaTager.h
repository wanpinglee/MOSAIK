#ifndef _ZATAGER_H_
#define _ZATAGER_H_

#include "AlignedRead.h"
#include "Alignment.h"
#include "CigarTager.h"
#include "MdTager.h"
#include "vector"

using namespace std;

class CZaTager {
	public:
		CZaTager();
		~CZaTager();
		//const char* GetZaTag( vector<Alignment>& ar1, vector<Alignment>& ar2 );
		const char* GetZaTag( const Alignment& query, const Alignment& mate, const bool& isFirstMate, const bool& isSingleton = false );
		// copy constructor
		CZaTager( CZaTager const & copy ) {
			bufferLen = copy.bufferLen;
			buffer    = new char [ bufferLen ];
			memcpy( buffer, copy.buffer, bufferLen );
			cigarTager = copy.cigarTager;
		};
		// assign operator
		CZaTager& operator=( CZaTager const & copy ) {
			 char * temp = new char [ copy.bufferLen ];
			 delete [] buffer;
			 bufferLen = copy.bufferLen;
			 buffer    = temp;
			 memcpy( buffer, copy.buffer, bufferLen );
			 cigarTager = copy.cigarTager;
			 return *this;
		};

	private:
		// our ZA tag buffer
		unsigned int bufferLen;
		char* buffer;

		// extend the buffer
		void ExtendBuffer( const unsigned int& length );

		// cigar tag
		CCigarTager cigarTager;

		// md tag
		CMdTager mdTager;
};

#endif
