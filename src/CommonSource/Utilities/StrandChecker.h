
#ifndef _STRANDCHECKER_H_
#define _STRANDCHECKER_H_

#include "stdio.h"
#include "Alignment.h"
#include "SequencingTechnologies.h"

inline bool isProperOrientation ( 
	const bool& reverseStrandMate1,
	const bool& reverseStrandMate2,
	const unsigned int& referenceBeginMate1,
	const unsigned int& referenceBeginMate2,
	const SequencingTechnologies& tech) {

	//bool proper = true;

	switch ( tech ) {
		
		case 1: // 454
			if ( reverseStrandMate1 != reverseStrandMate2 ) return false;
			if ( ( referenceBeginMate1 < referenceBeginMate2 ) && reverseStrandMate1 ) return false;
			if ( ( referenceBeginMate1 > referenceBeginMate2 ) && !reverseStrandMate1 ) return false;

			return true;
		break;


		case 16: // SOLiD
			if ( reverseStrandMate1 != reverseStrandMate2 ) return false;
			if ( ( referenceBeginMate1 < referenceBeginMate2 ) && !reverseStrandMate1 ) return false;
			if ( ( referenceBeginMate1 > referenceBeginMate2 ) && reverseStrandMate1 ) return false;

			return true;
		break;

		default:
			if ( reverseStrandMate1 == reverseStrandMate2 ) return false;
			if ( ( referenceBeginMate1 < referenceBeginMate2 ) && reverseStrandMate1 ) return false;
			if ( ( referenceBeginMate1 > referenceBeginMate2 ) && !reverseStrandMate1 ) return false;

			return true;
		break;
	}


	//fprintf(stderr, "%u; %u; %u; %s; %s; %s;\n", tech, referenceBeginMate1, referenceBeginMate2, (reverseStrandMate1?"true":"false"), (reverseStrandMate2?"true":"false"), (proper?"true":"false") );

	//return proper;

}

#endif
