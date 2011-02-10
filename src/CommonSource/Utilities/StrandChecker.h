
#ifndef _STRANDCHECKER_H_
#define _STRANDCHECKER_H_

#include "stdio.h"
#include "Alignment.h"
#include "SequencingTechnologies.h"

inline bool isProperOrientation ( 
	const bool& reverseStrandQuery,
	const bool& reverseStrandMate,
	const unsigned int& referenceBeginQuery,
	const unsigned int& referenceBeginMate,
	const bool& isQueryFirstMate,
	const SequencingTechnologies& tech) {

	//bool proper = true;
	//
	//fprintf(stderr, "%s; %s; %s; %u; %u\n", (isQueryFirstMate?"true":"false"), (reverseStrandQuery?"true":"false"), (reverseStrandMate?"true":"false"), referenceBeginQuery, referenceBeginMate );

	switch ( tech ) {
		
		case 1: // 454
			if ( reverseStrandQuery != reverseStrandMate ) return false;
			if ( isQueryFirstMate ) {
				if ( reverseStrandQuery && ( referenceBeginQuery < referenceBeginMate) ) return true;
				if ( !reverseStrandQuery && ( referenceBeginQuery > referenceBeginMate) ) return true;
			} else {
				if ( reverseStrandQuery && ( referenceBeginQuery > referenceBeginMate) ) return true;
				if ( !reverseStrandQuery && ( referenceBeginQuery < referenceBeginMate) ) return true;
			}

			return false;
		break;


		case 16: // SOLiD
			if ( reverseStrandQuery != reverseStrandMate ) return false;
			if ( isQueryFirstMate ) {
				if ( reverseStrandQuery && ( referenceBeginQuery > referenceBeginMate) ) return true;
				if ( !reverseStrandQuery && ( referenceBeginQuery < referenceBeginMate) ) return true;
			} else {
				if ( reverseStrandQuery && ( referenceBeginQuery < referenceBeginMate) ) return true;
				if ( !reverseStrandQuery && ( referenceBeginQuery > referenceBeginMate) ) return true;
			}

			return false;
		break;

		default:
			if ( reverseStrandQuery == reverseStrandMate ) return false;
			//if ( isQueryFirstMate ) {
				if ( reverseStrandQuery && ( referenceBeginQuery > referenceBeginMate) ) return true;
				if ( !reverseStrandQuery && ( referenceBeginQuery < referenceBeginMate) ) return true;
			//} else {
				//if ( reverseStrandQuery && ( referenceBeginQuery > referenceBeginMate) ) return true;
				//if ( !reverseStrandQuery && ( referenceBeginQuery < referenceBeginMate) ) return true;
			//}

			return false;
		break;
	}


	//fprintf(stderr, "%u; %u; %u; %s; %s; %s;\n", tech, referenceBeginMate1, referenceBeginMate2, (reverseStrandMate1?"true":"false"), (reverseStrandMate2?"true":"false"), (proper?"true":"false") );

	//return proper;

}

#endif
