#ifndef _BestNSecondBestSelection_H_
#define _BestNSecondBestSelection_H_


#include <limits.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "Alignment.h"
#include "SequencingTechnologies.h"


using namespace std;

namespace BestNSecondBestSelection {

	inline bool IsBetterPair ( 
		const Alignment& competitor_mate1, 
		const Alignment& competitor_mate2,
	        const unsigned int competitor_fragmentLength, 
		const Alignment& mate1,
		const Alignment& mate2, 
		const unsigned int fragmentLength,
		const unsigned int expectedFragmentLength,
		const SequencingTechnologies& tech,
		const unsigned int numMate1Bases,
		const unsigned int numMate2Bases);

	void Select ( 
		vector<Alignment>& mate1Set, 
		vector<Alignment>& mate2Set, 
		const unsigned int expectedFragmentLength,
		const SequencingTechnologies& tech,
		const unsigned int numMate1Bases,
		const unsigned int numMate2Bases,
		const bool& considerMate1 = true,
		const bool& considerMate2 = true );
		//const unsigned int highestSwScoreMate1 = 0,
		//const unsigned int highestSwScoreMate2 = 0);

	void CalculateFragmentLength( 
		const Alignment& al1, 
		const Alignment& al2, 
		const SequencingTechnologies& tech, 
		unsigned int& length );
}

#endif
