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
		Alignment& bestMate1,
		Alignment& bestMate2,
		vector<Alignment*>& mate1Set, 
		vector<Alignment*>& mate2Set, 
		const unsigned int& expectedFragmentLength,
		const SequencingTechnologies& tech,
		const unsigned int& numMate1Bases,
		const unsigned int& numMate2Bases,
		const bool& considerMate1,
		const bool& considerMate2,
		const bool& reserSwScore);

	void CalculateFragmentLength( 
		const Alignment& al1, 
		const Alignment& al2, 
		const SequencingTechnologies& tech, 
		unsigned int& length );
}

#endif
