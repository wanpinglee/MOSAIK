#ifndef _BestNSecondBestSelection_H_
#define _BestNSecondBestSelection_H_


#include <limits.h>

#include <algorithm>
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
		const SequencingTechnologies& tech);

	void Select ( vector<Alignment>& mate1Set, vector<Alignment>& mate2Set, const unsigned int expectedFragmentLength,
		const SequencingTechnologies& tech, const bool& considerMate1 = true, const bool& considerMate2 = true );
}

#endif
