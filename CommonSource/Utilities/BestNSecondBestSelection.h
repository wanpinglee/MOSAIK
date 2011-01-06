#ifndef _BestNSecondBestSelection_H_
#define _BestNSecondBestSelection_H_


#include <limits.h>

#include <algorithm>
#include <vector>

#include "Alignment.h"


using namespace std;

namespace BestNSecondBestSelection {

	inline bool LessThanMQ ( const Alignment& al1, const Alignment& al2);
	inline bool IsBetterPair ( 
		const Alignment& competitor_mate1, 
		const Alignment& competitor_mate2,
	        const unsigned int competitor_fragmentLength, 
		const Alignment& mate1,
		const Alignment& mate2, 
		const unsigned int fragmentLength,
		const unsigned int expectedFragmentLength );

	void Select ( vector<Alignment>& mate1Set, vector<Alignment>& mate2Set, const unsigned int expectedFragmentLength );
}

#endif
