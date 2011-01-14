#include "BestNSecondBestSelection.h"


// greater-than operator of mapping qualities
inline bool BestNSecondBestSelection::LessThanMQ ( const Alignment& al1, const Alignment& al2){
	return al1.Quality < al2.Quality;
}

// compare the given proper pairs
inline bool BestNSecondBestSelection::IsBetterPair ( 
	const Alignment& competitor_mate1, 
	const Alignment& competitor_mate2, 
	const unsigned int competitor_fragmentLength, 
	const Alignment& mate1, 
	const Alignment& mate2, 
	const unsigned int fragmentLength,
	const unsigned int expectedFragmentLength) {

	// rescured mate always wins
	if ( competitor_mate1.WasRescued ) return true;
	if ( competitor_mate2.WasRescued ) return true;
	// proper pair always wins improper pair
	bool competitor_model = ( competitor_mate1.IsReverseStrand != competitor_mate2.IsReverseStrand ) ? true : false;
	bool current_model    = ( mate1.IsReverseStrand != mate2.IsReverseStrand ) ? true : false;
	if ( competitor_model && !current_model ) return true;
	if ( !competitor_model && current_model ) return false;

	double competitor = ( competitor_mate1.Quality + competitor_mate2.Quality ) / 100.00;
	double current    = ( mate1.Quality + mate2.Quality ) / 100.00;

	unsigned int competitor_diff = ( expectedFragmentLength > competitor_fragmentLength ) ? 
		expectedFragmentLength - competitor_fragmentLength : 
		competitor_fragmentLength - expectedFragmentLength;
	
	unsigned int diff            = ( expectedFragmentLength > fragmentLength ) ? 
		expectedFragmentLength - fragmentLength : 
		fragmentLength - expectedFragmentLength;

	competitor += ( 200 - (int)competitor_diff ) / 200;
	current    += ( 200 - (int)diff ) / 200;

	if ( competitor > current ) return true;
	else return false;
}

// Select and only keep best and 2nd best
void BestNSecondBestSelection::Select ( 
	vector<Alignment>& mate1Set, 
	vector<Alignment>& mate2Set, 
	const unsigned int expectedFragmentLength,
	const bool& considerMate1,
	const bool& considerMate2) {
	
	vector<Alignment> newMate1Set;
	vector<Alignment> newMate2Set;
	// store the number of alignments
	unsigned int nMate1 = 0;
	unsigned int nMate2 = 0;

	Alignment bestMate1, bestMate2; 
	Alignment secondBestMate1, secondBestMate2;

	bool isMate1Aligned = !mate1Set.empty() && considerMate1;
	bool isMate2Aligned = !mate2Set.empty() && considerMate2;
	bool best           = false;
	bool secondBest     = false;

	unsigned int bestFl          = INT_MAX;
	unsigned int secondBestFl    = INT_MAX;
	
	if ( isMate1Aligned && isMate2Aligned ) {
		
		nMate1 = mate1Set.size();
		nMate2 = mate2Set.size();

		sort( mate1Set.begin(), mate1Set.end() );
		sort( mate2Set.begin(), mate2Set.end() );


		vector<Alignment>::iterator lastMinM2 = mate2Set.begin();
		for ( vector<Alignment>::iterator ite = mate1Set.begin(); ite != mate1Set.end(); ++ite ) {
			for ( vector<Alignment>::iterator ite2 = lastMinM2; ite2 != mate2Set.end(); ++ite2 ) {
				unsigned int length = ( ite->ReferenceBegin > ite2->ReferenceBegin) 
					? ite->ReferenceEnd - ite2->ReferenceBegin 
					: ite2->ReferenceEnd - ite->ReferenceBegin;
				
				if ( ( ite->ReferenceIndex == ite2->ReferenceIndex ) 
					&& ( length > ( 2 * expectedFragmentLength ) ) ) {
					if ( ite->ReferenceBegin > ite2->ReferenceBegin ) {
						lastMinM2 = ( ite->ReferenceBegin > ite2->ReferenceBegin) ? ite2 : lastMinM2;
						continue;
					} else {
						break;
					}
					
				// in the fragment length threshold
				} else {
					if ( IsBetterPair( *ite, *ite2, length, bestMate1, bestMate2, bestFl, expectedFragmentLength ) ) {
						// store the current best as second best
						if ( best ) {
							secondBest = true;
							secondBestMate1 = bestMate1;
							secondBestMate2 = bestMate2;
							secondBestFl    = bestFl;
						}
						best = true;
						bestMate1 = *ite;
						bestMate2 = *ite2;
						bestFl    = length;

					} else {
						if ( best && IsBetterPair( *ite, *ite2, length, secondBestMate1, secondBestMate2, secondBestFl, expectedFragmentLength ) ) {
							secondBest = true;
							secondBestMate1 = *ite;
							secondBestMate2 = *ite2;
							secondBestFl    = length;
						}
					}
				}
				
			}
		}

		
		if ( best ) {
			newMate1Set.push_back( bestMate1 );
			newMate2Set.push_back( bestMate2 );
		}
		else {
			// pick up mates having highest MQ
			sort ( mate1Set.begin(), mate1Set.end(), LessThanMQ );
			sort ( mate2Set.begin(), mate2Set.end(), LessThanMQ );
			newMate1Set.push_back( *mate1Set.rbegin() );
			newMate2Set.push_back( *mate2Set.rbegin() );
		}

		
		if ( secondBest ) {
			newMate1Set.begin()->NextBestQuality = secondBestMate1.Quality;
			newMate2Set.begin()->NextBestQuality = secondBestMate2.Quality;
		} else {
			if ( best ) {
				sort ( mate1Set.begin(), mate1Set.end(), LessThanMQ );
				sort ( mate2Set.begin(), mate2Set.end(), LessThanMQ );
			}

			if ( nMate1 == 1 )
				newMate1Set.begin()->NextBestQuality = 0;
			else
				newMate1Set.begin()->NextBestQuality = ( mate1Set.rbegin() + 1 )->Quality;

			if ( nMate2 == 1 )
				newMate2Set.begin()->NextBestQuality = 0;
			else
				newMate2Set.begin()->NextBestQuality = ( mate2Set.rbegin() + 1 )->Quality;
		}

		newMate1Set.begin()->NumMapped = nMate1;
		newMate2Set.begin()->NumMapped = nMate2;
		
		mate1Set.clear();
		mate2Set.clear();
		mate1Set = newMate1Set;
		mate2Set = newMate2Set;
		
		
	} else if ( isMate1Aligned ) {
		nMate1 = mate1Set.size();

		sort ( mate1Set.begin(), mate1Set.end(), LessThanMQ );
		// note: the size of mate1Set must be larger than one
		vector<Alignment>::reverse_iterator ite = mate1Set.rbegin();
		// the one having the highest MQ
		newMate1Set.push_back( *ite );
		ite++;
		newMate1Set.begin()->NextBestQuality = ite->Quality;
		newMate1Set.begin()->NumMapped = nMate1;

		mate1Set.clear();
		mate1Set = newMate1Set;

	} else if ( isMate2Aligned ) {
		nMate2 = mate2Set.size();
		
		sort ( mate2Set.begin(), mate2Set.end(), LessThanMQ );
		// note: the size of mate2Set must be larger than one
		vector<Alignment>::reverse_iterator ite = mate2Set.rbegin();
		// the one having the highest MQ
		newMate2Set.push_back( *ite );
		ite++;
		newMate2Set.begin()->NextBestQuality = ite->Quality;
		newMate2Set.begin()->NumMapped = nMate2;

		mate2Set.clear();
		mate2Set = newMate2Set;

	}
}


