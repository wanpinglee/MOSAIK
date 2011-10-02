#include "BestNSecondBestSelection.h"
#include "StrandChecker.h"

// compare the given proper pairs
inline bool BestNSecondBestSelection::IsBetterPair ( 
	const Alignment& competitor_mate1, 
	const Alignment& competitor_mate2, 
	const unsigned int competitor_fragmentLength, 
	const Alignment& mate1, 
	const Alignment& mate2,
	const unsigned int fragmentLength,
	const unsigned int expectedFragmentLength,
	const SequencingTechnologies& tech,
	const unsigned int numMate1Bases, 
	const unsigned int numMate2Bases) {

	// rescured mate always wins
	//bool competitor_wasRescued = competitor_mate1.WasRescued || competitor_mate2.WasRescued;
	//bool wasRescued            = mate1.WasRescued || mate2.WasRescued;
	//if ( competitor_wasRescued && !wasRescued ) return true;
	//if ( !competitor_wasRescued && wasRescued ) return false;

	//if ( competitor_mate1.WasRescued ) return true;
	//if ( competitor_mate2.WasRescued ) return true;
	//if ( mate1.WasRescued ) return false;
	//if ( mate2.WasRescued ) return false;

	// proper pair always wins improper pair
	//bool competitor_model = ( competitor_mate1.IsReverseStrand != competitor_mate2.IsReverseStrand ) ? true : false;
	//bool current_model    = ( mate1.IsReverseStrand != mate2.IsReverseStrand ) ? true : false;
	bool competitor_model = StrandChecker::isProperOrientation ( 
		competitor_mate1.IsReverseStrand, 
		competitor_mate2.IsReverseStrand, 
		competitor_mate1.ReferenceBegin, 
		competitor_mate2.ReferenceBegin,
		true,
		tech);
	
	bool current_model    = StrandChecker::isProperOrientation (
		mate1.IsReverseStrand,
		mate2.IsReverseStrand,
		mate1.ReferenceBegin,
		mate2.ReferenceBegin,
		true,
		tech);
	
	if ( competitor_model && !current_model ) return true;
	if ( !competitor_model && current_model ) return false;

	//double competitor = ( competitor_mate1.Quality + competitor_mate2.Quality ) / 100.00;
	//double current    = ( mate1.Quality + mate2.Quality ) / 100.00;

	unsigned int competitor_diff = ( expectedFragmentLength > competitor_fragmentLength ) ? 
		expectedFragmentLength - competitor_fragmentLength : 
		competitor_fragmentLength - expectedFragmentLength;
	
	unsigned int diff            = ( expectedFragmentLength > fragmentLength ) ? 
		expectedFragmentLength - fragmentLength : 
		fragmentLength - expectedFragmentLength;

	//competitor += ( 200 - (int)competitor_diff ) / 200;
	//current    += ( 200 - (int)diff ) / 200;
	
	//if ( competitor > current ) return true;
	//else return false;
	
	unsigned int diff_diff = ( competitor_diff < diff ) ? diff - competitor_diff : competitor_diff - diff;

	if ( diff_diff < ( expectedFragmentLength / 2 ) ) {
	//if ( competitor_diff < expectedFragmentLength ) {
		float competitor_swScore = ( competitor_mate1.SwScore + competitor_mate2.SwScore ) / (float)( ( numMate1Bases + numMate2Bases ) * 10 );
		float swScore            = ( mate1.SwScore + mate2.SwScore ) / (float)( ( numMate1Bases + numMate2Bases ) * 10 );

		competitor_swScore *= 1.5;
		swScore *= 1.5;

		float competitor_fragScore = ( expectedFragmentLength - competitor_diff ) / (float) ( expectedFragmentLength );
		float fragScore            = ( expectedFragmentLength - diff ) / (float) ( expectedFragmentLength );
		//float competitor_fragScore = 0;
		//float fragScore = 0;

		//float competitor_mqScore   = ( competitor_mate1.Quality + competitor_mate2.Quality ) / 200.0;
		//float mqScore              = ( mate1.Quality + mate2.Quality ) / 200.0;
		float competitor_mqScore   = 0;
		float mqScore              = 0;

		float competitor_finalScore = competitor_swScore + competitor_fragScore + competitor_mqScore;
		float finalScore            = swScore + fragScore + mqScore;

		return competitor_finalScore >= finalScore;

	} else {
		if ( competitor_diff < diff ) return true;
		else return false;
	}
	//return competitor_diff < diff;
}

void BestNSecondBestSelection::CalculateFragmentLength( const Alignment& al1, const Alignment& al2, const SequencingTechnologies& tech, unsigned int& length ) {
	
	bool strand1 = !al1.IsReverseStrand;
	//bool strand2 = !al2.IsReverseStrand;
	
	switch( tech ) {
		case ST_454:
			//if ( strand1 == strand2 ) {
			//	if ( al1.ReferenceIndex != al2.ReferenceIndex )
			//		okay = false;
			//	else {
			//		okay = true;
					length = strand1 ? (int64_t)al1.ReferenceEnd - (int64_t)al2.ReferenceBegin + 1 : (int64_t)al2.ReferenceEnd - (int64_t)al1.ReferenceBegin + 1;
			//	}
			//}
			break;
		case ST_SOLID:
			//if ( strand1 == strand2 ) {
			//	if ( al1.ReferenceIndex != al2.ReferenceIndex )
			//		okay = false;
			//	else {
			//		okay = true;
					length = strand1 ? (int64_t)al2.ReferenceEnd - (int64_t)al1.ReferenceBegin + 1 : (int64_t)al1.ReferenceEnd - (int64_t)al2.ReferenceBegin + 1;
			//	}
			//}
			break;

		case ST_ILLUMINA_LONG:
			//if ( strand1 != strand2 ) {
			//	if ( al1.ReferenceIndex != al2.ReferenceIndex )
			//		okay = false;
			//	else {
			//		okay = true;
					length = strand1 ? (int64_t)al1.ReferenceEnd - (int64_t)al2.ReferenceBegin + 1 : (int64_t)al2.ReferenceEnd - (int64_t)al1.ReferenceBegin + 1;
			//	}
			//}
		break;
		default:
			//if ( strand1 != strand2 ) {
			//	if ( al1.ReferenceIndex != al2.ReferenceIndex )
			//		okay = false;
			//	else {
			//		okay = true;
					length = strand1 ? (int64_t)al2.ReferenceEnd - (int64_t)al1.ReferenceBegin + 1 : (int64_t)al1.ReferenceEnd - (int64_t)al2.ReferenceBegin + 1;
			//	}
			//}
	}

}

// Select and only keep best and 2nd best
void BestNSecondBestSelection::Select ( 
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
	const bool& resetNextSwScore) {
	
	//vector<Alignment> newMate1Set;
	//vector<Alignment> newMate2Set;
	// store the number of alignments
	unsigned int nMate1 = 0;
	unsigned int nMate2 = 0;

	//Alignment bestMate1, bestMate2; 
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

		// sort by positions
		sort( mate1Set.begin(), mate1Set.end(), Alignment_LessThanPosition() );
		sort( mate2Set.begin(), mate2Set.end(), Alignment_LessThanPosition() );


		vector<Alignment*>::iterator lastMinM2 = mate2Set.begin();
		bestMate1 = **mate1Set.begin();
		bestMate2 = **mate2Set.begin();

		for ( vector<Alignment*>::iterator ite = mate1Set.begin(); ite != mate1Set.end(); ++ite ) {
			//if ( ( ite->SwScore / (float) highestSwScoreMate1 ) < 0.9 ) continue;

			for ( vector<Alignment*>::iterator ite2 = lastMinM2; ite2 != mate2Set.end(); ++ite2 ) {
				//if ( ( ite2->SwScore / (float) highestSwScoreMate2 ) < 0.9 ) continue;

				// fragment length
				unsigned int length = 0;
				CalculateFragmentLength( **ite, **ite2, tech, length );
				

				if ( (*ite)->ReferenceIndex == (*ite2)->ReferenceIndex ) {
					if ( length > ( 2 * expectedFragmentLength ) ) {
						if ( (*ite)->ReferenceBegin > (*ite2)->ReferenceBegin ) {
							lastMinM2 = ite2;
							continue;
						} else {
							break;
						}
					
				// in the fragment length threshold
					} else {
						if ( IsBetterPair(**ite, **ite2, length, bestMate1, bestMate2, 
						                  bestFl, expectedFragmentLength, tech, numMate1Bases, 
								  numMate2Bases ) ) {
							// store the current best as second best
							if ( best ) {
								secondBest = true;
								secondBestMate1 = bestMate1;
								secondBestMate2 = bestMate2;
								secondBestFl    = bestFl;
							}
							best = true;
							bestMate1 = **ite;
							bestMate2 = **ite2;
							bestFl    = length;
	
						} else {
							if (best && IsBetterPair(**ite, **ite2, length, secondBestMate1, 
							                         secondBestMate2, secondBestFl, expectedFragmentLength, 
										 tech, numMate1Bases, numMate2Bases ) ) {
								secondBest = true;
								secondBestMate1 = **ite;
								secondBestMate2 = **ite2;
								secondBestFl    = length;
							}
						}
					}

				// located at different chromosomes
				} else {
					if ( (*ite)->ReferenceIndex > (*ite2)->ReferenceIndex ) {
						lastMinM2 = ite2;
						continue;
					} else {
						break;
					}
				}
			}
		}

		
		if (!best) {
			// pick up mates having highest MQ
			sort ( mate1Set.begin(), mate1Set.end(), Alignment_LessThanMq() );
			sort ( mate2Set.begin(), mate2Set.end(), Alignment_LessThanMq() );
			bestMate1 = **(mate1Set.rbegin());
			bestMate2 = **(mate2Set.rbegin());
		}

		
		if (secondBest) {
			bestMate1.NextBestQuality = ( bestMate1 == secondBestMate1 ) ? 0 : secondBestMate1.Quality;
			bestMate2.NextBestQuality = ( bestMate2 == secondBestMate2 ) ? 0 : secondBestMate2.Quality;
			if (resetNextSwScore) {
				bestMate1.NextSwScore = ( bestMate1 == secondBestMate1 ) ? 0 : secondBestMate1.SwScore;
				bestMate2.NextSwScore = ( bestMate2 == secondBestMate2 ) ? 0 : secondBestMate2.SwScore;
			}
		} else {  // !secondBest
			if (best) {
				sort (mate1Set.begin(), mate1Set.end(), Alignment_LessThanMq());
				sort (mate2Set.begin(), mate2Set.end(), Alignment_LessThanMq());
			}

			if (nMate1 == 1) {
				bestMate1.NextBestQuality = 0;
				if (resetNextSwScore) bestMate1.NextSwScore = 0;
			} else {
				bestMate1.NextBestQuality = (*(mate1Set.rbegin() + 1))->Quality;
				if (resetNextSwScore) bestMate1.NextSwScore = (*(mate1Set.rbegin() + 1))->SwScore;
			}

			if (nMate2 == 1) {
				bestMate2.NextBestQuality = 0;
				if (resetNextSwScore) bestMate2.NextSwScore = 0;
			} else {
				bestMate2.NextBestQuality = (*(mate2Set.rbegin() + 1))->Quality;
				if (resetNextSwScore) bestMate2.NextSwScore = (*(mate2Set.rbegin() + 1))->SwScore;
			}
		}

		bestMate1.NumMapped = nMate1;
		bestMate2.NumMapped = nMate2;
	
	} else if ( isMate1Aligned ) {
		nMate1 = mate1Set.size();
		sort ( mate1Set.begin(), mate1Set.end(), Alignment_LessThanMq() );

		// note: the size of mate1Set must be larger than one
		vector<Alignment*>::reverse_iterator ite = mate1Set.rbegin();
		// the one having the highest MQ
		bestMate1 = **ite;
		ite++;
		bestMate1.NextBestQuality = (*ite)->Quality;
		bestMate1.NumMapped = nMate1;
		if (resetNextSwScore)
			bestMate1.NextSwScore = (*ite)->SwScore;

	} else if ( isMate2Aligned ) {
		nMate2 = mate2Set.size();
		
		sort ( mate2Set.begin(), mate2Set.end(), Alignment_LessThanMq() );
		// note: the size of mate2Set must be larger than one
		vector<Alignment*>::reverse_iterator ite = mate2Set.rbegin();
		// the one having the highest MQ
		bestMate2 = **ite;
		ite++;
		bestMate2.NextBestQuality = (*ite)->Quality;
		bestMate2.NumMapped = nMate2;
		if (resetNextSwScore)
			bestMate2.NextSwScore = (*ite)->SwScore;
	} // end if-else
}


