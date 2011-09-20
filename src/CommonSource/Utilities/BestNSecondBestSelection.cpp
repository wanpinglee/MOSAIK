#include "BestNSecondBestSelection.h"


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
	bool competitor_model = isProperOrientation ( 
		competitor_mate1.IsReverseStrand, 
		competitor_mate2.IsReverseStrand, 
		competitor_mate1.ReferenceBegin, 
		competitor_mate2.ReferenceBegin,
		true,
		tech);
	
	bool current_model    = isProperOrientation (
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

		//float competitor_fragScore = ( expectedFragmentLength - competitor_diff ) / (float) ( expectedFragmentLength );
		//float fragScore            = ( expectedFragmentLength - diff ) / (float) ( expectedFragmentLength );
		float competitor_fragScore = 0;
		float fragScore = 0;

		float competitor_mqScore   = ( competitor_mate1.Quality + competitor_mate2.Quality ) / 200.0;
		float mqScore              = ( mate1.Quality + mate2.Quality ) / 200.0;

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
	vector<Alignment*>& mate1Set, 
	vector<Alignment*>& mate2Set, 
	const unsigned int expectedFragmentLength,
	const SequencingTechnologies& tech,
	const unsigned int numMate1Bases,
	const unsigned int numMate2Bases,
	const bool& considerMate1,
	const bool& considerMate2) {
	//const unsigned int highestSwScoreMate1,
	//const unsigned int highestSwScoreMate2) {
	
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
						if ( IsBetterPair( **ite, **ite2, length, bestMate1, bestMate2, bestFl, expectedFragmentLength, tech, numMate1Bases, numMate2Bases ) ) {
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
							if ( best && IsBetterPair( **ite, **ite2, length, secondBestMate1, secondBestMate2, secondBestFl, expectedFragmentLength, tech, numMate1Bases, numMate2Bases ) ) {
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

		
		if ( best ) {
			newMate1Set.push_back( bestMate1 );
			newMate2Set.push_back( bestMate2 );
		}
		else {
			// pick up mates having highest MQ
			sort ( mate1Set.begin(), mate1Set.end(), Alignment_LessThanMq() );
			sort ( mate2Set.begin(), mate2Set.end(), Alignment_LessThanMq() );
			newMate1Set.push_back( **(mate1Set.rbegin()) );
			newMate2Set.push_back( **(mate2Set.rbegin()) );
		}

		
		if ( secondBest ) {
			newMate1Set.begin()->NextBestQuality = ( bestMate1 == secondBestMate1 ) ? 0 : secondBestMate1.Quality;
			newMate1Set.begin()->NextSwScore     = ( bestMate1 == secondBestMate1 ) ? 0 : secondBestMate1.SwScore;
			newMate2Set.begin()->NextBestQuality = ( bestMate2 == secondBestMate2 ) ? 0 : secondBestMate2.Quality;
			newMate2Set.begin()->NextSwScore     = ( bestMate2 == secondBestMate2 ) ? 0 : secondBestMate2.SwScore;
		} else {
			if ( best ) {
				sort ( mate1Set.begin(), mate1Set.end(), Alignment_LessThanMq() );
				sort ( mate2Set.begin(), mate2Set.end(), Alignment_LessThanMq() );
			}

			if ( nMate1 == 1 )
				newMate1Set.begin()->NextBestQuality = 0;
			else
				newMate1Set.begin()->NextBestQuality = (*(mate1Set.rbegin() + 1))->Quality;

			if ( nMate2 == 1 )
				newMate2Set.begin()->NextBestQuality = 0;
			else
				newMate2Set.begin()->NextBestQuality = (*(mate2Set.rbegin() + 1))->Quality;
		}

		newMate1Set.begin()->NumMapped = nMate1;
		newMate2Set.begin()->NumMapped = nMate2;
		
		//mate1Set.clear();
		//mate2Set.clear();
		for (unsigned int i = 0; i < newMate1Set.size(); ++i)
		  *mate1Set[i] = newMate1Set[i];
		
		for (unsigned int i = 0; i < newMate2Set.size(); ++i)
		  *mate2Set[i] = newMate2Set[i];
		
		
	} else if ( isMate1Aligned ) {
		nMate1 = mate1Set.size();

		sort ( mate1Set.begin(), mate1Set.end(), Alignment_LessThanMq() );
		// note: the size of mate1Set must be larger than one
		vector<Alignment*>::reverse_iterator ite = mate1Set.rbegin();
		// the one having the highest MQ
		newMate1Set.push_back( **ite );
		ite++;
		newMate1Set.begin()->NextBestQuality = (*ite)->Quality;
		newMate1Set.begin()->NumMapped = nMate1;

		//mate1Set.clear();
		//mate1Set = newMate1Set;
		for (unsigned int i = 0; i < newMate1Set.size(); ++i)
		  *mate1Set[i] = newMate1Set[i];

	} else if ( isMate2Aligned ) {
		nMate2 = mate2Set.size();
		
		sort ( mate2Set.begin(), mate2Set.end(), Alignment_LessThanMq() );
		// note: the size of mate2Set must be larger than one
		vector<Alignment*>::reverse_iterator ite = mate2Set.rbegin();
		// the one having the highest MQ
		newMate2Set.push_back( **ite );
		ite++;
		newMate2Set.begin()->NextBestQuality = (*ite)->Quality;
		newMate2Set.begin()->NumMapped = nMate2;

		//mate2Set.clear();
		//mate2Set = newMate2Set;
		for (unsigned int i = 0; i < newMate2Set.size(); ++i)
		  *mate2Set[i] = newMate2Set[i];
	}
}


