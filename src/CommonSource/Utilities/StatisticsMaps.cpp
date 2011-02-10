#include "StatisticsMaps.h"

CStatisticsMaps::CStatisticsMaps( void )
	: _fragmentLength(0)
	, _localSearchRadius(0)
	, _allowedMismatch(0)
	, _setExpectedStatistics(false)
	, nFragment(10000)
	, nFrangmentOver(0)
	, nFrangmentUnder(0)
	, minFragment(-99)
	, nReadLength(1000)
	, nReadLengthOver(0)
	, nReadLengthUnder(0)
	, nMultiplicity(1000)
	, nMultiplicityOver(0)
	, nMultiplicityUnder(0)
	, nMappingQuality(100)
	, nMappingQualityOver(0)
	, nMappingQualityUnder(0)
	, nMismatch(100)
	, nMismatchOver(0)
	, nMismatchUnder(0)
	, non_unique(0)
	, non_multiple(0)
	, unique_unique(0)
	, unique_multiple(0)
	, multiple_multiple(0)
{
	fragments        = new uint64_t [ nFragment ];
	readLengths      = new uint64_t [ nReadLength ];
	multiplicities   = new uint64_t [ nMultiplicity ];
	mappingQualities = new uint64_t [ nMappingQuality ];
	mismatches       = new uint64_t [ nMismatch ];
	
	memset ( fragments, 0, nFragment * sizeof(uint64_t) );
	memset ( readLengths, 0, nReadLength * sizeof(uint64_t) );
	memset ( multiplicities, 0, nMultiplicity * sizeof(uint64_t) );
	memset ( mappingQualities, 0, nMappingQuality * sizeof(uint64_t) );
	memset ( mismatches, 0, nMismatch * sizeof(uint64_t) );
}

CStatisticsMaps::~CStatisticsMaps( void ) {
	if ( fragments )        delete [] fragments;
	if ( readLengths )      delete [] readLengths;
	if ( multiplicities )   delete [] multiplicities;
	if ( mappingQualities ) delete [] mappingQualities;
	if ( mismatches )       delete [] mismatches;

	fragments        = NULL;
	readLengths      = NULL;
	multiplicities   = NULL;
	mappingQualities = NULL;
	mismatches       = NULL;
}

void CStatisticsMaps::SetExpectedStatistics( const uint32_t fragmentLength, const uint32_t localSearchRadius, const float allowedMismatch ) {
	_fragmentLength        = fragmentLength;
	_localSearchRadius     = localSearchRadius;
	_allowedMismatch       = allowedMismatch;
	_setExpectedStatistics = true;
}

void CStatisticsMaps::Reset( void ) {
	
	nFrangmentOver       = 0;
	nFrangmentUnder      = 0;
	nReadLengthOver      = 0;
	nReadLengthUnder     = 0;
	nMultiplicityOver    = 0;
	nMultiplicityUnder   = 0;
	nMappingQualityOver  = 0;
	nMappingQualityUnder = 0;
	nMismatchOver        = 0;
	nMismatchUnder       = 0;
	non_unique           = 0;
	non_multiple         = 0;
	unique_unique        = 0;
	unique_multiple      = 0;
	multiple_multiple    = 0;

	memset ( fragments, 0, nFragment * sizeof(uint64_t) );
	memset ( readLengths, 0, nReadLength * sizeof(uint64_t) );
	memset ( multiplicities, 0, nMultiplicity * sizeof(uint64_t) );
	memset ( mappingQualities, 0, nMappingQuality * sizeof(uint64_t) );
	memset ( mismatches, 0, nMismatch * sizeof(uint64_t) );

}

inline void CStatisticsMaps::SaveFragment( const Alignment& al1, const Alignment& al2, const SequencingTechnologies& tech ) {
	// true: +; false: -.
	bool strand1 = !al1.IsReverseStrand;
	bool strand2 = !al2.IsReverseStrand;
	bool okay    = false;
	int64_t length = 0;
	switch( tech ) {
		case ST_454:
			if ( strand1 == strand2 ) {
				okay = true;
				length = strand1 ? al1.ReferenceEnd - al2.ReferenceBegin + 1 : al2.ReferenceEnd - al1.ReferenceBegin + 1;
			}
			break;
		case ST_SOLID:
			if ( strand1 == strand2 ) {
				okay = true;
				length = strand1 ? al2.ReferenceEnd - al1.ReferenceBegin + 1 : al1.ReferenceEnd - al2.ReferenceBegin + 1;
			}
			break;
		default:
			if ( strand1 != strand2 ) {
				okay = true;

				if ( strand1 && (al2.ReferenceEnd < al1.ReferenceBegin) )
					fprintf(stderr, "%s\n", al1.Query.CData());
					//cerr << al1.Query.CData() << endl;
				else if ( !strand1 && (al1.ReferenceEnd < al2.ReferenceBegin) )
					fprintf(stderr, "%s\n", al1.Query.CData());
					//cerr << al1.Query.CData() << endl;

				length = strand1 ? (int64_t)al2.ReferenceEnd - (int64_t)al1.ReferenceBegin + 1 : (int64_t)al1.ReferenceEnd - (int64_t)al2.ReferenceBegin + 1;
			}
	}

	if ( okay ) {
		if ( length < minFragment )
			nFrangmentUnder++;
		else if ( length > ( minFragment + (int64_t)nFragment ) )
			nFrangmentOver++;
		else
			fragments[ length - minFragment ]++;
	}
}

inline void CStatisticsMaps::SaveReadLength( const unsigned int length ) {
	
	if ( length > nReadLength ) {
		nReadLengthOver++;
	} else {
		readLengths[ length ]++;
	}
}

inline void CStatisticsMaps::SaveMultiplicity( const unsigned int nAlignment ) {
	
	if ( nAlignment > nMultiplicity ) {
		nMappingQualityOver++;
	} else {
		multiplicities[ nAlignment ]++;
	}
}

inline void CStatisticsMaps::SaveMappingQuality( const unsigned char mq ) {
	
	//for ( vector<Alignment>::iterator ite = als.begin(); ite !=als.end(); ++ite ) {
		//unsigned char mq = ite->Quality;
		// mq shouldn't larger than nMappingQuality = 100
		if ( mq > nMappingQuality ) {
			nMappingQualityOver++;
		} else {
			mappingQualities[ mq ]++;
		}
	//}
}

inline void CStatisticsMaps::SaveMismatch( const unsigned short mm ) {
	
	//for ( vector<Alignment>::iterator ite = als.begin(); ite !=als.end(); ++ite ) {
		//unsigned short mm = ite->NumMismatches;
		if ( mm > nMismatch ) {
			nMismatchOver++;
		} else {
			mismatches[ mm ]++;
		}
	//}
}

inline void CStatisticsMaps::PrintMap( 
	FILE *const fOut, 
	const char* title,
	const uint64_t& size, 
	const uint64_t& over, 
	const uint64_t& under, 
	const uint64_t *const array, 
	const int64_t& start  ) {
	
	fprintf( fOut, "\n%s\n", title );

	// calculate sum and count
	int64_t sum   = 0;
	int64_t count = 0;
	for ( uint64_t i = 0; i < size; ++i ) {
		sum += array[i] * ( start + i );
		count += array[i];
	}

	double mean = sum / (double) count;
	long double std = 0;
	for ( uint64_t i = 0; i < size; ++i ) {
		double temp1 = ( start + i ) - mean;
		double temp2 = pow( temp1, 2.0 );
		double temp3 = temp2 * array[i];

		std += temp3;
	}
	//double std = sqrt( (( pow(sum,2.0)) -(( 1.0/count) * (pow(sum,2.0))))/ (count -1.0));
	std = sqrt( ( std/count ) );
	fprintf( fOut, "\tTOT\tMEAN\tSTD\tIN\tOVER\tUNDER\n" );
	fprintf( fOut, "\t%lu\t%10.3f\t%10.3Lf\t%lu\t%lu\t%lu\n", count + over + under, mean, std, count, over, under);
	fprintf( fOut, "\tbin\tx\tn\tcum\n");

	uint64_t cum = under;
	for ( uint64_t i = 0; i < size; ++i ) {
		if ( array[i] == 0 )
			continue;

		cum += array[i];
		fprintf( fOut, "\t%lu\t%lu\t%lu\t%lu\n", i, i + start, array[i], cum );
	}
}

void CStatisticsMaps::PrintMaps( const char* filename, const char* readGroupId ) {
	FILE* fOut;
	fOut = fopen( filename, "w" );

	// print read group ID
	fprintf( fOut, "RG:%s\n", readGroupId );
	
	if ( fOut != NULL ) {
		char buffer[1024];
		uint8_t n = sprintf( buffer, "LF fragment mapping length (-mfl: %u; -ls: %u)", _fragmentLength, _localSearchRadius);
		if ( n > 1024 ) {
			printf("ERROR: The buffer for LF title is insufficient.\n");
			exit(1);
		}
		PrintMap( fOut, buffer,                nFragment,       nFrangmentOver,      nFrangmentUnder,      fragments,        minFragment );
		PrintMap( fOut, "LR read length",               nReadLength,     nReadLengthOver,     nReadLengthUnder,     readLengths,      0 );
		PrintMap( fOut, "NA read mapping multiplicity", nMultiplicity,   nMultiplicityOver,   nMultiplicityUnder,   multiplicities,   0 );
		PrintMap( fOut, "RQ read map quality",          nMappingQuality, nMappingQualityOver, nMappingQualityUnder, mappingQualities, 0 );
		n = sprintf( buffer, "MM read map mismatch (-mm/-mmp: %4.2f)", _allowedMismatch );
		if ( n > 1024 ) {
			printf("ERROR: The buffer for MM title is insufficient.\n");
			exit(1);
		}
		PrintMap( fOut, buffer,                nMismatch,       nMismatchOver,       nMismatchUnder,       mismatches,       0 );

		// print pair combinations
		uint64_t total = non_unique + non_multiple + unique_unique + unique_multiple + multiple_multiple;
		fprintf( fOut, "\nRC pair multiplicity combinations\n");
		fprintf( fOut, "\tTOT\tMEAN\tSTD\tIN\tOVER\tUNDER\n" );
		fprintf( fOut, "\t%lu\t0\t0\t%lu\t0\t0\n", total, total );
		fprintf( fOut, "\tbin\tcombo\tn\tcum\tlabel\n");
		fprintf( fOut, "\t1\t1\t%lu\t%lu\t0-1\n", non_unique,        non_unique );
		fprintf( fOut, "\t2\t2\t%lu\t%lu\t0-N\n", non_multiple,      non_unique + non_multiple );
		fprintf( fOut, "\t4\t4\t%lu\t%lu\t1-1\n", unique_unique,     non_unique + non_multiple + non_multiple );
		fprintf( fOut, "\t5\t5\t%lu\t%lu\t1-N\n", unique_multiple,   non_unique + non_multiple + non_multiple + unique_multiple );
		fprintf( fOut, "\t8\t8\t%lu\t%lu\tN-N\n", multiple_multiple, total );
		fclose( fOut );
	} else {
		printf("ERROR: The statistics maps cannot be printed out.\n");
		exit(1);
	}
}

void CStatisticsMaps::SaveRecord( 
	const Alignment& al1, 
	const Alignment& al2, 
	const bool isPairedEnd, 
	const SequencingTechnologies& tech ) {
	
	//bool isUU = isPairedEnd && ( als1.size() == 1 ) && ( als2.size() == 1 );

	// calculate fragment length
	//if ( isUU ) {
	//	
	//}

	char* qPtr;
	char* rPtr;
	//SaveReadLength( al1.BaseQualities.Length() );
	if ( al1.IsMapped ) {
		qPtr = (char*)al1.Query.CData();
		rPtr = (char*)al1.Reference.CData();
		uint32_t mappedLength = 0;
		for ( uint32_t i = 0; i < al1.Query.Length(); ++i ) {
			if ( ( *qPtr != '-' ) && ( *rPtr != 'Z' ) )
				mappedLength++;

			++qPtr;
			++rPtr;
		}

		SaveReadLength( mappedLength );
		SaveMultiplicity( al1.NumMapped );
		SaveMappingQuality( al1.Quality );
		SaveMismatch( al1.NumMismatches );
	}

	if ( isPairedEnd ) {
		//SaveReadLength( al2.BaseQualities.Length() );
		if ( al2.IsMapped ) {
			qPtr = (char*)al2.Query.CData();
			rPtr = (char*)al2.Reference.CData();
			uint32_t mappedLength = 0;
			for ( uint32_t i = 0; i < al2.Query.Length(); ++i ) {
				if ( ( *qPtr != '-' ) && ( *rPtr != 'Z' ) )
					mappedLength++;

				++qPtr;
				++rPtr;
			}

			SaveReadLength( mappedLength );
			SaveMultiplicity( al2.NumMapped );
			SaveMappingQuality( al2.Quality );
			SaveMismatch( al2.NumMismatches );
		}

		if ( ( ( al1.NumMapped == 0 ) && ( al2.NumMapped == 1 ) ) || ( ( al1.NumMapped == 1 ) && ( al2.NumMapped == 0 ) ) )
			non_unique++;

		if ( ( ( al1.NumMapped == 0 ) && ( al2.NumMapped > 1 ) ) || ( ( al1.NumMapped > 1 ) && ( al2.NumMapped == 0 ) ) )
			non_multiple++;

		if ( ( al1.NumMapped == 1 ) && ( al2.NumMapped == 1 ) ) {
			unique_unique++;
			SaveFragment( al1, al2, tech );
		}

		if ( ( ( al1.NumMapped == 1 ) && ( al2.NumMapped > 1 ) ) || ( ( al1.NumMapped > 1 ) && ( al2.NumMapped == 1 ) ) )
			unique_multiple++;

		if ( ( al1.NumMapped > 0 ) && ( al2.NumMapped > 1 ) )
			multiple_multiple++;
	}
}
