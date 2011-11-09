#include "StatisticsMaps.h"

CStatisticsMaps::CStatisticsMaps( void )
	: _fragmentLength(0)
	, _localSearchRadius(0)
	, _allowedMismatch(0.0)
	, _setExpectedStatistics(false)
	, nFragment(10000)
	, nFrangmentOver(0)
	, nFrangmentUnder(0)
	, minFragment(-99)
	//, nIsize(10000)
	//, nIsizeOver(0)
	//, nIsizeUnder(0)
	, nReadLength(1000)
	, nReadLengthOver(0)
	, nReadLengthUnder(0)
	, nMultiplicity(5000)
	, nMultiplicityOver(0)
	, nMultiplicityUnder(0)
	, nMappingQuality(256)
	, nMappingQualityOver(0)
	, nMappingQualityUnder(0)
	, minMappingQuality(0)
	, nMismatch(100)
	, nMismatchOver(0)
	, nMismatchUnder(0)
	, U(0)
	, M(0)
	, X(0)
	, F(0)
	, UU(0)
	, UM(0)
	, MM(0)
	, UF(0)
	, MF(0)
	, UX(0)
	, MX(0)
	, FF(0)
	, FX(0)
	, XX(0)
	, UU_localRescue(0)
	, UU_localConsistance(0)
	, UM_localRescue(0)
	, UM_localConsistance(0)
	, MM_localRescue(0)
	, MM_localConsistance(0)
	, f1_f2(0)
	, f1_r2(0)
	, r1_f2(0)
	, r1_r2(0)
	, f2_f1(0)
	, f2_r1(0)
	, r2_f1(0)
	, r2_r1(0)
{
	fragments        = new uint64_t [ nFragment ];
	//isizes           = new uint64_t [ nIsize ];
	readLengths      = new uint64_t [ nReadLength ];
	multiplicities   = new uint64_t [ nMultiplicity ];
	mappingQualities = new uint64_t [ nMappingQuality ];
	mismatches       = new uint64_t [ nMismatch ];
	
	memset ( fragments,        0, nFragment       * sizeof(uint64_t) );
	//memset ( isizes,           0, nIsize          * sizeof(uint64_t) );
	memset ( readLengths,      0, nReadLength     * sizeof(uint64_t) );
	memset ( multiplicities,   0, nMultiplicity   * sizeof(uint64_t) );
	memset ( mappingQualities, 0, nMappingQuality * sizeof(uint64_t) );
	memset ( mismatches,       0, nMismatch       * sizeof(uint64_t) );
}

CStatisticsMaps::~CStatisticsMaps( void ) {
	if ( fragments )        delete [] fragments;
	//if ( isizes )           delete [] isizes;
	if ( readLengths )      delete [] readLengths;
	if ( multiplicities )   delete [] multiplicities;
	if ( mappingQualities ) delete [] mappingQualities;
	if ( mismatches )       delete [] mismatches;

	fragments        = NULL;
	//isizes           = NULL;
	readLengths      = NULL;
	multiplicities   = NULL;
	mappingQualities = NULL;
	mismatches       = NULL;
}

// Note: once starting to calculate the map, don't use this function to change minFragment
void CStatisticsMaps::SetLfMin( int64_t flMin ) {
	minFragment = flMin;
}

void CStatisticsMaps::SetMqMin( uint8_t mqMin ) {
	minMappingQuality = mqMin;
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
	//nIsizeOver           = 0;
	//nIsizeUnder          = 0;
	nReadLengthOver      = 0;
	nReadLengthUnder     = 0;
	nMultiplicityOver    = 0;
	nMultiplicityUnder   = 0;
	nMappingQualityOver  = 0;
	nMappingQualityUnder = 0;
	nMismatchOver        = 0;
	nMismatchUnder       = 0;
	U                    = 0;
	M                    = 0;
	X                    = 0;
	F                    = 0;
	UU                   = 0;
	UM                   = 0;
	MM                   = 0;
	UF                   = 0;
	MF                   = 0;
	UX                   = 0;
	MX                   = 0;
	FF                   = 0;
	FX                   = 0;
	XX                   = 0;
	UU_localRescue       = 0;
	UU_localConsistance  = 0;
	UM_localRescue       = 0;
	UM_localConsistance  = 0;
	MM_localRescue       = 0;
	MM_localConsistance  = 0;
	f1_f2                = 0;
	f1_r2                = 0;
	r1_f2                = 0;
	r1_r2                = 0;
	f2_f1                = 0;
	f2_r1                = 0;
	r2_f1                = 0;
	r2_r1                = 0;

	memset ( fragments,        0, nFragment * sizeof(uint64_t) );
	//memset ( isizes,           0, nIsize * sizeof(uint64_t) );
	memset ( readLengths,      0, nReadLength * sizeof(uint64_t) );
	memset ( multiplicities,   0, nMultiplicity * sizeof(uint64_t) );
	memset ( mappingQualities, 0, nMappingQuality * sizeof(uint64_t) );
	memset ( mismatches,       0, nMismatch * sizeof(uint64_t) );

}

inline void CStatisticsMaps::SaveModel( const Alignment& al1, const Alignment& al2 ) {

	if ( al1.ReferenceIndex != al2.ReferenceIndex )
		return;

	if ( !al1.IsReverseStrand && !al2.IsReverseStrand && ( al1.ReferenceBegin <= al2.ReferenceBegin ) )
		f1_f2++;
	
	if ( !al1.IsReverseStrand && al2.IsReverseStrand  && ( al1.ReferenceBegin <= al2.ReferenceBegin ) )
		f1_r2++;

	if ( al1.IsReverseStrand  && !al2.IsReverseStrand && ( al1.ReferenceBegin <= al2.ReferenceBegin ) )
		r1_f2++;
	
	if ( al1.IsReverseStrand  && al2.IsReverseStrand  && ( al1.ReferenceBegin <= al2.ReferenceBegin ) )
		r1_r2++;
	
	if ( !al1.IsReverseStrand && !al2.IsReverseStrand && ( al1.ReferenceBegin > al2.ReferenceBegin ) )
		f2_f1++;

	if ( al1.IsReverseStrand  && !al2.IsReverseStrand && ( al1.ReferenceBegin > al2.ReferenceBegin ) )
		f2_r1++;

	if ( !al1.IsReverseStrand && al2.IsReverseStrand  && ( al1.ReferenceBegin > al2.ReferenceBegin ) )
		r2_f1++;

	if ( al1.IsReverseStrand  && al2.IsReverseStrand  && ( al1.ReferenceBegin > al2.ReferenceBegin ) )
		r2_r1++;
}

void CStatisticsMaps::SaveFragment( const Alignment& al1, const Alignment& al2, const SequencingTechnologies& tech ) {
	// true: +; false: -.
	bool strand1 = !al1.IsReverseStrand;
	bool strand2 = !al2.IsReverseStrand;
	bool okay    = false;
	int64_t length = 0;
	switch( tech ) {
		case ST_454:
			if ( strand1 == strand2 ) {
				if ( al1.ReferenceIndex != al2.ReferenceIndex )
					okay = false;
				else {
					okay = true;
					length = strand1 ? (int64_t)al1.ReferenceEnd - (int64_t)al2.ReferenceBegin + 1 : (int64_t)al2.ReferenceEnd - (int64_t)al1.ReferenceBegin + 1;
				}
			}
			break;
		case ST_SOLID:
			if ( strand1 == strand2 ) {
				if ( al1.ReferenceIndex != al2.ReferenceIndex )
					okay = false;
				else {
					okay = true;
					length = strand1 ? (int64_t)al2.ReferenceEnd - (int64_t)al1.ReferenceBegin + 1 : (int64_t)al1.ReferenceEnd - (int64_t)al2.ReferenceBegin + 1;
				}
			}
			break;

		case ST_ILLUMINA_LONG:
			if ( strand1 != strand2 ) {
				if ( al1.ReferenceIndex != al2.ReferenceIndex )
					okay = false;
				else {
					okay = true;
					length = strand1 ? (int64_t)al1.ReferenceEnd - (int64_t)al2.ReferenceBegin + 1 : (int64_t)al2.ReferenceEnd - (int64_t)al1.ReferenceBegin + 1;
				}
			}
		break;
		default:
			if ( strand1 != strand2 ) {
				if ( al1.ReferenceIndex != al2.ReferenceIndex )
					okay = false;
				else {
					okay = true;
					length = strand1 ? (int64_t)al2.ReferenceEnd - (int64_t)al1.ReferenceBegin + 1 : (int64_t)al1.ReferenceEnd - (int64_t)al2.ReferenceBegin + 1;
				}
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

//inline void CStatisticsMaps::SaveIsize( const Alignment& al1, const Alignment& al2 ) {
//	unsigned int al1_5Prime = al1.IsReverseStrand ? al1.ReferenceEnd : al1.ReferenceBegin;
//	unsigned int al2_5Prime = al2.IsReverseStrand ? al2.ReferenceEnd : al2.ReferenceBegin;
//
//	unsigned int isize = ( al1_5Prime < al2_5Prime ) ? al2_5Prime - al1_5Prime : al1_5Prime - al2_5Prime;
//	if ( isize > nIsize ) {
//		nIsizeOver++;
//	} else {
//		isizes[ isize ]++;
//	}
//
//}

inline void CStatisticsMaps::SaveReadLength( const unsigned int length ) {
	
	if ( length > nReadLength ) {
		nReadLengthOver++;
	} else {
		readLengths[ length ]++;
	}
}

inline void CStatisticsMaps::SaveMultiplicity( const unsigned int nAlignment ) {
	
	if ( nAlignment > nMultiplicity ) {
		nMultiplicityOver++;
	} else {
		multiplicities[ nAlignment ]++;
	}
}

void CStatisticsMaps::SavePairMultiplicity( const unsigned int nMate1Alignments, const unsigned int nMate2Alignments, const bool mate1Rescued, const bool mate2Rescued, const bool mate1FilteredOut, const bool mate2FilteredOut, const bool isProperPair ) {
	
	// MM pairs
	if ( nMate1Alignments > 1 && nMate2Alignments > 1 ) {
		MM++;
		if ( mate1Rescued || mate2Rescued )
			MM_localRescue++;
		else if ( isProperPair )
			MM_localConsistance++;
	}
	
	// UU pairs
	if ( nMate1Alignments == 1 && nMate2Alignments == 1 ) {
		UU++;
		if ( mate1Rescued || mate2Rescued )
			UU_localRescue++;
		else if ( isProperPair )
			UU_localConsistance++;
	}

	// UM pairs
	if ( ( nMate1Alignments == 1 && nMate2Alignments > 1 ) || ( nMate1Alignments > 1 && nMate2Alignments == 1) ) {
		UM++;
		if ( mate1Rescued || mate2Rescued )
			UM_localRescue++;
		else if ( isProperPair )
			UM_localConsistance++;
	}

	// UF pairs
	if ( ( nMate1Alignments == 1 && ( ( nMate2Alignments == 0 ) && mate2FilteredOut ) ) || ( ( ( nMate1Alignments == 0 ) && mate1FilteredOut ) && nMate2Alignments == 1 ) ) 
		UF++;
	
	// MF pairs
	if ( ( nMate1Alignments > 1 && ( ( nMate2Alignments == 0 ) && mate2FilteredOut ) ) || ( ( ( nMate1Alignments == 0 ) && mate1FilteredOut ) && nMate2Alignments > 1 ) )
		MF++;
	
	// UX pairs
	if ( ( nMate1Alignments == 1 && ( ( nMate2Alignments == 0 ) && !mate2FilteredOut ) ) || ( ( ( nMate1Alignments == 0 ) && !mate1FilteredOut ) && nMate2Alignments == 1 ) )
		UX++;

	// MX pairs
	if ( ( nMate1Alignments > 1 && ( ( nMate2Alignments == 0 ) && !mate2FilteredOut ) ) || ( ( ( nMate1Alignments == 0 ) && !mate1FilteredOut ) && nMate2Alignments > 1 ) )
		MX++;

	// FF pairs
	if ( ( ( nMate1Alignments == 0 ) && mate1FilteredOut ) && ( ( nMate2Alignments == 0 ) && mate2FilteredOut ) )
		FF++;

	// FX pairs
	if ( ( ( ( nMate1Alignments == 0 ) && !mate1FilteredOut ) && ( ( nMate2Alignments == 0 ) && mate2FilteredOut ) )
		|| ( ( ( nMate1Alignments == 0 ) && mate1FilteredOut ) && ( ( nMate2Alignments == 0 ) && !mate2FilteredOut ) ) ) 
		FX++;

	// XX pairs
	if ( ( ( nMate1Alignments == 0 ) && !mate1FilteredOut ) && ( ( nMate2Alignments == 0 ) && !mate2FilteredOut ) )
		XX++;

}

inline void CStatisticsMaps::SaveMateMultiplicity( const unsigned int nAlignment, const bool isFilteredOut ) {
	
	if ( nAlignment == 0 ) {
		if ( isFilteredOut )
			F++;
		else
			X++;
	} else {
		if ( nAlignment == 1 )
			U++;
		else
			M++;
	}
}

inline void CStatisticsMaps::SaveMappingQuality( const unsigned char mq ) {
	
	//for ( vector<Alignment>::iterator ite = als.begin(); ite !=als.end(); ++ite ) {
		//unsigned char mq = ite->Quality;
		// mq shouldn't larger than nMappingQuality = 256
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

void CStatisticsMaps::PrintMap( 
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
	uint64_t count = 0;
	for ( uint64_t i = 0; i < size; ++i ) {
		sum += array[i] * ( start + i );
		count += array[i];
	}

	double mean = sum / (double) count;
	long double std = 0;
	for ( uint64_t i = 0; i < size; ++i ) {
		double temp1 = ( start + (int64_t)i ) - mean;
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
		fprintf( fOut, "\t%lu\t%ld\t%lu\t%lu\n", i, i + start, array[i], cum );
	}
}

void CStatisticsMaps::PrintMaps( const char* filename, const vector<MosaikReadFormat::ReadGroup>& readGroup, const unsigned char statMappingQuality ) {
	FILE* fOut;
	fOut = fopen( filename, "w" );

	// print the header
	for ( vector<MosaikReadFormat::ReadGroup>::const_iterator ite = readGroup.begin(); ite != readGroup.end(); ++ite ) {
		fprintf( fOut, "@RG\tID:%s\tSM:%s", ( ite->ReadGroupID.empty() ? "unknown" : ite->ReadGroupID.c_str() ), ( ite->SampleName.empty() ? "unknown" : ite->SampleName.c_str() ) );
		if ( !ite->LibraryName.empty() )  fprintf( fOut, "\tLB:%s", ite->LibraryName.c_str() );
		if ( !ite->Description.empty() )  fprintf( fOut, "\tDS:%s", ite->Description.c_str() );
		if ( !ite->PlatformUnit.empty() ) fprintf( fOut, "\tPU:%s", ite->PlatformUnit.c_str() );
		                                  fprintf( fOut, "\tPI:%u", ite->MedianFragmentLength );
		if ( !ite->CenterName.empty() )   fprintf( fOut, "\tCN:%s", ite->CenterName.c_str() );
		
		switch( ite->SequencingTechnology ) {
			case ST_454:
				fprintf( fOut, "\tPL:454" );
				break;
			case ST_HELICOS:
				fprintf( fOut, "\tPL:helicos" );
				break;
			case ST_ILLUMINA:
				fprintf( fOut, "\tPL:illumina" );
				break;
			case ST_ILLUMINA_LONG:
				fprintf( fOut, "\tPL:illumina long" );
				break;
			case ST_PACIFIC_BIOSCIENCES:
				fprintf( fOut, "\tPL:pacific biosciences" );
				break;
			case ST_SOLID:
				fprintf( fOut, "\tPL:solid" );
				break;
			case ST_SANGER:
				fprintf( fOut, "\tPL:sanger" );
				break;
			default:
				fprintf( fOut, "\tPL:unknown" );
		}
		fprintf( fOut, "\n" );
	}

	fprintf( fOut, "Mapping quality threshold:%u\n", statMappingQuality );
	
	if ( fOut != NULL ) {
		//PrintMap( fOut, "IS isize", nIsize, nIsizeOver, nIsizeUnder, isizes, 0);

		char buffer[1024];
		int n = sprintf( buffer, "LF fragment mapping length (-mfl: %u; -ls: %u)", _fragmentLength, _localSearchRadius);
		if ( n > 1024 ) {
			printf("ERROR: The buffer for LF title is insufficient.\n");
			exit(1);
		}
		PrintMap( fOut, buffer,                nFragment,       nFrangmentOver,      nFrangmentUnder,      fragments,        minFragment );
		
		PrintMap( fOut, "LR read mapping length",               nReadLength,     nReadLengthOver,     nReadLengthUnder,     readLengths,      0 );
		
		PrintMap( fOut, "NA read mapping multiplicity", nMultiplicity,   nMultiplicityOver,   nMultiplicityUnder,   multiplicities,   0 );
		
		PrintMap( fOut, "RQ read map quality",          nMappingQuality, nMappingQualityOver, nMappingQualityUnder, mappingQualities, 0 );
		
		n = sprintf( buffer, "MM read map mismatch (-mm/-mmp: %4.2f)", _allowedMismatch );
		if ( n > 1024 ) {
			printf("ERROR: The buffer for MM title is insufficient.\n");
			exit(1);
		}
		PrintMap( fOut, buffer,                nMismatch,       nMismatchOver,       nMismatchUnder,       mismatches,       0 );

		// print mate multiplicity combinations (MC)
		uint64_t total = U + M + F + X;
		fprintf( fOut, "\nMC mate multiplicity combinations\n");
		fprintf( fOut, "\tTOT\tMEAN\tSTD\tIN\tOVER\tUNDER\n" );
		fprintf( fOut, "\t%lu\t0\t0\t%lu\t0\t0\n", total, total );
		fprintf( fOut, "\tbin\tcombo\tn\tcum\tlabel\n");
		fprintf( fOut, "\t1\t1\t%lu\t%lu\tUniquely-aligned (U)\n", U, U );
		fprintf( fOut, "\t2\t2\t%lu\t%lu\tMultiply-aligned (M)\n", M, U + M );
		fprintf( fOut, "\t3\t3\t%lu\t%lu\tFiltered-out (F)\n",     F, U + M + F );
		fprintf( fOut, "\t4\t4\t%lu\t%lu\tUnmapped (X)\n",         X, U + M + F + X );
		
		// print pair multiplicity combinations (RC)
		total = UU + UM + MM + UF + MF + UX + MX + FF + FX + XX;
		fprintf( fOut, "\nRC pair multiplicity combinations\n");
		fprintf( fOut, "\tTOT\tMEAN\tSTD\tIN\tOVER\tUNDER\n" );
		fprintf( fOut, "\t%lu\t0\t0\t%lu\t0\t0\n", total, total );
		fprintf( fOut, "\tbin\tcombo\tn\tcum\tlabel\n");
		fprintf( fOut, "\t1\t1\t%lu\t%lu\tU-U\n", UU, UU );
		fprintf( fOut, "\t2\t2\t%lu\t%lu\tU-M\n", UM, UU + UM );
		fprintf( fOut, "\t3\t3\t%lu\t%lu\tM-M\n", MM, UU + UM + MM );
		fprintf( fOut, "\t4\t4\t%lu\t%lu\tU-F\n", UF, UU + UM + MM + UF );
		fprintf( fOut, "\t5\t5\t%lu\t%lu\tM-F\n", MF, UU + UM + MM + UF + MF );
		fprintf( fOut, "\t6\t6\t%lu\t%lu\tU-X\n", UX, UU + UM + MM + UF + MF + UX );
		fprintf( fOut, "\t7\t7\t%lu\t%lu\tM-X\n", MX, UU + UM + MM + UF + MF + UX + MX );
		fprintf( fOut, "\t8\t8\t%lu\t%lu\tF-F\n", FF, UU + UM + MM + UF + MF + UX + MX + FF );
		fprintf( fOut, "\t9\t9\t%lu\t%lu\tF-X\n", FX, UU + UM + MM + UF + MF + UX + MX + FF + FX );
		fprintf( fOut, "\t10\t10\t%lu\t%lu\tX-X\n", XX, total );

		// print local search rescues
		total = UU_localRescue + UU_localConsistance + UM_localRescue + UM_localConsistance + MM_localRescue + MM_localConsistance;
		fprintf( fOut, "\nLS local search\n");
		fprintf( fOut, "\tTOT\tMEAN\tSTD\tIN\tOVER\tUNDER\n" );
		fprintf( fOut, "\t%lu\t0\t0\t%lu\t0\t0\n", total, total );
		fprintf( fOut, "\tbin\tcombo\tn\tcum\tlabel\n");
		fprintf( fOut, "\t1\t1\t%lu\t%lu\tUU-LocalRescue\n", UU_localRescue, UU_localRescue );
		fprintf( fOut, "\t2\t2\t%lu\t%lu\tUM-LocalRescue\n", UM_localRescue, UU_localRescue + UM_localRescue );
		fprintf( fOut, "\t3\t3\t%lu\t%lu\tMM-LocalRescue\n", MM_localRescue, UU_localRescue + UM_localRescue + MM_localRescue );
		fprintf( fOut, "\t4\t4\t%lu\t%lu\tUU-localConsistance\n", UU_localConsistance, UU_localRescue + UM_localRescue + MM_localRescue + UU_localConsistance );
		fprintf( fOut, "\t5\t5\t%lu\t%lu\tUM-localConsistance\n", UM_localConsistance, UU_localRescue + UM_localRescue + MM_localRescue + UU_localConsistance + UM_localConsistance );
		fprintf( fOut, "\t6\t6\t%lu\t%lu\tMM-localConsistance\n", MM_localConsistance, total);

		// print orientation combinations (OC)
		total = f1_f2 + f1_r2 + r1_f2 + r1_r2 + f2_f1 + f2_r1 + r2_f1 + r2_r1;
		fprintf( fOut, "\nOC orientation combinations\n");
		fprintf( fOut, "\tTOT\tMEAN\tSTD\tIN\tOVER\tUNDER\n" );
		fprintf( fOut, "\t%lu\t0\t0\t%lu\t0\t0\n", total, total );
		fprintf( fOut, "\tbin\tcombo\tn\tcum\tlabel\n");
		fprintf( fOut, "\t1\t1\t%lu\t%lu\tf1-f2\n", f1_f2, f1_f2 );
		fprintf( fOut, "\t2\t2\t%lu\t%lu\tf1-r2\n", f1_r2, f1_f2 + f1_r2 );
		fprintf( fOut, "\t3\t3\t%lu\t%lu\tr1-f2\n", r1_f2, f1_f2 + f1_r2 + r1_f2 );
		fprintf( fOut, "\t4\t4\t%lu\t%lu\tr1-r2\n", r1_r2, f1_f2 + f1_r2 + r1_f2 + r1_r2 );
		fprintf( fOut, "\t5\t5\t%lu\t%lu\tf2-f1\n", f2_f1, f1_f2 + f1_r2 + r1_f2 + r1_r2 + f2_f1 );
		fprintf( fOut, "\t6\t6\t%lu\t%lu\tf2-r1\n", f2_r1, f1_f2 + f1_r2 + r1_f2 + r1_r2 + f2_f1 + f2_r1 );
		fprintf( fOut, "\t7\t7\t%lu\t%lu\tr2-f1\n", r2_f1, f1_f2 + f1_r2 + r1_f2 + r1_r2 + f2_f1 + f2_r1 + r2_f1 );
		fprintf( fOut, "\t8\t8\t%lu\t%lu\tr2-r1\n", r2_r1, f1_f2 + f1_r2 + r1_f2 + r1_r2 + f2_f1 + f2_r1 + r2_f1 + r2_r1 );
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

	//char* qPtr;
	//char* rPtr;
	//SaveReadLength( al1.BaseQualities.Length() );
	if ( al1.IsMapped ) {
		//qPtr = (char*)al1.Query.CData();
		//rPtr = (char*)al1.Reference.CData();
		//uint32_t mappedLength = 0;
		//for ( uint32_t i = 0; i < al1.Query.Length(); ++i ) {
		//	if ( ( *qPtr != '-' ) && ( *rPtr != 'Z' ) )
		//		mappedLength++;
		//
		//	++qPtr;
		//	++rPtr;
		//}

		SaveReadLength( al1.MappedLength );
		SaveMultiplicity( al1.NumMapped );
		SaveMappingQuality( al1.Quality );
		SaveMismatch( al1.NumMismatches );
	}
	SaveMateMultiplicity( al1.NumMapped, al1.IsFilteredOut );

	if ( isPairedEnd ) {
		//SaveReadLength( al2.BaseQualities.Length() );
		if ( al2.IsMapped ) {
			//qPtr = (char*)al2.Query.CData();
			//rPtr = (char*)al2.Reference.CData();
			//uint32_t mappedLength = 0;
			//for ( uint32_t i = 0; i < al2.Query.Length(); ++i ) {
			//	if ( ( *qPtr != '-' ) && ( *rPtr != 'Z' ) )
			//		mappedLength++;
			//
			//	++qPtr;
			//	++rPtr;
			//}

			SaveReadLength( al2.MappedLength );
			SaveMultiplicity( al2.NumMapped );
			SaveMappingQuality( al2.Quality );
			SaveMismatch( al2.NumMismatches );
		}
		SaveMateMultiplicity( al2.NumMapped, al2.IsFilteredOut );

		SavePairMultiplicity( al1.NumMapped, al2.NumMapped, al1.WasRescued, al2.WasRescued, al1.IsFilteredOut, al2.IsFilteredOut, ( al1.IsResolvedAsProperPair & al2.IsResolvedAsProperPair ) );
		//if ( ( ( al1.NumMapped == 0 ) && ( al2.NumMapped == 1 ) ) || ( ( al1.NumMapped == 1 ) && ( al2.NumMapped == 0 ) ) )
		//	non_unique++;

		//if ( ( ( al1.NumMapped == 0 ) && ( al2.NumMapped > 1 ) ) || ( ( al1.NumMapped > 1 ) && ( al2.NumMapped == 0 ) ) )
		//	non_multiple++;

		if ( ( al1.NumMapped == 1 ) && ( al2.NumMapped == 1 ) ) {
		//	unique_unique++;
			SaveFragment( al1, al2, tech );
			SaveModel( al1, al2 );
			//SaveIsize( al1, al2 );
		}

		//if ( ( ( al1.NumMapped == 1 ) && ( al2.NumMapped > 1 ) ) || ( ( al1.NumMapped > 1 ) && ( al2.NumMapped == 1 ) ) )
		//	unique_multiple++;

		//if ( ( al1.NumMapped > 0 ) && ( al2.NumMapped > 1 ) )
		//	multiple_multiple++;
	}
}
