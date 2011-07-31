#ifndef _STATISTICSMAPS_H_
#define _STATISTICSMAPS_H_

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <cmath>
#include <string>
#include <vector>

#include "Alignment.h"
#include "BamHeader.h"
#include "Read.h"
#include "ReadGroup.h"
#include "SequencingTechnologies.h"

using namespace std;

class CStatisticsMaps {
	public:
		CStatisticsMaps( void );
		CStatisticsMaps( int64_t mfl );
		~CStatisticsMaps( void );
		
		void SaveRecord( const Alignment& al1, const Alignment& al2, const bool isPairedEnd, const SequencingTechnologies& tech );
		void PrintMaps( const char* filename, const vector<MosaikReadFormat::ReadGroup>& readGroup, const unsigned char statMappingQuality );
		void SetExpectedStatistics( const uint32_t fragmentLength, const uint32_t localSearchRadius, const float allowedMismatch );
		void Reset( void );

		void SetLfMin( int64_t flMin );
		void SetMqMin( uint8_t mqMin );

	private:
		// copy constructor
		CStatisticsMaps( const CStatisticsMaps& copy );
		// assing operator
		CStatisticsMaps& operator= ( const CStatisticsMaps& copy );

		inline void SavePairMultiplicity( 
			  const unsigned int nMate1Alignments
			, const unsigned int nMate2Alignments
			, const bool mate1Rescued
			, const bool mate2Rescued
			, const bool mate1FilteredOut
			, const bool mate2FilteredOut
			, const bool isProperPair );
		inline void SaveMateMultiplicity( const unsigned int nAlignment, const bool isFilteredOut );
		inline void SaveModel( const Alignment& al1, const Alignment& al2 );
		inline void SaveFragment( const Alignment& al1, const Alignment& al2, const SequencingTechnologies& tech );
		inline void SaveIsize( const Alignment& al1, const Alignment& al2 );
		inline void SaveReadLength( const unsigned int length );
		inline void SaveMultiplicity( const unsigned int nAlignment );
		inline void SaveMappingQuality( const unsigned char mq );
		inline void SaveMismatch( const unsigned short mm );
		inline void PrintMap( 
			  FILE *const fOut 
			, const char* title 
			, const uint64_t& size
			, const uint64_t& over
			, const uint64_t& under
			, const uint64_t *const array
			, const int64_t& start  );

		// original setting
		uint32_t _fragmentLength;
		uint32_t _localSearchRadius;
		float   _allowedMismatch;
		bool    _setExpectedStatistics;

		// fragment length map
		const uint64_t nFragment;
		uint64_t nFrangmentOver;
		uint64_t nFrangmentUnder;
		int64_t  minFragment;
		uint64_t* fragments;
		// Isize map
		//const uint64_t nIsize;
		//uint64_t nIsizeOver;
		//uint64_t nIsizeUnder;
		//uint64_t* isizes;
		// mapped read length map
		const uint64_t nReadLength;
		uint64_t nReadLengthOver;
		uint64_t nReadLengthUnder;
		uint64_t* readLengths;
		// multiplicity map
		const uint64_t nMultiplicity;
		uint64_t nMultiplicityOver;
		uint64_t nMultiplicityUnder;
		uint64_t* multiplicities;
		// mapping quality map
		const uint16_t nMappingQuality;
		uint64_t nMappingQualityOver;
		uint64_t nMappingQualityUnder;
		uint8_t  minMappingQuality;
		uint64_t* mappingQualities;
		// mismatch map
		const uint64_t nMismatch;
		uint64_t nMismatchOver;
		uint64_t nMismatchUnder;
		uint64_t* mismatches;
		// mate map
		uint64_t U;  // uniquely aligned mates
		uint64_t M;  // multiply aligned mates
		uint64_t X;  // unmapped mates
		uint64_t F;  // filtered-out mates; mates are mapped but cannot pass the mismatch filters
		// pair map
		uint64_t UU; // uniquely-uniquely aligned pairs
		uint64_t UM; // uniquely-multiply aligned pairs
		uint64_t MM; // multiply-multiply aligned pairs
		uint64_t UF; // uniquely-filtered aligned pairs; in which filtered out mates are mapped but cannot pass the mismatch filters
		uint64_t MF; // multiply-filtered aligned pairs; in which filtered out mates are mapped but cannot pass the mismatch filters
		uint64_t UX; // uniquely-unmapped aligned pairs
		uint64_t MX; // multiply-unmapped aligned pairs
		uint64_t FF; // filtered-filtered aligned pairs
		uint64_t FX; // filtered-unmapped aligned pairs
		uint64_t XX; // unmapped-unmapped aligned pairs
		uint64_t UU_localRescue;
		uint64_t UU_localConsistance;
		uint64_t UM_localRescue;
		uint64_t UM_localConsistance;
		uint64_t MM_localRescue;
		uint64_t MM_localConsistance;
		// model map
		uint64_t f1_f2;
		uint64_t f1_r2;
		uint64_t r1_f2;
		uint64_t r1_r2;
		uint64_t f2_f1;
		uint64_t f2_r1;
		uint64_t r2_f1;
		uint64_t r2_r1;

};

#endif
