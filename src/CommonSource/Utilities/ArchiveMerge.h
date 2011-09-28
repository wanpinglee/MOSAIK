/*
 * =====================================================================================
 *
 *       Filename:  ArchiveMerge.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/25/2010 10:10:13 AM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Wan-Ping Lee
 *        Company:  Marth Lab., Biology, Boston College
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <limits.h>

#include "AlignedRead.h"
#include "Alignment.h"
#include "AlignmentReader.h"
#include "AlignmentWriter.h"
#include "BamWriter.h"
#include "BestNSecondBestSelection.h"
#include "Entropy.h"
#include "FileUtilities.h"
#include "MosaikString.h"
#include "QualityNeuralNetwork.h"
#include "ReadGroup.h"
#include "ReadStatus.h"
#include "SortNMergeUtilities.h"
#include "StatisticsMaps.h"
#include "ZaTager.h"

class CArchiveMerge
{
	public:
		struct StatisticsCounters {
	                // reads
			uint64_t UU;
	                uint64_t UM;
	                uint64_t UF;
	                uint64_t MM;
		        uint64_t MF;
			uint64_t UX;
			uint64_t MX;
			uint64_t FF;
			uint64_t FX;
			uint64_t XX;
			uint64_t UU_localRescue;
			uint64_t UU_localConsistance;
			uint64_t UM_localRescue;
			uint64_t UM_localConsistance;
			uint64_t MM_localRescue;
			uint64_t MM_localConsistance;
		        // mates
			uint64_t FilteredOutMates;
			uint64_t MultipleMates;
			uint64_t UniqueMates;
			uint64_t Unmapped;

			StatisticsCounters() 
				// reads
				: UU(0)
				, UM(0)
				, UF(0)
				, MM(0)
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
				// mates
				, FilteredOutMates(0)
				, MultipleMates(0)
				, UniqueMates(0)
				, Unmapped(0)
			{}

		};
	
	public:
		/* ====================  LIFECYCLE     ======================================= */
		//CArchiveMerge (vector < string > inputFilenames, string outputFilename, unsigned int *readNo);                             /* constructor */
		CArchiveMerge (
			const vector <string>& inputFilenames, 
			const string& outputFilename, 
			uint64_t            *readNo, 
			const bool&          isSolid,
			const string&        commandLine,
			const string&        paired_end_ann_file,
			const string&        single_end_ann_file,
			const unsigned int&  fragmentLength = 0,
			const unsigned int&  localAlignmentSearchRadius = 0,
			const bool&          hasSpecial = false,
			const unsigned char& statMappingQuality = 20 );


		void Merge();

		void GetStatisticsCounters ( StatisticsCounters& counter );
		void PrintStatisticsMaps( 
		    const string filename, 
		    const vector<MosaikReadFormat::ReadGroup>& readGroup, 
		    const uint8_t fragmentLength, 
		    const uint8_t localSearchRadius, 
		    const float allowedMismatch );
	private:
		/* ====================  DATA MEMBERS  ======================================= */
		vector < string > _inputFilenames;
		string _outputFilename;

		//unsigned int               _nMaxAlignment;
		uint64_t*                  _readNo;
		bool                       _isSolid;
		unsigned int               _expectedFragmentLength;
		unsigned int               _localAlignmentSearchRadius;
		vector< unsigned int >     _refIndex;
		vector<ReferenceSequence>           _referenceSequences;
		vector<ReferenceSequence>           _referenceSequencesWoSpecial;
		vector<MosaikReadFormat::ReadGroup> _readGroups;
		AlignmentStatus _alignmentStatus;
		bool _isPairedEnd;
		bool _hasSpecial;
		unsigned char _statMappingQuality;
		SequencingTechnologies _sequencingTechnologies;
		map<unsigned int, MosaikReadFormat::ReadGroup> _readGroupsMap;
		
		string _specialArchiveName;
		string _specialCode1;
		string _specialCode2;
		MosaikReadFormat::CAlignmentReader _specialReader;
		Mosaik::AlignedRead                _specialAl;
		bool                               _specialArchiveEmpty;
		vector<ReferenceSequence>          _specialReferenceSequences;

		StatisticsCounters _counters;
		CStatisticsMaps    _statisticsMaps;

		BamHeader _sHeader; // special reads
		BamHeader _rHeader; // regular bam

		CBamWriter _sBam; // multiply alignments
		CBamWriter _rBam;

		// ZA tagers
		CZaTager za1, za2;

		// Entropy
		Entropy _entropy;

		// neural-net
		QualityNeuralNetwork _mqCalculator;

		void UpdateReferenceIndex ( Mosaik::AlignedRead& mr, const unsigned int& owner );
		void CopyReferenceString( vector<ReferenceSequence>& refVec );
		void PrintReferenceSequence( vector<ReferenceSequence>& refVec );
		void CalculateStatisticsCounters( const Mosaik::AlignedRead& alignedRead );

		void WriteAlignment( Mosaik::AlignedRead& r );
		void SetAlignmentFlags(
		        Alignment& al,
		        const Alignment& mate,
		        const bool& isPair,
		        const bool& isProperPair,
		        const bool& isFirstMate,
		        const bool& isPairTech,
		        const bool& isItselfMapped,
		        const bool& isMateMapped,
		        const Mosaik::AlignedRead& r);
		CArchiveMerge (const CArchiveMerge&);
		CArchiveMerge& operator= (const CArchiveMerge&);

}; /* -----  end of class CArchiveMerge  ----- */

