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
#include "FileUtilities.h"
#include "MosaikString.h"
#include "ReadGroup.h"
#include "ReadStatus.h"
#include "SortNMergeUtilities.h"
#include "StatisticsMaps.h"
#include "ZaTager.h"


/*
 * =====================================================================================
 *        Class:  CArchiveMerge
 *  Description:  
 * =====================================================================================
 */
class CArchiveMerge
{
	public:

		struct StatisticsCounters {
	                // reads
			uint64_t AlignedReads;
	                uint64_t BothNonUniqueReads;
	                uint64_t BothUniqueReads;
	                uint64_t OneNonUniqueReads;
		        uint64_t OrphanedReads;
		        // mates
			uint64_t FilteredOutMates;
			uint64_t NonUniqueMates;
			uint64_t UniqueMates;

			StatisticsCounters() 
				// reads
				: AlignedReads(0)
				, BothNonUniqueReads(0)
				, BothUniqueReads(0)
				, OneNonUniqueReads(0)
				, OrphanedReads(0)
				// mates
				, FilteredOutMates(0)
				, NonUniqueMates(0)
				, UniqueMates(0)
			{}

		};
	
	public:
		/* ====================  LIFECYCLE     ======================================= */
		//CArchiveMerge (vector < string > inputFilenames, string outputFilename, unsigned int *readNo);                             /* constructor */
		CArchiveMerge (
			vector < string > inputFilenames, 
			string outputFilename, 
			unsigned int *readNo, 
			const unsigned int fragmentLength = 0, 
			const bool hasSpecial = false );

		void Merge();

		void GetStatisticsCounters ( StatisticsCounters& counter );
		void PrintStatisticsMaps( const string filename, const string readGroupId );
		/* ====================  ACCESSORS     ======================================= */

		/* ====================  MUTATORS      ======================================= */

		/* ====================  OPERATORS     ======================================= */

	protected:
		/* ====================  DATA MEMBERS  ======================================= */

	private:
		/* ====================  DATA MEMBERS  ======================================= */
		vector < string > _inputFilenames;
		string _outputFilename;

                //unsigned int               _nMaxAlignment;
		unsigned int*              _readNo;
		unsigned int               _expectedFragmentLength;
		vector< unsigned int >     _refIndex;
		vector<ReferenceSequence>           _referenceSequences;
		vector<ReferenceSequence>           _referenceSequencesWoSpecial;
		vector<MosaikReadFormat::ReadGroup> _readGroups;
		AlignmentStatus _alignmentStatus;
		bool _isPairedEnd;
		bool _hasSpecial;
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
		BamHeader _uHeader; // unaligned reads
		BamHeader _rHeader; // regular bam

		CBamWriter _sBam; // multiply alignments
		CBamWriter _uBam; // unaligned reads
		CBamWriter _rBam;

		// ZA tagers
		CZaTager za1, za2;
		
		
		inline void UpdateReferenceIndex ( Mosaik::AlignedRead& mr, const unsigned int& owner );
		void CopyReferenceString( vector<ReferenceSequence>& refVec );
		void PrintReferenceSequence( vector<ReferenceSequence>& refVec );
		void CalculateStatisticsCounters( const Mosaik::AlignedRead& alignedRead );

		void WriteAlignment( Mosaik::AlignedRead& r );
		inline void SetAlignmentFlags(
		        Alignment& al,
		        const Alignment& mate,
		        const bool& isPair,
		        const bool& isProperPair,
		        const bool& isFirstMate,
		        const bool& isPairTech,
		        const bool& isItselfMapped,
		        const bool& isMateMapped,
		        const Mosaik::AlignedRead& r);

}; /* -----  end of class CArchiveMerge  ----- */

