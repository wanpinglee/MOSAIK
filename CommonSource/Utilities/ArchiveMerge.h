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

#include "Alignment.h"
#include "SortNMergeUtilities.h"
#include "FileUtilities.h"
#include "AlignmentWriter.h"
#include "AlignmentReader.h"
#include "MosaikString.h"
#include "ReadGroup.h"
#include "ReadStatus.h"


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
		CArchiveMerge (vector < string > inputFilenames, string outputFilename, unsigned int *readNo);                             /* constructor */
		CArchiveMerge (vector < string > inputFilenames, string outputFilename, unsigned int nMaxAlignment, unsigned int *readNo);

		void Merge();

		void GetStatisticsCounters ( StatisticsCounters& counter );
		/* ====================  ACCESSORS     ======================================= */

		/* ====================  MUTATORS      ======================================= */

		/* ====================  OPERATORS     ======================================= */

	protected:
		/* ====================  DATA MEMBERS  ======================================= */

	private:
		/* ====================  DATA MEMBERS  ======================================= */
		vector < string > _inputFilenames;
		string _outputFilename;

                unsigned int               _nMaxAlignment;
		unsigned int*              _readNo;
		vector< unsigned int >     _refIndex;
		vector<ReferenceSequence>           _referenceSequences;
		vector<MosaikReadFormat::ReadGroup> _readGroups;
		AlignmentStatus _alignmentStatus;

		StatisticsCounters _counters;
		
		
		inline void UpdateReferenceIndex ( Mosaik::AlignedRead& mr, const unsigned int& owner );
		void CopyReferenceString( vector<ReferenceSequence>& refVec );
		void PrintReferenceSequence( vector<ReferenceSequence>& refVec );
		void CalculateStatisticsCounters( const Mosaik::AlignedRead& alignedRead );

}; /* -----  end of class CArchiveMerge  ----- */

