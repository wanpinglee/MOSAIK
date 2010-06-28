/*
 * =====================================================================================
 *
 *       Filename:  ArchiveSort.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/16/2010 04:17:31 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Wan-Ping Lee
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <list>
#include <algorithm>
#include <limits.h>

#include "Alignment.h"
#include "SortNMergeUtilities.h"
#include "FileUtilities.h"
#include "PosixThreads.h"
#include "ReferenceSequence.h"
#include "ReadGroup.h"
#include "ReadStatus.h"
#include "AlignmentWriter.h"
#include "AlignmentReader.h"
#include "MosaikString.h"

using namespace std;


/*
 * =====================================================================================
 *        Class:  AlignedArchiveMerge
 *  Description:  
 * =====================================================================================
 */

#ifndef CArchiveSort_H
#define CArchiveSort_H

class CArchiveSort {
	
	public:
		/* ====================  LIFECYCLE     ======================================= */
		/*  constructor */ 
		CArchiveSort ( string inputFilename, string outputFilename, unsigned int *readCounter, pthread_mutex_t *readCounterMutex, unsigned int medianFragmentLength );
		~CArchiveSort ();
		void Sort();

		/* ====================  ACCESSORS     ======================================= */

		/* ====================  MUTATORS      ======================================= */

		/* ====================  OPERATORS     ======================================= */

	protected:
		/* ====================  DATA MEMBERS  ======================================= */

	private:
		/* ====================  DATA MEMBERS  ======================================= */ 
		string _inputFilename;
		string _outputFilename;
		vector< string >          _tempSortedFiles;
		
		unsigned int _alignedReadCacheSize;
		unsigned int *_readCounter;
		pthread_mutex_t *_readCounterMutex;
		unsigned int _medianFragmentLength;
		unsigned int _extendedFragmentLength;

		vector<ReferenceSequence>*          _referenceSequences;
		vector<MosaikReadFormat::ReadGroup> _readGroups;
		AlignmentStatus _alignmentStatus;

		void SortNStoreCache( vector< string >& tempFiles, list<Mosaik::AlignedRead>& _alignedReadCache );
		void SortNStoreTemp ( vector< string >& tempFiles );


}; /* -----  end of class AlignedArchiveMerge  ----- */

#endif
