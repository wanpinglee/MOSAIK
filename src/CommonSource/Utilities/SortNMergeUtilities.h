/*
 * =====================================================================================
 *
 *       Filename:  SortNMergeUtilities.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/25/2010 12:58:20 PM
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
#include <iostream>

#include "Mosaik.h"
#include "AlignedRead.h"
#include "AlignmentReader.h"
#include "Alignment.h"


using namespace std;

#ifndef SortNMergeUtilities_H
#define SortNMergeUtilities_H

namespace SortNMergeUtilities {	
	struct AlignedReadPair {
		Mosaik::AlignedRead read;
		unsigned int owner;

		AlignedReadPair()
		    : read()
		    , owner(0)
		{}

		bool Clear(){
			read.Clear();
			owner = 0;
			return true;
		}

		bool operator<(AlignedReadPair& x) {
			if ( x.read.Name.empty() ) return false;
			if (   read.Name.empty() ) return true;
			return read.Name < x.read.Name;
		}

		bool operator>(AlignedReadPair& x) {
			if ( x.read.Name.empty() ) return false;
			if (   read.Name.empty() ) return true;
		}
	};



	// given a vector of AlignedReadPair, find the min element and put it in minIte
	void FindMinElement ( vector<AlignedReadPair>& input, vector<AlignedReadPair>::iterator& minIte );

	// given a vector of filenames, open them as AlignmentReader
	void OpenMosaikReader ( vector<MosaikReadFormat::CAlignmentReader*>& readers, const vector<string>& files );

	// close readers
	void CloseMosaikReader ( vector<MosaikReadFormat::CAlignmentReader*>& readers );

	bool LoadNextReadPair ( MosaikReadFormat::CAlignmentReader* reader, unsigned int readerNo, vector<AlignedReadPair>& reads );
  
	// keep the proper pair according to the given fragment length
	void KeepProperPair ( Mosaik::AlignedRead& mr, const unsigned int mfl );
	
}
#endif
