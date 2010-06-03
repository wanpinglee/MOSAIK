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
	inline void FindMinElement ( vector<AlignedReadPair>& input, vector<AlignedReadPair>::iterator& minIte ) {

        	minIte = input.begin();
	        for ( vector< AlignedReadPair >::iterator ite = input.begin() + 1; ite != input.end(); ite++ ){
		        if ( ite->read.Name.empty() )
				continue;
			                
			if ( minIte->read.Name.empty() )
                	        minIte = ite;
	                else if ( minIte->read.Name > ite->read.Name )
        	                minIte = ite;
		}

		if ( minIte->read.Name.empty() ) {
        	       	cout << "ERROR: The input vector is empty." << endl;
               		exit(1);
	       	}

	}


	// given a vector of filenames, open them as AlignmentReader
	inline void OpenMosaikReader ( vector<MosaikReadFormat::CAlignmentReader>& readers, const vector<string>& files ) {
        
		MosaikReadFormat::CAlignmentReader *reader;
	        readers.reserve(files.size());

	        for ( unsigned int i = 0; i < files.size(); i++ ) {
	                reader = new MosaikReadFormat::CAlignmentReader;
	                reader->Open( files[i] );
	                readers.push_back(*reader);
	        }
	}



	inline bool LoadNextReadPair ( MosaikReadFormat::CAlignmentReader& reader, unsigned int readerNo, vector<AlignedReadPair>& reads ) {


		if ( readerNo > reads.size() - 1 ) {
			cout << "ERROR: The reader number is worng." << endl;
			exit(1);
		}
		
        	Mosaik::AlignedRead mr;
	        AlignedReadPair     mrp;

        	mr.Clear();
	        reads[readerNo].Clear();
		
		if ( !reader.LoadNextRead(mr) ) 
        	        return false;
	        else{
        	        mrp.read         = mr;
	               	mrp.owner        = readerNo;
	                reads[readerNo]  = mrp;
        	        return true;
	        }
	}
  
	// keep the proper pair according to the given fragment length
	inline void KeepProperPair ( Mosaik::AlignedRead& mr, const unsigned int mfl ) {
		
		// sort alignments by their positions	
		mr.SortAlignment();
		
		vector<Alignment>::iterator lastMinM2 = mr.Mate2Alignments.begin();

		for ( vector<Alignment>::iterator ite = mr.Mate1Alignments.begin(); ite != mr.Mate1Alignments.end(); ite++ ) {
			for ( vector<Alignment>::iterator ite2 = lastMinM2; ite2 != mr.Mate2Alignments.end(); ite2++ ) {
				unsigned int length = ( ite->ReferenceBegin > ite2->ReferenceBegin) ? ite->ReferenceBegin - ite2->ReferenceBegin : ite2->ReferenceBegin - ite->ReferenceBegin;
				if ( length > mfl ) {
					lastMinM2 = ( ite->ReferenceBegin > ite2->ReferenceBegin) ? ite2 : lastMinM2;
					continue;
				}
				else {
					if ( ite->IsReverseStrand != ite2->IsReverseStrand ) {
						ite->Mark  = true;
						ite2->Mark = true;
					}

				}
			}

		}

		
		vector<Alignment> newM1;
		newM1.reserve(mr.Mate1Alignments.size());
		for ( unsigned int i = 0; i < mr.Mate1Alignments.size(); i++ ) {
			if ( mr.Mate1Alignments[i].Mark )
				newM1.push_back( mr.Mate1Alignments[i] );
				
		}

		// if mate1 becomes unique or empty, we have to append one mate back
		if ( (newM1.size() == 1 && mr.Mate1Alignments.size() > 1) || (newM1.empty() && mr.Mate1Alignments.size() == 1) )  {
			for ( unsigned int i = 0; i < mr.Mate1Alignments.size(); i++ ) {
				if ( !mr.Mate1Alignments[i].Mark ) {
					newM1.push_back( mr.Mate1Alignments[i] );
					break;
				}
			}
		} else if ( newM1.empty() && mr.Mate1Alignments.size() > 1 ) {
			bool found = false;
			for ( unsigned int i = 0; i < mr.Mate1Alignments.size(); i++ ) {
				if ( !mr.Mate1Alignments[i].Mark ) {
					newM1.push_back( mr.Mate1Alignments[i] );
					if ( found )
						break;
					else
						found = true;
				}
			}
		}

		mr.Mate1Alignments.clear();
		mr.Mate1Alignments = newM1;


		vector<Alignment> newM2;
		newM2.reserve(mr.Mate2Alignments.size());
		for ( unsigned int i = 0; i < mr.Mate2Alignments.size(); i++ ) {
			if ( mr.Mate2Alignments[i].Mark )
				newM2.push_back( mr.Mate2Alignments[i] );
				
		}
		// if mate2 becomes unique or empty, we have to append one mate back
		if ( (newM2.size() == 1 && mr.Mate2Alignments.size() > 1) || (newM2.empty() && mr.Mate2Alignments.size() == 1) )  {
			for ( unsigned int i = 0; i < mr.Mate2Alignments.size(); i++ ) {
				if ( !mr.Mate2Alignments[i].Mark ) {
					newM2.push_back( mr.Mate2Alignments[i] );
					break;
				}
			}
		} else if ( newM2.empty() && mr.Mate2Alignments.size() > 1 ) {
			bool found = false;
			for ( unsigned int i = 0; i < mr.Mate2Alignments.size(); i++ ) {
				if ( !mr.Mate2Alignments[i].Mark ) {
					newM2.push_back( mr.Mate2Alignments[i] );
					if ( found )
						break;
					else
						found = true;
				}
			}
		}

		mr.Mate2Alignments.clear();
		mr.Mate2Alignments = newM2;
	}
	
}
#endif
