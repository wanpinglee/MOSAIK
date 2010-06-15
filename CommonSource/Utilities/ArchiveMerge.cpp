/*
 * =====================================================================================
 *
 *       Filename:  ArchiveMerge.cpp
 *
 *    Description:  Given sorted archived, the program would merge them
 *
 *        Version:  1.0
 *        Created:  03/25/2010 10:07:12 AM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Wan-Ping Lee
 *        Company:  Marth Lab., Biology, Boston College
 *
 * =====================================================================================
 */

#include "ArchiveMerge.h"

CArchiveMerge::CArchiveMerge (vector < string > inputFilenames, string outputFilename, unsigned int *readNo)
	: _inputFilenames(inputFilenames)
	, _outputFilename(outputFilename)
	, _nMaxAlignment(1000)
	, _readNo(readNo)
{
	
	MosaikReadFormat::CAlignmentReader reader;
	reader.Open( _inputFilenames[0] );
	reader.GetReadGroups(_readGroups);
	_alignmentStatus = reader.GetStatus();
	reader.Close();

	_refIndex.resize( _inputFilenames.size(), 0 );

	for ( unsigned int i = 0; i < _inputFilenames.size(); i++ ) {
		reader.Open( _inputFilenames[i] );
		
		vector<ReferenceSequence>*  referenceSequences;
		referenceSequences = reader.GetReferenceSequences();
		CopyReferenceString( *referenceSequences );
		_referenceSequences.insert( _referenceSequences.end(), referenceSequences->begin(), referenceSequences->end() );
		_refIndex[i] = ( i == 0 ) ? referenceSequences->size() : referenceSequences->size() + _refIndex[i-1];

		reader.Close();
	}
	
}

CArchiveMerge::CArchiveMerge (vector < string > inputFilenames, string outputFilename, unsigned int nMaxAlignment, unsigned int *readNo)
	: _inputFilenames(inputFilenames)
	, _outputFilename(outputFilename)
	, _nMaxAlignment (nMaxAlignment)
	, _readNo(readNo)
{
	MosaikReadFormat::CAlignmentReader reader;
	reader.Open( _inputFilenames[0] );
	reader.GetReadGroups(_readGroups);
	_alignmentStatus = reader.GetStatus();
	reader.Close();

	_refIndex.resize( _inputFilenames.size(), 0 );


	for ( unsigned int i = 0; i < _inputFilenames.size(); i++ ) {
		reader.Open( _inputFilenames[i] );
		
		vector<ReferenceSequence> referenceSequences;
		reader.GetReferenceSequences(referenceSequences);
		CopyReferenceString( referenceSequences );
		_referenceSequences.insert( _referenceSequences.end(), referenceSequences.begin(), referenceSequences.end() );
		_refIndex[i] = ( i == 0 ) ? referenceSequences.size() : referenceSequences.size() + _refIndex[i-1];

		referenceSequences.clear();
		reader.Close();
	}

}


void CArchiveMerge::PrintReferenceSequence( vector<ReferenceSequence>& refVec ){

	for ( unsigned int i = 0; i < refVec.size(); i++ ) {
		cout << i << endl;
		cout << "\t" << "#Aligned   " << refVec[i].NumAligned << endl;
		cout << "\t" << "Begin      " << refVec[i].Begin << endl;
		cout << "\t" << "End        " << refVec[i].End << endl;
		cout << "\t" << "#Bases     " << refVec[i].NumBases << endl;
		cout << "\t" << "Name       " << refVec[i].Name << endl;
		cout << "\t" << "Bases      " << refVec[i].Bases << endl;
		cout << "\t" << "AssemblyId " << refVec[i].GenomeAssemblyID << endl;
		cout << "\t" << "Species    " << refVec[i].Species << endl;
		cout << "\t" << "MD5        " << refVec[i].MD5 << endl;
		cout << "\t" << "URI        " << refVec[i].URI << endl;
	}
	cout << endl;

}


void CArchiveMerge::CopyReferenceString( vector<ReferenceSequence>& refVec ){
	char temp[1024];

	for ( unsigned int i = 0; i < refVec.size(); i++ ) {
		
		if ( refVec[i].Name.size() > 1024 ) {
			cout << "ERROR: The length of the reference name is larger than 1024." << endl;
			exit(1);
		}

		memcpy( temp, refVec[i].Name.c_str(), refVec[i].Name.size() );
		temp[ refVec[i].Name.size() ] = 0;
		refVec[i].Name = temp;


		//if ( refVec[i].Bases.size() > 1024 ) {
		//	cout << "ERROR:  The length of the reference bases is larger than 1024." << endl;
		//	exit(1);
		//}

		//memcpy( temp, refVec[i].Bases.c_str(), refVec[i].Bases.size() );
		//temp[ refVec[i].Bases.size() ] = 0;
		//refVec[i].Bases = temp;
		
		if ( refVec[i].GenomeAssemblyID.size() > 1024 ) {
			cout << "ERROR:  The length of the reference genomeAssemblyID is larger than 1024." << endl;
			exit(1);
		}

		memcpy( temp, refVec[i].GenomeAssemblyID.c_str(), refVec[i].GenomeAssemblyID.size() );
		temp[ refVec[i].GenomeAssemblyID.size() ] = 0;
		refVec[i].GenomeAssemblyID = temp;


		if ( refVec[i].Species.size() > 1024 ) {
			cout << "ERROR:  The length of the reference species is larger than 1024." << endl;
			exit(1);
		}

		memcpy( temp, refVec[i].Species.c_str(), refVec[i].Species.size() );
		temp[ refVec[i].Species.size() ] = 0;
		refVec[i].Species = temp;


		if ( refVec[i].MD5.size() > 1024 ) {
			cout << "ERROR:  The length of the reference MD5 is larger than 1024." << endl;
			exit(1);
		}

		memcpy( temp, refVec[i].MD5.c_str(), refVec[i].MD5.size() );
		temp[ refVec[i].MD5.size() ] = 0;
		refVec[i].MD5 = temp;


		if ( refVec[i].URI.size() > 1024 ) {
			cout << "ERROR:  The length of the reference URI is larger than 1024." << endl;
			exit(1);
		}

		memcpy( temp, refVec[i].URI.c_str(), refVec[i].URI.size() );
		temp[ refVec[i].URI.size() ] = 0;
		refVec[i].URI = temp;

	}
}

// update the reference index
inline void CArchiveMerge::UpdateReferenceIndex ( Mosaik::AlignedRead& mr, const unsigned int& owner ) {

	if ( owner >= _refIndex.size() ) {
		cout << "ERROR: The ID of the temporary file " << owner << " is unexpected." << endl;
		exit(1);
        }

        unsigned int offset = ( owner == 0 ) ? 0 : _refIndex[ owner - 1 ];

        for ( unsigned int i = 0; i < mr.Mate1Alignments.size(); i++ )
		mr.Mate1Alignments[i].ReferenceIndex += offset;

        for ( unsigned int i = 0; i < mr.Mate2Alignments.size(); i++ )
		mr.Mate2Alignments[i].ReferenceIndex += offset;   
         
}

void CArchiveMerge::GetStatisticsCounters ( CArchiveMerge::StatisticsCounters& counter ) {
	counter = _counters;
}

void CArchiveMerge::CalculateStatisticsCounters( const Mosaik::AlignedRead& alignedRead ) {
	// reads
	_counters.AlignedReads++;
	
	if ( alignedRead.Mate1Alignments.size() > 1 && alignedRead.Mate2Alignments.size() > 1 )
		_counters.BothNonUniqueReads++;
	
	if ( alignedRead.Mate1Alignments.size() == 1 && alignedRead.Mate2Alignments.size() == 1 )
		_counters.BothUniqueReads++;

	if ( (alignedRead.Mate1Alignments.size() == 1 && alignedRead.Mate2Alignments.size() > 1) || (alignedRead.Mate1Alignments.size() > 1 && alignedRead.Mate2Alignments.size() == 1) )
		_counters.OneNonUniqueReads++;

	if ( (alignedRead.Mate1Alignments.size() != 0 && alignedRead.Mate2Alignments.size() == 0) || (alignedRead.Mate1Alignments.size() == 0 && alignedRead.Mate2Alignments.size() != 0) )
		_counters.OrphanedReads++;
	
	// mates
	if ( alignedRead.Mate1Alignments.size() == 0 )
		_counters.FilteredOutMates++;

	if ( alignedRead.Mate2Alignments.size() == 0 )
		_counters.FilteredOutMates++;

	if ( alignedRead.Mate1Alignments.size() > 1 )
		_counters.NonUniqueMates++;

	if ( alignedRead.Mate2Alignments.size() > 1 )
		_counters.NonUniqueMates++;

	if ( alignedRead.Mate1Alignments.size() == 1 )
		_counters.UniqueMates++;

	if ( alignedRead.Mate2Alignments.size() == 1 )
		_counters.UniqueMates++;
}


void CArchiveMerge::Merge() {
	
	//unsigned int nRead = 0;
	unsigned int nTemp = _inputFilenames.size();
	
	// initialize MOSAIK readers for all temp files
	vector< MosaikReadFormat::CAlignmentReader > readers;
	SortNMergeUtilities::OpenMosaikReader( readers, _inputFilenames );

	
	Mosaik::AlignedRead mr;
	unsigned int nDone = 0;
	vector< bool > done(nTemp, false);
	vector< SortNMergeUtilities::AlignedReadPair > reads(nTemp);
	// first load
	for ( unsigned int i = 0; i < nTemp; i++ ) {
		if ( !SortNMergeUtilities::LoadNextReadPair(readers[i], i, reads) ) {
			done[i] = true;
			nDone++;
		} else {
			UpdateReferenceIndex( reads[i].read, i );
			*_readNo = *_readNo + 1;
		}
	}

	
	// prepare MOSAIK writer
	MosaikReadFormat::CAlignmentWriter writer;
	writer.Open(_outputFilename, _referenceSequences, _readGroups, _alignmentStatus);
	writer.AdjustPartitionSize(1000);

	
	// pick the min one
	vector< SortNMergeUtilities::AlignedReadPair >::iterator ite;
	Mosaik::AlignedRead minRead;

	SortNMergeUtilities::FindMinElement(reads, ite);
	//ite = min_element( reads.begin(), reads.end(), CmpAlignedReadPair );
	minRead = ite->read;
	if ( !SortNMergeUtilities::LoadNextReadPair(readers[ite->owner], ite->owner, reads) ) {
		done[ite->owner] = true;
		nDone++;
	} else {
		UpdateReferenceIndex( reads[ite->owner].read, ite->owner );
		*_readNo = *_readNo + 1;
	}

	unsigned int tempNo = UINT_MAX;
	while ( nDone != nTemp - 1 ) {
		SortNMergeUtilities::FindMinElement(reads, ite);
		//ite = min_element( reads.begin(), reads.end(), CmpAlignedReadPair );
		
		if ( ite->read.Name > minRead.Name ) {
			//if ( ++nRead % 100000 == 0 )
			//	cout << "\t" << nRead << " have been merged." << endl;

			//cout << minRead.Mate1Alignments.size() << "\t" << minRead.Mate2Alignments.size() << endl;
			if ( minRead.Mate1Alignments.size() > _nMaxAlignment )
				minRead.Mate1Alignments.erase( minRead.Mate1Alignments.begin() + _nMaxAlignment, minRead.Mate1Alignments.end() );
			if ( minRead.Mate2Alignments.size() > _nMaxAlignment )
				minRead.Mate2Alignments.erase( minRead.Mate2Alignments.begin() + _nMaxAlignment, minRead.Mate2Alignments.end() );
			
			writer.SaveAlignedRead(minRead);
			CalculateStatisticsCounters(minRead);
			minRead.Clear();
			minRead = ite->read;
		} else {
			// merge their alignments
			bool isFull = ( ite->read.Mate1Alignments.size() + minRead.Mate1Alignments.size() > minRead.Mate1Alignments.max_size() )
			            ||( ite->read.Mate2Alignments.size() + minRead.Mate2Alignments.size() > minRead.Mate2Alignments.max_size() );
			if ( isFull ) {
				cout << "ERROR: Too many alignments waiting for writing." << endl;
				exit(1);
			}

			if ( ite->read.Mate1Alignments.size() > 0 )
				minRead.Mate1Alignments.insert( minRead.Mate1Alignments.begin(), ite->read.Mate1Alignments.begin(), ite->read.Mate1Alignments.end() );
			if ( ite->read.Mate2Alignments.size() > 0 )
				minRead.Mate2Alignments.insert( minRead.Mate2Alignments.begin(), ite->read.Mate2Alignments.begin(), ite->read.Mate2Alignments.end() );

			// accordant flag
			minRead.IsLongRead |= ite->read.IsLongRead;

			// TODO: think this more
			if ( minRead.Mate1Alignments.size() > _nMaxAlignment ) {
				random_shuffle(minRead.Mate1Alignments.begin(), minRead.Mate1Alignments.end());
				minRead.Mate1Alignments.erase( minRead.Mate1Alignments.begin() + _nMaxAlignment, minRead.Mate1Alignments.end() );
			}
			if ( minRead.Mate2Alignments.size() > _nMaxAlignment ) {
				random_shuffle(minRead.Mate2Alignments.begin(), minRead.Mate2Alignments.end());
				minRead.Mate2Alignments.erase( minRead.Mate2Alignments.begin() + _nMaxAlignment, minRead.Mate2Alignments.end() );
			}

		}
		
		tempNo = ite->owner;
		
		if ( tempNo >= nTemp ) {
			cout << "ERROR: Read ID is wrong." << endl;
			exit(1);
		}
		
		if ( !SortNMergeUtilities::LoadNextReadPair(readers[tempNo], tempNo, reads) ) {
			done[tempNo] = true;
			nDone++;
		} else {
			UpdateReferenceIndex( reads[tempNo].read, tempNo );
			*_readNo = *_readNo + 1;
		}
	}	
	

	unsigned int owner = UINT_MAX;
	for ( unsigned int i = 0; i < nTemp; i++ ) {
		if ( !done[i] ) {

			owner = reads[i].owner;

			if ( reads[i].read.Name > minRead.Name ) {
				//if ( ++nRead % 100000 == 0 )
				//	cout << "\t" << nRead << " have been merged." << endl;

				//cout << minRead.Mate1Alignments.size() << "\t" << minRead.Mate2Alignments.size() << endl;
				writer.SaveAlignedRead(minRead);
				CalculateStatisticsCounters(minRead);
				minRead.Clear();
				minRead = reads[i].read;
			} else {
				if ( reads[i].read.Mate1Alignments.size() > 0 )
					minRead.Mate1Alignments.insert( minRead.Mate1Alignments.end(), reads[i].read.Mate1Alignments.begin(), reads[i].read.Mate1Alignments.end() );
				if ( reads[i].read.Mate2Alignments.size() > 0 )
					minRead.Mate2Alignments.insert( minRead.Mate2Alignments.end(), reads[i].read.Mate2Alignments.begin(), reads[i].read.Mate2Alignments.end() );
				
				// accordant flag
				minRead.IsLongRead |= ite->read.IsLongRead;
			}

			while ( true ) {
				mr.Clear();
				if ( !readers[i].LoadNextRead(mr) ) 
					break;
				else {
					UpdateReferenceIndex( mr, owner );
					*_readNo = *_readNo + 1;
				}
				
				if ( mr.Name > minRead.Name ) {
					//UpdateReferenceIndex( mr, owner );
					//if ( ++nRead % 100000 == 0 )
					//	cout << "\t" << nRead << " have been merged." << endl;

					//cout << minRead.Mate1Alignments.size() << "\t" << minRead.Mate2Alignments.size() << endl;
					writer.SaveAlignedRead(minRead);
					CalculateStatisticsCounters(minRead);
					minRead.Clear();
					minRead = mr;
				} else {
					//UpdateReferenceIndex( mr, owner );
					if( mr.Mate1Alignments.size() > 0 )
						minRead.Mate1Alignments.insert( minRead.Mate1Alignments.end(), mr.Mate1Alignments.begin(), mr.Mate1Alignments.end() );
					if( mr.Mate2Alignments.size() > 0 )
						minRead.Mate2Alignments.insert( minRead.Mate2Alignments.end(), mr.Mate2Alignments.begin(), mr.Mate2Alignments.end() );
					
					// accordant flag
					minRead.IsLongRead |= ite->read.IsLongRead;
				}

			}

		}
	}


	//UpdateReferenceIndex(minRead, owner);
	//cout << minRead.Mate1Alignments.size() << "\t" << minRead.Mate2Alignments.size() << endl;
	writer.SaveAlignedRead(minRead);
	CalculateStatisticsCounters(minRead);


	writer.Close();

	
	// close readers
	for ( unsigned int i = 0; i < readers.size(); i++ )
		readers[i].Close();
}
