/*
 * =====================================================================================
 *
 *       Filename:  AlignedArchiveMerge.cpp
 *
 *    Description:  Merge a bunch of MOSAIK aligned archives
 *
 *        Version:  1.0
 *        Created:  03/16/2010 04:07:26 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Wan-Ping Lee
 *        Company:  Garth Lab, Biology, Boston College
 *
 * =====================================================================================
 */

#include "ArchiveSort.h"


CArchiveSort::CArchiveSort ( string inputFilename, string outputFilename, unsigned int *readCounter, pthread_mutex_t *readCounterMutex )
	:_inputFilename(inputFilename)
	,_outputFilename(outputFilename)
	,_readCounter(readCounter)
	,_readCounterMutex(readCounterMutex)
{
	_alignedReadCacheSize = 50000;
}

CArchiveSort::~CArchiveSort (){
	/*     
	for ( unsigned int i = 0; i < _tempFiles.size(); i++ ){
		for ( unsigned int j = 0; j < _tempFiles[i].size(); j++ )
			rm(_tempFiles[i][j].c_str());
	}
	*/

	
}

void CArchiveSort::SortNStoreCache( vector<string>& tempFiles, list<Mosaik::AlignedRead>& _alignedReadCache ) {
	
	//_alignedReadCache.sort(CmpAlignedRead);
	_alignedReadCache.sort();

	string tempFilename;
	CFileUtilities::GetTempFilename(tempFilename);
	tempFiles.push_back(tempFilename);

	MosaikReadFormat::CAlignmentWriter writer;
	writer.Open(tempFilename, *_referenceSequences, _readGroups, _alignmentStatus);
	writer.AdjustPartitionSize(1000);

	for ( list<Mosaik::AlignedRead>::iterator ite = _alignedReadCache.begin(); ite != _alignedReadCache.end(); ite++ ) 
		writer.SaveAlignedRead( *ite );
	
	
	writer.Close();
}


void CArchiveSort::SortNStoreTemp( vector<string>& tempFiles ){
	
	unsigned int nTemp = tempFiles.size();
	unsigned int nRead = 0;

	// initialize MOSAIK readers for all temp files
	vector< MosaikReadFormat::CAlignmentReader > readers;
	SortNMergeUtilities::OpenMosaikReader( readers, tempFiles );

	Mosaik::AlignedRead mr;
	unsigned int nDone = 0;
	vector<bool> done(nTemp, false);
	vector< SortNMergeUtilities::AlignedReadPair > reads(nTemp);

	// load the top element in each temp
	for ( unsigned int i = 0; i < nTemp; i++ ) {
		if ( !SortNMergeUtilities::LoadNextReadPair(readers[i], i, reads) ) {
			done[i] = true;
			nDone++;
		}
	}

	// prepare writer
	MosaikReadFormat::CAlignmentWriter writer;
	writer.Open(_outputFilename, *_referenceSequences, _readGroups, _alignmentStatus);
	// TODO: Consider what is the perfect number.
	writer.AdjustPartitionSize(1000);

	// pick the min one
	vector< SortNMergeUtilities::AlignedReadPair >::iterator ite;
	unsigned int tempNo = UINT_MAX;
	while ( nDone != ( nTemp - 1 ) ) {
		SortNMergeUtilities::FindMinElement(reads, ite);
		writer.SaveAlignedRead(ite->read);
		nRead++;

		pthread_mutex_lock(_readCounterMutex);
		*_readCounter = *_readCounter + 1;
		pthread_mutex_unlock(_readCounterMutex);
		
		tempNo = ite->owner;
		
		if ( tempNo >= nTemp ) {
			cout << "ERROR: Read temporary wrongly." << endl;
			exit(1);
		}
		
		if ( !SortNMergeUtilities::LoadNextReadPair(readers[tempNo], tempNo, reads) ) {
			if ( done[tempNo] ) {
				cout << "ERROR: The temporary file has been empty." << endl;
				exit(1);
			}
			done[tempNo] = true;
			nDone++;
		}
		
	}

	//cout << "\t\t-" << nRead << " reads are finished." << endl;
  
	// look for the remaining ones and store them
	for ( unsigned int i = 0; i < nTemp; i++ ) {
		if ( !done[i] ) {

			writer.SaveAlignedRead(reads[i].read);
			nRead++;
	
			pthread_mutex_lock(_readCounterMutex);
			*_readCounter = *_readCounter + 1;
			pthread_mutex_unlock(_readCounterMutex);

			while ( true ) {
				mr.Clear();
				if ( !readers[i].LoadNextRead(mr) ) 
					break;
				
				writer.SaveAlignedRead(mr);
				nRead++;

				pthread_mutex_lock(_readCounterMutex);
				*_readCounter = *_readCounter + 1;
				pthread_mutex_unlock(_readCounterMutex);
			}
		}
	}

	//cout << "\t\t-" << nRead << " reads are finished." << endl;

	writer.Close();
	
	// close readers
	for ( unsigned int i = 0; i < readers.size(); i++ ) {
		readers[i].Close();
	}
}

/*  
void CopyReferenceString( Mosaik::RefVector& refVec ){
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
	cout << endl;
}
*/

void CArchiveSort::Sort(){

	MosaikReadFormat::CAlignmentReader reader;
	reader.Open( _inputFilename );

		
	// set up archive header
	_referenceSequences = reader.GetReferenceSequences();
	reader.GetReadGroups(_readGroups);
	_alignmentStatus    = reader.GetStatus();
	
	unsigned int totalRead = reader.GetNumReads();
	unsigned int count     = 0;
	unsigned int nRead     = 0;
	vector<string> tempFiles;
		
	list<Mosaik::AlignedRead> _alignedReadCache;

	// allocate space for _alignedReadCache
	if ( totalRead < _alignedReadCacheSize )
		_alignedReadCache.resize( totalRead );
	else
		_alignedReadCache.resize( _alignedReadCacheSize );

	list<Mosaik::AlignedRead>::iterator ite;
	ite = _alignedReadCache.begin();

	
	Mosaik::AlignedRead mr;
	while ( true ) {
			
		mr.Clear();
		if ( !reader.LoadNextRead(mr) )
			break;

		// sort alignments by their positions
		bool isMm = ( mr.Mate1Alignments.size() > 1 ) && ( mr.Mate2Alignments.size() > 1 );
		bool isUm = (( mr.Mate1Alignments.size() == 1 ) && ( mr.Mate2Alignments.size() > 1 )) || (( mr.Mate1Alignments.size() > 1 ) && ( mr.Mate2Alignments.size() == 1 ));
		if ( isMm || isUm ) 
			SortNMergeUtilities::KeepProperPair(mr, 750);

		if ( count == _alignedReadCacheSize ) {
			// sort and store alignment
			SortNStoreCache(tempFiles, _alignedReadCache);
			//cout << "\t\t-" << nRead << " reads are sorted." << endl;
			pthread_mutex_lock(_readCounterMutex);
			*_readCounter = *_readCounter + _alignedReadCache.size();
			pthread_mutex_unlock(_readCounterMutex);

			_alignedReadCache.clear();
			// allocate space for _alignedReadCache
			unsigned int remaining = totalRead - nRead;
			if ( remaining < _alignedReadCacheSize )
				_alignedReadCache.resize( remaining );
			else
				_alignedReadCache.resize( _alignedReadCacheSize );
			
			// reset the pointer
			ite = _alignedReadCache.begin();
			count = 0;

		}

		*ite = mr;
		ite++;
		count++;
		nRead++;
	}
		
	reader.Close();

	// sort and store the last part
	if ( count > 0 ) {
		SortNStoreCache(tempFiles, _alignedReadCache);
		//cout << "\t\t-" << nRead << " reads are sorted." << endl;
		pthread_mutex_lock(_readCounterMutex);
		*_readCounter = *_readCounter + _alignedReadCache.size();
		pthread_mutex_unlock(_readCounterMutex);
	}
	_alignedReadCache.clear();

	// sort and store the current aligned archive
	//cout << "Sorting whole aligned archive." << endl;
	SortNStoreTemp(tempFiles);

	// delete the temperary files
	for ( unsigned int i = 0; i < tempFiles.size(); i++ )
		rm ( tempFiles[i].c_str() );

}


//void CArchiveSort::Sort( vector<string> filenames ){
	
	// initialize
	//Mosaik::CAlignmentReader reader;
	//reader.Open( filenames[0] );
	//_referenceSequences = reader.GetReferenceData();
	//_readGroups = reader.GetReadGroups();
	//_alignmentStatus = reader.GetStatus();

	//reader.Close();

	//_refIndex.resize( filenames.size(), 0 );
	//_nBase.resize   ( filenames.size(), 0 );
	//_nRead.resize   ( filenames.size(), 0 );

	//for ( unsigned int i = 0; i < filenames.size(); i++ ) {
	//	cout << "Sorting partial " << i+1 << " aligned archive." << endl;
		// open our MOSAIK alignment reader
	//	reader.Open( filenames[i] );

		// load reference data
	//	Mosaik::RefVector referenceSequences;
	//	referenceSequences = reader.GetReferenceData();
	//	CopyReferenceString( referenceSequences );
	//	_referenceSequences.insert( _referenceSequences.end(), referenceSequences.begin(), referenceSequences.end() );
	//	_refIndex[i] = ( i == 0 ) ? referenceSequences.size() : referenceSequences.size() + _refIndex[i-1];


	//	_nBase[i] = reader.GetNumBases();
	//	_nRead[i] = reader.GetNumReads();
	//	reader.Close();

	//	cout << "\t# of Reads: " << _nRead[i] << " in " << i+1 << " aligned archive." << endl;

	//	Sort( (void*) &filenames[i] );
	//}

//}
