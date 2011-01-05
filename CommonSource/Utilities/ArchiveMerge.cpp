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

/*
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
*/

CArchiveMerge::CArchiveMerge ( vector < string > inputFilenames, string outputFilename, unsigned int *readNo, const unsigned int fragmentLength )
	: _inputFilenames(inputFilenames)
	, _outputFilename(outputFilename)
	, _readNo(readNo)
	, _expectedFragmentLength( fragmentLength )
{
	MosaikReadFormat::CAlignmentReader reader;
	reader.Open( _inputFilenames[0] );
	reader.GetReadGroups(_readGroups);
	_alignmentStatus = reader.GetStatus();
	reader.Close();

	_isPairedEnd = ( _alignmentStatus && AS_PAIRED_END_READ != 0 ) ? true : false;

	_refIndex.resize( _inputFilenames.size(), 0 );

	for ( vector<MosaikReadFormat::ReadGroup>::iterator ite = _readGroups.begin(); ite != _readGroups.end(); ++ite ) {
		ite->ReadGroupCode = ite->GetCode( *ite );
		_readGroupsMap[ ite->ReadGroupCode ] = *ite ;
	}


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

	_mHeader.SortOrder = SORTORDER_UNSORTED;
	_uHeader.SortOrder = SORTORDER_UNSORTED;
	_rHeader.SortOrder = SORTORDER_UNSORTED;

	_mHeader.pReferenceSequences = &_referenceSequences;
	_uHeader.pReferenceSequences = &_referenceSequences;
	_rHeader.pReferenceSequences = &_referenceSequences;

	_mHeader.pReadGroups = &_readGroups;
	_uHeader.pReadGroups = &_readGroups;
	_rHeader.pReadGroups = &_readGroups;

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
	
	unsigned int nMate1Alignments = 0;
	unsigned int nMate2Alignments = 0;

	if ( !alignedRead.Mate1Alignments.empty() )
		nMate1Alignments = alignedRead.Mate1Alignments[0].NumMapped;
	if ( !alignedRead.Mate2Alignments.empty() )
		nMate2Alignments = alignedRead.Mate2Alignments[0].NumMapped;
	
	// reads
	if ( ( nMate1Alignments > 0 ) || ( nMate2Alignments > 0 ) )
		_counters.AlignedReads++;
	
	if ( nMate1Alignments > 1 && nMate2Alignments > 1 )
		_counters.BothNonUniqueReads++;
	
	if ( nMate1Alignments == 1 && nMate2Alignments == 1 )
		_counters.BothUniqueReads++;

	if ( ( nMate1Alignments == 1 && nMate2Alignments > 1 ) || ( nMate1Alignments > 1 && nMate2Alignments == 1) )
		_counters.OneNonUniqueReads++;

	if ( ( nMate1Alignments != 0 && nMate2Alignments == 0 ) || ( nMate1Alignments == 0 && nMate2Alignments != 0 ) )
		_counters.OrphanedReads++;
	
	// mates
	if ( nMate1Alignments == 0 )
		_counters.FilteredOutMates++;

	if ( nMate2Alignments == 0 )
		_counters.FilteredOutMates++;

	if ( nMate1Alignments > 1 )
		_counters.NonUniqueMates++;

	if ( nMate2Alignments > 1 )
		_counters.NonUniqueMates++;

	if ( nMate1Alignments == 1 )
		_counters.UniqueMates++;

	if ( nMate2Alignments == 1 )
		_counters.UniqueMates++;
}


void CArchiveMerge::WriteAlignment( Mosaik::AlignedRead& r ) {

	unsigned int nMate1Alignments = 0;
	unsigned int nMate2Alignments = 0;

	vector<Alignment> newMate1Set, newMate2Set;

	for ( vector<Alignment>::iterator ite = r.Mate1Alignments.begin(); ite != r.Mate1Alignments.end(); ++ite ) {
		nMate1Alignments += ite->NumMapped;
		if ( ite->IsMapped ) newMate1Set.push_back( *ite );
	}
	
	for ( vector<Alignment>::iterator ite = r.Mate2Alignments.begin(); ite != r.Mate2Alignments.end(); ++ite ) {
		nMate2Alignments += ite->NumMapped;
		if ( ite->IsMapped ) newMate2Set.push_back( *ite );
	}

	if ( nMate1Alignments > 0 ) {
		r.Mate1Alignments.clear();
		r.Mate1Alignments = newMate1Set;
	}

	if ( nMate2Alignments > 0 ) {
		r.Mate2Alignments.clear();
		r.Mate2Alignments = newMate2Set;
	}
		

	//r.Mate1Alignments.clear();
	//r.Mate2Alignments.clear();
	//r.Mate1Alignments = newMate1Set;
	//r.Mate2Alignments = newMate2Set;

	const bool isMate1Unique   = ( nMate1Alignments == 1 ) ? true : false;
	const bool isMate2Unique   = ( nMate2Alignments == 1 ) ? true : false;
	//const bool isMate1Aligned  = ( nMate1Alignments > 0 ) ? true : false;
	//const bool isMate2Aligned  = ( nMate2Alignments > 0 ) ? true : false;
	const bool isMate1Multiple = ( nMate1Alignments > 1 ) ? true : false;
	const bool isMate2Multiple = ( nMate2Alignments > 1 ) ? true : false;
	bool isMate1Empty    = ( nMate1Alignments == 0 ) ? true : false;
	bool isMate2Empty    = ( nMate2Alignments == 0 ) ? true : false;

//cout << nMate1Alignments << "\t" << nMate2Alignments << endl;

	// UU, UM, and MM pair
	if ( ( isMate1Unique && isMate2Unique )
		|| ( isMate1Unique && isMate2Multiple )
		|| ( isMate1Multiple && isMate2Unique )
		|| ( isMate1Multiple && isMate2Multiple ) ) {

//cout << "case1" << endl;
			
		if ( ( isMate1Unique && isMate2Multiple )
			|| ( isMate1Multiple && isMate2Unique )
			|| ( isMate1Multiple && isMate2Multiple ) )
				BestNSecondBestSelection::Select( r.Mate1Alignments, r.Mate2Alignments, _expectedFragmentLength );

			isMate1Empty = r.Mate1Alignments.empty();
			isMate2Empty = r.Mate2Alignments.empty();
			
			// sanity check
			if ( isMate1Empty | isMate2Empty ) {
				cout << "ERROR: One of mate sets is empty after apllying best and second best selection." << endl;
				exit(1);
			}

			// patch the information for reporting
			Alignment al1 = r.Mate1Alignments[0], al2 = r.Mate2Alignments[0];
			
			// TODO: handle fragment length for others sequencing techs
			int fl = ( al1.IsReverseStrand ) 
				? 0 - ( 2 * _expectedFragmentLength ) 
				: 2 * _expectedFragmentLength;
			bool properPair1 = false, properPair2 = false;
			properPair1 = al1.SetPairFlags( al2, fl,  !al1.IsReverseStrand );
			properPair2 = al2.SetPairFlags( al1, -fl, !al2.IsReverseStrand );

			if ( properPair1 != properPair2 ) {
				cout << "ERROR: An inconsistent proper pair is found." << endl;
				exit(1);
			}

			CZaTager za1, za2;
			const char* zaTag1 = za1.GetZaTag( al1, al2, true );
			const char* zaTag2 = za2.GetZaTag( al2, al1, false );

			SetAlignmentFlags( al1, al2, true, properPair1, true, _isPairedEnd, true, true, r );
			SetAlignmentFlags( al2, al1, true, properPair2, false, _isPairedEnd, true, true, r );

			_rBam.SaveAlignment( al1, zaTag1 );
			_rBam.SaveAlignment( al2, zaTag2 );
		
	// UX and MX pair
	} else if ( ( isMate1Empty || isMate2Empty )
		&&  !( isMate1Empty && isMate2Empty ) ) {
		
//cout << "case2" << endl;

		if ( isMate1Multiple || isMate2Multiple ) 
			BestNSecondBestSelection::Select( r.Mate1Alignments, r.Mate2Alignments, _expectedFragmentLength );

		isMate1Empty = r.Mate1Alignments.empty();
		isMate2Empty = r.Mate2Alignments.empty();
		
		bool isFirstMate;
		if ( !isMate1Empty ) {
			isFirstMate = true;
		} else if ( !isMate2Empty ) {
			isFirstMate = false;
		} else {
			cout << "ERROR: Both mates are empty after applying best and second best selection." << endl;
			exit(1);
		}
	
		// patch the information for reporting
		Alignment al = isFirstMate ? r.Mate1Alignments[0] : r.Mate2Alignments[0];
		Alignment unmappedAl = !isFirstMate ? r.Mate1Alignments[0] : r.Mate2Alignments[0];;

		SetAlignmentFlags( al, unmappedAl, false, false, isFirstMate, _isPairedEnd, true, false, r );
		SetAlignmentFlags( unmappedAl, al, true, false, !isFirstMate, _isPairedEnd, false, true, r );

		_rBam.SaveAlignment( al, 0 );
		_uBam.SaveAlignment( unmappedAl, 0, true );

	
	// XX
	} else if ( isMate1Empty && isMate2Empty ) {
		
//cout << "case3" << endl;

		Alignment unmappedAl1, unmappedAl2;
		unmappedAl1 = r.Mate1Alignments[0];
		unmappedAl2 = r.Mate2Alignments[0];

		SetAlignmentFlags( unmappedAl1, unmappedAl2, false, false, true, _isPairedEnd, false, false, r );
		SetAlignmentFlags( unmappedAl2, unmappedAl1, false, false, true, _isPairedEnd, false, false, r );

		_uBam.SaveAlignment( unmappedAl1, 0, true );
		_uBam.SaveAlignment( unmappedAl2, 0, true );
	
	} else {
		cout << "ERROR: Unknown pairs." << endl;
		cout << r.Mate1Alignments.size() << "\t" << r.Mate2Alignments.size() << endl;
		exit(1);
	}

}

inline void CArchiveMerge::SetAlignmentFlags( 
	Alignment& al,
	const Alignment& mate,
	const bool& isPair,
	const bool& isProperPair,
	const bool& isFirstMate,
	const bool& isPairTech,
	const bool& isItselfMapped,
	const bool& isMateMapped,
	const Mosaik::AlignedRead& r) {

	al.IsResolvedAsPair       = isPair;
	al.IsResolvedAsProperPair = isProperPair;
	al.IsFirstMate            = isFirstMate;
	al.IsPairedEnd            = isPairTech;
	al.Name                   = r.Name;
	al.IsMateReverseStrand    = mate.IsReverseStrand;
	al.MateReferenceIndex     = mate.ReferenceIndex;
	al.MateReferenceBegin     = mate.ReferenceBegin;
	al.IsMapped               = isItselfMapped;
	al.IsMateMapped           = isMateMapped;

	map<unsigned int, MosaikReadFormat::ReadGroup>::iterator rgIte;
	rgIte = _readGroupsMap.find( r.ReadGroupCode );
	// sanity check
	if ( rgIte == _readGroupsMap.end() ) {
		cout << "ERROR: ReadGroup cannot be found." << endl;
		exit(1);
	} else
		al.ReadGroup = rgIte->second.ReadGroupID;
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
	//MosaikReadFormat::CAlignmentWriter writer;
	//writer.Open(_outputFilename, _referenceSequences, _readGroups, _alignmentStatus, ALIGNER_SIGNATURE);
	//writer.AdjustPartitionSize(1000);
	
	// prepare BAM writers
	_mBam.Open( _outputFilename + ".multiple.bam", _mHeader );
	_uBam.Open( _outputFilename + ".unaligned.bam", _uHeader );
	_rBam.Open( _outputFilename + ".bam", _rHeader );

	
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

			//cout << minRead.Mate1Alignments.size() << "\t" << minRead.Mate2Alignments.size() << endl;
			//if ( minRead.Mate1Alignments.size() > _nMaxAlignment )
			//	minRead.Mate1Alignments.erase( minRead.Mate1Alignments.begin() + _nMaxAlignment, minRead.Mate1Alignments.end() );
			//if ( minRead.Mate2Alignments.size() > _nMaxAlignment )
			//	minRead.Mate2Alignments.erase( minRead.Mate2Alignments.begin() + _nMaxAlignment, minRead.Mate2Alignments.end() );
			
			//BestNSecondBestSelection::Select( minRead.Mate1Alignments, minRead.Mate2Alignments, _expectedFragmentLength );
			WriteAlignment( minRead );

			//writer.SaveAlignedRead(minRead);
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
			//if ( minRead.Mate1Alignments.size() > _nMaxAlignment ) {
			//	random_shuffle(minRead.Mate1Alignments.begin(), minRead.Mate1Alignments.end());
			//	minRead.Mate1Alignments.erase( minRead.Mate1Alignments.begin() + _nMaxAlignment, minRead.Mate1Alignments.end() );
			//}
			//if ( minRead.Mate2Alignments.size() > _nMaxAlignment ) {
			//	random_shuffle(minRead.Mate2Alignments.begin(), minRead.Mate2Alignments.end());
			//	minRead.Mate2Alignments.erase( minRead.Mate2Alignments.begin() + _nMaxAlignment, minRead.Mate2Alignments.end() );
			//}

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
				//writer.SaveAlignedRead(minRead);
				WriteAlignment( minRead );
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
					//writer.SaveAlignedRead(minRead);
					WriteAlignment( minRead );
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
	//writer.SaveAlignedRead(minRead);
	WriteAlignment( minRead );
	CalculateStatisticsCounters(minRead);


	//writer.Close();
	
	// Close BAMs
	_mBam.Close();
	_uBam.Close();
	_rBam.Close();

	
	// close readers
	for ( unsigned int i = 0; i < readers.size(); i++ )
		readers[i].Close();
}
