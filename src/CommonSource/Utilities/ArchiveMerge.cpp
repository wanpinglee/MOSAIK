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

CArchiveMerge::CArchiveMerge ( 
	const vector <string>& inputFilenames, 
	const string& outputFilename, 
	uint64_t           *readNo,
	const bool&          isSolid,
	const string&        commandLine,
	const string&        paired_end_ann_file,
	const string&        single_end_ann_file,
	const unsigned int&  fragmentLength,
	const unsigned int&  localAlignmentSearchRadius,
	const bool&          hasSpecial,
	const unsigned char& statMappingQuality)
	
	: _inputFilenames(inputFilenames)
	, _outputFilename(outputFilename)
	, _readNo                    (readNo)
	, _isSolid                   (isSolid)
	, _expectedFragmentLength    (fragmentLength)
	, _localAlignmentSearchRadius(localAlignmentSearchRadius)
	, _refIndex()
	, _referenceSequences()
	, _referenceSequencesWoSpecial()
	, _readGroups()
	, _alignmentStatus(0)
	, _isPairedEnd(false)
	, _hasSpecial                (hasSpecial)
	, _statMappingQuality        (statMappingQuality)
	, _sequencingTechnologies()
	, _readGroupsMap()
	, _specialArchiveName()
	, _specialCode1()
	, _specialCode2()
	, _specialReader()
	, _specialAl()
	, _specialArchiveEmpty(true)
	, _specialReferenceSequences()
	, _counters()
	, _statisticsMaps()
	, _sHeader()
	, _rHeader()
	, _sBam()
	, _rBam()
	, za1()
	, za2()
	, _entropy()
	, _mqCalculator()
{
	//_statisticsMaps = mStatisticsMaps;
	
	MosaikReadFormat::CAlignmentReader reader;
	reader.Open( _inputFilenames[0] );
	reader.GetReadGroups(_readGroups);
	_alignmentStatus = reader.GetStatus();
	reader.Close();

	_isPairedEnd = ( ( _alignmentStatus & AS_PAIRED_END_READ ) != 0 ) ? true : false;

	_refIndex.resize( _inputFilenames.size(), 0 );

	for ( vector<MosaikReadFormat::ReadGroup>::iterator ite = _readGroups.begin(); ite != _readGroups.end(); ++ite ) {
		ite->ReadGroupCode = ite->GetCode( *ite );
		_readGroupsMap[ ite->ReadGroupCode ] = *ite ;
	}

	_sequencingTechnologies = _readGroups[0].SequencingTechnology;

	// the last archive is the one containing alignments located at special references

	for ( unsigned int i = 0; i < _inputFilenames.size(); i++ ) {
		reader.Open( _inputFilenames[i] );
		
		// grab reference info from the archive
		vector<ReferenceSequence> referenceSequences;
		reader.GetReferenceSequences(referenceSequences);

		CopyReferenceString( referenceSequences );
		
		_referenceSequences.insert( _referenceSequences.end(), referenceSequences.begin(), referenceSequences.end() );
		_refIndex[i] = ( i == 0 ) ? referenceSequences.size() : referenceSequences.size() + _refIndex[i-1];

		// don't include the special references
		if ( ( _hasSpecial ) && ( i != ( _inputFilenames.size() - 1 ) ) )
			_referenceSequencesWoSpecial.insert( _referenceSequencesWoSpecial.end(), referenceSequences.begin(), referenceSequences.end() );

		// includes the special references
		if ( ( _hasSpecial ) && ( i == ( _inputFilenames.size() - 1 ) ) )
			_specialReferenceSequences.insert( _specialReferenceSequences.end(), referenceSequences.begin(), referenceSequences.end() );

		referenceSequences.clear();
		reader.Close();
	}

	if ( !_hasSpecial )
		_referenceSequencesWoSpecial = _referenceSequences;

	_sHeader.SortOrder = SORTORDER_UNSORTED;
	//_uHeader.SortOrder = SORTORDER_UNSORTED;
	_rHeader.SortOrder = SORTORDER_UNSORTED;

	_sHeader.pReferenceSequences = &_referenceSequences;
	//_uHeader.pReferenceSequences = &_referenceSequencesWoSpecial;
	_rHeader.pReferenceSequences = &_referenceSequencesWoSpecial;

	_sHeader.pReadGroups = &_readGroups;
	//_uHeader.pReadGroups = &_readGroups;
	_rHeader.pReadGroups = &_readGroups;

	ProgramGroup pg;
	pg.ID = "MosaikAligner";
	stringstream ss;
	ss << (int)MOSAIK_MAJOR_VERSION << "." << (int)MOSAIK_MINOR_VERSION << "." << (int)MOSAIK_BUILD_VERSION;
	pg.VN = ss.str();
	pg.CL = commandLine;

	_sHeader.pg.ID = "MosaikAligner";
	_rHeader.pg.ID = "MosaikAligner";
	_sHeader.pg.VN = ss.str();
	_rHeader.pg.VN = ss.str();
	_sHeader.pg.CL = commandLine;
	_rHeader.pg.CL = commandLine;

	_mqCalculator.Open(paired_end_ann_file, single_end_ann_file);

}


void CArchiveMerge::PrintStatisticsMaps( 
	const string filename, 
	const vector<MosaikReadFormat::ReadGroup>& readGroup, 
	const uint8_t fragmentLength, 
	const uint8_t localSearchRadius, 
	const float allowedMismatch ) {
	
	_statisticsMaps.SetExpectedStatistics( fragmentLength, localSearchRadius, allowedMismatch );
	_statisticsMaps.PrintMaps( filename.c_str(), readGroup, _statMappingQuality );
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
void CArchiveMerge::UpdateReferenceIndex ( Mosaik::AlignedRead& mr, const unsigned int& owner ) {

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


	bool mate1FilteredOut = false;
	bool mate2FilteredOut = false;
	bool mate1Rescued     = false;
	bool mate2Rescued     = false;
	bool isProperPair     = false;

	if ( !alignedRead.Mate1Alignments.empty() ) {
		nMate1Alignments = alignedRead.Mate1Alignments[0].NumMapped;
		mate1FilteredOut = alignedRead.Mate1Alignments[0].IsFilteredOut;
		mate1Rescued     = alignedRead.Mate1Alignments[0].WasRescued;
		isProperPair     = alignedRead.Mate1Alignments[0].IsResolvedAsProperPair;
	}

	if ( !alignedRead.Mate2Alignments.empty() ) {
		nMate2Alignments = alignedRead.Mate2Alignments[0].NumMapped;
		mate2FilteredOut = alignedRead.Mate2Alignments[0].IsFilteredOut;
		mate2Rescued     = alignedRead.Mate2Alignments[0].WasRescued;
		isProperPair    &= alignedRead.Mate2Alignments[0].IsResolvedAsProperPair;
	}
	
	// =====
	// reads
	// =====
	if ( _isPairedEnd ) {
	// MM pairs
	if ( nMate1Alignments > 1 && nMate2Alignments > 1 ) {
		_counters.MM++;
		if ( mate1Rescued || mate2Rescued )
			_counters.MM_localRescue++;
		else if ( isProperPair )
			_counters.MM_localConsistance++;
	}
	
	// UU pairs
	if ( nMate1Alignments == 1 && nMate2Alignments == 1 ) {
		_counters.UU++;
		if ( mate1Rescued || mate2Rescued )
			_counters.UU_localRescue++;
		else if ( isProperPair )
			_counters.UU_localConsistance++;
	}

	// UM pairs
	if ( ( nMate1Alignments == 1 && nMate2Alignments > 1 ) || ( nMate1Alignments > 1 && nMate2Alignments == 1) ) {
		_counters.UM++;
		if ( mate1Rescued || mate2Rescued )
			_counters.UM_localRescue++;
		else if ( isProperPair )
			_counters.UM_localConsistance++;
	}

	// UF pairs
	if ( ( nMate1Alignments == 1 && ( ( nMate2Alignments == 0 ) && mate2FilteredOut ) ) || ( ( ( nMate1Alignments == 0 ) && mate1FilteredOut ) && nMate2Alignments == 1 ) ) 
		_counters.UF++;
	
	// MF pairs
	if ( ( nMate1Alignments > 1 && ( ( nMate2Alignments == 0 ) && mate2FilteredOut ) ) || ( ( ( nMate1Alignments == 0 ) && mate1FilteredOut ) && nMate2Alignments > 1 ) )
		_counters.MF++;
	
	// UX pairs
	if ( ( nMate1Alignments == 1 && ( ( nMate2Alignments == 0 ) && !mate2FilteredOut ) ) || ( ( ( nMate1Alignments == 0 ) && !mate1FilteredOut ) && nMate2Alignments == 1 ) )
		_counters.UX++;

	// MX pairs
	if ( ( nMate1Alignments > 1 && ( ( nMate2Alignments == 0 ) && !mate2FilteredOut ) ) || ( ( ( nMate1Alignments == 0 ) && !mate1FilteredOut ) && nMate2Alignments > 1 ) )
		_counters.MX++;

	// FF pairs
	if ( ( ( nMate1Alignments == 0 ) && mate1FilteredOut ) && ( ( nMate2Alignments == 0 ) && mate2FilteredOut ) )
		_counters.FF++;
	
	// FX pairs
	if ( ( ( ( nMate1Alignments == 0 ) && mate1FilteredOut ) && ( ( nMate2Alignments == 0 ) && !mate2FilteredOut ) )
		|| ( ( ( nMate1Alignments == 0 ) && !mate1FilteredOut ) && ( ( nMate2Alignments == 0 ) && mate2FilteredOut ) ) )
		 _counters.FX++;

	// XX pairs
	if ( ( ( nMate1Alignments == 0 ) && !mate1FilteredOut ) && ( ( nMate2Alignments == 0 ) && !mate2FilteredOut ) )
		_counters.XX++;
	}


	// =====
	// mates
	// =====
	if ( nMate1Alignments == 0 && mate1FilteredOut )
		_counters.FilteredOutMates++;
	
	if ( nMate1Alignments == 0 && !mate1FilteredOut )
		_counters.Unmapped++;

	if ( nMate1Alignments > 1 )
		_counters.MultipleMates++;

	if ( nMate1Alignments == 1 )
		_counters.UniqueMates++;

	if ( _isPairedEnd ) {
		if ( nMate2Alignments == 0 && mate2FilteredOut )
			_counters.FilteredOutMates++;

		if ( nMate2Alignments == 0 && !mate2FilteredOut )
			_counters.Unmapped++;

		if ( nMate2Alignments > 1 )
			_counters.MultipleMates++;

		if ( nMate2Alignments == 1 )
			_counters.UniqueMates++;
	}

}


void CArchiveMerge::WriteAlignment( Mosaik::AlignedRead& r ) {

	_specialCode1.clear();
	_specialCode2.clear();

	bool isMate1Special = false;
	bool isMate2Special = false;

	// grab special tag
	if ( _hasSpecial ) {
		// compare their names
		while ( ( _specialAl < r ) && !_specialArchiveEmpty ) {
			_specialAl.Clear();
			_specialArchiveEmpty = !_specialReader.LoadNextRead( _specialAl );
			if ( !_specialArchiveEmpty ) *_readNo = *_readNo + 1;
		}

		if ( _specialAl.Name == r.Name ) {
			if ( r.Mate1Alignments.size() > 1 ) {
				if ( _specialAl.Mate1Alignments[0].IsMapped ) {
					unsigned int referenceIndex = _specialAl.Mate1Alignments[0].ReferenceIndex;
					_specialCode1 = _specialReferenceSequences[ referenceIndex ].Species;
					_specialCode1.resize(3);
					_specialCode1[2] = 0;
					isMate1Special = true;
				}
			}

			if ( r.Mate2Alignments.size() > 1 ) {
				if ( _specialAl.Mate2Alignments[0].IsMapped ) {
					unsigned int referenceIndex = _specialAl.Mate2Alignments[0].ReferenceIndex;
					_specialCode2 = _specialReferenceSequences[ referenceIndex ].Species;
					_specialCode2.resize(3);
					_specialCode2[2] = 0;
					isMate2Special = true;
				}
			}
		}
	}

	unsigned int nMate1Alignments = 0;
	unsigned int nMate2Alignments = 0;
	unsigned int nMate1Hashes     = 0;
	unsigned int nMate2Hashes     = 0;
	float mate1SwScore     = 0;
	float mate2SwScore     = 0;
	float mate1NextSwScore = 0;
	float mate2NextSwScore = 0;

	vector<Alignment*> newMate1Set, newMate2Set;
	Mosaik::Read read;

	string mate1Cs, mate1Cq, mate2Cs, mate2Cq;
	bool isMate1FilteredOut = false, isMate2FilteredOut = false;
/*
if (r.Name == "15_23564567_23565516_0:0:0_1:0:0_7552") {
  cerr << "mate1\t" << mate1SwScore << "\t" << mate1NextSwScore << endl;
  for ( vector<Alignment>::iterator ite = r.Mate1Alignments.begin(); ite != r.Mate1Alignments.end(); ++ite ) {
    cerr << (*ite).ReferenceIndex << "\t" 
         << (*ite).ReferenceBegin << "\t" 
	 << (*ite).SwScore << "\t" 
	 << (*ite).NextSwScore << "\t"
	 << (*ite).NumHash << "\t" 
	 << (*ite).NumMapped << endl;
  }
  cerr << "mate2\t" << mate2SwScore << "\t" << mate2NextSwScore << endl;
  for ( vector<Alignment>::iterator ite = r.Mate2Alignments.begin(); ite != r.Mate2Alignments.end(); ++ite ) {
    cerr << (*ite).ReferenceIndex << "\t" 
         << (*ite).ReferenceBegin << "\t" 
	 << (*ite).SwScore << "\t" 
	 << (*ite).NextSwScore << "\t"
	 << (*ite).NumHash << "\t" 
	 << (*ite).NumMapped << endl;
  }
}
*/
	for ( vector<Alignment>::iterator ite = r.Mate1Alignments.begin(); ite != r.Mate1Alignments.end(); ++ite ) {
		nMate1Alignments   += ite->NumMapped;
		nMate1Hashes       += ite->NumHash;
		isMate1FilteredOut |= ite->IsFilteredOut;
		ite->SpecialCode    = _specialCode1;
		if (ite->SwScore >= mate1SwScore) {
			mate1NextSwScore = mate1SwScore;
			mate1SwScore     = ite->SwScore;
		} else if (ite->SwScore >= mate1NextSwScore){
			mate1NextSwScore = ite->SwScore;
		}
		if (ite->IsMapped) newMate1Set.push_back( &*ite );
		// the record is in the first archive, and contains complete bases and base qualities.
		if ( ite->Owner == 0 ) {
			read.Mate1.Bases     = ite->Query;
			read.Mate1.Qualities = ite->BaseQualities;
			if ( read.Mate1.Bases.Length() > 0 ) read.Mate1.Bases.Remove('-');
			if ( _isSolid ) {
				mate1Cs = ite->CsQuery;
				mate1Cq = ite->CsBaseQualities;
			}
		}
	}
	
	for ( vector<Alignment>::iterator ite = r.Mate2Alignments.begin(); ite != r.Mate2Alignments.end(); ++ite ) {
		nMate2Alignments   += ite->NumMapped;
		nMate2Hashes       += ite->NumHash;
		isMate2FilteredOut |= ite->IsFilteredOut;
		ite->SpecialCode    = _specialCode2;
		if (ite->SwScore >= mate2SwScore) {
			mate2NextSwScore = mate2SwScore;
			mate2SwScore     = ite->SwScore;
		} else if (ite->SwScore >= mate2NextSwScore) {
			mate2NextSwScore = ite->SwScore;
		}
		if (ite->IsMapped) newMate2Set.push_back( &*ite );
		// the record is in the first archive, and contains complete bases and base qualities.
		if (ite->Owner == 0) {
			read.Mate2.Bases     = ite->Query;
			read.Mate2.Qualities = ite->BaseQualities;
			if ( read.Mate2.Bases.Length() > 0 ) read.Mate2.Bases.Remove('-');
			if ( _isSolid ) {
				mate2Cs = ite->CsQuery;
				mate2Cq = ite->CsBaseQualities;
			}
		}
	}

	const bool isMate1Unique   = ( newMate1Set.size() == 1 ) ? true : false;
	const bool isMate2Unique   = ( newMate2Set.size() == 1 ) ? true : false;
	const bool isMate1Multiple = ( newMate1Set.size() > 1 ) ? true : false;
	const bool isMate2Multiple = ( newMate2Set.size() > 1 ) ? true : false;
	bool isMate1Empty    = ( newMate1Set.size() == 0 ) ? true : false;
	bool isMate2Empty    = ( newMate2Set.size() == 0 ) ? true : false;

/*
if (r.Name == "11_67645641_67646650_1:0:0_4:0:0_3e62") {
  cerr << "mate1\t" << mate1SwScore << "\t" << mate1NextSwScore << endl;
  for ( vector<Alignment>::iterator ite = r.Mate1Alignments.begin(); ite != r.Mate1Alignments.end(); ++ite ) {
    cerr << (*ite).ReferenceIndex << "\t" 
         << (*ite).ReferenceBegin << "\t" 
	 << (*ite).SwScore << "\t" 
	 << (*ite).NextSwScore << "\t"
	 << (*ite).NumHash << "\t" 
	 << (*ite).NumMapped << endl;
  }
  cerr << "mate2\t" << mate2SwScore << "\t" << mate2NextSwScore << endl;
  for ( vector<Alignment>::iterator ite = r.Mate2Alignments.begin(); ite != r.Mate2Alignments.end(); ++ite ) {
    cerr << (*ite).ReferenceIndex << "\t" 
         << (*ite).ReferenceBegin << "\t" 
	 << (*ite).SwScore << "\t" 
	 << (*ite).NextSwScore << "\t"
	 << (*ite).NumHash << "\t" 
	 << (*ite).NumMapped << endl;
  }
}
*/
	// UU, UM, and MM pair
	if ( (isMate1Unique && isMate2Unique)
		|| (isMate1Unique && isMate2Multiple)
		|| (isMate1Multiple && isMate2Unique)
		|| (isMate1Multiple && isMate2Multiple) ) {

		Alignment al1 = *newMate1Set[0], al2 = *newMate2Set[0];
		if ( ( isMate1Unique && isMate2Multiple )
	          || ( isMate1Multiple && isMate2Unique )
	          || ( isMate1Multiple && isMate2Multiple ) )
			BestNSecondBestSelection::Select( al1, al2, newMate1Set, newMate2Set, _expectedFragmentLength, 
			    _sequencingTechnologies, read.Mate1.Bases.Length(), read.Mate2.Bases.Length(), true, true, false );

		// TODO: handle fragment length for others sequencing techs
		int minFl = _expectedFragmentLength - _localAlignmentSearchRadius;
		int maxFl = _expectedFragmentLength + _localAlignmentSearchRadius;
		bool properPair1 = false, properPair2 = false;
		al1.IsFirstMate = true;
		al2.IsFirstMate = false;

		// MM pair is always an improper pair
		//if ( !isMate1Multiple || !isMate2Multiple ) {
			properPair1 = al1.SetPairFlagsAndFragmentLength( al2, minFl, maxFl, _sequencingTechnologies );
			properPair2 = al2.SetPairFlagsAndFragmentLength( al1, minFl, maxFl, _sequencingTechnologies );
		//}

		if ( properPair1 != properPair2 ) {
			cout << "ERROR: An inconsistent proper pair is found." << endl;
			exit(1);
		}
		
		SetAlignmentFlags( al1, al2, true, properPair1, true, _isPairedEnd, true, true, r );
		SetAlignmentFlags( al2, al1, true, properPair2, false, _isPairedEnd, true, true, r );

		al1.CsQuery         = mate1Cs;
		al2.CsQuery         = mate2Cs;
		al1.CsBaseQualities = mate1Cq;
		al2.CsBaseQualities = mate2Cq;
		
		al1.NumMapped = nMate1Alignments;
		al2.NumMapped = nMate2Alignments;
		al1.NumHash   = nMate1Hashes;
		al2.NumHash   = nMate2Hashes;
		al1.QueryLength = al1.Query.Length();
		al2.QueryLength = al2.Query.Length();
		if (al1.SwScore == mate1SwScore) 
			al1.NextSwScore = (al1.NextSwScore > mate1NextSwScore) ? al1.NextSwScore: mate1NextSwScore;
		else 
			al1.NextSwScore = mate1SwScore;  // implies al1.SwScore < mate1SwScore
		if (al2.SwScore == mate2SwScore) 
			al2.NextSwScore = (al2.NextSwScore > mate2NextSwScore) ? al2.NextSwScore: mate2NextSwScore;
		else 
			al2.NextSwScore = mate2SwScore;  // implies al2.SwScore < mate2SwScore
		
		al1.Entropy = _entropy.shannon_H(al1.Query.Data(), al1.QueryLength);
		al2.Entropy = _entropy.shannon_H(al2.Query.Data(), al2.QueryLength);

		al1.RecalibratedQuality = GetMappingQuality(al1, al1.QueryLength, al2, al2.QueryLength);
		al2.RecalibratedQuality = GetMappingQuality(al1, al1.QueryLength, al2, al2.QueryLength);

		CZaTager za1, za2;
		const char* zaTag1 = za1.GetZaTag( al1, al2, true );
		const char* zaTag2 = za2.GetZaTag( al2, al1, false );

		//ostringstream zaTag1Stream, zaTag2Stream;
		//zaTag1Stream << "" << al1.SwScore << ";" << al1.NextSwScore << ";" << al1.NumLongestMatchs << ";" << al1.Entropy << ";" << al1.NumMapped << ";" << al1.NumHash;
		//zaTag2Stream << "" << al2.SwScore << ";" << al2.NextSwScore << ";" << al2.NumLongestMatchs << ";" << al2.Entropy << ";" << al2.NumMapped << ";" << al2.NumHash;

		//const char* zaTag1 = zaTag1Stream.str().c_str();
		//const char* zaTag2 = zaTag2Stream.str().c_str();

		if ( isMate2Special ) {
			Alignment genomicAl = al1;
			Alignment specialAl = _specialAl.Mate2Alignments[0];
			SetAlignmentFlags( specialAl, genomicAl, true, false, false, _isPairedEnd, true, true, r );

			const char* zas1Tag = za1.GetZaTag( genomicAl, specialAl, true );
			const char* zas2Tag = za2.GetZaTag( specialAl, genomicAl, false );

			specialAl.CsQuery         = mate2Cs;
			specialAl.CsBaseQualities = mate2Cq;

			_sBam.SaveAlignment( genomicAl, zas1Tag, false, false, _isSolid );
			_sBam.SaveAlignment( specialAl, zas2Tag, false, false, _isSolid );
		}

		if ( isMate1Special ) {
			Alignment genomicAl = al2;
			Alignment specialAl = _specialAl.Mate1Alignments[0];
			SetAlignmentFlags( specialAl, genomicAl, true, false, true, _isPairedEnd, true, true, r );

			//CZaTager zas1, zas2;

			const char* zas1Tag = za1.GetZaTag( genomicAl, specialAl, false );
			const char* zas2Tag = za2.GetZaTag( specialAl, genomicAl, true );

			specialAl.CsQuery         = mate1Cs;
			specialAl.CsBaseQualities = mate1Cq;

			_sBam.SaveAlignment( genomicAl, zas1Tag, false, false, _isSolid );
			_sBam.SaveAlignment( specialAl, zas2Tag, false, false, _isSolid );
		}

		// GetStatisticsCounters needs some information
		r.Mate1Alignments[0] = al1;
		r.Mate2Alignments[0] = al2;

		//_rBam.SaveAlignment( al1, zaTag1, false, false, _isSolid );
		//_rBam.SaveAlignment( al2, zaTag2, false, false, _isSolid );
		_rBam.SaveAlignment( al1, zaTag1, false, false, _isSolid );
		_rBam.SaveAlignment( al2, zaTag2, false, false, _isSolid );

		//if ( ( _statMappingQuality <= al1.Quality ) && ( _statMappingQuality <= al2.Quality ) )
			_statisticsMaps.SaveRecord( al1, al2, _isPairedEnd, _sequencingTechnologies );


	// UX and MX pair
	} else if ( (isMate1Empty || isMate2Empty)
		&& !(isMate1Empty && isMate2Empty) ) {
		
		Alignment al1, al2, unmappedAl;
		if (!isMate1Empty) al1 = *newMate1Set[0];
		if (!isMate2Empty) al2 = *newMate2Set[0];
		if ( isMate1Multiple || isMate2Multiple ) {
			if (r.Name == "9_67107412_67108404_1:0:0_0:0:0_266e") cerr << "bestselection" << endl;
			BestNSecondBestSelection::Select(al1, al2, newMate1Set, newMate2Set, _expectedFragmentLength, 
			    _sequencingTechnologies, read.Mate1.Bases.Length(), read.Mate2.Bases.Length(), 
			    (isMate1Empty ? false : true), (isMate2Empty ? false : true), false);
		}

		//isMate1Empty = r.Mate1Alignments.empty();
		//isMate2Empty = r.Mate2Alignments.empty();
		
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
		Alignment al         = isFirstMate ? al1 : al2;
		al.CsQuery           = isFirstMate ? mate1Cs : mate2Cs;
		al.CsBaseQualities   = isFirstMate ? mate1Cq : mate2Cq;

		//Alignment unmappedAl = !isFirstMate ? r.Mate1Alignments[0] : r.Mate2Alignments[0];

		if ( _isPairedEnd ) {
			unmappedAl = !isFirstMate ? r.Mate1Alignments[0]: r.Mate2Alignments[0];
			unmappedAl.Query           = isFirstMate ? read.Mate1.Bases : read.Mate2.Bases;
			unmappedAl.BaseQualities   = isFirstMate ? read.Mate1.Qualities : read.Mate2.Qualities;
			unmappedAl.ReferenceIndex  = al.ReferenceIndex;
			unmappedAl.ReferenceBegin  = al.ReferenceBegin;
			unmappedAl.CsQuery         = !isFirstMate ? mate1Cs : mate2Cs;
			unmappedAl.CsBaseQualities = !isFirstMate ? mate1Cq : mate2Cq;
		}

		SetAlignmentFlags( al, unmappedAl, false, false, isFirstMate, _isPairedEnd, true, false, r );
		al.NumMapped = isFirstMate ? nMate1Alignments : nMate2Alignments;
		al.NumHash   = isFirstMate ? nMate1Hashes : nMate2Hashes;
		al.QueryLength = al.Query.Length();
		al.Entropy = _entropy.shannon_H(al.Query.Data(), al.QueryLength);
		if (isFirstMate) {
			if (al.SwScore == mate1SwScore) 
				al.NextSwScore = (al.NextSwScore > mate1NextSwScore) ? al.NextSwScore: mate1NextSwScore;
			else 
				al.NextSwScore = mate1SwScore;
		} else {
			if (al.SwScore == mate2SwScore) 
				al.NextSwScore = (al.NextSwScore > mate2NextSwScore) ? al.NextSwScore: mate2NextSwScore;
			else 
				al.NextSwScore = mate2SwScore;
		}

		al.RecalibratedQuality = GetMappingQuality(al, al.QueryLength);

		SetAlignmentFlags( unmappedAl, al, true, false, !isFirstMate, _isPairedEnd, false, true, r );
		unmappedAl.NumMapped = 0;
		unmappedAl.NumHash   = 0;
		
		// show the original MQs in ZAs, and zeros in MQs fields of a BAM
		const char* zaTag1 = za1.GetZaTag( al, unmappedAl, isFirstMate, !_isPairedEnd, true );
		const char* zaTag2 = za2.GetZaTag( unmappedAl, al, !isFirstMate, !_isPairedEnd, false );

		/*
		ostringstream zaTag1Stream, zaTag2Stream;
		zaTag1Stream << "" << al.SwScore << ";" 
		             << al.NextSwScore << ";" 
			     << al.NumLongestMatchs << ";" 
			     << al.Entropy << ";" 
			     << al.NumMapped << ";" 
			     << al.NumHash;
		zaTag2Stream << "" << unmappedAl.SwScore << ";" 
		             << unmappedAl.NextSwScore << ";" 
			     << unmappedAl.NumLongestMatchs << ";" 
			     << unmappedAl.Entropy << ";" 
			     << unmappedAl.NumMapped << ";" 
			     << unmappedAl.NumHash;
		*/

		// store the alignment
		_rBam.SaveAlignment( al, zaTag1, false, false, _isSolid );

		// store mate1 special hits
		// NOTE: we consider mate2 special hits in the next block
		if ( isMate1Special ) {
			Alignment specialAl = _specialAl.Mate1Alignments[0];
			SetAlignmentFlags( specialAl, al, !isFirstMate, false, true, _isPairedEnd, true, !isFirstMate, r );

			//CZaTager zas1, zas2;

			const char* zas2Tag = isFirstMate 
				? za2.GetZaTag( specialAl, al, true, !_isPairedEnd, true ) 
				: za2.GetZaTag( specialAl, al, true );

			specialAl.CsQuery         = mate1Cs;
			specialAl.CsBaseQualities = mate1Cq;

			_sBam.SaveAlignment( specialAl, zas2Tag, false, false, _isSolid );
		}

		
		if ( _isPairedEnd ) {
			// store mate2 alignment in regular and unmapped bams
			_rBam.SaveAlignment( unmappedAl, zaTag2, true, false, _isSolid );
			//_uBam.SaveAlignment( unmappedAl, 0, true, false, _isSolid );

			// store special hits
			if ( isMate2Special ) {
				Alignment specialAl = _specialAl.Mate2Alignments[0];
				SetAlignmentFlags( specialAl, al, isFirstMate, false, false, _isPairedEnd, true, isFirstMate, r );

				//CZaTager zas1, zas2;

				const char* zas2Tag = !isFirstMate 
					? za2.GetZaTag( specialAl, al, false, !_isPairedEnd, true ) 
					: za2.GetZaTag( specialAl, al, false );

				specialAl.CsQuery         = mate2Cs;
				specialAl.CsBaseQualities = mate2Cq;

				_sBam.SaveAlignment( specialAl, zas2Tag, false, false, _isSolid );
			}
		}

		unmappedAl.IsFilteredOut = isFirstMate ? isMate1FilteredOut : isMate2FilteredOut;
		// GetStatisticsCounters needs some information
		r.Mate1Alignments[0] = isFirstMate ? al : unmappedAl;
		if ( _isPairedEnd )
			r.Mate2Alignments[0] = isFirstMate ? unmappedAl : al;
			
		//if ( _statMappingQuality <= al.Quality ) 
			_statisticsMaps.SaveRecord( ( isFirstMate ? al : unmappedAl ), ( !isFirstMate ? al : unmappedAl ), _isPairedEnd, _sequencingTechnologies );

	
	// XX
	} else if ( isMate1Empty && isMate2Empty ) {
		

		Alignment unmappedAl1, unmappedAl2;

		unmappedAl1 = r.Mate1Alignments[0];
		SetAlignmentFlags( unmappedAl1, unmappedAl2, false, false, true, _isPairedEnd, false, false, r );
		unmappedAl1.NumMapped = nMate1Alignments;
		
		unmappedAl1.Query           = read.Mate1.Bases;
		unmappedAl1.BaseQualities   = read.Mate1.Qualities;
		unmappedAl1.CsQuery         = mate1Cs;
		unmappedAl1.CsBaseQualities = mate1Cq;

		//_uBam.SaveAlignment( unmappedAl1, 0, true, false, _isSolid );
		_rBam.SaveAlignment( unmappedAl1, 0, true, false, _isSolid );
		

		// store special hits
		if ( isMate1Special ) {
			Alignment specialAl = _specialAl.Mate1Alignments[0];
			SetAlignmentFlags( specialAl, unmappedAl2, false, false, true, _isPairedEnd, true, false, r );

			//CZaTager zas1, zas2;

			const char* zas2Tag = za2.GetZaTag( specialAl, unmappedAl2, true, !_isPairedEnd, true );

			specialAl.CsQuery         = mate1Cs;
			specialAl.CsBaseQualities = mate1Cq;

			_sBam.SaveAlignment( specialAl, zas2Tag, false, false, _isSolid );
		}


		if ( _isPairedEnd ) {
			unmappedAl2 = r.Mate2Alignments[0];
			SetAlignmentFlags( unmappedAl2, unmappedAl1, false, false, false, _isPairedEnd, false, false, r );
			unmappedAl2.NumMapped = nMate2Alignments;
			
			unmappedAl2.Query           = read.Mate2.Bases;
			unmappedAl2.BaseQualities   = read.Mate2.Qualities;
			unmappedAl2.CsQuery         = mate2Cs;
			unmappedAl2.CsBaseQualities = mate2Cq;

			//_uBam.SaveAlignment( unmappedAl2, 0, true, false, _isSolid );
			_rBam.SaveAlignment( unmappedAl2, 0, true, false, _isSolid );

			// store special hits
			if ( isMate2Special ) {
				Alignment specialAl = _specialAl.Mate2Alignments[0];
				SetAlignmentFlags( specialAl, unmappedAl1, false, false, false, _isPairedEnd, true, false, r );

				//CZaTager zas1, zas2;

				const char* zas2Tag = za2.GetZaTag( specialAl, unmappedAl1, false, !_isPairedEnd, true );

				specialAl.CsQuery         = mate2Cs;
				specialAl.CsBaseQualities = mate2Cq;

				_sBam.SaveAlignment( specialAl, zas2Tag, false, false, _isSolid );
			}
		}

		// GetStatisticsCounters need some information
		unmappedAl1.IsFilteredOut = isMate1FilteredOut;
		unmappedAl2.IsFilteredOut = isMate2FilteredOut;

		r.Mate1Alignments[0] = unmappedAl1;
		if ( _isPairedEnd )
			r.Mate2Alignments[0] = unmappedAl2;

		_statisticsMaps.SaveRecord( unmappedAl1, unmappedAl2, _isPairedEnd, _sequencingTechnologies );
	
	} else {
		cout << "ERROR: Unknown pairs." << endl;
		cout << r.Mate1Alignments.size() << "\t" << r.Mate2Alignments.size() << endl;
		exit(1);
	}

}

unsigned char CArchiveMerge::GetMappingQuality (const Alignment& al, 
                                                const int& al_length) {

	QualityNeuralNetwork::FannInputs mate1Ann;

	mate1Ann.read_length   = al_length;
	mate1Ann.swScore       = al.SwScore;
	mate1Ann.nextSwScore   = al.NextSwScore;
	mate1Ann.longest_match = al.NumLongestMatchs;
	mate1Ann.entropy       = al.Entropy;
	mate1Ann.numMappings   = al.NumMapped;
	mate1Ann.numHashes     = al.NumHash;
	
	return _mqCalculator.GetQualitySe(mate1Ann);
}

unsigned char CArchiveMerge::GetMappingQuality (const Alignment& al1, 
                                                const int& al1_length, 
					        const Alignment& al2,
						const int& al2_length) {
 	//int fl = (al1.ReferenceBegin > al2.ReferenceBegin) ? (al1.ReferenceBegin - al2.ReferenceBegin) : (al2.ReferenceBegin - al1.ReferenceBegin);
	//fl += al1.Query.Length();
	//fl = abs(mSettings.MedianFragmentLength - fl);
	int flDiff = _expectedFragmentLength - abs(al1.FragmentLength);
	//if (al1.ReferenceIndex != al2.ReferenceIndex)
	//	flDiff = INT_MAX - 1;
	//else
		flDiff = abs(flDiff);

	QualityNeuralNetwork::FannInputs mate1Ann;
	QualityNeuralNetwork::FannInputs mate2Ann;
	
	mate1Ann.read_length   = al1_length;
	mate1Ann.swScore       = al1.SwScore;
	mate1Ann.nextSwScore   = al1.NextSwScore;
	mate1Ann.longest_match = al1.NumLongestMatchs;
	mate1Ann.entropy       = al1.Entropy;
	mate1Ann.numMappings   = al1.NumMapped;
	mate1Ann.numHashes     = al1.NumHash;

	mate2Ann.read_length   = al2_length;
	mate2Ann.swScore       = al2.SwScore;
	mate2Ann.nextSwScore   = al2.NextSwScore;
	mate2Ann.longest_match = al2.NumLongestMatchs;
	mate2Ann.entropy       = al2.Entropy;
	mate2Ann.numMappings   = al2.NumMapped;
	mate2Ann.numHashes     = al2.NumHash;
	return _mqCalculator.GetQualityPe(mate1Ann, mate2Ann, flDiff);
}

void CArchiveMerge::SetAlignmentFlags( 
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

	// GetFragmentAlignmentQuality
	/*
	if ( isProperPair ) {
		const bool isUU = ( al.NumMapped == 1 ) && ( mate.NumMapped == 1 );
		const bool isMM = ( al.NumMapped > 1 ) && ( mate.NumMapped > 1 );

		int aq = al.Quality;
		if ( isUU )      aq = (int) ( UU_COEFFICIENT * aq + UU_INTERCEPT );
                else if ( isMM ) aq = (int) ( MM_COEFFICIENT * aq + MM_INTERCEPT );
                else             aq = (int) ( UM_COEFFICIENT * aq + UM_INTERCEPT );

                if(aq < 0)       al.Quality = 0;
                else if(aq > 99) al.Quality = 99;
                else             al.Quality = aq;
	}
	*/

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
	
	if ( _hasSpecial ) {
		// the last one is special archive
		_specialArchiveName = *_inputFilenames.rbegin();
		vector < string >::iterator ite = _inputFilenames.end() - 1;
		_inputFilenames.erase( ite );
		_specialArchiveEmpty = false;
		_specialReader.Open( _specialArchiveName );

		if ( !_specialArchiveEmpty ) {
			_specialArchiveEmpty = !_specialReader.LoadNextRead( _specialAl );
			if ( !_specialArchiveEmpty ) *_readNo = *_readNo + 1;
		}
	}

	unsigned int nTemp = _inputFilenames.size();

	// sanity check
	if ( nTemp == 0 ) {
		printf("ERROR: There is no temporary archive.\n");
		exit(1);
	}
	
	// initialize MOSAIK readers for all temp files
	vector< MosaikReadFormat::CAlignmentReader* > readers;
	SortNMergeUtilities::OpenMosaikReader( readers, _inputFilenames );
	
	Mosaik::AlignedRead mr;
	unsigned int nDone = 0;
	vector< bool > done(nTemp, false);
	vector< SortNMergeUtilities::AlignedReadPair > reads(nTemp);

	// prepare BAM writers
	_sBam.Open( _outputFilename + ".special.bam", _sHeader );
	_rBam.Open( _outputFilename + ".bam", _rHeader );

	
	// pick the min one
	vector< SortNMergeUtilities::AlignedReadPair >::iterator ite;
	Mosaik::AlignedRead minRead;


	bool isDone = false;
	while ( !isDone ) {
		//SortNMergeUtilities::FindMinElement(reads, ite);
		
		minRead.Clear();

		for ( unsigned int i = 0; i < nTemp; i++ ) {
			if ( !SortNMergeUtilities::LoadNextReadPair(readers[i], i, reads) ) {
				done[i] = true;
				nDone++;
				isDone = true;
			} else {
				UpdateReferenceIndex( reads[i].read, i );
				*_readNo = *_readNo + 1;

				if ( i == 0 )
					minRead = reads[i].read;
				else {
					minRead.Mate1Alignments.push_back( reads[i].read.Mate1Alignments[0] );
					if ( _isPairedEnd )
					minRead.Mate2Alignments.push_back( reads[i].read.Mate2Alignments[0] );
				}
			}
		}

		// sanity check
		if ( isDone && ( nDone != nTemp ) ) {
			cout << "ERROR: The number of reads in archives do not match." << endl;
			exit(1);
		}

		if ( !isDone ) {
			WriteAlignment( minRead );
			CalculateStatisticsCounters( minRead );
		}

	}	
	
	
	// Close BAMs
	_sBam.Close();
	//_uBam.Close();
	_rBam.Close();

	
	// close readers
	SortNMergeUtilities::CloseMosaikReader( readers );

	if ( _hasSpecial )
		_specialReader.Close();
}
