// ***************************************************************************
// CMosaikText - exports alignments to various file formats.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "MosaikText.h"
        const char SORTED_POSITION[ ] = "coordinate";
	        const char SORTED_READNAME[ ] = "queryname";
		        const char SORTED_UNSORTED[ ] = "unsorted";


// constructor
CMosaikText::CMosaikText(void) {}

// destructor
CMosaikText::~CMosaikText(void) {}

// set sorting order
void CMosaikText::SetSortingOrder ( const unsigned short sortingModel ) {
	switch ( sortingModel ) {
		case 0: 
			mFlags.IsSortingByPosition  = true;
			mSettings.SortingModel      = SORTED_POSITION;
			break;
		case 1: 
			mFlags.IsSortingByPosition  = false;
			mSettings.SortingModel      = SORTED_READNAME;
			break;
		case 2: 
			mFlags.IsSortingByPosition  = false;
			mSettings.SortingModel      = SORTED_UNSORTED;
			break;
		default:
			cout << "The sorting model is invalid." << endl;
			exit(1);
	}
}

// input fastq file
void CMosaikText::ParseFastqFile ( const string& filename ) {
	mFlags.EnableFastqPatching = true;
	mSettings.inputFastqFilename = filename;
}

// input 2nd mate fastq file
void CMosaikText::ParseFastq2File ( const string& filename ) {
	mSettings.inputFastq2Filename = filename;
}

// enables AXT output
void CMosaikText::EnableAxtOutput(const string& filename) { 
	mFlags.IsAxtEnabled   = true;
	mSettings.AxtFilename = filename;
}

// enables BAM output
void CMosaikText::EnableBamOutput(const string& filename) {
	mFlags.IsBamEnabled   = true;
	mSettings.BamFilename = filename;
}

// enables BED output
void CMosaikText::EnableBedOutput(const string& filename) { 
	mFlags.IsBedEnabled   = true;
	mSettings.BedFilename = filename;
}

// enables Eland output
void CMosaikText::EnableElandOutput(const string& filename) { 
	mFlags.IsElandEnabled   = true;
	mSettings.ElandFilename = filename;
}

// enables FASTQ output
void CMosaikText::EnableFastqOutput(const string& filename) {
	mFlags.IsFastqEnabled   = true;
	mSettings.FastqFilename = filename;
}

// enables Psl output
void CMosaikText::EnablePslOutput(const string& filename) { 
	mFlags.IsPslEnabled   = true;
	mSettings.PslFilename = filename;
}

// enables the reference sequence filter
void CMosaikText::EnableReferenceFilter(const string& referenceName, const string& alignmentFilename) {
	mFlags.UseReferenceFilter = true;

	MosaikReadFormat::CAlignmentReader reader;
	reader.Open(alignmentFilename);

	// find the reference ID associated with the reference name
	bool foundReferenceIndex = false;

	unsigned int currentReferenceIndex = 0;
	uint64_t numAlignments             = 0;

	vector<ReferenceSequence>* pRefSeqs = reader.GetReferenceSequences();
	for(vector<ReferenceSequence>::const_iterator rsIter = pRefSeqs->begin(); rsIter != pRefSeqs->end(); ++rsIter, ++currentReferenceIndex) {
		if(rsIter->Name == referenceName) {
			foundReferenceIndex = true;
			numAlignments = rsIter->NumAligned;
			break;
		}		
	}

	reader.Close();

	// make sure the reference exists in the alignment archive
	if(!foundReferenceIndex) {
		printf("ERROR: The chosen reference sequence (%s) was not found in the alignment archive.\n", referenceName.c_str());
		exit(1);
	}

	// make sure we have some alignments to filter
	if(numAlignments == 0) {
		printf("ERROR: There were no alignments assigned to the chosen reference sequence (%s).\n", referenceName.c_str());
		exit(1);
	}

	mSettings.FilteredReferenceIndex    = currentReferenceIndex;
	mSettings.NumFilteredReferenceReads = numAlignments;
}

// enables SAM output
void CMosaikText::EnableSamOutput(const string& filename) {
	mFlags.IsSamEnabled   = true;
	mSettings.SamFilename = filename;

	const unsigned int filenameLength = filename.size();
	if((filenameLength >= 3) && (filename.substr(filenameLength - 3) != ".gz")) mSettings.SamFilename += ".gz";
}

// enables screen output
void CMosaikText::EnableScreenOutput(void) {
	mFlags.IsScreenEnabled = true;
}

// when triggered, the coverage calculation will only include unique reads
void CMosaikText::EvaluateUniqueReadsOnly(void) { 
	cout << "- evaluating unique reads only." << endl;
	mFlags.EvaluateUniqueReadsOnly = true;
}

// opens the output file stream for the AXT file
void CMosaikText::InitializeAxt(void) {
	mStreams.axt = fopen(mSettings.AxtFilename.c_str(), "wb");

	if(!mStreams.axt) {
		printf("ERROR: Unable to open the AXT file for writing.\n");
		exit(1);
	}
}

// opens the output file stream for the BAM file
void CMosaikText::InitializeBam(vector<ReferenceSequence>* pRefSeqs, vector<MosaikReadFormat::ReadGroup>& readGroups) {

	// set the sort order
	BamHeader header;
	//header.SortOrder = (isSortedByPosition ? SORTORDER_POSITION : SORTORDER_UNSORTED);
	if ( mSettings.SortingModel == SORTED_POSITION )
		header.SortOrder = SORTORDER_POSITION;
	else if ( mSettings.SortingModel == SORTED_READNAME )
		header.SortOrder = SORTORDER_READNAME;
	else
		header.SortOrder = SORTORDER_UNSORTED;

	// set the reference sequences and read groups
	header.pReferenceSequences = pRefSeqs;
	header.pReadGroups         = &readGroups;

	// open the bam file
	mStreams.bam.Open(mSettings.BamFilename, header);
}

// opens the output file stream for the BED file
void CMosaikText::InitializeBed(void) {
	mStreams.bed = fopen(mSettings.BedFilename.c_str(), "wb");

	if(!mStreams.bed) {
		printf("ERROR: Unable to open the BED file for writing.\n");
		exit(1);
	}
}

// opens the output file stream for the Eland file
void CMosaikText::InitializeEland(void) {
	mStreams.eland = fopen(mSettings.ElandFilename.c_str(), "wb");

	if(!mStreams.eland) {
		printf("ERROR: Unable to open the Eland file for writing.\n");
		exit(1);
	}
}

// opens the output file stream for the SAM file
void CMosaikText::InitializeSam(vector<ReferenceSequence>* pRefSeqs, vector<MosaikReadFormat::ReadGroup>& readGroups) {
	mStreams.sam = gzopen(mSettings.SamFilename.c_str(), "wb");

	if(!mStreams.sam) {
		printf("ERROR: Unable to open the SAM file for writing.\n");
		exit(1);
	}

	// store the header
	gzprintf(mStreams.sam, "@HD\tVN:1.0\tSO:");
	//if(isSortedByPosition) gzprintf(mStreams.sam, "coordinate\n");
	//else gzprintf(mStreams.sam, "unsorted\n");
	gzprintf( mStreams.sam, mSettings.SortingModel.c_str() );
	gzprintf( mStreams.sam, "\n");


	// store the sequence dictionary
	vector<ReferenceSequence>::const_iterator rsIter;
	for(rsIter = pRefSeqs->begin(); rsIter != pRefSeqs->end(); ++rsIter) {
		gzprintf(mStreams.sam, "@SQ\tSN:%s\tLN:%u", rsIter->Name.c_str(), rsIter->NumBases);
		if(!rsIter->GenomeAssemblyID.empty()) gzprintf(mStreams.sam, "\tAS:%s", rsIter->GenomeAssemblyID.c_str());
		if(!rsIter->MD5.empty())              gzprintf(mStreams.sam, "\tM5:%s", rsIter->MD5.c_str());
		if(!rsIter->URI.empty())              gzprintf(mStreams.sam, "\tUR:%s", rsIter->URI.c_str());
		if(!rsIter->Species.empty())          gzprintf(mStreams.sam, "\tSP:%s", rsIter->Species.c_str());
		gzprintf(mStreams.sam, "\n");
	}

	// store the read groups
	vector<MosaikReadFormat::ReadGroup>::const_iterator rgIter;
	for(rgIter = readGroups.begin(); rgIter != readGroups.end(); ++rgIter) {
		gzprintf(mStreams.sam, "@RG\tID:%s\tSM:%s", rgIter->ReadGroupID.c_str(), rgIter->SampleName.c_str());
		if(!rgIter->LibraryName.empty())      gzprintf(mStreams.sam, "\tLB:%s", rgIter->LibraryName.c_str());
		if(!rgIter->Description.empty())      gzprintf(mStreams.sam, "\tDS:%s", rgIter->Description.c_str());
		if(!rgIter->PlatformUnit.empty())     gzprintf(mStreams.sam, "\tPU:%s", rgIter->PlatformUnit.c_str());
		if(rgIter->MedianFragmentLength != 0) gzprintf(mStreams.sam, "\tPI:%u", rgIter->MedianFragmentLength);
		if(!rgIter->CenterName.empty())       gzprintf(mStreams.sam, "\tCN:%s", rgIter->CenterName.c_str());

		switch(rgIter->SequencingTechnology) {
			case ST_454:
				gzprintf(mStreams.sam, "\tPL:454\n");
				break;
			case ST_HELICOS:
				gzprintf(mStreams.sam, "\tPL:helicos\n");
				break;
			case ST_ILLUMINA:
				gzprintf(mStreams.sam, "\tPL:illumina\n");
				break;
			case ST_PACIFIC_BIOSCIENCES:
				gzprintf(mStreams.sam, "\tPL:pacific biosciences\n");
				break;
			case ST_SOLID:
				gzprintf(mStreams.sam, "\tPL:solid\n");
				break;
			case ST_SANGER:
				gzprintf(mStreams.sam, "\tPL:sanger\n");
				break;
			default:
				gzprintf(mStreams.sam, "\tPL:unknown\n");
				break;
		}
	}
}

// the lessthan operator used in MosaikText
inline bool CMosaikText::PositionLessThan( const Mosaik::AlignedRead& ar1, const Mosaik::AlignedRead& ar2 ) {
	
	unsigned int nMate1Ar1 = ar1.Mate1Alignments.size();
	unsigned int nMate1Ar2 = ar2.Mate1Alignments.size();

	if ( nMate1Ar1 == 0 ) return true;
	if ( nMate1Ar2 == 0 ) return false;

	unsigned int refIndexAr1 = ar1.Mate1Alignments.begin()->ReferenceIndex;
	unsigned int refIndexAr2 = ar2.Mate1Alignments.begin()->ReferenceIndex;
	unsigned int refBeginAr1 = ar1.Mate1Alignments.begin()->ReferenceBegin;
	unsigned int refBeginAr2 = ar2.Mate1Alignments.begin()->ReferenceBegin;

	if ( refIndexAr1 == refIndexAr2 )
		return refBeginAr1 < refBeginAr2;

	return refIndexAr1 < refIndexAr2;

}

// set the settings of input MOSAIK archive
void CMosaikText::SetArchiveSetting( const string& alignmentFilename ) {
	MosaikReadFormat::CAlignmentReader reader;
	reader.Open(alignmentFilename);
	// retrieve the read groups
	reader.GetReadGroups(mArchiveSetting.readGroups);
	// retrieve the reference sequence vector
	reader.GetReferenceSequences( mArchiveSetting.pReferenceSequences );
	// retrieve the alignment status (clear the sorting status)
	mArchiveSetting.as = reader.GetStatus();
	// retrieve the archive signature
	reader.GetSignature( mArchiveSetting.signature );
	reader.Close();
}

// sort FASTQ by read names
void CMosaikText::SortFastqByName( const string& inputFastqFilename, string& outputFastqFilename ) {
	CFileUtilities::GetTempFilename( outputFastqFilename );
	CFastq fastqReader;
	fastqReader.Open( inputFastqFilename );
	fastqReader.SortByName( outputFastqFilename );
	fastqReader.Close();
}

// initialize our patching buffers
void CMosaikText::InitializePatchingBuffer ( void ) {
	_bufferSize1 = 1024;
	_bufferSize2 = 1024;
	_originalReverseBase1    = new char [ _bufferSize1 ];
	_originalReverseBase2    = new char [ _bufferSize1 ];
	_originalReverseQuality1 = new char [ _bufferSize2 ];
	_originalReverseQuality2 = new char [ _bufferSize2 ];

	// reserve buffer for soft cliping notation
	// we're using "Z" indicating soft cliping in reference
	_clipSize   = 1024;
	_clipBuffer = new char [ _clipSize ];
	memset( _clipBuffer, 'Z', _clipSize );
}

// free patching buffers
void CMosaikText::FreePatchingBuffer ( void ) {
	_bufferSize1 = 0;
	_bufferSize2 = 0;
	if ( _originalReverseBase1 )    delete [] _originalReverseBase1;
	if ( _originalReverseBase2 )    delete [] _originalReverseBase2;
	if ( _originalReverseQuality1 ) delete [] _originalReverseQuality1;
	if ( _originalReverseQuality2 ) delete [] _originalReverseQuality2;

	_originalReverseBase1    = NULL;
	_originalReverseBase2    = NULL;
	_originalReverseQuality1 = NULL;
	_originalReverseQuality2 = NULL;


	_clipSize = 0;
	if ( _clipBuffer ) delete [] _clipBuffer;
	_clipBuffer = NULL;
}

// given a read name, search it in FASTQs
void CMosaikText::SearchReadInFastq ( const CMosaikString& readName, CFastq& fastqReader1, CFastq& fastqReader2, const bool hasFastq2 ) {

	// search the read
	while ( readName != _readName1 ) {
		bool loadFastq1 = false, loadFastq2 = false;
		loadFastq1 = fastqReader1.LoadNextMate( _readName1, _m1 );
		if ( hasFastq2 )
			loadFastq2 = fastqReader2.LoadNextMate( _readName2, _m2 );

		if ( !loadFastq1 || ( hasFastq2 && !loadFastq2 ) ) {
			cout << "ERROR: The targeted read (" << readName  << ") cannot be found in FASTQs." << endl;
			exit(1);
		}
	}

	// set patching buffers
	// Get Reverse Complement
	unsigned int queryLengthMate1 = _m1.Bases.Length();
	unsigned int queryLengthMate2 = 0;
	// prepare the larger buffer
	if ( queryLengthMate1 > _bufferSize1 ) {
		_bufferSize1 = queryLengthMate1;
		delete [] _originalReverseBase1;
		delete [] _originalReverseQuality1;
		_originalReverseBase1    = new char [ _bufferSize1 ];
		_originalReverseQuality1 = new char [ _bufferSize1 ];
	}
	memcpy( _originalReverseBase1, _m1.Bases.CData(), queryLengthMate1 );
	CSequenceUtilities::GetReverseComplement( _originalReverseBase1, queryLengthMate1 );

	// handle 2nd mate
	if ( hasFastq2 ) {
		queryLengthMate2 = _m2.Bases.Length();
		// prepare the larger buffer
		if ( queryLengthMate2 > _bufferSize2 ) {
			_bufferSize2 = queryLengthMate2;
			delete [] _originalReverseBase2;
			delete [] _originalReverseQuality2;
			_originalReverseBase2    = new char [ _bufferSize2 ];
			_originalReverseQuality2 = new char [ _bufferSize2 ];
		}

		memcpy( _originalReverseBase2, _m2.Bases.CData(), queryLengthMate2);
		CSequenceUtilities::GetReverseComplement( _originalReverseBase2, queryLengthMate2 );
	}

}

// given an alignedReadCache, sort them by positions and sorte them in a temp file
void CMosaikText::StoreReadCache ( CAlignedReadCache& cache ) {

	//cache.SortByPosition();
	
	// get temp filename
	string filename;
	CFileUtilities::GetTempFilename( filename );
	_tempFiles.push_back( filename );

	// prepare archive
	MosaikReadFormat::CAlignmentWriter aw;
	aw.Open(filename, mArchiveSetting.pReferenceSequences, mArchiveSetting.readGroups, mArchiveSetting.as, ALIGNER_SIGNATURE );
	Mosaik::AlignedRead sortedAr;
	
	cache.Rewind();
	while( cache.LoadNextAlignedRead( sortedAr ) ) {
		aw.SaveAlignedRead( sortedAr );
		sortedAr.Clear();
	}
	aw.Close();
	//cache.Reset();
}

// patchs trimmed infomation back from FASTQs
void CMosaikText::PatchInfo( const string& alignmentFilename, const string& inputFastqFilename, const string& inputFastq2Filename ) {
	//=================================================
	// Sort FASTQ by read names and load the first mate
	//=================================================
	string tempFilename1, tempFilename2;
	CFastq fastqReader1, fastqReader2;
	SortFastqByName( inputFastqFilename, tempFilename1 );
	fastqReader1.Open( tempFilename1 );
	fastqReader1.LoadNextMate( _readName1, _m1 );

	bool hasFastq2 = !inputFastq2Filename.empty();
	if ( hasFastq2 ) {
		SortFastqByName( inputFastq2Filename, tempFilename2 );
		fastqReader2.Open( tempFilename2 );
		fastqReader2.LoadNextMate( _readName2, _m2 );
	}

	//==========================================
	// Compare info in MOSAIK archive and FASTQs
	//==========================================
	// NOTE: the archive has been sorted by read names by MosaikSort
	// NOTE: several temp archives would be generated for storing alignments
	//       and they would be sorted by positions
	
	// open archive
	MosaikReadFormat::CAlignmentReader reader;
	reader.Open(alignmentFilename);

	// initialize buffers
	InitializePatchingBuffer();

	// retrieve AlignedReadCache
	unsigned int cacheSize = 50000;
	CAlignedReadCache cache( cacheSize );

	// read all alignments from archive and compare them with reads in FASTQs
	Mosaik::AlignedRead ar;
	while(reader.LoadNextRead(ar)) {
		// search the read in FASTQs and set FASTQs info in patching buffer
		// NOTE: _readName1 should be equal _readName2
		if ( ar.Name != _readName1 )
			SearchReadInFastq( ar.Name, fastqReader1, fastqReader2, hasFastq2 );

		bool isLongRead = false;
		for ( vector<Alignment>::iterator ite = ar.Mate1Alignments.begin(); ite != ar.Mate1Alignments.end(); ite++ ) {
			Mosaik::Mate currentMate;
			unsigned int queryLength;
			char* reverseMate;
			char* reverseQuality;
			if ( ite->IsFirstMate || !hasFastq2 ) {
				currentMate    = _m1;
				queryLength    = _m1.Bases.Length();
				reverseMate    = _originalReverseBase1;
				reverseQuality = _originalReverseQuality1;
			}
			else if ( !ite->IsFirstMate && hasFastq2 ) {
				currentMate    = _m2;
				queryLength    = _m2.Bases.Length();
				reverseMate    = _originalReverseBase2;
				reverseQuality = _originalReverseQuality2;
			}
			else {
				cout << "ERROR: The read(" << ar.Name << ") cannot be indicated as first or second mate." << endl;
				exit(1);
			}
				
			// get the length of the aligned bases
			CMosaikString query(ite->Query);
			query.Remove('-');
			unsigned int alignedLength = query.Length();
			// assign the original qualities to the alignment

			// compare qualities for MOSAIK mess
			CMosaikString tempq3 = ite->BaseQualities;
			if ( ite->IsReverseStrand )
				tempq3.Reverse();
			tempq3.Increment(33);
			string tempq2 = tempq3.Data();
			currentMate.Qualities.Increment(33);
			string tempq1 = currentMate.Qualities.Data();
			currentMate.Qualities.Decrement(33);
			size_t tempFound = tempq1.find( tempq2 );
			if ( tempFound == string::npos ) {
				cerr << "ERROR: The qualities cannot be found in the FASTQs." << endl;
				cerr << "       Read name:" << ar.Name << endl;
				cerr << "   Base qualites:" << tempq2 << endl;
				cerr << "           FASTQ:" << tempq1 << endl;
			}

			ite->BaseQualities.Copy( currentMate.Qualities.CData(), queryLength );
			ite->QueryBegin = 0;
			ite->QueryEnd   = queryLength - 1;

			if ( ite->IsReverseStrand )
				ite->BaseQualities.Reverse();
				
			// add soft chip
			if ( alignedLength != queryLength ) {
				char* originalBases;
				if ( ite->IsReverseStrand )
					originalBases = reverseMate;
					//ite->BaseQualities.Reverse();
				else
					originalBases = currentMate.Bases.Data();
					
				// find the aligned bases in the original bases
				string originalBasesStr = originalBases;
				size_t found;
				found = originalBasesStr.find(query.CData());
				
				if ( found == string::npos ) {
						cout << "ERROR: The trimmed bases cannot be found in the FASTQs." << endl;
						cout << "       Read name:" << ar.Name << endl;
						cout << "   Aligned bases:" << query.CData() << endl;
						if ( ite->IsReverseStrand )
							cout << "         Reverse: True" << endl;
						else
							cout << "         Reverse: False" << endl;
						exit(1);
				}
				else {
					// soft clip at the beginning
					if ( found != 0 ) {
						ite->Reference.Prepend( _clipBuffer, (unsigned int) found );
						ite->Query.Prepend( originalBases, (unsigned int) found );
					}

					// soft clip at the end
					unsigned int endFound = (unsigned int) found + alignedLength;
					if ( endFound < queryLength ) {
						char* suffix = originalBases + endFound;
						unsigned int suffixLength = queryLength - endFound;
						ite->Reference.Append( _clipBuffer, suffixLength );
						ite->Query.Append( suffix, suffixLength );
					}
				}
			}

			isLongRead = isLongRead || ( ite->Reference.Length() > 255 ) || ( ite->QueryEnd > 255 );
		}
		// store the aligned read to cache
		ar.IsLongRead = isLongRead;
		cache.Add( ar );

		if ( cache.isFull() ) {
			if ( mFlags.IsSortingByPosition )
				cache.SortByPosition();
			// store aligned read to a temp file
			// _tempFiles would collect all temp filenames
			StoreReadCache( cache );
			cache.Reset();
		}
	}
	
	if ( !cache.isEmpty() ) {
		if ( mFlags.IsSortingByPosition )
			cache.SortByPosition();
		// sort and store aligned read to a temp file
		// _tempFiles would collect all temp filenames
		StoreReadCache( cache );
		cache.Clear();
	}

	reader.Close();
	fastqReader1.Close();
	if ( hasFastq2 ) fastqReader2.Close();
		
	// delete temp FASTQ files
	rm( tempFilename1.c_str() );
	if ( hasFastq2 ) rm( tempFilename2.c_str() );

	// free buffers
	FreePatchingBuffer();
}

// merge a vector of partially sorted archive; sorting them by positions
void CMosaikText::MergeSortedArchive ( const vector <string>& filenames, const string& outArchiveFilename ) {
	// prepare archive readers
	unsigned int nTemp = filenames.size();
	
	vector<MosaikReadFormat::CAlignmentReader*> readers( nTemp );
	MosaikReadFormat::CAlignmentReader* readerPtr;
	for ( unsigned int i = 0; i < nTemp; i++ ) {
		readerPtr = new MosaikReadFormat::CAlignmentReader;
		readerPtr->Open( filenames[i] );
		readers[i] = readerPtr;
	}


	// prepare archive writer
	MosaikReadFormat::CAlignmentWriter aw;
	aw.Open( outArchiveFilename, mArchiveSetting.pReferenceSequences, mArchiveSetting.readGroups, mArchiveSetting.as, ALIGNER_SIGNATURE );


	// list for the top in each temp file
	list<Mosaik::AlignedRead> tops;
	unsigned int nDone = 0;
	vector<bool> dones( nTemp );
	Mosaik::AlignedRead ar;

	// load the first alignment in each temp archive
	for ( unsigned int i = 0; i < nTemp; i++ ) {
		ar.Clear();
		
		if ( readers[i]->LoadNextRead( ar ) ) {
			ar.Owner = i;
			tops.push_back( ar );
			dones[i] = false;
		}
		else {
			dones[i] = true;
			nDone++;
		}
	}

	unsigned int tempId;
	Mosaik::AlignedRead nextMin;
	bool isTempEmpty;
	// read all alignments and sort them
	while ( nDone != ( nTemp - 1 ) ) {
		tops.sort( PositionLessThan );
		tempId = tops.begin()->Owner;

		// save min
		aw.SaveAlignedRead( *tops.begin() );
		tops.pop_front();

		nextMin = *tops.begin();

		isTempEmpty = false;
		ar.Clear();
		if ( !dones[tempId] && readers[tempId]->LoadNextRead( ar ) )
			ar.Owner = tempId;
		else {
			nDone++;
			dones[tempId] = true;
			isTempEmpty = true;
			rm( filenames[tempId].c_str() );
		}


		// the file is empty
		if ( isTempEmpty ) continue;

		while ( PositionLessThan( ar, nextMin ) ) {
			aw.SaveAlignedRead( ar );
			ar.Clear();
			if ( !readers[tempId]->LoadNextRead( ar ) ) {
				nDone++;
				dones[tempId] = true;
				isTempEmpty = true;
				rm( _tempFiles[tempId].c_str() );
				break;
			}
		}

		ar.Owner = tempId;
		if ( !isTempEmpty )
			tops.push_back( ar );

	}

	// store the remaining aligned reads
	if ( tops.size() != 1 ) {
		cout << "ERROR: More than one aligned reads remain." << endl;
		exit(1);
	}

	tempId = tops.begin()->Owner;
	aw.SaveAlignedRead( *tops.begin() );
	ar.Clear();
	while ( readers[tempId]->LoadNextRead( ar ) ) {
		aw.SaveAlignedRead( ar );
		ar.Clear();
	}
	rm( filenames[tempId].c_str() );

	// close the reader and writer
	aw.Close();
	for ( vector<MosaikReadFormat::CAlignmentReader*>::iterator ite = readers.begin(); ite != readers.end(); ite++ )
		(*ite)->Close();
	readers.clear();

	_tempFiles.clear();
}

// sort aligned archive by positions
void CMosaikText::SortAlignmentByPosition( const string& inputArchive, const string& outputArchive ) {

	//=================
	//partially sorting
	//=================

	// open archive
	MosaikReadFormat::CAlignmentReader reader;
	reader.Open(inputArchive);
	
	// retrieve AlignedReadCache
	unsigned int cacheSize = 50000;
	CAlignedReadCache cache( cacheSize );
	
	_tempFiles.clear();
	// read all alignments from archive and compare them with reads in FASTQs
	Mosaik::AlignedRead ar;
	while(reader.LoadNextRead(ar)) {
		// store the aligned read to cache
		cache.Add( ar );

		if ( cache.isFull() ) {
			if ( mFlags.IsSortingByPosition )
				cache.SortByPosition();
			// sort and store aligned read to a temp file
			// _tempFiles would collect all temp filenames
			StoreReadCache( cache );
			cache.Reset();
		}
	}

	if ( !cache.isEmpty() ) {
		if ( mFlags.IsSortingByPosition )
			cache.SortByPosition();
		// sort and store aligned read to a temp file
		// _tempFiles would collect all temp filenames
		StoreReadCache( cache );
		cache.Clear();
	}

	//===========================================
	// globally sorting and merging as an archive
	//===========================================
	MergeSortedArchive( _tempFiles, outputArchive );
}

// parses the specified MOSAIK alignment file and matching anchors file
void CMosaikText::ParseMosaikAlignmentFile ( const string& alignmentFilename ) {

	// ============================
	// patch trimmed information
	// ============================
	
	string filename;
	bool isNewSorted = ( strcmp(mArchiveSetting.signature, SORT_SIGNATURE) == 0 ) ? true : false;
	if ( mFlags.EnableFastqPatching ) {
		// patch the trimmed info back and store alignments in temp files
		// note that the alignments in temp files would be sorted by positions
		PatchInfo( alignmentFilename, mSettings.inputFastqFilename, mSettings.inputFastq2Filename );

		// merge the temp archives which are generated by PatchInfo
		CFileUtilities::GetTempFilename( filename );
		MergeSortedArchive( _tempFiles, filename );
		_tempFiles.clear();
	}
	else {
		if ( isNewSorted && mFlags.IsSortingByPosition ) {
			// note that MosaikSort sorts alignments by names
			CFileUtilities::GetTempFilename( filename );
			SortAlignmentByPosition( alignmentFilename, filename );
		} else 
			filename = alignmentFilename;
	}

	
	// open the alignment archive
	MosaikReadFormat::CAlignmentReader reader;
	reader.Open( filename );

	// retrieve the alignment archive status
	//const AlignmentStatus as = reader.GetStatus();
	//const bool isSortedByPosition = ((mArchiveSetting.as & AS_SORTED_ALIGNMENT) != 0 ? true : false);

	// retrieve the reference sequences
	//vector<ReferenceSequence>* pRefSeqs = reader.GetReferenceSequences();

	// retrieve the read groups
	//vector<MosaikReadFormat::ReadGroup> readGroups;	
	//reader.GetReadGroups(readGroups);

	// retrieve the total number of reads
	uint64_t numReads = reader.GetNumReads();

	// jump to the desired reference index
	//if(mFlags.UseReferenceFilter && isSortedByPosition) {
	if(mFlags.UseReferenceFilter && mFlags.IsSortingByPosition) {
		reader.Jump(mSettings.FilteredReferenceIndex, 0);
		numReads = mSettings.NumFilteredReferenceReads;
	}

	// open our output file streams
	if(mFlags.IsAxtEnabled)   InitializeAxt();
	if(mFlags.IsBamEnabled)   InitializeBam(&mArchiveSetting.pReferenceSequences, mArchiveSetting.readGroups);
	if(mFlags.IsBedEnabled)   InitializeBed();
	if(mFlags.IsElandEnabled) InitializeEland();
	if(mFlags.IsSamEnabled)   InitializeSam(&mArchiveSetting.pReferenceSequences, mArchiveSetting.readGroups);

	// initialize
	mCurrentRead      = 0;
	mCurrentAlignment = 0;

	if(!mFlags.IsScreenEnabled) {
		CConsole::Heading(); printf("Converting alignment archive:\n"); CConsole::Reset();
		//CProgressBar<uint64_t>::StartThread(&mCurrentRead, 0, numReads, (isSortedByPosition ? "alignments" : "reads"));
		CProgressBar<uint64_t>::StartThread(&mCurrentRead, 0, numReads, "alignments");
	}

	// retrieve all reads from the alignment reader
	Mosaik::AlignedRead ar;
	while(reader.LoadNextRead(ar)) {
		
		// stop processing reads if we're already past the current reference sequence
		if(mFlags.UseReferenceFilter && mFlags.IsSortingByPosition && 
			(ar.Mate1Alignments.begin()->ReferenceIndex > mSettings.FilteredReferenceIndex)) break;

		const unsigned int numMate1Alignments = (unsigned int)ar.Mate1Alignments.size();
		const unsigned int numMate2Alignments = (unsigned int)ar.Mate2Alignments.size();

		// convert the read group code to a read group ID string
		MosaikReadFormat::ReadGroup rg = reader.GetReadGroupFromCode(ar.ReadGroupCode);
		const bool isColorspace = (rg.SequencingTechnology == ST_SOLID ? true : false);

		// dump the mate 1 alignments
		if(numMate1Alignments > 0) {
			if((mFlags.EvaluateUniqueReadsOnly && (numMate1Alignments == 1)) || !mFlags.EvaluateUniqueReadsOnly) {
				// prepare ZA tag
				if ( numMate2Alignments > 0 ) {
					bool hasLessPosition = ( ar.Mate1Alignments[0].ReferenceBegin < ar.Mate2Alignments[0].ReferenceBegin ) ? true : false;
					if ( hasLessPosition && isNewSorted ) 
						zaString = (char*) zaTager.GetZaTag( ar.Mate1Alignments, ar.Mate2Alignments );
					else 
						zaString = 0;
				}

				ProcessAlignments(1, ar.Name, isColorspace, ar.Mate1Alignments, rg.ReadGroupID);
			}
		}
		zaString = 0;

		// we shouldn't use any alignment in ar.Mate2Alignments if it is the sorted archive 
		// since the information there are just used for ZA tags
		// dump the mate 2 alignments
		if( !isNewSorted && numMate2Alignments > 0) {
			if((mFlags.EvaluateUniqueReadsOnly && (numMate2Alignments == 1)) || !mFlags.EvaluateUniqueReadsOnly) {
				ProcessAlignments(2, ar.Name, isColorspace, ar.Mate2Alignments, rg.ReadGroupID);
			}
		}

		// increment the read counter
		++mCurrentRead;
		ar.Clear();
	}
	
	if ( mFlags.EnableFastqPatching || ( isNewSorted && mFlags.IsSortingByPosition ) )
		rm(filename.c_str());

	// wait for the progress bar to finish
	if(!mFlags.IsScreenEnabled) CProgressBar<uint64_t>::WaitThread();

	// close our file streams
	reader.Close();
	if(mFlags.IsAxtEnabled)   fclose(mStreams.axt);
	if(mFlags.IsBamEnabled)   mStreams.bam.Close();
	if(mFlags.IsBedEnabled)   fclose(mStreams.bed);
	if(mFlags.IsElandEnabled) fclose(mStreams.eland);
	if(mFlags.IsSamEnabled)   gzclose(mStreams.sam);
}

// processes the alignments according to the chosen file format
void CMosaikText::ProcessAlignments(const unsigned char mateNum, const CMosaikString& readName, const bool isColorspace, vector<Alignment>& alignments, const string& readGroupID) {

	// initialize
	vector<Alignment>::iterator alIter;
	const bool isUnique = (alignments.size() == 1 ? true : false);

	// ===================================
	// first pass: collect additional data
	// ===================================

	if(mFlags.IsElandEnabled) {

		// initialize
		unsigned int mismatchCounts[3];
		uninitialized_fill(mismatchCounts, mismatchCounts + 3, 0);

		unsigned short lowestMismatchCount = MAX_SHORT;

		unsigned char bestAlignmentQuality = 0;
		vector<Alignment>::const_iterator bestAlIter = alignments.begin();

		// record mismatch counts and identify the best alignment
		for(alIter = alignments.begin(); alIter != alignments.end(); ++alIter) {

			unsigned short numMismatches = alIter->NumMismatches;

			if(numMismatches < lowestMismatchCount) lowestMismatchCount = numMismatches;
			if(numMismatches <= 2) ++mismatchCounts[numMismatches];

			if(alIter->Quality > bestAlignmentQuality) {
				bestAlignmentQuality = alIter->Quality;
				bestAlIter           = alIter;
			}
		}

		// write the best alignment
		fprintf(mStreams.eland, "%llu\t%s\t%s\t%c\t%u\t%u\t%u\t%u\t", (unsigned long long)(mCurrentRead + 1),
			readName.CData(), bestAlIter->Query.CData(), (isUnique ? 'U' : 'R'), lowestMismatchCount,
			mismatchCounts[0], mismatchCounts[1], mismatchCounts[2]);
	}

	// ===========================
	// second pass: writing output
	// ===========================

	for(alIter = alignments.begin(); alIter != alignments.end(); ++alIter) {

		// skip the alignment if filtering is activated
		// TODO: of course this is useless for the ELAND output
		if(mFlags.UseReferenceFilter && (alIter->ReferenceIndex != mSettings.FilteredReferenceIndex)) continue;

		if(mFlags.IsAxtEnabled) {
			fprintf(mStreams.axt, "%llu %s %u %u %s %u %u %c %u\n", (unsigned long long)(mCurrentAlignment + 1),
				alIter->ReferenceName, alIter->ReferenceBegin + 1, alIter->ReferenceEnd + 1, readName.CData(),
				alIter->QueryBegin + 1, alIter->QueryEnd + 1, (alIter->IsReverseStrand ? '-' : '+'), alIter->Quality);
			fprintf(mStreams.axt, "%s\n%s\n", alIter->Reference.CData(), alIter->Query.CData());
		}

		
		if(mFlags.IsBamEnabled) mStreams.bam.SaveAlignment(readName, readGroupID, alIter, zaString);

		
		// BED is in 0-based format.
		if(mFlags.IsBedEnabled) {
			fprintf(mStreams.bed, "%s %u %u %s 1 %c %u %u %s\n", alIter->ReferenceName, alIter->ReferenceBegin,
				alIter->ReferenceEnd + 1, readName.CData(), (alIter->IsReverseStrand ? '-' : '+'), 
				alIter->ReferenceBegin, alIter->ReferenceEnd + 1, (alIter->IsReverseStrand ? "0,0,255" : "255,0,0"));
		}

		if(mFlags.IsElandEnabled && isUnique) {
			fprintf(mStreams.eland, "\t%s\t%u\t%c\t..", alIter->ReferenceName, alIter->ReferenceBegin + 1,
				(alIter->IsReverseStrand ? 'R' : 'F'));

			const char* pReference = alIter->Reference.CData();
			const char* pQuery     = alIter->Query.CData();

			// TODO: update this - it probably doesn't work correctly when alignment is on reverse strand
			for(unsigned short j = 0; j < (unsigned short)alIter->Reference.Length(); ++j) {
				if(pReference[j] != pQuery[j]) fprintf(mStreams.eland, "\t%u", j + alIter->QueryBegin + 1);
			}
		}

		// SAM
		if(mFlags.IsSamEnabled) WriteSamEntry(readName, readGroupID, alIter);

		if(mFlags.IsScreenEnabled) {
			CMosaikString bq = alIter->BaseQualities;
			bq.Increment(33);

			printf("%llu %s %u %u %s %u %u %c %u %s %s%s %s%u\n", (unsigned long long)(mCurrentAlignment + 1), alIter->ReferenceName, 
				alIter->ReferenceBegin + 1, alIter->ReferenceEnd + 1, readName.CData(), alIter->QueryBegin + 1,
				alIter->QueryEnd + 1, (alIter->IsReverseStrand ? '-' : '+'), alIter->Quality, readGroupID.c_str(),
				( alIter->IsFirstMate ? "mate1" : "mate2" ), (alIter->WasRescued ? " [rescued]" : ""), (alIter->IsMapped ? " [mapped]" : "" ), alIter->NumMapped);
			printf("%s\n%s\n%s\n\n", alIter->Reference.CData(), alIter->Query.CData(), bq.CData());
		}

		// increment our alignment counter
		++mCurrentAlignment;
	}

	// ============
	// final output
	// ============

	if(mFlags.IsElandEnabled) fprintf(mStreams.eland, "\n");
}

// parses the specified MOSAIK read file
void CMosaikText::ParseMosaikReadFile(const string& readFilename) {

	MosaikReadFormat::CReadReader reader;
	reader.Open(readFilename);

	// retrieve the sequencing technology
	const MosaikReadFormat::ReadGroup rg = reader.GetReadGroup();
	const bool isColorspace = (rg.SequencingTechnology == ST_SOLID ? true : false);

	// retrieve the read status
	const ReadStatus rs = reader.GetStatus();
	const bool isPairedEnd = ((rs & RS_PAIRED_END_READ) != 0 ? true : false);

	// open our output file streams
	if(mFlags.IsFastqEnabled) {
		mStreams.fastq = fopen(mSettings.FastqFilename.c_str(), "wb");

		if(!mStreams.fastq) {
			printf("ERROR: Unable to open the FASTQ file for writing.\n");
			exit(1);
		}
	}

	// parse all of the reads
	if(!mFlags.IsScreenEnabled) {
		printf("- parsing reads... ");
		fflush(stdout);
	}

	// retrieve all reads from the read reader
	Mosaik::Read mr;
	while(reader.LoadNextRead(mr)) {
		if(mr.Mate1.Bases.Length() != 0) ProcessMate(1, mr.Name, isColorspace, mr.Mate1, isPairedEnd);
		if(mr.Mate2.Bases.Length() != 0) ProcessMate(2, mr.Name, isColorspace, mr.Mate2, isPairedEnd);
	}

	if(!mFlags.IsScreenEnabled) printf("finished.\n");

	// clean up
	reader.Close();
	if(mFlags.IsFastqEnabled) fclose(mStreams.fastq);
}

// processes the mates according to the chosen file format
void CMosaikText::ProcessMate(const unsigned char mateNum, const CMosaikString& readName, const bool isColorspace, Mosaik::Mate& mate, const bool isPairedEnd) {

	// increment the base qualities
	mate.Qualities.Increment(33);

	if(isColorspace) {
		mCS.ConvertReadPseudoColorspaceToColorspace(mate.Bases);
		mate.Bases.Prepend(mate.SolidPrefixTransition, 2);
		if(mFlags.IsFastqEnabled) mate.Qualities.Prepend("!?", 2);
		else mate.Qualities.Prepend("?", 1);
	}

	if(mFlags.IsFastqEnabled) {
		fprintf(mStreams.fastq, "@%s", readName.CData());
		if(isPairedEnd) fprintf(mStreams.fastq, " (mate %u)", mateNum);
		fprintf(mStreams.fastq, "\n%s\n+\n%s\n", mate.Bases.CData(), mate.Qualities.CData());
	}

	if(mFlags.IsScreenEnabled) {
		printf("%s", readName.CData());
		if(isPairedEnd) printf(" (mate %u)", mateNum);
		printf("\n%s\n%s\n\n", mate.Bases.CData(), mate.Qualities.CData());
	}
}

// writes the current alignment to the SAM output file
void CMosaikText::WriteSamEntry(const CMosaikString& readName, const string& readGroupID, const vector<Alignment>::iterator& alIter ) {

	// =================
	// set the SAM flags
	// =================

	unsigned int flag                = 0;
	unsigned int queryPosition5Prime = 0;
	unsigned int matePosition5Prime  = 0;
	int insertSize                   = 0;

	if(alIter->IsPairedEnd) {

		flag |= BAM_SEQUENCED_AS_PAIRS;

		// first or second mate?
		flag |= (alIter->IsFirstMate ? BAM_QUERY_FIRST_MATE : BAM_QUERY_SECOND_MATE);

		if(alIter->IsResolvedAsPair) {

			flag |= BAM_PROPER_PAIR;
			if(alIter->IsMateReverseStrand) flag |= BAM_MATE_REVERSE_COMPLEMENT;

			// sanity check
			if(alIter->ReferenceIndex != alIter->MateReferenceIndex) {
				printf("ERROR: The resolved paired-end reads occur on different reference sequences.\n");
				exit(1);
			}

			// set the 5' coordinates
			queryPosition5Prime = (alIter->IsReverseStrand     ? alIter->ReferenceEnd     : alIter->ReferenceBegin);
			matePosition5Prime  = (alIter->IsMateReverseStrand ? alIter->MateReferenceEnd : alIter->MateReferenceBegin);

			// calculate the insert size
			insertSize = matePosition5Prime - queryPosition5Prime;

		} else flag |= BAM_MATE_UNMAPPED;
	}

	if(alIter->IsReverseStrand) flag |= BAM_QUERY_REVERSE_COMPLEMENT;

	const char* pCigar = cigarTager.GetCigarTag( alIter->Reference.CData(), alIter->Query.CData(), alIter->Reference.Length() );
	const char* pMd    = mdTager.GetMdTag( alIter->Reference.CData(), alIter->Query.CData(), alIter->Reference.Length() );
	
	// sanity check
	// ===================
	// write the alignment
	// ===================

	// B7_591:4:96:693:509	73	seq1	1	99	36M	*	0	0	CACTAGTGGCTCATTGTAAATGTGTGGTTTAACTCG	<<<<<<<<<<<<<<<;<<<<<<<<<5<<<<<;:<;7
	// <QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> <TAGS>

	// remove the gaps from the read
	CMosaikString query(alIter->Query);
	query.Remove('-');

	// shift the base qualities 
	CMosaikString bq(alIter->BaseQualities);
	bq.Increment(33);

	// sanity check
	alIter->BaseQualities.CheckQuality();

	gzprintf(mStreams.sam, "%s\t%u\t%s\t%u\t%u\t%s\t", readName.CData(), flag, alIter->ReferenceName, 
		alIter->ReferenceBegin + 1, alIter->Quality, pCigar);

	if(alIter->IsResolvedAsPair) {
		// N.B. we already checked that the reference indexes were identical
		gzprintf(mStreams.sam, "=\t%u\t%d\t", alIter->MateReferenceBegin + 1, insertSize);
	} else gzprintf(mStreams.sam, "*\t0\t0\t");

	if ( zaString != 0 )
		gzprintf(mStreams.sam, "%s\t%s\tRG:Z:%s\tNM:i:%u\tMD:Z:%s\tZA:Z:%s\n", query.CData(), bq.CData(), readGroupID.c_str(), alIter->NumMismatches, pMd, zaString );
	else
		gzprintf(mStreams.sam, "%s\t%s\tRG:Z:%s\tNM:i:%u\tMD:Z:%s\n", query.CData(), bq.CData(), readGroupID.c_str(), alIter->NumMismatches, pMd);
}
