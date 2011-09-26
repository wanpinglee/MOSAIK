// ***************************************************************************
// CMosaikAligner - delegates the read alignments to worker threads and
//                  displays the final statistics.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "MosaikAligner.h"

#include "QualityNeuralNetwork.h"

// constructor
CMosaikAligner::CMosaikAligner(unsigned char hashSize, CAlignmentThread::AlignerAlgorithmType algorithmType, CAlignmentThread::AlignerModeType algorithmMode, unsigned char numThreads, const string inputCommandLine )
	: mAlgorithm(algorithmType)
	, mMode(algorithmMode)
	, mReference(NULL)
	, mReferenceLength(0)
	, mpDNAHash(NULL)
	, commandLine( inputCommandLine )
{
	// initialization
	mSettings.HashSize            = hashSize;
	mSettings.AllocatedReadLength = 0;
	mSettings.NumThreads          = numThreads;
}

// deconstructor
CMosaikAligner::~CMosaikAligner(void) {

	// delete the hash table
	// this cause some problems with double freeing of memory
	//delete mpDNAHash;
}

void CMosaikAligner::AlignReadArchiveLowMemory(void) {
	// ==============
	// initialization
	// ==============

	string inputReadArchiveFilename  = mSettings.InputReadArchiveFilename;
		
	// define our read format reader and writer
	MosaikReadFormat::CReadReader in;
	in.Open(inputReadArchiveFilename);
	MosaikReadFormat::ReadGroup readGroup = in.GetReadGroup();
	ReadStatus readStatus          = in.GetStatus();
	mSettings.SequencingTechnology = readGroup.SequencingTechnology;
	mSettings.MedianFragmentLength = readGroup.MedianFragmentLength;

	//vector<MosaikReadFormat::ReadGroup> readGroups;
	readGroups.clear();
	readGroups.push_back(readGroup);
	for ( vector<MosaikReadFormat::ReadGroup>::iterator ite = readGroups.begin(); ite != readGroups.end(); ++ite ) {
		readGroupsMap[ ite->ReadGroupCode ] = *ite;
	}

	// close open file streams
	in.Close();

        // retrieve the concatenated reference sequence length
	vector<ReferenceSequence> referenceSequences;
	vector<ReferenceSequence> referenceSequencesBs;

	// get reference information
	MosaikReadFormat::CReferenceSequenceReader refseq;
	refseq.Open(mSettings.ReferenceFilename);
	refseq.GetReferenceSequences(referenceSequences);
	mReferenceLength = refseq.GetReferenceSequenceLength();
	const unsigned int numRefSeqs = refseq.GetNumReferenceSequences();
	refseq.Close();

	// set the bottom lines of LF and MQ in stat map
	mStatisticsMaps.SetLfMin( 0 - (int)mSettings.MedianFragmentLength );
	//mStatisticsMaps.SetMqMin( mStatisticsCounters.StatMappingQuality );

	// retrieve the basespace reference filenames
	if(mFlags.EnableColorspace) {
		
		MosaikReadFormat::CReferenceSequenceReader bsRefSeq;
		bsRefSeq.Open(mSettings.BasespaceReferenceFilename);

		if(!bsRefSeq.HasSameReferenceSequences(referenceSequences)) {
			printf("ERROR: The basespace and colorspace reference sequence archives do not seem to represent the same FASTA file.\n"); 
			exit(1);
		}

		bsRefSeq.GetReferenceSequences(referenceSequencesBs);

		bsRefSeq.Close();
	}

	// check special reference
	if ( mSReference.enable ) {
		// the special references should be appended after normal ones
		for ( vector<ReferenceSequence>::reverse_iterator rit = referenceSequences.rbegin(); rit != referenceSequences.rend(); ++rit ) {
			size_t found = rit->Name.find( mSReference.prefix );
			if ( found != string::npos ) {
				mSReference.found = true;
				mSReference.nReference++;
				mSReference.begin = rit->Begin;
				rit->Species = rit->Name.substr( found + mSReference.prefix.size() + 1, 2 );
				rit->Special = true;
			}
			else
				break;
		}
	}


	// this reference vector is for regular, multiply, and unmapped bams
	vector<ReferenceSequence> referenceSequencesWoSpecial;
	// sanity check
	if ( mSReference.nReference > referenceSequences.size() ) {
		cout << "ERROR: The number of detected special references are larger than the one of input references." << endl;
		exit(1);
	}
	const unsigned int nReference = referenceSequences.size() - mSReference.nReference;
	vector<ReferenceSequence>::const_iterator ite1;
	ite1 = mFlags.EnableColorspace ? referenceSequencesBs.begin() : referenceSequences.begin();
	for ( unsigned int i = 0; i < nReference; ++i ) {
		referenceSequencesWoSpecial.push_back( *ite1 );
		++ite1;
	}
	
	// both of full- and low-memory MOSAIK need multiply-mapped bam in AlignmentThread.cpp
	mBams.mHeader.SortOrder           = SORTORDER_UNSORTED;
	mBams.mHeader.pReferenceSequences = &referenceSequencesWoSpecial;
	mBams.mHeader.pReadGroups         = &readGroups;
	mBams.mBam.Open( mSettings.OutputReadArchiveFilename + ".multiple.bam", mBams.mHeader);

	// ===================
	// full-memory version
	// ===================
	if ( !mFlags.UseLowMemory ) {
		
		// ==============================
		// set the headers of bam writers
		// ==============================
		//mBams.mHeader.SortOrder = SORTORDER_UNSORTED;
		mBams.sHeader.SortOrder = SORTORDER_UNSORTED;
		//mBams.uHeader.SortOrder = SORTORDER_UNSORTED;
		mBams.rHeader.SortOrder = SORTORDER_UNSORTED;

		//mBams.mHeader.pReferenceSequences = &referenceSequencesWoSpecial;
		mBams.sHeader.pReferenceSequences = mFlags.EnableColorspace ? &referenceSequencesBs : &referenceSequences;
		//mBams.uHeader.pReferenceSequences = &referenceSequencesWoSpecial;
		mBams.rHeader.pReferenceSequences = &referenceSequencesWoSpecial;
	
		//mBams.mHeader.pReadGroups = &readGroups;
		mBams.sHeader.pReadGroups = &readGroups;
		//mBams.uHeader.pReadGroups = &readGroups;
		mBams.rHeader.pReadGroups = &readGroups;

		ProgramGroup pg;
		pg.ID = "MosaikAligner";
		stringstream ss;
		ss << (int)MOSAIK_MAJOR_VERSION << "." << (int)MOSAIK_MINOR_VERSION << "." << (int)MOSAIK_BUILD_VERSION;
		pg.VN = ss.str();
		pg.CL = commandLine;

		mBams.sHeader.pg.ID = "MosaikAligner";
		//mBams.uHeader.pg.ID = "MosaikAligner";
		mBams.rHeader.pg.ID = "MosaikAligner";
		mBams.sHeader.pg.VN = ss.str();
		//mBams.uHeader.pg.VN = ss.str();
		mBams.rHeader.pg.VN = ss.str();
		mBams.sHeader.pg.CL = commandLine;
		//mBams.uHeader.pg.CL = commandLine;
		mBams.rHeader.pg.CL = commandLine;

		//mBams.mBam.Open( mSettings.OutputReadArchiveFilename + ".multiple.bam", mBams.mHeader);
		mBams.sBam.Open( mSettings.OutputReadArchiveFilename + ".special.bam", mBams.sHeader);
		//mBams.uBam.Open( mSettings.OutputReadArchiveFilename + ".unaligned.bam", mBams.uHeader);
		mBams.rBam.Open( mSettings.OutputReadArchiveFilename + ".bam", mBams.rHeader);
	

		// prepare BS reference sequence for SOLiD data
		char** pBsRefSeqs = NULL;
		if(mFlags.EnableColorspace) {

			cout << "- loading basespace reference sequences... ";
			cout.flush();

			MosaikReadFormat::CReferenceSequenceReader bsRefSeq;
			bsRefSeq.Open(mSettings.BasespaceReferenceFilename);


			bsRefSeq.CopyReferenceSequences(pBsRefSeqs);
			bsRefSeq.Close();

			cout << "finished." << endl;
		}

		// prepare reference sequence
		refseq.Open(mSettings.ReferenceFilename);
		cout << "- loading reference sequence... ";
		cout.flush();
		refseq.LoadConcatenatedSequence(mReference);
		cout << "finished." << endl;
		refseq.Close();
		
		unsigned int* pRefBegin = new unsigned int[numRefSeqs];
		unsigned int* pRefEnd   = new unsigned int[numRefSeqs];
		char** pRefSpecies      = new char*[numRefSeqs];
		bool* pRefSpecial       = new bool[numRefSeqs];
		for(unsigned int j = 0; j < numRefSeqs; j++) {
			pRefBegin[j]   = referenceSequences[j].Begin;
			pRefEnd[j]     = referenceSequences[j].End;
			pRefSpecial[j] = referenceSequences[j].Special;
			if ( pRefSpecial[j] ) {
				pRefSpecies[j] = new char [3];
				memcpy( pRefSpecies[j], referenceSequences[j].Species.data(), 2 );
				pRefSpecies[j][2] = 0;
			} else {
				pRefSpecies[j] = new char [1];
				pRefSpecies[j][0] = 0;
			}
		}
		
		// initialize our hash tables
		if(!mFlags.IsUsingJumpDB) {
			InitializeHashTables(CalculateHashTableSize(mReferenceLength, mSettings.HashSize), 0, 0, 0, mFlags.UseLowMemory, 0, false);
			HashReferenceSequence(refseq);
		}
		else {
			InitializeHashTables(CalculateHashTableSize(mReferenceLength, mSettings.HashSize), pRefBegin[0], pRefEnd[numRefSeqs - 1], 0, mFlags.UseLowMemory, 0, mSReference.found);
			mpDNAHash->LoadKeysNPositions();
		}

		// set the hash positions threshold
		if(mFlags.IsUsingHashPositionThreshold && (mAlgorithm == CAlignmentThread::AlignerAlgorithm_ALL))
			mpDNAHash->RandomizeAndTrimHashPositions(mSettings.HashPositionThreshold);

		// localize the read archive filenames
		string outputReadArchiveFilename = mSettings.OutputReadArchiveFilename;

		// define our read format reader and writer
		MosaikReadFormat::CReadReader inn;
		inn.Open(inputReadArchiveFilename);

		// set the alignment status flags
		AlignmentStatus alignmentStatus = AS_UNSORTED_READ | readStatus;
		if(mMode == CAlignmentThread::AlignerMode_ALL) alignmentStatus |= AS_ALL_MODE;
		else alignmentStatus |= AS_UNIQUE_MODE;

		// define our output file
		mFlags.UseBamOutput     = true;
		mFlags.UseArchiveOutput = false;
		
		MosaikReadFormat::CAlignmentWriter out;
		//out.Open(mSettings.OutputReadArchiveFilename.c_str(), referenceSequences, readGroups, alignmentStatus, ALIGNER_SIGNATURE);

		mFlags.SaveMultiplyBam = true;
		AlignReadArchive(inn, out, pRefBegin, pRefEnd, pRefSpecies, pRefSpecial, pBsRefSeqs, 0);

		// close open file streams
		inn.Close();

		// solid references should be one-base longer after converting back to basespace
		//if(mFlags.EnableColorspace) out.AdjustSolidReferenceBases();
		//out.Close();

		// free memory
		if(mFlags.IsUsingJumpDB) mpDNAHash->FreeMemory();
		if(pRefBegin)   delete [] pRefBegin;
		if(pRefEnd)     delete [] pRefEnd;
		if(mReference)  delete [] mReference;
		if(pRefSpecial) delete [] pRefSpecial;
		if(pBsRefSeqs) {
			for(unsigned int i = 0; i < numRefSeqs; ++i) delete [] pBsRefSeqs[i];
			delete [] pBsRefSeqs;
		}
		if(pRefSpecies) {
			for(unsigned int i = 0; i < numRefSeqs; ++i) delete [] pRefSpecies[i];
			delete [] pRefSpecies;
		}
		pRefBegin   = NULL;
		pRefEnd     = NULL;
		mReference  = NULL;
		pRefSpecial = NULL;
		pBsRefSeqs  = NULL;
		pRefSpecies = NULL;


		// close mBams
		//mBams.mBam.Close();
		mBams.sBam.Close();
		//mBams.uBam.Close();
		mBams.rBam.Close();
	
	}
	else {
/*
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/1");
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/2");
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/3");
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/4");
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/5");
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/6");
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/7");
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/8");
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/9");
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/10");
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/11");
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/12");
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/13");
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/14");
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/15");
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/16");
            outputFilenames.push_back("/home/wanping/tools/seqan-trunk/build/Release/reads/illumina_readLength100/bak/test/17");
            //outputFilenames.push_back("/home/wanping/bugTest/Al/buffer/tmp/SE/18");
            MergeArchives();

	//string tempDir1;
	//CFileUtilities::GetTempDirectory( tempDir1 );
	//if ( CFileUtilities::DirExists( tempDir1.c_str() ) ) {
	//        cout << endl << "- cleaning up temp files...";
	//	CFileUtilities::DeleteDir( tempDir1 );
	//	cout << "finished." << endl;
	//}
            PrintStatistics();
            exit(1);
*/
		
		// grouping reference and store information in referenceGroups vector
		GroupReferences( referenceSequences );
		
		// get hash statistics for adjusting mhp for each reference group and reserve memory
		vector< unsigned int > nHashs;             // the numbers of hash positions in each reference group
		vector< unsigned int > expectedMemories;   // the numbers of hashs in each reference group
		uint64_t nTotalHash;
		GetHashStatistics( nHashs, expectedMemories, nTotalHash, referenceSequences );
		

		// align reads again per chromosome group
		for ( unsigned int i = 0; i < referenceGroups.size(); ++i ) {
	        	unsigned int startRef = referenceGroups[i].first;
			unsigned int endRef   = referenceGroups[i].first + referenceGroups[i].second - 1;

			CConsole::Heading();
		        if ( referenceGroups[i].second > 1 )
				cout << endl << "Aligning chromosome " << startRef + 1 << "-" << endRef + 1 << " (of " << numRefSeqs << "):" << endl;
			else
				cout << endl << "Aligning chromosome " << startRef + 1 << " (of " << numRefSeqs << "):" << endl;
		        CConsole::Reset();

			// initialize our hash tables
			// calculate expected memories for jump data
			unsigned int expectedMemory = nHashs[i] + expectedMemories[i];
			// reserve 3% more memory for unexpected usage
			// TODO: for small mhp, 3% more allowed space is not enough
			expectedMemory =  expectedMemory * 1.03;

			InitializeHashTables(
				0, 
				referenceSequences[startRef].Begin, 
				referenceSequences[endRef].End, 
				referenceSequences[startRef].Begin, 
				mFlags.UseLowMemory, expectedMemory, false);

			// set the hash positions threshold
			if(mFlags.IsUsingHashPositionThreshold && (mAlgorithm == CAlignmentThread::AlignerAlgorithm_ALL)) { 
				double ratio = nHashs[i] / (double)nTotalHash;
				unsigned int positionThreshold = ( mSReference.found && ( i == referenceGroups.size() - 1 ) ) 
					? mSReference.count 
					: ceil(ratio * (double)mSettings.HashPositionThreshold);
				mpDNAHash->RandomizeAndTrimHashPositions(positionThreshold);
			}

			// load jump data
			mpDNAHash->LoadKeysNPositions();

			// set reference information
			unsigned int* pRefBegin = new unsigned int[ referenceGroups[i].second ];
			unsigned int* pRefEnd   = new unsigned int[ referenceGroups[i].second ];
			char** pRefSpecies      = new char* [ referenceGroups[i].second ];
			bool* pRefSpecial       = new bool [ referenceGroups[i].second ];
			for ( unsigned int j = 0; j < referenceGroups[i].second; j++ ){
				pRefBegin[j]    = referenceSequences[ startRef + j ].Begin - referenceSequences[ startRef ].Begin;
				pRefEnd[j]      = referenceSequences[ startRef + j ].End   - referenceSequences[ startRef ].Begin;
				pRefSpecial[j]  = referenceSequences[ startRef + j ].Special;

				if ( pRefSpecial[j] ) {
					pRefSpecies[j]  = new char [3];
					memcpy( pRefSpecies[j], referenceSequences[ startRef + j ].Species.data(), 2 );
					pRefSpecies[j][2] = 0;
				} else {
					pRefSpecies[j]  = new char [1];
					pRefSpecies[j][0] = 0;
				}
			}

			// prepare BS reference sequence for SOLiD data
			char** pBsRefSeqs = NULL;
			if(mFlags.EnableColorspace) {
	
				cout << "- loading basespace reference sequences... ";
				cout.flush();

				MosaikReadFormat::CReferenceSequenceReader bsRefSeq;
				bsRefSeq.Open(mSettings.BasespaceReferenceFilename);


				bsRefSeq.CopyReferenceSequences(pBsRefSeqs, startRef, referenceGroups[i].second);
				bsRefSeq.Close();

				cout << "finished." << endl;
			}

			// prepare reference sequence
			refseq.Open(mSettings.ReferenceFilename);
			cout << "- loading reference sequence... ";
			cout.flush();
			//refseq.LoadConcatenatedSequence(mReference);
			refseq.LoadConcatenatedSequence(mReference, startRef, referenceGroups[i].second);
			refseq.Close();

			cout << "finished." << endl;
			
			
			// localize the read archive filenames
			// get a temporary file name
			string tempFilename;
			CFileUtilities::GetTempFilename(tempFilename);
			outputFilenames.push_back(tempFilename);

			// define our read format reader and writer
			MosaikReadFormat::CReadReader inn;
			inn.Open(inputReadArchiveFilename);

			// set the alignment status flags
			AlignmentStatus alignmentStatus = AS_UNSORTED_READ | readStatus;
			if(mMode == CAlignmentThread::AlignerMode_ALL) alignmentStatus |= AS_ALL_MODE;
			else alignmentStatus |= AS_UNIQUE_MODE;

			// prepare a new vector for the current chromosome for opening out archive
			vector<ReferenceSequence> smallReferenceSequences;
			for ( unsigned int j = 0; j < referenceGroups[i].second; j++ ){
				if ( mFlags.EnableColorspace )
					smallReferenceSequences.push_back(referenceSequencesBs[startRef+j]);
				else
					smallReferenceSequences.push_back(referenceSequences[startRef+j]);

			}

			// define our output file
			mFlags.UseBamOutput     = false;
			mFlags.UseArchiveOutput = true;

			MosaikReadFormat::CAlignmentWriter out;
			out.Open(tempFilename.c_str(), smallReferenceSequences, readGroups, alignmentStatus, ALIGNER_SIGNATURE);
			//out.AdjustPartitionSize(20000/referenceGroups.size());
			out.AdjustPartitionSize(1000);

			mFlags.SaveMultiplyBam = ( mSReference.found && ( i == referenceGroups.size() - 1 ) ) ? false : true;
			mFlags.SaveUnmappedBasesInArchive = ( i == 0 ) ? true : false; 
			AlignReadArchive(inn, out, pRefBegin, pRefEnd, pRefSpecies, pRefSpecial, pBsRefSeqs, referenceGroups[i].first );

			// close open file streams
			inn.Close();

			// solid references should be one-base longer after converting back to basespace
			//if(mFlags.EnableColorspace) out.AdjustSolidReferenceBases();
			out.Close();

			// free memory
			if(mFlags.IsUsingJumpDB) mpDNAHash->FreeMemory();
			if(pRefBegin)   delete [] pRefBegin;
			if(pRefEnd)     delete [] pRefEnd;
			if(mReference)  delete [] mReference;
			if(pRefSpecial) delete [] pRefSpecial;
			if(pBsRefSeqs) {
				for(unsigned int j = 0; j < referenceGroups[i].second; j++) delete [] pBsRefSeqs[j];
				delete [] pBsRefSeqs;
			}
			if(pRefSpecies) {
				for(unsigned int j = 0; j < referenceGroups[i].second; j++) delete [] pRefSpecies[j];
				delete [] pRefSpecies;
			}
			pRefBegin   = NULL;
			pRefEnd     = NULL;
			mReference  = NULL;
			pRefSpecial = NULL;
			pBsRefSeqs  = NULL;
			pRefSpecies = NULL;
		}
	}

	mBams.mBam.Close();

	if ( mFlags.UseLowMemory )
		MergeArchives();

	// clean up temp files
	string tempDir;
	CFileUtilities::GetTempDirectory( tempDir );
	if ( CFileUtilities::DirExists( tempDir.c_str() ) ) {
	        cout << endl << "- cleaning up temp files...";
		CFileUtilities::DeleteDir( tempDir );
		cout << "finished." << endl;
	}

	if ( !mSReference.found ) {
		string sBamName = mSettings.OutputReadArchiveFilename + ".special.bam";
		rm ( sBamName.c_str() );
	}


	PrintStatistics();
}


void CMosaikAligner::GetHashStatistics( 
	vector<unsigned int>& nHashs, 
	vector<unsigned int>& expectedMemories, 
	uint64_t& nTotalHash, 
	const vector<ReferenceSequence>& referenceSequences ) {

	//unsigned int length = referenceSequences.size();
	unsigned int length = referenceGroups.size();
	unsigned int begin  = referenceSequences[0].Begin;
	unsigned int end    = referenceSequences[ referenceSequences.size() - 1].End;
	unsigned int offset = 0;

	// initial JumpDnaHash
	CJumpDnaHash hash(
		mSettings.HashSize, 
		mSettings.JumpFilenameStub, 
		0, 
		mFlags.KeepJumpKeysInMemory, 
		mFlags.KeepJumpPositionsInMemory, 
		mSettings.NumCachedHashes, 
		begin, end, offset, 
		0, false, false, 0, 0.0);

	// set mhp number to JumpDnaHash
	if(mFlags.IsUsingHashPositionThreshold && (mAlgorithm == CAlignmentThread::AlignerAlgorithm_ALL))
		hash.RandomizeAndTrimHashPositions(mSettings.HashPositionThreshold);

	vector< pair<unsigned int, unsigned int> > references;
	for ( unsigned int i = 0; i < length; i++ ) {
		unsigned int startRef = referenceGroups[i].first;
		unsigned int endRef   = referenceGroups[i].first + referenceGroups[i].second - 1;
		pair<unsigned int, unsigned int> temp(referenceSequences[startRef].Begin, referenceSequences[endRef].End);
//cerr << referenceGroups[i].first << "\t" << referenceGroups[i].second << "\t" << referenceSequences[endRef].End - referenceSequences[startRef].Begin + 1 << endl;
		references.push_back(temp);
	}
	//exit(1);

	nHashs.resize(length, 0);
	expectedMemories.resize(length, 0);
	hash.GetHashStatistics(references, nHashs, expectedMemories, mSReference.found, mSReference.begin, mSReference.count );

	nTotalHash = 0;
	for ( unsigned int i = 0; i < nHashs.size(); i++ ) {
		nTotalHash += nHashs[i];
	}

}


void CMosaikAligner::GroupReferences( const vector<ReferenceSequence>& referenceSequences ) {
	// find the largest reference
	unsigned int max = 0;
	unsigned int length = 0;
	for ( unsigned int i = 0; i < referenceSequences.size(); i++ ) {
		length = referenceSequences[i].End - referenceSequences[i].Begin + 1;
		max = ( length > max ) ? length : max;
	}

	unsigned int start = 0;
	unsigned int accLength = referenceSequences[0].End - referenceSequences[0].Begin + 1;
	
	unsigned int nNormalReferences = mSReference.found ? referenceSequences.size() - mSReference.nReference : referenceSequences.size() ;
	//unsigned int nNormalReferences = referenceSequences.size();
	
	for ( unsigned int i = 1; i < nNormalReferences; i++ ) {
		length = referenceSequences[i].End - referenceSequences[i].Begin + 1;
		if ( ( accLength + length ) > max ) {
			pair<unsigned int, unsigned int> tmp (start, i - start);
			referenceGroups.push_back(tmp);
			start = i;
			accLength = length;
		}
		else
			accLength += length;

	}

	if ( accLength != 0 ) {
		pair<unsigned int, unsigned int> tmp (start, nNormalReferences - start);
		referenceGroups.push_back(tmp);
	}


	// handle special references
	if ( mSReference.found ) {
		start = referenceSequences.size() - mSReference.nReference;
		pair<unsigned int, unsigned int> tmp (start, mSReference.nReference);
		referenceGroups.push_back(tmp);
	}


	//for ( unsigned int i = 0; i < referenceGroups.size(); i++) {
	//	cout << referenceGroups[i].first << "\t" << referenceGroups[i].second << endl;
	//}
	

}


void CMosaikAligner::MergeArchives(void) {

	// set active threads
	unsigned int nThread = ( mSettings.NumThreads < outputFilenames.size() ) ? mSettings.NumThreads : outputFilenames.size();
	
	vector< string > temporaryFiles;
	temporaryFiles.resize( outputFilenames.size() );
	for ( unsigned int i = 0; i < outputFilenames.size(); i++ ) {
		string tempFilename;
		CFileUtilities::GetTempFilename(tempFilename);
		temporaryFiles[i] = tempFilename;
	}

        // calculate total # of reads
        uint64_t nReads = 0;
	//unsigned int nOutputFilenames = mSReference.found ? outputFilenames.size() - 1 : outputFilenames.size() ;
        for ( unsigned int i = 0 ; i < outputFilenames.size(); i++ ) {
	        MosaikReadFormat::CAlignmentReader reader;
                reader.Open( outputFilenames[i] );
                nReads += reader.GetNumReads();
                reader.Close();
        }

	// at most, total 8,000,000 alignments will be put in memory
	const unsigned int alignmentInMemory        = 8000000;
	const uint64_t     avgReads                 = nReads / outputFilenames.size();
	unsigned int       alignmentCacheEachSorter = alignmentInMemory / nThread;
	unsigned int       fileEachSorter           = avgReads / alignmentCacheEachSorter;

	while ( ( nThread > 1 ) && ( fileEachSorter * nThread ) > 1000 ) {
		--nThread;
		alignmentCacheEachSorter = alignmentInMemory / nThread;
		fileEachSorter           = avgReads / alignmentCacheEachSorter;
	}

	CConsole::Heading();
	cout << endl << "Sorting alignment archive:" << endl;
	CConsole::Reset();
	SortThread sThread ( outputFilenames, temporaryFiles, nThread, nReads, mSettings.MedianFragmentLength, alignmentCacheEachSorter );
	if ( mFlags.IsQuietMode )
		sThread.SetQuietMode();
	sThread.Start();

	for ( unsigned int i = 0; i < outputFilenames.size(); i++ )
		//cerr << outputFilenames[i] << endl;
		rm(outputFilenames[i].c_str());
	//cerr << endl;


	CConsole::Heading();
	cout << "Merging alignment archive:" << endl;
	CConsole::Reset();


        uint64_t readNo        = 0;
	if ( !mFlags.IsQuietMode )
        	CProgressBar<uint64_t>::StartThread(&readNo, 0, nReads, "reads");

        CArchiveMerge merger( 
		temporaryFiles, mSettings.OutputReadArchiveFilename, &readNo, mFlags.EnableColorspace, commandLine, 
		mSettings.MedianFragmentLength, mSettings.LocalAlignmentSearchRadius, mSReference.found, mStatisticsCounters.StatMappingQuality );

	merger.Merge();

	if ( !mFlags.IsQuietMode )
	        CProgressBar<uint64_t>::WaitThread();

	
	for ( unsigned int i = 0; i < temporaryFiles.size(); i++ )
		//cerr << temporaryFiles[i] << endl;
		rm(temporaryFiles[i].c_str());

	// get statistics information
	float allowedMm = mFilters.UseMismatchFilter ? (float)mFilters.MaxNumMismatches :  mFilters.MaxMismatchPercent;
	string mapFile = mSettings.OutputReadArchiveFilename + ".stat";
	merger.PrintStatisticsMaps( mapFile, readGroups, mSettings.MedianFragmentLength, mSettings.LocalAlignmentSearchRadius, allowedMm  );


	CArchiveMerge::StatisticsCounters mergeCounters;
	merger.GetStatisticsCounters( mergeCounters );


	// mates
	mStatisticsCounters.FilteredOutMates   = mergeCounters.FilteredOutMates;
	mStatisticsCounters.MultipleMates      = mergeCounters.MultipleMates;
	mStatisticsCounters.UniqueMates        = mergeCounters.UniqueMates;
	mStatisticsCounters.Unmapped           = mergeCounters.Unmapped;

	// pairs
	mStatisticsCounters.UU                   = mergeCounters.UU;
	mStatisticsCounters.UM                   = mergeCounters.UM;
	mStatisticsCounters.UF                   = mergeCounters.UF;
	mStatisticsCounters.MM                   = mergeCounters.MM;
	mStatisticsCounters.MF                   = mergeCounters.MF;
	mStatisticsCounters.UX                   = mergeCounters.UX;
	mStatisticsCounters.MX                   = mergeCounters.MX;
	mStatisticsCounters.FF                   = mergeCounters.FF;
	mStatisticsCounters.FX                   = mergeCounters.FX;
	mStatisticsCounters.XX                   = mergeCounters.XX;
	mStatisticsCounters.UU_localRescue       = mergeCounters.UU_localRescue;
	mStatisticsCounters.UU_localConsistance  = mergeCounters.UU_localConsistance;
	mStatisticsCounters.UM_localRescue       = mergeCounters.UM_localRescue;
	mStatisticsCounters.UM_localConsistance  = mergeCounters.UM_localConsistance;
	mStatisticsCounters.MM_localRescue       = mergeCounters.MM_localRescue;
	mStatisticsCounters.MM_localConsistance  = mergeCounters.MM_localConsistance;

}

// aligns the read archive
void CMosaikAligner::AlignReadArchive(
	MosaikReadFormat::CReadReader& in, 
	MosaikReadFormat::CAlignmentWriter& out, 
	unsigned int* pRefBegin, 
	unsigned int* pRefEnd, 
	char** pRefSpecies, 
	bool*  pRefSpecial, 
	char** pBsRefSeqs,
	const unsigned int referenceOffset) {

	ReadStatus readStatus          = in.GetStatus();
	
	const bool isPairedEnd = (readStatus == RS_PAIRED_END_READ ? true : false);

	// open the unaligned read report file
	//FILE* unalignedStream = NULL;
	//if(mFlags.IsReportingUnalignedReads) {
	//	if(fopen_s(&unalignedStream, mSettings.UnalignedReadReportFilename.c_str(), "wb") != 0) {
	//		cout << "ERROR: Unable to open the unaligned read FASTQ file for output." << endl;
	//		exit(1);
	//	}
	//}
	


	// localize our read and reference counts. Initialize our statistical counters
	uint64_t numReadArchiveReads = in.GetNumReads();
	uint64_t readCounter = 0;

	// initialize our threads
	pthread_t* activeThreads = new pthread_t[mSettings.NumThreads];

	CAlignmentThread::ThreadData td;
	td.Algorithm           = mAlgorithm;
	td.ReferenceLen        = mReferenceLength;
	td.Filters             = mFilters;
	td.Flags               = mFlags;
	td.Mode                = mMode;
	td.pReference          = mReference;
	td.pCounters           = &mStatisticsCounters;
	td.pMaps               = &mStatisticsMaps;
	td.pDnaHash            = mpDNAHash;
	td.pIn                 = &in;
	td.pOut                = &out;
	//td.pUnalignedStream    = unalignedStream;
	td.pRefBegin           = pRefBegin;
	td.pRefEnd             = pRefEnd;
	td.pRefSpecies         = pRefSpecies;
	td.pRefSpecial         = pRefSpecial;
	td.Settings            = mSettings;
	td.pReadCounter        = &readCounter;
	td.IsPairedEnd         = isPairedEnd;
	td.pBsRefSeqs          = pBsRefSeqs;
	td.pBams               = &mBams;
	td.SpecialReference    = mSReference;
	td.pReadGroups         = &readGroupsMap;
	td.ReferenceOffset     = referenceOffset;
	td.paired_end_ann_file = mPeNeuralNetworkFilename;
	td.single_end_ann_file = mSeNeuralNetworkFilename;


	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	pthread_mutex_init(&CAlignmentThread::mGetReadMutex,              NULL);
	//pthread_mutex_init(&CAlignmentThread::mReportUnalignedMate1Mutex, NULL);
	//pthread_mutex_init(&CAlignmentThread::mReportUnalignedMate2Mutex, NULL);
	pthread_mutex_init(&CAlignmentThread::mSaveReadMutex,             NULL);
	pthread_mutex_init(&CAlignmentThread::mStatisticsMutex,           NULL);
	pthread_mutex_init(&CAlignmentThread::mStatisticsMapsMutex,       NULL);
	pthread_mutex_init(&CAlignmentThread::mSaveMultipleBamMutex,      NULL);
	pthread_mutex_init(&CAlignmentThread::mSaveSpecialBamMutex,       NULL);
	pthread_mutex_init(&CAlignmentThread::mSaveUnmappedBamMutex,      NULL);

	pthread_mutex_init(&CAbstractDnaHash::mJumpCacheMutex,            NULL);
	pthread_mutex_init(&CAbstractDnaHash::mJumpKeyMutex,              NULL);
	pthread_mutex_init(&CAbstractDnaHash::mJumpPositionMutex,         NULL);

	// ===========================
	// start our alignment threads
	// ===========================

	// initialize our progress bar
	
	if ( !mFlags.UseLowMemory ) {
		CConsole::Heading();
		cout << endl;
	}
	cout << "Aligning read library (" << numReadArchiveReads << "):" << endl;
	if ( !mFlags.UseLowMemory )
	CConsole::Reset();

	if ( !mFlags.IsQuietMode )
		CProgressBar<uint64_t>::StartThread(&readCounter, 0, numReadArchiveReads, "reads");

	// create our threads
	for(unsigned int i = 0; i < mSettings.NumThreads; i++)
		pthread_create(&activeThreads[i], &attr, CAlignmentThread::StartThread, (void*)&td);

	pthread_attr_destroy(&attr);

	CBenchmark alignmentBench;
	alignmentBench.Start();

	// wait for the threads to complete
	void* status = NULL;
	for(unsigned int i = 0; i < mSettings.NumThreads; i++) 
		pthread_join(activeThreads[i], &status);

	// wait for the progress bar to finish
	if ( !mFlags.IsQuietMode )
		CProgressBar<uint64_t>::WaitThread();

	alignmentBench.Stop();

	// free up some memory
	delete [] activeThreads;
	activeThreads = NULL;

	//if(mFlags.IsReportingUnalignedReads) fclose(unalignedStream);
	
}

// print our statistics
void CMosaikAligner::PrintStatistics () {

	// for low-memory version, the map is printed when merging archive.
	if ( !mFlags.UseLowMemory ) {
		float allowedMm = mFilters.UseMismatchFilter ? (float)mFilters.MaxNumMismatches : mFilters.MaxMismatchPercent;
		mStatisticsMaps.SetExpectedStatistics( mSettings.MedianFragmentLength, mSettings.LocalAlignmentSearchRadius, allowedMm );
		string mapFile = mSettings.OutputReadArchiveFilename + ".stat";
		mStatisticsMaps.PrintMaps( mapFile.c_str(), readGroups, mStatisticsCounters.StatMappingQuality );
	}
	
	MosaikReadFormat::CReadReader in;
        string inputReadArchiveFilename  = mSettings.InputReadArchiveFilename;
	in.Open(inputReadArchiveFilename);

	const uint64_t numReadArchiveReads = in.GetNumReads();
	ReadStatus readStatus        = in.GetStatus();
	const bool isPairedEnd = (readStatus == RS_PAIRED_END_READ ? true : false);
	const uint64_t totalMates = isPairedEnd ? numReadArchiveReads * 2 : numReadArchiveReads;
	
	if ( mFlags.UseLowMemory ) {
		mStatisticsCounters.ShortMates       = 0;
		mStatisticsCounters.FailedHashMates  = 0;
		mStatisticsCounters.TooManyNsMates   = 0;
	}

	const uint64_t totalAlignedMates = mStatisticsCounters.UniqueMates + mStatisticsCounters.MultipleMates;
	const uint64_t totalAlignedReads = mStatisticsCounters.UU + mStatisticsCounters.UM + mStatisticsCounters.UX + mStatisticsCounters.MM + mStatisticsCounters.MX + mStatisticsCounters.UF + mStatisticsCounters.MF;;

	// print our alignment statistics (mates) if don't enable low-memory algorithm

	
	printf("\n");
	CConsole::Heading(); printf("Alignment statistics (mates):\n"); CConsole::Reset();
	printf("================================================\n");

	uint64_t totalUnaligned = mFlags.UseLowMemory ? mStatisticsCounters.Unmapped
		: mStatisticsCounters.ShortMates + mStatisticsCounters.TooManyNsMates + mStatisticsCounters.FailedHashMates;
	if ( totalUnaligned > 0 ) {
		if ( !mFlags.UseLowMemory ) {
			if ( mStatisticsCounters.ShortMates > 0 )
				printf("   # too short:              %9llu (%5.1f %%)\n", (unsigned long long) mStatisticsCounters.ShortMates, ( mStatisticsCounters.ShortMates / (double)totalMates) * 100.0);
			if ( mStatisticsCounters.TooManyNsMates > 0 )
				printf("   # too many N's:           %9llu (%5.1f %%)\n", (unsigned long long) mStatisticsCounters.TooManyNsMates, ( mStatisticsCounters.TooManyNsMates / (double)totalMates) * 100.0);
			if ( mStatisticsCounters.FailedHashMates > 0 )
				printf("   # failed hash:            %9llu (%5.1f %%)\n", (unsigned long long) mStatisticsCounters.FailedHashMates, ( mStatisticsCounters.FailedHashMates/ (double)totalMates) * 100.0);

			printf("------------------------------------------------\n");
		}
		printf("# unaligned mates(");
		CConsole::Bold();
		printf("X");
		CConsole::Reset();
		printf("):        %9llu (%5.1f %%)\n", (unsigned long long) totalUnaligned, ( totalUnaligned / (double)totalMates) * 100.0);
	}

	if ( mStatisticsCounters.FilteredOutMates > 0 ) {
		printf("# filtered out(");
		CConsole::Bold();
		printf("F");
		CConsole::Reset();
		printf("):           %9llu (%5.1f %%)\n", (unsigned long long) mStatisticsCounters.FilteredOutMates, ( mStatisticsCounters.FilteredOutMates / (double)totalMates) * 100.0);
	}

	
	if(mStatisticsCounters.UniqueMates > 0) {
		printf("# uniquely aligned mates(");
		CConsole::Bold();
		printf("U");
		CConsole::Reset();
		printf("): %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.UniqueMates,      (mStatisticsCounters.UniqueMates      / (double)totalMates) * 100.0);
	}

	if(mStatisticsCounters.MultipleMates > 0) {
		printf("# multiply aligned mates(");
		CConsole::Bold();
		printf("M");
		CConsole::Reset();
		printf("): %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.MultipleMates,   (mStatisticsCounters.MultipleMates   / (double)totalMates) * 100.0);
	}

	printf("================================================\n");
	printf("total aligned:               ");
	CConsole::Bold(); printf("%9llu", (unsigned long long)totalAlignedMates); CConsole::Reset();
	printf(" (");
	CConsole::Bold(); printf("%5.1f %%", (totalAlignedMates / (double)totalMates) * 100.0); CConsole::Reset();
	printf(")\n");
	printf("total:                       %9llu\n", (unsigned long long)totalMates);

	// print our local alignment search statistics
	// we don't print out local alignment information when the low-memory approach is enabled.
	//if( !mFlags.UseLowMemory && mFlags.UseLocalAlignmentSearch ) {
	//	printf("\n");
	//	CConsole::Heading(); printf("Local alignment search statistics:\n"); CConsole::Reset();
	//	printf("===================================\n");

	//	double rescuedAlignmentsPercent = mStatisticsCounters.AdditionalLocalMates / (double)totalMates * 100.0;
	//	printf("rescued mates:  %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.AdditionalLocalMates, rescuedAlignmentsPercent);
	//}

	// print our alignment statistics (reads)

	uint64_t totalRescues       = mStatisticsCounters.UU_localRescue + mStatisticsCounters.UM_localRescue + mStatisticsCounters.MM_localRescue;
	uint64_t totalConsistencies = mStatisticsCounters.UU_localConsistance + mStatisticsCounters.UM_localConsistance + mStatisticsCounters.MM_localConsistance;
	uint64_t totalPartially     = mStatisticsCounters.UX + mStatisticsCounters.MX + mStatisticsCounters.UF + mStatisticsCounters.MF;
	uint64_t totalCompletely    = mStatisticsCounters.UU + mStatisticsCounters.UM + mStatisticsCounters.MM;
	uint64_t totalMissing       = mStatisticsCounters.XX + mStatisticsCounters.FF + mStatisticsCounters.FX;
	if( ( totalPartially + totalCompletely + totalMissing ) > 0 ) {
		printf("\n");
		CConsole::Heading(); printf("Alignment statistics (pairs):\n"); CConsole::Reset();
		printf("===================================================================\n");
		CConsole::Bold(); 
		printf("                                  Local rescues   Frag. consistency\n");
		CConsole::Reset();
		
		if ( totalCompletely > 0 ) {
			CConsole::Bold();
			printf("Completely aligned pairs\n");
			CConsole::Reset();
			printf("-------------------------------------------------------------------\n");
			if ( mStatisticsCounters.UU > 0 )
			printf("# U-U pairs:  %9llu (%5.1f %%) %13llu   %17llu\n", (unsigned long long)mStatisticsCounters.UU, (mStatisticsCounters.UU / (double)numReadArchiveReads) * 100.0, (unsigned long long)mStatisticsCounters.UU_localRescue, (unsigned long long)mStatisticsCounters.UU_localConsistance );
		
			if ( mStatisticsCounters.UM > 0 )
			printf("# U-M pairs:  %9llu (%5.1f %%) %13llu   %17llu\n", (unsigned long long)mStatisticsCounters.UM, (mStatisticsCounters.UM / (double)numReadArchiveReads) * 100.0, (unsigned long long)mStatisticsCounters.UM_localRescue, (unsigned long long)mStatisticsCounters.UM_localConsistance);
		
			if ( mStatisticsCounters.MM > 0 )
			printf("# M-M pairs:  %9llu (%5.1f %%) %13llu   %17llu\n", (unsigned long long)mStatisticsCounters.MM, (mStatisticsCounters.MM / (double)numReadArchiveReads) * 100.0, (unsigned long long)mStatisticsCounters.MM_localRescue, (unsigned long long)mStatisticsCounters.MM_localConsistance);
		}
		
		if ( totalPartially > 0 ) {
			CConsole::Bold();
			printf("\nPartially aligned pairs\n");
			CConsole::Reset();
			printf("-------------------------------------------------------------------\n");
			
			if ( mStatisticsCounters.UF > 0 )
			printf("# U-F pairs:  %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.UF, (mStatisticsCounters.UF / (double)numReadArchiveReads) * 100.0);
			
			if ( mStatisticsCounters.UX > 0 )
			printf("# U-X pairs:  %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.UX, (mStatisticsCounters.UX / (double)numReadArchiveReads) * 100.0);

			if ( mStatisticsCounters.MF > 0 )
			printf("# M-F pairs:  %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.MF, (mStatisticsCounters.MF / (double)numReadArchiveReads) * 100.0);
			
			if ( mStatisticsCounters.MX > 0 )
			printf("# M-X pairs:  %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.MX, (mStatisticsCounters.MX / (double)numReadArchiveReads) * 100.0);
		}

		if ( totalMissing > 0 ) {
			CConsole::Bold();
			printf("\nBoth ends unaliged pairs\n");
			CConsole::Reset();
			printf("-------------------------------------------------------------------\n");
			
			if ( mStatisticsCounters.FF > 0 )
			printf("# F-F pairs:  %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.FF, (mStatisticsCounters.FF / (double)numReadArchiveReads) * 100.0);
			
			if ( mStatisticsCounters.FF > 0 )
			printf("# F-X pairs:  %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.FX, (mStatisticsCounters.FX / (double)numReadArchiveReads) * 100.0);

			if ( mStatisticsCounters.XX > 0 )
			printf("# X-X pairs:  %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.XX, (mStatisticsCounters.XX / (double)numReadArchiveReads) * 100.0);
		}

		printf("===================================================================\n");
		printf("total aligned:");
		CConsole::Bold(); printf("%9llu", (unsigned long long)totalAlignedReads); CConsole::Reset();
		printf(" (");
		CConsole::Bold(); printf("%5.1f %%", (totalAlignedReads / (double)numReadArchiveReads) * 100.0); CConsole::Reset();
		printf(") %13llu   %17llu\n", (unsigned long long)totalRescues, (unsigned long long)totalConsistencies);
		printf("total:        ");
		printf("%9llu", (unsigned long long) numReadArchiveReads);
		printf("\n");
	}

	// print our jump cache statistics
	if( !mFlags.UseLowMemory && mFlags.IsUsingJumpDB && (mSettings.NumCachedHashes > 0)) {
		printf("\n");
		CConsole::Heading(); printf("Jump database cache statistics:\n"); CConsole::Reset();
		printf("====================================\n");

		uint64_t cacheHits = 0, cacheMisses = 0, cacheTotal = 0;
		CJumpDnaHash* pJump = (CJumpDnaHash*)mpDNAHash;
		pJump->GetCacheStatistics(cacheHits, cacheMisses);

		cacheTotal = cacheHits + cacheMisses;
		double cacheHitsPercent = cacheHits / (double)cacheTotal * 100.0;

		printf("cache hits:   %10llu (%5.1f %%)\n", (unsigned long long)cacheHits, cacheHitsPercent);
		printf("cache misses: %10llu\n", (unsigned long long)cacheMisses);
	}

	fflush(stdout);

}


// estimates the appropriate hash table size
unsigned char CMosaikAligner::CalculateHashTableSize(const unsigned int referenceLength, const unsigned char hashSize) {

	// define our regression constants
	const double TOP_X_INTERCEPT   =  9.75487079E-01;
	const double TOP_X             = -3.23233529E-10;
	const double TOP_X2            =  9.75872364E-20;

	const double HALF_X_INTERCEPT  = -5.17783119E-01;
	const double HALF_LNX          =  7.65128834E-01;

	const double SLOPE_X_INTERCEPT =  8.64727395E-01;
	const double SLOPE_X           =  6.77627395E-10;
	const double SLOPE_X2          = -2.26963233E-19;

	// calculate the top, slope, and halfway point
	double refLen   = (double)referenceLength;
	double refLen2  = refLen * refLen;
	double lnRefLen = log(refLen);

	double top   = refLen2 * TOP_X2   + refLen * TOP_X   + TOP_X_INTERCEPT;
	double slope = refLen2 * SLOPE_X2 + refLen * SLOPE_X + SLOPE_X_INTERCEPT;
	double half  = lnRefLen * HALF_LNX + HALF_X_INTERCEPT;

	// calculate the bit size (Boltzmann sigmoid function)
	double numHashes80 = top / (1 + exp((half - hashSize) / slope)) * refLen / 0.8;
	unsigned char bitSize  = (unsigned char)ceil(log(numHashes80) / log(2.0));

	// always allocate at least 2 million entries (roughly 25 MB memory for fast) in our hash table
	if(bitSize < 21) bitSize = 21;

	return bitSize;
}

void CMosaikAligner::SetPeNeuralNetworkFilename( const string& neuralNetworkFilename ) {
	mPeNeuralNetworkFilename = neuralNetworkFilename;
}

void CMosaikAligner::SetSeNeuralNetworkFilename( const string& neuralNetworkFilename ) {
	mSeNeuralNetworkFilename = neuralNetworkFilename;
}


// enables special references checker
void CMosaikAligner::EnableSpecialReference ( const string referencePrefix ) {
	mSReference.enable = true;
	mSReference.prefix = referencePrefix;
}
// sets special hashes percentage
void CMosaikAligner::SetSpecialHashCount ( const unsigned int count ) {
	mSReference.count = count;
}

// outputs multiply mapped alignments
void CMosaikAligner::OutputMultiply( void ) {
	mFlags.OutputMultiply = true;
}

// sets quiet mode
void CMosaikAligner::SetQuietMode( void ) {
	mFlags.IsQuietMode = true;
}

// enables the alignment candidate threshold
void CMosaikAligner::EnableAlignmentCandidateThreshold(const unsigned short alignmentCandidateThreshold) {
	mFlags.IsUsingAlignmentCandidateThreshold   = true;
	mSettings.AlignmentCandidateThreshold = alignmentCandidateThreshold;
}

// enables the banded Smith-Waterman algorithm
void CMosaikAligner::EnableBandedSmithWaterman(const unsigned int bandwidth) {
	mFlags.UseBandedSmithWaterman = true;
	mSettings.Bandwidth           = bandwidth;
}

// enable the low-memory algorithm
void CMosaikAligner::EnableLowMemory(void) {
	mFlags.UseLowMemory  = true;
}

// Enables SOLiD colorspace translation
void CMosaikAligner::EnableColorspace(const string& basespaceReferenceFilename) {
	mFlags.EnableColorspace              = true;
	mSettings.BasespaceReferenceFilename = basespaceReferenceFilename;
}

// enables the hash position threshold
void CMosaikAligner::EnableHashPositionThreshold(const unsigned short hashPositionThreshold) {
	mFlags.IsUsingHashPositionThreshold = true;
	mSettings.HashPositionThreshold = hashPositionThreshold;
}

// enable the hash region threshold
void CMosaikAligner::EnableHashRegionThreshold(const unsigned short hashRegionThreshold) {
	mFlags.IsUsingHashRegionThreshold = true;
	mSettings.HashRegionThreshold = hashRegionThreshold;
}

// enables the use of the jump database
void CMosaikAligner::EnableJumpDB(const string& filenameStub, const unsigned int numCachedHashes, const bool keepKeysInMemory, const bool keepPositionsInMemory) {
	mFlags.IsUsingJumpDB             = true;
	mFlags.KeepJumpKeysInMemory      = keepKeysInMemory;
	mFlags.KeepJumpPositionsInMemory = keepPositionsInMemory;
	mSettings.JumpFilenameStub       = filenameStub;
	mSettings.NumCachedHashes        = numCachedHashes;
}

// enables the local alignment search
void CMosaikAligner::EnableLocalAlignmentSearch(const unsigned int radius) {
	mFlags.UseLocalAlignmentSearch       = true;
	mSettings.LocalAlignmentSearchRadius = radius;
}

// enables paired-end read output
void CMosaikAligner::EnablePairedEndOutput(void) {
	mFlags.UsePairedEndOutput = true;
}

// sets mapping quality threshold for stat map
void CMosaikAligner::SetStatMappingQuality( const unsigned char mq ) {
	mStatisticsCounters.StatMappingQuality = mq;
}

// enables reporting of unaligned reads
//void CMosaikAligner::EnableUnalignedReadReporting(const string& unalignedReadReportFilename) {
//	mSettings.UnalignedReadReportFilename = unalignedReadReportFilename;
//	mFlags.IsReportingUnalignedReads = true;
//}

// hashes the reference sequence
void CMosaikAligner::HashReferenceSequence(MosaikReadFormat::CReferenceSequenceReader& refseq) {

	// retrieve the 2-bit reference sequence and associated masking sequence
	char* twoBitConcatenatedSequence = NULL;
	char* maskSequence               = NULL;
	unsigned int numMaskedPositions;

	refseq.Load2BitConcatenatedSequence(twoBitConcatenatedSequence, maskSequence, numMaskedPositions);
	unsigned int numBases = refseq.GetReferenceSequenceLength();

	// initialization
	unsigned char rightMasks[4], rightShifts[4];
	unsigned char leftMasks[4] = { 0xff, 0x3f, 0x0f, 0x03 };

	char* pAnchor = twoBitConcatenatedSequence;

	// DEBUG by WAN-Ping on 20100920
	//unsigned int counter = 1;
	//for ( unsigned int i = 0; i < numBases / 4; i++ ) {
	//	char temp = *pAnchor;
	//	for ( unsigned int j = 0; j < 4; j++ ) {
	//		char temp1 = temp;
	//		unsigned short shiftBit = 6 - (j * 2);
	//		temp1 = temp1 >> shiftBit;
	//		temp1 &= 0x03;
	//		switch(temp1) {
	//			case 0: cout << 'A'; break;
	//			case 1: cout << 'C'; break;
	//			case 2: cout << 'G'; break;
	//			case 3: cout << 'T'; break;
	//			default: cout << 'N'; break;
	//		}
	//		if ( (counter % 70) == 0 )
	//			cout << endl;
	//		counter++;
	//
	//	}
	//	pAnchor++;
	//}

	//exit(1);
	// END od DEBUG

	uint64_t key, tKey;
	unsigned char hashSize = mSettings.HashSize;
	unsigned char hashBits = hashSize * 2;

	unsigned char observedBytes[4];
	for(unsigned char i = 0; i < 4; i++) {
		unsigned char sumBits = (i * 2) + hashBits;
		observedBytes[i]      = (unsigned char)ceil(sumBits / 8.0);
	}

	unsigned char unusedRightBits = observedBytes[0] * 8 - hashBits;

	switch(unusedRightBits) {
	case 0:
		rightMasks[0]  = 0xff; // 1111 1111
		rightMasks[1]  = 0xc0; // 1100 0000
		rightMasks[2]  = 0xf0; // 1111 0000
		rightMasks[3]  = 0xfc; // 1111 1100
		break;
	case 2:
		rightMasks[0]  = 0xfc; // 1111 1100
		rightMasks[1]  = 0xff; // 1111 1111
		rightMasks[2]  = 0xc0; // 1100 0000
		rightMasks[3]  = 0xf0; // 1111 0000
		break;
	case 4:
		rightMasks[0]  = 0xf0; // 1111 0000
		rightMasks[1]  = 0xfc; // 1111 1100
		rightMasks[2]  = 0xff; // 1111 1111
		rightMasks[3]  = 0xc0; // 1100 0000
		break;
	case 6:
		rightMasks[0]  = 0xc0; // 1100 0000
		rightMasks[1]  = 0xf0; // 1111 0000
		rightMasks[2]  = 0xfc; // 1111 1100
		rightMasks[3]  = 0xff; // 1111 1111
		break;
	default:
		cout << "ERROR: Unknown unused right bit combination." << endl;
		exit(1);
		break;
	}

	short rightShift = unusedRightBits;
	for(unsigned char i = 0; i < 4; i++) {
		rightShifts[i] = (unsigned char)rightShift;
		rightShift -= 2;
		if(rightShift < 0) rightShift = 6;
	}

	unsigned int currentByte;
	unsigned char cycle = 0;
	unsigned char shift;
	unsigned int numCompletedBytes = 0;

	// initialize masking variables
	unsigned int* pMask = (unsigned int*)maskSequence;
	unsigned int maskBegin = 0xffffffff;
	unsigned int maskEnd   = 0xffffffff;
	unsigned int maskIndex = 0;
	unsigned int maskMaxIndex = numMaskedPositions * 2;

	if(numMaskedPositions > 0) {
		maskBegin = pMask[maskIndex++];
		maskEnd   = pMask[maskIndex++];
	}

	unsigned int j = 0;
	unsigned int maxPositions = (unsigned int)(numBases - hashSize + 1);

	CConsole::Heading();
	cout << endl << "Hashing reference sequence:" << endl;
	CConsole::Reset();
	if ( !mFlags.IsQuietMode )
		CProgressBar<unsigned int>::StartThread(&j, 0, maxPositions, "ref bases");

	//unsigned int counter = 1;
	for(; j < maxPositions; j++) {

		// update mask
		if(j > maskEnd) {
			if(maskIndex == maskMaxIndex) {
				maskBegin = 0xffffffff;
				maskEnd   = 0xffffffff;
			} else {
				maskBegin = pMask[maskIndex++];
				maskEnd   = pMask[maskIndex++];
			}
		}

		// initialize
		currentByte = numCompletedBytes + observedBytes[cycle] - 1;
		shift       = -rightShifts[cycle];

		// add the rightmost byte
		tKey = pAnchor[currentByte--] & rightMasks[cycle];
		key  = tKey >> -shift;
		shift += 8;

		// add the internal bytes
		while(currentByte > numCompletedBytes) {
			tKey = pAnchor[currentByte--] & 0xff;
			key |= tKey << shift;
			shift += 8;
		}

		// add the leftmost byte
		tKey = pAnchor[currentByte] & leftMasks[cycle++];
		key |= tKey << shift;

		// DEBUG
		// show the current hash
		//for ( unsigned int k = 0; k < hashSize; k++ ) {
		//	uint64_t tempKey = key;
		//	//unsigned short shiftBit = hashSize * 2 - ( k + 1 ) * 2;
		//	//char currentBase = tempKey >> shiftBit;
		//	char currentBase = tempKey >> (hashSize * 2 - 2);
		//	currentBase &= 0x03;
		//	switch(currentBase) {
		//		case 0: cout << 'A'; break;
		//		case 1: cout << 'C'; break;
		//		case 2: cout << 'G'; break;
		//		case 3: cout << 'T'; break;
		//		default: cout << 'N'; break;
		//	}
		//}
		//if ( ( counter % 70) == 0 )
		//	cout << endl;
		//counter++;
		//cout << endl;
		//END of DEBUG

		// increment our offset cycle and 2-bit position
		if(cycle > 3) {
			cycle = 0;
			numCompletedBytes++;
		}

		// check if we should mask this
		bool maskPosition = false;
		unsigned int jEnd = j + hashSize - 1;
		if((j >= maskBegin) && (j <= maskEnd))      maskPosition = true;
		if((maskBegin >= j) && (maskBegin <= jEnd)) maskPosition = true;
		if((maskBegin >= j) && (maskEnd <= jEnd))   maskPosition = true;
		if((j >= maskBegin) && (jEnd <= maskEnd))   maskPosition = true;

		if(maskPosition) continue;

		mpDNAHash->Add(key, j);
	}

	if ( !mFlags.IsQuietMode )
		CProgressBar<unsigned int>::WaitThread();
	cout << endl;

	// clean up
	delete [] twoBitConcatenatedSequence;
	delete [] maskSequence;
}

// initializes the hash tables
void CMosaikAligner::InitializeHashTables(const unsigned char bitSize, const unsigned int begin, const unsigned int end, const unsigned int offset, const bool useLowMemory, const unsigned int expectedMemory, const bool bubbleSpecialHashes) {

	// decide which DNA hash table to use
	switch(mAlgorithm) {
	case CAlignmentThread::AlignerAlgorithm_FAST:
	case CAlignmentThread::AlignerAlgorithm_SINGLE:
		if(mFlags.IsUsingJumpDB) {
			mpDNAHash = new CJumpDnaHash(mSettings.HashSize, mSettings.JumpFilenameStub, 1, mFlags.KeepJumpKeysInMemory, mFlags.KeepJumpPositionsInMemory, mSettings.NumCachedHashes, begin, end, offset, expectedMemory, useLowMemory, bubbleSpecialHashes, mSReference.begin, mSReference.count);
		} else mpDNAHash = new CDnaHash(bitSize, mSettings.HashSize);
		break;
	case CAlignmentThread::AlignerAlgorithm_MULTI:
		if(mFlags.IsUsingJumpDB) {
			mpDNAHash = new CJumpDnaHash(mSettings.HashSize, mSettings.JumpFilenameStub, 9, mFlags.KeepJumpKeysInMemory, mFlags.KeepJumpPositionsInMemory, mSettings.NumCachedHashes, begin, end, offset, expectedMemory, useLowMemory, bubbleSpecialHashes, mSReference.begin, mSReference.count);
		} else mpDNAHash = new CMultiDnaHash(bitSize, mSettings.HashSize);
		break;
	case CAlignmentThread::AlignerAlgorithm_ALL:
		if(mFlags.IsUsingJumpDB) {
			mpDNAHash = new CJumpDnaHash(mSettings.HashSize, mSettings.JumpFilenameStub, 0, mFlags.KeepJumpKeysInMemory, mFlags.KeepJumpPositionsInMemory, mSettings.NumCachedHashes, begin, end, offset, expectedMemory, useLowMemory, bubbleSpecialHashes, mSReference.begin, mSReference.count);
		} else mpDNAHash = new CUbiqDnaHash(bitSize, mSettings.HashSize);
		break;
	default:
		cout << "ERROR: Unknown alignment algorithm specified." << endl;
		exit(1);
		break;
	}
}

// sets the filenames used by the aligner
void CMosaikAligner::SetFilenames(const string& inputReadArchiveFilename, const string& outputReadArchiveFilename, const string& referenceFilename) {
	mSettings.InputReadArchiveFilename  = inputReadArchiveFilename;
	mSettings.OutputReadArchiveFilename = outputReadArchiveFilename;
	mSettings.ReferenceFilename         = referenceFilename;
}

// enables the use of the entire read length when calculating mismatches
void CMosaikAligner::UseAlignedReadLengthForMismatchCalculation(void) {
	mFlags.UseAlignedReadLengthForMismatchCalculation = true;
}

void CMosaikAligner::SetLocalAlignmentSearchMqThreshold ( const unsigned char LocalAlignmentSearchHighMqThreshold, const unsigned char LocalAlignmentSearchLowMqThreshold ) {
	mFilters.LocalAlignmentSearchHighMqThreshold = LocalAlignmentSearchHighMqThreshold;
	mFilters.LocalAlignmentSearchLowMqThreshold  = LocalAlignmentSearchLowMqThreshold;
}
