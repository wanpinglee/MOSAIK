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


// constructor
CMosaikAligner::CMosaikAligner(unsigned char hashSize, CAlignmentThread::AlignerAlgorithmType algorithmType, CAlignmentThread::AlignerModeType algorithmMode, unsigned char numThreads)
	: mAlgorithm(algorithmType)
	, mMode(algorithmMode)
	, mReference(NULL)
	, mReferenceLength(0)
	, mpDNAHash(NULL)
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

	// get reference information
	MosaikReadFormat::CReferenceSequenceReader refseq;
	refseq.Open(mSettings.ReferenceFilename);
	refseq.GetReferenceSequences(referenceSequences);
	mReferenceLength = refseq.GetReferenceSequenceLength();
	const unsigned int numRefSeqs = refseq.GetNumReferenceSequences();
	refseq.Close();

	// retrieve the basespace reference filenames
	if(mFlags.EnableColorspace) {
		
		MosaikReadFormat::CReferenceSequenceReader bsRefSeq;
		bsRefSeq.Open(mSettings.BasespaceReferenceFilename);

		if(!bsRefSeq.HasSameReferenceSequences(referenceSequences)) {
			printf("ERROR: The basespace and colorspace reference sequence archives do not seem to represent the same FASTA file.\n"); 
			exit(1);
		}

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
	vector<ReferenceSequence>::iterator ite1 = referenceSequences.begin();
	for ( unsigned int i = 0; i < nReference; ++i ) {
		referenceSequencesWoSpecial.push_back( *ite1 );
		++ite1;
	}
	
	// both of full- and low-memory MOSAIK need multiply-mapped bam in AlignmentThread.cpp
	mBams.mHeader.SortOrder = SORTORDER_UNSORTED;
	mBams.mHeader.pReferenceSequences = &referenceSequencesWoSpecial;
	mBams.mHeader.pReadGroups = &readGroups;
	mBams.mBam.Open( mSettings.OutputReadArchiveFilename + ".multiple.bam", mBams.mHeader);

	if ( !mFlags.UseLowMemory ) {
		
		// ==============================
		// set the headers of bam writers
		// ==============================
		//mBams.mHeader.SortOrder = SORTORDER_UNSORTED;
		mBams.sHeader.SortOrder = SORTORDER_UNSORTED;
		mBams.uHeader.SortOrder = SORTORDER_UNSORTED;
		mBams.rHeader.SortOrder = SORTORDER_UNSORTED;

		//mBams.mHeader.pReferenceSequences = &referenceSequencesWoSpecial;
		mBams.sHeader.pReferenceSequences = &referenceSequences;
		mBams.uHeader.pReferenceSequences = &referenceSequencesWoSpecial;
		mBams.rHeader.pReferenceSequences = &referenceSequencesWoSpecial;
	
		//mBams.mHeader.pReadGroups = &readGroups;
		mBams.sHeader.pReadGroups = &readGroups;
		mBams.uHeader.pReadGroups = &readGroups;
		mBams.rHeader.pReadGroups = &readGroups;

		//mBams.mBam.Open( mSettings.OutputReadArchiveFilename + ".multiple.bam", mBams.mHeader);
		mBams.sBam.Open( mSettings.OutputReadArchiveFilename + ".special.bam", mBams.sHeader);
		mBams.uBam.Open( mSettings.OutputReadArchiveFilename + ".unaligned.bam", mBams.uHeader);
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
		mBams.uBam.Close();
		mBams.rBam.Close();
	
	}
	else {

//***************** DEBUG ********************************
/*
	vector< string > temporaryFiles;
	temporaryFiles.push_back("tmp/ao504h7jbxfxr3hfmzspx6jug5ljcs5s.tmp");
	temporaryFiles.push_back("tmp/iuhqvma8ejonntt5ogpqsg7onkr91yc9.tmp");
	temporaryFiles.push_back("tmp/2vhgjba0ptm42bh8vbrev76ujy0s3oub.tmp");
	temporaryFiles.push_back("tmp/dofn2ecxd5oaizadu9bxnktheg35xib1.tmp");
	temporaryFiles.push_back("tmp/vrxf86h014fd9ky232etxbsk86b0qt4s.tmp");
	temporaryFiles.push_back("tmp/edccaeyanitd1x40qcsmu6j3r7bm8q5c.tmp");
	temporaryFiles.push_back("tmp/q1151muwpjti0scghwm9jvooypy49pd5.tmp");
	temporaryFiles.push_back("tmp/eh243s5v89r0mpz4oo9ef1n720uiz393.tmp");
	temporaryFiles.push_back("tmp/djdip5c7xkj72u33hpxsi96m7zm7pja9.tmp");
	temporaryFiles.push_back("tmp/0qzaje4lo4kg1z8hfsgwc8fgzur8tfmj.tmp");
	temporaryFiles.push_back("tmp/f6jz0n7o310ccernuc7cct55g3ecu90z.tmp");
	temporaryFiles.push_back("tmp/jwir9ysea08s978mc3tt2dgt4yp1sfth.tmp");
	temporaryFiles.push_back("tmp/3khy8mq7pplkja45gorpqw7iozajys07.tmp");
	temporaryFiles.push_back("tmp/qzstdph28dneek3yaqnvya18ek8ivter.tmp");
	temporaryFiles.push_back("tmp/f4vbvstfzladceyi8s3pjibnnl1wxk43.tmp");
	temporaryFiles.push_back("tmp/fbi2wpf71z4937b094ibwjpdeek2hwj2.tmp");
	temporaryFiles.push_back("tmp/3emt1cfob2wcmqrw595w6y20qeg5fu2o.tmp");
	temporaryFiles.push_back("tmp/jg7spl5m1uwky36vwrv5s63iwxf6ii47.tmp");
        
        // calculate total # of reads
        unsigned int nReads = 0;
        for ( unsigned int i = 0 ; i < temporaryFiles.size(); i++ ) {
	        MosaikReadFormat::CAlignmentReader reader;
                reader.Open( temporaryFiles[i] );
                nReads += reader.GetNumReads();
                reader.Close();
        }
	
	CConsole::Heading();
	cout << "Merging alignment archive:" << endl;
	CConsole::Reset();
	unsigned int readNo        = 0;
	//unsigned int nMaxAlignment = 1000;
        CProgressBar<unsigned int>::StartThread(&readNo, 0, nReads, "reads");
        CArchiveMerge merger( temporaryFiles, mSettings.OutputReadArchiveFilename, &readNo, mSettings.MedianFragmentLength, true );
        merger.Merge();
        CProgressBar<unsigned int>::WaitThread();
	exit(1);
*/
//********************************************************

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
				//cout << positionThreshold << endl;
				mpDNAHash->RandomizeAndTrimHashPositions(positionThreshold);
//cout << "p: " << positionThreshold << endl;
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
				smallReferenceSequences.push_back(referenceSequences[startRef+j]);
			}

			// define our output file
			mFlags.UseBamOutput     = false;
			mFlags.UseArchiveOutput = true;

			MosaikReadFormat::CAlignmentWriter out;
			out.Open(tempFilename.c_str(), smallReferenceSequences, readGroups, alignmentStatus, ALIGNER_SIGNATURE);
			out.AdjustPartitionSize(20000/referenceGroups.size());

			mFlags.SaveMultiplyBam = ( mSReference.found && ( i == referenceGroups.size() - 1 ) ) ? false : true;
			mFlags.SaveUnmappedBasesInArchive = ( i == 0 ) ? true : false; 
			AlignReadArchive(inn, out, pRefBegin, pRefEnd, pRefSpecies, pRefSpecial, pBsRefSeqs, referenceGroups[i].first );

			// close open file streams
			inn.Close();

			// solid references should be one-base longer after converting back to basespace
			if(mFlags.EnableColorspace) out.AdjustSolidReferenceBases();
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
		references.push_back(temp);
	}

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
        unsigned int nReads = 0;
	//unsigned int nOutputFilenames = mSReference.found ? outputFilenames.size() - 1 : outputFilenames.size() ;
        for ( unsigned int i = 0 ; i < outputFilenames.size(); i++ ) {
	        MosaikReadFormat::CAlignmentReader reader;
                reader.Open( outputFilenames[i] );
                nReads += reader.GetNumReads();
                reader.Close();
        }




	// if nThread is too large, it'll open too many files at the same time.
	// Then, we'll get an error since system doesn't allow us to open any file.
	if ( nThread > 7 )
		nThread = 7;

	// prepare reference offset vector for SOLiD
	//vector<unsigned int> refOffsets;
	//refOffsets.resize(referenceGroups.size());
	//for ( unsigned int i = 0; i < referenceGroups.size(); i++ ) {
	//	unsigned int startRef = referenceGroups[i].first;
	//	refOffsets[i] = referenceSequences[startRef].Begin;
	//}

	CConsole::Heading();
	cout << endl << "Sorting alignment archive:" << endl;
	CConsole::Reset();
	SortThread sThread ( outputFilenames, temporaryFiles, nThread, nReads, mSettings.MedianFragmentLength );
	sThread.Start();

	for ( unsigned int i = 0; i < outputFilenames.size(); i++ )
		rm(outputFilenames[i].c_str());


	CConsole::Heading();
	cout << "Merging alignment archive:" << endl;
	CConsole::Reset();

	//mStatisticsMaps.Reset();

        unsigned int readNo        = 0;
	//unsigned int nMaxAlignment = 1000;
        CProgressBar<unsigned int>::StartThread(&readNo, 0, nReads, "reads");
        CArchiveMerge merger( temporaryFiles, mSettings.OutputReadArchiveFilename, &readNo, mSettings.MedianFragmentLength, mSReference.found );
        merger.Merge();
        CProgressBar<unsigned int>::WaitThread();

	
	for ( unsigned int i = 0; i < temporaryFiles.size(); i++ )
		rm(temporaryFiles[i].c_str());

	// get statistics information
	string mapFile = mSettings.OutputReadArchiveFilename + ".stat";
	merger.PrintStatisticsMaps( mapFile, readGroups[0].ReadGroupID );

	CArchiveMerge::StatisticsCounters mergeCounters;
	merger.GetStatisticsCounters( mergeCounters );

	mStatisticsCounters.AlignedReads       = mergeCounters.AlignedReads;
	mStatisticsCounters.BothNonUniqueReads = mergeCounters.BothNonUniqueReads;
	mStatisticsCounters.BothUniqueReads    = mergeCounters.BothUniqueReads;
	mStatisticsCounters.OneNonUniqueReads  = mergeCounters.OneNonUniqueReads;
	mStatisticsCounters.OrphanedReads      = mergeCounters.OrphanedReads;
	mStatisticsCounters.FilteredOutMates   = mergeCounters.FilteredOutMates;
	mStatisticsCounters.NonUniqueMates     = mergeCounters.NonUniqueMates;
	mStatisticsCounters.UniqueMates        = mergeCounters.UniqueMates;

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
	FILE* unalignedStream = NULL;
	if(mFlags.IsReportingUnalignedReads) {
		if(fopen_s(&unalignedStream, mSettings.UnalignedReadReportFilename.c_str(), "wb") != 0) {
			cout << "ERROR: Unable to open the unaligned read FASTQ file for output." << endl;
			exit(1);
		}
	}
	


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
	td.pUnalignedStream    = unalignedStream;
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


	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	pthread_mutex_init(&CAlignmentThread::mGetReadMutex,              NULL);
	pthread_mutex_init(&CAlignmentThread::mReportUnalignedMate1Mutex, NULL);
	pthread_mutex_init(&CAlignmentThread::mReportUnalignedMate2Mutex, NULL);
	pthread_mutex_init(&CAlignmentThread::mSaveReadMutex,             NULL);
	pthread_mutex_init(&CAlignmentThread::mStatisticsMutex,           NULL);
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
	CProgressBar<uint64_t>::WaitThread();

	alignmentBench.Stop();

	// free up some memory
	delete [] activeThreads;
	activeThreads = NULL;

	if(mFlags.IsReportingUnalignedReads) fclose(unalignedStream);
	
}

// print our statistics
void CMosaikAligner::PrintStatistics () {

	// for low-memory version, the map is printed when merging archive.
	if ( !mFlags.UseLowMemory ) {
		string mapFile = mSettings.OutputReadArchiveFilename + ".stat";
		mStatisticsMaps.PrintMaps( mapFile.c_str(), readGroups[0].ReadGroupID.c_str() );
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
		mStatisticsCounters.UnalignedReads   = numReadArchiveReads - mStatisticsCounters.AlignedReads;
		mStatisticsCounters.FilteredOutMates = mStatisticsCounters.UnalignedReads;
	}

	const uint64_t totalAlignedMates = mStatisticsCounters.UniqueMates + mStatisticsCounters.NonUniqueMates;
	const uint64_t totalAlignedReads = mStatisticsCounters.AlignedReads;

	// print our alignment statistics (mates) if don't enable low-memory algorithm

	printf("\n");
	CConsole::Heading(); printf("Alignment statistics (mates):\n"); CConsole::Reset();
	printf("===================================\n");

	if(mStatisticsCounters.ShortMates > 0)
		printf("# too short:    %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.ShortMates,       (mStatisticsCounters.ShortMates       / (double)totalMates) * 100.0);

	if(mStatisticsCounters.FailedHashMates > 0)
		printf("# failed hash:  %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.FailedHashMates,  (mStatisticsCounters.FailedHashMates  / (double)totalMates) * 100.0);

	if(mStatisticsCounters.FilteredOutMates > 0)
		printf("# filtered out: %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.FilteredOutMates, (mStatisticsCounters.FilteredOutMates / (double)totalMates) * 100.0);

	if(mStatisticsCounters.UniqueMates > 0)
		printf("# unique:       %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.UniqueMates,      (mStatisticsCounters.UniqueMates      / (double)totalMates) * 100.0);

	if(mStatisticsCounters.NonUniqueMates > 0)
		printf("# non-unique:   %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.NonUniqueMates,   (mStatisticsCounters.NonUniqueMates   / (double)totalMates) * 100.0);

	printf("-----------------------------------\n");
	printf("total:          %9llu\n", (unsigned long long)totalMates);
	printf("total aligned:  ");
	CConsole::Bold(); printf("%9llu", (unsigned long long)totalAlignedMates); CConsole::Reset();
	printf(" (");
	CConsole::Bold(); printf("%5.1f %%", (totalAlignedMates / (double)totalMates) * 100.0); CConsole::Reset();
	printf(")\n");

	// print our local alignment search statistics
	// we don't print out local alignment information when the low-memory approach is enabled.
	if( !mFlags.UseLowMemory && mFlags.UseLocalAlignmentSearch ) {
		printf("\n");
		CConsole::Heading(); printf("Local alignment search statistics:\n"); CConsole::Reset();
		printf("===================================\n");

		double rescuedAlignmentsPercent = mStatisticsCounters.AdditionalLocalMates / (double)totalMates * 100.0;
		printf("rescued mates:  %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.AdditionalLocalMates, rescuedAlignmentsPercent);
	}

	// print our alignment statistics (reads)
	if(isPairedEnd) {
		printf("\n");
		CConsole::Heading(); printf("Alignment statistics (reads):\n"); CConsole::Reset();
		printf("============================================\n");
		printf("# unaligned:             %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.UnalignedReads,     (mStatisticsCounters.UnalignedReads     / (double)numReadArchiveReads) * 100.0);
		printf("# orphaned:              %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.OrphanedReads,      (mStatisticsCounters.OrphanedReads      / (double)numReadArchiveReads) * 100.0);
		printf("# both mates unique:     %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.BothUniqueReads,    (mStatisticsCounters.BothUniqueReads    / (double)numReadArchiveReads) * 100.0);
		printf("# one mate non-unique:   %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.OneNonUniqueReads,  (mStatisticsCounters.OneNonUniqueReads  / (double)numReadArchiveReads) * 100.0);
		printf("# both mates non-unique: %9llu (%5.1f %%)\n", (unsigned long long)mStatisticsCounters.BothNonUniqueReads, (mStatisticsCounters.BothNonUniqueReads / (double)numReadArchiveReads) * 100.0);
		printf("--------------------------------------------\n");
		printf("total reads:             ");
		CConsole::Bold(); printf("%9llu", (unsigned long long)numReadArchiveReads); CConsole::Reset();
		printf("\n");
		printf("total reads aligned:     ");
		CConsole::Bold(); printf("%9llu", (unsigned long long)totalAlignedReads); CConsole::Reset();
		printf(" (");
		CConsole::Bold(); printf("%5.1f %%", (totalAlignedReads / (double)numReadArchiveReads) * 100.0); CConsole::Reset();
		printf(")\n");
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

	//if ( !mFlags.UseLowMemory ) {
	//	printf("\n");
	//	CConsole::Heading(); printf("Miscellaneous statistics:\n"); CConsole::Reset();
	//	printf("==================================\n");
	//	printf("aligned mate bp:        %10llu\n", (unsigned long long)mStatisticsCounters.MateBasesAligned);
	//	printf("alignment candidates/s: %10.1f\n", mStatisticsCounters.AlignmentCandidates / alignmentBench.GetElapsedWallTime());
	//}
	
	
	
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
// enables special references checker
void CMosaikAligner::EnableSpecialReference ( const string referencePrefix ) {
	mSReference.enable = true;
	mSReference.prefix = referencePrefix;
}
// sets special hashes percentage
void CMosaikAligner::SetSpecialHashCount ( const unsigned int count ) {
	mSReference.count = count;
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

// enables reporting of unaligned reads
void CMosaikAligner::EnableUnalignedReadReporting(const string& unalignedReadReportFilename) {
	mSettings.UnalignedReadReportFilename = unalignedReadReportFilename;
	mFlags.IsReportingUnalignedReads = true;
}

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
