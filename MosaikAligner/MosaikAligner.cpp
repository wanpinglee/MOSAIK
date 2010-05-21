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

// aligns the read archive
void CMosaikAligner::AlignReadArchive(void) {

	// ==============
	// initialization
	// ==============

	// retrieve the concatenated reference sequence length
	vector<ReferenceSequence> referenceSequences;

	MosaikReadFormat::CReferenceSequenceReader refseq;
	refseq.Open(mSettings.ReferenceFilename);
	refseq.GetReferenceSequences(referenceSequences);
	mReferenceLength = refseq.GetReferenceSequenceLength();
	const unsigned int numRefSeqs = refseq.GetNumReferenceSequences();

	// retrieve the basespace reference filenames
	char** pBsRefSeqs = NULL;
	if(mFlags.EnableColorspace) {

		cout << "- loading basespace reference sequences... ";
		cout.flush();

		MosaikReadFormat::CReferenceSequenceReader bsRefSeq;
		bsRefSeq.Open(mSettings.BasespaceReferenceFilename);

		if(!bsRefSeq.HasSameReferenceSequences(referenceSequences)) {
			printf("ERROR: The basespace and colorspace reference sequence archives do not seem to represent the same FASTA file.\n"); 
			exit(1);
		}

		bsRefSeq.CopyReferenceSequences(pBsRefSeqs);
		bsRefSeq.Close();

		cout << "finished." << endl;
	}

	// initialize our hash tables
	InitializeHashTables(CalculateHashTableSize(mReferenceLength, mSettings.HashSize));

	// hash the concatenated reference sequence
	if(!mFlags.IsUsingJumpDB) HashReferenceSequence(refseq);

	cout << "- loading reference sequence... ";
	cout.flush();
	refseq.LoadConcatenatedSequence(mReference);
	cout << "finished." << endl;

	refseq.Close();

	// create our reference sequence LUTs
	unsigned int* pRefBegin = new unsigned int[numRefSeqs];
	unsigned int* pRefEnd   = new unsigned int[numRefSeqs];

	for(unsigned int j = 0; j < numRefSeqs; j++) {
		pRefBegin[j] = referenceSequences[j].Begin;
		pRefEnd[j]   = referenceSequences[j].End;
	}

	// set the hash positions threshold
	if(mFlags.IsUsingHashPositionThreshold && (mAlgorithm == CAlignmentThread::AlignerAlgorithm_ALL)) 
		mpDNAHash->RandomizeAndTrimHashPositions(mSettings.HashPositionThreshold);

	// localize the read archive filenames
	string inputReadArchiveFilename  = mSettings.InputReadArchiveFilename;
	string outputReadArchiveFilename = mSettings.OutputReadArchiveFilename;

	// define our read format reader and writer
	MosaikReadFormat::CReadReader in;
	in.Open(inputReadArchiveFilename);
	MosaikReadFormat::ReadGroup readGroup = in.GetReadGroup();
	ReadStatus readStatus          = in.GetStatus();
	mSettings.SequencingTechnology = readGroup.SequencingTechnology;
	mSettings.MedianFragmentLength = readGroup.MedianFragmentLength;

	const bool isPairedEnd = (readStatus == RS_PAIRED_END_READ ? true : false);

	vector<MosaikReadFormat::ReadGroup> readGroups;
	readGroups.push_back(readGroup);

	// set the alignment status flags
	AlignmentStatus alignmentStatus = AS_UNSORTED_READ | readStatus;
	if(mMode == CAlignmentThread::AlignerMode_ALL) alignmentStatus |= AS_ALL_MODE;
	else alignmentStatus |= AS_UNIQUE_MODE;

	MosaikReadFormat::CAlignmentWriter out;
	out.Open(mSettings.OutputReadArchiveFilename.c_str(), referenceSequences, readGroups, alignmentStatus);

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
	td.pDnaHash            = mpDNAHash;
	td.pIn                 = &in;
	td.pOut                = &out;
	td.pUnalignedStream    = unalignedStream;
	td.pRefBegin           = pRefBegin;
	td.pRefEnd             = pRefEnd;
	td.Settings            = mSettings;
	td.pReadCounter        = &readCounter;
	td.IsPairedEnd         = isPairedEnd;
	td.pBsRefSeqs          = pBsRefSeqs;

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
	CConsole::Heading();
	cout << endl << "Aligning read library (" << numReadArchiveReads << "):" << endl;
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
	delete [] mReference;
	delete [] activeThreads;
	if(pRefBegin) delete [] pRefBegin;
	if(pRefEnd)   delete [] pRefEnd;

	if(pBsRefSeqs) {
		for(unsigned int i = 0; i < numRefSeqs; ++i) delete [] pBsRefSeqs[i];
		delete [] pBsRefSeqs;
	}

	// close open file streams
	in.Close();
	
	// solid references should be one-base longer after converting back to basespace
	if(mFlags.EnableColorspace) out.AdjustSolidReferenceBases();
	out.Close();

	if(mFlags.IsReportingUnalignedReads) fclose(unalignedStream);
	if(mFlags.IsUsingJumpDB) mpDNAHash->FreeMemory();

	// ====================
	// print our statistics
	// ====================

	const uint64_t totalMates = mStatisticsCounters.ShortMates +
		mStatisticsCounters.FailedHashMates +
		mStatisticsCounters.UniqueMates +
		mStatisticsCounters.NonUniqueMates +
		mStatisticsCounters.FilteredOutMates;

	const uint64_t totalAlignedMates = mStatisticsCounters.UniqueMates + mStatisticsCounters.NonUniqueMates;
	const uint64_t totalAlignedReads = mStatisticsCounters.AlignedReads;

	// print our alignment statistics (mates)
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
	if(mFlags.UseLocalAlignmentSearch) {
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
	if(mFlags.IsUsingJumpDB && (mSettings.NumCachedHashes > 0)) {
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

	printf("\n");
	CConsole::Heading(); printf("Miscellaneous statistics:\n"); CConsole::Reset();
	printf("==================================\n");
	printf("aligned mate bp:        %10llu\n", (unsigned long long)mStatisticsCounters.MateBasesAligned);
	printf("alignment candidates/s: %10.1f\n", mStatisticsCounters.AlignmentCandidates / alignmentBench.GetElapsedWallTime());
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
		rightMasks[0]  = 0xff;
		rightMasks[1]  = 0xc0;
		rightMasks[2]  = 0xf0;
		rightMasks[3]  = 0xfc;
		break;
	case 2:
		rightMasks[0]  = 0xfc;
		rightMasks[1]  = 0xff;
		rightMasks[2]  = 0xc0;
		rightMasks[3]  = 0xf0;
		break;
	case 4:
		rightMasks[0]  = 0xf0;
		rightMasks[1]  = 0xfc;
		rightMasks[2]  = 0xff;
		rightMasks[3]  = 0xc0;
		break;
	case 6:
		rightMasks[0]  = 0xc0;
		rightMasks[1]  = 0xf0;
		rightMasks[2]  = 0xfc;
		rightMasks[3]  = 0xff;
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
void CMosaikAligner::InitializeHashTables(const unsigned char bitSize) {

	// decide which DNA hash table to use
	switch(mAlgorithm) {
	case CAlignmentThread::AlignerAlgorithm_FAST:
	case CAlignmentThread::AlignerAlgorithm_SINGLE:
		if(mFlags.IsUsingJumpDB) {
			mpDNAHash = new CJumpDnaHash(mSettings.HashSize, mSettings.JumpFilenameStub, 1, mFlags.KeepJumpKeysInMemory, mFlags.KeepJumpPositionsInMemory, mSettings.NumCachedHashes);
		} else mpDNAHash = new CDnaHash(bitSize, mSettings.HashSize);
		break;
	case CAlignmentThread::AlignerAlgorithm_MULTI:
		if(mFlags.IsUsingJumpDB) {
			mpDNAHash = new CJumpDnaHash(mSettings.HashSize, mSettings.JumpFilenameStub, 9, mFlags.KeepJumpKeysInMemory, mFlags.KeepJumpPositionsInMemory, mSettings.NumCachedHashes);
		} else mpDNAHash = new CMultiDnaHash(bitSize, mSettings.HashSize);
		break;
	case CAlignmentThread::AlignerAlgorithm_ALL:
		if(mFlags.IsUsingJumpDB) {
			mpDNAHash = new CJumpDnaHash(mSettings.HashSize, mSettings.JumpFilenameStub, 0, mFlags.KeepJumpKeysInMemory, mFlags.KeepJumpPositionsInMemory, mSettings.NumCachedHashes);
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
