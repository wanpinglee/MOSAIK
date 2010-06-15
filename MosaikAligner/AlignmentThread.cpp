// ***************************************************************************
// CAlignmentThread - aligns all of the reads within a worker thread.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg & Wan-Ping Lee
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "AlignmentThread.h"

// register our thread mutexes
pthread_mutex_t CAlignmentThread::mGetReadMutex;
pthread_mutex_t CAlignmentThread::mReportUnalignedMate1Mutex;
pthread_mutex_t CAlignmentThread::mReportUnalignedMate2Mutex;
pthread_mutex_t CAlignmentThread::mSaveReadMutex;
pthread_mutex_t CAlignmentThread::mStatisticsMutex;

// define our constants
const double CAlignmentThread::P_ERR_REF  = pow(10.0, -REFERENCE_SEQUENCE_QUALITY / 10.0);
const double CAlignmentThread::P_CORR_REF = 1.0 - CAlignmentThread::P_ERR_REF;
const double CAlignmentThread::ONE_THIRD  = 1.0 / 3.0;
const double CAlignmentThread::TWO_NINTHS = 2.0 / 9.0;

// constructor
CAlignmentThread::CAlignmentThread(AlignerAlgorithmType& algorithmType, FilterSettings& filters, FlagData& flags, AlignerModeType& algorithmMode, char* pAnchor, unsigned int referenceLen, CAbstractDnaHash* pDnaHash, AlignerSettings& settings, unsigned int* pRefBegin, unsigned int* pRefEnd, char** pBsRefSeqs)
	: mAlgorithm(algorithmType)
	, mMode(algorithmMode)
	, mSettings(settings)
	, mFilters(filters)
	, mFlags(flags)
	, mReference(pAnchor)
	, mForwardRead(NULL)
	, mReverseRead(NULL)
	, mReferenceLength(referenceLen)
	, mpDNAHash(pDnaHash)
	, mSW(CPairwiseUtilities::MatchScore, CPairwiseUtilities::MismatchScore, CPairwiseUtilities::GapOpenPenalty, CPairwiseUtilities::GapExtendPenalty)
	, mBSW(CPairwiseUtilities::MatchScore, CPairwiseUtilities::MismatchScore, CPairwiseUtilities::GapOpenPenalty, CPairwiseUtilities::GapExtendPenalty, settings.Bandwidth)
	, mReferenceBegin(pRefBegin)
	, mReferenceEnd(pRefEnd)
{
	// calculate our base quality LUT
	for(unsigned char i = 0; i < 100; i++) mBaseQualityLUT[i] = pow(10.0, -i / 10.0);

	// set our flags
	if(algorithmMode == AlignerMode_ALL) mFlags.IsAligningAllReads = true;

	// assign the reference sequences to the colorspace utilities object
	mCS.SetReferenceSequences(pBsRefSeqs);
}

// destructor
CAlignmentThread::~CAlignmentThread(void) {
	if(mForwardRead) delete [] mForwardRead;
	if(mReverseRead) delete [] mReverseRead;
}

// activates the current alignment thread
void* CAlignmentThread::StartThread(void* arg) {
	ThreadData* pTD = (ThreadData*)arg;

	// align reads
	CAlignmentThread at(pTD->Algorithm, pTD->Filters, pTD->Flags, pTD->Mode, pTD->pReference, pTD->ReferenceLen, pTD->pDnaHash, pTD->Settings, pTD->pRefBegin, pTD->pRefEnd, pTD->pBsRefSeqs);
	at.AlignReadArchive(pTD->pIn, pTD->pOut, pTD->pUnalignedStream, pTD->pReadCounter, pTD->IsPairedEnd);

	vector<ReferenceSequence>::iterator refIter;

	// synchronize the statistics
	pthread_mutex_lock(&mStatisticsMutex);
	pTD->pCounters->AdditionalLocalMates += at.mStatisticsCounters.AdditionalLocalMates;
	pTD->pCounters->AlignedReads         += at.mStatisticsCounters.AlignedReads;
	pTD->pCounters->AlignmentCandidates  += at.mStatisticsCounters.AlignmentCandidates;
	pTD->pCounters->BothNonUniqueReads   += at.mStatisticsCounters.BothNonUniqueReads;
	pTD->pCounters->BothUniqueReads      += at.mStatisticsCounters.BothUniqueReads;
	pTD->pCounters->FailedHashMates      += at.mStatisticsCounters.FailedHashMates;
	pTD->pCounters->FilteredOutMates     += at.mStatisticsCounters.FilteredOutMates;
	pTD->pCounters->MateBasesAligned     += at.mStatisticsCounters.MateBasesAligned;
	pTD->pCounters->NonUniqueMates       += at.mStatisticsCounters.NonUniqueMates;
	pTD->pCounters->OneNonUniqueReads    += at.mStatisticsCounters.OneNonUniqueReads;
	pTD->pCounters->OrphanedReads        += at.mStatisticsCounters.OrphanedReads;
	pTD->pCounters->ShortMates           += at.mStatisticsCounters.ShortMates;
	pTD->pCounters->UnalignedReads       += at.mStatisticsCounters.UnalignedReads;
	pTD->pCounters->UniqueMates          += at.mStatisticsCounters.UniqueMates;
	pthread_mutex_unlock(&mStatisticsMutex);

	// exit the thread
	pthread_exit((void*)0);
	return 0;
}

// aligns the read archive
void CAlignmentThread::AlignReadArchive(MosaikReadFormat::CReadReader* pIn, MosaikReadFormat::CAlignmentWriter* pOut, FILE* pUnalignedStream, uint64_t* pReadCounter, bool isPairedEnd) {

	// create our local alignment models
	const bool isUsing454      = (mSettings.SequencingTechnology == ST_454      ? true : false);
	const bool isUsingIllumina = (mSettings.SequencingTechnology == ST_ILLUMINA ? true : false);
	const bool isUsingSOLiD    = (mSettings.SequencingTechnology == ST_SOLID    ? true : false);

	// catch unsupported local alignment search sequencing technologies
	if(mFlags.UseLocalAlignmentSearch && (!isUsing454 && !isUsingIllumina)) {
		cout << "ERROR: This sequencing technology is not currently supported for local alignment search." << endl;
		exit(1);
	}

	// initialize our status variables
	AlignmentStatusType mate1Status, mate2Status;

	// derive the minimum span length
	unsigned int minSpanLength = mSettings.HashSize;
	if(mFlags.IsUsingAlignmentCandidateThreshold && (mSettings.AlignmentCandidateThreshold > minSpanLength)) 
		minSpanLength = mSettings.AlignmentCandidateThreshold;

	// decide if we need to calculate the correction coefficient
	const bool calculateCorrectionCoefficient = mFlags.IsUsingJumpDB && mFlags.IsUsingHashPositionThreshold;

	// keep reading until no reads remain
	Mosaik::Read mr;
	CNaiveAlignmentSet mate1Alignments(mReferenceLength, (isUsingIllumina || isUsingSOLiD)), mate2Alignments(mReferenceLength, (isUsingIllumina || isUsingSOLiD));

	while(true) {

		mr.clear();
		pthread_mutex_lock(&mGetReadMutex);
		bool hasMoreReads = pIn->LoadNextRead(mr);
		if(hasMoreReads) *pReadCounter = *pReadCounter + 1;
		pthread_mutex_unlock(&mGetReadMutex);

		// quit if we've processed all of the reads
		if(!hasMoreReads) break;

		// specify if this is a paired-end read
		const unsigned short numMate1Bases = (unsigned short)mr.Mate1.Bases.Length();
		const unsigned short numMate2Bases = (unsigned short)mr.Mate2.Bases.Length();
		const bool areBothMatesPresent = (((numMate1Bases != 0) && (numMate2Bases != 0)) ? true : false);

		// ====================
		// align the first mate
		// ====================

		bool isMate1Aligned = false;
		mate1Alignments.Clear();
		if(numMate1Bases != 0) {

			// align the read
			if(AlignRead(mate1Alignments, mr.Mate1.Bases.CData(), mr.Mate1.Qualities.CData(), numMate1Bases, mate1Status)) {

				// calculate the alignment qualities
				mate1Alignments.CalculateAlignmentQualities(calculateCorrectionCoefficient, minSpanLength);
				isMate1Aligned = true;

			} else {

				// write the unaligned read if specified
				if(mFlags.IsReportingUnalignedReads) {
					pthread_mutex_lock(&mReportUnalignedMate1Mutex);
					mr.Mate1.Qualities.Increment(33);

					if(isUsingSOLiD) {
						mCS.ConvertReadPseudoColorspaceToColorspace(mr.Mate1.Bases);
						mr.Mate1.Bases.Prepend(mr.Mate1.SolidPrefixTransition, 2);
						mr.Mate1.Qualities.Prepend("!?", 2);
					}

					fprintf(pUnalignedStream, "@%s (mate 1, length=%u)\n%s\n+\n%s\n", mr.Name.CData(), 
						(isUsingSOLiD ? numMate1Bases + 1 : numMate1Bases), mr.Mate1.Bases.CData(), 
						mr.Mate1.Qualities.CData());

					pthread_mutex_unlock(&mReportUnalignedMate1Mutex);
				}
			}
		}

		// =====================
		// align the second mate
		// =====================

		bool isMate2Aligned = false;
		mate2Alignments.Clear();
		if(numMate2Bases != 0) {

			// align the read
			if(AlignRead(mate2Alignments, mr.Mate2.Bases.CData(), mr.Mate2.Qualities.CData(), numMate2Bases, mate2Status)) {

				// calculate the alignment qualities
				mate2Alignments.CalculateAlignmentQualities(calculateCorrectionCoefficient, minSpanLength);
				isMate2Aligned = true;

			} else {

				// write the unaligned read if specified
				if(mFlags.IsReportingUnalignedReads) {
					pthread_mutex_lock(&mReportUnalignedMate2Mutex);
					mr.Mate2.Qualities.Increment(33);

					if(isUsingSOLiD) {
						mCS.ConvertReadPseudoColorspaceToColorspace(mr.Mate2.Bases);
						mr.Mate2.Bases.Prepend(mr.Mate2.SolidPrefixTransition, 2);
						mr.Mate2.Qualities.Prepend("!?", 2);
					}

					fprintf(pUnalignedStream, "@%s (mate 2, length=%u)\n%s\n+\n%s\n", mr.Name.CData(), 
						(isUsingSOLiD ? numMate2Bases + 1 : numMate2Bases), mr.Mate2.Bases.CData(), 
						mr.Mate2.Qualities.CData());

					pthread_mutex_unlock(&mReportUnalignedMate2Mutex);
				}
			}
		}

		// ======================
		// local alignment search
		// ======================

		const bool isMate1Unique = mate1Alignments.IsUnique();
		const bool isMate2Unique = mate2Alignments.IsUnique();

		// we can only perform a local alignment search if both mates are present
		if(areBothMatesPresent) {

			// perform local alignment search WHEN MATE1 IS UNIQUE
			if(mFlags.UseLocalAlignmentSearch && isMate1Unique) {

				// extract the unique begin and end coordinates
				AlignmentSet::const_iterator uniqueIter = mate1Alignments.GetSet()->begin();
				const unsigned int refIndex    = uniqueIter->ReferenceIndex;
				const unsigned int uniqueBegin = mReferenceBegin[refIndex] + uniqueIter->ReferenceBegin;
				const unsigned int uniqueEnd   = mReferenceBegin[refIndex] + uniqueIter->ReferenceEnd;

				// create the appropriate local alignment search model
				LocalAlignmentModel lam;
				if(uniqueIter->IsReverseStrand) {
					if(isUsingIllumina) lam.IsTargetBeforeUniqueMate = true;
					else if(isUsing454) lam.IsTargetReverseStrand    = true;
					// NB: we caught other technologies during initialization
				} else {
					if(isUsingIllumina) lam.IsTargetReverseStrand    = true;
					else if(isUsing454) lam.IsTargetBeforeUniqueMate = true;
					// NB: we caught other technologies during initialization
				}

				Alignment al;
				if(RescueMate(lam, mr.Mate2.Bases, uniqueBegin, uniqueEnd, refIndex, al)) {

					const char* pQualities = mr.Mate2.Qualities.CData();

					// add the alignment to the alignment set if it passes the filters
					if(ApplyReadFilters(al, pQualities, mr.Mate2.Bases.Length())) {
						al.WasRescued = true;
						if(mate2Alignments.Add(al)) {
							mStatisticsCounters.AdditionalLocalMates++;
							mate2Status    = ALIGNMENTSTATUS_GOOD;
							isMate2Aligned = true;
						}
					}

					// increment our candidates counter
					mStatisticsCounters.AlignmentCandidates++;
				}
			}

			// perform local alignment search WHEN MATE2 IS UNIQUE
			if(mFlags.UseLocalAlignmentSearch && isMate2Unique) {

				// extract the unique begin and end coordinates
				AlignmentSet::const_iterator uniqueIter = mate2Alignments.GetSet()->begin();
				const unsigned int refIndex    = uniqueIter->ReferenceIndex;
				const unsigned int uniqueBegin = mReferenceBegin[refIndex] + uniqueIter->ReferenceBegin;
				const unsigned int uniqueEnd   = mReferenceBegin[refIndex] + uniqueIter->ReferenceEnd;

				// create the appropriate local alignment search model
				LocalAlignmentModel lam;
				if(uniqueIter->IsReverseStrand) {
					if(isUsingIllumina) lam.IsTargetBeforeUniqueMate = true;
					else if(isUsing454) {
						lam.IsTargetReverseStrand    = true;
						lam.IsTargetBeforeUniqueMate = true;
					}
					// NB: we caught other technologies during initialization
				} else {
					if(isUsingIllumina) lam.IsTargetReverseStrand = true;
					// NB: defaults are appropriate for 454
					// NB: we caught other technologies during initialization
				}

				Alignment al;
				if(RescueMate(lam, mr.Mate1.Bases, uniqueBegin, uniqueEnd, refIndex, al)) {

					const char* pQualities = mr.Mate1.Qualities.CData();

					// add the alignment to the alignment set if it passes the filters
					if(ApplyReadFilters(al, pQualities, mr.Mate1.Bases.Length())) {
						al.WasRescued = true;
						if(mate1Alignments.Add(al)) {
							mStatisticsCounters.AdditionalLocalMates++;
							mate1Status    = ALIGNMENTSTATUS_GOOD;
							isMate1Aligned = true;
						}
					}

					// increment our candidates counter
					mStatisticsCounters.AlignmentCandidates++;
				}
			}
		}

		// =================
		// update statistics
		// =================

		// update the alignment status statistics for mate1
		if(numMate1Bases != 0) {
			switch(mate1Status) {
			case ALIGNMENTSTATUS_TOOSHORT:
				mStatisticsCounters.ShortMates++;
				break;
			case ALIGNMENTSTATUS_FAILEDHASH:
				mStatisticsCounters.FailedHashMates++;
				break;
			case ALIGNMENTSTATUS_FILTEREDOUT:
				mStatisticsCounters.FilteredOutMates++;
				break;
			case ALIGNMENTSTATUS_GOOD:
				if(mate1Alignments.IsUnique()) mStatisticsCounters.UniqueMates++;
				else mStatisticsCounters.NonUniqueMates++;
				break;
			}
		}

		// update the alignment status statistics for mate2
		if(numMate2Bases != 0) {
			switch(mate2Status) {
			case ALIGNMENTSTATUS_TOOSHORT:
				mStatisticsCounters.ShortMates++;
				break;
			case ALIGNMENTSTATUS_FAILEDHASH:
				mStatisticsCounters.FailedHashMates++;
				break;
			case ALIGNMENTSTATUS_FILTEREDOUT:
				mStatisticsCounters.FilteredOutMates++;
				break;
			case ALIGNMENTSTATUS_GOOD:
				if(mate2Alignments.IsUnique()) mStatisticsCounters.UniqueMates++;
				else mStatisticsCounters.NonUniqueMates++;
				break;
			}
		}

		// update the paired-end read statistics
		if(isPairedEnd) {
			if(isMate1Aligned || isMate2Aligned) {
				if(isMate1Aligned && isMate2Aligned) {
					const bool isM1Unique = mate1Alignments.IsUnique();
					const bool isM2Unique = mate2Alignments.IsUnique();
					if(isM1Unique && isM2Unique) mStatisticsCounters.BothUniqueReads++;
					else if(!isM1Unique && !isM2Unique) mStatisticsCounters.BothNonUniqueReads++;
					else mStatisticsCounters.OneNonUniqueReads++;
				} else mStatisticsCounters.OrphanedReads++;
			} else mStatisticsCounters.UnalignedReads++;
		}

		// update the mate bases aligned
		if(isMate1Aligned) mStatisticsCounters.MateBasesAligned += numMate1Bases;
		if(isMate2Aligned) mStatisticsCounters.MateBasesAligned += numMate2Bases;

		// =================
		// save aligned read
		// =================

		// if any of the two mates aligned, save the read
		if(isMate1Aligned || isMate2Aligned) {
			mStatisticsCounters.AlignedReads++;
			pthread_mutex_lock(&mSaveReadMutex);
			pOut->SaveRead(mr, mate1Alignments, mate2Alignments);
			pthread_mutex_unlock(&mSaveReadMutex);
		}
	}
}

// aligns the read against the reference sequence and returns true if the read was aligned
bool CAlignmentThread::AlignRead(CNaiveAlignmentSet& alignments, const char* query, const char* qualities, const unsigned int queryLength, AlignmentStatusType& status) {

	// set the alignment status to good
	status = ALIGNMENTSTATUS_GOOD;

	// return if the read is smaller than the hash size
	if(queryLength < mSettings.HashSize) {
		status = ALIGNMENTSTATUS_TOOSHORT;
		return false;
	}

	// resize the forward and reverse reads if required
	if(queryLength >= mSettings.AllocatedReadLength) {

		// clean up
		if(mForwardRead) delete [] mForwardRead;
		if(mReverseRead) delete [] mReverseRead;

		// create a larger allocated read length
		mSettings.AllocatedReadLength = queryLength + ALLOCATION_EXTENSION;

		try {
			mForwardRead = new char[mSettings.AllocatedReadLength];
			mReverseRead = new char[mSettings.AllocatedReadLength];
		} catch(bad_alloc) {
			cout << "ERROR: Unable to allocate enough memory for the forward and reverse read buffers." << endl;
			exit(1);
		}
	}

	// calculate the number of bases to extend past the hash regions during alignment
	unsigned int numExtensionBases = 0;
	if(mFilters.UseMismatchFilter)        numExtensionBases = mFilters.MaxNumMismatches;
	if(mFilters.UseMismatchPercentFilter) numExtensionBases = (unsigned int)(queryLength * mFilters.MaxMismatchPercent);
	if(numExtensionBases < 2)             numExtensionBases = 2;

	// statistics variables
	int64_t hashRegionLength = 0;

	// control variables
	bool evaluateReverseReads = true;

	const bool alignAllReads    = mFlags.IsAligningAllReads;
	const bool useFastAlgorithm = (mAlgorithm == AlignerAlgorithm_FAST ? true : false);

	bool ret = true; // assume we will align the read

	try {

		//
		// copy the query to the forward and reverse reads
		//

		// N.B. using global reads prevents unnecessary initialization, but
		// only speeds up FAST algorithm performance by 1.4 %
		strncpy_s(mForwardRead, mSettings.AllocatedReadLength, query, queryLength);
		strncpy_s(mReverseRead, mSettings.AllocatedReadLength, query, queryLength);

		mForwardRead[queryLength] = 0;
		mReverseRead[queryLength] = 0;

		if(mFlags.EnableColorspace) CSequenceUtilities::ReverseSequence(mReverseRead, queryLength);
		else CSequenceUtilities::GetReverseComplement(mReverseRead, queryLength);

		// ======================================
		// consolidate the hashes for our regions
		// ======================================

		// used for all algorithms except fast
		vector<HashRegion> forwardRegions, reverseRegions;

		// used for fast algorithm
		HashRegion fastHashRegion;
		bool isFastHashRegionReverseStrand = false;
		char* fastHashRead = NULL;

		if(useFastAlgorithm) {

			// find the largest hash region in either direction and then pick
			// the largest of the two directions
			HashRegion forwardHashRegion, reverseHashRegion;
			int64_t forwardHashRegionLength = 0, reverseHashRegionLength = 0;
			int64_t* pHashRegionLength = NULL;

			GetFastReadCandidate(forwardHashRegion, mForwardRead, queryLength, alignments.GetFwdMhpOccupancyList());
			GetFastReadCandidate(reverseHashRegion, mReverseRead, queryLength, alignments.GetRevMhpOccupancyList());

			// detect failed hashes
			if((forwardHashRegion.End == 0) && (reverseHashRegion.End == 0)) {
				status = ALIGNMENTSTATUS_FAILEDHASH;
				return false;
			}

			forwardHashRegionLength = forwardHashRegion.End - forwardHashRegion.Begin + 1;
			reverseHashRegionLength = reverseHashRegion.End - reverseHashRegion.Begin + 1;

			if(reverseHashRegionLength > forwardHashRegionLength) {
				fastHashRegion    = reverseHashRegion;
				fastHashRead      = mReverseRead;
				pHashRegionLength = &reverseHashRegionLength;
				isFastHashRegionReverseStrand = true;
			} else {
				fastHashRegion    = forwardHashRegion;
				fastHashRead      = mForwardRead;
				pHashRegionLength = &forwardHashRegionLength;
			}

			// enforce double-hits
			if(*pHashRegionLength <= mSettings.HashSize) {
				status = ALIGNMENTSTATUS_FAILEDHASH;
				return false;
			}

		} else {

			GetReadCandidates(forwardRegions, mForwardRead, queryLength, alignments.GetFwdMhpOccupancyList());
			GetReadCandidates(reverseRegions, mReverseRead, queryLength, alignments.GetRevMhpOccupancyList());

			// detect failed hashes
			if(forwardRegions.empty() && reverseRegions.empty()) {
				status = ALIGNMENTSTATUS_FAILEDHASH;
				return false;
			}
		}

		// =======================
		// evaluate the fast reads
		// =======================

		if(useFastAlgorithm) {

			// create a new alignment data structure
			Alignment al;
			al.IsReverseStrand = isFastHashRegionReverseStrand;

			// perform a Smith-Waterman alignment
			AlignRegion(fastHashRegion, al, fastHashRead, queryLength, numExtensionBases);

			// add the alignment to the vector if it passes the filters
			if(ApplyReadFilters(al, qualities, queryLength)) alignments.Add(al);

			// increment our candidates counter
			mStatisticsCounters.AlignmentCandidates++;

		} else {

			for(unsigned int i = 0; i < (unsigned int)forwardRegions.size(); i++) {

				// enforce alignment candidate thresholds
				if(mFlags.IsUsingAlignmentCandidateThreshold) {
					hashRegionLength = forwardRegions[i].End - forwardRegions[i].Begin + 1;
					if(hashRegionLength < mSettings.AlignmentCandidateThreshold) continue;
				}

				// create a new alignment data structure
				Alignment al;
				//al.IsReverseStrand = false;

				// perform a Smith-Waterman alignment
				AlignRegion(forwardRegions[i], al, mForwardRead, queryLength, numExtensionBases);

				// add the alignment to the alignments vector
				if(ApplyReadFilters(al, qualities, queryLength)) alignments.Add(al);

				// increment our candidates counter
				mStatisticsCounters.AlignmentCandidates++;

				// check if we can prematurely stop
				if(!alignAllReads && (alignments.GetCount() > 1)) {
					evaluateReverseReads = false;				
					break;
				}
			}

			// ==========================
			// evaluate the reverse reads
			// ==========================

			if(alignAllReads) evaluateReverseReads = true;

			if(evaluateReverseReads) {

				for(unsigned int i = 0; i < (unsigned int)reverseRegions.size(); i++) {

					// enforce alignment candidate thresholds
					if(mFlags.IsUsingAlignmentCandidateThreshold) {
						hashRegionLength = reverseRegions[i].End - reverseRegions[i].Begin + 1;
						if(hashRegionLength < mSettings.AlignmentCandidateThreshold) continue;
					}

					// create a new alignment data structure
					Alignment al;
					al.IsReverseStrand = true;

					// perform a Smith-Waterman alignment
					AlignRegion(reverseRegions[i], al, mReverseRead, queryLength, numExtensionBases);

					// add the alignment to the alignments vector
					if(ApplyReadFilters(al, qualities, queryLength)) alignments.Add(al);

					// increment our candidates counter
					mStatisticsCounters.AlignmentCandidates++;

					// check if we can prematurely stop
					if(!alignAllReads && (alignments.GetCount() > 1)) break;
				}
			}
		}

		// no alignments because of filtering
		if(alignments.IsEmpty()) {
			status = ALIGNMENTSTATUS_FILTEREDOUT;
			ret = false;
		}

	} catch(bad_alloc &ba) {

		cout << "ERROR: Could not allocate enough memory to create forward and reverse aligned sequences: " << ba.what() << endl;
		exit(1);
	}

	return ret;
}

// aligns the read against a specified hash region using Smith-Waterman-Gotoh
void CAlignmentThread::AlignRegion(const HashRegion& r, Alignment& alignment, char* query, unsigned int queryLength, unsigned int extensionBases) {

	// define the begin coordinate of our alignment region
	unsigned int begin = r.End;
	if(begin >= queryLength) begin -= queryLength - 1;
	else begin = 0;

	if(r.Begin < begin) begin = r.Begin;

	if(begin >= extensionBases) begin -= extensionBases;
	else begin = 0;

	// define the end coordinate of our alignment region
	// N.B. built-in extension of 1
	unsigned int end  = r.Begin + queryLength - 1;
	if(r.End > end) end = r.End;

	end += extensionBases;


	// make sure the endpoints are within the reference sequence
	unsigned int referenceIndex = 0;
	while(r.Begin > mReferenceEnd[referenceIndex]) referenceIndex++;


	const unsigned int refBegin = mReferenceBegin[referenceIndex];
	const unsigned int refEnd   = mReferenceEnd[referenceIndex];

	if(begin < refBegin) begin = refBegin;
	if(end   < refBegin) end   = refBegin;
	if(begin > refEnd)   begin = refEnd;
	if(end   > refEnd)   end   = refEnd;

	// adjust the begin and end positions if the reference is masked
	while(mReference[begin] == 'X') begin++;
	while(mReference[end]   == 'X') end--;

	// perform a Smith-Waterman alignment on our region
	char* pAnchor = mReference + begin;

	
	if ( begin > mReferenceLength ) {
		cout << "ERROR: The hash region excceds the reference region." << endl;
		exit(1);
	}

	// determine if the specified bandwidth is enough to accurately align using the banded algorithm
	bool hasEnoughBandwidth = false;
	HashRegion diagonalRegion = r;

	if(mFlags.UseBandedSmithWaterman) {

		diagonalRegion.Begin -= begin;
		diagonalRegion.End   -= begin;

		unsigned int rowStart = min(diagonalRegion.Begin, (unsigned int)diagonalRegion.QueryBegin);

		diagonalRegion.Begin      -= rowStart;
		diagonalRegion.QueryBegin -= rowStart;

		hasEnoughBandwidth = (queryLength - diagonalRegion.QueryBegin) > mSettings.Bandwidth;
		hasEnoughBandwidth = hasEnoughBandwidth && (((end - begin + 1) - diagonalRegion.Begin) > mSettings.Bandwidth / 2);
	}

	if(mFlags.UseBandedSmithWaterman && hasEnoughBandwidth) {
		mBSW.Align(alignment, pAnchor, (end - begin + 1), query, queryLength, diagonalRegion);
	} else {
		mSW.Align(alignment, pAnchor, (end - begin + 1), query, queryLength);
	}

	// adjust the reference start positions
	//if ( !mFlags.UseLowMemory )
		alignment.ReferenceIndex = referenceIndex;
	//else
	//	alignment.ReferenceIndex = 0;
	alignment.ReferenceBegin += begin - refBegin;
	alignment.ReferenceEnd   += begin - refBegin;
}

// returns true if the alignment passes all of the user-specified filters
bool CAlignmentThread::ApplyReadFilters(Alignment& al, const char* qualities, const unsigned int queryLength) {

	// assuming this is a good read
	bool ret = true;

	// copy the base qualities
	// TODO: there may be a significant penalty involved in performing this on all aligned reads
	al.BaseQualities.Copy(qualities + al.QueryBegin, al.QueryEnd - al.QueryBegin + 1);
	if(al.IsReverseStrand) al.BaseQualities.Reverse();

	unsigned short numNonAlignedBases = queryLength - al.QueryLength;

	// convert from colorspace to basespace
	if(mFlags.EnableColorspace) {
		mCS.ConvertAlignmentToBasespace(al);		
		numNonAlignedBases = (queryLength + 1) - al.QueryLength;
	}

	// calculate the total number of mismatches
	const unsigned short numTotalMismatches = al.NumMismatches + (mFlags.UseAlignedReadLengthForMismatchCalculation ? 0 : numNonAlignedBases);

	// check to see if this alignment meets the maximum mismatch threshold
	if(mFilters.UseMismatchFilter && (numTotalMismatches > mFilters.MaxNumMismatches)) ret = false;

	// check to see if this alignment meets the maximum mismatch threshold
	if(mFilters.UseMismatchPercentFilter) {
		double percentMismatch = (double)numTotalMismatches / (double)al.QueryLength;
		if(percentMismatch > mFilters.MaxMismatchPercent) ret = false;
	}

	// check to see if this alignment meets the minimum percentage alignment threshold
	if(mFilters.UseMinAlignmentPercentFilter) {
		double percentageAligned = (double)al.QueryLength / (double)queryLength;
		if(percentageAligned < mFilters.MinPercentAlignment) ret = false;
	}

	// check to see if this alignment meets the minimum alignment threshold
	if(mFilters.UseMinAlignmentFilter && (al.QueryLength < mFilters.MinAlignment)) ret = false;

	return ret;
}

// creates the hash for a supplied fragment
void CAlignmentThread::CreateHash(const char* fragment, const unsigned char fragmentLen, uint64_t& key) {

	// set the key to zero
	key = 0;

	if(fragmentLen > 32) {
		cout << "ERROR: This hash table can only handle fragments smaller or equal to 32 bases." << endl;
		exit(1);
	}

	// human genome 36.2
	//
	// A: 843953565
	// C: 584268578
	// G: 584621685
	// M: 1
	// N: 129484
	// R: 2
	// T: 845168978
	char translation[26] = { 0, 3, 1, 3, -1, -1, 2, 3, -1, -1, 3, -1, 0, 3, -1, -1, -1, 0, 2, 3, -1, 0, 3, 1, 3, -1 };

	// convert each nucleotide to its 2-bit representation
	for(unsigned char i = 0; i < fragmentLen; i++) {

		// convert [A,C,G,T] to [0,1,2,3]
		char tValue = translation[fragment[i] - 'A'];

		// catch any unrecognized nucleotides
		if(tValue < 0) {
			cout << "ERROR: Unrecognized nucleotide in hash table: " << fragment[i] << endl;
			cout << "- fragment: ";
			for(unsigned char j = 0; j < fragmentLen; j++) cout << std::hex << (int)fragment[j];
			cout << endl;
			exit(1);
		}

		// shift the key and add the new value
		key = key << 2 | tValue;
	}
}

// consolidates hash hits into a read candidate (fast algorithm)
void CAlignmentThread::GetFastReadCandidate(HashRegion& region, char* query, const unsigned int queryLength, MhpOccupancyList* pMhpOccupancyList) {

	// localize the hash size
	unsigned char hashSize = mSettings.HashSize;

	// initialize the mhp occupancy list
	const unsigned int numHashes = queryLength - hashSize + 1;
	pMhpOccupancyList->resize(numHashes);
	MhpOccupancyList::iterator mhpIter = pMhpOccupancyList->begin();

	// get hash hits from the hash region tree
	AVLTree::CHashRegionTree hrt(queryLength, hashSize);
	uint64_t key;
	char* pQuery = query;

	for(unsigned int i = 0; i < numHashes; ++i, ++pQuery, ++mhpIter) {
		CreateHash(pQuery, hashSize, key);
		mhpIter->Begin = i;
		mhpIter->End   = i + hashSize - 1;
		mpDNAHash->Get(key, i, hrt, mhpIter->Occupancy);
	}

	// find the largest region
	unsigned int regionLength, largestRegionLength = 0;
	hrt.GotoFirstEntry();
	while(HashRegion* r = hrt.GetTraversalHashRegion()) {
		if(!hrt.GetNextEntry()) break;

		regionLength = r->End - r->Begin + 1;
		if(regionLength > largestRegionLength) {
			largestRegionLength = regionLength;
			region = *r;
		}
	}
}

// consolidates hash hits into read candidates
void CAlignmentThread::GetReadCandidates(vector<HashRegion>& regions, char* query, const unsigned int queryLength, MhpOccupancyList* pMhpOccupancyList) {

	// localize the hash size
	unsigned char hashSize = mSettings.HashSize;

	// initialize the mhp occupancy list
	const unsigned int numHashes = queryLength - hashSize + 1;
	pMhpOccupancyList->resize(numHashes);
	MhpOccupancyList::iterator mhpIter = pMhpOccupancyList->begin();

	// get hash hits from the hash region tree
	AVLTree::CHashRegionTree hrt(queryLength, hashSize);
	uint64_t key;
	char* pQuery = query;

	for(unsigned int i = 0; i < numHashes; ++i, ++pQuery, ++mhpIter) {
		CreateHash(pQuery, hashSize, key);
		mhpIter->Begin = i;
		mhpIter->End   = i + hashSize - 1;
		mpDNAHash->Get(key, i, hrt, mhpIter->Occupancy);
	}

	// add the consolidated regions
	regions.resize(hrt.GetCount());
	vector<HashRegion>::iterator hrIter = regions.begin();

	hrt.GotoFirstEntry();
	while(HashRegion* r = hrt.GetTraversalHashRegion()) {
		if(!hrt.GetNextEntry()) break;
		*hrIter = *r;
		hrIter++;
	}

	// sort the hash regions according to length (descending)
	if(!mFlags.IsAligningAllReads) sort(regions.begin(), regions.end(), SortHashRegionByLength());
}

// attempts to rescue the mate paired with a unique mate
bool CAlignmentThread::RescueMate(const LocalAlignmentModel& lam, const CMosaikString& bases, const unsigned int uniqueBegin, const unsigned int uniqueEnd, const unsigned int refIndex, Alignment& al) {

	// calculate the target regions using the local alignment models
	const unsigned int refBegin = mReferenceBegin[refIndex];
	const unsigned int refEnd   = mReferenceEnd[refIndex];
	unsigned int begin = uniqueBegin;
	unsigned int end   = uniqueEnd;

	if(lam.IsTargetBeforeUniqueMate) {

		if(begin >= mSettings.MedianFragmentLength)       begin -= mSettings.MedianFragmentLength;
		else begin = 0;

		if(begin >= mSettings.LocalAlignmentSearchRadius) begin -= mSettings.LocalAlignmentSearchRadius;
		else begin = 0;

		if(end   >= mSettings.MedianFragmentLength)       end   -= mSettings.MedianFragmentLength;
		else end = 0;

		end += mSettings.LocalAlignmentSearchRadius;

	} else {

		begin += mSettings.MedianFragmentLength;

		if(begin >= mSettings.LocalAlignmentSearchRadius) begin -= mSettings.LocalAlignmentSearchRadius;
		else begin = 0;

		end  += mSettings.MedianFragmentLength + mSettings.LocalAlignmentSearchRadius;
	}

	// make sure the endpoints are within the reference sequence
	if(begin < refBegin) begin = refBegin;
	if(end   < refBegin) end   = refBegin;
	if(begin > refEnd)   begin = refEnd;
	if(end   > refEnd)   end   = refEnd;

	// adjust the start position if the reference starts with a J nucleotide
	while(mReference[begin] == 'X') begin++;

	// adjust the stop position if the reference ends with a J nucleotide
	while(mReference[end] == 'X')   end--;

	// quit if we don't have a region to align against
	if(begin == end) return false;

	// prepare for alignment (borrow the forward read buffer)
	const char* query              = bases.CData();
	const unsigned int queryLength = bases.Length();

	strncpy_s(mForwardRead, mSettings.AllocatedReadLength, query, queryLength);
	mForwardRead[queryLength] = 0;

	if(lam.IsTargetReverseStrand) {
		if(mFlags.EnableColorspace) CSequenceUtilities::ReverseSequence(mForwardRead, queryLength);
		else CSequenceUtilities::GetReverseComplement(mForwardRead, queryLength);
	}

	// align according to the model
	al.IsReverseStrand = (lam.IsTargetReverseStrand ? true : false);
	char* pAnchor = mReference + begin;
	mSW.Align(al, pAnchor, (end - begin + 1), mForwardRead, queryLength);

	// adjust the reference start positions
	al.ReferenceIndex = refIndex;
	al.ReferenceBegin += begin - refBegin;
	al.ReferenceEnd   += begin - refBegin;

	// an alignment was performed
	return true;
}
