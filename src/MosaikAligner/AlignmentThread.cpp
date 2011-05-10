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
//pthread_mutex_t CAlignmentThread::mReportUnalignedMate1Mutex;
//pthread_mutex_t CAlignmentThread::mReportUnalignedMate2Mutex;
pthread_mutex_t CAlignmentThread::mSaveReadMutex;
pthread_mutex_t CAlignmentThread::mStatisticsMutex;
pthread_mutex_t CAlignmentThread::mStatisticsMapsMutex;
pthread_mutex_t CAlignmentThread::mSaveMultipleBamMutex;
pthread_mutex_t CAlignmentThread::mSaveSpecialBamMutex;
pthread_mutex_t CAlignmentThread::mSaveUnmappedBamMutex;

// define our constants
const double CAlignmentThread::P_ERR_REF  = pow(10.0, -REFERENCE_SEQUENCE_QUALITY / 10.0);
const double CAlignmentThread::P_CORR_REF = 1.0 - CAlignmentThread::P_ERR_REF;
const double CAlignmentThread::ONE_THIRD  = 1.0 / 3.0;
const double CAlignmentThread::TWO_NINTHS = 2.0 / 9.0;
//const double CAlignmentThread::UU_COEFFICIENT = 0.486796848;
//const double CAlignmentThread::UU_INTERCEPT   = 35.45967112;
//const double CAlignmentThread::UM_COEFFICIENT = 0.426395518;
//const double CAlignmentThread::UM_INTERCEPT   = 19.29236958;
//const double CAlignmentThread::MM_COEFFICIENT = 0.327358673;
//const double CAlignmentThread::MM_INTERCEPT   = 4.350331532;

// constructor
CAlignmentThread::CAlignmentThread(
	const AlignerAlgorithmType& algorithmType, 
	const FilterSettings&       filters, 
	const FlagData&             flags, 
	const AlignerModeType&      algorithmMode, 
	char*                       pAnchor, 
	const unsigned int          referenceLen, 
	CAbstractDnaHash*           pDnaHash, 
	const AlignerSettings&      settings, 
	unsigned int*               pRefBegin, 
	unsigned int*               pRefEnd, 
	char**                      pRefSpecies, 
	bool*                       pRefSpecial, 
	char**                      pBsRefSeqs, 
	const SReference&           SpecialReference, 
	map<unsigned int, MosaikReadFormat::ReadGroup>* pReadGroupsMap,
	const unsigned int          referenceOffset)
	: mAlgorithm(algorithmType)
	, mMode(algorithmMode)
	, mSettings(settings)
	, mFilters(filters)
	, mFlags(flags)
	, mReference(pAnchor)
	, mSReference(SpecialReference)
	, mForwardRead(NULL)
	, mReverseRead(NULL)
	, mReferenceLength(referenceLen)
	, mpDNAHash(pDnaHash)
	, mSW(CPairwiseUtilities::MatchScore, CPairwiseUtilities::MismatchScore, CPairwiseUtilities::GapOpenPenalty, CPairwiseUtilities::GapExtendPenalty)
	, mBSW(CPairwiseUtilities::MatchScore, CPairwiseUtilities::MismatchScore, CPairwiseUtilities::GapOpenPenalty, CPairwiseUtilities::GapExtendPenalty, settings.Bandwidth)
	, mReferenceBegin(pRefBegin)
	, mReferenceEnd(pRefEnd)
	, mReferenceSpecies(pRefSpecies)
	, mReferenceSpecial(pRefSpecial)
//	, softClippedIdentifierLength(2048)
	, mReadGroupsMap(pReadGroupsMap)
	, mReferenceOffset( referenceOffset )
{
	// calculate our base quality LUT
	for(unsigned char i = 0; i < 100; i++) mBaseQualityLUT[i] = pow(10.0, -i / 10.0);

	// set our flags
	if(algorithmMode == AlignerMode_ALL) mFlags.IsAligningAllReads = true;

	// assign the reference sequences to the colorspace utilities object
	mCS.SetReferenceSequences(pBsRefSeqs);

	// initialize our soft clip identifier
//	softClippedIdentifier = new char [ softClippedIdentifierLength + 1 ];
//	memset(softClippedIdentifier, 'Z', softClippedIdentifierLength);
//	softClippedIdentifier[ softClippedIdentifierLength ] = 0;

}

// destructor
CAlignmentThread::~CAlignmentThread(void) {
	if ( mForwardRead ) {
		delete [] mForwardRead; 
		mForwardRead = NULL;
	}

	if ( mReverseRead ) {
		delete [] mReverseRead; 
		mReverseRead = NULL;
	}

//	if ( softClippedIdentifier ) {
//		delete [] softClippedIdentifier; 
//		softClippedIdentifier = NULL;
//	}
}

// activates the current alignment thread
void* CAlignmentThread::StartThread(void* arg) {
	ThreadData* pTD = (ThreadData*)arg;

	// align reads
	CAlignmentThread at(
		pTD->Algorithm, 
		pTD->Filters, 
		pTD->Flags, 
		pTD->Mode, 
		pTD->pReference, 
		pTD->ReferenceLen, 
		pTD->pDnaHash, 
		pTD->Settings, 
		pTD->pRefBegin, 
		pTD->pRefEnd, 
		pTD->pRefSpecies, 
		pTD->pRefSpecial, 
		pTD->pBsRefSeqs, 
		pTD->SpecialReference, 
		pTD->pReadGroups,
		pTD->ReferenceOffset);

	at.AlignReadArchive(
		pTD->pIn, 
		pTD->pOut,
		//pTD->pUnalignedStream, 
		pTD->pReadCounter, 
		pTD->IsPairedEnd, 
		pTD->pMaps, 
		pTD->pBams,
		pTD->pCounters->StatMappingQuality);

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
void CAlignmentThread::AlignReadArchive(
	MosaikReadFormat::CReadReader* pIn, 
	MosaikReadFormat::CAlignmentWriter* pOut, 
	//FILE*     pUnalignedStream, 
	uint64_t* pReadCounter, 
	bool      isPairedEnd, 
	CStatisticsMaps* pMaps, 
	BamWriters*      pBams,
	unsigned char statMappingQuality) {

	// create our local alignment models
	const bool isUsing454          = (mSettings.SequencingTechnology == ST_454      ? true : false);
	const bool isUsingIllumina     = (mSettings.SequencingTechnology == ST_ILLUMINA ? true : false);
	const bool isUsingSOLiD        = (mSettings.SequencingTechnology == ST_SOLID    ? true : false);
	const bool isUsingIlluminaLong = (mSettings.SequencingTechnology == ST_ILLUMINA_LONG    ? true : false);

	// catch unsupported local alignment search sequencing technologies
	if(mFlags.UseLocalAlignmentSearch && (!isUsing454 && !isUsingIllumina && !isUsingIlluminaLong && !isUsingSOLiD )) {
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
	CNaiveAlignmentSet mate1Alignments(mReferenceLength, (isUsingIllumina || isUsingSOLiD || isUsingIlluminaLong)), mate2Alignments(mReferenceLength, (isUsingIllumina || isUsingSOLiD || isUsingIlluminaLong));

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

		// specify the percent of N's in the read
		unsigned short numMate1NBases = 0;
		//unsigned short numMate1ABases = 0;
		//unsigned short numMate1TBases = 0;
		//unsigned short numMate1CBases = 0;
		//unsigned short numMate1GBases = 0;
		char* basePtr = mr.Mate1.Bases.Data();
		for ( unsigned short i = 0; i < numMate1Bases; ++i ) {
			if ( *basePtr == 'N' )
				numMate1NBases++;
		//	if ( *basePtr == 'A' )
		//		numMate1ABases++;
		//	if ( *basePtr == 'T' )
		//		numMate1TBases++;
		//	if ( *basePtr == 'C' )
		//		numMate1CBases++;
		//	if ( *basePtr == 'G' )
		//		numMate1GBases++;
			++basePtr;
		}

		unsigned short numMate2NBases = 0;
		//unsigned short numMate2ABases = 0;
		//unsigned short numMate2TBases = 0;
		//unsigned short numMate2CBases = 0;
		//unsigned short numMate2GBases = 0;
		basePtr = mr.Mate2.Bases.Data();
		for ( unsigned short i = 0; i < numMate2Bases; ++i ) {
			if ( *basePtr == 'N' )
				numMate2NBases++;
		//	if ( *basePtr == 'A' )
		//		numMate2ABases++;
		//	if ( *basePtr == 'T' )
		//		numMate2TBases++;
		//	if ( *basePtr == 'C' )
		//		numMate2CBases++;
		//	if ( *basePtr == 'G' )
		//		numMate2GBases++;
			++basePtr;
		}

		const bool isTooManyNMate1 = ( ( numMate1NBases / (double) numMate1Bases ) > 0.7 ) 
		//	|| ( ( numMate1ABases / (double) numMate1Bases ) > 0.6 ) )
			//|| ( ( numMate1TBases / (double) numMate1Bases ) > 0.5 ) )
			//|| ( ( numMate1CBases / (double) numMate1Bases ) > 0.5 )
			//|| ( ( numMate1GBases / (double) numMate1Bases ) > 0.5 ) ) 
			? true : false;
		const bool isTooManyNMate2 = ( ( numMate2NBases / (double) numMate2Bases ) > 0.7 ) 
		//	|| ( ( numMate2ABases / (double) numMate1Bases ) > 0.6 ) )
			//|| ( ( numMate2TBases / (double) numMate1Bases ) > 0.5 ) )
			//|| ( ( numMate2CBases / (double) numMate1Bases ) > 0.5 )
			//|| ( ( numMate2GBases / (double) numMate1Bases ) > 0.5 ) ) 
			? true : false;

		//const bool isTooManyNMate1 = false;
		//const bool isTooManyNMate2 = false;

		const bool areBothMatesPresent = (((numMate1Bases != 0) && (numMate2Bases != 0)) ? true : false);

		// ====================
		// align the first mate
		// ====================

		bool isMate1Aligned = false;
		mate1Alignments.Clear();
		if(numMate1Bases != 0 && !isTooManyNMate1 ) {

			// align the read
			if(AlignRead(mate1Alignments, mr.Mate1.Bases.CData(), mr.Mate1.Qualities.CData(), numMate1Bases, mate1Status)) {
				// calculate the alignment qualities
				mate1Alignments.CalculateAlignmentQualities(calculateCorrectionCoefficient, minSpanLength);
				isMate1Aligned = true;

			} 
		}

		// =====================
		// align the second mate
		// =====================

		bool isMate2Aligned = false;
		mate2Alignments.Clear();
		if(numMate2Bases != 0 && !isTooManyNMate2 ) {

			// align the read
			if(AlignRead(mate2Alignments, mr.Mate2.Bases.CData(), mr.Mate2.Qualities.CData(), numMate2Bases, mate2Status)) {

				// calculate the alignment qualities
				mate2Alignments.CalculateAlignmentQualities(calculateCorrectionCoefficient, minSpanLength);
				isMate2Aligned = true;

			} 
		}

		// ======================
		// local alignment search
		// ======================

		bool isMate1Unique = mate1Alignments.IsUnique();
		bool isMate2Unique = mate2Alignments.IsUnique();

		// we can only perform a local alignment search if both mates are present
		if(areBothMatesPresent) {
			// perform local alignment search WHEN MATE1 IS UNIQUE
			if(mFlags.UseLocalAlignmentSearch && isMate1Unique && !isTooManyNMate2 ) {

				// extract the unique begin and end coordinates
				AlignmentSet::const_iterator uniqueIter = mate1Alignments.GetSet()->begin();
				const unsigned int refIndex    = uniqueIter->ReferenceIndex;
				const unsigned int uniqueBegin = mReferenceBegin[refIndex] + uniqueIter->ReferenceBegin;
				const unsigned int uniqueEnd   = mReferenceBegin[refIndex] + uniqueIter->ReferenceEnd;

				// create the appropriate local alignment search model
				LocalAlignmentModel lam;
				if(uniqueIter->IsReverseStrand) {
					if( isUsingIllumina ) 
						lam.IsTargetBeforeUniqueMate = true;
					else if( isUsing454 ) 
						lam.IsTargetReverseStrand    = true;
					else if ( isUsingSOLiD ) {
						lam.IsTargetBeforeUniqueMate = true;
						lam.IsTargetReverseStrand    = true;
					}
						
					// NB: we caught other technologies during initialization
				} else {
					if( isUsingIllumina ) 
						lam.IsTargetReverseStrand    = true;
					else if( isUsing454 ) 
						lam.IsTargetBeforeUniqueMate = true;
					else if ( isUsingIlluminaLong ) {
						lam.IsTargetReverseStrand    = true;
						lam.IsTargetBeforeUniqueMate = true;
					}
					// NB: we caught other technologies during initialization
				}

				bool settleLocalSearchWindow = false, isAlExisting = false;
				unsigned int localSearchBegin, localSearchEnd;
				settleLocalSearchWindow = SettleLocalSearchRegion( lam, refIndex, uniqueBegin, uniqueEnd, localSearchBegin, localSearchEnd );
				
				if ( settleLocalSearchWindow )
					// check if there are alignments already sitting in the region
					isAlExisting = mate2Alignments.CheckExistence( refIndex, localSearchBegin - mReferenceBegin[refIndex], localSearchEnd - mReferenceBegin[refIndex] );
				
				if ( !isAlExisting ) {

					Alignment al;
					
					if( RescueMate( lam, mr.Mate2.Bases, localSearchBegin, localSearchEnd, refIndex, al ) ) {

						const char* pQualities = mr.Mate2.Qualities.CData();
						const char* pBases = mr.Mate2.Bases.CData();

						// add the alignment to the alignment set if it passes the filters
						if( ApplyReadFilters( al, pBases, pQualities, mr.Mate2.Bases.Length() ) ) {
							al.WasRescued = true;
							if( mate2Alignments.Add( al ) ) {
								mStatisticsCounters.AdditionalLocalMates++;
								mate2Status    = ALIGNMENTSTATUS_GOOD;
								isMate2Aligned = true;
							}
						}

						// increment our candidates counter
						mStatisticsCounters.AlignmentCandidates++;
					}
					

				}
			}

			// perform local alignment search WHEN MATE2 IS UNIQUE
			if(mFlags.UseLocalAlignmentSearch && isMate2Unique && !isTooManyNMate1 ) {

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

				bool settleLocalSearchWindow = false, isAlExisting = false;
				unsigned int localSearchBegin, localSearchEnd;
				settleLocalSearchWindow = SettleLocalSearchRegion( lam, refIndex, uniqueBegin, uniqueEnd, localSearchBegin, localSearchEnd );

				if ( settleLocalSearchWindow )
					// check if there are alignments already sitting in the region
					isAlExisting = mate1Alignments.CheckExistence( refIndex, localSearchBegin - mReferenceBegin[refIndex], localSearchEnd - mReferenceBegin[refIndex]);

				if ( !isAlExisting ) {

					Alignment al;
					
					if(RescueMate(lam, mr.Mate1.Bases, localSearchBegin, localSearchEnd, refIndex, al)) {

						const char* pQualities = mr.Mate1.Qualities.CData();
						const char* pBases = mr.Mate1.Bases.Data();

						// add the alignment to the alignment set if it passes the filters
						if( ApplyReadFilters( al, pBases, pQualities, mr.Mate1.Bases.Length() ) ) {
							al.WasRescued = true;
							if( mate1Alignments.Add( al ) ) {
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


		// process alignments mapped in special references and delete them in vectors
		vector<Alignment> mate1Set = *mate1Alignments.GetSet();
		vector<Alignment> mate2Set = *mate2Alignments.GetSet();
		bool isLongRead = mate1Alignments.HasLongAlignment() || mate2Alignments.HasLongAlignment();
		mate1Alignments.Clear();
		mate2Alignments.Clear();


		Alignment mate1SpecialAl, mate2SpecialAl;
		bool isMate1Special = false, isMate2Special = false;

		// we don't remove special alignment here,
		// the special alignments will be considered when merging archives
		if ( mSReference.found && !mFlags.UseArchiveOutput )
			ProcessSpecialAlignment( mate1Set, mate2Set, mate1SpecialAl, mate2SpecialAl, isMate1Special, isMate2Special );
		
		// deleting special alignments may let mate1Alignments or mate2Alignments become empty
		// so we have to check them again
		isMate1Aligned = !mate1Set.empty();
		isMate2Aligned = !mate2Set.empty();

		isMate1Unique = ( mate1Set.size() == 1 ) ? true: false;
		isMate2Unique = ( mate2Set.size() == 1 ) ? true: false;

		bool isMate1Multiple = ( mate1Set.size() > 1 ) ? true: false;
		bool isMate2Multiple = ( mate2Set.size() > 1 ) ? true: false;

		bool isMate1Empty = mate1Set.empty();
		bool isMate2Empty = mate2Set.empty();
		
		// update the paired-end read statistics
		if(isPairedEnd) {
			if(isMate1Aligned || isMate2Aligned) {
				if(isMate1Aligned && isMate2Aligned) {
					if(isMate1Unique && isMate2Unique) mStatisticsCounters.BothUniqueReads++;
					else if(!isMate1Unique && !isMate2Unique) mStatisticsCounters.BothNonUniqueReads++;
					else mStatisticsCounters.OneNonUniqueReads++;
				} else mStatisticsCounters.OrphanedReads++;
			} else mStatisticsCounters.UnalignedReads++;
		}

		// update the mate bases aligned
		if(isMate1Aligned) mStatisticsCounters.MateBasesAligned += numMate1Bases;
		if(isMate2Aligned) mStatisticsCounters.MateBasesAligned += numMate2Bases;


		// save chromosomes and positions of multiple alignments in bam
		if ( mFlags.SaveMultiplyBam ) {
			// -om is enabled
			if ( mFlags.OutputMultiply ) {
				if ( isMate1Multiple ) {
					Alignment mateAl;
					if ( !isMate2Empty ) mateAl = mate2Set[0];

					vector<Alignment> mate1SetTemp = mate1Set;
					for(vector<Alignment>::iterator alIter = mate1SetTemp.begin(); alIter != mate1SetTemp.end(); ++alIter) {
						if ( !isMate2Empty )
							alIter->SetPairFlagsAndFragmentLength(mateAl, 0, 0, mSettings.SequencingTechnology);
						alIter->RecalibratedQuality = alIter->Quality;
						SetRequiredInfo( *alIter, mateAl, mr.Mate1, mr, !isMate2Empty, false, true, isPairedEnd, true, !isMate2Empty);
					}
					
					pthread_mutex_lock(&mSaveMultipleBamMutex);
					for(vector<Alignment>::iterator alIter = mate1SetTemp.begin(); alIter != mate1SetTemp.end(); ++alIter)
						pBams->mBam.SaveAlignment( *alIter, 0, false, false, mFlags.EnableColorspace );
					pthread_mutex_unlock(&mSaveMultipleBamMutex);
				}
				if ( isPairedEnd && isMate2Multiple ) {
					Alignment mateAl;
					if ( !isMate1Empty ) mateAl = mate1Set[0];

					vector<Alignment> mate2SetTemp = mate2Set;
					for(vector<Alignment>::iterator alIter = mate2SetTemp.begin(); alIter != mate2SetTemp.end(); ++alIter) {
						if ( !isMate1Empty )
							alIter->SetPairFlagsAndFragmentLength(mateAl, 0, 0, mSettings.SequencingTechnology);
						alIter->RecalibratedQuality = alIter->Quality;
						SetRequiredInfo( *alIter, mateAl, mr.Mate2, mr, !isMate1Empty, false, false, isPairedEnd, true, !isMate1Empty);
					}
					
					pthread_mutex_lock(&mSaveMultipleBamMutex);
					for(vector<Alignment>::iterator alIter = mate2SetTemp.begin(); alIter != mate2SetTemp.end(); ++alIter)
						pBams->mBam.SaveAlignment( *alIter, 0, false, false, mFlags.EnableColorspace );
					pthread_mutex_unlock(&mSaveMultipleBamMutex);
				}
			}
			else {
				if ( isMate1Multiple ) {
					pthread_mutex_lock(&mSaveMultipleBamMutex);
					for(vector<Alignment>::iterator alIter = mate1Set.begin(); alIter != mate1Set.end(); ++alIter) {
						pBams->mBam.SaveReferencePosition( alIter->ReferenceIndex + mReferenceOffset, alIter->ReferenceBegin, alIter->ReferenceEnd );
					}
					pthread_mutex_unlock(&mSaveMultipleBamMutex);
				}
				if ( isPairedEnd && isMate2Multiple ) {
					pthread_mutex_lock(&mSaveMultipleBamMutex);
					for(vector<Alignment>::iterator alIter = mate2Set.begin(); alIter != mate2Set.end(); ++alIter)
						pBams->mBam.SaveReferencePosition( alIter->ReferenceIndex + mReferenceOffset, alIter->ReferenceBegin, alIter->ReferenceEnd );
					pthread_mutex_unlock(&mSaveMultipleBamMutex);
				}
			}
		}
		

		// ===================================
		// Save alignments to BAMs or archives
		// ===================================
		
		// UU, UM, and MM pair
		if ( ( isMate1Unique && isMate2Unique )
			|| ( isMate1Unique && isMate2Multiple )
			|| ( isMate1Multiple && isMate2Unique )
			|| ( isMate1Multiple && isMate2Multiple ) ) {
		
			if ( ( isMate1Unique && isMate2Multiple )
				|| ( isMate1Multiple && isMate2Unique )
				|| ( isMate1Multiple && isMate2Multiple ) )
				BestNSecondBestSelection::Select( mate1Set, mate2Set, mSettings.MedianFragmentLength, mSettings.SequencingTechnology );

			isMate1Empty = mate1Set.empty();
			isMate2Empty = mate2Set.empty();
			
			// sanity check
			if ( isMate1Empty | isMate2Empty ) {
				cout << "ERROR: One of mate sets is empty after apllying best and second best selection." << endl;
				exit(1);
			}

			// patch the information for reporting
			Alignment al1 = mate1Set[0], al2 = mate2Set[0];
			
			// TODO: handle fragment length for others sequencing techs
			int minFl = mSettings.MedianFragmentLength - mSettings.LocalAlignmentSearchRadius;
			int maxFl = mSettings.MedianFragmentLength + mSettings.LocalAlignmentSearchRadius;
			bool properPair1 = false, properPair2 = false;
			al1.IsFirstMate = true;
			al2.IsFirstMate = false;
			// MM pair is always an improper pair
			if ( !isMate1Multiple || !isMate2Multiple ) {
				properPair1 = al1.SetPairFlagsAndFragmentLength( al2, minFl, maxFl, mSettings.SequencingTechnology );
				properPair2 = al2.SetPairFlagsAndFragmentLength( al1, minFl, maxFl, mSettings.SequencingTechnology );
			}

			// sanity checker
			if ( properPair1 != properPair2 ) {
				cout << "ERROR: An inconsistent proper pair is found." << endl;
				exit(1);
			}

			const bool isUU = isMate1Unique && isMate2Unique;
			const bool isMM = isMate1Multiple && isMate2Multiple;


			SetRequiredInfo( al1, al2, mr.Mate1, mr, true, properPair1, true, isPairedEnd, true, true );
			SetRequiredInfo( al2, al1, mr.Mate2, mr, true, properPair2, false, isPairedEnd, true, true );
			

			if ( mFlags.UseArchiveOutput ) {
				//bool isLongRead = mate1Alignments.HasLongAlignment() || mate2Alignments.HasLongAlignment();
				isLongRead |= ( ( al1.CsQuery.size() > 255 ) || ( al2.CsQuery.size() > 255 ) );
				pthread_mutex_lock(&mSaveReadMutex);
				pOut->SaveRead( mr, al1, al2, isLongRead, true, true, mFlags.SaveUnmappedBasesInArchive );
				pthread_mutex_unlock(&mSaveReadMutex);
			
			} else {

				//Note: the following function will set RecalibratedQuality and RecalibratedQuality will be shown in the bam
				al1.RecalibrateQuality(isUU, isMM);
				al2.RecalibrateQuality(isUU, isMM);
				
				if (  isMate2Special  ) {
					Alignment genomicAl = al1;
					Alignment specialAl = mate2SpecialAl;

					SetRequiredInfo( specialAl, genomicAl, mr.Mate2, mr, true, false, false, isPairedEnd, true, true );
				
					//CZaTager zas1, zas2;

					const char *zas1Tag = za1.GetZaTag( genomicAl, al2, true );
					const char *zas2Tag = za2.GetZaTag( specialAl, genomicAl, false );
					pthread_mutex_lock(&mSaveSpecialBamMutex);
					pBams->sBam.SaveAlignment( genomicAl, zas1Tag, false, false, mFlags.EnableColorspace );
					pBams->sBam.SaveAlignment( specialAl, zas2Tag, false, false, mFlags.EnableColorspace );
					pthread_mutex_unlock(&mSaveSpecialBamMutex);
				}

				if (  isMate1Special  ) {
					Alignment genomicAl = al2;
					Alignment specialAl = mate1SpecialAl;
					SetRequiredInfo( specialAl, genomicAl, mr.Mate1, mr, true, false, true, isPairedEnd, true, true );
	
					//CZaTager zas1, zas2;

					const char *zas1Tag = za1.GetZaTag( genomicAl, al1, false );
					const char *zas2Tag = za2.GetZaTag( specialAl, genomicAl, true );
					pthread_mutex_lock(&mSaveSpecialBamMutex);
					pBams->sBam.SaveAlignment( genomicAl, zas1Tag, false, false, mFlags.EnableColorspace );
					pBams->sBam.SaveAlignment( specialAl, zas2Tag, false, false, mFlags.EnableColorspace );
					pthread_mutex_unlock(&mSaveSpecialBamMutex);
				}


				//CZaTager za1, za2;
				const char* zaTag1 = za1.GetZaTag( al1, al2, true );
				const char* zaTag2 = za2.GetZaTag( al2, al1, false );
				//if ( isMate1Multiple ) al1.Quality = 0;
				//if ( isMate2Multiple ) al2.Quality = 0;
				pthread_mutex_lock(&mSaveReadMutex);
				pBams->rBam.SaveAlignment( al1, zaTag1, false, false, mFlags.EnableColorspace );
				pBams->rBam.SaveAlignment( al2, zaTag2, false, false, mFlags.EnableColorspace );
				pthread_mutex_unlock(&mSaveReadMutex);

				if ( ( statMappingQuality <= al1.Quality ) && ( statMappingQuality <= al2.Quality ) ) {
					pthread_mutex_lock(&mStatisticsMapsMutex);
					pMaps->SaveRecord( al1, al2, isPairedEnd, mSettings.SequencingTechnology );
					pthread_mutex_unlock(&mStatisticsMapsMutex);
				}

			}

			mStatisticsCounters.AlignedReads++;

		// UX and MX pair
		} else if ( ( isMate1Empty || isMate2Empty )
			&&  !( isMate1Empty && isMate2Empty ) ) {

			if ( isMate1Multiple || isMate2Multiple ) 
				BestNSecondBestSelection::Select( mate1Set, mate2Set, mSettings.MedianFragmentLength, mSettings.SequencingTechnology );

			isMate1Empty = mate1Set.empty();
			isMate2Empty = mate2Set.empty();
			
			bool isFirstMate;
			if ( !isMate1Empty ) {
				isFirstMate = true;
			} else if ( !isMate2Empty ) {
				isFirstMate = false;
				if ( !isPairedEnd ) {
					cout << "ERROR: The sequence technology is single-end, but second mate is aligned." << endl;
					exit(1);
				}
			} else {
				cout << "ERROR: Both mates are empty after applying best and second best selection." << endl;
				exit(1);
			}
		
			// patch the information for reporting
			Alignment al = isFirstMate ? mate1Set[0] : mate2Set[0];
			Alignment unmappedAl;

			SetRequiredInfo( al, unmappedAl, ( isFirstMate ? mr.Mate1 : mr.Mate2 ), mr, false, false, isFirstMate, isPairedEnd, true, false );
			SetRequiredInfo( unmappedAl, al, ( isFirstMate ? mr.Mate2 : mr.Mate1 ), mr, true, false, !isFirstMate, isPairedEnd, false, true );


			if ( mFlags.UseArchiveOutput ) {
				//bool isLongRead = mate1Alignments.HasLongAlignment() || mate2Alignments.HasLongAlignment();
				isLongRead |= ( ( al.CsQuery.size() > 255 ) || ( unmappedAl.CsQuery.size() > 255 ) );
				pthread_mutex_lock(&mSaveReadMutex);
				pOut->SaveRead( mr, ( isFirstMate ? al : unmappedAl ), ( isFirstMate ? unmappedAl : al ), isLongRead, true, isPairedEnd, mFlags.SaveUnmappedBasesInArchive );
				pthread_mutex_unlock(&mSaveReadMutex);

			} else {
			
				al.RecalibratedQuality = al.Quality;
				
				// show the original MQs in ZAs, and zeros in MQs fields of a BAM
				const char* zaTag1 = za1.GetZaTag( al, unmappedAl, isFirstMate, !isPairedEnd, true );
				const char* zaTag2 = za2.GetZaTag( unmappedAl, al, !isFirstMate, !isPairedEnd, false );

				// Note: RecalibratedQuality will be shown in the bam
				if ( isFirstMate && isMate1Multiple )
					al.RecalibratedQuality = 0;
				else if ( !isFirstMate && isMate2Multiple )
					al.RecalibratedQuality = 0;


				if ( isPairedEnd ) {
					
					if ( isMate1Special ) {
						
						Alignment specialAl = mate1SpecialAl;
						SetRequiredInfo( specialAl, al, mr.Mate1, mr, !isFirstMate, false, true, isPairedEnd, true, !isFirstMate );
	
						//CZaTager zas1, zas2;
						const char *zas2Tag = isFirstMate ? za2.GetZaTag( specialAl, al, true, !isPairedEnd, true ) : za2.GetZaTag( specialAl, al, true );

						pthread_mutex_lock(&mSaveSpecialBamMutex);
						pBams->sBam.SaveAlignment( specialAl, zas2Tag, false, false, mFlags.EnableColorspace );
						pthread_mutex_unlock(&mSaveSpecialBamMutex);
					}
					
					if ( isMate2Special ) {
						
						Alignment specialAl = mate2SpecialAl ;
						SetRequiredInfo( specialAl, al, mr.Mate2, mr, isFirstMate, false, false, isPairedEnd, true, isFirstMate );
	
						//CZaTager zas1, zas2;

						const char *zas2Tag = !isFirstMate ? za2.GetZaTag( specialAl, al, false, !isPairedEnd, true ) : za2.GetZaTag( specialAl, al, false );

						pthread_mutex_lock(&mSaveSpecialBamMutex);
						pBams->sBam.SaveAlignment( specialAl, zas2Tag, false, false, mFlags.EnableColorspace );
						pthread_mutex_unlock(&mSaveSpecialBamMutex);
					}
					
					unmappedAl.ReferenceBegin = al.ReferenceBegin;
					unmappedAl.ReferenceIndex = al.ReferenceIndex;

					pthread_mutex_lock(&mSaveReadMutex);
					pBams->rBam.SaveAlignment( al, zaTag1, false, false, mFlags.EnableColorspace );
					pBams->rBam.SaveAlignment( unmappedAl, zaTag2, true, false, mFlags.EnableColorspace );  // noCigarMdNm
					pthread_mutex_unlock(&mSaveReadMutex);

					
					pthread_mutex_lock(&mSaveUnmappedBamMutex);
					pBams->uBam.SaveAlignment( unmappedAl, 0, true, false, mFlags.EnableColorspace );
					pthread_mutex_unlock(&mSaveUnmappedBamMutex);
				}
				// single end
				else {
					
					// store special hits
					if ( isMate1Special ) {
						Alignment specialAl = mate1SpecialAl;
						SetRequiredInfo( specialAl, al, mr.Mate1, mr, false, false, true, isPairedEnd, true, false );
						
						const char *zas2Tag = za2.GetZaTag( specialAl, al, true, !isPairedEnd, true );
						pthread_mutex_lock(&mSaveSpecialBamMutex);
						pBams->sBam.SaveAlignment( specialAl, zas2Tag, false, false, mFlags.EnableColorspace );
						pthread_mutex_unlock(&mSaveSpecialBamMutex);
					}
					
					pthread_mutex_lock(&mSaveReadMutex);
					pBams->rBam.SaveAlignment( al, zaTag1, false, false, mFlags.EnableColorspace );
					pthread_mutex_unlock(&mSaveReadMutex);
				}
				
				if ( statMappingQuality <= al.Quality ) {
					pthread_mutex_lock(&mStatisticsMapsMutex);
					pMaps->SaveRecord( ( isFirstMate ? al : unmappedAl ), ( isFirstMate ? unmappedAl : al), isPairedEnd, mSettings.SequencingTechnology );
					pthread_mutex_unlock(&mStatisticsMapsMutex);
				}
			}
			
			mStatisticsCounters.AlignedReads++;
		
		// XX
		} else if ( isMate1Empty && isMate2Empty ) {
			
			Alignment unmappedAl1, unmappedAl2;
			SetRequiredInfo( unmappedAl1, unmappedAl2, mr.Mate1, mr, false, false, true, isPairedEnd, false, false );
			SetRequiredInfo( unmappedAl2, unmappedAl1, mr.Mate2, mr, false, false, false, isPairedEnd, false, false );

			if ( mFlags.UseArchiveOutput ) {

				bool isLongReadXX = ( ( unmappedAl1.QueryEnd > 255 ) || ( unmappedAl2.QueryEnd > 255 ) ) ? true : false;
				isLongReadXX |= ( ( unmappedAl1.CsQuery.size() > 255 ) || ( unmappedAl2.CsQuery.size() > 255 ) );
			
				pthread_mutex_lock(&mSaveReadMutex);
				pOut->SaveRead( mr, unmappedAl1, unmappedAl2, isLongReadXX, true, isPairedEnd, mFlags.SaveUnmappedBasesInArchive );
				pthread_mutex_unlock(&mSaveReadMutex);

			} else {
				
				if ( isPairedEnd ) {

					// store special hits
					if ( isMate1Special ) {
						
						Alignment specialAl = mate1SpecialAl;
						SetRequiredInfo( specialAl, unmappedAl2, mr.Mate1, mr, false, false, true, isPairedEnd, true, false );
	
						//CZaTager zas1, zas2;
						const char *zas2Tag = za2.GetZaTag( specialAl, unmappedAl2, true, !isPairedEnd, true );

						pthread_mutex_lock(&mSaveSpecialBamMutex);
						pBams->sBam.SaveAlignment( specialAl, zas2Tag, false, false, mFlags.EnableColorspace );
						pthread_mutex_unlock(&mSaveSpecialBamMutex);
					}
					
					if ( isMate2Special ) {
						
						Alignment specialAl = mate2SpecialAl;
						SetRequiredInfo( specialAl, unmappedAl1, mr.Mate2, mr, false, false, false, isPairedEnd, true, false );
	
						//CZaTager zas1, zas2;

						const char *zas2Tag = za2.GetZaTag( specialAl, unmappedAl1, false, !isPairedEnd, true );

						pthread_mutex_lock(&mSaveSpecialBamMutex);
						pBams->sBam.SaveAlignment( specialAl, zas2Tag, false, false, mFlags.EnableColorspace );
						pthread_mutex_unlock(&mSaveSpecialBamMutex);
					}
					

					pthread_mutex_lock(&mSaveUnmappedBamMutex);
					pBams->uBam.SaveAlignment( unmappedAl1, 0, true, false, mFlags.EnableColorspace );
					pBams->uBam.SaveAlignment( unmappedAl2, 0, true, false, mFlags.EnableColorspace );
					pthread_mutex_unlock(&mSaveUnmappedBamMutex);
				} else {

					// store special hits
					if ( isMate1Special ) {
						Alignment specialAl = mate1SpecialAl;
						SetRequiredInfo( specialAl, unmappedAl2, mr.Mate1, mr, false, false, true, isPairedEnd, true, false );
						
						const char *zas2Tag = za2.GetZaTag( specialAl, unmappedAl2, true, !isPairedEnd, true );
						pthread_mutex_lock(&mSaveSpecialBamMutex);
						pBams->sBam.SaveAlignment( specialAl, zas2Tag, false, false, mFlags.EnableColorspace );
						pthread_mutex_unlock(&mSaveSpecialBamMutex);
					}
					
					pthread_mutex_lock(&mSaveUnmappedBamMutex);
					pBams->uBam.SaveAlignment( unmappedAl1, 0, true, false, mFlags.EnableColorspace );
					pthread_mutex_unlock(&mSaveUnmappedBamMutex);
				}

				//pthread_mutex_lock(&mStatisticsMapsMutex);
				//pMaps->SaveRecord( unmappedAl1, unmappedAl2, isPairedEnd, mSettings.SequencingTechnology );
				//pthread_mutex_unlock(&mStatisticsMapsMutex);

			}

		} else {
			cout << "ERROR: Unknown pairs." << endl;
			//cout << mate1Alignments.GetCount() << "\t" << mate2Alignments.GetCount() << endl;
			exit(1);
		}

	}
}

// adjust alignment quality
//inline int CAlignmentThread::AdjustMappingQuality ( const unsigned char& mq, const bool& isUU, const bool& isMM ) {
//	int aq = mq;
//
//	if ( isUU )      aq = (int) ( UU_COEFFICIENT * aq + UU_INTERCEPT );
//	else if ( isMM ) aq = (int) ( MM_COEFFICIENT * aq + MM_INTERCEPT );
//	else             aq = (int) ( UM_COEFFICIENT * aq + UM_INTERCEPT );
//
//	if(aq < 0)       aq = 0;
//	else if(aq > 99) aq = 99;
//	
//
//	return aq;
//}

// Set the required information and flag for alignments
void CAlignmentThread::SetRequiredInfo (
	Alignment& al,
	const Alignment& mate,
	const Mosaik::Mate& m,
	const Mosaik::Read& r,
	const bool& isPair,
	const bool& isProperPair,
	const bool& isFirstMate,
	const bool& isPairTech,
	const bool& isItselfMapped,
	const bool& isMateMapped) {

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
	rgIte = mReadGroupsMap->find( r.ReadGroupCode );
	// sanity check
	if ( rgIte == mReadGroupsMap->end() ) {
		cout << "ERROR: ReadGroup cannot be found." << endl;
		exit(1);
	}	
	else 
		al.ReadGroup = rgIte->second.ReadGroupID;

	// fill out the alignment
	// not SOLiD
	if ( !mFlags.EnableColorspace ) {
		// copy qualities
		al.BaseQualities    = m.Qualities;

		// handle bases
		CMosaikString patchBases   = m.Bases;
		unsigned int patchStartLen = al.QueryBegin;
		unsigned int patchEndLen   = patchBases.Length() - al.QueryEnd - 1;

		if ( al.IsReverseStrand ) {
			// reverse qualities
			al.BaseQualities.Reverse();
			// reverse complement bases
			patchBases.ReverseComplement();
			// re-calculate patching start and end
			unsigned int temp;
			temp = patchStartLen;
			patchStartLen = patchEndLen;
			patchEndLen = temp;
		}

		if ( !isItselfMapped ) {
			al.NumMapped = 0;
			if ( ( !mFlags.UseArchiveOutput ) || ( mFlags.UseArchiveOutput && mFlags.SaveUnmappedBasesInArchive ) ) {
				al.Query = m.Bases;
				al.Reference.Copy( 'Z', al.Query.Length() );
				al.IsJunk = false;
			} else {
				al.IsJunk = true;
			}
		}
		else {	
			// patch bases
			if ( patchStartLen > 0 ) {
				al.Query.Prepend    ( patchBases.CData(), patchStartLen );
				al.Reference.Prepend( 'Z', patchStartLen );
			}

			if ( patchEndLen > 0 ) {
				const unsigned int length = patchBases.Length();
				const unsigned int start  = length - patchEndLen;
				const char* startPoint    = patchBases.CData() + start;
				// sanity check
				if ( length > patchBases.Length() ) {
					cout << "ERROR: The soft chip position is wrong." << endl;
					exit(1);
				}
				al.Query.Append    ( startPoint, patchEndLen );
				al.Reference.Append( 'Z', patchEndLen );
			}
		}

		al.QueryBegin = 0;
		al.QueryEnd   = m.Bases.Length() - 1;
	}
	else {
		// fill out Colorspace raw bases and qualites
		al.CsQuery.clear();
		al.CsBaseQualities.clear();
		if ( ( !mFlags.UseArchiveOutput ) || ( mFlags.UseArchiveOutput && mFlags.SaveUnmappedBasesInArchive ) ) {
			// raw sequence
			CMosaikString rawCS = m.Bases;
			mCS.ConvertReadPseudoColorspaceToColorspace( rawCS );
			al.CsQuery.insert( 0, m.SolidPrefixTransition, SOLID_PREFIX_LENGTH );
			al.CsQuery += rawCS.CData();

			// raw base qualities
			// Note: if the first quality base is not '!'
			CMosaikString rawCQ = m.Qualities;
			rawCQ.Increment(33);
			char prefix = '!';
			rawCQ.Prepend( &prefix, 1);
			al.CsBaseQualities = rawCQ.CData();
		}

		if ( !isItselfMapped ) {
			al.NumMapped = 0;
			al.IsJunk = true;
		} else {
			// patch N's
			--al.QueryEnd;
			const unsigned int readLength    = m.Bases.Length();
			const unsigned int patchStartLen = al.IsReverseStrand ? readLength - al.QueryEnd - 1 : al.QueryBegin;
			const unsigned int patchEndLen   = al.IsReverseStrand ? al.QueryBegin : readLength - al.QueryEnd - 1;

			// sanity checker
			if ( al.QueryEnd > ( readLength - 1 ) ) {
				cout << "ERROR: The aligned length is larger than the read length." << endl;
				exit(1);
			}

			if ( patchStartLen == 0 ) {
				al.Query.TrimBegin(1);
				al.BaseQualities.TrimBegin(1);
				al.Reference.TrimBegin(1);
			}
			else if ( patchStartLen > 1 ) {
				al.Query.Prepend( 'N', patchStartLen - 1 );
				al.BaseQualities.Prepend( (char)0, patchStartLen - 1 );
				al.Reference.Prepend( 'Z', patchStartLen - 1 );
			}

			if ( patchEndLen > 0 ) {
				al.Query.Append( 'N', patchEndLen );
				al.BaseQualities.Append( (char)0, patchEndLen );
				al.Reference.Append( 'Z', patchEndLen );
			}

			al.QueryBegin  = 0;
			al.QueryEnd    = readLength - 1;
			al.QueryLength = al.QueryEnd - al.QueryBegin + 1;
		}
	}

}

// handle and then delete special alignments
// also assign two special characters for the general alignments
void CAlignmentThread::ProcessSpecialAlignment ( vector<Alignment>& mate1Set, vector<Alignment>& mate2Set, 
	Alignment& mate1SpecialAl, Alignment& mate2SpecialAl,
	bool& isMate1Special, bool& isMate2Special ) {

	unsigned int nMobAl = 0;
	string specialCode;
	specialCode.resize(3);
	for ( vector<Alignment>::iterator ite = mate1Set.begin(); ite != mate1Set.end(); ++ite ) {

		if ( ite->IsMappedSpecialReference ) {
			
			nMobAl++;
			char* tempCode = mReferenceSpecies[ ite->ReferenceIndex ];
			specialCode[0] = *tempCode;
			tempCode++;
			specialCode[1] = *tempCode;
			specialCode[2] = 0;

		}

	}

	if ( ( nMobAl == mate1Set.size() ) && ( mate1Set.size() != 0 ) ) {

		isMate1Special = true;
		sort ( mate1Set.begin(), mate1Set.end(), LessThanMQ );
		mate1SpecialAl = *mate1Set.rbegin();
		mate1SpecialAl.SpecialCode = specialCode;
		mate1SpecialAl.NumMapped   = nMobAl;
		mate1Set.clear();

	} else if ( nMobAl > 0 ) {
		vector<Alignment> newMate1Set;
		for ( vector<Alignment>::iterator ite = mate1Set.begin(); ite != mate1Set.end(); ++ite ) {

			if ( !ite->IsMappedSpecialReference ) {
				ite->CanBeMappedToSpecialReference = true;
				ite->SpecialCode = specialCode;
				newMate1Set.push_back( *ite );
			} else {
				isMate1Special = true;
				mate1SpecialAl = ( ite->Quality >= mate1SpecialAl.Quality ) ? *ite : mate1SpecialAl;
				mate1SpecialAl.SpecialCode = specialCode;
				mate1SpecialAl.NumMapped   = nMobAl;
			}
		
		}
		mate1Set.clear();
		mate1Set = newMate1Set;
	}

	unsigned int nMobAl2 = 0;
	string specialCode2;
	specialCode2.resize(3);

	for ( vector<Alignment>::iterator ite = mate2Set.begin(); ite != mate2Set.end(); ++ite ) {
		
		if ( ite->IsMappedSpecialReference ) {
			
			nMobAl2++;
			char* tempCode = mReferenceSpecies[ ite->ReferenceIndex ];
			specialCode2[0] = *tempCode;
			tempCode++;
			specialCode2[1] = *tempCode;
			specialCode2[2] = 0;

		}
	}

	if ( ( nMobAl2 == mate2Set.size() ) && ( mate2Set.size() != 0 ) ) {

		isMate2Special = true;
		sort ( mate2Set.begin(), mate2Set.end(), LessThanMQ );
		mate2SpecialAl = *mate2Set.rbegin();
		mate2SpecialAl.SpecialCode = specialCode2;
		mate2SpecialAl.NumMapped   = nMobAl2;
		mate2Set.clear();
	
	} else if ( nMobAl2 > 0 ) {
		vector<Alignment> newMate2Set;
		for ( vector<Alignment>::iterator ite = mate2Set.begin(); ite != mate2Set.end(); ++ite ) {
			
			if ( !ite->IsMappedSpecialReference ) {
				ite->CanBeMappedToSpecialReference = true;
				ite->SpecialCode = specialCode2;
				newMate2Set.push_back( *ite );
			} else {
				isMate2Special = true;
				mate2SpecialAl = ( ite->Quality >= mate2SpecialAl.Quality ) ? *ite : mate2SpecialAl;
				mate2SpecialAl.SpecialCode = specialCode2;
				mate2SpecialAl.NumMapped   = nMobAl2;
			}

		}
		mate2Set.clear();
		mate2Set = newMate2Set;
	}
}


// greater-than operator of mapping qualities
//inline bool CAlignmentThread::LessThanMQ ( const Alignment& al1, const Alignment& al2){
//	return al1.Quality < al2.Quality;
//}


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
		if(mForwardRead) {
			delete [] mForwardRead; 
			mForwardRead = NULL;
		}

		if(mReverseRead) {
			delete [] mReverseRead; 
			mReverseRead = NULL;
		}

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
	//if ( mSettings.SequencingTechnology == ST_454 ) {
		if(mFilters.UseMismatchFilter)        numExtensionBases = mFilters.MaxNumMismatches;
		if(mFilters.UseMismatchPercentFilter) numExtensionBases = (unsigned int)(queryLength * mFilters.MaxMismatchPercent);
		if(numExtensionBases < 2)             numExtensionBases = 2;
	//} else {
	//	numExtensionBases = queryLength;
	//}

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
			//al.IsReverseStrand = isFastHashRegionReverseStrand;
			//al.BaseQualities.Copy( qualities, queryLength );
			//if ( al.IsReverseStrand ) al.BaseQualities.Reverse();
	
			// perform a Smith-Waterman alignment
			AlignRegion(fastHashRegion, al, fastHashRead, queryLength, numExtensionBases);

			// add the alignment to the vector if it passes the filters
			if( ApplyReadFilters( al, query, qualities, queryLength ) ) {
				// the base qualities of SOLiD reads are attached in ApplyReadFilters
				//if( mFlags.EnableColorspace )
				//	al.BaseQualities.Copy(qualities, queryLength);
				alignments.Add(al);
			}

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
				if( ApplyReadFilters( al, query, qualities, queryLength ) ) {	
					// the base qualities of SOLiD reads are attached in ApplyReadFilters
					//if( mFlags.EnableColorspace )
					//	al.BaseQualities.Copy(qualities, queryLength);
					alignments.Add(al);
				}

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
					if( ApplyReadFilters(al, query, qualities, queryLength ) ) {
						// the base qualities of SOLiD reads are attached in ApplyReadFilters
						//if( mFlags.EnableColorspace ) {
						//	al.BaseQualities.Copy( qualities, queryLength);
						//	al.BaseQualities.Reverse();
						//}
						alignments.Add(al);
					}

					// increment our candidates counter
					mStatisticsCounters.AlignmentCandidates++;

					// check if we can prematurely stop
					if(!alignAllReads && (alignments.GetCount() > 1)) break;
				}
			}
		}

		// no alignments because of filtering
		if( alignments.IsEmpty() ) {
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
bool CAlignmentThread::ApplyReadFilters(Alignment& al, const char* bases, const char* qualities, const unsigned int queryLength) {

	unsigned int queryLength1 = mFlags.EnableColorspace ? queryLength + 1 : queryLength;
	
	// assuming this is a good read
	bool ret = true;

	unsigned short numNonAlignedBases = queryLength1 - al.QueryLength;
	al.BaseQualities.Copy( qualities + al.QueryBegin, al.QueryEnd - al.QueryBegin + 1 );

	// convert from colorspace to basespace
	if( mFlags.EnableColorspace ) {
		//cerr << al.NumMismatches;
		//al.BaseQualities.Copy( qualities + al.QueryBegin, al.QueryEnd - al.QueryBegin + 1 );
		if( al.IsReverseStrand ) al.BaseQualities.Reverse();
		if ( !mCS.ConvertAlignmentToBasespace( al ) ) ret = false;
		numNonAlignedBases = queryLength1 - al.QueryLength;
		//cerr << " " << al.NumMismatches << " " << al.QueryLength << endl;

	} else {
		//al.BaseQualities.Copy( qualities + al.QueryBegin, al.QueryEnd - al.QueryBegin + 1 );
	
		// don't count leading and lagging N's as mismatches
		unsigned int pos = 0;
		while( ( bases[pos] == 'N' ) && ( pos < queryLength1 ) ) {
			numNonAlignedBases--;
			pos++;
		}
		pos = queryLength1 - 1;
		while( ( bases[pos] == 'N' ) && ( pos >= 0 ) ) {
			numNonAlignedBases--;
			pos--;
		}
	}

	// calculate the total number of mismatches
	//const unsigned short numTotalMismatches = al.NumMismatches + (mFlags.UseAlignedReadLengthForMismatchCalculation ? 0 : numNonAlignedBases);
	const unsigned short numTotalMismatches = al.NumMismatches + ( (mFilters.UseMinAlignmentFilter || mFilters.UseMinAlignmentPercentFilter ) ? 0 : numNonAlignedBases);

	// check to see if this alignment meets the maximum mismatch threshold
	if(mFilters.UseMismatchFilter && (numTotalMismatches > mFilters.MaxNumMismatches)) ret = false;

	// check to see if this alignment meets the maximum mismatch threshold
	if(mFilters.UseMismatchPercentFilter) {
		double percentMismatch = (double)numTotalMismatches / (double)queryLength1;
		if(percentMismatch > mFilters.MaxMismatchPercent) ret = false; 
	}

	// check to see if this alignment meets the minimum percentage alignment threshold
	if(mFilters.UseMinAlignmentPercentFilter) {
		double percentageAligned = (double)al.QueryLength / (double)queryLength1;
		if(percentageAligned < mFilters.MinPercentAlignment) ret = false;
	}

	// check to see if this alignment meets the minimum alignment threshold
	if(mFilters.UseMinAlignmentFilter && (al.QueryLength < mFilters.MinAlignment)) ret = false;

	// aligned to special references
	if ( ret && mSReference.found )
		al.IsMappedSpecialReference = mReferenceSpecial[ al.ReferenceIndex ];

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
	if( !mFlags.IsAligningAllReads) sort(regions.begin(), regions.end(), SortHashRegionByLength() );
}

// settles the local Smith-Waterman window
bool CAlignmentThread::SettleLocalSearchRegion( const LocalAlignmentModel& lam, const unsigned int refIndex, const unsigned int uniqueBegin, const unsigned int uniqueEnd, unsigned int& localSearchBegin, unsigned int& localSearchEnd ) {

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
	if(begin == end) {
		localSearchBegin = 0;
		localSearchEnd   = 0;
		return false;
	} else {
		localSearchBegin = begin;
		localSearchEnd   = end;
		return true;
	}

}

// attempts to rescue the mate paired with a unique mate
bool CAlignmentThread::RescueMate(const LocalAlignmentModel& lam, const CMosaikString& bases, const unsigned int begin, const unsigned int end, const unsigned int refIndex, Alignment& al) {

/*
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
*/
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
	al.ReferenceBegin += begin - mReferenceBegin[refIndex];
	al.ReferenceEnd   += begin - mReferenceBegin[refIndex];

	if ( al.QueryBegin >= al.QueryEnd )
		return false;
	else
		// an alignment was performed
		return true;
}
