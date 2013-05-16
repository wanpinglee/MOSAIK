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
pthread_mutex_t CAlignmentThread::mSaveReadMutex;
pthread_mutex_t CAlignmentThread::mStatisticsMutex;
pthread_mutex_t CAlignmentThread::mStatisticsMapsMutex;
pthread_mutex_t CAlignmentThread::mSaveMultipleBamMutex;
pthread_mutex_t CAlignmentThread::mSaveSpecialBamMutex;
pthread_mutex_t CAlignmentThread::mSaveUnmappedBamMutex;

bool FilterMateOut ( const unsigned int length, char* basePtr ) {

	unsigned int count = 0;
	for ( unsigned int i = 0; i < length; ++i ) {
		if ( *basePtr == 'N' )
			++count;

		++basePtr;
	}

	return ( ( count / (float) length ) > 0.7 );
}

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
	const unsigned int          referenceOffset,
	string                      i_paired_end_ann_file,
	string                      i_single_end_ann_file)
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
	, mSW(CPairwiseUtilities::MatchScore, CPairwiseUtilities::MismatchScore, CPairwiseUtilities::GapOpenPenalty, CPairwiseUtilities::GapExtendPenalty, flags.NotCountGapAsMismatch)
	, mBSW(CPairwiseUtilities::MatchScore, CPairwiseUtilities::MismatchScore, CPairwiseUtilities::GapOpenPenalty, CPairwiseUtilities::GapExtendPenalty, settings.Bandwidth, flags.NotCountGapAsMismatch)
	, mSSW(10, 9, 15, 1)
	, mReferenceBegin(pRefBegin)
	, mReferenceEnd(pRefEnd)
	, mReferenceSpecies(pRefSpecies)
	, mReferenceSpecial(pRefSpecial)
	, mReadGroupsMap(pReadGroupsMap)
	, mReferenceOffset(referenceOffset)
	, paired_end_ann_file(i_paired_end_ann_file)
	, single_end_ann_file(i_single_end_ann_file)
{
	// set our flags
	if(algorithmMode == AlignerMode_ALL) mFlags.IsAligningAllReads = true;

	// assign the reference sequences to the colorspace utilities object
	mCS.SetReferenceSequences(pBsRefSeqs);

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
		pTD->ReferenceOffset,
		pTD->paired_end_ann_file,
		pTD->single_end_ann_file);

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
	pTD->pCounters->FailedHashMates         += at.mStatisticsCounters.FailedHashMates;
	pTD->pCounters->FailedHashMates_Rescue  += at.mStatisticsCounters.FailedHashMates_Rescue;
	pTD->pCounters->FilteredOutMates        += at.mStatisticsCounters.FilteredOutMates;
	pTD->pCounters->FilteredOutMates_Rescue += at.mStatisticsCounters.FilteredOutMates_Rescue;
	pTD->pCounters->ShortMates              += at.mStatisticsCounters.ShortMates;
	pTD->pCounters->ShortMates_Rescue       += at.mStatisticsCounters.ShortMates_Rescue;
	pTD->pCounters->TooManyNsMates          += at.mStatisticsCounters.TooManyNsMates;
	pTD->pCounters->TooManyNsMates_Rescue   += at.mStatisticsCounters.TooManyNsMates_Rescue;
	pTD->pCounters->MultipleMates           += at.mStatisticsCounters.MultipleMates;
	pTD->pCounters->MultipleMates_Rescue    += at.mStatisticsCounters.MultipleMates_Rescue;
	pTD->pCounters->UniqueMates             += at.mStatisticsCounters.UniqueMates;
	pTD->pCounters->UniqueMates_Rescue      += at.mStatisticsCounters.UniqueMates_Rescue;
	pTD->pCounters->UU                      += at.mStatisticsCounters.UU;
	pTD->pCounters->UM                      += at.mStatisticsCounters.UM;
	pTD->pCounters->UF                      += at.mStatisticsCounters.UF;
	pTD->pCounters->MM                      += at.mStatisticsCounters.MM;
	pTD->pCounters->MF                      += at.mStatisticsCounters.MF;
	pTD->pCounters->UX                      += at.mStatisticsCounters.UX;
	pTD->pCounters->MX                      += at.mStatisticsCounters.MX;
	pTD->pCounters->FF                      += at.mStatisticsCounters.FF;
	pTD->pCounters->FX                      += at.mStatisticsCounters.FX;
	pTD->pCounters->XX                      += at.mStatisticsCounters.XX;
	pTD->pCounters->UU_localRescue          += at.mStatisticsCounters.UU_localRescue;
	pTD->pCounters->UU_localConsistance     += at.mStatisticsCounters.UU_localConsistance;
	pTD->pCounters->UM_localRescue          += at.mStatisticsCounters.UM_localRescue;
	pTD->pCounters->UM_localConsistance     += at.mStatisticsCounters.UM_localConsistance;
	pTD->pCounters->MM_localRescue          += at.mStatisticsCounters.MM_localRescue;
	pTD->pCounters->MM_localConsistance     += at.mStatisticsCounters.MM_localConsistance;
	pthread_mutex_unlock(&mStatisticsMutex);

	// exit the thread
	pthread_exit((void*)0);
	return 0;
}

void CAlignmentThread::SearchLocalRegion(
    const vector<Alignment*>& anchorVector,
    const Mosaik::Mate& mate,
    CNaiveAlignmentSet* mateVector) {
	
	for (vector<Alignment*>::const_iterator uniqueIter = anchorVector.begin();
	    uniqueIter != anchorVector.end(); ++uniqueIter) {
	// extract the unique begin and end coordinates
	// Note that for multiply mappings, 
	// the best one alignment is at the end of the vector
	//vector<Alignment*>::const_iterator uniqueIter = anchorVector.end() - 1;

	const unsigned int refIndex    = (*uniqueIter)->ReferenceIndex;
	const unsigned int uniqueBegin = mReferenceBegin[refIndex] + (*uniqueIter)->ReferenceBegin;
	const unsigned int uniqueEnd   = mReferenceBegin[refIndex] + (*uniqueIter)->ReferenceEnd;

	// create the appropriate local alignment search model
	// NB: we caught other technologies during initialization
	LocalAlignmentModel lam;
	if((*uniqueIter)->IsReverseStrand) {
		if( alInfo.isUsingIllumina )
			lam.IsTargetBeforeUniqueMate = true;
		else if ( alInfo.isUsing454 )
			lam.IsTargetReverseStrand    = true;
		else if ( alInfo.isUsingSOLiD ) {
			lam.IsTargetBeforeUniqueMate = true;
			lam.IsTargetReverseStrand    = true;
		}
	} else {
		if( alInfo.isUsingIllumina )
			lam.IsTargetReverseStrand    = true;
		else if ( alInfo.isUsing454 )
			lam.IsTargetBeforeUniqueMate = true;
		else if ( alInfo.isUsingIlluminaLong ) {
			lam.IsTargetReverseStrand    = true;
			lam.IsTargetBeforeUniqueMate = true;
		}
	}

	unsigned int localSearchBegin = 0;
	unsigned int localSearchEnd   = 0;

	// settle the local search region
	bool settleLocalSearchWindow = SettleLocalSearchRegion( lam, refIndex, uniqueBegin, uniqueEnd, localSearchBegin, localSearchEnd );
	
	// check if the mate is already sitting in the local region
	bool isAlExisting = false;
	if ( settleLocalSearchWindow )
		isAlExisting = mateVector->CheckExistence(refIndex, localSearchBegin - mReferenceBegin[refIndex], localSearchEnd - mReferenceBegin[refIndex]);

	if ( !isAlExisting ) {
		Alignment al;

		if( RescueMate( lam, mate.Bases, localSearchBegin, localSearchEnd, refIndex, al ) ) {

			const char* pQualities = mate.Qualities.CData();
			const char* pBases = mate.Bases.CData();

			// add the alignment to the alignment set if it passes the filters
			if( ApplyReadFilters( al, pBases, pQualities, mate.Bases.Length() ) ) {
				al.WasRescued = true;
				mateVector->Add( al ); 
			}
		}
	}
	} // end for


}

void CAlignmentThread::UpdateSeStatistics ( const enum AlignmentStatusType& mateStatus, const Alignment& al ) {
	
	bool isMateRescued = al.WasRescued;
	if ( !al.IsMapped ) {
		switch( mateStatus ) {
			case ALIGNMENTSTATUS_FAILEDHASH:
				mStatisticsCounters.FailedHashMates++;
				if ( isMateRescued ) mStatisticsCounters.FailedHashMates_Rescue++;
			break;
			case ALIGNMENTSTATUS_TOOSHORT:
				mStatisticsCounters.ShortMates++;
				if ( isMateRescued ) mStatisticsCounters.ShortMates_Rescue++;
			break;
			case ALIGNMENTSTATUS_TOOMANYNS:
				mStatisticsCounters.TooManyNsMates++;
				if ( isMateRescued ) mStatisticsCounters.TooManyNsMates_Rescue++;
			break;
			case ALIGNMENTSTATUS_FILTEREDOUT:
				mStatisticsCounters.FilteredOutMates++;
				if ( isMateRescued ) mStatisticsCounters.FilteredOutMates_Rescue++;
			break;
			default:
			break;
		}
	} else {
		if ( al.NumMapped == 1 ) {
			mStatisticsCounters.UniqueMates++;
			if ( al.WasRescued )
				mStatisticsCounters.UniqueMates_Rescue++;
		} else {
			mStatisticsCounters.MultipleMates++;
			if ( al.WasRescued )
				mStatisticsCounters.MultipleMates_Rescue++;
		}
	}
}

// update statistics
void CAlignmentThread::UpdateStatistics ( 
	  const enum AlignmentStatusType& mate1Status
	, const enum AlignmentStatusType& mate2Status
	, const Alignment& al1
	, const Alignment& al2
	, const bool& isProperPair) {
	
	UpdateSeStatistics( mate1Status, al1 );
	if ( alInfo.isPairedEnd ) {
		UpdateSeStatistics( mate2Status, al2 );

		// update paired-end statistics
		bool UU = ( al1.NumMapped == 1 ) && ( al2.NumMapped == 1 );
		bool MM = ( al1.NumMapped > 1 )  && ( al2.NumMapped > 1 );
		bool UM = ( ( al1.NumMapped == 1 ) && ( al2.NumMapped > 1 ) )
			||( ( al1.NumMapped > 1 )  && ( al2.NumMapped == 1 ) );
		bool UF = ( ( al1.NumMapped == 1 ) && ( !al2.IsMapped && ( mate2Status == ALIGNMENTSTATUS_FILTEREDOUT ) ) )
			||( ( al2.NumMapped == 1 ) && ( !al1.IsMapped && ( mate1Status == ALIGNMENTSTATUS_FILTEREDOUT ) ) );
		bool UX = ( ( al1.NumMapped == 1 ) && ( !al2.IsMapped && ( mate2Status != ALIGNMENTSTATUS_FILTEREDOUT ) ) )
			||( ( al2.NumMapped == 1 ) && ( !al1.IsMapped && ( mate1Status != ALIGNMENTSTATUS_FILTEREDOUT ) ) );
		bool MF = ( ( al1.NumMapped > 1 )  && ( !al2.IsMapped && ( mate2Status == ALIGNMENTSTATUS_FILTEREDOUT ) ) )
			||( ( al2.NumMapped > 1 )  && ( !al1.IsMapped && ( mate1Status == ALIGNMENTSTATUS_FILTEREDOUT ) ) );
		bool MX = ( ( al1.NumMapped > 1 )  && ( !al2.IsMapped && ( mate2Status != ALIGNMENTSTATUS_FILTEREDOUT ) ) )
			||( ( al2.NumMapped > 1 )  && ( !al1.IsMapped && ( mate1Status != ALIGNMENTSTATUS_FILTEREDOUT ) ) );
		bool FF = ( !al1.IsMapped && ( mate1Status == ALIGNMENTSTATUS_FILTEREDOUT ) ) && ( !al2.IsMapped && ( mate2Status == ALIGNMENTSTATUS_FILTEREDOUT ) );
		bool FX = ( ( !al1.IsMapped && ( mate1Status == ALIGNMENTSTATUS_FILTEREDOUT ) ) && ( !al2.IsMapped && ( mate2Status == !ALIGNMENTSTATUS_FILTEREDOUT ) ) )
			||( ( !al1.IsMapped && ( mate1Status == !ALIGNMENTSTATUS_FILTEREDOUT ) ) && ( !al2.IsMapped && ( mate2Status == ALIGNMENTSTATUS_FILTEREDOUT ) ) );
		bool XX = ( !al1.IsMapped && ( mate1Status != ALIGNMENTSTATUS_FILTEREDOUT ) ) && ( !al2.IsMapped && ( mate2Status != ALIGNMENTSTATUS_FILTEREDOUT ) );
		
		if ( UU ) {
			mStatisticsCounters.UU++;

			if ( al1.WasRescued || al2.WasRescued )
				mStatisticsCounters.UU_localRescue++;
			else if ( isProperPair )
				mStatisticsCounters.UU_localConsistance++;
		}

		if ( UM ) {
			mStatisticsCounters.UM++;

			if ( al1.WasRescued || al2.WasRescued )
				mStatisticsCounters.UM_localRescue++;
			else if ( isProperPair )
				mStatisticsCounters.UM_localConsistance++;
		}

		if ( MM ) {
			mStatisticsCounters.MM++;

			if ( al1.WasRescued || al2.WasRescued )
				mStatisticsCounters.MM_localRescue++;
			else if ( isProperPair )
				mStatisticsCounters.MM_localConsistance++;
		}

		if ( UX )
			mStatisticsCounters.UX++;
		
		if ( MX )
			mStatisticsCounters.MX++;

		if ( XX )
			mStatisticsCounters.XX++;

		if ( UF )
			mStatisticsCounters.UF++;

		if ( MF )
			mStatisticsCounters.MF++;

		if ( FF )
			mStatisticsCounters.FF++;

		if ( FX )
			mStatisticsCounters.FX++;
	}


}

// save multiply alignment in buffer
void CAlignmentThread::SaveMultiplyAlignment(
    const vector<Alignment*>& mate1Set, 
    const vector<Alignment*>& mate2Set, 
    const Mosaik::Read& mr, 
    BamWriters* const pBams, 
    CStatisticsMaps* const pMaps ) {
	
	bool isMate1Multiple = ( mate1Set.size() > 1 ) ? true : false;
	bool isMate2Multiple = ( mate2Set.size() > 1 ) ? true : false;
	bool isMate1Empty    = ( mate1Set.size() == 0 ) ? true : false;
	bool isMate2Empty    = ( mate2Set.size() == 0 ) ? true : false;

	AlignmentStatusType mate1Status, mate2Status;
	// -om is enabled
	if (mFlags.OutputMultiplyComplete) {
		if (isMate1Multiple) {
			Alignment mateAl;
			if ( !isMate2Empty ) {
				mateAl = *(mate2Set[0]);
				mateAl.ReferenceIndex += mReferenceOffset;
			}
			vector<Alignment*> mate1SetTemp = mate1Set;
			for(vector<Alignment*>::iterator alIter = mate1SetTemp.begin(); alIter != mate1SetTemp.end(); ++alIter) {
				
				if ( !isMate2Empty )
					(*alIter)->SetPairFlagsAndFragmentLength(mateAl, 0, 0, mSettings.SequencingTechnology);

				(*alIter)->ReferenceIndex += mReferenceOffset;
				(*alIter)->RecalibratedQuality = (*alIter)->Quality;
				SetRequiredInfo( **alIter, mate1Status, mateAl, mr.Mate1, mr, !isMate2Empty, false, true, alInfo.isPairedEnd, true, !isMate2Empty);

				AlignmentBamBuffer buffer;
				buffer.al = **alIter;
				buffer.zaString        = (char) 0;
				buffer.noCigarMdNm     = false;
				buffer.notShowRnamePos = false;

				bamMultiplyBuffer.push( buffer );
			}
		}
		if (alInfo.isPairedEnd && isMate2Multiple) {
			Alignment mateAl;
			if (!isMate1Empty) {
				mateAl = *(mate1Set[0]);
				mateAl.ReferenceIndex += mReferenceOffset;
			}

			vector<Alignment*> mate2SetTemp = mate2Set;
			for(vector<Alignment*>::iterator alIter = mate2SetTemp.begin(); alIter != mate2SetTemp.end(); ++alIter) {
				if ( !isMate1Empty )
					(*alIter)->SetPairFlagsAndFragmentLength(mateAl, 0, 0, mSettings.SequencingTechnology);
				
				(*alIter)->ReferenceIndex += mReferenceOffset;
				(*alIter)->RecalibratedQuality = (*alIter)->Quality;
				SetRequiredInfo( **alIter, mate2Status, mateAl, mr.Mate2, mr, !isMate1Empty, false, false, alInfo.isPairedEnd, true, !isMate1Empty);

				AlignmentBamBuffer buffer;
				buffer.al = **alIter;
				buffer.zaString        = (char) 0;
				buffer.noCigarMdNm     = false;
				buffer.notShowRnamePos = false;

				bamMultiplyBuffer.push( buffer );
			}
		}
		// buffer is full; save and clear it
		if (bamMultiplyBuffer.size() > _bufferSize) {
			//AlignmentBamBuffer buffer;
			pthread_mutex_lock(&mSaveMultipleBamMutex);
			while(!bamMultiplyBuffer.empty()) {
				//buffer = bamMultiplyBuffer.front();
				pBams->mBam.SaveAlignment(bamMultiplyBuffer.front().al, 
				                          bamMultiplyBuffer.front().zaString.c_str(), 
							  bamMultiplyBuffer.front().noCigarMdNm, 
							  bamMultiplyBuffer.front().notShowRnamePos, 
							  alInfo.isUsingSOLiD);
				bamMultiplyBuffer.pop();
			}
			pthread_mutex_unlock(&mSaveMultipleBamMutex);
		}
	} else if (mFlags.OutputMultiplyIncomplete) {
		if (isMate1Multiple) {
			SimpleBamRecordBuffer buffer;
			for(vector<Alignment*>::const_iterator alIter = mate1Set.begin(); alIter != mate1Set.end(); ++alIter) {
				buffer.refIndex = (*alIter)->ReferenceIndex + mReferenceOffset;
				buffer.refBegin = (*alIter)->ReferenceBegin;
				buffer.refEnd   = (*alIter)->ReferenceEnd;
				bamMultiplySimpleBuffer.push( buffer );
			}
		}
		if (alInfo.isPairedEnd && isMate2Multiple) {
			SimpleBamRecordBuffer buffer;
			for(vector<Alignment*>::const_iterator alIter = mate2Set.begin(); alIter != mate2Set.end(); ++alIter) {
				buffer.refIndex = (*alIter)->ReferenceIndex + mReferenceOffset;
				buffer.refBegin = (*alIter)->ReferenceBegin;
				buffer.refEnd   = (*alIter)->ReferenceEnd;
				bamMultiplySimpleBuffer.push( buffer );
			}
		}
		// buffer is full; save and clear it
		if (bamMultiplySimpleBuffer.size() > _bufferSize) {
			//SimpleBamRecordBuffer buffer;
			pthread_mutex_lock(&mSaveMultipleBamMutex);
			while( !bamMultiplySimpleBuffer.empty() ) {
				//buffer = bamMultiplySimpleBuffer.front();
				pBams->mBam.SaveReferencePosition( bamMultiplySimpleBuffer.front().refIndex, 
				                                   bamMultiplySimpleBuffer.front().refBegin, 
								   bamMultiplySimpleBuffer.front().refEnd );
				bamMultiplySimpleBuffer.pop();
			}
			pthread_mutex_unlock(&mSaveMultipleBamMutex);
		}
	} else {
		// nothing
	}

}

void CAlignmentThread::SaveNClearBuffers( BamWriters* const pBams, CStatisticsMaps* const pMaps, MosaikReadFormat::CAlignmentWriter* const pOut ) {
	if ( !bamBuffer.empty() )
		WriteAlignmentBufferToFile( pBams, pMaps, pOut );

	if ( !archiveBuffer.empty() )
		WriteAlignmentBufferToFile( pBams, pMaps, pOut );

	if ( !bamMultiplyBuffer.empty() ) {
		//AlignmentBamBuffer buffer;
		pthread_mutex_lock(&mSaveMultipleBamMutex);
		while( !bamMultiplyBuffer.empty() ) {
			//buffer = bamMultiplyBuffer.front();
			pBams->mBam.SaveAlignment( bamMultiplyBuffer.front().al, 
			                           bamMultiplyBuffer.front().zaString.c_str(), 
						   bamMultiplyBuffer.front().noCigarMdNm, 
						   bamMultiplyBuffer.front().notShowRnamePos, 
						   alInfo.isUsingSOLiD );
			bamMultiplyBuffer.pop();
		}
		pthread_mutex_unlock(&mSaveMultipleBamMutex);
	}

	if ( !bamMultiplySimpleBuffer.empty() ) {
		SimpleBamRecordBuffer buffer;
		pthread_mutex_lock(&mSaveMultipleBamMutex);
		while( !bamMultiplySimpleBuffer.empty() ) {
			buffer = bamMultiplySimpleBuffer.front();
			pBams->mBam.SaveReferencePosition( buffer.refIndex, buffer.refBegin, buffer.refEnd );
			bamMultiplySimpleBuffer.pop();
		}
		pthread_mutex_unlock(&mSaveMultipleBamMutex);
	}

	if ( !bamSpecialBuffer.empty() ) {
		const bool processedBamData = true;
		//AlignmentBamBuffer buffer;
		pthread_mutex_lock(&mSaveSpecialBamMutex);
		while( !bamSpecialBuffer.empty() ) {
			//buffer = bamSpecialBuffer.front();
			//bamSpecialBuffer.pop();
			pBams->sBam.SaveAlignment( bamSpecialBuffer.front().al, 
			                           bamSpecialBuffer.front().zaString.c_str(), 
						   bamSpecialBuffer.front().noCigarMdNm, 
						   bamSpecialBuffer.front().notShowRnamePos, 
						   mFlags.EnableColorspace, 
						   processedBamData);
			bamSpecialBuffer.pop();
		}
		pthread_mutex_unlock(&mSaveSpecialBamMutex);
	}

}

// Save alignment in buffer
void CAlignmentThread::SaveBamAlignment( const Alignment& al, const char* zaString, const bool& noCigarMdNm, const bool& notShowRnamePos, const bool& isSpecial ) {
	AlignmentBamBuffer buffer;
	buffer.al              = al;
	buffer.noCigarMdNm     = noCigarMdNm;
	buffer.notShowRnamePos = notShowRnamePos;
	if ( zaString == (char)0 ) 
		buffer.zaString.clear();
	else
		buffer.zaString = zaString;
	
	if ( !noCigarMdNm ) {
		bamMisc.CreatePackedCigar( buffer.al, buffer.al.PackedCigar, buffer.al.NumCigarOperation, alInfo.isUsingSOLiD );
		buffer.al.MdString = mdTager.GetMdTag( buffer.al.Reference.CData(), buffer.al.Query.CData(),  buffer.al.Reference.Length() );
	}
	buffer.al.Query.Remove('-');

	if ( buffer.al.Query.Length() != 0 )
		bamMisc.EncodeQuerySequence( buffer.al.Query.CData(), buffer.al.EncodedQuery );

	// after this, Reference and Query should not be used again
	//   since they are already free.
	buffer.al.Reference.clearMemory();
	buffer.al.Query.clearMemory();

	if (isSpecial)
		bamSpecialBuffer.push( buffer );
	else
		bamBuffer.push( buffer );
}

inline void CAlignmentThread::SaveArchiveAlignment ( const Mosaik::Read& mr, const Alignment& al1, const Alignment& al2, const bool& isLongRead ){
	AlignmentArchiveBuffer buffer;
	buffer.mr  = mr;
	buffer.al1 = al1;
	buffer.al2 = al2;
	buffer.isLongRead = isLongRead;

	archiveBuffer.push( buffer );
}

// write special buffer to bam
void CAlignmentThread::WriteSpecialAlignmentBufferToFile( BamWriters* const pBams ) {
	//AlignmentBamBuffer buffer;
	const bool processedBamData = true;
	pthread_mutex_lock(&mSaveSpecialBamMutex);
	while( !bamSpecialBuffer.empty() ) {
		//buffer = bamSpecialBuffer.front();
		//bamSpecialBuffer.pop();
		const char* za = bamSpecialBuffer.front().zaString.empty() ? 0 : bamSpecialBuffer.front().zaString.c_str();
		pBams->sBam.SaveAlignment( bamSpecialBuffer.front().al, 
		                           za, 
					   bamSpecialBuffer.front().noCigarMdNm, 
					   bamSpecialBuffer.front().notShowRnamePos, 
					   mFlags.EnableColorspace, 
					   processedBamData);
		bamSpecialBuffer.pop();
	}
	pthread_mutex_unlock(&mSaveSpecialBamMutex);
}

// write alignment buffer to bam/archive
void CAlignmentThread::WriteAlignmentBufferToFile( BamWriters* const pBams, CStatisticsMaps* const pMaps, MosaikReadFormat::CAlignmentWriter* const pOut ) {
	// bam output
	if (!alInfo.isUsingLowMemory) {
		AlignmentBamBuffer buffer1, buffer2;
		Alignment dumpAl;

		const bool processedBamData = true;

		if ( !alInfo.isPairedEnd ) {
			pthread_mutex_lock(&mSaveReadMutex);
			while ( !bamBuffer.empty() ) {
				//buffer1 = bamBuffer.front();
				const char* za = bamBuffer.front().zaString.empty() ? 0 : bamBuffer.front().zaString.c_str();
				pBams->rBam.SaveAlignment(bamBuffer.front().al, 
				                          za, 
							  bamBuffer.front().noCigarMdNm, 
							  bamBuffer.front().notShowRnamePos, 
							  mFlags.EnableColorspace, 
							  processedBamData, 
							  mFlags.ReportZnTag);
				pMaps->SaveRecord(bamBuffer.front().al, dumpAl, alInfo.isPairedEnd, mSettings.SequencingTechnology);
				bamBuffer.pop();
			}
			pthread_mutex_unlock(&mSaveReadMutex);
		} else {
			pthread_mutex_lock(&mSaveReadMutex);
			while ( !bamBuffer.empty() ) {
				buffer1 = bamBuffer.front();
				bamBuffer.pop();
				buffer2 = bamBuffer.front();
				bamBuffer.pop();
				const char* za1 = buffer1.zaString.empty() ? 0 : buffer1.zaString.c_str();
				pBams->rBam.SaveAlignment(buffer1.al, za1, buffer1.noCigarMdNm, buffer1.notShowRnamePos, mFlags.EnableColorspace, processedBamData, mFlags.ReportZnTag);
				const char* za2 = buffer2.zaString.empty() ? 0 : buffer2.zaString.c_str();
				pBams->rBam.SaveAlignment(buffer2.al, za2, buffer2.noCigarMdNm, buffer2.notShowRnamePos, mFlags.EnableColorspace, processedBamData, mFlags.ReportZnTag);
				bool buffer1IsMate1 = buffer1.al.IsFirstMate;
				pMaps->SaveRecord((buffer1IsMate1 ? buffer1.al : buffer2.al), 
				                  (buffer1IsMate1 ? buffer2.al : buffer1.al), 
						  alInfo.isPairedEnd, 
						  mSettings.SequencingTechnology );
			}
			pthread_mutex_unlock(&mSaveReadMutex);
		}
	} else {
		//AlignmentArchiveBuffer buffer;
		pthread_mutex_lock(&mSaveReadMutex);
		while ( !archiveBuffer.empty() ) {
			//buffer = archiveBuffer.front();
			pOut->SaveRead(archiveBuffer.front().mr, 
			               archiveBuffer.front().al1, 
				       archiveBuffer.front().al2, 
				       archiveBuffer.front().isLongRead, 
				       true, 
				       alInfo.isPairedEnd, mFlags.SaveUnmappedBasesInArchive);
			archiveBuffer.pop();
		}
		pthread_mutex_unlock(&mSaveReadMutex);
	}

}

// aligns the read archive
void CAlignmentThread::AlignReadArchive(
	MosaikReadFormat::CReadReader* pIn, 
	MosaikReadFormat::CAlignmentWriter* pOut, 
	uint64_t* pReadCounter, 
	bool      isPairedEnd, 
	CStatisticsMaps* pMaps, 
	BamWriters*      pBams,
	unsigned char statMappingQuality) {

	alInfo.isUsing454          = (mSettings.SequencingTechnology == ST_454      ? true : false);
	alInfo.isUsingIllumina     = (mSettings.SequencingTechnology == ST_ILLUMINA ? true : false);
	alInfo.isUsingSOLiD        = (mSettings.SequencingTechnology == ST_SOLID    ? true : false);
	alInfo.isUsingIlluminaLong = (mSettings.SequencingTechnology == ST_ILLUMINA_LONG    ? true : false);
	alInfo.isPairedEnd         = isPairedEnd;
	alInfo.isUsingLowMemory    = mFlags.UseArchiveOutput;

	_bufferSize = 1000;

	// catch unsupported local alignment search sequencing technologies
	if( mFlags.UseLocalAlignmentSearch && ( !alInfo.isUsing454 && !alInfo.isUsingIllumina && !alInfo.isUsingIlluminaLong && !alInfo.isUsingSOLiD ) ) {
		cout << "ERROR: This sequencing technology is not currently supported for local alignment search." << endl;
		exit(1);
	}

	// initialize our status variables
	enum AlignmentStatusType mate1Status, mate2Status;

	// derive the minimum span length
	unsigned int minSpanLength = mSettings.HashSize;
	if(mFlags.IsUsingAlignmentCandidateThreshold && (mSettings.AlignmentCandidateThreshold > minSpanLength)) 
		minSpanLength = mSettings.AlignmentCandidateThreshold;

	// create neural networks
	if (!alInfo.isUsingLowMemory)
		mqCalculator.Open(paired_end_ann_file, single_end_ann_file);

	// keep reading until no reads remain
	Mosaik::Read mr;
	CNaiveAlignmentSet mate1Alignments(mReferenceLength, (alInfo.isUsingIllumina || alInfo.isUsingSOLiD || alInfo.isUsingIlluminaLong ) );
	CNaiveAlignmentSet mate2Alignments(mReferenceLength, (alInfo.isUsingIllumina || alInfo.isUsingSOLiD || alInfo.isUsingIlluminaLong ) );

	unsigned int alignmentBufferSize = 1000;

	// fragment length window
	int minFl = mSettings.MedianFragmentLength - mSettings.LocalAlignmentSearchRadius;
	int maxFl = mSettings.MedianFragmentLength + mSettings.LocalAlignmentSearchRadius;

	while(true) {

		// =============================
		// load reads from input archive
		// =============================
		bool hasMoreReads = true;
		if ( readBuffer.empty() ) {
			unsigned int nReadBuffer = 0;
			mr.clear();
			pthread_mutex_lock(&mGetReadMutex);
			hasMoreReads = pIn->LoadNextRead(mr);
			while( hasMoreReads && ( nReadBuffer < alignmentBufferSize ) ) {
				readBuffer.push( mr );
				++nReadBuffer;
				mr.clear();
				*pReadCounter = *pReadCounter + 1;
				hasMoreReads = pIn->LoadNextRead(mr);
			}
			if ( hasMoreReads ) {
				readBuffer.push( mr );
				*pReadCounter = *pReadCounter + 1;
			}
			pthread_mutex_unlock(&mGetReadMutex);
		}

		// quit if we've processed all of the reads
		if ( readBuffer.empty() && !hasMoreReads ) {
			SaveNClearBuffers( pBams, pMaps, pOut );
			break;
		}
		
		// ===========================
		// take a read from the buffer
		// ===========================
		mr.clear();
		mr = readBuffer.front();
		readBuffer.pop();

		#ifdef VERBOSE_DEBUG
		cerr << "Read: " << mr.Name.CData() << endl;
		#endif

		// specify if this is a paired-end read
		const unsigned short numMate1Bases = (unsigned short)mr.Mate1.Bases.Length();
		const unsigned short numMate2Bases = (unsigned short)mr.Mate2.Bases.Length();

		const bool isTooManyNMate1 = FilterMateOut( numMate1Bases, mr.Mate1.Bases.Data() );
		const bool isTooManyNMate2 = FilterMateOut( numMate2Bases, mr.Mate2.Bases.Data() );

		if ( isTooManyNMate1 )
			mate1Status = ALIGNMENTSTATUS_TOOMANYNS;
		if ( isTooManyNMate2 )
			mate2Status = ALIGNMENTSTATUS_TOOMANYNS;

		const bool areBothMatesPresent = ( ( ( numMate1Bases != 0 ) && ( numMate2Bases != 0 ) ) ? true : false );

		// ====================
		// align the first mate
		// ====================

		bool isMate1Aligned = false;
		int numMate1Hashes  = 0;
		mate1Alignments.Clear();
		if( numMate1Bases != 0 && !isTooManyNMate1 ) {
			#ifdef VERBOSE_DEBUG
				cerr << "=== Align mate1 ===" << endl;
			#endif
			// align the read
			if (AlignRead(mate1Alignments, mr.Mate1.Bases.CData(), mr.Mate1.Qualities.CData(), numMate1Bases, mate1Status, &numMate1Hashes)) 
				isMate1Aligned = true;
			#ifdef VERBOSE_DEBUG
				if (isMate1Aligned) cerr << "mate1 is mapped." << endl;
				else cerr << "mate1 is NOT mapped!" << endl;
			#endif
		}

		// =====================
		// align the second mate
		// =====================

		bool isMate2Aligned = false;
		int numMate2Hashes  = 0;
		mate2Alignments.Clear();
		if( numMate2Bases != 0 && !isTooManyNMate2 ) {
			#ifdef VERBOSE_DEBUG
				cerr << "=== Align mate2 ===" << endl;
			#endif
			// align the read
			if (AlignRead(mate2Alignments, mr.Mate2.Bases.CData(), mr.Mate2.Qualities.CData(), numMate2Bases, mate2Status, &numMate2Hashes)) 
				isMate2Aligned = true;
			#ifdef VERBOSE_DEBUG
				if (isMate2Aligned) cerr << "mate2 is mapped." << endl;
				else cerr << "mate2 is NOT mapped!" << endl;
			#endif
		}

		// ======================
		// local alignment search
		// ======================

		vector<Alignment*> mate1Set, mate2Set;

		// we can only perform a local alignment search if both mates are present
		if(areBothMatesPresent) {
			// search local region
			if (mFlags.UseLocalAlignmentSearch && !isTooManyNMate2) {
				#ifdef VERBOSE_DEBUG
					cerr << "=== Local Search mate2 ===" << endl;
				#endif
				mate1Alignments.GetSet(&mate1Set);
				if (mate1Alignments.IsUnique()) {
					SearchLocalRegion(mate1Set, mr.Mate2, &mate2Alignments);
				} else if (mate1Alignments.IsMultiple()){ // do local search for some good multiply mappings
					bool considerMate1Unique = TreatBestAsUnique(&mate1Set, numMate1Bases);
					if (considerMate1Unique) // TreatBestAsUnique puts candidates in mate1Set
						SearchLocalRegion(mate1Set, mr.Mate2, &mate2Alignments);
				} // end if-else
			}

			// search local region
			if (mFlags.UseLocalAlignmentSearch && !isTooManyNMate1) {
				#ifdef VERBOSE_DEBUG
					cerr << "=== Local Search mate1 ===" << endl;
				#endif
				mate2Alignments.GetSet(&mate2Set);
				if (mate2Alignments.IsUnique()) {
					SearchLocalRegion(mate2Set, mr.Mate1, &mate1Alignments);
				} else if (mate2Alignments.IsMultiple()){ // do local search for some good multiply mappings
					bool considerMate2Unique = TreatBestAsUnique(&mate2Set, numMate2Bases);
					if (considerMate2Unique) // TreatBestAsUnique puts candidates in mate1Set
						SearchLocalRegion(mate2Set, mr.Mate1, &mate1Alignments);
				} // end if-else
			}
		}

		// process alignments mapped in special references and delete them in vectors
		mate1Alignments.GetSet(&mate1Set);
		mate2Alignments.GetSet(&mate2Set);
		//bool isLongRead = mate1Alignments.HasLongAlignment() || mate2Alignments.HasLongAlignment();


		// For low-memory, we don't remove special alignment here,
		// the special alignments will be considered when merging archives
		Alignment mate1SpecialAl, mate2SpecialAl;
		bool isMate1Special = false, isMate2Special = false;
		if (mSReference.found && !alInfo.isUsingLowMemory) {
			ProcessSpecialAlignment(&mate1Set, &mate1SpecialAl, &isMate1Special);
			if (isPairedEnd)
			  ProcessSpecialAlignment(&mate2Set, &mate2SpecialAl, &isMate2Special);
		}
		
		// deleting special alignments may let mate1Alignments or mate2Alignments become empty
		// so we have to check them again
		isMate1Aligned = !mate1Set.empty();
		isMate2Aligned = !mate2Set.empty();

		const bool isMate1Unique = ( mate1Set.size() == 1 ) ? true: false;
		const bool isMate2Unique = ( mate2Set.size() == 1 ) ? true: false;
		const bool isMate1Multiple = ( mate1Set.size() > 1 ) ? true: false;
		const bool isMate2Multiple = ( mate2Set.size() > 1 ) ? true: false;
		const bool isMate1Empty = mate1Set.empty();
		const bool isMate2Empty = mate2Set.empty();


		// ===================================
		// Save alignments to BAMs or archives
		// ===================================

		// save chromosomes and positions of multiple alignments in bam
		if (mFlags.SaveMultiplyBam && (mFlags.OutputMultiplyComplete || mFlags.OutputMultiplyIncomplete)) {
			SaveMultiplyAlignment( mate1Set, mate2Set, mr, pBams, pMaps );
		}
		
		// UU, UM, and MM pair
		if ( ( isMate1Unique && isMate2Unique )
			|| ( isMate1Unique && isMate2Multiple )
			|| ( isMate1Multiple && isMate2Unique )
			|| ( isMate1Multiple && isMate2Multiple ) ) {
	
			Alignment al1 = *mate1Set[0], al2 = *mate2Set[0];
			if ( ( isMate1Unique && isMate2Multiple )
				|| ( isMate1Multiple && isMate2Unique )
				|| ( isMate1Multiple && isMate2Multiple ) )
				// After selecting, only best and second best (if there) will be kept in mate1Set and mate2Set
				// The function also sets NumMapped of alignments
				BestNSecondBestSelection::Select(al1, al2, mate1Set, mate2Set, mSettings.MedianFragmentLength, 
				    mSettings.SequencingTechnology, numMate1Bases, numMate2Bases, true, true, true);

			bool properPair1 = false, properPair2 = false;
			al1.IsFirstMate = true;
			al2.IsFirstMate = false;
			// MM pair is always an improper pair
			properPair1 = al1.SetPairFlagsAndFragmentLength( al2, minFl, maxFl, mSettings.SequencingTechnology );
			properPair2 = al2.SetPairFlagsAndFragmentLength( al1, minFl, maxFl, mSettings.SequencingTechnology );

			// sanity checker
			if ( properPair1 != properPair2 ) {
				cout << "ERROR: An inconsistent proper pair is found." << endl;
				exit(1);
			}

			al1.NumHash = numMate1Hashes;
			al2.NumHash = numMate2Hashes;
			SetRequiredInfo(al1, mate1Status, al2, mr.Mate1, mr, true, properPair1, true, isPairedEnd, true, true);
			SetRequiredInfo(al2, mate2Status, al1, mr.Mate2, mr, true, properPair2, false, isPairedEnd, true, true);

			if (!alInfo.isUsingLowMemory) {
				al1.RecalibratedQuality = GetMappingQuality(al1, al1.QueryLength, al2, al2.QueryLength);
				al2.RecalibratedQuality = GetMappingQuality(al2, al2.QueryLength, al1, al1.QueryLength);
			}

			// Since Reference Begin may be changed, applying the following function to reset fragment length is necessary.
			if ( mFlags.EnableColorspace && ( !isMate1Multiple || !isMate2Multiple ) ) {
				al1.SetPairFlagsAndFragmentLength( al2, minFl, maxFl, mSettings.SequencingTechnology );
				al2.SetPairFlagsAndFragmentLength( al1, minFl, maxFl, mSettings.SequencingTechnology );
			}

			if (alInfo.isUsingLowMemory) {
				//bool isLongRead = mate1Alignments.HasLongAlignment() || mate2Alignments.HasLongAlignment();
				bool isLongRead = ( ( al1.QueryEnd > 255 ) || ( al2.QueryEnd > 255 ) ) ? true : false;
				isLongRead |= ((al1.Reference.Length() > 255) || (al2.Reference.Length() > 255));
				isLongRead |= ( ( al1.CsQuery.size() > 255 ) || ( al2.CsQuery.size() > 255 ) );
				SaveArchiveAlignment( mr, al1, al2, isLongRead );
			} else {
				if (isMate2Special) {
					Alignment genomicAl = al1;
					Alignment specialAl = mate2SpecialAl;
					specialAl.NumMapped = al2.NumMapped;

					SetRequiredInfo(specialAl, mate2Status, genomicAl, mr.Mate2, mr, true, false, false, isPairedEnd, true, true);
				
					const char *zas1Tag = za1.GetZaTag(genomicAl, specialAl, true);
					SaveBamAlignment(genomicAl, zas1Tag, false, false, true);
					const char *zas2Tag = za2.GetZaTag(specialAl, genomicAl, false);
					SaveBamAlignment(specialAl, zas2Tag, false, false, true);
				}
				if (isMate1Special) {
					Alignment genomicAl = al2;
					Alignment specialAl = mate1SpecialAl;
					specialAl.NumMapped = al1.NumMapped;

					SetRequiredInfo(specialAl, mate1Status, genomicAl, mr.Mate1, mr, true, false, true, isPairedEnd, true, true);
	
					const char *zas1Tag = za1.GetZaTag(genomicAl, specialAl, false);
					SaveBamAlignment(genomicAl, zas1Tag, false, false, true);
					const char *zas2Tag = za2.GetZaTag(specialAl, genomicAl, true);
					SaveBamAlignment(specialAl, zas2Tag, false, false, true);
				}

				const char* zaTag1 = za1.GetZaTag(al1, al2, true);
				SaveBamAlignment(al1, zaTag1, false, false, false);
				const char* zaTag2 = za2.GetZaTag(al2, al1, false);
				SaveBamAlignment(al2, zaTag2, false, false, false);
			}

			UpdateStatistics( mate1Status, mate2Status, al1, al2, properPair1 );

		// UX and MX pair
		} else if ( ( isMate1Empty || isMate2Empty )
			&&  !( isMate1Empty && isMate2Empty ) ) {

			Alignment al1, al2, unmappedAl;
			if ( isMate1Multiple || isMate2Multiple ) 
				BestNSecondBestSelection::Select( al1, al2, mate1Set, mate2Set, mSettings.MedianFragmentLength, 
				    mSettings.SequencingTechnology, numMate1Bases, numMate2Bases, true, true, true);

			bool isFirstMate;
			Alignment al;
			if ( !mate1Set.empty() ) {
				isFirstMate = true;
				al = isMate1Multiple ? al1 : *mate1Set[0];
			} else if ( !mate2Set.empty() ) {
				isFirstMate = false;
				al = isMate2Multiple ? al2 : *mate2Set[0];
				if ( !isPairedEnd ) {
					cout << "ERROR: The sequence technology is single-end, but second mate is aligned." << endl;
					exit(1);
				}
			} else {
				cout << "ERROR: Both mates are empty after applying best and second best selection." << endl;
				exit(1);
			}

			al.NumHash         = isFirstMate ? numMate1Hashes : numMate2Hashes;
			unmappedAl.NumHash = !isFirstMate ? numMate1Hashes : numMate2Hashes;
			SetRequiredInfo( al, (isFirstMate ? mate1Status : mate2Status), unmappedAl, (isFirstMate ? mr.Mate1 : mr.Mate2)
				, mr, false, false, isFirstMate, isPairedEnd, true, false);
			if (isPairedEnd)
				SetRequiredInfo(unmappedAl, (isFirstMate ? mate2Status : mate1Status ),
				    al, (isFirstMate ? mr.Mate2 : mr.Mate1), mr, true, false, !isFirstMate, isPairedEnd, false, true );

			if (!alInfo.isUsingLowMemory) {
				al.RecalibratedQuality = GetMappingQuality(al, al.QueryLength);
			}

			if (alInfo.isUsingLowMemory) {
				bool isLongRead = ( ( al.QueryEnd > 255 ) || ( unmappedAl.QueryEnd > 255 ) ) ? true : false;
				isLongRead |= ((al.Reference.Length() > 255) || (unmappedAl.Reference.Length() > 255));
				isLongRead |= ( ( al.CsQuery.size() > 255 ) || ( unmappedAl.CsQuery.size() > 255 ) );
				SaveArchiveAlignment( mr, ( isFirstMate ? al : unmappedAl ), ( isFirstMate ? unmappedAl : al ), isLongRead );
			} else {
				// show the original MQs in ZAs, and zeros in MQs fields of a BAM
				if (isPairedEnd) {
					unmappedAl.ReferenceBegin = al.ReferenceBegin;
					unmappedAl.ReferenceIndex = al.ReferenceIndex;
					const char* zaTag1 = za1.GetZaTag(al, unmappedAl, isFirstMate, !isPairedEnd, true);
					SaveBamAlignment(al, zaTag1, false, false, false);
					const char* zaTag2 = za2.GetZaTag(unmappedAl, al, !isFirstMate, !isPairedEnd, false);
					SaveBamAlignment(unmappedAl, zaTag2, true, false, false);
				} else {
					const char* zaTag1 = za1.GetZaTag(al, unmappedAl, isFirstMate, !isPairedEnd, true);
					SaveBamAlignment(al, zaTag1, false, false, false);
				}

				// store special hits
				if (isMate1Special) {
					Alignment specialAl = mate1SpecialAl;
					SetRequiredInfo( specialAl, mate1Status, al, mr.Mate1, mr, !isFirstMate, false, true, isPairedEnd, true, !isFirstMate );
					if (isFirstMate) { // the other mate is missing
					  specialAl.NumMapped = al.NumMapped;
					  const char *zas2Tag = za2.GetZaTag(specialAl, unmappedAl, true, !isPairedEnd, true);
					  SaveBamAlignment(specialAl, zas2Tag, false, false, true);
					} else if (isPairedEnd){ // the mate is mapped; myself has special alignments only
					  const char *zas1Tag = za1.GetZaTag(al, specialAl, false, !isPairedEnd, false);
					  SaveBamAlignment(al, zas1Tag, false, false, true);
					  const char *zas2Tag = za2.GetZaTag(specialAl, al, true, !isPairedEnd, false);
					  SaveBamAlignment(specialAl, zas2Tag, false, false, true);
					}
				}
				
				if (isPairedEnd && isMate2Special) {
					// store special hits
					Alignment specialAl = mate2SpecialAl;
					SetRequiredInfo( specialAl, mate2Status, al, mr.Mate2, mr, isFirstMate, false, false, isPairedEnd, true, isFirstMate );
					if (!isFirstMate) { // the other mate is missing
					  specialAl.NumMapped = al.NumMapped;
					  const char *zas2Tag = za2.GetZaTag(specialAl, unmappedAl, false, !isPairedEnd, true);
					  SaveBamAlignment(specialAl, zas2Tag, false, false, true);
					} else { // the mate is mapped; myself has special alignments only
					  const char *zas1Tag = za1.GetZaTag(al, specialAl, true, !isPairedEnd, false);
					  SaveBamAlignment(al, zas1Tag, false, false, true);
					  const char *zas2Tag = za2.GetZaTag(specialAl, al, false, !isPairedEnd, false);
					  SaveBamAlignment(specialAl, zas2Tag, false, false, true);
					}
				}
			}

			UpdateStatistics( ( isFirstMate ? mate1Status : mate2Status ) , ( isFirstMate ? mate2Status : mate1Status ), al, unmappedAl, false );
		// XX
		} else if ( isMate1Empty && isMate2Empty ) {
			
			Alignment unmappedAl1, unmappedAl2;
			unmappedAl1.NumHash = numMate1Hashes;
			unmappedAl2.NumHash = numMate2Hashes;
			SetRequiredInfo( unmappedAl1, mate1Status, unmappedAl2, mr.Mate1, mr, false, false, true, isPairedEnd, false, false );
			SetRequiredInfo( unmappedAl2, mate2Status, unmappedAl1, mr.Mate2, mr, false, false, false, isPairedEnd, false, false );

			if (alInfo.isUsingLowMemory) {

				bool isLongRead = ( ( unmappedAl1.QueryEnd > 255 ) || ( unmappedAl2.QueryEnd > 255 ) ) ? true : false;
				isLongRead |= ((unmappedAl1.Reference.Length() > 255) || (unmappedAl2.Reference.Length() > 255));
				isLongRead |= ( ( unmappedAl1.CsQuery.size() > 255 ) || ( unmappedAl2.CsQuery.size() > 255 ) );
				SaveArchiveAlignment( mr, unmappedAl1, unmappedAl2, isLongRead );
			} else {
				// store special hits
				if ( isMate1Special ) {
					Alignment specialAl = mate1SpecialAl;
					SetRequiredInfo( specialAl, mate1Status, unmappedAl2, mr.Mate1, mr, false, false, true, isPairedEnd, true, false );
					const char *zas2Tag = za2.GetZaTag( specialAl, unmappedAl2, true, !isPairedEnd, true );
					SaveBamAlignment(specialAl, zas2Tag, false, false, true);
				}

				if ( isPairedEnd ) {
					if ( isMate2Special ) {
						Alignment specialAl = mate2SpecialAl;
						SetRequiredInfo( specialAl, mate2Status, unmappedAl1, mr.Mate2, mr, false, false, false, isPairedEnd, true, false );
						const char *zas2Tag = za2.GetZaTag( specialAl, unmappedAl1, false, !isPairedEnd, true );
						SaveBamAlignment(specialAl, zas2Tag, false, false, true);
					}
					SaveBamAlignment(unmappedAl1, 0, true, false, false);
					SaveBamAlignment(unmappedAl2, 0, true, false, false);

				} else {
					SaveBamAlignment(unmappedAl1, 0, true, false, false);
				}

			}

			UpdateStatistics( mate1Status, mate2Status, unmappedAl1, unmappedAl2, false );

		} else {
			cout << "ERROR: Unknown pairs." << endl;
			//cout << mate1Alignments.GetCount() << "\t" << mate2Alignments.GetCount() << endl;
			exit(1);
		}

		if (!alInfo.isUsingLowMemory) {
			if ( bamBuffer.size() > alignmentBufferSize )
				WriteAlignmentBufferToFile( pBams, pMaps, pOut );
		}
		else { 
			if ( archiveBuffer.size() > alignmentBufferSize )
				WriteAlignmentBufferToFile( pBams, pMaps, pOut );
		}
		

		if ( bamSpecialBuffer.size() > alignmentBufferSize ) {
			WriteSpecialAlignmentBufferToFile( pBams );
		}

	} // end while
	
	SaveNClearBuffers(pBams, pMaps, pOut);
}

unsigned char CAlignmentThread::GetMappingQuality (const Alignment& al, 
                                                   const int& al_length) {

	mate1Ann.read_length   = al_length;
	mate1Ann.swScore       = al.SwScore;
	mate1Ann.nextSwScore   = al.NextSwScore;
	mate1Ann.longest_match = al.NumLongestMatchs;
	mate1Ann.entropy       = al.Entropy;
	mate1Ann.numMappings   = al.NumMapped;
	mate1Ann.numHashes     = al.NumHash;
	
	return mqCalculator.GetQualitySe(mate1Ann);
}

unsigned char CAlignmentThread::GetMappingQuality (const Alignment& al1, 
                                                   const int& al1_length, 
						   const Alignment& al2,
						   const int& al2_length) {
 	//int fl = (al1.ReferenceBegin > al2.ReferenceBegin) ? (al1.ReferenceBegin - al2.ReferenceBegin) : (al2.ReferenceBegin - al1.ReferenceBegin);
	//fl += al1.Query.Length();
	//fl = abs(mSettings.MedianFragmentLength - fl);
	int flDiff = mSettings.MedianFragmentLength - abs(al1.FragmentLength);
	//if (al1.ReferenceIndex != al2.ReferenceIndex)
	//	flDiff = INT_MAX - 1;
	//else
		flDiff = abs(flDiff);

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
	return mqCalculator.GetQualityPe(mate1Ann, mate2Ann, flDiff);
}

// treat the best alignment as an unique mapping and than turn on local search
bool CAlignmentThread::TreatBestAsUnique (vector<Alignment*>* mateSet, const unsigned int& readLength) {
	sort(mateSet->begin(), mateSet->end(), Alignment_LessThanMq());

	// Note that there are at least two alignments
	vector<Alignment*>::reverse_iterator rit = mateSet->rbegin();
	Alignment* bestAl  = *rit;
	unsigned short mq1 = (*rit)->Quality;
	float swScore1     = (*rit)->SwScore;
	++rit;
	Alignment* secondBestAl = *rit;
	unsigned short mq2      = (*rit)->Quality;
	float swScore2          = (*rit)->SwScore;

	mateSet->clear();

	if ( swScore1 > ( readLength * 9 ) ) {
		if (swScore1 == swScore2) {
			mateSet->push_back(bestAl);
			mateSet->push_back(secondBestAl);
		} else {
			mateSet->push_back(bestAl);
		}
		return true;
	}

	if ( ( mq1 > mFilters.LocalAlignmentSearchHighMqThreshold ) && ( mq2 < mFilters.LocalAlignmentSearchLowMqThreshold ) ) {
		mateSet->push_back(bestAl);
		return true;
	} else {
		return false;
	}
}

// Set the required information and flag for alignments
// Note: Don't apply this function for SOLiD alignments more than once
void CAlignmentThread::SetRequiredInfo (
	Alignment& al,
	const AlignmentStatusType& status,
	Alignment& mate,
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
	al.IsFilteredOut          = ( status == ALIGNMENTSTATUS_FILTEREDOUT );

	map<unsigned int, MosaikReadFormat::ReadGroup>::iterator rgIte;
	rgIte = mReadGroupsMap->find( r.ReadGroupCode );
	// sanity check
	if ( rgIte == mReadGroupsMap->end() ) {
		cout << "ERROR: ReadGroup cannot be found." << endl;
		exit(1);
	}	
	else 
		al.ReadGroup = rgIte->second.ReadGroupID;

	// calculate entropy
	if (isItselfMapped && !alInfo.isUsingLowMemory) 
		al.Entropy = entropy_.shannon_H((char*) m.Bases.CData(), m.Bases.Length());

	// fill out the alignment
	// not SOLiD
	if ( !mFlags.EnableColorspace ) {
		// copy qualities
		al.BaseQualities = m.Qualities;

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
			al.MappedLength = 0;
			if ((!alInfo.isUsingLowMemory) 
			    || (alInfo.isUsingLowMemory && mFlags.SaveUnmappedBasesInArchive)) {
				al.Query = m.Bases;
				al.Reference.Copy( 'Z', al.Query.Length() );
				al.IsJunk = false;
			} else {
				al.IsJunk = true;
			}
		}
		else {	
			string queryString, referenceString;
			// patch bases
			if ( patchStartLen > 0 ) {
				queryString.append( patchBases.CData(), patchStartLen );
				referenceString.append( patchStartLen, 'Z' );
				//al.Query.Prepend    ( patchBases.CData(), patchStartLen );
				//al.Reference.Prepend( 'Z', patchStartLen );
			}

			queryString.append( al.Query.CData() );
			referenceString.append( al.Reference.CData() );

			if ( patchEndLen > 0 ) {
				const unsigned int length = patchBases.Length();
				const unsigned int start  = length - patchEndLen;
				const char* startPoint    = patchBases.CData() + start;
				// sanity check
				if ( length > patchBases.Length() ) {
					cout << "ERROR: The soft chip position is wrong." << endl;
					exit(1);
				}

				queryString.append( startPoint, patchEndLen );
				referenceString.append( patchEndLen, 'Z' );
				//al.Query.Append    ( startPoint, patchEndLen );
				//al.Reference.Append( 'Z', patchEndLen );
			}

			al.Query = queryString.c_str();
			al.Reference = referenceString.c_str();

			char* qPtr = al.Query.Data();;
			char* rPtr = al.Reference.Data();
			al.MappedLength = 0;
			for ( uint32_t i = 0; i < al.Query.Length(); ++i ) {
				if ( ( *qPtr != '-' ) && ( *rPtr != 'Z' ) )
					++al.MappedLength;
				++qPtr;
				++rPtr;
			}

		}

		al.QueryBegin = 0;
		al.QueryEnd   = m.Bases.Length() - 1;
		al.QueryLength = al.QueryEnd - al.QueryBegin + 1;

	}
	else {
		// fill out Colorspace raw bases and qualites
		al.CsQuery.clear();
		al.CsBaseQualities.clear();
		if ( (!alInfo.isUsingLowMemory) 
		    || (alInfo.isUsingLowMemory && mFlags.SaveUnmappedBasesInArchive)) {
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
			al.MappedLength = 0;
		} else {
			
			//string queryString, referenceString, qualityString;
			
			// patch N's
			--al.QueryEnd;  //CS query end
			const unsigned int readLength    = m.Bases.Length();
			const unsigned int patchStartLen = al.IsReverseStrand ? readLength - al.QueryEnd - 1 : al.QueryBegin;
			const unsigned int patchEndLen   = al.IsReverseStrand ? al.QueryBegin : readLength - al.QueryEnd - 1;
			// sanity checker
			if ( al.QueryEnd > ( readLength - 1 ) ) {
				cout << "ERROR: The aligned length is larger than the read length." << endl;
				cout << "       Query end: " << al.QueryEnd << "\tread length: " << readLength << "\t# of mappings" << al.NumMapped << endl;
				cout << "       Read name: " << r.Name.CData() << endl;
				cout << "       Aligned bases :" << al.Query.CData() << endl;
				exit(1);
			}
			if ( patchStartLen == 0 ) {
				al.Query.TrimBegin(1);
				al.BaseQualities.TrimBegin(1);
				al.Reference.TrimBegin(1);
				al.ReferenceBegin++;
				mate.MateReferenceBegin++;
			}else if ( patchStartLen > 1 ) {
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

			char* qPtr = al.Query.Data();;
			char* rPtr = al.Reference.Data();
			al.MappedLength = 0;
			for ( uint32_t i = 0; i < al.Query.Length(); ++i ) {
				if ( ( *qPtr != '-' ) && ( *rPtr != 'Z' ) )
					++al.MappedLength;
				++qPtr;
				++rPtr;
			}
		}
	}

}

// handle and then delete special alignments
// also assign two special characters for the general alignments
void CAlignmentThread::ProcessSpecialAlignment ( 
    vector<Alignment*>* mate1Set, 
    Alignment* mate1SpecialAl, 
    bool* isMate1Special ) {

  unsigned int nMobAl = 0;
  string specialCode;
  specialCode.resize(3);
  for (vector<Alignment*>::const_iterator ite = mate1Set->begin(); ite != mate1Set->end(); ++ite) {
    if ((*ite)->IsMappedSpecialReference) {
      ++nMobAl;
      char* tempCode = mReferenceSpecies[(*ite)->ReferenceIndex];
      specialCode[0] = *tempCode;
      ++tempCode;
      specialCode[1] = *tempCode;
      specialCode[2] = 0;
    } // end if
  } // end for

  if ((nMobAl == mate1Set->size()) && (mate1Set->size() != 0)) { // all alignments are special
    *isMate1Special = true;
    sort(mate1Set->begin(), mate1Set->end(), Alignment_LessThanMq());
    *mate1SpecialAl = **(mate1Set->rbegin());
    mate1SpecialAl->SpecialCode = specialCode;
    mate1SpecialAl->NumMapped   = nMobAl;
    mate1Set->clear();

  } else if ( nMobAl > 0 ) { // some alignments are special
    vector<Alignment*> newMate1Set;
    for ( vector<Alignment*>::iterator ite = mate1Set->begin(); ite != mate1Set->end(); ++ite ) {
      if ( !(*ite)->IsMappedSpecialReference ) { 
	(*ite)->CanBeMappedToSpecialReference = true;
	(*ite)->SpecialCode = specialCode;
	newMate1Set.push_back( *ite );
      } else {
        *isMate1Special = true;
	if ((*ite)->Quality >= mate1SpecialAl->Quality)
	  *mate1SpecialAl = *(*ite);
	mate1SpecialAl->SpecialCode = specialCode;
	mate1SpecialAl->NumMapped   = nMobAl;
      } // end if-else
    } // end for
    mate1Set->clear();
    *mate1Set = newMate1Set;
  } // end if-else-if

}

// aligns the read against the reference sequence and returns true if the read was aligned
bool CAlignmentThread::AlignRead(CNaiveAlignmentSet& alignments, 
                                 const char* query, 
				 const char* qualities, 
				 const unsigned int& queryLength,
				 AlignmentStatusType& status,
				 int* numHash) {

	// set the alignment status to be INITIAL
	status = ALIGNMENTSTATUS_INITIAL;

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
	//int64_t hashRegionLength = 0;

	// control variables
	bool evaluateReverseReads = true;

	const bool alignAllReads    = mFlags.IsAligningAllReads;
	const bool useFastAlgorithm = (mAlgorithm == AlignerAlgorithm_FAST ? true : false);

	bool ret = true; // assume we will align the read
	//unsigned int numHash = 0;

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

			*numHash = forwardRegions.size() + reverseRegions.size();

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
				alignments.Add(al);
			}
		} else {
			for(unsigned int i = 0; i < (unsigned int)forwardRegions.size(); ++i) {
				
				// enforce alignment candidate thresholds
				// @Wan-Ping Lee: The following check is moved to GetReadCandidates function
				//if(mFlags.IsUsingAlignmentCandidateThreshold) {
				//	hashRegionLength = forwardRegions[i].End - forwardRegions[i].Begin + 1;
				//	if(hashRegionLength < mSettings.AlignmentCandidateThreshold) continue;
				//}

				// create a new alignment data structure
				Alignment al;
				//al.IsReverseStrand = false;

				// perform a Smith-Waterman alignment
				AlignRegion(forwardRegions[i], al, mForwardRead, queryLength, numExtensionBases);

				// add the alignment to the alignments vector
				if( ApplyReadFilters( al, query, qualities, queryLength ) ) {	
					alignments.Add(al);
				}

				// check if we can prematurely stop
				if( !alignAllReads && ( alignments.GetCount() > 1 ) ) {
					evaluateReverseReads = false;				
					break;
				}

				if ( mFlags.IsUsingHashRegionThreshold && ( i > mSettings.HashRegionThreshold ) )
					break;
			}

			// ==========================
			// evaluate the reverse reads
			// ==========================

			if(alignAllReads) evaluateReverseReads = true;

			if(evaluateReverseReads) {
				for(unsigned int i = 0; i < (unsigned int)reverseRegions.size(); ++i) {
					
					// enforce alignment candidate thresholds
					// @Wan-Ping Lee: The following check is moved to GetReadCandidates function 
					//if(mFlags.IsUsingAlignmentCandidateThreshold) {
					//	hashRegionLength = reverseRegions[i].End - reverseRegions[i].Begin + 1;
					//	if(hashRegionLength < mSettings.AlignmentCandidateThreshold) continue;
					//}

					// create a new alignment data structure
					Alignment al;
					al.IsReverseStrand = true;

					// perform a Smith-Waterman alignment
					AlignRegion(reverseRegions[i], al, mReverseRead, queryLength, numExtensionBases);

					// add the alignment to the alignments vector
					if( ApplyReadFilters(al, query, qualities, queryLength ) ) {
						alignments.Add(al);
					}

					// check if we can prematurely stop
					if (!alignAllReads && (alignments.GetCount() > 1))
						break;
					if ( mFlags.IsUsingHashRegionThreshold && ( i > mSettings.HashRegionThreshold ) )
						break;
				}
			}
		}

		// no alignments because of filtering
		if( alignments.IsEmpty() ) {
			ret = false;
			status = ALIGNMENTSTATUS_FILTEREDOUT;
		}

	} catch(bad_alloc &ba) {
		cout << "ERROR: Could not allocate enough memory to create forward and reverse aligned sequences: " << ba.what() << endl;
		exit(1);
	}

	return ret;
}

// aligns the read against a specified hash region using Smith-Waterman-Gotoh
void CAlignmentThread::AlignRegion(const HashRegion& r, Alignment& alignment, char* query, const unsigned int& queryLength, unsigned int& extensionBases) {

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
	while(r.Begin > mReferenceEnd[referenceIndex]) ++referenceIndex;


	const unsigned int refBegin = mReferenceBegin[referenceIndex];
	const unsigned int refEnd   = mReferenceEnd[referenceIndex];


	if(begin < refBegin) begin = refBegin;
	if(end   < refBegin) end   = refBegin;
	if(begin > refEnd)   begin = refEnd;
	if(end   > refEnd)   end   = refEnd;

	// adjust the begin and end positions if the reference is masked
	while(mReference[begin] == 'X') ++begin;
	while(mReference[end]   == 'X') --end;

	// perform a Smith-Waterman alignment on our region
	char* pAnchor = mReference + begin;

	
	if ( begin > mReferenceLength ) {
		cout << "ERROR: The hash region excceds the reference region." << endl;
		exit(1);
	}


	// determine if the specified bandwidth is enough to accurately align using the banded algorithm
	//bool hasEnoughBandwidth = false;
	//HashRegion diagonalRegion = r;

	/*
	if(mFlags.UseBandedSmithWaterman) {

		diagonalRegion.Begin -= begin;
		diagonalRegion.End   -= begin;

		unsigned int rowStart = min(diagonalRegion.Begin, (unsigned int)diagonalRegion.QueryBegin);

		diagonalRegion.Begin      -= rowStart;
		diagonalRegion.QueryBegin -= rowStart;

		hasEnoughBandwidth = (queryLength - diagonalRegion.QueryBegin) > mSettings.Bandwidth;
		hasEnoughBandwidth = hasEnoughBandwidth && (((end - begin + 1) - diagonalRegion.Begin) > mSettings.Bandwidth / 2);
	}
	*/

	//if(mFlags.UseBandedSmithWaterman && hasEnoughBandwidth) {
	//	mBSW.Align(alignment, pAnchor, (end - begin + 1), query, queryLength, diagonalRegion);
	//} else {
	//	mSW.Align(alignment, pAnchor, (end - begin + 1), query, queryLength);
	//}
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment ssw_alignment;
	mSSW.Align(query, pAnchor, (end - begin + 1), filter, &ssw_alignment);
	ConvertSswToAlignment(ssw_alignment, pAnchor, query, &alignment);

	// adjust the reference start positions
	//if ( !mFlags.UseLowMemory )
		alignment.ReferenceIndex = referenceIndex;
	//else
	//	alignment.ReferenceIndex = 0;
	alignment.ReferenceBegin += begin - refBegin;
	alignment.ReferenceEnd   += begin - refBegin;
}

// returns true if the alignment passes all of the user-specified filters
bool CAlignmentThread::ApplyReadFilters(Alignment& al, const char* bases, const char* qualities, const unsigned int& queryLength) {

	unsigned int queryLength1 = mFlags.EnableColorspace ? queryLength + 1 : queryLength;
	// assuming this is a good read
	bool ret = true;

	unsigned short numNonAlignedBases = queryLength1 - al.QueryLength;
	al.BaseQualities.Copy( qualities + al.QueryBegin, al.QueryEnd - al.QueryBegin + 1 );

	// convert from colorspace to basespace
	if( mFlags.EnableColorspace ) {
		if( al.IsReverseStrand ) al.BaseQualities.Reverse();
		if ( !mCS.ConvertAlignmentToBasespace( al ) ) ret = false;
		numNonAlignedBases = queryLength1 - al.QueryLength;
	} else {
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
void CAlignmentThread::CreateHash(const char* fragment, const unsigned char& fragmentLen, uint64_t& key) {

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
void CAlignmentThread::GetFastReadCandidate(HashRegion& region, char* query, const unsigned int& queryLength, MhpOccupancyList* pMhpOccupancyList) {

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
void CAlignmentThread::GetReadCandidates(vector<HashRegion>& regions, char* query, const unsigned int& queryLength, MhpOccupancyList* pMhpOccupancyList) {

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
	regions.reserve(hrt.GetCount());
	//vector<HashRegion>::iterator hrIter = regions.begin();

	hrt.GotoFirstEntry();
	while(HashRegion* r = hrt.GetTraversalHashRegion()) {
		if(!hrt.GetNextEntry()) break;
		//*hrIter = *r;
		//hrIter++;
		if(mFlags.IsUsingAlignmentCandidateThreshold) {
		  int length = r->End - r->Begin + 1;
		  if (length >= mSettings.AlignmentCandidateThreshold)
		    regions.push_back(*r);
		}
	}

	// sort the hash regions according to length (descending)
	if ((!mFlags.IsAligningAllReads)
	   || (mFlags.IsUsingHashRegionThreshold && (regions.size() > mSettings.HashRegionThreshold))) {
	  sort(regions.begin(), regions.end(), SortHashRegionByLength());
	}
	  
}

// settles the local Smith-Waterman window
bool CAlignmentThread::SettleLocalSearchRegion(
    const LocalAlignmentModel& lam, 
    const unsigned int& refIndex, 
    const unsigned int& uniqueBegin, 
    const unsigned int& uniqueEnd, 
    unsigned int& localSearchBegin, 
    unsigned int& localSearchEnd ) {

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
	while(mReference[begin] == 'X') ++begin;

	// adjust the stop position if the reference ends with a J nucleotide
	while(mReference[end] == 'X')   --end;


	// quit if we don't have a region to align against
	if (begin >= end) {
		localSearchBegin = 0;
		localSearchEnd   = 0;
		return false;
	} else if ((end - begin) > (mSettings.LocalAlignmentSearchRadius * 10)) {
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
bool CAlignmentThread::RescueMate(const LocalAlignmentModel& lam, 
                                  const CMosaikString& bases, 
				  const unsigned int& begin, 
				  const unsigned int& end, 
				  const unsigned int& refIndex, 
				  Alignment& al) {
	// prepare for alignment (borrow the forward read buffer)
	const char* query              = bases.CData();
	const unsigned int queryLength = bases.Length();

	strncpy_s(mForwardRead, mSettings.AllocatedReadLength, query, queryLength);
	mForwardRead[queryLength] = 0;

	if(lam.IsTargetReverseStrand) {
		if(mFlags.EnableColorspace) 
			CSequenceUtilities::ReverseSequence(mForwardRead, queryLength);
		else 
			CSequenceUtilities::GetReverseComplement(mForwardRead, queryLength);
	}

	// align according to the model
	al.IsReverseStrand = (lam.IsTargetReverseStrand ? true : false);
	char* pAnchor = mReference + begin;
	//mSW.Align(al, pAnchor, (end - begin + 1), mForwardRead, queryLength);
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment ssw_alignment;
	mSSW.Align(mForwardRead, pAnchor, (end - begin + 1), filter, &ssw_alignment);
	ConvertSswToAlignment(ssw_alignment, pAnchor, query, &al);

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
