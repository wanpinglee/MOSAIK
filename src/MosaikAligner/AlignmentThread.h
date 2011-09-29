// ***************************************************************************
// CAlignmentThread - aligns all of the reads within a worker thread.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg & Wan-Ping Lee
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <limits.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <string>
#include <vector>
#include "AbstractDnaHash.h"
#include "AlignmentWriter.h"
#include "BandedSmithWaterman.h"
#include "BamWriter.h"
#include "BestNSecondBestSelection.h"
#include "CigarTager.h"
#include "ColorspaceUtilities.h"
#include "Entropy.h"
#include "MdTager.h"
#include "NaiveAlignmentSet.h"
#include "PairwiseUtilities.h"
#include "PosixThreads.h"
#include "QualityNeuralNetwork.h"
#include "ReadReader.h"
#include "ReferenceSequence.h"
#include "SequenceUtilities.h"
#include "SmithWatermanGotoh.h"
#include "StatisticsMaps.h"
#include "ZaTager.h"

using namespace std;

#define ALLOCATION_EXTENSION 10

class CAlignmentThread {
public:
	// our enumerated alignment algorithms
	enum AlignerAlgorithmType {
		AlignerAlgorithm_FAST,
		AlignerAlgorithm_SINGLE,
		AlignerAlgorithm_MULTI,
		AlignerAlgorithm_ALL
	};
	// our enumerated alignment modes
	enum AlignerModeType {
		AlignerMode_UNIQUE,
		AlignerMode_ALL
	};
	// define our alignment configuration structure
	struct AlignerSettings {
		string AlignedReadReportFilename;
		string BasespaceReferenceFilename;
		string InputReadArchiveFilename;
		string JumpFilenameStub;
		string OutputReadArchiveFilename;
		string ReferenceFilename;
		//string UnalignedReadReportFilename;
		unsigned int AllocatedReadLength;
		unsigned int Bandwidth;
		unsigned int LocalAlignmentSearchRadius;
		unsigned int MedianFragmentLength;
		unsigned int NumCachedHashes;
		unsigned short AlignmentCandidateThreshold;
		unsigned short HashPositionThreshold;
		unsigned short HashRegionThreshold;
		unsigned char HashSize;
		unsigned char NumThreads;
		SequencingTechnologies SequencingTechnology;
	};
	// stores the filter settings
	struct FilterSettings {
		bool UseMinAlignmentFilter;
		bool UseMinAlignmentPercentFilter;
		bool UseMismatchFilter;
		bool UseMismatchPercentFilter;

		double MinPercentAlignment;
		double MaxMismatchPercent;
		unsigned int MinAlignment;
		unsigned int MaxNumMismatches;
		unsigned char LocalAlignmentSearchHighMqThreshold;
		unsigned char LocalAlignmentSearchLowMqThreshold;

		FilterSettings()
			: UseMinAlignmentFilter(CPairwiseUtilities::UseMinAlignmentFilter)
			, UseMinAlignmentPercentFilter(CPairwiseUtilities::UseMinAlignmentPercentFilter)
			, UseMismatchFilter(CPairwiseUtilities::UseMismatchFilter)
			, UseMismatchPercentFilter(CPairwiseUtilities::UseMismatchPercentFilter)
			, MinPercentAlignment(CPairwiseUtilities::MinPercentAlignment)
			, MaxMismatchPercent(CPairwiseUtilities::MaxMismatchPercent)
			, MinAlignment(CPairwiseUtilities::MinAlignment)
			, MaxNumMismatches(CPairwiseUtilities::MaxNumMismatches)
			, LocalAlignmentSearchHighMqThreshold(30)
			, LocalAlignmentSearchLowMqThreshold(10)
		{}
	};
	// define our boolean flags structure
	struct FlagData {
		bool EnableColorspace;
		bool IsAligningAllReads;
		bool IsQuietMode;
		//bool IsReportingUnalignedReads;
		bool IsUsingAlignmentCandidateThreshold;
		bool IsUsingHashPositionThreshold;
		bool IsUsingHashRegionThreshold;
		bool IsUsingJumpDB;
		bool KeepJumpKeysInMemory;
		bool KeepJumpPositionsInMemory;
		bool OutputMultiply;
		bool UseAlignedReadLengthForMismatchCalculation;
		bool UseBandedSmithWaterman;
		bool UseLocalAlignmentSearch;
		bool UsePairedEndOutput;
		bool UseLowMemory;
		bool UseBamOutput;
		bool UseArchiveOutput;
		bool SaveMultiplyBam;
		bool SaveUnmappedBasesInArchive;

		FlagData()
			: EnableColorspace(false)
			, IsAligningAllReads(false)
			, IsQuietMode(false)
			//, IsReportingUnalignedReads(false)
			, IsUsingAlignmentCandidateThreshold(false)
			, IsUsingHashPositionThreshold(false)
			, IsUsingHashRegionThreshold(false)
			, IsUsingJumpDB(false)
			, KeepJumpKeysInMemory(false)
			, KeepJumpPositionsInMemory(false)
			, OutputMultiply(false)
			, UseAlignedReadLengthForMismatchCalculation(false)
			, UseBandedSmithWaterman(false)
			, UseLocalAlignmentSearch(false)
			, UsePairedEndOutput(false)
			, UseLowMemory(false)
			, UseBamOutput(false)
			, UseArchiveOutput(false)
			, SaveMultiplyBam(false)
			, SaveUnmappedBasesInArchive(false)
		{}
	};
	// stores the statistical counters
	struct StatisticsCounters {
		// single-end
		uint64_t FailedHashMates;
		uint64_t FailedHashMates_Rescue;
		uint64_t FilteredOutMates;   // try to aligned but passing filters fails.
		uint64_t FilteredOutMates_Rescue;
		uint64_t ShortMates;
		uint64_t ShortMates_Rescue;
		uint64_t TooManyNsMates;
		uint64_t TooManyNsMates_Rescue;
		uint64_t MultipleMates;
		uint64_t MultipleMates_Rescue;
		uint64_t UniqueMates;
		uint64_t UniqueMates_Rescue;
		uint64_t Unmapped;
		// paired-end
		uint64_t UU;  // uniquely-uniquely aligned pairs
		uint64_t UM;  // uniquely-multiply aligned pairs
		uint64_t UF;  // uniquely-filtered aligned pairs; in which filtered out mates are mapped but cannot pass the mismatch filters
		uint64_t MM;  // multiply-multiply aligned pairs
		uint64_t MF;  // multiply-filtered aligned pairs; in which filtered out mates are mapped but cannot pass the mismatch filters
		uint64_t UX;  // uniquely-unmapped aligned pairs
		uint64_t MX;  // multiply-unmapped aligned pairs
		uint64_t FF;  // filtered-filtered aligned pairs
		uint64_t FX;  // filtered-unmapped aligned pairs
		uint64_t XX;  // unmapped-unmapped aligned pairs
		uint64_t UU_localRescue;
		uint64_t UU_localConsistance;
		uint64_t UM_localRescue;
		uint64_t UM_localConsistance;
		uint64_t MM_localRescue;
		uint64_t MM_localConsistance;
		unsigned char StatMappingQuality;

		StatisticsCounters() 
			: FailedHashMates(0)
			, FailedHashMates_Rescue(0)
			, FilteredOutMates(0)
			, FilteredOutMates_Rescue(0)
			, ShortMates(0)
			, ShortMates_Rescue(0)
			, TooManyNsMates(0)
			, TooManyNsMates_Rescue(0)
			, MultipleMates(0)
			, MultipleMates_Rescue(0)
			, UniqueMates(0)
			, UniqueMates_Rescue(0)
			, Unmapped(0)
			, UU(0)
			, UM(0)
			, UF(0)
			, MM(0)
			, MF(0)
			, UX(0)
			, MX(0)
			, FF(0)
			, FX(0)
			, XX(0)
			, UU_localRescue(0)
			, UU_localConsistance(0)
			, UM_localRescue(0)
			, UM_localConsistance(0)
			, MM_localRescue(0)
			, MM_localConsistance(0)
			, StatMappingQuality(20)
		{}
	};
	// bam writers
	struct BamWriters {
		BamHeader mHeader; // multiply alignments
		BamHeader sHeader; // special alignments
		//BamHeader uHeader; // unaligned reads
		BamHeader rHeader; // regular bam

		CBamWriter mBam;   // multiply alignments
		CBamWriter sBam;   // special alignments
		//CBamWriter uBam;   // unaligned reads
		CBamWriter rBam;
	};
	// special reference 
	struct SReference {
		bool     enable;
		bool     found;
		uint64_t begin;
		uint64_t nReference;
		double   count;
		string   prefix;

		SReference()
			: enable(false)
			, found(false)
			, begin(0)
			, nReference(0)
			, count(0)
		{}
	};
	// constructor
	CAlignmentThread(
		const AlignerAlgorithmType& algorithmType, 
		const FilterSettings&       filters, 
		const FlagData&             flags, 
		const AlignerModeType&      algorithmMode, 
		char*                       pReference, 
		const unsigned int          referenceLen, 
		CAbstractDnaHash*           pDnaHash, 
		const AlignerSettings&      settings, 
		unsigned int*               pRefBegin, 
		unsigned int*               pRefEnd, 
		char**                      pRefSpecies, 
		bool*                       pRefSpecial,
		char**                      pBsRefSeqs, 
		const SReference&           SpecialReference,
		map <unsigned int, MosaikReadFormat::ReadGroup>* pReadGroupsMap,
		const unsigned int          referenceOffset,
		string                      i_paired_end_ann_file,
		string                      i_single_end_ann_file
	);

	// destructor
	~CAlignmentThread(void);
	// define our thread data structure
	struct ThreadData {
		AlignerAlgorithmType Algorithm;
		AlignerModeType      Mode;
		AlignerSettings      Settings;
		FilterSettings       Filters;
		FlagData             Flags;
		StatisticsCounters*  pCounters;
		CStatisticsMaps*     pMaps;
		CAbstractDnaHash*    pDnaHash;
		MosaikReadFormat::CReadReader*      pIn;
		MosaikReadFormat::CAlignmentWriter* pOut;
		unsigned int  ReferenceLen;
		char*         pReference;
		unsigned int* pRefBegin;
		unsigned int* pRefEnd;
		char**        pRefSpecies;
		bool*         pRefSpecial;
		uint64_t*     pReadCounter;
		bool          IsPairedEnd;
		char**        pBsRefSeqs;
		BamWriters*   pBams;
		SReference    SpecialReference;
		map< unsigned int, MosaikReadFormat::ReadGroup >* pReadGroups;
		unsigned int  ReferenceOffset;
		string        paired_end_ann_file;
		string        single_end_ann_file;
	};
	// aligns the read archive
	void AlignReadArchive(
		MosaikReadFormat::CReadReader*      pIn, 
		MosaikReadFormat::CAlignmentWriter* pOut, 
		//FILE*            pUnalignedStream,
		uint64_t*        pReadCounter,
		bool             isPairedEnd, 
		CStatisticsMaps* pMaps, 
		BamWriters*      pBams,
		unsigned char    statMappingQuality
	);
	// activates the current alignment thread
	static void* StartThread(void* arg);
	// register our thread mutexes
	static pthread_mutex_t mGetReadMutex;
	static pthread_mutex_t mSaveReadMutex;
	static pthread_mutex_t mStatisticsMutex;
	static pthread_mutex_t mStatisticsMapsMutex;
	static pthread_mutex_t mSaveMultipleBamMutex;
	static pthread_mutex_t mSaveUnmappedBamMutex;
	static pthread_mutex_t mSaveSpecialBamMutex;
	// stores the statistical counters
	StatisticsCounters mStatisticsCounters;
private:

	// ===========
	// data struct
	// ===========
	enum AlignmentStatusType { 
		ALIGNMENTSTATUS_FAILEDHASH,
		ALIGNMENTSTATUS_TOOSHORT, 
		ALIGNMENTSTATUS_TOOMANYNS, 
		ALIGNMENTSTATUS_FILTEREDOUT,
		ALIGNMENTSTATUS_INITIAL };
	// define a comparison function for sorting our hash regions
	struct SortHashRegionByLength {
		bool operator()(const HashRegion& hr1, const HashRegion& hr2) {
			return (hr2.End - hr2.Begin) < (hr1.End - hr1.Begin);
		}
	};
	// our local alignment model data structure used in mate rescue
	struct LocalAlignmentModel {
		bool IsTargetBeforeUniqueMate;
		bool IsTargetReverseStrand;

		LocalAlignmentModel(void)
			: IsTargetBeforeUniqueMate(false)
			, IsTargetReverseStrand(false)
		{}
	};
	// data structure for alignment buffer
	struct AlignmentBamBuffer {
		Alignment al;
		string    zaString;
		bool      noCigarMdNm;
		bool      notShowRnamePos;

		AlignmentBamBuffer( void )
			: noCigarMdNm(false)
			, notShowRnamePos(false)
		{}
	};
	struct AlignmentArchiveBuffer {
		Mosaik::Read mr;
		Alignment al1;
		Alignment al2;
		bool isLongRead;

		AlignmentArchiveBuffer( void )
			: isLongRead(false)
		{}
	};
	// data structure for multiply alignment buffer
	struct SimpleBamRecordBuffer {
		unsigned short refIndex;
		unsigned int   refBegin;
		unsigned int   refEnd;
	};
	// common info of the alignment
	struct AlignmentInfo {
		bool isUsing454;
		bool isUsingIllumina;
		bool isUsingSOLiD;
		bool isUsingIlluminaLong;
		bool isPairedEnd;
		bool isUsingLowMemory;
	};

	// =========
	// functions
	// =========

	// aligns the read against the reference sequence and returns true if the read was aligned
	bool AlignRead(CNaiveAlignmentSet& alignments, const char* query, const char* qualities, const unsigned int queryLength, AlignmentStatusType& status);
	// aligns the read against a specified hash region using Smith-Waterman-Gotoh
	void AlignRegion(const HashRegion& r, Alignment& alignment, char* query, unsigned int queryLength, unsigned int extensionBases);
	// returns true if the alignment passes all of the user-specified filters
	bool ApplyReadFilters(Alignment& al, const char* bases, const char* qualities, const unsigned int queryLength);
	// creates the hash for a supplied fragment
	void CreateHash(const char* fragment, const unsigned char fragmentLen, uint64_t& key);
	// consolidates hash hits into a read candidate (fast algorithm)
	void GetFastReadCandidate(HashRegion& region, char* query, const unsigned int queryLength, MhpOccupancyList* pMhpOccupancyList);
	// consolidates hash hits into read candidates
	void GetReadCandidates(vector<HashRegion>& regions, char* query, const unsigned int queryLength, MhpOccupancyList* pMhpOccupancyList);
	// settles the local Smith-Waterman window
	bool SettleLocalSearchRegion( const LocalAlignmentModel& lam, const unsigned int refIndex
		, const unsigned int uniqueBegin, const unsigned int uniqueEnd, unsigned int& localSearchBegin, unsigned int& localSearchEnd );
	// attempts to rescue the mate paired with a unique mate
	bool RescueMate(const LocalAlignmentModel& lam, const CMosaikString& bases, const unsigned int& begin
		, const unsigned int& end, const unsigned int& refIndex, Alignment& al);
	// Prepare bam required info
	void SetRequiredInfo ( Alignment& al, const AlignmentStatusType& status, Alignment& mate, const Mosaik::Mate& m, const Mosaik::Read& r
		, const bool& isPair, const bool& isProperPair, const bool& isFirstMate, const bool& isPairTech, const bool& isItselfMapped, const bool& isMateMapped);
	void ProcessSpecialAlignment (vector<Alignment*>* mate1Set 
		, Alignment* mate1SpecialAl, bool* mate1Special);
	// treat the best alignment as an unique mapping and than turn on local search
	bool TreatBestAsUnique(vector<Alignment*>* mateSet, const unsigned int readLength);
	// update statistics
	void UpdateStatistics ( const enum AlignmentStatusType& mate1Status, const enum AlignmentStatusType& mate2Status
		, const Alignment al1, const Alignment al2, const bool isProperPair );
	void UpdateSeStatistics ( const enum AlignmentStatusType& mateStatus, const Alignment al );
	void SearchLocalRegion(const vector<Alignment*>& anchorVector, const Mosaik::Mate& mate, CNaiveAlignmentSet* mateVector);
	inline void SaveBamAlignment( const Alignment& al, const char* zaString, const bool noCigarMdNm, const bool notShowRnamePos, const bool isSpecial );
	inline void SaveArchiveAlignment ( const Mosaik::Read& mr, const Alignment& al1, const Alignment& al2, const bool isLongRead );
	void WriteAlignmentBufferToFile( BamWriters* const pBams, CStatisticsMaps* const pMaps, MosaikReadFormat::CAlignmentWriter* const pOut );
	void WriteSpecialAlignmentBufferToFile( BamWriters* const pBams );
	void SaveMultiplyAlignment(const vector<Alignment*>& mate1Set, const vector<Alignment*>& mate2Set, const Mosaik::Read& mr
		, BamWriters* const pBams, CStatisticsMaps* const pMaps);
	void SaveNClearBuffers( BamWriters* const pBams, CStatisticsMaps* const pMaps, MosaikReadFormat::CAlignmentWriter* const pOut );
	unsigned char GetMappingQuality(const Alignment& al, const int& al_length);
	unsigned char GetMappingQuality(
	    const Alignment& al1, 
	    const int& al1_length, 
	    const Alignment& al2, 
	    const int& al2_length);

	// ====
	// data
	// ====
	
	// denotes the active alignment algorithm
	AlignerAlgorithmType mAlgorithm;
	// denotes the active alignment mode
	AlignerModeType mMode;
	// stores the alignment configuration
	AlignerSettings mSettings;
	// stores the filter configuration
	FilterSettings mFilters;
	// stores our boolean flags
	FlagData mFlags;
	// the reference sequence
	char* mReference;
	// the sepcial references
	SReference mSReference;
	// our forward and reverse complement copy of the read
	char* mForwardRead;
	char* mReverseRead;
	// the length of the reference sequence
	unsigned int mReferenceLength;
	// the hash-table associated with the specified alignment algorithm
	CAbstractDnaHash* mpDNAHash;
	// our Smith-Waterman-Gotoh local alignment algorithms
	CSmithWatermanGotoh mSW;
	CBandedSmithWaterman mBSW;
	// our reference sequence LUTs
	unsigned int* mReferenceBegin;
	unsigned int* mReferenceEnd;
	char** mReferenceSpecies;
	bool*  mReferenceSpecial;
	// for soft clips
	//const unsigned int softClippedIdentifierLength;
	char* softClippedIdentifier;
	// our colorspace to basespace converter
	CColorspaceUtilities mCS;
	vector<ReferenceSequence> mpBsRefSeqs;
	// read groups map
	map<unsigned int, MosaikReadFormat::ReadGroup>* mReadGroupsMap;
	// reference offset used for low-memory multiply-mapped bam
	unsigned int mReferenceOffset;
	// ZA tagers
	CZaTager za1, za2;
	// MD tager
	CMdTager mdTager;
	CBamWriter bamMisc;
	unsigned int _bufferSize;
	// alignment buffer used for speedup
	queue<AlignmentBamBuffer> bamBuffer;         // for bam output; full-memory version
	queue<AlignmentBamBuffer> bamMultiplyBuffer;
	queue<SimpleBamRecordBuffer> bamMultiplySimpleBuffer;
	queue<AlignmentArchiveBuffer> archiveBuffer; // for archive output; low-memory version
	queue<AlignmentBamBuffer> bamSpecialBuffer;
	// reads/mates buffer
	queue<Mosaik::Read> readBuffer;
	AlignmentInfo alInfo;

	// neural-net
	QualityNeuralNetwork mqCalculator;
	QualityNeuralNetwork::FannInputs mate1Ann;
	QualityNeuralNetwork::FannInputs mate2Ann;
	string paired_end_ann_file;
	string single_end_ann_file;
	
	//Entropy
	Entropy entropy_;
};
