// ***************************************************************************
// CPairedEndSort - resolves paired-end reads and creates an alignment archive
//                  sorted by reference sequence position.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iostream>
#include <list>
#include <set>
#include <map>
#include <string>
#include "sqlite3.h"
#include "AlignmentReader.h"
#include "AlignmentStatus.h"
#include "AlignmentWriter.h"
#include "AlignedRead.h"
#include "ConsoleUtilities.h"
#include "ReadGroup.h"
#include "Mosaik.h"
#include "ProgressBar.h"
#include "ProgressCounter.h"
#include "SequencingTechnologies.h"
#include "UnorderedMap.h"

using namespace std;

// here we assume that we'll need space for Sanger length reads (reference bases,
// query bases, query base qualities)
static double DEFAULT_CONFIDENCE_INTERVAL = 0.9973;
#define GAP '-'

// define our serialization status flags
#define PE_IS_LONG_READ           1
#define PE_IS_REVERSE_STRAND      2
#define PE_IS_MATE_REVERSE_STRAND 4
#define PE_IS_FIRST_MATE          8
#define PE_IS_RESOLVED_AS_PAIR    16
#define PE_WAS_RESCUED            32

#define SQL_BUFFER_SIZE           10240
#define DUMMY_MODEL               100
#define MODEL_COUNT_THRESHOLD     0.1

#define UU_COEFFICIENT  0.486796848
#define UU_INTERCEPT   35.45967112
#define UM_COEFFICIENT  0.426395518
#define UM_INTERCEPT   19.29236958
#define MM_COEFFICIENT  0.327358673
#define MM_INTERCEPT    4.350331532

// define our model count data structure
struct ModelType {
	unsigned char ID;
	unsigned int Count;

	// constructor
	ModelType(void)
		: ID(0)
		, Count(0)
	{}

	// our less-than operator
	bool operator<(const ModelType& mt) const {
		return mt.Count < Count;
	}
};

class CPairedEndSort {
public:
	// constructor
	CPairedEndSort(const unsigned int numCachedReads);
	// destructor
	~CPairedEndSort(void);
	// configures which read pair types should be resolved
	void ConfigureResolution(const bool uo, const bool uu, const bool um, const bool mm);
	// disables fragment alignment quality calculation
	void DisableFragmentAlignmentQuality(void);
	// allows any fragment length when evaluating unique mate-pairs
	void EnableAllUniqueFragmentLengths(void);
	// enables consed renaming
	void EnableConsedRenaming(void);
	// enables duplicate read filtering
	void EnableDuplicateFiltering(const string& duplicateDirectory);
	// enables the sampling of all read pairs
	void EnableFullFragmentLengthSampling(void);
	// resolves the paired-end reads found in the specified input file
	void ResolvePairedEndReads(const string& inputFilename, const string& outputFilename);
	// sets the desired confidence interval
	void SetConfidenceInterval(const double& percent);
	// patch the original fastq information
	void PatchFastq(void);

private:
	// define our sort configuration structure
	struct SortSettings {
		unsigned char AlignmentModel1;
		unsigned char AlignmentModel2;
		string DuplicateDirectory;
		string UnresolvedFilename;
		double ConfidenceInterval;
		unsigned int NumCachedReads;

		SortSettings() 
			: ConfidenceInterval(DEFAULT_CONFIDENCE_INTERVAL)
			, NumCachedReads(0)
		{}
	} mSettings;
	// define our boolean flags structure
	struct FlagData {
		bool AllowAllUniqueFragmentLengths;
		bool RemoveDuplicates;
		bool RenameMates;
		bool ResolveMM;
		bool ResolveUM;
		bool ResolveUO;
		bool ResolveUU;
		bool SampleAllFragmentLengths;
		bool UseFragmentAlignmentQuality;
		bool PatchFastq;
		//bool SortByName;

		FlagData()
			: AllowAllUniqueFragmentLengths(false)
			, RemoveDuplicates(false)
			, RenameMates(false)
			, ResolveMM(false)
			, ResolveUM(false)
			, ResolveUO(false)
			, ResolveUU(false)
			, SampleAllFragmentLengths(false)
			, UseFragmentAlignmentQuality(true)
			//, PatchFastq(false)
			//, SortByName(false)
		{}
	} mFlags;
	// retrieves an alignment from the specified temporary file and adds it to the specified list
	void AddAlignment(FILE* tempFile, const unsigned int owner, list<Alignment>& alignments);
	// returns the current alignment model based on the order and orientation of the mates
	static inline unsigned char GetCurrentModel(unsigned int m1Begin, bool m1IsReverseStrand, unsigned int m2Begin, bool m2IsReverseStrand);
	// calculates the fragment alignment quality based on the fragment class
	static inline unsigned char GetFragmentAlignmentQuality(int aq, const bool isUU, const bool isMM);
	static inline unsigned char GetFragmentAlignmentQuality(int aq, const bool isUU, const bool isMM, unsigned int fragmentLen);
	// retrieves an alignment from the specified temporary file
	bool GetAlignment(FILE* tempFile, const unsigned int owner, Alignment& al);
	// records the observed gaps in the specified reference 
	void RecordReferenceGaps(Alignment& al);
	// serializes the specified vector to a temporary file
	uint64_t Serialize(list<Alignment>& alignmentCache, const unsigned int numEntries);
	// our temporary file vector
	vector<string> mTempFiles;
	// our output buffer
	unsigned char* mBuffer;
	unsigned int mBufferLen;
	// our reference gap hash map vector and associated iterator
	vector<unordered_map<unsigned int, unsigned short> > mRefGapVector;
	unordered_map<unsigned int, unsigned short>::iterator mRefGapIter;
	// sort alignments by their names
	static inline bool NameLessThan(const Alignment& al1, const Alignment& al2);
};

// returns the current alignment model based on the order and orientation of the mates
inline unsigned char CPairedEndSort::GetCurrentModel(unsigned int m1Begin, bool m1IsReverseStrand, unsigned int m2Begin, bool m2IsReverseStrand) {

	unsigned char currentModel = DUMMY_MODEL;

	if(m1Begin < m2Begin) {

		if(!m1IsReverseStrand && !m2IsReverseStrand) currentModel = 0;
		if(!m1IsReverseStrand && m2IsReverseStrand)  currentModel = 1;
		if(m1IsReverseStrand  && !m2IsReverseStrand) currentModel = 2;
		if(m1IsReverseStrand  && m2IsReverseStrand)  currentModel = 3;

	} else {

		if(!m2IsReverseStrand && !m1IsReverseStrand) currentModel = 4;
		if(!m2IsReverseStrand && m1IsReverseStrand)  currentModel = 5;
		if(m2IsReverseStrand  && !m1IsReverseStrand) currentModel = 6;
		if(m2IsReverseStrand  && m1IsReverseStrand)  currentModel = 7;
	}

	return currentModel;
}

// calculates the fragment alignment quality based on the fragment class
inline unsigned char CPairedEndSort::GetFragmentAlignmentQuality(int aq, const bool isUU, const bool isMM) {

	if(isUU)      aq = (int)(UU_COEFFICIENT * aq + UU_INTERCEPT);
	else if(isMM) aq = (int)(MM_COEFFICIENT * aq + MM_INTERCEPT);
	else          aq = (int)(UM_COEFFICIENT * aq + UM_INTERCEPT);

	if(aq < 0)  aq = 0;
	if(aq > 99) aq = 99;
	return (unsigned char)aq;
}

// sort alignments by their names
inline bool CPairedEndSort::NameLessThan(const Alignment& al1, const Alignment& al2){
	return al1.Name < al2.Name;
}
