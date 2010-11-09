// ***************************************************************************
// CMosaikText - exports alignments to various file formats.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <zlib.h>
#include "AlignedRead.h"
#include "AlignedReadCache.h"
#include "AlignmentReader.h"
#include "AlignmentWriter.h"
#include "ColorspaceUtilities.h"
#include "BamWriter.h"
#include "Fastq.h"
#include "FileUtilities.h"
#include "Mosaik.h"
#include "MosaikString.h"
#include "ProgressBar.h"
#include "Read.h"
#include "ReadReader.h"
#include "SequenceUtilities.h"

using namespace std;

#define CIGAR_BUFFER_SIZE 4096

// define some SAM/BAM FLAGS
#define BAM_SEQUENCED_AS_PAIRS         1
#define BAM_PROPER_PAIR                2
#define BAM_QUERY_UNMAPPED             4
#define BAM_MATE_UNMAPPED              8
#define BAM_QUERY_REVERSE_COMPLEMENT  16
#define BAM_MATE_REVERSE_COMPLEMENT   32
#define BAM_QUERY_FIRST_MATE          64
#define BAM_QUERY_SECOND_MATE        128
#define BAM_SECONDARY_ALIGNMENT      256
#define BAM_FAILS_PLATFORM_QUALITY   512
#define BAM_PCR_OPTICAL_DUPLICATE   1024

class CMosaikText {
public:
	// constructor
	CMosaikText(void);
	// destructor
	~CMosaikText(void);
	// enables AXT output
	void EnableAxtOutput(const string& filename);
	// enables BAM output
	void EnableBamOutput(const string& filename);
	// enables BED output
	void EnableBedOutput(const string& filename);
	// enables Eland output
	void EnableElandOutput(const string& filename);
	// enables FASTQ output
	void EnableFastqOutput(const string& filename);
	// enables Psl output
	void EnablePslOutput(const string& filename);
	// enables the reference sequence filter
	void EnableReferenceFilter(const string& referenceName, const string& alignmentFilename);
	// enables SAM output
	void EnableSamOutput(const string& filename);
	// enables screen output
	void EnableScreenOutput(void);
	// when triggered, the coverage calculation will only include unique reads
	void EvaluateUniqueReadsOnly(void);
	// parses the specified MOSAIK alignment file
	void ParseMosaikAlignmentFile(const string& alignmentFilename);
	// set the settings of input MOSAIK archive
	void SetArchiveSetting(const string& alignmentFilename);
	// parses the specified MOSAIK read file
	void ParseMosaikReadFile(const string& readFilename);
	// parse the fastq file
	void ParseFastqFile(const string& readFilename);
	// parse the 2nd mate fastq file
	void ParseFastq2File(const string& readFilename);
	// set sorting order
	void SetSortingOrder ( const unsigned short sortingModel );

private:
	struct PslBlock {
		unsigned int ReferenceStart;
		unsigned int QueryStart;
		unsigned int Length;
	};
	// our file stream structure
	struct Streams {
		FILE* axt;
		FILE* bed;
		FILE* eland;
		FILE* fastq;
		FILE* psl;
		gzFile sam;
		CBamWriter bam;
	} mStreams;
	// our flags data structure
	struct Flags {
		bool EvaluateUniqueReadsOnly;
		bool IsAxtEnabled;
		bool IsBamEnabled;
		bool IsBedEnabled;
		bool IsElandEnabled;
		bool IsFastqEnabled;
		bool IsPslEnabled;
		bool IsSamEnabled;
		bool IsScreenEnabled;
		bool UseReferenceFilter;
		bool EnableFastqPatching;
		bool IsSortingByPosition; // false: dose not change the order in the input archive

		Flags(void)
			: EvaluateUniqueReadsOnly(false)
			, IsAxtEnabled(false)
			, IsBamEnabled(false)
			, IsBedEnabled(false)
			, IsElandEnabled(false)
			, IsFastqEnabled(false)
			, IsPslEnabled(false)
			, IsSamEnabled(false)
			, IsScreenEnabled(false)
			, UseReferenceFilter(false)
			, EnableFastqPatching(false)
			, IsSortingByPosition(true)
		{}
	} mFlags;
	// our settings data structure
	struct Settings {
		string AxtFilename;
		string BamFilename;
		string BedFilename;
		string ElandFilename;
		string FastqFilename;
		string PslFilename;
		string SamFilename;
		string inputFastqFilename;
		string inputFastq2Filename;
		string SortingModel;
		unsigned int FilteredReferenceIndex;
		uint64_t NumFilteredReferenceReads;

	} mSettings;
	// settings of input MOSAIK archive
	struct ArchiveSetting {
		vector<MosaikReadFormat::ReadGroup> readGroups;
		vector<ReferenceSequence> pReferenceSequences;
		AlignmentStatus as;
		char* signature;

		ArchiveSetting(void)
			:signature(NULL)
		{}
	} mArchiveSetting;

	
	// opens the output file stream for the AXT file
	void InitializeAxt(void);
	// opens the output file stream for the BAM file
	void InitializeBam(vector<ReferenceSequence>* pRefSeqs, vector<MosaikReadFormat::ReadGroup>& readGroups);
	// opens the output file stream for the BED file
	void InitializeBed(void);
	// opens the output file stream for the Eland file
	void InitializeEland(void);
	// opens the output file stream for the SAM file
	void InitializeSam(vector<ReferenceSequence>* pRefSeqs, vector<MosaikReadFormat::ReadGroup>& readGroups);
	// processes the alignments according to the chosen file format
	void ProcessAlignments(const unsigned char mateNum, const CMosaikString& readName, const bool isColorspace, vector<Alignment>& alignments, const string& readGroupID);
	// processes the mates according to the chosen file format
	void ProcessMate(const unsigned char mateNum, const CMosaikString& readName, const bool isColorspace, Mosaik::Mate& mate, const bool isPairedEnd);
	// writes the current alignment to the SAM output file
	void WriteSamEntry(const CMosaikString& readName, const string& readGroupID, const vector<Alignment>::iterator& alIter);
	// patchs trimmed infomation back from FASTQs
	void PatchInfo( const string& alignmentFilename, const string& inputFastqFilename, const string& inputFastq2Filename );
	// given an alignedReadCache, sort them by positions and sorte them in a temp file
	void StoreReadCache ( CAlignedReadCache& cache );
	// given a read name, search it in FASTQs
	void SearchReadInFastq ( const CMosaikString& readName, CFastq& fastqReader1, CFastq& fastqReader2, const bool hasFastq2 );
	// initialize our patching buffers
	void InitializePatchingBuffer ( void );
	// free patching buffers
	void FreePatchingBuffer ( void );
	// sort FASTQ by read names
	void SortFastqByName( const string& inputFastqFilename, string& outputFastqFilename );
	// merge a vector of partially sorted archive; sorting them by positions
	void MergeSortedArchive ( const vector <string>& filenames, const string& outArchiveFilename );
	// sort aligned archive by positions
	void SortAlignmentByPosition( const string& inputArchive, const string& outputArchive );
	// cigar buffer
	char mCigarBuffer[CIGAR_BUFFER_SIZE];
	// our current read and alignment counters
	uint64_t mCurrentAlignment;
	uint64_t mCurrentRead;
	// our buffer for patching function
	unsigned int _bufferSize1;
	unsigned int _bufferSize2;
	unsigned int _clipSize;
	char* _originalReverseBase1;
	char* _originalReverseBase2;
	char* _originalReverseQuality1;
	char* _originalReverseQuality2;
	char* _clipBuffer;
	CMosaikString  _readName1, _readName2;
	Mosaik::Mate   _m1, _m2;
	vector<string> _tempFiles;
	// our colorspace to basespace converter
	CColorspaceUtilities mCS;
	// the lessthan operator used in MosaikText
	static inline bool PositionLessThan( const Mosaik::AlignedRead& ar1, const Mosaik::AlignedRead& ar );
};
