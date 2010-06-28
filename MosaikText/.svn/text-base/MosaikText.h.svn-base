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
#include "AlignmentReader.h"
#include "ColorspaceUtilities.h"
#include "BamWriter.h"
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
	// parses the specified MOSAIK read file
	void ParseMosaikReadFile(const string& readFilename);

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
		unsigned int FilteredReferenceIndex;
		uint64_t NumFilteredReferenceReads;
	} mSettings;
	// opens the output file stream for the AXT file
	void InitializeAxt(void);
	// opens the output file stream for the BAM file
	void InitializeBam(const bool isSortedByPosition, vector<ReferenceSequence>* pRefSeqs, vector<MosaikReadFormat::ReadGroup>& readGroups);
	// opens the output file stream for the BED file
	void InitializeBed(void);
	// opens the output file stream for the Eland file
	void InitializeEland(void);
	// opens the output file stream for the SAM file
	void InitializeSam(const bool isSortedByPosition, vector<ReferenceSequence>* pRefSeqs, vector<MosaikReadFormat::ReadGroup>& readGroups);
	// processes the alignments according to the chosen file format
	void ProcessAlignments(const unsigned char mateNum, const CMosaikString& readName, const bool isColorspace, vector<Alignment>& alignments, const string& readGroupID);
	// processes the mates according to the chosen file format
	void ProcessMate(const unsigned char mateNum, const CMosaikString& readName, const bool isColorspace, Mosaik::Mate& mate, const bool isPairedEnd);
	// writes the current alignment to the SAM output file
	void WriteSamEntry(const CMosaikString& readName, const string& readGroupID, const vector<Alignment>::iterator& alIter);
	// cigar buffer
	char mCigarBuffer[CIGAR_BUFFER_SIZE];
	// our current read and alignment counters
	uint64_t mCurrentAlignment;
	uint64_t mCurrentRead;
	// our colorspace to basespace converter
	CColorspaceUtilities mCS;
};
