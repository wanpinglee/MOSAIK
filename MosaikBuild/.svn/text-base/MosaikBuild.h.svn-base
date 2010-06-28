// ***************************************************************************
// BuildMain.cpp - imports the reads and reference sequences into MOSAIK.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#ifdef WIN32
#include <regex>
using namespace std::tr1;
#endif
#include <set>
#include <sstream>
#include "ColorspaceUtilities.h"
#include "FastLZIO.h"
#include "Mosaik.h"
#include "ConsoleUtilities.h"
#include "ConversionUtilities.h"
#include "Fasta.h"
#include "Fastq.h"
#include "MD5.h"
#include "MosaikString.h"
#include "ReadWriter.h"
#include "ProgressBar.h"
#include "ProgressCounter.h"
#include "ReadStatus.h"
#include "ReferenceSequence.h"
#include "ReferenceSequenceStatus.h"
#include "RegexUtilities.h"
#include "SequenceUtilities.h"
#include "SequencingTechnologies.h"
#include "SHA1.h"
#include "SRF.h"

using namespace std;

#define NUM_N_BASES_ALLOWED  4
#define MIN_READ_LENGTH     20

#define NORMAL_FASTQ_OFFSET 33
#define SOLEXA_FASTQ_OFFSET 64

#define NUM_REFERENCE_DIVIDER_BASES 500

class CMosaikBuild {
	friend class CMosaikBuildTests;
public:
	// constructor
	CMosaikBuild(const MosaikReadFormat::ReadGroup& md);
	// destructor
	~CMosaikBuild(void);
	// creates a read group ID from the current time
	static void CreateReadGroupID(string& readGroupID);
	// creates a MOSAIK reference archive
	void CreateReferenceArchive(const string& fastaFilename, const string& archiveFilename);
	// Enables the processing of base qualities
	void EnableBaseQualities(const string& filename);
	// Enables the processing of base qualities for the 2nd mate
	void EnableBaseQualities2(const string& filename);
	// Enables trimming of bases and qualities
	void EnableBaseTrimming(const unsigned short prefixTrim, const unsigned short suffixTrim);
	// Enables SOLiD colorspace translation
	void EnableColorspace(void);
	// Enables Helicos processing
	void EnableHelicosProcessing(void);
	// Enables instrument info removal
	void EnableInstrumentInfoRemoval(void);
	// Enables trimming the first bases from the read name
	void EnableReadNameTrimming(const unsigned char prefixTrim, const unsigned char suffixTrim);
	// Enables the addition of a user specified read name prefix
	void EnableReadNamePrefix(const string& prefix);
	// Enables a limit on the number of reads written to the read archive
	void EnableReadLimit(const uint64_t readLimit);
	// Parses an Illumina Bustard directory
	void ParseBustard(const string& directory, const string& lanes, const string& outputFilename, const bool splitReads);
	// Parses the sequence and quality FASTA files while writing to our read archive
	void ParseFasta(const string& readFastaFilename, const string& outputFilename);
	// Parses the sequence and quality paired-end FASTA files while writing to our read archive
	void ParsePEFasta(string& readFastaFilename, string& readFastaFilename2, const string& outputFilename);
	// Parses the reads and base qualities from a FASTQ file
	void ParseFastq(vector<string>& fastqFiles, const string& outputFilename);
	// Parses the reads and base qualities from a paired-end FASTQ file
	void ParsePEFastq(vector<string>& mate1Files, vector<string>& mate2Files, const string& outputFilename);
	// Parses an Illumina Gerald directory
	void ParseGerald(const string& directory, const string& lanes, const string& outputFilename);
	// Parses the SRF archive
	void ParseSRF(vector<string>& srfFiles, const string& outputFilename);
	// Sets the default base quality when a data set lacks BQ data
	void SetAssignedBaseQuality(unsigned char bq);
	// Sets the Genome Assembly ID [used when creating reference archives]
	void SetGenomeAssemblyID(const string& id);
	// Sets the maximum number of N's allowed
	void SetNumNBasesAllowed(const unsigned char numNBasesAllowed);
	// Sets the species name [used when creating reference archives]
	void SetSpecies(const string& name);
	// Sets the uniform resource identifier
	void SetURI(const string& uri);

private:
	// Stores settings pertinent to the building process
	struct BuildSettings {
		string BaseQualityFastaFilename;
		string BaseQualityFasta2Filename;
		string GenomeAssemblyID;
		string ReadFastaFilename;
		string SrfFilename;
		string Species;
		string UniformResourceIdentifier;
		unsigned char AssignedBaseQuality;
	} mSettings;
	// stores the colorspace triplet
	struct ColorspaceName {
		unsigned short first;
		unsigned short second;
		unsigned short third;

		// our less-than operator
		bool operator<(const ColorspaceName& cn) const {
			if((first == cn.first) && (second == cn.second)) return third < cn.third;
			if(first == cn.first) return second < cn.second;
			return first < cn.first;
		}
	};
	// stores the conversion statistics
	struct Statistics {
		uint64_t NumTotalMates;
		uint64_t NumReadsWritten;
		uint64_t NumBasesWritten;
		uint64_t NumMate1Orphaned;
		uint64_t NumMate2Orphaned;
		bool IsPairedEnd;

		Statistics(void)
			: NumTotalMates(0)
			, NumReadsWritten(0)
			, NumBasesWritten(0)
			, NumMate1Orphaned(0)
			, NumMate2Orphaned(0)
			, IsPairedEnd(false)
		{}
	};
	// stores the endpoints for each masked section
	struct MaskedPosition {
		unsigned int Begin;
		unsigned int End;

		MaskedPosition(unsigned int pos)
			: Begin(pos)
			, End(pos)
		{}
	};
	// our metadata object
	MosaikReadFormat::ReadGroup mReadGroup;
	// activates the specified Illumina lanes
	void ActivateIlluminaLanes(const string& lanes);
	// returns the colorspace name for the given read name
	static void GetColorspaceName(const CMosaikString& readName, ColorspaceName& cn);
	// trims the mate
	void ProcessMate(Mosaik::Mate& mate);
	// trims the read name and adds a read name prefix
	void ProcessReadName(CMosaikString& readName);
	// returns true if a swap had to performed to guarantee that filename1 is the F3 read (filename2 = R3)
	static void ReorderSolidFastaFilenames(string& filename1, string& filename2);
	// shows the conversion statistics
	void ShowStatistics(const Statistics& s);
	// toggles if SOLiD colorspace translation should be performed
	bool mEnableColorspace;
	// Denotes the presence of base qualities data
	bool mHasBaseQualities;
	bool mHasBaseQualities2;
	// toggles the availability of a read name prefix
	bool mHasReadNamePrefix;
	// toggles the use of the read limit
	bool mHasReadLimit;
	// toggles the trimming of reads
	bool mTrimReads;
	// toggles the trimming of read names
	bool mTrimReadNames;
	// toggles the removal of instrument info
	bool mRemoveInstrumentInfo;
	// toggles the trimming of reads with N's
	unsigned char mNumNBasesAllowed;
	unsigned int mNumLeadingNsTrimmed;
	unsigned int mNumLaggingNsTrimmed;
	unsigned int mNumMatesDeleted;
	unsigned short mMinimumReadLength;
	// our output buffer
	unsigned char* mBuffer;
	// the output buffer size
	unsigned int mBufferLen;
	// the number of prefix bases to trim from the read
	unsigned short mReadPrefixTrim;
	// the number of suffix bases to trim from the read
	unsigned short mReadSuffixTrim;
	// the number of prefix bases to trim from the read name
	unsigned char mReadNamePrefixTrim;
	unsigned char mReadNameSuffixTrim;
	// the maximum number of reads to process
	uint64_t mReadLimit;
	// the read name prefix
	CMosaikString mReadNamePrefix;
	// specifies the lanes that are allowed
	bool mAllowedLanes[8];
};
