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
//

#pragma once


#include <iostream>
#include <fstream>
#include <cmath>
#include <utility>
#include "AlignmentThread.h"
#include "AlignmentWriter.h"
#include "AlignmentReader.h"
#include "ArchiveMerge.h"
#include "BamWriter.h"
#include "Benchmark.h"
#include "ConsoleUtilities.h"
#include "DnaHash.h"
#include "JumpDnaHash.h"
#include "ReadReader.h"
#include "MultiDnaHash.h"
#include "PosixThreads.h"
#include "ProgressBar.h"
#include "ReferenceSequence.h"
#include "ReferenceSequenceReader.h"
#include "UbiqDnaHash.h"
#include "FileUtilities.h"
#include "SortThread.h"
#include "SortNMergeUtilities.h"
#include "StatisticsMaps.h"

using namespace std;


class CMosaikAligner {
public:
	// constructor
	CMosaikAligner(
	    unsigned char hashSize, 
	    CAlignmentThread::AlignerAlgorithmType algorithmType, 
	    CAlignmentThread::AlignerModeType algorithmMode, 
	    unsigned char numThreads, 
	    const string inputCommandLine);
	// destructor
	~CMosaikAligner(void);
	// aligns the read archive chromosome by chromosome
	void AlignReadArchiveLowMemory(void);
	// enables the alignment candidate threshold
	void EnableAlignmentCandidateThreshold(
	    const unsigned short alignmentCandidateThreshold);
	// enables the banded Smith-Waterman algorithm
	void EnableBandedSmithWaterman(const unsigned int bandwidth);
	// enables low-memory algorithm
	void EnableLowMemory(void);
	// enables SOLiD colorspace translation
	void EnableColorspace(const string& basespaceReferenceFilename);
	// enables the hash position threshold
	void EnableHashPositionThreshold(const unsigned short hashPositionThreshold);
	// enable the hash region threshold
	void EnableHashRegionThreshold(const unsigned short hashRegionThreshold);
	// enables the use of the jump database
	void EnableJumpDB(
	    const string& filenameStub, 
	    const unsigned int cacheSizeMB, 
	    const bool keepKeysInMemory, 
	    const bool keepPositionsInMemory);
	// enables the local alignment search
	void EnableLocalAlignmentSearch(const unsigned int radius);
	// enables paired-end read output
	void EnablePairedEndOutput(void);
	// sets mapping quality threshold for stat map
	void SetStatMappingQuality(const unsigned char mq);
	void SetLocalAlignmentSearchMqThreshold (
	    const unsigned char LocalAlignmentSearchHighMqThreshold, 
	    const unsigned char LocalAlignmentSearchLowMqThreshold );
	// enables reporting of unaligned reads
	//void EnableUnalignedReadReporting(const string& unalignedReadReportFilename);
	// enables special references checker
	void EnableSpecialReference (const string referencePrefix);
	// sets special hashes percentage
	void SetSpecialHashCount (const unsigned int count);
	// sets the filenames used by the aligner
	void SetFilenames(const string& inputReadArchiveFilename, 
	    const string& outputReadArchiveFilename, 
	    const string& referenceFilename);
	// enables the use of the aligned read length when calculating mismatches
	void UseAlignedReadLengthForMismatchCalculation(void);
	// sets quiet mode
	void SetQuietMode(void);
	void SetIsNoLog(void);
	// outputs multiply mapped alignments
	void OutputStdout(void);
	void OutputMultiply(void);
	void SetPeNeuralNetworkFilename(const string& neuralNetworkFilename);
	void SetSeNeuralNetworkFilename(const string& neuralNetworkFilename);
	void EnableZnTag();

private:
	// copy constructor
	CMosaikAligner(const CMosaikAligner& copy);
	// assign operator
	CMosaikAligner& operator=(const CMosaikAligner& copy);
	// denotes the active alignment algorithm
	CAlignmentThread::AlignerAlgorithmType mAlgorithm;
	// denotes the active alignment mode
	CAlignmentThread::AlignerModeType mMode;
	// stores the alignment configuration
	CAlignmentThread::AlignerSettings mSettings;
	// stores the filter configuration
	CAlignmentThread::FilterSettings mFilters;
	// stores our boolean flags
	CAlignmentThread::FlagData mFlags;
	// stores the statistical counters
	CAlignmentThread::StatisticsCounters mStatisticsCounters;
	// stores the statistical maps
	CStatisticsMaps mStatisticsMaps;
	// bam writers
	CAlignmentThread::BamWriters mBams;
	// special reference
	CAlignmentThread::SReference mSReference;
	// estimates the appropriate hash table size
	static unsigned char CalculateHashTableSize(const unsigned int referenceLength, const unsigned char hashSize);
	// hashes the reference sequence
	void HashReferenceSequence(MosaikReadFormat::CReferenceSequenceReader& refseq);
	// initializes the hash tables
	void InitializeHashTables(
	    const unsigned char bitSize, 
	    const unsigned int  begin, 
	    const unsigned int  end, 
	    const unsigned int  offset, 
	    const bool          useLowMemory, 
	    const unsigned int  expectedMemory, 
	    const bool          bubbleSpecialHashes);
	// the reference sequence
	char* mReference;
	// the length of the reference sequence
	unsigned int mReferenceLength;
	// the hash-table associated with the specified alignment algortihm
	CAbstractDnaHash* mpDNAHash;
	// temporary output files for chromosome-by-chromosome alignment
	vector<string> outputFilenames;
	// reference groups for low-memory algorithm
	//            start_ref,    length
	vector<pair <unsigned int, unsigned int> > referenceGroups;
	// command line
	string commandLine;
	// neuralNetworkFilename
	string mPeNeuralNetworkFilename;
	string mSeNeuralNetworkFilename;
	// read groups
	vector<MosaikReadFormat::ReadGroup> readGroups;
	map<unsigned int, MosaikReadFormat::ReadGroup> readGroupsMap;
	// merge the aligned archives generated by chromosome-by-chromosome alignment
	void MergeArchives(void);
	void PrintStatistics(void);
	void GroupReferences(const vector<ReferenceSequence>& referenceSequences);
	void GetHashStatistics( 
	    vector<unsigned int>& nHashs, 
	    vector<unsigned int>& expectedMemories, 
	    uint64_t& nTotalHash, 
	    const vector<ReferenceSequence>& referenceSequences);
	// aligns the read archive
	void AlignReadArchive(
		MosaikReadFormat::CReadReader&      in, 
		MosaikReadFormat::CAlignmentWriter& out, 
		unsigned int*                       pRefBegin, 
		unsigned int*                       pRefEnd, 
		char**                              pRefSpecies, 
		bool*                               pRefSpecial, 
		char**                              pBsRefSeqs,
		const unsigned int                  referenceOffset);
};

