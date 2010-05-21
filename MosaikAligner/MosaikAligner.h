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

#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include "AlignmentThread.h"
#include "AlignmentWriter.h"
#include "ReferenceSequence.h"
#include "Benchmark.h"
#include "ConsoleUtilities.h"
#include "DnaHash.h"
#include "JumpDnaHash.h"
#include "ReadReader.h"
#include "MultiDnaHash.h"
#include "PosixThreads.h"
#include "ProgressBar.h"
#include "ReferenceSequenceReader.h"
#include "UbiqDnaHash.h"

using namespace std;

class CMosaikAligner {
public:
	// constructor
	CMosaikAligner(unsigned char hashSize, CAlignmentThread::AlignerAlgorithmType algorithmType, 
		CAlignmentThread::AlignerModeType algorithmMode, unsigned char numThreads);
	// destructor
	~CMosaikAligner(void);
	// aligns the read archive
	void AlignReadArchive(void);
	// enables the alignment candidate threshold
	void EnableAlignmentCandidateThreshold(const unsigned short alignmentCandidateThreshold);
	// enables the banded Smith-Waterman algorithm
	void EnableBandedSmithWaterman(const unsigned int bandwidth);
	// enables SOLiD colorspace translation
	void EnableColorspace(const string& basespaceReferenceFilename);
	// enables the hash position threshold
	void EnableHashPositionThreshold(const unsigned short hashPositionThreshold);
	// enables the use of the jump database
	void EnableJumpDB(const string& filenameStub, const unsigned int cacheSizeMB, const bool keepKeysInMemory, const bool keepPositionsInMemory);
	// enables the local alignment search
	void EnableLocalAlignmentSearch(const unsigned int radius);
	// enables paired-end read output
	void EnablePairedEndOutput(void);
	// enables reporting of unaligned reads
	void EnableUnalignedReadReporting(const string& unalignedReadReportFilename);
	// sets the filenames used by the aligner
	void SetFilenames(const string& inputReadArchiveFilename, const string& outputReadArchiveFilename, const string& referenceFilename);
	// enables the use of the aligned read length when calculating mismatches
	void UseAlignedReadLengthForMismatchCalculation(void);
private:
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
	// estimates the appropriate hash table size
	static unsigned char CalculateHashTableSize(const unsigned int referenceLength, const unsigned char hashSize);
	// hashes the reference sequence
	void HashReferenceSequence(MosaikReadFormat::CReferenceSequenceReader& refseq);
	// initializes the hash tables
	void InitializeHashTables(const unsigned char bitSize);
	// the reference sequence
	char* mReference;
	// the length of the reference sequence
	unsigned int mReferenceLength;
	// the hash-table associated with the specified alignment algortihm
	CAbstractDnaHash* mpDNAHash;
};
