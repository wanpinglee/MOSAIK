// ***************************************************************************
// AlignerMain.cpp - collects command-line parameters for MosaikAligner.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include <iostream>
#include <sstream>
#include "AlignmentThread.h"
#include "Benchmark.h"
#include "ConsoleUtilities.h"
#include "FileUtilities.h"
#include "Mosaik.h"
#include "MosaikAligner.h"
#include "Options.h"
#include "PairwiseUtilities.h"
#include "ReferenceSequenceStatus.h"
#include "SequenceUtilities.h"

using namespace std;

// define our default values
string DEFAULT_ALGORITHM = "all";
string DEFAULT_MODE      = "all";

unsigned char DEFAULT_HASH_SIZE      = 15;
unsigned char DEFAULT_NUM_MISMATCHES = 4;
unsigned int DEFAULT_BANDWIDTH       = 9;
unsigned int DEFAULT_NUM_THREADS     = 1;

#define MIN_HASH_SIZE     4
#define MAX_HASH_SIZE     32

// create a configuration variable struct
struct ConfigurationSettings {

	// flags
	bool CheckAlignmentQuality;
	bool CheckMinAlignment;
	bool CheckMinAlignmentPercent;
	bool CheckMismatchPercent;
	bool CheckNumMismatches;
	bool EnableAlignmentCandidateThreshold;
	bool EnableColorspace;
	bool EnableDoubleHashHits;
	bool HasAlgorithm;
	bool HasAlignmentsFilename;
	bool HasBandwidth;
	bool HasBasespaceReferencesFilename;
	bool HasGapExtendPenalty;
	bool HasGapOpenPenalty;
	bool HasHashSize;
	bool HasHomoPolymerGapOpenPenalty;
	bool HasJumpCacheMemory;
	bool HasLocalAlignmentSearchRadius;
	bool HasMatchScore;
	bool HasMismatchScore;
	bool HasMode;
	bool HasNumThreads;
	bool HasReadsFilename;
	bool HasReferencesFilename;
	bool KeepJumpKeysOnDisk;
	bool KeepJumpPositionsOnDisk;
	bool LimitHashPositions;
	bool RecordUnalignedReads;
	bool UseAlignedLengthForMismatches;
	bool UseJumpDB;
	bool UseLowMemory;

	// filenames
	string AlignmentsFilename;
	string BasespaceReferencesFilename;
	string JumpFilenameStub;
	string ReadsFilename;
	string ReferencesFilename;
	string UnalignedReadsFilename;

	// parameters
	double MinimumAlignmentPercentage;
	double MismatchPercent;
	float GapExtendPenalty;
	float GapOpenPenalty;
	float HomoPolymerGapOpenPenalty;
	float MatchScore;
	float MismatchScore;
	string Algorithm;
	string Mode;
	unsigned char AlignmentCandidateThreshold;
	unsigned char AlignmentQualityThreshold;
	unsigned char HashSize;
	unsigned int Bandwidth;
	unsigned int HashPositionThreshold;
	unsigned int JumpCacheMemory;
	unsigned int LocalAlignmentSearchRadius;
	unsigned int MinimumAlignment;
	unsigned int NumMismatches;
	unsigned int NumThreads;

	// constructor
	ConfigurationSettings()
		: CheckAlignmentQuality(false)
		, CheckMinAlignment(false)
		, CheckMinAlignmentPercent(false)
		, CheckMismatchPercent(false)
		, CheckNumMismatches(false)
		, EnableAlignmentCandidateThreshold(false)
		, EnableColorspace(false)
		, EnableDoubleHashHits(false)
		, HasAlgorithm(false)
		, HasAlignmentsFilename(false)
		, HasBandwidth(false)
		, HasBasespaceReferencesFilename(false)
		, HasGapExtendPenalty(false)
		, HasGapOpenPenalty(false)
		, HasHashSize(false)
		, HasHomoPolymerGapOpenPenalty(false)
		, HasJumpCacheMemory(false)
		, HasLocalAlignmentSearchRadius(false)
		, HasMatchScore(false)
		, HasMismatchScore(false)
		, HasMode(false)
		, HasNumThreads(false)
		, HasReadsFilename(false)
		, HasReferencesFilename(false)
		, KeepJumpKeysOnDisk(false)
		, KeepJumpPositionsOnDisk(false)
		, LimitHashPositions(false)
		, RecordUnalignedReads(false)
		, UseAlignedLengthForMismatches(false)
		, UseJumpDB(false)
		, UseLowMemory(false)
		, Algorithm(DEFAULT_ALGORITHM)
		, Mode(DEFAULT_MODE)
		, HashSize(DEFAULT_HASH_SIZE)
		, JumpCacheMemory(0)
		, NumMismatches(DEFAULT_NUM_MISMATCHES)
		, NumThreads(DEFAULT_NUM_THREADS)
	{}
};

int main(int argc, char* argv[]) {

	CConsole::Initialize();
	ConfigurationSettings settings;

	printf("------------------------------------------------------------------------------\n");
	printf("Mosaik"); CConsole::Red(); printf("Aligner"); CConsole::Reset();
	printf(" %u.%u.%04u                                              %s\n", 
		MOSAIK_MAJOR_VERSION, MOSAIK_MINOR_VERSION, MOSAIK_BUILD_VERSION, MOSAIK_VERSION_DATE);
	printf("Michael Stromberg & Wan-Ping Lee  Marth Lab, Boston College Biology Department\n");
	printf("------------------------------------------------------------------------------\n\n");

	// =================================
	// configure the command line parser
	// =================================

	// set general info about the program
	COptions::SetProgramInfo("MosaikAligner", "pairwise aligns a MOSAIK read file", "-in <filename> -out <filename> -ia <filename>");

	// add the input/output options
	OptionGroup* pIoOpts = COptions::CreateOptionGroup("Input/output: (required)");
	COptions::AddValueOption("-ia",  "MOSAIK reference filename", "the input reference file",  "An input MOSAIK reference file",  settings.HasReferencesFilename, settings.ReferencesFilename, pIoOpts);
	COptions::AddValueOption("-in",  "MOSAIK read filename",      "the input read file",       "An input MOSAIK read file",       settings.HasReadsFilename,      settings.ReadsFilename,      pIoOpts);
	COptions::AddValueOption("-out", "MOSAIK alignment filename", "the output alignment file", "An output MOSAIK alignment file", settings.HasAlignmentsFilename, settings.AlignmentsFilename, pIoOpts);
	COptions::AddValueOption("-ibs", "MOSAIK reference filename", "enables colorspace to basespace conversion using the supplied BASESPACE reference archive",  "",  settings.HasBasespaceReferencesFilename, settings.BasespaceReferencesFilename, pIoOpts);

	// add the essential options
	OptionGroup* pEssentialOpts = COptions::CreateOptionGroup("Essential parameters");
	COptions::AddValueOption("-a",  "algorithm",  "alignment algorithm: fast, single, multi, or all",  "", settings.HasAlgorithm, settings.Algorithm, pEssentialOpts, DEFAULT_ALGORITHM);
	COptions::AddValueOption("-m",  "mode",       "alignment mode: unique or all",                     "", settings.HasMode,      settings.Mode,      pEssentialOpts, DEFAULT_MODE);
	COptions::AddValueOption("-hs", "hash size",  "hash size [4 - 32]",                                "", settings.HasHashSize,  settings.HashSize,  pEssentialOpts, DEFAULT_HASH_SIZE);

	// add the filtering options
	OptionGroup* pFilterOpts = COptions::CreateOptionGroup("Filtering");
	COptions::AddValueOption("-act",  "threshold",      "the alignment candidate threshold (length)", "", settings.EnableAlignmentCandidateThreshold, settings.AlignmentCandidateThreshold, pFilterOpts);
	COptions::AddOption("-dh", "require at least two hash hits",                                          settings.EnableDoubleHashHits,                                                    pFilterOpts);
	COptions::AddValueOption("-ls",   "radius",         "enable local alignment search for PE reads", "", settings.HasLocalAlignmentSearchRadius,     settings.LocalAlignmentSearchRadius,  pFilterOpts);
	COptions::AddValueOption("-mhp",  "hash positions", "enable an alignment quality threshold",      "", settings.LimitHashPositions,                settings.HashPositionThreshold,       pFilterOpts);
	COptions::AddValueOption("-min",  "# nucleotides",  "enable an alignment quality threshold",      "", settings.CheckMinAlignment,                 settings.MinimumAlignment,            pFilterOpts);
	COptions::AddValueOption("-minp", "percent",        "the # of mismatches allowed",                "", settings.CheckMinAlignmentPercent,          settings.MinimumAlignmentPercentage,  pFilterOpts);
	COptions::AddValueOption("-mm",   "mismatches",     "the # of mismatches allowed",                "", settings.CheckNumMismatches,                settings.NumMismatches,               pFilterOpts);
	COptions::AddValueOption("-mmp",  "threshold",      "enable an alignment quality threshold",      "", settings.CheckMismatchPercent,              settings.MismatchPercent,             pFilterOpts);
	COptions::AddOption("-mmal", "when enabled, unaligned portions of the read will not count as a mismatch", settings.UseAlignedLengthForMismatches,                                       pFilterOpts);

	// TODO: we need to move the alignment quality calculation up to ApplyReadFilters in order to make this option useable
	//COptions::AddValueOption("-aq",   "threshold",      "enable an alignment quality threshold",      "", settings.CheckAlignmentQuality,             settings.AlignmentQualityThreshold,   pFilterOpts);

	// add the performance options
	OptionGroup* pPerformanceOpts = COptions::CreateOptionGroup("Performance");
	COptions::AddValueOption("-p",  "processors", "uses the specified number of processors", "", settings.HasNumThreads, settings.NumThreads, pPerformanceOpts);
	COptions::AddValueOption("-bw", "bandwidth",  "specifies the Smith-Waterman bandwidth", "",  settings.HasBandwidth,  settings.Bandwidth,  pPerformanceOpts, DEFAULT_BANDWIDTH);
	COptions::AddOption("-lm",                    "keeps the keys file on disk",              settings.UseLowMemory,                                 pPerformanceOpts);

	// add the jump database options
	OptionGroup* pJumpOpts = COptions::CreateOptionGroup("Jump database");
	COptions::AddValueOption("-j",  "filename stub", "uses the specified jump database",     "", settings.UseJumpDB,                 settings.JumpFilenameStub, pJumpOpts);
	//COptions::AddValueOption("-jc", "# of hashes",   "caches the most recently used hashes", "", settings.HasJumpCacheMemory,        settings.JumpCacheMemory,  pJumpOpts);
	COptions::AddOption("-kd",                       "keeps the keys file on disk",              settings.KeepJumpKeysOnDisk,                                 pJumpOpts);
	COptions::AddOption("-pd",                       "keeps the positions file on disk",         settings.KeepJumpPositionsOnDisk,                            pJumpOpts);

	// add the reporting options
	OptionGroup* pReportingOpts = COptions::CreateOptionGroup("Reporting");
	COptions::AddValueOption("-rur", "FASTQ filename", "stores unaligned reads in a FASTQ file", "", settings.RecordUnalignedReads, settings.UnalignedReadsFilename, pReportingOpts);

	// add the pairwise alignment scoring options
	OptionGroup* pPairwiseOpts = COptions::CreateOptionGroup("Pairwise Alignment Scores");
	COptions::AddValueOption("-ms",   "match score",        "the match score",             "", settings.HasMatchScore,                settings.MatchScore,                pPairwiseOpts, CPairwiseUtilities::MatchScore);
	COptions::AddValueOption("-mms",  "mismatch score",     "the mismatch score",          "", settings.HasMismatchScore,             settings.MismatchScore,             pPairwiseOpts, CPairwiseUtilities::MismatchScore);
	COptions::AddValueOption("-gop",  "gap open penalty",   "the gap open penalty",        "", settings.HasGapOpenPenalty,            settings.GapOpenPenalty,            pPairwiseOpts, CPairwiseUtilities::GapOpenPenalty);
	COptions::AddValueOption("-gep",  "gap extend penalty", "the gap extend penalty",      "", settings.HasGapExtendPenalty,          settings.GapExtendPenalty,          pPairwiseOpts, CPairwiseUtilities::GapExtendPenalty);
	COptions::AddValueOption("-hgop", "gap open penalty",   "enables the homopolymer gop", "", settings.HasHomoPolymerGapOpenPenalty, settings.HomoPolymerGapOpenPenalty, pPairwiseOpts, CPairwiseUtilities::HomoPolymerGapOpenPenalty);

	// parse the current command line
	COptions::Parse(argc, argv);

	// =============================
	// check for missing information
	// =============================

	bool foundError = false;
	ostringstream errorBuilder;
	const string ERROR_SPACER(7, ' ');

	if(settings.EnableAlignmentCandidateThreshold && settings.EnableDoubleHashHits) {
		errorBuilder << ERROR_SPACER << "Please specify either an alignment candidate threshold (-act) or double-hash hits (-dh). Double-hash hits are equivalent to '-act <hash size + 1>." << endl;
		foundError = true;
	}

	if((settings.HasJumpCacheMemory || settings.KeepJumpKeysOnDisk || settings.KeepJumpPositionsOnDisk) && !settings.UseJumpDB) {
		errorBuilder << ERROR_SPACER << "Jump database settings were specified, but the jump database was not explicitly chosen. Please use the -j parameter." << endl;
		foundError = true;
	}

	if(settings.UseJumpDB) {
		string keysFilename      = settings.JumpFilenameStub + "_keys.jmp";
		string metaFilename      = settings.JumpFilenameStub + "_meta.jmp";
		string positionsFilename = settings.JumpFilenameStub + "_positions.jmp";

		CFileUtilities::CheckFile(keysFilename.c_str(), true);
		CFileUtilities::CheckFile(metaFilename.c_str(), true);
		CFileUtilities::CheckFile(positionsFilename.c_str(), true);

		if(!settings.KeepJumpKeysOnDisk && !settings.KeepJumpPositionsOnDisk && settings.HasJumpCacheMemory) 
			settings.HasJumpCacheMemory = false;
	}

	if(!settings.CheckNumMismatches && !settings.CheckMismatchPercent && !settings.CheckAlignmentQuality) {
		settings.CheckNumMismatches = true;
	}

	// figure out which algorithm to use
	CSequenceUtilities::LowercaseSequence(settings.Algorithm);
	CAlignmentThread::AlignerAlgorithmType algorithmType  = CAlignmentThread::AlignerAlgorithm_ALL;

	if(settings.Algorithm == "fast")        algorithmType = CAlignmentThread::AlignerAlgorithm_FAST;
	else if(settings.Algorithm == "single") algorithmType = CAlignmentThread::AlignerAlgorithm_SINGLE;
	else if(settings.Algorithm == "multi")  algorithmType = CAlignmentThread::AlignerAlgorithm_MULTI;
	else if(settings.Algorithm == "all")    algorithmType = CAlignmentThread::AlignerAlgorithm_ALL;
	else {
		errorBuilder << ERROR_SPACER << "Unknown algorithm type. Please choose between 'fast', 'single', 'multi', or 'all'. The default value is '" << DEFAULT_ALGORITHM << "'." << endl;
		foundError = true;
	}

	// set the hash positions threshold
	if(settings.LimitHashPositions) {

		// make sure we're using the all algorithm
		if(algorithmType != CAlignmentThread::AlignerAlgorithm_ALL) {
			errorBuilder << ERROR_SPACER << "Setting the hash positions threshold is only applicable when using the 'all' algorithm. This can be set by using the '-a all' parameter." << endl;
			foundError = true;
		}

		if(settings.HashPositionThreshold == 0) {
			errorBuilder << ERROR_SPACER << "The hash position threshold should be larger than 0. Use the -mhp parameter to change the hash position threshold." << endl;
			foundError = true;
		}
	}

	// figure out which alignment mode to use
	CSequenceUtilities::LowercaseSequence(settings.Algorithm);
	CAlignmentThread::AlignerModeType modeType = CAlignmentThread::AlignerMode_ALL;

	if(settings.Mode == "unique")     modeType = CAlignmentThread::AlignerMode_UNIQUE;
	else if(settings.Mode == "all")   modeType = CAlignmentThread::AlignerMode_ALL;
	else {
		errorBuilder << ERROR_SPACER << "Unknown mode type. Please choose between 'unique' or 'all'. The default value is '" << DEFAULT_MODE << "'." << endl;
		foundError = true;
	}

	// set the maximum mismatch percentage
	if(settings.CheckMismatchPercent) {

		// make sure our value is within bounds
		if((settings.MismatchPercent < 0.0) || (settings.MismatchPercent > 1.0)) {
			errorBuilder << ERROR_SPACER << "The maximum mismatch percentage should be between 0.0 and 1.0." << endl;
			foundError = true;
		}

		CPairwiseUtilities::MaxMismatchPercent       = settings.MismatchPercent;
		CPairwiseUtilities::UseMismatchPercentFilter = true;
	}

	// set the minimum aligned percentage
	if(settings.CheckMinAlignmentPercent) {

		// make sure our value is within bounds
		if((settings.MinimumAlignmentPercentage < 0.0) || (settings.MinimumAlignmentPercentage > 1.0)) {
			errorBuilder << ERROR_SPACER << "The minimum alignment percentage should be between 0.0 and 1.0." << endl;
			foundError = true;
		}

		// assign the minimum percentage alignment
		CPairwiseUtilities::MinPercentAlignment          = settings.MinimumAlignmentPercentage;
		CPairwiseUtilities::UseMinAlignmentPercentFilter = true;
	}

	// check the minimum alignment quality
	if(settings.CheckAlignmentQuality && (settings.AlignmentQualityThreshold > 99)) {
		errorBuilder << ERROR_SPACER << "The alignment quality threshold should be between 0 and 99." << endl;
		foundError = true;
	}

	// set the hash size
	if(settings.HasHashSize && ((settings.HashSize < MIN_HASH_SIZE) || (settings.HashSize > MAX_HASH_SIZE))) {
		errorBuilder << ERROR_SPACER << "The hash size should be between " << MIN_HASH_SIZE << " and " << MAX_HASH_SIZE << ". The default value is " << DEFAULT_HASH_SIZE << "." << endl;
		foundError = true;
	}

	// set the number of threads
	if(settings.HasNumThreads && (settings.NumThreads < 1)) {
		errorBuilder << ERROR_SPACER << "At least one processor should be specified. Use the -p parameter to change the number of desired processors." << endl;
		foundError = true;
	}

	// set the Smith-Waterman bandwidth
	if(settings.HasBandwidth && ((settings.Bandwidth % 2) != 1)) {
		errorBuilder << ERROR_SPACER << "The bandwidth must be an odd number. Use the -bw parameter to change the bandwidth." << endl;
		foundError = true;
	}

	// test if the specified input files exist and are in the right format
	SequencingTechnologies seqTech;
	ReadStatus readStatus;
	MosaikReadFormat::CReadReader::CheckFile(settings.ReadsFilename, seqTech, readStatus, true);
	MosaikReadFormat::CReferenceSequenceReader::CheckFile(settings.ReferencesFilename, true);

	switch(seqTech) {
		case ST_454:
			if(!settings.HasHomoPolymerGapOpenPenalty) {
				CPairwiseUtilities::UseHomoPolymerGapOpenPenalty = true;
				settings.HasHomoPolymerGapOpenPenalty            = true;
				settings.HomoPolymerGapOpenPenalty               = CPairwiseUtilities::HomoPolymerGapOpenPenalty;
			}
			break;
		case ST_SOLID:
			settings.EnableColorspace = true;
			break;
	}

	// check if we have a supplied basespace reference archive when aligning a SOLiD read archive
	if(seqTech == ST_SOLID) {

		// force -ibs to be basespace and -ia to be colorspace if aligning SOLiD reads
		if(settings.HasBasespaceReferencesFilename && settings.HasReferencesFilename) {
			MosaikReadFormat::CReferenceSequenceReader csRef, bsRef;

			// retrieve the colorspace status
			csRef.Open(settings.ReferencesFilename);
			ReferenceSequenceStatus csStatus = csRef.GetStatus();
			csRef.Close();

			// retrieve the basespace status
			bsRef.Open(settings.BasespaceReferencesFilename);
			ReferenceSequenceStatus bsStatus = bsRef.GetStatus();
			bsRef.Close();

			if(csStatus != REF_COLORSPACE) {
				errorBuilder << ERROR_SPACER << "Expected to find a colorspace reference sequence archive (" << settings.ReferencesFilename << ") with the -ia parameter, but found a basespace reference sequence archive." << endl;
				foundError = true;
			}

			if(bsStatus != REF_UNKNOWN) {
				errorBuilder << ERROR_SPACER << "Expected to find a basespace reference sequence archive (" << settings.BasespaceReferencesFilename << ") with the -ibs parameter, but found a colorspace reference sequence archive." << endl;
				foundError = true;
			}
		} 

		if(!settings.HasBasespaceReferencesFilename)  { 
			errorBuilder << ERROR_SPACER << "When aligning SOLiD read archives, both a colorspace reference archive AND a basespace reference archive are required. Use the -ibs parameter to supply the basespace reference archive filename." << endl;
			foundError = true;
		}

		if( settings.UseLowMemory ) {
			errorBuilder << ERROR_SPACER << "The low-memory algorithm does not support for SOLiD reads yet." << endl;
			foundError = true;
		}
	}

	// print the errors if any were found
	if(foundError) {

		CConsole::Red();
		printf("ERROR: Some problems were encountered when parsing the command line options:\n");
		CConsole::Reset();

		printf("%s\n", errorBuilder.str().c_str());
		printf("For a complete list of command line options, type \"%s -h\"\n", argv[0]);
		exit(1);
	}

	// ===================================================
	// Parse configuration strings and set class variables
	// ===================================================

	// set the minimum alignment quality
	if(settings.CheckAlignmentQuality) {
		CPairwiseUtilities::MinAlignmentQuality          = settings.AlignmentQualityThreshold;
		CPairwiseUtilities::UseMinAlignmentQualityFilter = true;
	}

	// set the maximum number of mismatches
	if(settings.CheckNumMismatches) {
		CPairwiseUtilities::MaxNumMismatches  = settings.NumMismatches;
		CPairwiseUtilities::UseMismatchFilter = true;
	}

	// set the minimum number of aligned nucleotides
	if(settings.CheckMinAlignment) {
		CPairwiseUtilities::MinAlignment          = settings.MinimumAlignment;
		CPairwiseUtilities::UseMinAlignmentFilter = true;
	}

	// set the Smith-Waterman scores
	if(settings.HasMatchScore || settings.HasMismatchScore || settings.HasGapOpenPenalty || settings.HasGapExtendPenalty) {
		if(settings.HasMatchScore)       CPairwiseUtilities::MatchScore       = settings.MatchScore;
		if(settings.HasMismatchScore)    CPairwiseUtilities::MismatchScore    = settings.MismatchScore;
		if(settings.HasGapOpenPenalty)   CPairwiseUtilities::GapOpenPenalty   = settings.GapOpenPenalty;
		if(settings.HasGapExtendPenalty) CPairwiseUtilities::GapExtendPenalty = settings.GapExtendPenalty;
	}

	// set the Smith-Waterman homo-polymer gap open penalty
	if(settings.HasHomoPolymerGapOpenPenalty) {
		CPairwiseUtilities::HomoPolymerGapOpenPenalty    = settings.HomoPolymerGapOpenPenalty;
		CPairwiseUtilities::UseHomoPolymerGapOpenPenalty = true;
	}

	// show warning message about unique alignments
	if(((readStatus & RS_PAIRED_END_READ) != 0) && (modeType != CAlignmentThread::AlignerMode_ALL)) {
		cout << "WARNING: A paired-end read archive was detected and the aligner mode (-m parameter) was not set to 'all'. Paired-end resolution in MosaikSort will be limited to unique vs unique reads.\n" << endl << endl;
	}

	// show warning messages dealing with the local alignment search radius
	if(settings.HasLocalAlignmentSearchRadius) {

		// show the warning message if we have a SE read archive
		if((readStatus & RS_SINGLE_END_READ) != 0) {
			cout << "WARNING: A single-end read archive was detected and the local alignment search was enabled. Local alignment search only works with paired-end reads.\n" << endl << endl;
			settings.HasLocalAlignmentSearchRadius = false;
		} else { 

			// show the warning message if we have a PE read archive with no mean fragment length
			MosaikReadFormat::CReadReader in;
			in.Open(settings.ReadsFilename);
			MosaikReadFormat::ReadGroup readGroup = in.GetReadGroup();
			in.Close();

			if(readGroup.MedianFragmentLength == 0) {
				cout << "WARNING: Local alignment search only works when the median fragment length (-mfl parameter) has been specified in MosaikBuild.\n" << endl << endl;
				settings.HasLocalAlignmentSearchRadius = false;
			}
		}
	}

	// show warning message about using the local alignment search with SE read archives
	if(((readStatus & RS_SINGLE_END_READ) != 0) && settings.HasLocalAlignmentSearchRadius) {
		cout << "WARNING: A single-end read archive was detected and the local alignment search was enabled. Local alignment search only works with paired-end reads.\n" << endl << endl;
		settings.HasLocalAlignmentSearchRadius = false;
	}

	// start benchmarking
	CBenchmark bench;
	bench.Start();

	// create our aligner
	CMosaikAligner ma(settings.HashSize, algorithmType, modeType, settings.NumThreads);

	// ===============
	// enable features
	// ===============

	// enable the hash positions threshold
	if(settings.LimitHashPositions) ma.EnableHashPositionThreshold(settings.HashPositionThreshold);

	// enable the alignment candidate threshold
	if(settings.EnableAlignmentCandidateThreshold) ma.EnableAlignmentCandidateThreshold(settings.AlignmentCandidateThreshold);

	// enable unaligned read reporting if specified
	if(settings.RecordUnalignedReads) ma.EnableUnalignedReadReporting(settings.UnalignedReadsFilename);

	// enable colorspace (SOLiD)
	if(settings.EnableColorspace) ma.EnableColorspace(settings.BasespaceReferencesFilename);

	// enable double-hash hits
	if(settings.EnableDoubleHashHits) ma.EnableAlignmentCandidateThreshold(settings.HashSize + 1);

	// enable entire read length mismatch checking
	if(settings.UseAlignedLengthForMismatches) ma.UseAlignedReadLengthForMismatchCalculation();

	// enable the jump database
	if(settings.UseJumpDB) ma.EnableJumpDB(settings.JumpFilenameStub, settings.JumpCacheMemory, !settings.KeepJumpKeysOnDisk, !settings.KeepJumpPositionsOnDisk);

	// enable the local alignment search
	if(settings.HasLocalAlignmentSearchRadius) ma.EnableLocalAlignmentSearch(settings.LocalAlignmentSearchRadius);

	// set the Smith-Waterman bandwidth
	if(settings.HasBandwidth) ma.EnableBandedSmithWaterman(settings.Bandwidth);

	// enable low-memory algorithm
	if(settings.UseLowMemory) ma.EnableLowMemory();

	// =============
	// set filenames
	// =============

	ma.SetFilenames(settings.ReadsFilename, settings.AlignmentsFilename, settings.ReferencesFilename);

	// ====================
	// echo enabled options
	// ====================

	cout << "- Using the following alignment algorithm: ";
	switch(algorithmType) {
		case CAlignmentThread::AlignerAlgorithm_FAST:
			cout << "single position (fast)" << endl;
			break;
		case CAlignmentThread::AlignerAlgorithm_SINGLE:
			cout << "single position" << endl;
			break;
		case CAlignmentThread::AlignerAlgorithm_MULTI:
			cout << "multiple position" << endl;
			break;
		case CAlignmentThread::AlignerAlgorithm_ALL:
			cout << "all positions" << endl;
			break;
		default:
			cout << "ERROR: Unknown alignment algorithm specified." << endl;
			exit(1);
			break;
	}

	cout << "- Using the following alignment mode: ";
	switch(modeType) {
		case CAlignmentThread::AlignerMode_ALL:
			cout << "aligning reads to all possible locations" << endl;
			break;
		case CAlignmentThread::AlignerMode_UNIQUE:
			cout << "only aligning unique reads" << endl;
			break;
		default:
			cout << "ERROR: Unknown alignment mode specified." << endl;
			exit(1);
			break;
	}

	if(settings.CheckAlignmentQuality)    cout << "- Using an alignment quality threshold of " << CPairwiseUtilities::MinAlignmentQuality << endl;
	if(settings.CheckNumMismatches)       cout << "- Using a maximum mismatch threshold of " << CPairwiseUtilities::MaxNumMismatches << endl;
	if(settings.CheckMismatchPercent)     cout << "- Using a maximum mismatch percent threshold of " << CPairwiseUtilities::MaxMismatchPercent << endl;
	if(settings.CheckMinAlignment)        cout << "- Using a minimum alignment threshold of " << CPairwiseUtilities::MinAlignment << endl;
	if(settings.CheckMinAlignmentPercent) cout << "- Using a minimum percent alignment threshold of " << CPairwiseUtilities::MinPercentAlignment << endl;
	if(settings.HasHashSize)              cout << "- Using a hash size of " << (unsigned int)settings.HashSize << endl;
	if(settings.EnableDoubleHashHits)     cout << "- Using double-hash hits" << endl;
	if(settings.EnableColorspace)         cout << "- Aligning in colorspace (SOLiD)" << endl;
	if(settings.HasNumThreads)            cout << "- Using " << (short)settings.NumThreads << (settings.NumThreads > 1 ? " processors" : " processor") << endl;
	if(settings.HasBandwidth)             cout << "- Using a Smith-Waterman bandwidth of " << settings.Bandwidth << endl;

	if(settings.EnableAlignmentCandidateThreshold) 
		cout << "- Using an alignment candidate threshold of " << (unsigned short)settings.AlignmentCandidateThreshold << "bp." << endl;

	if(settings.HasLocalAlignmentSearchRadius)
		cout << "- Using a local alignment search radius of " << settings.LocalAlignmentSearchRadius << "bp." << endl;

	if(settings.LimitHashPositions) cout << "- Setting hash position threshold to " << settings.HashPositionThreshold << endl;

	if(settings.UseJumpDB) {
		cout << "- Using a jump database for hashing";
		if(settings.HasJumpCacheMemory) cout << " with a " << settings.JumpCacheMemory << " element cache";
		cout << ".";

		if(!settings.KeepJumpKeysOnDisk && !settings.KeepJumpPositionsOnDisk)     cout << " Storing keys & positions in memory.";
		else if(!settings.KeepJumpKeysOnDisk && settings.KeepJumpPositionsOnDisk) cout << " Storing keys in memory.";
		else if(settings.KeepJumpKeysOnDisk && !settings.KeepJumpPositionsOnDisk) cout << " Storing positions in memory.";
		cout << endl;
	}

	if(settings.RecordUnalignedReads)
		cout << "- Reporting all unaligned reads to " << settings.UnalignedReadsFilename << "." << endl;

	if(settings.HasHomoPolymerGapOpenPenalty) 
		cout << "- Using a homo-polymer gap open penalty of " << CPairwiseUtilities::HomoPolymerGapOpenPenalty << endl;

	if(settings.HasMatchScore || settings.HasMismatchScore || settings.HasGapOpenPenalty || settings.HasGapExtendPenalty) 
		cout << "- Updating Smith-Waterman scoring scheme (match, mismatch, gap open, gap extend): (" <<
		CPairwiseUtilities::MatchScore << ", " << CPairwiseUtilities::MismatchScore << ", " <<
		CPairwiseUtilities::GapOpenPenalty << ", " << CPairwiseUtilities::GapExtendPenalty << ")" << endl;

	// ==============
	// Start aligning
	// ==============

	ma.AlignReadArchiveLowMemory();
	//ma.AlignReadArchive();

	// ==================
	// Show total runtime
	// ==================

	// stop benchmarking
	bench.Stop();

	// show the benchmarking results
	cout << endl;
	bench.DisplayTime("MosaikAligner");

	return 0;
}
