// ***************************************************************************
// SortMain.cpp - collects command-line parameters for MosaikSort.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include <iomanip>
#include <iostream>
#include <sstream>
#include "AlignmentReader.h"
#include "Benchmark.h"
#include "ConsoleUtilities.h"
#include "Mosaik.h"
#include "Options.h"
#include "PairedEndSort.h"
#include "SingleEndSort.h"

using namespace std;

// define our default values
const unsigned int DEFAULT_CACHE_SIZE = 6000000;

// create a configuration variable struct
struct ConfigurationSettings {

	// flags
	bool AllowAllUniqueFragmentLengths;
	bool DisableFragmentAlignmentQuality;
	bool HasCacheSize;
	bool HasConfidenceInterval;
	bool HasDuplicateDirectory;
	bool HasInputMosaikAlignmentFilename;
	bool HasOutputMosaikAlignmentFilename;
	bool IgnoreUM;
	bool IgnoreUU;
	bool IgnoreUniqueOrphans;
	bool ResolveMM;
	bool SampleAllFragmentLengths;
	bool UseConsedRenaming;
	bool UseNonUniqueReads;
	bool HasInputFastqFilename;
	bool HasInputFastq2Filename;

	// filenames
	string DuplicateDirectory;
	string InputMosaikAlignmentFilename;
	string OutputMosaikAlignmentFilename;
	string UnresolvedMosaikAlignmentFilename;
	string InputFastqFilename;
	string InputFastq2Filename;

	// parameters
	double ConfidenceInterval;
	unsigned int CacheSize;

	// constructor
	ConfigurationSettings()
		: AllowAllUniqueFragmentLengths(false)
		, DisableFragmentAlignmentQuality(false)
		, HasCacheSize(false)
		, HasConfidenceInterval(false)
		, HasDuplicateDirectory(false)
		, HasInputMosaikAlignmentFilename(false)
		, HasOutputMosaikAlignmentFilename(false)
		, IgnoreUM(false)
		, IgnoreUU(false)
		, IgnoreUniqueOrphans(false)
		, ResolveMM(false)
		, SampleAllFragmentLengths(false)
		, UseConsedRenaming(false)
		, UseNonUniqueReads(false)
		, HasInputFastqFilename(false)
		, HasInputFastq2Filename(false)
		, ConfidenceInterval(DEFAULT_CONFIDENCE_INTERVAL)
		, CacheSize(DEFAULT_CACHE_SIZE)
	{}
};

int main(int argc, char* argv[]) {

	CConsole::Initialize();
	ConfigurationSettings settings;

	printf("------------------------------------------------------------------------------\n");
	printf("Mosaik"); CConsole::Red(); printf("Sort"); CConsole::Reset();
	printf(" %u.%u.%04u                                                 %s\n", 
		MOSAIK_MAJOR_VERSION, MOSAIK_MINOR_VERSION, MOSAIK_BUILD_VERSION, MOSAIK_VERSION_DATE);
	printf("Michael Stromberg & Wan-Ping Lee  Marth Lab, Boston College Biology Department\n");
	printf("------------------------------------------------------------------------------\n\n");

	// =================================
	// configure the command line parser
	// =================================

	// set general info about the program
	COptions::SetProgramInfo("MosaikSort", "sorts the alignments by position", "-in <filename> -out <filename>");

	// add the input/output options
	OptionGroup* pIOOpts = COptions::CreateOptionGroup("Input & Output");
	COptions::AddOption("-consed",                                "appends a number to read names for consed compatibility",                     settings.UseConsedRenaming,                                                        pIOOpts);
	COptions::AddValueOption("-dup", "directory",                 "enables duplicate filtering with databases in the specified directory", "",   settings.HasDuplicateDirectory,            settings.DuplicateDirectory,            pIOOpts);
	COptions::AddValueOption("-in",  "MOSAIK alignment filename", "the input MOSAIK alignment file",      "An input MOSAIK alignment filename",  settings.HasInputMosaikAlignmentFilename,  settings.InputMosaikAlignmentFilename,  pIOOpts);
	//COptions::AddValueOption("-q",   "FASTQ filename",            "the original FASTQ filename",          "An input FASTQ filename",             settings.HasInputFastqFilename,  settings.InputFastqFilename,  pIOOpts);
	//COptions::AddValueOption("-q2",  "FASTQ filename",            "the original 2nd mate FASTQ filename", "An input FASTQ filename",             settings.HasInputFastq2Filename,  settings.InputFastq2Filename,  pIOOpts);
	COptions::AddValueOption("-out", "MOSAIK alignment filename", "the output MOSAIK alignment file",     "An output MOSAIK alignment filename", settings.HasOutputMosaikAlignmentFilename, settings.OutputMosaikAlignmentFilename, pIOOpts);
	COptions::AddValueOption("-mem", "# of alignments in memory", "sets the sorting cache size",                                           "",   settings.HasCacheSize,                     settings.CacheSize,                     pIOOpts, DEFAULT_CACHE_SIZE);

	// add the single-end options
	OptionGroup* pSingleOpts = COptions::CreateOptionGroup("Single-end Options");
	COptions::AddOption("-nu", "include non-unique reads", settings.UseNonUniqueReads, pSingleOpts);

	// add the paired-end options
	OptionGroup* pPairedOpts = COptions::CreateOptionGroup("Paired-end Options");
	COptions::AddOption("-afl",                             "allows all fragment lengths when resolving unique vs unique read pairs", settings.AllowAllUniqueFragmentLengths,                      pPairedOpts);
	COptions::AddValueOption("-ci", "confidence interval",  "sets the desired fragment length confidence interval",  "",              settings.HasConfidenceInterval, settings.ConfidenceInterval, pPairedOpts, DEFAULT_CONFIDENCE_INTERVAL);
	COptions::AddOption("-sa",                              "sample all unique vs unique fragment lengths",                           settings.SampleAllFragmentLengths,                           pPairedOpts);
	COptions::AddOption("-dfaq",                            "disables fragment alignment quality calculation",                        settings.DisableFragmentAlignmentQuality,                    pPairedOpts);

	// add the paired-end resolution options
	OptionGroup* pResolveOpts = COptions::CreateOptionGroup("Paired-end Resolution");
	COptions::AddOption("-iuo", "ignore unique orphaned reads",            settings.IgnoreUniqueOrphans, pResolveOpts);
	COptions::AddOption("-iuu", "ignore unique vs unique read pairs",      settings.IgnoreUU,            pResolveOpts);
	COptions::AddOption("-ium", "ignore unique vs multiple read pairs",    settings.IgnoreUM,            pResolveOpts);
	COptions::AddOption("-rmm", "resolve multiple vs multiple read pairs", settings.ResolveMM,           pResolveOpts);

	// parse the current command line
	COptions::Parse(argc, argv);

	// test to see if the specified input files exist
	SequencingTechnologies seqTech;
	AlignmentStatus alignmentStatus;
	MosaikReadFormat::CAlignmentReader::CheckFile(settings.InputMosaikAlignmentFilename, seqTech, alignmentStatus, true);

	// activate paired-end sorting mode
	bool isAlignmentArchivePairedEnd = false;
	if((alignmentStatus & AS_PAIRED_END_READ) != 0) isAlignmentArchivePairedEnd = true;

	// =============================
	// check for missing information
	// =============================
	
	bool foundError = false;
	ostringstream errorBuilder;
	const string ERROR_SPACER(7, ' ');

	// convert the number of standard deviations
	if(settings.HasConfidenceInterval && ((settings.ConfidenceInterval <= 0.0) || (settings.ConfidenceInterval > 1.0))) {
		errorBuilder << ERROR_SPACER << "The confidence interval should be between 0.0 and 1.0." << endl;
		foundError = true;
	}

	// October 19th, 2010
	// patching paired-end archives needs two the original FASTQ files
	/*
	if ( isAlignmentArchivePairedEnd ) {
		bool noneFastqFiles = !settings.HasInputFastqFilename && !settings.HasInputFastq2Filename;
		bool bothFastqFiles = settings.HasInputFastqFilename && settings.HasInputFastq2Filename;
		if ( !noneFastqFiles && !bothFastqFiles ) {
			errorBuilder << ERROR_SPACER << "For paired-end data, both two original FASTQ files are needed, -q and -q2." << endl;
			foundError = true;
		}
	} else {
		if ( settings.HasInputFastq2Filename ) {
			errorBuilder << ERROR_SPACER << " For single-end data, please use -q." << endl;
			foundError = true;
		}
	}
	*/

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

	// check if the duplicate directory exists
	if(settings.HasDuplicateDirectory) {
		CFileUtilities::CheckDirectory(settings.DuplicateDirectory);

		const unsigned int directoryNameLen = (unsigned int)settings.DuplicateDirectory.size();
		if(settings.DuplicateDirectory[directoryNameLen - 1] != OS_DIRECTORY_SEPARATOR) settings.DuplicateDirectory += OS_DIRECTORY_SEPARATOR;
	}

	// test to see if the specified input files exist
	//SequencingTechnologies seqTech;
	//AlignmentStatus alignmentStatus;
	//MosaikReadFormat::CAlignmentReader::CheckFile(settings.InputMosaikAlignmentFilename, seqTech, alignmentStatus, true);

	// activate paired-end sorting mode
	//bool isAlignmentArchivePairedEnd = false;
	//if((alignmentStatus & AS_PAIRED_END_READ) != 0) isAlignmentArchivePairedEnd = true;

	// start benchmarking
	CBenchmark bench;
	bench.Start();

	if(isAlignmentArchivePairedEnd) {

		// resolve paired-end reads
		CPairedEndSort pes(settings.CacheSize);

		// display which types are being resolved
		printf("- resolving the following types of read pairs: ");
		if(!settings.IgnoreUniqueOrphans) printf("[unique orphans] ");
		if(!settings.IgnoreUU)            printf("[unique vs unique] ");
		if(!settings.IgnoreUM)            printf("[unique vs multiple] ");
		if(settings.ResolveMM)            printf("[multiple vs multiple] ");
		printf("\n");

		pes.ConfigureResolution(!settings.IgnoreUniqueOrphans, !settings.IgnoreUU, !settings.IgnoreUM, settings.ResolveMM);

		// enforce unique read pair constraints
		if(settings.AllowAllUniqueFragmentLengths) {
			printf("- allowing all fragment lengths for reads with unique mates\n");
			pes.EnableAllUniqueFragmentLengths();
		}

		// sample all read pairs
		if(settings.SampleAllFragmentLengths) {
			printf("- using entire data set for fragment length confidence interval calculation\n");
			pes.EnableFullFragmentLengthSampling();
		}

		// set the desired confidence interval
		if(settings.HasConfidenceInterval) {
			printf("- setting the confidence interval to %0.4f\n", settings.ConfidenceInterval);
			pes.SetConfidenceInterval(settings.ConfidenceInterval);
		}

		// enable read filtering
		if(settings.HasDuplicateDirectory) {
			printf("- enabling duplicate filtering\n");
			pes.EnableDuplicateFiltering(settings.DuplicateDirectory);
		}

		// enable consed renaming
		if(settings.UseConsedRenaming) {
			printf("- enabling consed renaming: mate number will be appended to the read name\n");
			pes.EnableConsedRenaming();
		}

		// disable fragment alignment quality calculation
		if(settings.DisableFragmentAlignmentQuality) {
			printf("- disabling fragment alignment quality calculation.\n");
			pes.DisableFragmentAlignmentQuality();
		}

		// patch the original fastq information
		//if(settings.HasInputFastqFilename) {
		//	printf("- patching the original FASTQ information.\n");
		//	pes.PatchFastq();
		//}

		// resolve the paired-end reads
		pes.ResolvePairedEndReads(settings.InputMosaikAlignmentFilename, settings.OutputMosaikAlignmentFilename);

	} else { 

		// sort single-end reads
		CSingleEndSort ses(settings.CacheSize);

		// enable unique resolution if the file was aligned in unique mode
		if(settings.UseNonUniqueReads && ((alignmentStatus & AS_UNIQUE_MODE) != 0)) {
			printf("- deactivating inclusion of non-unique reads since the archive was aligned in 'unique' mode\n");
			settings.UseNonUniqueReads = false;
		}

		// enable unique resolution if specified
		if(settings.UseNonUniqueReads) {
			printf("- using non-uniquely aligned reads\n");
			ses.EnableNonUniqueMode();
		}

		// enable duplicate read filtering
		if(settings.HasDuplicateDirectory) {
			printf("- enabling duplicate filtering\n");
			ses.EnableDuplicateFiltering(settings.DuplicateDirectory);
		}

		// enable consed renaming
		if(settings.UseConsedRenaming && settings.UseNonUniqueReads) {
			printf("- enabling consed renaming: alignment counter will be appended to the read name\n");
			ses.EnableConsedRenaming();
		}

		ses.SaveAlignmentsOrderedByPosition(settings.InputMosaikAlignmentFilename, settings.OutputMosaikAlignmentFilename);
	}

	// ==================
	// Show total runtime
	// ==================

	// stop benchmarking
	bench.Stop();

	// show the benchmarking results
	bench.DisplayTime("MosaikSort");

	return 0;
}
