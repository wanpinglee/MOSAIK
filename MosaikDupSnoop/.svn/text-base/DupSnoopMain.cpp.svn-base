// ***************************************************************************
// DupSnoopMain.cpp - collects command-line parameters for MosaikDupSnoop.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include <iostream>
#include <vector>
#include <cstdio>
#include "AlignmentReader.h"
#include "AlignmentStatus.h"
#include "Benchmark.h"
#include "ConsoleUtilities.h"
#include "FileUtilities.h"
#include "Mosaik.h"
#include "MosaikDupSnoop.h"
#include "Options.h"
#include "SequencingTechnologies.h"

using namespace std;

// create a configuration variable struct
struct ConfigurationSettings {

	// flags
	bool AllowAllUniqueFragmentLengths;
	bool HasConfidenceInterval;
	bool HasInputFilename;
	bool HasOutputDirectory;
	bool IgnoreUM;
	bool IgnoreUU;
	bool IgnoreUniqueOrphans;
	bool ResolveMM;

	// filenames
	string OutputDirectory;
	vector<string> InputFiles;

	// parameters
	double ConfidenceInterval;

	// constructor
	ConfigurationSettings()
		: AllowAllUniqueFragmentLengths(false)
		, HasConfidenceInterval(false)
		, HasInputFilename(false)
		, HasOutputDirectory(false)
		, IgnoreUM(false)
		, IgnoreUU(false)
		, IgnoreUniqueOrphans(false)
		, ResolveMM(false)
		, ConfidenceInterval(DEFAULT_CONFIDENCE_INTERVAL)
	{}
};

int main(int argc, char* argv[]) {

	CConsole::Initialize();
	ConfigurationSettings settings;

	printf("------------------------------------------------------------------------------\n");
	printf("Mosaik"); CConsole::Red(); printf("DupSnoop"); CConsole::Reset();
	printf(" %u.%u.%04u                                             %s\n", 
		MOSAIK_MAJOR_VERSION, MOSAIK_MINOR_VERSION, MOSAIK_BUILD_VERSION, MOSAIK_VERSION_DATE);
	printf("Michael Stromberg                 Marth Lab, Boston College Biology Department\n");
	printf("------------------------------------------------------------------------------\n\n");

	// =================================
	// configure the command line parser
	// =================================

	// set general info about the program
	COptions::SetProgramInfo("MosaikDupSnoop", "stores duplicate fragment data in library-specific databases", "-in <filename> -od <directory>");

	// add the input/output options
	OptionGroup* pIoOpts = COptions::CreateOptionGroup("Input/output: (required)");
	COptions::AddValueOption("-in", "filename or directory", "the input MOSAIK alignment filename", "An unsorted MOSAIK alignment filename", settings.HasInputFilename,   settings.InputFiles,      pIoOpts);
	COptions::AddValueOption("-od", "directory",             "the output directory",                "An output directory",                   settings.HasOutputDirectory, settings.OutputDirectory, pIoOpts);

	// add the paired-end options
	OptionGroup* pPairedEndOpts = COptions::CreateOptionGroup("Paired-end Options");
	COptions::AddOption("-afl",                             "allows all fragment lengths when resolving unique vs unique read pairs", settings.AllowAllUniqueFragmentLengths,                      pPairedEndOpts);
	COptions::AddValueOption("-ci", "confidence interval",  "sets the desired fragment length confidence interval",  "",              settings.HasConfidenceInterval, settings.ConfidenceInterval, pPairedEndOpts, DEFAULT_CONFIDENCE_INTERVAL);
	COptions::AddOption("-iuo",                             "ignore unique orphaned reads",                                           settings.IgnoreUniqueOrphans,                                pPairedEndOpts);
	COptions::AddOption("-iuu",                             "ignore unique vs unique read pairs",                                     settings.IgnoreUU,                                           pPairedEndOpts);
	COptions::AddOption("-ium",                             "ignore unique vs multiple read pairs",                                   settings.IgnoreUM,                                           pPairedEndOpts);
	COptions::AddOption("-rmm",                             "resolve multiple vs multiple read pairs",                                settings.ResolveMM,                                          pPairedEndOpts);

	// parse the current command line
	COptions::Parse(argc, argv);

	// =============================
	// check for missing information
	// =============================

	bool foundError = false;
	ostringstream errorBuilder;
	const string ERROR_SPACER(7, ' ');

	if(settings.InputFiles.empty()) {
		errorBuilder << ERROR_SPACER << "No unsorted MOSAIK alignment filenames were found. Please use the -in parameter." << endl;
		foundError = true;
	}

	// convert the number of standard deviations
	if(settings.HasConfidenceInterval && ((settings.ConfidenceInterval <= 0.0) || (settings.ConfidenceInterval > 1.0))) {
		errorBuilder << ERROR_SPACER << "The confidence interval should be between 0.0 and 1.0." << endl;
		foundError = true;
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

	// check that the output directory exists
	CFileUtilities::CheckDirectory(settings.OutputDirectory);

	CMosaikDupSnoop ds;

	// display which types are being resolved
	printf("- resolving the following types of read pairs: ");
	if(!settings.IgnoreUniqueOrphans) printf("[unique orphans] ");
	if(!settings.IgnoreUU)            printf("[unique vs unique] ");
	if(!settings.IgnoreUM)            printf("[unique vs multiple] ");
	if(settings.ResolveMM)            printf("[multiple vs multiple] ");
	printf("\n");

	ds.ConfigureResolution(!settings.IgnoreUniqueOrphans, !settings.IgnoreUU, !settings.IgnoreUM, settings.ResolveMM);

	// enforce unique read pair constraints
	if(settings.AllowAllUniqueFragmentLengths) {
		printf("- allowing all fragment lengths for reads with unique mates\n");
		ds.EnableAllUniqueFragmentLengths();
	}

	// set the desired confidence interval
	if(settings.HasConfidenceInterval) {
		printf("- setting the confidence interval to %0.4f\n", settings.ConfidenceInterval);
		ds.SetConfidenceInterval(settings.ConfidenceInterval);
	}

	// add a directory separator to the output directory
	if(settings.OutputDirectory[settings.OutputDirectory.size() - 1] != OS_DIRECTORY_SEPARATOR) {
		settings.OutputDirectory += OS_DIRECTORY_SEPARATOR;
	}

	printf("\n");

	// ================
	// Record fragments
	// ================

	// start benchmarking
	CBenchmark bench;
	bench.Start();

	// records the fragments found in the input files
	ds.RecordFragments(settings.InputFiles, settings.OutputDirectory);

	// stop benchmarking
	bench.Stop();
	bench.DisplayTime("MosaikDupSnoop");

	return 0;
}
