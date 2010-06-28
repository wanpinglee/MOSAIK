// ***************************************************************************
// CoverageMain.cpp - collects command-line parameters for MosaikCoverage.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include <iomanip>
#include <iostream>
#include "AlignmentStatus.h"
#include "Benchmark.h"
#include "ConsoleUtilities.h"
#include "FileUtilities.h"
#include "Mosaik.h"
#include "MosaikCoverage.h"
#include "Options.h"
#include "SequencingTechnologies.h"

using namespace std;

unsigned char MIN_COVERAGE = 1;

// create a configuration variable struct
struct ConfigurationSettings {

	// flags
	bool CreatePostscriptGraphs;
	bool HasInputAlignmentsFilename;
	bool HasInputReferencesFilename;
        bool HasMinAlignmentQuality;
        bool HasMinCoverage;
	bool HasOutputDirectory;
	bool EvaluateUniqueReadsOnly;
	bool SkipOutput;

	// filenames
	vector<string> InputFiles;
	string InputReferencesFilename;
	string OutputDirectory;

        // parameters
        unsigned char MinAlignmentQuality;
        unsigned char MinCoverage;

	// constructor
	ConfigurationSettings()
		: CreatePostscriptGraphs(false)
		, HasInputAlignmentsFilename(false)
		, HasInputReferencesFilename(false)
                , HasMinAlignmentQuality(false)
		, HasMinCoverage(false)
		, HasOutputDirectory(false)
		, EvaluateUniqueReadsOnly(false)
		, SkipOutput(false)
		, MinCoverage(MIN_COVERAGE)
	{}
};

int main(int argc, char* argv[]) {

	CConsole::Initialize();
	ConfigurationSettings settings;

	printf("------------------------------------------------------------------------------\n");
	printf("Mosaik"); CConsole::Red(); printf("Coverage"); CConsole::Reset();
	printf(" %u.%u.%04u                                             %s\n", 
		MOSAIK_MAJOR_VERSION, MOSAIK_MINOR_VERSION, MOSAIK_BUILD_VERSION, MOSAIK_VERSION_DATE);
	printf("Michael Stromberg                 Marth Lab, Boston College Biology Department\n");
	printf("------------------------------------------------------------------------------\n\n");

	// =================================
	// configure the command line parser
	// =================================

	// set general info about the program
	COptions::SetProgramInfo("MosaikCoverage", "produces coverage files given a MOSAIK alignment file", "-in <filename> [-in <filename2>] -ia <filename>");

	// add the input options
	OptionGroup* pInputOpts = COptions::CreateOptionGroup("Input Options");
	COptions::AddValueOption("-ia",  "filename",              "the input reference file",        "An input MOSAIK reference file", settings.HasInputReferencesFilename, settings.InputReferencesFilename, pInputOpts);
	COptions::AddValueOption("-in",  "filename or directory", "the alignment file or directory", "An input MOSAIK alignment file", settings.HasInputAlignmentsFilename, settings.InputFiles,              pInputOpts);
	COptions::AddOption("-u",                                 "limit coverage analysis to unique reads",                           settings.EvaluateUniqueReadsOnly,                                      pInputOpts);

	// add the output options
	OptionGroup* pOutputOpts = COptions::CreateOptionGroup("Output Options");
	COptions::AddValueOption("-od", "directory",         "the output directory", "",                       settings.HasOutputDirectory,     settings.OutputDirectory,     pOutputOpts);
	COptions::AddValueOption("-aq", "alignment quality", "sets the minimum alignment quality", "",         settings.HasMinAlignmentQuality, settings.MinAlignmentQuality, pOutputOpts);
	COptions::AddOption("-cg",                           "creates coverage graphs using gnuplot",          settings.CreatePostscriptGraphs,                               pOutputOpts);
	COptions::AddValueOption("-mc", "coverage",          "coverage below this value will be set to 0", "", settings.HasMinCoverage,         settings.MinCoverage,         pOutputOpts, MIN_COVERAGE);
	COptions::AddOption("-ngc",                          "skips making graphs and coverage files",         settings.SkipOutput,                                           pOutputOpts);
	
	
	// parse the current command line
	COptions::Parse(argc, argv);

	// ================================================================
	// Expand the input file list with sorted MOSAIK alignment archives
	// ================================================================

	vector<string> expandedFileList;
	vector<string>::const_iterator sIter;
	AlignmentStatus alignmentStatus;
	SequencingTechnologies seqTech;

	for(sIter = settings.InputFiles.begin(); sIter != settings.InputFiles.end(); ++sIter) {
		if(CFileUtilities::DirExists(sIter->c_str())) {

			vector<string> dirFiles;
			CFileUtilities::SearchDirectory(dirFiles, sIter->c_str());

			for(unsigned int i = 0; i < (unsigned int)dirFiles.size(); i++) {
				if(MosaikReadFormat::CAlignmentReader::CheckFile(dirFiles[i], seqTech, alignmentStatus, false)) 
					expandedFileList.push_back(dirFiles[i]);
			}

		} else expandedFileList.push_back(*sIter);
	}

	// ===================================================
	// Parse configuration strings and set class variables
	// ===================================================

	// test to see if the specified input files exist
	MosaikReadFormat::CReferenceSequenceReader::CheckFile(settings.InputReferencesFilename, true);

	if(settings.HasOutputDirectory) {
		CFileUtilities::CheckDirectory(settings.OutputDirectory.c_str());

		if(settings.OutputDirectory.size() > 1)
			if(settings.OutputDirectory[settings.OutputDirectory.size() - 1] == OS_DIRECTORY_SEPARATOR)
				settings.OutputDirectory = settings.OutputDirectory.substr(0, settings.OutputDirectory.size() - 1);
	} else settings.OutputDirectory = "";

	// start benchmarking
	CBenchmark bench;
	bench.Start();

	CMosaikCoverage mc;

	if(settings.EvaluateUniqueReadsOnly) mc.EvaluateUniqueReadsOnly();
	if(settings.CreatePostscriptGraphs)  mc.EnableGraphCreation();
	if(settings.SkipOutput)              mc.DisableOutput();
	if(settings.HasMinAlignmentQuality)  mc.SetMinAlignmentQuality(settings.MinAlignmentQuality);
	printf("- statistics will use %ux coverage as the threshold.\n", settings.MinCoverage);
	mc.SetOutputDirectory(settings.OutputDirectory);
	mc.ParseMosaikAlignmentFile(expandedFileList, settings.InputReferencesFilename);
	mc.SaveCoverage(settings.MinCoverage);

	// ==================
	// Show total runtime
	// ==================

	// stop benchmarking
	bench.Stop();

	// show the benchmarking results
	bench.DisplayTime("\nMosaikCoverage");

	return 0;
}
