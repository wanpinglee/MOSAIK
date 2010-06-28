// ***************************************************************************
// MergeMain.cpp - collects command-line parameters for MosaikMerge.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include <cstdio>
#include <iostream>
#include <vector>
#include "AlignmentReader.h"
#include "AlignmentStatus.h"
#include "Benchmark.h"
#include "ConsoleUtilities.h"
#include "FileUtilities.h"
#include "MosaikMerge.h"
#include "Options.h"

using namespace std;

// define our default values
unsigned int DEFAULT_CACHE_SIZE = 6000000;

// create a configuration variable struct
struct ConfigurationSettings {

	// flags
	bool HasCacheSize;
	bool HasInputFilename;
	bool HasOutputFilename;

	// filenames
	string OutputFilename;
	vector<string> InputFiles;

	// parameters
	unsigned int CacheSize;

	// constructor
	ConfigurationSettings()
		: HasCacheSize(false)
		, HasInputFilename(false)
		, HasOutputFilename(false)
		, CacheSize(DEFAULT_CACHE_SIZE)
	{}
};

int main(int argc, char* argv[]) {

	CConsole::Initialize();
	ConfigurationSettings settings;

	printf("------------------------------------------------------------------------------\n");
	printf("Mosaik"); CConsole::Red(); printf("Merge"); CConsole::Reset();
	printf(" %u.%u.%04u                                                %s\n", 
		MOSAIK_MAJOR_VERSION, MOSAIK_MINOR_VERSION, MOSAIK_BUILD_VERSION, MOSAIK_VERSION_DATE);
	printf("Michael Stromberg                 Marth Lab, Boston College Biology Department\n");
	printf("------------------------------------------------------------------------------\n\n");

	// =================================
	// configure the command line parser
	// =================================

	// set general info about the program
	COptions::SetProgramInfo("MosaikMerge", "merges several sorted MOSAIK alignment files", "-in <filename> [-in <filename2>] -out <filename>");

	// add the input/output options
	OptionGroup* pOpts = COptions::CreateOptionGroup("Options");
	COptions::AddValueOption("-in",  "filename|directory",        "any number of MOSAIK alignment files or directories", "A sorted input MOSAIK alignment filename", settings.HasInputFilename,  settings.InputFiles,     pOpts);
	COptions::AddValueOption("-out", "filename",                  "the output MOSAIK alignment filename",                "An output MOSAIK alignment filename",      settings.HasOutputFilename, settings.OutputFilename, pOpts);
	COptions::AddValueOption("-mem", "# of alignments in memory", "sets the sorting cache size",                         "",                                         settings.HasCacheSize,      settings.CacheSize,      pOpts, DEFAULT_CACHE_SIZE);

	// parse the current command line
	COptions::Parse(argc, argv);

	// =============================
	// check for missing information
	// =============================

	bool foundError = false;
	ostringstream errorBuilder;
	const string ERROR_SPACER(7, ' ');

	// check if we have any files to process
	if(settings.InputFiles.empty()) {
		errorBuilder << ERROR_SPACER << "None of the files seemed to be sorted alignment archives. Please use MosaikSort after aligning with MosaikAligner." << endl;
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
				MosaikReadFormat::CAlignmentReader::CheckFile(dirFiles[i], seqTech, alignmentStatus, false);
				if((alignmentStatus & AS_SORTED_ALIGNMENT) != 0) expandedFileList.push_back(dirFiles[i]);
			}

		} else expandedFileList.push_back(*sIter);
	}

	// ===============
	// Merge the files
	// ===============

	// start benchmarking
	CBenchmark bench;
	bench.Start();

	CMosaikMerge mm(settings.CacheSize);
	mm.MergeFiles(expandedFileList, settings.OutputFilename);

	// ==================
	// Show total runtime
	// ==================

	// stop benchmarking
	bench.Stop();
	bench.DisplayTime("\nMosaikMerge");

	return 0;
}
