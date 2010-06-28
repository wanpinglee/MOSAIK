// ***************************************************************************
// TextMain.cpp - collects command-line parameters for MosaikText.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include <iomanip>
#include <iostream>
#include "AlignmentReader.h"
#include "AlignmentStatus.h"
#include "Benchmark.h"
#include "ConsoleUtilities.h"
#include "FileUtilities.h"
#include "Mosaik.h"
#include "MosaikText.h"
#include "Options.h"
#include "ReadReader.h"
#include "ReadStatus.h"
#include "SequencingTechnologies.h"

using namespace std;

// define our default values

// create a configuration variable struct
struct ConfigurationSettings {

	// flags
	bool EnableScreenOutput;
	bool HasAxtFilename;
	bool HasBamFilename;
	bool HasBedFilename;
	bool HasElandFilename;
	bool HasFastqFilename;
	bool HasInputAlignmentsFilename;
	bool HasInputReadsFilename;
	bool HasSamFilename;
	bool EvaluateUniqueReadsOnly;
	bool UseReferenceFilter;

	// filenames
	string AxtFilename;
	string BamFilename;
	string BedFilename;
	string ElandFilename;
	string FastqFilename;
	string InputAlignmentsFilename;
	string InputReadsFilename;
	string SamFilename;

	// parameters
	string FilteredReferenceName;

	// constructor
	ConfigurationSettings()
		: EnableScreenOutput(false)
		, HasAxtFilename(false)
		, HasBamFilename(false)
		, HasBedFilename(false)
		, HasElandFilename(false)
		, HasFastqFilename(false)
		, HasInputAlignmentsFilename(false)
		, HasInputReadsFilename(false)
		, HasSamFilename(false)
		, EvaluateUniqueReadsOnly(false)
		, UseReferenceFilter(false)
	{}
};

int main(int argc, char* argv[]) {

	CConsole::Initialize();
	ConfigurationSettings settings;

	printf("------------------------------------------------------------------------------\n");
	printf("Mosaik"); CConsole::Red(); printf("Text"); CConsole::Reset();
	printf(" %u.%u.%04u                                                 %s\n", 
		MOSAIK_MAJOR_VERSION, MOSAIK_MINOR_VERSION, MOSAIK_BUILD_VERSION, MOSAIK_VERSION_DATE);
	printf("Michael Stromberg                 Marth Lab, Boston College Biology Department\n");
	printf("------------------------------------------------------------------------------\n\n");

	// =================================
	// configure the command line parser
	// =================================

	// set general info about the program
	COptions::SetProgramInfo("MosaikText", "exports reads and alignments from MOSAIK", "-in <filename> -out <filename> -ia <filename>");

	// add the read archive options
	OptionGroup* pReadArchiveOpts = COptions::CreateOptionGroup("Read Archive Options");
	COptions::AddValueOption("-ir",    "MOSAIK read filename", "the input read file",             "", settings.HasInputReadsFilename, settings.InputReadsFilename, pReadArchiveOpts);
	COptions::AddValueOption("-fastq", "FASTQ filename",       "stores the data in a FASTQ file", "", settings.HasFastqFilename,      settings.FastqFilename,      pReadArchiveOpts);
	COptions::AddOption("-screen",                             "displays the reads on the screen",    settings.EnableScreenOutput,                                 pReadArchiveOpts);

	// add the alignment archive options
	OptionGroup* pAlignmentArchiveOpts = COptions::CreateOptionGroup("Alignment Archive Options");
	COptions::AddValueOption("-in",    "MOSAIK alignment filename", "the input alignment file",                 "", settings.HasInputAlignmentsFilename, settings.InputAlignmentsFilename, pAlignmentArchiveOpts);
	COptions::AddValueOption("-axt",   "axt filename",              "stores the data in an AXT file",           "", settings.HasAxtFilename,             settings.AxtFilename,             pAlignmentArchiveOpts);
	COptions::AddValueOption("-bam",   "bam filename",              "stores the data in a BAM file",            "", settings.HasBamFilename,             settings.BamFilename,             pAlignmentArchiveOpts);
	COptions::AddValueOption("-bed",   "bed filename",              "stores the data in a BED file",            "", settings.HasBedFilename,             settings.BedFilename,             pAlignmentArchiveOpts);
	COptions::AddValueOption("-eland", "eland filename",            "stores the data in an Eland file",         "", settings.HasElandFilename,           settings.ElandFilename,           pAlignmentArchiveOpts);
	COptions::AddValueOption("-ref",   "reference sequence name",   "displays output for a specific reference", "", settings.UseReferenceFilter,         settings.FilteredReferenceName,   pAlignmentArchiveOpts);
	COptions::AddValueOption("-sam",   "sam filename",              "stores the data in a SAM file",            "", settings.HasSamFilename,             settings.SamFilename,             pAlignmentArchiveOpts);
	COptions::AddOption("-screen",                                  "displays the alignments on the screen",        settings.EnableScreenOutput,                                           pAlignmentArchiveOpts);
	COptions::AddOption("-u",                                       "limit output to unique reads",                 settings.EvaluateUniqueReadsOnly,                                      pAlignmentArchiveOpts);

	// parse the current command line
	COptions::Parse(argc, argv);

	// =============================
	// check for missing information
	// =============================

	bool foundError = false;
	ostringstream errorBuilder;
	const string ERROR_SPACER(7, ' ');

	if(!settings.HasInputAlignmentsFilename && !settings.HasInputReadsFilename) {
		errorBuilder << ERROR_SPACER << "A MOSAIK input file was not specified. Please use either the -in parameter (for alignments) or the -ir parameter (for reads)." << endl;
		foundError = true;
	}

	if(settings.HasInputAlignmentsFilename && settings.HasInputReadsFilename) {
		errorBuilder << ERROR_SPACER << "Both MOSAIK input options were specified. Please use either the -in parameter (for alignments) or the -ir parameter (for reads)." << endl;
		foundError = true;
	}

	if(settings.HasInputAlignmentsFilename) {
		if(!settings.EnableScreenOutput && !settings.HasAxtFilename && !settings.HasBamFilename && !settings.HasBedFilename && !settings.HasElandFilename && !settings.HasSamFilename) {
			errorBuilder << ERROR_SPACER << "Please specify an output file format. AXT (-axt), BAM (-bam), BED (-bed), Eland (-eland), SAM (-sam), or screen (-screen)." << endl;
			foundError = true;
		}
	} else {
		if(!settings.EnableScreenOutput && !settings.HasFastqFilename) {
			errorBuilder << ERROR_SPACER << "Please specify an output file format. FASTQ (-fastq) or screen (-screen)." << endl;
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

	// test to see if the specified input files exist
	SequencingTechnologies seqTech;
	AlignmentStatus alignmentStatus;
	ReadStatus readStatus;

	if(settings.HasInputAlignmentsFilename) MosaikReadFormat::CAlignmentReader::CheckFile(settings.InputAlignmentsFilename, seqTech, alignmentStatus, true);
	if(settings.HasInputReadsFilename)      MosaikReadFormat::CReadReader::CheckFile(settings.InputReadsFilename, seqTech, readStatus, true);

	// start benchmarking
	CBenchmark bench;
	bench.Start();

	CMosaikText mt;
	if(settings.EnableScreenOutput) mt.EnableScreenOutput();

	if(settings.HasInputAlignmentsFilename) {

		if(settings.EvaluateUniqueReadsOnly) mt.EvaluateUniqueReadsOnly();
		if(settings.HasAxtFilename)          mt.EnableAxtOutput(settings.AxtFilename);
		if(settings.HasBamFilename)          mt.EnableBamOutput(settings.BamFilename);
		if(settings.HasBedFilename)          mt.EnableBedOutput(settings.BedFilename);
		if(settings.HasElandFilename)        mt.EnableElandOutput(settings.ElandFilename);
		if(settings.HasSamFilename)          mt.EnableSamOutput(settings.SamFilename);
		if(settings.UseReferenceFilter)      mt.EnableReferenceFilter(settings.FilteredReferenceName, settings.InputAlignmentsFilename);

		if(!settings.EnableScreenOutput) {
			printf("- converting the alignment archive to the following formats:");
			if(settings.HasAxtFilename)          printf(" AXT");
			if(settings.HasBamFilename)          printf(" BAM");
			if(settings.HasBedFilename)          printf(" BED");
			if(settings.HasElandFilename)        printf(" Eland");
			if(settings.HasSamFilename)          printf(" SAM");
			printf("\n\n");
		}

		mt.ParseMosaikAlignmentFile(settings.InputAlignmentsFilename);

	} else {

		if(settings.HasFastqFilename) mt.EnableFastqOutput(settings.FastqFilename);

		mt.ParseMosaikReadFile(settings.InputReadsFilename);
	}

	// ==================
	// Show total runtime
	// ==================

	// stop benchmarking
	bench.Stop();

	// show the benchmarking results
	bench.DisplayTime("\nMosaikText");

	return 0;
}
