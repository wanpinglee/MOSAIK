// ***************************************************************************
// AssemblerMain.cpp - collects command-line parameters for MosaikAssembler.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include <iostream>
#include "AlignmentReader.h"
#include "AlignmentStatus.h"
#include "Benchmark.h"
#include "ConsoleUtilities.h"
#include "FileUtilities.h"
#include "MosaikAssembler.h"
#include "Options.h"
#include "ReferenceSequenceReader.h"
#include "RegexUtilities.h"
#include "SequencingTechnologies.h"

using namespace std;

// define our default values
const CMosaikAssembler::AssemblyFormatType DEFAULT_ASSEMBLY_FORMAT = CMosaikAssembler::AssemblyFormat_ACE;
string DEFAULT_ASSEMBLY_FORMAT_STRING                              = "ace";
unsigned char DEFAULT_REFERENCE_SEQUENCE_BASE_QUALITY              = 40;

// create a configuration variable struct
struct ConfigurationSettings {

	// flags
	bool EnableFastaFileCreation;
	bool HasAlignmentFilename;
	bool HasAssemblyFormat;
	bool HasReferenceFilename;
	bool HasOutputFilenameStub;
	bool HasReferenceBaseQuality;
	bool HasRegionOfInterest;

	// filenames
	string AlignmentFilename;
	string ReferenceFilename;
	string OutputFilenameStub;

	// parameters
	string AssemblyFormat;
	unsigned char ReferenceBaseQuality;
	string RegionOfInterest;

	// constructor
	ConfigurationSettings()
		: EnableFastaFileCreation(false)
		, HasAlignmentFilename(false)
		, HasAssemblyFormat(false)
		, HasOutputFilenameStub(false)
		, HasReferenceBaseQuality(false)
		, HasRegionOfInterest(false)
		, ReferenceBaseQuality(DEFAULT_REFERENCE_SEQUENCE_BASE_QUALITY)
	{}
};



int main(int argc, char* argv[]) {

	CConsole::Initialize();
	ConfigurationSettings settings;

	printf("------------------------------------------------------------------------------\n");
	printf("Mosaik"); CConsole::Red(); printf("Assembler"); CConsole::Reset();
	printf(" %u.%u.%04u                                            %s\n", 
		MOSAIK_MAJOR_VERSION, MOSAIK_MINOR_VERSION, MOSAIK_BUILD_VERSION, MOSAIK_VERSION_DATE);
	printf("Michael Stromberg                 Marth Lab, Boston College Biology Department\n");
	printf("------------------------------------------------------------------------------\n\n");

	// =================================
	// configure the command line parser
	// =================================

	// set general info about the program
	COptions::SetProgramInfo("MosaikAssembler", "produces an assembly file from a MOSAIK alignment file", "-in <filename> -out <file stub> -ia <filename>");

	// add the input/output options
	OptionGroup* pIoOpts = COptions::CreateOptionGroup("Input/output: (required)");
	COptions::AddValueOption("-ia",  "MOSAIK reference filename",  "the input reference file",               "An MOSAIK reference file",       settings.HasReferenceFilename,  settings.ReferenceFilename,  pIoOpts);
	COptions::AddValueOption("-in",  "MOSAIK alignment filename",  "the sorted input MOSAIK alignment file", "A sorted MOSAIK alignment file", settings.HasAlignmentFilename,  settings.AlignmentFilename,  pIoOpts);
	COptions::AddValueOption("-out", "MOSAIK alignment filename",  "the assembly output filename stub",      "An output filename stub",        settings.HasOutputFilenameStub, settings.OutputFilenameStub, pIoOpts);

	// add the options
	OptionGroup* pOpts = COptions::CreateOptionGroup("Options");
	COptions::AddOption("-cf",                               "creates FASTA read files for gigaBayes",      settings.EnableFastaFileCreation,                                pOpts);
	COptions::AddValueOption("-f",   "assembly file format", "assembly format: 'ace' or 'gig'",         "", settings.HasAssemblyFormat,       settings.AssemblyFormat,       pOpts, DEFAULT_ASSEMBLY_FORMAT_STRING);
	COptions::AddValueOption("-rbq", "base quality",         "sets the default reference base quality", "", settings.HasReferenceBaseQuality, settings.ReferenceBaseQuality, pOpts);
	COptions::AddValueOption("-roi", "refseq name",          "assembles a specific reference sequence", "", settings.HasRegionOfInterest,     settings.RegionOfInterest,     pOpts);

	// parse the current command line
	COptions::Parse(argc, argv);

	// =============================
	// check for missing information
	// =============================

	bool foundError = false;
	ostringstream errorBuilder;
	const string ERROR_SPACER(7, ' ');

	// check the reference base quality
	if(settings.HasReferenceBaseQuality && (settings.ReferenceBaseQuality > 99)) {
		errorBuilder << ERROR_SPACER << "The reference sequence base quality should be between 0 and 99." << endl;
		foundError = true;
	}

	// convert the assembly format
	CMosaikAssembler::AssemblyFormatType assemblyFormat = DEFAULT_ASSEMBLY_FORMAT;
	if(settings.HasAssemblyFormat) {
		if(settings.AssemblyFormat == "ace") {
			assemblyFormat = CMosaikAssembler::AssemblyFormat_ACE;
		} else if(settings.AssemblyFormat == "gig") {
			assemblyFormat = CMosaikAssembler::AssemblyFormat_GIG;
		} else {
			errorBuilder << ERROR_SPACER << "Unknown assembly format (" << assemblyFormat << ") specified. Please choose between 'ace' or 'gig'. The default value is 'ace'." << endl;
			foundError = true;
		}
	}

	// test to see if the specified input files exist
	SequencingTechnologies seqTech;
	AlignmentStatus alignmentStatus;
	MosaikReadFormat::CAlignmentReader::CheckFile(settings.AlignmentFilename, seqTech, alignmentStatus, true);

	if((alignmentStatus & AS_SORTED_ALIGNMENT) == 0) {
		errorBuilder << ERROR_SPACER << "Unable to assemble this alignment file because it hasn't been sorted yet. Please use MosaikSort to remedy this." << endl;
		foundError = true;
	}

	// make sure the region of interest exists
	if(settings.HasRegionOfInterest) {
		MosaikReadFormat::CAlignmentReader reader;
		reader.Open(settings.AlignmentFilename);

		unsigned int regionOfInterestIndex = 0;
		if(!reader.GetReferenceSequenceIndex(settings.RegionOfInterest, regionOfInterestIndex)) {
			errorBuilder << ERROR_SPACER << "Unable to find the desired reference sequence in this alignment archive: [" << settings.RegionOfInterest << "]." << endl;
			foundError = true;
		}

		reader.Close();
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

	MosaikReadFormat::CReferenceSequenceReader::CheckFile(settings.ReferenceFilename, true);	

	// start benchmarking
	CBenchmark bench;
	bench.Start();

	// create our assembler
	CMosaikAssembler ma(assemblyFormat, settings.ReferenceBaseQuality);

	// ===============
	// enable features
	// ===============

	// select a specific region of interest
	if(settings.HasRegionOfInterest) ma.SelectRegionOfInterest(settings.RegionOfInterest);

	// enable FASTA file creation
	if(settings.EnableFastaFileCreation) ma.EnableFastaFileCreation();

	// ===========================
	// assemble the alignment file
	// ===========================

	ma.Assemble(settings.ReferenceFilename, settings.AlignmentFilename, settings.OutputFilenameStub);

	// ==================
	// Show total runtime
	// ==================

	// show the benchmarking results
	bench.Stop();
	bench.DisplayTime("\nMosaikAssembler");

	return 0;
}
