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
	bool HasInputFastqFilename;
	bool HasInputFastq2Filename;
	bool HasSortingOrder;
	bool IsQuietMode;

	// filenames
	string AxtFilename;
	string BamFilename;
	string BedFilename;
	string ElandFilename;
	string FastqFilename;
	string InputAlignmentsFilename;
	string InputReadsFilename;
	string SamFilename;
	string InputFastqFilename;
	string InputFastq2Filename;

	// parameters
	string FilteredReferenceName;
	string SortingOrder;
	unsigned short SortingModel; //0: coordinate; 1: queryname; 2: unsorted


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
		, HasInputFastqFilename(false)
		, HasInputFastq2Filename(false)
		, HasSortingOrder(false)
		, IsQuietMode(false)
	{}
};

int main(int argc, char* argv[]) {

	CConsole::Initialize();
	ConfigurationSettings settings;

	printf("------------------------------------------------------------------------------\n");
	printf("Mosaik"); CConsole::Red(); printf("Text"); CConsole::Reset();
	printf(" %u.%u.%04u                                                 %s\n", 
		MOSAIK_MAJOR_VERSION, MOSAIK_MINOR_VERSION, MOSAIK_BUILD_VERSION, MOSAIK_VERSION_DATE);
	printf("Michael Stromberg & Wan-Ping Lee  Marth Lab, Boston College Biology Department\n");
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
	COptions::AddValueOption("-q",     "FASTQ filename",            "the original FASTQ file",                  "", settings.HasInputFastqFilename,      settings.InputFastqFilename,      pAlignmentArchiveOpts);
	COptions::AddValueOption("-q2",    "FASTQ filename",            "the original 2nd mate FASTQ file",         "", settings.HasInputFastq2Filename,     settings.InputFastq2Filename,     pAlignmentArchiveOpts);
	COptions::AddValueOption("-axt",   "axt filename",              "stores the data in an AXT file",           "", settings.HasAxtFilename,             settings.AxtFilename,             pAlignmentArchiveOpts);
	COptions::AddValueOption("-bam",   "bam filename",              "stores the data in a BAM file",            "", settings.HasBamFilename,             settings.BamFilename,             pAlignmentArchiveOpts);
	COptions::AddValueOption("-bed",   "bed filename",              "stores the data in a BED file",            "", settings.HasBedFilename,             settings.BedFilename,             pAlignmentArchiveOpts);
	COptions::AddValueOption("-eland", "eland filename",            "stores the data in an Eland file",         "", settings.HasElandFilename,           settings.ElandFilename,           pAlignmentArchiveOpts);
	COptions::AddValueOption("-ref",   "reference sequence name",   "displays output for a specific reference", "", settings.UseReferenceFilter,         settings.FilteredReferenceName,   pAlignmentArchiveOpts);
	COptions::AddValueOption("-sam",   "sam filename",              "stores the data in a SAM file",            "", settings.HasSamFilename,             settings.SamFilename,             pAlignmentArchiveOpts);
	COptions::AddOption("-screen",                                  "displays the alignments on the screen",        settings.EnableScreenOutput,                                           pAlignmentArchiveOpts);
	COptions::AddOption("-u",                                       "limit output to unique reads",                 settings.EvaluateUniqueReadsOnly,                                      pAlignmentArchiveOpts);
	COptions::AddValueOption("-sort",  "order",                     "sorting order: [queryname, coordinate(default)]", "", settings.HasSortingOrder,     settings.SortingOrder,            pAlignmentArchiveOpts);

        // add the interface options
	OptionGroup* pInterface = COptions::CreateOptionGroup("Interface Options");
	COptions::AddOption("-quiet",  "enable progress bars and counters", settings.IsQuietMode, pInterface);

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

	CSequenceUtilities::LowercaseSequence(settings.SortingOrder);
	if(settings.HasSortingOrder) {
		if ( ( settings.SortingOrder != "queryname" ) && ( settings.SortingOrder != "coordinate" ) ) {
			errorBuilder << ERROR_SPACER << "Unknown sorting order. Please choose between \"queryname\" or \"coordinate\". The default sorting order is \"coordinate\"." << endl;
			foundError = true;
		}
	}

	if(settings.HasSortingOrder && !settings.HasInputAlignmentsFilename ) {
		errorBuilder << ERROR_SPACER << "A MOSAIK archive was not specified. Please use the -in parameter." << endl;
		foundError = true;
	}
        
	if(settings.HasInputAlignmentsFilename) {
		// test to see if the specified input files exist
		SequencingTechnologies seqTech;
        	AlignmentStatus alignmentStatus;
	        MosaikReadFormat::CAlignmentReader::CheckFile(settings.InputAlignmentsFilename, seqTech, alignmentStatus, true);
		
		// activate paired-end sorting mode
	        bool isAlignmentArchivePairedEnd = false;
        	if((alignmentStatus & AS_PAIRED_END_READ) != 0) isAlignmentArchivePairedEnd = true;

		bool noneFastqFiles = !settings.HasInputFastqFilename && !settings.HasInputFastq2Filename;
		bool bothFastqFiles = settings.HasInputFastqFilename && settings.HasInputFastq2Filename;
		if ( isAlignmentArchivePairedEnd ) {
			if ( !noneFastqFiles && !bothFastqFiles ) {
				errorBuilder << ERROR_SPACER << "For paired-end data, both two original FASTQ files are needed, -q and -q2." << endl;
				foundError = true;
			}
		}
		else {
			if ( settings.HasInputFastq2Filename ){
				errorBuilder << ERROR_SPACER << "For single-end data, please use -q." << endl;
				foundError = true;
			}
		}

		if ( ( seqTech == ST_SOLID ) && !noneFastqFiles ) {
			errorBuilder << ERROR_SPACER << "Soft clipping does not support for SOLiD technology." << endl;
			foundError = true;
		}

		if ( ( settings.HasBedFilename || settings.HasElandFilename ) && !noneFastqFiles ) {
			errorBuilder << ERROR_SPACER << "Soft clipping only supports for SAM/BAM file formats." << endl;
			foundError = true;
		}

		MosaikReadFormat::CAlignmentReader reader;
		reader.Open( settings.InputAlignmentsFilename );
		char* signature;
		reader.GetSignature( signature );
		reader.Close();
		bool newSortedArchive = ( strcmp( signature, SORT_SIGNATURE ) == 0 ) ? true : false;
		bool isSorted         = ( ( alignmentStatus & AS_SORTED_ALIGNMENT ) != 0 ) ? true : false;

		if ( settings.HasSortingOrder && ( !newSortedArchive || !isSorted ) ) {
			errorBuilder << ERROR_SPACER << "ERROR: The input archive should be processed by MosaikSort when sorting alignments by positions or querynames." << endl;
			foundError = true;
		}

		if ( !newSortedArchive && !noneFastqFiles ) {
			errorBuilder << ERROR_SPACER << "ERROR: The input archive should be processed by MosaikSort when reporting soft-clipped alignments." << endl;
			foundError = true;
		}

		// must be isSorted if newSortedArchive 
		if ( newSortedArchive ) {
			if ( settings.HasSortingOrder && ( settings.SortingOrder == "queryname" ) ) {
				settings.SortingModel = 1; // sort by queryname
				settings.SortingOrder = "queryname";
			} else {
				settings.SortingModel = 0; // sort by coordinate
				settings.SortingOrder = "coordinate";
			}
		} else {
			if ( isSorted ) {
				settings.SortingModel = 0; // sort by coordinate
				settings.SortingOrder = "coordinate";
			} else
				settings.SortingModel = 2; // unsorted
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
	if(settings.IsQuietMode)         mt.SetQuietMode();

	if(settings.HasInputAlignmentsFilename) {

		mt.SetSortingOrder( settings.SortingModel );
		if ( ( settings.SortingModel == 0 ) || ( settings.SortingModel == 1 ))
			printf("- sorting alignments by %s\n", settings.SortingOrder.c_str() );

		if(settings.EvaluateUniqueReadsOnly) mt.EvaluateUniqueReadsOnly();
		if(settings.HasAxtFilename)          mt.EnableAxtOutput(settings.AxtFilename);
		if(settings.HasBamFilename)          mt.EnableBamOutput(settings.BamFilename);
		if(settings.HasBedFilename)          mt.EnableBedOutput(settings.BedFilename);
		if(settings.HasElandFilename)        mt.EnableElandOutput(settings.ElandFilename);
		if(settings.HasSamFilename)          mt.EnableSamOutput(settings.SamFilename);
		if(settings.UseReferenceFilter)      mt.EnableReferenceFilter(settings.FilteredReferenceName, settings.InputAlignmentsFilename);

		if(settings.HasInputFastqFilename) {
			mt.ParseFastqFile(settings.InputFastqFilename);
			printf("- Reporting soft-clipped alignments.\n");
		}
		if(settings.HasInputFastq2Filename)  mt.ParseFastq2File(settings.InputFastq2Filename);

		if(!settings.EnableScreenOutput) {
			printf("- converting the alignment archive to the following formats:");
			if(settings.HasAxtFilename)          printf(" AXT");
			if(settings.HasBamFilename)          printf(" BAM");
			if(settings.HasBedFilename)          printf(" BED");
			if(settings.HasElandFilename)        printf(" Eland");
			if(settings.HasSamFilename)          printf(" SAM");
			printf("\n\n");
		}

		mt.SetArchiveSetting(settings.InputAlignmentsFilename);
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
