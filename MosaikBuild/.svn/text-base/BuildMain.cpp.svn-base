// ***************************************************************************
// BuildMain.cpp - collects command-line parameters for MosaikBuild.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include <iostream>
#include "Benchmark.h"
#include "ConsoleUtilities.h"
#include "FileUtilities.h"
#include "MosaikBuild.h"
#include "Options.h"
#include "ReadGroup.h"
#include "SequencingTechnologies.h"

using namespace std;

// create a configuration variable struct
struct ConfigurationSettings {

	// flags
	bool EnableColorspace;
	bool HasBaseQualityFasta2Filename;
	bool HasBaseQualityFastaFilename;
	bool HasBustardDirectory;
	bool HasCenterName;
	bool HasDescription;
	bool HasFastq2Filename;
	bool HasFastqFilename;
	bool HasGenomeAssemblyID;
	bool HasGeraldDirectory;
	bool HasIlluminaLanes;
	bool HasLibraryName;
	bool HasMedianFragmentLength;
	bool HasOutputReadsFilename;
	bool HasOutputReferenceFilename;
	bool HasPlatformUnit;
	bool HasReadFasta2Filename;
	bool HasReadFastaFilename;
	bool HasReadGroupID;
	bool HasReadLimit;
	bool HasReadNamePrefix;
	bool HasSampleName;
	bool HasSequenceTechnology;
	bool HasSpeciesName;
	bool HasSrfFilename;
	bool HasTrimPrefixBases;
	bool HasTrimPrefixName;
	bool HasTrimSuffixBases;
	bool HasTrimSuffixName;
	bool HasUniformResourceIdentifier;
	bool SetNumNBasesAllowed;
	bool SplitBustardReads;
	bool UseAssignedBQ;

	// filenames
	string BaseQualityFasta2Filename;
	string BaseQualityFastaFilename;
	string BustardDirectory;
	string Fastq2Filename;
	string FastqFilename;
	string GeraldDirectory;
	string OutputReadsFilename;
	string OutputReferenceFilename;
	string ReadFasta2Filename;
	string ReadFastaFilename;
	string SrfFilename;

	// parameters
	string CenterName;
	string Description;
	string GenomeAssemblyID;
	string IlluminaLanesString;
	string LibraryName;
	string PlatformUnit;
	string ReadGroupID;
	string ReadNamePrefix;
	string SampleName;
	string SequenceTechnologyString;
	string SpeciesName;
	string UniformResourceIdentifier;
	uint64_t ReadLimit;
	unsigned char AssignedBQ;
	unsigned int MedianFragmentLength;
	unsigned int NumNBasesAllowed;
	unsigned int NumTrimPrefixBases;
	unsigned int NumTrimPrefixName;
	unsigned int NumTrimSuffixBases;
	unsigned int NumTrimSuffixName;

	// constructor
	ConfigurationSettings()
		: EnableColorspace(false)
		, HasBaseQualityFasta2Filename(false)
		, HasBaseQualityFastaFilename(false)
		, HasBustardDirectory(false)
		, HasCenterName(false)
		, HasDescription(false)
		, HasFastq2Filename(false)
		, HasFastqFilename(false)
		, HasGenomeAssemblyID(false)
		, HasGeraldDirectory(false)
		, HasIlluminaLanes(false)
		, HasLibraryName(false)
		, HasMedianFragmentLength(false)
		, HasOutputReadsFilename(false)
		, HasOutputReferenceFilename(false)
		, HasPlatformUnit(false)
		, HasReadFasta2Filename(false)
		, HasReadFastaFilename(false)
		, HasReadGroupID(false)
		, HasReadLimit(false)
		, HasReadNamePrefix(false)
		, HasSampleName(false)
		, HasSequenceTechnology(false)
		, HasSpeciesName(false)
		, HasSrfFilename(false)
		, HasTrimPrefixBases(false)
		, HasTrimPrefixName(false)
		, HasTrimSuffixBases(false)
		, HasTrimSuffixName(false)
		, HasUniformResourceIdentifier(false)
		, SetNumNBasesAllowed(false)
		, SplitBustardReads(false)
		, UseAssignedBQ(false)
		, NumNBasesAllowed(4)
	{}
};

int main(int argc, char* argv[]) {

	CConsole::Initialize();
	ConfigurationSettings settings;

	printf("------------------------------------------------------------------------------\n");
	printf("Mosaik"); CConsole::Red(); printf("Build"); CConsole::Reset();
	printf(" %u.%u.%04u                                                %s\n", 
		MOSAIK_MAJOR_VERSION, MOSAIK_MINOR_VERSION, MOSAIK_BUILD_VERSION, MOSAIK_VERSION_DATE);
	printf("Michael Stromberg                 Marth Lab, Boston College Biology Department\n");
	printf("------------------------------------------------------------------------------\n\n");

	// =================================
	// configure the command line parser
	// =================================

	// set general info about the program
	COptions::SetProgramInfo("MosaikBuild", "converts external read formats to native MOSAIK formats", "[OPTIONS] [-out|-oa] <filename>");

	// add the reference sequence options
	OptionGroup* pRefOpts = COptions::CreateOptionGroup("Conversion (Reference Sequence)");
	COptions::AddOption("-cs",  "translate reference to colorspace", settings.EnableColorspace, pRefOpts);
	COptions::AddValueOption("-fr",  "FASTA reference filename",  "the FASTA reference sequences file",      "", settings.HasReadFastaFilename,         settings.ReadFastaFilename,         pRefOpts);
	COptions::AddValueOption("-ga",  "genome assembly ID",        "the genome assembly ID. e.g. HG18",       "", settings.HasGenomeAssemblyID,          settings.GenomeAssemblyID,          pRefOpts);
	COptions::AddValueOption("-oa",  "MOSAIK reference filename", "the output reference file",               "", settings.HasOutputReferenceFilename,   settings.OutputReferenceFilename,   pRefOpts);
	COptions::AddValueOption("-sn",  "species name",              "the species name. e.g. \"Homo sapiens\"", "", settings.HasSpeciesName,               settings.SpeciesName,               pRefOpts);
	COptions::AddValueOption("-uri", "uniform resource ID",       "the URI (e.g. URL or URN)",               "", settings.HasUniformResourceIdentifier, settings.UniformResourceIdentifier, pRefOpts);

	// add the FASTA options
	OptionGroup* pFastaOpts = COptions::CreateOptionGroup("Conversion (FASTA)");
	COptions::AddValueOption("-fr",         "FASTA read filename",    "the FASTA reads file",            "", settings.HasReadFastaFilename,         settings.ReadFastaFilename,         pFastaOpts);
	COptions::AddValueOption("-fq",         "FASTA quality filename", "the FASTA base qualities file",   "", settings.HasBaseQualityFastaFilename,  settings.BaseQualityFastaFilename,  pFastaOpts);
	COptions::AddValueOption("-fr2",        "FASTA read filename",    "the FASTA 2nd mate",              "", settings.HasReadFasta2Filename,        settings.ReadFasta2Filename,        pFastaOpts);
	COptions::AddValueOption("-fq2",        "FASTA quality filename", "the FASTA BQ 2nd mate",           "", settings.HasBaseQualityFasta2Filename, settings.BaseQualityFasta2Filename, pFastaOpts);
	COptions::AddValueOption("-assignQual", "base quality",           "assigns a quality for each base", "", settings.UseAssignedBQ,                settings.AssignedBQ,                pFastaOpts);

	// add the FASTQ options
	OptionGroup* pFastqOpts = COptions::CreateOptionGroup("Conversion (FASTQ)");
	COptions::AddValueOption("-q",  "FASTQ filename or directory", "the FASTQ file or directory", "", settings.HasFastqFilename,  settings.FastqFilename,   pFastqOpts);
	COptions::AddValueOption("-q2", "FASTQ filename or directory", "the FASTQ 2nd mate",          "", settings.HasFastq2Filename, settings.Fastq2Filename,  pFastqOpts);

	// add the SRF options
	OptionGroup* pSrfOpts = COptions::CreateOptionGroup("Conversion (Short Read Format)");
	COptions::AddValueOption("-srf", "SRF filename or directory", "the SRF file or directory", "", settings.HasSrfFilename,  settings.SrfFilename, pSrfOpts);

	// add the Bustard options
	OptionGroup* pBustardOpts = COptions::CreateOptionGroup("Conversion (Illumina Bustard)");
	COptions::AddValueOption("-bd", "Bustard directory", "the Illumina Bustard directory", "", settings.HasBustardDirectory, settings.BustardDirectory,    pBustardOpts);
	COptions::AddValueOption("-il", "lanes", "the desired lanes e.g 5678 for lanes 5-8",   "", settings.HasIlluminaLanes,    settings.IlluminaLanesString, pBustardOpts);
	COptions::AddOption("-split",  "splits the read into two mates", settings.SplitBustardReads, pBustardOpts);

	// add the Gerald options
	OptionGroup* pGeraldOpts = COptions::CreateOptionGroup("Conversion (Illumina Gerald)");
	COptions::AddValueOption("-gd", "Gerald directory", "the Illumina Gerald directory",            "", settings.HasGeraldDirectory, settings.GeraldDirectory,     pGeraldOpts);
	COptions::AddValueOption("-il", "lanes",            "the desired lanes e.g 5678 for lanes 5-8", "", settings.HasIlluminaLanes,   settings.IlluminaLanesString, pGeraldOpts);

	// add the read archive read group options
	OptionGroup* pMetadataOpts = COptions::CreateOptionGroup("Read Archive Metadata");
	COptions::AddValueOption("-cn", "center name",             "sequencing center name. e.g. broad",                                                "", settings.HasCenterName,           settings.CenterName,               pMetadataOpts);
	COptions::AddValueOption("-ds", "description",             "read group description",                                                            "", settings.HasDescription,          settings.Description,              pMetadataOpts);
	COptions::AddValueOption("-id", "identifier",              "read group ID. e.g. SRR009060",                                                     "", settings.HasReadGroupID,          settings.ReadGroupID,              pMetadataOpts);
	COptions::AddValueOption("-ln", "library name",            "library name. e.g. g1k-sc-NA18944-JPT-1",                                           "", settings.HasLibraryName,          settings.LibraryName,              pMetadataOpts);
	COptions::AddValueOption("-mfl", "median fragment length", "median fragment length. e.g. 150",                                                  "", settings.HasMedianFragmentLength, settings.MedianFragmentLength,     pMetadataOpts);
	COptions::AddValueOption("-pu", "run name & lane",         "the platform unit. e.g. IL12_490_5",                                                "", settings.HasPlatformUnit,         settings.PlatformUnit,             pMetadataOpts);
	COptions::AddValueOption("-sam", "sample name",            "sample name. e.g. NA12878",                                                         "", settings.HasSampleName,           settings.SampleName,               pMetadataOpts);
	COptions::AddValueOption("-st", "sequencing technology",   "sets the sequencing technology: '454', 'helicos', 'illumina', 'sanger' or 'solid'", "", settings.HasSequenceTechnology,   settings.SequenceTechnologyString, pMetadataOpts);

	// add the read archive options
	OptionGroup* pReadArchiveOpts = COptions::CreateOptionGroup("Read Archive Options");
	COptions::AddValueOption("-out", "MOSAIK read filename", "the output read file",                     "", settings.HasOutputReadsFilename, settings.OutputReadsFilename, pReadArchiveOpts);
	COptions::AddValueOption("-p",   "read name prefix",     "adds the prefix to each read name",        "", settings.HasReadNamePrefix,      settings.ReadNamePrefix,      pReadArchiveOpts);
	COptions::AddValueOption("-rl",  "# of reads",           "limits the # of reads processed",          "", settings.HasReadLimit,           settings.ReadLimit,           pReadArchiveOpts);
	COptions::AddValueOption("-tn",  "# of characters",      "sets the max # of internal Ns allowed",    "", settings.SetNumNBasesAllowed,    settings.NumNBasesAllowed,    pReadArchiveOpts);
	COptions::AddValueOption("-tp",  "# of beginning bases", "trims the first # of bases",               "", settings.HasTrimPrefixBases,     settings.NumTrimPrefixBases,  pReadArchiveOpts);
	COptions::AddValueOption("-ts",  "# of end bases",       "trims the last # of bases",                "", settings.HasTrimSuffixBases,     settings.NumTrimSuffixBases,  pReadArchiveOpts);
	COptions::AddValueOption("-tpr", "# of characters",      "trims the first characters from the name", "", settings.HasTrimPrefixName,      settings.NumTrimPrefixName,   pReadArchiveOpts);
	COptions::AddValueOption("-tsr", "# of characters",      "trims the last characters from the name",  "", settings.HasTrimSuffixName,      settings.NumTrimSuffixName,   pReadArchiveOpts);

	// parse the current command line
	COptions::Parse(argc, argv);

	// =============================
	// check for missing information
	// =============================

	bool foundError = false;
	ostringstream errorBuilder;
	const string ERROR_SPACER(7, ' ');

	if(!settings.HasReadFastaFilename && !settings.HasFastqFilename && !settings.HasSrfFilename && !settings.HasBustardDirectory && !settings.HasGeraldDirectory) {
		errorBuilder << ERROR_SPACER << "The input read file not specified. Please select either the -bd, -fr, -gd, -q, or -srf parameters." << endl;
		foundError = true;
	}

	if((settings.HasBustardDirectory && !settings.HasIlluminaLanes) || (settings.HasGeraldDirectory && !settings.HasIlluminaLanes)) {
		errorBuilder << ERROR_SPACER << "The desired Illumina lanes were not specified. Please use the -il parameter. e.g. \"-il 578\" means use lanes 5, 7, and 8." << endl;
		foundError = true;
	}

	if(!settings.HasOutputReferenceFilename && settings.HasReadFastaFilename) {

		if(settings.HasBaseQualityFastaFilename && settings.UseAssignedBQ) settings.UseAssignedBQ = false;

		if(!settings.HasBaseQualityFastaFilename && !settings.UseAssignedBQ) {
			errorBuilder << ERROR_SPACER << "Base quality parameters were not specified. Please use either the -fq or -assignQual parameter." << endl;
			foundError = true;
		}
	}

	if(settings.HasOutputReferenceFilename && (settings.HasFastqFilename || settings.HasSrfFilename)) {
		errorBuilder << ERROR_SPACER << "MOSAIK reference file creation is only supported with FASTA input files (-fr)." << endl;
		foundError = true;
	}

	if(settings.HasOutputReferenceFilename && (settings.ReadFastaFilename == settings.OutputReadsFilename)) {
		errorBuilder << ERROR_SPACER << "The input filename is the same as the output filename. Please select a different output filename." << endl;
		foundError = true;
	}

	if(!settings.HasOutputReadsFilename && !settings.HasOutputReferenceFilename) {
		errorBuilder << ERROR_SPACER << "An output read filename or reference filename was not specified. Please select either the -out or -oa parameter." << endl;
		foundError = true;
	}

	// handle the obvious choices before we complain about missing sequencing technology
	if(settings.HasBustardDirectory || settings.HasGeraldDirectory) settings.SequenceTechnologyString = "illumina";
	if(settings.EnableColorspace)                                   settings.SequenceTechnologyString = "solid";

	// check the base quality
	if(settings.UseAssignedBQ && ((settings.AssignedBQ < 1) || (settings.AssignedBQ > 99))) {
		errorBuilder << ERROR_SPACER << "The assigned base quality should be between 1 and 99." << endl;
		foundError = true;
	}

	// check the sequencing technology
	SequencingTechnologies seqTech = ST_UNKNOWN;

	if(!settings.HasOutputReferenceFilename) {
		CSequenceUtilities::LowercaseSequence(settings.SequenceTechnologyString);

		if(settings.SequenceTechnologyString == "454") {
			seqTech = ST_454;
		} else if(settings.SequenceTechnologyString == "helicos") {
			seqTech = ST_HELICOS;
		} else if(settings.SequenceTechnologyString == "illumina") {
			seqTech = ST_ILLUMINA;
		} else if(settings.SequenceTechnologyString == "sanger") {
			seqTech = ST_SANGER;
		} else if(settings.SequenceTechnologyString == "solid") {
			seqTech = ST_SOLID;
			settings.EnableColorspace = true;
		} else {
			errorBuilder << ERROR_SPACER << "Unknown sequencing technology. Please choose between '454', 'helicos', 'illumina', 'sanger', or 'solid'." << endl;
			foundError = true;
		}
	}

	// enforce metadata string length limitations
	if(settings.CenterName.size() > 255) {
		errorBuilder << ERROR_SPACER << "The 'center name' must be shorter than 256 characters." << endl;
		foundError = true;
	}

	if(settings.Description.size() > 65535) {
		errorBuilder << ERROR_SPACER << "The 'description' must be shorter than 65536 characters." << endl;
		foundError = true;
	}

	if(settings.LibraryName.size() > 255) {
		errorBuilder << ERROR_SPACER << "The 'library name' must be shorter than 256 characters." << endl;
		foundError = true;
	}

	if(settings.PlatformUnit.size() > 255) {
		errorBuilder << ERROR_SPACER << "The 'platform unit' must be shorter than 256 characters." << endl;
		foundError = true;
	}

	if(settings.ReadGroupID.size() > 255) {
		errorBuilder << ERROR_SPACER << "The 'read group ID' must be shorter than 256 characters." << endl;
		foundError = true;
	}

	if(settings.SampleName.size() > 255) {
		errorBuilder << ERROR_SPACER << "The 'sample name' must be shorter than 256 characters." << endl;
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

	// create a read group ID if one was not chosen
	if(!settings.HasReadGroupID) {
		settings.HasReadGroupID = true;
		CMosaikBuild::CreateReadGroupID(settings.ReadGroupID);
	}

	// create an unknown sample name
	if(!settings.HasSampleName) {
		settings.HasSampleName = true;
		settings.SampleName    = "unknown";
	}

	// build our metadata object
	MosaikReadFormat::ReadGroup rg;
	rg.CenterName           = settings.CenterName;
	rg.Description          = settings.Description;
	rg.LibraryName          = settings.LibraryName;
	rg.MedianFragmentLength = settings.MedianFragmentLength;
	rg.PlatformUnit         = settings.PlatformUnit;
	rg.ReadGroupID          = settings.ReadGroupID;
	rg.SampleName           = settings.SampleName;
	rg.SequencingTechnology = seqTech;

	// test to see if the specified input files exist
	vector<string> mate1Files, mate2Files;

	if(settings.HasReadFastaFilename) {
		CFasta::CheckFile(settings.ReadFastaFilename, true);
	} else if(settings.HasFastqFilename) {

		// ================
		// check first mate
		// ================

		{
			// check if this is a directory or file
			bool parseDirectory = false;
			if(CFileUtilities::DirExists(settings.FastqFilename.c_str())) parseDirectory = true;
			else CFastq::CheckFile(settings.FastqFilename, true);

			// populate our file vector
			if(parseDirectory) {

				vector<string> dirFiles;
				CFileUtilities::SearchDirectory(dirFiles, settings.FastqFilename.c_str());

				// look for FASTQ files
				for(unsigned int i = 0; i < (unsigned int)dirFiles.size(); i++)
					if(CFastq::CheckFile(dirFiles[i], false)) 
						mate1Files.push_back(dirFiles[i]);

			} else mate1Files.push_back(settings.FastqFilename);
		}

		// =================
		// check second mate
		// =================

		if(settings.HasFastq2Filename) {

			// check if this is a directory or file
			bool parseDirectory = false;
			if(CFileUtilities::DirExists(settings.Fastq2Filename.c_str())) parseDirectory = true;
			else CFastq::CheckFile(settings.Fastq2Filename, true);

			// populate our file vector
			if(parseDirectory) {

				vector<string> dirFiles;
				CFileUtilities::SearchDirectory(dirFiles, settings.Fastq2Filename.c_str());

				// look for FASTQ files
				for(unsigned int i = 0; i < (unsigned int)dirFiles.size(); i++)
					if(CFastq::CheckFile(dirFiles[i], false)) 
						mate2Files.push_back(dirFiles[i]);

			} else mate2Files.push_back(settings.Fastq2Filename);
		}

	} else if(settings.HasSrfFilename) {

		// check if this is a directory or file
		bool parseDirectory = false;
		if(CFileUtilities::DirExists(settings.SrfFilename.c_str())) parseDirectory = true;
		else CSRF::CheckFile(settings.SrfFilename, true);

		// populate our file vector
		if(parseDirectory) {

			vector<string> dirFiles;
			CFileUtilities::SearchDirectory(dirFiles, settings.SrfFilename.c_str());

			// look for SRF files
			for(unsigned int i = 0; i < (unsigned int)dirFiles.size(); i++)
				if(CSRF::CheckFile(dirFiles[i], false)) 
					mate1Files.push_back(dirFiles[i]);

		} else mate1Files.push_back(settings.SrfFilename);

	}

	if(settings.HasBaseQualityFastaFilename)  CFileUtilities::CheckFile(settings.BaseQualityFastaFilename.c_str(), true);
	if(settings.HasBaseQualityFasta2Filename) CFileUtilities::CheckFile(settings.BaseQualityFasta2Filename.c_str(), true);

	if(!mate2Files.empty() && settings.HasTrimSuffixName) {
		settings.HasTrimSuffixName = false;
		settings.NumTrimSuffixName = 0;
		cout << "ERROR: Suffix trimming will not performed on the read name because it may interfere with the automatic paired-end/mate-pair suffix trimming." << endl << endl;
	}

	// start benchmarking
	CBenchmark bench;
	bench.Start();

	// time to create a new library
	CMosaikBuild mb(rg);

	// output the metadata information
	if(!settings.HasOutputReferenceFilename) {
		if(settings.HasCenterName)               cout << "- setting center name to: " << rg.CenterName << endl;
		if(settings.HasDescription)              cout << "- setting description to: " << rg.Description << endl;
		if(settings.HasReadGroupID)              cout << "- setting read group ID to: " << rg.ReadGroupID << endl;
		if(settings.HasLibraryName)              cout << "- setting library name to: " << rg.LibraryName << endl;
		if(settings.HasMedianFragmentLength)     cout << "- setting median fragment length to: " << rg.MedianFragmentLength << endl;
		if(settings.HasPlatformUnit)             cout << "- setting platform unit to: " << rg.PlatformUnit << endl;
		if(settings.HasSampleName)               cout << "- setting sample name to: " << rg.SampleName << endl;
		cout << "- setting sequencing technology to: " << settings.SequenceTechnologyString << endl;
	}

	// enable colorspace handling
	if(settings.EnableColorspace) {
		if(settings.HasOutputReferenceFilename) cout << "- enabling AB SOLiD colorspace conversion" << endl;	
		mb.EnableColorspace();
	}

	// set the max number of N's allowed
	if(!settings.HasOutputReferenceFilename) {
		cout << "- trimming leading and lagging N's. ";
		if(settings.NumNBasesAllowed == 0) cout << "Mates with interior N's will not be deleted." << endl;	
		else cout << "Mates with >" << settings.NumNBasesAllowed << " interior N's will be deleted." << endl;
		if(settings.SetNumNBasesAllowed) mb.SetNumNBasesAllowed(settings.NumNBasesAllowed);
	} else {

		if(settings.HasGenomeAssemblyID) {
			cout << "- setting genome assembly ID to \"" << settings.GenomeAssemblyID << "\"" << endl;
			mb.SetGenomeAssemblyID(settings.GenomeAssemblyID);
		}

		if(settings.HasUniformResourceIdentifier) {
			cout << "- setting URI to \"" << settings.UniformResourceIdentifier << "\"" << endl;
			mb.SetURI(settings.UniformResourceIdentifier);
		}

		if(settings.HasSpeciesName) {
			cout << "- setting species name to \"" << settings.SpeciesName << "\"" << endl;
			mb.SetSpecies(settings.SpeciesName);
		}

		cout << "- converting " << settings.ReadFastaFilename << " to a reference sequence archive." << endl;
	}

	// enable read and read name trimming
	if(settings.HasTrimPrefixBases || settings.HasTrimSuffixBases) {
		cout << "- trimming the first " << settings.NumTrimPrefixBases << " and the last " << settings.NumTrimSuffixBases << " bases" << endl;
		mb.EnableBaseTrimming(settings.NumTrimPrefixBases, settings.NumTrimSuffixBases);		
	}

	if(settings.HasTrimPrefixName || settings.HasTrimSuffixName) {
		cout << "- trimming the first " << settings.NumTrimPrefixName << " and the last " << settings.NumTrimSuffixName << " characters of the read name" << endl;	
		mb.EnableReadNameTrimming(settings.NumTrimPrefixName, settings.NumTrimSuffixName);
	}

	// enable the addition of a user specified read name prefix
	if(settings.HasReadNamePrefix) {
		cout << "- prepending all read names with \"" << settings.ReadNamePrefix << "\"" << endl;	
		mb.EnableReadNamePrefix(settings.ReadNamePrefix);
	}

	// enable the read limit
	if(settings.HasReadLimit) {
		cout << "- limiting read archive to " << settings.ReadLimit << " reads" << endl;	
		mb.EnableReadLimit(settings.ReadLimit);
	}

	cout << endl;

	// ================
	// parse read files
	// ================

	if(settings.HasOutputReferenceFilename) {

		mb.CreateReferenceArchive(settings.ReadFastaFilename, settings.OutputReferenceFilename);

	} else if(settings.HasReadFastaFilename) {

		// enable base quality parsing
		if(settings.UseAssignedBQ) {
			mb.SetAssignedBaseQuality(settings.AssignedBQ);
		} else {
			if(settings.HasBaseQualityFastaFilename)  mb.EnableBaseQualities(settings.BaseQualityFastaFilename);
			if(settings.HasBaseQualityFasta2Filename) mb.EnableBaseQualities2(settings.BaseQualityFasta2Filename);
		}

		// parse the FASTA files
		if(settings.HasReadFasta2Filename) {
			mb.ParsePEFasta(settings.ReadFastaFilename, settings.ReadFasta2Filename, settings.OutputReadsFilename);
		} else {
			mb.ParseFasta(settings.ReadFastaFilename, settings.OutputReadsFilename);
		}

	} else if(settings.HasFastqFilename) {

		// parse the FASTQ file
		if(!mate2Files.empty()) {
			mb.ParsePEFastq(mate1Files, mate2Files, settings.OutputReadsFilename);
		} else {
			mb.ParseFastq(mate1Files, settings.OutputReadsFilename);
		}

	} else if(settings.HasSrfFilename) {

		// parse the SRF file
		mb.ParseSRF(mate1Files, settings.OutputReadsFilename);

	} else if(settings.HasBustardDirectory) {

		// parse the Illumina Bustard directory
		mb.ParseBustard(settings.BustardDirectory, settings.IlluminaLanesString, settings.OutputReadsFilename, settings.SplitBustardReads);

	} else if(settings.HasGeraldDirectory) {

		// parse the Illumina Gerald directory
		mb.ParseGerald(settings.GeraldDirectory, settings.IlluminaLanesString, settings.OutputReadsFilename);

	} else {
		cout << "ERROR: None of the parsing routines were picked. Nothing to do." << endl;
		exit(1);
	}

	// ==================
	// Show total runtime
	// ==================

	// stop benchmarking
	bench.Stop();

	// show the benchmarking results
	cout << endl;
	bench.DisplayTime("MosaikBuild");

	return 0;
}
