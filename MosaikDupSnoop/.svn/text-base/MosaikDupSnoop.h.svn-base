// ***************************************************************************
// CMosaikDupSnoop - records duplicate fragments in sequencing libraries.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <algorithm>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "sqlite3.h"
#include "Alignment.h"
#include "AlignmentReader.h"
#include "Benchmark.h"
#include "ConsoleUtilities.h"
#include "ConversionUtilities.h"
#include "FileUtilities.h"
#include "ProgressBar.h"
#include "ProgressCounter.h"

using namespace std;

static double DEFAULT_CONFIDENCE_INTERVAL = 0.9973;
#define SQL_BUFFER_SIZE             10240

#ifdef WIN32
#define snprintf _snprintf
#endif

class CMosaikDupSnoop {
public:
	// constructor
	CMosaikDupSnoop(void);
	// destructor
	~CMosaikDupSnoop(void);
	// configures which read pair types should be resolved
	void ConfigureResolution(const bool uo, const bool uu, const bool um, const bool mm);
	// allows any fragment length when evaluating unique mate-pairs
	void EnableAllUniqueFragmentLengths(void);
	// records the fragments found in the input files
	void RecordFragments(const vector<string>& inputFiles, const string& outputDirectory);
	// sets the desired confidence interval
	void SetConfidenceInterval(const double& percent);

private:
	// define our sort configuration structure
	struct SortSettings {
		double ConfidenceInterval;

		SortSettings() 
			: ConfidenceInterval(DEFAULT_CONFIDENCE_INTERVAL)
		{}
	} mSettings;
	// define our boolean flags structure
	struct FlagData {
		bool AllowAllUniqueFragmentLengths;
		bool ResolveMM;
		bool ResolveUM;
		bool ResolveUO;
		bool ResolveUU;

		FlagData()
			: AllowAllUniqueFragmentLengths(false)
			, ResolveMM(false)
			, ResolveUM(false)
			, ResolveUO(false)
			, ResolveUU(false)
		{}
	} mFlags;
	// calculates the aggregate quality
	unsigned int CalculateAggregateQuality(vector<Alignment>::iterator& alIter);
	// consolidates the PairedFragments table taking into account variable endpoints
	void ConsolidatePairedEndFragments(const string& outputDirectory, const set<string>& libraryNames);
	// creates databases for each library if they don't already exist
	void CreateLibraryDatabases(const string& outputDirectory, const set<string>& libraryNames);
	// populates the library databases with alignment data
	void PopulateLibraryDatabases(const vector<string>& inputFiles, const string& outputDirectory);
	// retrieves the library names used in the supplied alignment archives
	void RetrieveLibraryNames(const vector<string>& inputFiles, set<string>& libraryNames);
	// our SQL data buffer
	char mSqlBuffer[SQL_BUFFER_SIZE];
};
