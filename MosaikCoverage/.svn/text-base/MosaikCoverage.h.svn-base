// ***************************************************************************
// CMosaikCoverage - exports coverage graphs using alignment archives.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "AlignmentReader.h"
#include "ConsoleUtilities.h"
#include "FileUtilities.h"
#include "ReferenceSequenceReader.h"

using namespace std;

class CMosaikCoverage {
public:
	// constructor
	CMosaikCoverage(void);
	// destructor
	~CMosaikCoverage(void);
	// toggles the creation of output files
	void DisableOutput(void);
	// toggles the creation of graphs via gnuplot
	void EnableGraphCreation(void);
	// when triggered, the coverage calculation will only include unique reads
	void EvaluateUniqueReadsOnly(void);
	// parses the specified MOSAIK alignment file and matching anchors file
	void ParseMosaikAlignmentFile(const vector<string>& alignmentFilenames, const string& anchorsFilename);
	// saves the coverage files to the specified output directory
	void SaveCoverage(const unsigned char minCoverage);
        // sets the minimum alignment quality threshold
	void SetMinAlignmentQuality(const unsigned char minAQ);
	// sets the output directory
	void SetOutputDirectory(const string& directory);

private:
	// our coverage array
	unsigned int** mpCoverageArray;
	// toggles the creation of output files
	bool mDisableOutput;
	// toggles unique read evaluation
	bool mEvaluateUniqueReadsOnly;
	// toggles postscript graph creation
	bool mCreatePostscriptGraphs;
	// toggles the alignment quality filter
	bool mUseAlignmentQualityFilter;
	unsigned char mAlignmentQualityThreshold;
	// the number of anchors in the vector
	unsigned int mNumRefSeqs;
	// specifies the output directory
	string mOutputDirectory;
	// specifies the reference sequence vector
	vector<ReferenceSequence> mReferenceSequences;
};
