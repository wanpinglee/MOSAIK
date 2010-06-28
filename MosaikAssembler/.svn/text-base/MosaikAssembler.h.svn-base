// ***************************************************************************
// CMosaikAssembler - consolidates gaps in the underlying alignments and
//                    exports the data set into an assembly format.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iostream>
#include <string>
#include "Ace.h"
#include "AbstractAssemblyFormat.h"
#include "AlignmentReader.h"
#include "Benchmark.h"
#include "ConsoleUtilities.h"
#include "GapInfo.h"
#include "GigaBayesFormat.h"
#include "ProgressBar.h"
#include "ReferenceSequenceReader.h"
#include "UnorderedMap.h"

using namespace std;

class CMosaikAssembler {
public:
	// our enumerated alignment algorithms
	enum AssemblyFormatType {
		AssemblyFormat_ACE,
		AssemblyFormat_GIG
	};
	// constructor
	CMosaikAssembler(AssemblyFormatType format, unsigned char referenceBaseQuality);
	// destructor
	~CMosaikAssembler(void);
	// assembles the specified alignment file
	void Assemble(const string& referenceFilename, const string& alignmentFilename, const string& outputFilenameStub);
	// enables FASTA file creation for gigaBayes
	void EnableFastaFileCreation(void);
	// selects the specified region of interest
	void SelectRegionOfInterest(const string& referenceName);

private:
	// define our assembly configuration structure
	struct AssemblySettings {
		unsigned char ReferenceSequenceBaseQuality;
		string RoiName;
		AssemblyFormatType AssemblyFormat;

		AssemblySettings()
			: AssemblyFormat(AssemblyFormat_ACE)
		{}
	} mSettings;
	// define our flags structure
	struct FlagData {
		bool EnableFastaFileCreation;
		bool HasRegionOfInterest;

		FlagData()
			: EnableFastaFileCreation(false)
			, HasRegionOfInterest(false)
		{}
	} mFlags;
	// define our reference pad structure
	struct ReferencePadInfo {
		unsigned int Position;
		unsigned short Length;
		char Nucleotide;
	};
	// define our read pad structure
	struct ReadPadInfo {
		unsigned int Position;
		unsigned short Length;
		char ReferenceNucleotide;
		char QueryNucleotide;
	};
	// inserts observed gaps to the specified read
	void AddGapsToRead(const vector<GapInfo>& gaps, const Alignment& al, CMosaikString& gappedRead);
	// inserts observed gaps to the reference sequence
	void AddGapsToReferenceSequence(const vector<GapInfo>& gaps, string& referenceSequence, CMosaikString& gappedReferenceSequence);
	// re-allocated the read pad buffer if more space is needed
	void CheckReadPadBuffer(const unsigned int readLength);
	// builds a ungapped to gapped conversion vector from the specified reference sequence
	void PopulateUngappedToGappedVector(const CMosaikString& reference, const unsigned int ungappedLength);
	// create our ungapped to gapped vector
	unsigned int* mpUngap2Gap;
	unsigned int mUngap2GapLen;
	// our read pad buffer
	ReadPadInfo* mReadPadBuffer;
	unsigned int mReadPadBufferLen;
	// our gap iterators
	vector<GapInfo>::const_iterator mGapIter;
};
