// ***************************************************************************
// ReadGroup.h - stores data related to a read group.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <string>
#include "SequencingTechnologies.h"
#include "SequenceUtilities.h"
#include "SHA1.h"

using namespace std;

namespace MosaikReadFormat {
	struct ReadGroup {
		unsigned int MedianFragmentLength;
		unsigned int ReadGroupCode;
		SequencingTechnologies SequencingTechnology;
		string CenterName;
		string Description;
		string LibraryName;
		string PlatformUnit;
		string ReadGroupID;
		string SampleName;

		// constructor
		ReadGroup(void)
			: MedianFragmentLength(0)
			, ReadGroupCode(0)
			, SequencingTechnology(ST_UNKNOWN)
		{}

		// create a 32-bit identifier for the read group code
		static unsigned int GetCode(const ReadGroup& readGroup) {
			const unsigned int readGroupCode = CSHA1::GenerateReadGroupCode(readGroup.ReadGroupID, readGroup.SampleName);
			return readGroupCode;
		}
	};
}
