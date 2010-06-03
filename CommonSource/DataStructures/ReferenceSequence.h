// ***************************************************************************
// ReferenceSequence.h - stores everything related to reference sequences.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <string>
#include "Mosaik.h"
#include "LargeFileSupport.h"

using namespace std;

struct ReferenceSequence {
	off_type BasesOffset;
	uint64_t NumAligned;
	unsigned int Begin;
	unsigned int End;
	unsigned int NumBases;
	string Name;
	string Bases;
	string GenomeAssemblyID;
	string Species;
	string MD5;
	string URI;

	// constructor
	ReferenceSequence()
		: BasesOffset(0)
		, NumAligned(0)
		, Begin(0)
		, End(0)
		, NumBases(0)
	{}
};
