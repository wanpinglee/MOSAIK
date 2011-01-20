// ***************************************************************************
// CAbstractAssemblyFormat - superclass to all of the exported assembly
//                           formats.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// This code is dual licenced under the GNU General Public License 2.0+ or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <string>
#include "Alignment.h"
#include "Mosaik.h"
#include "MosaikString.h"
#include "UnorderedMap.h"

using namespace std;

class CAbstractAssemblyFormat {
public:
	// constructor
	CAbstractAssemblyFormat(void);
	// destructor
	virtual ~CAbstractAssemblyFormat(void) = 0;
	// opens the ace file
	virtual void Open(const CMosaikString& filename) = 0;
	// closes the ace file
	virtual void Close(void) = 0;
	// saves the specified reference sequence to the header
	virtual void SaveHeader(CMosaikString& reference, const string& referenceName, const unsigned int ungappedRefLength, const uint64_t& numSequences) = 0;
	// saves the specified read to the index and reads files
	virtual void SaveRead(const Alignment& al, CMosaikString& gappedRead) = 0;
	// sets the ungapped to gapped vector
	void SetUngappedToGappedVector(unsigned int* pVector);

protected:
	// our file state
	bool mIsOpen;
	// the ungapped to gapped conversion vector
	unsigned int* mpUngap2Gap;
	// our gapped reference length
	unsigned int mGappedRefLength;
};
