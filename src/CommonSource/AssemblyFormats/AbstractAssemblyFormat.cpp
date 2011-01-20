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

#include "AbstractAssemblyFormat.h"

// constructor
CAbstractAssemblyFormat::CAbstractAssemblyFormat(void)
: mIsOpen(false)
{}

// destructor
CAbstractAssemblyFormat::~CAbstractAssemblyFormat(void) {}

// sets the ungapped to gapped vector
void CAbstractAssemblyFormat::SetUngappedToGappedVector(unsigned int* pVector) {
	mpUngap2Gap = pVector;
}
