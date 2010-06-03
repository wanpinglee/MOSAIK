// ***************************************************************************
// GapInfo - stores the reference gap locations. (Used in CMosaikAssembler)
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

struct GapInfo {
	unsigned int Position;
	unsigned int Length;

	bool operator<(const GapInfo& gi) const { 
		return Position < gi.Position; 
	}
};
