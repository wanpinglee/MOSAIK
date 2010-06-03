// ***************************************************************************
// Read.h - stores everything related to reads. (pre-alignment)
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include "MosaikString.h"

#define SOLID_PREFIX_LENGTH 2

namespace Mosaik {

	struct Mate {
		CMosaikString Bases;
		CMosaikString Qualities;
		char SolidPrefixTransition[SOLID_PREFIX_LENGTH];
	};

	struct Read {
		unsigned int ReadGroupCode;
		CMosaikString Name;
		Mate Mate1;
		Mate Mate2;
	};
}
