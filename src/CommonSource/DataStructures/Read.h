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

#define SOLID_PREFIX_LENGTH 1

namespace Mosaik {

	struct Mate {
		CMosaikString Bases;
		CMosaikString Qualities;
		char SolidPrefixTransition[SOLID_PREFIX_LENGTH];

		bool clear(void) {
			Bases.clear();
			Qualities.clear();

			return true;
		}
	};

	struct Read {
		unsigned int ReadGroupCode;
		unsigned int Owner; // the temporary file that contains the read
		CMosaikString Name;
		Mate Mate1;
		Mate Mate2;
		
		bool clear(void) {
			ReadGroupCode = 0;
			Owner         = 0;
			Name.clear();
			Mate1.clear();
			Mate2.clear();

			return true;
		}

		bool operator<(const Read& i) const{
			return Name < i.Name;
		}
	};
}
