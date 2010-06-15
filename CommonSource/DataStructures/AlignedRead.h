// ***************************************************************************
// AlignedRead.h - stores a collection of mate 1 and mate 2 alignments.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <vector>
#include <algorithm>
#include "Alignment.h"
#include "MosaikString.h"

using namespace std;

namespace Mosaik {
	struct AlignedRead {
		unsigned int ReadGroupCode;
		CMosaikString Name;
		vector<Alignment> Mate1Alignments;
		vector<Alignment> Mate2Alignments;
		bool IsLongRead;
		bool IsPairedEnd;

		// constructor
		AlignedRead()
			: ReadGroupCode(0)
			, IsLongRead(false)
			, IsPairedEnd(false)
		{}

                bool Clear(){
                        ReadGroupCode   = 0;
                        Name.clear();
                        Mate1Alignments.clear();
                        Mate2Alignments.clear();
                        IsLongRead = false;
                        IsPairedEnd= false;

			if ( !Mate1Alignments.empty() ) {
				cout << "ERROR: Clearing AlignedRead is failed." << endl;
				exit(1);
			}

			if ( !Mate2Alignments.empty() ) {
				cout << "ERROR: Clearing AlignedRead is failed." << endl;
				exit(1);
			}

                        return true;
		}

                bool SortAlignment() {
                        sort(Mate1Alignments.begin(), Mate1Alignments.end() );
                        sort(Mate2Alignments.begin(), Mate2Alignments.end() );

                        return true;
                }

                bool operator<( AlignedRead& x ) {
                        return Name < x.Name;
                }



	};
}
