// ***************************************************************************
// CSequenceUtilities - handles basic sequence manipulation.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#ifndef SEQUENCEUTILITIES_H_
#define SEQUENCEUTILITIES_H_

#include <cstdlib>
#include <cstring>
#include <map>
#include <string>

using namespace std;

namespace CSequenceUtilities {
	// Performs an in-place reverse complement conversion
	void GetReverseComplement(char* seqBases, const unsigned int seqLength);
	// Performs an in-place sequence reversal using a C string
	void ReverseSequence(char* seqBases, const unsigned int seqLength);
	// Converts an STL string to uppercase
	void UppercaseSequence(string& s);
	// Converts an STL string to lowercase
	void LowercaseSequence(string& s);
	// Trims the carriage return at the end of a string
	void Chomp(char* s);
	void ChompQuality(char* s);
};

#endif // SEQUENCEUTILITIES_H_
