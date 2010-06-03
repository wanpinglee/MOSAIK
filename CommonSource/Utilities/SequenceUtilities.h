// ***************************************************************************
// CSequenceUtilities - handles basic sequence manipulation.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <cstdlib>
#include <cstring>
#include <map>
#include <string>

using namespace std;

class CSequenceUtilities {
public:
	// Performs an in-place reverse complement conversion
	static void GetReverseComplement(char* seqBases, const unsigned int seqLength);
	// Performs an in-place sequence reversal using a C string
	static inline void ReverseSequence(char* seqBases, const unsigned int seqLength);
	// Converts an STL string to uppercase
	static inline void UppercaseSequence(string& s);
	// Converts an STL string to lowercase
	static inline void LowercaseSequence(string& s);
	// Trims the carriage return at the end of a string
	static inline void Chomp(char* s);
};

// Trims the carriage return at the end of a string
inline void CSequenceUtilities::Chomp(char* s) {
	size_t sLen = strlen(s);
	if(sLen == 0) return;
	sLen--;

	while((s[sLen] == 10) || (s[sLen] == 13)) {
		s[sLen] = 0;
		sLen--;
		if(sLen < 0) break;
	}
}

// Converts a STL string to uppercase
inline void CSequenceUtilities::UppercaseSequence(string& s) {
	char* sc = (char*)s.data();
	for(unsigned int i = 0; i < (unsigned int)s.size(); i++) sc[i] = toupper(sc[i]);
}

// Converts an STL string to lowercase
inline void CSequenceUtilities::LowercaseSequence(string& s) {
	char* sc = (char*)s.data();
	for(unsigned int i = 0; i < (unsigned int)s.size(); i++) sc[i] = tolower(sc[i]);
}

// Performs an in-place sequence reversal using a C string
inline void CSequenceUtilities::ReverseSequence(char* seqBases, const unsigned int seqLength) {
	for(unsigned int i = seqLength; i >= (seqLength / 2) + 1; i--) 
		swap(seqBases[seqLength - i], seqBases[i - 1]);
}
