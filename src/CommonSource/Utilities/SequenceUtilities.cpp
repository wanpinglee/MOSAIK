// ***************************************************************************
// CSequenceUtilities - handles basic sequence manipulation.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "SequenceUtilities.h"

// Performs an in-place reverse complement conversion
void CSequenceUtilities::GetReverseComplement(char* seqBases, const unsigned int seqLength) {

	// reverse the sequence
	ReverseSequence(seqBases, seqLength);

	// swap the bases
	for(unsigned int i = 0; i < seqLength; i++) {
		switch(seqBases[i]) {
			case 'A':
				seqBases[i] = 'T';
				break;
			case 'C':
				seqBases[i] = 'G';
				break;
			case 'G':
				seqBases[i] = 'C';
				break;
			case 'T':
				seqBases[i] = 'A';
				break;
			default:
				break;
		}
	}
}


// Trims the carriage return at the end of a string
void CSequenceUtilities::Chomp(char* s) {
	size_t sLen = strlen(s);
	if(sLen == 0) return;
	sLen--;

	while((s[sLen] == 10) || (s[sLen] == 13)) {
		s[sLen] = 0;
		sLen--;
		if(sLen < 0) break;
	}
}

// Trims the carriage return at the end of a string
void CSequenceUtilities::ChompQuality(char* s) {
	size_t sLen = strlen(s);
	if(sLen == 0) return;
	sLen--;
	
	while((s[sLen] == 10) || (s[sLen] == 13)) {
		s[sLen] = ' ';
		sLen--;
		if(sLen < 0) break;
	}
}

// Converts a STL string to uppercase
void CSequenceUtilities::UppercaseSequence(string& s) {
	char* sc = (char*)s.data();
	for(unsigned int i = 0; i < (unsigned int)s.size(); i++) sc[i] = toupper(sc[i]);
}

// Converts an STL string to lowercase
void CSequenceUtilities::LowercaseSequence(string& s) {
	char* sc = (char*)s.data();
	for(unsigned int i = 0; i < (unsigned int)s.size(); i++) sc[i] = tolower(sc[i]);
}

// Performs an in-place sequence reversal using a C string
void CSequenceUtilities::ReverseSequence(char* seqBases, const unsigned int seqLength) {
	for(unsigned int i = seqLength; i >= (seqLength / 2) + 1; i--) 
		swap(seqBases[seqLength - i], seqBases[i - 1]);
}

