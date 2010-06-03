// ***************************************************************************
// CRegexUtilities - performs regular expression-like tasks.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iostream>
#ifdef WIN32
#include <regex>
using namespace std::tr1;
#endif
#include <string>
#include <vector>
#include <cstdlib>
#include "ConversionUtilities.h"
#include "MosaikString.h"

using namespace std;

class CRegexUtilities {
public:
	// converts a space separated quality string into a compressed quality string
	static void ConvertQualities(string& qualities, CMosaikString& compQualities);
	// extracts the genome assembly ID from a FASTA/FASTQ header
	static void ExtractGenomeAssemblyID(const string& line, CMosaikString& genomeAssemblyID);
	// extracts the sequence name from a FASTA/FASTQ header
	static void ExtractSequenceName(const string& line, CMosaikString& name);
	// extracts the species name from a FASTA/FASTQ header
	static void ExtractSpecies(const string& line, CMosaikString& species);
	// extracts the URI from a FASTA/FASTQ header
	static void ExtractURI(const string& line, CMosaikString& uri);

private:
	// Trims the carriage return at the end of a string
	static inline void Chomp(char* s);
#ifdef WIN32
	// specifies our genome assembly ID regular expression
	static regex mGenomeAssemblyIDRegex;
	// specifies our sequence name regular expression
	static regex mSequenceNameRegex;
	// specifies our species name regular expression
	static regex mSpeciesRegex;
	// specifies our URI regular expression
	static regex mUriRegex;
#endif
};

// splits a string according to the supplied delimiter
template<typename I>
void SplitString(I& inserter, const string& delimiter, const string& s) {

	string::size_type lpos = 0;
	string::size_type pos  = s.find_first_of(delimiter, lpos);

	while(lpos != string::npos)	{
		if(lpos != pos) {
			if(pos != string::npos) *inserter = s.substr(lpos,pos - lpos);
			else *inserter = s.substr(lpos);			
			inserter++;
		}

		lpos = (pos == string::npos ? string::npos : pos + 1);
		pos  = s.find_first_of(delimiter, lpos);
	}
}

// splits a string according to the supplied delimiter
template<typename I>
void SplitStringEmpty(I& inserter, const string& delimiter, const string& s) {

	string::size_type lpos = 0;
	string::size_type pos  = s.find_first_of(delimiter, lpos);

	while(lpos != string::npos)	{
		if(pos != string::npos) *inserter = s.substr(lpos,pos - lpos);
		else *inserter = s.substr(lpos);			
		inserter++;

		lpos = (pos == string::npos ? string::npos : pos + 1);
		pos  = s.find_first_of(delimiter, lpos);
	}
}

// Trims the carriage return at the end of a string
inline void CRegexUtilities::Chomp(char* s) {
	size_t sLen = strlen(s);
	if(sLen == 0) return;
	sLen--;

	while((s[sLen] == 10) || (s[sLen] == 13)) {
		s[sLen] = 0;
		sLen--;
		if(sLen < 0) break;
	}
}
