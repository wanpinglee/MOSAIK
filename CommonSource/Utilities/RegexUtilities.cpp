// ***************************************************************************
// CRegexUtilities - performs regular expression-like tasks.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "RegexUtilities.h"

#ifdef WIN32

// specifies our genome assembly ID regular expression
regex CRegexUtilities::mGenomeAssemblyIDRegex("GA\\((.+?)\\)");

// specifies our sequence name regular expression
regex CRegexUtilities::mSequenceNameRegex("^[+>@]\\s*(\\S+)");

// specifies our species name regular expression
regex CRegexUtilities::mSpeciesRegex("SN\\((.+?)\\)");

// specifies our URI regular expression
regex CRegexUtilities::mUriRegex("URI\\((.+?)\\)");

#endif

// converts a space separated quality string into a compressed quality string
// NOTE: this function has horrible amounts of overhead, but lean and mean code that I had before
//       failed some of the unit tests.
void CRegexUtilities::ConvertQualities(string& qualities, CMosaikString& compQualities) {

	vector<string> columns;
	vector<string>::const_iterator sIter;

	char* pQualities = (char*)qualities.c_str();
	Chomp(pQualities);

	back_insert_iterator<vector<string> > backiter(columns);
	SplitString(backiter, " ", pQualities);
	const unsigned int numQualities = (unsigned int)columns.size();

	compQualities.Reserve(numQualities);
	compQualities.SetLength(numQualities);

	unsigned char* pCompQualities = (unsigned char*)compQualities.Data();

	for(sIter = columns.begin(); sIter != columns.end(); ++sIter, ++pCompQualities) {
		if(sIter->empty()) continue;
		*pCompQualities = GetUnsignedChar((char*)sIter->c_str());
	}
}

// extracts the genome assembly ID from a FASTA/FASTQ header
void CRegexUtilities::ExtractGenomeAssemblyID(const string& line, CMosaikString& genomeAssemblyID) {
#ifdef WIN32

	cmatch results;
	if(!regex_search(line.c_str(), results, mGenomeAssemblyIDRegex)) {
		genomeAssemblyID.SetLength(0);
		return;
	}
	genomeAssemblyID = results[1].str().c_str();

#else

	// TODO: replace this with the TR1 regex above when it finally works in gcc. It doesn't work in gcc 4.3.3

	// find the GA tag
	const string gaTag = "GA(";
	string::size_type gaPos = line.find(gaTag.c_str());

	if(gaPos == string::npos) {
		genomeAssemblyID.SetLength(0);
		return;
	}

	// find the matching end parenthesis
	const unsigned int start = gaPos + gaTag.size();
	unsigned int stop = start;

	const char* pBuffer = line.data();
	unsigned int lineLen = line.size();

	if(stop < lineLen) {
		while(pBuffer[stop] != ')') {
			stop++;
			if(stop == lineLen) break;
		}
	}

	if(start == stop) {
		cout << "ERROR: could not parse genome assembly ID from FASTA header." << endl;
		cout << "       " << line << endl;
		exit(1);
	}

	genomeAssemblyID = line.substr(start, stop - start).c_str();

#endif
}

// extracts the sequence name from a header
void CRegexUtilities::ExtractSequenceName(const string& line, CMosaikString& name) {
#ifdef WIN32

	cmatch results;
	if(!regex_search(line.c_str(), results, mSequenceNameRegex)) {
		printf("ERROR: The sequence name was not captured by our regular expression. Sequence name: [%s]\n", line.c_str());
		exit(1);
	}
	name = results[1].str().c_str();

#else

	// TODO: replace this with the TR1 regex above when it finally works in gcc. It doesn't work in gcc 4.3.3

	// get rid of the leading greater than sign
	string s = line.substr(1);

	// extract the first non-whitespace segment
	char* pName = (char*)s.data();
	unsigned int nameLen = (unsigned int)s.size();

	unsigned int start = 0;
	while((pName[start] == 32) || (pName[start] == 9) || (pName[start] == 10) || (pName[start] == 13)) {
		start++;
		if(start == nameLen) break;
	}

	unsigned int stop  = start;
	if(stop < nameLen) {
		while((pName[stop] != 32) && (pName[stop] != 9) && (pName[stop] != 10) && (pName[stop] != 13)) {
			stop++;
			if(stop == nameLen) break;
		}
	}

	if(start == stop) {
		cout << "ERROR: could not parse read name from FASTA header." << endl;
		cout << "       " << line << endl;
		exit(1);
	}

	name = s.substr(start, stop - start).c_str();

#endif
}

// extracts the species name from a FASTA/FASTQ header
void CRegexUtilities::ExtractSpecies(const string& line, CMosaikString& species) {
#ifdef WIN32

	cmatch results;
	if(!regex_search(line.c_str(), results, mSpeciesRegex)) {
		species.SetLength(0);
		return;
	}
	species = results[1].str().c_str();

#else

	// TODO: replace this with the TR1 regex above when it finally works in gcc. It doesn't work in gcc 4.3.3

	// find the SN tag
	const string snTag = "SN(";
	string::size_type snPos = line.find(snTag.c_str());

	if(snPos == string::npos) {
		species.SetLength(0);
		return;
	}

	// find the matching end parenthesis
	const unsigned int start = snPos + snTag.size();
	unsigned int stop = start;

	const char* pBuffer = line.data();
	unsigned int lineLen = line.size();

	if(stop < lineLen) {
		while(pBuffer[stop] != ')') {
			stop++;
			if(stop == lineLen) break;
		}
	}

	if(start == stop) {
		cout << "ERROR: could not parse genome assembly ID from FASTA header." << endl;
		cout << "       " << line << endl;
		exit(1);
	}

	species = line.substr(start, stop - start).c_str();

#endif
}

// extracts the URI from a FASTA/FASTQ header
void CRegexUtilities::ExtractURI(const string& line, CMosaikString& uri) {
#ifdef WIN32

	cmatch results;
	if(!regex_search(line.c_str(), results, mUriRegex)) {
		uri.SetLength(0);
		return;
	}
	uri = results[1].str().c_str();

#else

	// TODO: replace this with the TR1 regex above when it finally works in gcc. It doesn't work in gcc 4.3.3

	// find the URI tag
	const string uriTag = "URI(";
	string::size_type uriPos = line.find(uriTag.c_str());

	if(uriPos == string::npos) {
		uri.SetLength(0);
		return;
	}

	// find the matching end parenthesis
	const unsigned int start = uriPos + uriTag.size();
	unsigned int stop = start;

	const char* pBuffer = line.data();
	unsigned int lineLen = line.size();

	if(stop < lineLen) {
		while(pBuffer[stop] != ')') {
			stop++;
			if(stop == lineLen) break;
		}
	}

	if(start == stop) {
		cout << "ERROR: could not parse genome assembly ID from FASTA header." << endl;
		cout << "       " << line << endl;
		exit(1);
	}

	uri = line.substr(start, stop - start).c_str();

#endif
}
