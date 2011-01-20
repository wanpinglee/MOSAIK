// ***************************************************************************
// COptions - parses command-line arguments and creates the help menu.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "ConsoleUtilities.h"
#include "Variant.h"

#ifndef WIN32
#include <stdint.h>
#endif

using namespace std;

#define ARGUMENT_LENGTH       35
#define DESC_LENGTH_FIRST_ROW 42
#define DESC_LENGTH           39
#define MAX_LINE_LENGTH       78

#ifdef WIN32
#define snprintf _snprintf
typedef __int64          int64_t;
typedef unsigned __int64 uint64_t;
#define strtoui64 _strtoui64
#else
#define strtoui64 strtoull
#endif

struct Option {
	string Argument;
	string ValueDescription;
	string Description;
	bool StoreValue;
	bool HasDefaultValue;
	variant_t DefaultValue;

	// constructor
	Option(void)
		: StoreValue(true)
		, HasDefaultValue(false)
	{}
};

struct OptionValue {
	bool* pFoundArgument;
	void* pValue;
	string ValueTypeDescription;
	bool UseVector;
	bool StoreValue;
	bool IsRequired;
	variant_t VariantValue;

	// constructor
	OptionValue(void)
		: pFoundArgument(NULL)
		, pValue(NULL)
		, UseVector(false)
		, StoreValue(true)
		, IsRequired(false)
	{}

};

struct OptionGroup {
	string Name;
	vector<Option> Options;
};

class COptions {
public:
	// adds a simple option to the parser
	static void AddOption(const string& argument, const string& optionDescription, bool& foundArgument, OptionGroup* pGroup);
	// adds a value option to the parser
	template<typename T>
	static void AddValueOption(const string& argument, const string& valueDescription, const string& optionDescription, const string& valueTypeDescription, bool& foundArgument, T& val, OptionGroup* pGroup);
	// adds a value option to the parser (default value)
	template<typename T, typename D>
	static void AddValueOption(const string& argument, const string& valueDescription, const string& optionDescription, const string& valueTypeDescription, bool& foundArgument, T& val, OptionGroup* pGroup, D& defaultValue);
	// creates an option group
	static OptionGroup* CreateOptionGroup(const string& groupName);
	// displays the help menu
	static void DisplayHelp(void);
	// parses the command line
	static void Parse(int argc, char* argv[]);
	// sets the program info
	static void SetProgramInfo(const string& programName, const string& description, const string& arguments);

private:
	// the program name
	static string mProgramName;
	// the main description
	static string mDescription;
	// the example arguments
	static string mExampleArguments;
	// stores the option groups
	static vector<OptionGroup> mOptionGroups;
	// stores the options in a map
	static map<string, OptionValue> mOptionsMap;
};


// adds a value option to the parser
template<typename T>
void COptions::AddValueOption(const string& argument, const string& valueDescription, const string& optionDescription, const string& valueTypeDescription, bool& foundArgument, T& val, OptionGroup* pGroup) {

	Option o;
	o.Argument         = argument;
	o.ValueDescription = valueDescription;
	o.Description      = optionDescription;
	pGroup->Options.push_back(o);

	OptionValue ov;
	ov.pFoundArgument       = &foundArgument;
	ov.pValue               = (void*)&val;
	ov.VariantValue         = val;
	ov.IsRequired           = (valueTypeDescription.empty() ? false : true);
	ov.ValueTypeDescription = valueTypeDescription;
	mOptionsMap[argument] = ov;
}

// adds a value option to the parser (default value)
template<typename T, typename D>
void COptions::AddValueOption(const string& argument, const string& valueDescription, const string& optionDescription, const string& valueTypeDescription, bool& foundArgument, T& val, OptionGroup* pGroup, D& defaultValue) {

	Option o;
	o.Argument         = argument;
	o.ValueDescription = valueDescription;
	o.Description      = optionDescription;
	o.DefaultValue     = defaultValue;
	o.HasDefaultValue  = true;
	pGroup->Options.push_back(o);

	OptionValue ov;
	ov.pFoundArgument       = &foundArgument;
	ov.pValue               = (void*)&val;
	ov.VariantValue         = val;
	ov.IsRequired           = (valueTypeDescription.empty() ? false : true);
	ov.ValueTypeDescription = valueTypeDescription;
	mOptionsMap[argument] = ov;
}
