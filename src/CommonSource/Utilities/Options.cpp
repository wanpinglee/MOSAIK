// ***************************************************************************
// COptions - parses command-line arguments and creates the help menu.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "Options.h"

// the program name
string COptions::mProgramName;

// the main description
string COptions::mDescription;

// the example arguments
string COptions::mExampleArguments;

// stores the option groups
vector<OptionGroup> COptions::mOptionGroups;

// stores the options in a map
map<string, OptionValue> COptions::mOptionsMap;

// adds a simple option to the parser
void COptions::AddOption(const string& argument, const string& optionDescription, bool& foundArgument, OptionGroup* pGroup) {

	Option o;
	o.Argument    = argument;
	o.Description = optionDescription;
	o.StoreValue  = false;
	pGroup->Options.push_back(o);

	OptionValue ov;
	ov.pFoundArgument = &foundArgument;
	ov.StoreValue     = false;

	mOptionsMap[argument] = ov;
}

// creates an option group
OptionGroup* COptions::CreateOptionGroup(const string& groupName) {
	OptionGroup og;
	og.Name = groupName;
	mOptionGroups.push_back(og);
	return &mOptionGroups[mOptionGroups.size() - 1];
}

// displays the help menu
void COptions::DisplayHelp(void) {

	// initialize
	char argumentBuffer[ARGUMENT_LENGTH + 1];
	ostringstream sb;

	char indentBuffer[MAX_LINE_LENGTH - DESC_LENGTH + 1];
	memset(indentBuffer, ' ', MAX_LINE_LENGTH - DESC_LENGTH);
	indentBuffer[MAX_LINE_LENGTH - DESC_LENGTH] = 0;

	// display the menu
	printf("Description: %s.\n\n", mDescription.c_str());

	printf("Usage: ");
	CConsole::Bold(); printf("%s", mProgramName.c_str()); CConsole::Reset();
	printf(" %s\n\n", mExampleArguments.c_str());

	vector<Option>::const_iterator optionIter;
	vector<OptionGroup>::const_iterator groupIter;
	for(groupIter = mOptionGroups.begin(); groupIter != mOptionGroups.end(); ++groupIter) {
		CConsole::Heading(); printf("%s:\n", groupIter->Name.c_str()); CConsole::Reset();

		for(optionIter = groupIter->Options.begin(); optionIter != groupIter->Options.end(); ++optionIter) {

			if(optionIter->StoreValue) snprintf(argumentBuffer, ARGUMENT_LENGTH + 1, "  %s <%s>", optionIter->Argument.c_str(), optionIter->ValueDescription.c_str());
			else snprintf(argumentBuffer, ARGUMENT_LENGTH + 1, "  %s", optionIter->Argument.c_str());
			printf("%-35s ", argumentBuffer);

			string description = optionIter->Description;

			// handle default values
			if(optionIter->HasDefaultValue) {
				sb.str("");
				sb << description << ". def: ";

				if(optionIter->DefaultValue.is_type<unsigned int>()) {
					sb << (unsigned int)optionIter->DefaultValue;
				} else if(optionIter->DefaultValue.is_type<unsigned char>()) {
					sb << (unsigned short)(unsigned char)optionIter->DefaultValue;
				} else if(optionIter->DefaultValue.is_type<float>()) {
					sb << fixed << setprecision(2) << (float)optionIter->DefaultValue;
				} else if(optionIter->DefaultValue.is_type<double>()) {
					sb << fixed << setprecision(4) << (double)optionIter->DefaultValue;
				} else if(optionIter->DefaultValue.is_type<string>()) {
					const string stringValue = optionIter->DefaultValue;
					sb << stringValue;
				} else {
					printf("ERROR: Found an unsupported data type for argument %s when casting the default value.\n", optionIter->Argument.c_str());
					exit(1);
				}

				description = sb.str(); 
			}

			if(description.size() <= DESC_LENGTH_FIRST_ROW) {

				printf("%s\n", description.c_str());

			} else {

				// handle the first row
				const char* pDescription = description.data();
				unsigned int cutIndex = DESC_LENGTH_FIRST_ROW;
				while(pDescription[cutIndex] != ' ') cutIndex--;
				printf("%s\n", description.substr(0, cutIndex).c_str());
				description = description.substr(cutIndex + 1);

				// handle subsequent rows
				while(description.size() > DESC_LENGTH) {
					pDescription = description.data();
					cutIndex = DESC_LENGTH;
					while(pDescription[cutIndex] != ' ') cutIndex--;
					printf("%s%s\n", indentBuffer, description.substr(0, cutIndex).c_str());
					description = description.substr(cutIndex + 1);
				}

				// handle last row
				printf("%s%s\n", indentBuffer, description.c_str());
			}			
		}

		printf("\n");
	}

	CConsole::Heading(); printf("Help:\n"); CConsole::Reset(); 
	printf("  --help, -h                        shows this help text\n");
	exit(1);
}

// parses the command line
void COptions::Parse(int argc, char* argv[]) {

	// initialize
	CConsole::Initialize();
	map<string, OptionValue>::const_iterator ovMapIter, checkMapIter;
	const int LAST_INDEX = argc - 1;
	ostringstream errorBuilder;
	bool foundError = false;
	char* end_ptr = NULL;

	const string ERROR_SPACER(7, ' ');

	// check if we should show the help menu
	bool showHelpMenu = false;
	if(argc > 1) {
		for(int i = 1; i < argc; i++) {
			const string argument = argv[i];
			if((argument == "-h") || (argument == "--help")) showHelpMenu = true;
		}
	} else showHelpMenu = true;

	if(showHelpMenu) DisplayHelp();

	// check each argument
	for(int i = 1; i < argc; i++) {
		const string argument = argv[i];
		ovMapIter = mOptionsMap.find(argument);

		if(ovMapIter == mOptionsMap.end()) {
			errorBuilder << ERROR_SPACER << "An unrecognized argument was found: " << argument << endl;
			foundError = true;

		} else {

			*ovMapIter->second.pFoundArgument = true;

			// grab the value
			if(ovMapIter->second.StoreValue) {

				if(i < LAST_INDEX) {

					// check if the next argument is really a command line option
					const string val = argv[i + 1]; 
					checkMapIter = mOptionsMap.find(val);

					if(checkMapIter == mOptionsMap.end()) {
						++i;
						if(ovMapIter->second.VariantValue.is_type<unsigned int>()) {

							const unsigned int uint32 = (unsigned int)strtoul(val.c_str(), &end_ptr, 10);
							unsigned int* varValue = (unsigned int*)ovMapIter->second.pValue;
							*varValue = uint32;

						} else if(ovMapIter->second.VariantValue.is_type<unsigned char>()) {

							const unsigned char uint8 = (unsigned char)strtoul(val.c_str(), &end_ptr, 10);
							unsigned char* varValue = (unsigned char*)ovMapIter->second.pValue;
							*varValue = uint8;

						} else if(ovMapIter->second.VariantValue.is_type<uint64_t>()) {

							const uint64_t uint64 = strtoui64(val.c_str(), &end_ptr, 10);
							uint64_t* varValue = (uint64_t*)ovMapIter->second.pValue;
							*varValue = uint64;

						} else if(ovMapIter->second.VariantValue.is_type<double>()) {

							const double d = strtod(val.c_str(), &end_ptr);
							double* varValue = (double*)ovMapIter->second.pValue;
							*varValue = d;

						} else if(ovMapIter->second.VariantValue.is_type<float>()) {

							const float f = (float)strtod(val.c_str(), &end_ptr);
							float* varValue = (float*)ovMapIter->second.pValue;
							*varValue = f;

						} else if(ovMapIter->second.VariantValue.is_type<string>()) {

							string* pStringValue = (string*)ovMapIter->second.pValue;
							*pStringValue = val;

						} else if(ovMapIter->second.VariantValue.is_type<vector<string> >()) {

							vector<string>* pVectorValue = (vector<string>*)ovMapIter->second.pValue;
							pVectorValue->push_back(val);

						} else {

							printf("ERROR: Found an unsupported data type for argument %s when parsing the arguments.\n", argument.c_str());
							exit(1);
						}

					} else {
						errorBuilder << ERROR_SPACER << "The argument (" << argument << ") expects a value, but none was found." << endl;
						foundError = true;
					}

				} else {
					errorBuilder << ERROR_SPACER << "The argument (" << argument << ") expects a value, but none was found." << endl;
					foundError = true;
				}
			}
		}
	}

	// check if we missed any required parameters
	for(ovMapIter = mOptionsMap.begin(); ovMapIter != mOptionsMap.end(); ++ovMapIter) {
		if(ovMapIter->second.IsRequired && !*ovMapIter->second.pFoundArgument) {
			errorBuilder << ERROR_SPACER << ovMapIter->second.ValueTypeDescription << " was not specified. Please use the " << ovMapIter->first << " parameter." << endl;
			foundError = true;
		}
	}

	// print the errors if any were found
	if(foundError) {

		CConsole::Red();
		printf("ERROR: Some problems were encountered when parsing the command line options:\n");
		CConsole::Reset();

		printf("%s\n", errorBuilder.str().c_str());
		printf("For a complete list of command line options, type \"%s -h\"\n", argv[0]);
		exit(1);
	}
}

// sets the program info
void COptions::SetProgramInfo(const string& programName, const string& description, const string& arguments) {
	mProgramName      = programName;
	mDescription      = description;
	mExampleArguments = arguments;
}
