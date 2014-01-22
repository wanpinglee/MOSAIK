#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>

#include <iostream>
#include <string>

#include "string_converter.h"
#include "parameter_parser_float.h"

using std::cout;
using std::endl;
using std::string;

bool CheckParameters(Parameters* param);
void PrintHelp(const string& program);

void ParseArgumentsOrDie(const int argc, char* const * argv, 
    Parameters* param) {

	if (argc == 1) { // no argument
		PrintHelp(argv[0]);
		exit(1);
	}

	// record command line
	param->command_line = argv[0];
	for ( int i = 1; i < argc; ++i ) {
		param->command_line += " ";
		param->command_line += argv[i];
	}

	const char *short_option = "hi:o:f:l:n:e:m:p";

	const struct option long_option[] = {
		{ "help", no_argument, NULL, 'h'},
		{ "input", required_argument, NULL, 'i'},
		{ "output", required_argument, NULL, 'o'},

		{ "fragment-length", no_argument, NULL, 'f'},
		{ "network-layer", no_argument, NULL, 'l'},
		{ "hidden-neurons", no_argument, NULL, 'n'},
		{ "desired-error", no_argument, NULL, 'e'},
		{ "max-epochs", no_argument, NULL, 'm'},
		{ "paired-end", no_argument, NULL, 'p'},

		{ 0, 0, 0, 0 }
	};

	int c = 0;
	bool help = false;
	while (true) {
		int optionIndex = 0;
		c = getopt_long( argc, argv, short_option, long_option, &optionIndex );
		
		if ( c == -1 ) // end of options
			break;
		
		switch ( c ) {
			// help
			case 'h':
				help = true;
				break;
			// i/o parameters
			case 'i':
				param->input_sam =  optarg;
				break;
			case 'o':
				param->output_ann = optarg;
				break;
			// operation parameters
			case 'f':
				if (!convert_from_string( optarg, param->fragment_length))
					cout << "WARNING: Cannot parse -f --fragment-length." << endl;
				break;
			case 'l':
				if (!convert_from_string( optarg, param->network_layer))
					cout << "WARNING: Cannot parse -l --network-layer." << endl;
				break;
			case 'n':
				if (!convert_from_string( optarg, param->hidden_neurons))
					cout << "WARNING: Cannot parse -n --hidden-neurons." << endl;
				break;

			case 'e':
				if (!convert_from_string( optarg, param->desired_error))
					cout << "WARNING: Cannot parse -e --desired-error." << endl;
				break;

			case 'm':
				if (!convert_from_string( optarg, param->max_epochs))
					cout << "WARNING: Cannot parse -m --max-epochs." << endl;
				break;
			case 'p':
				param->paired_end = true;
				break;
			default:
				break;
		}

	}

	if (help) {
		PrintHelp(argv[0]);
		exit(1);
	}

	if (!CheckParameters(param)) exit(1);
}

// true: passing the checker
bool CheckParameters(Parameters* param) {
	
	bool errorFound = false;
	// necessary parameters
	if ( param->input_sam.empty() ) {
		cout << "ERROR: Please specific an input file, -i." << endl;
		errorFound = true;
	}
	
	if ( param->output_ann.empty() ) {
		cout << "ERROR: Please specific an output file, -o." << endl;
		errorFound = true;
	}
	
	// unnecessary parameters
	if ( param->network_layer == 0 ) {
		cout << "WARNING: -l should not be zero. Set it to default, 3." << endl;
		param->network_layer = 3;
	}

	if ( param->hidden_neurons == 0 ) {
		cout << "WARNING: -n should not be zero. Set it to default, 100." << endl;
		param->hidden_neurons = 100;
	}

	if ( ( param->desired_error < 0.0 ) || ( param->desired_error > 1.0 ) ) {
		cout << "WARNING: -e should be in [0.0 - 1.0]. Set it to default, 0.01." << endl;
		param->desired_error = 0.01;
	}

	if ( param->max_epochs == 0 ) {
		cout << "WARNING: -m should not be zero. Set it to default, 10000." << endl;
		param->max_epochs = 10000;
	}

	if ( param->paired_end && (param->fragment_length == 0 )) {
		cout << "ERROR: Please specific fragment length, -f." << endl;
		errorFound = true;
	}

	if ( !param->paired_end && (param->fragment_length > 0 )) {
		cout << "WARNING: Unused setting -f." << endl;
		errorFound = true;
	}

	return !errorFound;

}

void PrintHelp(const string& program) {
	cout
		<< endl
		<< "usage: " << program << " [OPTIONS] -i <FILE> -o <FILE>"
		<< endl
		<< endl
		<< "Help:" << endl
		<< endl
		<< "   -h --help             Print this help dialog." << endl
		<< endl
		<< "Input & Output:" << endl
		<< endl
		<< "   -i --input <FILE>     Input SAM file that includes XC tags." << endl
		<< "   -o --output <FILE>    Output ANN file." << endl
		<< endl
		<< "Operations" << endl
		<< endl
		<< "   -f --fragment-length <INT>" << endl
		<< "                         Fragment length." << endl
		<< "   -l --network-layer <INT>" << endl
		<< "                         Network layer." << endl
		<< "                         Default: 3" << endl
		<< "   -n --hidden-neurons <INT>"  << endl
		<< "                         Hidden neurons." << endl
		<< "                         Default: 100" << endl
		<< "   -e --desired-error <FLOAT>"  << endl
		<< "                         Desired error [0.0 - 1.0]." << endl
		<< "                         Default: 0.01" << endl
		<< "   -m --max-epochs <INT>"  << endl
		<< "                         Max epoch." << endl
		<< "                         Default: 10000" << endl
		<< "   -p --paired-end       Default: false"  << endl
		<< endl;
}

