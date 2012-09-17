#ifndef ParameterParser_H_
#define ParameterParser_H_

#include <string>

using std::string;

struct Parameters {
  // i/o parameters
  string input_sam;             // -i  --input
  string output_ann;            // -o  --output

  string reference_filename;
  string hash_filename;

  // operation parameters
  int   fragment_length;  // -f --fragmenr-length
  int   network_layer;    // -l --network-layer
  int   hidden_neurons;   // -n --hidden-neurons
  float desired_error;    // -e --desired-error
  int   max_epochs;       // -m --max-epochs
  bool  paired_end;       // -p --paired-end
	
  // command line
  string command_line;

  // default values
  Parameters()
      : input_sam()
      , output_ann()
      , fragment_length(0)
      , network_layer(3)
      , hidden_neurons(100)
      , desired_error(0.01)
      , max_epochs(10000)
      , paired_end(false)
  {}
};

void ParseArgumentsOrDie(const int argc, char* const * argv, Parameters* param);

#endif
