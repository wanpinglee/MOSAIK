#include "fann.h"
#include "convert.h"
#include <fstream>
#include <iostream>
#include <getopt.h>
#include <limits.h>
#include <string>
#include <iostream>
#include <vector>

#include "sam_parser_float.h"
#include "parameter_parser_float.h"

using namespace std;
//using namespace vcf;

template<typename T>
bool convert_from_string(const std::string& s, T& r) {
        istringstream iss(s);
        iss >> r;

        return (iss.fail() || ((std::size_t) iss.tellg()) != s.size()) ? false : true;
}

struct fann_train_data *read_from_array(vector<float>& din, vector<float>& dout, unsigned int num_data, unsigned int num_input, unsigned int num_output) {
  unsigned int i, j;
  fann_type *data_input, *data_output;
  struct fann_train_data *data =
    (struct fann_train_data *) malloc(sizeof(struct fann_train_data));
  if(data == NULL) {
    fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
    return NULL;
  }

  fann_init_error_data((struct fann_error *) data);

  data->num_data = num_data;
  data->num_input = num_input;
  data->num_output = num_output;
  data->input = (float **) calloc(num_data, sizeof(float *));
  if(data->input == NULL) {
    fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
    fann_destroy_train(data);
    return NULL;
  }

  data->output = (float **) calloc(num_data, sizeof(float *));
  if(data->output == NULL) {
    fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
    fann_destroy_train(data);
    return NULL;
  }

  data_input = (float *) calloc(num_input * num_data, sizeof(float));
  if(data_input == NULL) {
    fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
    fann_destroy_train(data);
    return NULL;
  }

  data_output = (float *) calloc(num_output * num_data, sizeof(float));
  if(data_output == NULL) {
    fann_error(NULL, FANN_E_CANT_ALLOCATE_MEM);
    fann_destroy_train(data);
    return NULL;
  }

  for(i = 0; i != num_data; i++) {
    data->input[i] = data_input;
    data_input += num_input;
   
    for(j = 0; j != num_input; j++) {
      data->input[i][j] = din[i*num_input+j];
    }
   
   
    data->output[i] = data_output;
    data_output += num_output;
   
    for(j = 0; j != num_output; j++) {
      data->output[i][j] = dout[i*num_output+j];
    }
  }
  return data;
}

void convertSamToFann(const NeuralElement& element, vector<float>* din) {

  //din->push_back(element.smith_waterman2); // best_sw
  din->push_back(element.smith_waterman);  // (best-next_best) / read_len
  din->push_back(element.longest_match);
  din->push_back(element.entropy);
  din->push_back(log10(element.num_mappings + 1));
  din->push_back(log10(element.num_hashes + 1));
}

int main(int argc, char** argv)
{   
    Parameters param;
    ParseArgumentsOrDie(argc, argv, &param);

    unsigned int num_output = 1;
    unsigned int num_input;

    if (param.paired_end)
      num_input = 11;
    else
      num_input = 5;

    vector<float> din, dout;

    //unsigned int num_layers = 3;
    //unsigned int num_neurons_hidden = 100;
    //float desired_error = (float) 0.005;
    //unsigned int max_epochs = 10000;
    unsigned int epochs_between_reports = 10;

    string annFile = param.output_ann;

    SamParser sam_parser;
    if (!sam_parser.Open(param.input_sam.c_str())) {
    	cout << "ERROR: Cannot open input files." << endl;
	return 1;
    }

    Record mate1, mate2;
    NeuralElement mate1_element, mate2_element;
    bool sam_okay = false;
    sam_okay = sam_parser.LoadNextNeuralElement(&mate1_element, &mate1);
    if (param.paired_end)
      sam_okay = sam_parser.LoadNextNeuralElement(&mate2_element, &mate2);
    while (sam_okay) {
	// first mate
	// al info
	convertSamToFann(mate1_element, &din);
	if (param.paired_end) {
	  // mate of al info
	  convertSamToFann(mate2_element, &din);
	  int fl = abs(mate1.isize);
	  fl = abs(param.fragment_length - fl);
	  din.push_back(log10(static_cast<float>(fl + 1)));
	}

	float output = mate1_element.correct ? 1.0 : -1.0;
	dout.push_back(output);

	if (param.paired_end) {
	  // second mate
	  // al info
	  convertSamToFann(mate2_element, &din);
	  // mate of al info
	  convertSamToFann(mate1_element, &din);
	  int fl = abs(mate2.isize);
	  fl = abs(param.fragment_length - fl);
	  din.push_back(log10(static_cast<float>(fl + 1)));

	  output = mate2_element.correct ? 1.0 : -1.0;
	  dout.push_back(output);
	}

	// load next
	sam_okay = sam_parser.LoadNextNeuralElement(&mate1_element, &mate1);
	if (param.paired_end)
	  sam_okay = sam_parser.LoadNextNeuralElement(&mate2_element, &mate2);
    }


    cout << "Eating sam is done." << endl;

    struct fann *ann = fann_create_standard(param.network_layer, num_input, param.hidden_neurons, num_output);

    fann_set_activation_function_hidden(ann, FANN_SIGMOID_SYMMETRIC);
    fann_set_activation_function_output(ann, FANN_SIGMOID_SYMMETRIC);

    unsigned int num_data = din.size() / num_input;
    //if (din.size() % num_input != 0) cerr << "what???" << endl;
    fann_train_data *data = read_from_array(din, dout, num_data, num_input, num_output);

    //fann_train_on_file(ann, "test.data", max_epochs, epochs_between_reports, desired_error);
    fann_train_on_data(ann, data, param.max_epochs, epochs_between_reports, param.desired_error);

    cerr << "saving neural net to " << param.output_ann << endl;
    fann_save(ann, param.output_ann.c_str());

    fann_destroy(ann);
    sam_parser.Close();
    return 0;
}
