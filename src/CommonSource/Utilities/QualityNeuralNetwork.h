#ifndef COMMONSOURCE_UTILITIES_QUALITYNEURALNETWORK_H_ 
#define COMMONSOURCE_UTILITIES_QUALITYNEURALNETWORK_H_

#include <string>
#include <vector>

#include "fann.h"

using std::string;
using std::vector;

class QualityNeuralNetwork {
 public:
  struct FannInputs {
    int read_length;
    float swScore;
    float nextSwScore;
    int longest_match;
    float entropy;
    int numMappings;
    int numHashes;
  };
  QualityNeuralNetwork()
      : pe_ann(NULL)
      , se_ann(NULL)
      , ann_open(false)
      , fann_inputs()
      , calc_out(NULL)
      {};
  ~QualityNeuralNetwork();
  void Open(const string& pe_file, const string& se_file);
  unsigned char GetQualitySe(const FannInputs& annInputs);
  unsigned char GetQualityPe(const FannInputs& annInputs1,
                             const FannInputs& annInputs2,
			     const int& fragment_length_diff);
 private:
  // variables
  fann*  pe_ann;
  fann*  se_ann;
  bool   ann_open;
  vector<fann_type> fann_inputs;
  fann_type *calc_out;
  // end of variables
  
  // Functions
  QualityNeuralNetwork (const QualityNeuralNetwork&);
  QualityNeuralNetwork& operator= (const QualityNeuralNetwork&);
  // end of Functions
};

#endif // COMMONSOURCE_UTILITIES_QUALITYNEURALNETWORK_H_
