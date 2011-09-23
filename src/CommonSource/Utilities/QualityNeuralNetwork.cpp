#include "QualityNeuralNetwork.h"

const int SwMatchScore = 10;
const unsigned char PHRED_MAX = 255;

inline unsigned char float2phred(long double prob) {
    if (prob == 1)
    return PHRED_MAX;  // guards against "-0"
    long double p = -10 * (long double) log10(prob);
    if (p < 0 || p > PHRED_MAX) // int overflow guard
      return 0;
    else
      return floor(p + 0.5);

}

QualityNeuralNetwork::~QualityNeuralNetwork() {
  ann_open = false;
  fann_destroy(pe_ann);
  fann_destroy(se_ann);
}

void QualityNeuralNetwork::Open(const string& pe_file, 
                                const string& se_file) {
  pe_ann = fann_create_from_file(pe_file.c_str());
  se_ann = fann_create_from_file(se_file.c_str());
  ann_open = true;
}

unsigned char QualityNeuralNetwork::GetQualitySe(const FannInputs& annInputs) {
  if (!ann_open) return 0;

  fann_inputs.clear();
  
  int swDiff = annInputs.swScore - annInputs.nextSwScore;
  fann_inputs.push_back(swDiff / (float)(annInputs.read_length * SwMatchScore));
  fann_inputs.push_back((float)annInputs.longest_match);
  fann_inputs.push_back(annInputs.entropy);
  fann_inputs.push_back(log10((float)(annInputs.numMappings + 1)));
  fann_inputs.push_back(log10((float)(annInputs.numHashes + 1)));
  
  calc_out = fann_run(se_ann, &fann_inputs[0]);
  return float2phred(1 - (1 + calc_out[0]) / 2);

}

unsigned char QualityNeuralNetwork::GetQualityPe(const FannInputs& annInputs1,
                                                 const FannInputs& annInputs2,
						 const int& fragment_length_diff) {
  if (!ann_open) return 0;

  fann_inputs.clear();

  int swDiff1 = annInputs1.swScore - annInputs1.nextSwScore;
  fann_inputs.push_back(swDiff1 / (float)(annInputs1.read_length * SwMatchScore));
  fann_inputs.push_back(annInputs1.longest_match / (float)annInputs1.read_length);
  fann_inputs.push_back(annInputs1.entropy);
  if (annInputs1.numHashes == 0) { // the mate is rescued by mate2
    fann_inputs.push_back(log10((float)(annInputs2.numMappings + 1)));
    fann_inputs.push_back(log10((float)(annInputs2.numHashes + 1)));
  } else {
    fann_inputs.push_back(log10((float)(annInputs1.numMappings + 1)));
    fann_inputs.push_back(log10((float)(annInputs1.numHashes + 1)));
  }

  int swDiff2 = annInputs2.swScore - annInputs2.nextSwScore;
  fann_inputs.push_back(swDiff2 / (float)(annInputs2.read_length * SwMatchScore));
  fann_inputs.push_back(annInputs2.longest_match / (float)annInputs2.read_length);
  fann_inputs.push_back(annInputs2.entropy);
  if (annInputs2.numHashes == 0) { // the mate is rescued by mate1
    fann_inputs.push_back(log10((float)(annInputs1.numMappings + 1)));
    fann_inputs.push_back(log10((float)(annInputs1.numHashes + 1)));
  } else {
    fann_inputs.push_back(log10((float)(annInputs2.numMappings + 1)));
    fann_inputs.push_back(log10((float)(annInputs2.numHashes + 1)));
  }

  fann_inputs.push_back(log10((float)(fragment_length_diff + 1)));

  calc_out = fann_run(pe_ann, &fann_inputs[0]);
  return float2phred(1 - (1 + calc_out[0]) / 2);
}

