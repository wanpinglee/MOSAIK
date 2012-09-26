#ifndef SAM_PARSER_H_
#define SAM_PARSER_H_

#include <fstream>
#include <string>

using std::ifstream;
using std::string;

struct Record {
  string name, rname, cigar, mrnm, seq, qual, rg, nm, md, za, zn, xc;
  unsigned int pos, mpos;
  int isize;
  unsigned short flag, mapq;
};

struct NeuralElement {
  float smith_waterman;
  float smith_waterman2;
  float smith_waterman3;
  float longest_match;
  float entropy;
  float num_mappings;
  float num_hashes;
  bool  correct;
};

class SamParser {
 public:
  bool LoadNextNeuralElement(NeuralElement* element, Record* record);


  // inlines
  inline bool Open(const char* filename) {
    filename = filename;
    input_sam.open(filename);
    return input_sam.good();
  }
  inline void Close() {
    input_sam.close();
  }
 private:
  string filename;
  ifstream input_sam;
};

#endif
