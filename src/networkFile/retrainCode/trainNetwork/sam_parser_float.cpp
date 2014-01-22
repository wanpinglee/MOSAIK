#include "sam_parser_float.h"

#include <iostream>
#include <sstream>
//#include <string>

using namespace std;

template<typename T>
bool convert_from_string(const std::string& s, T& r) {
        std::istringstream iss(s);
        iss >> r;

        return (iss.fail() || ((std::size_t) iss.tellg()) != s.size()) ? false : true;
}

bool IsPositionCorrect(
    const string& mapped_chr, 
    const int& mapped_position, 
    const string& correct_position) {
  size_t found1, found2;
  found1 = correct_position.find(':');
  found1 = correct_position.find( ':', found1 + 1 );
  found2 = correct_position.find( ';', found1 + 1 );
  string chr = correct_position.substr(found1 + 1, found2 - found1 - 1);
  int pos;
  convert_from_string(correct_position.substr(found2 + 1), pos);

  if (mapped_chr != chr) return false;
  
  if ( ( ( pos - 20 ) < mapped_position ) && ( mapped_position < ( pos + 20 ) ) )
    return true;
  else
    return false;

}

bool SamParser::LoadNextNeuralElement(NeuralElement* element, Record* record) {
  string line;
  getline (input_sam, line);
  
  if (input_sam.eof()) return false;
  
  std::istringstream buffer (line);

  buffer >> record->name;
  buffer >> record->flag;
  buffer >> record->rname;
  buffer >> record->pos;
  buffer >> record->mapq;
  buffer >> record->cigar;
  buffer >> record->mrnm;
  buffer >> record->mpos;
  buffer >> record->isize;
  buffer >> record->seq;
  buffer >> record->qual;
  buffer >> record->rg;
  buffer >> record->nm;
  buffer >> record->md;
  buffer >> record->za;
  buffer >> record->zn;
  buffer >> record->xc;
/*
cout << record->name << "\t"
     << record->flag << "\t"
     << record->rname << "\t"
     << record->pos << "\t"
     << record->mapq << "\t"
     << record->cigar << "\t"
     << record->mrnm << "\t"
     << record->mpos << "\t"
     << record->isize << "\t"
     << record->seq << endl;
*/  
  int read_length = record->seq.size();
  if (read_length == 0) return false;
  
  size_t found1, found2;
  float best_sw, next_sw;
  found1 = record->zn.find(':');
  found1 = record->zn.find(':', found1 + 1);
  found2 = record->zn.find(';');
  convert_from_string(record->zn.substr(found1 + 1, found2 - found1 - 1), best_sw);

  found1 = found2;
  found2 = record->zn.find(';', found1 + 1);
  convert_from_string(record->zn.substr(found1 + 1, found2 - found1 - 1), next_sw);
  element->smith_waterman = (best_sw - next_sw) / static_cast<float>((read_length * 10));
  element->smith_waterman2 = best_sw;
  element->smith_waterman3 = next_sw;

  found1 = found2;
  found2 = record->zn.find(';', found1 + 1);
  convert_from_string(record->zn.substr(found1 + 1, found2 - found1 - 1), element->longest_match);
  element->longest_match = element->longest_match / static_cast<float>(read_length);

  found1 = found2;
  found2 = record->zn.find(';', found1 + 1 );
  convert_from_string(record->zn.substr( found1 + 1, found2 - found1 - 1 ), element->entropy);

  found1 = found2;
  found2 = record->zn.find(';', found1 + 1 );
  convert_from_string(record->zn.substr( found1 + 1, found2 - found1 - 1 ), element->num_mappings);

  //if (element->num_mappings == 1)
  //  element->smith_waterman = 1.0;
  //else if (element->num_mappings == 0)
  //  element->smith_waterman = -1.0;

  found1 = found2;
  convert_from_string(record->zn.substr(found1 + 1), element->num_hashes);

  element->correct = IsPositionCorrect(record->rname, record->pos, record->xc);

  return true;
}
