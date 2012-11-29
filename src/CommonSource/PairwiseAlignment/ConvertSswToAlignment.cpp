#include <ConvertSswToAlignment.h>

bool ConvertSswToAlignment(
    const StripedSmithWaterman::Alignment& ssw_al,
    const char* ref,
    const char* query,
    Alignment* al) {
  al->QueryBegin  = ssw_al.query_begin;
  al->QueryEnd    = ssw_al.query_end;
  al->QueryLength = ssw_al.query_end - ssw_al.query_begin + 1;

  al->ReferenceBegin = ssw_al.ref_begin;
  al->ReferenceEnd   = ssw_al.ref_end;
  al->NumMismatches  = ssw_al.mismatches;
  al->SwScore        = ssw_al.sw_score;

  al->NumLongestMatchs = 0;

  char* ref_ptr   = (char*)ref + al->ReferenceBegin;
  char* query_ptr = (char*)query + al->QueryBegin;
  
  for (unsigned int i = 0; i < ssw_al.cigar.size(); ++i) {
    int op  = ssw_al.cigar[i] & 0x0000000f;
    int len = ssw_al.cigar[i] >> 4;
    switch (op) {
      case 0: //M
        al->Reference.Append(ref_ptr, len);
	al->Query.Append(query_ptr, len);
	ref_ptr   += len;
	query_ptr += len;
	break;
      case 1: //I
        al->Query.Append(query_ptr, len);
	query_ptr += len;
	al->Reference.Append('-', len);
        break;
      case 2: //D
        al->Reference.Append(ref_ptr, len);
	ref_ptr   += len;
	al->Query.Append('-', len);
        break;
    }
  }

  return true;

}

