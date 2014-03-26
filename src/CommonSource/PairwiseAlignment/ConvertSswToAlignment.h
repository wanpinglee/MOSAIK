#include <Alignment.h>
#include <ssw_cpp.h>

bool ConvertSswToAlignment(
    const StripedSmithWaterman::Alignment& ssw_al,
    const char* ref,
    const char* query,
    const int& queryLength,
    Alignment* al);
