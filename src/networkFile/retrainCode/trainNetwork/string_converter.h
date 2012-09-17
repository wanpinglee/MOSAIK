
#ifndef UTILITIES_MISCELLANEOUS_STRING_CONVERTER_H_
#define UTILITIES_MISCELLANEOUS_STRING_CONVERTER_H_

#include <sstream>

// converts the string into the specified type, setting r to the converted
// value and returning true/false on success or failure
template<typename T>
bool convert_from_string(const std::string& s, T& r) {
	std::istringstream iss(s);
	iss >> r;
	
	return (iss.fail() || ((std::size_t) iss.tellg()) != s.size()) ? false : true;
}


//template<typename T>
//std::string convert(const T& r) {
//	std::ostringstream iss;
//	iss << r;
//	
//	return iss.str();
//}

#endif // UTILITIES_STRING_CONVERTER_H_; utilities/string_converter.h
