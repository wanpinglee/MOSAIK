// ***************************************************************************
// ConversionUtilities.h - converts strings into numerical data types.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <cstdio>

inline unsigned long long GetUInt64(char* s) {
	char* end_ptr = NULL;
	unsigned long long ull = (unsigned long long)strtoul(s, &end_ptr, 10);

	if(s == end_ptr) {
		printf("ERROR: Could not convert the string (%s) to an unsigned long long.\n", s);
		exit(1);
	}

	return ull;
}

inline unsigned int GetUnsignedInt(char* s) {
	char* end_ptr = NULL;
	unsigned int ui = (unsigned int)strtoul(s, &end_ptr, 10);

	if(s == end_ptr) {
		printf("ERROR: Could not convert the string (%s) to an unsigned integer.\n", s);
		exit(1);
	}

	return ui;
}

inline unsigned short GetUnsignedShort(char* s) {
	char* end_ptr = NULL;
	unsigned short us = (unsigned short)strtoul(s, &end_ptr, 10);

	if(s == end_ptr) {
		printf("ERROR: Could not convert the string (%s) to an unsigned short.\n", s);
		exit(1);
	}

	return us;
}

inline unsigned char GetUnsignedChar(char* s) {
	char* end_ptr = NULL;
	unsigned char uc = (unsigned char)strtoul(s, &end_ptr, 10);

	if(s == end_ptr) {
		printf("ERROR: Could not convert the string (%s) to an unsigned character.\n", s);
		exit(1);
	}

	return uc;
}

inline double GetDouble(char* s) {
	char* end_ptr = NULL;
	double d = strtod(s, &end_ptr);

	if(s == end_ptr) {
		printf("ERROR: Could not convert the string (%s) to a double.\n", s);
		exit(1);
	}

	return d;
}

inline float GetFloat(char* s) {
	char* end_ptr = NULL;
	float f = (float)strtod(s, &end_ptr);

	if(s == end_ptr) {
		printf("ERROR: Could not convert the string (%s) to a float.\n", s);
		exit(1);
	}

	return f;
}
