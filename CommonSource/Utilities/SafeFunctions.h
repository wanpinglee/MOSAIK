// ***************************************************************************
// SafeFunctions.h - provides platform independence while using the "safe
//                   function" syntax from WIN32.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#ifndef WIN32

#include <cstdio>
#include <cstdarg>
#include <cstring>

typedef int errno_t;
#define EINVAL 22

inline int sprintf_s(char *buffer, size_t sizeOfBuffer, const char *format, ...) {
	va_list argp;
	va_start(argp, format);
	unsigned int err = vsprintf(buffer, format, argp);
	va_end(argp);
	return err;
}

inline errno_t strncpy_s(char *strDest, size_t sizeInBytes, const char *strSource, size_t count) {
	strncpy(strDest, strSource, count);
	return 0;
}

inline errno_t strcat_s(char *strDestination, size_t sizeInBytes, const char *strSource) {
	strcat(strDestination, strSource);
	return 0;
}

inline errno_t memmove_s(void *dest, size_t numberOfElements, const void *src, size_t count) {
	memmove(dest, src, count);
	return 0;
}

inline errno_t fopen_s(FILE** pFile, const char *filename, const char *mode) {
	*pFile = fopen(filename, mode);
	if(*pFile) return 0;
	return EINVAL;
}

inline errno_t localtime_s(struct tm* _tm, const time_t *time) {
	memcpy(_tm, localtime(time), sizeof(*_tm));
	return 0;
}

inline errno_t gmtime_s(struct tm* _tm, const time_t *time) {
	memcpy(_tm, gmtime(time), sizeof(*_tm));
	return 0;
}

inline errno_t asctime_s(char* buffer, size_t sizeInBytes, const struct tm *_tm) {
	char* time = asctime(_tm);
	unsigned int timeLen = strlen(time);
	strncpy(buffer, time, timeLen);
	buffer[timeLen] = 0;
	return 0;
}

inline char* strtok_s(char* strToken, const char *strDelimit, char **context) {
	return strtok(strToken, strDelimit);
}

#endif
