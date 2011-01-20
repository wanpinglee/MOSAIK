// ***************************************************************************
// CTimeSupport - handles high-resolution timekeeping tasks.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#ifdef WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#endif
#include "Mosaik.h"

using namespace std;

class CTimeSupport {
public:
	// converts the supplied string into a 64-bit time object
	static uint64_t ConvertStringToTime(const string& time);
	// returns the current time
	static uint64_t GetSystemTime(void);
	// displays the date associated with the supplied parameter
	static string ConvertTimeToString(const uint64_t& currentTime);
private:
	// get the offset between local time and UTC
	static time_t GetUtcOffset(void);
	// converts a time_t variable to our 64-bit notation
	static uint64_t ConvertTimeT(const time_t& time);
	// stores a static list of the weekdays
	static const char* WEEKDAYS[];
	// stores a static list of the months
	static const char* MONTHS[];
	//
	static const uint64_t SECS_BETWEEN_EPOCHS;
	//
	static const uint64_t SECS_TO_100NS;
};
