// ***************************************************************************
// CTimeSupport - handles high-resolution timekeeping tasks.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "TimeSupport.h"

const char* CTimeSupport::WEEKDAYS[] = { "Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat" };
const char* CTimeSupport::MONTHS[]   = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };

#ifdef WIN32
const uint64_t CTimeSupport::SECS_BETWEEN_EPOCHS = 11644473600;
const uint64_t CTimeSupport::SECS_TO_100NS       = 10000000;
#else
const uint64_t CTimeSupport::SECS_BETWEEN_EPOCHS = 11644473600llu;
const uint64_t CTimeSupport::SECS_TO_100NS       = 10000000llu;
#endif

// returns the current time (UTC)
uint64_t CTimeSupport::GetSystemTime(void) {

	uint64_t currentTime = 0;

#ifdef WIN32

	FILETIME ft;
	ULARGE_INTEGER ul;
	GetSystemTimeAsFileTime(&ft);
	ul.HighPart = ft.dwHighDateTime;
	ul.LowPart  = ft.dwLowDateTime;
	currentTime = ul.QuadPart;

#else

	timeval ft;
	gettimeofday(&ft, NULL);
	currentTime = (ft.tv_sec + SECS_BETWEEN_EPOCHS) * SECS_TO_100NS + ft.tv_usec * 10;

#endif

	return currentTime;
}

// converts the supplied string into a 64-bit time object (UTC)
uint64_t CTimeSupport::ConvertStringToTime(const string& timeString) {

	char buffer[26];
	struct tm time_tm;
	char* end_ptr = NULL;

	// reset the structure
	time_t convertTimeT = time(NULL);
	gmtime_s(&time_tm, &convertTimeT);

	char* pTime           = (char*)timeString.c_str();
	unsigned char timeLen = (unsigned char)timeString.size();

	// establish the month
	bool foundMonth = false;
	for(unsigned int i = 0; i < 12; i++) {
		if(strncmp(pTime + 4, MONTHS[i], 3) == 0) {
			time_tm.tm_mon = i;
			foundMonth     = true;
			break;
		}
	}

	if(!foundMonth) {
		cout << "ERROR: Unable to convert the month when parsing the timestamp." << endl;
		exit(1);
	}

	// establish the day of the month
	unsigned char startPos = 8;
	while(pTime[startPos] == ' ') startPos++;

	unsigned char endPos = startPos;
	while(pTime[endPos] != ' ') endPos++;
	unsigned char monthDayLen = endPos - startPos;

	bool hasLeadingZero = false;
	if((buffer[startPos] == '0') && (monthDayLen > 1)) hasLeadingZero = true;

	memcpy(buffer, pTime + startPos, monthDayLen);
	buffer[monthDayLen] = 0;

	time_tm.tm_mday = (int)strtol(buffer, &end_ptr, 10);
	if(end_ptr == buffer) {
		cout << "ERROR: Could not convert the day of the month string to an integer." << endl;
		exit(1);
	}

	// establish the hour
	startPos = endPos + 1;
	endPos = startPos;
	while(pTime[endPos] != ':') endPos++;
	unsigned int hourLen = endPos - startPos;

	memcpy(buffer, pTime + startPos, hourLen);
	buffer[hourLen] = 0;

	time_tm.tm_hour = (int)strtol(buffer, &end_ptr, 10);
	if(end_ptr == buffer) {
		cout << "ERROR: Could not convert the hour string to an integer." << endl;
		exit(1);
	}

	// establish the minutes
	memcpy(buffer, pTime + startPos + hourLen + 1, 2);
	buffer[2] = 0;

	time_tm.tm_min = (int)strtol(buffer, &end_ptr, 10);
	if(end_ptr == buffer) {
		cout << "ERROR: Could not convert the minute string to an integer." << endl;
		exit(1);
	}

	// establish the seconds
	memcpy(buffer, pTime + startPos + hourLen + 4, 2);
	buffer[2] = 0;

	time_tm.tm_sec = (int)strtol(buffer, &end_ptr, 10);
	if(end_ptr == buffer) {
		cout << "ERROR: Could not convert the second string to an integer." << endl;
		exit(1);
	}

	// establish the year
	memcpy(buffer, pTime + timeLen - 4, 4);
	buffer[4] = 0;

	time_tm.tm_year = (int)strtol(buffer, &end_ptr, 10) - 1900;
	if(end_ptr == buffer) {
		cout << "ERROR: Could not convert the year string to an integer." << endl;
		exit(1);
	}

	// UTC has no DST
	time_tm.tm_isdst = 0;

	// automatically guess the rest
	time_tm.tm_yday  = -1;
	time_tm.tm_wday  = -1;

	// convert the tm structure into epoch time
	return ConvertTimeT(mktime(&time_tm) - GetUtcOffset());
}

// displays the date associated with the supplied parameter
string CTimeSupport::ConvertTimeToString(const uint64_t& currentTime) {

	string ret;

#ifdef WIN32

	FILETIME ft;
	SYSTEMTIME st;
	ULARGE_INTEGER ul;

	ul.QuadPart       = currentTime;
	ft.dwHighDateTime = ul.HighPart;
	ft.dwLowDateTime  = ul.LowPart;

	FileTimeToSystemTime(&ft, &st);

	ostringstream sb;

	// "Fri Feb 9 12:15:48 2007"
	sb << setfill('0') << WEEKDAYS[st.wDayOfWeek] << " " << MONTHS[st.wMonth - 1]
	<< " " << st.wDay << " " << st.wHour << ":" << setw(2) << st.wMinute 
		<< ":" << setw(2) << st.wSecond << " " << st.wYear;

	ret = sb.str();

#else

	struct tm localTime;
	char timebuf[26];

	time_t currentTimeT = (time_t)((currentTime - 116444736000000000llu) / 10000000.0);
	localtime_s(&localTime, &currentTimeT);
	asctime_s(timebuf, 26, &localTime);

	ret = timebuf;

#endif

	return ret;
}

// converts a time_t variable to our 64-bit notation
uint64_t CTimeSupport::ConvertTimeT(const time_t& timeT) {
	return (timeT + SECS_BETWEEN_EPOCHS) * SECS_TO_100NS;
}

// get the offset between local time and UTC
time_t CTimeSupport::GetUtcOffset(void) {

	time_t currentTime = time(NULL);
	struct tm local_tm, utc_tm;
	time_t local_seconds, utc_seconds;

	localtime_s(&local_tm, &currentTime);
	gmtime_s(&utc_tm, &currentTime);

	local_seconds = mktime(&local_tm);
	utc_seconds   = mktime(&utc_tm);

	return utc_seconds - local_seconds;
}
