// ***************************************************************************
// CProgressCounter - displays a counter depicting the current progress of a 
//                    particular task. (operates in its own thread)
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <memory>
#include <cmath>
#include <ctime>
#include "SafeFunctions.h"

#ifdef WIN32
#include <windows.h>
#else
#include <sys/ioctl.h>
#endif

using namespace std;

typedef unsigned long long progress_counter_time_t;

// define our sleep interval (0.5 s)
#define PROGRESS_COUNTER_SLEEP_MS         500
#define PROGRESS_COUNTER_THREAD_SUCCESS   100

// define our multi-platform sleep function
#ifdef WIN32
#define ProgressCounterSleep(a) Sleep(a)
#else
#define ProgressCounterSleep(a) usleep(a * 1000)
#endif

template<class K>
class CProgressCounter {
public:
	CProgressCounter(K* pCurrentValue, bool* pRunning, const string units);
	~CProgressCounter();
	// starts the main thread
	static void StartThread(K* pCurrentValue, bool* pRunning, const string units);
	// wait for the main thread to finish
	static void WaitThread();
private:
	// displays the progress counter
	void Display(void);
	// returns the current time in milliseconds
	progress_counter_time_t GetSystemMilliseconds(void);
	// gets the current twirl character
	char GetTwirlCharacter(void);
	// returns the number with digit grouping
	static string AddDigitGrouping(const unsigned int value);
	// returns the number with digit grouping
	static string AddDigitGroupingFloat(const float value);
	// the current value
	K* mpCurrentValue;
	// the running state
	bool* mpRunning;
	// the rotating twirl orientation
	unsigned char mTwirlState;
	// the start time in milliseconds
	progress_counter_time_t mStartTime;
	// the units string
	string mUnits;
	// our thread parameters
	struct ThreadParams {
		K* pCurrentValue;
		bool* pRunning;
		string Units;
	};
#ifdef WIN32
	// our main thread
	static HANDLE mThread;
	// monitors the value parameter
	static DWORD WINAPI Monitor(PVOID arg);
#else
	// our main thread
	static pthread_t mThread;
	// monitors the value parameter
	static void* Monitor(void* arg);
#endif
};

template<class K>
#ifdef WIN32
HANDLE CProgressCounter<K>::mThread = 0;
#else
pthread_t CProgressCounter<K>::mThread = 0;
#endif

template<class K>
CProgressCounter<K>::CProgressCounter(K* pCurrentValue, bool* pRunning, const string units) 
: mpCurrentValue(pCurrentValue)
, mpRunning(pRunning)
, mTwirlState(0)
, mUnits(units)
{
	mStartTime = GetSystemMilliseconds();
}

template<class K>
CProgressCounter<K>::~CProgressCounter() {}

// starts the main thread
template<class K>
void CProgressCounter<K>::StartThread(K* pCurrentValue, bool* pRunning, const string units) {

	// initialize our thread parameters
	// N.B. it is the responsibility of the thread function to free the memory
	ThreadParams* pTP = new ThreadParams();
	pTP->pCurrentValue = pCurrentValue;
	pTP->pRunning      = pRunning;
	pTP->Units         = units;

#ifdef WIN32

	DWORD threadId;
	mThread = CreateThread(NULL, 0, CProgressCounter::Monitor, (PVOID)pTP, 0, &threadId);

	if(mThread == NULL) {
		cout << "ERROR: Unable to create thread for progress counter." << endl;
		exit(1);
	}

#else

	int res = pthread_create(&mThread, NULL, CProgressCounter::Monitor, (void*)pTP);

	if(res != 0) {
		cout << "ERROR: Unable to create thread for progress counter." << endl;
		exit(1);
	}

#endif
}

// wait for the main thread to finish
template<class K>
void CProgressCounter<K>::WaitThread() {
#ifdef WIN32

	// wait for the thread to finish
	if(WaitForSingleObject(mThread, INFINITE) != WAIT_OBJECT_0) {
		cout << "ERROR: Unable to join thread in progress counter." << endl;
		exit(1);
	}

	mThread = NULL;

#else

	void* threadResult;
	int res = pthread_join(mThread, &threadResult);

	if(res != 0) {
		cout << "ERROR: Unable to join thread in progress counter." << endl;
		exit(1);
	}

#endif
}

// Gets the current twirl character
template<class K>
char CProgressCounter<K>::GetTwirlCharacter(void) {

	char twirlCh;

	switch(mTwirlState) {
		case 0:
		case 4:
			twirlCh = '|';
			break;
		case 1:
		case 5:
			twirlCh = '/';
			break;
		case 2:
		case 6:
			twirlCh = '-';
			break;
		case 3:
		case 7:
			twirlCh = '\\';
			break;
		default:
			twirlCh = '*';
			break;
	}

	// increment the twirl state
	mTwirlState++;

	// reset the twirl state
	if(mTwirlState > 7) mTwirlState = 0;

	// return the twirl character
	return twirlCh;
}

// monitors the value parameter
template<class K>
#ifdef WIN32
DWORD WINAPI CProgressCounter<K>::Monitor(PVOID arg) {
#else
void* CProgressCounter<K>::Monitor(void* arg) {
#endif

	// initialize a new progress counter with the thread parameters
	ThreadParams* pTP = (ThreadParams*)arg;
	CProgressCounter<K> counter(pTP->pCurrentValue, pTP->pRunning, pTP->Units);

	// poll the current value
	while(*pTP->pRunning) {
		counter.Display();
		ProgressCounterSleep(PROGRESS_COUNTER_SLEEP_MS);
	}

	// display the final progress counter
	counter.Display();

	// cleanup
	if(pTP) delete pTP;

	// exit the thread
#ifdef WIN32
	return PROGRESS_COUNTER_THREAD_SUCCESS;
#else
	pthread_exit((void*)"SUCCESS");
#endif
}

// displays the progress counter
template<class K>
void CProgressCounter<K>::Display(void) {

	// localize the currentvalue
	K currentValue = *mpCurrentValue;
	bool finalPass = !*mpRunning;

	// display the progress counter
	if(finalPass) {

		// get the elapsed time
		progress_counter_time_t elapsedTime = GetSystemMilliseconds() - mStartTime;
		float rate = currentValue / (elapsedTime / 1000.0f);

		printf("\r%s: %s (%s %s/s)\n", mUnits.c_str(), AddDigitGrouping(currentValue).c_str(), AddDigitGroupingFloat(rate).c_str(), mUnits.c_str());

	} else printf("\r%s: %s %c", mUnits.c_str(), AddDigitGrouping(currentValue).c_str(), GetTwirlCharacter());

	fflush(stdout);
}

// returns the current time in milliseconds
template<class K>
progress_counter_time_t CProgressCounter<K>::GetSystemMilliseconds(void) {

	progress_counter_time_t currentTime = 0;

#ifdef WIN32

	FILETIME ft;
	ULARGE_INTEGER ul;
	GetSystemTimeAsFileTime(&ft);
	ul.HighPart = ft.dwHighDateTime;
	ul.LowPart  = ft.dwLowDateTime;
	currentTime = ul.QuadPart / 10000;

#else

	struct timeval tv;
	struct timezone tz;
	gettimeofday(&tv, &tz);
	currentTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;

#endif

	return currentTime;
}

// returns the number with digit grouping
template<class K>
string CProgressCounter<K>::AddDigitGrouping(const unsigned int value) {

	// convert the value to a string
	ostringstream sb;
	sb << value;
	string temp = sb.str();
	sb.str("");

	// no further processing required if the value is less than 1000
	if(value < 1000) return temp;

	// calculate the string length and ignore the sign
	unsigned char tempLen = (unsigned char)temp.size();

	// calculate how many leading digits before inserting the comma
	// N.B. this will never be a negative value
	unsigned char numLeadingDigits = tempLen % 3;

	bool printComma = false;
	if(numLeadingDigits > 0) printComma = true;

	// insert the leading digits
	unsigned char i = 0;
	for(; i < numLeadingDigits; i++) sb << temp[i];

	// insert the remaining digits
	for(; i < tempLen; i += 3) {
		if(printComma) sb << ",";
		for(unsigned char j = 0; j < 3; j++) sb << temp[i + j];
		printComma = true;
	}

	// return our string
	return sb.str();
}

// returns the number with digit grouping
template<class K>
string CProgressCounter<K>::AddDigitGroupingFloat(const float value) {

	// figure out which precision to use
	streamsize precision = 1;
	if(value < 1.0f)             precision = 4;
	else if(value < 10.0f)       precision = 2;
	else if(value >= 1000000.0f) precision = 0;

	// convert the value to a string
	ostringstream sb;
	sb << fixed << setprecision(precision) << value;
	string temp = sb.str();
	sb.str("");

	// no further processing required if the value is less than 1000
	if(value < 1000.0f) return temp;

	// calculate the string length and ignore the sign
	unsigned char decimalChars = 0;
	if(precision > 0) decimalChars = (unsigned char)precision + 1;

	unsigned char tempLen = (unsigned char)temp.size() - decimalChars;

	// calculate how many leading digits before inserting the comma
	// N.B. this will never be a negative value
	unsigned char numLeadingDigits = tempLen % 3;

	bool printComma = false;
	if(numLeadingDigits > 0) printComma = true;

	// insert the leading digits
	unsigned char i = 0;
	for(; i < numLeadingDigits; i++) sb << temp[i];

	// insert the remaining digits
	for(; i < tempLen; i += 3) {
		if(printComma) sb << ",";
		for(unsigned char j = 0; j < 3; j++) sb << temp[i + j];
		printComma = true;
	}

	// add the decimal portion
	if(precision > 0) sb << temp.substr(tempLen);

	// return our string
	return sb.str();
}
