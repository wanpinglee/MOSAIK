// ***************************************************************************
// CProgressBar - displays a bar depicting the current progress of a 
//                particular task. (operates in its own thread)
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
#include "ConsoleUtilities.h"
#include "PosixThreads.h"
#include "SafeFunctions.h"

#ifndef WIN32
#include <sys/ioctl.h>
#endif

using namespace std;

// define our sleep interval (0.5 s)
#define PROGRESS_BAR_SLEEP_MS         500
#define PROGRESS_BAR_THREAD_SUCCESS   100
#define PROGRESS_BAR_HISTORY_LENGTH   30	// 15 s history
#define PROGRESS_BAR_MS_TO_100NS      10000
#define PROGRESS_BAR_INITIAL_WAIT     3000

#define PROGRESS_BAR_PERCENT_LENGTH   4
#define PROGRESS_BAR_ENDPOINTS_LENGTH 2
#define PROGRESS_BAR_RATE_LENGTH      14
#define PROGRESS_BAR_ETA_LENGTH       15
#define PROGRESS_BAR_TWIRL_LENGTH     2

typedef unsigned long long progress_bar_time_t;

// define our multi-platform sleep function
#ifdef WIN32
#define ProgressBarSleep(a) Sleep(a)
#else
#define ProgressBarSleep(a) usleep(a * 1000)
#endif

template<class K>
class CProgressBar {
public:
	// constructor
	CProgressBar(K* pCurrentValue, const K minValue, const K maxValue, const string units);
	// destructor
	~CProgressBar();
	// starts the main thread
	static void StartThread(K* pCurrentValue, const K minValue, const K maxValue, const string units);
	// wait for the main thread to finish
	static void WaitThread();
private:
	// returns the completion rate
	float CalculateRate(const progress_bar_time_t& elapsedTime);
	// displays the progress bar
	void Display(const bool isFinalPass);
	// gets the current twirl character
	char GetTwirlCharacter(void);
	// returns the screen size
	static unsigned short GetScreenWidth(void);
	// returns the number with digit grouping
	static string AddDigitGrouping(const float value);
	// returns the current time in milliseconds
	static progress_bar_time_t GetSystemMilliseconds(void);
	// saves the average rate in the history
	void SaveAverageRate(float averageRate);
	// resets the rate history
	void ResetHistory(void);
	// the current value
	K* mpCurrentValue;
	// the minimum value
	K mMinValue;
	// the maximum value
	K mMaxValue;
	// the rotating twirl orientation
	unsigned char mTwirlState;
	// the screen width
	unsigned short mScreenWidth;
	// the progress bar length
	unsigned short mBarLength;
	// the units string
	string mUnits;
	// specifies the buffer used when writing the progress bar
	char* mBuffer;
	// specifies the current buffer size
	unsigned int mBufferLen;
	// specifies the current rate calculation wait time
	progress_bar_time_t mCalculationWaitTime;
	// our history
	struct HistorySettings {
		progress_bar_time_t StartTime;
		progress_bar_time_t LastTime;
		progress_bar_time_t LastElapsedTime;
		K                   LastValue;
		float               LastWeightedRate;
		float               Rates[PROGRESS_BAR_HISTORY_LENGTH];
		unsigned char       Position;
		unsigned char       NumRatesUsed;
		unsigned int        RemainingSeconds;
	} mHistory;
	// our thread parameters
	struct ThreadParams {
		K* pCurrentValue;
		K MinValue;
		K MaxValue;
		string Units;
	};
	// our main thread
	static pthread_t mThread;
	// monitors the value parameter
	static void* Monitor(void* arg);
};

template<class K>
pthread_t CProgressBar<K>::mThread;

template<class K>
CProgressBar<K>::CProgressBar(K* pCurrentValue, const K minValue, const K maxValue, const string units) 
: mpCurrentValue(pCurrentValue)
, mMinValue(minValue)
, mMaxValue(maxValue)
, mTwirlState(0)
, mUnits(units)
, mBuffer(NULL)
{
	// get the current screen width
	mScreenWidth = GetScreenWidth();

	// calculate the progress bar length according to:
	mBarLength = (mScreenWidth - 1) - PROGRESS_BAR_PERCENT_LENGTH - PROGRESS_BAR_ENDPOINTS_LENGTH - 
		PROGRESS_BAR_RATE_LENGTH - PROGRESS_BAR_ETA_LENGTH - (unsigned short)units.size() -
		PROGRESS_BAR_TWIRL_LENGTH;

	// allocate memory for our buffer
	mBufferLen = mScreenWidth + 1000;

	// get the start time
	mHistory.StartTime        = GetSystemMilliseconds();
	mHistory.LastTime         = 0;
	mHistory.LastElapsedTime  = 0;
	mHistory.LastValue        = 0;
	mHistory.LastWeightedRate = 0.0f;
	mHistory.RemainingSeconds = 0;

	// set the default calculation wait time
	mCalculationWaitTime = 1000;

	// reset the history
	ResetHistory();

	try {
		mBuffer = new char[mBufferLen];
	} catch(bad_alloc) {
		cout << "ERROR: Unable to allocate the buffer used by the progress bar." << endl;
		exit(1);
	}
}

template<class K>
CProgressBar<K>::~CProgressBar() {
	if(mBuffer) delete [] mBuffer;
}

// starts the main thread
template<class K>
void CProgressBar<K>::StartThread(K* pCurrentValue, const K minValue, const K maxValue, const string units) {

	// initialize our thread parameters
	// N.B. it is the responsibility of the thread function to free the memory
	ThreadParams* pTP = new ThreadParams();
	pTP->pCurrentValue = pCurrentValue;
	pTP->MinValue      = minValue;
	pTP->MaxValue      = maxValue;
	pTP->Units         = units;

	int res = pthread_create(&mThread, NULL, CProgressBar::Monitor, (void*)pTP);

	if(res != 0) {
		cout << "ERROR: Unable to create thread for progress bar." << endl;
		exit(1);
	}
}

// wait for the main thread to finish
template<class K>
void CProgressBar<K>::WaitThread() {

	void* threadResult;
	int res = pthread_join(mThread, &threadResult);

	if(res != 0) {
		cout << "ERROR: Unable to join thread in progress bar." << endl;
		exit(1);
	}
}

// Gets the current twirl character
template<class K>
char CProgressBar<K>::GetTwirlCharacter(void) {

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
	++mTwirlState;

	// reset the twirl state
	if(mTwirlState > 7) mTwirlState = 0;

	// return the twirl character
	return twirlCh;
}

// monitors the value parameter
template<class K>
void* CProgressBar<K>::Monitor(void* arg) {

	// initialize a new progress bar with the thread parameters
	ThreadParams* pTP = (ThreadParams*)arg;

	// increase the scope to force the deconstructors
	{
		CProgressBar<K> bar(pTP->pCurrentValue, pTP->MinValue, pTP->MaxValue, pTP->Units);

		// poll the current value
		while(*pTP->pCurrentValue < pTP->MaxValue) {
			bar.Display(false);
			ProgressBarSleep(PROGRESS_BAR_SLEEP_MS);
		}

		// display the final progress bar
		bar.Display(true);
	}

	// cleanup
	if(pTP) delete pTP;

	// exit the thread
	pthread_exit((void*)0);
	return 0;
}

// displays the progress bar
template<class K>
void CProgressBar<K>::Display(const bool isFinalPass) {

	// localize the currentvalue
	K currentValue = *mpCurrentValue;

	// add a logic helper
	unsigned char percentageCompleted = (unsigned char)(currentValue / (float)mMaxValue * 100.0f);

	// get the current time and calculate the elapsed time
	progress_bar_time_t currentTime = GetSystemMilliseconds();
	progress_bar_time_t elapsedTime = currentTime - mHistory.StartTime;
	if(isFinalPass) {
		elapsedTime = mHistory.LastElapsedTime;
		if(elapsedTime < 1000) elapsedTime = 1000;
	}

	if(currentValue > mMaxValue) {
		cout << "ERROR: The current value is greater than maximum value in the progress bar." << endl;
		cout << "       current value: " << currentValue << ", maximum value: " << mMaxValue << endl;
		exit(1);
	}

	if(currentValue < mMinValue) {
		cout << "ERROR: The current value is less than the minimum value in the progress bar." << endl;
		cout << "       current value: " << currentValue << ", minimum value: " << mMinValue << endl;
		exit(1);
	}

	// initialize our buffer variables
	char* pBuffer = mBuffer;
	unsigned short bufferRemaining = mBufferLen;

	// calculate and print the percentage completed
	if(percentageCompleted < 100) sprintf_s(pBuffer, bufferRemaining, "%2d%% [", percentageCompleted);
	else sprintf_s(pBuffer, bufferRemaining, "100%%[");
	pBuffer         += 5;
	bufferRemaining -= 5;

	*pBuffer = 0; ++pBuffer;

	// create the progress bar
	unsigned short displayedArrowLength = (unsigned short)(mBarLength * (currentValue / (float)mMaxValue));
	unsigned short whiteSpaceLength     = mBarLength - displayedArrowLength;

	memset(pBuffer, '=', displayedArrowLength);
	if((displayedArrowLength < mBarLength) && (displayedArrowLength > 0)) pBuffer[displayedArrowLength - 1] = '>';
	pBuffer         += displayedArrowLength;
	bufferRemaining -= displayedArrowLength;

	if(whiteSpaceLength > 0) {
		memset(pBuffer, ' ', whiteSpaceLength);
		pBuffer         += whiteSpaceLength;
		bufferRemaining -= whiteSpaceLength;
	}

	*pBuffer = 0; ++pBuffer;

	// calculate the rate
	float rate = mHistory.LastWeightedRate;
	progress_bar_time_t intervalMilliseconds = elapsedTime - mHistory.LastTime;

	if(!isFinalPass) {
		if(intervalMilliseconds > mCalculationWaitTime)
			rate = CalculateRate(elapsedTime);
	} else rate = (mMaxValue - mMinValue) / (elapsedTime / 1000.0f);
	//cout << "rate: " << rate << ", elapsed time: " << elapsedTime << endl;

	unsigned short rateLength = PROGRESS_BAR_RATE_LENGTH + (unsigned short)mUnits.size() + 1;
	if(isFinalPass || (elapsedTime > PROGRESS_BAR_INITIAL_WAIT)) {

		sprintf_s(pBuffer, bufferRemaining, "] %9s %s/s ", AddDigitGrouping(rate).c_str(), mUnits.c_str());

	} else {

		memset(pBuffer, ' ', rateLength);
		pBuffer[rateLength] = 0;
		pBuffer[0] = ']';
	} 

	pBuffer         += rateLength;
	bufferRemaining -= rateLength;

	if(!isFinalPass) {

		// calculate the ETA after some initial time has passed
		if(elapsedTime > PROGRESS_BAR_INITIAL_WAIT) {

			// calculate the number of remaining seconds
			K remainingValues = mMaxValue - currentValue;
			unsigned int remainingSeconds = (unsigned int)(remainingValues / mHistory.LastWeightedRate);

			// calculate hours
			unsigned int hours = remainingSeconds / 3600; 
			remainingSeconds %= 3600;

			// calculate minutes
			unsigned int minutes = remainingSeconds / 60; 
			remainingSeconds %= 60;

			// calculate seconds
			unsigned int seconds = remainingSeconds;

			//ETA 365.3 days 
			// ETA 123:33:12
			//     ETA 33:12
			//      ETA 12 s 

			if(hours > 999) {        // show days

				float days = remainingValues / rate / 86400.0f;
				sprintf_s(pBuffer, bufferRemaining, "ETA %3.1f days ", days);

			} else if(hours > 0) {   // show hours

				sprintf_s(pBuffer, bufferRemaining, " ETA %3d:%02d:%02d ", hours, minutes, seconds);

			} else if(minutes > 0) { // show minutes

				sprintf_s(pBuffer, bufferRemaining, "     ETA %02d:%02d ", minutes, seconds);

			} else {                 // show seconds

				sprintf_s(pBuffer, bufferRemaining, "      ETA %2d s ", seconds);
			}

		} else {

			// fill the ETA field with blanks during the initial wait
			memset(pBuffer, ' ', PROGRESS_BAR_ETA_LENGTH);
			pBuffer[PROGRESS_BAR_ETA_LENGTH] = 0;
		}

	} else { 

		unsigned int elapsedSeconds = (unsigned int)(elapsedTime / 1000.0f);

		// calculate hours
		unsigned int hours = elapsedSeconds / 3600; 
		elapsedSeconds %= 3600;

		// calculate minutes
		unsigned int minutes = elapsedSeconds / 60; 
		elapsedSeconds %= 60;

		// calculate seconds
		unsigned int seconds = elapsedSeconds;

		if(hours > 999) {        // show days

			float days = elapsedTime / 86400.0f;
			sprintf_s(pBuffer, bufferRemaining, " in %3.1f days ", days);

		} else if(hours > 0) {   // show hours

			sprintf_s(pBuffer, bufferRemaining, "  in %3d:%02d:%02d ", hours, minutes, seconds);

		} else if(minutes > 0) { // show minutes

			sprintf_s(pBuffer, bufferRemaining, "      in %02d:%02d ", minutes, seconds);

		} else {                 // show seconds

			sprintf_s(pBuffer, bufferRemaining, "       in %2d s ", seconds);
		}
	}

	pBuffer += PROGRESS_BAR_ETA_LENGTH;

	// add the twirl
	if(!isFinalPass) { 
		*pBuffer = GetTwirlCharacter();
	} else {
		*pBuffer = ' ';
		++pBuffer;
		*pBuffer = '\n';
	}
	++pBuffer;

	// add the null terminator
	*pBuffer = 0;

	// display the progress bar
	printf("\r%s", mBuffer);
	CConsole::ProgressBar(); printf("%s", mBuffer + 6); CConsole::Reset();
	printf("%s", mBuffer + mBarLength + 7);
	fflush(stdout);

	mHistory.LastElapsedTime = elapsedTime;
}

// returns the screen size
template<class K>
unsigned short CProgressBar<K>::GetScreenWidth(void) {

	unsigned short screenWidth = 80;

#ifdef WIN32

	CONSOLE_SCREEN_BUFFER_INFO csbiInfo;
	if(GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbiInfo)) 
		screenWidth = csbiInfo.dwSize.X;

#else

	struct winsize ws;
	if(ioctl(fileno(stdout), TIOCGWINSZ, &ws) == 0) 
		screenWidth = ws.ws_col;

#endif

	return screenWidth;
}

// returns the number with digit grouping
template<class K>
string CProgressBar<K>::AddDigitGrouping(const float value) {

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
	for(; i < numLeadingDigits; ++i) sb << temp[i];

	// insert the remaining digits
	for(; i < tempLen; i += 3) {
		if(printComma) sb << ",";
		for(unsigned char j = 0; j < 3; ++j) sb << temp[i + j];
		printComma = true;
	}

	// add the decimal portion
	if(precision > 0) sb << temp.substr(tempLen);

	// return our string
	return sb.str();
}

// returns the completion rate
template<class K>
float CProgressBar<K>::CalculateRate(const progress_bar_time_t& elapsedTime) {

	// return if no time has elapsed yet
	if(elapsedTime == 0) return 0.0f;

	// calculate the interval time and values
	K currentValue = *mpCurrentValue;
	float intervalSeconds = (elapsedTime - mHistory.LastTime) / 1000.0f;
	K intervalValue       = currentValue - mHistory.LastValue;

	// return if we're seeing less than one value per interval
	SaveAverageRate(intervalValue / intervalSeconds);

	// calculate the weights for recent rates
	float tempSum = 0.0f;

	unsigned char pos = mHistory.Position + 1;
	if(pos >= PROGRESS_BAR_HISTORY_LENGTH) pos = 0;

	unsigned int currentWeight = 1;
	for(unsigned int i = 0; i < PROGRESS_BAR_HISTORY_LENGTH; ++i, ++currentWeight) {
		tempSum += currentWeight * mHistory.Rates[pos++];
		if(pos >= PROGRESS_BAR_HISTORY_LENGTH) pos = 0;
	}

	// calculated the added sum based on our history length
	unsigned int weightedSum   = 0;
	currentWeight = PROGRESS_BAR_HISTORY_LENGTH;
	for(unsigned int i = 0; i < mHistory.NumRatesUsed; ++i, --currentWeight) {
		weightedSum += currentWeight;
	}

	float weightedRecentAverage = tempSum / (float)weightedSum;

	// update the last time, value, and weighted average
	mHistory.LastTime         = elapsedTime;
	mHistory.LastValue        = currentValue;
	mHistory.LastWeightedRate = weightedRecentAverage;

	// update the calculation waiting time
	mCalculationWaitTime = (progress_bar_time_t)(weightedSum / (float)tempSum * 10000.0f);
	if(mCalculationWaitTime < 1000) mCalculationWaitTime = 1000;

	return weightedRecentAverage;
}

// returns the current time in milliseconds
template<class K>
progress_bar_time_t CProgressBar<K>::GetSystemMilliseconds(void) {

	progress_bar_time_t currentTime = 0;

#ifdef WIN32

	FILETIME ft;
	ULARGE_INTEGER ul;
	GetSystemTimeAsFileTime(&ft);
	ul.HighPart = ft.dwHighDateTime;
	ul.LowPart  = ft.dwLowDateTime;
	currentTime = ul.QuadPart / PROGRESS_BAR_MS_TO_100NS;

#else

	struct timeval tv;
	struct timezone tz;
	gettimeofday(&tv, &tz);
	currentTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;

#endif

	return currentTime;
}

// saves the average rate in the history
template<class K>
void CProgressBar<K>::SaveAverageRate(float averageRate) {

	++mHistory.Position;
	++mHistory.NumRatesUsed;

	if(mHistory.Position >= PROGRESS_BAR_HISTORY_LENGTH) mHistory.Position = 0;

	if(mHistory.NumRatesUsed > PROGRESS_BAR_HISTORY_LENGTH) 
		mHistory.NumRatesUsed = PROGRESS_BAR_HISTORY_LENGTH;

	mHistory.Rates[mHistory.Position] = averageRate;
}

// resets the rate history
template<class K>
void CProgressBar<K>::ResetHistory(void) {
	mHistory.NumRatesUsed     = 0;
	mHistory.Position         = 0;
	uninitialized_fill(mHistory.Rates, mHistory.Rates + PROGRESS_BAR_HISTORY_LENGTH, 0.0f);
}
