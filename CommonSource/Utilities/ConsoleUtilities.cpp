// ***************************************************************************
// CConsole - adds color to command-line programs.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "ConsoleUtilities.h"

#ifdef WIN32
HANDLE CConsole::mStandardOutput;
CONSOLE_SCREEN_BUFFER_INFO CConsole::mCSBI;
#else
bool CConsole::mShowColors;
#endif

void CConsole::Bold(void) {
#ifdef WIN32
	SetConsoleTextAttribute(mStandardOutput, FOREGROUND_INTENSITY | FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
#else
	printf(ANSI_BOLD);
#endif
}

void CConsole::Heading(void) {
#ifdef WIN32
	SetConsoleTextAttribute(mStandardOutput, FOREGROUND_INTENSITY | FOREGROUND_RED | FOREGROUND_GREEN);
#else
	printf((mShowColors ? ANSI_YELLOW : ANSI_RED)); 
#endif
}

void CConsole::Initialize(void) {
#ifdef WIN32
	mStandardOutput = GetStdHandle(STD_OUTPUT_HANDLE);
	GetConsoleScreenBufferInfo(mStandardOutput, &mCSBI);
#else
	mShowColors = false;
	char* terminal = getenv("TERM");
	if(terminal) {
		if(strncmp(terminal, "xterm", 5) == 0) mShowColors = true;
	}
#endif
}

void CConsole::ProgressBar(void) {
#ifdef WIN32
	SetConsoleTextAttribute(mStandardOutput, FOREGROUND_INTENSITY | FOREGROUND_GREEN);
#else
	if(mShowColors) printf(ANSI_GREEN);
#endif
}

void CConsole::Red(void) {
#ifdef WIN32
	SetConsoleTextAttribute(mStandardOutput, FOREGROUND_INTENSITY | FOREGROUND_RED);
#else
	printf(ANSI_RED);
#endif
}

void CConsole::Reset(void) {
#ifdef WIN32
	SetConsoleTextAttribute(mStandardOutput, mCSBI.wAttributes);
#else
	printf(ANSI_RESET);
#endif
}
