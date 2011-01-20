// ***************************************************************************
// CConsole - adds color to command-line programs.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#ifdef WIN32
#include <windows.h>
#else
#include <cstdio>
#include <cstdlib>
#include <cstring>

// our ANSI escape codes
#define ANSI_BOLD   "\033[1m"
#define ANSI_YELLOW "\033[1;33m"
#define ANSI_GREEN  "\033[1;32m"
#define ANSI_RED    "\033[1;31m"
#define ANSI_RESET  "\033[0m"

#endif

class CConsole {
public:
	static void Bold(void);
	static void Heading(void);
	static void Initialize(void);
	static void ProgressBar(void);
	static void Red(void);
	static void Reset(void);
private:
#ifdef WIN32
	static HANDLE mStandardOutput;
	static CONSOLE_SCREEN_BUFFER_INFO mCSBI;
#else
	static bool mShowColors;
#endif
};
