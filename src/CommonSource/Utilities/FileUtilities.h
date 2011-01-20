// ***************************************************************************
// CFileUtilities - performs basic file operations.
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
#include <sys/time.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#endif

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include "Mosaik.h"
#include "LargeFileSupport.h"

#define COPY_BUFFER_SIZE     409600
#define TEMP_FILENAME_LENGTH 32

#ifdef WIN32
#define rm(a) _unlink(a)
#else
#define rm(a) unlink(a)
#endif

using namespace std;

class CFileUtilities {
public:
	// checks if a file exists, exits otherwise
	static bool CheckFile(const char* filename, bool showError);
	// checks if a directory exists, exits otherwise
	static void CheckDirectory(const string& directory);
	// checks if a directory exists, creates it otherwise
	static void CreateDir(const char* directory);
	// calculates the file size for the given filename
	static void GetFileSize(const string& filename, uint64_t& fileSize);
	// searches a directory for filenames that contain the specified string
	static void SearchDirectory(vector<string>& filenames, const char* directory);
	// moves the specified file to the specified directory
	static void MoveFile(const char* filename, const char* directory);
	// moves the specified file to the specified directory
	static void CopyFile(const char* filename, const char* directory);
	// returns the temp directory for the appropriate platform
	static void GetTempDirectory(string& tempDirectory);
	// generates a random filename in the temp directory
	static void GetTempFilename(string& tempFilename);
	// returns the file size for the specified filename
	static off_type GetFileSize(const char* filename);
	// returns true if a directory exists, false otherwise
	static bool DirExists(const char* directory);
};
