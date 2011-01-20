// ***************************************************************************
// CFileUtilities - performs basic file operations.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "FileUtilities.h"

#ifdef WIN32
#pragma warning (disable:6309)
#pragma warning (disable:6387)
#endif

// checks if a file exists, exits otherwise
bool CFileUtilities::CheckFile(const char* filename, bool showError) {

	bool fileExists = false;

	// open the file
	ifstream in(filename);

	if(in.fail()) {
		if(showError) {
			cout << "ERROR: Could not open " << filename << " for reading." << endl;
			exit(1);
		}
	} else fileExists = true;

	// close the file
	in.close();

	return fileExists;
}

// checks if a directory exists, exits otherwise
void CFileUtilities::CheckDirectory(const string& directory) {

	if(!DirExists(directory.c_str())) {
		cout << "ERROR: Could not find directory: " << directory << endl;
		exit(1);
	}
}

// checks if a directory exists, creates it otherwise
void CFileUtilities::CreateDir(const char* directory) {

	// return if the directory already exists
	if(DirExists(directory)) return;

#ifdef WIN32

	// convert the directory from multi-byte to wide characters
	int numWideChar = MultiByteToWideChar(CP_ACP, 0, directory, (int)strlen(directory) + 1, NULL, 0);
	LPWSTR wDirectory = new WCHAR[numWideChar];
	MultiByteToWideChar(CP_ACP, 0, directory, (int)strlen(directory) + 1, wDirectory, numWideChar);

	// N.B. Make sure that you have a trailing slash before calling this function
	if(!CreateDirectory(wDirectory, NULL)) {
		cout << "ERROR: Unable to create directory: " << directory << endl;
		exit(1);
	}

	// cleanup
	delete [] wDirectory;

#else

	if(mkdir(directory, S_IRWXU | S_IRGRP | S_IXGRP) != 0) {
		cout << "ERROR: Unable to create directory: " << directory << endl;
		exit(1);
	}

#endif
}

// returns true if a directory exists, false otherwise
bool CFileUtilities::DirExists(const char* directory) {

	bool foundDirectory = false;

#ifdef WIN32

	// remove 
	int dirLen = (int)strlen(directory);

	// convert the directory from multi-byte to wide characters
	int numWideChar = MultiByteToWideChar(CP_ACP, 0, directory, dirLen + 1, NULL, 0);
	LPWSTR wDirectory = new WCHAR[numWideChar];
	MultiByteToWideChar(CP_ACP, 0, directory, dirLen + 1, wDirectory, numWideChar);

	// check the directory
	DWORD dwAttr = GetFileAttributes(wDirectory);
	if(dwAttr == FILE_ATTRIBUTE_DIRECTORY) foundDirectory = true;

	// perform some clean up
	delete [] wDirectory;

#else

	DIR *pDirectory = opendir(directory);
	if(pDirectory != NULL) foundDirectory = true;

#endif

	return foundDirectory;
}

// calculates the file size for the given filename
void CFileUtilities::GetFileSize(const string& filename, uint64_t& fileSize) {

#ifdef WIN32

	HANDLE hFile = CreateFileA(filename.c_str(), GENERIC_READ, 0, NULL, OPEN_EXISTING, 0, NULL);
	LARGE_INTEGER tempSize;
	GetFileSizeEx(hFile, &tempSize);
	CloseHandle(hFile);

	// convert the file size to something usable
	memcpy((char*)&fileSize, (char*)&tempSize.QuadPart, SIZEOF_UINT64);

#else 

	FILE* hFile = NULL;
	fopen_s(&hFile, filename.c_str(), "rb");

	struct stat64 statBuffer;
	if(fstat64(fileno(hFile), &statBuffer) == -1) {
		cout << "ERROR: Unable to get the file size for " << filename << endl;
		exit(1);
	}

	fclose(hFile);
	fileSize = statBuffer.st_size;

#endif
}

// searches a directory for filenames that contain the specified string
void CFileUtilities::SearchDirectory(vector<string>& filenames, const char* directory) {

	// Check if we have a trailing slash
	unsigned int directoryNameLen = (unsigned int)strlen(directory);
	char tDirectory[DIRECTORY_NAME_LENGTH];

	if(directory[directoryNameLen - 1] != OS_DIRECTORY_SEPARATOR) {
		sprintf_s(tDirectory, DIRECTORY_NAME_LENGTH * SIZEOF_CHAR, "%s%c", directory, OS_DIRECTORY_SEPARATOR);
	} else sprintf_s(tDirectory, DIRECTORY_NAME_LENGTH * SIZEOF_CHAR, "%s", directory);

#ifdef WIN32

	WIN32_FIND_DATA FindFileData;
	char filename[DIRECTORY_NAME_LENGTH];
	wchar_t findPattern[DIRECTORY_NAME_LENGTH];

	// Get the first file
	swprintf_s(findPattern, DIRECTORY_NAME_LENGTH, L"%S*.*", tDirectory); 
	HANDLE hList = FindFirstFile(findPattern, &FindFileData);

	bool isFinished = false;
	while(!isFinished) {

		// create our filename
		sprintf_s(filename, DIRECTORY_NAME_LENGTH, "%s%S", tDirectory, FindFileData.cFileName);

		// Add the file to the vector if it is a normal file
		if(FindFileData.dwFileAttributes != FILE_ATTRIBUTE_DIRECTORY) filenames.push_back(filename);

		// find the next file
		if(!FindNextFile(hList, &FindFileData)) 
			if(GetLastError() == ERROR_NO_MORE_FILES) isFinished = true;
	}

	FindClose(hList);

#else

	DIR* dirp = opendir(tDirectory);
	char filename[255];

	while(dirp) {
		dirent* dp;

		if((dp = readdir(dirp)) != NULL) {
			if((strcmp(dp->d_name, ".") != 0) && (strcmp(dp->d_name, "..") != 0)) {
				sprintf(filename, "%s%s", tDirectory, dp->d_name);
				filenames.push_back(filename);
			}
		} else break;
	}
#endif
}

// moves the specified file to the specified directory
void CFileUtilities::MoveFile(const char* filename, const char* directory) {
	CopyFile(filename, directory);
	if(rm(filename) != 0) {
		cout << "ERROR: Unable to erase file (" << filename << ") during move operation." << endl;
		exit(1);
	}
}

// copies the specified file to the specified directory
// TODO: replace with platform specific commands
void CFileUtilities::CopyFile(const char* filename, const char* directory) {

	// extract the base filename
	string baseFilename = filename;
	string::size_type delimiterPosition = baseFilename.rfind(OS_DIRECTORY_SEPARATOR);

	if(delimiterPosition != string::npos)
		baseFilename = baseFilename.substr(delimiterPosition + 1);

	// re-construct our output filename
	string outFilename = directory + baseFilename;

	// open our files
	ifstream in(filename, ios::binary);
	ofstream out(outFilename.c_str(), ios::binary);

	if(!in) {
		cout << "ERROR: Unable to open (" << filename << ") for copying (input)." << endl;
		exit(1);
	}

	if(!out) {
		cout << "ERROR: Unable to open (" << outFilename << ") for copying (output)." << endl;
		exit(1);
	}

	// define our copy buffer
	char* buffer = new char[COPY_BUFFER_SIZE];

	// copy the file
	while(in) {
		in.read(buffer, COPY_BUFFER_SIZE);
		std::streamsize numBytesRead = in.gcount();
		out.write(buffer, numBytesRead);
	}

	// close our files
	out.close();
	in.close();

	// cleanup
	delete [] buffer;
}

// returns the temp directory for the appropriate platform
void CFileUtilities::GetTempDirectory(string& tempDirectory) {

	const char* MOSAIK_TMP = "MOSAIK_TMP";

#ifdef WIN32

	// allocate the memory needed to store the larger variable
	size_t requiredTmpSize       = 0;
	size_t requiredMosaikTmpSize = 0;

	getenv_s(&requiredTmpSize,       NULL, 0, "TMP");
	getenv_s(&requiredMosaikTmpSize, NULL, 0, MOSAIK_TMP);

	size_t requiredSize = max(requiredTmpSize, requiredMosaikTmpSize);

	char* tmpDir = new char[requiredSize];

	// get the environment variables
	bool foundProblem = false;
	if(requiredMosaikTmpSize != 0) {
		getenv_s(&requiredMosaikTmpSize, tmpDir, requiredMosaikTmpSize, MOSAIK_TMP);
	} 

	getenv_s(&requiredTmpSize,       tmpDir, requiredTmpSize,       "TMP");

	tempDirectory = tmpDir;
	tempDirectory += '\\';

	// clean up
	delete [] tmpDir;
#else
	char* tmpDir = getenv(MOSAIK_TMP);
	if(tmpDir) {
		tempDirectory = tmpDir;
		if(tempDirectory[tempDirectory.size() - 1] != OS_DIRECTORY_SEPARATOR)
			tempDirectory += OS_DIRECTORY_SEPARATOR;
	} else tempDirectory = "/tmp/";
#endif
}

// generates a random filename in the temp directory
void CFileUtilities::GetTempFilename(string& tempFilename) {

	unsigned int seed;

#ifdef WIN32
	FILETIME ft;
	GetSystemTimeAsFileTime(&ft);
	seed = ft.dwLowDateTime;
#else
	timeval ft;
	gettimeofday(&ft, NULL);
	seed = (unsigned int)ft.tv_usec;
#endif

	// seed the random number generator with supplied seed and time
	srand(seed);

	// define our set of random characters
	string randomChars = "abcdefghijklmnopqrstuvwxyz0123456789";
	unsigned int numRandomChars = (unsigned int)randomChars.size();

	// get our temporary directory
	string tempDirectory;
	GetTempDirectory(tempDirectory);

	// define a stringbuilder
	ostringstream sb;

	bool filenameExists = true;
	while(filenameExists) {

		sb << tempDirectory;

		// build our random filename
		for(unsigned int i = 0; i < TEMP_FILENAME_LENGTH; i++)
			sb << randomChars[rand() % numRandomChars];

		// add our file extension
		sb << ".tmp";
		tempFilename = sb.str();
		sb.str("");

		// check if the file exists
		filenameExists = CheckFile(tempFilename.c_str(), false);
	}
}

// returns the file size for the specified filename
off_type CFileUtilities::GetFileSize(const char* filename) {

	FILE* FILEHANDLE  = NULL;
	off_type fileSize = 0;


	if(fopen_s(&FILEHANDLE, filename, "rb") != 0) {
		cout << "ERROR: Unable to open file (" << filename << ") when getting file size." << endl;
		exit(1);
	}

	if(FILEHANDLE) {

		if(fseek64(FILEHANDLE, 0, SEEK_END) != 0) {
			cout << "ERROR: Unable to go to the end of the file (" << filename << ")" << endl;
			exit(1);
		}

		fileSize = ftell64(FILEHANDLE);
		fclose(FILEHANDLE);
	}

	return fileSize;
}
