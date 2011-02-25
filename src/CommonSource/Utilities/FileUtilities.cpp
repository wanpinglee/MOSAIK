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

	bool fileExists = true;

#ifdef WIN32
	// open the file
        int fd = _open(filename, O_RDONLY);
#else
	int fd = open(filename, O_RDONLY);
#endif

	if (fd < 0)
        {
            if(showError) 
            {
                cout << "ERROR: Could not open " << filename << " for reading." << endl;
                exit(1);
            }
	} 
        else 
            fileExists = false;

	// close the file
        close(fd);

	return fileExists;
}

// checks if a file exists, exits otherwise
bool CFileUtilities::CheckTempFile(const char* filename, bool showError) {

	bool fileExists = true;

	// open the file
#ifdef WIN32
        int fd = _open(filename, O_RDWR|O_CREAT|O_EXCL, 0600);
#else
	int fd = open(filename, O_RDWR|O_CREAT|O_EXCL, 0600);
#endif
	if (fd < 0)
        {
            if(showError) 
            {
                cout << "ERROR: Could not open " << filename << " for reading." << endl;
                exit(1);
            }
	} 
        else 
            fileExists = false;

	// close the file
        close(fd);

	return fileExists;
}

// checks if a directory exists, exits otherwise
//void CFileUtilities::CheckDirectory(const string& directory) {
//
//	if( !DirExists( directory.c_str() ) ) {
//		cout << "ERROR: Could not find directory: " << directory << endl;
//		exit(1);
//	}
//}

// delete the directory
bool CFileUtilities::DeleteDir ( string directory ) {
	
	if( !DirExists( directory.c_str() ) ) return false;
	uint64_t directoryLen = directory.size();

#ifdef WIN32
	if ( directory[ directoryLen - 1 ] != '\\' )
		directory += '\\';
#else
	if ( directory[ directoryLen - 1 ] != OS_DIRECTORY_SEPARATOR )
		directory += OS_DIRECTORY_SEPARATOR;
#endif

	// first off, we need to create a pointer to a directory
	DIR *pdir = NULL; // remember, it's good practice to initialise a pointer to NULL!
	pdir = opendir ( directory.c_str() );
	if (pdir == NULL) // if pdir wasn't initialised correctly
		return false;

	string file;
	struct dirent *pent = NULL;
	while ( ( pent = readdir ( pdir ) ) != NULL ) { // while there is still something in the directory to list
		if ( pent == NULL ) // if pent has not been initialised correctly
			return false; // we couldn't do it
			
		file = directory + pent->d_name; // concatenate the strings to get the complete path
		remove( file.c_str() );
	}

	// finally, let's clean up
	closedir ( pdir ); // close the directory
	if (!rmdir( directory.c_str() ) ) return false; // delete the directory

	return true;
}

// checks if a directory exists, creates it otherwise
void CFileUtilities::CreateDir(const char* directory) {

	// return if the directory already exists
	if( DirExists( directory ) ) return;

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
	struct stat st;

	if ( stat( directory, &st ) == 0 )
		foundDirectory = true;
	
	
	//DIR *pDirectory = opendir(directory);
	//if ( pDirectory != NULL ) foundDirectory = true;

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
	// get hostname
	char hostname[256];
	uint32_t isHostname;
	isHostname = gethostname( hostname, sizeof(hostname) );


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

	if ( isHostname == 0 ) {
		tempDirectory += hostname;
		tempDirectory += '.';
	}

	char pidChar[ 16 ];
	memset( pidChar, 0, sizeof(char) * 16 );
	sprintf( pidChar, "%u", getpid() );
	tempDirectory += pidChar;
	tempDirectory += '\\';

	// checks if a directory exists, creates it otherwise
	CreateDir( tempDirectory.c_str() );


	// clean up
	delete [] tmpDir;
#else
	char* tmpDir = getenv( MOSAIK_TMP );
	if( tmpDir ) {
		tempDirectory = tmpDir;
		if( tempDirectory[ tempDirectory.size() - 1 ] != OS_DIRECTORY_SEPARATOR)
			tempDirectory += OS_DIRECTORY_SEPARATOR;
	} else {
		tempDirectory = "/tmp/";
	}

	if ( isHostname == 0 ) {
		tempDirectory += hostname;
		tempDirectory += '.';
	}

	char pidChar[ 16 ];
	memset( pidChar, 0, sizeof(char) * 16 );
	sprintf( pidChar, "%u", getpid() );
	tempDirectory += pidChar;
	tempDirectory += OS_DIRECTORY_SEPARATOR;

	// checks if a directory exists, creates it otherwise
	CreateDir( tempDirectory.c_str() );

#endif
}

// generates a random filename in the temp directory
void CFileUtilities::GetTempFilename(string& tempFilename) {

        uint64_t value = 0;
        //char hostname[256];

#ifdef WIN32
        value += _getpid();
#else
	timeval ft;
	gettimeofday(&ft, NULL);
        value += ((uint64_t) ft.tv_usec << 16) ^ ft.tv_sec ^ getpid();
#endif

	// seed the random number generator with supplied seed and time
	// define our set of random characters
        static const char letters[]= "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
	//static const char letters[]= "abcdefghijklmnopqrstuvwxyz0123456789";

	// get our temporary directory
	string tempDirectory;
	GetTempDirectory( tempDirectory );

        // get host name
        //int ghnVal = gethostname( hostname, sizeof(hostname) );
        //if (ghnVal == 0)
        //{
        //    tempFilename = tempDirectory + hostname + "XXXXXX";
        //}
        //else
        //{
        //    tempFilename = tempDirectory + "XXXXXX";
        //}
	tempFilename = tempDirectory + "XXXXXX";

	bool filenameExists = true;
        unsigned int count  = 0;
        unsigned int length = tempFilename.size();
	while(count++ < MAX_TMP_TRYING_TIME) 
        {
            uint64_t tempValue = value;

            // build our random filename
            for(unsigned int i = 0; i != 6; ++i)
            {
                tempFilename[length - 6 + i] = letters[tempValue % 62];
                tempValue /= 62;
            }

            // check if the file exists
            filenameExists = CheckTempFile(tempFilename.c_str(), false);

            if (filenameExists) {
#ifdef WIN32
		value += 7777;
#else
		gettimeofday(&ft, NULL);
		value += ((uint64_t) ft.tv_usec << 16) ^ ft.tv_sec ^ getpid();
#endif
            }else 
                return;
	}

        // exceed the maximum trying time
        // report an error
        cout << "Can not create a temporary file: " << tempFilename << endl;
        exit(1);
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
