// ***************************************************************************
// CMemoryUtilities - used in monitoring memory usage and dynamic buffer
//                    reallocation.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <iomanip>

using namespace std;

#ifdef WIN32

#include <windows.h>
#include "psapi.h"

typedef DWORD process_id_t;
typedef SIZE_T mem_size_t;

#else

#include <sys/types.h>
#include <unistd.h>

typedef pid_t process_id_t;
typedef unsigned long long mem_size_t;
#define MEM_USAGE_BUFFER_LEN 4096
#define BYTES_PER_MEMORY_PAGE 4096

#endif

class CMemoryUtilities {
public:
	// returns the current process ID
	static process_id_t GetProcessID(void);
	// returns the number of bytes used by the current process
	static mem_size_t GetMemoryUsage(process_id_t processID);
	// checks if the buffer is large enough to accomodate the requested size
	static void CheckBufferSize(char* &pBuffer, unsigned int& bufferLen, const unsigned int requestedBytes);
	// checks if the buffer is large enough to accomodate the requested size
	static void CheckBufferSize(unsigned char* &pBuffer, unsigned int& bufferLen, const unsigned int requestedBytes);
};
