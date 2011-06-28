// ***************************************************************************
// Mosaik.h - handles global definitions throughout the MOSAIK suite.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#ifndef _MOSAIK_H_
#define _MOSAIK_H_

#ifndef WIN32
#include <stdint.h>
#include "SafeFunctions.h"
#endif

// ==============
// MOSAIK version
// ==============

#define MOSAIK_VERSION_DATE "2011-06-28"

// adopt a major.minor.build version number [1].[1].[3]
const unsigned char  MOSAIK_MAJOR_VERSION = 2;
const unsigned char  MOSAIK_MINOR_VERSION = 0;
const unsigned short MOSAIK_BUILD_VERSION = 122;

// ================================
// Platform specific variable sizes
// ================================

// Windows Vista 32-bit
// Fedora Core 7 32-bit
// Fedora Core 6 64-bit
// Itanium2      64-bit
#define SIZEOF_CHAR          1
#define SIZEOF_WCHAR         2
#define SIZEOF_SHORT         2
#define SIZEOF_INT           4
#define SIZEOF_FLOAT         4
#define SIZEOF_DOUBLE        8
#define SIZEOF_UINT64        8
#define MOSAIK_LITTLE_ENDIAN 1

#ifdef WIN32
typedef signed long long    int64_t;
typedef unsigned long long uint64_t;
#endif

#define NEGATIVE_ONE_INT     0xffffffff
#define NEGATIVE_TWO_INT     0xfffffffe
#define NEGATIVE_THREE_INT   0xfffffffd
#define NEGATIVE_FOUR_INT    0xfffffffc
#define MAX_SHORT            0xffff

// ==========================
// Platform specific file I/O 
// ==========================

#ifdef WIN32
const char OS_DIRECTORY_SEPARATOR = '\\';
#else
const char OS_DIRECTORY_SEPARATOR = '/';
#endif

#define DIRECTORY_NAME_LENGTH    255

// ====================================
// Enable unit test diagnostic messages
// ====================================

#ifdef UNITTEST
#define SILENTMODE if(0)
#else
#define SILENTMODE
#endif

// =================
// Aligner constants
// =================

#define NUM_REFERENCE_DIVIDER_BASES 500
const double HASH_REGION_EXTENSION_PERCENT     = 0.025;
const unsigned char REFERENCE_SEQUENCE_QUALITY = 40;

// ==================
// Archive signatures
// ==================

const unsigned char SIGNATURE_LENGTH = 6;
const char ALIGNER_SIGNATURE[ ]  = "MSKAA\4";
const char ALIGNER_SIGNATURE5[ ]  = "MSKAA\5";
const char SORT_SIGNATURE[ ]     = "MSKSA\1";

// ===========================
// Mapping Qualities constants
// ===========================
const double UU_COEFFICIENT = 0.486796848;
const double UU_INTERCEPT   = 35.45967112;
const double UM_COEFFICIENT = 0.426395518;
const double UM_INTERCEPT   = 19.29236958;
const double MM_COEFFICIENT = 0.327358673;
const double MM_INTERCEPT   = 4.350331532;


#endif
