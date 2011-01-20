// ***************************************************************************
// CAbstractDnaHash - superclass to the genome hash maps.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "AbstractDnaHash.h"

// defines the code for the empty dna hash code
const uint64_t CAbstractDnaHash::DNA_HASH_EMPTY_KEY = -1;

// specifies the largest resizeable hash table size
const unsigned int CAbstractDnaHash::LargestResizeableSize = 1 << 31;

// register our thread mutexes
pthread_mutex_t CAbstractDnaHash::mJumpCacheMutex;
pthread_mutex_t CAbstractDnaHash::mJumpKeyMutex;
pthread_mutex_t CAbstractDnaHash::mJumpPositionMutex;

CAbstractDnaHash::CAbstractDnaHash() 
: mHashes(NULL)
, mMemoryAllocated(false)
{}

CAbstractDnaHash::~CAbstractDnaHash() {}
