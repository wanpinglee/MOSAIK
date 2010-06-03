// ***************************************************************************
// LargeFileSupport.h - alleviates some cross-platform headaches when using
//                      large files. (64-bit addressing)
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#define SIZEOF_OFF_TYPE    8

#ifdef WIN32
#define ftell64(a)     _ftelli64(a)
#define fseek64(a,b,c) _fseeki64(a,b,c)
typedef __int64 off_type;
#elif SPARC
#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
typedef off64_t off_type;
#else

#ifdef __FreeBSD__
#define fstat64(a,b)   fstat(a,b)
#define stat64         stat
#endif

#define ftell64(a)     ftello(a)
#define fseek64(a,b,c) fseeko(a,b,c)
typedef off_t off_type;
#endif
