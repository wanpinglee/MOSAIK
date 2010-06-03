// ***************************************************************************
// UnorderedMap.h - handles all of the complications concerning hash maps on
//                  various platforms.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#ifdef __APPLE__

#include <ext/hash_map>

#ifndef CXX
namespace __gnu_cxx {

	template<>
	struct hash<std::string> {
		size_t operator()(const std::string& x) const {
			return hash<const char*>()( x.c_str() );
		}
	};

	template<>
	struct hash<uint64_t> {
		size_t operator()(const uint64_t r) const {
			return (size_t)r;
		}
	};
}
#define CXX
#endif

using namespace __gnu_cxx;

#define unordered_map hash_map

#else // all decent C++ compilers

#ifdef WIN32
#include <unordered_map>
#else // Linux
#include <tr1/unordered_map>
#endif

using namespace std::tr1;

#endif
