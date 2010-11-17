// ***************************************************************************
// CMosaikString - a fast and lightweight string class with some improvements
//                 for handling sequence data.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <ostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "SafeFunctions.h"

#define GROWTH_FACTOR 10
#define PACK_MASK     15

typedef unsigned short uint16_t;
#define swapByte(a, b) { char c = (a); (a) = (b); (b) = c; }
#define get16bits(d) (*((const uint16_t*)(d)))

class CMosaikString {
public:
	// constructor
	CMosaikString(void);
	// copy constructor
	CMosaikString(const CMosaikString& ms);
	// copy constructor
	CMosaikString(const char* c_str);
	// copy constructor
	CMosaikString(const char* c_str, const unsigned int len);
	// destructor
	~CMosaikString(void);
	// assignment operator
	CMosaikString& operator=(const CMosaikString& ms);
	// assignment operator
	CMosaikString& operator=(const char* c_str);
	// larger than operator
	bool operator>(const CMosaikString& c) const;
	// less than operator
	bool operator<(const CMosaikString& c) const;
	// not equal operator
	bool operator!=(const CMosaikString& c) const;
	// equal operator
	bool operator==(const CMosaikString& c) const;
	// element operator
	char& operator[](const int index);
	// appends the specified string to the current string
	void Append(const char* s);
	// appends the specified string to the current string
	void Append(const char* s, const unsigned int sLen);
	// returns a const pointer to the data
	const char* CData(void) const;
	// copies the specified c-style string
	void Copy(const char* string, const unsigned int numBytes);
	// returns a pointer to the data
	char* Data(void);
	// decrements each character in the string by the specified amount
	void Decrement(const char amount);
	// fills the string with numBytes copies of the ch
	void Fill(const char ch, const unsigned int numBytes);
	// returns the hash value
	size_t GetHash(void) const;
	// increments each character in the string by the specified amount
	void Increment(const char amount);
	// joins two strings (used by the Smith-Waterman caching algorithm)
	void Join(const char* s1, const unsigned int s1Length, const char* s2, const unsigned int s2Length);
	// returns the size of the data
	unsigned int Length(void) const;
	// packs both the original bases and the supplied bases into 4-bit notation
	void Pack(const CMosaikString& ms);
	// prepends the specified string to the current string
	void Prepend(const CMosaikString& ms);
	// prepends the specified string to the current string
	void Prepend(const char* s, const unsigned int sLen);
	// replaces all occurrences of the first parameter with the second parameter
	void Replace(const char oldCh, const char newCh);
	// removes all occurrences of the specified character
	void Remove(const char ch);
	// reserve the specified number of bytes (destructive)
	void Reserve(const unsigned int numBytes);
	// reverses the contents of the string
	void Reverse(void);
	// sets the length to the specified size
	void SetLength(const unsigned int length);
	// retrieves the specified substring
	std::string Substring(const unsigned int position, const unsigned int length) const;
	// trims the first specified number of bytes
	void TrimBegin(unsigned int numBytes);
	// trims the last specified number of bytes
	void TrimEnd(unsigned int numBytes);
	// unpacks both the original bases and the supplied bases from 4-bit notation
	void Unpack(CMosaikString& ms);
	// converts the string to uppercase
	void Uppercase(void);
	// output operators
	friend std::ostream& operator<<(std::ostream& out, const CMosaikString& ms);
	// is the string empty?
	bool empty(void);
	// clear the string
	bool clear(void);
	// checks values of qualities which shouldn't be larger than 60
	bool CheckQuality( void );

private:
	// our underlying data
	char* mData;
	// the allocated length of the string
	unsigned int mAllocatedLength;
	// the string length
	unsigned int mLength;
	// the packing vector
	static const char* FOUR_BIT_PACKING;
	// the unpacking vector
	static const char* FOUR_BIT_UNPACKING;
};
