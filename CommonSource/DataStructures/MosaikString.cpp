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

#include "MosaikString.h"

// our packing and unpacking vectors
// ASCII 0-90
//                                                                                                                                                                                       -                                                           A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
const char* CMosaikString::FOUR_BIT_PACKING   = "\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xc\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\xf\x0\xf\x1\xf\xf\xf\x2\xf\xf\xf\x9\xf\x4\xa\xf\xf\xf\x5\x7\x3\xf\xf\x6\xb\x8\xd";
const char* CMosaikString::FOUR_BIT_UNPACKING = "ACGTMRWSYKNX-ZXX";
// Z is used for determing soft-clip locations

// constructor
CMosaikString::CMosaikString(void)
: mData(NULL)
, mAllocatedLength(0)
, mLength(0)
{}

// destructor
CMosaikString::~CMosaikString(void) {
	if( mAllocatedLength > 1 ) delete [] mData;
	else if ( mAllocatedLength == 1 ) delete mData;
	mData = NULL;

	//delete [] FOUR_BIT_PACKING;
	//delete [] FOUR_BIT_UNPACKING;
}

// copy constructor
CMosaikString::CMosaikString(const CMosaikString& ms) 
: mData(NULL)
, mAllocatedLength(ms.mLength + GROWTH_FACTOR)
, mLength(ms.mLength)
{
	mData = new char[mAllocatedLength];
	memcpy(mData, ms.mData, mLength);
	mData[mLength] = 0;
}

// copy constructor
CMosaikString::CMosaikString(const char* c_str)
: mData(NULL)
, mAllocatedLength(0)
, mLength(0)
{
	// don't assign a null character or self assigned string
	if(!c_str || (mData == c_str)) return;

	// copy the string
	mLength = (unsigned int)strlen(c_str);
	mAllocatedLength = mLength + GROWTH_FACTOR;
	mData = new char[mAllocatedLength];
	memcpy(mData, c_str, mLength);
	mData[mLength] = 0;
}

// copy constructor
CMosaikString::CMosaikString(const char* c_str, const unsigned int len)
: mData(NULL)
, mAllocatedLength(0)
, mLength(0)
{
	// don't assign a null character or self assigned string
	if(!c_str || (mData == c_str)) return;

	// copy the string
	mLength = len;
	mAllocatedLength = len + GROWTH_FACTOR;
	mData = new char[mAllocatedLength];
	memcpy(mData, c_str, len);
	mData[len] = 0;
}

// assignment operator
CMosaikString& CMosaikString::operator=(const CMosaikString& ms) {

	// return if this is a self assignment
	if(this == &ms) return *this;

	// skip strings that are zero in length
	if(ms.mLength == 0) {
		mLength = 0;
		return *this;
	}

	// copy the string
	Reserve(ms.mLength);
	mLength = ms.mLength;
	memcpy(mData, ms.mData, mLength);
	mData[mLength] = 0;

	return *this;
}

// assignment operator
CMosaikString& CMosaikString::operator=(const char* c_str) {

	// don't assign a null character or self assigned string
	if(!c_str || (mData == c_str)) return *this;

	// copy the string
	const unsigned int numBytes = (unsigned int)strlen(c_str);
	Reserve(numBytes);
	mLength = numBytes;
	memcpy(mData, c_str, numBytes);
	mData[numBytes] = 0;

	return *this;
}

// larger than operator
bool CMosaikString::operator>(const CMosaikString& ms) const {
	if(mLength == 0)                        return false;
	if(ms.mLength == 0)                     return true;
	if(strcmp(mData, ms.mData) > 0)         return true;
	return false;
}

// less than operator
bool CMosaikString::operator<(const CMosaikString& ms) const {
	//if((mLength == 0) || (ms.mLength == 0)) return true;
	
	if(mLength == 0)                        return true;
	if(ms.mLength == 0)                     return false;
	if(strcmp(mData, ms.mData) < 0)         return true;
	return false;
}

// not equal operator
bool CMosaikString::operator!=(const CMosaikString& ms) const {
	
	if((mLength == 0) || (ms.mLength == 0)) return true;
	if(strcmp(mData, ms.mData) != 0)        return true;
	return false;
}

// equal operator
bool CMosaikString::operator==(const CMosaikString& ms) const {
	if((mLength == 0) || (ms.mLength == 0)) return false;
	if(strcmp(mData, ms.mData) == 0)        return true;
	return false;
}

// element operator
char& CMosaikString::operator[](const int index) {
	return mData[index];
}

std::ostream& operator<<(std::ostream& out, const CMosaikString& ms) {
	if(ms.mLength == 0) return out;
	return out << ms.mData;
}

bool CMosaikString::empty(void) {
	if ( mLength == 0 )
		return true;

	return false;
}

bool CMosaikString::clear(void) {
	//if ( mData ) delete [] mData;
	//mData = NULL;
	memset(mData, 0, mAllocatedLength);
	mLength = 0;
	//mAllocatedLength = 0;

	return true;
}

// appends the specified string to the current string
void CMosaikString::Append(const char* s) {

	// check the allocated room
	const unsigned int suffixLength  = (unsigned int)strlen(s);
	const unsigned int currentLength = mLength;
	unsigned int newLength = suffixLength + currentLength;

	if((newLength + 1) > mAllocatedLength) {

		// save the old string
		char* newData = new char[currentLength + 1];
		memcpy(newData, mData, currentLength);
		newData[currentLength] = 0;

		// copy the old string
		Reserve(newLength);
		memcpy(mData, newData, currentLength);

		// clean up
		delete [] newData;

	}

	// copy the suffix
	memcpy(mData + currentLength, s, suffixLength);

	mLength = newLength;
	mData[newLength] = 0;
}

// append the specified string to the current string
void CMosaikString::Append(const char* s, const unsigned int sLen) {

	// check the allocated room
	const unsigned int suffixLength  = sLen;
	const unsigned int currentLength = mLength;
	unsigned int newLength = suffixLength + currentLength;

	if((newLength + 1) > mAllocatedLength) {

		// save the old string
		char* newData = new char[currentLength + 1];
		memcpy(newData, mData, currentLength);
		newData[currentLength] = 0;

		// copy the old string
		Reserve(newLength);
		memcpy(mData, newData, currentLength);

		// clean up
		delete [] newData;
	}
	
	// copy the suffix
	memcpy(mData + currentLength, s, suffixLength);

	mLength = newLength;
	mData[newLength] = 0;

}

// copies the specified c-style string
void CMosaikString::Copy(const char* string, const unsigned int numBytes) {
	Reserve(numBytes);
	memcpy(mData, string, numBytes);
	mData[numBytes] = 0;
	mLength = numBytes;
}

// returns a pointer to the data
char* CMosaikString::Data(void) {
	return mData;
}

// returns a const pointer to the data
const char* CMosaikString::CData(void) const {
	return (const char*)mData;
}

// decrements each character in the string by the specified amount
void CMosaikString::Decrement(const char amount) {
	for(unsigned int i = 0; i < mLength; i++) mData[i] -= amount;
}

// fills the string with numBytes copies of the ch
void CMosaikString::Fill(const char ch, const unsigned int numBytes) {
	Reserve(numBytes);
	memset(mData, ch, numBytes);
	mData[numBytes] = 0;
	mLength = numBytes;
}

// returns the hash value
size_t CMosaikString::GetHash(void) const {

	// N.B. caching the hash value doesn't help
	// Hsieh hash: slightly faster than Java string hash
	unsigned int hash = 0;
	if(mLength == 0) return hash;

	const char* data = mData;
	unsigned int len = mLength;

	unsigned int rem = len & 3;
	len >>= 2;

	unsigned int tmp;
	for(; len > 0; len--) {
		hash  += get16bits(data);
		tmp    = (get16bits(data + 2) << 11) ^ hash;
		hash   = (hash << 16) ^ tmp;
		data  += 4;
		hash  += hash >> 11;
	}

	switch(rem) {
		case 3:
			hash += get16bits(data);
			hash ^= hash << 16;
			hash ^= data[2] << 18;
			hash += hash >> 11;
			break;
		case 2:
			hash += get16bits(data);
			hash ^= hash << 11;
			hash += hash >> 17;
			break;
		case 1:
			hash += *data;
			hash ^= hash << 10;
			hash += hash >> 1;
	}

	hash ^= hash << 3;
	hash += hash >> 5;
	hash ^= hash << 4;
	hash += hash >> 17;
	hash ^= hash << 25;
	hash += hash >> 6;

	return (size_t)hash;
}

// checks values of qualities which shouldn't be larger than 60
void CMosaikString::CheckQuality( void ) {
	for(unsigned int i = 0; i < mLength; i++) {
		if ( mData[i] > 60 ) {
			printf("ERROR: The base quality is larger than 60.\n");
			exit(1);
		}
	}
}

// increments each character in the string by the specified amount
void CMosaikString::Increment(const char amount) {
	for(unsigned int i = 0; i < mLength; i++) mData[i] += amount;
}

// joins two strings (used by the Smith-Waterman caching algorithm)
void CMosaikString::Join(const char* s1, const unsigned int s1Length, const char* s2, const unsigned int s2Length) {

	const unsigned int JOINED_LENGTH = s1Length + s2Length + 1;
	Reserve(JOINED_LENGTH);

	unsigned int bufferOffset = s1Length;
	memcpy(mData, s1, s1Length);
	mData[bufferOffset++] = '+';
	memcpy(mData + bufferOffset, s2, s2Length);

	mLength = JOINED_LENGTH;
	mData[mLength] = 0;
}

// returns the size of the data
unsigned int CMosaikString::Length(void) const {
	return mLength;
}

// packs both the original bases and the supplied bases into 4-bit notation
void CMosaikString::Pack(const CMosaikString& ms) {

	// make sure that both strings are the same length
	if(mLength != ms.mLength) {
		printf("ERROR: Both strings must be the same length for 4-bit packing to succeed\n");
		exit(1);
	}

	// OBS: FOUR_BIT_PACKING translations have been checked

	// pack the current string
	for(unsigned int i = 0; i < mLength; i++) {
		mData[i] = FOUR_BIT_PACKING[mData[i]] | (FOUR_BIT_PACKING[ms.mData[i]] << 4);
	}
}

// prepends the specified string to the current string
void CMosaikString::Prepend(const CMosaikString& ms) {
	Prepend(ms.mData, ms.mLength);
}

// prepends the specified string to the current string
void CMosaikString::Prepend(const char* s, const unsigned int sLen) {

	// check the allocated room
	const unsigned int prefixLength  = sLen;
	const unsigned int currentLength = mLength;
	unsigned int newLength = prefixLength + currentLength;

	if((newLength + 1) > mAllocatedLength) {

		// save the old string
		char* newData = new char[currentLength + 1];
		memcpy(newData, mData, currentLength);
		newData[currentLength] = 0;

		// copy the prefix
		Reserve(newLength);
		memcpy(mData, s, prefixLength);

		// copy the old string
		memcpy(mData + prefixLength, newData, currentLength);

		// clean up
		delete [] newData;

	} else {

		memmove_s(mData + prefixLength, mAllocatedLength - prefixLength, mData, currentLength);
		memcpy(mData, s, prefixLength);
	}

	mLength = newLength;
	mData[newLength] = 0;
}

// removes all occurrences of the specified character
void CMosaikString::Remove(const char ch) {

	unsigned int currentRemovePosition = 0;
	for(unsigned int i = 0; i < mLength; i++) {

		// copy the characters if the coordinates are different
		if(i != currentRemovePosition) mData[currentRemovePosition] = mData[i];

		// increment the remove position
		if(mData[i] != ch) currentRemovePosition++;
	}

	mLength = currentRemovePosition;
	mData[mLength] = 0;
}

// replaces all occurrences of the first parameter with the second parameter
void CMosaikString::Replace(const char oldCh, const char newCh) {
	for(unsigned int i = 0; i < mLength; i++) {
		if(mData[i] == oldCh) mData[i] = newCh;
	}
}

// reserve the specified number of bytes
void CMosaikString::Reserve(const unsigned int numBytes) {

	if((numBytes + 1) > mAllocatedLength) {
		mAllocatedLength = numBytes + GROWTH_FACTOR + 1;
		if(mData) delete [] mData;
		mData = new char[mAllocatedLength];
	}

	// reset the data
	mLength = 0;
	//mData[0] = 0;
	memset(mData, 0, mAllocatedLength);
}

// reverses the contents of the string
void CMosaikString::Reverse(void) {
	for(unsigned int y = mLength; y >= (mLength / 2) + 1; y--)
		swapByte(mData[mLength - y], mData[y - 1]);
}

// Performs an in-place reverse complement conversion
void CMosaikString::ReverseComplement(void) {

	for(unsigned int y = mLength; y >= (mLength / 2) + 1; y--)
		swapByte(mData[mLength - y], mData[y - 1]);

	for(unsigned int i = 0; i < mLength; i++) {
		switch(mData[i]) {
			case 'A':
				mData[i] = 'T';
				break;
			case 'C':
				mData[i] = 'G';
				break;
			case 'G':
				mData[i] = 'C';
				break;
			case 'T':
				mData[i] = 'A';
				break;
			default:
				break;
		}
	}
}

// sets the length to the specified size
void CMosaikString::SetLength(const unsigned int length) {

	if(mAllocatedLength == 0) {
		mLength = 0;
		return;
	}

	if(length < mAllocatedLength) mLength = length;
	else mLength = mAllocatedLength - 1;
	mData[mLength] = 0;
}

// retrieves the specified substring
// TODO: add some bounds-checking here
std::string CMosaikString::Substring(const unsigned int position, const unsigned int length) const {
	std::string s;
	s.resize(length);
	char* pS = (char*)s.data();
	memcpy(pS, mData + position, length);
	return s;
}

// trims the first specified number of bytes
void CMosaikString::TrimBegin(unsigned int numBytes) {

	if(numBytes > mLength) {
		mLength = 0;
		mData[0] = 0;
		return;
	}

	unsigned int newLength = mLength - numBytes;
	memmove_s(mData, mAllocatedLength, mData + numBytes, newLength);
	mLength = newLength;
	mData[mLength] = 0;
}

// trims the last specified number of bytes
void CMosaikString::TrimEnd(unsigned int numBytes) {

	if(numBytes > mLength) {
		mLength = 0;
		mData[0] = 0;
		return;
	}

	mLength = mLength - numBytes;
	mData[mLength] = 0;
}

// unpacks both the original bases and the supplied bases from 4-bit notation
void CMosaikString::Unpack(CMosaikString& ms) {

	// make sure that the supplied string is the same length
	if(mLength != ms.mLength) {
		ms.Reserve(mLength);
		ms.SetLength(mLength);
	}

	// OBS: FOUR_BIT_UNPACKING translations have been checked

	// unpack the current string
	// PACK_MASK = 15 = 0x1111
	for(unsigned int i = 0; i < mLength; i++) {		
		ms.mData[i] = FOUR_BIT_UNPACKING[(mData[i] >> 4) & PACK_MASK];
		mData[i]    = FOUR_BIT_UNPACKING[mData[i] & PACK_MASK];
	}
}

// converts the string to uppercase
void CMosaikString::Uppercase(void) {
	for(unsigned int i = 0; i < mLength; i++) mData[i] = toupper(mData[i]);
}
