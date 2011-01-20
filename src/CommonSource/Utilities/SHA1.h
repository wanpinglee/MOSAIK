// ***************************************************************************
// CSHA1 - A C++ implementation of the SHA1 hashing algorithm (FIPS 180-1)
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

using namespace std;

#define NUM_SHA_IV_WORDS     5
#define PGP_SHA_BLOCKBYTES  64
#define NUM_SHA_KEY_WORDS   16
#define SHA_HASH_LENGTH     20
#define SHA_KEY_LENGTH      32
#define MAX_STRING_LENGTH  512

// The SHA f()-functions
#define f1(x,y,z) (z ^ (x & (y ^ z)))			// Rounds  0 - 19
#define f2(x,y,z) (x ^ y ^ z)					// Rounds 20 - 39
#define f3(x,y,z) ((x & y) + (z & (x ^ y) ))	// Rounds 40 - 59
#define f4(x,y,z) (x ^ y ^ z)					// Rounds 60 - 79

// The SHA Mysterious Constants.
#define K2	0x5A827999L	// Rounds  0-19 - floor(sqrt(2) * 2^30)
#define K3	0x6ED9EBA1L	// Rounds 20-39 - floor(sqrt(3) * 2^30)
#define K5	0x8F1BBCDCL	// Rounds 40-59 - floor(sqrt(5) * 2^30)
#define K10	0xCA62C1D6L	// Rounds 60-79 - floor(sqrt(10) * 2^30)

// 32-bit rotate left
#define ROTL(n,X) ((X << n) | (X >> (32-n)))

// The prototype SHA sub-round
#define subRound(a, b, c, d, e, f, k, data) \
	(e += ROTL(5,a) + f(b, c, d) + k + data, b = ROTL(30, b))

// The initial expanding function
#define expandx(W,i) (t = W[i&15] ^ W[(i-14)&15] ^ W[(i-8)&15] ^ W[(i-3)&15], ROTL(1, t))
#define expand(W,i) (W[i&15] = expandx(W,i))

#define swap_int32(x) \
	(((x & 0x000000ff) << 24) + \
	((x & 0x0000ff00) <<  8) + \
	((x & 0x00ff0000) >>  8) + \
	((x & 0xff000000) >> 24))

class CSHA1 {
public:
	// constructor
	CSHA1(void);
	// destructor
	~CSHA1(void);
	// generates the read group code based on the supplied text
	static unsigned int GenerateReadGroupCode(const string& readGroupID, const string& sampleName);

private:
	// define our SHA context
	struct SHA1_Data {
		unsigned int key[NUM_SHA_KEY_WORDS];
		unsigned int iv[NUM_SHA_IV_WORDS];
		unsigned int bytes;
	};
	// final wrapup - pad to 64-byte boundary
	static void FinalizeSHA1(SHA1_Data& data, unsigned int* key, const unsigned char numWords);
	// initialize the SHA values
	static void InitializeSHA1(SHA1_Data& data);
	// Shuffle the bytes into big-endian order within words, as per the SHA spec.
	static void SwapBytes(unsigned int* dest, unsigned char const *src, unsigned char words);
	// perform the SHA transformation
	static void TransformSHA1(unsigned int* block, unsigned int* key);
	// update SHA for a block of data
	static void UpdateSHA1(SHA1_Data& data, void const *bufIn, size_t len);
	// our passphrase key
	unsigned char mKey[SHA_KEY_LENGTH];
};
