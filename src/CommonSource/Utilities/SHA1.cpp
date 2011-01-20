// ***************************************************************************
// CSHA1 - A C++ implementation of the SHA1 hashing algorithm (FIPS 180-1)
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "SHA1.h"

// constructor
CSHA1::CSHA1(void) {}

// destructor
CSHA1::~CSHA1(void) {}

// generates the read group code based on the supplied text
unsigned int CSHA1::GenerateReadGroupCode(const string& readGroupID, const string& sampleName) {

	// concatenate the two strings
	char concatenatedString[MAX_STRING_LENGTH];
	sprintf(concatenatedString, "%s%s", readGroupID.c_str(), sampleName.c_str());
	const unsigned int concatenatedStringLen = (unsigned int)strlen(concatenatedString);

	// initialize our SHA1 object
	CSHA1 sha1;

	// initialize our contexts
	SHA1_Data contexts[2];
	sha1.InitializeSHA1(contexts[0]);
	sha1.InitializeSHA1(contexts[1]);
	contexts[1].bytes = 1;

	// perform a SHA1 hash on the hash string
	sha1.UpdateSHA1(contexts[0], concatenatedString, concatenatedStringLen);
	sha1.UpdateSHA1(contexts[1], concatenatedString, concatenatedStringLen);

	// extract the final combined string from an array of hash private buffers
	unsigned int* pKey = (unsigned int*)sha1.mKey;
	sha1.FinalizeSHA1(contexts[0],     pKey, 5);
	sha1.FinalizeSHA1(contexts[1], pKey + 5, 3);

	//printf("SHA1: { ");
	//for(unsigned char m = 0; m < SHA_HASH_LENGTH; m++) printf("%02X ", sha1.mKey[m]);
	//printf("}\n");

	// copy the first 4 bytes from the hash
	unsigned int readGroupCode;
	memcpy((char*)&readGroupCode, sha1.mKey, 4);

	return readGroupCode;
}

// initialize the SHA values
void CSHA1::InitializeSHA1(SHA1_Data& data) {

	// Set the h-vars to their initial values
	data.iv[0] = 0x67452301;
	data.iv[1] = 0xEFCDAB89;
	data.iv[2] = 0x98BADCFE;
	data.iv[3] = 0x10325476;
	data.iv[4] = 0xC3D2E1F0;

	// Initialise bit count
	data.bytes = 0;

	memset((char*)&data.key, 0, PGP_SHA_BLOCKBYTES);
}

// perform the SHA transformation
void CSHA1::TransformSHA1(unsigned int* block, unsigned int* key) {

	// Set up first buffer
	register unsigned int A = block[0];
	register unsigned int B = block[1];
	register unsigned int C = block[2];
	register unsigned int D = block[3];
	register unsigned int E = block[4];
	register unsigned int t;

	// Heavy mangling, in 4 sub-rounds of 20 interations each.
	subRound( A, B, C, D, E, f1, K2, key[ 0] );
	subRound( E, A, B, C, D, f1, K2, key[ 1] );
	subRound( D, E, A, B, C, f1, K2, key[ 2] );
	subRound( C, D, E, A, B, f1, K2, key[ 3] );
	subRound( B, C, D, E, A, f1, K2, key[ 4] );
	subRound( A, B, C, D, E, f1, K2, key[ 5] );
	subRound( E, A, B, C, D, f1, K2, key[ 6] );
	subRound( D, E, A, B, C, f1, K2, key[ 7] );
	subRound( C, D, E, A, B, f1, K2, key[ 8] );
	subRound( B, C, D, E, A, f1, K2, key[ 9] );
	subRound( A, B, C, D, E, f1, K2, key[10] );
	subRound( E, A, B, C, D, f1, K2, key[11] );
	subRound( D, E, A, B, C, f1, K2, key[12] );
	subRound( C, D, E, A, B, f1, K2, key[13] );
	subRound( B, C, D, E, A, f1, K2, key[14] );
	subRound( A, B, C, D, E, f1, K2, key[15] );
	subRound( E, A, B, C, D, f1, K2, expand(key, 16) );
	subRound( D, E, A, B, C, f1, K2, expand(key, 17) );
	subRound( C, D, E, A, B, f1, K2, expand(key, 18) );
	subRound( B, C, D, E, A, f1, K2, expand(key, 19) );

	subRound( A, B, C, D, E, f2, K3, expand(key, 20) );
	subRound( E, A, B, C, D, f2, K3, expand(key, 21) );
	subRound( D, E, A, B, C, f2, K3, expand(key, 22) );
	subRound( C, D, E, A, B, f2, K3, expand(key, 23) );
	subRound( B, C, D, E, A, f2, K3, expand(key, 24) );
	subRound( A, B, C, D, E, f2, K3, expand(key, 25) );
	subRound( E, A, B, C, D, f2, K3, expand(key, 26) );
	subRound( D, E, A, B, C, f2, K3, expand(key, 27) );
	subRound( C, D, E, A, B, f2, K3, expand(key, 28) );
	subRound( B, C, D, E, A, f2, K3, expand(key, 29) );
	subRound( A, B, C, D, E, f2, K3, expand(key, 30) );
	subRound( E, A, B, C, D, f2, K3, expand(key, 31) );
	subRound( D, E, A, B, C, f2, K3, expand(key, 32) );
	subRound( C, D, E, A, B, f2, K3, expand(key, 33) );
	subRound( B, C, D, E, A, f2, K3, expand(key, 34) );
	subRound( A, B, C, D, E, f2, K3, expand(key, 35) );
	subRound( E, A, B, C, D, f2, K3, expand(key, 36) );
	subRound( D, E, A, B, C, f2, K3, expand(key, 37) );
	subRound( C, D, E, A, B, f2, K3, expand(key, 38) );
	subRound( B, C, D, E, A, f2, K3, expand(key, 39) );

	subRound( A, B, C, D, E, f3, K5, expand(key, 40) );
	subRound( E, A, B, C, D, f3, K5, expand(key, 41) );
	subRound( D, E, A, B, C, f3, K5, expand(key, 42) );
	subRound( C, D, E, A, B, f3, K5, expand(key, 43) );
	subRound( B, C, D, E, A, f3, K5, expand(key, 44) );
	subRound( A, B, C, D, E, f3, K5, expand(key, 45) );
	subRound( E, A, B, C, D, f3, K5, expand(key, 46) );
	subRound( D, E, A, B, C, f3, K5, expand(key, 47) );
	subRound( C, D, E, A, B, f3, K5, expand(key, 48) );
	subRound( B, C, D, E, A, f3, K5, expand(key, 49) );
	subRound( A, B, C, D, E, f3, K5, expand(key, 50) );
	subRound( E, A, B, C, D, f3, K5, expand(key, 51) );
	subRound( D, E, A, B, C, f3, K5, expand(key, 52) );
	subRound( C, D, E, A, B, f3, K5, expand(key, 53) );
	subRound( B, C, D, E, A, f3, K5, expand(key, 54) );
	subRound( A, B, C, D, E, f3, K5, expand(key, 55) );
	subRound( E, A, B, C, D, f3, K5, expand(key, 56) );
	subRound( D, E, A, B, C, f3, K5, expand(key, 57) );
	subRound( C, D, E, A, B, f3, K5, expand(key, 58) );
	subRound( B, C, D, E, A, f3, K5, expand(key, 59) );

	subRound( A, B, C, D, E, f4, K10, expand(key, 60) );
	subRound( E, A, B, C, D, f4, K10, expand(key, 61) );
	subRound( D, E, A, B, C, f4, K10, expand(key, 62) );
	subRound( C, D, E, A, B, f4, K10, expand(key, 63) );
	subRound( B, C, D, E, A, f4, K10, expand(key, 64) );
	subRound( A, B, C, D, E, f4, K10, expand(key, 65) );
	subRound( E, A, B, C, D, f4, K10, expand(key, 66) );
	subRound( D, E, A, B, C, f4, K10, expand(key, 67) );
	subRound( C, D, E, A, B, f4, K10, expand(key, 68) );
	subRound( B, C, D, E, A, f4, K10, expand(key, 69) );
	subRound( A, B, C, D, E, f4, K10, expand(key, 70) );
	subRound( E, A, B, C, D, f4, K10, expand(key, 71) );
	subRound( D, E, A, B, C, f4, K10, expand(key, 72) );
	subRound( C, D, E, A, B, f4, K10, expand(key, 73) );
	subRound( B, C, D, E, A, f4, K10, expand(key, 74) );
	subRound( A, B, C, D, E, f4, K10, expand(key, 75) );
	subRound( E, A, B, C, D, f4, K10, expand(key, 76) );
	subRound( D, E, A, B, C, f4, K10, expandx(key, 77) );
	subRound( C, D, E, A, B, f4, K10, expandx(key, 78) );
	subRound( B, C, D, E, A, f4, K10, expandx(key, 79) );

	// Build message digest
	block[0] += A;
	block[1] += B;
	block[2] += C;
	block[3] += D;
	block[4] += E;
}

// update SHA for a block of data
void CSHA1::UpdateSHA1(SHA1_Data& data, void const *bufIn, size_t len) {

	unsigned char *buf = (unsigned char *)bufIn;

	// Update bitcount
	unsigned int i = (unsigned int)data.bytes % PGP_SHA_BLOCKBYTES;
	data.bytes += (unsigned int)len;

	// i is always less than PGP_SHA_BLOCKBYTES.
	if(PGP_SHA_BLOCKBYTES-i > len) {
		memmove((unsigned char *)data.key + i, buf, len);
		return;
	}

	if(i) {	// First chunk is an odd size
		memcpy((unsigned char *)data.key + i, buf, PGP_SHA_BLOCKBYTES - i);
		SwapBytes(data.key, (unsigned char *)data.key, NUM_SHA_KEY_WORDS);
		TransformSHA1(data.iv, data.key);
		buf += PGP_SHA_BLOCKBYTES-i;
		len -= PGP_SHA_BLOCKBYTES-i;
	}

	// Process data in 64-byte chunks
	while(len >= PGP_SHA_BLOCKBYTES) {
		SwapBytes(data.key, buf, NUM_SHA_KEY_WORDS);
		TransformSHA1(data.iv, data.key);
		buf += PGP_SHA_BLOCKBYTES;
		len -= PGP_SHA_BLOCKBYTES;
	}

	// Handle any remaining bytes of data.
	if(len) memmove(data.key, buf, len);
}

// final wrapup - pad to 64-byte boundary
void CSHA1::FinalizeSHA1(SHA1_Data& data, unsigned int* key, const unsigned char numWords) {

	unsigned int i = (unsigned int)data.bytes % PGP_SHA_BLOCKBYTES;
	unsigned char *p = (unsigned char *)data.key + i; // First unused byte

	// Set the first char of padding to 0x80.  There is always room.
	*p++ = 0x80;

	// Bytes of padding needed to make 64 bytes (0..63)
	i = PGP_SHA_BLOCKBYTES - 1 - i;

	if(i < 8) {	// Padding forces an extra block
		memset(p, 0, i);
		SwapBytes(data.key, (unsigned char *)data.key, 16);
		TransformSHA1(data.iv, data.key);
		p = (unsigned char *)data.key;
		i = 64;
	}

	memset(p, 0, i - 8);
	SwapBytes(data.key, (unsigned char *)data.key, 14);

	// Append length in bits and transform
	data.key[14] = (unsigned int)(data.bytes >> 29);
	data.key[15] = (unsigned int)data.bytes << 3;
	TransformSHA1(data.iv, data.key);

	for(i = 0; i < NUM_SHA_IV_WORDS; i++) data.iv[i] = swap_int32(data.iv[i]);
	memcpy(key, data.iv, numWords * 4);
}

// Shuffle the bytes into big-endian order within words, as per the SHA spec.
void CSHA1::SwapBytes(unsigned int* dest, unsigned char const *src, unsigned char words) {
	do {
		*dest++ = (unsigned int)((unsigned)src[0] << 8 | src[1]) << 16 | ((unsigned)src[2] << 8 | src[3]);
		src += 4;
	} while(--words);
}
