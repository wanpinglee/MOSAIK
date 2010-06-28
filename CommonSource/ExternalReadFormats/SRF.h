// ***************************************************************************
// CSRF - imports reads from the SRF file format.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg [Adapted from io_lib: James Bonfield]
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <assert.h>
#include <cmath>
#include <errno.h>
#include "LargeFileSupport.h"
#include "Read.h"
#include "SafeFunctions.h"

using namespace std;

#define ZTR_TYPE_BASE	0x42415345
#define ZTR_TYPE_CNF4	0x434e4634
#define ZTR_TYPE_HUFF   0x48554646

#define SRFB_CONTAINER    'S'
#define SRFB_XML          'X'
#define SRFB_TRACE_HEADER 'H'
#define SRFB_TRACE_BODY   'R'
#define SRFB_INDEX        'I'

#define SRF_MAGIC   "SSRF"
#define SRF_VERSION "1.3"

#define MF_READ    1
#define MF_WRITE   2
#define MF_APP     4
#define MF_BINARY  8
#define MF_TRUNC  16

#define CODE_USER	128
#define SYM_EOF     256
#define NAMELEN     512
#define CHECK_BYTES  12

#define ZTR_MAGIC         "\256ZTR\r\n\032\n"
#define ZTR_VERSION_MAJOR 1
#define ZTR_VERSION_MINOR 2

#define ZTR_FORM_RAW     0
#define ZTR_FORM_XRLE2   4
#define ZTR_FORM_STHUFF 77
#define ZTR_FORM_QSHIFT 79

#define iswap_int4(x) \
	(((x & 0x000000ff) << 24) + \
	((x & 0x0000ff00) <<  8) + \
	((x & 0x00ff0000) >>  8) + \
	((x & 0xff000000) >> 24))

#ifdef SP_BIG_ENDIAN
#define be_int8(x) (x)
#define be_int4(x) (x)
#define be_int2(x) (x)
#define be_int1(x) (x)

#define le_int8(x) iswap_int8((x))
#define le_int4(x) iswap_int4((x))
#define le_int2(x) iswap_int2((x))
#define le_int1(x) (x)
#else
#define be_int8(x) iswap_int8((x))
#define be_int4(x) iswap_int4((x))
#define be_int2(x) iswap_int2((x))
#define be_int1(x) (x)

#define le_int8(x) (x)
#define le_int4(x) (x)
#define le_int2(x) (x)
#define le_int1(x) (x)
#endif

class CSRF {
public:
	// constructor
	CSRF(void);
	// destructor
	~CSRF(void);
	// validates the supplied SRF file
	static bool CheckFile(const string& filename, const bool showError);
	// closes the SRF file
	void Close(void);
	// gets the next read from the SRF file
	bool GetRead(Mosaik::Read& mr);
	// opens the SRF file
	void Open(const string& filename);
private:
	// container header
	struct ContainerHeader {
		char block_type;
		char version[256];
		char container_type;
		char base_caller[256];
		char base_caller_version[256];
	};
	// trace header
	struct TraceHeader {
		char block_type; 
		char read_prefix_type;
		char id_prefix[256];
		unsigned int trace_hdr_size;
		unsigned char *trace_hdr;

		TraceHeader()
			: trace_hdr(NULL)
		{}
	};
	// trace body
	struct TraceBody {
		char block_type;
		char read_id[256];
		unsigned char flags;
		unsigned int trace_size;
		unsigned char *trace;

		TraceBody()
			: trace(NULL)
		{}
	};
	// the ZTR header
	struct ZtrHeader {
		unsigned char  magic[8];	  // 0xae5a54520d0a1a0a (be)
		unsigned char  version_major; // ZTR_VERSION_MAJOR
		unsigned char  version_minor; // ZTR_VERSION_MINOR
	};
	// the ZTR chunk
	struct ZtrChunk {
		unsigned int type;		// chunk type (be)
		unsigned int mdlength;	// length of meta data field (be)
		char *mdata;			// meta data
		unsigned int dlength;	// length of data field (be)
		char *data;				// a format byte and the data itself
		int ztr_owns;			// boolean: true if we can free (meta)data

		ZtrChunk()
			: mdata(NULL)
			, data(NULL)
		{}
	};
	// A single symbol and it's encoding
	struct huffman_code_t {
		signed int symbol; // 0-255 character, 256 = exception code, 257 = EOF
		int nbits;
		unsigned int code;
		int freq;
	};
	// A collection of huffman_code_t along with decoding optimisations
	struct huffman_codes_t {
		huffman_code_t *codes;
		int ncodes;
		int codes_static;
		huffman_code_t lookup[258]; // Mapping of symbol character to code
		int max_code_len;

		huffman_codes_t()
			: codes(NULL)
		{}
	};
	// Use for store_bits() and GetBits()
	struct block_t {
		unsigned char *data;
		size_t alloc;
		size_t byte;
		int bit;

		block_t()
			: data(NULL)
		{}
	};
	// Byte-wise jumping table
	struct h_jump4_t {	
		unsigned short jump;
		unsigned char symbol[4];
		unsigned char nsymbols;
		unsigned char top_bit;   // bit 9 of symbol[]
	};
	// Tree and jump-table data structures used for fast decoding.
	struct htree_t {
		unsigned short c[2]; // child node
		signed short l[2]; // symbol to emit on transition. -1 => none
	};
	// A collection of huffman_codes_t, for use with the multi-code codec
	struct huffman_codeset_t {
		huffman_codes_t **codes;
		int ncodes;
		int code_set; // (128-255) The user specified number for this encoding

		// Cached binary version of codeset, assumes last block
		block_t *blk;
		int      bit_num; // if 1st block, which bit will stored codes end on

		// Cache DecodeHuffmanStream parameters
		h_jump4_t (*decode_J4)[16];
		htree_t *decode_t;

		huffman_codeset_t()
			: codes(NULL)
			, blk(NULL)
			, decode_t(NULL)
		{}
	};
	// the ZTR Huffman codeset
	struct ztr_hcode_t {
		int ztr_owns; // true is ZTR is to free the data later
		huffman_codeset_t* codes;

		ztr_hcode_t()
			: codes(NULL)
		{}
	};
	// the main ZTR structure
	struct ZTR_t {
		// General bits to do with the ZTR file format
		ZtrHeader header;
		ZtrChunk* chunk;
		int nchunks;

		// Cached huffman encoding/decoding tables for STHUFF format
		ztr_hcode_t* hcodes;
		int nhcodes;
		int hcodes_checked;

		ZTR_t()
			: chunk(NULL)
			, hcodes(NULL)
			, hcodes_checked(0)
		{}
	};
	// the extended file
	struct mFILE {
		FILE *fp;
		char *data;
		size_t alloced;
		int eof;
		int mode; /* open mode in MF_?? define bit pattern */
		size_t size;
		size_t offset;
		size_t flush_pos;

		mFILE()
			: fp(NULL)
			, data(NULL)
		{}
	};
	// master SRF object
	struct SRF_t {
		FILE *fp;

		// Cached copies of each of the most recent chunk types loaded
		ContainerHeader ch;
		TraceHeader     th;
		TraceBody       tb;

		// Private: cached data for use by FetchNextTrace
		ZTR_t *ztr;
		mFILE *mf;
		long mf_pos, mf_end;

		SRF_t()
			: fp(NULL)
			, ztr(NULL)
			, mf(NULL)
		{}

	} mSrfData;

	// ====================
	// compression routines
	// ====================

	// reorders quality data from an interleaved 4-byte aligned format to its RAW format 
	char* DeinterleaveQualityData(char* qold, int qlen, int* new_len);
	// reverses multi-byte run length encoding
	char* ExpandMultiByteRLE(char* comp, int comp_len, int* uncomp_len);
	// implements decompression using a set of static huffman codes stored using the Deflate algorithm
	char* InflateStaticHuffman(ZTR_t* ztr, char* comp, int comp_len, int* uncomp_len);

	// ================
	// Huffman routines
	// ================

	// allocates and returns a new block_t struct of a specified default size
	block_t* CreateBlock(unsigned char* data, size_t size);
	// decode a huffman stream from 'block' using huffman codes 'c'
	block_t* DecodeHuffmanStream(block_t* in, huffman_codeset_t* cs);
	// deallocates memory created by CreateBlock()
	void DestroyBlock(block_t* blk, int keep_data);
	// generates canonical huffman codes given a set of symbol bit lengths
	int GenerateCanonicalHuffmanCodes(huffman_codes_t* c);
	// reads up to 24-bits worth of data and returns
	static signed int GetBits(block_t* block, int nbits);
	// gets Huffman codes from the stream
	huffman_codeset_t* GetHuffmanCodes(block_t* block, int* bfinal);
	// it restores huffman_codes_t structs from the a serialised data stream
	huffman_codes_t* GetHuffmanCodesSingle(block_t* block);
	// A slow version of the above huffman_decode function
	int GetNextSymbol(block_t* in, int* htab);
	// initialize the huffman tables for decoding
	int InitializeDecodeHuffmanTables(huffman_codeset_t* cs);
	// ensures a block_t holds at least 'size' bytes
	int ResizeBlock(block_t* blk, size_t size);
	// reverses the order of bits in the bottom nbits of val
	unsigned int ReverseBitOrder(unsigned int val, int nbits);
	// stores nbytes bytes, padding to align on the next byte boundary
	void StoreBytes(block_t* block, unsigned char* val, int nbytes);

	// ====================
	// memory file routines
	// ====================

	// for creating existing mFILE pointers directly from memory buffers
	mFILE* mfcreate(char* data, int size);
	// memory fread
	size_t mfread(void* ptr, size_t size, size_t nmemb, mFILE* mf);
	// memory rewind
	void mfrewind(mFILE* mf);
	// memory fseek
	int mfseek(mFILE* mf, long offset, int whence);
	// memory ftell
	long mftell(mFILE* mf);
	// memory ftruncate
	void mftruncate(mFILE* mf, long offset);
	// memory fwrite
	size_t mfwrite(void* ptr, size_t size, size_t nmemb, mFILE* mf);

	// ============
	// SRF routines
	// ============

	// decodes a partial ZTR file consisting of data in 'mf'
	ZTR_t* DecodePartialZtr(SRF_t* srf, mFILE* mf, ZTR_t* z);
	// fetches the next trace from an SRF container as a ZTR object
	ZTR_t* FetchNextTrace(SRF_t* srf, char* name);
	// returns the type of the next block
	int GetBlockType(SRF_t* srf);
	// reads a container header and stores the result in 'ch'
	int ReadContainerHeader(SRF_t* srf, ContainerHeader* ch);
	// reads a pascal-style string from the srf file
	int ReadPascalString(SRF_t* srf, char *str);
	// reads a trace header + trace 'blob' and stores the result in 'th'
	int ReadTraceBody(SRF_t* srf, TraceBody* tb, int no_trace);
	// reads a data header and stores the result in 'th'
	int ReadTraceHeader(SRF_t* srf, TraceHeader* th);
	// read unsigned 32-bit values in big-endian format
	int ReadUInt32(SRF_t* srf, unsigned int* val);

	// ============
	// ZTR routines
	// ============

	// adds a user-defined huffman_codeset_t code-set to the available code sets used by huffman_decode
	ztr_hcode_t* AddUserDefinedHuffmanCodes(ZTR_t *ztr, huffman_codeset_t* codes, int ztr_owns);
	// delete ZTR
	void DeleteZtr(ZTR_t* ztr);
	// creates a copy of ZTR_t 'src'
	ZTR_t* DuplicateZtr(ZTR_t* src);
	// searches for chunks of a specific type
	ZtrChunk** FindZtrChunks(ZTR_t* ztr, unsigned int type, int* nchunks_p);
	// allocates and initialises a ZTR_t structure
	ZTR_t* NewZtr(void);
	// reads a ZTR chunk header and metadata, but not the main data segment
	ZtrChunk* ReadZtrChunkHeader(mFILE* mf);
	// reads a ZTR file header
	int ReadZtrHeader(mFILE* mf, ZtrHeader* h);
	// searches through the cached huffman_codeset_t tables looking for a stored huffman code of type 'code_set'
	ztr_hcode_t* SearchCachedHuffmanTables(ZTR_t* ztr, int code_set);
	// uncompresses an individual chunk from all levels of compression
	int UncompressChunk(ZTR_t* ztr, ZtrChunk* chunk);
	// uncompresses a ztr (in memory)
	int UncompressZtr(ZTR_t* ztr);

	// ========================
	// private member variables
	// ========================

	// denotes the status of the input stream
	bool mIsOpen;
	// our pre-computed base quality LUT
	char mBaseQualityLUT[256];
};
