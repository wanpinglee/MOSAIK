// ***************************************************************************
// CBamWriter - exports alignment data into the BAM file format.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <sstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <zlib.h>
#include "Alignment.h"
#include "ReadGroup.h"
#include "ReferenceSequence.h"
#include "SequencingTechnologies.h"

using namespace std;

// define our sorting types
typedef unsigned char SortOrderType;
const SortOrderType SORTORDER_UNSORTED = 0;
const SortOrderType SORTORDER_READNAME = 10;
const SortOrderType SORTORDER_POSITION = 20;

// define our BAM header structure
struct BamHeader {
	SortOrderType SortOrder;
	string Version;
	vector<MosaikReadFormat::ReadGroup>* pReadGroups;
	vector<ReferenceSequence>* pReferenceSequences;

	// constructor
	BamHeader(void)
		: SortOrder(SORTORDER_UNSORTED)
		, Version("1.0")
		, pReadGroups(NULL)
		, pReferenceSequences(NULL)
	{}
};

// our zlib constants
#define GZIP_ID1             31
#define GZIP_ID2            139
#define CM_DEFLATE            8
#define FLG_FEXTRA            4
#define OS_UNKNOWN          255
#define BGZF_XLEN             6
#define BGZF_ID1             66
#define BGZF_ID2             67
#define BGZF_LEN              2
#define GZIP_WINDOW_BITS    -15
#define Z_DEFAULT_MEM_LEVEL   8

// our BZGF constants
#define BLOCK_HEADER_LENGTH    18
#define BLOCK_FOOTER_LENGTH     8
#define MAX_BLOCK_SIZE      65536
#define DEFAULT_BLOCK_SIZE  65536

// our BAM constants
#define BAM_CORE_SIZE  32
#define BAM_CMATCH      0
#define BAM_CINS        1
#define BAM_CDEL        2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD        6
#define BAM_CIGAR_SHIFT 4

#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)

// define some SAM/BAM FLAGS
#define BAM_SEQUENCED_AS_PAIRS         1
#define BAM_PROPER_PAIR                2
#define BAM_QUERY_UNMAPPED             4
#define BAM_MATE_UNMAPPED              8
#define BAM_QUERY_REVERSE_COMPLEMENT  16
#define BAM_MATE_REVERSE_COMPLEMENT   32
#define BAM_QUERY_FIRST_MATE          64
#define BAM_QUERY_SECOND_MATE        128
#define BAM_SECONDARY_ALIGNMENT      256
#define BAM_FAILS_PLATFORM_QUALITY   512
#define BAM_PCR_OPTICAL_DUPLICATE   1024

// define some tag lengths
#define MISMATCH_TAG_LEN 7

// define our BZGF structure
struct BGZF {
	unsigned int UncompressedBlockSize;
	unsigned int CompressedBlockSize;
	unsigned int BlockLength;
	unsigned int BlockOffset;
	uint64_t BlockAddress;
	bool IsOpen;
	FILE* Stream;
	char* UncompressedBlock;
	char* CompressedBlock;

	// constructor
	BGZF(void)
		: UncompressedBlockSize(DEFAULT_BLOCK_SIZE)
		, CompressedBlockSize(MAX_BLOCK_SIZE)
		, BlockLength(0)
		, BlockOffset(0)
		, BlockAddress(0)
		, IsOpen(false)
		, Stream(NULL)
		, UncompressedBlock(NULL)
		, CompressedBlock(NULL)
	{
		try {
			CompressedBlock   = new char[CompressedBlockSize];
			UncompressedBlock = new char[UncompressedBlockSize];
		} catch(bad_alloc&) {
			printf("ERROR: Unable to allocate memory for our BGZF object.\n");
			exit(1);
		}
	}

	// destructor
	~BGZF(void) {
		if(CompressedBlock)   delete [] CompressedBlock;
		if(UncompressedBlock) delete [] UncompressedBlock;
	}
};

class CBamWriter {
public:
	// constructor
	CBamWriter(void);
	// destructor
	~CBamWriter(void);
	// closes the alignment archive
	void Close(void);
	// opens the alignment archive
	void Open(const string& filename, const BamHeader& header);
	// saves the alignment to the alignment archive
	void SaveAlignment(const CMosaikString& readName, const string& readGroupID, const vector<Alignment>::iterator& alIter);
private:
	// closes the BAM file
	void BgzfClose(void);
	// compresses the current block
	int BgzfDeflateBlock(void);
	// flushes the data in the BGZF block
	void BgzfFlushBlock(void);
	// opens the BAM file for writing
	void BgzfOpen(const string& filename);
	// packs an unsigned integer into the specified buffer
	static inline void BgzfPackUnsignedInt(char* buffer, unsigned int value);
	// packs an unsigned short into the specified buffer
	static inline void BgzfPackUnsignedShort(char* buffer, unsigned short value);
	// writes the supplied data into the BGZF buffer
	unsigned int BgzfWrite(const char* data, const unsigned int dataLen);
	// calculates the minimum bin that contains a region [begin, end)
	static inline unsigned int CalculateMinimumBin(unsigned int begin, unsigned int end);
	// creates a packed cigar string from the supplied alignment
	static void CreatePackedCigar(const Alignment& al, string& packedCigar, unsigned int& numCigarOperations);
	// encodes the supplied query sequence into 4-bit notation
	static void EncodeQuerySequence(const CMosaikString& query, string& encodedQuery);
	// our BGZF output object
	BGZF mBGZF;
};

// packs an unsigned integer into the specified buffer
inline void CBamWriter::BgzfPackUnsignedInt(char* buffer, unsigned int value) {
	buffer[0] = (char)value;
	buffer[1] = (char)(value >> 8);
	buffer[2] = (char)(value >> 16);
	buffer[3] = (char)(value >> 24);
}

// packs an unsigned short into the specified buffer
inline void CBamWriter::BgzfPackUnsignedShort(char* buffer, unsigned short value) {
	buffer[0] = (char)value;
	buffer[1] = (char)(value >> 8);
}

// calculates the minimum bin that contains a region [begin, end)
inline unsigned int CBamWriter::CalculateMinimumBin(unsigned int begin, unsigned int end) {
	--end;
	if((begin >> 14) == (end >> 14)) return 4681 + (begin >> 14);
	if((begin >> 17) == (end >> 17)) return  585 + (begin >> 17);
	if((begin >> 20) == (end >> 20)) return   73 + (begin >> 20);
	if((begin >> 23) == (end >> 23)) return    9 + (begin >> 23);
	if((begin >> 26) == (end >> 26)) return    1 + (begin >> 26);
	return 0;
}
