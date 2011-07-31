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
#include "BamHeader.h"
#include "MdTager.h"
#include "ReadGroup.h"
#include "ReferenceSequence.h"
#include "SequencingTechnologies.h"

using namespace std;

// define our sorting types
typedef unsigned char SortOrderType;
const SortOrderType SORTORDER_UNSORTED = 0;
const SortOrderType SORTORDER_READNAME = 10;
const SortOrderType SORTORDER_POSITION = 20;

// program group in header
struct ProgramGroup {
	string ID;  // program name
	string VN;  // program version
	string CL;  // command line
};

// define our BAM header structure
struct BamHeader {
	SortOrderType SortOrder;
	string Version;
	vector<MosaikReadFormat::ReadGroup>* pReadGroups;
	vector<ReferenceSequence>* pReferenceSequences;
	ProgramGroup pg;

	// constructor
	BamHeader(void)
		: SortOrder(SORTORDER_UNSORTED)
		, Version("1.0")
		, pReadGroups(NULL)
		, pReferenceSequences(NULL)
	{}
};

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

	// copy constructor
	BGZF ( const BGZF & copy ) {
		CompressedBlockSize   = copy.CompressedBlockSize;
		UncompressedBlockSize = copy.UncompressedBlockSize;
		CompressedBlock   = new char[ CompressedBlockSize ];
		UncompressedBlock = new char[ UncompressedBlockSize ];
		memcpy( CompressedBlock, copy.CompressedBlock, CompressedBlockSize );
		memcpy( UncompressedBlock, copy.UncompressedBlock, UncompressedBlockSize );

	}

	// assign operator
	BGZF& operator=( const BGZF & copy ) {
		CompressedBlockSize    = copy.CompressedBlockSize;
		UncompressedBlockSize  = copy.UncompressedBlockSize;
		char* temp_CompressedBlock   = new char[ CompressedBlockSize ];
		char* temp_UncompressedBlock = new char[ UncompressedBlockSize ];
		memcpy( temp_CompressedBlock, copy.CompressedBlock, CompressedBlockSize );
		memcpy( temp_UncompressedBlock, copy.UncompressedBlock, UncompressedBlockSize );
		delete [] CompressedBlock;
		delete [] UncompressedBlock;
		CompressedBlock   = temp_CompressedBlock;
		UncompressedBlock = temp_UncompressedBlock;

		// unsafe
		Stream = copy.Stream;

		return *this;

	}

};

class CBamWriter {
public:
	// constructor
	//CBamWriter(void);
	// destructor
	~CBamWriter(void);
	// closes the alignment archive
	void Close(void);
	// opens the alignment archive
	void Open(const string& filename, const BamHeader& header);
	// saves the alignment to the alignment archive
	void SaveAlignment(const Alignment al, const char* zaString, const bool& noCigarMdNm, const bool& notShowRnamePos, const bool& isSolid, const bool processedBamData = false );
	// saves the reference and position of an alignment to the alignment archive
	void SaveReferencePosition( const unsigned int refIndex, const unsigned int refBegin, const unsigned int refEnd );
	// creates a packed cigar string from the supplied alignment
	void CreatePackedCigar(const Alignment& al, string& packedCigar, unsigned short& numCigarOperations, const bool isSolid );
	// encodes the supplied query sequence into 4-bit notation
	void EncodeQuerySequence(const CMosaikString& query, string& encodedQuery);
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
	//static void CreatePackedCigar(const Alignment& al, string& packedCigar, unsigned int& numCigarOperations, const bool& isSolid );
	// encodes the supplied query sequence into 4-bit notation
	//static void EncodeQuerySequence(const CMosaikString& query, string& encodedQuery);
	// MD tager
	CMdTager mdTager;
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
