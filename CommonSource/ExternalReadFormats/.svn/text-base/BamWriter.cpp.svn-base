// ***************************************************************************
// CBamWriter - exports alignment data into the BAM file format.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "BamWriter.h"

// constructor
CBamWriter::CBamWriter(void)
{}

// destructor
CBamWriter::~CBamWriter(void) {
	if(mBGZF.IsOpen) BgzfClose();
}

// closes the BAM file
void CBamWriter::BgzfClose(void) {

	mBGZF.IsOpen = false;

	// flush the BGZF block
	BgzfFlushBlock();

	// add an empty block
	mBGZF.BlockOffset = 0;
	int blockLength = BgzfDeflateBlock();
	fwrite(mBGZF.CompressedBlock, 1, blockLength, mBGZF.Stream);

	// flush and close
	fflush(mBGZF.Stream);
	fclose(mBGZF.Stream);
}

// compresses the current block
int CBamWriter::BgzfDeflateBlock(void) {

	// initialize the gzip header
	char* buffer = mBGZF.CompressedBlock;
	unsigned int bufferSize = mBGZF.CompressedBlockSize;

	memset(buffer, 0, 18);
	buffer[0]  = GZIP_ID1;
	buffer[1]  = (char)GZIP_ID2;
	buffer[2]  = CM_DEFLATE;
	buffer[3]  = FLG_FEXTRA;
	buffer[9]  = (char)OS_UNKNOWN;
	buffer[10] = BGZF_XLEN;
	buffer[12] = BGZF_ID1;
	buffer[13] = BGZF_ID2;
	buffer[14] = BGZF_LEN;

	// loop to retry for blocks that do not compress enough
	int inputLength = mBGZF.BlockOffset;
	int compressedLength = 0;

	while(true) {

		z_stream zs;
		zs.zalloc    = NULL;
		zs.zfree     = NULL;
		zs.next_in   = (Bytef*)mBGZF.UncompressedBlock;
		zs.avail_in  = inputLength;
		zs.next_out  = (Bytef*)&buffer[BLOCK_HEADER_LENGTH];
		zs.avail_out = bufferSize - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

		// initialize the zlib compression algorithm
		if(deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY) != Z_OK) {
			printf("ERROR: zlib deflate initialization failed.\n");
			exit(1);
		}

		// compress the data
		int status = deflate(&zs, Z_FINISH);
		if(status != Z_STREAM_END) {
			deflateEnd(&zs);

			// reduce the input length and try again
			if(status == Z_OK) {
				inputLength -= 1024;
				if(inputLength < 0) {
					printf("ERROR: input reduction failed.\n");
					exit(1);
				}
				continue;
			}

			printf("ERROR: zlib deflate failed.\n");
			exit(1);
		}

		// finalize the compression routine
		if(deflateEnd(&zs) != Z_OK) {
			printf("ERROR: deflate end failed.\n");
			exit(1);
		}

		compressedLength = zs.total_out;
		compressedLength += BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;

		if(compressedLength > MAX_BLOCK_SIZE) {
			printf("ERROR: deflate overflow.\n");
			exit(1);
		}

		break;
	}

	// store the compressed length
	BgzfPackUnsignedShort(&buffer[16], (unsigned short)(compressedLength - 1));

	// store the CRC32 checksum
	unsigned int crc = crc32(0, NULL, 0);
	crc = crc32(crc, (Bytef*)mBGZF.UncompressedBlock, inputLength);
	BgzfPackUnsignedInt(&buffer[compressedLength - 8], crc);
	BgzfPackUnsignedInt(&buffer[compressedLength - 4], inputLength);

	// ensure that we have less than a block of data left
	int remaining = mBGZF.BlockOffset - inputLength;
	if(remaining > 0) {
		if(remaining > inputLength) {
			printf("ERROR: remainder too large.\n");
			exit(1);
		}

		memcpy(mBGZF.UncompressedBlock, mBGZF.UncompressedBlock + inputLength, remaining);
	}

	mBGZF.BlockOffset = remaining;
	return compressedLength;
}

// flushes the data in the BGZF block
void CBamWriter::BgzfFlushBlock(void) {

	// flush all of the remaining blocks
	while(mBGZF.BlockOffset > 0) {

		// compress the data block
		int blockLength = BgzfDeflateBlock();

		// flush the data to our output stream
		int numBytesWritten = fwrite(mBGZF.CompressedBlock, 1, blockLength, mBGZF.Stream);

		if(numBytesWritten != blockLength) {
			printf("ERROR: Expected to write %u bytes during flushing, but wrote %u bytes.\n", blockLength, numBytesWritten);
			exit(1);
		}

		mBGZF.BlockAddress += blockLength;
	}
}

// opens the BAM file for writing
void CBamWriter::BgzfOpen(const string& filename) {

	mBGZF.Stream = fopen(filename.c_str(), "wb");

	if(!mBGZF.Stream) {
		printf("ERROR: Unable to open the BAM file (%s) for writing.\n", filename.c_str());
		exit(1);
	}

	mBGZF.IsOpen = true;
}

// writes the supplied data into the BGZF buffer
unsigned int CBamWriter::BgzfWrite(const char* data, const unsigned int dataLen) {

	// initialize
	unsigned int numBytesWritten = 0;
	const char* input = data;
	unsigned int blockLength = mBGZF.UncompressedBlockSize;

	// copy the data to the buffer
	while(numBytesWritten < dataLen) {
		unsigned int copyLength = min(blockLength - mBGZF.BlockOffset, dataLen - numBytesWritten);
		char* buffer = mBGZF.UncompressedBlock;
		memcpy(buffer + mBGZF.BlockOffset, input, copyLength);

		mBGZF.BlockOffset += copyLength;
		input             += copyLength;
		numBytesWritten   += copyLength;

		if(mBGZF.BlockOffset == blockLength) BgzfFlushBlock();
	}

	return numBytesWritten;
}

// closes the alignment archive
void CBamWriter::Close(void) {
	if(mBGZF.IsOpen) BgzfClose();
}

// creates a cigar string from the supplied alignment
void CBamWriter::CreatePackedCigar(const Alignment& al, string& packedCigar, unsigned int& numCigarOperations) {

	// initialize
	const char* pReference = al.Reference.CData();
	const char* pQuery     = al.Query.CData();

	const unsigned int numBases = al.Reference.Length();
	unsigned int currentPos = 0;

	packedCigar.resize(numBases * SIZEOF_INT);
	unsigned int* pPackedCigar = (unsigned int*)packedCigar.data();
	numCigarOperations = 0;

	// create the cigar string by parsing the reference and query strings 
	while(currentPos < numBases) {

		unsigned int testPos         = currentPos;
		unsigned int operationLength = 0;

		if((pReference[currentPos] != '-') && (pQuery[currentPos] != '-')) {

			// find the matches or mismatches
			while((pReference[testPos] != '-') && (pQuery[testPos] != '-') && (testPos < numBases)) {
				++testPos;					
				++operationLength;
			}

			*pPackedCigar = operationLength << BAM_CIGAR_SHIFT | BAM_CMATCH;

		} else if(pReference[currentPos] == '-') {

			// find the insertions
			while((pReference[testPos] == '-') && (testPos < numBases)) {
				++testPos;					
				++operationLength;
			}

			*pPackedCigar = operationLength << BAM_CIGAR_SHIFT | BAM_CINS;

		} else if(pQuery[currentPos] == '-') {

			// find the deletions
			while((pQuery[testPos] == '-') && (testPos < numBases)) {
				++testPos;					
				++operationLength;
			}

			*pPackedCigar = operationLength << BAM_CIGAR_SHIFT | BAM_CDEL;

		} else {
			printf("ERROR: CIGAR string generation failed.\n");
			exit(1);
		}

		// increment our position
		++pPackedCigar;
		++numCigarOperations;
		currentPos += operationLength;

		// make sure aren't creating a buffer overflow
		if(numCigarOperations == numBases) {
			printf("ERROR: buffer overflow detected when creating the packed cigar string.\n");
			exit(1);
		}
	}

	// resize the packed cigar string
	packedCigar.resize(numCigarOperations * SIZEOF_INT);
}

// encodes the supplied query sequence into 4-bit notation
void CBamWriter::EncodeQuerySequence(const CMosaikString& query, string& encodedQuery) {

	// prepare the encoded query string
	const unsigned int queryLen = query.Length();
	const unsigned int encodedQueryLen = (unsigned int)((queryLen / 2.0) + 0.5);
	encodedQuery.resize(encodedQueryLen);
	char* pEncodedQuery = (char*)encodedQuery.data();
	const char* pQuery = (const char*)query.CData();

	unsigned char nucleotideCode;
	bool useHighWord = true;

	while(*pQuery) {

		switch(*pQuery) {
			case '=':
				nucleotideCode = 0;
				break;
			case 'A':
				nucleotideCode = 1;
				break;
			case 'C':
				nucleotideCode = 2;
				break;
			case 'G':
				nucleotideCode = 4;
				break;
			case 'T':
				nucleotideCode = 8;
				break;
			case 'N':
				nucleotideCode = 15;
				break;
			default:
				printf("ERROR: Only the following bases are supported in the BAM format: {=, A, C, G, T, N}. Found [%c]\n", *pQuery);
				exit(1);
		}

		// pack the nucleotide code
		if(useHighWord) {
			*pEncodedQuery = nucleotideCode << 4;
			useHighWord = false;
		} else {
			*pEncodedQuery |= nucleotideCode;
			++pEncodedQuery;
			useHighWord = true;
		}

		// increment the query position
		++pQuery;
	}
}

// opens the alignment archive
void CBamWriter::Open(const string& filename, const BamHeader& header) {

	// open the BGZF file for writing
	BgzfOpen(filename);

	// ====================
	// write the SAM header
	// ====================

	// build header tag
	ostringstream sb;
	sb << "@HD\tVN:" << header.Version << "\tSO:";

	switch(header.SortOrder) {
		case SORTORDER_POSITION:
			sb << "coordinate" << endl;
			break;
		case SORTORDER_READNAME:
			sb << "queryname" << endl;
			break;
		default:
			sb << "unsorted" << endl;
	}

	// build the sequence dictionary
	const unsigned int numReferenceSequences = header.pReferenceSequences->size();
	vector<ReferenceSequence>::const_iterator rsIter;
	for(rsIter = header.pReferenceSequences->begin(); rsIter != header.pReferenceSequences->end(); ++rsIter) {
		sb << "@SQ\tSN:" << rsIter->Name << "\tLN:" << rsIter->NumBases;
		if(!rsIter->GenomeAssemblyID.empty()) sb << "\tAS:" << rsIter->GenomeAssemblyID;
		if(!rsIter->MD5.empty())              sb << "\tM5:" << rsIter->MD5;
		if(!rsIter->URI.empty())              sb << "\tUR:" << rsIter->URI;
		if(!rsIter->Species.empty())          sb << "\tSP:" << rsIter->Species;
		sb << endl;
	}

	// build the read groups
	vector<MosaikReadFormat::ReadGroup>::const_iterator rgIter;
	for(rgIter = header.pReadGroups->begin(); rgIter != header.pReadGroups->end(); ++rgIter) {
		sb << "@RG\tID:" << rgIter->ReadGroupID << "\tSM:" << rgIter->SampleName;
		if(!rgIter->LibraryName.empty())      sb << "\tLB:" << rgIter->LibraryName;
		if(!rgIter->Description.empty())      sb << "\tDS:" << rgIter->Description;
		if(!rgIter->PlatformUnit.empty())     sb << "\tPU:" << rgIter->PlatformUnit;
		if(rgIter->MedianFragmentLength != 0) sb << "\tPI:" << rgIter->MedianFragmentLength;
		if(!rgIter->CenterName.empty())       sb << "\tCN:" << rgIter->CenterName;

		switch(rgIter->SequencingTechnology) {
			case ST_454:
				sb << "\tPL:454" << endl;
				break;
			case ST_HELICOS:
				sb << "\tPL:helicos" << endl;
				break;
			case ST_ILLUMINA:
				sb << "\tPL:illumina" << endl;
				break;
			case ST_PACIFIC_BIOSCIENCES:
				sb << "\tPL:pacific biosciences" << endl;
				break;
			case ST_SOLID:
				sb << "\tPL:solid" << endl;
				break;
			case ST_SANGER:
				sb << "\tPL:sanger" << endl;
				break;
			default:
				sb << "\tPL:unknown" << endl;
		}
	}

	// distill the header text
	string samHeader = sb.str();
	sb.str("");

	// ================
	// write the header
	// ================

	// write the BAM signature
	const unsigned char SIGNATURE_LENGTH = 4;
	const char* BAM_SIGNATURE = "BAM\1";
	BgzfWrite(BAM_SIGNATURE, SIGNATURE_LENGTH);

	// write the SAM header text length
	const unsigned int samHeaderLen = samHeader.size();
	BgzfWrite((char*)&samHeaderLen, SIZEOF_INT);

	//printf("samHeaderLen: %u\n%s\n", samHeaderLen, samHeader.c_str());

	// write the SAM header text
	if(samHeaderLen > 0) BgzfWrite(samHeader.data(), samHeaderLen);

	// write the number of reference sequences
	BgzfWrite((char*)&numReferenceSequences, SIZEOF_INT);

	// =============================
	// write the sequence dictionary
	// =============================

	for(rsIter = header.pReferenceSequences->begin(); rsIter != header.pReferenceSequences->end(); ++rsIter) {

		// write the reference sequence name length
		const unsigned int referenceSequenceNameLen = rsIter->Name.size() + 1;
		BgzfWrite((char*)&referenceSequenceNameLen, SIZEOF_INT);

		// write the reference sequence name
		BgzfWrite(rsIter->Name.c_str(), referenceSequenceNameLen);

		// write the reference sequence length
		BgzfWrite((char*)&rsIter->NumBases, SIZEOF_INT);
	}
}

// saves the alignment to the alignment archive
void CBamWriter::SaveAlignment(const CMosaikString& readName, const string& readGroupID, const vector<Alignment>::iterator& alIter) {

	// =================
	// set the BAM flags
	// =================

	// define our flags
	unsigned int flag                = 0;
	unsigned int queryPosition5Prime = 0;
	unsigned int matePosition5Prime  = 0;
	int insertSize                   = 0;

	if(alIter->IsPairedEnd) {

		flag |= BAM_SEQUENCED_AS_PAIRS;

		// first or second mate?
		flag |= (alIter->IsFirstMate ? BAM_QUERY_FIRST_MATE : BAM_QUERY_SECOND_MATE);

		if(alIter->IsResolvedAsPair) {

			flag |= BAM_PROPER_PAIR;
			if(alIter->IsMateReverseStrand) flag |= BAM_MATE_REVERSE_COMPLEMENT;

			// sanity check
			if(alIter->ReferenceIndex != alIter->MateReferenceIndex) {
				printf("ERROR: The resolved paired-end reads occur on different reference sequences.\n");
				exit(1);
			}

			// set the 5' coordinates
			queryPosition5Prime = (alIter->IsReverseStrand     ? alIter->ReferenceEnd     : alIter->ReferenceBegin);
			matePosition5Prime  = (alIter->IsMateReverseStrand ? alIter->MateReferenceEnd : alIter->MateReferenceBegin);

			// calculate the insert size
			insertSize = matePosition5Prime - queryPosition5Prime;

		} else flag |= BAM_MATE_UNMAPPED;
	}

	if(alIter->IsReverseStrand) flag |= BAM_QUERY_REVERSE_COMPLEMENT;

	// ==========================
	// construct the cigar string
	// ==========================

	string packedCigar;
	unsigned int numCigarOperations;
	CreatePackedCigar(*alIter, packedCigar, numCigarOperations);
	const unsigned int packedCigarLen = packedCigar.size();

	// ===================
	// write the alignment
	// ===================

	// remove the gaps from the read
	CMosaikString query(alIter->Query);
	query.Remove('-');

	// initialize
	const unsigned int nameLen  = readName.Length() + 1;
	const unsigned int queryLen = query.Length();

	// encode the query
	string encodedQuery;
	EncodeQuerySequence(query, encodedQuery);
	const unsigned int encodedQueryLen = encodedQuery.size();

	// create our read group tag
	string readGroupTag;
	const unsigned int readGroupTagLen = 3 + readGroupID.size() + 1;
	readGroupTag.resize(readGroupTagLen);
	char* pReadGroupTag = (char*)readGroupTag.data();
	sprintf(pReadGroupTag, "RGZ%s", readGroupID.c_str());

	// create our mismatch tag
	string mismatchTag = "NMi";
	mismatchTag.resize(MISMATCH_TAG_LEN);
	const unsigned int numMismatches = alIter->NumMismatches;
	memcpy((char*)mismatchTag.data() + 3, (char*)&numMismatches, SIZEOF_INT);

	// retrieve our bin
	unsigned int bin = CalculateMinimumBin(alIter->ReferenceBegin, alIter->ReferenceEnd);

	// assign the BAM core data
	unsigned int buffer[8];
	buffer[0] = alIter->ReferenceIndex;
	buffer[1] = alIter->ReferenceBegin;
	buffer[2] = (bin << 16) | (alIter->Quality << 8) | nameLen;
	buffer[3] = (flag << 16) | numCigarOperations;
	buffer[4] = queryLen;

	if(alIter->IsResolvedAsPair) {
		buffer[5] = alIter->MateReferenceIndex;
		buffer[6] = alIter->MateReferenceBegin;
		buffer[7] = insertSize;
	} else {
		buffer[5] = 0xffffffff;
		buffer[6] = 0xffffffff;
		buffer[7] = 0;
	}

	// write the block size
	const unsigned int dataBlockSize = nameLen + packedCigarLen + encodedQueryLen + queryLen + readGroupTagLen + MISMATCH_TAG_LEN;
	const unsigned int blockSize = BAM_CORE_SIZE + dataBlockSize;
	BgzfWrite((char*)&blockSize, SIZEOF_INT);

	// write the BAM core
	BgzfWrite((char*)&buffer, BAM_CORE_SIZE);

	// write the query name
	BgzfWrite(readName.CData(), nameLen);

	// write the packed cigar
	BgzfWrite(packedCigar.data(), packedCigarLen);

	// write the encoded query sequence
	BgzfWrite(encodedQuery.data(), encodedQueryLen);

	// write the base qualities
	BgzfWrite(alIter->BaseQualities.CData(), queryLen);

	// write the read group tag
	BgzfWrite(readGroupTag.data(), readGroupTagLen);

	// write the mismatch tag
	BgzfWrite(mismatchTag.data(), MISMATCH_TAG_LEN);
}
