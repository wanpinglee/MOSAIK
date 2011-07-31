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
//CBamWriter::CBamWriter(void)
//{}

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

/*
void CBamWriter::TranslateCigarToPackCigar ( const string cigar, string packCigar ) {

	packCigar.resize( cigar.size() * SIZEOF_INT );
	for ( unsigned int i = 0; i < cigar.size(); ++i ) {
		switch( cigar[i] ) {
		}
	}
}
*/


// creates a cigar string from the supplied alignment
void CBamWriter::CreatePackedCigar( const Alignment& al, string& packedCigar, unsigned short& numCigarOperations, const bool isSolid ) {

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

		if((pReference[currentPos] != '-') && (pQuery[currentPos] != '-') && (pReference[currentPos] != 'Z')) {

			// find the matches or mismatches
			while((pReference[testPos] != '-') && (pQuery[testPos] != '-') && (testPos < numBases) && (pReference[testPos] != 'Z') && (testPos < numBases)) {
				++testPos;					
				++operationLength;
			}

			*pPackedCigar = operationLength << BAM_CIGAR_SHIFT | BAM_CMATCH;

		} else if ( pReference[currentPos] == 'Z' ) {
			while( ( pReference[testPos] == 'Z' ) && ( testPos < numBases ) ){
				++testPos;
				++operationLength;
			}

			*pPackedCigar = operationLength << BAM_CIGAR_SHIFT | BAM_CSOFT_CLIP;

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
			case ST_ILLUMINA_LONG:
				sb << "\tPL:illumina long" << endl;
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

	if ( !header.pg.ID.empty() ) {
		sb << "@PG\tID:" << header.pg.ID;
		if ( !header.pg.VN.empty() )
			sb << "\tVN:" << header.pg.VN;
		if ( !header.pg.CL.empty() )
			sb << "\tCL:" << header.pg.CL;
		sb << endl;
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
// saves the reference and position of an alignment to the alignment archive
void CBamWriter::SaveReferencePosition( const unsigned int refIndex, const unsigned int refBegin, const unsigned int refEnd ) {
	// =================
	// set the BAM flags
	// =================
	
	unsigned int flag                = 0;

	// retrieve our bin
	unsigned int bin = CalculateMinimumBin(refBegin, refEnd);
	const unsigned int quality = 0;
	const unsigned int nameLen = 0;
	const unsigned int numCigarOperations = 0;
	const unsigned int packedCigarLen = 0;
	const unsigned int queryLen = 0;
	const unsigned int encodedQueryLen = 0;

	// assign the BAM core data
	unsigned int buffer[8];
	buffer[0] = refIndex;
	buffer[1] = refBegin;
	buffer[2] = (bin << 16) | (quality << 8) | nameLen;
	buffer[3] = (flag << 16) | numCigarOperations;
	buffer[4] = queryLen;    // read_len
	buffer[5] = 0xffffffff;  // mate_rID
	buffer[6] = 0xffffffff;  // mate_pos
	buffer[7] = 0;           // ins_size

	const char* startChar = '\0';
	
	// write the block size
	const unsigned int dataBlockSize = nameLen + packedCigarLen + encodedQueryLen + queryLen;
	//const unsigned int dataBlockSize = nameLen;
	const unsigned int blockSize = BAM_CORE_SIZE + dataBlockSize;
	BgzfWrite((char*)&blockSize, SIZEOF_INT);

	// write the BAM core
	BgzfWrite((char*)&buffer, BAM_CORE_SIZE);

	// write the query name
	BgzfWrite(startChar, nameLen);

	// write the packed cigar
	BgzfWrite(startChar, packedCigarLen);

	// write the encoded query sequence
	BgzfWrite(startChar, encodedQueryLen);

	// write the base qualities
	BgzfWrite(startChar, queryLen);

}


// saves the alignment to the alignment archive
void CBamWriter::SaveAlignment(const Alignment al, const char* zaString, const bool& noCigarMdNm, const bool& notShowRnamePos, const bool& isSolid, const bool processedBamData ) {

	// =================
	// set the BAM flags
	// =================

	// define our flags
	unsigned int flag                = 0;
	//unsigned int queryPosition5Prime = 0;
	//unsigned int matePosition5Prime  = 0;
	int insertSize                   = 0;

	if(al.IsPairedEnd) {

		flag |= BAM_SEQUENCED_AS_PAIRS;

		// first or second mate?
		flag |= (al.IsFirstMate ? BAM_QUERY_FIRST_MATE : BAM_QUERY_SECOND_MATE);

		if(al.IsResolvedAsPair) {

			if ( al.IsResolvedAsProperPair ) 
				flag |= BAM_PROPER_PAIR;
				
			if(al.IsMateReverseStrand) flag |= BAM_MATE_REVERSE_COMPLEMENT;

			// sanity check
			//if(alIter->ReferenceIndex != alIter->MateReferenceIndex) {
			//	printf("ERROR: The resolved paired-end reads occur on different reference sequences.\n");
			//	exit(1);
			//}

			// set the 5' coordinates
			//queryPosition5Prime = (alIter->IsReverseStrand     ? alIter->ReferenceEnd     : alIter->ReferenceBegin);
			//matePosition5Prime  = (alIter->IsMateReverseStrand ? alIter->MateReferenceEnd : alIter->MateReferenceBegin);

			// calculate the insert size
			//insertSize = matePosition5Prime - queryPosition5Prime;
			insertSize = al.FragmentLength;

		}

		if ( !al.IsMapped )
			flag |= BAM_QUERY_UNMAPPED;
		if ( !al.IsMateMapped )
			flag |= BAM_MATE_UNMAPPED;
			
	}

	if(al.IsReverseStrand) flag |= BAM_QUERY_REVERSE_COMPLEMENT;

	// ==========================
	// construct the cigar string
	// ==========================

	string packedCigar;
	unsigned short numCigarOperations = 0;
	if ( !noCigarMdNm ) {
		if ( !processedBamData )
			CreatePackedCigar( al, packedCigar, numCigarOperations, isSolid );
		else {
			packedCigar = al.PackedCigar;
			numCigarOperations = al.NumCigarOperation;
		}
	}
	else
		packedCigar = "\0";
	const unsigned int packedCigarLen = !noCigarMdNm ? packedCigar.size() : 0;

	// ===================
	// write the alignment
	// ===================

	// remove the gaps from the read
	CMosaikString query;
	if ( !processedBamData ) {
		query = al.Query.CData();
		query.Remove('-');
	}

	// initialize
	const unsigned int nameLen  = al.Name.Length() + 1;
	const unsigned int queryLen = processedBamData ? al.QueryLength : query.Length();

	// sanity check
	//al.BaseQualities.CheckQuality();
	//if ( queryLen != alIter->BaseQualities.Length() ) {
	//        printf("ERROR: The lengths of bases(%u) and qualities(%u) of Read (%s) didn't match.\n", queryLen, alIter->BaseQualities.Length(), readName.CData());
        //        exit(1);
        //}
	
	// encode the query
	string encodedQuery;
	if ( !processedBamData )
		EncodeQuerySequence(query, encodedQuery);
	else
		encodedQuery = al.EncodedQuery;
	const unsigned int encodedQueryLen = encodedQuery.size();

	// create our read group tag
	string readGroupTag;
	const unsigned int readGroupTagLen = 3 + al.ReadGroup.size() + 1;
	readGroupTag.resize(readGroupTagLen);
	char* pReadGroupTag = (char*)readGroupTag.data();
	sprintf(pReadGroupTag, "RGZ%s", al.ReadGroup.c_str());

	// create our mismatch tag
	string mismatchTag;
	unsigned int numMismatches = 0;
	unsigned int nmTagLen = 0;
	if ( !noCigarMdNm ) {
		mismatchTag = "NMi";
		mismatchTag.resize(MISMATCH_TAG_LEN);
		nmTagLen = MISMATCH_TAG_LEN;
		numMismatches = al.NumMismatches;
		memcpy((char*)mismatchTag.data() + 3, (char*)&numMismatches, SIZEOF_INT);
	}

	// create our MD tag
	string mdTag;
	char* pMd = 0;
	unsigned int mdTagLen = 0;
	char* pMdTag;
	if ( !noCigarMdNm ) {
		if ( !processedBamData ) 
			pMd = (char*) mdTager.GetMdTag( al.Reference.CData(), al.Query.CData(), al.Reference.Length() );
		else
			pMd = (char*) al.MdString.c_str();

		mdTagLen = 3 + strlen( pMd ) + 1;
		mdTag.resize( mdTagLen );
		pMdTag = (char*)mdTag.data();
		sprintf(pMdTag, "MDZ%s", pMd);
	}

	// create our za tag
	unsigned int zaTagLen = 0;
	string zaTag;
	char* pZaTag;
	if ( zaString != 0 ) {
		zaTagLen = 3 + strlen( zaString ) + 1;
		zaTag.resize( zaTagLen );
		pZaTag = (char*)zaTag.data();
		sprintf(pZaTag, "ZAZ%s",zaString);
	}

	// create our cs tag
	unsigned int csTagLen = 0;
	string csTag;
	char* pCsTag;
	if ( isSolid ) {
		csTagLen = 3 + strlen ( al.CsQuery.c_str() ) + 1;
		csTag.resize( csTagLen );
		pCsTag = (char*)csTag.data();
		sprintf( pCsTag, "CSZ%s", al.CsQuery.c_str() );
	}

	// create our cq tag
	unsigned int cqTagLen = 0;
	string cqTag;
	char* pCqTag;
	if ( isSolid ) {
		cqTagLen = 3 + strlen ( al.CsBaseQualities.c_str() ) + 1;
		cqTag.resize( cqTagLen );
		pCqTag = (char*)cqTag.data();
		sprintf( pCqTag, "CQZ%s", al.CsBaseQualities.c_str() );
	}

	// retrieve our bin
	unsigned int bin = CalculateMinimumBin(al.ReferenceBegin, al.ReferenceEnd);

	// assign the BAM core data
	unsigned int buffer[8] = {0};
	buffer[0] = notShowRnamePos ? 0xffffffff : al.ReferenceIndex;
	buffer[1] = notShowRnamePos ? 0xffffffff : al.ReferenceBegin;
	buffer[2] = (bin << 16) | (al.RecalibratedQuality << 8) | nameLen;
	buffer[3] = (flag << 16) | numCigarOperations;
	buffer[4] = queryLen;

	if(al.IsResolvedAsPair) {
		buffer[5] = al.IsMateMapped ? al.MateReferenceIndex : 0xffffffff;
		buffer[6] = al.IsMateMapped ? al.MateReferenceBegin : 0xffffffff;
		buffer[7] = insertSize;
	} else {
		buffer[5] = 0xffffffff;
		buffer[6] = 0xffffffff;
		buffer[7] = 0;
	}

	// write the block size
	const unsigned int dataBlockSize = nameLen + packedCigarLen + encodedQueryLen + queryLen + readGroupTagLen + nmTagLen + mdTagLen + zaTagLen + csTagLen + cqTagLen;
	const unsigned int blockSize = BAM_CORE_SIZE + dataBlockSize;
	BgzfWrite((char*)&blockSize, SIZEOF_INT);

	// write the BAM core
	BgzfWrite((char*)&buffer, BAM_CORE_SIZE);

	// write the query name
	BgzfWrite(al.Name.CData(), nameLen);

	// write the packed cigar
	BgzfWrite(packedCigar.data(), packedCigarLen);

	// write the encoded query sequence
	BgzfWrite(encodedQuery.data(), encodedQueryLen);

	// write the base qualities
	BgzfWrite(al.BaseQualities.CData(), queryLen);

	// write the read group tag
	BgzfWrite(readGroupTag.data(), readGroupTagLen);

	// write the mismatch tag
	if ( !noCigarMdNm )
		BgzfWrite(mismatchTag.data(), MISMATCH_TAG_LEN);

	// write the MD tag
	if ( !noCigarMdNm )
	BgzfWrite(mdTag.data(), mdTagLen);

	// write the ZA tag
	if ( zaString != 0 )
		BgzfWrite(zaTag.data(), zaTagLen);

	// write the cs tag
	if ( isSolid )
		BgzfWrite(csTag.data(), csTagLen);

	// write the cq tag
	if ( isSolid )
		BgzfWrite(cqTag.data(), cqTagLen);
}
