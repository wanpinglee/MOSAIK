// ***************************************************************************
// CGigaBayesFormat - exports MOSAIK assemblies into the GigaBayes (GIG)
//                    format.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "GigaBayesFormat.h"

// constructor
CGigaBayesFormat::CGigaBayesFormat(unsigned char referenceBaseQuality) 
: mBuffer(NULL)
, mBufferLen(0)
, mBQBuffer(NULL)
, mBQBufferLen(0)
, mReferenceSequenceBaseQuality(referenceBaseQuality)
{
	// set the file open flag
	mIsOpen = false;
}

// destructor
CGigaBayesFormat::~CGigaBayesFormat(void) {
	if(mIsOpen)   Close();
	if(mBuffer)   delete [] mBuffer;
	if(mBQBuffer) delete [] mBQBuffer;
}

// checks if the base quality buffer is large enough to accomodate the requested size
void CGigaBayesFormat::CheckBQBufferSize(const unsigned int requestedLength) {
	try {
		if(requestedLength > mBQBufferLen) {
			mBQBufferLen = requestedLength + 10;
			delete [] mBQBuffer;
			mBQBuffer = new unsigned short[mBQBufferLen];
		}
	} catch(bad_alloc) {
		cout << "ERROR: Out of memory when allocating " << requestedLength << " unsigned shorts." << endl;
		exit(1);
	}
}

// opens the GIG file
void CGigaBayesFormat::Open(const CMosaikString& filename) {

	// open the main file stream
	if(fopen_s(&mOutputStream, filename.CData(), "wb") != 0) {
		printf("ERROR: Unable to open the GIG file stream for writing.\n");
		exit(1);
	}

	// open the temporary read index
	CFileUtilities::GetTempFilename(mTempReadIndexFilename);
	if(fopen_s(&mReadIndexStream, mTempReadIndexFilename.c_str(), "w+b") != 0) {
		cout << "ERROR: Unable to open the GIG read index file stream for writing." << endl;
		exit(1);
	}

	mIsOpen = true;
	mGappedRefLength  = 0;

	mContigIndexOffset = 0;
}

// closes the GIG file
void CGigaBayesFormat::Close(void) {

	mIsOpen = false;

	// close the GIG file
	//cout << endl << "- joining temporary files... ";
	//cout.flush();

	// retrieve the read index file offset
	mReadIndexOffset = ftell64(mOutputStream);

	// create a large buffer
	const unsigned int requestedBytes = 52428800; // 50 MB
	CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, requestedBytes);

	// append to the output stream
	printf("\n- appending read index to read data... ");
	fflush(stdout);

	rewind(mReadIndexStream);

	size_t numBytesRead = 0, numBytesWritten = 0;
	while(true) {
		numBytesRead = fread(mBuffer, 1, requestedBytes, mReadIndexStream);
		//printf("bytes read: %u\n", (unsigned int)numBytesRead);
		if(feof(mReadIndexStream)) break;
		numBytesWritten = fwrite(mBuffer, 1, numBytesRead, mOutputStream);
		//printf("bytes written: %u\n", (unsigned int)numBytesWritten);
	}

	numBytesWritten = fwrite(mBuffer, 1, numBytesRead, mOutputStream);
	//printf("bytes written: %u\n", (unsigned int)numBytesWritten);

	printf("finished.\n");

	// ======================
	// write the contig index
	// ======================

	mContigIndexOffset = ftell64(mOutputStream);

	const unsigned int CONTIG_NAME_LENGTH = mContigName.size();
	const unsigned int requestedBytesContigIndex = SIZEOF_INT + CONTIG_NAME_LENGTH + SIZEOF_OFF_TYPE + SIZEOF_CHAR;
	CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, requestedBytesContigIndex);	
	unsigned int bufferOffset = 0;

	// store the name length
	memcpy(mBuffer + bufferOffset, (char*)&CONTIG_NAME_LENGTH, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the read name
	memcpy(mBuffer + bufferOffset, mContigName.c_str(), CONTIG_NAME_LENGTH);
	bufferOffset += CONTIG_NAME_LENGTH;
	mBuffer[bufferOffset++] = 0;

	// store the read offset
	//printf("contig file offset: %llu\n", (unsigned long long)mContigFileOffset);
	memcpy(mBuffer + bufferOffset, (char*)&mContigFileOffset, SIZEOF_OFF_TYPE);
	bufferOffset += SIZEOF_OFF_TYPE;

	//printf("buffer offset: %u, requested bytes: %u\n", bufferOffset, requestedBytesContigIndex);

	// sanity check
	if(bufferOffset > mBufferLen) {
		printf("ERROR: Buffer overflow detected when writing the reference sequence into the read index.\n");
		exit(1);
	}

	// write the contig index
	fwrite(mBuffer, bufferOffset, 1, mOutputStream);

	// adjust the index positions in the header
	fseek64(mOutputStream, mContigIndexOffsetLocus, SEEK_SET);
	fwrite((char*)&mContigIndexOffset, SIZEOF_OFF_TYPE, 1, mOutputStream);

	fseek64(mOutputStream, mReadIndexOffsetLocus, SEEK_SET);
	fwrite((char*)&mReadIndexOffset, SIZEOF_OFF_TYPE, 1, mOutputStream);

	// close our files
	fclose(mOutputStream);
	fclose(mReadIndexStream);

	// delete our temporary files
	rm(mTempReadIndexFilename.c_str());

	//cout << "finished." << endl;
}

// saves the specified reference sequence to the header
void CGigaBayesFormat::SaveHeader(CMosaikString& reference, const string& referenceName, const unsigned int ungappedRefLength, const uint64_t& numSequences) {

	// write the GIG file signature
	const unsigned int ASSEMBLY_NAME_LENGTH  = 8;
	const unsigned int SIGNATURE_LENGTH      = 7;
	const unsigned int REFERENCE_NAME_LENGTH = referenceName.size();

	mGappedRefLength = reference.Length();
	mContigName = referenceName;

	const char* BASE_SEGMENT_NAME = ".MosaikReference";
	const unsigned int BASE_SEGMENT_NAME_LENGTH = strlen(BASE_SEGMENT_NAME);

	const unsigned int requestedBytesHeader      = SIGNATURE_LENGTH + ASSEMBLY_NAME_LENGTH + 2 * SIZEOF_OFF_TYPE + 4 * SIZEOF_INT + 2 * SIZEOF_CHAR;
	const unsigned int requestedBytesContig      = REFERENCE_NAME_LENGTH + mGappedRefLength + 2 * SIZEOF_OFF_TYPE + 4 * SIZEOF_INT + mGappedRefLength * SIZEOF_SHORT + 2 * SIZEOF_CHAR;
	const unsigned int requestedBytesBaseSegment = 3 * SIZEOF_INT + BASE_SEGMENT_NAME_LENGTH + SIZEOF_CHAR;
	const unsigned int requestedBytesReference   = mGappedRefLength + BASE_SEGMENT_NAME_LENGTH + 8 * SIZEOF_INT + 3 * SIZEOF_CHAR + mGappedRefLength * SIZEOF_SHORT;
	const unsigned int requestedBytes            = requestedBytesHeader + requestedBytesContig + requestedBytesBaseSegment + requestedBytesReference;

	CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, requestedBytes);	
	unsigned int bufferOffset = 0;

	// store the signature length
	memcpy(mBuffer + bufferOffset, (char*)&SIGNATURE_LENGTH, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the signature
	memcpy(mBuffer + bufferOffset, "GIG-0.0", SIGNATURE_LENGTH);
	bufferOffset += SIGNATURE_LENGTH;
	mBuffer[bufferOffset++] = 0;

	// store the number of contigs
	const unsigned int numContigs = 1;
	memcpy(mBuffer + bufferOffset, (char*)&numContigs, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the contig index offset
	mContigIndexOffsetLocus = SIZEOF_INT + SIGNATURE_LENGTH + SIZEOF_CHAR + SIZEOF_INT;
	const off_type FILE_OFFSET = 0;
	memcpy(mBuffer + bufferOffset, (char*)&FILE_OFFSET, SIZEOF_OFF_TYPE);
	bufferOffset += SIZEOF_OFF_TYPE;

	// store the contig block offset
	const off_type CONTIG_BLOCK_OFFSET = requestedBytesHeader;
	memcpy(mBuffer + bufferOffset, (char*)&CONTIG_BLOCK_OFFSET, SIZEOF_OFF_TYPE);
	bufferOffset += SIZEOF_OFF_TYPE;

	// store the number of assembled reads
	// TODO: change this to uint64_t - Gabor has this as an int
	int numSequencesInt = (int)(numSequences + 1);
	memcpy(mBuffer + bufferOffset, (char*)&numSequencesInt, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the length of the assembly name
	memcpy(mBuffer + bufferOffset, (char*)&ASSEMBLY_NAME_LENGTH, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the assembly name
	memcpy(mBuffer + bufferOffset, "ASSEMBLY", ASSEMBLY_NAME_LENGTH);
	bufferOffset += ASSEMBLY_NAME_LENGTH;
	mBuffer[bufferOffset++] = 0;

	// ------------------------------------------------------------------

	// store the current location
	mContigFileOffset = requestedBytesHeader;

	// store the number of contig reads
	// TODO: change this to uint64_t - Gabor has this as an int
	memcpy(mBuffer + bufferOffset, (char*)&numSequencesInt, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the read index offset
	// TODO: find the read index offset
	mReadIndexOffsetLocus = requestedBytesHeader + SIZEOF_INT;
	memcpy(mBuffer + bufferOffset, (char*)&FILE_OFFSET, SIZEOF_OFF_TYPE);
	bufferOffset += SIZEOF_OFF_TYPE;

	// store the read block offset
	const off_type READ_BLOCK_OFFSET = requestedBytesHeader + requestedBytesContig + requestedBytesBaseSegment;
	memcpy(mBuffer + bufferOffset, (char*)&READ_BLOCK_OFFSET, SIZEOF_OFF_TYPE);
	bufferOffset += SIZEOF_OFF_TYPE;

	// store the contig name length
	memcpy(mBuffer + bufferOffset, (char*)&REFERENCE_NAME_LENGTH, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the padded contig length
	memcpy(mBuffer + bufferOffset, (char*)&mGappedRefLength, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the number of base segments
	const unsigned int NUM_BASE_SEGMENTS = 1;
	memcpy(mBuffer + bufferOffset, (char*)&NUM_BASE_SEGMENTS, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the contig name
	memcpy(mBuffer + bufferOffset, referenceName.c_str(), REFERENCE_NAME_LENGTH);
	bufferOffset += REFERENCE_NAME_LENGTH;
	mBuffer[bufferOffset++] = 0;

	// store the contig bases
	memcpy(mBuffer + bufferOffset, reference.CData(), mGappedRefLength);
	bufferOffset += mGappedRefLength;
	mBuffer[bufferOffset++] = 0;

	// store the contig base qualities
	// TODO: change this to char - Gabor has this as a short
	const char* pReference = reference.CData();

	const unsigned short ZERO_BQ      = 0;
	const unsigned short REFERENCE_BQ = mReferenceSequenceBaseQuality;

	for(unsigned int i = 0; i < mGappedRefLength; i++, pReference++) {
		if(*pReference == '-') memcpy(mBuffer + bufferOffset, (char*)&ZERO_BQ, SIZEOF_SHORT);
		else memcpy(mBuffer + bufferOffset, (char*)&REFERENCE_BQ, SIZEOF_SHORT);
		bufferOffset += SIZEOF_SHORT;
	}

	// ------------------------------------------------------------------

	// store the base segment name length
	memcpy(mBuffer + bufferOffset, (char*)&BASE_SEGMENT_NAME_LENGTH, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the base segment start position
	const unsigned int BASE_SEGMENT_START = 1;
	memcpy(mBuffer + bufferOffset, (char*)&BASE_SEGMENT_START, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the base segment end position
	memcpy(mBuffer + bufferOffset, (char*)&mGappedRefLength, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the base segment name
	memcpy(mBuffer + bufferOffset, BASE_SEGMENT_NAME, BASE_SEGMENT_NAME_LENGTH);
	bufferOffset += BASE_SEGMENT_NAME_LENGTH;
	mBuffer[bufferOffset++] = 0;

	// =======================================
	// write the reference sequence read entry
	// =======================================

	// retrieve the read offset
	const off_type REFERENCE_FILE_OFFSET = ftell64(mOutputStream);

	// store the name length
	memcpy(mBuffer + bufferOffset, (char*)&BASE_SEGMENT_NAME_LENGTH, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the read length
	memcpy(mBuffer + bufferOffset, (char*)&mGappedRefLength, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the begin coordinate
	memcpy(mBuffer + bufferOffset, (char*)&BASE_SEGMENT_START, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the end coordinate
	memcpy(mBuffer + bufferOffset, (char*)&mGappedRefLength, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the quality trim begin endpoint
	const unsigned int TRIM_BEGIN = 1;
	memcpy(mBuffer + bufferOffset, (char*)&TRIM_BEGIN, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the quality trim end endpoint
	memcpy(mBuffer + bufferOffset, (char*)&mGappedRefLength, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the alignment trim begin endpoint
	memcpy(mBuffer + bufferOffset, (char*)&TRIM_BEGIN, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the alignment trim end endpoint
	memcpy(mBuffer + bufferOffset, (char*)&mGappedRefLength, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the orientation
	const bool IsReverseStrand = false;
	memcpy(mBuffer + bufferOffset, (char*)&IsReverseStrand, SIZEOF_CHAR);
	bufferOffset += SIZEOF_CHAR;

	// store the read name
	memcpy(mBuffer + bufferOffset, BASE_SEGMENT_NAME, BASE_SEGMENT_NAME_LENGTH);
	bufferOffset += BASE_SEGMENT_NAME_LENGTH;
	mBuffer[bufferOffset++] = 0;

	// store the read bases
	memcpy(mBuffer + bufferOffset, reference.CData(), mGappedRefLength);
	bufferOffset += mGappedRefLength;
	mBuffer[bufferOffset++] = 0;

	// store the read qualities
	for(unsigned int i = 0; i < mGappedRefLength; i++, pReference++) {
		memcpy(mBuffer + bufferOffset, (char*)&ZERO_BQ, SIZEOF_SHORT);
		bufferOffset += SIZEOF_SHORT;
	}

	//printf("buffer offset: %u, requested bytes: %u\n", bufferOffset, requestedBytes);

	// sanity check
	if(bufferOffset > mBufferLen) {
		printf("ERROR: Buffer overflow detected when writing the header.\n");
		exit(1);
	}

	// write the header
	fwrite(mBuffer, bufferOffset, 1, mOutputStream);

	// =============================================
	// write the reference sequence read index entry
	// =============================================

	const unsigned int requestedBytesReadIndex = SIZEOF_INT + BASE_SEGMENT_NAME_LENGTH + SIZEOF_OFF_TYPE + SIZEOF_CHAR;
	CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, requestedBytesReadIndex);	
	bufferOffset = 0;

	// store the name length
	memcpy(mBuffer + bufferOffset, (char*)&BASE_SEGMENT_NAME_LENGTH, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the read name
	memcpy(mBuffer + bufferOffset, BASE_SEGMENT_NAME, BASE_SEGMENT_NAME_LENGTH);
	bufferOffset += BASE_SEGMENT_NAME_LENGTH;
	mBuffer[bufferOffset++] = 0;

	// store the read offset
	memcpy(mBuffer + bufferOffset, (char*)&REFERENCE_FILE_OFFSET, SIZEOF_OFF_TYPE);
	bufferOffset += SIZEOF_OFF_TYPE;

	//printf("buffer offset: %u, requested bytes: %u\n", bufferOffset, requestedBytesReadIndex);

	// sanity check
	if(bufferOffset > mBufferLen) {
		printf("ERROR: Buffer overflow detected when writing the reference sequence into the read index.\n");
		exit(1);
	}

	// write the header
	fwrite(mBuffer, bufferOffset, 1, mReadIndexStream);
}

// saves the specified read to the index and reads files
void CGigaBayesFormat::SaveRead(const Alignment& al, CMosaikString& gappedRead) {

	// ====================
	// write the read entry
	// ====================

	const unsigned int NAME_LENGTH = al.Name.Length();
	const unsigned int READ_LENGTH = gappedRead.Length();

	const unsigned int requestedBytes = READ_LENGTH + NAME_LENGTH + 8 * SIZEOF_INT + 3 * SIZEOF_CHAR + READ_LENGTH * SIZEOF_SHORT;
	CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, requestedBytes);
	unsigned int bufferOffset = 0;

	// retrieve the read offset
	//const off_type READ_FILE_OFFSET = ftell64(mOutputStream); // never used

	// store the name length
	memcpy(mBuffer + bufferOffset, (char*)&NAME_LENGTH, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the read length
	memcpy(mBuffer + bufferOffset, (char*)&READ_LENGTH, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the begin coordinate
	const int BEGIN_COORDINATE = (int)mpUngap2Gap[al.ReferenceBegin] + 1;
	memcpy(mBuffer + bufferOffset, (char*)&BEGIN_COORDINATE, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the end coordinate
	const int END_COORDINATE = BEGIN_COORDINATE + READ_LENGTH - 1;
	memcpy(mBuffer + bufferOffset, (char*)&END_COORDINATE, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the quality trim begin endpoint
	const unsigned int TRIM_BEGIN = 1;
	memcpy(mBuffer + bufferOffset, (char*)&TRIM_BEGIN, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the quality trim end endpoint
	memcpy(mBuffer + bufferOffset, (char*)&READ_LENGTH, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the alignment trim begin endpoint
	memcpy(mBuffer + bufferOffset, (char*)&TRIM_BEGIN, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the alignment trim end endpoint
	memcpy(mBuffer + bufferOffset, (char*)&READ_LENGTH, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the orientation
	memcpy(mBuffer + bufferOffset, (char*)&al.IsReverseStrand, SIZEOF_CHAR);
	bufferOffset += SIZEOF_CHAR;

	// store the read name
	memcpy(mBuffer + bufferOffset, al.Name.CData(), NAME_LENGTH);
	bufferOffset += NAME_LENGTH;
	mBuffer[bufferOffset++] = 0;

	// store the read bases
	memcpy(mBuffer + bufferOffset, gappedRead.CData(), READ_LENGTH);
	bufferOffset += READ_LENGTH;
	mBuffer[bufferOffset++] = 0;

	// store the read qualities
	//const unsigned short ZERO_BQ = 0; // never used
	//unsigned int ungappedPos = 0;     // never used

	// reverse the base qualities if required
	// TODO: check with Gabor about this... MOSAIK already has these reversed
	CMosaikString baseQualities = al.BaseQualities;
	if(al.IsReverseStrand) baseQualities.Reverse();

	const char* pBaseQualities = baseQualities.CData();
	const char* pRead          = gappedRead.CData();

	CheckBQBufferSize(READ_LENGTH);

	// copy the base qualities (use previous base quality for INDELs)
	bool containsIndel = false;
	unsigned short previousBaseQuality = 0;
	for(unsigned int i = 0; i < READ_LENGTH; i++) {
		if(pRead[i] == '-') {
			mBQBuffer[i] = previousBaseQuality;
			containsIndel = true;
		} else {
			mBQBuffer[i] = *pBaseQualities;
			previousBaseQuality = *pBaseQualities;
			pBaseQualities++;
		}
	}

	// copy the base qualities (use downstream base quality for INDELs if lower)
	if(containsIndel) {
		unsigned short previousBaseQuality = 0;
		for(int i = READ_LENGTH - 1; i >= 0; i--) {
			if((pRead[i] == '-') && (mBQBuffer[i] > previousBaseQuality)) {
				mBQBuffer[i] = previousBaseQuality;
			} else previousBaseQuality = mBQBuffer[i];
		}
	}

	// store the base qualities to the file stream
	for(unsigned int i = 0; i < READ_LENGTH; i++) {
		memcpy(mBuffer + bufferOffset, (char*)&mBQBuffer[i], SIZEOF_SHORT);
		bufferOffset += SIZEOF_SHORT;
	}

	//printf("buffer offset: %u, requested bytes: %u\n", bufferOffset, requestedBytes);

	// sanity check
	if(bufferOffset > mBufferLen) {
		printf("ERROR: Buffer overflow detected when writing read data.\n");
		exit(1);
	}

	// write the index entry
	fwrite(mBuffer, bufferOffset, 1, mOutputStream);

	// =====================
	// write the index entry
	// =====================

	//const unsigned int requestedBytesReadIndex = SIZEOF_INT + NAME_LENGTH + SIZEOF_OFF_TYPE + SIZEOF_CHAR;
	CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, requestedBytes);
	bufferOffset = 0;

	// store the name length
	memcpy(mBuffer + bufferOffset, (char*)&NAME_LENGTH, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// store the read name
	memcpy(mBuffer + bufferOffset, al.Name.CData(), NAME_LENGTH);
	bufferOffset += NAME_LENGTH;
	mBuffer[bufferOffset++] = 0;

	// store the read offset
	// TODO: Ask Gabor why the read offsets aren't saved
	const off_type ZERO_OFFSET = 0;
	memcpy(mBuffer + bufferOffset, (char*)&ZERO_OFFSET, SIZEOF_OFF_TYPE);
	bufferOffset += SIZEOF_OFF_TYPE;

	//printf("buffer offset: %u, requested bytes: %u\n", bufferOffset, requestedBytesReadIndex);

	// sanity check
	if(bufferOffset > mBufferLen) {
		printf("ERROR: Buffer overflow detected when writing read index data.\n");
		exit(1);
	}

	// write the index entry
	fwrite(mBuffer, bufferOffset, 1, mReadIndexStream);

	//const unsigned int gappedReadLength = gappedRead.Length();
	//fprintf(mReadStream, "RD %s %u 0 0\n", al.Name.CData(), gappedReadLength);
	//gappedRead.Replace('-', '*');
	//ChopSequence(mReadStream, gappedRead);
	//fprintf(mReadStream, "\nQA 1 %u 1 %u\nDS CHROMAT_FILE: %s.scf PHD_FILE: %s.scf.phd.1 TIME: %s\n\n", gappedReadLength, gappedReadLength, al.Name.CData(), al.Name.CData(), mTimeString.c_str());
}
