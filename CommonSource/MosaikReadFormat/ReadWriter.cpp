// ***************************************************************************
// CReadWriter - stores reads in a MOSAIK read archive.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "ReadWriter.h"

namespace MosaikReadFormat {

	// constructor
	CReadWriter::CReadWriter(void)
		: mIsOpen(false)
		, mOutStream(NULL)
		, mNumReads(0)
		, mNumBases(0)
		, mBuffer(NULL)
		, mBufferLen(10485760)
		, mBufferPosition(0)
		, mCompressionBuffer(NULL)
		, mCompressionBufferLen(0)
		, mPartitionSize(20000)
		, mPartitionMembers(0)
		, mIsSOLiD(false)
	{
		// set the buffer threshold
		mBufferThreshold = mBufferLen - MEMORY_BUFFER_SIZE;

		// initialize the read and index buffer
		try {
			mBuffer = new unsigned char[mBufferLen];
		} catch(bad_alloc) {
			cout << "ERROR: Unable to allocate enough memory for the read archive output buffer." << endl;
			exit(1);
		}
	}

	// destructor
	CReadWriter::~CReadWriter(void) {
		if(mIsOpen)            Close();
		if(mBuffer)            delete [] mBuffer;
		if(mCompressionBuffer) delete [] mCompressionBuffer;
	}

	// checks the buffer
	void CReadWriter::AdjustBuffer(void) {

		// allocate a new buffer
		unsigned int newBufferLen = mBufferLen << 1;
		unsigned char* newBuffer = NULL;

		// DEBUG
		//cout << "changing buffer size to: " << newBufferLen << " bytes." << endl;

		try {
			newBuffer = new unsigned char[newBufferLen];
		} catch(bad_alloc) {
			cout << "ERROR: Unable to reallocate enough memory for the read archive output buffer." << endl;
			exit(1);
		}

		// copy the old data and destroy the old buffer
		memcpy(newBuffer, mBuffer, mBufferLen);
		delete [] mBuffer;

		// repoint the new buffer
		mBuffer          = newBuffer;
		mBufferLen       = newBufferLen;
		mBufferThreshold = newBufferLen - MEMORY_BUFFER_SIZE;
	}

	// closes the read archive
	void CReadWriter::Close(void) {

		// prevent the archive from being updated elsewhere
		mIsOpen = false;

		// flush the buffer
		if(mPartitionMembers > 0) WritePartition();

		// =================
		// update the header
		// =================

		// update the number of reads in the archive
		fseek64(mOutStream, UPDATE_HEADER_OFFSET, SEEK_SET);
		fwrite((char*)&mNumReads, SIZEOF_UINT64, 1, mOutStream);

		// update the number of bases in the archive
		fwrite((char*)&mNumBases, SIZEOF_UINT64, 1, mOutStream);

		// close the file stream
		fclose(mOutStream);
	}

	// retrieves the number of bases written
	uint64_t CReadWriter::GetNumBases(void) const {
		return mNumBases;
	}

	// retrieves the number of reads written
	uint64_t CReadWriter::GetNumReads(void) const {
		return mNumReads;
	}

	// opens the read archive
	void CReadWriter::Open(const string& filename, const ReadStatus rs, const ReadGroup& readGroup) {

		if(mIsOpen) {
			cout << "ERROR: An attempt was made to open an already open read archive." << endl;
			exit(1);
		}

		mOutputFilename = filename;

		if(fopen_s(&mOutStream, filename.c_str(), "wb") != 0) {
			cout << "ERROR: Could not open the compressed read archive (" << mOutputFilename << ") for writing." << endl;
			exit(1);
		}

		mIsOpen = true;

		// initialization
		mBufferPosition     = 0;
		mPartitionMembers   = 0;

		// ================
		// write the header
		// ================

		// MOSAIK_SIGNATURE[6]       0  -  5
		// STATUS[1]                 6  -  6
		// SEQUENCING_TECHNOLOGY[1]  7  -  7
		// ARCHIVE_DATE[8]           8  - 15
		// NUM_READS[8]              16 - 23
		// NUM_BASES[8]              24 - 31
		// MEDIAN_FRAGMENT_LENGTH[4] 32 - 35
		// CENTER_NAME_LEN[1]        36 - 36
		// LIBRARY_NAME_LEN[1]       37 - 37
		// PLATFORM_UNIT_LEN[1]      38 - 38
		// READ_GROUP_ID_LEN[1]      39 - 39
		// SAMPLE_NAME_LEN[1]        40 - 40
		// DESCRIPTION_LEN[2]        41 - 42
		// RESERVED[8]               43 - 50
		// CENTER_NAME[*]            51
		// DESCRIPTION[*]
		// LIBRARY_NAME[*]
		// PLATFORM_UNIT[*]
		// READ_GROUP_ID[*]
		// SAMPLE_NAME[*]

		// write the MOSAIK signature
		const unsigned char SIGNATURE_LENGTH = 6;
		const char* MOSAIK_SIGNATURE = "MSKRA\1";
		fwrite(MOSAIK_SIGNATURE, SIGNATURE_LENGTH, 1, mOutStream);

		// write the read status (currently single end or paired end)
		fputc((unsigned char)rs, mOutStream);

		// write the sequencing technology
		fputc((unsigned char)readGroup.SequencingTechnology, mOutStream);
		if(readGroup.SequencingTechnology == ST_SOLID) mIsSOLiD = true;

		// write the archive date
		uint64_t currentTime = CTimeSupport::GetSystemTime();
		fwrite((char*)&currentTime, SIZEOF_UINT64, 1, mOutStream);

		// skip the number of reads and bases
		fseek64(mOutStream, 2 * SIZEOF_UINT64, SEEK_CUR);

		// write the median fragment length
		fwrite((char*)&readGroup.MedianFragmentLength, SIZEOF_INT, 1, mOutStream);

		// write the metadata string lengths: the lengths are checked in BuildMain.cpp
		const unsigned char centerNameLen   = (unsigned char)readGroup.CenterName.size();
		const unsigned char libraryNameLen  = (unsigned char)readGroup.LibraryName.size();
		const unsigned char platformUnitLen = (unsigned char)readGroup.PlatformUnit.size();
		const unsigned char readGroupIDLen  = (unsigned char)readGroup.ReadGroupID.size();
		const unsigned char sampleNameLen   = (unsigned char)readGroup.SampleName.size();
		const unsigned short descriptionLen = (unsigned short)readGroup.Description.size();

		fputc(centerNameLen,   mOutStream);
		fputc(libraryNameLen,  mOutStream);
		fputc(platformUnitLen, mOutStream);
		fputc(readGroupIDLen,  mOutStream);
		fputc(sampleNameLen,   mOutStream);
		fwrite((char*)&descriptionLen, SIZEOF_SHORT, 1, mOutStream);

		// write the reserved bytes
		const uint64_t reserved = 0;
		fwrite((char*)&reserved, SIZEOF_UINT64, 1, mOutStream);

		// convert the center name to lowercase
		string centerName = readGroup.CenterName;
		CSequenceUtilities::LowercaseSequence(centerName);

		// write the metadata strings
		fwrite(centerName.c_str(),      centerNameLen,   1, mOutStream);
		fwrite(readGroup.Description.c_str(),  descriptionLen,  1, mOutStream);
		fwrite(readGroup.LibraryName.c_str(),  libraryNameLen,  1, mOutStream);
		fwrite(readGroup.PlatformUnit.c_str(), platformUnitLen, 1, mOutStream);
		fwrite(readGroup.ReadGroupID.c_str(),  readGroupIDLen,  1, mOutStream);
		fwrite(readGroup.SampleName.c_str(),   sampleNameLen,   1, mOutStream);
	}

	// saves the read to the read archive
	void CReadWriter::SaveRead(const Mosaik::Read& mr) {

		// initialize
		unsigned char readNameLen    = (unsigned char)mr.Name.Length();
		unsigned short numMate1Bases = (unsigned short)mr.Mate1.Bases.Length();
		unsigned short numMate2Bases = (unsigned short)mr.Mate2.Bases.Length();

		// return if both mates are empty
		if((numMate1Bases == 0) && (numMate2Bases == 0)) return;

		mNumBases += numMate1Bases + numMate2Bases;

		bool isPairedEnd = false;
		if(numMate2Bases > 0) isPairedEnd = true;

		// calculate the entry size
		unsigned short entrySize   = 2 * numMate1Bases + readNameLen + SIZEOF_SHORT + 2;
		if(isPairedEnd) entrySize += 2 * numMate2Bases + SIZEOF_SHORT;

		if(mIsSOLiD) {
			entrySize += SOLID_PREFIX_LENGTH;
			if(isPairedEnd) entrySize += SOLID_PREFIX_LENGTH;
		}

		// check the memory buffer
		if(mBufferPosition >= (mBufferLen - entrySize)) AdjustBuffer();

		// ============================
		// serialize data to our buffer
		// ============================

		// no need to leave space for entry size
		unsigned int bufferOffset = mBufferPosition;

		// store the read type
		mBuffer[bufferOffset++] = (isPairedEnd ? 1 : 0);

		// store the read name
		mBuffer[bufferOffset++] = readNameLen;

		memcpy(mBuffer + bufferOffset, mr.Name.CData(), readNameLen);
		bufferOffset += readNameLen;

		// ============
		// store mate 1
		// ============

		// store the read length
		memcpy(mBuffer + bufferOffset, (char*)&numMate1Bases, SIZEOF_SHORT);
		bufferOffset += SIZEOF_SHORT;

		// store the bases
		memcpy(mBuffer + bufferOffset, mr.Mate1.Bases.CData(), numMate1Bases);
		bufferOffset += numMate1Bases;

		if(mIsSOLiD) {
			memcpy(mBuffer + bufferOffset, mr.Mate1.SolidPrefixTransition, SOLID_PREFIX_LENGTH);
			bufferOffset += SOLID_PREFIX_LENGTH;
		}

		// store the qualities
		memcpy(mBuffer + bufferOffset, mr.Mate1.Qualities.CData(), numMate1Bases);
		bufferOffset += numMate1Bases;

		// ============
		// store mate 2
		// ============

		if(isPairedEnd) {

			// store the read length
			memcpy(mBuffer + bufferOffset, (char*)&numMate2Bases, SIZEOF_SHORT);
			bufferOffset += SIZEOF_SHORT;

			// store the bases
			memcpy(mBuffer + bufferOffset, mr.Mate2.Bases.CData(), numMate2Bases);
			bufferOffset += numMate2Bases;

			if(mIsSOLiD) {
				memcpy(mBuffer + bufferOffset, mr.Mate2.SolidPrefixTransition, SOLID_PREFIX_LENGTH);
				bufferOffset += SOLID_PREFIX_LENGTH;
			}

			// store the qualities
			memcpy(mBuffer + bufferOffset, mr.Mate2.Qualities.CData(), numMate2Bases);
			bufferOffset += numMate2Bases;
		}

		// check the buffer
		if(bufferOffset >= mBufferLen) {
			cout << endl << "ERROR: Buffer overrun detected when saving read. Used " << bufferOffset << " bytes, but allocated " << mBufferLen << " bytes." << endl;
			exit(1);
		}

		// update the partition variables and buffer position
		mPartitionMembers++;
		mBufferPosition = bufferOffset;

		// flush the buffer
		if(mPartitionMembers >= mPartitionSize) WritePartition();

		// increment the read counter
		mNumReads++;
	}

	// write partition to disk
	void CReadWriter::WritePartition(void) {

		// check the compression buffer size
		unsigned int requestedSize = (unsigned int)(mBufferPosition * 1.05);
		CMemoryUtilities::CheckBufferSize(mCompressionBuffer, mCompressionBufferLen, requestedSize);

		// compress the partition
		int compressedSize = fastlz_compress_level(FASTLZ_BETTER_COMPRESSION, mBuffer, mBufferPosition, mCompressionBuffer);

		// write the uncompressed partition entry size
		fwrite((char*)&mBufferPosition, SIZEOF_INT, 1, mOutStream);

		// write the compressed partition entry size
		fwrite((char*)&compressedSize, SIZEOF_INT, 1, mOutStream);

		// write the partition member size
		fwrite((char*)&mPartitionMembers, SIZEOF_SHORT, 1, mOutStream);

		// write the partition
		fwrite(mCompressionBuffer, compressedSize, 1, mOutStream);

		mPartitionMembers = 0;
		mBufferPosition   = 0;
	}
}
