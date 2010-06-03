// ***************************************************************************
// CReadReader - loads reads from the MOSAIK read archive.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "ReadReader.h"

namespace MosaikReadFormat {

	// constructor
	CReadReader::CReadReader(void)
		: mIsOpen(false)
		, mInStream(NULL)
		, mNumReads(0)
		, mNumBases(0)
		, mCurrentRead(0)
		, mBuffer(NULL)
		, mBufferPtr(NULL)
		, mBufferLen(0)
		, mCompressionBuffer(NULL)
		, mCompressionBufferLen(0)
		, mPartitionSize(0)
		, mPartitionMembers(0)
		, mIsSOLiD(false)
	{}

	// destructor
	CReadReader::~CReadReader(void) {
		if(mIsOpen)            Close();
		if(mBuffer)            delete [] mBuffer;
		if(mCompressionBuffer) delete [] mCompressionBuffer;
	}

	// validates the supplied read archive file
	bool CReadReader::CheckFile(const string& filename, SequencingTechnologies& st, ReadStatus& rs, const bool showError) {

		// read in the first 6 characters
		char signature[7];
		signature[6] = 0;
		bool foundError = false;

		const char* MOSAIK_SIGNATURE = "MSKRA\1";

		// open the MOSAIK read archive
		FILE* checkStream = NULL;
		if(fopen_s(&checkStream, filename.c_str(), "rb") != 0) {
			if(showError) {
				printf("ERROR: Could not open %s when validating the read archive.\n", filename.c_str());
				exit(1);
			}

			foundError = true;
		}

		// retrieve the MOSAIK read archive signature
		if(!foundError) {

			// check if we were able to read 6 bytes
			if(fread(signature, 1, 6, checkStream) < 6) {
				if(showError) {
					printf("ERROR: It seems that the input file (%s) is not in the MOSAIK read format.\n", filename.c_str());
					exit(1);
				}

				foundError = true;
			}

			// check if the read signatures match
			if(!foundError && (strncmp(signature, MOSAIK_SIGNATURE, 5) != 0)) {
				if(showError) {
					printf("ERROR: It seems that the input file (%s) is not in the MOSAIK read format.\n", filename.c_str());
					exit(1);
				}

				foundError = true;
			}

			// check if the file format is from another version
			if(!foundError && (MOSAIK_SIGNATURE[5] != signature[5])) {
				if(showError) {
					printf("ERROR: It seems that the input file (%s) was created in another version of MosaikBuild. A new read archive is required.\n", filename.c_str());
					printf("       file version: %hu, expected version: %hu\n", signature[5], MOSAIK_SIGNATURE[5]);
					exit(1);
				}

				foundError = true;
			}
		}

		// grab the read status and sequencing technology
		rs = RS_UNKNOWN;
		st = ST_UNKNOWN;
		if(!foundError) {
			rs = (ReadStatus)fgetc(checkStream);
			st = (SequencingTechnologies)fgetc(checkStream);
		}

		// close the file
		if(checkStream) fclose(checkStream);

		// return the appropriate values
		if(foundError) return false;
		return true;
	}

	// closes the read archive
	void CReadReader::Close(void) {
		mIsOpen = false;
		fclose(mInStream);
	}

	// gets the metadata object
	ReadGroup CReadReader::GetReadGroup(void) const {
		return mReadGroup;
	}

	// gets the archive read count
	uint64_t CReadReader::GetNumReads(void) const {
		if(!mIsOpen) return 0;
		return mNumReads;
	}

	// gets the read archive sequencing technology
	SequencingTechnologies CReadReader::GetSequencingTechnology(void) const {
		return mReadGroup.SequencingTechnology;
	}

	// gets the read archive status
	ReadStatus CReadReader::GetStatus(void) const {
		return mStatus;
	}

	// loads the next read from the read archive
	bool CReadReader::LoadNextRead(Mosaik::Read& mr) {

		if(!mIsOpen) {
			cout << "ERROR: An attempt was made to get reads from a read archive that hasn't been opened yet." << endl;
			exit(1);
		}

		// check if there are any more reads
		if(mCurrentRead >= mNumReads) return false;

		// ==================
		// read the partition
		// ==================

		if(mPartitionMembers == mPartitionSize) {

			// read the uncompressed partition entry size
			unsigned int uncompressedSize = 0;
			fread((char*)&uncompressedSize, SIZEOF_INT, 1, mInStream);

			if(feof(mInStream)) return false;

			// read the compressed partition entry size
			int compressedSize = 0;
			fread((char*)&compressedSize, SIZEOF_INT, 1, mInStream);

			// read the partition member size
			mPartitionMembers = 0;
			fread((char*)&mPartitionSize, SIZEOF_SHORT, 1, mInStream);

			// check the compression buffer size
			CMemoryUtilities::CheckBufferSize(mCompressionBuffer, mCompressionBufferLen, compressedSize);
			CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, uncompressedSize);

			// read and uncompress the partition
			int numBytesRead = (int)fread(mCompressionBuffer, 1, compressedSize, mInStream);

			if(numBytesRead != compressedSize) {
				cout << "ERROR: Tried to read " << compressedSize << " bytes, but received only " << numBytesRead << " bytes." << endl;
				exit(1);
			}

			int result = fastlz_decompress(mCompressionBuffer, compressedSize, mBuffer, mBufferLen);

			if(result == 0) {
				cout << "ERROR: Unable to properly uncompress the current data partition." << endl;
				exit(1);
			}

			// set the buffer pointer
			mBufferPtr = mBuffer;
		}

		// retrieve the read type
		const bool isPairedEnd = (*mBufferPtr == 0 ? false : true);
		mBufferPtr++;

		// retrieve the read name
		unsigned char readNameLen = *mBufferPtr;
		mBufferPtr++;

		mr.Name.Copy((const char*)mBufferPtr, readNameLen);
		mBufferPtr += readNameLen;

		// set the read group code
		mr.ReadGroupCode = mReadGroup.ReadGroupCode;

		// ===============
		// retrieve mate 1
		// ===============

		// retrieve the read length
		unsigned short numMate1Bases = 0;
		memcpy((char*)&numMate1Bases, mBufferPtr, SIZEOF_SHORT);
		mBufferPtr += SIZEOF_SHORT;

		// retrieve the bases
		mr.Mate1.Bases.Copy((const char*)mBufferPtr, numMate1Bases);
		mBufferPtr += numMate1Bases;

		if(mIsSOLiD) {
			memcpy(mr.Mate1.SolidPrefixTransition, mBufferPtr, SOLID_PREFIX_LENGTH);
			mBufferPtr += SOLID_PREFIX_LENGTH;
		}

		// retrieve the qualities
		mr.Mate1.Qualities.Copy((const char*)mBufferPtr, numMate1Bases);
		mBufferPtr += numMate1Bases;

		// ===============
		// retrieve mate 2
		// ===============

		if(isPairedEnd) {

			// retrieve the read length
			unsigned short numMate2Bases = 0;
			memcpy((char*)&numMate2Bases, mBufferPtr, SIZEOF_SHORT);
			mBufferPtr += SIZEOF_SHORT;

			// retrieve the bases
			mr.Mate2.Bases.Copy((const char*)mBufferPtr, numMate2Bases);
			mBufferPtr += numMate2Bases;

			if(mIsSOLiD) {
				memcpy(mr.Mate2.SolidPrefixTransition, mBufferPtr, SOLID_PREFIX_LENGTH);
				mBufferPtr += SOLID_PREFIX_LENGTH;
			}

			// retrieve the qualities
			mr.Mate2.Qualities.Copy((const char*)mBufferPtr, numMate2Bases);
			mBufferPtr += numMate2Bases;

		} else {

			// reset mate2
			mr.Mate2.Bases.SetLength(0);
			mr.Mate2.Qualities.SetLength(0);
		}

		// increment the read counter
		mCurrentRead++;
		mPartitionMembers++;

		return true;
	}

	// opens the read archive
	void CReadReader::Open(const string& filename) {

		if(mIsOpen) {
			cout << "ERROR: An attempt was made to open an already open read archive." << endl;
			exit(1);
		}

		mInputFilename = filename;

		if(fopen_s(&mInStream, filename.c_str(), "rb") != 0) {
			cout << "ERROR: Could not open the compressed read archive (" << mInputFilename << ") for reading." << endl;
			exit(1);
		}

		mIsOpen = true;

		// ===============
		// read the header
		// ===============

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

		// skip the MOSAIK signature
		const unsigned char SIGNATURE_LENGTH = 6;
		fseek64(mInStream, SIGNATURE_LENGTH, SEEK_SET);

		// read the read status (currently single end or paired end)
		mStatus = (ReadStatus)fgetc(mInStream);

		// read the sequencing technology
		mReadGroup.SequencingTechnology = (SequencingTechnologies)fgetc(mInStream);
		if(mReadGroup.SequencingTechnology == ST_SOLID) mIsSOLiD = true;

		// skip the archive date
		fseek64(mInStream, SIZEOF_UINT64, SEEK_CUR);

		// retrieve the number of reads
		fread((char*)&mNumReads, SIZEOF_UINT64, 1, mInStream);

		// retrieve the number of bases
		fread((char*)&mNumBases, SIZEOF_UINT64, 1, mInStream);

		// read the median fragment length
		fread((char*)&mReadGroup.MedianFragmentLength, SIZEOF_INT, 1, mInStream);

		// read the metadata string lengths
		const unsigned char centerNameLen   = (unsigned char)fgetc(mInStream);
		const unsigned char libraryNameLen  = (unsigned char)fgetc(mInStream);
		const unsigned char platformUnitLen = (unsigned char)fgetc(mInStream);
		const unsigned char readGroupIDLen  = (unsigned char)fgetc(mInStream);
		const unsigned char sampleNameLen   = (unsigned char)fgetc(mInStream);

		unsigned short descriptionLen = 0;
		fread((char*)&descriptionLen, SIZEOF_SHORT, 1, mInStream);

		mReadGroup.CenterName.resize(centerNameLen);
		mReadGroup.LibraryName.resize(libraryNameLen);
		mReadGroup.PlatformUnit.resize(platformUnitLen);
		mReadGroup.ReadGroupID.resize(readGroupIDLen);
		mReadGroup.SampleName.resize(sampleNameLen);
		mReadGroup.Description.resize(descriptionLen);

		// skip the reserved bytes
		fseek64(mInStream, SIZEOF_UINT64, SEEK_CUR);

		// read the metadata strings
		if(centerNameLen > 0)   fread((void*)mReadGroup.CenterName.data(),   centerNameLen,   1, mInStream);
		if(descriptionLen > 0)  fread((void*)mReadGroup.Description.data(),  descriptionLen,  1, mInStream);
		if(libraryNameLen > 0)  fread((void*)mReadGroup.LibraryName.data(),  libraryNameLen,  1, mInStream);
		if(platformUnitLen > 0) fread((void*)mReadGroup.PlatformUnit.data(), platformUnitLen, 1, mInStream);

		fread((void*)mReadGroup.ReadGroupID.data(),  readGroupIDLen,  1, mInStream);
		fread((void*)mReadGroup.SampleName.data(),   sampleNameLen,   1, mInStream);

		// create our read group code
		mReadGroup.ReadGroupCode = ReadGroup::GetCode(mReadGroup);

		// get the reads offset
		mReadsOffset = ftell64(mInStream);

		//// DEBUG
		//cout << "# reads:                " << mNumReads << endl;
		//cout << "# bases:                " << mNumBases << endl;
		//cout << "median fragment length: " << mReadGroup.MedianFragmentLength << endl;
		//cout << "center name:            " << mReadGroup.CenterName << endl;
		//cout << "library name:           " << mReadGroup.LibraryName << endl;
		//cout << "platform unit:          " << mReadGroup.PlatformUnit << endl;
		//cout << "read group ID:          " << mReadGroup.ReadGroupID << endl;
		//cout << "sample name:            " << mReadGroup.SampleName << endl;
		//cout << "description:            " << mReadGroup.Description << endl;
		//exit(1);
	}

	// sets the file pointer to the beginning of the read data
	void CReadReader::Rewind(void) {
		fseek64(mInStream, mReadsOffset, SEEK_SET);
		mCurrentRead      = 0;
		mPartitionMembers = 0;
		mPartitionSize    = 0;
	}
}
