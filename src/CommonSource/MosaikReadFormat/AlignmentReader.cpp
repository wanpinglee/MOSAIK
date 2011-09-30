// ***************************************************************************
// CAlignmentReader - loads alignments from the MOSAIK alignment archive.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "AlignmentReader.h"

namespace MosaikReadFormat {

	// define our MOSAIK file signature
	//const char* CAlignmentReader::MOSAIK_SIGNATURE         = "MSKAA\4";
	//const unsigned char CAlignmentReader::SIGNATURE_LENGTH = 6;

	// constructor
	CAlignmentReader::CAlignmentReader(void)
		: mIsOpen(false)
		, mInStream(NULL)
		, mNumReads(0)
		, mNumBases(0)
		, mCurrentRead(0)
		, mReadsOffset(0)
		, mReferenceGapOffset(0)
		, mIndexOffset(0)
		, mBuffer(NULL)
		, mBufferPtr(NULL)
		, mBufferLen(0)
		, mCompressionBuffer(NULL)
		, mCompressionBufferLen(0)
		, mPartitionSize(0)
		, mPartitionMembers(0)
		//, mRefSeqLUT(NULL)
		, mStatus(AS_UNKNOWN)
		, mSeqTech(ST_UNKNOWN)
		, MosaikSignature(NULL)
	{}

	// destructor
	CAlignmentReader::~CAlignmentReader(void) {
		if(mIsOpen)            Close();
		//if(mBuffer)            delete mBuffer;
		//if(mCompressionBuffer) delete mCompressionBuffer;
		//if(MosaikSignature)    delete [] MosaikSignature;

		// delete the reference sequence LUT
		// Mark this will cause memory leakage; however, it causes segmentation fault
		//for(unsigned short i = 0; i < mNumRefSeqs; ++i) delete [] mRefSeqLUT[i];
		//if ( mRefSeqLUT ) delete mRefSeqLUT;

		mBuffer            = NULL;
		mCompressionBuffer = NULL;
		MosaikSignature    = NULL;
		//mRefSeqLUT         = NULL;
	}

	// checks to see if this is truly an MOSAIK alignment archive
	bool CAlignmentReader::CheckFile(const string& filename, SequencingTechnologies& st, AlignmentStatus& as, const bool showError) {

		// read in the first 6 characters
		char signature[SIGNATURE_LENGTH + 1];
		signature[SIGNATURE_LENGTH] = 0;
		bool foundError = false;

		// open the MOSAIK alignment archive
		FILE* checkStream = NULL;
		if(fopen_s(&checkStream, filename.c_str(), "rb") != 0) {
			if(showError) {
				printf("ERROR: Could not open %s when validating the alignment archive.\n", filename.c_str());
				exit(1);
			}

			foundError = true;
		}

		// retrieve the MOSAIK alignment archive signature
		if(!foundError) {

			// check if we were able to read 6 bytes
			if(fread(signature, 1, SIGNATURE_LENGTH, checkStream) < 6) {
				if(showError) {
					printf("ERROR: It seems that the input file (%s) is not in the MOSAIK read format.\n", filename.c_str());
					exit(1);
				}

				foundError = true;
			}

			// check if the read signatures match
			//if(!foundError && (strncmp(signature, MOSAIK_SIGNATURE, 5) != 0)) {
			if ( ( strncmp( signature, ALIGNER_SIGNATURE, SIGNATURE_LENGTH - 1 ) != 0 )
			   &&( strncmp( signature, SORT_SIGNATURE, SIGNATURE_LENGTH - 1 ) != 0 ) ) {
			//if(!foundError && (strncmp(signature, ALIGNER_SIGNATURE, 5) != 0)) {
				if(showError) {
					printf("ERROR: It seems that the input file (%s) is not in the MOSAIK alignment format.\n", filename.c_str());
					exit(1);
				}

				foundError = true;
			}

			// check if the file format is from another version
			if( !foundError && ( ( signature[5] != ALIGNER_SIGNATURE[5] ) && ( signature[5] != ALIGNER_SIGNATURE5[5] ) && ( signature[5] != SORT_SIGNATURE[5] ) ) ) {
			//if( !foundError && ( ( signature[5] != ALIGNER_SIGNATURE[5] ) && ( signature[5] != ALIGNER_SIGNATURE5[5] ) ) ) {
				if(showError) {
					//char version = ( strncmp( signature, ALIGNER_SIGNATURE, SIGNATURE_LENGTH - 1 ) == 0 ) ? ALIGNER_SIGNATURE[5] : SORT_SIGNATURE[5];

					//printf("ERROR: It seems that the input file (%s) was created in another version of MosaikAligner. This version of MOSAIK expected to find an alignment archive using version: %hu, but the alignment archive uses version: %hu. A new alignment archive is required.\n", filename.c_str(), version, signature[5]);
					printf("ERROR: It seems that the input file (%s) was created in another version of MosaikAligner."
					       " This version of MOSAIK expected to find an alignment archive using version: 4 or 5,"
					       " but the alignment archive uses version: %hu. A new alignment archive is required.\n", filename.c_str(), signature[5]);
					exit(1);
				}

				foundError = true;
			}
		}

		// grab the alignment status and sequencing technology
		as = AS_UNKNOWN;
		st = ST_UNKNOWN;
		if(!foundError) {
			as = (AlignmentStatus)fgetc(checkStream);
			fread((char*)&st, SIZEOF_SHORT, 1, checkStream);
		}

		// close the file
		if( checkStream != NULL ) fclose(checkStream);

		// return the appropriate values
		if(foundError) return false;
		return true;
	}

	// closes the alignment archive
	void CAlignmentReader::Close(void) {
		mIsOpen = false;
		mRefSeqLUT.clear();
		mReferenceSequences.clear();
		mRefSeqGaps.clear();
		mReadGroups.clear();
		mHeaderTags.clear();
		mReadGroupLUT.clear();

		if(mBuffer)            delete mBuffer;
		if(mCompressionBuffer) delete mCompressionBuffer;

		mBuffer            = NULL;

		if ( MosaikSignature ) delete [] MosaikSignature;
		MosaikSignature = NULL;
		fclose(mInStream);
	}

	// returns the a pointer to the header tags map
	map<unsigned char, Tag>* CAlignmentReader::GetHeaderTags(void) {
		return &mHeaderTags;
	}

	// returns the number of reads in the archive
	uint64_t CAlignmentReader::GetNumBases(void) const {
		if(!mIsOpen) return 0;
		return mNumBases;
	}

	// returns the number of reads in the archive
	uint64_t CAlignmentReader::GetNumReads(void) const {
		if(!mIsOpen) return 0;
		return mNumReads;
	}

	// retrieves the read group given a read group code
	ReadGroup CAlignmentReader::GetReadGroupFromCode(const unsigned int code) {
		map<unsigned int, ReadGroup>::const_iterator rgIter = mReadGroupLUT.find(code);

		if(rgIter == mReadGroupLUT.end()) {
			printf("ERROR: The following read group ID was not found in the lookup table: %u\n", code);
			exit(1);
		}

		return rgIter->second;
	}

	// retrieves the read groups vector
	void CAlignmentReader::GetReadGroups(vector<ReadGroup>& readGroups) const {
		readGroups.resize(mReadGroups.size());
		vector<ReadGroup>::const_iterator rgIter;
		vector<ReadGroup>::iterator crgIter = readGroups.begin();

		// force a deep copy
		for(rgIter = mReadGroups.begin(); rgIter != mReadGroups.end(); ++rgIter, ++crgIter) {
			*crgIter = *rgIter;
		}
	}

	// retrieves the reference sequence gaps
	vector<vector<GapInfo> >* CAlignmentReader::GetReferenceSequenceGaps(void) {
		return &mRefSeqGaps;
	}

	// if the reference name is found, the referenceIndex will be set and the function will return true
	bool CAlignmentReader::GetReferenceSequenceIndex(const string& referenceName, unsigned int& referenceIndex) const {
		
		// initialize
		referenceIndex = 0;
		bool foundReferenceName = false;
		vector<ReferenceSequence>::const_iterator rsIter;
		
		// search for the correct reference index
		for(rsIter = mReferenceSequences.begin(); rsIter != mReferenceSequences.end(); ++rsIter, ++referenceIndex) {
			if(rsIter->Name == referenceName) {
				foundReferenceName = true;
				break;
			}
		}

		// return true if we found the reference name
		if(!foundReferenceName) return false;
		return true;
	}

	// retrieves the reference sequence data
	//vector<ReferenceSequence>* CAlignmentReader::GetReferenceSequences(void) {
	//	return &mReferenceSequences;
	//}

	// retrieves the reference sequence data
	void CAlignmentReader::GetReferenceSequences( vector<ReferenceSequence>& refVec) {
		refVec = mReferenceSequences;
	}

	// gets the alignment archive sequencing technology
	SequencingTechnologies CAlignmentReader::GetSequencingTechnology(void) const {
		return mSeqTech;
	}

	// retrieves the signature
	void CAlignmentReader::GetSignature ( char*& signature ) {
		if ( signature ) delete [] signature;

		signature = new char [ SIGNATURE_LENGTH + 1 ];
		memcpy( signature, MosaikSignature, SIGNATURE_LENGTH );
		signature[ SIGNATURE_LENGTH ] = 0;
	}

	// retrieves the file status
	AlignmentStatus CAlignmentReader::GetStatus(void) const {
		return mStatus;
	}

	// jumps to the block containing the specified reference index and position
	void CAlignmentReader::Jump(const unsigned int referenceIndex, const unsigned int referencePosition) {

		// ===============
		// parse the index
		// ===============

		if(mIndexOffset == 0) {
			cout << "ERROR: Cannot jump to the desired compressed block because the index offset was not set." << endl;
			exit(1);
		}

		// jump to the index offset and read the number of entries
		fseek64(mInStream, mIndexOffset, SEEK_SET);

		unsigned int numIndexEntries = 0;
		fread((char*)&numIndexEntries, SIZEOF_INT, 1, mInStream);

		// load the index
		CFastLZIO fio;
		char* pBuffer = mBuffer;
		fio.Read(pBuffer, mBufferLen, mInStream);
		mBuffer = pBuffer;

		// find the block containing the specified reference index and position
		unsigned int bufferOffset = 0;

		unsigned int index    = 0;
		unsigned int position = 0;
		off_type offset       = 0;

		bool foundBlock = false;
		for(unsigned int i = 0; i < numIndexEntries; ++i) {

			// retrieve the reference index
			memcpy((char*)&index, mBuffer + bufferOffset, SIZEOF_INT);
			bufferOffset += SIZEOF_INT;

			// store the reference position
			memcpy((char*)&position, mBuffer + bufferOffset, SIZEOF_INT);
			bufferOffset += SIZEOF_INT;

			// store the file offset
			memcpy((char*)&offset, mBuffer + bufferOffset, SIZEOF_UINT64);
			bufferOffset += SIZEOF_UINT64;

			// keep going until we find a compression block that is past our desired index and position
			if(index > referenceIndex) foundBlock = true;
			if((index == referenceIndex) && (position >= referencePosition)) foundBlock = true;
			if(foundBlock) break;
		}

		if(!foundBlock) {
			cout << "ERROR: A suitable compression block was not found in the index." << endl;
			exit(1);
		}

		fseek64(mInStream, offset, SEEK_SET);
		mCurrentRead      = 0;
		mPartitionMembers = 0;
		mPartitionSize    = 0;
	}

/*
	// loads the next alignment from the alignment archive
	bool CAlignmentReader::LoadNextAlignment(Alignment& al) {

		if(!mIsOpen) {
			cout << "ERROR: An attempt was made to get reads from an alignment archive that hasn't been opened yet." << endl;
			exit(1);
		}

		// check if we have already processed all of the reads
		if(mCurrentRead >= mNumReads) return false;

		// read the partition
		if(mPartitionMembers == mPartitionSize) {
			if(!ReadPartition()) return false;
		}

		// initialize
		unsigned char readStatus = RF_UNKNOWN;

		unsigned int numMate1Alignments = 0;
		unsigned int numMate2Alignments = 0;
		unsigned int numMate1OriginalAlignments = 0;
		unsigned int numMate2OriginalAlignments = 0;

		// load the read header
		LoadReadHeader(al.Name, al.ReadGroupCode, readStatus, numMate1Alignments, numMate2Alignments, numMate1OriginalAlignments, numMate2OriginalAlignments);

		// interpret the read status
		const bool isLongRead           = ((readStatus & RF_IS_LONG_READ)            != 0 ? true : false);
		const bool isPairedInSequencing = ((readStatus & RF_IS_PAIRED_IN_SEQUENCING) != 0 ? true : false);
		const bool isResolvedAsPair     = ((readStatus & RF_RESOLVED_AS_PAIR)        != 0 ? true : false);

		// deserialize the alignment
		ReadAlignment(al, isLongRead, isPairedInSequencing, isResolvedAsPair, numMate1OriginalAlignments, numMate2OriginalAlignments);

		// increment the read counter
		++mCurrentRead;
		++mPartitionMembers;

		return true;
	}
*/

	// loads the next read from the alignment archive
	bool CAlignmentReader::LoadNextRead(Mosaik::AlignedRead& ar) {

		if(!mIsOpen) {
			cout << "ERROR: An attempt was made to get reads from an alignment archive that hasn't been opened yet." << endl;
			exit(1);
		}

		// check if we have already processed all of the reads
		if(mCurrentRead >= mNumReads) return false;

		// read the partition
		if(mPartitionMembers == mPartitionSize) {
			if(!ReadPartition()) return false;
		}

		// initialize
		unsigned char readStatus = RF_UNKNOWN;

		int numMate1Alignments = 0;
		int numMate2Alignments = 0;
		int numMate1OriginalAlignments = 0;
		int numMate2OriginalAlignments = 0;
		int numMate1Hash = 0;
		int numMate2Hash = 0;
		
		// load the read header
		LoadReadHeader(ar.Name, ar.ReadGroupCode, readStatus, numMate1Alignments, numMate2Alignments, 
		    numMate1OriginalAlignments, numMate2OriginalAlignments, numMate1Hash, numMate2Hash);

		// interpret the read status
		const bool haveMate1        = ((readStatus & RF_HAVE_MATE1)              != 0 ? true : false);
		const bool haveMate2        = ((readStatus & RF_HAVE_MATE2)              != 0 ? true : false);
		const bool isResolvedAsPair = ((readStatus & RF_RESOLVED_AS_PAIR)        != 0 ? true : false);
		ar.hasCsString              = ((readStatus & RF_HAS_CS_STRING)           != 0 ? true : false);
		ar.IsLongRead               = ((readStatus & RF_IS_LONG_READ)            != 0 ? true : false);
		ar.IsPairedEnd              = ((readStatus & RF_IS_PAIRED_IN_SEQUENCING) != 0 ? true : false);
		ar.IsResolvedAsPair         = ((readStatus & RF_RESOLVED_AS_PAIR)        != 0 ? true : false);


		// =================================
		// deserialize each mate 1 alignment
		// =================================

		ar.Mate1Alignments.resize(numMate1Alignments);
		if (haveMate1) 
			ReadAlignments(ar.Mate1Alignments, ar.IsLongRead, ar.IsPairedEnd, 
				isResolvedAsPair, ar.ReadGroupCode, numMate1OriginalAlignments,
				numMate2OriginalAlignments, numMate1Hash, numMate2Hash, ar.hasCsString);

		// =================================
		// deserialize each mate 2 alignment
		// =================================

		ar.Mate2Alignments.resize(numMate2Alignments);
		if (haveMate2) 
			ReadAlignments(ar.Mate2Alignments, ar.IsLongRead, ar.IsPairedEnd, 
				isResolvedAsPair, ar.ReadGroupCode, numMate1OriginalAlignments, 
				numMate2OriginalAlignments, numMate1Hash, numMate2Hash, ar.hasCsString);

		// increment the read counter
		++mCurrentRead;
		++mPartitionMembers;

		return true;
	}

	// load the read header from disk
	void CAlignmentReader::LoadReadHeader(
		CMosaikString& readName, 
		unsigned int&  readGroupCode, 
		unsigned char& readStatus, 
		int&  numMate1Alignments, 
		int&  numMate2Alignments,
		int&  numMate1OriginalAlignments,
		int&  numMate2OriginalAlignments,
		int&  numMate1Hashes,
		int&  numMate2Hashes) {

		// get the read name
		const unsigned char readNameLength = (unsigned char)*mBufferPtr;
		++mBufferPtr;

		readName.Copy((const char*)mBufferPtr, readNameLength);
		mBufferPtr += readNameLength;

		// get the read group code
		memcpy((char*)&readGroupCode, mBufferPtr, SIZEOF_INT);
		mBufferPtr += SIZEOF_INT;

		// get the read status
		readStatus = (unsigned char)*mBufferPtr;
		++mBufferPtr;

		const bool haveMate1 = ((readStatus & RF_HAVE_MATE1) != 0 ? true : false);
		const bool haveMate2 = ((readStatus & RF_HAVE_MATE2) != 0 ? true : false);

		// get the number of mate 1 alignments
		if(haveMate1) {
			memcpy((char*)&numMate1Alignments, mBufferPtr, SIZEOF_INT);
			mBufferPtr += SIZEOF_INT;
			memcpy((char*)&numMate1OriginalAlignments, mBufferPtr, SIZEOF_INT);
			mBufferPtr += SIZEOF_INT;
			memcpy((char*)&numMate1Hashes, mBufferPtr, SIZEOF_INT);
			mBufferPtr += SIZEOF_INT;
		}

		// get the number of mate 2 alignments
		if(haveMate2) {
			memcpy((char*)&numMate2Alignments, mBufferPtr, SIZEOF_INT);
			mBufferPtr += SIZEOF_INT;
			memcpy((char*)&numMate2OriginalAlignments, mBufferPtr, SIZEOF_INT);
			mBufferPtr += SIZEOF_INT;
			memcpy((char*)&numMate2Hashes, mBufferPtr, SIZEOF_INT);
			mBufferPtr += SIZEOF_INT;
		}
	}

	// opens the alignment archive
	void CAlignmentReader::Open(const string& filename) {

		if(mIsOpen) {
			cout << "ERROR: An attempt was made to open an already open alignment archive." << endl;
			exit(1);
		}

		mInputFilename = filename;

		mInStream = NULL;
		if(fopen_s(&mInStream, filename.c_str(), "rb") != 0) {
			cout << "ERROR: Could not open the compressed alignment archive (" << mInputFilename << ") for reading." << endl;
			exit(1);
		}

		mIsOpen = true;
		
		// ===============
		// read the header
		// ===============

		// MOSAIK_SIGNATURE[6]	   0  -  5
		// STATUS[1]               6  -  6
		// SEQUENCE_TECHNOLOGY[2]  7  -  8
		// ARCHIVE_DATE[8]		   9  - 16
		// NUM_REFERENCE_SEQS[4]   17 - 20
		// NUM_READ_GROUPS[4]      21 - 24
		// NUM_READS[8]            25 - 32
		// NUM_BASES[8]            33 - 40
		// REFERENCES_OFFSET[8]    41 - 48
		// REFERENCE_GAP_OFFSET[8] 49 - 57
		// INDEX_OFFSET[8]         58 - 63
		// NUM_READ_GROUP_TAGS[1]  64 - 64
		// READ_GROUPS[*]

		// check the MOSAIK signature
		char signature[SIGNATURE_LENGTH + 1];
		signature[SIGNATURE_LENGTH] = 0;
		fread( signature, SIGNATURE_LENGTH, 1, mInStream );

		// check if the read signatures match
		//if(strncmp(signature, MOSAIK_SIGNATURE, 5) != 0) {
		if ( ( strncmp( signature, ALIGNER_SIGNATURE, SIGNATURE_LENGTH - 1 ) != 0 ) && ( strncmp( signature, SORT_SIGNATURE, SIGNATURE_LENGTH - 1 ) != 0 ) ) {
		//if(strncmp(signature, ALIGNER_SIGNATURE, 5) != 0) {
			printf("ERROR: It seems that the input file (%s) is not in the MOSAIK alignment format.\n", 
				filename.c_str());
			exit(1);
		}

		//if(MOSAIK_SIGNATURE[5] != signature[5]) {
		if ( ( signature[5] != ALIGNER_SIGNATURE[5] ) && ( signature[5] != ALIGNER_SIGNATURE5[5] ) && ( signature[5] != SORT_SIGNATURE[5] ) ) {
		//if ( ( signature[5] != ALIGNER_SIGNATURE[5] ) && ( signature[5] != ALIGNER_SIGNATURE5[5] ) ) {
			//char version = ( strncmp( signature, ALIGNER_SIGNATURE, SIGNATURE_LENGTH - 1 ) == 0 ) ? ALIGNER_SIGNATURE[5] : SORT_SIGNATURE[5];
			//printf("ERROR: It seems that the input file (%s) was created in another version of MosaikAligner. "
			//	"This version of MOSAIK expected to find an alignment archive using version: %hu, but the "
			//	"alignment archive uses version: %hu. A new alignment archive is required.\n", 
			//	filename.c_str(), version, signature[5]);
			
			printf("ERROR: It seems that the input file (%s) was created in another version of MosaikAligner. "
				"This version of MOSAIK expected to find an alignment archive using version: 4 or 5, but the "
				"alignment archive uses version: %hu. A new alignment archive is required.\n", 
				filename.c_str(), signature[5]);
			exit(1);
		}

		MosaikSignature = new char [ SIGNATURE_LENGTH + 1 ];
		memcpy( MosaikSignature, signature, SIGNATURE_LENGTH );
		MosaikSignature[ SIGNATURE_LENGTH ] = 0;

		// retrieve the alignment file status
		mStatus = (AlignmentStatus)fgetc(mInStream);

		// retrieve the sequencing technology
		fread((char*)&mSeqTech, SIZEOF_SHORT, 1, mInStream);

		// skip the archive date
		fseek64(mInStream, SIZEOF_UINT64, SEEK_CUR);

		// retrieve the number of reference sequences
		fread((char*)&mNumRefSeqs, SIZEOF_INT, 1, mInStream);

		// retrieve the number of read groups
		unsigned int numReadGroups;
		fread((char*)&numReadGroups, SIZEOF_INT, 1, mInStream);

		// retrieve the number of reads
		fread((char*)&mNumReads, SIZEOF_UINT64, 1, mInStream);

		if(mNumReads == 0) {
			printf("ERROR: The alignment archive header indicates that no reads are contained in\n");
			printf("       this file. This might happen when the file was not closed properly -\n");
			printf("       usually from a killed process or a crash. Your only recourse is to\n");
			printf("       realign this data set.\n");
			printf("       filename: [%s]\n", filename.c_str());
			exit(1);
		}

		// retrieve the number of bases
		fread((char*)&mNumBases, SIZEOF_UINT64, 1, mInStream);

		// retrieve the references offset
		off_type referencesOffset = 0;
		fread((char*)&referencesOffset, SIZEOF_OFF_TYPE, 1, mInStream);

		// retrieve the reference gaps offset
		fread((char*)&mReferenceGapOffset, SIZEOF_OFF_TYPE, 1, mInStream);

		// retrieve the index offset
		fread((char*)&mIndexOffset, SIZEOF_OFF_TYPE, 1, mInStream);

		// retrieve the number of header tags
		const unsigned char numHeaderTags = (unsigned char)fgetc(mInStream);

		if(numHeaderTags != 0) {
			for(unsigned char j = 0; j < numHeaderTags; j++) {
				Tag tag;
				ReadTag(tag);
				mHeaderTags[tag.ID] = tag;
			}
		}

		// DEBUG
		//cout << "mStatus:             " << (short)mStatus << endl;
		//cout << "mSeqTech:            " << mSeqTech << endl;
		//cout << "mNumRefSeqs:         " << mNumRefSeqs << endl;
		//cout << "numReadGroups:       " << numReadGroups << endl;
		//cout << "mNumReads:           " << mNumReads << endl;
		//cout << "mNumBases:           " << mNumBases << endl;
		//cout << "referencesOffset:    " << referencesOffset << endl;
		//cout << "mReferenceGapOffset: " << mReferenceGapOffset << endl;
		//cout << "mIndexOffset:        " << mIndexOffset << endl;
		//cout << "numHeaderTags:       " << (unsigned short)numHeaderTags << endl << endl;

		// retrieve the read groups
		mReadGroups.resize(numReadGroups);

		vector<ReadGroup>::iterator rgIter;
		for(rgIter = mReadGroups.begin(); rgIter != mReadGroups.end(); ++rgIter) {

			// read the metadata string lengths
			const unsigned char centerNameLen   = (unsigned char)fgetc(mInStream);
			const unsigned char libraryNameLen  = (unsigned char)fgetc(mInStream);
			const unsigned char platformUnitLen = (unsigned char)fgetc(mInStream);
			const unsigned char readGroupIDLen  = (unsigned char)fgetc(mInStream);
			const unsigned char sampleNameLen   = (unsigned char)fgetc(mInStream);

			unsigned short descriptionLen = 0;
			fread((char*)&descriptionLen, SIZEOF_SHORT, 1, mInStream);
			fread((char*)&rgIter->SequencingTechnology, SIZEOF_SHORT, 1, mInStream);
			fread((char*)&rgIter->MedianFragmentLength, SIZEOF_INT, 1, mInStream);

			rgIter->CenterName.resize(centerNameLen);
			rgIter->LibraryName.resize(libraryNameLen);
			rgIter->PlatformUnit.resize(platformUnitLen);
			rgIter->ReadGroupID.resize(readGroupIDLen);
			rgIter->SampleName.resize(sampleNameLen);
			rgIter->Description.resize(descriptionLen);

			// read the metadata strings
			fread((void*)rgIter->CenterName.data(),   centerNameLen,   1, mInStream);
			fread((void*)rgIter->Description.data(),  descriptionLen,  1, mInStream);
			fread((void*)rgIter->LibraryName.data(),  libraryNameLen,  1, mInStream);
			fread((void*)rgIter->PlatformUnit.data(), platformUnitLen, 1, mInStream);
			fread((void*)rgIter->ReadGroupID.data(),  readGroupIDLen,  1, mInStream);
			fread((void*)rgIter->SampleName.data(),   sampleNameLen,   1, mInStream);
			
			// set the read group code
			rgIter->ReadGroupCode = ReadGroup::GetCode(*rgIter);

			// add the read group to our LUT
			mReadGroupLUT[rgIter->ReadGroupCode] = *rgIter;

			// retrieve the number of read group tags
			const unsigned char numReadGroupTags = (unsigned char)fgetc(mInStream);

			if(numReadGroupTags != 0) {
				printf("ERROR: Found %u read group tags, but support for read group tags has not been implemented yet.\n", numReadGroupTags);
				exit(1);
			}

			//// DEBUG
			//cout << "center name:            " << rgIter->CenterName << endl;
			//cout << "description:            " << rgIter->Description << endl;
			//cout << "library name:           " << rgIter->LibraryName << endl;
			//cout << "platform unit:          " << rgIter->PlatformUnit << endl;
			//cout << "read group ID:          " << rgIter->ReadGroupID << endl;
			//cout << "sample name:            " << rgIter->SampleName << endl;
			//cout << "sequencing technology:  " << rgIter->SequencingTechnology << endl;
			//cout << "median fragment length: " << rgIter->MedianFragmentLength << endl << endl;
		}

		// store the reads offset
		mReadsOffset = ftell64(mInStream);

		// ============================
		// read the reference sequences
		// ============================

		// jump to the reference sequence section
		fseek64(mInStream, referencesOffset, SEEK_SET);

		mReferenceSequences.resize(mNumRefSeqs);
		//mRefSeqLUT = new char*[mNumRefSeqs];
		mRefSeqLUT.resize( mNumRefSeqs );

		unsigned int currentRefSeq = 0;
		vector<ReferenceSequence>::iterator rsIter;
		for(rsIter = mReferenceSequences.begin(); rsIter != mReferenceSequences.end(); ++rsIter, ++currentRefSeq) {

			// REFERENCE_SEQ_NAME_LEN[1]                0 -  0 
			// REFERENCE_SEQ_SPECIES_LEN[1]             1 -  1
			// REFERENCE_SEQ_GENOME_ASSEMBLY_ID_LEN[1]  2 -  2
			// REFERENCE_SEQ_URI_LEN[1]                 3 -  3
			// REFERENCE_SEQ_NUM_BASES[4]               4 -  7
			// REFERENCE_SEQ_SEQ_OFFSET[8]              8 - 15
			// REFERENCE_SEQ_MD5[16]                   16 - 31
			// REFERENCE_SEQ_NAME[X]                   32 - XX
			// REFERENCE_SEQ_SPECIES[X]
			// REFERENCE_SEQ_GENOME_ASSEMBLY_ID[X]
			// REFERENCE_SEQ_URI[X]

			// read the name length
			const unsigned char nameLen = fgetc(mInStream);

			// read the species length
			const unsigned char speciesLen = fgetc(mInStream);

			// read the genome assembly id length
			const unsigned char genomeAssemblyIDLen = fgetc(mInStream);

			// read the uri length
			const unsigned char uriLen = fgetc(mInStream);

			// read the number of bases
			fread((char*)&rsIter->NumBases, SIZEOF_INT, 1, mInStream);

			// write the number of aligned reads
			fread((char*)&rsIter->NumAligned, SIZEOF_UINT64, 1, mInStream);

			// read the MD5 checksum
			rsIter->MD5.resize(32);
			char* pBuffer = (char*)rsIter->MD5.data();
			fread(pBuffer, 32, 1, mInStream);

			// read the reference name
			rsIter->Name.resize(nameLen);
			pBuffer = (char*)rsIter->Name.data();
			fread(pBuffer, nameLen, 1, mInStream);

			//mRefSeqLUT[currentRefSeq] = new char[nameLen + 1];
			//mRefSeqLUT[currentRefSeq].resize( nameLen + 1 );
			//memcpy(mRefSeqLUT[currentRefSeq], pBuffer, nameLen);
			mRefSeqLUT[currentRefSeq].insert( 0, pBuffer, nameLen );
			//mRefSeqLUT[currentRefSeq][nameLen] = 0;
			mRefSeqLUT[currentRefSeq].push_back(0);

			// read the species name
			if(speciesLen > 0) {
				rsIter->Species.resize(speciesLen);
				pBuffer = (char*)rsIter->Species.data();
				fread(pBuffer, speciesLen, 1, mInStream);
			}

			// read the genome assembly ID
			if(genomeAssemblyIDLen > 0) {
				rsIter->GenomeAssemblyID.resize(genomeAssemblyIDLen);
				pBuffer = (char*)rsIter->GenomeAssemblyID.data();
				fread(pBuffer, genomeAssemblyIDLen, 1, mInStream);
			}

			// read the URI
			if(uriLen > 0) {
				rsIter->URI.resize(uriLen);
				pBuffer = (char*)rsIter->URI.data();
				fread(pBuffer, uriLen, 1, mInStream);
			}

			// retrieve the number of reference sequence tags
			const unsigned char numReferenceSequenceTags = (unsigned char)fgetc(mInStream);

			if(numReferenceSequenceTags != 0) {
				printf("ERROR: Found reference sequence tags, but support for reference sequence tags has not been implemented yet.\n");
				exit(1);
			}

			//// DEBUG
			//cout << "# bases:                " << rsIter->NumBases << endl;
			//cout << "md5:                    " << rsIter->MD5 << endl;
			//cout << "name:                   " << rsIter->Name << endl;
			//cout << "species:                " << rsIter->Species << endl;
			//cout << "genome assembly ID:     " << rsIter->GenomeAssemblyID << endl;
			//cout << "URI:                    " << rsIter->URI << endl;
		}

		// ================================
		// read the reference sequence gaps
		// ================================

		CFastLZIO fio;
		if(mReferenceGapOffset != 0) {

			// jump to the reference gap location
			fseek64(mInStream, mReferenceGapOffset, SEEK_SET);

			// read the reference gaps vector
			fio.Read(mBuffer, mBufferLen, mInStream);

			unsigned int bufferOffset = 0;
			vector<GapInfo>::iterator gvIter;
			vector<vector<GapInfo> >::iterator rsgIter;

			mRefSeqGaps.resize(mNumRefSeqs);
			for(rsgIter = mRefSeqGaps.begin(); rsgIter != mRefSeqGaps.end(); ++rsgIter) {

				// retrieve the number of gaps for this reference sequence
				unsigned int numGaps = 0;
				memcpy((char*)&numGaps, mBuffer + bufferOffset, SIZEOF_INT);
				bufferOffset += SIZEOF_INT;

				// pre-allocate the reference gap vector
				rsgIter->resize(numGaps);

				for(gvIter = rsgIter->begin(); gvIter != rsgIter->end(); ++gvIter) {

					// retrieve the reference gap position
					memcpy((char*)&gvIter->Position, mBuffer + bufferOffset, SIZEOF_INT);
					bufferOffset += SIZEOF_INT;

					// retrieve the reference gap length
					memcpy((char*)&gvIter->Length, mBuffer + bufferOffset, SIZEOF_SHORT);
					bufferOffset += SIZEOF_SHORT;
				}
			}
		}

		// restore our file position
		Rewind();
	}

	// deserializes each alignment and stores them in the supplied vector
	void CAlignmentReader::ReadAlignments(
		vector<Alignment>&   alignments, 
		const bool&          isLongRead, 
		const bool&          isPairedInSequencing, 
		const bool&          isResolvedAsPair, 
		const unsigned int&  readGroupCode,
		const int&           numMate1OriginalAlignments,
		const int&           numMate2OriginalAlignments,
		const int&           numMate1Hashes,
		const int&           numMate2Hashes,
		const bool&          hasCsString) {
		
		vector<Alignment>::iterator alIter;
		for(alIter = alignments.begin(); alIter != alignments.end(); ++alIter) {
			ReadAlignment(*alIter, isLongRead, isPairedInSequencing, isResolvedAsPair,
			    numMate1OriginalAlignments, numMate2OriginalAlignments, 
			    numMate1Hashes, numMate2Hashes, hasCsString);
			alIter->ReadGroupCode = readGroupCode;
		}
	}

	// deserialize the alignment
	void CAlignmentReader::ReadAlignment(
		Alignment& al, 
		const bool& isLongRead, 
		const bool& isPairedInSequencing, 
		const bool& isResolvedAsPair,
		const int&  numMate1OriginalAlignments,
		const int&  numMate2OriginalAlignments,
		const int&  numMate1Hashes,
		const int&  numMate2Hashes,
		const bool& hasCsString) {

		// get the reference sequence start position
		memcpy((char*)&al.ReferenceBegin, mBufferPtr, SIZEOF_INT);
		mBufferPtr += SIZEOF_INT;

		// get the reference sequence end position
		memcpy((char*)&al.ReferenceEnd, mBufferPtr, SIZEOF_INT);
		mBufferPtr += SIZEOF_INT;

		// get the reference sequence index
		memcpy((char*)&al.ReferenceIndex, mBufferPtr, SIZEOF_INT);
		mBufferPtr += SIZEOF_INT;
		al.ReferenceName = (char*)mRefSeqLUT[al.ReferenceIndex].c_str();
		//strcpy( al.ReferenceName, mRefSeqLUT[al.ReferenceIndex].c_str() );

		// get the alignment's best and 2nd best Smith-Waterman scores
		memcpy((char*)&al.SwScore, mBufferPtr, SIZEOF_FLOAT);
		mBufferPtr += SIZEOF_FLOAT;
		memcpy((char*)&al.NextSwScore, mBufferPtr, SIZEOF_FLOAT);
		mBufferPtr += SIZEOF_FLOAT;

		// get the alignment's longest perfect match
		memcpy((char*)&al.NumLongestMatchs, mBufferPtr, SIZEOF_SHORT);
		mBufferPtr += SIZEOF_SHORT;

		// get the alignment status flag
		const unsigned char status = (unsigned char)*mBufferPtr;
		++mBufferPtr;

		al.IsFirstMate         = false;
		al.IsReverseStrand     = false;
		al.IsMateReverseStrand = false;
		al.WasRescued          = false;
		al.IsMapped            = true;
		al.IsFilteredOut       = false;

		al.IsPairedEnd      = isPairedInSequencing;
		al.IsResolvedAsPair = isResolvedAsPair;

		if(isPairedInSequencing) {
			if((status & AF_IS_FIRST_MATE) != 0) al.IsFirstMate = true;
			if(isResolvedAsPair && ((status & AF_IS_MATE_REVERSE_STRAND) != 0)) al.IsMateReverseStrand = true;
		}

		if((status & AF_IS_REVERSE_STRAND) != 0) al.IsReverseStrand = true;
		if((status & AF_WAS_RESCUED)       != 0) al.WasRescued      = true;
		if((status & AF_IS_UNMAPPED)       != 0) al.IsMapped        = false;
		if((status & AF_IS_JUNK)           != 0) al.IsJunk          = true;
		if((status & AF_IS_FILTEREDOUT)    != 0) al.IsFilteredOut   = true;

		if ( ( isPairedInSequencing && al.IsFirstMate ) || !isPairedInSequencing ) {
			al.NumMapped = numMate1OriginalAlignments;
			al.NumHash   = numMate1Hashes;
		}
		else {
			al.NumMapped = numMate2OriginalAlignments;
			al.NumHash   = numMate2Hashes;
		}

		if ( hasCsString ) {
			unsigned short csLen = 0;
			if( isLongRead ) {
				memcpy((char*)&csLen, mBufferPtr, SIZEOF_SHORT);
				mBufferPtr += SIZEOF_SHORT;
			} else {
				csLen = (unsigned char)*mBufferPtr;
				++mBufferPtr;
			}

                        // retrieve the colorspace raw sequence
			al.CsQuery.insert(0, (const char*)mBufferPtr, csLen);
			mBufferPtr += csLen;
			// retrieve the colorspace raw base qualities
			al.CsBaseQualities.insert(0, (const char*)mBufferPtr, csLen);
			mBufferPtr += csLen;
		}

		if ( !al.IsJunk ) {

			// get the number of mismatches
			memcpy((char*)&al.NumMismatches, mBufferPtr, SIZEOF_SHORT);
			mBufferPtr += SIZEOF_SHORT;

			// get mate pair information
			if(isResolvedAsPair) {

				// get the mate reference sequence start position
				memcpy((char*)&al.MateReferenceBegin, mBufferPtr, SIZEOF_INT);
				mBufferPtr += SIZEOF_INT;

				// get the mate reference sequence end position
				memcpy((char*)&al.MateReferenceEnd, mBufferPtr, SIZEOF_INT);
				mBufferPtr += SIZEOF_INT;

				// get the mate reference sequence index
				memcpy((char*)&al.MateReferenceIndex, mBufferPtr, SIZEOF_INT);
				mBufferPtr += SIZEOF_INT;

			} else {

				al.MateReferenceBegin = 0;
				al.MateReferenceEnd   = 0;
				al.MateReferenceIndex = 0;
			}

			unsigned short pairwiseLength = 0;

			if(isLongRead) {

				// get the pairwise length
				memcpy((char*)&pairwiseLength, mBufferPtr, SIZEOF_SHORT);
				mBufferPtr += SIZEOF_SHORT;

				// get the query begin
				memcpy((char*)&al.QueryBegin, mBufferPtr, SIZEOF_SHORT);
				mBufferPtr += SIZEOF_SHORT;

				// get the query end
				memcpy((char*)&al.QueryEnd, mBufferPtr, SIZEOF_SHORT);
				mBufferPtr += SIZEOF_SHORT;

			} else {

				// get the pairwise length
				pairwiseLength = (unsigned char)*mBufferPtr;
				++mBufferPtr;

				// get the query begin
				al.QueryBegin = (unsigned char)*mBufferPtr;
				++mBufferPtr;

				// get the query end
				al.QueryEnd = (unsigned char)*mBufferPtr;
				++mBufferPtr;
			}

			// retrieve the packed pairwise alignment
			al.Reference.Copy((const char*)mBufferPtr, pairwiseLength);
			mBufferPtr += pairwiseLength;

			// unpack the pairwise query bases
			al.Reference.Unpack(al.Query);

			// get the pairwise query base qualities
			const unsigned short bqLength = al.QueryEnd - al.QueryBegin + 1;
			al.BaseQualities.Copy((const char*)mBufferPtr, bqLength);
			mBufferPtr += bqLength;
		}

		//al.BaseQualities.Increment(33);
		//cout << al.Reference.CData() << endl << al.Query.CData() << endl << al.BaseQualities.CData() << endl;
		
		// read the number of tags present in this alignment
		const unsigned char numTags = (unsigned char)*mBufferPtr;
		++mBufferPtr;

		if(numTags != 0) {
			cout << "ERROR: Tags have not been implemented yet." << endl;
			exit(1);
		}

		// DEBUG
		//cout << "reference index: " << al.ReferenceIndex << ", begin: " << al.ReferenceBegin << ", end: " << al.ReferenceEnd << endl;
		//cout << "query begin: " << al.QueryBegin << ", end: " << al.QueryEnd << ", pairwise length: " << pairwiseLength << endl << endl;
	}

	// reads a new compressed partition (returns false if EOF occurs)
	bool CAlignmentReader::ReadPartition(void) {

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
			cout << "ERROR: Tried to read " << compressedSize << " bytes, but received only " << numBytesRead << " bytes (" << mInputFilename << ") [read the partition: LoadNextRead]" << endl;
			exit(1);
		}

		int result = fastlz_decompress(mCompressionBuffer, compressedSize, mBuffer, mBufferLen);

		if(result == 0) {
			cout << "ERROR: Unable to properly uncompress the current data partition." << endl;
			exit(1);
		}

		// set the buffer pointer
		mBufferPtr = mBuffer;

		return true;
	}

	// reads the tag from disk
	void CAlignmentReader::ReadTag(Tag& tag) {

		// read the tag ID
		tag.ID = fgetc(mInStream);

		// read the tag type
		TagType tagType = fgetc(mInStream);
		tag.Type = tagType;

		// read the data
		unsigned short stringLength;

		switch(tagType) {
			case TT_CHAR:
				tag.Char = fgetc(mInStream);
				break;
			case TT_DOUBLE:
				fread((char*)&tag.Double, SIZEOF_DOUBLE, 1, mInStream);
				break;
			case TT_FLOAT:
				fread((char*)&tag.Float, SIZEOF_FLOAT, 1, mInStream);
				break;
			case TT_INT16:
				fread((char*)&tag.Int16, SIZEOF_SHORT, 1, mInStream);
				break;
			case TT_INT32:
				fread((char*)&tag.Int32, SIZEOF_INT, 1, mInStream);
				break;
			case TT_INT64:
				fread((char*)&tag.Int64, SIZEOF_UINT64, 1, mInStream);
				break;
			case TT_STRING:
				fread((char*)&stringLength, SIZEOF_SHORT, 1, mInStream);
				fread(tag.String, stringLength, 1, mInStream);
				break;
			case TT_UCHAR:
				tag.UChar = fgetc(mInStream);
				break;
			case TT_UINT16:
				fread((char*)&tag.UInt16, SIZEOF_SHORT, 1, mInStream);
				break;
			case TT_UINT32:
				fread((char*)&tag.UInt32, SIZEOF_INT, 1, mInStream);
				break;
			case TT_UINT64:
				fread((char*)&tag.UInt64, SIZEOF_UINT64, 1, mInStream);
				break;
			default:
				printf("ERROR: Unknown tag storage type found: %c [%u]\n", tagType, tagType);
				exit(1);
		}
	}

	// sets the file pointer to the beginning of the read data
	void CAlignmentReader::Rewind(void) {
		fseek64(mInStream, mReadsOffset, SEEK_SET);
		mCurrentRead      = 0;
		mPartitionMembers = 0;
		mPartitionSize    = 0;
	}
}
