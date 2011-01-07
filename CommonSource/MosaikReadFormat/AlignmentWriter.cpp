// ***************************************************************************
// CAlignmentWriter - stores reads in a MOSAIK alignment archive.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "AlignmentWriter.h"

namespace MosaikReadFormat {

	// constructor
	CAlignmentWriter::CAlignmentWriter(void)
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
		, mpRefGapVector(NULL)
		, mStatus(AS_UNKNOWN)
		, mIsPairedEndArchive(false)
		, mLastReferenceIndex(0)
		, mLastReferencePosition(0)
		, mStoreIndex(false)
		, MosaikSignature(NULL)
	{
		// set the buffer threshold
		mBufferThreshold = mBufferLen - MEMORY_BUFFER_SIZE;

		// initialize the read and index buffer
		try {
			mBuffer = new unsigned char[mBufferLen];
		} catch(bad_alloc) {
			cout << "ERROR: Unable to allocate enough memory for the alignment output buffer." << endl;
			exit(1);
		}
	}

	// destructor
	CAlignmentWriter::~CAlignmentWriter(void) {
		if(mIsOpen)            Close();
		if(mBuffer)            delete [] mBuffer;
		if(mCompressionBuffer) delete mCompressionBuffer;
		if(MosaikSignature)    delete [] MosaikSignature;
	}

	// adds a header tag
	void CAlignmentWriter::AddHeaderTag(const Tag& tag) {

		// we can only add header tags before opening the file
		if(mIsOpen) {
			printf("ERROR: Header tags can only be added before opening the alignment archive.\n");
			exit(1);
		}

		map<unsigned char, Tag>::iterator htIter = mHeaderTags.find(tag.ID);
		if(htIter == mHeaderTags.end()) {
			mHeaderTags[tag.ID] = tag;
		} else htIter->second = tag;
	}

	// let the references be increased by 1
	void CAlignmentWriter::AdjustSolidReferenceBases(void) {
		// we can only adjust references when the file is open
		if(!mIsOpen)
			printf("Warning: Unable to adjust the numbers of references bases before opeing the alignment archive.");
		else
			for ( vector<ReferenceSequence>::iterator rsIte = mReferenceSequences.begin(); rsIte != mReferenceSequences.end(); rsIte++ )
				rsIte->NumBases++;
		
	}

	// checks the buffer
	void CAlignmentWriter::AdjustBuffer(void) {

		// allocate a new buffer
		unsigned int newBufferLen = mBufferLen << 1;
		unsigned char* newBuffer = NULL;

		try {
			newBuffer = new unsigned char[newBufferLen];
		} catch(bad_alloc) {
			cout << "ERROR: Unable to reallocate enough memory for the alignment output buffer." << endl;
			exit(1);
		}

		// copy the old data and destroy the old buffer
		memcpy(newBuffer, mBuffer, mBufferLen);
		if ( mBuffer ) delete [] mBuffer;

		// repoint the new buffer
		mBuffer          = newBuffer;
		mBufferLen       = newBufferLen;
		mBufferThreshold = newBufferLen - MEMORY_BUFFER_SIZE;
	}

	// adjust the size of partition; the default is 20000
	void CAlignmentWriter::AdjustPartitionSize(unsigned short size) {
                mPartitionSize = size;
	}

	// closes the alignment archive
	void CAlignmentWriter::Close(void) {

		// the archive is not open.
		if ( !mIsOpen )
			return;
		
		// prevent the archive from being updated elsewhere
		mIsOpen = false;

		// flush the buffer
		if(mPartitionMembers > 0) WritePartition();

		// =======================================
		// save the reference sequence information
		// =======================================

		const off_type referenceOffset = ftell64(mOutStream);

		// write the reference sequence dictionary
		vector<ReferenceSequence>::const_iterator rsIter;
		for(rsIter = mReferenceSequences.begin(); rsIter != mReferenceSequences.end(); rsIter++) {

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

			// get the string lengths
			unsigned int nameLen             = (unsigned int)rsIter->Name.size();
			unsigned int speciesLen          = (unsigned int)rsIter->Species.size();
			unsigned int genomeAssemblyIDLen = (unsigned int)rsIter->GenomeAssemblyID.size();
			unsigned int uriLen              = (unsigned int)rsIter->URI.size();

			// write the name length
			fputc((unsigned char)nameLen, mOutStream);

			// write the species length
			fputc((unsigned char)speciesLen, mOutStream);

			// write the genome assembly id length
			fputc((unsigned char)genomeAssemblyIDLen, mOutStream);

			// write the URI length
			fputc((unsigned char)uriLen, mOutStream);

			// write the number of bases
			fwrite((char*)&rsIter->NumBases, SIZEOF_INT, 1, mOutStream);

			// write the number of aligned reads
			fwrite((char*)&rsIter->NumAligned, SIZEOF_UINT64, 1, mOutStream);

			// write the MD5 checksum
			fwrite(rsIter->MD5.data(), 32, 1, mOutStream);

			// write the reference name
			fwrite(rsIter->Name.data(), nameLen, 1, mOutStream);

			// write the species name
			if(speciesLen > 0) fwrite(rsIter->Species.data(), speciesLen, 1, mOutStream);

			// write the genome assembly ID
			if(genomeAssemblyIDLen > 0) fwrite(rsIter->GenomeAssemblyID.data(), genomeAssemblyIDLen, 1, mOutStream);

			// write the URI
			if(uriLen > 0) fwrite(rsIter->URI.data(), uriLen, 1, mOutStream);

			// write the number of reference sequence tags (hard coded as 0 for now)
			fputc(0, mOutStream);
		}

		// ================================
		// save the reference sequence gaps
		// ================================

		CFastLZIO fio;
		off_type referenceGapFileOffset = 0;
		if(mpRefGapVector) {

			// get the current file offset
			referenceGapFileOffset = ftell64(mOutStream);

			// calculate how much buffer space we'll need
			unsigned int requestedBytes = 0;
			for(unsigned int i = 0; i < mNumRefSeqs; i++) requestedBytes += SIZEOF_INT + (unsigned int)mpRefGapVector->at(i).size() * (SIZEOF_INT + SIZEOF_SHORT);
			CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, requestedBytes);

			// create a sorting vector
			vector<GapInfo> refGaps;
			vector<GapInfo>::iterator refGapsIter;
			unordered_map<unsigned int, unsigned short>::const_iterator refHashMapIter;

			// fill our buffer
			unsigned int bufferOffset = 0;
			vector<unordered_map<unsigned int, unsigned short> >::const_iterator rgIter;
			for(rgIter = mpRefGapVector->begin(); rgIter != mpRefGapVector->end(); rgIter++) {

				// store the number of reference gaps
				const unsigned int numReferenceGaps = (unsigned int)rgIter->size();
				memcpy(mBuffer + bufferOffset, (char*)&numReferenceGaps, SIZEOF_INT);
				bufferOffset += SIZEOF_INT;

				// store the reference gap data if available
				if(numReferenceGaps > 0) {

					// populate the sorting vector
					refGaps.resize(numReferenceGaps);

					refGapsIter = refGaps.begin();
					for(refHashMapIter = rgIter->begin(); refHashMapIter != rgIter->end(); refHashMapIter++, refGapsIter++) {
						refGapsIter->Position = refHashMapIter->first;
						refGapsIter->Length   = refHashMapIter->second;
					}

					// sort the reference gaps
					sort(refGaps.begin(), refGaps.end());

					// store the gaps for this reference sequence
					for(refGapsIter = refGaps.begin(); refGapsIter != refGaps.end(); refGapsIter++) {

						// store the reference gap position
						memcpy(mBuffer + bufferOffset, (char*)&refGapsIter->Position, SIZEOF_INT);
						bufferOffset += SIZEOF_INT;

						// store the reference gap length
						memcpy(mBuffer + bufferOffset, (char*)&refGapsIter->Length, SIZEOF_SHORT);
						bufferOffset += SIZEOF_SHORT;
					}
				}
			}

			// write the reference sequence gaps
			fio.Write((char*)mBuffer, bufferOffset, mOutStream);
		}

		// ==============
		// save the index
		// ==============

		off_type indexFileOffset = 0;
		if(mStoreIndex) {

			// get the current file offset
			indexFileOffset = ftell64(mOutStream);

			// calculate how much buffer space we'll need
			const unsigned int numIndexEntries = (unsigned int)mIndex.size();
			const unsigned int requestedBytes = numIndexEntries * (SIZEOF_UINT64 + 2 * SIZEOF_INT);
			CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, requestedBytes);

			// store the number of index entries
			fwrite((char*)&numIndexEntries, SIZEOF_INT, 1, mOutStream);

			unsigned int bufferOffset = 0;
			vector<IndexEntry>::const_iterator indexIter;
			for(indexIter = mIndex.begin(); indexIter != mIndex.end(); indexIter++) {

				// store the reference index
				memcpy(mBuffer + bufferOffset, (char*)&indexIter->ReferenceIndex, SIZEOF_INT);
				bufferOffset += SIZEOF_INT;

				// store the reference position
				memcpy(mBuffer + bufferOffset, (char*)&indexIter->Position, SIZEOF_INT);
				bufferOffset += SIZEOF_INT;

				// store the file offset
				memcpy(mBuffer + bufferOffset, (char*)&indexIter->Offset, SIZEOF_UINT64);
				bufferOffset += SIZEOF_UINT64;
			}

			// write the index
			fio.Write((char*)mBuffer, bufferOffset, mOutStream);
		}

		// =================
		// update the header
		// =================

		// update the number of reads in the archive
		fseek64(mOutStream, NUM_READS_OFFSET, SEEK_SET);
		fwrite((char*)&mNumReads, SIZEOF_UINT64, 1, mOutStream);

		// update the number of bases in the archive
		fwrite((char*)&mNumBases, SIZEOF_UINT64, 1, mOutStream);

		// update the references offset in the archive
		fwrite((char*)&referenceOffset, SIZEOF_UINT64, 1, mOutStream);

		// update the reference gap offset in the archive
		fwrite((char*)&referenceGapFileOffset, SIZEOF_UINT64, 1, mOutStream);

		// update the index offset in the archive		
		fwrite((char*)&indexFileOffset, SIZEOF_UINT64, 1, mOutStream);

		// DEBUG
		//cout << endl;
		//cout << "mNumReads:              " << mNumReads << endl;
		//cout << "mNumBases:              " << mNumBases << endl;
		//cout << "referenceOffset:        " << referenceOffset << endl;
		//cout << "referenceGapFileOffset: " << referenceGapFileOffset << endl;
		//cout << "indexFileOffset:        " << indexFileOffset << endl << endl;

		// close the file stream
		fclose(mOutStream);
	}

	// retrieves the number of bases written
	uint64_t CAlignmentWriter::GetNumBases(void) const {
		return mNumBases;
	}

	// retrieves the number of reads written
	uint64_t CAlignmentWriter::GetNumReads(void) const {
		return mNumReads;
	}

	// opens the alignment archive
	void CAlignmentWriter::Open(const string& filename, const vector<ReferenceSequence>& referenceSequences, const vector<ReadGroup>& readGroups, const AlignmentStatus as, const string& signature) {

		if(mIsOpen) {
			cout << "ERROR: An attempt was made to open an already open alignment archive." << endl;
			exit(1);
		}

		mOutputFilename = filename;

		if(fopen_s(&mOutStream, filename.c_str(), "wb") != 0) {
			cout << "ERROR: Could not open the compressed alignment archive (" << mOutputFilename << ") for writing." << endl;
			exit(1);
		}

		mIsOpen = true;

		// initialization
		mBufferPosition   = 0;
		mPartitionMembers = 0;
		mStatus           = as;

		// copy the reference sequence statistics
		// N.B. we copy these because g++ was making shallow copies before, this is only a temporary fix
		mNumRefSeqs = (unsigned int)referenceSequences.size();
		mReferenceSequences.resize(mNumRefSeqs);

		vector<ReferenceSequence>::const_iterator rsIter = referenceSequences.begin();
		vector<ReferenceSequence>::iterator crsIter;
		for(crsIter = mReferenceSequences.begin(); crsIter != mReferenceSequences.end(); crsIter++, rsIter++) {
			crsIter->GenomeAssemblyID = rsIter->GenomeAssemblyID;
			crsIter->MD5              = rsIter->MD5;
			crsIter->Name             = rsIter->Name;
			crsIter->NumAligned       = 0;
			crsIter->NumBases         = rsIter->NumBases;
			crsIter->Species          = rsIter->Species;
			crsIter->URI              = rsIter->URI;
		}

		// create our composite sequencing technology
		unsigned int numReadGroups = (unsigned int)readGroups.size();
		SequencingTechnologies st = ST_UNKNOWN;

		vector<ReadGroup>::const_iterator rgIter;
		for(rgIter = readGroups.begin(); rgIter != readGroups.end(); rgIter++) st |= rgIter->SequencingTechnology;

		// ================
		// write the header
		// ================

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

		// NB: the following blocks occur at the end of the file
		// REFERENCE_SEQS[*]
		// REFERENCE_GAPS[*]
		// INDEX[*]

		// write the MOSAIK signature
		// the definations of MOSAIK signatures are in Mosaik.h
		//const unsigned char SIGNATURE_LENGTH = 6;
		//const char* MOSAIK_SIGNATURE = "MSKAA\4";
		//fwrite(MOSAIK_SIGNATURE, SIGNATURE_LENGTH, 1, mOutStream);
		if ( ( signature != ALIGNER_SIGNATURE ) && ( signature != ALIGNER_SIGNATURE5 ) && ( signature != SORT_SIGNATURE ) ) {
			cout << "ERROR: The signature for MOSAIK archive is invalid." << endl;
			exit(1);
		}
		fwrite( signature.c_str(), SIGNATURE_LENGTH, 1, mOutStream );

		// write the alignment status
		fputc((unsigned char)as, mOutStream);
		if((mStatus & AS_SORTED_ALIGNMENT) != 0) mStoreIndex         = true;
		if((mStatus & AS_PAIRED_END_READ)  != 0) mIsPairedEndArchive = true;

		// write the sequencing technology
		fwrite((char*)&st, SIZEOF_SHORT, 1, mOutStream);

		// write the archive date
		uint64_t currentTime = CTimeSupport::GetSystemTime();
		fwrite((char*)&currentTime, SIZEOF_UINT64, 1, mOutStream);

		// write the number of reference sequences
		fwrite((char*)&mNumRefSeqs, SIZEOF_INT, 1, mOutStream);

		// write the number of read groups
		fwrite((char*)&numReadGroups, SIZEOF_INT, 1, mOutStream);

		// skip the # of reads, # of bases, references offset, reference gap offset, index offset
		fseek64(mOutStream, 2 * SIZEOF_UINT64 + 3 * SIZEOF_OFF_TYPE, SEEK_CUR);

		// write the number of header tags (hard coded as 0 for now)
		const unsigned char numHeaderTags = (unsigned char)mHeaderTags.size();
		fputc(numHeaderTags, mOutStream);

		// write the header tags
		if(numHeaderTags != 0) {
			map<unsigned char, Tag>::const_iterator htIter;
			for(htIter = mHeaderTags.begin(); htIter != mHeaderTags.end(); htIter++) {
				WriteTag(htIter);
			}
		}

		// write the read groups
		for(rgIter = readGroups.begin(); rgIter != readGroups.end(); rgIter++) {

			// write the metadata string lengths
			const unsigned char centerNameLen   = (unsigned char)rgIter->CenterName.size();
			const unsigned char libraryNameLen  = (unsigned char)rgIter->LibraryName.size();
			const unsigned char platformUnitLen = (unsigned char)rgIter->PlatformUnit.size();
			const unsigned char readGroupIDLen  = (unsigned char)rgIter->ReadGroupID.size();
			const unsigned char sampleNameLen   = (unsigned char)rgIter->SampleName.size();
			const unsigned short descriptionLen = (unsigned short)rgIter->Description.size();

			fputc(centerNameLen,   mOutStream);
			fputc(libraryNameLen,  mOutStream);
			fputc(platformUnitLen, mOutStream);
			fputc(readGroupIDLen,  mOutStream);
			fputc(sampleNameLen,   mOutStream);
			fwrite((char*)&descriptionLen, SIZEOF_SHORT, 1, mOutStream);
			fwrite((char*)&rgIter->SequencingTechnology, SIZEOF_SHORT, 1, mOutStream);
			fwrite((char*)&rgIter->MedianFragmentLength, SIZEOF_INT, 1, mOutStream);

			// write the metadata strings
			fwrite(rgIter->CenterName.c_str(),   centerNameLen,   1, mOutStream);
			fwrite(rgIter->Description.c_str(),  descriptionLen,  1, mOutStream);
			fwrite(rgIter->LibraryName.c_str(),  libraryNameLen,  1, mOutStream);
			fwrite(rgIter->PlatformUnit.c_str(), platformUnitLen, 1, mOutStream);
			fwrite(rgIter->ReadGroupID.c_str(),  readGroupIDLen,  1, mOutStream);
			fwrite(rgIter->SampleName.c_str(),   sampleNameLen,   1, mOutStream);

			// write the number of read group tags (hard coded as 0 for now)
			fputc(0, mOutStream);
		}
	}

	// saves the read to the alignment archive
	void CAlignmentWriter::SaveAlignedRead(const Mosaik::AlignedRead& ar ) {

		// check the memory buffer
		if(mBufferPosition > mBufferThreshold) AdjustBuffer();

		// initialize
		const unsigned int numMate1Alignments = (unsigned int)ar.Mate1Alignments.size();
		const unsigned int numMate2Alignments = (unsigned int)ar.Mate2Alignments.size();

		const bool haveMate1 = (numMate1Alignments != 0 ? true : false);
		const bool haveMate2 = (numMate2Alignments != 0 ? true : false);

		unsigned int numMate1OriginalAlignments = 0;
		unsigned int numMate2OriginalAlignments = 0;

		if ( haveMate1 )
			numMate1OriginalAlignments = ar.Mate1Alignments[0].NumMapped;
		if ( haveMate2 )
			numMate2OriginalAlignments = ar.Mate2Alignments[0].NumMapped;
		
		// check if this is a long read
		const bool isLongRead = ar.IsLongRead;

		// derive our read status
		unsigned char readStatus = RF_UNKNOWN;

		if(haveMate1)      readStatus |= RF_HAVE_MATE1;
		if(haveMate2)      readStatus |= RF_HAVE_MATE2;
		if(isLongRead)     readStatus |= RF_IS_LONG_READ;
		if(ar.IsPairedEnd) readStatus |= RF_IS_PAIRED_IN_SEQUENCING;
		if(ar.IsResolvedAsPair) readStatus |= RF_RESOLVED_AS_PAIR;

		// write the read header
		WriteReadHeader(ar.Name, ar.ReadGroupCode, readStatus, numMate1Alignments, numMate2Alignments, numMate1OriginalAlignments, numMate2OriginalAlignments );

		// ===============================
		// serialize each mate 1 alignment
		// ===============================

		vector<Alignment>::const_iterator alIter;
		const Alignment *pAl = NULL, *pAlBegin = NULL;

		if(haveMate1) {
			pAlBegin = &ar.Mate1Alignments[0];
			for(alIter = ar.Mate1Alignments.begin(); alIter != ar.Mate1Alignments.end(); ++alIter) {
				pAl = pAlBegin + (alIter - ar.Mate1Alignments.begin());
				WriteAlignment(pAl, isLongRead, ar.IsPairedEnd, alIter->IsFirstMate, ar.IsResolvedAsPair);
				//WriteAlignment(pAl, isLongRead, ar.IsPairedEnd, true, false);
			}
		}

		// ===============================
		// serialize each mate 2 alignment
		// ===============================

		if(haveMate2) {
			pAlBegin = &ar.Mate2Alignments[0];
			for(alIter = ar.Mate2Alignments.begin(); alIter != ar.Mate2Alignments.end(); ++alIter) {
				pAl = pAlBegin + (alIter - ar.Mate2Alignments.begin());
				WriteAlignment(pAl, isLongRead, ar.IsPairedEnd, alIter->IsFirstMate, ar.IsResolvedAsPair);
				//WriteAlignment(pAl, isLongRead, ar.IsPairedEnd, false, false);
			}
		}

		// flush the buffer
		mPartitionMembers++;
		if(mPartitionMembers >= mPartitionSize) WritePartition();

		// increment the read counter
		mNumReads++;
	}

	// saves the alignment to the alignment archive
	void CAlignmentWriter::SaveAlignment(Alignment* pAl) {

		// DEBUG
		//printf("Name: %s, reference: %2u, orientation: %c, alignment quality: %2u\n", pAl->Name.CData(), pAl->ReferenceIndex, (pAl->IsReverseStrand ? 'R' : 'F'), pAl->Quality);
		//printf("%8u %s %u\n", pAl->ReferenceBegin, pAl->Anchor.CData(), pAl->ReferenceEnd);
		//printf("%8u %s %u\n", pAl->QueryBegin, pAl->Query.CData(), pAl->QueryEnd);
		//printf("         %s\n\n", pAl->BaseQualities.CData());

		// check the memory buffer
		if(mBufferPosition > mBufferThreshold) AdjustBuffer();

		// check if this is a long read
		const unsigned short pairwiseLength = (unsigned short)pAl->Reference.Length();

		bool isLongRead = false;
		if((pAl->QueryEnd > 255) || (pairwiseLength > 255)) isLongRead = true;

		// derive our read status
		unsigned char readStatus = RF_HAVE_MATE1;

		if(isLongRead)            readStatus |= RF_IS_LONG_READ;
		if(pAl->IsPairedEnd)      readStatus |= RF_IS_PAIRED_IN_SEQUENCING;
		if(pAl->IsResolvedAsPair) readStatus |= RF_RESOLVED_AS_PAIR;

		// write the read header
		// has one first mate and doesn't have any second mate
		WriteReadHeader(pAl->Name, pAl->ReadGroupCode, readStatus, 1, 0);

		// =======================
		// serialize our alignment
		// =======================

		WriteAlignment(pAl, isLongRead, pAl->IsPairedEnd, pAl->IsFirstMate, pAl->IsResolvedAsPair);

		// flush the buffer
		mPartitionMembers++;
		if(mPartitionMembers >= mPartitionSize) WritePartition();

		// increment the read counter
		mNumReads++;
	}

	// saves the paired-end read to the alignment archive
	void CAlignmentWriter::SaveRead(const Mosaik::Read& mr, CNaiveAlignmentSet& mate1Alignments, CNaiveAlignmentSet& mate2Alignments) {

		// check the memory buffer
		if(mBufferPosition > mBufferThreshold) AdjustBuffer();

		// initialize
		const unsigned int numMate1Alignments = mate1Alignments.GetCount();
		const unsigned int numMate2Alignments = mate2Alignments.GetCount();

		const bool haveMate1 = (numMate1Alignments != 0 ? true : false);
		const bool haveMate2 = (numMate2Alignments != 0 ? true : false);

		// check if this is a long read
		bool isLongRead = false;
		if(mate1Alignments.HasLongAlignment() || mate2Alignments.HasLongAlignment()) isLongRead = true;

		// derive our read status
		unsigned char readStatus = RF_UNKNOWN;

		if(haveMate1)           readStatus |= RF_HAVE_MATE1;
		if(haveMate2)           readStatus |= RF_HAVE_MATE2;
		if(isLongRead)          readStatus |= RF_IS_LONG_READ;
		if(mIsPairedEndArchive) readStatus |= RF_IS_PAIRED_IN_SEQUENCING;

		// write the read header
		WriteReadHeader(mr.Name, mr.ReadGroupCode, readStatus, numMate1Alignments, numMate2Alignments);

		// ===============================
		// serialize each mate 1 alignment
		// ===============================

		if(haveMate1) {
			AlignmentSet* pMateSet = mate1Alignments.GetSet();
			for(AlignmentSet::iterator alIter = pMateSet->begin(); alIter != pMateSet->end(); ++alIter) {
				WriteAlignment(&(*alIter), isLongRead, mIsPairedEndArchive, true, false);
			}
		}

		// ===============================
		// serialize each mate 2 alignment
		// ===============================

		if(haveMate2) {
			AlignmentSet* pMateSet = mate2Alignments.GetSet();
			for(AlignmentSet::iterator alIter = pMateSet->begin(); alIter != pMateSet->end(); ++alIter) {
				WriteAlignment(&(*alIter), isLongRead, mIsPairedEndArchive, false, false);
			}
		}

		// flush the buffer
		mPartitionMembers++;
		if(mPartitionMembers >= mPartitionSize) WritePartition();

		// increment the read counter
		mNumReads++;
	}

	// saves the paired-end read to the alignment archive
	void CAlignmentWriter::SaveRead(
		const Mosaik::Read& mr, 
		const Alignment & mate1Alignment, 
		const Alignment & mate2Alignment, 
		const bool & isLongRead,
		const bool & isSaveMate1,
		const bool & isSaveMate2) {

		// check the memory buffer
		if(mBufferPosition > mBufferThreshold) AdjustBuffer();

		// initialize
		const unsigned int numMate1Alignments = isSaveMate1 ? 1 : 0;
		const unsigned int numMate2Alignments = isSaveMate2 ? 1 : 0;
		unsigned int numMate1OriginalAlignments = 0;
		unsigned int numMate2OriginalAlignments = 0;

		const bool haveMate1 = isSaveMate1;
		const bool haveMate2 = isSaveMate2;

		if ( haveMate1 )
			numMate1OriginalAlignments = mate1Alignment.NumMapped;
		if ( haveMate2 )
			numMate2OriginalAlignments = mate2Alignment.NumMapped;

		// check if this is a long read
		//bool isLongRead = false;
		//if(mate1Alignments.HasLongAlignment() || mate2Alignments.HasLongAlignment()) isLongRead = true;

		// derive our read status
		unsigned char readStatus = RF_UNKNOWN;

		if(haveMate1)           readStatus |= RF_HAVE_MATE1;
		if(haveMate2)           readStatus |= RF_HAVE_MATE2;
		if(isLongRead)          readStatus |= RF_IS_LONG_READ;
		if(mIsPairedEndArchive) readStatus |= RF_IS_PAIRED_IN_SEQUENCING;

		// write the read header
		WriteReadHeader( mr.Name, mr.ReadGroupCode, readStatus, numMate1Alignments, numMate2Alignments, numMate1OriginalAlignments, numMate2OriginalAlignments );

		// ===============================
		// serialize each mate 1 alignment
		// ===============================

		if(haveMate1) {
			//AlignmentSet* pMateSet = mate1Alignments.GetSet();
			//for(vector<Alignment>::const_iterator alIter = mate1Alignments.begin(); alIter != mate1Alignments.end(); ++alIter) {
			//	WriteAlignment(&(*alIter), isLongRead, mIsPairedEndArchive, true, false);
			//}
			WriteAlignment( &mate1Alignment, isLongRead, mIsPairedEndArchive, true, false );
		}

		// ===============================
		// serialize each mate 2 alignment
		// ===============================

		if(haveMate2) {
			//AlignmentSet* pMateSet = mate2Alignments.GetSet();
			//for(vector<Alignment>::const_iterator alIter = mate2Alignments.begin(); alIter != mate2Alignments.end(); ++alIter) {
			//	WriteAlignment(&(*alIter), isLongRead, mIsPairedEndArchive, false, false);
			//}
			WriteAlignment( &mate2Alignment, isLongRead, mIsPairedEndArchive, false, false );
		}

		// flush the buffer
		mPartitionMembers++;
		if(mPartitionMembers >= mPartitionSize) WritePartition();

		// increment the read counter
		mNumReads++;
	}

	// serializes the specified alignment
	void CAlignmentWriter::WriteAlignment(const Alignment* pAl, const bool isLongRead, const bool isPairedEnd, const bool isFirstMate, const bool isResolvedAsPair) {

		// check the memory buffer
		if(mBufferPosition > mBufferThreshold) AdjustBuffer();

		// store the reference sequence start position
		mLastReferencePosition = pAl->ReferenceBegin;
		memcpy(mBuffer + mBufferPosition, (char*)&pAl->ReferenceBegin, SIZEOF_INT);
		mBufferPosition += SIZEOF_INT;

		// store the reference sequence end position
		memcpy(mBuffer + mBufferPosition, (char*)&pAl->ReferenceEnd, SIZEOF_INT);
		mBufferPosition += SIZEOF_INT;

		// store the reference sequence index
		if ( pAl->ReferenceIndex > mNumRefSeqs ) {
			cout << "ERROR: The reference index is out the range when saving alignments." << endl;
			exit(1);
		}
		mLastReferenceIndex = pAl->ReferenceIndex;
		mReferenceSequences[pAl->ReferenceIndex].NumAligned++;
		memcpy(mBuffer + mBufferPosition, (char*)&pAl->ReferenceIndex, SIZEOF_INT);
		mBufferPosition += SIZEOF_INT;

		// store the alignment quality
		mBuffer[mBufferPosition++] = pAl->Quality;

		// store the alignment status flag
		unsigned char status = AF_UNKNOWN;

		if(isPairedEnd && isFirstMate)                   status |= AF_IS_FIRST_MATE;
		if(isPairedEnd && !isFirstMate)                  status |= AF_IS_SECOND_MATE;
		if(pAl->IsReverseStrand)                         status |= AF_IS_REVERSE_STRAND;
		if(isResolvedAsPair && pAl->IsMateReverseStrand) status |= AF_IS_MATE_REVERSE_STRAND;
		if(pAl->WasRescued)                              status |= AF_WAS_RESCUED;
		if(!pAl->IsMapped)                               status |= AF_IS_UNMAPPED;
		//if(pAl->IsJunk)                                  status |= AF_IS_JUNK;


		// not really sure how this applies to single-end and paired-end reads. Disabling the flag for now.
		//if(!isPrimaryAlignment)                              status |= AF_IS_NOT_PRIMARY;

		mBuffer[mBufferPosition++] = status;

		//if ( pAl->IsJunk )
		//	mBuffer[mBufferPosition++] = 0;
		//else {	
		
		// store the number of mismatches
		memcpy(mBuffer + mBufferPosition, (char*)&pAl->NumMismatches, SIZEOF_SHORT);
		mBufferPosition += SIZEOF_SHORT;

		// get the pairwise length
		const unsigned short pairwiseLength = (unsigned short)pAl->Reference.Length();

		// add mate pair information
		if(isResolvedAsPair) {

			// store the mate reference sequence start position
			memcpy(mBuffer + mBufferPosition, (char*)&pAl->MateReferenceBegin, SIZEOF_INT);
			mBufferPosition += SIZEOF_INT;

			// store the mate reference sequence end position
			memcpy(mBuffer + mBufferPosition, (char*)&pAl->MateReferenceEnd, SIZEOF_INT);
			mBufferPosition += SIZEOF_INT;

			// store the mate reference sequence index
			memcpy(mBuffer + mBufferPosition, (char*)&pAl->MateReferenceIndex, SIZEOF_INT);
			mBufferPosition += SIZEOF_INT;
		}

		if(isLongRead) {

			// store the pairwise length
			memcpy(mBuffer + mBufferPosition, (char*)&pairwiseLength, SIZEOF_SHORT);
			mBufferPosition += SIZEOF_SHORT;

			// store the query begin
			memcpy(mBuffer + mBufferPosition, (char*)&pAl->QueryBegin, SIZEOF_SHORT);
			mBufferPosition += SIZEOF_SHORT;

			// store the query end
			memcpy(mBuffer + mBufferPosition, (char*)&pAl->QueryEnd, SIZEOF_SHORT);
			mBufferPosition += SIZEOF_SHORT;

		} else {

			// store the pairwise length
			mBuffer[mBufferPosition++] = (unsigned char)pairwiseLength;

			// store the query begin
			mBuffer[mBufferPosition++] = (unsigned char)pAl->QueryBegin;

			// store the query end
			mBuffer[mBufferPosition++] = (unsigned char)pAl->QueryEnd;
		}

		// pack the pairwise-alignment string
		CMosaikString packString(pAl->Reference);
		packString.Pack(pAl->Query);

		// store the packed pairwise alignment
		memcpy(mBuffer + mBufferPosition, packString.CData(), pairwiseLength);
		mBufferPosition += pairwiseLength;

		// DEBUG
		//printf("\nReference: %s\n", pAl->Reference.CData());
		//printf("Query:     %s\n", pAl->Query.CData());
		////printf("Packed data: ");
		////for(unsigned short k = 0; k < pairwiseLength; ++k) {
		////	printf("%02x ", packString[k]);
		////}
		//const char* pRef = pAl->Reference.CData();
		//const char* pQry = pAl->Query.CData();
		//for(unsigned short k = 0; k < pairwiseLength; ++k) {
		//	if(pRef[k] != pQry[k]) {
		//		printf("ref: %c, query: %c, pack string: %02x\n", pRef[k], pQry[k], packString[k]);
		//	}
		//	//printf("%02x ", packString[k]);
		//}
		//printf("\n\n");
		////exit(1);

		// store the pairwise query base qualities
		const unsigned int bqLength = pAl->BaseQualities.Length();
		memcpy(mBuffer + mBufferPosition, pAl->BaseQualities.CData(), bqLength);
		mBufferPosition += bqLength;

		// write the number of tags present in this alignment (hard coded as 0 for now)
		mBuffer[mBufferPosition++] = 0;

		// update our statistics
		mNumBases += bqLength;
		//}
		
		// check the buffer
		if(mBufferPosition >= mBufferLen) {
			cout << endl << "ERROR: Buffer overrun detected when saving read. Used " << mBufferPosition << " bytes, but allocated " << mBufferLen << " bytes." << endl;
			exit(1);
		}
	}

	// write the read header to disk
	void CAlignmentWriter::WriteReadHeader(
		const CMosaikString& readName, 
		const unsigned int   readGroupCode, 
		const unsigned char  readStatus, 
		const unsigned int   numMate1Alignments, 
		const unsigned int   numMate2Alignments, 
		const unsigned int   numMate1OriginalAlignments,
		const unsigned int   numMate2OriginalAlignments
		) {

		// store the read name
		const unsigned char readNameLen = (unsigned char)readName.Length();
		mBuffer[mBufferPosition++] = readNameLen;
		memcpy(mBuffer + mBufferPosition, readName.CData(), readNameLen);
		mBufferPosition += readNameLen;

		// store the read group code
		memcpy(mBuffer + mBufferPosition, (char*)&readGroupCode, SIZEOF_INT);
		mBufferPosition += SIZEOF_INT;

		// store the read status flag
		mBuffer[mBufferPosition++] = readStatus;

		// store the number of mate 1 alignments
		if(numMate1Alignments != 0) {
			memcpy(mBuffer + mBufferPosition, (char*)&numMate1Alignments, SIZEOF_INT);
			mBufferPosition += SIZEOF_INT;
			memcpy(mBuffer + mBufferPosition, (char*)&numMate1OriginalAlignments, SIZEOF_INT);
			mBufferPosition += SIZEOF_INT;
		}

		// store the number of mate 2 alignments
		if(numMate2Alignments != 0) {
			memcpy(mBuffer + mBufferPosition, (char*)&numMate2Alignments, SIZEOF_INT);
			mBufferPosition += SIZEOF_INT;
			memcpy(mBuffer + mBufferPosition, (char*)&numMate2OriginalAlignments, SIZEOF_INT);
			mBufferPosition += SIZEOF_INT;
		}
	}

	// set reference gaps vector
	void CAlignmentWriter::SetReferenceGaps(vector<unordered_map<unsigned int, unsigned short> >* pRefGapVector) {
		mpRefGapVector = pRefGapVector;
	}

	// write partition to disk
	void CAlignmentWriter::WritePartition(void) {

		// store the partition index entry
		if(mStoreIndex) {
			IndexEntry ie;
			ie.Offset         = ftell64(mOutStream);
			ie.ReferenceIndex = mLastReferenceIndex;
			ie.Position       = mLastReferencePosition;
			mIndex.push_back(ie);
		}

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
		memset(mBuffer, 0, mBufferLen);
	}

	// writes the tag to disk
	void CAlignmentWriter::WriteTag(const map<unsigned char, Tag>::const_iterator& htIter) {

		// localize the tag type
		TagType tagType = htIter->second.Type;

		// write the tag ID
		fputc(htIter->second.ID, mOutStream);

		// write the tag type
		fputc(tagType, mOutStream);

		// write the data
		unsigned short stringLength;

		switch(tagType) {
			case TT_CHAR:
				fputc(htIter->second.Char, mOutStream);
				break;
			case TT_DOUBLE:
				fwrite((char*)&htIter->second.Double, SIZEOF_DOUBLE, 1, mOutStream);
				break;
			case TT_FLOAT:
				fwrite((char*)&htIter->second.Float, SIZEOF_FLOAT, 1, mOutStream);
				break;
			case TT_INT16:
				fwrite((char*)&htIter->second.Int16, SIZEOF_SHORT, 1, mOutStream);
				break;
			case TT_INT32:
				fwrite((char*)&htIter->second.Int32, SIZEOF_INT, 1, mOutStream);
				break;
			case TT_INT64:
				fwrite((char*)&htIter->second.Int64, SIZEOF_UINT64, 1, mOutStream);
				break;
			case TT_STRING:
				stringLength = (unsigned short)strlen(htIter->second.String);
				fwrite((char*)&stringLength, SIZEOF_SHORT, 1, mOutStream);
				fwrite(htIter->second.String, stringLength, 1, mOutStream);
				break;
			case TT_UCHAR:
				fputc(htIter->second.UChar, mOutStream);
				break;
			case TT_UINT16:
				fwrite((char*)&htIter->second.UInt16, SIZEOF_SHORT, 1, mOutStream);
				break;
			case TT_UINT32:
				fwrite((char*)&htIter->second.UInt32, SIZEOF_INT, 1, mOutStream);
				break;
			case TT_UINT64:
				fwrite((char*)&htIter->second.UInt64, SIZEOF_UINT64, 1, mOutStream);
				break;
			default:
				printf("ERROR: Unknown tag storage type found: %c\n", tagType);
				exit(1);
		}
	}
}
