// ***************************************************************************
// CReferenceSequenceReader - loads reference sequences from the MOSAIK 
//                            reference sequence archive.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "ReferenceSequenceReader.h"

namespace MosaikReadFormat {

	// constructor
	CReferenceSequenceReader::CReferenceSequenceReader(void)
		: mIsOpen(false)
		, mInStream(NULL)
		, mConcatenatedOffset(0)
		, mConcatenated2bOffset(0)
		, mIndexOffset(0)
		, mReferenceBasesOffset(0)
		, mMaskedRegionsOffset(0)
		, mNumReferenceSequences(0)
		, mConcatenatedLen(0)
		, mConcatenated2bLen(0)
		, mStatus(0)
	{}

	// destructor
	CReferenceSequenceReader::~CReferenceSequenceReader(void) {
		if(mIsOpen) Close();
	}

	// checks to see if this is truly a MOSAIK reference sequence archive
	bool CReferenceSequenceReader::CheckFile(const string& filename, const bool showError) {

		// read in the first 6 characters
		char signature[7];
		signature[6] = 0;
		bool foundError = false;

		const char* MOSAIK_SIGNATURE = "MSKRS\2";

		// open the MOSAIK reference sequence archive
		FILE* checkStream = fopen(filename.c_str(), "rb");

		if(!checkStream) {
			if(showError) {
				printf("ERROR: Could not open %s when validating the reference sequence archive.\n", filename.c_str());
				exit(1);
			}

			foundError = true;
		}

		// retrieve the MOSAIK reference sequence archive signature
		if(!foundError) {

			// check if we were able to read 6 bytes
			if(fread(signature, 1, 6, checkStream) < 6) {
				if(showError) {
					printf("ERROR: It seems that the input file (%s) is not in the MOSAIK reference sequence archive format.\n", filename.c_str());
					exit(1);
				}

				foundError = true;
			}

			// check if the signatures match
			if(!foundError && (strncmp(signature, MOSAIK_SIGNATURE, 5) != 0)) {
				if(showError) {
					printf("ERROR: It seems that the input file (%s) is not in the MOSAIK reference sequence archive format.\n", filename.c_str());
					exit(1);
				}

				foundError = true;
			}

			// check if the file format is an older version
			if(!foundError && (MOSAIK_SIGNATURE[5] > signature[5])) {
				if(showError) {
					printf("ERROR: It seems that the input file (%s) was created in an older version of MosaikBuild. A new reference sequence archive is required.\n", filename.c_str());
					printf("       file version: %hu, expected version: %hu\n", signature[5], MOSAIK_SIGNATURE[5]);
					exit(1);
				}

				foundError = true;
			}
		}

		// close the file
		if(checkStream) fclose(checkStream);

		// return the appropriate values
		if(foundError) return false;
		return true;
	}

	// closes the reference sequence archive
	void CReferenceSequenceReader::Close(void) {
		if(mIsOpen) {
			mIsOpen = false;
			fclose(mInStream);
		}
	}

	// copies the reference sequences from this archive into the supplied character array
	// NOTE: caller frees the memory
	void CReferenceSequenceReader::CopyReferenceSequences(char** &pSeqs) {

		// get a sorted vector of the reference sequences
		vector<ReferenceSequence> refSeqs;
		GetReferenceSequences(refSeqs);
		pSeqs = new char*[mNumReferenceSequences];

		// grab the appropriate sequences
		string bases;
		vector<ReferenceSequence>::const_iterator rsIter;
		unsigned int currentReferenceIndex = 0;
		for(rsIter = refSeqs.begin(); rsIter != refSeqs.end(); ++rsIter, ++currentReferenceIndex) {

			// retrieve the bases for the current sequence
			GetReferenceSequence(rsIter->Name, bases);
			const unsigned int numBases = (unsigned int)bases.size();

			// copy the bases to our character array
			pSeqs[currentReferenceIndex] = new char[numBases + 1];
			memcpy(pSeqs[currentReferenceIndex], bases.data(), numBases);
			pSeqs[currentReferenceIndex][numBases] = 0;
		}
	}

	// returns the number of reference sequences in this archive
	unsigned int CReferenceSequenceReader::GetNumReferenceSequences(void) const {
		return mNumReferenceSequences;
	}

	// returns the reference sequence length
	unsigned int CReferenceSequenceReader::GetReferenceSequenceLength(void) const {
		return mConcatenatedLen;
	}

	// retrieves the desired reference sequence and places it in the specified string
	void CReferenceSequenceReader::GetReferenceSequence(const string& name, string& bases) {

		// check if the reference sequence is in the index unordered map
		unordered_map<string, ReferenceSequence>::const_iterator indexIter;

		indexIter = mIndex.find(name);
		if(indexIter == mIndex.end()) {
			cout << "ERROR: Unable to find the specified reference sequence (" << name << ") in the index." << endl;
			exit(1);
		}

		// jump to the reference sequence
		fseek64(mInStream, indexIter->second.BasesOffset, SEEK_SET);

		// read the reference sequence
		mFIO.Read(bases, mInStream);
	}

	// adds the reference sequences to the supplied vector
	void CReferenceSequenceReader::GetReferenceSequences(vector<ReferenceSequence>& referenceSequences) {

		// clear the vector
		referenceSequences.clear();

		// add each reference sequence to the vector
		unordered_map<string, ReferenceSequence>::const_iterator indexIter;
		for(indexIter = mIndex.begin(); indexIter != mIndex.end(); indexIter++)		
			referenceSequences.push_back(indexIter->second);

		// sort the reference regions
		sort(referenceSequences.begin(), referenceSequences.end(), SortReferenceSequencesByBeginAsc());
	}

	// returns the reference sequence status
	ReferenceSequenceStatus CReferenceSequenceReader::GetStatus(void) const {
		return mStatus;
	}

	// returns true if the reference sequences in this archive match those from the supplied vector
	bool CReferenceSequenceReader::HasSameReferenceSequences(vector<ReferenceSequence>& otherSeqs) {

		// get a sorted vector of the reference sequences
		vector<ReferenceSequence> refSeqs;
		GetReferenceSequences(refSeqs);

		// check for the same number of reference sequences
		if(refSeqs.size() != otherSeqs.size()) return false;

		// check the names of the references sequences
		bool foundError = false;
		vector<ReferenceSequence>::const_iterator rsIter;
		vector<ReferenceSequence>::const_iterator otherIter = otherSeqs.begin();

		for(rsIter = refSeqs.begin(); rsIter != refSeqs.end(); ++rsIter, ++otherIter) {
			if(rsIter->Name != otherIter->Name) {
				foundError = true;
				break;
			}
		}

		return !foundError;
	}

	// initializes the supplied pointer with the concatenated reference sequence
	void CReferenceSequenceReader::LoadConcatenatedSequence(char* &referenceSequence) {

		// jump to the reference sequence
		fseek64(mInStream, mConcatenatedOffset, SEEK_SET);

		// read the reference sequence
		unsigned int referenceSequenceLen = 0;
		mFIO.Read(referenceSequence, referenceSequenceLen, mInStream);

		// sanity checking
		if(referenceSequenceLen != mConcatenatedLen) {
			printf("ERROR: Found a mismatch between the uncompressed concatenated sequence length and the header concatenated sequence length.\n");
			exit(1);
		}
	}

	//initializes the supplied pointer with the 2-bit concatenated reference sequence
	void CReferenceSequenceReader::Load2BitConcatenatedSequence(char* &referenceSequence, char* &maskSequence, unsigned int& numMaskedPositions) {

		// jump to the reference sequence
		fseek64(mInStream, mConcatenated2bOffset, SEEK_SET);

		// read the reference sequence
		unsigned int concatenated2bLen = 0;
		mFIO.Read(referenceSequence, concatenated2bLen, mInStream);

		// sanity checking
		if(concatenated2bLen != mConcatenated2bLen) {
			printf("ERROR: Found a mismatch between the uncompressed concatenated 2-bit sequence length and the header concatenated 2-bit sequence length.\n");
			exit(1);
		}

		// read the number of masked positions
		fread((char*)&numMaskedPositions, SIZEOF_INT, 1, mInStream);

		// read the masked sequence
		if(numMaskedPositions > 0) {

			const unsigned int expectedMaskedSequenceLength = numMaskedPositions * 2 * SIZEOF_INT;
			unsigned int maskedSequenceLength = 0;
			mFIO.Read(maskSequence, maskedSequenceLength, mInStream);

			// sanity checking
			if(maskedSequenceLength != expectedMaskedSequenceLength) {
				printf("ERROR: Found a mismatch between the uncompressed mask sequence length and the header mask sequence length.\n");
				exit(1);
			}
		}
	}

	// opens the reference sequence archive
	void CReferenceSequenceReader::Open(const string& filename) {

		if(mIsOpen) {
			printf("ERROR: An attempt was made to open an already open reference sequence archive.\n");
			exit(1);
		}

		mInStream = fopen(filename.c_str(), "rb");

		if(!mInStream) {
			printf("ERROR: Could not open the reference sequence archive (%s) for reading.\n", filename.c_str());
			exit(1);
		}

		// ===============
		// read the header
		// ===============

		// MOSAIK_SIGNATURE[6]	       0  -  5
		// STATUS[1]                   6  -  6
		// ARCHIVE_DATE[8]		       7  - 14
		// NUM_REFERENCES[4]           15 - 18
		// CONCATENATED_LENGTH[4]      19 - 22
		// CONCATENATED_OFFSET[8]      23 - 30
		// CONCATENATED_2BIT_LENGTH[4] 31 - 34
		// CONCATENATED_2BIT_OFFSET[8] 35 - 42
		// INDEX_OFFSET[8]             43 - 50
		// REFERENCE_BASES_OFFSET[8]   51 - 58
		// MASKED_REGIONS_OFFSET[8]    59 - 66
		// RESERVED[8]                 67 - 74

		// skip the signature
		fseek64(mInStream, 6, SEEK_SET);

		// read the status
		mStatus = fgetc(mInStream);

		// skip the archive date
		fseek64(mInStream, SIZEOF_UINT64, SEEK_CUR);

		// read the number of reference sequences
		fread((char*)&mNumReferenceSequences, SIZEOF_INT, 1, mInStream);

		// read the concatenated reference sequence length
		fread((char*)&mConcatenatedLen, SIZEOF_INT, 1, mInStream);

		// read the concatenated reference offset
		fread((char*)&mConcatenatedOffset, SIZEOF_OFF_TYPE, 1, mInStream);

		// read the concatenated 2-bit reference sequence length
		fread((char*)&mConcatenated2bLen, SIZEOF_INT, 1, mInStream);

		// read the the concatenated 2-bit reference offset
		fread((char*)&mConcatenated2bOffset, SIZEOF_OFF_TYPE, 1, mInStream);

		// read the index offset
		fread((char*)&mIndexOffset, SIZEOF_OFF_TYPE, 1, mInStream);

		// read the reference bases offset
		fread((char*)&mReferenceBasesOffset, SIZEOF_OFF_TYPE, 1, mInStream);

		// read the masked regions offset
		fread((char*)&mMaskedRegionsOffset, SIZEOF_OFF_TYPE, 1, mInStream);

		// skip the reserved field
		fseek64(mInStream, SIZEOF_UINT64, SEEK_CUR);

		// ==============
		// read the index
		// ==============

		// clear the unordered map
		mIndex.clear();

		// jump to the beginning of the index
		fseek64(mInStream, mIndexOffset, SEEK_SET);

		// read the index
		for(unsigned int i = 0; i < mNumReferenceSequences; i++) {

			// REFERENCE_SEQ_NAME_LEN[1]                0 -  0 
			// REFERENCE_SEQ_SPECIES_LEN[1]             1 -  1
			// REFERENCE_SEQ_GENOME_ASSEMBLY_ID_LEN[1]  2 -  2
			// REFERENCE_SEQ_URI_LEN[1]                 3 -  3
			// REFERENCE_SEQ_NUM_BASES[4]               4 -  7
			// REFERENCE_SEQ_BEGIN[4]                   8 - 11
			// REFERENCE_SEQ_END[4]                    12 - 15
			// REFERENCE_SEQ_SEQ_OFFSET[8]             16 - 23
			// REFERENCE_SEQ_MD5[16]                   24 - 39
			// REFERENCE_SEQ_NAME[X]                   40 - XX
			// REFERENCE_SEQ_SPECIES[X]
			// REFERENCE_SEQ_GENOME_ASSEMBLY_ID[X]
			// REFERENCE_SEQ_URI[X]

			ReferenceSequence rs;

			// read the name length
			const unsigned char nameLen = fgetc(mInStream);

			// read the species length
			const unsigned char speciesLen = fgetc(mInStream);

			// read the genome assembly id length
			const unsigned char genomeAssemblyIDLen = fgetc(mInStream);

			// read the URI length
			const unsigned char uriLen = fgetc(mInStream);

			// read the number of bases
			fread((char*)&rs.NumBases, SIZEOF_INT, 1, mInStream);

			// read the concatenated begin coordinate
			fread((char*)&rs.Begin, SIZEOF_INT, 1, mInStream);

			// read the concatenated end coordinate
			fread((char*)&rs.End, SIZEOF_INT, 1, mInStream);

			// read the bases offset
			fread((char*)&rs.BasesOffset, SIZEOF_OFF_TYPE, 1, mInStream);

			// read the MD5 checksum
			rs.MD5.resize(32);
			char* pBuffer = (char*)rs.MD5.data();
			fread(pBuffer, 32, 1, mInStream);

			// read the reference name
			rs.Name.resize(nameLen);
			pBuffer = (char*)rs.Name.data();
			fread(pBuffer, nameLen, 1, mInStream);

			// read the species name
			if(speciesLen > 0) {
				rs.Species.resize(speciesLen);
				pBuffer = (char*)rs.Species.data();
				fread(pBuffer, speciesLen, 1, mInStream);
			}

			// read the genome assembly ID
			if(genomeAssemblyIDLen > 0) {
				rs.GenomeAssemblyID.resize(genomeAssemblyIDLen);
				pBuffer = (char*)rs.GenomeAssemblyID.data();
				fread(pBuffer, genomeAssemblyIDLen, 1, mInStream);
			}

			// read the URI
			if(uriLen > 0) {
				rs.URI.resize(uriLen);
				pBuffer = (char*)rs.URI.data();
				fread(pBuffer, uriLen, 1, mInStream);
			}

			//// DEBUG
			//cout << "# bases:                " << rs.NumBases << endl;
			//cout << "md5:                    " << rs.MD5 << endl;
			//cout << "name:                   " << rs.Name << endl;
			//cout << "species:                " << rs.Species << endl;
			//cout << "genome assembly ID:     " << rs.GenomeAssemblyID << endl;
			//cout << "URI:                    " << rs.URI << endl;

			// add the reference sequence to the index
			mIndex[rs.Name] = rs;
		}
	}
}
