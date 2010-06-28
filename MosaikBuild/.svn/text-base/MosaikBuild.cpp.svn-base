// ***************************************************************************
// BuildMain.cpp - imports the reads and reference sequences into MOSAIK.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "MosaikBuild.h"

// constructor
CMosaikBuild::CMosaikBuild(const MosaikReadFormat::ReadGroup& md)
: mReadGroup(md)
, mEnableColorspace(false)
, mHasBaseQualities(false)
, mHasReadNamePrefix(false)
, mHasReadLimit(false)
, mTrimReads(false)
, mTrimReadNames(false)
, mRemoveInstrumentInfo(false)
, mNumNBasesAllowed(NUM_N_BASES_ALLOWED)
, mNumLeadingNsTrimmed(0)
, mNumLaggingNsTrimmed(0)
, mNumMatesDeleted(0)
, mMinimumReadLength(MIN_READ_LENGTH)
, mBuffer(NULL)
, mBufferLen(256)
, mReadPrefixTrim(0)
, mReadSuffixTrim(0)
, mReadNamePrefixTrim(0)
, mReadNameSuffixTrim(0)
{
	// initialize the read and index buffer
	try {
		mBuffer = new unsigned char[mBufferLen];
	} catch(bad_alloc) {
		cout << "ERROR: Unable to allocate enough memory for the I/O buffer." << endl;
		exit(1);
	}
}

// destructor
CMosaikBuild::~CMosaikBuild(void) {
	mBufferLen = 0;
	if(mBuffer) {
		delete [] mBuffer;
		mBuffer = NULL;
	}
}

// activates the specified Illumina lanes
void CMosaikBuild::ActivateIlluminaLanes(const string& lanes) {

	// clear the allowed lanes
	for(unsigned char i = 0; i < 8; i++) mAllowedLanes[i] = false;

	// activate the necessary lanes
	unsigned int numLanes = lanes.size();
	const char* s = lanes.data();
	for(unsigned int i = 0; i < numLanes; i++) {
		unsigned char currentLane = s[i] - '1';
		mAllowedLanes[currentLane] = true;
	}
}

// creates a read group ID from the current time
void CMosaikBuild::CreateReadGroupID(string& readGroupID) {

	const char* symbols = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	char readGroup[256];
	char* pReadGroup = readGroup;
	const unsigned char outputBase = 36;

	// retrieve the current high-resolution time
	uint64_t hash = CTimeSupport::GetSystemTime();

	unsigned int count = 0;
	while((hash > 0) && (count++ < 256)) {
		*pReadGroup = symbols[hash % outputBase];
		pReadGroup++;
		hash = hash / outputBase;
	}

	*pReadGroup = 0;

	readGroupID.resize(count);
	memcpy((void*)readGroupID.data(), readGroup, count);

	// reverse the read group ID
	CSequenceUtilities::ReverseSequence((char*)readGroupID.data(), count);
}

// creates a MOSAIK reference archive
void CMosaikBuild::CreateReferenceArchive(const string& fastaFilename, const string& archiveFilename) {

	// initialize
	vector<ReferenceSequence> references;
	unsigned int concatenatedLength = 0;

	off_type concatenatedOffset   = 0;
	off_type concatenated2bOffset = 0;
	off_type indexOffset          = 0;
	off_type referenceBasesOffset = 0;
	off_type maskedRegionsOffset  = 0;

	CColorspaceUtilities csu;

	// load the reference sequences into the vector
	// TODO: fix this - we're using more memory than we should because of the transition from Mate to ReferenceSequence
	{
		unsigned int referenceBegin = 0;
		CMosaikString referenceName;
		Mosaik::Mate m;

		bool foundDuplicateReferenceName = false;
		set<string> referenceNames;
		set<string> duplicateReferenceNames;
		set<string>::const_iterator rnIter;
		string rnString;

		CConsole::Heading(); printf("- parsing reference sequences:\n"); CConsole::Reset();

		// initialize our counter
		unsigned int numRefSeqsParsed = 0;
		bool isRunning = true;
		CProgressCounter<unsigned int>::StartThread(&numRefSeqsParsed, &isRunning, "ref seqs");

		// open the FASTA file
		CFasta fasta;
		fasta.Open(fastaFilename);
		fasta.SetAssignedBaseQuality(mSettings.AssignedBaseQuality);

		// parse the FASTA file
		FastaTags ft;
		while(fasta.LoadNextMate(ft, m)) {

			// convert the reference sequence to colorspace if necessary
			if(mEnableColorspace) csu.ConvertReadBasespaceToPseudoColorspace(m.Bases);

			// sanity check: make sure we don't have the same reference name
			rnString = ft.Name.CData();
			rnIter = referenceNames.find(rnString);

			if(rnIter != referenceNames.end()) {
				duplicateReferenceNames.insert(rnString);
				foundDuplicateReferenceName = true;
			}

			referenceNames.insert(rnString);

			// assign the reference sequence data
			ReferenceSequence rs;
			rs.Name             = ft.Name.CData();
			rs.Bases            = m.Bases.CData();
			rs.NumBases         = rs.Bases.size();
			rs.Begin            = referenceBegin;
			rs.End              = referenceBegin + rs.NumBases - 1;

			// set the genome assembly ID
			if(ft.GenomeAssemblyID.Length() == 0) rs.GenomeAssemblyID = mSettings.GenomeAssemblyID;
			else rs.GenomeAssemblyID = ft.GenomeAssemblyID.CData();

			// set the species name
			if(ft.Species.Length() == 0) rs.Species = mSettings.Species;
			else rs.Species = ft.Species.CData();

			// set the URI
			if(ft.URI.Length() == 0) rs.URI = mSettings.UniformResourceIdentifier;
			else rs.URI = ft.URI.CData();

			concatenatedLength += rs.NumBases;
			referenceBegin = rs.End  + NUM_REFERENCE_DIVIDER_BASES + 1;

			if(!mEnableColorspace) CSequenceUtilities::UppercaseSequence(rs.Bases);
			references.push_back(rs);
			numRefSeqsParsed++;
		}

		// close the FASTA file
		fasta.Close();

		// stop the progress counter
		isRunning = false;
		CProgressCounter<unsigned int>::WaitThread();

		// stop processing the reference sequences if we found a duplicate reference sequence
		if(foundDuplicateReferenceName) {
			printf("\nERROR: Found duplicate reference sequence names:\n");
			for(rnIter = duplicateReferenceNames.begin(); rnIter != duplicateReferenceNames.end(); rnIter++) printf("- %s\n", rnIter->c_str());
			exit(1);
		}
	}

	// check if we have data
	if(references.empty()) {
		printf("ERROR: Unable to create the reference archive- could not find any reference sequences in the FASTA file.\n");
		exit(1);
	}

	// add the divider bases to the concatenated length
	if(references.size() > 1) concatenatedLength += (references.size() - 1) * NUM_REFERENCE_DIVIDER_BASES;

	// open the reference archive
	FILE* refStream = fopen(archiveFilename.c_str(), "wb");

	if(!refStream) {
		printf("ERROR: Unable to open the reference archive (%s) for writing.\n", archiveFilename.c_str());
		exit(1);
	}

	// ==============================
	// write to the reference archive
	// ==============================

	CFastLZIO fio;

	// initialize a reference spacer sequence (used between sequences in the concatenated sequence)
	char referenceDivider[NUM_REFERENCE_DIVIDER_BASES];
	uninitialized_fill(referenceDivider, referenceDivider + NUM_REFERENCE_DIVIDER_BASES, 'J');

	// ================
	// write the header
	// ================

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

	// write the signature
	fwrite("MSKRS\2", 6, 1, refStream);

	// write the status
	ReferenceSequenceStatus status = REF_UNKNOWN;
	if(mEnableColorspace) status |= REF_COLORSPACE;
	fputc(status, refStream);

	// write the archive date
	uint64_t currentTime = CTimeSupport::GetSystemTime();
	fwrite((char*)&currentTime, SIZEOF_UINT64, 1, refStream);

	// write the number of reference sequences
	const unsigned int numReferenceSequences = references.size();
	fwrite((char*)&numReferenceSequences, SIZEOF_INT, 1, refStream);

	// write the concatenated reference sequence length
	fwrite((char*)&concatenatedLength, SIZEOF_INT, 1, refStream);

	// write the concatenated reference offset [placeholder]
	fwrite((char*)&concatenatedOffset, SIZEOF_OFF_TYPE, 1, refStream);

	// write the concatenated 2-bit reference sequence length
	const unsigned int concatenated2bLength = (unsigned int)ceil(concatenatedLength / 4.0);
	fwrite((char*)&concatenated2bLength, SIZEOF_INT, 1, refStream);

	// write the concatenated 2-bit reference offset [placeholder]
	fwrite((char*)&concatenated2bOffset, SIZEOF_OFF_TYPE, 1, refStream);

	// write the index offset [placeholder]
	fwrite((char*)&indexOffset, SIZEOF_OFF_TYPE, 1, refStream);

	// write the reference bases offset [placeholder]
	fwrite((char*)&referenceBasesOffset, SIZEOF_OFF_TYPE, 1, refStream);

	// write the masked regions offset [placeholder]
	fwrite((char*)&maskedRegionsOffset, SIZEOF_OFF_TYPE, 1, refStream);

	// write the reserved field [placeholder]
	uint64_t reserved = 0;
	fwrite((char*)&reserved, SIZEOF_UINT64, 1, refStream);

	// ==================================
	// write the reference sequence bases
	// ==================================

	// get the current offset
	referenceBasesOffset = ftell64(refStream);

	CConsole::Heading(); printf("\n- writing reference sequences:\n"); CConsole::Reset();

	unsigned int currentRefSeq = 0;
	CProgressBar<unsigned int>::StartThread(&currentRefSeq, 0, numReferenceSequences, "ref seqs");

	// write the reference sequences
	vector<ReferenceSequence>::iterator rsIter;
	for(rsIter = references.begin(); rsIter != references.end(); rsIter++) {
		rsIter->BasesOffset = ftell64(refStream);
		fio.Write(rsIter->Bases.data(), rsIter->NumBases, refStream);

		currentRefSeq++;
	}

	CProgressBar<unsigned int>::WaitThread();

	// =========================
	// calculating MD5 checksums
	// =========================

	{
		CConsole::Heading(); printf("\n- calculating MD5 checksums:\n"); CConsole::Reset();

		unsigned char MD5[16];
		currentRefSeq = 0;
		CProgressBar<unsigned int>::StartThread(&currentRefSeq, 0, numReferenceSequences, "ref seqs");
		for(rsIter = references.begin(); rsIter != references.end(); rsIter++) {

			// generate the MD5 checksum
			memset(MD5, 0, 16);
			MD5_CTX context;
			MD5Init(&context); 
			MD5Update(&context, (unsigned char*)rsIter->Bases.data(), rsIter->NumBases);
			MD5Final(MD5, &context);

			// copy the checksum
			rsIter->MD5.resize(32);
			char* pMD5String = (char*)rsIter->MD5.data();
			for(unsigned int i = 0; i < 16; i++) {
				sprintf(pMD5String, "%02X", MD5[i]);
				pMD5String += 2;
			}

			currentRefSeq++;
		}

		CProgressBar<unsigned int>::WaitThread();
	}

	// ==================================
	// write the reference sequence index
	// ==================================

	// get the current offset
	indexOffset = ftell64(refStream);

	CConsole::Heading(); printf("\n- writing reference sequence index:\n"); CConsole::Reset();

	currentRefSeq = 0;
	CProgressBar<unsigned int>::StartThread(&currentRefSeq, 0, numReferenceSequences, "ref seqs");

	// write the index
	for(rsIter = references.begin(); rsIter != references.end(); rsIter++) {

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

		// enforce the maximum reference sequence name length
		unsigned int nameLen = (unsigned int)rsIter->Name.size();
		string referenceName = rsIter->Name;

		if(nameLen > 255) {
			nameLen = 255;
			referenceName = referenceName.substr(0, nameLen);
		}

		// enforce the maximum species name length
		unsigned int speciesLen = (unsigned int)rsIter->Species.size();
		string species = rsIter->Species;

		if(speciesLen > 255) {
			speciesLen = 255;
			species = species.substr(0, speciesLen);
		}

		// enforce the maximum genome assembly ID length
		unsigned int genomeAssemblyIDLen = (unsigned int)rsIter->GenomeAssemblyID.size();
		string genomeAssemblyID = rsIter->GenomeAssemblyID;

		if(genomeAssemblyIDLen > 255) {
			genomeAssemblyIDLen = 255;
			genomeAssemblyID = genomeAssemblyID.substr(0, genomeAssemblyIDLen);
		}

		// enforce the maximum URI length
		unsigned int uriLen = (unsigned int)rsIter->URI.size();
		string uri = rsIter->URI;

		if(uriLen > 255) {
			uriLen = 255;
			uri = uri.substr(0, uriLen);
		}

		// write the name length
		fputc((unsigned char)nameLen, refStream);

		// write the species length
		fputc((unsigned char)speciesLen, refStream);

		// write the genome assembly id length
		fputc((unsigned char)genomeAssemblyIDLen, refStream);

		// write the URI length
		fputc((unsigned char)uriLen, refStream);

		// write the number of bases
		fwrite((char*)&rsIter->NumBases, SIZEOF_INT, 1, refStream);

		// write the concatenated begin coordinate
		fwrite((char*)&rsIter->Begin, SIZEOF_INT, 1, refStream);

		// write the concatenated end coordinate
		fwrite((char*)&rsIter->End, SIZEOF_INT, 1, refStream);

		// write the bases offset
		fwrite((char*)&rsIter->BasesOffset, SIZEOF_OFF_TYPE, 1, refStream);

		// write the MD5 checksum
		fwrite(rsIter->MD5.data(), 32, 1, refStream);

		// write the reference name
		fwrite(referenceName.data(), nameLen, 1, refStream);

		// write the species name
		if(speciesLen > 0) fwrite(species.data(), speciesLen, 1, refStream);

		// write the genome assembly ID
		if(genomeAssemblyIDLen > 0) fwrite(genomeAssemblyID.data(), genomeAssemblyIDLen, 1, refStream);

		// write the genome assembly ID
		if(uriLen > 0) fwrite(uri.data(), uriLen, 1, refStream);

		currentRefSeq++;
	}

	CProgressBar<unsigned int>::WaitThread();

	// ========================================
	// creating concatenated reference sequence
	// ========================================

	char* concatenatedReference = NULL;

	try {
		concatenatedReference = new char[concatenatedLength + 1];
	} catch(bad_alloc) {
		printf("ERROR: Unable to allocate enough memory (%u bytes) to create the concatenated reference sequence.\n", concatenatedLength + 1);
		exit(1);
	}

	CConsole::Heading(); printf("\n- creating concatenated reference sequence:\n"); CConsole::Reset();

	currentRefSeq = 0;
	CProgressBar<unsigned int>::StartThread(&currentRefSeq, 0, numReferenceSequences, "ref seqs");

	// build the concatenated reference sequence
	char* pBuffer = concatenatedReference;
	for(rsIter = references.begin(); rsIter != references.end(); rsIter++) {

		// add the divider
		if(rsIter != references.begin()) {
			memcpy(pBuffer, referenceDivider, NUM_REFERENCE_DIVIDER_BASES);
			pBuffer += NUM_REFERENCE_DIVIDER_BASES;
		}

		// copy the reference sequence
		memcpy(pBuffer, rsIter->Bases.data(), rsIter->NumBases);
		pBuffer += rsIter->NumBases;

		currentRefSeq++;
	}

	CProgressBar<unsigned int>::WaitThread();

	// =======================================
	// writing concatenated reference sequence
	// =======================================

	// get the current offset
	concatenatedOffset = ftell64(refStream);

	printf("\n- writing concatenated reference sequence...        ");
	fflush(stdout);

	// write the concatenated reference sequence
	fio.Write(concatenatedReference, concatenatedLength, refStream);

	printf("finished.\n");

	// ==============================================
	// creating concatenated 2-bit reference sequence
	// ==============================================

	printf("- creating concatenated 2-bit reference sequence... ");
	fflush(stdout);

	char* concatenated2bReference = NULL;

	try {
		concatenated2bReference = new char[concatenated2bLength + 1];
	} catch(bad_alloc) {
		printf("ERROR: Unable to allocate enough memory (%u bytes) to create the concatenated 2-bit reference sequence.\n", concatenated2bLength + 1);
		exit(1);
	}

	// human genome 36.2
	//
	// A: 843953565
	// C: 584268578
	// G: 584621685
	// M: 1
	// N: 129484
	// R: 2
	// T: 845168978
	//                        A  B  C  D   E   F  G  H   I   J  K   L  M  N   O   P   Q  R  S  T   U  V  W   X  Y   Z
	char translation[26]  = { 0, 3, 1, 3, -1, -1, 2, 3, -1, -1, 3, -1, 0, 3, -1, -1, -1, 0, 2, 3, -1, 0, 3, -1, 3, -1 };

	unsigned int currentBase = 0;
	unsigned int offset = 0;
	char shift = 6;

	vector<MaskedPosition> maskedPositions;
	int maskIndex = -1;
	unsigned int lastMaskedBase = 0xfffffffe;

	while(currentBase < concatenatedLength) {

		// get the translated base
		char base = concatenatedReference[currentBase];
		char tBase = -1;
		if((base >= 'A') && (base <= 'Z')) tBase = translation[base - 'A'];

		// mask the base if it's not A, C, G, or T
		unsigned char twoBit  = 0;
		if((tBase >= 0) && (tBase <= 3)) {

			twoBit = tBase;

		} else {

			if(currentBase == (lastMaskedBase + 1)) {
				maskedPositions[maskIndex].End++;
			} else {
				MaskedPosition mp(currentBase);
				maskedPositions.push_back(mp);
				maskIndex++;
			}

			lastMaskedBase = currentBase;
		}

		// store the packed data
		concatenated2bReference[offset] |= twoBit  << shift;

		// update the bit shifts
		shift -= 2;
		if(shift < 0) {
			offset++;
			shift = 6;
		}

		currentBase++;
	}

	printf("finished.\n");

	// =============================================
	// writing concatenated 2-bit reference sequence
	// =============================================

	// get the current offset
	concatenated2bOffset = ftell64(refStream);

	printf("- writing concatenated 2-bit reference sequence...  ");
	fflush(stdout);

	// write the concatenated reference sequence
	fio.Write(concatenated2bReference, concatenated2bLength, refStream);

	printf("finished.\n");

	// ======================
	// writing masking vector
	// ======================

	// get the current offset
	maskedRegionsOffset = ftell64(refStream);

	// write the number of masked regions
	const unsigned int numMaskedRegions = maskedPositions.size();
	fwrite((char*)&numMaskedRegions, SIZEOF_INT, 1, refStream);

	if(numMaskedRegions > 0) {

		printf("- writing masking vector...                         ");
		fflush(stdout);

		const unsigned int numBytesWritten = numMaskedRegions * SIZEOF_INT * 2;

		char* maskedBuffer = NULL;
		try {
			maskedBuffer = new char[numBytesWritten];
		} catch(bad_alloc) {
			printf("ERROR: Unable to allocate memory for the masked region buffer.\n");
			exit(1);
		}

		pBuffer = maskedBuffer;
		vector<MaskedPosition>::const_iterator mpIter;
		for(mpIter = maskedPositions.begin(); mpIter != maskedPositions.end(); mpIter++) {
			memcpy(pBuffer, (char*)&mpIter->Begin, SIZEOF_INT);
			pBuffer += SIZEOF_INT;
			memcpy(pBuffer, (char*)&mpIter->End, SIZEOF_INT);
			pBuffer += SIZEOF_INT;
		}

		fio.Write(maskedBuffer, numBytesWritten, refStream);
		printf("finished.\n");

		// clean up
		delete [] maskedBuffer;
	}

	// clean up
	fio.Clear();
	delete [] concatenatedReference;
	delete [] concatenated2bReference;

	// =================
	// update the header
	// =================

	fseek64(refStream, 23, SEEK_SET);

	// write the concatenated reference offset
	fwrite((char*)&concatenatedOffset, SIZEOF_OFF_TYPE, 1, refStream);
	fseek64(refStream, SIZEOF_INT, SEEK_CUR);

	// write the concatenated 2-bit reference offset
	fwrite((char*)&concatenated2bOffset, SIZEOF_OFF_TYPE, 1, refStream);

	// write the index offset
	fwrite((char*)&indexOffset, SIZEOF_OFF_TYPE, 1, refStream);

	// write the reference bases offset
	fwrite((char*)&referenceBasesOffset, SIZEOF_OFF_TYPE, 1, refStream);

	// write the masked regions offset
	fwrite((char*)&maskedRegionsOffset, SIZEOF_OFF_TYPE, 1, refStream);

	// close the reference sequence archive
	fclose(refStream);
}

// Enables the processing of base qualities
void CMosaikBuild::EnableBaseQualities(const string& filename) {
	mHasBaseQualities                  = true;
	mSettings.BaseQualityFastaFilename = filename;
}

// Enables the processing of base qualities for the 2nd mate
void CMosaikBuild::EnableBaseQualities2(const string& filename) {
	mHasBaseQualities2                  = true;
	mSettings.BaseQualityFasta2Filename = filename;
}

// Enables trimming of bases and qualities
void CMosaikBuild::EnableBaseTrimming(const unsigned short prefixTrim, const unsigned short suffixTrim) {
	mTrimReads      = true;
	mReadPrefixTrim = prefixTrim;
	mReadSuffixTrim = suffixTrim;
}

// Enables SOLiD colorspace translation
void CMosaikBuild::EnableColorspace(void) {
	mEnableColorspace = true;
}

// Enables instrument info removal
void CMosaikBuild::EnableInstrumentInfoRemoval(void) {
	mRemoveInstrumentInfo = true;
}

// Enables trimming the first bases from the read name
void CMosaikBuild::EnableReadNameTrimming(const unsigned char prefixTrim, const unsigned char suffixTrim) {
	mTrimReadNames      = true;
	mReadNamePrefixTrim = prefixTrim;
	mReadNameSuffixTrim = suffixTrim;
}

// Enables the addition of a user specified read name prefix
void CMosaikBuild::EnableReadNamePrefix(const string& prefix) {
	mHasReadNamePrefix = true;
	mReadNamePrefix    = prefix.c_str();
}

// Enables a limit on the number of reads written to the read archive
void CMosaikBuild::EnableReadLimit(const uint64_t readLimit) {
	mHasReadLimit = true;
	mReadLimit    = readLimit;
}

// returns the colorspace name for the given read name
void CMosaikBuild::GetColorspaceName(const CMosaikString& readName, ColorspaceName& cn) {
	vector<string> columns;
	back_insert_iterator<vector<string> > backiter(columns);
	SplitString(backiter, "_", readName.CData());
	cn.first  = GetUnsignedShort((char*)columns[0].c_str());
	cn.second = GetUnsignedShort((char*)columns[1].c_str());
	cn.third  = GetUnsignedShort((char*)columns[2].c_str());
}

// Parses an Illumina Bustard directory
void CMosaikBuild::ParseBustard(const string& directory, const string& lanes, const string& outputFilename, const bool splitReads) {

	// activate our Illumina lanes
	ActivateIlluminaLanes(lanes);

	// create our Illumina base quality LUT
	char baseQualityLUT[256];
	int bq;
	for(int i = -128; i < 128; i++) {
		bq = (int)(-10.0 * log10(1.0 / (1.0 + pow(10.0, ((double)i / 10.0)))) + 0.5);
		if(bq < 1) bq = 1;
		baseQualityLUT[i + 128] = bq;
	}

	// find the sequence files
	vector<string> files, bustardFiles;
	CFileUtilities::SearchDirectory(files, directory.c_str());

	// Bustard files look like this: s_4_0200_seq.txt & s_4_0200_prb.txt
	for(unsigned int i = 0; i < (unsigned int)files.size(); i++) {

		// get the filepath
		string filepath = files[i];

		// extract the filename
		string filename = filepath;
		string::size_type slashPos = filepath.rfind(OS_DIRECTORY_SEPARATOR);
		if(slashPos != string::npos) filename = filepath.substr(slashPos + 1); 
		unsigned int filenameLen = filename.size();

		if(filenameLen != 16) continue;
		if(filename.substr(8) != "_seq.txt") continue;

		// extract the current lane
		unsigned char currentLane = filename[2] - '1';
		if(mAllowedLanes[currentLane]) bustardFiles.push_back(filepath.substr(0, filepath.size() - 8));
	}

	if(bustardFiles.empty()) {
		cout << "ERROR: could not find any Bustard files in the given directory." << endl;
		exit(1);
	}

	unsigned int numBustardFiles = (unsigned int)bustardFiles.size();
	cout << "- found " << numBustardFiles << " files" << endl;

	CConsole::Heading();
	SILENTMODE cout << endl << "- parsing Bustard file" << ((numBustardFiles > 1) ? "s" : "") << ":" << endl;
	CConsole::Reset();

	// initialize our writer
	MosaikReadFormat::CReadWriter writer;
	writer.Open(outputFilename, (splitReads ? RS_PAIRED_END_READ : RS_SINGLE_END_READ), mReadGroup);

	unsigned int fBufferSize = 4096;
	char* fBuffer = new char[fBufferSize];
	unsigned int lane, tile, xcoord, ycoord;
	int aBQ, cBQ, gBQ, tBQ, illuminaBQ;

	bool isRunning = true;
	unsigned int numReadsParsed = 0;
	CProgressCounter<unsigned int>::StartThread(&numReadsParsed, &isRunning, "reads");

	for(unsigned int currentFile = 0; currentFile < numBustardFiles; currentFile++) {

		// open our sequence file
		string seqFilename = bustardFiles[currentFile] + "_seq.txt";

		FILE* seq = NULL;
		fopen_s(&seq, seqFilename.c_str(), "rb");

		if(!seq) {
			cout << "ERROR: Could not open the Illumina Bustard sequence file (" << seqFilename << ") for reading." << endl;
			exit(1);
		}

		// open our probability file
		string prbFilename = bustardFiles[currentFile] + "_prb.txt";

		FILE* prb = NULL;
		fopen_s(&prb, prbFilename.c_str(), "rb");

		if(!prb) {
			cout << "ERROR: Could not open the Illumina Bustard probability file (" << prbFilename << ") for reading." << endl;
			exit(1);
		}

		// read the entire file
		while(true) {

			// get the next sequence line
			fscanf(seq, "%d\t%d\t%d\t%d\t%s", &lane, &tile, &xcoord, &ycoord, fBuffer);

			// stop if we found the EOF
			if(feof(seq)) break;

			// define a new read
			Mosaik::Read mr;
			mr.Mate1.Bases = fBuffer;
			unsigned int numBases = mr.Mate1.Bases.Length();
			mr.Mate1.Bases.Replace('.', 'N');
			mr.Mate1.Qualities.Reserve(numBases);

			// assign the read name
			sprintf_s(fBuffer, fBufferSize, "%d_%d_%d_%d", lane, tile, xcoord, ycoord);
			mr.Name = fBuffer;

			// retrieve the base qualities
			const char* bases = mr.Mate1.Bases.Data();
			char* qualities = (char*)mr.Mate1.Qualities.Data();
			for(unsigned int i = 0; i < numBases; i++) {

				fscanf(prb, "%d %d %d %d", &aBQ, &cBQ, &gBQ, &tBQ);

				switch(bases[i]) {
					case 'A':
						illuminaBQ = aBQ;
						break;
					case 'C':
						illuminaBQ = cBQ;
						break;
					case 'G':
						illuminaBQ = gBQ;
						break;
					case 'T':
						illuminaBQ = tBQ;
						break;
					default:
						illuminaBQ = aBQ;
						break;
				}

				// convert the Illumina base quality to the phred BQ
				qualities[i] = baseQualityLUT[illuminaBQ + 128];
			}

			if(splitReads) {

				unsigned char readLength = numBases / 2;

				// assign the second mate
				mr.Mate2.Bases     = mr.Mate1.Bases;
				mr.Mate2.Qualities = mr.Mate1.Bases;

				mr.Mate2.Bases.TrimBegin(readLength);
				mr.Mate2.Qualities.TrimBegin(readLength);

				// assign the first mate
				mr.Mate1.Bases.TrimEnd(readLength);
				mr.Mate1.Qualities.TrimEnd(readLength);

				// save both mates
				ProcessReadName(mr.Name);
				ProcessMate(mr.Mate1);
				ProcessMate(mr.Mate2);
				writer.SaveRead(mr);
				numReadsParsed++;

			} else {

				// save our read
				ProcessReadName(mr.Name);
				ProcessMate(mr.Mate1);
				writer.SaveRead(mr);
				numReadsParsed++;
			}

			// limit the number of reads being processed
			if(mHasReadLimit && (numReadsParsed >= mReadLimit)) break;
		}

		// close our files
		fclose(seq);
		fclose(prb);

		// limit the number of reads being processed
		if(mHasReadLimit && (numReadsParsed >= mReadLimit)) break;
	}

	// stop the progress counter
	isRunning = false;
	CProgressCounter<unsigned int>::WaitThread();

	// close files
	writer.Close();

	// cleanup
	delete [] fBuffer;

	// show some statistics
	Statistics s;
	s.IsPairedEnd     = splitReads;
	s.NumBasesWritten = writer.GetNumBases();
	s.NumReadsWritten = writer.GetNumReads();
	s.NumTotalMates   = (splitReads ? numReadsParsed * 2 : numReadsParsed);
}

// Parses the sequence and quality FASTA files while writing to our read archive
void CMosaikBuild::ParseFasta(const string& readFastaFilename, const string& outputFilename) {

	CConsole::Heading();
	SILENTMODE cout << "- parsing FASTA files:" << endl;
	CConsole::Reset();

	// initialize our writer
	MosaikReadFormat::CReadWriter writer;
	writer.Open(outputFilename, RS_SINGLE_END_READ, mReadGroup);

	// initialize our reader
	CFasta reader;
	reader.SetAssignedBaseQuality(mSettings.AssignedBaseQuality);
	if(mHasBaseQualities) reader.EnableBaseQualityFile(mSettings.BaseQualityFastaFilename);
	reader.Open(readFastaFilename);

	// initialize our counter
	unsigned int numReadsParsed = 0;
	bool isRunning = true;
	CProgressCounter<unsigned int>::StartThread(&numReadsParsed, &isRunning, "reads");

	CColorspaceUtilities csu;

	Mosaik::Read r;
	FastaTags ft;
	while(reader.LoadNextMate(ft, r.Mate1)) {

		if(mEnableColorspace) {
			const char* bases = r.Mate1.Bases.CData();
			memcpy((char*)&r.Mate1.SolidPrefixTransition, bases, 2);
			r.Mate1.Bases.TrimBegin(2);
			r.Mate1.Qualities.TrimBegin(1);
			csu.ConvertReadColorspaceToPseudoColorspace(r.Mate1.Bases);
		}

		r.Name = ft.Name;
		ProcessReadName(r.Name);
		ProcessMate(r.Mate1);
		writer.SaveRead(r);
		numReadsParsed++;
		if(mHasReadLimit && (numReadsParsed >= mReadLimit)) break;
	}

	// stop the progress counter
	isRunning = false;
	CProgressCounter<unsigned int>::WaitThread();

	// close the input and output files
	reader.Close();
	writer.Close();

	// show some statistics
	Statistics s;
	s.NumBasesWritten = writer.GetNumBases();
	s.NumReadsWritten = writer.GetNumReads();
	s.NumTotalMates   = numReadsParsed;
	ShowStatistics(s);
}

// Parses the sequence and quality paired-end FASTA files while writing to our read archive
void CMosaikBuild::ParsePEFasta(string& readFastaFilename, string& readFastaFilename2, const string& outputFilename) {

	CConsole::Heading();
	SILENTMODE cout << "- parsing paired-end/mate-pair FASTA files:" << endl;
	CConsole::Reset();

	// initialize our writer
	MosaikReadFormat::CReadWriter writer;
	writer.Open(outputFilename, RS_PAIRED_END_READ, mReadGroup);

	bool removedMateSuffix = false;
	unsigned short numSuffixCharactersRemoved = 0;

	// re-arrange the FASTA filenames if we are parsing a SOLiD csfasta file
	// the F3 read will always be mate1 and the R3 read will always be mate2
	if(mEnableColorspace) {
		ReorderSolidFastaFilenames(readFastaFilename, readFastaFilename2);
		if(mHasBaseQualities) ReorderSolidFastaFilenames(mSettings.BaseQualityFastaFilename, mSettings.BaseQualityFasta2Filename);
	}

	// initialize our mate 1 reader
	CFasta reader;
	reader.SetAssignedBaseQuality(mSettings.AssignedBaseQuality);
	if(mHasBaseQualities) reader.EnableBaseQualityFile(mSettings.BaseQualityFastaFilename);
	reader.Open(readFastaFilename);

	// initialize our mate 2 reader
	CFasta reader2;
	reader2.SetAssignedBaseQuality(mSettings.AssignedBaseQuality);
	if(mHasBaseQualities2) reader2.EnableBaseQualityFile(mSettings.BaseQualityFasta2Filename);
	reader2.Open(readFastaFilename2);

	// initialize our counter
	unsigned int numReadsParsed = 0;
	bool isRunning = true;
	CProgressCounter<unsigned int>::StartThread(&numReadsParsed, &isRunning, "reads");

	if(mEnableColorspace) {
		removedMateSuffix = true;
		numSuffixCharactersRemoved = 3;
	}

	uint64_t numMate1Orphaned = 0;
	uint64_t numMate2Orphaned = 0;
	uint64_t numMate1Written  = 0;
	uint64_t numMate2Written  = 0;
	unsigned int numMate1Read = 0;
	unsigned int numMate2Read = 0;

	CColorspaceUtilities csu;

	FastaTags ftags1, ftags2;
	ColorspaceName cn1, cn2;
	CMosaikString mate2Name;
	Mosaik::Read r;
	while(true) {

		bool ret1 = reader.LoadNextMate(ftags1, r.Mate1);
		bool ret2 = reader2.LoadNextMate(ftags2, r.Mate2);
		numMate1Read = 1;
		numMate2Read = 1;

		if(!ret1 || !ret2) break;

		// trim the read names
		if(!removedMateSuffix) {

			numSuffixCharactersRemoved = 0;
			while((ftags1.Name != ftags2.Name) && (ftags1.Name.Length() > 0)) {
				ftags1.Name.TrimEnd(1);
				ftags2.Name.TrimEnd(1);
				numSuffixCharactersRemoved++;
			}

			char lastChar = ftags1.Name[ftags1.Name.Length() - 1];
			if((lastChar == '/') || (lastChar == '_') || (lastChar == '.') || (lastChar == '|')) {
				ftags1.Name.TrimEnd(1);
				ftags2.Name.TrimEnd(1);
				numSuffixCharactersRemoved++;
			}

			removedMateSuffix = true;

		} else {

			ftags1.Name.TrimEnd(numSuffixCharactersRemoved);
			ftags2.Name.TrimEnd(numSuffixCharactersRemoved);
		}

		if(ftags1.Name != ftags2.Name) {

			if(mEnableColorspace) {
				GetColorspaceName(ftags1.Name, cn1);
				GetColorspaceName(ftags2.Name, cn2);

				while(ftags1.Name != ftags2.Name) {
					if(cn1 < cn2) {
						if(!reader.LoadNextMate(ftags1, r.Mate1)) break;
						numMate1Read++;
						ftags1.Name.TrimEnd(numSuffixCharactersRemoved);
						GetColorspaceName(ftags1.Name, cn1);
					} else if(cn2 < cn1) {
						if(!reader2.LoadNextMate(ftags2, r.Mate2)) break;
						numMate2Read++;
						ftags2.Name.TrimEnd(numSuffixCharactersRemoved);
						GetColorspaceName(ftags2.Name, cn2);
					}
				}

				if(ftags1.Name != ftags2.Name) {
					cout << "ERROR: Resynchronization colorspace support failed." << endl;
					exit(1);
				}

			} else {
				cout << "ERROR: The mate1 read name did not match the mate2 read name. Resynchronization support needs to be implemented." << endl;
				cout << "mate 1 name: " << ftags1.Name << endl;
				cout << "mate 2 name: " << ftags2.Name << endl;
				exit(1);
			}
		}

		if(mEnableColorspace) {
			const char* bases = r.Mate1.Bases.CData();
			memcpy((char*)&r.Mate1.SolidPrefixTransition, bases, 2);
			r.Mate1.Bases.TrimBegin(2);
			r.Mate1.Qualities.TrimBegin(1);

			const char* bases2 = r.Mate2.Bases.CData();
			memcpy((char*)&r.Mate2.SolidPrefixTransition, bases2, 2);
			r.Mate2.Bases.TrimBegin(2);
			r.Mate2.Qualities.TrimBegin(1);

			csu.ConvertReadColorspaceToPseudoColorspace(r.Mate1.Bases);
			csu.ConvertReadColorspaceToPseudoColorspace(r.Mate2.Bases);
		}

		numMate1Orphaned += numMate1Read - 1;
		numMate2Orphaned += numMate2Read - 1;
		numMate1Written++;
		numMate2Written++;

		r.Name = ftags1.Name;
		ProcessReadName(r.Name);
		ProcessMate(r.Mate1);
		ProcessMate(r.Mate2);
		writer.SaveRead(r);
		numReadsParsed++;
		if(mHasReadLimit && (numReadsParsed >= mReadLimit)) break;
	}

	// stop the progress counter
	isRunning = false;
	CProgressCounter<unsigned int>::WaitThread();

	// close the input and output files
	reader.Close();
	reader2.Close();
	writer.Close();

	// show some statistics
	Statistics s;
	s.IsPairedEnd      = true;
	s.NumBasesWritten  = writer.GetNumBases();
	s.NumReadsWritten  = writer.GetNumReads();
	s.NumTotalMates    = numReadsParsed;
	s.NumMate1Orphaned = numMate1Orphaned;
	s.NumMate2Orphaned = numMate2Orphaned;
	ShowStatistics(s);
}

// Parses the reads and base qualities from a FASTQ file
void CMosaikBuild::ParseFastq(vector<string>& fastqFiles, const string& outputFilename) {

	CConsole::Heading();
	SILENTMODE cout << "- parsing FASTQ file" << ((fastqFiles.size() > 1) ? "s" : "") << ":" << endl;
	CConsole::Reset();

	// initialize our writer
	MosaikReadFormat::CReadWriter writer;
	writer.Open(outputFilename, RS_SINGLE_END_READ, mReadGroup);

	CColorspaceUtilities csu;

	bool isRunning = true;
	unsigned int numReadsParsed = 0;
	CProgressCounter<unsigned int>::StartThread(&numReadsParsed, &isRunning, "reads");

	// process each FASTQ file
	const unsigned int numFastqFiles = (unsigned int)fastqFiles.size();
	for(unsigned int currentFile = 0; currentFile < numFastqFiles; currentFile++) {

		// initialize our reader
		CFastq reader;
		reader.Open(fastqFiles[currentFile]);

		Mosaik::Read r;
		while(reader.LoadNextMate(r.Name, r.Mate1)) {

			if(mEnableColorspace) {
				const char* bases = r.Mate1.Bases.CData();
				memcpy((char*)&r.Mate1.SolidPrefixTransition, bases, 2);
				r.Mate1.Bases.TrimBegin(2);
				r.Mate1.Qualities.TrimBegin(2);
				csu.ConvertReadColorspaceToPseudoColorspace(r.Mate1.Bases);
			}

			ProcessReadName(r.Name);
			ProcessMate(r.Mate1);
			writer.SaveRead(r);
			numReadsParsed++;

			if(mHasReadLimit && (numReadsParsed >= mReadLimit)) break;
		}

		// close our reader
		reader.Close();

		// limit the number of reads being processed
		if(mHasReadLimit && (numReadsParsed >= mReadLimit)) break;
	}

	// stop the progress counter
	isRunning = false;
	CProgressCounter<unsigned int>::WaitThread();

	// close files
	writer.Close();

	// show some statistics
	Statistics s;
	s.NumBasesWritten = writer.GetNumBases();
	s.NumReadsWritten = writer.GetNumReads();
	s.NumTotalMates   = numReadsParsed;
	ShowStatistics(s);
}

// Parses the reads and base qualities from a paired-end FASTQ file
void CMosaikBuild::ParsePEFastq(vector<string>& mate1Files, vector<string>& mate2Files, const string& outputFilename) {

	// check that we have the same number of files for each mate
	const unsigned int numMate1Files = mate1Files.size();
	const unsigned int numMate2Files = mate2Files.size();

	if(numMate1Files != numMate2Files) {
		cout << "ERROR: Found a different number of files for mate1 (" << numMate1Files << ") and mate2 (" << numMate2Files << ")." << endl;
		exit(1);
	}

	// sort the filenames
	sort(mate1Files.begin(), mate1Files.end());
	sort(mate2Files.begin(), mate2Files.end());

	CConsole::Heading();
	SILENTMODE cout << "- parsing paired-end/mate-pair FASTQ files:" << endl;
	CConsole::Reset();

	// initialize our writer
	MosaikReadFormat::CReadWriter writer;
	writer.Open(outputFilename, RS_PAIRED_END_READ, mReadGroup);

	CColorspaceUtilities csu;

	bool removedMateSuffix = false;
	unsigned short numSuffixCharactersRemoved = 0;

	bool isRunning = true;
	unsigned int numReadsParsed = 0;
	CProgressCounter<unsigned int>::StartThread(&numReadsParsed, &isRunning, "reads");

	// process each FASTQ file
	for(unsigned int currentFile = 0; currentFile < numMate1Files; currentFile++) {

		// initialize our reader
		CFastq reader1, reader2;
		reader1.Open(mate1Files[currentFile]);
		reader2.Open(mate2Files[currentFile]);

		CMosaikString mate2Name;
		Mosaik::Read r;
		while(true) {
			bool ret1 = reader1.LoadNextMate(r.Name, r.Mate1);
			bool ret2 = reader2.LoadNextMate(mate2Name, r.Mate2);

			if(!ret1 && !ret2) break;
			if(!ret1 || !ret2) {
				cout << "ERROR: One of the FASTQ parsers reached the end of the file before the other." << endl;
				exit(1);
			}

			// trim the read names
			if(!removedMateSuffix) {
				numSuffixCharactersRemoved = 0;
				while((r.Name != mate2Name) && (r.Name.Length() > 0)) {
					r.Name.TrimEnd(1);
					mate2Name.TrimEnd(1);
					numSuffixCharactersRemoved++;
				}

				char lastChar = r.Name[r.Name.Length()-1];
				if((lastChar == '/') || (lastChar == '_') || (lastChar == '.') || (lastChar == '|')) {
					r.Name.TrimEnd(1);
					mate2Name.TrimEnd(1);
					numSuffixCharactersRemoved++;
				}

				removedMateSuffix = true;

			} else {

				r.Name.TrimEnd(numSuffixCharactersRemoved);
				mate2Name.TrimEnd(numSuffixCharactersRemoved);
			}

			if(r.Name != mate2Name) {
				cout << "ERROR: The mate1 read name did not match the mate2 read name. Resynchronization support needs to be implemented." << endl;
				exit(1);
			}

			if(mEnableColorspace) {
				const char* bases = r.Mate1.Bases.CData();
				memcpy((char*)&r.Mate1.SolidPrefixTransition, bases, 2);
				r.Mate1.Bases.TrimBegin(2);
				r.Mate1.Qualities.TrimBegin(2);
				csu.ConvertReadColorspaceToPseudoColorspace(r.Mate1.Bases);

				const char* bases2 = r.Mate2.Bases.CData();
				memcpy((char*)&r.Mate2.SolidPrefixTransition, bases2, 2);
				r.Mate2.Bases.TrimBegin(2);
				r.Mate2.Qualities.TrimBegin(2);
				csu.ConvertReadColorspaceToPseudoColorspace(r.Mate2.Bases);
			}

			ProcessReadName(r.Name);
			ProcessMate(r.Mate1);
			ProcessMate(r.Mate2);
			writer.SaveRead(r);
			numReadsParsed++;
			if(mHasReadLimit && (numReadsParsed >= mReadLimit)) break;
		}

		// close our reader
		reader1.Close();
		reader2.Close();

		// limit the number of reads being processed
		if(mHasReadLimit && (numReadsParsed >= mReadLimit)) break;
	}

	// stop the progress counter
	isRunning = false;
	CProgressCounter<unsigned int>::WaitThread();

	// close files
	writer.Close();

	// show some statistics
	Statistics s;
	s.IsPairedEnd     = true;
	s.NumBasesWritten = writer.GetNumBases();
	s.NumReadsWritten = writer.GetNumReads();
	s.NumTotalMates   = numReadsParsed * 2;
	ShowStatistics(s);
}

// Parses an Illumina Gerald directory
void CMosaikBuild::ParseGerald(const string& directory, const string& lanes, const string& outputFilename) {

	// activate our Illumina lanes
	ActivateIlluminaLanes(lanes);

	// create our Illumina base quality LUT
	char baseQualityLUT[256];
	int bq;
	for(int i = -128; i < 128; i++) {
		bq = (int)(-10.0 * log10(1.0 / (1.0 + pow(10.0, ((double)i / 10.0)))) + 0.5);
		if(bq < 1) bq = 1;
		baseQualityLUT[i + 128] = bq;
	}

	// find the sequence files
	vector<string> files, geraldFiles;
	CFileUtilities::SearchDirectory(files, directory.c_str());

	// Gerald files look like this: s_1_sequence.txt
	for(unsigned int i = 0; i < (unsigned int)files.size(); i++) {

		// get the filepath
		string filepath = files[i];

		// extract the filename
		string filename = filepath;
		string::size_type slashPos = filepath.rfind(OS_DIRECTORY_SEPARATOR);
		if(slashPos != string::npos) filename = filepath.substr(slashPos + 1); 
		unsigned int filenameLen = filename.size();

		if(filenameLen != 16) continue;
		if(filename.substr(3) != "_sequence.txt") continue;

		// extract the current lane
		unsigned char currentLane = filename[2] - '1';
		if(mAllowedLanes[currentLane]) geraldFiles.push_back(filepath);
	}

	if(geraldFiles.empty()) {
		cout << "ERROR: could not find any Gerald files in the given directory." << endl;
		exit(1);
	}

	const unsigned int numGeraldFiles = (unsigned int)geraldFiles.size();
	cout << "- found " << numGeraldFiles << " files" << endl;

	CConsole::Heading();
	SILENTMODE cout << "- parsing Gerald file" << ((numGeraldFiles > 1) ? "s" : "") << ":" << endl;
	CConsole::Reset();

	// initialize our writer
	MosaikReadFormat::CReadWriter writer;
	writer.Open(outputFilename, RS_SINGLE_END_READ, mReadGroup);

	bool isRunning = true;
	unsigned int numReadsParsed = 0;
	CProgressCounter<unsigned int>::StartThread(&numReadsParsed, &isRunning, "reads");

	// process each FASTQ file
	for(unsigned int currentFile = 0; currentFile < numGeraldFiles; currentFile++) {

		// initialize our reader
		CFastq reader;
		reader.Open(geraldFiles[currentFile]);
		reader.SetOffset(ILLUMINA_FASTQ_OFFSET);

		Mosaik::Read r;
		while(reader.LoadNextMate(r.Name, r.Mate1)) {
			ProcessReadName(r.Name);
			ProcessMate(r.Mate1);
			writer.SaveRead(r);
			numReadsParsed++;
			if(mHasReadLimit && (numReadsParsed >= mReadLimit)) break;
		}

		// close our reader
		reader.Close();

		// limit the number of reads being processed
		if(mHasReadLimit && (numReadsParsed >= mReadLimit)) break;
	}

	// stop the progress counter
	isRunning = false;
	CProgressCounter<unsigned int>::WaitThread();

	// close files
	writer.Close();

	// show some statistics
	Statistics s;
	s.NumBasesWritten = writer.GetNumBases();
	s.NumReadsWritten = writer.GetNumReads();
	s.NumTotalMates   = numReadsParsed;
	ShowStatistics(s);
}

// Parses the SRF archive
void CMosaikBuild::ParseSRF(vector<string>& srfFiles, const string& outputFilename) {

	CConsole::Heading();
	SILENTMODE cout << "- parsing SRF file" << ((srfFiles.size() > 1) ? "s" : "") << ":" << endl;
	CConsole::Reset();

	// initialize our writer
	MosaikReadFormat::CReadWriter writer;
	writer.Open(outputFilename, RS_SINGLE_END_READ, mReadGroup);

	bool isRunning = true;
	unsigned int numReadsParsed = 0;
	CProgressCounter<unsigned int>::StartThread(&numReadsParsed, &isRunning, "reads");

	// process each SRF file
	unsigned int numSrfFiles = (unsigned int)srfFiles.size();
	for(unsigned int currentFile = 0; currentFile < numSrfFiles; currentFile++) {

		// initialize our reader
		CSRF reader;
		reader.Open(srfFiles[currentFile]);

		Mosaik::Read mr;
		while(reader.GetRead(mr)) {
			mr.Mate1.Bases.Uppercase();
			ProcessReadName(mr.Name);
			ProcessMate(mr.Mate1);
			writer.SaveRead(mr);
			numReadsParsed++;

			// limit the number of reads being processed
			if(mHasReadLimit && (numReadsParsed >= mReadLimit)) break;
		}

		// close our reader
		reader.Close();

		// limit the number of reads being processed
		if(mHasReadLimit && (numReadsParsed >= mReadLimit)) break;
	}

	// stop the progress counter
	isRunning = false;
	CProgressCounter<unsigned int>::WaitThread();

	// close files
	writer.Close();

	// show some statistics
	Statistics s;
	s.NumBasesWritten = writer.GetNumBases();
	s.NumReadsWritten = writer.GetNumReads();
	s.NumTotalMates   = numReadsParsed;
	ShowStatistics(s);
}

// trims the mate
void CMosaikBuild::ProcessMate(Mosaik::Mate& mate) {

	if(mTrimReads) {

		if(mReadSuffixTrim > 0) {
			mate.Bases.TrimEnd(mReadSuffixTrim);
			mate.Qualities.TrimEnd(mReadSuffixTrim);
		}

		if(mReadPrefixTrim > 0) {
			mate.Bases.TrimBegin(mReadPrefixTrim);
			mate.Qualities.TrimBegin(mReadPrefixTrim);
		}
	}

	// trim the leading and lagging N's
	const char* pBases = mate.Bases.CData();
	int pos = 0;
	unsigned short numLeadingNs = 0;
	while((pBases[pos] == 'N') && (pos < (int)mate.Bases.Length())) {
		numLeadingNs++;
		pos++;
	}

	mate.Bases.TrimBegin(numLeadingNs);
	mate.Qualities.TrimBegin(numLeadingNs);

	pBases = mate.Bases.CData();
	pos = mate.Bases.Length() - 1;
	unsigned short numLaggingNs = 0;
	while((pBases[pos] == 'N') && (pos >= 0)) {
		numLaggingNs++;
		pos--;
	}

	mate.Bases.TrimEnd(numLaggingNs);
	mate.Qualities.TrimEnd(numLaggingNs);

	// determine if the read has too many remaining N's
	unsigned short numRemainingNs = 0;
	const unsigned int mateLength = mate.Bases.Length();
	pBases = mate.Bases.CData();
	for(unsigned short i = 0; i < mateLength; i++) if(pBases[i] == 'N') numRemainingNs++;

	bool deleteMate = false;
	if((mNumNBasesAllowed != 0) && (numRemainingNs > mNumNBasesAllowed)) deleteMate = true;
	if(mateLength < mMinimumReadLength)                                  deleteMate = true;

	// delete the mate
	if(deleteMate) {
		mate.Bases.SetLength(0);
		mate.Qualities.SetLength(0);
		mNumMatesDeleted++;
	} else {
		mNumLeadingNsTrimmed += numLeadingNs;
		mNumLaggingNsTrimmed += numLaggingNs;
	}
}

// trims the read name and adds a read name prefix
void CMosaikBuild::ProcessReadName(CMosaikString& readName) {
	if(mTrimReadNames && (mReadNamePrefixTrim > 0)) readName.TrimBegin(mReadNamePrefixTrim);
	if(mHasReadNamePrefix)                          readName.Prepend(mReadNamePrefix);
	if(mTrimReadNames && (mReadNameSuffixTrim > 0)) readName.TrimEnd(mReadNameSuffixTrim);
}

// returns true if a swap had to performed to guarantee that filename1 is the F3 read (filename2 = R3)
void CMosaikBuild::ReorderSolidFastaFilenames(string& filename1, string& filename2) {

	// initialize
	char buffer[2048];
#ifdef WIN32
	regex tagRegex("\\s--tag=(\\S+)\\s"); // " --tag=F3 "
	cmatch results;
#endif

	// ======================================
	// retrieve the first line from filename1
	// ======================================

	gzFile f1Stream = gzopen(filename1.c_str(), "rb");

	if(f1Stream == NULL) {
		printf("ERROR: Could not open FASTA file (%s) when reordering filenames.\n", filename1.c_str());
		exit(1);
	}

	gzgets(f1Stream, buffer, 2048);
	gzclose(f1Stream);

	string f1Header = buffer;

	// ======================================
	// retrieve the first line from filename2
	// ======================================

	gzFile f2Stream = gzopen(filename2.c_str(), "rb");

	if(f2Stream == NULL) {
		printf("ERROR: Could not open FASTA file (%s) when reordering filenames.\n", filename2.c_str());
		exit(1);
	}

	gzgets(f2Stream, buffer, 2048);
	gzclose(f2Stream);

	string f2Header = buffer;

	// sanity checks
	if(f1Header[0] != '#') {
		printf("ERROR: filename1 (%s) didn't contain a valid csfasta header.\n", filename1.c_str());
		exit(1);
	}

	if(f2Header[0] != '#') {
		printf("ERROR: filename2 (%s) didn't contain a valid csfasta header.\n", filename2.c_str());
		exit(1);
	}

	// ===========================
	// find the tag from filename1
	// ===========================

#ifdef WIN32

	if(!regex_search(f1Header.c_str(), results, tagRegex)) {
		printf("ERROR: Unable to find the tag keyword in the SOLiD FASTA header.\n");
		exit(1);
	}

	string tag1 = results[1].str().c_str();

#else

	// TODO: replace this with the TR1 regex above when it finally works in gcc. It doesn't work in gcc 4.3.3
	// find the tag element
	const unsigned int f1HeaderLen = f1Header.size();
	const char* pHeader = f1Header.c_str();

	const string TAG = "--tag=";
	string::size_type tagPos = f1Header.find(TAG.c_str());

	if(tagPos == string::npos) {
		printf("ERROR: Unable to find the tag keyword in the SOLiD FASTA header.\n");
		exit(1);
	}

	unsigned int start = tagPos + TAG.size();

	unsigned int stop = start;
	if(stop < f1HeaderLen) {
		while((pHeader[stop] != 32) && (pHeader[stop] != 9) && (pHeader[stop] != 10) && (pHeader[stop] != 13)) {
			stop++;
			if(stop == f1HeaderLen) break;
		}
	}

	if(start == stop) {
		printf("ERROR: Unable to parse tag from the SOLiD FASTA header.\n");
		exit(1);
	}

	string tag1 = f1Header.substr(start, stop - start);

#endif

	// sanity check
	if((tag1 != "F3") && (tag1 != "R3")) {
		printf("ERROR: The tag in filename1 (%s) was neither 'F3' nor 'R3'.\n", filename1.c_str());
		exit(1);
	}

	// ===========================
	// find the tag from filename2
	// ===========================

#ifdef WIN32

	if(!regex_search(f2Header.c_str(), results, tagRegex)) {
		printf("ERROR: Unable to find the tag keyword in the SOLiD FASTA header.\n");
		exit(1);
	}

	string tag2 = results[1].str().c_str();

#else

	// TODO: replace this with the TR1 regex above when it finally works in gcc. It doesn't work in gcc 4.3.3
	// find the tag element
	const unsigned int f2HeaderLen = f2Header.size();
	pHeader = f2Header.c_str();

	tagPos = f2Header.find(TAG.c_str());

	if(tagPos == string::npos) {
		printf("ERROR: Unable to find the tag keyword in the SOLiD FASTA header.\n");
		exit(1);
	}

	start = tagPos + TAG.size();

	stop = start;
	if(stop < f2HeaderLen) {
		while((pHeader[stop] != 32) && (pHeader[stop] != 9) && (pHeader[stop] != 10) && (pHeader[stop] != 13)) {
			stop++;
			if(stop == f2HeaderLen) break;
		}
	}

	if(start == stop) {
		printf("ERROR: Unable to parse tag from the SOLiD FASTA header.\n");
		exit(1);
	}

	string tag2 = f2Header.substr(start, stop - start);

#endif

	// sanity checks
	if((tag2 != "F3") && (tag2 != "R3")) {
		printf("ERROR: The tag in filename2 (%s) was neither 'F3' nor 'R3'.\n", filename2.c_str());
		exit(1);
	}

	if(tag1 == tag2) {
		printf("ERROR: The tags in filename1 (%s) and filename2 (%s) are the same: '%s'.\n", filename1.c_str(), filename2.c_str(), tag1.c_str());
		exit(1);
	}

	// ==========================
	// reorder the tags if needed
	// ==========================

	if(tag1 == "R3") {
		string t = filename1;
		filename1 = filename2;
		filename2 = t;
	}
}

// Sets the default base quality when a data set lacks BQ data
void CMosaikBuild::SetAssignedBaseQuality(unsigned char baseQuality) {
	mSettings.AssignedBaseQuality = baseQuality;
}

// Sets the Genome Assembly ID [used when creating reference archives]
void CMosaikBuild::SetGenomeAssemblyID(const string& id) {
	mSettings.GenomeAssemblyID = id;
}

// Sets the maximum number of N's allowed
void CMosaikBuild::SetNumNBasesAllowed(const unsigned char numNBasesAllowed) {
	mNumNBasesAllowed = numNBasesAllowed;
}

// Sets the species name [used when creating reference archives]
void CMosaikBuild::SetSpecies(const string& name) {
	mSettings.Species = name;
}

// Sets the uniform resource identifier
void CMosaikBuild::SetURI(const string& uri) {
	mSettings.UniformResourceIdentifier = uri;
}

// shows the conversion statistics
void CMosaikBuild::ShowStatistics(const Statistics& s) {

	cout << endl;
	CConsole::Heading(); cout << "Filtering statistics:" << endl; CConsole::Reset();
	cout << "============================================" << endl;

	// show trimming statistics
	if(mNumMatesDeleted > 0)
		cout << "# " << (s.IsPairedEnd ? "mates" : "reads") << " deleted:       " << setw(11) << mNumMatesDeleted << " (" << setw(5) << fixed << setprecision(1) << (mNumMatesDeleted / (double)s.NumTotalMates) * 100.0 << " %)" << endl;

	if(mNumLeadingNsTrimmed > 0)
		cout << "# leading N's trimmed: " << setw(11) << mNumLeadingNsTrimmed << endl;

	if(mNumLaggingNsTrimmed > 0)
		cout << "# lagging N's trimmed: " << setw(11) << mNumLaggingNsTrimmed << endl;

	if((mNumMatesDeleted > 0) || (mNumLeadingNsTrimmed > 0) || (mNumLaggingNsTrimmed > 0)) cout << "--------------------------------------------" << endl;

	// show orphan statistics
	if(s.NumMate1Orphaned > 0)
		cout << "# orphaned mate1:      " << setw(11) << s.NumMate1Orphaned << endl;

	if(s.NumMate2Orphaned > 0)
		cout << "# orphaned mate2:      " << setw(11) << s.NumMate2Orphaned << endl;

	if((s.NumMate1Orphaned > 0) || (s.NumMate2Orphaned > 0)) cout << "--------------------------------------------" << endl;

	// show general statistics
	if(s.NumReadsWritten > 0)
		cout << "# reads written:       " << setw(11) << s.NumReadsWritten << endl;

	if(s.NumBasesWritten > 0)
		cout << "# bases written:       " << setw(11) << s.NumBasesWritten << endl;
}
