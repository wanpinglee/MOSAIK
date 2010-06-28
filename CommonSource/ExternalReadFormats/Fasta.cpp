// ***************************************************************************
// CFasta - imports reads from the FASTA file format.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "Fasta.h"

// define our uncompressed/compressed macros
#define READ_FILE_CLOSE          (mAreBasesCompressed ? gzclose(mInZStream)                               : fclose(mInStream))
#define READ_FILE_EOF            (mAreBasesCompressed ? gzeof(mInZStream)                                 : feof(mInStream))
#define READ_FILE_GETC           (mAreBasesCompressed ? gzgetc(mInZStream)                                : fgetc(mInStream))
#define READ_FILE_GETS           (mAreBasesCompressed ? gzgets(mInZStream, mBaseBuffer, mBaseBufferLen)   : fgets(mBaseBuffer, mBaseBufferLen, mInStream))
#define READ_FILE_OFFSET         (mAreBasesCompressed ? gztell(mInZStream)                                : ftell64(mInStream))
#define READ_FILE_OPEN(filename) (mAreBasesCompressed ? mInZStream = gzopen(filename, "rb")               : mInStream = fopen(filename, "rb"))
#define READ_FILE_REWIND         (mAreBasesCompressed ? gzseek(mInZStream, mReadDataBaseOffset, SEEK_SET) : fseek64(mInStream, mReadDataBaseOffset, SEEK_SET))
#define READ_FILE_UNGETC(ch)     (mAreBasesCompressed ? gzungetc(ch, mInZStream)                          : ungetc(ch, mInStream))
#define READ_STREAM              (mAreBasesCompressed ? mInZStream                                        : mInStream)

#define QUALITY_FILE_CLOSE          (mAreBasesCompressed         ? gzclose(mInQualityZStream)                                   : fclose(mInQualityStream))
#define QUALITY_FILE_EOF            (mAreBaseQualitiesCompressed ? gzeof(mInQualityZStream)                                     : feof(mInQualityStream))
#define QUALITY_FILE_GETC           (mAreBaseQualitiesCompressed ? gzgetc(mInQualityZStream)                                    : fgetc(mInQualityStream))
#define QUALITY_FILE_GETS           (mAreBaseQualitiesCompressed ? gzgets(mInQualityZStream, mQualityBuffer, mQualityBufferLen) : fgets(mQualityBuffer, mQualityBufferLen, mInQualityStream))
#define QUALITY_FILE_OFFSET         (mAreBaseQualitiesCompressed ? gztell(mInQualityZStream)                                    : ftell64(mInQualityStream))
#define QUALITY_FILE_OPEN(filename) (mAreBaseQualitiesCompressed ? mInQualityZStream = gzopen(filename, "rb")                   : mInQualityStream = fopen(filename, "rb"))
#define QUALITY_FILE_REWIND         (mAreBaseQualitiesCompressed ? gzseek(mInQualityZStream, mReadDataQualityOffset, SEEK_SET)  : fseek64(mInQualityStream, mReadDataQualityOffset, SEEK_SET))
#define QUALITY_FILE_UNGETC(ch)     (mAreBaseQualitiesCompressed ? gzungetc(ch, mInQualityZStream)                              : ungetc(ch, mInQualityStream))
#define QUALITY_STREAM              (mAreBaseQualitiesCompressed ? mInQualityZStream                                            : mInQualityStream)

// constructor
CFasta::CFasta(void)
: mIsOpen(false)
, mAreBasesCompressed(false)
, mAreBaseQualitiesCompressed(false)
, mHasBaseQualityFile(false)
, mInStream(NULL)
, mInQualityStream(NULL)
, mInZStream(NULL)
, mInQualityZStream(NULL)
, mBaseBuffer(NULL)
, mBaseBufferLen(0)
, mQualityBuffer(NULL)
, mQualityBufferLen(0)
, mReadDataBaseOffset(0)
, mReadDataQualityOffset(0)
, mAssignedBaseQuality(30)
{}

// destructor
CFasta::~CFasta(void) {
	if(mIsOpen)        Close();
	if(mBaseBuffer)    delete [] mBaseBuffer;
	if(mQualityBuffer) delete [] mQualityBuffer;
}

// Checks if a file is truly a FASTA file
bool CFasta::CheckFile(const string& filename, const bool showError) {

	// open the FASTA file
	bool foundError = false;
	gzFile checkStream = gzopen(filename.c_str(), "rb");
	if(checkStream == NULL) {
		if(showError) {
			cout << "ERROR: Could not open FASTA file (" << filename << ") when performing integrity check." << endl;
			exit(1);
		}

		foundError = true;
	}

	// retrieve the first character
	const unsigned int BUFFERLEN = 1024;
	char* buffer = new char[BUFFERLEN];
	buffer[0] = '#';
	if(!foundError) {

		while(buffer[0] == '#' && !gzeof(checkStream)) gzgets(checkStream, buffer, BUFFERLEN);

		// check if the FASTA file starts with the header
		if(buffer[0] != '>') {
			if(showError) {
				cout << "ERROR: It seems that the input file (" << filename << ") is not in FASTA format." << endl;
				cout << "       buffer: " << buffer << endl;
				exit(1);
			}

			foundError = true;
		}
	}
	delete [] buffer;

	// close the file
	gzclose(checkStream);

	// return the appropriate values
	if(foundError) return false;
	return true;
}

// closes the FASTA file
void CFasta::Close(void) {
	mIsOpen = false;
	READ_FILE_CLOSE;
	if(mHasBaseQualityFile) QUALITY_FILE_CLOSE;
}

// enables parsing of the base quality file
void CFasta::EnableBaseQualityFile(const string& filename) {
	mHasBaseQualityFile  = true;
	mBaseQualityFilename = filename;
}

// loads the next read from the FASTA file
bool CFasta::LoadNextMate(FastaTags& ft, Mosaik::Mate& m) {

	if(!mIsOpen) {
		printf("ERROR: LoadNextMate was called before the FASTA files were opened.\n");
		exit(1);
	}

	// ======================
	// retrieve the read name
	// ======================

	READ_FILE_GETS;
	if(READ_FILE_EOF) return false;

	// sanity check
	if(mBaseBuffer[0] != '>') {
		printf("ERROR: Expected a '>' in the FASTA header, found '%c'.\n", mBaseBuffer[0]);
		exit(1);
	}

	CRegexUtilities::ExtractSequenceName(mBaseBuffer, ft.Name);
	CRegexUtilities::ExtractSpecies(mBaseBuffer, ft.Species);
	CRegexUtilities::ExtractGenomeAssemblyID(mBaseBuffer, ft.GenomeAssemblyID);
	CRegexUtilities::ExtractURI(mBaseBuffer, ft.URI);

	if(mHasBaseQualityFile) {
		QUALITY_FILE_GETS;
		if(QUALITY_FILE_EOF) return false;

		// sanity check
		if(mQualityBuffer[0] != '>') {
			printf("ERROR: Expected a '>' in the FASTA header, found '%c'.\n", mQualityBuffer[0]);
			exit(1);
		}

		CMosaikString qualityReadName;
		CRegexUtilities::ExtractSequenceName(mQualityBuffer, qualityReadName);

		// sanity check
		if(ft.Name != qualityReadName) {
			printf("ERROR: The read names in the read file and base quality file didn't match.\n");
			exit(1);
		}
	}

	// ==================
	// retrieve the bases
	// ==================

	ostringstream sb;
	while(true) {
		char ch = READ_FILE_GETC;
		READ_FILE_UNGETC(ch);
		if((ch == '>') || READ_FILE_EOF) break;		
		READ_FILE_GETS;
		CSequenceUtilities::Chomp(mBaseBuffer);
		sb << mBaseBuffer;
	}

	m.Bases = sb.str().c_str();
	sb.str("");

	// ======================
	// retrieve the qualities
	// ======================

	if(mHasBaseQualityFile) {
		while(true) {
			char ch = QUALITY_FILE_GETC;
			QUALITY_FILE_UNGETC(ch);
			if((ch == '>') || QUALITY_FILE_EOF) break;
			QUALITY_FILE_GETS;
			CSequenceUtilities::Chomp(mQualityBuffer);
			sb << mQualityBuffer;
		}

		string qualities = sb.str();
		CRegexUtilities::ConvertQualities(qualities, m.Qualities);

	} else m.Qualities.Fill(mAssignedBaseQuality, m.Bases.Length());

	m.Bases.Uppercase();

	return true;
}

// opens the FASTA file
void CFasta::Open(const string& filename) {

	// ===================================
	// check the compression state (bases)
	// ===================================

	mAreBasesCompressed = false;

	FILE* checkStream = NULL;
	if(fopen_s(&checkStream, filename.c_str(), "rb") != 0) {
		printf("ERROR: Unable to open the read FASTA file.\n");
		exit(1);
	}

	const unsigned short GZIP_MAGIC_NUMBER = 0x8b1f;
	unsigned short magicNumber = 0;
	fread((char*)&magicNumber, SIZEOF_SHORT, 1, checkStream);
	fclose(checkStream);

	if(magicNumber == GZIP_MAGIC_NUMBER) mAreBasesCompressed = true;

	// =======================================
	// check the compression state (qualities)
	// =======================================

	mAreBaseQualitiesCompressed = false;

	if(mHasBaseQualityFile) {
		if(fopen_s(&checkStream, mBaseQualityFilename.c_str(), "rb") != 0) {
			printf("ERROR: Unable to open the base quality FASTA file.\n");
			exit(1);
		}

		fread((char*)&magicNumber, SIZEOF_SHORT, 1, checkStream);
		fclose(checkStream);

		if(magicNumber == GZIP_MAGIC_NUMBER) mAreBaseQualitiesCompressed = true;
	}

	// ====================
	// open the FASTA files
	// ====================

	READ_FILE_OPEN(filename.c_str());

	if(!READ_STREAM) {
		printf("ERROR: Unable to open the FASTA read file (%s) for reading.\n", filename.c_str());
		exit(1);
	}

	if(mHasBaseQualityFile) {
		QUALITY_FILE_OPEN(mBaseQualityFilename.c_str());

		if(!QUALITY_STREAM) {
			printf("ERROR: Unable to open the FASTA base quality file (%s) for reading.\n", mBaseQualityFilename.c_str());
			exit(1);
		}
	}

	mIsOpen = true;

	// ============================================
	// skip over comments (found in .csfasta files)
	// ============================================

	// allocate our buffers
	CMemoryUtilities::CheckBufferSize(mBaseBuffer,    mBaseBufferLen,    1024);
	CMemoryUtilities::CheckBufferSize(mQualityBuffer, mQualityBufferLen, 1024);

	// skip over the read comments
	while(true) {
		char ch = READ_FILE_GETC;
		if(READ_FILE_EOF) break;
		READ_FILE_UNGETC(ch);
		if(ch == '#') READ_FILE_GETS;
		else break;
	}

	// skip over the base quality comments
	if(mHasBaseQualityFile) {
		while(true) {
			char ch = QUALITY_FILE_GETC;
			if(QUALITY_FILE_EOF) break;
			QUALITY_FILE_UNGETC(ch);
			if(ch == '#') QUALITY_FILE_GETS;
			else break;
		}
	}

	// record the data offsets
	mReadDataBaseOffset = READ_FILE_OFFSET;
	if(mHasBaseQualityFile) mReadDataQualityOffset = QUALITY_FILE_OFFSET;
}

// sets the file pointer to the beginning of the read data
void CFasta::Rewind(void) {
	READ_FILE_REWIND;
	if(mHasBaseQualityFile) QUALITY_FILE_REWIND;
}

// sets the assigned base quality if a base quality file is not specified
void CFasta::SetAssignedBaseQuality(const unsigned char baseQuality) {
	mAssignedBaseQuality = baseQuality;
}
