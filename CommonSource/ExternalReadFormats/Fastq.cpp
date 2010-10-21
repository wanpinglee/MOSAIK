// ***************************************************************************
// CFastq - imports reads from the FASTQ file format.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "Fastq.h"

// define our uncompressed/compressed macros
#define READ_FILE_CLOSE          (mAreBasesCompressed ? gzclose(mInZStream)                           : fclose(mInStream))
#define READ_FILE_EOF            (mAreBasesCompressed ? gzeof(mInZStream)                             : feof(mInStream))
#define READ_FILE_GETC           (mAreBasesCompressed ? gzgetc(mInZStream)                            : fgetc(mInStream))
#define READ_FILE_GETS           (mAreBasesCompressed ? gzgets(mInZStream, mBuffer, mBufferLen)       : fgets(mBuffer, mBufferLen, mInStream))
#define READ_FILE_OFFSET         (mAreBasesCompressed ? gztell(mInZStream)                            : ftell64(mInStream))
#define READ_FILE_OPEN(filename) (mAreBasesCompressed ? mInZStream = gzopen(filename, "rb")           : mInStream = fopen(filename, "rb"))
#define READ_FILE_REWIND         (mAreBasesCompressed ? gzseek(mInZStream, mReadDataOffset, SEEK_SET) : fseek64(mInStream, mReadDataOffset, SEEK_SET))
#define READ_FILE_UNGETC(ch)     (mAreBasesCompressed ? gzungetc(ch, mInZStream)                      : ungetc(ch, mInStream))
#define READ_STREAM              (mAreBasesCompressed ? mInZStream                                    : mInStream)

// constructor
CFastq::CFastq(void) 
: mIsOpen(false)
, mAreBasesCompressed(false)
, mIsFastqStyleKnown(false)
, mUsingIlluminaStyle(false)
, mFastqOffset(NORMAL_FASTQ_OFFSET)
, mInStream(NULL)
, mInZStream(NULL)
, mBuffer(NULL)
, mBufferLen(0)
, mReadDataOffset(0)
{
	// calculate our base quality LUT
	for(short i = -128; i < 128; i++) {
		int bq = (int)(-10.0 * log10(1.0 / (1.0 + pow(10.0, (double)i / 10.0))));
		if(bq < 1) bq = 1;
		mIlluminaToPhredLUT[i + 128] = bq;
	}
}

// destructor
CFastq::~CFastq(void) {
	if(mIsOpen) Close();
	if(mBuffer) delete [] mBuffer;
}

// checks to see if this is truly a FASTQ file
bool CFastq::CheckFile(const string& filename, const bool showError) {

	// open the FASTQ file
	bool foundError = false;
	gzFile checkStream = gzopen(filename.c_str(), "rb");
	if(checkStream == NULL) {
		if(showError) {
			cout << "ERROR: Could not open FASTQ file (" << filename << ") when performing integrity check." << endl;
			exit(1);
		}

		foundError = true;
	}

	// retrieve the first character
	const unsigned int BUFFERLEN = 1024;
	char* buffer = new char[BUFFERLEN];

	if(!foundError) {

		// check if the FASTQ file starts with the header
		gzgets(checkStream, buffer, BUFFERLEN);
		if(buffer[0] != '@') {
			if(showError) {
				cout << "ERROR: It seems that the input file (" << filename << ") is not in FASTQ format." << endl;
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

// closes the alignment archive
void CFastq::Close(void) {
	mIsOpen = false;
	READ_FILE_CLOSE;
}

// create a temporary FASTQ which is sorted by read names
void CFastq::SortByName(const string filename) {

	if(!mIsOpen) {
		printf("ERROR: SortByName was called before the FASTQ file was opened.\n");
		exit(1);
	}

	this->Rewind();

	CMosaikString readName;
	Mosaik::Mate m;

	unsigned int cacheSize = 10000000;

	vector<string> tempFastqs;
	unsigned int readCount = 0;
	list<Mosaik::Read> cache;
	
	// partially sort the file and store them in several temp files
	while ( this->LoadNextMate(readName, m) ) {
		
		readCount++;

		Mosaik::Read r;
		r.Name  = readName;
		r.Mate1 = m;

		cache.push_back(r);
		
		// full
		if ( readCount == cacheSize ) {
			// sort by read names
			cache.sort();
			
			// retrieve a temporary filename
			string tempFilename;
			CFileUtilities::GetTempFilename(tempFilename);
			tempFastqs.push_back(tempFilename);

			ofstream cacheFile;
			cacheFile.open(tempFilename.c_str(), ios::out);
			for ( list<Mosaik::Read>::iterator ite = cache.begin(); ite != cache.end(); ite++ ) {
				cacheFile << '@' << ite->Name << endl;
				cacheFile << ite->Mate1.Bases << endl;
				ite->Mate1.Qualities.Increment(mFastqOffset);
				cacheFile << '+' << ite->Name << endl;
				cacheFile << ite->Mate1.Qualities << endl;
			}
			cacheFile.close();

			readCount = 0;
		}

	}
	// sort the last list
	if ( readCount > 0 ) {
		// sort by read names
		cache.sort();

		// retrieve a temporary filename
		string tempFilename;
		CFileUtilities::GetTempFilename(tempFilename);
		tempFastqs.push_back(tempFilename);

		ofstream cacheFile;
		cacheFile.open(tempFilename.c_str(), ios::out);
		for ( list<Mosaik::Read>::iterator ite = cache.begin(); ite != cache.end(); ite++ ) {
			cacheFile << '@' << ite->Name << endl;
			cacheFile << ite->Mate1.Bases << endl;
			ite->Mate1.Qualities.Increment(mFastqOffset);
			cacheFile << '+' << ite->Name << endl;
			cacheFile << ite->Mate1.Qualities << endl;
		}
		cacheFile.close();
	}


	// globally sort the file
	unsigned int nTemp = tempFastqs.size();
	vector<CFastq> tempReaders;
	tempReaders.resize( nTemp );
	for ( unsigned int i = 0; i < nTemp; i++ )
		tempReaders[i].Open(tempFastqs[i].c_str());

	// load the top element in each temp
	vector<Mosaik::Read> tops;
	vector<bool> dones;
	unsigned int nDone = 0;
	for ( vector<CFastq>::iterator ite = tempReaders.begin(); ite != tempReaders.end(); ite++ ) {
		Mosaik::Read r;
		if ( ite->LoadNextMate( readName, m ) ) {
			r.Name = readName;
			r.Mate1 = m;
			tops.push_back(r);
			dones.push_back(false);
		}
		else {
			r.clear();
			tops.push_back(r);
			dones.push_back(true);
		}
	}

	ofstream file;
	file.open(filename.c_str(), ios::out);

	// pick the min one
	while ( nDone != ( nTemp - 1 ) ) {
		unsigned int minId;
		minId = FindMinRead(tops);

		if ( dones[minId] ) {
			cout << "ERROR: The targeted temporary has been empty." << endl;
			exit(1);
		}

		file << '@' << tops[minId].Name << endl;
		file << tops[minId].Mate1.Bases << endl;
		tops[minId].Mate1.Qualities.Increment(mFastqOffset);
		file << '+' << tops[minId].Name << endl;
		file << tops[minId].Mate1.Qualities << endl;


		if ( tempReaders[minId].LoadNextMate( readName, m ) ) {
			tops[minId].clear();
			tops[minId].Name = readName;
			tops[minId].Mate1 = m;
		}
		else {
			tops[minId].clear();
			dones[minId] = true;
			nDone++;
			//delete temp file
			rm( tempFastqs[minId].c_str() );
		}
	}

	for ( unsigned int i = 0; i < dones.size(); i++ ) {
		if ( !dones[i] ) {
			// dump the last one on tops
			file << '@' << tops[i].Name << endl;
			file << tops[i].Mate1.Bases << endl;
			tops[i].Mate1.Qualities.Increment(mFastqOffset);
			file << '+' << tops[i].Name << endl;
			file << tops[i].Mate1.Qualities << endl;

			// dump the remaining reads in temp file
			while( tempReaders[i].LoadNextMate( readName, m ) ) {
				file << '@' << readName.CData() << endl;
				file << m.Bases << endl;
				m.Qualities.Increment(mFastqOffset);
				file << '+' << readName.CData() << endl;
				file << m.Qualities << endl;
			}
			// delete temp file
			rm( tempFastqs[i].c_str() );
		}
	}

	file.close();

}


// given a vector containing reads, find the min read and return the vector id
inline unsigned int CFastq::FindMinRead ( vector<Mosaik::Read>& tops ) {
	unsigned int minId = 0;
	for ( unsigned int i = 1; i < tops.size(); i++ ) {
		if ( tops[i].Name.empty() )
			continue;
		else {
			if ( tops[minId].Name.empty() )
				minId = i;
			else {
				if ( tops[i] < tops[minId] )
					minId = i;
			}
		}
	}

	return minId;

}

// loads the next read from the FASTQ file
bool CFastq::LoadNextMate(CMosaikString& readName, Mosaik::Mate& m) {

	if(!mIsOpen) {
		printf("ERROR: LoadNextMate was called before the FASTQ file was opened.\n");
		exit(1);
	}

	// ==============================
	// retrieve the read name (bases)
	// ==============================

	READ_FILE_GETS;
	if(READ_FILE_EOF) return false;

	// sanity check
	if(mBuffer[0] != '@') {
		printf("ERROR: Expected a '@' in the FASTQ header, found '%c'.\n", mBuffer[0]);
		exit(1);
	}

	CRegexUtilities::ExtractSequenceName(mBuffer, readName);

	// ==================
	// retrieve the bases
	// ==================

	ostringstream sb;
	while(true) {
		char ch = READ_FILE_GETC;
		READ_FILE_UNGETC(ch);
		if((ch == '+') || READ_FILE_EOF) break;
		READ_FILE_GETS;
		CSequenceUtilities::Chomp(mBuffer);
		sb << mBuffer;
	}

	m.Bases = sb.str().c_str();
	const unsigned int numBases = m.Bases.Length();
	sb.str("");

	// ==================================
	// retrieve the read name (qualities)
	// ==================================

	READ_FILE_GETS;

	// ======================
	// retrieve the qualities
	// ======================

	unsigned int bufferLen = 0;
	while(true) {
		char ch = READ_FILE_GETC;
		READ_FILE_UNGETC(ch);
		if(READ_FILE_EOF) break;
		READ_FILE_GETS;
		CSequenceUtilities::Chomp(mBuffer);
		bufferLen += strlen(mBuffer);
		sb << mBuffer;
		if(bufferLen >= numBases) break;
	}

	m.Qualities = sb.str().c_str();

	// sanity check
	if(m.Qualities.Length() != m.Bases.Length()) {
		printf("ERROR: The number of qualities (%u) do not match the number of bases (%u) in %s.\n", m.Qualities.Length(), m.Bases.Length(), readName.CData());
		exit(1);
	}

	// determine the FASTQ style
	// what does this mean? October 19th, 2010
	//if(!mIsFastqStyleKnown) {
	//	char firstBase = m.Qualities[0];
	//
	//	if(firstBase <= 72) {
	//		mFastqOffset = NORMAL_FASTQ_OFFSET;
	//		mUsingIlluminaStyle = false;
	//	} else { 
	//		mFastqOffset = ILLUMINA_FASTQ_OFFSET;
	//		mUsingIlluminaStyle = true;
	//	}
	//
	//	mIsFastqStyleKnown = true;
	//}
	mFastqOffset = NORMAL_FASTQ_OFFSET;
	
	
	m.Bases.Uppercase();
	m.Qualities.Decrement(mFastqOffset);

	// convert the base qualities
	if(mUsingIlluminaStyle) {
		for(unsigned int i = 0; i < m.Qualities.Length(); i++) 
			m.Qualities[i] = mIlluminaToPhredLUT[m.Qualities[i] + 128];
	}

	return true;
}

// opens the alignment archive
void CFastq::Open(const string& filename) {

	// ===========================
	// check the compression state
	// ===========================

	mAreBasesCompressed = false;

	FILE* checkStream = NULL;
	if(fopen_s(&checkStream, filename.c_str(), "rb") != 0) {
		printf("ERROR: Unable to open the read FASTQ file, filename:%s.\n", filename.c_str());
		exit(1);
	}

	const unsigned short GZIP_MAGIC_NUMBER = 0x8b1f;
	unsigned short magicNumber = 0;
	fread((char*)&magicNumber, SIZEOF_SHORT, 1, checkStream);
	fclose(checkStream);

	if(magicNumber == GZIP_MAGIC_NUMBER) mAreBasesCompressed = true;

	// ===================
	// open the FASTQ file
	// ===================

	READ_FILE_OPEN(filename.c_str());

	if(!READ_STREAM) {
		printf("ERROR: Unable to open the FASTQ file (%s) for reading.\n", filename.c_str());
		exit(1);
	}

	mIsOpen = true;

	// ==============
	// initialization
	// ==============

	CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, 5120);

	// record the data offsets
	mReadDataOffset = READ_FILE_OFFSET;

	// reset the FASTQ style
	mIsFastqStyleKnown  = false;
	mUsingIlluminaStyle = false;
}

// sets the file pointer to the beginning of the read data
void CFastq::Rewind(void) {
	READ_FILE_REWIND;
}

// sets the BQ offset
void CFastq::SetOffset(const unsigned char offset) {
	mFastqOffset       = offset;
	mIsFastqStyleKnown = true;

	if(mFastqOffset == ILLUMINA_FASTQ_OFFSET) mUsingIlluminaStyle = true;
	else mUsingIlluminaStyle = false;
}
