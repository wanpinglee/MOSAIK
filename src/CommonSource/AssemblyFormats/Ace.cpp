// ***************************************************************************
// CAce - exports MOSAIK assemblies into the ACE format.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// This code is dual licenced under the GNU General Public License 2.0+ or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "Ace.h"

// constructor
CAce::CAce(unsigned char referenceBaseQuality)
: mReferenceSequenceBaseQuality(referenceBaseQuality)
{

	// set the file open flag
	mIsOpen = false;

	// get the current time
	GetTime(mTimeString);
}

// destructor
CAce::~CAce(void) {
	if(mIsOpen) Close();
}

// write a sequence to the specified file stream in rows of ACE_LINE_LENGTH bases
void CAce::ChopSequence(FILE* out, const CMosaikString& sequence) {

	// initialize
	const char* pSequence        = sequence.CData();
	const unsigned int seqLength = sequence.Length();

	const char* pSeqPtr = pSequence;
	const unsigned int numRows = (unsigned int)(seqLength / (double)ACE_LINE_LENGTH);

	for(unsigned int i = 0; i < numRows; i++) {
		fwrite(pSeqPtr, ACE_LINE_LENGTH, 1, out);
		pSeqPtr += ACE_LINE_LENGTH;
		fprintf(out, "\n");
	}

	fprintf(out, "%s\n", pSeqPtr);
}

// closes the ace file
void CAce::Close(void) {

	mIsOpen = false;

	// add the byte segment to our header stream
	fprintf(mHeaderStream, "BS 1 %u .MosaikReference\n\n", mGappedRefLength);

	// create a large buffer
	const unsigned int requestedBytes = 52428800; // 50 MB
	char* buffer = new char[requestedBytes];

	// append to the output stream
	printf("\n- appending read data to header data... ");
	fflush(stdout);

	rewind(mReadStream);

	size_t numBytesRead = 0, numBytesWritten = 0;
	while(true) {
		numBytesRead = fread(buffer, 1, requestedBytes, mReadStream);
		if(feof(mReadStream)) break;
		numBytesWritten = fwrite(buffer, 1, numBytesRead, mHeaderStream);
	}

	numBytesWritten = fwrite(buffer, 1, numBytesRead, mHeaderStream);

	printf("finished.\n");

	//// append the read stream
	//char buffer[256];
	//fseek64(mReadStream, 0, SEEK_SET);

	//while(true) {
	//	fgets(buffer, 256, mReadStream);
	//	if(feof(mReadStream)) break;
	//	fprintf(mHeaderStream, "%s", buffer);
	//}

	// clean up
	delete [] buffer;

	// close our files
	fclose(mHeaderStream);
	fclose(mReadStream);

	// delete our temporary files
	rm(mTempReadFilename.c_str());
}

// retrieves the current time (asctime but not modified by regional settings)
#define ACE_TIME_BUFFER_LEN 25
void CAce::GetTime(string& timeString) {

	// initialize
	struct tm local_tm;
	char timeBuffer[ACE_TIME_BUFFER_LEN];

	const char* WEEKDAYS[] = { "Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat" };
	const char* MONTHS[]   = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };

	// get the current time
	time_t local;
	time(&local);

	// get the local time
	localtime_s(&local_tm, &local);

	// Thu Jul 27 15:33:47 2000
	sprintf_s(timeBuffer, ACE_TIME_BUFFER_LEN, "%s %s %2u %02u:%02u:%02u %4u", WEEKDAYS[local_tm.tm_wday], MONTHS[local_tm.tm_mon], 
		local_tm.tm_mday, local_tm.tm_hour, local_tm.tm_min, local_tm.tm_sec, local_tm.tm_year + 1900);

	timeString = timeBuffer;
}

// opens the ace file
void CAce::Open(const CMosaikString& filename) {

	// create two separate temporary files (header, reads)

	// open the header file
	//CFileUtilities::GetTempFilename(mTempHeaderFilename, processID++);
	if(fopen_s(&mHeaderStream, filename.CData(), "wb") != 0) {
		cout << "ERROR: Unable to open the ace header stream for writing." << endl;
		exit(1);
	}

	// open the reads file
	CFileUtilities::GetTempFilename(mTempReadFilename);
	if(fopen_s(&mReadStream, mTempReadFilename.c_str(), "w+b") != 0) {
		cout << "ERROR: Unable to open the ace read stream for writing." << endl;
		exit(1);
	}

	mIsOpen = true;
	mGappedRefLength  = 0;
}

// saves the specified reference sequence to the header
void CAce::SaveHeader(CMosaikString& reference, const string& referenceName, const unsigned int ungappedRefLength, const uint64_t& numSequences) {

	// write the assembly header
	// N.B. the extra sequence is the MOSAIK reference sequence
	fprintf(mHeaderStream, "AS 1 %llu\n\n", (unsigned long long)(numSequences + 1));

	// write the contig header
	mGappedRefLength = reference.Length();
	fprintf(mHeaderStream, "CO %s %u %llu 1 U\n", referenceName.c_str(), mGappedRefLength, (unsigned long long)(numSequences + 1));

	// write the reference sequence
	reference.Replace('-', '*');
	ChopSequence(mHeaderStream, reference);

	// ===========================================
	// write the reference sequence base qualities
	// ===========================================

	fprintf(mHeaderStream, "\n\nBQ\n");

	unsigned int remainingBuffer = 4 * ACE_LINE_LENGTH + 1;
	char* qualityLineBuffer = new char[remainingBuffer];
	char* pBufferPtr = qualityLineBuffer;

	for(unsigned int i = 0; i < ACE_LINE_LENGTH; i++) {
		int numBytesWritten = sprintf_s(pBufferPtr, remainingBuffer, " %u", mReferenceSequenceBaseQuality);
		pBufferPtr      += numBytesWritten;
		remainingBuffer -= numBytesWritten;
	}
	*pBufferPtr = 0;

	const unsigned int numRows = (unsigned int)(ungappedRefLength / (double)ACE_LINE_LENGTH);
	const unsigned int remainingQualities = ungappedRefLength - (numRows * ACE_LINE_LENGTH);

	for(unsigned int i = 0; i < numRows; i++) fprintf(mHeaderStream, "%s\n", qualityLineBuffer);
	for(unsigned int i = 0; i < remainingQualities; i++) fprintf(mHeaderStream, " %u", mReferenceSequenceBaseQuality);
	fprintf(mHeaderStream, "\n\n");

	// clean up
	if(qualityLineBuffer) delete [] qualityLineBuffer;

	// =======================================
	// write the reference sequence read entry
	// =======================================

	fprintf(mReadStream, "RD .MosaikReference %u 0 0\n", mGappedRefLength);
	ChopSequence(mReadStream, reference);
	fprintf(mReadStream, "\nQA 1 %u 1 %u\nDS CHROMAT_FILE: MosaikReference.scf PHD_FILE: MosaikReference.scf.phd.1 TIME: %s\n\n\n", mGappedRefLength, mGappedRefLength, mTimeString.c_str());

	// ===========================
	// write the first index entry
	// ===========================

	fprintf(mHeaderStream, "AF .MosaikReference U 1\n");
}

// saves the specified read to the index and reads files
void CAce::SaveRead(const Alignment& al, CMosaikString& gappedRead) {

	// =====================
	// write the index entry
	// =====================

	int offset = (int)mpUngap2Gap[al.ReferenceBegin] + 1;
	fprintf(mHeaderStream, "AF %s %c %d\n", al.Name.CData(), (al.IsReverseStrand ? 'C' : 'U'), offset);

	// ====================
	// write the read entry
	// ====================

	const unsigned int gappedReadLength = gappedRead.Length();
	fprintf(mReadStream, "RD %s %u 0 0\n", al.Name.CData(), gappedReadLength);
	gappedRead.Replace('-', '*');
	ChopSequence(mReadStream, gappedRead);
	fprintf(mReadStream, "\nQA 1 %u 1 %u\nDS CHROMAT_FILE: %s.scf PHD_FILE: %s.scf.phd.1 TIME: %s\n\n", gappedReadLength, gappedReadLength, al.Name.CData(), al.Name.CData(), mTimeString.c_str());
}
