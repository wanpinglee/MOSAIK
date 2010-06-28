// ***************************************************************************
// CSingleEndSort - sorts single-end reads by reference sequence position.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "SingleEndSort.h"

// constructor
CSingleEndSort::CSingleEndSort(const unsigned int numCachedAlignments)
: mNumCachedAlignments(numCachedAlignments)
, mBuffer(NULL)
, mBufferLen(0)
, mSortNonUniqueMates(false)
, mRemoveDuplicates(false)
, mRenameReads(false)
{
}

// destructor
CSingleEndSort::~CSingleEndSort(void) {}

// retrieves an alignment from the specified temporary file and adds it to the specified vector
void CSingleEndSort::AddAlignment(FILE* tempFile, const unsigned int owner, list<Alignment>& alignments) {
	Alignment al;
	if(GetAlignment(tempFile, owner, al)) alignments.push_back(al);
}

// enables consed renaming
void CSingleEndSort::EnableConsedRenaming(void) {
	mRenameReads = true;
}

// enables duplicate read filtering
void CSingleEndSort::EnableDuplicateFiltering(const string& duplicateDirectory) {
	mRemoveDuplicates   = true;
	mDuplicateDirectory = duplicateDirectory;
}

// processes multiply aligned reads
void CSingleEndSort::EnableNonUniqueMode(void) {
	mSortNonUniqueMates = true;
}

// retrieves an alignment from the specified temporary file
bool CSingleEndSort::GetAlignment(FILE* tempFile, const unsigned int owner, Alignment& al) {

	// assign the owner
	al.Owner = owner;

	// get the entry size
	unsigned short entrySize = 0;
	fread((char*)&entrySize, SIZEOF_SHORT, 1, tempFile);
	if(feof(tempFile)) return false;

	// get the entry
	CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, entrySize);
	fread(mBuffer, entrySize, 1, tempFile);
	if(feof(tempFile)) return false;
	unsigned char* pBuffer = mBuffer;

	// get the alignment name
	const unsigned char readNameLen = *pBuffer;
	pBuffer++;

	al.Name.Copy((const char*)pBuffer, readNameLen);
	pBuffer += readNameLen;

	// get the read group code
	memcpy((char*)&al.ReadGroupCode, pBuffer, SIZEOF_INT);
	pBuffer += SIZEOF_INT;

	// get the status
	const unsigned char status = *pBuffer;
	pBuffer++;

	const bool isLongRead = ((status & SE_IS_LONG_READ)      != 0 ? true : false);
	al.IsReverseStrand    = ((status & SE_IS_REVERSE_STRAND) != 0 ? true : false);

	// get the reference sequence start position
	memcpy((char*)&al.NumMismatches, pBuffer, SIZEOF_SHORT);
	pBuffer += SIZEOF_SHORT;

	// get the reference sequence start position
	memcpy((char*)&al.ReferenceBegin, pBuffer, SIZEOF_INT);
	pBuffer += SIZEOF_INT;

	// get the reference sequence end position
	memcpy((char*)&al.ReferenceEnd, pBuffer, SIZEOF_INT);
	pBuffer += SIZEOF_INT;

	// get the reference sequence index
	memcpy((char*)&al.ReferenceIndex, pBuffer, SIZEOF_INT);
	pBuffer += SIZEOF_INT;

	// get the alignment quality
	al.Quality = *pBuffer;
	pBuffer++;

	unsigned short pairwiseLength = 0;
	if(isLongRead) {

		// get the pairwise length
		memcpy((char*)&pairwiseLength, pBuffer, SIZEOF_SHORT);
		pBuffer += SIZEOF_SHORT;

		// get the query begin
		memcpy((char*)&al.QueryBegin, pBuffer, SIZEOF_SHORT);
		pBuffer += SIZEOF_SHORT;

		// get the query end
		memcpy((char*)&al.QueryEnd, pBuffer, SIZEOF_SHORT);
		pBuffer += SIZEOF_SHORT;

	} else {

		// get the pairwise length
		pairwiseLength = *pBuffer;
		pBuffer++;

		// get the query begin
		al.QueryBegin = *pBuffer;
		pBuffer++;

		// get the query end
		al.QueryEnd = *pBuffer;
		pBuffer++;
	}

	// get the pairwise reference bases
	al.Reference.Copy((const char*)pBuffer, pairwiseLength);
	pBuffer += pairwiseLength;

	// get the pairwise query bases
	al.Query.Copy((const char*)pBuffer, pairwiseLength);
	pBuffer += pairwiseLength;
	RecordReferenceGaps(al);

	// get the pairwise query base qualities
	const unsigned short bqLength = al.QueryEnd - al.QueryBegin + 1;
	al.BaseQualities.Copy((const char*)pBuffer, bqLength);
	pBuffer += bqLength;

	//// DEBUG
	//al.BaseQualities.Increment(33);
	//printf("Name: %s, reference index: %2u, orientation: %c, alignment quality: %2u\n", 
	//	al.Name.CData(), al.ReferenceIndex, (al.IsReverseStrand ? 'R' : 'F'), al.Quality);
	//printf("alternate quality: %2u, read group code: %u\n", al.AlternateQuality, al.ReadGroupCode);
	//printf("%8u %s %u\n", al.ReferenceBegin, al.Reference.CData(), al.ReferenceEnd);
	//printf("%8u %s %u\n", al.QueryBegin, al.Query.CData(), al.QueryEnd);
	//printf("         %s\n\n", al.BaseQualities.CData());
	//exit(1);

	return true;
}

// records the observed gaps in the specified reference 
void CSingleEndSort::RecordReferenceGaps(Alignment& al) {

	//if(al.IsReverseStrand) CorrectHomopolymerGapOrder(al);

	// localize some data
	const char* pReference         = al.Reference.CData();
	const unsigned int refLength   = al.Reference.Length();

	unordered_map<unsigned int, unsigned short>* pHashMap = &mRefGapVector[al.ReferenceIndex];

	// initialize
	unsigned short gapLength   = 0;
	unsigned int checkPosition = 0;
	unsigned int gapPosition   = 0;

	// find the gaps
	unsigned int ungappedRefPos = al.ReferenceBegin;
	for(unsigned int i = 0; i < refLength; i++) {

		// check if we have a gap
		if(pReference[i] == GAP) {

			// find the gap length
			gapPosition   = ungappedRefPos - 1;
			gapLength     = 0;
			checkPosition = i;
			while((checkPosition < refLength) && (pReference[checkPosition++] == GAP)) gapLength++;

			// find the gap position in the hash map
			mRefGapIter = pHashMap->find(gapPosition);

			// add or modify the gap length in the hash map
			if(mRefGapIter != pHashMap->end()) {
				const unsigned short oldGapLength = mRefGapIter->second;
				mRefGapIter->second = max(oldGapLength, gapLength);
			} else (*pHashMap)[gapPosition] = gapLength;

			// jump over the gap
			if(gapLength > 1) i += gapLength - 1;

		} else { // a gap was not found

			// increment our ungapped reference position
			ungappedRefPos++;
		}
	}
}

// sorts the input alignments and saves them to the output file
void CSingleEndSort::SaveAlignmentsOrderedByPosition(const string& inputFilename, const string& outputFilename) {

	// -----------------------------------------
	// sort and serialize the alignments to disk
	// -----------------------------------------

	// initialize
	vector<Alignment>::const_iterator mate1Iter;
	uint64_t currentRead                             = 0;
	uint64_t numSerializedAlignments                 = 0;
	uint64_t numParsedAlignments                     = 0;
	uint64_t numSerializedAlignmentsDuplicateRemoved = 0;

	uint64_t numNonUniqueAlignments                  = 0;
	uint64_t numNonUniqueReads                       = 0;
	uint64_t numUniqueReads                          = 0;

	list<Alignment> alignmentCache;
	alignmentCache.resize(mNumCachedAlignments);
	list<Alignment>::iterator acIter = alignmentCache.begin();
	unsigned int numCachedEntries    = 0;

	char readNameBuffer[READ_NAME_BUFFER_SIZE];

	MosaikReadFormat::CAlignmentReader reader;
	reader.Open(inputFilename);

	vector<MosaikReadFormat::ReadGroup> readGroups;
	reader.GetReadGroups(readGroups);

	vector<ReferenceSequence>* pReferenceSequences = reader.GetReferenceSequences();
	mRefGapVector.resize(pReferenceSequences->size());
	const uint64_t numReads = reader.GetNumReads();

	// ================================================================================
	// gather all fragments that belonging to this read group from the fragment library
	// ================================================================================

	set<string> prunedReadNames;
	set<string>::const_iterator prnIter;

	if(mRemoveDuplicates) {

		// derive the database filename
		string databasePath = mDuplicateDirectory + readGroups[0].LibraryName + ".db";
		CFileUtilities::CheckFile(databasePath.c_str(), true);

		printf("- retrieving read names from the duplicate database... ");
		fflush(stdout);

		// open the database
		char* errorMessage = NULL;
		char** sqlResults  = NULL;
		sqlite3* db        = NULL;
		int numRows        = 0;
		int numColumns     = 0;

		char* sqlBuffer = new char[SQL_BUFFER_SIZE];

		if(sqlite3_open(databasePath.c_str(), &db)) {
			printf("ERROR: Unable to open the database (%s) for writing.\n", databasePath.c_str());
			printf("       error message: %s\n", sqlite3_errmsg(db));
			exit(1);
		}

		// set exclusive locking mode
		if(sqlite3_exec(db, "PRAGMA locking_mode = EXCLUSIVE;", NULL, 0, &errorMessage) != SQLITE_OK) {
			printf("ERROR: Unable to set exclusive locking mode on the database (%s).\n", databasePath.c_str());
			printf("       error message: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		// disable the synchronous mode
		if(sqlite3_exec(db, "PRAGMA synchronous = OFF;", NULL, 0, &errorMessage) != SQLITE_OK) {
			printf("ERROR: Unable to disable synchronous mode on the database (%s).\n", databasePath.c_str());
			printf("       error message: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		// use memory journaling
		if(sqlite3_exec(db, "PRAGMA journal_mode = MEMORY;", NULL, 0, &errorMessage) != SQLITE_OK) {
			printf("ERROR: Unable to use memory journaling on the database (%s).\n", databasePath.c_str());
			printf("       error message: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		// =============================================
		// look up the database ID for our read group ID
		// =============================================

		const string readGroupID = readGroups[0].ReadGroupID;

		sprintf_s(sqlBuffer, SQL_BUFFER_SIZE, "SELECT ID FROM ReadGroups WHERE Name='%s';", readGroupID.c_str());

		if(sqlite3_get_table(db, sqlBuffer, &sqlResults, &numRows, &numColumns, &errorMessage) != SQLITE_OK) {
			printf("ERROR: The SQL query resulted in the following error: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		// sanity check: ensure that we have an answer with 1 column
		if(numColumns != 1) {
			printf("ERROR: Expected one column when receiving the read group database code.\n");
			printf("rows: %u, columns: %u\n", numRows, numColumns);
			exit(1);
		}

		const string readGroupCode = sqlResults[1];

		sqlite3_free_table(sqlResults);

		// ============================
		// retrieve the best read names
		// ============================

		sprintf_s(sqlBuffer, SQL_BUFFER_SIZE, "SELECT BestReadName FROM SingleFragments WHERE BestReadGroupID=%s;", readGroupCode.c_str());

		if(sqlite3_get_table(db, sqlBuffer, &sqlResults, &numRows, &numColumns, &errorMessage) != SQLITE_OK) {
			printf("ERROR: The SQL query resulted in the following error: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		// sanity check: ensure that we have an answer with 1 column
		if(numColumns != 1) {
			printf("ERROR: Expected one columns when receiving best read names.\n");
			printf("rows: %u, columns: %u\n", numRows, numColumns);
			exit(1);
		}

		// store the read names
		for(int i = 1; i <= numRows; i++) prunedReadNames.insert(sqlResults[i]);

		// clean up
		sqlite3_free_table(sqlResults);
		sqlite3_close(db);

		printf("%u read names loaded.\n\n", (unsigned int)prunedReadNames.size());
	}

	// ===============================
	// sort and serialize unique reads
	// ===============================

	CConsole::Heading();
	printf("- phase 1 of 2: serialize alignments:\n");
	CConsole::Reset();

	CProgressBar<uint64_t>::StartThread(&currentRead, 0, numReads, "reads");

	Mosaik::AlignedRead ar;	
	while(reader.LoadNextRead(ar)) {

		numParsedAlignments += ar.Mate1Alignments.size();

		// figure out if our mate is unique
		const bool isMate1Unique = (ar.Mate1Alignments.size() == 1 ? true : false);

		if(!isMate1Unique) {
			numNonUniqueReads++;
			numNonUniqueAlignments += ar.Mate1Alignments.size();
		} else numUniqueReads++;

		if(!mSortNonUniqueMates && !isMate1Unique) {
			currentRead++;
			continue;
		}

		// skip the read if using the duplicate database
		if(mRemoveDuplicates) {
			const string readName = ar.Name.CData();
			prnIter = prunedReadNames.find(readName);
			if(prnIter == prunedReadNames.end()) {
				numSerializedAlignmentsDuplicateRemoved++;
				currentRead++;
				continue;
			}
		}

		// process each alignment
		unsigned int alignmentNum = 1;
		for(mate1Iter = ar.Mate1Alignments.begin(); mate1Iter != ar.Mate1Alignments.end(); mate1Iter++, alignmentNum++) {

			// copy the mate1 alignment
			*acIter = *mate1Iter;
			acIter->Name = ar.Name;
			acIter->ReadGroupCode = ar.ReadGroupCode;
			if(mRenameReads && mSortNonUniqueMates) {
				sprintf_s(readNameBuffer, READ_NAME_BUFFER_SIZE, ".%u", alignmentNum);
				acIter->Name.Append(readNameBuffer);
			}
			numCachedEntries++;
			acIter++;

			// serialize the alignment cache
			if(acIter == alignmentCache.end()) {
				numSerializedAlignments += Serialize(alignmentCache, numCachedEntries);
				acIter = alignmentCache.begin();
				numCachedEntries = 0;
			}
		}

		// increment the read counter
		currentRead++;
	}

	// wait for the progress bar to end
	CProgressBar<uint64_t>::WaitThread();

	// close the input alignment archive
	reader.Close();

	// ====================================================
	// serialize the remaining reads in the alignment cache
	// ====================================================

	if(acIter != alignmentCache.begin()) numSerializedAlignments += Serialize(alignmentCache, numCachedEntries);

	// ------------------------------------------------------------
	// consolidate the sorted files and create a new alignment file
	// ------------------------------------------------------------

	CConsole::Heading();
	printf("\n- phase 2 of 2: restitch serialized alignments:\n");
	CConsole::Reset();

	// open our output file
	MosaikReadFormat::CAlignmentWriter aw;
	aw.Open(outputFilename, *pReferenceSequences, readGroups, AS_SORTED_ALIGNMENT);

	// allocate the file stream array
	const unsigned int numTempFiles = (unsigned int)mTempFiles.size();
	FILE** tempFile = new FILE*[numTempFiles];

	// sanity check: check the number of temp files
	if(numTempFiles > 65535) {
		printf("ERROR: More than 65535 temporary files were produced during sorting. The owner variable needs to be redimensioned.\n");
		exit(1);
	}

	// allocate the alignments list
	list<Alignment>::iterator bestIter, nextBestIter;
	list<Alignment> alignments;
	Alignment al, *pBestAlignment = NULL;
	uint64_t numSavedAlignments = 0;

	// open the file streams
	for(unsigned int i = 0; i < mTempFiles.size(); i++) {
		if(fopen_s(&tempFile[i], mTempFiles[i].c_str(), "rb") != 0) {
			cout << "ERROR: Unable to open temporary file (" << mTempFiles[i] << ") for reading." << endl;
			exit(1);
		}
	}

	// fill the alignment list
	for(unsigned int i = 0; i < numTempFiles; i++) AddAlignment(tempFile[i], i, alignments);

	// show the progress bar
	CProgressBar<uint64_t>::StartThread(&numSavedAlignments, 0, numSerializedAlignments, "alignments");

	// keep processing until only one active file remains
	unsigned short bestOwner = 0;
	while(alignments.size() > 1) {

		// sort the alignment list
		alignments.sort();

		// grab the two best alignments
		bestIter     = alignments.begin();
		nextBestIter = alignments.begin();
		nextBestIter++;

		// save the best alignment
		pBestAlignment = &(*bestIter);
		aw.SaveAlignment(pBestAlignment);
		numSavedAlignments++;

		// remove the best alignment
		bestOwner = bestIter->Owner;
		alignments.pop_front();

		// grab another alignment from the best alignment's file
		if(!GetAlignment(tempFile[bestOwner], bestOwner, al)) continue;

		// save these alignments as long as they are better than the next best
		bool isFileEmpty = false;
		while(al < *nextBestIter) {
			aw.SaveAlignment(&al);
			numSavedAlignments++;
			isFileEmpty = !GetAlignment(tempFile[bestOwner], bestOwner, al);
			if(isFileEmpty) break;
		}

		// add the new alignment to the sorting list
		if(!isFileEmpty) alignments.push_back(al);
	}

	// save the last alignment in our vector
	pBestAlignment = &(*alignments.begin());
	bestOwner = pBestAlignment->Owner;
	aw.SaveAlignment(pBestAlignment);
	numSavedAlignments++;

	// save the alignments in the remaining file
	while(GetAlignment(tempFile[bestOwner], bestOwner, al)) {
		aw.SaveAlignment(&al);
		numSavedAlignments++;
	}

	// wait for the progress bar to end
	CProgressBar<uint64_t>::WaitThread();

	// close the file streams
	aw.SetReferenceGaps(&mRefGapVector);
	aw.Close();
	for(unsigned int i = 0; i < numTempFiles; i++) fclose(tempFile[i]);

	// clean up
	if(tempFile) delete [] tempFile;

	// ======================
	// display our statistics
	// ======================

	const uint64_t numTotalReads              = numUniqueReads + numNonUniqueReads;
	const uint64_t numTotalAlignments         = numUniqueReads + numNonUniqueAlignments;
	const uint64_t numTotalResolvedAlignments = numSerializedAlignments;

	CConsole::Heading(); printf("\nSingle-end read statistics:\n"); CConsole::Reset();

	printf("======================================================\n");
	printf("                     reads              alignments   \n");
	printf("------------------------------------------------------\n");

	printf("# non-unique: %9llu (%4.1f %%)   %9llu (%4.1f %%)\n", (unsigned long long)numNonUniqueReads, (numNonUniqueReads / (double)numTotalReads) * 100.0, (unsigned long long)numNonUniqueAlignments, (numNonUniqueAlignments / (double)numTotalAlignments) * 100.0);
	printf("# unique:     %9llu (%4.1f %%)   %9llu (%4.1f %%)\n", (unsigned long long)numUniqueReads,    (numUniqueReads    / (double)numTotalReads) * 100.0, (unsigned long long)numUniqueReads,         (numUniqueReads         / (double)numTotalAlignments) * 100.0);

	printf("------------------------------------------------------\n");
	printf("total:        %9llu            %9llu\n\n", (unsigned long long)numTotalReads, (unsigned long long)numTotalAlignments);

	if(mRemoveDuplicates) {
		printf("removed duplicates:            %9llu (%4.1f %%)\n",   (unsigned long long)numSerializedAlignmentsDuplicateRemoved, (numSerializedAlignmentsDuplicateRemoved / (double)numTotalAlignments) * 100.0);
		printf("total written:                 %9llu (%4.1f %%)\n\n", (unsigned long long)numTotalResolvedAlignments,              (numTotalResolvedAlignments              / (double)numTotalAlignments) * 100.0);
	}
}

// serializes the specified list to a temporary file
uint64_t CSingleEndSort::Serialize(list<Alignment>& alignmentCache, const unsigned int numEntries) {

	uint64_t numSerializedAlignments = 0;

	// return if we have nothing to serialize
	if(numEntries == 0) return numSerializedAlignments;

	// remove reads if the alignment cache is not full
	if(numEntries != mNumCachedAlignments) alignmentCache.resize(numEntries);

	// sort if have more than one read
	if(numEntries > 1) alignmentCache.sort();

	// retrieve a temporary filename
	string tempFilename;
	CFileUtilities::GetTempFilename(tempFilename);
	mTempFiles.push_back(tempFilename);

	// open the temporary file
	FILE* temp = NULL;
	fopen_s(&temp, tempFilename.c_str(), "wb");

	if(!temp) {
		cout << "ERROR: Unable to open temporary file (" << tempFilename << ") for writing." << endl;
		exit(1);
	}

	// ============================
	// serialize data to our buffer
	// ============================

	list<Alignment>::const_iterator alIter;
	for(alIter = alignmentCache.begin(); alIter != alignmentCache.end(); alIter++) {

		// calculate the entry size
		const unsigned char readNameLen     = (unsigned char)alIter->Name.Length();
		const unsigned short pairwiseLength = (unsigned short)alIter->Reference.Length();
		const unsigned short bqLength       = alIter->QueryEnd - alIter->QueryBegin + 1;

		const bool isLongRead = ((alIter->QueryEnd > 255) || (pairwiseLength > 255) ? true : false);
		const unsigned int requestedBytes = (SIZEOF_INT * 4) + (SIZEOF_SHORT * 2) + 6 + readNameLen + pairwiseLength + pairwiseLength + bqLength + (isLongRead ? 3 : 0);

		// check our memory allocation
		if(requestedBytes > mBufferLen) CMemoryUtilities::CheckBufferSize(mBuffer, mBufferLen, requestedBytes);

		// store the entry size
		unsigned int bufferOffset = 0;
		const unsigned short entrySize = requestedBytes - SIZEOF_SHORT;
		memcpy(mBuffer + bufferOffset, (char*)&entrySize, SIZEOF_SHORT);
		bufferOffset += SIZEOF_SHORT;

		// store the read name
		mBuffer[bufferOffset++] = readNameLen;
		memcpy(mBuffer + bufferOffset, alIter->Name.CData(), readNameLen);
		bufferOffset += readNameLen;

		// store the read group code
		memcpy(mBuffer + bufferOffset, (char*)&alIter->ReadGroupCode, SIZEOF_INT);
		bufferOffset += SIZEOF_INT;

		// store the status
		unsigned char status = 0;
		if(isLongRead)              status |= SE_IS_LONG_READ;
		if(alIter->IsReverseStrand) status |= SE_IS_REVERSE_STRAND;

		mBuffer[bufferOffset++] = status;

		// store the reference sequence start position
		memcpy(mBuffer + bufferOffset, (char*)&alIter->NumMismatches, SIZEOF_SHORT);
		bufferOffset += SIZEOF_SHORT;

		// store the reference sequence start position
		memcpy(mBuffer + bufferOffset, (char*)&alIter->ReferenceBegin, SIZEOF_INT);
		bufferOffset += SIZEOF_INT;

		// store the reference sequence end position
		memcpy(mBuffer + bufferOffset, (char*)&alIter->ReferenceEnd, SIZEOF_INT);
		bufferOffset += SIZEOF_INT;

		// store the reference sequence index
		memcpy(mBuffer + bufferOffset, (char*)&alIter->ReferenceIndex, SIZEOF_INT);
		bufferOffset += SIZEOF_INT;

		// store the alignment quality
		mBuffer[bufferOffset++] = alIter->Quality;

		if(isLongRead) {

			// store the pairwise length
			memcpy(mBuffer + bufferOffset, (char*)&pairwiseLength, SIZEOF_SHORT);
			bufferOffset += SIZEOF_SHORT;

			// store the query begin
			memcpy(mBuffer + bufferOffset, (char*)&alIter->QueryBegin, SIZEOF_SHORT);
			bufferOffset += SIZEOF_SHORT;

			// store the query end
			memcpy(mBuffer + bufferOffset, (char*)&alIter->QueryEnd, SIZEOF_SHORT);
			bufferOffset += SIZEOF_SHORT;

		} else {

			// store the pairwise length
			mBuffer[bufferOffset++] = (unsigned char)pairwiseLength;

			// store the query begin
			mBuffer[bufferOffset++] = (unsigned char)alIter->QueryBegin;

			// store the query end
			mBuffer[bufferOffset++] = (unsigned char)alIter->QueryEnd;
		}

		// store the pairwise reference bases
		memcpy(mBuffer + bufferOffset, alIter->Reference.CData(), pairwiseLength);
		bufferOffset += pairwiseLength;

		// store the pairwise query bases
		memcpy(mBuffer + bufferOffset, alIter->Query.CData(), pairwiseLength);
		bufferOffset += pairwiseLength;

		// store the pairwise query base qualities
		memcpy(mBuffer + bufferOffset, alIter->BaseQualities.CData(), bqLength);
		bufferOffset += bqLength;

		// sanity check
		if(bufferOffset != requestedBytes) {
			printf("ERROR: Mismatch found between the allocated buffer size and the used buffer size during serialization.\n");
			printf("       allocated buffer size: %u, used buffer size: %u\n", requestedBytes, bufferOffset);
			exit(1);
		}

		// write the buffer
		fwrite(mBuffer, bufferOffset, 1, temp);

		// increment our serialized alignments counter
		numSerializedAlignments++;
	}

	// close the temp file
	fclose(temp);

	return numSerializedAlignments;
}
