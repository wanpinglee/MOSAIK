// ***************************************************************************
// CMosaikMerge - merges two or more sorted alignment archives.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "MosaikMerge.h"

// constructor
CMosaikMerge::CMosaikMerge(const unsigned int numCachedAlignments)
: mNumCachedAlignments(numCachedAlignments)
, mBuffer(NULL)
, mBufferLen(0)
{
}

// destructor
CMosaikMerge::~CMosaikMerge(void) {
	if(mBuffer) delete [] mBuffer;

	// delete our temporary files
	for(unsigned int i = 0; i < mTempFiles.size(); i++) rm(mTempFiles[i].c_str());
}

// retrieves an alignment from the specified temporary file and adds it to the specified list
void CMosaikMerge::AddAlignment(FILE* tempFile, const unsigned int owner, list<Alignment>& alignments) {
	Alignment al;
	if(GetAlignment(tempFile, owner, al)) alignments.push_back(al);
}

// retrieves a read from the specified temporary file
bool CMosaikMerge::GetAlignment(FILE* tempFile, const unsigned int owner, Alignment& al) {

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

	const bool isLongRead  = ((status & MERGE_IS_LONG_READ)               != 0 ? true : false);
	al.IsReverseStrand     = ((status & MERGE_IS_REVERSE_COMPLEMENT)      != 0 ? true : false);
	al.IsMateReverseStrand = ((status & MERGE_IS_MATE_REVERSE_COMPLEMENT) != 0 ? true : false);

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

	// get the mate reference sequence start position
	memcpy((char*)&al.MateReferenceBegin, pBuffer, SIZEOF_INT);
	pBuffer += SIZEOF_INT;

	// get the mate reference sequence end position
	memcpy((char*)&al.MateReferenceEnd, pBuffer, SIZEOF_INT);
	pBuffer += SIZEOF_INT;

	// get the mate reference sequence index
	memcpy((char*)&al.MateReferenceIndex, pBuffer, SIZEOF_INT);
	pBuffer += SIZEOF_INT;

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
	RecordReferenceGaps(al.ReferenceIndex, al.ReferenceBegin, al.Reference);

	// get the pairwise query bases
	al.Query.Copy((const char*)pBuffer, pairwiseLength);
	pBuffer += pairwiseLength;

	// get the pairwise query base qualities
	const unsigned short bqLength = al.QueryEnd - al.QueryBegin + 1;
	al.BaseQualities.Copy((const char*)pBuffer, bqLength);
	pBuffer += bqLength;

	// DEBUG
	//al.BaseQualities.Increment(33);
	//printf("Name: %s, reference: %2u, orientation: %c, alignment quality: %2u\n", al.Name.CData(), al.ReferenceIndex, (al.IsReverseStrand ? 'R' : 'F'), al.Quality);
	//printf("%8u %s %u\n", al.ReferenceBegin, al.Reference.CData(), al.ReferenceEnd);
	//printf("%8u %s %u\n", al.QueryBegin, al.Query.CData(), al.QueryEnd);
	//printf("         %s\n\n", al.BaseQualities.CData());

	//printf("%20s %2u %9u %9u\n", al.Name.CData(), al.ReferenceIndex, al.ReferenceBegin, al.ReferenceEnd);

	return true;
}

// merges the files contained in the file vector and stores them in a specified output file
void CMosaikMerge::MergeFiles(vector<string>& fileVector, string& outputFilename) {

	// ------------------------------------------------------------------------------------------------------
	// calculate the total number of alignments, grab the read groups, and check reference sequence integrity
	// ------------------------------------------------------------------------------------------------------

	vector<MosaikReadFormat::ReadGroup> readGroups;
	vector<MosaikReadFormat::ReadGroup>::const_iterator rgIter;

	vector<ReferenceSequence> referenceSequences;
	vector<ReferenceSequence>::iterator rsIter;

	uint64_t numTotalAlignments = 0;
	vector<string>::const_iterator svIter;

	for(svIter = fileVector.begin(); svIter != fileVector.end(); svIter++) {

		MosaikReadFormat::CAlignmentReader reader;
		reader.Open(*svIter);

		// copy the read groups
		vector<MosaikReadFormat::ReadGroup> tempReadGroups;
		reader.GetReadGroups(tempReadGroups);
		for(rgIter = tempReadGroups.begin(); rgIter != tempReadGroups.end(); rgIter++) {
			readGroups.push_back(*rgIter);
		}

		// copy the reference sequences
		vector<ReferenceSequence>::const_iterator trsIter;
		vector<ReferenceSequence>* pReferenceSequences = reader.GetReferenceSequences();

		if(referenceSequences.empty()) {
			for(trsIter = pReferenceSequences->begin(); trsIter != pReferenceSequences->end(); trsIter++) {
				referenceSequences.push_back(*trsIter);
			}

		} else { // check the reference sequences

			const unsigned int numReferenceSequences     = referenceSequences.size();
			const unsigned int numTempReferenceSequences = pReferenceSequences->size();

			// check the common set of reference sequences
			const unsigned int numCommonReferenceSequences = min(numReferenceSequences, numTempReferenceSequences);

			trsIter = pReferenceSequences->begin();
			rsIter  = referenceSequences.begin();
			for(unsigned int i = 0; i < numCommonReferenceSequences; ++i, ++trsIter, ++rsIter) {
				if(trsIter->MD5 != rsIter->MD5) {
					printf("ERROR: Found an inconsistency in the reference sequences. Merged alignment archives must share the same subset of reference sequences.\n");
					exit(1);
				}
			}

			// add any additional reference sequences
			if(numTempReferenceSequences > numReferenceSequences) {
				const unsigned int numAdditionalReferenceSequences = numTempReferenceSequences - numReferenceSequences;
				for(unsigned int i = 0; i < numAdditionalReferenceSequences; ++i, ++trsIter) {
					referenceSequences.push_back(*trsIter);
				}
			}
		}

		numTotalAlignments += reader.GetNumReads();
		reader.Close();
	}

	// -----------------------------------------
	// sort and serialize the alignments to disk
	// -----------------------------------------

	CConsole::Heading();
	printf("- phase 1 of 2: serialize alignments:\n");
	CConsole::Reset();

	// initialize our alignment cache
	list<Alignment> alignmentCache;
	alignmentCache.resize(mNumCachedAlignments);
	list<Alignment>::iterator acIter = alignmentCache.begin();

	unsigned int numCachedEntries    = 0;
	uint64_t numSerializedAlignments = 0;
	uint64_t currentAlignment        = 0;

	AlignmentStatus as = AS_SORTED_ALIGNMENT;

	CProgressBar<uint64_t>::StartThread(&currentAlignment, 0, numTotalAlignments, "alignments");

	// open each alignment file
	for(svIter = fileVector.begin(); svIter != fileVector.end(); svIter++) {

		MosaikReadFormat::CAlignmentReader reader;
		reader.Open(*svIter);

		// get the alignment status
		as |= (reader.GetStatus() & 0xf3);

		// process all of the alignments
		while(reader.LoadNextAlignment(*acIter)) {
			currentAlignment++;
			numCachedEntries++;
			acIter++;

			if(acIter == alignmentCache.end()) {
				numSerializedAlignments += Serialize(alignmentCache, numCachedEntries);
				acIter = alignmentCache.begin();
				numCachedEntries = 0;
			}
		}

		reader.Close();
	}

	// wait for the progress bar to finish
	CProgressBar<uint64_t>::WaitThread();

	if(acIter != alignmentCache.begin()) numSerializedAlignments += Serialize(alignmentCache, numCachedEntries);
	mRefGapVector.resize(referenceSequences.size());

	// ------------------------------------------------------------
	// sort the serialized reads
	// ------------------------------------------------------------

	CConsole::Heading();
	printf("\n- phase 2 of 2: restitch serialized alignments:\n");
	CConsole::Reset();

	// open our output file
	MosaikReadFormat::CAlignmentWriter aw;
	aw.Open(outputFilename, referenceSequences, readGroups, as);

	// allocate the file stream array
	const unsigned int numTempFiles = mTempFiles.size();
	FILE** tempFile = new FILE*[numTempFiles];

	// perform a sanity check
	if(numTempFiles == 0) {
		cout << "ERROR: Unable to continue merging files. There are no stored temporary file names." << endl;
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
}

// records the observed gaps in the specified reference 
void CMosaikMerge::RecordReferenceGaps(const unsigned short refIndex, const unsigned int refBegin, const CMosaikString& refSeq) {

	// localize some data
	const char* pReference       = refSeq.CData();
	const unsigned int refLength = refSeq.Length();
	unordered_map<unsigned int, unsigned short>* pHashMap = &mRefGapVector[refIndex];

	// initialize
	unsigned short gapLength   = 0;
	unsigned int checkPosition = 0;
	unsigned int gapPosition   = 0;

	// find the gaps
	unsigned int ungappedRefPos = refBegin;
	for(unsigned int i = 0; i < refLength; i++) {

		// check if we have a gap
		if(pReference[i] == '-') {

			// find the gap length
			gapLength     = 0;
			checkPosition = i;
			while((checkPosition < refLength) && (pReference[checkPosition++] == '-')) gapLength++;

			// find the gap position in the hash map
			gapPosition = ungappedRefPos - 1;
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

// serializes the specified list to a temporary file
uint64_t CMosaikMerge::Serialize(list<Alignment>& alignmentCache, const unsigned int numEntries) {

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
		const unsigned int requestedBytes = (SIZEOF_INT * 7) + SIZEOF_SHORT + 6 + readNameLen + pairwiseLength + pairwiseLength + bqLength + (isLongRead ? 3 : 0);

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
		if(isLongRead)                  status |= MERGE_IS_LONG_READ;
		if(alIter->IsReverseStrand)     status |= MERGE_IS_REVERSE_COMPLEMENT;
		if(alIter->IsMateReverseStrand) status |= MERGE_IS_MATE_REVERSE_COMPLEMENT;

		mBuffer[bufferOffset++] = status;

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

		// store the mate reference sequence start position
		memcpy(mBuffer + bufferOffset, (char*)&alIter->MateReferenceBegin, SIZEOF_INT);
		bufferOffset += SIZEOF_INT;

		// store the mate reference sequence end position
		memcpy(mBuffer + bufferOffset, (char*)&alIter->MateReferenceEnd, SIZEOF_INT);
		bufferOffset += SIZEOF_INT;

		// store the mate reference sequence index
		memcpy(mBuffer + bufferOffset, (char*)&alIter->MateReferenceIndex, SIZEOF_INT);
		bufferOffset += SIZEOF_INT;

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
