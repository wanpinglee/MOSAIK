// ***************************************************************************
// CPairedEndSort - resolves paired-end reads and creates an alignment archive
//                  sorted by reference sequence position.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "PairedEndSort.h"

// constructor
CPairedEndSort::CPairedEndSort(const unsigned int numCachedReads)
: mBuffer(NULL)
, mBufferLen(0)
{
	// set the cache size
	mSettings.NumCachedReads = numCachedReads;
}

// destructor
CPairedEndSort::~CPairedEndSort(void) {
	if(mBuffer) delete [] mBuffer;

	// delete our temporary files
	for(unsigned int i = 0; i < mTempFiles.size(); i++) rm(mTempFiles[i].c_str());
}

// retrieves an alignment from the specified temporary file and adds it to the specified list
void CPairedEndSort::AddAlignment(FILE* tempFile, const unsigned int owner, list<Alignment>& alignments) {
	Alignment al;
	if(GetAlignment(tempFile, owner, al)) alignments.push_back(al);
}

// configures which read pair types should be resolved
void CPairedEndSort::ConfigureResolution(const bool uo, const bool uu, const bool um, const bool mm) {
	mFlags.ResolveUO = uo;
	mFlags.ResolveUU = uu;
	mFlags.ResolveUM = um;
	mFlags.ResolveMM = mm;
}

// disables fragment alignment quality calculation
void CPairedEndSort::DisableFragmentAlignmentQuality(void) {
	mFlags.UseFragmentAlignmentQuality = false;
}

// allows any fragment length when evaluating unique mate-pairs
void CPairedEndSort::EnableAllUniqueFragmentLengths(void) {
	mFlags.AllowAllUniqueFragmentLengths = true;
}

// enables closest mate selection in um read pairs
void CPairedEndSort::EnableClosestMultipleMateSelection(void) {
	mFlags.FindClosestMultipleMate = true;
}

// enables consed renaming
void CPairedEndSort::EnableConsedRenaming(void) {
	mFlags.RenameMates = true;
}

// enables duplicate read filtering
void CPairedEndSort::EnableDuplicateFiltering(const string& duplicateDirectory) {
	mFlags.RemoveDuplicates      = true;
	mSettings.DuplicateDirectory = duplicateDirectory;
}

// enables the sampling of all read pairs
void CPairedEndSort::EnableFullFragmentLengthSampling(void) {
	mFlags.SampleAllFragmentLengths = true;
}

// retrieves a read from the specified temporary file
bool CPairedEndSort::GetAlignment(FILE* tempFile, const unsigned int owner, Alignment& al) {

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

	const bool isLongRead  = ((status & PE_IS_LONG_READ)           != 0 ? true : false);
	al.IsReverseStrand     = ((status & PE_IS_REVERSE_STRAND)      != 0 ? true : false);
	al.IsMateReverseStrand = ((status & PE_IS_MATE_REVERSE_STRAND) != 0 ? true : false);
	al.IsFirstMate         = ((status & PE_IS_FIRST_MATE)          != 0 ? true : false);
	al.IsResolvedAsPair    = ((status & PE_IS_RESOLVED_AS_PAIR)    != 0 ? true : false);
	al.WasRescued          = ((status & PE_WAS_RESCUED)            != 0 ? true : false);
	al.IsPairedEnd         = true;

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

	if(al.IsResolvedAsPair) {

		// get the mate reference sequence start position
		memcpy((char*)&al.MateReferenceBegin, pBuffer, SIZEOF_INT);
		pBuffer += SIZEOF_INT;

		// get the mate reference sequence end position
		memcpy((char*)&al.MateReferenceEnd, pBuffer, SIZEOF_INT);
		pBuffer += SIZEOF_INT;

		// get the mate reference sequence index
		memcpy((char*)&al.MateReferenceIndex, pBuffer, SIZEOF_INT);
		pBuffer += SIZEOF_INT;

	} /*else al.MateReferenceIndex = ALIGNMENT_NO_MATE_INFO*/;

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
void CPairedEndSort::RecordReferenceGaps(Alignment& al) {

	// localize some data
	const char* pReference       = al.Reference.CData();
	const unsigned int refLength = al.Reference.Length();

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

// resolves the paired-end reads found in the specified input file
void CPairedEndSort::ResolvePairedEndReads(const string& inputFilename, const string& outputFilename) {

	// open the alignment reader
	MosaikReadFormat::CAlignmentReader reader;
	reader.Open(inputFilename);

	// retrieve the read groups
	vector<MosaikReadFormat::ReadGroup> readGroups;
	reader.GetReadGroups(readGroups);

	// retrieve the alignment status (clear the sorting status)
	AlignmentStatus as = (reader.GetStatus() & 0xf3) | AS_SORTED_ALIGNMENT;

	// retrieve the reference sequence vector
	vector<ReferenceSequence>* pReferenceSequences = reader.GetReferenceSequences();
	mRefGapVector.resize(pReferenceSequences->size());

	// retrieve the number of reads in the archive
	const uint64_t numReads = reader.GetNumReads();

	// statistical counters
	uint64_t numBothUnique                        = 0;
	uint64_t numBothNonUnique                     = 0;
	uint64_t numOneNonUnique                      = 0;

	uint64_t numOrphaned                          = 0;
	uint64_t numUniqueOrphansResolved             = 0;
	uint64_t numUniqueOrphansDuplicateRemoved     = 0;

	uint64_t numBothUniqueResolved                = 0;
	uint64_t numBothNonUniqueResolved             = 0;
	uint64_t numOneNonUniqueResolved              = 0;

	uint64_t numBothUniqueDuplicateRemoved        = 0;
	uint64_t numBothNonUniqueDuplicateRemoved     = 0;
	uint64_t numOneNonUniqueDuplicateRemoved      = 0;

	// sanity check: we should only have one read group
	if(readGroups.size() != 1) {
		printf("ERROR: Expected one read group in the alignment archive. Found %u read groups.\n", (unsigned int)readGroups.size());
		exit(1);
	}

	// figure out if we are using paired-end or mate-pair sequencing
	// TODO: this should be automatically detected, but for now this is technology-based
	//bool usingMatePair  = false;
	//switch(readGroups[0].SequencingTechnology) {
	//case ST_454:
	//	usingMatePair = true;
	//	break;
	//case ST_ILLUMINA:
	//	break;
	//case ST_SOLID:
	//	usingMatePair = true;
	//	break;
	//default:
	//	printf("ERROR: Cannot figure out if we expect paired-end or mate-pair resolution with this sequencing technology. Assuming paired-end resolution.\n");
	//	break;
	//}

	//if(usingMatePair) printf("- resolving mate-pair alignments\n");
	//else printf("- resolving paired-end alignments\n");

	// ================================================================================
	// gather all fragments that belonging to this read group from the fragment library
	// ================================================================================

	set<string> prunedReadNames;
	set<string> prunedOrphanReadNames;
	set<string>::const_iterator prnIter;

	if(mFlags.RemoveDuplicates) {

		// derive the database filename
		string databasePath = mSettings.DuplicateDirectory + readGroups[0].LibraryName + ".db";
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

		sprintf_s(sqlBuffer, SQL_BUFFER_SIZE, "SELECT BestReadName FROM PairedFragments WHERE BestReadGroupID=%s;", readGroupCode.c_str());

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
		sqlite3_free_table(sqlResults);

		printf("%u read names loaded.\n", (unsigned int)prunedReadNames.size());

		// ===================================
		// retrieve the best orphan read names
		// ===================================

		printf("- retrieving orphan read names from the duplicate database... ");
		fflush(stdout);

		sprintf_s(sqlBuffer, SQL_BUFFER_SIZE, "SELECT BestReadName FROM OrphanFragments WHERE BestReadGroupID=%s;", readGroupCode.c_str());

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
		for(int i = 1; i <= numRows; i++) prunedOrphanReadNames.insert(sqlResults[i]);

		// clean up
		sqlite3_free_table(sqlResults);
		sqlite3_close(db);

		printf("%u orphan read names loaded.\n", (unsigned int)prunedOrphanReadNames.size());
	}

	// ===========================
	// gather fragment length data
	// ===========================

	// initialize
	unsigned int minFragmentLength    = 0;
	unsigned int medianFragmentLength = 0;
	unsigned int maxFragmentLength    = 0;

	Mosaik::AlignedRead ar;

	{
		ModelType models[8];
		for(unsigned char j = 0; j < 8; j++) models[j].ID = j + 1;

		unsigned int numFragmentLengthsCollected = 0;
		unsigned int numFragmentLengthsDesired   = 1000000;

		if(mFlags.SampleAllFragmentLengths || (numReads < numFragmentLengthsDesired)) numFragmentLengthsDesired = (unsigned int)numReads;

		vector<vector<unsigned int> > fragmentVectorPerModel;
		fragmentVectorPerModel.resize(8);
		for(unsigned char j = 0; j < 8; j++) fragmentVectorPerModel[j].reserve(numFragmentLengthsDesired);

		CConsole::Heading();
		printf("\n- phase 1 of 3: building fragment length distribution:\n");
		CConsole::Reset();

		bool gatheringFragmentLengths = true;
		CProgressCounter<unsigned int>::StartThread(&numFragmentLengthsCollected, &gatheringFragmentLengths, "samples");

		for(; numFragmentLengthsCollected < numFragmentLengthsDesired; numFragmentLengthsCollected++) {

			// get the next read
			if(!reader.LoadNextRead(ar)) break;

			// figure out which mates are unique
			const bool isMate1Unique = (ar.Mate1Alignments.size() == 1 ? true : false);
			const bool isMate2Unique = (ar.Mate2Alignments.size() == 1 ? true : false);

			// process if both are unique
			if(isMate1Unique && isMate2Unique) {

				// retrieve the start coordinates, read orientation, and reference index
				const unsigned int mate1ReferenceIndex = ar.Mate1Alignments[0].ReferenceIndex;
				const unsigned int mate2ReferenceIndex = ar.Mate2Alignments[0].ReferenceIndex;

				const unsigned int mate1RefBegin = ar.Mate1Alignments[0].ReferenceBegin;
				const unsigned int mate2RefBegin = ar.Mate2Alignments[0].ReferenceBegin;

				const unsigned int mate1RefEnd = ar.Mate1Alignments[0].ReferenceEnd;
				const unsigned int mate2RefEnd = ar.Mate2Alignments[0].ReferenceEnd;

				const bool mate1IsReverseStrand = ar.Mate1Alignments[0].IsReverseStrand;
				const bool mate2IsReverseStrand = ar.Mate2Alignments[0].IsReverseStrand;

				// skip collecting statistics if the reference indexes are not identical
				if(mate1ReferenceIndex != mate2ReferenceIndex) continue;

				// collect the model count statistics
				const unsigned int fragmentLength = (mate1RefBegin < mate2RefBegin ? mate2RefEnd - mate1RefBegin + 1 : mate1RefEnd - mate2RefBegin + 1);
				const unsigned char currentModel  = GetCurrentModel(mate1RefBegin, mate1IsReverseStrand, mate2RefBegin, mate2IsReverseStrand);

				if(fragmentLength < 10000) fragmentVectorPerModel[currentModel].push_back(fragmentLength);
				models[currentModel].Count++;
			}
		}

		// wait for the thread to end
		gatheringFragmentLengths = false;
		CProgressCounter<unsigned int>::WaitThread();
		printf("\n");

		// sanity check: make sure there is a large difference between the active
		// model counts and the unused model counts.
		sort(models, models + 8);

		const unsigned int activeModelCountSum = models[0].Count + models[1].Count;
		const unsigned int unusedModelCountSum = models[2].Count + models[3].Count + models[4].Count + models[5].Count + models[6].Count + models[7].Count;
		const double unusedPercentage = (double)unusedModelCountSum / (double)activeModelCountSum;

		mSettings.AlignmentModel1 = models[0].ID - 1;
		mSettings.AlignmentModel2 = models[1].ID - 1;

		// check the unused percentage
		if(unusedPercentage > MODEL_COUNT_THRESHOLD) {

			printf("ERROR: When determining whether to apply mate-pair or paired-end constraints, an irregularity in the alignment model counts was discovered.\n\n");

			printf("       Normal mate-pair data sets have the highest counts for alignment models:  4 & 5.\n");
			printf("       Normal paired-end data sets have the highest counts for alignment models: 2 & 6.\n\n");

			printf("       We expect that the ratio of the 6 lowest counts to the 2 highest counts to be no larger than %.2f, but in this data set the ratio was %.2f\n\n", MODEL_COUNT_THRESHOLD, unusedPercentage);

			for(unsigned char i = 0; i < 8; i++) printf("- alignment model %u: %9u hits\n", models[i].ID, models[i].Count);
			exit(1);
		}

		// emit a warning if the best alignment models are non-standard
		const bool isModel2Top = (models[0].ID == 2) || (models[1].ID == 2);
		const bool isModel4Top = (models[0].ID == 4) || (models[1].ID == 4);
		const bool isModel5Top = (models[0].ID == 5) || (models[1].ID == 5);
		const bool isModel6Top = (models[0].ID == 6) || (models[1].ID == 6);

		bool isMatePair  = (isModel4Top && isModel5Top ? true : false);
		bool isPairedEnd = (isModel2Top && isModel6Top ? true : false);

		if(isMatePair) {
			printf("- resolving mate-pair alignments\n");
		} else if(isPairedEnd) {
			printf("- resolving paired-end alignments\n");
		} else {
			printf("- WARNING: Found a non-standard alignment model configuration. Using alignment models %u & %u\n", models[0].ID, models[1].ID);
		}

		// add the fragments from the best alignment models
		vector<unsigned int> fragmentVector;
		fragmentVector.insert(fragmentVector.end(), fragmentVectorPerModel[mSettings.AlignmentModel1].begin(), fragmentVectorPerModel[mSettings.AlignmentModel1].end());
		fragmentVector.insert(fragmentVector.end(), fragmentVectorPerModel[mSettings.AlignmentModel2].begin(), fragmentVectorPerModel[mSettings.AlignmentModel2].end());

		// =========================================
		// calculate the min and max fragment length
		// =========================================

		if(fragmentVector.empty()) {
			printf("ERROR: Cannot calculate the min and max fragment length because no fragment length samples were collected.\n");
			exit(1);
		}

		// identify the min and max points on our confidence interval
		const double halfNonConfidenceInterval = (1.0 - mSettings.ConfidenceInterval) / 2.0;
		const unsigned int fragmentVectorLength = (unsigned int)fragmentVector.size();
		sort(fragmentVector.begin(), fragmentVector.end());

		minFragmentLength    = fragmentVector[(unsigned int)(fragmentVectorLength * halfNonConfidenceInterval)];
		medianFragmentLength = fragmentVector[(unsigned int)(fragmentVectorLength * 0.5)];
		maxFragmentLength    = fragmentVector[(unsigned int)(fragmentVectorLength * (1.0 - halfNonConfidenceInterval))];
	}

	// ==================
	// resolve read pairs
	// ==================

	CConsole::Heading();
	printf("\n- phase 2 of 3: resolve read pairs:\n");
	CConsole::Reset();

	// initialize
	set<string>::const_iterator dsIter;
	vector<Alignment>::const_iterator mate1Iter, mate2Iter;
	uint64_t currentRead             = 0;
	uint64_t numSerializedAlignments = 0;

	list<Alignment> alignmentCache;
	alignmentCache.resize(mSettings.NumCachedReads);
	list<Alignment>::iterator acIter = alignmentCache.begin();
	unsigned int numCachedEntries = 0;

	Alignment *pMate1Al = NULL, *pMate2Al = NULL, *pMate1UMAl = NULL, *pMate2UMAl = NULL;
	int shortestDistanceFromMedian = 0;
	bool foundClosestMate = false;

	CProgressBar<uint64_t>::StartThread(&currentRead, 0, numReads, "reads");

	// rewind the alignment reader
	reader.Rewind();

	while(reader.LoadNextRead(ar)) {

		// localize the read name
		const string readName = ar.Name.CData();

		// figure out which mates are unique
		const bool isMate1Unique = (ar.Mate1Alignments.size() == 1 ? true : false);
		const bool isMate2Unique = (ar.Mate2Alignments.size() == 1 ? true : false);

		// determine the paired-end class
		const bool isUO = ar.Mate1Alignments.empty() || ar.Mate2Alignments.empty();
		const bool isUU =  isMate1Unique &&  isMate2Unique;
		const bool isMM = !isMate1Unique && !isMate2Unique;
		const bool isUM = !isUO && !isUU && !isMM;

		// handle orphans
		if(isUO) {

			// resolve unique orphans
			if(mFlags.ResolveUO && (isMate1Unique || isMate2Unique)) {

				bool skipOrphan = false;
				if(mFlags.RemoveDuplicates) {
					prnIter = prunedOrphanReadNames.find(readName);
					if(prnIter == prunedOrphanReadNames.end()) {
						skipOrphan = true;
						numUniqueOrphansDuplicateRemoved++;
					}
				}

				if(!skipOrphan) {

					if(isMate1Unique) pMate1Al = &ar.Mate1Alignments[0];
					else pMate1Al = &ar.Mate2Alignments[0];

					// copy the mate1 alignment
					*acIter = *pMate1Al;
					acIter->Name               = ar.Name;
					acIter->ReadGroupCode      = ar.ReadGroupCode;
					if(mFlags.UseFragmentAlignmentQuality) acIter->Quality = GetFragmentAlignmentQuality(acIter->Quality, acIter->Query.Length(), isUO, isUU, isMM);
					//acIter->MateReferenceIndex = ALIGNMENT_NO_MATE_INFO;
					numCachedEntries++;
					acIter++;

					// serialize the alignment cache
					if(acIter == alignmentCache.end()) {
						numSerializedAlignments += Serialize(alignmentCache, numCachedEntries);
						acIter = alignmentCache.begin();
						numCachedEntries = 0;
					}

					numUniqueOrphansResolved++;
				}
			} 

			numOrphaned++;
			currentRead++;
			continue;
		}

		// increment the mate-pair type counters
		bool skipReadPair = false;

		if(isUU) {
			numBothUnique++;
			if(!mFlags.ResolveUU) skipReadPair = true;
		} else if(isMM) {
			numBothNonUnique++;
			if(!mFlags.ResolveMM) skipReadPair = true;
		} else {
			numOneNonUnique++;
			if(!mFlags.ResolveUM) skipReadPair = true;
		}

		// skip the read if necessary
		if(skipReadPair) {
			currentRead++;		
			continue;
		}

		pMate1Al   = NULL;
		pMate2Al   = NULL;
		pMate1UMAl = NULL;
		pMate2UMAl = NULL;

		shortestDistanceFromMedian = DEFAULT_SHORTEST_MEDIAN_DISTANCE;
		foundClosestMate = false;

		unsigned int numMatches = 0;
		for(unsigned int m1Index = 0; m1Index < ar.Mate1Alignments.size(); m1Index++) {

			// localize the mate 1 variables
			const unsigned int mate1ReferenceIndex = ar.Mate1Alignments[m1Index].ReferenceIndex;
			const unsigned int mate1RefBegin       = ar.Mate1Alignments[m1Index].ReferenceBegin;
			const unsigned int mate1RefEnd         = ar.Mate1Alignments[m1Index].ReferenceEnd;
			const bool mate1IsReverseStrand        = ar.Mate1Alignments[m1Index].IsReverseStrand;			

			for(unsigned int m2Index = 0; m2Index < ar.Mate2Alignments.size(); m2Index++) {

				// localize the mate 2 variables
				const unsigned int mate2ReferenceIndex = ar.Mate2Alignments[m2Index].ReferenceIndex;
				const unsigned int mate2RefBegin       = ar.Mate2Alignments[m2Index].ReferenceBegin;
				const unsigned int mate2RefEnd         = ar.Mate2Alignments[m2Index].ReferenceEnd;
				const bool mate2IsReverseStrand        = ar.Mate2Alignments[m2Index].IsReverseStrand;

				// check that the mates appear on the same reference
				const bool isOnSameRef = (mate1ReferenceIndex == mate2ReferenceIndex);

				// check for the proper orientation and order
				bool isProperlyOrdered = false;
				const unsigned char currentModel = GetCurrentModel(mate1RefBegin, mate1IsReverseStrand, mate2RefBegin, mate2IsReverseStrand);
				if((currentModel == mSettings.AlignmentModel1) || (currentModel == mSettings.AlignmentModel2)) isProperlyOrdered = true;

				// check for the proper fragment length
				const unsigned int fragmentLength = (mate1RefBegin < mate2RefBegin ? mate2RefEnd - mate1RefBegin + 1 : mate1RefEnd - mate2RefBegin + 1);

				bool isFragmentLengthOK = false;
				if((fragmentLength >= minFragmentLength) && (fragmentLength <= maxFragmentLength)) isFragmentLengthOK = true;
				if(mFlags.AllowAllUniqueFragmentLengths && isUU) isFragmentLengthOK = true;

				// check if this is the closest multiple mate
				if(mFlags.FindClosestMultipleMate && isUM && isProperlyOrdered && isOnSameRef) {
					const int distanceFromMedian = abs((int)fragmentLength - (int)medianFragmentLength); 

					if(distanceFromMedian < shortestDistanceFromMedian) {
						shortestDistanceFromMedian = distanceFromMedian;
						pMate1UMAl = &ar.Mate1Alignments[m1Index];
						pMate2UMAl = &ar.Mate2Alignments[m2Index];
						ExchangeMateInfo(pMate1UMAl, pMate2UMAl);
						foundClosestMate = true;
					}
				}

				// set the pointer to the properly resolved pair
				if(isProperlyOrdered && isFragmentLengthOK && isOnSameRef) {
					pMate1Al = &ar.Mate1Alignments[m1Index];
					pMate2Al = &ar.Mate2Alignments[m2Index];
					ExchangeMateInfo(pMate1Al, pMate2Al);
					numMatches++;
				}
			}
		}

		// handle duplicates
		if(mFlags.RemoveDuplicates && (numMatches == 1)) {

			prnIter = prunedReadNames.find(readName);
			if(prnIter == prunedReadNames.end()) {
				numMatches = 0;

				// increment the mate-pair type resolved counters
				if(isUU)      numBothUniqueDuplicateRemoved++;
				else if(isMM) numBothNonUniqueDuplicateRemoved++;
				else numOneNonUniqueDuplicateRemoved++;
			}
		}

		// save the read
		if((numMatches == 1) || foundClosestMate) {

			// increment the mate-pair type resolved counters
			if(isUU)      numBothUniqueResolved++;
			else if(isMM) numBothNonUniqueResolved++;
			else          numOneNonUniqueResolved++;

			// change the pointers if we a closer mate to choose from
			if(foundClosestMate && (numMatches != 1)) {
				pMate1Al = pMate1UMAl;
				pMate2Al = pMate2UMAl;
			}

			// store the summed alignment quality
			// (the SE alignment qualities are capped at 99 and thus no overflow for the PE alignment qualities)
			//const unsigned char fragmentAlignmentQuality = pMate1Al->Quality + pMate2Al->Quality;

			// copy the mate1 alignment
			*acIter = *pMate1Al;
			acIter->Name = ar.Name;
			if(mFlags.RenameMates) acIter->Name.Append("/1");

			if(mFlags.UseFragmentAlignmentQuality) 
				acIter->Quality = GetFragmentAlignmentQuality(acIter->Quality, acIter->Query.Length(), isUO, isUU, isMM);

			acIter->ReadGroupCode = ar.ReadGroupCode;
			numCachedEntries++;
			acIter++;

			// serialize the alignment cache
			if(acIter == alignmentCache.end()) {
				numSerializedAlignments += Serialize(alignmentCache, numCachedEntries);
				acIter = alignmentCache.begin();
				numCachedEntries = 0;
			}

			// copy the mate2 alignment
			*acIter = *pMate2Al;
			acIter->Name = ar.Name;
			if(mFlags.RenameMates) acIter->Name.Append("/2");

			if(mFlags.UseFragmentAlignmentQuality) 
				acIter->Quality = GetFragmentAlignmentQuality(acIter->Quality, acIter->Query.Length(), isUO, isUU, isMM);

			acIter->ReadGroupCode = ar.ReadGroupCode;
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

	// wait for the progress bar to finish
	CProgressBar<uint64_t>::WaitThread();

	// close our files
	reader.Close();

	// ====================================================
	// serialize the remaining reads in the alignment cache
	// ====================================================

	if(acIter != alignmentCache.begin()) numSerializedAlignments += Serialize(alignmentCache, numCachedEntries);

	// =======================
	// sort the resolved reads
	// =======================

	CConsole::Heading();
	printf("\n- phase 3 of 3: sort resolved read pairs:\n");
	CConsole::Reset();

	// open our output file
	MosaikReadFormat::CAlignmentWriter aw;
	aw.Open(outputFilename, *pReferenceSequences, readGroups, as);

	// allocate the file stream array
	const unsigned int numTempFiles = (unsigned int)mTempFiles.size();
	FILE** tempFile = new FILE*[numTempFiles];

	// sanity check: check the number of temp files
	if(numTempFiles > 65535) {
		printf("ERROR: More than 65535 temporary files were produced during paired-end resolution. The owner variable needs to be redimensioned.\n");
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

	const uint64_t numTotalReads                  = numBothUnique + numOneNonUnique + numBothNonUnique + numOrphaned;
	const uint64_t numTotalDuplicatesRemovedReads = numBothUniqueDuplicateRemoved + numOneNonUniqueDuplicateRemoved + numBothNonUniqueDuplicateRemoved + numUniqueOrphansDuplicateRemoved;
	const uint64_t numTotalResolvedReads          = numBothUniqueResolved + numOneNonUniqueResolved + numBothNonUniqueResolved + numUniqueOrphansResolved;

	CConsole::Heading(); printf("\nPaired-end read statistics:\n"); CConsole::Reset();

	if(mFlags.RemoveDuplicates) {

		printf("=====================================================================================\n");
		printf("                                original      removed duplicates          resolved   \n");
		printf("-------------------------------------------------------------------------------------\n");

		printf("# orphaned:              %9llu (%4.1f %%)   %9llu (%4.1f %%)   %9llu (%4.1f %%)\n", (unsigned long long)numOrphaned,      (numOrphaned      / (double)numTotalReads) * 100.0, (unsigned long long)numUniqueOrphansDuplicateRemoved, (numUniqueOrphansDuplicateRemoved / (double)numTotalReads) * 100.0, (unsigned long long)numUniqueOrphansResolved, (numUniqueOrphansResolved / (double)numTotalReads) * 100.0);
		printf("# both mates unique:     %9llu (%4.1f %%)   %9llu (%4.1f %%)   %9llu (%4.1f %%)\n", (unsigned long long)numBothUnique,    (numBothUnique    / (double)numTotalReads) * 100.0, (unsigned long long)numBothUniqueDuplicateRemoved,    (numBothUniqueDuplicateRemoved    / (double)numTotalReads) * 100.0, (unsigned long long)numBothUniqueResolved,    (numBothUniqueResolved    / (double)numTotalReads) * 100.0);
		printf("# one mate non-unique:   %9llu (%4.1f %%)   %9llu (%4.1f %%)   %9llu (%4.1f %%)\n", (unsigned long long)numOneNonUnique,  (numOneNonUnique  / (double)numTotalReads) * 100.0, (unsigned long long)numOneNonUniqueDuplicateRemoved,  (numOneNonUniqueDuplicateRemoved  / (double)numTotalReads) * 100.0, (unsigned long long)numOneNonUniqueResolved,  (numOneNonUniqueResolved  / (double)numTotalReads) * 100.0);
		printf("# both mates non-unique: %9llu (%4.1f %%)   %9llu (%4.1f %%)   %9llu (%4.1f %%)\n", (unsigned long long)numBothNonUnique, (numBothNonUnique / (double)numTotalReads) * 100.0, (unsigned long long)numBothNonUniqueDuplicateRemoved, (numBothNonUniqueDuplicateRemoved / (double)numTotalReads) * 100.0, (unsigned long long)numBothNonUniqueResolved, (numBothNonUniqueResolved / (double)numTotalReads) * 100.0);

		printf("-------------------------------------------------------------------------------------\n");
		printf("total:                   %9llu            %9llu (%4.1f %%)   %9llu (%4.1f %%)\n\n", (unsigned long long)numTotalReads, (unsigned long long)numTotalDuplicatesRemovedReads, (numTotalDuplicatesRemovedReads / (double)numTotalReads) * 100.0, (unsigned long long)numTotalResolvedReads, (numTotalResolvedReads / (double)numTotalReads) * 100.0);

	} else {

		printf("================================================================\n");
		printf("                                original             resolved   \n");
		printf("----------------------------------------------------------------\n");

		printf("# orphaned:              %9llu (%4.1f %%)   %9llu (%4.1f %%)\n", (unsigned long long)numOrphaned,      (numOrphaned      / (double)numTotalReads) * 100.0, (unsigned long long)numUniqueOrphansResolved, (numUniqueOrphansResolved / (double)numTotalReads) * 100.0);
		printf("# both mates unique:     %9llu (%4.1f %%)   %9llu (%4.1f %%)\n", (unsigned long long)numBothUnique,    (numBothUnique    / (double)numTotalReads) * 100.0, (unsigned long long)numBothUniqueResolved,    (numBothUniqueResolved    / (double)numTotalReads) * 100.0);
		printf("# one mate non-unique:   %9llu (%4.1f %%)   %9llu (%4.1f %%)\n", (unsigned long long)numOneNonUnique,  (numOneNonUnique  / (double)numTotalReads) * 100.0, (unsigned long long)numOneNonUniqueResolved,  (numOneNonUniqueResolved  / (double)numTotalReads) * 100.0);
		printf("# both mates non-unique: %9llu (%4.1f %%)   %9llu (%4.1f %%)\n", (unsigned long long)numBothNonUnique, (numBothNonUnique / (double)numTotalReads) * 100.0, (unsigned long long)numBothNonUniqueResolved, (numBothNonUniqueResolved / (double)numTotalReads) * 100.0);

		printf("----------------------------------------------------------------\n");
		printf("total:                   %9llu            %9llu (%4.1f %%)\n\n", (unsigned long long)numTotalReads, (unsigned long long)numTotalResolvedReads, (numTotalResolvedReads / (double)numTotalReads) * 100.0);
	}

	CConsole::Heading(); printf("Fragment statistics:\n"); CConsole::Reset();
	printf("================================\n");
	printf("min target frag len:     %7u\n", minFragmentLength);
	printf("median target frag len:  %7u\n", medianFragmentLength);
	printf("max target frag len:     %7u\n\n", maxFragmentLength);
}

// serializes the specified list to a temporary file
uint64_t CPairedEndSort::Serialize(list<Alignment>& alignmentCache, const unsigned int numEntries) {

	uint64_t numSerializedAlignments = 0;

	// return if we have nothing to serialize
	if(numEntries == 0) return numSerializedAlignments;

	// remove reads if the alignment cache is not full
	if(numEntries != mSettings.NumCachedReads) alignmentCache.resize(numEntries);

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
		const unsigned int requestedBytes = (SIZEOF_INT * 4) + (SIZEOF_SHORT * 2) + 6 + readNameLen + pairwiseLength + pairwiseLength + bqLength + (isLongRead ? 3 : 0) + (alIter->IsResolvedAsPair ? SIZEOF_INT * 3 : 0);

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
		if(isLongRead)                  status |= PE_IS_LONG_READ;
		if(alIter->IsReverseStrand)     status |= PE_IS_REVERSE_STRAND;
		if(alIter->IsMateReverseStrand) status |= PE_IS_MATE_REVERSE_STRAND;
		if(alIter->IsFirstMate)         status |= PE_IS_FIRST_MATE;
		if(alIter->IsResolvedAsPair)    status |= PE_IS_RESOLVED_AS_PAIR;
		if(alIter->WasRescued)          status |= PE_WAS_RESCUED;

		mBuffer[bufferOffset++] = status;

		// store the number of mismatches
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

		if(alIter->IsResolvedAsPair) {

			// store the mate reference sequence start position
			memcpy(mBuffer + bufferOffset, (char*)&alIter->MateReferenceBegin, SIZEOF_INT);
			bufferOffset += SIZEOF_INT;

			// store the mate reference sequence end position
			memcpy(mBuffer + bufferOffset, (char*)&alIter->MateReferenceEnd, SIZEOF_INT);
			bufferOffset += SIZEOF_INT;

			// store the mate reference sequence index
			memcpy(mBuffer + bufferOffset, (char*)&alIter->MateReferenceIndex, SIZEOF_INT);
			bufferOffset += SIZEOF_INT;
		}

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

// sets the desired confidence interval
void CPairedEndSort::SetConfidenceInterval(const double& percent) {
	mSettings.ConfidenceInterval = percent;
}
