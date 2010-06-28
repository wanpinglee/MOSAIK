// ***************************************************************************
// CMosaikDupSnoop - records duplicate fragments in sequencing libraries.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "MosaikDupSnoop.h"

// constructor
CMosaikDupSnoop::CMosaikDupSnoop(void) {}

// destructor
CMosaikDupSnoop::~CMosaikDupSnoop(void) {}

// calculates the aggregate quality
unsigned int CMosaikDupSnoop::CalculateAggregateQuality(vector<Alignment>::iterator& alIter) {

	unsigned int aggregateQuality = 0;
	const char* pBQ = alIter->BaseQualities.CData();
	for(unsigned int i = alIter->QueryBegin; i <= alIter->QueryEnd; i++) {
		aggregateQuality += pBQ[i];
	}

	// add the alignment quality
	aggregateQuality |= (alIter->Quality << 24);

	return aggregateQuality;
}

// configures which read pair types should be resolved
void CMosaikDupSnoop::ConfigureResolution(const bool uo, const bool uu, const bool um, const bool mm) {
	mFlags.ResolveUO = uo;
	mFlags.ResolveUU = uu;
	mFlags.ResolveUM = um;
	mFlags.ResolveMM = mm;
}

// consolidates the PairedFragments table taking into account variable endpoints
void CMosaikDupSnoop::ConsolidatePairedEndFragments(const string& outputDirectory, const set<string>& libraryNames) {

	// initialize
	ostringstream sb;
	char* errorMessage = NULL;
	char** sqlResults  = NULL;
	sqlite3* db        = NULL;
	int numRows        = 0;
	int numColumns     = 0;

	set<string>::const_iterator lnIter;
	for(lnIter = libraryNames.begin(); lnIter != libraryNames.end(); lnIter++) {

		CConsole::Heading(); printf("Consolidating fragments from library: %s\n", lnIter->c_str()); CConsole::Reset();

		// create the database name
		sb << outputDirectory << *lnIter << ".db";
		const string libraryFilename = sb.str();
		sb.str("");

		// open the database
		if(sqlite3_open(libraryFilename.c_str(), &db)) {
			printf("ERROR: Unable to open the database (%s) for writing.\n", libraryFilename.c_str());
			printf("       error message: %s\n", sqlite3_errmsg(db));
			exit(1);
		}

		// set exclusive locking mode
		if(sqlite3_exec(db, "PRAGMA locking_mode = EXCLUSIVE;", NULL, 0, &errorMessage) != SQLITE_OK) {
			printf("ERROR: Unable to set exclusive locking mode on the database (%s).\n", libraryFilename.c_str());
			printf("       error message: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		// disable the synchronous mode
		if(sqlite3_exec(db, "PRAGMA synchronous = OFF;", NULL, 0, &errorMessage) != SQLITE_OK) {
			printf("ERROR: Unable to disable synchronous mode on the database (%s).\n", libraryFilename.c_str());
			printf("       error message: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		// use memory journaling
		if(sqlite3_exec(db, "PRAGMA journal_mode = MEMORY;", NULL, 0, &errorMessage) != SQLITE_OK) {
			printf("ERROR: Unable to use memory journaling on the database (%s).\n", libraryFilename.c_str());
			printf("       error message: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		// =====================================
		// check if we have any paired fragments
		// =====================================

		sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "SELECT COUNT(*) FROM PairedFragments;");

		if(sqlite3_get_table(db, mSqlBuffer, &sqlResults, &numRows, &numColumns, &errorMessage) != SQLITE_OK) {
			printf("ERROR: The SQL query resulted in the following error: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		if((numRows != 1) || (numColumns != 1)) {
			printf("ERROR: Expected one column when receiving the paired fragment count.\n");
			printf("rows: %u, columns: %u\n", numRows, numColumns);
			exit(1);
		}

		const unsigned int numPairedFragmentRows = GetUnsignedInt(sqlResults[1]);

		sqlite3_free_table(sqlResults);

		if(numPairedFragmentRows == 0) {
			cout << "- no paired-end fragments found. skipping library." << endl << endl;
			sqlite3_close(db);
			continue;
		}

		// ================================================================================
		// retrieve the fragments with a fixed begin coordinate and variable end coordinate
		// ================================================================================

		sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "SELECT a.ROWID, a.Count, a.BestReadQuality, a.BestReadName, a.BestReadGroupID, b.ROWID, b.Count, b.BestReadQuality, b.BestReadName, b.BestReadGroupID FROM PairedFragments AS a, PairedFragments AS b WHERE (a.ReferenceIndex = b.ReferenceIndex) AND (a.Begin = b.Begin) AND (a.End != b.End) AND (a.End BETWEEN b.End-2 AND b.End+2);");

		cout << "- performing first consolidation SQL query... ";
		cout.flush();

		CBenchmark bench1;
		bench1.Start();

		if(sqlite3_get_table(db, mSqlBuffer, &sqlResults, &numRows, &numColumns, &errorMessage) != SQLITE_OK) {
			printf("ERROR: The SQL query resulted in the following error: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		bench1.Stop();
		double wallTime = bench1.GetElapsedWallTime();
		cout << wallTime << " s." << endl;

		// ===================================================================================
		// consolidate the fragments with a fixed begin coordinate and variable end coordinate
		// ===================================================================================

		if(numRows > 0) {

			// sanity check: ensure that we have an answer with 10 columns
			if(numColumns != 10) {
				printf("ERROR: Expected ten columns when receiving the consolidation results.\n");
				printf("rows: %u, columns: %u\n", numRows, numColumns);
				exit(1);
			}

			char* bestReadName    = NULL;
			char* bestReadGroupID = NULL;
			char* bestReadQuality = NULL;

			char* aROWID          = NULL;
			char* bROWID          = NULL;

			cout << "- consolidating results... ";
			cout.flush();

			CBenchmark bench2;
			bench2.Start();

			for(int i = 1; i <= numRows; i++) {

				const unsigned int currentIndex = numColumns * i;

				const unsigned int aCount = GetUnsignedInt(sqlResults[currentIndex + 1]);
				const unsigned int bCount = GetUnsignedInt(sqlResults[currentIndex + 6]);
				const unsigned int newCount = aCount + bCount;

				const unsigned int aBestReadQuality = GetUnsignedInt(sqlResults[currentIndex + 2]);
				const unsigned int bBestReadQuality = GetUnsignedInt(sqlResults[currentIndex + 7]);

				aROWID = sqlResults[currentIndex];
				bROWID = sqlResults[currentIndex + 5];

				// figure out which is the best read
				if(aBestReadQuality > bBestReadQuality) {
					bestReadQuality = sqlResults[currentIndex + 2];
					bestReadName    = sqlResults[currentIndex + 3];
					bestReadGroupID = sqlResults[currentIndex + 4];
				} else {
					bestReadQuality = sqlResults[currentIndex + 7];
					bestReadName    = sqlResults[currentIndex + 8];
					bestReadGroupID = sqlResults[currentIndex + 9];
				}

				if(aCount > bCount) {

					// ===============================
					// merge the results from b into a
					// ===============================

					sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "UPDATE PairedFragments SET Count=%u, BestReadName='%s', BestReadQuality=%s, BestReadGroupID=%s WHERE ROWID=%s;", 
						newCount, bestReadName, bestReadQuality, bestReadGroupID, aROWID);

					if(sqlite3_exec(db, mSqlBuffer, NULL, 0, &errorMessage) != SQLITE_OK) {
						printf("ERROR: Unable to update the fragment in the database (%s).\n", libraryFilename.c_str());
						printf("       error message: %s\n", errorMessage);
						sqlite3_free(errorMessage);
						exit(1);
					}

					// delete b
					sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "DELETE FROM PairedFragments WHERE ROWID=%s;", bROWID);

					if(sqlite3_exec(db, mSqlBuffer, NULL, 0, &errorMessage) != SQLITE_OK) {
						printf("ERROR: Unable to update the fragment in the database (%s).\n", libraryFilename.c_str());
						printf("       error message: %s\n", errorMessage);
						sqlite3_free(errorMessage);
						exit(1);
					}

				} else {

					// ===============================
					// merge the results from a into b
					// ===============================

					sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "UPDATE PairedFragments SET Count=%u, BestReadName='%s', BestReadQuality=%s, BestReadGroupID=%s WHERE ROWID=%s;", 
						newCount, bestReadName, bestReadQuality, bestReadGroupID, bROWID);

					if(sqlite3_exec(db, mSqlBuffer, NULL, 0, &errorMessage) != SQLITE_OK) {
						printf("ERROR: Unable to update the fragment in the database (%s).\n", libraryFilename.c_str());
						printf("       error message: %s\n", errorMessage);
						sqlite3_free(errorMessage);
						exit(1);
					}

					// delete b
					sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "DELETE FROM PairedFragments WHERE ROWID=%s;", aROWID);

					if(sqlite3_exec(db, mSqlBuffer, NULL, 0, &errorMessage) != SQLITE_OK) {
						printf("ERROR: Unable to update the fragment in the database (%s).\n", libraryFilename.c_str());
						printf("       error message: %s\n", errorMessage);
						sqlite3_free(errorMessage);
						exit(1);
					}
				}
			}

			bench2.Stop();
			wallTime = bench2.GetElapsedWallTime();
			cout << wallTime << " s." << endl;
		}

		sqlite3_free_table(sqlResults);

		// ================================================================================
		// retrieve the fragments with a variable begin coordinate and fixed end coordinate
		// ================================================================================

		sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "SELECT a.ROWID, a.Count, a.BestReadQuality, a.BestReadName, a.BestReadGroupID, b.ROWID, b.Count, b.BestReadQuality, b.BestReadName, b.BestReadGroupID FROM PairedFragments AS a, PairedFragments AS b WHERE (a.ReferenceIndex = b.ReferenceIndex) AND (a.End = b.End) AND (a.Begin != b.Begin) AND (a.Begin BETWEEN b.Begin-2 AND b.Begin+2);");

		cout << "- performing second consolidation SQL query... ";
		cout.flush();

		CBenchmark bench3;
		bench3.Start();

		if(sqlite3_get_table(db, mSqlBuffer, &sqlResults, &numRows, &numColumns, &errorMessage) != SQLITE_OK) {
			printf("ERROR: The SQL query resulted in the following error: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		bench3.Stop();
		wallTime = bench3.GetElapsedWallTime();
		cout << wallTime << " s." << endl;

		// ===================================================================================
		// consolidate the fragments with a variable begin coordinate and fixed end coordinate
		// ===================================================================================

		if(numRows > 0) {

			// sanity check: ensure that we have an answer with 10 columns
			if(numColumns != 10) {
				printf("ERROR: Expected ten columns when receiving the consolidation results.\n");
				printf("rows: %u, columns: %u\n", numRows, numColumns);
				exit(1);
			}

			char* bestReadName    = NULL;
			char* bestReadGroupID = NULL;
			char* bestReadQuality = NULL;

			char* aROWID          = NULL;
			char* bROWID          = NULL;

			cout << "- consolidating results... ";
			cout.flush();

			CBenchmark bench4;
			bench4.Start();

			for(int i = 1; i <= numRows; i++) {

				const unsigned int currentIndex = numColumns * i;

				const unsigned int aCount = GetUnsignedInt(sqlResults[currentIndex + 1]);
				const unsigned int bCount = GetUnsignedInt(sqlResults[currentIndex + 6]);
				const unsigned int newCount = aCount + bCount;

				const unsigned int aBestReadQuality = GetUnsignedInt(sqlResults[currentIndex + 2]);
				const unsigned int bBestReadQuality = GetUnsignedInt(sqlResults[currentIndex + 7]);

				aROWID = sqlResults[currentIndex];
				bROWID = sqlResults[currentIndex + 5];

				// figure out which is the best read
				if(aBestReadQuality > bBestReadQuality) {
					bestReadQuality = sqlResults[currentIndex + 2];
					bestReadName    = sqlResults[currentIndex + 3];
					bestReadGroupID = sqlResults[currentIndex + 4];
				} else {
					bestReadQuality = sqlResults[currentIndex + 7];
					bestReadName    = sqlResults[currentIndex + 8];
					bestReadGroupID = sqlResults[currentIndex + 9];
				}

				if(aCount > bCount) {

					// ===============================
					// merge the results from b into a
					// ===============================

					sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "UPDATE PairedFragments SET Count=%u, BestReadName='%s', BestReadQuality=%s, BestReadGroupID=%s WHERE ROWID=%s;", 
						newCount, bestReadName, bestReadQuality, bestReadGroupID, aROWID);

					if(sqlite3_exec(db, mSqlBuffer, NULL, 0, &errorMessage) != SQLITE_OK) {
						printf("ERROR: Unable to update the fragment in the database (%s).\n", libraryFilename.c_str());
						printf("       error message: %s\n", errorMessage);
						sqlite3_free(errorMessage);
						exit(1);
					}

					// delete b
					sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "DELETE FROM PairedFragments WHERE ROWID=%s;", bROWID);

					if(sqlite3_exec(db, mSqlBuffer, NULL, 0, &errorMessage) != SQLITE_OK) {
						printf("ERROR: Unable to update the fragment in the database (%s).\n", libraryFilename.c_str());
						printf("       error message: %s\n", errorMessage);
						sqlite3_free(errorMessage);
						exit(1);
					}

				} else {

					// ===============================
					// merge the results from a into b
					// ===============================

					sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "UPDATE PairedFragments SET Count=%u, BestReadName='%s', BestReadQuality=%s, BestReadGroupID=%s WHERE ROWID=%s;", 
						newCount, bestReadName, bestReadQuality, bestReadGroupID, bROWID);

					if(sqlite3_exec(db, mSqlBuffer, NULL, 0, &errorMessage) != SQLITE_OK) {
						printf("ERROR: Unable to update the fragment in the database (%s).\n", libraryFilename.c_str());
						printf("       error message: %s\n", errorMessage);
						sqlite3_free(errorMessage);
						exit(1);
					}

					// delete b
					sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "DELETE FROM PairedFragments WHERE ROWID=%s;", aROWID);

					if(sqlite3_exec(db, mSqlBuffer, NULL, 0, &errorMessage) != SQLITE_OK) {
						printf("ERROR: Unable to update the fragment in the database (%s).\n", libraryFilename.c_str());
						printf("       error message: %s\n", errorMessage);
						sqlite3_free(errorMessage);
						exit(1);
					}
				}
			}

			bench4.Stop();
			wallTime = bench4.GetElapsedWallTime();
			cout << wallTime << " s." << endl << endl;
		}

		// clean up
		sqlite3_free_table(sqlResults);
		sqlite3_close(db);
	}
}

// creates databases for each library if they don't already exist
void CMosaikDupSnoop::CreateLibraryDatabases(const string& outputDirectory, const set<string>& libraryNames) {

	// create the tables for each database
	char* errorMessage = NULL;
	ostringstream sb;

	set<string>::const_iterator lnIter;
	for(lnIter = libraryNames.begin(); lnIter != libraryNames.end(); ++lnIter) {

		// create the database name
		sb << outputDirectory << *lnIter << ".db";
		const string libraryName = sb.str();
		sb.str("");

		// open the database
		sqlite3* db = NULL;

		if(sqlite3_open(libraryName.c_str(), &db)) {
			printf("ERROR: Unable to open the database (%s) for writing.\n", libraryName.c_str());
			printf("       error message: %s\n", sqlite3_errmsg(db));
			exit(1);
		}

		// set the page size
		if(sqlite3_exec(db, "PRAGMA page_size = 4096;", NULL, 0, &errorMessage) != SQLITE_OK) {
			printf("ERROR: Unable to set the page size on the database (%s).\n", libraryName.c_str());
			printf("       error message: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		// set the default cache size
		if(sqlite3_exec(db, "PRAGMA default_cache_size = 262144;", NULL, 0, &errorMessage) != SQLITE_OK) {
			printf("ERROR: Unable to set the default cache size on the database (%s).\n", libraryName.c_str());
			printf("       error message: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		// create the tables
		/*
		BEGIN;

		CREATE TABLE IF NOT EXISTS PairedFragments (
		ReferenceIndex INTEGER,
		Begin INTEGER,
		End INTEGER,
		Count INTEGER,
		BestReadName TEXT,
		BestReadQuality INTEGER,
		BestReadGroupID INTEGER
		);

		CREATE TABLE IF NOT EXISTS OrphanFragments (
		ReferenceIndex INTEGER,
		Begin INTEGER,
		End INTEGER,
		Count INTEGER,
		BestReadName TEXT,
		BestReadQuality INTEGER,
		BestReadGroupID INTEGER
		);

		CREATE TABLE IF NOT EXISTS SingleFragments (
		ReferenceIndex INTEGER,
		Begin INTEGER,
		End INTEGER,
		Count INTEGER,
		BestReadName TEXT,
		BestReadQuality INTEGER,
		BestReadGroupID INTEGER
		);

		CREATE TABLE IF NOT EXISTS ReadGroups (
		ID INTEGER PRIMARY KEY,
		Name TEXT
		);

		CREATE INDEX IF NOT EXISTS BestReadGroupIDPairedFragmentsIndex ON PairedFragments(BestReadGroupID);
		CREATE INDEX IF NOT EXISTS BeginPairedFragmentsIndex ON PairedFragments(Begin, End, ReferenceIndex);
		CREATE INDEX IF NOT EXISTS EndPairedFragmentsIndex ON PairedFragments(End, Begin, ReferenceIndex);

		CREATE INDEX IF NOT EXISTS BestReadGroupIDSingleFragmentsIndex ON SingleFragments(BestReadGroupID);
		CREATE INDEX IF NOT EXISTS BeginSingleFragmentsIndex ON SingleFragments(Begin, End, ReferenceIndex);
		CREATE INDEX IF NOT EXISTS EndSingleFragmentsIndex ON SingleFragments(End, Begin, ReferenceIndex);

		CREATE INDEX IF NOT EXISTS BestReadGroupIDOrphanFragmentsIndex ON OrphanFragments(BestReadGroupID);
		CREATE INDEX IF NOT EXISTS BeginOrphanFragmentsIndex ON OrphanFragments(Begin, End, ReferenceIndex);
		CREATE INDEX IF NOT EXISTS EndOrphanFragmentsIndex ON OrphanFragments(End, Begin, ReferenceIndex);
		*/

		const string createTablesSQL = "BEGIN; CREATE TABLE IF NOT EXISTS PairedFragments (ReferenceIndex INTEGER, Begin INTEGER, End INTEGER, Count INTEGER, BestReadName TEXT, BestReadQuality INTEGER, BestReadGroupID INTEGER); CREATE TABLE IF NOT EXISTS OrphanFragments (ReferenceIndex INTEGER, Begin INTEGER, End INTEGER, Count INTEGER, BestReadName TEXT, BestReadQuality INTEGER, BestReadGroupID INTEGER); CREATE TABLE IF NOT EXISTS SingleFragments (ReferenceIndex INTEGER, Begin INTEGER, End INTEGER, Count INTEGER, BestReadName TEXT, BestReadQuality INTEGER, BestReadGroupID INTEGER); CREATE TABLE IF NOT EXISTS ReadGroups (ID INTEGER PRIMARY KEY, Name TEXT); CREATE INDEX IF NOT EXISTS BestReadGroupIDPairedFragmentsIndex ON PairedFragments(BestReadGroupID); CREATE INDEX IF NOT EXISTS BeginPairedFragmentsIndex ON PairedFragments(Begin, End, ReferenceIndex); CREATE INDEX IF NOT EXISTS EndPairedFragmentsIndex ON PairedFragments(End, Begin, ReferenceIndex); CREATE INDEX IF NOT EXISTS BestReadGroupIDSingleFragmentsIndex ON SingleFragments(BestReadGroupID); CREATE INDEX IF NOT EXISTS BeginSingleFragmentsIndex ON SingleFragments(Begin, End, ReferenceIndex); CREATE INDEX IF NOT EXISTS EndSingleFragmentsIndex ON SingleFragments(End, Begin, ReferenceIndex); CREATE INDEX IF NOT EXISTS BestReadGroupIDOrphanFragmentsIndex ON OrphanFragments(BestReadGroupID); CREATE INDEX IF NOT EXISTS BeginOrphanFragmentsIndex ON OrphanFragments(Begin, End, ReferenceIndex); CREATE INDEX IF NOT EXISTS EndOrphanFragmentsIndex ON OrphanFragments(End, Begin, ReferenceIndex);";

		if(sqlite3_exec(db, createTablesSQL.c_str(), NULL, 0, &errorMessage) != SQLITE_OK) {
			printf("ERROR: Unable to create tables in the database (%s).\n", libraryName.c_str());
			printf("       error message: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		// commit the database
		if(sqlite3_exec(db, "COMMIT;", NULL, 0, &errorMessage) != SQLITE_OK) {
			printf("ERROR: Unable to commit the database (%s).\n", libraryName.c_str());
			printf("       error message: %s\n", errorMessage);
			sqlite3_free(errorMessage);
		}

		// close the database
		sqlite3_close(db);
	}
}

// allows any fragment length when evaluating unique mate-pairs
void CMosaikDupSnoop::EnableAllUniqueFragmentLengths(void) {
	mFlags.AllowAllUniqueFragmentLengths = true;
}

// populates the library databases with alignment data
void CMosaikDupSnoop::PopulateLibraryDatabases(const vector<string>& inputFiles, const string& outputDirectory) {

	// initialize
	ostringstream sb;
	char* errorMessage = NULL;
	char** sqlResults  = NULL;
	sqlite3* db        = NULL;

	vector<string>::const_iterator sIter;
	for(sIter = inputFiles.begin(); sIter != inputFiles.end(); ++sIter) {

		// open the file and store the library name
		MosaikReadFormat::CAlignmentReader reader;
		reader.Open(*sIter);
		AlignmentStatus as = reader.GetStatus();

		vector<MosaikReadFormat::ReadGroup> readGroups;
		vector<MosaikReadFormat::ReadGroup>::const_iterator rgIter;
		reader.GetReadGroups(readGroups);

		// retrieve the number of reads in the archive
		const uint64_t numReads = reader.GetNumReads();

		// sanity check: we only support one read group per unsorted alignment file for now
		if(readGroups.size() != 1) {
			printf("ERROR: Currently only one read group per unsorted alignment archive is supported.\n");
			exit(1);
		}

		// create the database name
		sb << outputDirectory << readGroups.begin()->LibraryName << ".db";
		const string libraryFilename = sb.str();
		sb.str("");

		// open the database
		if(sqlite3_open(libraryFilename.c_str(), &db)) {
			printf("ERROR: Unable to open the database (%s) for writing.\n", libraryFilename.c_str());
			printf("       error message: %s\n", sqlite3_errmsg(db));
			exit(1);
		}

		// set exclusive locking mode
		if(sqlite3_exec(db, "PRAGMA locking_mode = EXCLUSIVE;", NULL, 0, &errorMessage) != SQLITE_OK) {
			printf("ERROR: Unable to set exclusive locking mode on the database (%s).\n", libraryFilename.c_str());
			printf("       error message: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		// disable the synchronous mode
		if(sqlite3_exec(db, "PRAGMA synchronous = OFF;", NULL, 0, &errorMessage) != SQLITE_OK) {
			printf("ERROR: Unable to disable synchronous mode on the database (%s).\n", libraryFilename.c_str());
			printf("       error message: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		//// use memory journaling
		//if(sqlite3_exec(db, "PRAGMA journal_mode = MEMORY;", NULL, 0, &errorMessage) != SQLITE_OK) {
		//	printf("ERROR: Unable to use memory journaling on the database (%s).\n", libraryFilename.c_str());
		//	printf("       error message: %s\n", errorMessage);
		//	sqlite3_free(errorMessage);
		//	exit(1);
		//}

		//if(sqlite3_exec(db, "PRAGMA temp_store = MEMORY;", NULL, 0, &errorMessage) != SQLITE_OK) {
		//	printf("ERROR: Unable to use memory temporary storage on the database (%s).\n", libraryFilename.c_str());
		//	printf("       error message: %s\n", errorMessage);
		//	sqlite3_free(errorMessage);
		//	exit(1);
		//}

		// add the read group to the database
		snprintf(mSqlBuffer, SQL_BUFFER_SIZE, "BEGIN; INSERT INTO ReadGroups (Name) VALUES ('%s');", readGroups.begin()->ReadGroupID.c_str());

		if(sqlite3_exec(db, mSqlBuffer, NULL, 0, &errorMessage) != SQLITE_OK) {
			printf("ERROR: Unable to add the read group the database (%s).\n", libraryFilename.c_str());
			printf("       error message: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		const unsigned int readGroupID = (unsigned int)sqlite3_last_insert_rowid(db);

		CConsole::Heading(); printf("Parsing %s:\n", sIter->c_str()); CConsole::Reset();

		const bool isPairedEndArchive = (((as & AS_PAIRED_END_READ) != 0) ? true : false);

		if(isPairedEndArchive) {

			// ===========================
			// gather fragment length data
			// ===========================

			// initialize
			unsigned int numFragmentLengthsCollected = 0;
			unsigned int numFragmentLengthsDesired   = 1000000;

			if(numReads < numFragmentLengthsDesired) numFragmentLengthsDesired = (unsigned int)numReads;

			vector<unsigned int> fragmentVector;
			fragmentVector.reserve(numFragmentLengthsDesired);

			// set the appropriate sequencing technology flags
			bool isUsing454 = false;
			switch(readGroups.begin()->SequencingTechnology) {
				case ST_454:
					isUsing454 = true;
					break;
				case ST_ILLUMINA:
					break;
				default:
					printf("ERROR: The sequencing technology specified in the alignment format is currently not supported.\n");
					break;
			}

			printf("- phase 1 of 2: building fragment length distribution:\n");

			bool gatheringFragmentLengths = true;
			CProgressCounter<unsigned int>::StartThread(&numFragmentLengthsCollected, &gatheringFragmentLengths, "samples");

			Mosaik::AlignedRead ar;
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

					// enforce the ordering criteria
					bool isProperlyOrdered = false;

					if(isUsing454) { // 454
						if(mate1RefBegin < mate2RefBegin) {
							if(mate1IsReverseStrand && mate2IsReverseStrand)      isProperlyOrdered = true;
						} else if(!mate2IsReverseStrand && !mate1IsReverseStrand) isProperlyOrdered = true;
					} else {         // Illumina
						if(mate1RefBegin < mate2RefBegin) {
							if(!mate1IsReverseStrand && mate2IsReverseStrand)    isProperlyOrdered = true;
						} else if(!mate2IsReverseStrand && mate1IsReverseStrand) isProperlyOrdered = true;
					}

					// store the fragment length
					if(isProperlyOrdered) {
						if(mate1ReferenceIndex == mate2ReferenceIndex) {
							const unsigned int fragmentLength = (mate1RefBegin < mate2RefBegin ? mate2RefEnd - mate1RefBegin + 1 : mate1RefEnd - mate2RefBegin + 1);
							if(fragmentLength < 10000) fragmentVector.push_back(fragmentLength);
						}
					}
				}
			}

			// wait for the thread to end
			gatheringFragmentLengths = false;
			CProgressCounter<unsigned int>::WaitThread();

			printf("\n");

			// =========================================
			// calculate the min and max fragment length
			// =========================================

			if(fragmentVector.empty()) {
				printf("ERROR: Cannot calculate the min and max fragment length because no fragment length samples were collected.\n");
				exit(1);
			}

			// identify the min and max points on our confidence interval
			const double halfNonConfidenceInterval = (1.0 - DEFAULT_CONFIDENCE_INTERVAL) / 2.0;
			const unsigned int fragmentVectorLength = fragmentVector.size();
			sort(fragmentVector.begin(), fragmentVector.end());

			unsigned int minFragmentLength = fragmentVector[(unsigned int)(fragmentVectorLength * halfNonConfidenceInterval)];
			unsigned int maxFragmentLength = fragmentVector[(unsigned int)(fragmentVectorLength * (1.0 - halfNonConfidenceInterval))];

			fragmentVector.clear();

			// =======================
			// record fragment lengths
			// =======================

			printf("- phase 2 of 2: record fragment lengths: (%u - %u bp)\n", minFragmentLength, maxFragmentLength);

			// initialize
			set<string>::const_iterator dsIter;
			vector<Alignment>::iterator alIter, mate1Iter, mate2Iter, resolvedMate1Iter, resolvedMate2Iter;
			uint64_t currentRead = 0;
			int numSqlRows = 0, numSqlColumns = 0;
			unsigned int numAlignments = 0;

			CProgressBar<uint64_t>::StartThread(&currentRead, 0, numReads, "reads");

			// rewind the alignment reader
			reader.Rewind();

			while(reader.LoadNextRead(ar)) {

				// figure out which mates are unique
				const bool isMate1Unique = (ar.Mate1Alignments.size() == 1 ? true : false);
				const bool isMate2Unique = (ar.Mate2Alignments.size() == 1 ? true : false);

				// handle orphans
				if(ar.Mate1Alignments.empty() || ar.Mate2Alignments.empty()) {

					// resolve unique orphans
					if(mFlags.ResolveUO && (isMate1Unique || isMate2Unique)) {

						// assign our iterator
						if(ar.Mate1Alignments.empty()) {
							alIter        = ar.Mate2Alignments.begin();
							numAlignments = ar.Mate2Alignments.size();
						} else {
							alIter        = ar.Mate1Alignments.begin();
							numAlignments = ar.Mate1Alignments.size();
						}

						// skip the orphan if non-unique
						if(numAlignments != 1) {
							currentRead++;
							continue;
						}

						// calculate the aggregate quality
						const unsigned int aggregateQuality = CalculateAggregateQuality(alIter);

						// best read info from the database
						sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "SELECT ROWID, Count, BestReadQuality FROM OrphanFragments WHERE ReferenceIndex=%u AND Begin=%u AND End=%u;", alIter->ReferenceIndex, alIter->ReferenceBegin, alIter->ReferenceEnd);

						if(sqlite3_get_table(db, mSqlBuffer, &sqlResults, &numSqlRows, &numSqlColumns, &errorMessage) != SQLITE_OK) {
							printf("ERROR: Unable to grab the fragment data from the database (%s).\n", libraryFilename.c_str());
							printf("       error message: %s\n", errorMessage);
							sqlite3_free(errorMessage);
							exit(1);
						}

						// add the fragment to the database
						if((numSqlRows == 0) && (numSqlColumns == 0)) {

							sqlite3_free_table(sqlResults);

							sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "INSERT INTO OrphanFragments VALUES (%u, %u, %u, 1, '%s', %u, %u);", alIter->ReferenceIndex, alIter->ReferenceBegin, alIter->ReferenceEnd, ar.Name.CData(), aggregateQuality, readGroupID);

							if(sqlite3_exec(db, mSqlBuffer, NULL, 0, &errorMessage) != SQLITE_OK) {
								printf("ERROR: Unable to add the fragment to the database (%s).\n", libraryFilename.c_str());
								printf("       error message: %s\n", errorMessage);
								sqlite3_free(errorMessage);
								exit(1);
							}

						} else { // update the entry if needed

							// sanity checking
							if((numSqlColumns != 3) || (numSqlRows != 1)) {
								printf("ERROR: Expected 1 row and 3 columns to be returned from the database.\n");
								printf("       rows: %u, columns: %u\n", numSqlRows, numSqlColumns);
								exit(1);
							}

							const char* rowID                       = sqlResults[3];
							const unsigned int newCount             = GetUnsignedInt(sqlResults[4]) + 1;
							const unsigned int bestAggregateQuality = GetUnsignedInt(sqlResults[5]);

							if(aggregateQuality > bestAggregateQuality) {

								sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "UPDATE OrphanFragments SET Count=%u, BestReadName='%s', BestReadQuality=%u, BestReadGroupID=%u WHERE ROWID=%s;", 
									newCount, ar.Name.CData(), aggregateQuality, readGroupID, rowID);

								if(sqlite3_exec(db, mSqlBuffer, NULL, 0, &errorMessage) != SQLITE_OK) {
									printf("ERROR: Unable to update the fragment in the database (%s).\n", libraryFilename.c_str());
									printf("       error message: %s\n", errorMessage);
									printf("       SQL query: [%s]\n", mSqlBuffer);
									sqlite3_free(errorMessage);
									exit(1);
								}
							}

							sqlite3_free_table(sqlResults);
						}
					}

					currentRead++;
					continue;
				}

				// increment the mate-pair type counters
				bool skipReadPair = false;

				if(isMate1Unique && isMate2Unique) {
					if(!mFlags.ResolveUU) skipReadPair = true;
				} else if(!isMate1Unique && !isMate2Unique) {
					if(!mFlags.ResolveMM) skipReadPair = true;
				} else {
					if(!mFlags.ResolveUM) skipReadPair = true;
				}

				// skip the read if necessary
				if(skipReadPair) {
					currentRead++;		
					continue;
				}

				unsigned int numMatches = 0;
				for(mate1Iter = ar.Mate1Alignments.begin(); mate1Iter != ar.Mate1Alignments.end(); ++mate1Iter) {

					// localize the mate 1 variables
					const unsigned int mate1ReferenceIndex = mate1Iter->ReferenceIndex;
					const unsigned int mate1RefBegin       = mate1Iter->ReferenceBegin;
					const unsigned int mate1RefEnd         = mate1Iter->ReferenceEnd;
					const bool mate1IsReverseStrand    = mate1Iter->IsReverseStrand;			

					for(mate2Iter = ar.Mate2Alignments.begin(); mate2Iter != ar.Mate2Alignments.end(); ++mate2Iter) {

						// localize the mate 2 variables
						const unsigned int mate2ReferenceIndex = mate2Iter->ReferenceIndex;
						const unsigned int mate2RefBegin       = mate2Iter->ReferenceBegin;
						const unsigned int mate2RefEnd         = mate2Iter->ReferenceEnd;
						const bool mate2IsReverseStrand    = mate2Iter->IsReverseStrand;	

						bool isProperlyOrdered = false;

						if(isUsing454) { // 454
							if(mate1RefBegin < mate2RefBegin) {
								if(mate1IsReverseStrand && mate2IsReverseStrand)      isProperlyOrdered = true;
							} else if(!mate2IsReverseStrand && !mate1IsReverseStrand) isProperlyOrdered = true;
						} else {         // Illumina
							if(mate1RefBegin < mate2RefBegin) {
								if(!mate1IsReverseStrand && mate2IsReverseStrand)     isProperlyOrdered = true;
							} else if(!mate2IsReverseStrand && mate1IsReverseStrand)  isProperlyOrdered = true;
						}

						// decide if this combination fits the constraints
						const unsigned int fragmentLength = (mate1RefBegin < mate2RefBegin ? mate2RefEnd - mate1RefBegin + 1 : mate1RefEnd - mate2RefBegin + 1);

						bool isFragmentLengthOK = false;
						if((fragmentLength >= minFragmentLength) && (fragmentLength <= maxFragmentLength)) isFragmentLengthOK = true;
						if(mFlags.AllowAllUniqueFragmentLengths && isMate1Unique && isMate2Unique)         isFragmentLengthOK = true;

						if(isProperlyOrdered && isFragmentLengthOK && (mate1ReferenceIndex == mate2ReferenceIndex)) {
							numMatches++;
							resolvedMate1Iter = mate1Iter;
							resolvedMate2Iter = mate2Iter;

							resolvedMate1Iter->MateReferenceIndex      = resolvedMate2Iter->ReferenceIndex;
							resolvedMate1Iter->MateReferenceBegin      = resolvedMate2Iter->ReferenceBegin;
							resolvedMate1Iter->MateReferenceEnd        = resolvedMate2Iter->ReferenceEnd;
							resolvedMate1Iter->IsMateReverseStrand     = resolvedMate2Iter->IsReverseStrand;

							resolvedMate2Iter->MateReferenceIndex      = resolvedMate1Iter->ReferenceIndex;
							resolvedMate2Iter->MateReferenceBegin      = resolvedMate1Iter->ReferenceBegin;
							resolvedMate2Iter->MateReferenceEnd        = resolvedMate1Iter->ReferenceEnd;
							resolvedMate2Iter->IsMateReverseStrand     = resolvedMate1Iter->IsReverseStrand;
						}
					}
				}

				// store the fragment length information
				if(numMatches == 1) {

					// calculate the fragment endpoints
					unsigned int fragmentBegin = 0, fragmentEnd = 0;
					if(resolvedMate1Iter->ReferenceBegin < resolvedMate2Iter->ReferenceBegin) {
						fragmentBegin = resolvedMate1Iter->ReferenceBegin;
						fragmentEnd   = resolvedMate2Iter->ReferenceEnd;
					} else {
						fragmentBegin = resolvedMate2Iter->ReferenceBegin;
						fragmentEnd   = resolvedMate1Iter->ReferenceEnd;
					}

					// calculate the aggregate quality
					const unsigned int aggregateQuality = CalculateAggregateQuality(resolvedMate1Iter) + CalculateAggregateQuality(resolvedMate2Iter);

					// best read info from the database
					sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "SELECT ROWID, Count, BestReadQuality FROM PairedFragments WHERE ReferenceIndex=%u AND Begin=%u AND End=%u;", resolvedMate1Iter->ReferenceIndex, fragmentBegin, fragmentEnd);

					if(sqlite3_get_table(db, mSqlBuffer, &sqlResults, &numSqlRows, &numSqlColumns, &errorMessage) != SQLITE_OK) {
						printf("ERROR: Unable to grab the fragment data from the database (%s).\n", libraryFilename.c_str());
						printf("       error message: %s\n", errorMessage);
						sqlite3_free(errorMessage);
						exit(1);
					}

					// add the fragment to the database
					if((numSqlRows == 0) && (numSqlColumns == 0)) {

						sqlite3_free_table(sqlResults);

						sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "INSERT INTO PairedFragments VALUES (%u, %u, %u, 1, '%s', %u, %u);", resolvedMate1Iter->ReferenceIndex, fragmentBegin, fragmentEnd, ar.Name.CData(), aggregateQuality, readGroupID);

						if(sqlite3_exec(db, mSqlBuffer, NULL, 0, &errorMessage) != SQLITE_OK) {
							printf("ERROR: Unable to add the fragment to the database (%s).\n", libraryFilename.c_str());
							printf("       error message: %s\n", errorMessage);
							sqlite3_free(errorMessage);
							exit(1);
						}

					} else { // update the entry if needed

						// sanity checking
						if((numSqlColumns != 3) || (numSqlRows != 1)) {
							printf("ERROR: Expected 1 row and 3 columns to be returned from the database.\n");
							printf("       rows: %u, columns: %u\n", numSqlRows, numSqlColumns);
							exit(1);
						}

						const char* rowID                       = sqlResults[3];
						const unsigned int newCount             = GetUnsignedInt(sqlResults[4]) + 1;
						const unsigned int bestAggregateQuality = GetUnsignedInt(sqlResults[5]);

						if(aggregateQuality > bestAggregateQuality) {

							sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "UPDATE PairedFragments SET Count=%u, BestReadName='%s', BestReadQuality=%u, BestReadGroupID=%u WHERE ROWID=%s;", 
								newCount, ar.Name.CData(), aggregateQuality, readGroupID, rowID);

							if(sqlite3_exec(db, mSqlBuffer, NULL, 0, &errorMessage) != SQLITE_OK) {
								printf("ERROR: Unable to update the fragment in the database (%s).\n", libraryFilename.c_str());
								printf("       error message: %s\n", errorMessage);
								sqlite3_free(errorMessage);
								exit(1);
							}
						}

						sqlite3_free_table(sqlResults);
					}
				}

				// increment the read counter
				currentRead++;
			}

			// wait for the progress bar to finish
			CProgressBar<uint64_t>::WaitThread();

			printf("\n");

		} else {

			// ======================
			// gather single-end data
			// ======================

			printf("- recording unique read lengths:\n");

			int numSqlRows = 0, numSqlColumns = 0;
			vector<Alignment>::iterator alIter;
			uint64_t currentRead = 0;

			CProgressBar<uint64_t>::StartThread(&currentRead, 0, numReads, "reads");

			Mosaik::AlignedRead ar;
			while(reader.LoadNextRead(ar)) {

				// process only the unique reads
				const bool isMate1Unique = (ar.Mate1Alignments.size() == 1 ? true : false);
				if(!isMate1Unique) {
					currentRead++;				
					continue;
				}

				// calculate the aggregate quality
				alIter = ar.Mate1Alignments.begin();
				const unsigned int aggregateQuality = CalculateAggregateQuality(alIter);

				// best read info from the database
				sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "SELECT ROWID, Count, BestReadQuality FROM SingleFragments WHERE ReferenceIndex=%u AND Begin=%u AND End=%u;", alIter->ReferenceIndex, alIter->ReferenceBegin, alIter->ReferenceEnd);

				if(sqlite3_get_table(db, mSqlBuffer, &sqlResults, &numSqlRows, &numSqlColumns, &errorMessage) != SQLITE_OK) {
					printf("ERROR: Unable to grab the fragment data from the database (%s).\n", libraryFilename.c_str());
					printf("       error message: %s\n", errorMessage);
					sqlite3_free(errorMessage);
					exit(1);
				}

				// add the fragment to the database
				if((numSqlRows == 0) && (numSqlColumns == 0)) {

					sqlite3_free_table(sqlResults);

					sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "INSERT INTO SingleFragments VALUES (%u, %u, %u, 1, '%s', %u, %u);", alIter->ReferenceIndex, alIter->ReferenceBegin, alIter->ReferenceEnd, ar.Name.CData(), aggregateQuality, readGroupID);

					if(sqlite3_exec(db, mSqlBuffer, NULL, 0, &errorMessage) != SQLITE_OK) {
						printf("ERROR: Unable to add the fragment to the database (%s).\n", libraryFilename.c_str());
						printf("       error message: %s\n", errorMessage);
						sqlite3_free(errorMessage);
						exit(1);
					}

				} else { // update the entry if needed

					// sanity checking
					if((numSqlColumns != 3) || (numSqlRows != 1)) {
						printf("ERROR: Expected 1 row and 3 columns to be returned from the database.\n");
						printf("       rows: %u, columns: %u\n", numSqlRows, numSqlColumns);
						exit(1);
					}

					const char* rowID                       = sqlResults[3];
					const unsigned int newCount             = GetUnsignedInt(sqlResults[4]) + 1;
					const unsigned int bestAggregateQuality = GetUnsignedInt(sqlResults[5]);

					if(aggregateQuality > bestAggregateQuality) {

						sprintf_s(mSqlBuffer, SQL_BUFFER_SIZE, "UPDATE SingleFragments SET Count=%u, BestReadName='%s', BestReadQuality=%u, BestReadGroupID=%u WHERE ROWID=%s;", 
							newCount, ar.Name.CData(), aggregateQuality, readGroupID, rowID);

						if(sqlite3_exec(db, mSqlBuffer, NULL, 0, &errorMessage) != SQLITE_OK) {
							printf("ERROR: Unable to update the fragment in the database (%s).\n", libraryFilename.c_str());
							printf("       error message: %s\n", errorMessage);
							sqlite3_free(errorMessage);
							exit(1);
						}
					}

					sqlite3_free_table(sqlResults);
				}

				// increment the read counter
				currentRead++;
			}

			// wait for the progress bar to finish
			CProgressBar<uint64_t>::WaitThread();

			printf("\n");
		}

		// close our files
		reader.Close();

		if(sqlite3_exec(db, "COMMIT;", NULL, 0, &errorMessage) != SQLITE_OK) {
			printf("ERROR: Unable to commit transaction in the database (%s).\n", libraryFilename.c_str());
			printf("       error message: %s\n", errorMessage);
			sqlite3_free(errorMessage);
			exit(1);
		}

		// close the database
		sqlite3_close(db);
	}
}

// records the fragments found in the input files
void CMosaikDupSnoop::RecordFragments(const vector<string>& inputFiles, const string& outputDirectory) {

	// =====================================================
	// extract the library names from the alignment archives
	// =====================================================

	set<string> libraryNames;
	set<string>::const_iterator lnIter;

	RetrieveLibraryNames(inputFiles, libraryNames);

	// display the libraries that were found
	CConsole::Heading(); printf("Databases for the following libraries will be created:\n"); CConsole::Reset();

	for(lnIter = libraryNames.begin(); lnIter != libraryNames.end(); lnIter++) {
		printf("- %s\n", lnIter->c_str());
	}

	printf("\n");

	// ====================
	// create the databases
	// ====================

	printf("Creating databases... ");
	fflush(stdout);

	CreateLibraryDatabases(outputDirectory, libraryNames);

	printf("finished.\n\n");

	// ===================================================
	// parse each alignment file and populate the database
	// ===================================================

	PopulateLibraryDatabases(inputFiles, outputDirectory);

	// ================================
	// consolidate paired-end fragments
	// ================================

	ConsolidatePairedEndFragments(outputDirectory, libraryNames);
}

// retrieves the library names used in the supplied alignment archives
void CMosaikDupSnoop::RetrieveLibraryNames(const vector<string>& inputFiles, set<string>& libraryNames) {

	// display the files that will be analyzed
	CConsole::Heading(); printf("Scanning the following alignment archives:\n"); CConsole::Reset();

	vector<string>::const_iterator sIter;
	for(sIter = inputFiles.begin(); sIter != inputFiles.end(); ++sIter) {
		printf("- %s\n", sIter->c_str());

		// open the file and store the library name
		MosaikReadFormat::CAlignmentReader reader;
		reader.Open(*sIter);

		vector<MosaikReadFormat::ReadGroup> readGroups;
		vector<MosaikReadFormat::ReadGroup>::const_iterator rgIter;
		reader.GetReadGroups(readGroups);
		for(rgIter = readGroups.begin(); rgIter != readGroups.end(); ++rgIter) {
			libraryNames.insert(rgIter->LibraryName);
		}

		reader.Close();
	}

	printf("\n");

	// sanity check: make sure we have some libraries
	if(libraryNames.empty()) {
		printf("ERROR: No libraries were identified when parsing alignment archives.\n");
		exit(1);
	}
}

// sets the desired confidence interval
void CMosaikDupSnoop::SetConfidenceInterval(const double& percent) {
	mSettings.ConfidenceInterval = percent;
}
