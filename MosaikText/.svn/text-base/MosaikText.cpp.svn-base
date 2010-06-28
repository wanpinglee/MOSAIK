// ***************************************************************************
// CMosaikText - exports alignments to various file formats.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "MosaikText.h"

// constructor
CMosaikText::CMosaikText(void) {}

// destructor
CMosaikText::~CMosaikText(void) {}

// enables AXT output
void CMosaikText::EnableAxtOutput(const string& filename) { 
	mFlags.IsAxtEnabled   = true;
	mSettings.AxtFilename = filename;
}

// enables BAM output
void CMosaikText::EnableBamOutput(const string& filename) {
	mFlags.IsBamEnabled   = true;
	mSettings.BamFilename = filename;
}

// enables BED output
void CMosaikText::EnableBedOutput(const string& filename) { 
	mFlags.IsBedEnabled   = true;
	mSettings.BedFilename = filename;
}

// enables Eland output
void CMosaikText::EnableElandOutput(const string& filename) { 
	mFlags.IsElandEnabled   = true;
	mSettings.ElandFilename = filename;
}

// enables FASTQ output
void CMosaikText::EnableFastqOutput(const string& filename) {
	mFlags.IsFastqEnabled   = true;
	mSettings.FastqFilename = filename;
}

// enables Psl output
void CMosaikText::EnablePslOutput(const string& filename) { 
	mFlags.IsPslEnabled   = true;
	mSettings.PslFilename = filename;
}

// enables the reference sequence filter
void CMosaikText::EnableReferenceFilter(const string& referenceName, const string& alignmentFilename) {
	mFlags.UseReferenceFilter = true;

	MosaikReadFormat::CAlignmentReader reader;
	reader.Open(alignmentFilename);

	// find the reference ID associated with the reference name
	bool foundReferenceIndex = false;

	unsigned int currentReferenceIndex = 0;
	uint64_t numAlignments             = 0;

	vector<ReferenceSequence>* pRefSeqs = reader.GetReferenceSequences();
	for(vector<ReferenceSequence>::const_iterator rsIter = pRefSeqs->begin(); rsIter != pRefSeqs->end(); ++rsIter, ++currentReferenceIndex) {
		if(rsIter->Name == referenceName) {
			foundReferenceIndex = true;
			numAlignments = rsIter->NumAligned;
			break;
		}		
	}

	reader.Close();

	// make sure the reference exists in the alignment archive
	if(!foundReferenceIndex) {
		printf("ERROR: The chosen reference sequence (%s) was not found in the alignment archive.\n", referenceName.c_str());
		exit(1);
	}

	// make sure we have some alignments to filter
	if(numAlignments == 0) {
		printf("ERROR: There were no alignments assigned to the chosen reference sequence (%s).\n", referenceName.c_str());
		exit(1);
	}

	mSettings.FilteredReferenceIndex    = currentReferenceIndex;
	mSettings.NumFilteredReferenceReads = numAlignments;
}

// enables SAM output
void CMosaikText::EnableSamOutput(const string& filename) {
	mFlags.IsSamEnabled   = true;
	mSettings.SamFilename = filename;

	const unsigned int filenameLength = filename.size();
	if((filenameLength >= 3) && (filename.substr(filenameLength - 3) != ".gz")) mSettings.SamFilename += ".gz";
}

// enables screen output
void CMosaikText::EnableScreenOutput(void) {
	mFlags.IsScreenEnabled = true;
}

// when triggered, the coverage calculation will only include unique reads
void CMosaikText::EvaluateUniqueReadsOnly(void) { 
	cout << "- evaluating unique reads only." << endl;
	mFlags.EvaluateUniqueReadsOnly = true;
}

// opens the output file stream for the AXT file
void CMosaikText::InitializeAxt(void) {
	mStreams.axt = fopen(mSettings.AxtFilename.c_str(), "wb");

	if(!mStreams.axt) {
		printf("ERROR: Unable to open the AXT file for writing.\n");
		exit(1);
	}
}

// opens the output file stream for the BAM file
void CMosaikText::InitializeBam(const bool isSortedByPosition, vector<ReferenceSequence>* pRefSeqs, vector<MosaikReadFormat::ReadGroup>& readGroups) {

	// set the sort order
	BamHeader header;
	header.SortOrder = (isSortedByPosition ? SORTORDER_POSITION : SORTORDER_UNSORTED);

	// set the reference sequences and read groups
	header.pReferenceSequences = pRefSeqs;
	header.pReadGroups         = &readGroups;

	// open the bam file
	mStreams.bam.Open(mSettings.BamFilename, header);
}

// opens the output file stream for the BED file
void CMosaikText::InitializeBed(void) {
	mStreams.bed = fopen(mSettings.BedFilename.c_str(), "wb");

	if(!mStreams.bed) {
		printf("ERROR: Unable to open the BED file for writing.\n");
		exit(1);
	}
}

// opens the output file stream for the Eland file
void CMosaikText::InitializeEland(void) {
	mStreams.eland = fopen(mSettings.ElandFilename.c_str(), "wb");

	if(!mStreams.eland) {
		printf("ERROR: Unable to open the Eland file for writing.\n");
		exit(1);
	}
}

// opens the output file stream for the SAM file
void CMosaikText::InitializeSam(const bool isSortedByPosition, vector<ReferenceSequence>* pRefSeqs, vector<MosaikReadFormat::ReadGroup>& readGroups) {
	mStreams.sam = gzopen(mSettings.SamFilename.c_str(), "wb");

	if(!mStreams.sam) {
		printf("ERROR: Unable to open the SAM file for writing.\n");
		exit(1);
	}

	// store the header
	gzprintf(mStreams.sam, "@HD\tVN:1.0\tSO:");
	if(isSortedByPosition) gzprintf(mStreams.sam, "coordinate\n");
	else gzprintf(mStreams.sam, "unsorted\n");

	// store the sequence dictionary
	vector<ReferenceSequence>::const_iterator rsIter;
	for(rsIter = pRefSeqs->begin(); rsIter != pRefSeqs->end(); ++rsIter) {
		gzprintf(mStreams.sam, "@SQ\tSN:%s\tLN:%u", rsIter->Name.c_str(), rsIter->NumBases);
		if(!rsIter->GenomeAssemblyID.empty()) gzprintf(mStreams.sam, "\tAS:%s", rsIter->GenomeAssemblyID.c_str());
		if(!rsIter->MD5.empty())              gzprintf(mStreams.sam, "\tM5:%s", rsIter->MD5.c_str());
		if(!rsIter->URI.empty())              gzprintf(mStreams.sam, "\tUR:%s", rsIter->URI.c_str());
		if(!rsIter->Species.empty())          gzprintf(mStreams.sam, "\tSP:%s", rsIter->Species.c_str());
		gzprintf(mStreams.sam, "\n");
	}

	// store the read groups
	vector<MosaikReadFormat::ReadGroup>::const_iterator rgIter;
	for(rgIter = readGroups.begin(); rgIter != readGroups.end(); ++rgIter) {
		gzprintf(mStreams.sam, "@RG\tID:%s\tSM:%s", rgIter->ReadGroupID.c_str(), rgIter->SampleName.c_str());
		if(!rgIter->LibraryName.empty())      gzprintf(mStreams.sam, "\tLB:%s", rgIter->LibraryName.c_str());
		if(!rgIter->Description.empty())      gzprintf(mStreams.sam, "\tDS:%s", rgIter->Description.c_str());
		if(!rgIter->PlatformUnit.empty())     gzprintf(mStreams.sam, "\tPU:%s", rgIter->PlatformUnit.c_str());
		if(rgIter->MedianFragmentLength != 0) gzprintf(mStreams.sam, "\tPI:%u", rgIter->MedianFragmentLength);
		if(!rgIter->CenterName.empty())       gzprintf(mStreams.sam, "\tCN:%s", rgIter->CenterName.c_str());

		switch(rgIter->SequencingTechnology) {
			case ST_454:
				gzprintf(mStreams.sam, "\tPL:454\n");
				break;
			case ST_HELICOS:
				gzprintf(mStreams.sam, "\tPL:helicos\n");
				break;
			case ST_ILLUMINA:
				gzprintf(mStreams.sam, "\tPL:illumina\n");
				break;
			case ST_PACIFIC_BIOSCIENCES:
				gzprintf(mStreams.sam, "\tPL:pacific biosciences\n");
				break;
			case ST_SOLID:
				gzprintf(mStreams.sam, "\tPL:solid\n");
				break;
			case ST_SANGER:
				gzprintf(mStreams.sam, "\tPL:sanger\n");
				break;
			default:
				gzprintf(mStreams.sam, "\tPL:unknown\n");
				break;
		}
	}
}

// parses the specified MOSAIK alignment file and matching anchors file
void CMosaikText::ParseMosaikAlignmentFile(const string& alignmentFilename) {

	// open the alignment archive
	MosaikReadFormat::CAlignmentReader reader;
	reader.Open(alignmentFilename);

	// retrieve the alignment archive status
	const AlignmentStatus as = reader.GetStatus();
	const bool isSortedByPosition = ((as & AS_SORTED_ALIGNMENT) != 0 ? true : false);

	// retrieve the reference sequences
	vector<ReferenceSequence>* pRefSeqs = reader.GetReferenceSequences();

	// retrieve the read groups
	vector<MosaikReadFormat::ReadGroup> readGroups;	
	reader.GetReadGroups(readGroups);

	// retrieve the total number of reads
	uint64_t numReads = reader.GetNumReads();

	// jump to the desired reference index
	if(mFlags.UseReferenceFilter && isSortedByPosition) {
		reader.Jump(mSettings.FilteredReferenceIndex, 0);
		numReads = mSettings.NumFilteredReferenceReads;
	}

	// open our output file streams
	if(mFlags.IsAxtEnabled)   InitializeAxt();
	if(mFlags.IsBamEnabled)   InitializeBam(isSortedByPosition, pRefSeqs, readGroups);
	if(mFlags.IsBedEnabled)   InitializeBed();
	if(mFlags.IsElandEnabled) InitializeEland();
	if(mFlags.IsSamEnabled)   InitializeSam(isSortedByPosition, pRefSeqs, readGroups);

	// initialize
	mCurrentRead      = 0;
	mCurrentAlignment = 0;

	if(!mFlags.IsScreenEnabled) {
		CConsole::Heading(); printf("Converting alignment archive:\n"); CConsole::Reset();
		CProgressBar<uint64_t>::StartThread(&mCurrentRead, 0, numReads, (isSortedByPosition ? "alignments" : "reads"));
	}

	// retrieve all reads from the alignment reader
	Mosaik::AlignedRead ar;
	while(reader.LoadNextRead(ar)) {

		// stop processing reads if we're already past the current reference sequence
		if(mFlags.UseReferenceFilter && isSortedByPosition && 
			(ar.Mate1Alignments.begin()->ReferenceIndex > mSettings.FilteredReferenceIndex)) break;

		const unsigned int numMate1Alignments = (unsigned int)ar.Mate1Alignments.size();
		const unsigned int numMate2Alignments = (unsigned int)ar.Mate2Alignments.size();

		// convert the read group code to a read group ID string
		MosaikReadFormat::ReadGroup rg = reader.GetReadGroupFromCode(ar.ReadGroupCode);
		const bool isColorspace = (rg.SequencingTechnology == ST_SOLID ? true : false);

		// dump the mate 1 alignments
		if(numMate1Alignments > 0) {
			if((mFlags.EvaluateUniqueReadsOnly && (numMate1Alignments == 1)) || !mFlags.EvaluateUniqueReadsOnly) {
				ProcessAlignments(1, ar.Name, isColorspace, ar.Mate1Alignments, rg.ReadGroupID);
			}
		}

		// dump the mate 2 alignments
		if(numMate2Alignments > 0) {
			if((mFlags.EvaluateUniqueReadsOnly && (numMate2Alignments == 1)) || !mFlags.EvaluateUniqueReadsOnly) {
				ProcessAlignments(2, ar.Name, isColorspace, ar.Mate2Alignments, rg.ReadGroupID);
			}
		}

		// increment the read counter
		++mCurrentRead;
	}

	// wait for the progress bar to finish
	if(!mFlags.IsScreenEnabled) CProgressBar<uint64_t>::WaitThread();

	// close our file streams
	reader.Close();
	if(mFlags.IsAxtEnabled)   fclose(mStreams.axt);
	if(mFlags.IsBamEnabled)   mStreams.bam.Close();
	if(mFlags.IsBedEnabled)   fclose(mStreams.bed);
	if(mFlags.IsElandEnabled) fclose(mStreams.eland);
	if(mFlags.IsSamEnabled)   gzclose(mStreams.sam);
}

// processes the alignments according to the chosen file format
void CMosaikText::ProcessAlignments(const unsigned char mateNum, const CMosaikString& readName, const bool isColorspace, vector<Alignment>& alignments, const string& readGroupID) {

	// initialize
	vector<Alignment>::iterator alIter;
	const bool isUnique = (alignments.size() == 1 ? true : false);

	// ===================================
	// first pass: collect additional data
	// ===================================

	if(mFlags.IsElandEnabled) {

		// initialize
		unsigned int mismatchCounts[3];
		uninitialized_fill(mismatchCounts, mismatchCounts + 3, 0);

		unsigned short lowestMismatchCount = MAX_SHORT;

		unsigned char bestAlignmentQuality = 0;
		vector<Alignment>::const_iterator bestAlIter = alignments.begin();

		// record mismatch counts and identify the best alignment
		for(alIter = alignments.begin(); alIter != alignments.end(); ++alIter) {

			unsigned short numMismatches = alIter->NumMismatches;

			if(numMismatches < lowestMismatchCount) lowestMismatchCount = numMismatches;
			if(numMismatches <= 2) ++mismatchCounts[numMismatches];

			if(alIter->Quality > bestAlignmentQuality) {
				bestAlignmentQuality = alIter->Quality;
				bestAlIter           = alIter;
			}
		}

		// write the best alignment
		fprintf(mStreams.eland, "%llu\t%s\t%s\t%c\t%u\t%u\t%u\t%u\t", (unsigned long long)(mCurrentRead + 1),
			readName.CData(), bestAlIter->Query.CData(), (isUnique ? 'U' : 'R'), lowestMismatchCount,
			mismatchCounts[0], mismatchCounts[1], mismatchCounts[2]);
	}

	// ===========================
	// second pass: writing output
	// ===========================

	for(alIter = alignments.begin(); alIter != alignments.end(); ++alIter) {

		// skip the alignment if filtering is activated
		// TODO: of course this is useless for the ELAND output
		if(mFlags.UseReferenceFilter && (alIter->ReferenceIndex != mSettings.FilteredReferenceIndex)) continue;

		if(mFlags.IsAxtEnabled) {
			fprintf(mStreams.axt, "%llu %s %u %u %s %u %u %c %u\n", (unsigned long long)(mCurrentAlignment + 1),
				alIter->ReferenceName, alIter->ReferenceBegin + 1, alIter->ReferenceEnd + 1, readName.CData(),
				alIter->QueryBegin + 1, alIter->QueryEnd + 1, (alIter->IsReverseStrand ? '-' : '+'), alIter->Quality);
			fprintf(mStreams.axt, "%s\n%s\n", alIter->Reference.CData(), alIter->Query.CData());
		}

		if(mFlags.IsBamEnabled) mStreams.bam.SaveAlignment(readName, readGroupID, alIter);

		if(mFlags.IsBedEnabled) {
			fprintf(mStreams.bed, "%s %u %u %s 1 %c %u %u %s\n", alIter->ReferenceName, alIter->ReferenceBegin + 1,
				alIter->ReferenceEnd + 1, readName.CData(), (alIter->IsReverseStrand ? '-' : '+'), 
				alIter->ReferenceBegin + 1, alIter->ReferenceEnd + 1, (alIter->IsReverseStrand ? "0,0,255" : "255,0,0"));
		}

		if(mFlags.IsElandEnabled && isUnique) {
			fprintf(mStreams.eland, "\t%s\t%u\t%c\t..", alIter->ReferenceName, alIter->ReferenceBegin + 1,
				(alIter->IsReverseStrand ? 'R' : 'F'));

			const char* pReference = alIter->Reference.CData();
			const char* pQuery     = alIter->Query.CData();

			// TODO: update this - it probably doesn't work correctly when alignment is on reverse strand
			for(unsigned short j = 0; j < (unsigned short)alIter->Reference.Length(); ++j) {
				if(pReference[j] != pQuery[j]) fprintf(mStreams.eland, "\t%u", j + alIter->QueryBegin + 1);
			}
		}

		if(mFlags.IsSamEnabled) WriteSamEntry(readName, readGroupID, alIter);

		if(mFlags.IsScreenEnabled) {
			CMosaikString bq = alIter->BaseQualities;
			bq.Increment(33);

			printf("%llu %s %u %u %s %u %u %c %u %s%s\n", (unsigned long long)(mCurrentAlignment + 1), alIter->ReferenceName, 
				alIter->ReferenceBegin + 1, alIter->ReferenceEnd + 1, readName.CData(), alIter->QueryBegin + 1,
				alIter->QueryEnd + 1, (alIter->IsReverseStrand ? '-' : '+'), alIter->Quality, readGroupID.c_str(),
				(alIter->WasRescued ? " [rescued]" : ""));
			printf("%s\n%s\n%s\n\n", alIter->Reference.CData(), alIter->Query.CData(), bq.CData());
		}

		// increment our alignment counter
		++mCurrentAlignment;
	}

	// ============
	// final output
	// ============

	if(mFlags.IsElandEnabled) fprintf(mStreams.eland, "\n");
}

// parses the specified MOSAIK read file
void CMosaikText::ParseMosaikReadFile(const string& readFilename) {

	MosaikReadFormat::CReadReader reader;
	reader.Open(readFilename);

	// retrieve the sequencing technology
	const MosaikReadFormat::ReadGroup rg = reader.GetReadGroup();
	const bool isColorspace = (rg.SequencingTechnology == ST_SOLID ? true : false);

	// retrieve the read status
	const ReadStatus rs = reader.GetStatus();
	const bool isPairedEnd = ((rs & RS_PAIRED_END_READ) != 0 ? true : false);

	// open our output file streams
	if(mFlags.IsFastqEnabled) {
		mStreams.fastq = fopen(mSettings.FastqFilename.c_str(), "wb");

		if(!mStreams.fastq) {
			printf("ERROR: Unable to open the FASTQ file for writing.\n");
			exit(1);
		}
	}

	// parse all of the reads
	if(!mFlags.IsScreenEnabled) {
		printf("- parsing reads... ");
		fflush(stdout);
	}

	// retrieve all reads from the read reader
	Mosaik::Read mr;
	while(reader.LoadNextRead(mr)) {
		if(mr.Mate1.Bases.Length() != 0) ProcessMate(1, mr.Name, isColorspace, mr.Mate1, isPairedEnd);
		if(mr.Mate2.Bases.Length() != 0) ProcessMate(2, mr.Name, isColorspace, mr.Mate2, isPairedEnd);
	}

	if(!mFlags.IsScreenEnabled) printf("finished.\n");

	// clean up
	reader.Close();
	if(mFlags.IsFastqEnabled) fclose(mStreams.fastq);
}

// processes the mates according to the chosen file format
void CMosaikText::ProcessMate(const unsigned char mateNum, const CMosaikString& readName, const bool isColorspace, Mosaik::Mate& mate, const bool isPairedEnd) {

	// increment the base qualities
	mate.Qualities.Increment(33);

	if(isColorspace) {
		mCS.ConvertReadPseudoColorspaceToColorspace(mate.Bases);
		mate.Bases.Prepend(mate.SolidPrefixTransition, 2);
		if(mFlags.IsFastqEnabled) mate.Qualities.Prepend("!?", 2);
		else mate.Qualities.Prepend("?", 1);
	}

	if(mFlags.IsFastqEnabled) {
		fprintf(mStreams.fastq, "@%s", readName.CData());
		if(isPairedEnd) fprintf(mStreams.fastq, " (mate %u)", mateNum);
		fprintf(mStreams.fastq, "\n%s\n+\n%s\n", mate.Bases.CData(), mate.Qualities.CData());
	}

	if(mFlags.IsScreenEnabled) {
		printf("%s", readName.CData());
		if(isPairedEnd) printf(" (mate %u)", mateNum);
		printf("\n%s\n%s\n\n", mate.Bases.CData(), mate.Qualities.CData());
	}
}

// writes the current alignment to the SAM output file
void CMosaikText::WriteSamEntry(const CMosaikString& readName, const string& readGroupID, const vector<Alignment>::iterator& alIter) {

	// =================
	// set the SAM flags
	// =================

	unsigned int flag                = 0;
	unsigned int queryPosition5Prime = 0;
	unsigned int matePosition5Prime  = 0;
	int insertSize                   = 0;

	if(alIter->IsPairedEnd) {

		flag |= BAM_SEQUENCED_AS_PAIRS;

		// first or second mate?
		flag |= (alIter->IsFirstMate ? BAM_QUERY_FIRST_MATE : BAM_QUERY_SECOND_MATE);

		if(alIter->IsResolvedAsPair) {

			flag |= BAM_PROPER_PAIR;
			if(alIter->IsMateReverseStrand) flag |= BAM_MATE_REVERSE_COMPLEMENT;

			// sanity check
			if(alIter->ReferenceIndex != alIter->MateReferenceIndex) {
				printf("ERROR: The resolved paired-end reads occur on different reference sequences.\n");
				exit(1);
			}

			// set the 5' coordinates
			queryPosition5Prime = (alIter->IsReverseStrand     ? alIter->ReferenceEnd     : alIter->ReferenceBegin);
			matePosition5Prime  = (alIter->IsMateReverseStrand ? alIter->MateReferenceEnd : alIter->MateReferenceBegin);

			// calculate the insert size
			insertSize = matePosition5Prime - queryPosition5Prime;

		} else flag |= BAM_MATE_UNMAPPED;
	}

	if(alIter->IsReverseStrand) flag |= BAM_QUERY_REVERSE_COMPLEMENT;

	// ==========================
	// construct the cigar string
	// ==========================

	char* pCigar = mCigarBuffer;
	const char* pReference = alIter->Reference.CData();
	const char* pQuery     = alIter->Query.CData();

	const unsigned short numBases = alIter->Reference.Length();
	unsigned short currentPos    = 0;
	unsigned int numBufferBytes  = 0;

	while(currentPos < numBases) {

		unsigned short testPos = currentPos;
		unsigned short operationLength = 0;
		int numWritten = 0;

		if((pReference[currentPos] != '-') && (pQuery[currentPos] != '-')) {

			while((pReference[testPos] != '-') && (pQuery[testPos] != '-') && (testPos < numBases)) {
				++testPos;					
				++operationLength;
			}

			numWritten = sprintf_s(pCigar, CIGAR_BUFFER_SIZE, "%uM", operationLength);

		} else if(pReference[currentPos] == '-') {

			while((pReference[testPos] == '-') && (testPos < numBases)) {
				++testPos;					
				++operationLength;
			}

			numWritten = sprintf_s(pCigar, CIGAR_BUFFER_SIZE, "%uI", operationLength);

		} else if(pQuery[currentPos] == '-') {

			while((pQuery[testPos] == '-') && (testPos < numBases)) {
				++testPos;					
				++operationLength;
			}

			numWritten = sprintf_s(pCigar, CIGAR_BUFFER_SIZE, "%uD", operationLength);

		} else {
			cout << "ERROR: CIGAR string generation failed." << endl;
			exit(1);
		}

		// increment our position
		pCigar         += numWritten;
		numBufferBytes += numWritten;
		currentPos     += operationLength;

		// make sure aren't creating a buffer overflow
		if(numBufferBytes >= CIGAR_BUFFER_SIZE) {
			printf("ERROR: buffer overflow detected when creating the cigar string.\n");
			exit(1);
		}
	}

	*pCigar = 0;

	// ===================
	// write the alignment
	// ===================

	// B7_591:4:96:693:509	73	seq1	1	99	36M	*	0	0	CACTAGTGGCTCATTGTAAATGTGTGGTTTAACTCG	<<<<<<<<<<<<<<<;<<<<<<<<<5<<<<<;:<;7
	// <QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL>

	// remove the gaps from the read
	CMosaikString query(alIter->Query);
	query.Remove('-');

	// shift the base qualities 
	CMosaikString bq(alIter->BaseQualities);
	bq.Increment(33);

	gzprintf(mStreams.sam, "%s\t%u\t%s\t%u\t%u\t%s\t", readName.CData(), flag, alIter->ReferenceName, 
		alIter->ReferenceBegin + 1, alIter->Quality, mCigarBuffer);

	if(alIter->IsResolvedAsPair) {
		// N.B. we already checked that the reference indexes were identical
		gzprintf(mStreams.sam, "=\t%u\t%u\t", alIter->MateReferenceBegin + 1, insertSize);
	} else gzprintf(mStreams.sam, "*\t*\t*\t");

	gzprintf(mStreams.sam, "%s\t%s\tRG:Z:%s\tNM:i:%u\n", query.CData(), bq.CData(), readGroupID.c_str(), alIter->NumMismatches);	
}
