// ***************************************************************************
// CMosaikAssembler - consolidates gaps in the underlying alignments and
//                    exports the data set into an assembly format.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "MosaikAssembler.h"

// constructor
CMosaikAssembler::CMosaikAssembler(AssemblyFormatType format, unsigned char referenceBaseQuality)
: mpUngap2Gap(NULL)
, mUngap2GapLen(0)
, mReadPadBuffer(NULL)
, mReadPadBufferLen(0)
{
	mSettings.AssemblyFormat               = format;
	mSettings.ReferenceSequenceBaseQuality = referenceBaseQuality;
}

// destructor
CMosaikAssembler::~CMosaikAssembler(void) {
	if(mReadPadBuffer) delete [] mReadPadBuffer;
	if(mpUngap2Gap)    delete [] mpUngap2Gap;
}

// inserts observed gaps to the specified read
void CMosaikAssembler::AddGapsToRead(const vector<GapInfo>& gaps, const Alignment& al, CMosaikString& gappedRead) {

	//printf("\nDEBUG: NEW READ: %s\n", al.Name.CData());
	//printf("-------------------------------------------\n");

	// return if there are no gaps in the reference sequence
	if(gaps.empty()) {
		gappedRead = al.Query;
		return;
	}

	// retrieve the first and last gap positions 
	const unsigned int firstGapPosition = gaps.front().Position;
	const unsigned int lastGapPosition  = gaps.back().Position;

	// handle cases where all gaps occur before or after read
	if((al.ReferenceEnd < firstGapPosition) || (al.ReferenceBegin > lastGapPosition)) {
		gappedRead = al.Query;
		return;
	}

	// check our read pad buffer
	const unsigned int numPairwiseBases = al.Query.Length();
	CheckReadPadBuffer(numPairwiseBases);

	// ===============================
	// convert the read to a pad array
	// ===============================

	const char* pReference = al.Reference.CData();
	const char* pQuery     = al.Query.CData();

	unsigned int ungappedRefPos = al.ReferenceBegin;

	for(unsigned int i = 0; i < numPairwiseBases; i++) {
		mReadPadBuffer[i].ReferenceNucleotide = pReference[i];
		mReadPadBuffer[i].QueryNucleotide     = pQuery[i];
		mReadPadBuffer[i].Position            = ungappedRefPos;
		mReadPadBuffer[i].Length              = 0;

		// increment the reference position
		if(pReference[i + 1] != '-') ungappedRefPos++;
	}

	// skip all gaps before the reference start
	vector<GapInfo>::const_iterator gvIter = mGapIter;
	while(gvIter->Position < al.ReferenceBegin) gvIter++;
	mGapIter = gvIter;

	// ====================================
	// add all gaps until the reference end
	// ====================================

	unsigned int addedGaps         = 0;
	unsigned int currentRefPos     = 0;

	while((gvIter != gaps.end()) && (gvIter->Position < al.ReferenceEnd)) {

		//printf("\nDEBUG: handling gap position: %u, length: %u\n", gvIter->Position, gvIter->Length);

		// find the insertion point
		while(mReadPadBuffer[currentRefPos].Position < gvIter->Position) currentRefPos++;

		//printf("DEBUG: insertion point positions: %u (%u), current reference position: %u, alleles: %c (ref: %c)\n", mReadPadBuffer[currentRefPos].Position, gvIter->Position, currentRefPos, mReadPadBuffer[currentRefPos].QueryNucleotide, mReadPadBuffer[currentRefPos].ReferenceNucleotide);

		// count the number of observed gaps
		unsigned int numObservedGaps = 0;
		unsigned int testRefPos = currentRefPos + 1;
		while((testRefPos < numPairwiseBases) && (mReadPadBuffer[testRefPos].ReferenceNucleotide == '-')) {
			numObservedGaps++;
			testRefPos++;
		}

		// insert the appropriate number of gaps
		const unsigned int numInsertedGaps = gvIter->Length - numObservedGaps;

		if((numInsertedGaps > 0) && (currentRefPos < numPairwiseBases)) {
			mReadPadBuffer[currentRefPos].Length = numInsertedGaps;
			addedGaps += numInsertedGaps;
			//printf("DEBUG: added %u gaps at position %u. Allele before gap: %c (ref: %c)\n", numInsertedGaps, currentRefPos, mReadPadBuffer[currentRefPos].QueryNucleotide, mReadPadBuffer[currentRefPos].ReferenceNucleotide);
		}

		// continue with the next gap
		currentRefPos += numObservedGaps;

		gvIter++;
	}

	// ======================================
	// convert our read pad array to a string
	// ======================================

	const unsigned int numGappedReadBases = numPairwiseBases + addedGaps;
	gappedRead.Reserve(numGappedReadBases);
	gappedRead.SetLength(numGappedReadBases);

	char* pGappedRef = gappedRead.Data();
	char* pUpdateRef = pGappedRef;

	for(unsigned int i = 0; i < numPairwiseBases; i++) {
		*pUpdateRef = mReadPadBuffer[i].QueryNucleotide;
		pUpdateRef++;

		const unsigned short gapLength = mReadPadBuffer[i].Length;
		if(gapLength > 0) {
			memset(pUpdateRef, '-', gapLength);
			pUpdateRef += gapLength;
		}
	}
}

// inserts observed gaps to the reference sequence
void CMosaikAssembler::AddGapsToReferenceSequence(const vector<GapInfo>& gaps, string& referenceSequence, CMosaikString& gappedReferenceSequence) {

	//printf("\nDEBUG: NEW REFERENCE:\n");

	// return if there are no gaps in the reference sequence
	const unsigned int numGaps = (unsigned int)gaps.size();
	if(numGaps == 0) {
		gappedReferenceSequence = referenceSequence.c_str();	
		return;
	}

	// initialize
	const char* pReference = referenceSequence.data();
	const unsigned int numReferenceBases = (unsigned int)referenceSequence.size();
	ReferencePadInfo* refPads = NULL;

	try {
		refPads = new ReferencePadInfo[numReferenceBases];
	} catch(bad_alloc) {
		printf("ERROR: Unable to allocate memory for %u bases when adding gaps to the reference sequence.\n", numReferenceBases);
		exit(1);
	}

	for(unsigned int i = 0; i < numReferenceBases; i++) {
		refPads[i].Nucleotide = pReference[i];
		refPads[i].Position   = i;
		refPads[i].Length     = 0;
	}

	// add the gaps to our reference sequence
	unsigned int addedGaps = 0;
	vector<GapInfo>::const_iterator gvIter;
	for(gvIter = gaps.begin(); gvIter != gaps.end(); gvIter++) {
		//printf("DEBUG: gap position: %u, length: %u\n", gvIter->Position, gvIter->Length);
		refPads[gvIter->Position].Length = gvIter->Length;
		addedGaps += gvIter->Length;
	}

	// convert our reference pad array to a string
	const unsigned int numGappedReferenceBases = numReferenceBases + addedGaps;
	gappedReferenceSequence.Reserve(numGappedReferenceBases);
	gappedReferenceSequence.SetLength(numGappedReferenceBases);
	char* pGappedRef = gappedReferenceSequence.Data();

	char* pUpdateRef = pGappedRef;
	for(unsigned int i = 0; i < numReferenceBases; i++) {
		*pUpdateRef = refPads[i].Nucleotide;
		pUpdateRef++;

		const unsigned short gapLength = refPads[i].Length;
		if(gapLength > 0) {
			memset(pUpdateRef, '-', gapLength);
			pUpdateRef += gapLength;
		}
	}

	// DEBUG
	//printf("DEBUG: gapped reference sequence:\n%s\n", gappedReferenceSequence.CData());

	// clean up
	if(refPads) delete [] refPads;
}

// assembles the specified alignment file
void CMosaikAssembler::Assemble(const string& referenceFilename, const string& alignmentFilename, const string& outputFilenameStub) {

	// =====================================================
	// find out which reference sequences have aligned reads
	// =====================================================

	// open our alignment archive
	MosaikReadFormat::CAlignmentReader reader;
	reader.Open(alignmentFilename);

	// check if we're creating ACE files
	const bool usingACE = (mSettings.AssemblyFormat == AssemblyFormat_ACE ? true : false);

	// load the first alignment
	Alignment al;
	reader.LoadNextAlignment(al);

	vector<ReferenceSequence>* pReferenceSequences   = reader.GetReferenceSequences();
	vector<vector<GapInfo> >* pReferenceSequenceGaps = reader.GetReferenceSequenceGaps();

	vector<ReferenceSequence>::const_iterator arIter;
	vector<vector<GapInfo> >::const_iterator rsgIter;

	// show the user which reference sequences will be processed
	cout << "===============================================================" << endl;
	cout << "alignment count   reference sequence" << endl;
	cout << "---------------------------------------------------------------" << endl;

	bool foundAlignedReads = false;
	for(arIter = pReferenceSequences->begin(); arIter != pReferenceSequences->end(); arIter++) {

		// skip if not our region of interest
		if(mFlags.HasRegionOfInterest && (arIter->Name != mSettings.RoiName)) continue;

		if(arIter->NumAligned > 0) {
			printf("%15llu   %s\n", (unsigned long long)arIter->NumAligned, arIter->Name.c_str());
			foundAlignedReads = true;
		}
	}

	// quit if there is nothing to assemble
	if(!foundAlignedReads) {
		cout << "(no aligned reads were found)" << endl;
		return;
	}

	// ===============================
	// process each reference sequence
	// ===============================

	// open our reference sequence archive
	MosaikReadFormat::CReferenceSequenceReader rsr;
	rsr.Open(referenceFilename);

	CMosaikString gappedPairwiseQuery;
	string referenceSequence;
	CAbstractAssemblyFormat* af = NULL;

	// DEBUG
	//cout << "DEBUG: reference sequence gaps size: " << pReferenceSequenceGaps->size() << endl;

	// TODO: make sure we have as many reference sequence gaps as we reference sequences
	//       maybe this is guaranteed by the alignment writer?

	unsigned short currentReferenceIndex = 0;
	rsgIter = pReferenceSequenceGaps->begin();
	for(arIter = pReferenceSequences->begin(); arIter != pReferenceSequences->end(); arIter++, rsgIter++, currentReferenceIndex++) {

		// skip empty reference sequences
		if(arIter->NumAligned == 0) continue;

		// skip if not our region of interest
		if(mFlags.HasRegionOfInterest && (arIter->Name != mSettings.RoiName)) continue;

		printf("\n");
		CConsole::Heading();
		printf("Processing reference sequence %s:\n", arIter->Name.c_str());
		CConsole::Reset();

		// get the reference sequence
		rsr.GetReferenceSequence(arIter->Name, referenceSequence);
		const unsigned int ungappedRefLength = referenceSequence.size();

		// insert gaps into the reference
		printf("- inserting gaps into reference sequence... ");
		fflush(stdout);
		CMosaikString gappedReferenceSequence;
		AddGapsToReferenceSequence(*rsgIter, referenceSequence, gappedReferenceSequence);
		printf("finished.\n");

		// create the ungapped to gapped conversion table
		printf("- creating ungapped to gapped conversion table... ");
		fflush(stdout);

		PopulateUngappedToGappedVector(gappedReferenceSequence, ungappedRefLength);

		printf("finished.\n");

		// derive our assembly and FASTA filenames
		CMosaikString assemblyFilename = outputFilenameStub.c_str();
		assemblyFilename.Append("_");
		assemblyFilename.Append(arIter->Name.c_str());
		CMosaikString fasta = assemblyFilename;

		if(usingACE) assemblyFilename.Append(".ace");
		else assemblyFilename.Append(".gig");

		fasta.Append(".fasta");
		CMosaikString fastaQual = fasta;
		fastaQual.Append(".qual");

		// initialize the correct assembly format
		if(usingACE) af = new CAce(mSettings.ReferenceSequenceBaseQuality);
		else af = new CGigaBayesFormat(mSettings.ReferenceSequenceBaseQuality);

		// create our assembly file and write the header
		printf("- writing assembly header... ");
		fflush(stdout);

		af->Open(assemblyFilename);
		af->SetUngappedToGappedVector(mpUngap2Gap);
		af->SaveHeader(gappedReferenceSequence, arIter->Name, ungappedRefLength, arIter->NumAligned);

		printf("finished.\n");

		printf("- locating first read for this reference sequence... ");
		fflush(stdout);

		if(mFlags.HasRegionOfInterest) reader.Jump(currentReferenceIndex, 0);

		bool haveMoreReads = true;
		uint64_t numSkippedReads = 0;

		do {
			if(al.ReferenceIndex == currentReferenceIndex) break;
			numSkippedReads++;
		} while((haveMoreReads = reader.LoadNextAlignment(al)));

		// technically this should never happen
		if(!haveMoreReads) {
			cout << "ERROR: Expected more reads to assemble." << endl;
			exit(1);
		}

		printf("finished.\n");

		// ================================
		// add gaps and save the alignments
		// ================================

		uint64_t numSavedAlignments = 0;
		const uint64_t numTotalAlignments = arIter->NumAligned;

		// initialize our iterators
		mGapIter = rsgIter->begin();

		// open our FASTA files
		FILE *fastaStream = NULL, *fastaQualStream = NULL;

		if(mFlags.EnableFastaFileCreation) {
			if(fopen_s(&fastaStream, fasta.CData(), "wb") != 0) {
				cout << "ERROR: Could not open the FASTA reads file for writing." << endl;
				exit(1);
			}

			if(fopen_s(&fastaQualStream, fastaQual.CData(), "wb") != 0) {
				cout << "ERROR: Could not open the FASTA base qualities file for writing." << endl;
				exit(1);
			}
		}

		// show the progress bar
		CConsole::Heading(); 
		cout << endl << "- saving alignments from " << arIter->Name << ":" << endl;
		CConsole::Reset(); 
		CProgressBar<uint64_t>::StartThread(&numSavedAlignments, 0, numTotalAlignments, "alignments");

		for(; numSavedAlignments < numTotalAlignments; numSavedAlignments++) {

			// add gaps and save the read
			AddGapsToRead(*rsgIter, al, gappedPairwiseQuery);
			af->SaveRead(al, gappedPairwiseQuery);

			// save the FASTA files for gigaBayes
			if(mFlags.EnableFastaFileCreation) {
				CMosaikString ungappedRead = al.Query;
				CMosaikString ungappedBaseQualities = al.BaseQualities;
				ungappedRead.Remove('-');
				if(al.IsReverseStrand) {
					CSequenceUtilities::GetReverseComplement(ungappedRead.Data(), ungappedRead.Length());
					ungappedBaseQualities.Reverse();
				}

				fprintf(fastaStream, ">%s\n%s\n", al.Name.CData(), ungappedRead.CData());
				fprintf(fastaQualStream, ">%s\n", al.Name.CData());
				for(unsigned int k = 0; k < ungappedBaseQualities.Length(); k++) fprintf(fastaQualStream, "%u ", ungappedBaseQualities[k]);
				fprintf(fastaQualStream, "\n");
			}

			// get the next alignment (al)
			haveMoreReads = reader.LoadNextAlignment(al);
		}

		// wait for the progress bar to end
		CProgressBar<uint64_t>::WaitThread();

		// close and delete the assembly file format
		af->Close();
		delete af;

		// close the FASTA files
		if(mFlags.EnableFastaFileCreation) {
			fclose(fastaStream);
			fclose(fastaQualStream);
		}
	}

	// close the MOSAIK archives
	reader.Close();
	rsr.Close();
}

// re-allocated the read pad buffer if more space is needed
void CMosaikAssembler::CheckReadPadBuffer(const unsigned int readLength) {
	if(readLength > mReadPadBufferLen) {
		unsigned int newLength = readLength + 10;
		if(mReadPadBuffer) delete [] mReadPadBuffer;
		mReadPadBuffer = new ReadPadInfo[newLength];
		mReadPadBufferLen = newLength;
	}
}

// enables FASTA file creation for gigaBayes
void CMosaikAssembler::EnableFastaFileCreation(void) {
	mFlags.EnableFastaFileCreation = true;
}

// builds a ungapped to gapped conversion vector from the specified reference sequence
void CMosaikAssembler::PopulateUngappedToGappedVector(const CMosaikString& reference, const unsigned int ungappedLength) {

	// redimension our vector (if needed)
	if(ungappedLength > mUngap2GapLen) {

		mUngap2GapLen = ungappedLength;
		if(mpUngap2Gap) delete [] mpUngap2Gap;

		try {
			mpUngap2Gap = new unsigned int[mUngap2GapLen];
		} catch(bad_alloc) {
			printf("ERROR: Unable to allocate %u unsigned integers for the ungapped-to-gapped coordinate vector.\n", mUngap2GapLen);
			exit(1);
		}
	}

	// initialize
	unsigned int gappedLength = reference.Length();	
	unsigned int ungappedPosition = 0;
	unsigned int gappedPosition   = 0;

	const char* pReference = reference.CData();

	// build the map
	bool reachedEnd = false;
	while(gappedPosition < gappedLength) {

		// skip over the gaps
		while(pReference[gappedPosition] == '-') {
			gappedPosition++;
			if(gappedPosition >= gappedLength) {
				reachedEnd = true;
				break;
			}
		}

		// quit if we have reached the end of the reference
		if(reachedEnd) break;

		//printf("DEBUG: gapped: %u, ungapped: %u, ref: %c\n", gappedPosition, ungappedPosition, pReference[gappedPosition]);

		mpUngap2Gap[ungappedPosition] = gappedPosition;

		// increment our positions
		gappedPosition++;
		ungappedPosition++;
	}
}

// selects the specified region of interest
void CMosaikAssembler::SelectRegionOfInterest(const string& referenceName) {
	mFlags.HasRegionOfInterest = true;
	mSettings.RoiName = referenceName;
}
