// ***************************************************************************
// CMosaikCoverage - exports coverage graphs using alignment archives.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "MosaikCoverage.h"

// constructor
CMosaikCoverage::CMosaikCoverage(void)
: mpCoverageArray(NULL)
, mDisableOutput(false)
, mEvaluateUniqueReadsOnly(false)
, mCreatePostscriptGraphs(false)
, mUseAlignmentQualityFilter(false)
, mNumRefSeqs(0)
{
}

// destructor
CMosaikCoverage::~CMosaikCoverage(void) {
	for(unsigned int i = 0; i < mNumRefSeqs; i++) if(mpCoverageArray[i]) delete [] mpCoverageArray[i];
	if(mpCoverageArray) delete [] mpCoverageArray;
}

// toggles the creation of output files
void CMosaikCoverage::DisableOutput(void) {
	cout << "- disabling output of coverage files or graphs." << endl;
	mDisableOutput = true;
}

// toggles the creation of graphs via gnuplot
void CMosaikCoverage::EnableGraphCreation(void) {
	cout << "- enabling graph creation." << endl;
	mCreatePostscriptGraphs = true;
}

// when triggered, the coverage calculation will only include unique reads
void CMosaikCoverage::EvaluateUniqueReadsOnly(void) { 
	cout << "- evaluating unique reads only." << endl;
	mEvaluateUniqueReadsOnly = true;
}

// parses the specified MOSAIK alignment file and matching anchors file
void CMosaikCoverage::ParseMosaikAlignmentFile(const vector<string>& alignmentFilenames, const string& anchorsFilename) {

	// retrieve the reference sequence metadata
	MosaikReadFormat::CReferenceSequenceReader refseq;
	refseq.Open(anchorsFilename);
	refseq.GetReferenceSequences(mReferenceSequences);
	mNumRefSeqs = mReferenceSequences.size();

	// initialize the coverage array
	CConsole::Heading(); printf("- initializing coverage array... "); CConsole::Reset();
	fflush(stdout);

	mpCoverageArray = new unsigned int*[mNumRefSeqs];
	unsigned int* allocatedLength = new unsigned int[mNumRefSeqs];

	vector<ReferenceSequence>::iterator refIter = mReferenceSequences.begin();
	for(unsigned int i = 0; i < mNumRefSeqs; i++, refIter++) {
		const unsigned int numBases = refIter->NumBases;
		mpCoverageArray[i] = new unsigned int[numBases];
		allocatedLength[i] = numBases;
		uninitialized_fill(mpCoverageArray[i], mpCoverageArray[i] + numBases, 0);

		// count the number of non-masked bases
		{
			string tempBases;
			refseq.GetReferenceSequence(refIter->Name, tempBases);
			char* s = (char*)tempBases.data();
			unsigned int numBases = 0;
			for(unsigned int j = 0; j < (unsigned int)tempBases.size(); j++) {
				if((s[j] != 'X') && (s[j] != 'N')) numBases++;
			}
			refIter->NumBases = numBases;
		}
	}

	printf("finished.\n");

	// parse all of the reads
	CConsole::Heading(); printf("- parsing read archive(s):\n"); CConsole::Reset();

	// ===================================
	// parse all of the alignment archives
	// ===================================

	vector<string>::const_iterator fnIter;
	for(fnIter = alignmentFilenames.begin(); fnIter != alignmentFilenames.end(); ++fnIter) {

		printf("  %s... ", fnIter->c_str());
		fflush(stdout);

		// open the alignment reader
		MosaikReadFormat::CAlignmentReader ar;
		ar.Open(*fnIter);

		// update the aligned read totals
		vector<ReferenceSequence>* pReferenceSequences = ar.GetReferenceSequences();
		for(unsigned int i = 0; i < mNumRefSeqs; i++) {
			mReferenceSequences[i].NumAligned = pReferenceSequences->at(i).NumAligned;
		}

		vector<Alignment>::const_iterator alIter;
		Mosaik::AlignedRead r;
		while(ar.LoadNextRead(r)) {

			// skip non-unique reads if requested
			const unsigned int numMate1Alignments = r.Mate1Alignments.size();
			const unsigned int numMate2Alignments = r.Mate2Alignments.size();

			if(mEvaluateUniqueReadsOnly && ((numMate1Alignments > 1) || (numMate2Alignments > 1))) continue;

			// handle the mate 1 alignments
			for(alIter = r.Mate1Alignments.begin(); alIter != r.Mate1Alignments.end(); alIter++) {
				for(unsigned int j = alIter->ReferenceBegin; j <= alIter->ReferenceEnd; j++) {
					if(j >= allocatedLength[alIter->ReferenceIndex]) {
						cout << "ERROR: Tried to write past the end of the coverage array for reference sequence " 
							<< alIter->ReferenceIndex << ": reference begin: " << alIter->ReferenceBegin << ", end: " 
							<< alIter->ReferenceEnd << " (mate 1: " << r.Name << ")" << endl;
						exit(1);
					}

					if(mUseAlignmentQualityFilter) {
					  if(alIter->Quality >= mAlignmentQualityThreshold) mpCoverageArray[alIter->ReferenceIndex][j]++;
					} else mpCoverageArray[alIter->ReferenceIndex][j]++;
				}
			}

			// handle the mate 2 alignments
			for(alIter = r.Mate2Alignments.begin(); alIter != r.Mate2Alignments.end(); alIter++) {
				for(unsigned int j = alIter->ReferenceBegin; j <= alIter->ReferenceEnd; j++) {
					if(j >= allocatedLength[alIter->ReferenceIndex]) {
						cout << "ERROR: Tried to write past the end of the coverage array for reference sequence " 
							<< alIter->ReferenceIndex << ": reference begin: " << alIter->ReferenceBegin << ", end: " 
							<< alIter->ReferenceEnd << " (mate 2: " << r.Name << ")" << endl;
						exit(1);
					}

					if(mUseAlignmentQualityFilter) {
					  if(alIter->Quality >= mAlignmentQualityThreshold) mpCoverageArray[alIter->ReferenceIndex][j]++;
					} else mpCoverageArray[alIter->ReferenceIndex][j]++;
				}
			}
		}

		// clean up
		ar.Close();
		printf("finished.\n");
	}

	// save the coverage files
	delete [] allocatedLength;
}

// saves the coverage files to the specified output directory
void CMosaikCoverage::SaveCoverage(const unsigned char minCoverage) {

	printf("\n");
	if(!mDisableOutput) {
		CConsole::Heading(); printf("- writing coverage files:"); CConsole::Reset();
	} else {
		CConsole::Heading(); printf("- calculating coverage statistics:"); CConsole::Reset();
	}
	printf("\n");

	// create our coverage files
	ostringstream sb;
	vector<ReferenceSequence>::const_iterator refIter;
	unsigned short currentRefSeq = 0;

	uint64_t numTotalCoveredBases = 0;
	uint64_t numTotalBases        = 0;
	uint64_t totalCoverageCount   = 0;

	for(refIter = mReferenceSequences.begin(); refIter != mReferenceSequences.end(); refIter++, currentRefSeq++) {

		// continue on to the next reference sequence if this one doesn't have
		// coverage information
		const unsigned int anchorLen = refIter->End - refIter->Begin + 1;

		unsigned int numCoveredBases = 0;
		uint64_t coverageCount       = 0;

		if(refIter->NumAligned > 0) {

			// convert our anchor name
			string anchorName = refIter->Name;
			char* pAnchorName = (char*)anchorName.c_str();
			for(unsigned int k = 0; k < (unsigned int)anchorName.size(); k++)
				if((pAnchorName[k] < 33) || (pAnchorName[k] > 122)) pAnchorName[k] = '_';

			// open the coverage file
			sb.str("");
			if(mOutputDirectory.empty()) sb << anchorName << ".cov";
			else sb << mOutputDirectory << OS_DIRECTORY_SEPARATOR << anchorName << ".cov";

			string outputFilename = sb.str();

			FILE* outStream = NULL;
			if(!mDisableOutput) {
				if(fopen_s(&outStream, outputFilename.c_str(), "wb") != 0) {
					cout << "ERROR: Could not open the coverage file (" << outputFilename << ") for writing." << endl;
					exit(1);
				}

				// write the header
				fprintf(outStream, "# MOSAIK coverage file for %s\n", refIter->Name.c_str());
				fprintf(outStream, "# reference-pos coverage\n");
			}

			// collect statistics and write coverage data to file
			for(unsigned int j = 0; j < anchorLen; j++) {
			        const unsigned int count = mpCoverageArray[currentRefSeq][j];
				
				if(!mDisableOutput) fprintf(outStream, "%u %u\n", j + 1, count);

				if(count >= minCoverage) {
					numCoveredBases++;
					coverageCount += count;
				}
			}

			// close the coverage file
			if(!mDisableOutput) {

				fclose(outStream);

				// create postscript graphs with gnuplot
				if(mCreatePostscriptGraphs) {

					sb.str("");
					sb << anchorName << ".cov";
					string outputFilenameStub = sb.str();

					// escape the anchor title
					string anchorTitle = refIter->Name;

					ofstream gnuplot("MosaikCoverage.plt");

					if(!mOutputDirectory.empty()) gnuplot << "cd \"" << mOutputDirectory << "\"" << endl;
					gnuplot << "set term postscript color" << endl;
					gnuplot << "set output \"" << anchorName << ".ps\"" << endl;
					gnuplot << "plot [1:" << anchorLen << "] \"" << outputFilenameStub << "\" title \"" << anchorTitle << "\"" << endl;
					gnuplot << "set term postscript color" << endl;
					gnuplot << "set output \"" << anchorName << ".ps\"" << endl;
					gnuplot << "set title \"Reference sequence coverage\"" << endl;
					gnuplot << "set xlabel \"reference sequence position (bp)\"" << endl;
					gnuplot << "set ylabel \"coverage\"" << endl;
					gnuplot << "set style data lines" << endl;
					gnuplot << "replot" << endl;
					gnuplot.close();
					system("gnuplot MosaikCoverage.plt");
					rm("MosaikCoverage.plt");

#ifndef WIN32
					sb.str("");
					sb << "ps2pdf ";
					if(mOutputDirectory.empty()) sb << anchorName << ".ps " << anchorName << ".pdf";
					else sb << mOutputDirectory << OS_DIRECTORY_SEPARATOR << anchorName << ".ps " << mOutputDirectory 
						<< OS_DIRECTORY_SEPARATOR << anchorName << ".pdf";
					string pdfConvertCommand = sb.str();
					system(pdfConvertCommand.c_str());
#endif
				}
			}
		}

		// ======================
		// display our statistics
		// ======================

		numTotalBases        += refIter->NumBases;
		numTotalCoveredBases += numCoveredBases;
		totalCoverageCount   += coverageCount;
		const double meanCoverage = coverageCount / (double)refIter->NumBases;

		printf("* coverage statistics for %s (%u bp): %u bp (%4.1f %%), mean: %.1fx\n", refIter->Name.c_str(), refIter->NumBases, numCoveredBases, (numCoveredBases / (double)refIter->NumBases) * 100.0, meanCoverage);
	}

	const double meanTotalCoverage = totalCoverageCount / (double)numTotalBases;

	printf("\n");
	CConsole::Heading(); printf("- total coverage (%llu bp): ", (unsigned long long)numTotalBases); CConsole::Reset();
	printf("%llu bp (%4.1f %%), mean: %.1fx\n", (unsigned long long)numTotalCoveredBases, (numTotalCoveredBases / (double)numTotalBases) * 100.0, meanTotalCoverage);
}

// sets the minimum alignment quality threshold
void CMosaikCoverage::SetMinAlignmentQuality(const unsigned char minAQ) {
  cout << "- setting the minimum alignment quality to: " << (int)minAQ << endl;
        mUseAlignmentQualityFilter = true;
	mAlignmentQualityThreshold = minAQ;
}

// sets the output directory
void CMosaikCoverage::SetOutputDirectory(const string& directory) {
	mOutputDirectory = directory;
}
