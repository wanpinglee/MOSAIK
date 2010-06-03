// ***************************************************************************
// CNaiveAlignmentSet - essentially a very naive mechanism for ensuring that
//                      duplicates alignments are not reported.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "NaiveAlignmentSet.h"

// constructor
CNaiveAlignmentSet::CNaiveAlignmentSet(unsigned int refLen, bool usingIllumina) 
: mHasLongAlignment(false)
, mAlignmentQuality(usingIllumina, refLen)
{}

// destructor
CNaiveAlignmentSet::~CNaiveAlignmentSet() {}

// adds an alignment to the set
bool CNaiveAlignmentSet::Add(Alignment& al) {

	// check if this is a long alignment
	const unsigned short pairwiseLength = (unsigned short)al.Reference.Length();
	if((al.QueryEnd > 255) || (pairwiseLength > 255)) mHasLongAlignment = true;

	// check to see if any of the other entries are similar
	bool foundSubset = false;
	AlignmentSet::iterator setIter;

	// TODO: update this with an interval tree
	if(!mAlignments.empty()) {
		for(setIter = mAlignments.begin(); setIter != mAlignments.end(); setIter++) {
			if(CheckOverlap(al, setIter)) {
				foundSubset = true;
				break;
			}
		}
	}

	// add the alignment to the alignment set
	if(!foundSubset) {
		mAlignments.push_back(al);
		return true;
	}

	// handle the subset: choose the bigger alignment
	unsigned int alLen = al.ReferenceEnd       - al.ReferenceBegin       + 1;
	unsigned int axLen = setIter->ReferenceEnd - setIter->ReferenceBegin + 1;
	if(alLen > axLen) *setIter = al;

	return false;
}

// calculates the alignment qualities for each alignment in the set
void CNaiveAlignmentSet::CalculateAlignmentQualities(const bool calculateCorrectionCoefficient, const unsigned short minSpanLength) {

	// sort the mhp occupancy lists
	if(calculateCorrectionCoefficient) {
		mFwdMhpOccupancyList.sort();
		mRevMhpOccupancyList.sort();
	}

	// calculate the alignment qualities
	AlignmentSet::iterator setIter;
	for(setIter = mAlignments.begin(); setIter != mAlignments.end(); ++setIter) {
		mAlignmentQuality.CalculateQuality(setIter);

		if(calculateCorrectionCoefficient) {
			const double correctionCoefficient = CalculateCorrectionCoefficient(setIter->QueryBegin, setIter->QueryEnd, 
				(setIter->IsReverseStrand ? mRevMhpOccupancyList : mFwdMhpOccupancyList), minSpanLength);

			// modify the alignment quality if the correction coefficient kicks in
			if(correctionCoefficient < 1.0) {
				const double Pcorrect = 1.0 - pow(10.0, -setIter->Quality / 10.0);
				setIter->Quality = (unsigned char)(-10.0 * log10(1.0 - correctionCoefficient * Pcorrect) + 0.5);
			}
		}
	}
}

// calculates the correction coefficient
double CNaiveAlignmentSet::CalculateCorrectionCoefficient(const unsigned short queryBegin, const unsigned short queryEnd, const MhpOccupancyList& mhpOccupancyList, const unsigned short minSpanLength) {

	// process the first mhp occupancy position that is within [queryBegin, queryEnd]
	MhpOccupancyList::const_iterator mhpIter = mhpOccupancyList.begin();
	while((mhpIter->Begin < queryBegin) || (mhpIter->End > queryEnd)) ++mhpIter;

	MhpOccupancyRegion mor(mhpIter->Begin, mhpIter->End, mhpIter->Occupancy);
	++mhpIter;

	// process the remaining mhp occupancy positions within [queryBegin, queryEnd]
	for(; mhpIter != mhpOccupancyList.end(); ++mhpIter) {

		// check if we have a long enough region
		if((mor.End - mor.Begin + 1) >= minSpanLength) break;

		// skip positions that lie outside of [queryBegin, queryEnd]
		if((mhpIter->Begin < queryBegin) || (mhpIter->End > queryEnd)) continue;

		// expand the region
		if(mhpIter->Begin < mor.Begin) {
			mor.Begin          = mhpIter->Begin;
			mor.BeginOccupancy = mhpIter->Occupancy;
		}

		if(mhpIter->End > mor.End) {
			mor.End          = mhpIter->End;
			mor.EndOccupancy = mhpIter->Occupancy;
		}
	}

	// return the correction coefficient
	return mor.BeginOccupancy * mor.EndOccupancy;
}

// checks if the current alignment is a subset of another alignment
bool CNaiveAlignmentSet::CheckOverlap(const Alignment& al1, AlignmentSet::iterator& al2) {

	bool observedOverlap = true;

	// CASE 1: detect if al1 and al2 occur on different reference sequences
	if(al1.ReferenceIndex != al2->ReferenceIndex) return false;

	// CASE 2: detect if the coordinates are the same
	if((al1.ReferenceBegin == al2->ReferenceBegin) && (al1.ReferenceEnd == al2->ReferenceEnd) &&
		(al1.QueryBegin == al2->QueryBegin) && (al1.QueryEnd == al2->QueryEnd)) return true;

	// CASE 3: detect if al1 occurs before al2
	if(al1.ReferenceEnd < al2->ReferenceBegin) observedOverlap = false;

	// CASE 4: detect if al2 occurs before al1
	if(al2->ReferenceEnd < al1.ReferenceBegin) observedOverlap = false;

	// CASE 5: ignore overlap if the phase is skewed
	if(observedOverlap) {

		long long al1BeginDiff, al1EndDiff, al2BeginDiff, al2EndDiff;

		// calculate the phase difference for al1
		if(al1.IsReverseStrand) {
			al1BeginDiff = al1.ReferenceBegin + al1.QueryEnd;
			al1EndDiff   = al1.ReferenceEnd   + al1.QueryBegin;
		} else {
			al1BeginDiff = al1.ReferenceBegin - al1.QueryBegin;
			al1EndDiff   = al1.ReferenceEnd   - al1.QueryEnd;
		}

		// calculate the phase difference for al2
		if(al2->IsReverseStrand) {
			al2BeginDiff = al2->ReferenceBegin + al2->QueryEnd;
			al2EndDiff   = al2->ReferenceEnd   + al2->QueryBegin;
		} else {
			al2BeginDiff = al2->ReferenceBegin - al2->QueryBegin;
			al2EndDiff   = al2->ReferenceEnd   - al2->QueryEnd;
		}

		// return true if the phase differences are equal
		if((al1BeginDiff == al2BeginDiff) && (al1EndDiff == al2EndDiff)) return true;
	}

	// not overlapping or not in phase
	return false;
}

// resets the counter and stored alignments
void CNaiveAlignmentSet::Clear(void) {
	mAlignments.clear();
	mHasLongAlignment = false;
}

// dumps the contents of the alignment set to standard output
void CNaiveAlignmentSet::Dump(void) const {

	cout << "Naive alignment set contents" << endl;
	cout << "============================" << endl;

	unsigned int count = 1;
	for(AlignmentSet::const_iterator setIter = mAlignments.begin(); setIter != mAlignments.end(); setIter++, count++) {
		printf("%3u, ref index: %3u, begin: %9u, end: %9u, query begin: %2u, end: %2u, orientation: %c, quality: %2u, mm: %u\n", 
			count, setIter->ReferenceIndex, setIter->ReferenceBegin, setIter->ReferenceEnd, setIter->QueryBegin, setIter->QueryEnd,
			(setIter->IsReverseStrand ? 'R' : 'F'), setIter->Quality, setIter->NumMismatches);
		printf("   %s\n", setIter->Reference.CData());
		printf("   %s\n", setIter->Query.CData());
	}

	cout << endl << "alignment set count: " << mAlignments.size() << endl;
}

// returns the number of alignments in the set
unsigned int CNaiveAlignmentSet::GetCount(void) const {
	return (unsigned int)mAlignments.size();
}

// retrieves the mhp occupancy list for the forward read
MhpOccupancyList* CNaiveAlignmentSet::GetFwdMhpOccupancyList(void) {
	return &mFwdMhpOccupancyList;
}

// retrieves the mhp occupancy list for the reverse read
MhpOccupancyList* CNaiveAlignmentSet::GetRevMhpOccupancyList(void) {
	return &mRevMhpOccupancyList;
}

// retrieves the alignment set
AlignmentSet* CNaiveAlignmentSet::GetSet(void) {
	return &mAlignments;
}

// retrieves the long alignment flag
bool CNaiveAlignmentSet::HasLongAlignment(void) const {
	return mHasLongAlignment;
}

// returns true if the alignment set is empty
bool CNaiveAlignmentSet::IsEmpty(void) const {
	return mAlignments.empty();
}

// returns true if the alignment set contains only one entry
bool CNaiveAlignmentSet::IsUnique(void) const {
	return (mAlignments.size() == 1 ? true : false);
}
