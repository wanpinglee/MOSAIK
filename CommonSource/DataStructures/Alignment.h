// ***************************************************************************
// Alignment.h - stores everything related to alignments.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <string>

#include "MosaikString.h"

using namespace std;

// our alignment structure [all members identified as (temp) are not serialized to disk]
struct Alignment {

	unsigned int MateReferenceBegin;   // required for SAM/BAM
	unsigned int MateReferenceEnd;     // required for SAM/BAM
	unsigned int MateReferenceIndex;   // required for SAM/BAM
	unsigned int ReferenceBegin;
	unsigned int ReferenceEnd;
	unsigned int ReferenceIndex;
	unsigned int Owner;                // the temporary file that contains the alignment
	unsigned int ReadGroupCode;        // the read group code (temp)
	unsigned int NumMapped;            // the total number of mapped alignments
	int FragmentLength;                // the fragment length with its pair
	unsigned short QueryLength;        // used during filtering (temp)
	unsigned short NumMismatches;      // number of mismatches
	unsigned short QueryBegin;
	unsigned short QueryEnd;
	unsigned char Quality;             // alignment quality
	unsigned char NextBestQuality;     // the next best alignment quality
	bool CanBeMappedToSpecialReference;// can the sequence be mapped to special references?`
	bool IsFirstMate;                  // is this alignment from the first mate in a paired-end read
	bool IsJunk;                       // are the fileds in this alignment used for other propose, e.g. counting total numbers of alignments?
	bool IsMateReverseStrand;          // read orientation for the mate
	bool IsPairedEnd;                  // is the read sequenced as a paired-end read
	bool IsResolvedAsPair;             // is the alignment part of resolved paired-end read
	bool IsResolvedAsProperPair;       // is the alignment resolved as proper pair
	bool IsReverseStrand;              // read orientation
	bool IsMappedSpecialReference;     // is this alignment mapped to the special references which is defined by "-sref"? 
	bool IsMapped;                     // is this alignment mapped?
	bool IsMateMapped;                 // is its mate mapped?
	bool WasRescued;                   // was the alignment rescued during local alignment search
	char* ReferenceName;               // only filled via CAlignmentReader (temp)
	CMosaikString Reference;
	CMosaikString Query;
	CMosaikString BaseQualities;
	CMosaikString Name;                // the read name
	string Cigar;
	string ReadGroup;                  // the read group string
	string SpecialCode;                // 2 letters to indicate the sequence can be mapped to which special reference
	bool Mark;

	// constructors
	Alignment(void)
		: MateReferenceBegin(0)
		, MateReferenceEnd(0)
		, MateReferenceIndex(0)
		, ReferenceBegin(0)
		, ReferenceEnd(0)
		, ReferenceIndex(0)
		, ReadGroupCode(0)
		, NumMapped(1)
		, FragmentLength(0)
		, QueryBegin(0)
		, QueryEnd(0)
		, Quality(0)
		, NextBestQuality(0)
		, CanBeMappedToSpecialReference(false)
		, IsFirstMate(false)
		, IsJunk(false)
		, IsMateReverseStrand(false)
		, IsPairedEnd(false)
		, IsResolvedAsPair(false)
		, IsReverseStrand(false)
		, IsMappedSpecialReference(false)
		, IsMapped(true)
		, IsMateMapped(false)
		, WasRescued(false)
		, ReferenceName(NULL)
		, Mark(false)
	{}

	// our less-than operator
	bool operator<(const Alignment& al) const {
		if(ReferenceIndex == al.ReferenceIndex) return ReferenceBegin < al.ReferenceBegin;
		return ReferenceIndex < al.ReferenceIndex;
	}

	
	bool SetPairFlags ( const Alignment& pairMate, const int& allowedFragmentLength, const bool& expectedMateStrand ) {
		unsigned int queryPosition5Prime = ( IsReverseStrand ) ? ReferenceEnd : ReferenceBegin;
		unsigned int matePosition5Prime  = ( pairMate.IsReverseStrand ) ? pairMate.ReferenceEnd : pairMate.ReferenceBegin;
		FragmentLength = ( ReferenceIndex != pairMate.ReferenceIndex ) ? 0 : matePosition5Prime - queryPosition5Prime;
		
		if ( expectedMateStrand != pairMate.IsReverseStrand )
			IsResolvedAsProperPair = false;
		else if ( ReferenceIndex != pairMate.ReferenceIndex )
			IsResolvedAsProperPair = false;
		else {
			if ( ( allowedFragmentLength >= 0 ) && ( FragmentLength <= allowedFragmentLength ) )
				IsResolvedAsProperPair = true;
			else if ( ( allowedFragmentLength < 0 ) && ( FragmentLength >= allowedFragmentLength ) )
				IsResolvedAsProperPair = true;
			else
				IsResolvedAsProperPair = false;
		}

		//IsMateReverseStrand = pairMate.IsReverseStrand;
		//MateReferenceIndex = pairMate.ReferenceIndex;
		//MateReferenceBegin = pairMate.ReferenceBegin;

		return IsResolvedAsProperPair;
	}
};
