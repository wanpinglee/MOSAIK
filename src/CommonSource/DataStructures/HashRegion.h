// ***************************************************************************
// HashRegion - defines where a hash matches in relation to both the read and
//              the reference sequence.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

struct HashRegion {
	unsigned int Begin;
	unsigned int End;
	unsigned short QueryBegin;
	unsigned short QueryEnd;
	unsigned short NumMismatches;

	HashRegion()
		: Begin(0)
		, End(0)
		, QueryBegin(0)
		, QueryEnd(0)
		, NumMismatches(0)
	{}

	bool operator<(const HashRegion& r) const {
		if(Begin      != r.Begin)      return Begin      < r.Begin;
		//if(End        != r.End)        return End        < r.End;
		if(QueryBegin != r.QueryBegin) return QueryBegin < r.QueryBegin;		
		if(End        != r.End)        return End        < r.End;
		if(QueryEnd   != r.QueryEnd)   return QueryEnd   < r.QueryEnd;
		return false;
	}
};
