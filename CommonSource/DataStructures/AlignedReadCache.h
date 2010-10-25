
#ifndef _AlignedReadCache_H_
#define _AlignedReadCache_H_

#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include "AlignedRead.h"
#include "Alignment.h"
#include "MosaikString.h"

using namespace std;

class CAlignedReadCache {
	public:

		CAlignedReadCache( unsigned int cacheSize );
		// add a new aligned read in _cache
		bool Add ( const Mosaik::AlignedRead& ar );
		// sort _cache by the position of the first mate in Mate1Alignments
		void SortByPosition ( void );
		// store _cache in the given filename by CAlignmentWriter
		//bool StoreCacheInFile ( const string& filename, const vector<ReferenceSequence>* pReferenceSequences, const vector<MosaikReadFormat::ReadGroup>& readGrou, const AlignmentStatus& as );
		// reset _cache
		bool Reset ( void );
		// return _full
		bool isFull ( void ) { return _full; };
		// load the next aligned read
		bool LoadNextAlignedRead( Mosaik::AlignedRead& ar );
		// rewind the load pointer
		void Rewind ( void ) { 
			_loadIte = _cache.begin();
			_loadNo = 0;
		};
		

	private:
		list<Mosaik::AlignedRead> _cache;
		list<Mosaik::AlignedRead>::iterator _currentIte;
		list<Mosaik::AlignedRead>::iterator _loadIte;
		unsigned int _cacheSize;
		unsigned int _currentNo;
		unsigned int _loadNo;
		bool _full;

		
		static inline bool SortBy1stMatePosition( const Mosaik::AlignedRead& ar1, const Mosaik::AlignedRead& ar2 );
		/*
		struct Sort {
			bool operator()( const Mosaik::AlignedRead ar1, const Mosaik::AlignedRead ar2 ) {
				unsigned int nMate1Ar1 = ar1.Mate1Alignments.size();
				unsigned int nMate1Ar2 = ar2.Mate1Alignments.size();

				if ( nMate1Ar1 == 0 ) return true;
				if ( nMate1Ar2 == 0 ) return false;
				
				unsigned int refIndexAr1 = ar1.Mate1Alignments.begin()->ReferenceIndex;
				unsigned int refIndexAr2 = ar2.Mate1Alignments.begin()->ReferenceIndex;
				unsigned int refBeginAr1 = ar1.Mate1Alignments.begin()->ReferenceBegin;
				unsigned int refBeginAr2 = ar2.Mate1Alignments.begin()->ReferenceBegin;
	
				if ( refIndexAr1 == refIndexAr2 )
					return refBeginAr1 < refBeginAr2;
	
				return refIndexAr1 < refIndexAr2;
			}
		};
		*/
};

#endif 
