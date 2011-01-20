
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
		~CAlignedReadCache();
		// add a new aligned read in _cache
		bool Add ( const Mosaik::AlignedRead& ar );
		// sort _cache by the position of the first mate in Mate1Alignments
		void SortByPosition ( void );
		// sort _cache by read names
		void SortByName ( void );
		// store _cache in the given filename by CAlignmentWriter
		//bool StoreCacheInFile ( const string& filename, const vector<ReferenceSequence>* pReferenceSequences, const vector<MosaikReadFormat::ReadGroup>& readGrou, const AlignmentStatus& as );
		// reset _cache
		bool Reset ( void );
		// return _full
		bool isFull ( void ) { return _full; };
		// return
		bool isEmpty ( void );
		// clear _cache
		void Clear( void );
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
		static inline bool SortByReadName( const Mosaik::AlignedRead& ar1, const Mosaik::AlignedRead& ar2 );
};

#endif 
