// ***************************************************************************
// CMruCache - keeps the most recently used hash/genome position pairs.
//             Essentially a LRU caching algorithm.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iostream>
#include "DoubleLinkedList.h"
#include "Mosaik.h"
#include "PosixThreads.h"
#include "UnorderedMap.h"

using namespace std;

template<class K, class V>
class CMruCache {
public:
	// constructor
	CMruCache(const unsigned int maxSize);
	// destructor
	~CMruCache(void);
	// gets the element from the hash map
	bool Get(const K& key, V& val);
	// retrieves the cache statistics
	void GetStatistics(uint64_t& cacheHits, uint64_t& cacheMisses) const;
	// our insert operator
	void Insert(const K& key, const V& val);
	// our cache mutex
	static pthread_mutex_t mMruCacheMutex;

private:
	// our most recently used list
	CDoubleLinkedList<K,V> mMruList;
	// our most recently used map
	unordered_map<K, LinkedListNode<K,V>*> mMruMap;
	// our maximum size
	unsigned int mMaxSize;
	// our statistical counters
	uint64_t mCacheHits;
	uint64_t mCacheMisses;
};

// register our mutex
template<class K, class V>
pthread_mutex_t CMruCache<K,V>::mMruCacheMutex;

// constructor
template<class K, class V>
CMruCache<K,V>::CMruCache(const unsigned int maxSize)
: mMaxSize(maxSize)
, mCacheHits(0)
, mCacheMisses(0)
{}

// destructor
template<class K, class V>
CMruCache<K,V>::~CMruCache(void) {}

// gets the element from the hash map
template<class K, class V>
bool CMruCache<K,V>::Get(const K& key, V& val) {

	bool ret = true;
	pthread_mutex_lock(&mMruCacheMutex);

	// search for our element
	typename unordered_map<K, LinkedListNode<K,V>*>::iterator hashIter = mMruMap.find(key);

	if(hashIter != mMruMap.end()) {
		LinkedListNode<K,V>* pNode = hashIter->second; 
		val = pNode->Value;
		mMruList.MoveToHead(pNode);
		mCacheHits++;
	} else {
		ret = false;
		mCacheMisses++;
	}

	pthread_mutex_unlock(&mMruCacheMutex);
	return ret;
}

// retrieves the cache statistics
template<class K, class V>
void CMruCache<K,V>::GetStatistics(uint64_t& cacheHits, uint64_t& cacheMisses) const {
	cacheHits   = mCacheHits;
	cacheMisses = mCacheMisses;
}

// our insert operator
template<class K, class V>
void CMruCache<K,V>::Insert(const K& key, const V& val) {

	pthread_mutex_lock(&mMruCacheMutex);
	typename unordered_map<K, LinkedListNode<K,V>*>::iterator hashIter;

	// remove the least recently used value
	if(mMruList.GetSize() == mMaxSize) {

		// delete the tail key from the hash map
		K tailKey = mMruList.GetTailKey();
		hashIter = mMruMap.find(tailKey);

		// show an error message if we can't find the element
		if(hashIter == mMruMap.end()) {
			cout << "ERROR: Could not find the tail element (" << tailKey << ") in the hash map." << endl;
			exit(1);
		}

		mMruMap.erase(hashIter);
		mMruList.DeleteTail();
	}

	hashIter = mMruMap.find(key);

	if(hashIter == mMruMap.end()) {
		LinkedListNode<K,V>* pNode = mMruList.Insert(key, val);
		mMruMap[key] = pNode;
	}

	pthread_mutex_unlock(&mMruCacheMutex);
}
