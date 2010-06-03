// ***************************************************************************
// CDoubleLinkedList - a fast and lightweight double-linked list.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iostream>

template <class K, class V>
struct LinkedListNode {
	K Key;
	V Value;
	LinkedListNode* pNext;
	LinkedListNode* pPrev;

	// constructor
	LinkedListNode()
		: pNext(NULL)
		, pPrev(NULL)
	{}

	// data constructor
	LinkedListNode(const K& k, const V& v) 
		: Key(k)
		, Value(v)
		, pNext(NULL)
		, pPrev(NULL)
	{}
};

template <class K, class V>
class CDoubleLinkedList {
public:
	// constructor
	CDoubleLinkedList(void);
	// destructor
	~CDoubleLinkedList(void);
	// deletes the tail node
	void DeleteTail(void);
	// dumps the contents of the list
	void Dump(void);
	// returns the size of the list
	unsigned int GetSize(void);
	// returns the key for the tail element
	K& GetTailKey(void);
	// add a node to the front
	LinkedListNode<K,V>* Insert(const K& k, const V& v);
	// moves the node to the front
	void MoveToHead(LinkedListNode<K,V>* pNode);

private:
	// our head node
	LinkedListNode<K,V>* mpHead;
	// our tail node
	LinkedListNode<K,V>* mpTail;
	// our current size
	unsigned int mSize;
};

// constructor
template <class K, class V>
CDoubleLinkedList<K,V>::CDoubleLinkedList(void)
: mpHead(NULL)
, mpTail(NULL)
, mSize(0)
{}

// destructor
template <class K, class V>
CDoubleLinkedList<K,V>::~CDoubleLinkedList(void) {
	LinkedListNode<K,V>* pNode = mpHead;

	while(pNode) {
		LinkedListNode<K,V>* pNextNode = pNode->pNext;
		delete pNode;
		pNode = pNextNode;
	}
}

// deletes the tail node
template <class K, class V>
void CDoubleLinkedList<K,V>::DeleteTail(void) {

	// handle empty lists
	if(!mpTail) return;

	// handle lists with one element
	if(mSize == 1) {
		delete mpHead;
		mpHead = NULL;
		mpTail = NULL;
		mSize = 0;
		return;
	}

	LinkedListNode<K,V>* pNewTail = mpTail->pPrev;
	pNewTail->pNext = NULL;
	delete mpTail;
	mpTail = pNewTail;
	mSize--;
}


// dumps the contents of the list
template <class K, class V>
void CDoubleLinkedList<K,V>::Dump(void) {

	std::cout << "list contains (" << mSize << "): ";

	LinkedListNode<K,V>* pNode = mpHead;
	while(pNode != NULL) {
		std::cout << "[" << pNode->Key << "] ";
		pNode = pNode->pNext;
	}

	std::cout << std::endl;
}

// returns the size of the list
template <class K, class V>
unsigned int CDoubleLinkedList<K,V>::GetSize(void) {
	return mSize;
}

// returns the key for the tail element
template <class K, class V>
K& CDoubleLinkedList<K,V>::GetTailKey(void) {
	return mpTail->Key;
}

// add a node to the front
template <class K, class V>
LinkedListNode<K,V>* CDoubleLinkedList<K,V>::Insert(const K& k, const V& v) {

	LinkedListNode<K,V>* pNode = new LinkedListNode<K,V>(k,v);

	if(!mpHead) {
		mpHead = pNode;
		mpTail = pNode;
	} else {
		mpHead->pPrev = pNode;
		pNode->pNext = mpHead;
		mpHead = pNode;
	}

	mSize++;

	return mpHead;
}

// moves the node to the front
template <class K, class V>
void CDoubleLinkedList<K,V>::MoveToHead(LinkedListNode<K,V>* pNode) {

	// handle null nodes, empty lists, and lists with one element
	if(!pNode || !mpHead || (mSize == 1) || (pNode == mpHead)) return;

	// restitch the extraction point
	LinkedListNode<K,V>* pPrev = pNode->pPrev;
	LinkedListNode<K,V>* pNext = pNode->pNext;

	pPrev->pNext = pNext;
	if(pNext) pNext->pPrev = pPrev;
	else mpTail = pPrev;

	// move the node to the front
	pNode->pPrev  = NULL;
	pNode->pNext  = mpHead;
	mpHead->pPrev = pNode;
	mpHead        = pNode;
}
