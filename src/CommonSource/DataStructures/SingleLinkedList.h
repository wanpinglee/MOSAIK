// ***************************************************************************
// CSingleLinkedList - a fast and lightweight single-linked list.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iostream>

template <class V>
struct LinkedListNode {
	V Value;
	LinkedListNode* pNext;

	// constructor
	LinkedListNode()
		: pNext(NULL)
	{}

	// data constructor
	LinkedListNode(const V& v) 
		: Value(v)
		, pNext(NULL)
	{}
};

template <class V>
class CSingleLinkedList {
public:
	// constructor
	CSingleLinkedList(void);
	// destructor
	~CSingleLinkedList(void);
	// deletes the head node
	void DeleteHead(void);
	// dumps the contents of the list
	void Dump(void);
	// returns the head node
	LinkedListNode<V>* GetHead(void);
	// returns the size of the list
	unsigned int GetSize(void);
	// returns the tail node
	LinkedListNode<V>* GetTail(void);
	// add a node to the back
	void InsertBack(const V& v);
	// add a node to the front
	void InsertFront(const V& v);

private:
	// our head node
	LinkedListNode<V>* mpHead;
	// our tail node
	LinkedListNode<V>* mpTail;
	// our current size
	unsigned int mSize;
};

// constructor
template <class V>
CSingleLinkedList<V>::CSingleLinkedList(void)
: mpHead(NULL)
, mpTail(NULL)
, mSize(0)
{}

// destructor
template <class V>
CSingleLinkedList<V>::~CSingleLinkedList(void) {
	LinkedListNode<V>* pNode = mpHead;

	while(pNode) {
		LinkedListNode<V>* pNextNode = pNode->pNext;
		delete pNode;
		pNode = pNextNode;
	}
}

// deletes the head node
template <class V>
void CSingleLinkedList<V>::DeleteHead(void) {

	// handle empty lists
	if(!mpHead) return;

	// handle lists with one element
	if(mSize == 1) {
		delete mpHead;
		mpHead = NULL;
		mpTail = NULL;
		mSize = 0;
		return;
	}

	LinkedListNode<V>* pNext = mpHead->pNext;
	delete mpHead;	
	mpHead = pNext;
	mSize--;
}

// dumps the contents of the list
template <class V>
void CSingleLinkedList<V>::Dump(void) {

	std::cout << "list contains (" << mSize << "): ";

	LinkedListNode<V>* pNode = mpHead;
	while(pNode != NULL) {
		std::cout << "[" << pNode->Value << "] ";
		pNode = pNode->pNext;
	}

	std::cout << std::endl;
}

// returns the key for the tail element
template <class V>
LinkedListNode<V>* CSingleLinkedList<V>::GetHead(void) {
	return mpHead;
}

// returns the size of the list
template <class V>
unsigned int CSingleLinkedList<V>::GetSize(void) {
	return mSize;
}

// returns the key for the tail element
template <class V>
LinkedListNode<V>* CSingleLinkedList<V>::GetTail(void) {
	return mpTail;
}

// add a node to the front
template <class V>
void CSingleLinkedList<V>::InsertFront(const V& v) {

	LinkedListNode<V>* pNode = new LinkedListNode<V>(v);

	if(!mpHead) {
		mpHead = pNode;
		mpTail = pNode;
	} else {
		pNode->pNext = mpHead;
		mpHead = pNode;
	}

	mSize++;
}

// add a node to the back
template <class V>
void CSingleLinkedList<V>::InsertBack(const V& v) {

	LinkedListNode<V>* pNode = new LinkedListNode<V>(v);

	if(!mpHead) {
		mpHead = pNode;
		mpTail = pNode;
	} else {
		mpTail->pNext = pNode;
		mpTail = pNode;
	}

	mSize++;
}
