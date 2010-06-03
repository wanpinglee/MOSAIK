// ***************************************************************************
// CAvlTree - an implementation of an AVL tree (G.M. Adelson-Velskii and 
//            E.M. Landis). Faster than a red-black tree (lower constant).
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iostream>

using namespace std;

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

namespace AVLTree {

	template<class K>
	struct AvlNode {
		K Data;
		AvlNode* Left;
		AvlNode* Right;
		AvlNode* Parent;
		int Balance;

		AvlNode(K& data, AvlNode* parent)
		: Data(data)
		, Left(NULL)
		, Right(NULL)
		, Parent(parent)
		, Balance(0)
		{}
	};

	template<class K>
	class CAvlTree {
	public:
		CAvlTree();
		~CAvlTree();
		// the current element count
		unsigned int mCount;
		// clears the tree
		void Clear();
		// dumps the tree in sorted order
		void DumpTree();
		// finds the key in the tree
		bool Find(K& key);
		// go to the first entry
		void GotoFirstEntry();
		// go to the last entry
		void GotoLastEntry();
		// moves on to the next element in the tree
		bool GetNextEntry(K& key);
		// moves on to the previous element in the tree
		bool GetPreviousEntry(K& key);
		// insert into the tree; duplicates are ignored
		void Insert(K& key);
	private:
		// our tree root
		AvlNode<K>* mRoot;
		// our traversal pointer root
		AvlNode<K>* mTraverse;
		// our first entry
		AvlNode<K>* mFirst;
		// our last entry
		AvlNode<K>* mLast;
		// rotates a given node left
		void LeftRotate(AvlNode<K>* n);
		// rotates a given node right
		void RightRotate(AvlNode<K>* n);
	};

	template<class K>
	CAvlTree<K>::CAvlTree()
	: mRoot(NULL)
	, mTraverse(NULL)
	, mFirst(NULL)
	, mLast(NULL)
	, mCount(0)
	{}

	template<class K>
	CAvlTree<K>::~CAvlTree() {
		Clear();
	}

	// clears the tree
	template<class K>
	void CAvlTree<K>::Clear() {

		// nothing to traverse
		if(mRoot == NULL) return;

		AvlNode<K>* n = mRoot;
		AvlNode<K>* d = NULL;

		// find the minimum leftmost level
		while(n->Left  != NULL) n = n->Left;
		while(n->Right != NULL) n = n->Right;

		do {

			// create a pointer to our delete node
			d = n;

			if(n->Parent == NULL) {
				delete d;
				break;

			// Next in order will be our parent's right child's leftmost child
			} else if((n->Parent->Right != NULL) && (n->Parent->Right != n)) {

				n = n->Parent->Right;
				while(n->Left  != NULL) n = n->Left;
				while(n->Right != NULL) n = n->Right;
				delete d;
			
			} else {
				n = n->Parent;
				delete d;
			}

		} while(n != NULL);

		// reset the counter and root pointer
		mRoot = NULL;
		mCount = 0;
	}

	// dumps the tree in sorted order
	template<class K>
	void CAvlTree<K>::DumpTree() {

		// DEBUG
		unsigned int level = 0;

		// nothing to traverse
		if(mRoot == NULL) return;

		AvlNode<K>* n = mRoot;

		// find the minimum node
		while(n->Left != NULL) {
			n = n->Left;
			level++;
		}

		do {

			cout << n->Data << ", level: " << level;
			if(n->Parent != NULL) {
				cout << ", parent: " << n->Parent->Data;
				if(n->Parent->Left == n) cout << " LEFT";
					else cout << " RIGHT";
			}
			cout << endl;

			// Next in order will be our right child's leftmost child (if null, our right child)
			if(n->Right != NULL) {

				n = n->Right;
				level++;
				while(n->Left != NULL) {
					n = n->Left;
					level++;
				}

			} else {

				// Next in order will be our closest ancestor that hasn't been visited yet
				while(true) {
				
					if(n->Parent == NULL) {
						n = NULL;
						break;
					}

					bool stop = false;
					if(n->Parent->Left == n) stop = true;

					n = n->Parent;
					level--;

					// If node is its parents left child, then its parent hasn't been visited yet
					if(stop) break;
				}
			}

		} while(n != NULL);
	}

	// finds the key in the tree
	template<class K>
	bool CAvlTree<K>::Find(K& key) {

		// nothing to traverse
		if(mRoot == NULL) return false;

		AvlNode<K>* n = mRoot;

		//unsigned int numComparisons = 1;

		// look for a NULL pointer
		while(n != NULL) {

			//cout << "# comparisons: " << numComparisons << endl;

			if(key < n->Data) n = n->Left;
			else if(key > n->Data) n = n->Right;
			else return true;

			//numComparisons++;
		}

		return false;
	}

	// go to the first entry
	template<class K>
	void CAvlTree<K>::GotoFirstEntry() {

		// empty tree
		if(mRoot == NULL) return;

		// set our traverse pointer to the first element
		mTraverse = mFirst;
	}

	// go to the last entry
	template<class K>
	void CAvlTree<K>::GotoLastEntry() {

		// empty tree
		if(mRoot == NULL) return;

		// set our traverse pointer to the last element
		mTraverse = mLast;
	}

	// find the next entry
	template<class K>
	bool CAvlTree<K>::GetNextEntry(K& key) {
		
		// return if our traversal pointer is NULL
		if(mTraverse == NULL) return false;

		// retrieve the key
		key = mTraverse->Data;

		// Next in order will be our right child's leftmost child (if null, our right child)
		if(mTraverse->Right != NULL) {

			mTraverse = mTraverse->Right;
			while(mTraverse->Left != NULL) mTraverse = mTraverse->Left;

		} else {

			// Next in order will be our closest ancestor that hasn't been visited yet
			while(true) {

				if(mTraverse->Parent == NULL) {
					mTraverse = NULL;
					break;
				}

				bool stop = false;
				if(mTraverse->Parent->Left == mTraverse) stop = true;

				mTraverse = mTraverse->Parent;

				// If node is its parents left child, then its parent hasn't been visited yet
				if(stop) break;
			}
		}

		return true;
	}

	// find the previous entry
	template<class K>
	bool CAvlTree<K>::GetPreviousEntry(K& key) {
		
		// return if our traversal pointer is NULL
		if(mTraverse == NULL) return false;

		// retrieve the key
		key = mTraverse->Data;

		// Next in order will be our left child's rightmost child (if null, our left child)
		if(mTraverse->Left != NULL) {

			mTraverse = mTraverse->Left;
			while(mTraverse->Right != NULL) mTraverse = mTraverse->Right;

		} else {

			// Next in order will be our closest ancestor that hasn't been visited yet
			while(true) {

				if(mTraverse->Parent == NULL) {
					mTraverse = NULL;
					break;
				}

				bool stop = false;
				if(mTraverse->Parent->Right == mTraverse) stop = true;

				mTraverse = mTraverse->Parent;

				// If node is its parents right child, then its parent hasn't been visited yet
				if(stop) break;
			}
		}

		return true;
	}

	// insert into the tree; duplicates are ignored
	template<class K>
	void CAvlTree<K>::Insert(K& key) {

		// A1. Initialization
		AvlNode<K>* n = mRoot;
		AvlNode<K>* p = NULL;
		AvlNode<K>* q = NULL;

		// A2. Find insertion point
		if(n != NULL) {

			// look for a NULL pointer
			while(n != NULL) {

				p = n;
				if(p->Balance != 0) q = p;

				if(key < n->Data) n = n->Left;
					else if(key > n->Data) n = n->Right;
					else return;
			}
		}

		// A3. Insert
		mCount++;
		n = new AvlNode<K>(key, p);

		// update our first and last pointers
		if(mFirst) {
			if(key < mFirst->Data) mFirst = n;
		} else mFirst = n;

		if(mLast) {
			if(key > mLast->Data) mLast = n;
		} else mLast = n;

		if(p != NULL) {
		
			if(key < p->Data) p->Left = n;
				else p->Right = n;

			// A4. Adjust balance factors
			while(p != q) {
				if(p->Left == n) p->Balance = -1;
					else p->Balance = 1;

				n = p;
				p = p->Parent;
			}

			// A5. Check for imbalance
			if(q != NULL) {

				if(q->Left == n) {
					--q->Balance;
				
					if(q->Balance == -2) {

						// A6. Left imbalance
						if(q->Left->Balance > 0) LeftRotate(q->Left);
						RightRotate(q);
					}
				}

				if(q->Right == n) {
					++q->Balance;

					if(q->Balance == 2) {
						
						// A7. Right imbalance
						if(q->Right->Balance < 0) RightRotate(q->Right);
						LeftRotate(q);	
					}
				}
			}

		} else mRoot = n; 

		// A8 All done.
	}

	// rotates a given node left
	template<class K>
	void CAvlTree<K>::LeftRotate(AvlNode<K>* n) {
		
		// L1. Do the rotation
		AvlNode<K>* r = n->Right;
		n->Right = r->Left;
		
		if(r->Left != NULL) r->Left->Parent = n;
		
		AvlNode<K>* p = n->Parent;
		r->Parent = p;

		if(p != NULL) {
			if(p->Left == n) p->Left = r;
			else p->Right = r;
		} else mRoot = r;

		r->Left = n;
		n->Parent = r;

		// L2. Recompute the balance factors
		n->Balance = n->Balance - (1 + MAX(r->Balance, 0));
		r->Balance = r->Balance - (1 - MIN(n->Balance, 0));
	}

	// rotates a given node right
	template<class K>
	void CAvlTree<K>::RightRotate(AvlNode<K>* n) {
		
		// R1. Do the rotation
		AvlNode<K>* l = n->Left;
		n->Left = l->Right;
		
		if(l->Right != NULL) l->Right->Parent = n;
		
		AvlNode<K>* p = n->Parent;
		l->Parent = p;

		if(p != NULL) {
			if(p->Left == n) p->Left = l;
			else p->Right = l;
		} else mRoot = l;

		l->Right = n;
		n->Parent = l;

		// R2. Recompute the balance factors
		n->Balance = n->Balance + (1 - MIN(l->Balance, 0));
		l->Balance = l->Balance + (1 + MAX(n->Balance, 0));
	}
}
