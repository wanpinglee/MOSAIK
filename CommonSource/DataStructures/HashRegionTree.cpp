// ***************************************************************************
// CHashRegionTree - aggregates hash positions and forms alignment candidates.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "HashRegionTree.h"

namespace AVLTree {

	CHashRegionTree::CHashRegionTree(unsigned int queryLen, unsigned char hashSize)
	: mRoot(NULL)
	, mTraverse(NULL)
	//, mFirst(NULL)
	//, mLast(NULL)
	, mQueryLength(queryLen)
	, mHashSize(hashSize)
	, mIndelSize(3)
	, mCount(0)
	{}

	CHashRegionTree::~CHashRegionTree() {
		Clear();
	}

	// clears the tree
	void CHashRegionTree::Clear() {

		// nothing to traverse
		if(mRoot == NULL) return;

		HashRegionAvlNode* n = mRoot;
		HashRegionAvlNode* d = NULL;

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

	// overload the cout operator
	ostream& operator<<(ostream& os, const HashRegion& hr) {
		return os << "reference begin: " << hr.Begin << ", end: " << hr.End << 
			", query begin: " << hr.QueryBegin << ", end: " << hr.QueryEnd;
	}

	// dumps the tree in sorted order
	void CHashRegionTree::DumpTree() {

		// DEBUG
		unsigned int level = 0;

		// nothing to traverse
		if(mRoot == NULL) return;

		HashRegionAvlNode* n = mRoot;

		// find the minimum node
		while(n->Left != NULL) {
			n = n->Left;
			level++;
		}

		do {

			cout << n->Data << ", level: " << level << endl;
			if(n->Parent != NULL) {
				cout << "- parent: " << n->Parent->Data;
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

	// returns the current size of the tree
	unsigned int CHashRegionTree::GetCount() {
		return mCount;
	}

	// gets the current hash region at the traversal pointer
	HashRegion* CHashRegionTree::GetTraversalHashRegion() {
		return &mTraverse->Data;
	}

	// go to the first entry
	void CHashRegionTree::GotoFirstEntry() {

		// empty tree
		if(mRoot == NULL) return;

		// set our traverse pointer to the first element
		//mTraverse = mFirst;

		mTraverse = mRoot;
		while(mTraverse->Left) mTraverse = mTraverse->Left;
	}

	// go to the last entry
	void CHashRegionTree::GotoLastEntry() {

		// empty tree
		if(mRoot == NULL) return;

		// set our traverse pointer to the last element
		//mTraverse = mLast;

		mTraverse = mRoot;
		while(mTraverse->Right) mTraverse = mTraverse->Right;
	}

	// find the next entry
	bool CHashRegionTree::GetNextEntry() {
		
		// return if our traversal pointer is NULL
		if(!mTraverse) return false;

		//cout << "------------" << endl;
		//if(!mTraverse->Left)  cout << "LEFT NULL" << endl;
		//if(!mTraverse->Right) cout << "RIGHT NULL" << endl;
		//cout << "------------" << endl;

		// Next in order will be our right child's leftmost child (if null, our right child)
		if(mTraverse->Right) {

			mTraverse = mTraverse->Right;
			while(mTraverse->Left) mTraverse = mTraverse->Left;

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
	bool CHashRegionTree::GetPreviousEntry(HashRegion& key) {
		
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

	// find the previous entry
	void CHashRegionTree::MoveToPreviousEntry() {
		
		// return if our traversal pointer is NULL
		if(mTraverse == NULL) return;

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
	}

	// update the tree
	void CHashRegionTree::Insert(HashRegion& key) {

		// A1. Initialization
		HashRegionAvlNode* n = mRoot;
		HashRegionAvlNode* p = NULL;
		HashRegionAvlNode* q = NULL;

		// A2. Find insertion point
		if(n) {

			// look for a NULL pointer
			while(n) {

				p = n;
				if(p->Balance != 0) q = p;

				if(key < n->Data) n = n->Left;
					else if(n->Data < key) n = n->Right;
					else break;
			}
		}

		// A3. Insert
		mCount++;
		n = new HashRegionAvlNode(key, p);

		if(p) {
		
			// add a temporary link to the node, but do not re-balance
			bool isLeftChild = false;
			if(key < p->Data) {
				p->Left = n;
				isLeftChild = true;
			} else p->Right = n;

			// point our traversal pointer to the current node
			mTraverse = n;

			// move forward a little bit
			HashRegion* phr = NULL;

			// keep looking at previous entries as long as an expected
			// query length can exist between the current position and
			// the previous entry's start position

			bool foundCandidate = false;

			//
			// check previous entries
			//

			MoveToPreviousEntry();
			while(mTraverse) {

				// localize our previous hash region
				phr = &mTraverse->Data;

				// stop checking previous entries if an expected query length
				// cannot exist prior to the sequence
				if((phr->Begin + mQueryLength - 1) < key.Begin) break;

				// define our sequence differences
				int diffAnchors = key.Begin      - phr->End;
				int diffQueries = key.QueryBegin - phr->QueryEnd;

				// check if the sequence occurs within an island
				if((key.Begin >= phr->Begin) && (key.End <= phr->End) && (key.QueryBegin >= phr->QueryBegin) 
					&& (key.QueryEnd <= phr->QueryEnd)) {
					foundCandidate = true;
					break;
				}

				// in phase
				if((!foundCandidate) && (diffAnchors == diffQueries) && (diffQueries >= -mHashSize) 
					&& (diffQueries <= mHashSize)) {
					foundCandidate = true;
				}

				// insertion
				if((!foundCandidate) && (diffAnchors == 1) && (diffQueries > 1) && (diffQueries <= (mIndelSize + 1))) {				
					foundCandidate = true;
					//phr->HasIndel  = true;
				}

				// deletion
				if((!foundCandidate) && (diffQueries == 1) && (diffAnchors > 1) && (diffAnchors <= (mIndelSize + 1))) {				
					foundCandidate = true;
					//phr->HasIndel  = true;
				}

				// adjust the endpoints if the result of an extension or indel
				if(foundCandidate) {
					phr->End      = key.End;
					phr->QueryEnd = key.QueryEnd;
					break;
				}

				MoveToPreviousEntry();
			}

			// update our first and last pointers
			//if(mFirst) {
			//	if(key < mFirst->Data) mFirst = n;
			//} else mFirst = n;

			//if(mLast) {
			//	if(mLast->Data < key) mLast = n;
			//} else mLast = n;

			// delete the previous node if we found a candidate
			if(foundCandidate) {

				// reset the parent's links
				if(isLeftChild) p->Left = NULL;
					else p->Right = NULL;

				// delete the candidate node
				delete n;

				// decrement the count
				mCount--;

			} else {

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
			}

		} else mRoot = n;
	}

	// rotates a given node left
	void CHashRegionTree::LeftRotate(HashRegionAvlNode* n) {
		
		// L1. Do the rotation
		HashRegionAvlNode* r = n->Right;
		n->Right = r->Left;
		
		if(r->Left != NULL) r->Left->Parent = n;
		
		HashRegionAvlNode* p = n->Parent;
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
	void CHashRegionTree::RightRotate(HashRegionAvlNode* n) {
		
		// R1. Do the rotation
		HashRegionAvlNode* l = n->Left;
		n->Left = l->Right;
		
		if(l->Right != NULL) l->Right->Parent = n;
		
		HashRegionAvlNode* p = n->Parent;
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

	// sets the expected query length
	void CHashRegionTree::SetExpectedQueryLength(unsigned int queryLen) {
		mQueryLength = queryLen;
	}
}
