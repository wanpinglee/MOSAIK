/*
 * =====================================================================================
 *
 *       Filename:  SortThread.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/25/2010 03:19:05 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "PosixThreads.h"
#include "ArchiveSort.h"
#include "ProgressCounter.h"
#include "ProgressBar.h"
#include "Alignment.h"

#ifndef SortThread_H
#define SortThread_H

using namespace std;

namespace SortThreadData {
  

//class SortThread {

//	private: 

	struct ThreadData {
		vector < string > _input;
		vector < string > _output;
		unsigned int _nArchive;
		unsigned int _archiveNo;
		uint64_t     _nRead;
		uint64_t     _readNo;
		unsigned int _medianFragmentLength;
		unsigned int _sorterCacheSize;

		pthread_attr_t  _attr;
		pthread_mutex_t _mutex;
		pthread_mutex_t _readCounterMutex;

		ThreadData()
			: _input()
			, _output()
			, _nArchive(0)
			, _archiveNo(0)
			, _nRead(0)
			, _readNo(0)
			, _medianFragmentLength(0)
			, _sorterCacheSize(0)
			, _attr()
			, _mutex()
			, _readCounterMutex()
			{}
	};
/*
	void* StartThread( void* arg ){
		//SortThreadData::ThreadData *td = (SortThreadData::ThreadData *)arg;
		ThreadData *td = (ThreadData *)arg;

		string input;
		string output;
		//bool gatheringCounter = true;

		//CProgressBar<uint64_t>::StartThread(&readCounter, 0, td->nRead * 2, "reads");
		//CProgressCounter<unsigned int>::StartThread(&readCounter, &gatheringCounter, "aligned reads");

			pthread_mutex_lock(&td->_mutex);
		while( td->_archiveNo < td->_nArchive ) {
			input  = td->_input [td->_archiveNo];
			output = td->_output[td->_archiveNo];
			td->_archiveNo++;
			pthread_mutex_unlock(&td->_mutex);

			CArchiveSort sorter( input, output, &td->_readNo, &td->_readCounterMutex );
			sorter.Sort();
		}

		pthread_exit(NULL);

		//CProgressBar<uint64_t>::WaitThread();
		//CProgressCounter<unsigned int>::WaitThread();
	}
*/
}

class SortThread {
	
	public:
		SortThread( const vector< string > input, const vector< string > output, const unsigned int nThread, uint64_t nRead, const unsigned int medianFragmentLength, const unsigned int sorterCacheSize )
			: _td()
			, _nThread(nThread)
			, _threads(NULL)
			, IsQuietMode(false)
		{
			_td._input      = input;
			_td._output     = output;
			_td._nArchive   = input.size();
			_td._archiveNo  = 0;
			_td._nRead      = nRead;
			_td._medianFragmentLength = medianFragmentLength;
			_td._sorterCacheSize      = sorterCacheSize;

			pthread_attr_init (&_td._attr);
			pthread_attr_setdetachstate(&_td._attr, PTHREAD_CREATE_JOINABLE);
			pthread_mutex_init(&_td._mutex, NULL);
			pthread_mutex_init(&_td._readCounterMutex, NULL);
			_threads = new pthread_t[nThread];
		}

		~SortThread(){
			pthread_attr_destroy (&_td._attr);
			pthread_mutex_destroy(&_td._mutex);
			pthread_mutex_destroy(&_td._readCounterMutex);
		}

		void Start();
		void SetQuietMode( void );
	
	private:
		SortThreadData::ThreadData _td;
		//ThreadData _td;
		
		unsigned int _nThread;
		pthread_t *_threads;
		bool IsQuietMode;

		SortThread (const SortThread&);
		SortThread& operator=(const SortThread&);
		
};

/*  
void SortThread::Start() {
	
	
	//bool gatheringCounter = true;
	//CProgressCounter<unsigned int>::StartThread(&_td._readNo, &gatheringCounter, "aligned reads");
	CProgressBar<unsigned int>::StartThread(&_td._readNo, 0, _td._nRead * 2, "reads");

	for ( unsigned int i = 0 ; i < _nThread; i++ ) {
		int rc = pthread_create(&_threads[i], &_td._attr, (void*)&StartThread, (void*)&_td);
		if ( rc ) {
			cout << "ERROR: Return code from pthread_create is " << rc << endl;
			exit(1);
		}
	}

	void* status;
	for ( unsigned int i = 0 ; i < _nThread; i++ ) {
		int rc = pthread_join(_threads[i], &status);
		if ( rc ) {
			cout << "ERROR: Return code from pthread_join is " << rc << endl;
			exit(1);
		}
	}

	//CProgressCounter<unsigned int>::WaitThread();
	CProgressBar<unsigned int>::WaitThread();
}
*/
#endif
