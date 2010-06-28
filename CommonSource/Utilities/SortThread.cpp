#include "SortThread.h"


	
	void* StartThread( void* arg ){
		SortThreadData::ThreadData *td = (SortThreadData::ThreadData *)arg;
		//ThreadData *td = (ThreadData *)arg;

		string input;
		string output;


			pthread_mutex_lock(&td->_mutex);
		
		while( td->_archiveNo < td->_nArchive ) {
			input  = td->_input [td->_archiveNo];
			output = td->_output[td->_archiveNo];
			td->_archiveNo++;

			pthread_mutex_unlock(&td->_mutex);

			CArchiveSort sorter( input, output, &td->_readNo, &td->_readCounterMutex, td->_medianFragmentLength );
			sorter.Sort();
		}

		pthread_exit(NULL);

	}
//}

//class SortThread {


void SortThread::Start() {
	
	
	CProgressBar<unsigned int>::StartThread(&_td._readNo, 0, _td._nRead * 2, "reads");

	for ( unsigned int i = 0 ; i < _nThread; i++ ) {
		int rc = pthread_create(&_threads[i], &_td._attr, StartThread, (void*)&_td);
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

