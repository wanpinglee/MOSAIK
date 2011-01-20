// ***************************************************************************
// PosixThreads.h - provides a wrapper for WIN32 threading routines so that
//                  POSIX threads calls can be used.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#ifdef WIN32

#include <windows.h>

#define PTHREAD_CANCEL_ASYNCHRONOUS 1
#define PTHREAD_CANCEL_ENABLE       2
#define PTHREAD_CANCEL_DEFERRED     3
#define PTHREAD_CANCEL_DISABLE      4
#define PTHREAD_CANCELED            5
#define PTHREAD_COND_INITIALIZER    6
#define PTHREAD_CREATE_DETACHED     7
#define PTHREAD_CREATE_JOINABLE     8
#define PTHREAD_EXPLICIT_SCHED      9
#define PTHREAD_INHERIT_SCHED       10
#define PTHREAD_MUTEX_DEFAULT       11
#define PTHREAD_MUTEX_ERRORCHECK    12
#define PTHREAD_MUTEX_NORMAL        13
#define PTHREAD_MUTEX_INITIALIZER   14
#define PTHREAD_MUTEX_RECURSIVE     15
#define PTHREAD_ONCE_INIT           16
#define PTHREAD_PRIO_INHERIT		17
#define PTHREAD_PRIO_NONE			18
#define PTHREAD_PRIO_PROTECT		19
#define PTHREAD_PROCESS_SHARED		20	
#define PTHREAD_PROCESS_PRIVATE		21
#define PTHREAD_RWLOCK_INITIALIZER	22
#define PTHREAD_SCOPE_PROCESS		23
#define PTHREAD_SCOPE_SYSTEM		24

typedef struct {
	HANDLE handle;
	unsigned int tid;
} pthread_t;

typedef struct {
	LPSECURITY_ATTRIBUTES threadAttributes;
	SIZE_T stackSize;
	void * stackAddr;
	DWORD creationFlags;
	int detachState;
	int contentionScope;
	int policy; /*supported values: SCHED_FIFO, SCHED_RR, and SCHED_OTHER*/
	int inheritSched;
	int detach;
} pthread_attr_t;

typedef struct {
	HANDLE mutex;
	int destroyed;
	int init;
	int lockedOrReferenced;
} pthread_mutex_t;

typedef struct {
	int protocol;
	int pShared;
	int prioCeiling;
	int type;
} pthread_mutexattr_t;

#define	EPERM 1
#define	EBUSY 16
#define	EINVAL 22

#ifdef __cplusplus
extern "C" {  
#endif	

	int pthread_create(pthread_t *thread, const pthread_attr_t *attr, void *(*startup)(void *), void *params);
	int pthread_join(pthread_t thread, void **value_ptr);
	void pthread_exit(void *value_ptr);
	int pthread_mutex_init(pthread_mutex_t *mutex, const pthread_mutexattr_t *attr);
	int pthread_mutex_lock(pthread_mutex_t *mutex);
	int pthread_mutex_unlock(pthread_mutex_t *mutex);
	int pthread_mutex_destroy(pthread_mutex_t *mutex);
	int pthread_mutexattr_init(pthread_mutexattr_t *attr);
	int pthread_mutexattr_destroy(pthread_mutexattr_t *attr);
	int pthread_attr_init(pthread_attr_t *attr);
	int pthread_attr_destroy(pthread_attr_t *attr);
	int pthread_attr_setdetachstate(pthread_attr_t *attr, int detachstate);

#ifdef __cplusplus
}
#endif

#else // LINUX
#include <pthread.h>
#endif
