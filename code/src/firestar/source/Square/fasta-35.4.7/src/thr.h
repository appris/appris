
/***************************************/
/* thread global variable declarations */
/***************************************/

/* $Id: thr.h 28 2008-06-30 16:31:45Z pearson $ */
/* $Revision$  */

#ifndef MAX_WORKERS
#define MAX_WORKERS 2
#endif

#ifndef XTERNAL
struct buf_head **worker_buf;  /* pointers to full buffers */
struct buf_head **reader_buf;  /* pointers to empty buffers */

/* protected by worker_mutex/woker_cond_var */
int worker_buf_workp, worker_buf_readp; /* indices into full-buffers ptrs */
int num_worker_bufs, reader_done;

/* protected by reader_mutex/reader_cond var */
int reader_buf_workp, reader_buf_readp; /* indices into empty-buffers ptrs */
int num_reader_bufs;
int reader_wait;

/* protected by start_mutex/start_cont_var */
int start_thread=1;        /* start-up predicate, 0 starts */
#else
extern struct buf_head **worker_buf;
extern struct buf_head **reader_buf;
extern int num_worker_bufs, reader_done;
extern int  num_reader_bufs, reader_wait;
extern int worker_buf_workp, worker_buf_readp;
extern int reader_buf_workp, reader_buf_readp;

extern int start_thread;
#endif

