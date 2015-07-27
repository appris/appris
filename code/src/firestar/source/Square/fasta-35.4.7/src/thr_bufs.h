
/* thr_bufs.h - structures for passing buffers of sequences to threads */

/* copyright (c) 2007 William R. Pearson and the U. of Virginia */

/* $Id: thr_bufs.h 28 2008-06-30 16:31:45Z pearson $ */
/* $Revision$  */

#include <sys/types.h>

struct thr_str {
  int worker;
  void *status;
  int max_work_buf;
  int qframe;
  struct pstruct *ppst;
  int qshuffle;
  int nrshuffle;
  unsigned char *aa0;
  int n0;
  int nm0;
  int max_tot;
  char info_lib_range[MAX_SSTR];
};


/* this structure passes library sequences to the worker threads
   and returns scores */

struct buf2_str {
  struct seq_record *seq;
  struct seq_record *seq_b;
  struct rstruct rst;
  int qframe;
  int frame;
  int r_score, qr_score;
  double r_escore, qr_escore;
  struct beststr *best_save;
  int have_ares;
  struct a_res_str a_res;
};

#define BUF2_DOWORK 0x1
#define BUF2_DOSHUF 0x2
#define BUF2_DOOPT 0x4
#define BUF2_DOALIGN 0x8

struct buf_head {
  int buf2_cnt;
  int buf2_type;	/* specify the type of computation - search, shuffle, etc */
  int have_data;
  int have_results;
  int have_best_save;
  unsigned char *start;
  struct buf2_str *buf2;
};
