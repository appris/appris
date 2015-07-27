/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/*  $Id: work_thr2.c 27 2008-06-30 16:27:31Z pearson $ */
/* $Revision: 28 $  */

/* work_thr.c - threaded worker */

/* modified 21-Oct-1998 to work with reverse complement for DNA */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <signal.h>

#include "defs.h"		/* various constants */
#include "mw.h"			/* defines beststr */
#include "structs.h"
#include "param.h"		/* pstruct rstruct */
#include "thr_bufs.h"

/***************************************/
/* thread global variable declarations */
/***************************************/

#define XTERNAL
#include "thr.h"
#undef XTERNAL

void alloc_pam (int, int, struct pstruct *);
int **alloc_pam2p(int **,int, int);
void revcomp(unsigned char *seq, int n, int *c_nt);
#ifdef WIN32
void pthread_exit(void *);
#else
void THR_EXIT(void *);
#endif

#ifdef DEBUG
extern struct buf_head *lib_buf2_list;
#endif

/* functions getting/sending buffers to threads (thr_sub.c) */
extern void wait_thr(void);
extern int get_wbuf(struct buf_head **cur_buf, int max_work_buf);
extern void put_wbuf(struct buf_head *cur_buf, int max_work_buf);

/* dropxx.c functions */
extern void init_work (unsigned char *aa0, int n0,
		       struct pstruct *ppst, void **f_arg);

extern void do_work (unsigned char *aa0, int n0, unsigned char *aa1, int n1,
		     int frame,
		     struct pstruct *ppst, void *f_str, int qr_flg,
		     struct rstruct *rst);

extern void close_work (unsigned char *, int, struct pstruct *, void **);

extern void irand(int);
extern int shuffle(unsigned char *, unsigned char *, int);
extern int wshuffle(unsigned char *, unsigned char *, int, int);
extern void qshuffle(unsigned char *aa0, int n0, int nm0);
extern void free_pam2p(int **);

extern void
buf_do_work(unsigned char **aa0, int n0, struct buf_head *lib_bhead_p,
	    struct pstruct *ppst, void **f_str);
extern void
buf_qshuf_work(unsigned char *aa0s, int n0, struct buf_head *lib_bhead_p,
	       struct pstruct *ppst, void *qf_str, int score_ix);
extern void
buf_shuf_work(unsigned char **aa0, int n0, unsigned char *aa1s, int zs_win,
	      struct buf_head *lib_bhead_p, struct pstruct *ppst,
	      void **f_str, int score_ix);
void
work_thread (struct thr_str *work_info)
{
  struct buf_head *cur_buf;
  struct buf2_str *cur_buf_p;
  struct buf2_str *p_rbuf;
  char info_lib_range[MAX_SSTR];
  unsigned char *aa1s=NULL;
  int cur_cnt;
  int my_worker;
  int i, j, npam, n0, nm0;
  int debug_lib, do_shuffle;
  int frame;
  int icnt;

  struct rstruct rrst;
  struct pstruct my_pst, *my_ppst;
  unsigned char *aa0[6], *aa0s;
  void *f_str[6], *qf_str;

  my_worker = work_info->worker;

  wait_thr();	/* wait for start_thread predicate to drop to  0 */

/* **************************************************************** */
/* let each thread have its own copy of the query */

  n0 = work_info->n0;
  nm0 = work_info->nm0;

  if ((aa0[0]=(unsigned char *)calloc((size_t)n0+2,sizeof(unsigned char)))
      ==NULL) {
    fprintf(stderr," cannot allocate aa00[%d] for worker %d\n",
	    n0, my_worker);
    exit(1);
  }
  *aa0[0]='\0';
  aa0[0]++;
  memcpy(aa0[0],work_info->aa0,n0+1);

  /* make certain that all but 0 have their own copy of pst */
  if (my_worker) {
    my_ppst = &my_pst;
    memcpy(my_ppst,work_info->ppst,sizeof(struct pstruct));
    my_ppst->pam2p[0] = my_ppst->pam2p[1] = NULL;

    alloc_pam(MAXSQ, MAXSQ, my_ppst);

    npam = (my_pst.ext_sq_set) ? my_pst.nsqx : my_pst.nsq;

    /* allocate local copy of pam2[][] */
    for (i=0; i<=npam; i++) {
      for (j=0; j<=npam; j++) {
	my_pst.pam2[0][i][j] = work_info->ppst->pam2[0][i][j];
	my_pst.pam2[1][i][j] = work_info->ppst->pam2[1][i][j];
      }
    }

    /* if this is a pssm search, allocate local copy of pam2p[][]*/
    if (work_info->ppst->pam_pssm && work_info->ppst->pam2p[0]) {
      my_ppst->pam2p[0] = alloc_pam2p(my_ppst->pam2p[0],n0,npam);
      my_ppst->pam2p[1] = alloc_pam2p(my_ppst->pam2p[1],n0,npam);

      for (i=0; i<n0; i++) {
	for (j=0; j <= npam; j++) {
	  my_pst.pam2p[0][i][j] = work_info->ppst->pam2p[0][i][j];
	  my_pst.pam2p[1][i][j] = work_info->ppst->pam2p[1][i][j];
	}
      }
    }
  }
  else my_ppst=work_info->ppst;

  /* note that aa[5,4,3,2] are never used, but are provided so that frame
     can range from 0 .. 5; likewise for f_str[5..2] */

  aa0[5] = aa0[4] = aa0[3] = aa0[2] = aa0[1] = aa0[0];

  init_work (aa0[0], n0, my_ppst, &f_str[0]);

  if (my_worker == 0) {
    /* label library size limits */
    if (my_ppst->n1_low > 0 && my_ppst->n1_high < BIGNUM) 
      sprintf(info_lib_range," (range: %d-%d)",my_ppst->n1_low,my_ppst->n1_high);
    else if (my_ppst->n1_low > 0) 
      sprintf(info_lib_range," (range: >%d)",my_ppst->n1_low);
    else if (my_ppst->n1_high < BIGNUM)
      sprintf(info_lib_range," (range: <%d)",my_ppst->n1_high);
    else
      info_lib_range[0]='\0';
    info_lib_range[sizeof(info_lib_range)-1]='\0';
    strncpy(work_info->info_lib_range,info_lib_range,MAX_SSTR);
  }

  f_str[5] = f_str[4] = f_str[3] = f_str[2] = f_str[1] = f_str[0];

  if (work_info->qframe == 2) {
    if ((aa0[1]=(unsigned char *)calloc((size_t)n0+2,sizeof(unsigned char)))==NULL) {
      fprintf(stderr," cannot allocate aa01[%d] for worker %d\n",
	    n0, my_worker);
    }
    *aa0[1]='\0';
    aa0[1]++;
    memcpy(aa0[1],work_info->aa0,n0+1);
    revcomp(aa0[1],n0,my_ppst->c_nt);
    init_work (aa0[1], n0, my_ppst, &f_str[1]);
  }

  if (work_info->qshuffle) {
    if ((aa0s=(unsigned char *)calloc(n0+2,sizeof(char)))==NULL) {
      fprintf(stderr,"cannot allocate aa0s[%d]\n",n0+2);
      exit(1);
    }
    *aa0s='\0';
    aa0s++;
    memcpy(aa0s,aa0[0],n0);
    qshuffle(aa0s,n0,nm0);
    init_work (aa0s, n0, my_ppst, &qf_str);
  }

  /* always allocate shuffle space */
  if((aa1s=calloc(work_info->max_tot+1,sizeof(char))) == NULL) {
    fprintf(stderr,"unable to allocate shuffled library sequence\n");
  }
  else {
    *aa1s=0;
    aa1s++;
    irand(0);
  }

  if (my_ppst->zsflag >= 10) { do_shuffle = 1; }
  else { do_shuffle = 0; }

/* **************************************************************** */
/* main work loop */

  while (get_wbuf(&cur_buf,work_info->max_work_buf)) {

    if (cur_buf->buf2_cnt <= 0) break;

    if (cur_buf->buf2_type & BUF2_DOWORK) {
      buf_do_work(aa0, n0, cur_buf, my_ppst, f_str);

      if (work_info->qshuffle) {
	buf_qshuf_work(aa0s, n0, cur_buf, my_ppst, qf_str, my_ppst->score_ix);
      }
    }

    if (cur_buf->buf2_type & BUF2_DOSHUF) {
      buf_shuf_work(aa0, n0, aa1s, my_ppst->zs_win,
		    cur_buf, my_ppst, f_str, my_ppst->score_ix);
    }

    /* 
    if (cur_buf->buf2_type & BUF2_DOOPT) {
      buf_do_opt(aa0, n0, cur_buf, my_ppst, f_str);
    }
    */
    /*
    if (cur_buf->buf2_type & BUF2_DOALIGN) {
      buf_do_align(aa0, n0, cur_buf, my_ppst, f_str);
    }
    */
    cur_buf->have_results = 1;
    cur_buf->have_data = 0;

    put_wbuf(cur_buf,work_info->max_work_buf);

  } /* end main while */

/* **************************************************************** */
/* all done - clean-up */

  close_work(aa0[0], n0, my_ppst, &f_str[0]);
  free(aa0[0]-1);
  if (work_info->qframe == 2) {
    close_work(aa0[1], n0, my_ppst, &f_str[1]);
    free(aa0[1]-1);
  }

  if (work_info->qshuffle) {
    close_work(aa0s, n0, my_ppst, &qf_str);
    free(aa0s-1);
  }

  free(aa1s-1);

  if (my_worker) {
    free(my_pst.pam2[1][0]);
    free(my_pst.pam2[0][0]);
    free(my_pst.pam2[1]);
    free(my_pst.pam2[0]);
  }

  if (my_worker && my_pst.pam_pssm) {
    free_pam2p(my_pst.pam2p[0]);
    free_pam2p(my_pst.pam2p[1]);
  }

/* **************************************************************** */
/* and exit */

#ifdef WIN32
  pthread_exit(&work_info->status);
#else
  THR_EXIT(&work_info->status);
#endif

}  /* end work_thread */

