
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/*  $Id: last_tat.c 27 2008-06-30 16:27:31Z pearson $ */
/* $Revision: 28 $  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "mm_file.h"

#include "structs.h"
#include "param.h"

#ifndef PCOMPLIB
#include "mw.h"
#else
#include "msg.h"
#include "p_mw.h"

void do_stage2(struct beststr **bptr, int nbest, const struct mngmsg *m_msg,
	       int s_func, struct qmng_str *qm_msp);
#endif


extern int (*ranlib) (char *str, int cnt,
	       fseek_t libpos, char *libstr,
	       struct lmf_str *lm_fd);

#define RANLIB (m_fptr->ranlib)

#define MAX_BLINE 200

int
re_getlib(unsigned char *, unsigned char **,int, int, int, int, int, long *, long *, 
	  struct lmf_str *m_fptr);

void
do_work(unsigned char *aa0, int n0, unsigned char *aa1, int n1, int frame, 
	struct pstruct *ppst, void *f_str, int qr_flg, struct rstruct *rst);

extern void
do_opt (unsigned char *aa0, int n0, unsigned char *aa1, int n1,
	int frame, struct pstruct *pst, void *f_str,
	struct rstruct *rst);

struct lmf_str *re_openlib(struct lmf_str *, int outtty);

void sortbestz (struct beststr **bptr, int nbest);

double zs_to_E(double zs,int n1, int isdna, long entries, struct db_str db);

double scale_one_score(int ipos, double escore, struct db_str db, void *rs_str);

void sortbests (struct beststr **bptr, int nbest)
{
    int gap, i, j;
    struct beststr *tmp;

    for (gap = nbest/2; gap > 0; gap /= 2)
	for (i = gap; i < nbest; i++)
	    for (j = i - gap; j >= 0; j-= gap) {
	      if (bptr[j]->rst.score[0] >= bptr[j + gap]->rst.score[0]) break;
	      tmp = bptr[j];
	      bptr[j] = bptr[j + gap];
	      bptr[j + gap] = tmp;
	    }
}

int
last_calc(
#ifndef PCOMPLIB
	  unsigned char **aa0, unsigned char *aa1save, int maxn,
#endif	  
	  struct beststr **bptr, int nbest, 
	  const struct mngmsg *m_msg, struct pstruct *ppst
#ifdef PCOMPLIB
	  , struct qmng_str *qm_msp
#else
	  , void **f_str
#endif
	  , void *rstat_str)
{
  unsigned char *aa1;
  int nopt, ib;
  struct beststr *bbp;
  long loffset, l_off;
  int n0, n1;
  struct rstruct rst;
  struct lmf_str *m_fptr;
  char bline[60];
  int tat_samp, tat_inc, loop_cnt, i;
  double min_escore, ess;

  n0 = m_msg->n0;

  sortbestz(bptr,nbest);

  tat_inc = 500;
/*
  if (zs_to_E(bptr[0]->zscore,bptr[0]->n1,0,ppst->zdb_size,m_msg->db)/ 
      zs_to_E(bptr[nbest-1]->zscore,bptr[nbest-1]->n1,0,ppst->zdb_size,m_msg->db)
	< 1e-20) { tat_inc /= 4 ;}
*/

/* || (zs_to_E(bptr[0]->zscore,bptr[0]->n1,0,ppst->zdb_size,m_msg->db)< 1e-5); */

  ib = tat_samp = 0;
  for (loop_cnt = 0; loop_cnt < 5; loop_cnt++) {
    tat_samp += tat_inc;
    nopt = min(nbest,tat_samp);
    min_escore = 1000000.0;
#ifndef PCOMPLIB
    for ( ; ib<nopt; ib++) {
      bbp = bptr[ib];

      if (bbp->rst.score[0] < 0) break;

      if ((m_fptr=re_openlib(bbp->seq->m_file_p,!m_msg->quiet))==NULL) {
	fprintf(stderr,"*** cannot re-open %s\n",bbp->seq->m_file_p->lb_name);
	exit(1);
      }
      RANLIB(bline,sizeof(bline),bbp->seq->lseek,bbp->seq->libstr,m_fptr);

      if (bbp->seq->aa1b == NULL) {
	n1 = re_getlib(aa1save,NULL,maxn,m_msg->maxt3,
		       m_msg->l_overlap,bbp->seq->cont,m_msg->term_code,
		       &loffset,&l_off,bbp->seq->m_file_p);
	aa1 = aa1save;
      }
      else {
	n1 = bbp->seq->n1;
	aa1 = bbp->seq->aa1b;
	loffset = bbp->seq->l_offset;
	l_off = bbp->seq->l_off;
      }

      do_opt(aa0[bbp->frame],m_msg->n0,aa1,n1,bbp->frame,ppst,
	     f_str[bbp->frame],&rst);
      memcpy(&(bbp->rst),&rst,sizeof(struct rstruct));

      if ((ess=scale_one_score(ib, bbp->rst.escore, m_msg->db, rstat_str)) < 
	min_escore) { min_escore = ess;}
      /*
      fprintf(stderr,"%d: %4d %2d %3d %.4g %.4g\n",
	      ib, bbp->rst.score[0], bbp->segnum,bbp->seglen,bbp->escore, ess);
      */
    }
#else
    do_stage2(&bptr[ib], nopt-ib, m_msg, DO_CALC_FLG, qm_msp);

    for ( ; ib < nopt; ib++) {
      if ((ess=scale_one_score(ib, bptr[ib]->rst.escore, m_msg->db, rstat_str)) < 
	min_escore) { min_escore = ess;}
      /*
      fprintf(stderr, "%d: %4d %2d %3d %.4g %.4g\n",
	      ib,bptr[ib]->rst.score[0],bptr[ib]->segnum,bptr[ib]->seglen,bptr[ib]->rst.escore,ess);
      */
    }
#endif


    if (min_escore > m_msg->e_cut) goto done;
  }
 done:
  return ib;
}
