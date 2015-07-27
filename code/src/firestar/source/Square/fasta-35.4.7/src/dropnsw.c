/* copyright (c) 1994, 1995, 1996 William R. Pearson */

/*  $Id: dropnsw.c 27 2008-06-30 16:27:31Z pearson $ */
/* $Revision: 87 $  */

/*
  this is a slower version of dropgsw.c that implements the Smith-Waterman
  algorithm.  It lacks the shortcuts in dropgsw.c that prevent scores less
  than the penalty for the first residue in a gap from being generated.

  Thus, dropnsw.c should be used for tests with very large gap penalties,
  and is more appropriate for programs like prss3, which are interested
  in accurate low scores.
*/

/* the do_walign() code in this file is not thread_safe */
/* init_work(), do_work(), are thread safe */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "defs.h"
#include "param.h"

static char *verstr="3.5 Sept 2006";

struct swstr { int H, E;};

struct f_struct {
  struct swstr *ss;
  struct swstr *f_ss;
  struct swstr *r_ss;
  int *waa_s, *waa_a;
  int **pam2p[2];
  int max_res;
  double aa0_f[MAXSQ];
  double *kar_p;
};

#define DROP_INTERN
#include "drop_func.h"

extern int do_karlin(const unsigned char *aa1, int n1,
		     int **pam2, const struct pstruct *ppst,
		     double *aa0_f, double *kar_p, double *lambda, double *H);
extern void aancpy(char *to, char *from, int count, struct pstruct *ppst);
int ALIGN(const unsigned char *A, const unsigned char *B, int M, int N,
	  int **W, int IW, int G, int H, int *S, int *NC,
	  struct f_struct *f_str);

/* initialize for Smith-Waterman optimal score */

void init_work (unsigned char *aa0, int n0,
		struct pstruct *ppst,
		struct f_struct **f_arg)
{
   int maxn0;
   int *pwaa_s, *pwaa_a;
   int e, f, i, j, q;
   int *res;
   struct f_struct *f_str;
   int **pam2p;
   struct swstr *ss, *f_ss, *r_ss;
   int nsq, ip;

   if (ppst->ext_sq_set) {
     nsq = ppst->nsqx; ip = 1;
   }
   else {
     nsq = ppst->nsq; ip = 0;
   }

   f_str = (struct f_struct *)calloc(1,sizeof(struct f_struct));

   /* allocate space for the scoring arrays */
   maxn0 = n0 + 2;
   if ((ss = (struct swstr *) calloc (maxn0, sizeof (struct swstr)))
	 == NULL) {
     fprintf (stderr, "cannot allocate ss array %3d\n", n0);
     exit (1);
   }
   ss++;
   f_str->ss = ss;

   if ((f_ss = (struct swstr *) calloc (maxn0, sizeof (struct swstr)))
       == NULL) {
     fprintf (stderr, "cannot allocate f_ss array %3d\n", n0);
     exit (1);
   }
   f_ss++;
   f_str->f_ss = f_ss;

   if ((r_ss = (struct swstr *) calloc (n0+2, sizeof (struct swstr)))
       == NULL) {
     fprintf (stderr, "cannot allocate r_ss array %3d\n", n0);
     exit (1);
   }
   r_ss++;
   f_str->r_ss = r_ss;

   /* initialize variable (-S) pam matrix */
   if ((f_str->waa_s= (int *)calloc((nsq+1)*(n0+1),sizeof(int))) == NULL) {
     fprintf(stderr,"cannot allocate waa_s array %3d\n",nsq*n0);
     exit(1);
   }

   if ((f_str->pam2p[1]= (int **)calloc((n0+1),sizeof(int *))) == NULL) {
     fprintf(stderr,"cannot allocate pam2p[1] array %3d\n",n0);
     exit(1);
   }

   pam2p = f_str->pam2p[1];
   if ((pam2p[0]=(int *)calloc((nsq+1)*(n0+1),sizeof(int))) == NULL) {
     fprintf(stderr,"cannot allocate pam2p[1][] array %3d\n",nsq*n0);
     exit(1);
   }

   for (i=1; i<n0; i++) {
     pam2p[i]= pam2p[0] + (i*(nsq+1));
   }

   /* initialize universal (alignment) matrix */
   if ((f_str->waa_a= (int *)calloc((nsq+1)*(n0+1),sizeof(int))) == NULL) {
     fprintf(stderr,"cannot allocate waa_a struct %3d\n",nsq*n0);
     exit(1);
   }

   if ((f_str->pam2p[0]= (int **)calloc((n0+1),sizeof(int *))) == NULL) {
     fprintf(stderr,"cannot allocate pam2p[1] array %3d\n",n0);
     exit(1);
   }

   pam2p = f_str->pam2p[0];
   if ((pam2p[0]=(int *)calloc((nsq+1)*(n0+1),sizeof(int))) == NULL) {
     fprintf(stderr,"cannot allocate pam2p[1][] array %3d\n",nsq*n0);
     exit(1);
   }

   for (i=1; i<n0; i++) {
     pam2p[i]= pam2p[0] + (i*(nsq+1));
   }

   /* 
      pwaa effectively has a sequence profile --
       pwaa[0..n0-1] has pam score for residue 0 (-BIGNUM)
       pwaa[n0..2n0-1] has pam scores for residue 1 (A)
       pwaa[2n0..3n-1] has pam scores for residue 2 (R), ...

       thus: pwaa = f_str->waa_s + (*aa1p++)*n0; sets up pwaa so that
       *pwaa++ rapidly moves though the scores of the aa1p[] position
       without further indexing

       For a real sequence profile, pwaa[0..n0-1] vs ['A'] could have
       a different score in each position.
   */

   if (ppst->pam_pssm) {
     pwaa_s = f_str->waa_s;
     pwaa_a = f_str->waa_a;
     for (e = 0; e <=nsq; e++)	{	/* for each residue in the alphabet */
       for (f = 0; f < n0; f++) {	/* for each position in aa0 */
	 *pwaa_s++ = f_str->pam2p[ip][f][e] = ppst->pam2p[ip][f][e];
	 *pwaa_a++ = f_str->pam2p[0][f][e]  = ppst->pam2p[0][f][e];
       }
     }
   }
   else {	/* initialize scanning matrix */
     pwaa_s = f_str->waa_s;
     pwaa_a = f_str->waa_a;
     for (e = 0; e <=nsq; e++)	/* for each residue in the alphabet */
       for (f = 0; f < n0; f++)	{	/* for each position in aa0 */
	 *pwaa_s++ = f_str->pam2p[ip][f][e]= ppst->pam2[ip][e][aa0[f]];
	 *pwaa_a++ = f_str->pam2p[0][f][e] = ppst->pam2[0][e][aa0[f]];
       }
   }

   /* minimum allocation for alignment */
   f_str->max_res =  max(3*n0/2,MIN_RES);

   *f_arg = f_str;
}

void close_work (const unsigned char *aa0, int n0,
		 struct pstruct *ppst, struct f_struct **f_arg)
{
  struct f_struct *f_str;

  f_str = *f_arg;

  if (f_str != NULL) {
    if (f_str->kar_p !=NULL) free(f_str->kar_p);
    f_str->ss--;
    free(f_str->ss);
    free(f_str->waa_a);
    free(f_str->pam2p[0][0]);
    free(f_str->pam2p[0]);
    free(f_str->waa_s);
    free(f_str->pam2p[1][0]);
    free(f_str->pam2p[1]);

    free(f_str);
    *f_arg = NULL;
  }
}

/* pstring1 is a message to the manager, currently 512 */
/*void get_param(struct pstruct *pstr,char **pstring1)*/
void    get_param (struct pstruct *ppst, char **pstring1, char *pstring2)
{
  char psi_str[120];

  char *pg_str="Smith-Waterman";

  if (ppst->pam_pssm) { strncpy(psi_str,"-PSI",sizeof(psi_str));}
  else { psi_str[0]='\0';}

  sprintf (pstring1[0], " %s (%s)", pg_str, verstr);
  sprintf (pstring1[1],
#ifdef OLD_FASTA_GAP
	   "%s matrix%s (%d:%d)%s, gap-penalty: %d/%d",
#else
	   "%s matrix%s (%d:%d)%s, open/ext: %d/%d",
#endif
	   ppst->pamfile, psi_str, ppst->pam_h,ppst->pam_l, 
	   (ppst->ext_sq_set)?"xS":"\0", ppst->gdelval, ppst->ggapval);

   if (pstring2 != NULL) {
     sprintf(pstring2,
#ifdef OLD_FASTA_GAP
	     "; pg_name_alg: %s\n; pg_ver_rel: %s\n; pg_matrix: %s (%d:%d)%s\n; pg_gap-pen: %d %d\n",
#else
	     "; pg_name_alg: %s\n; pg_ver_rel: %s\n; pg_matrix: %s (%d:%d)%s\n; pg_open-ext: %d %d\n",
#endif
	     pg_str,verstr,psi_str,ppst->pam_h,ppst->pam_l, 
	     (ppst->ext_sq_set)?"xS":"\0",ppst->gdelval,ppst->ggapval);
   }
}


void do_work (const unsigned char *aa0, int n0,
	      const unsigned char *aa1, int n1,
	      int frame,
	      const struct pstruct *ppst, struct f_struct *f_str,
	      int qr_flg,
	      struct rstruct *rst)
{
   const unsigned char *aa0p, *aa1p;
   register struct swstr *ssj;
   struct swstr *ss, *f_ss, *r_ss;
   register int *pwaa;
   int *waa;
   register int i, j;
   int     e, f, h, p;
   int     q, r, m;
   int     score;

   double lambda, H, K;

   rst->escore = 1.0;
   rst->segnum = rst->seglen = 1;

   waa = f_str->waa_s;
   ss = f_str->ss;
   f_ss = f_str->f_ss;
   r_ss = f_str->r_ss;

#ifdef OLD_FASTA_GAP
   q = -(ppst->gdelval - ppst->ggapval);
#else
   q = -ppst->gdelval;
#endif
   r = -ppst->ggapval;
   m = q + r;

   /* initialize 0th row */
   for (ssj=ss; ssj<&ss[n0]; ssj++) {
     ssj->H = 0;
     ssj->E = -q;
   }

   score = 0;
   aa1p = aa1;
   while (*aa1p) {
     h = p = 0;
     f = -q;
     pwaa = waa + (*aa1p++ * n0);
     for (ssj = ss, aa0p = aa0; ssj < ss+n0; ssj++) {
       if ((h =   h     - m) > (f =   f     - r)) f = h;
       if ((h = ssj->H - m) > (e = ssj->E - r)) e = h;
       h = p + *pwaa++;
       if (h < 0 ) h = 0;
       if (h < f ) h = f;
       if (h < e ) h = e;
       p = ssj->H;
       ssj->H = h;
       ssj->E = e;
       if (h > score) score = h;
     }
   }				/* done with forward pass */

   rst->score[0] = score;

   if(ppst->zsflag == 6 || ppst->zsflag == 16 && 
     (do_karlin(aa1, n1, ppst->pam2[0], ppst,f_str->aa0_f, 
		f_str->kar_p, &lambda, &H)>0)) {
     rst->comp = 1.0/lambda;
     rst->H = H;
   }
  else {rst->comp = rst->H = -1.0;}
}				/* here we should be all done */

void    do_opt (const unsigned char *aa0, int n0,
		const unsigned char *aa1, int n1,
		int frame,
		struct pstruct *pst, struct f_struct *f_str,
		struct rstruct *rstr)
{
}

void
do_walign (const unsigned char *aa0, int n0,
	   const unsigned char *aa1, int n1,
	   int frame,
	   struct pstruct *ppst, 
	   struct f_struct *f_str, 
	   struct a_res_str *a_res,
	   int *have_ares )
{
  *have_ares = 0x3;	/* set 0x2 bit to indicate local copy */
  
  /* now we need alignment storage - get it */
  if ((a_res->res = (int *)calloc((size_t)f_str->max_res,sizeof(int)))==NULL) {
    fprintf(stderr," *** cannot allocate alignment results array %d\n",f_str->max_res);
    exit(1);
   }

  a_res->next = NULL;

  a_res->sw_score = sw_walign(f_str->pam2p[0], n0, aa1, n1, 
#ifdef OLD_FASTA_GAP
			      -(ppst->gdelval - ppst->ggapval),
#else
			      -ppst->gdelval,
#endif
			      -ppst->ggapval,

			      f_str->ss, a_res);
}

/*
#define XTERNAL
#include "upam.h"

void
print_seq_prof(unsigned char *A, int M,
	       unsigned char *B, int N,
	       int **w, int iw, int dir) {
  char c_max;
  int i_max, j_max, i,j;

  char *c_dir="LRlr";

  for (i=1; i<=min(60,M); i++) {
    fprintf(stderr,"%c",aa[A[i]]);
  }
  fprintf(stderr, - %d\n,M);

  for (i=0; i<min(60,M); i++) {
    i_max = -1;
    for (j=1; j<21; j++) {
      if (w[iw+i][j]> i_max) {
	i_max = w[iw+i][j]; 
	j_max = j;
      }
    }
    fprintf(stderr,"%c",aa[j_max]);
  }
  fputc(':',stderr);
  for (i=1; i<=min(60,N); i++) {
    fprintf(stderr,"%c",aa[B[i]]);
  }
  fprintf(stderr," -%c: %d,%d\n",c_dir[dir],M,N);
}
*/

void
pre_cons(const unsigned char *aa1, int n1, int frame, struct f_struct *f_str) {}

/* aln_func_vals - set up aln.qlfact, qlrev, llfact, llmult, frame, llrev */
/* call from calcons, calc_id, calc_code */
void 
aln_func_vals(int frame, struct a_struct *aln) {

  aln->llfact = aln->llmult = aln->qlfact = 1;
  aln->qlrev = aln->llrev = 0;
  aln->frame = 0;
}

#ifdef PCOMPLIB
#include "p_mw.h"
void
update_params(struct qmng_str *qm_msg, struct pstruct *ppst)
{
  ppst->n0 = qm_msg->n0;
}
#endif
