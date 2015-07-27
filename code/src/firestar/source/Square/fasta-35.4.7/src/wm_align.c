
/* $Id: wm_align.c 27 2008-06-30 16:27:31Z pearson $  */
/* $Revision: 28 $  */

#include <stdio.h>
#include <stdlib.h>

#include "defs.h"

#include "aln_structs.h"

struct swstr {int H, E;};

int 
NW_ALIGN(int IW, const unsigned char *B,
	 int M, int N,
	 int **W, int G, int H, int *S, int *NC);

static int
CHECK_SCORE(int IW, const unsigned char *B,
	    int M, int N,
	    int *S, int **W, int G, int H, int *nres);

int
sw_walign (int **pam2p, int n0,
	   const unsigned char *aa1, int n1,
	   int q, int r,
	   struct swstr *ss,
	   struct a_res_str *a_res
	   )
{
   const unsigned char *aa1p;
   register int i, j;
   register struct swstr *ssj;
   int e, f, h, p;
   int qr;
   int     score;
   int cost, I, J, K, L;

   qr = q + r;

   /* initialize 0th row */
   for (ssj=ss; ssj<ss+n0; ssj++) {
     ssj->H = 0;
     ssj->E = -q;
   }

   /* I = saved position in aa1
      J = saved position in aa0
   */
   score = I = J = 0;
   aa1p = aa1;
   i = 0;
   while (*aa1p) {
     h = p = 0;
     f = -q;
     /* pwaa = waa + (*aa1p++ * n0); */
     for (ssj = ss, j=0; j < n0; ssj++, j++) {
       if ((h =   h     - qr) > (f =   f     - r)) f = h;
       if ((h = ssj->H - qr) > (e = ssj->E - r)) e = h;
       /* h = p + *pwaa++; */
       h = p + pam2p[j][*aa1p];
       if (h < 0 ) h = 0;
       if (h < f ) h = f;
       if (h < e ) h = e;
       p = ssj->H;
       ssj->H = h;
       ssj->E = e;
       if (h > score) {
	 score = h;
	 I = i;
	 /* J = (int)(ssj-ss);  */
	 J = j;
       }
     }
     i++;
     aa1p++;
   }				/* done with forward pass */
   if (score <= 0) return 0;  

  /* to get the start point, go backwards */
  
   /* K = begin in aa1
      L = begin in aa0
   */
  cost = K = L = 0;
  for (ssj=ss+J; ssj>=ss; ssj--) {
    ssj->H=ssj->E= -1;
  }
  
  for (i=I,aa1p=aa1+I; i>=0; i--) {
    h = f = -1;
    p = (i == I) ? 0 : -1;
    for (ssj = ss+J, j=J; ssj>=ss; ssj--,j--) {
      f = max (f,h-q)-r;
      ssj->E=max(ssj->E,ssj->H-q)-r;
      h = max(max(ssj->E,f), p+pam2p[j][aa1[i]]);
      p = ssj->H;
      ssj->H=h;
      if (h > cost) {
	cost = h;
	K = i;
	L = (int)(ssj-ss);
	if (cost >= score) goto found;
      }
    }
  }
  
found:	

  /*  printf(" %d: L: %3d-%3d/%3d; K: %3d-%3d/%3d\n",score,L,J,n0,K,I,n1); */

  a_res->max0 = J+1; a_res->min0 = L; a_res->max1 = I+1; a_res->min1 = K;
  
  NW_ALIGN(L,&aa1[K-1],J-L+1,I-K+1,pam2p,q,r,a_res->res,&a_res->nres);

  return score;
}

#define gap(k)  ((k) <= 0 ? 0 : q+r*(k))	/* k-symbol indel cost */

/* Append "Delete k" op */
#define DEL(k)				\
{ if (*last < 0)			\
    *last = (*sapp)[-1] -= (k);		\
  else {				\
    *last = (*sapp)[0] = -(k);		\
    (*sapp)++;				\
  }					\
}

/* Append "Insert k" op */
#define INS(k)				\
{ if (*last > 0)			\
    *last = (*sapp)[-1] += (k);		\
  else {				\
    *last = (*sapp)[0] = (k);		\
    (*sapp)++;				\
  }					\
}

/* align(A,B,M,N,tb,te,last) returns the cost of an optimum conversion between
   A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */

static int 
nw_align(int iw,	/* beginning of alignment in pam2p profile */
	 const unsigned char *B,	/* second sequence aa1 */
	 int M, int N,			/* length of profile, aa1 */
	 int tb, int te,
	 int **w, int q, int r, 	/* pam2p profile, open, ext */
	 struct swstr *f_ss,		/* forward, reverse row matrix */
	 struct swstr *r_ss, 
	 int dir,			/* dir [0..3] is not currently used */
	 int **sapp, int *last)
{

  int midi, midj, type;	/* Midpoint, type, and cost */
  int midc;
  int c1, c2;

  register int   i, j;
  register int c, e, d, s;
  int qr, t, *wa;

/*   print_seq_prof(A,M,B,N,w,iw); */

/*  m = g + h;  */
  qr = q + r;

/* Boundary cases: M <= 1 or N == 0 */

  if (N <= 0) {
    if (M > 0) {DEL(M)}
    return -gap(M);
  }

  if (M <= 1) {
    if (M <= 0) { 
      INS(N)
      return -gap(N);
    }

    if (tb < te) tb = te;
    midc = (tb-r) - gap(N);
    midj = 0;
/*  wa = w[A[1]]; */
    wa = w[iw];
    for (j = 1; j <= N; j++) {
      c = -gap(j-1) + wa[B[j]] - gap(N-j);
      if (c > midc) { midc = c; midj = j;}
    }
    if (midj == 0) { DEL(1) INS(N) }
    else  {
      if (midj > 1) { INS(midj-1)}
      *last = (*sapp)[0] = 0;
      (*sapp)++;
      if (midj < N) { INS(N-midj)}
    }
    return midc;
  }

/* Divide: Find optimum midpoint (midi,midj) of cost midc */

  midi = M/2;		/* Forward phase:                          */
  f_ss[0].H = 0;	/*   Compute H(M/2,k) & E(M/2,k) for all k */
  f_ss[0].E = t = -q;
  for (j = 1; j <= N; j++) {
    f_ss[j].H = t = t-r;
    f_ss[j].E = t-q;
  }
  t = tb;
  for (i = 1; i <= midi; i++) {
    s = f_ss[0].H;
    f_ss[0].H = c = t = t-r;
    e = t-q;
/*    wa = w[A[i]]; */
    wa = w[iw+i-1];
    for (j = 1; j <= N; j++) {
      if ((c =   c   - qr) > (e =   e   - r)) e = c;
      if ((c = f_ss[j].H - qr) > (d = f_ss[j].E - r)) d = c;
      c = s + wa[B[j]];
      if (e > c) c = e;
      if (d > c) c = d;
      s = f_ss[j].H;
      f_ss[j].H = c;
      f_ss[j].E = d;
    }
  }
  f_ss[0].E = f_ss[0].H;

  r_ss[N].H = 0;		/* Reverse phase:                  */
  t = -q;			/*   Compute R(M/2,k) & S(M/2,k) for all k */

  for (j = N-1; j >= 0; j--) {
    r_ss[j].H = t = t-r;
    r_ss[j].E = t-q;
  }

  t = te;
  for (i = M-1; i >= midi; i--) {
    s = r_ss[N].H;
    r_ss[N].H = c = t = t-r;
    e = t-q;
/*    wa = w[A[i+1]]; */
    wa = w[iw+i];
    for (j = N-1; j >= 0; j--) {
      if ((c =   c   - qr) > (e =   e   - r)) { e = c; }
      if ((c = r_ss[j].H - qr) > (d = r_ss[j].E - r)) { d = c; }
      c = s + wa[B[j+1]];
      if (e > c) c = e;
      if (d > c) c = d;
      s = r_ss[j].H;
      r_ss[j].H = c;
      r_ss[j].E = d;
    }
  }
  r_ss[N].E = r_ss[N].H;

  midc = f_ss[0].H+r_ss[0].H;		/* Find optimal midpoint */
  midj = 0;
  type = 1;

  for (j = 0; j <= N; j++) {
    if ((c = f_ss[j].H + r_ss[j].H) >= midc) {
      if (c > midc || (f_ss[j].H != f_ss[j].E && r_ss[j].H == r_ss[j].E)) {
	midc = c;
	midj = j;
      }
    }
  }

  for (j = N; j >= 0; j--) {
    if ((c = f_ss[j].E + r_ss[j].E + q) > midc) {
      midc = c;
      midj = j;
      type = 2;
    }
  }

/* Conquer: recursively around midpoint */

  if (type == 1) {
    c1 = nw_align(iw,B,midi,midj,tb,-q,w,q,r,f_ss, r_ss,0,sapp,last);
    c2 = nw_align(iw+midi,B+midj,M-midi,N-midj,-q,te,w,q,r,f_ss, r_ss,1,sapp,last);
  }
  else {
    nw_align(iw,B,midi-1,midj,tb,0,w,q,r,f_ss, r_ss,2,sapp,last);
    DEL(2);
    nw_align(iw+midi+1,B+midj,M-midi-1,N-midj,0,te,w,q,r,f_ss,r_ss,3,sapp,last);
  }
  return midc;
}

/* Interface and top level of comparator */

int 
NW_ALIGN(int IW, const unsigned char *B,
	 int M, int N,
	 int **W, int G, int H, int *S, int *NC)
{ 
  struct swstr *f_ss, *r_ss;
  int *sapp, last;
  int c, ck;

  sapp = S;
  last = 0;

   if ((f_ss = (struct swstr *) calloc (N+2, sizeof (struct swstr)))
       == NULL) {
     fprintf (stderr, " *** cannot allocate f_ss array %3d\n", N+2);
     exit (1);
   }
   f_ss++;

   if ((r_ss = (struct swstr *) calloc (N+2, sizeof (struct swstr)))
       == NULL) {
     fprintf (stderr, " *** cannot allocate r_ss array %3d\n", N+2);
     exit (1);
   }
   r_ss++;

  /*   print_seq_prof(A,M,W,IW); */
  c = nw_align(IW,B,M,N,-G,-G,W,G,H,f_ss, r_ss,0,&sapp,&last);	/* OK, do it */

  ck = CHECK_SCORE(IW,B,M,N,S,W,G,H,NC);
  if (c != ck) {
    fprintf(stderr," *** Check_score error. %d != %d ***\n",c,ck);
  }

  f_ss--; r_ss--;
  free(r_ss); free(f_ss);

  return c;
}

/* CHECK_SCORE - return the score of the alignment stored in S */

static int CHECK_SCORE(int iw, const unsigned char *B,
		       int M, int N,
		       int *S, int **w,
		       int g, int h, int *NC)
{ 
  register int   i,  j, op, nc;
  int score;

  /*  print_seq_prof(A,M,w,iw); */

  score = i = j = op = nc = 0;
  while (i < M || j < N) {
    op = *S++;
    if (op == 0) {
      score += w[iw+i][B[++j]];
      i++;
      nc++;
    }
    else if (op > 0) {
      score = score - (g+op*h);
      j += op;
      nc += op;
    } else {
      score = score - (g-op*h);
      i -= op;
      nc -= op;
    }
  }
  *NC = nc;
  return score;
}

