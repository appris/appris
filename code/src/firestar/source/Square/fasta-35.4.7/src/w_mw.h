
/* $Id: w_mw.h 28 2008-06-30 16:31:45Z pearson $ */
/* $Revision$  */

/* 21-July-2000 - changes for p2_complib/p2_workcomp:
   there are now two sequence numbers; the old (worker) seqnm,
   and a new manager (master) sequence number, m_seqnm
*/

#ifndef BFR
#define BFR	300
#endif
#ifndef BFR2
#define BFR2	100
#endif

#define MAXSQL	125000
#define MMAXSQL 2000000
#ifndef MAXWRKR
#define MAXWRKR	64
#endif
#define MAXLSEQ	50000
#define DESLIN	60
#define NDES	100

struct qmng_str
{
  int n0;	/* query sequence length */
  int nm0;	/* number of segments */
  int escore_flg;	/* use escores */
  int qshuffle;		/* query shuffle */
  int pam_pssm;	/* flag for pssm/profile search */
  int s_func;	/* for p_workcomp: func==0>simple comparison, ==1>alignments */
  int slist;	/* number of alignments to do */
  int seqnm;	/* query sequence number - used for identity searches */
  int have_ann;	/* have an annotation string for this query */
  char libstr[MAX_FN];
};

struct comstr
{
  int m_seqnm;		/* sequence number */
  int seqnm;		/* sequence number */
  struct rstruct rst;
  int frame;
  int r_score, qr_score;
  double r_escore, qr_escore;
};

struct comstr2
{
  int m_seqnm;		/* sequence number */
  int seqnm;		/* sequence number */
  struct rstruct rst;
  int sw_score;

  /* int a_len; */	/* consensus alignment length */
  /* int min0, max0, min1, max1;
  int nident, ngap_q, ngap_l; */ /* number of identities, gaps in q, l */
  
  struct a_struct aln_d;
  float percent, gpercent;
  int aln_code_n;
  int ann_code_n;
};

/* The message structure */

struct wrkmsg
{
    char lname [80];	/* name of the library */
    char libenv[80];	/* directory in which library resides */
    int lb_off;		/* offset in the library */
    int lb_stop;	/* stop position in library */
    int lb_code;	/* continue code */
    int lb_size;	/* library size */
    int p_size;		/* parcel size */
    int libfn;		/* current library being searched */
    int stage;		/* current stage number */
};

struct sqs
{
  int n1;		/* size of library sequence */
  unsigned char *aa1; 	/* sequence data */
};

#include "aln_structs.h"

struct sqs2
{
  int n1;		/* size of library sequence */
  int m_seqnm;		/* location in master list */
  unsigned char *aa1;
  int walign_dflg[2];
  int sw_score[2];
  struct a_res_str a_res[2];	/* need a_res for each frame */
};

struct stage2_str {
  int m_seqnm;	/* manager sequence number */
  int seqnm;	/* worker sequence number */
  int frame;	/* query frame */
};
