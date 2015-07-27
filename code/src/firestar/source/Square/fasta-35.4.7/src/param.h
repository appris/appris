
/* $Id: param.h 67 2008-10-27 15:55:11Z pearson $ */
/* $Revision$  */

#include <sys/types.h>

#ifndef P_STRUCT
#define P_STRUCT

#define MAXSQ 56

/* Concurrent read version */

struct fastr {
  int ktup;
  int cgap;
  int pgap;
  int pamfact;
  int scfact;
  int bestoff;
  int bestscale;
  int bkfact;
  int bktup;
  int bestmax;
  int altflag;
  int optflag;
  int iniflag;
  int optcut;
  int optcut_set;
  int optwid;
};

struct prostr {
    int gopen;
    int gextend;
    int width;
};

struct pstruct		/* parameters */
{
  int n0;	/* length of query sequence, used for statistics */
  int gdelval;	/* value gap open (-10) */
  int ggapval;	/* value for additional residues in gap (-2) */
  int gshift;	/* frameshift for fastx, fasty */
  int gsubs;	/* nt substitution in fasty */
  int p_d_mat;	/* dna match penalty */
  int p_d_mis;	/* dna mismatch penalty */
  int p_d_set;	/* using match/mismatch */
  int n1_low;
  int n1_high;	/* sequence length limits */
  int score_ix;	/* index to sorted score */
  int show_ident;	/* flag - show identical lalign alignment */
  int nseq;	/* number of different sequences (for lalign) */
  int zsflag;	/* use scalebest() */
  int zsflag_f;	/* use scalebest() */
  int zs_win;
  int histint;		/* histogram interval */
  char sq[MAXSQ+1];
  int hsq[MAXSQ+1];
  int nsq;		/* length of normal sq */
  int ext_sq_set;	/* flag for using extended alphabet */
  char sqx[MAXSQ];
  int hsqx[MAXSQ+1];
  int c_nt[MAXSQ+1];
  int nsqx;	/* length of extended sq */
  int nsq_e;	/* effective nsq */
  int dnaseq;	/* -1 = not set (protein); 0 = protein; 1 = DNA; 2 = other, 3 RNA */
  int nt_align;	/* DNA/RNA alignment = 1 */
  int debug_lib;
  int tr_type;	/* codon table */
  int sw_flag;
  char pamfile[120];	/* pam file type */
  char pgpfile[120];
  int pgpfile_type;
  float pamscale;
  int pam_pssm;
  int pam_set;
  int have_pam2;
  int **pam2[2];
  int **pam2p[2];
  int pamoff;	/* offset for pam values */
  int pam_l, pam_h, pam_xx, pam_xm;	/* lowest, highest pam value */
  int pam_x_set;
  int pam_ms;		/* use a Mass Spec pam matrix */
  int maxlen;
  int max_repeat;	/* used for repeat count in ssearch34/lalign */
  int repeat_thresh;
  char *other_info;
  double e_cut;		/* cutoff for scores */
  double e_cut_r; 	/* cutoff for multiple local alignments */
  long zdb_size; 	/* force database size */
  int zdb_size_set;	/* flag for user -Z */
  int pgm_id;
  int pseudocts;
  int shuff_node;
  union {
    struct fastr fa;
    struct prostr pr;
  } param_u;
};

#include "rstruct.h"

/* the seq_record has all the invariant data about a sequence -
   sequence length, libstr, sequence itself, etc.
   it does not have the results information
   we can have 1, 2, or 6 (obsolete tfasta) results records for a sequence,
   but there will still be only one sequence record.
*/

struct seq_record {
  int n1;
  int *n1tot_p;
  unsigned char *aa1b;		/* sequence buffer */
  unsigned char *aa1_ann;	/* annotation string */
#ifdef USE_FSEEKO
  off_t lseek;
#else
  long lseek;
#endif
#ifndef PCOMPLIB
  struct lmf_str *m_file_p;
#else
  int m_seqnm;
  int seqnm;
  int wrkr;
  char *bline;
#endif
  int cont;
  long l_offset;
  long l_off;
  char libstr[MAX_UID];
#ifdef SUPERFAMNUM
  int nsfnum;
  int sfnum[10];
#endif
};

#endif	/* P_STRUCT */

#include "aln_structs.h"
