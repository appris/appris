/* Concurrent read version */

/* $Id: mw.h 28 2008-06-30 16:31:45Z pearson $ */
/* $Revision$  */

#include <sys/types.h>

#include "param.h"
#include "aln_structs.h"

#ifndef FSEEK_T_DEF
#ifndef USE_FSEEKO
typedef long fseek_t;
#else
typedef off_t fseek_t;
#endif
#endif

struct beststr {
  struct seq_record *seq;	/* sequence info */
  struct beststr *bbp_link;
  struct rstruct rst;		/* results info */

  int n1;		/* duplicate of seq.n1, used for error checking/debugging */
  int frame;		/* in buf2_str */
  double zscore;	/* the z-score mostly exists for sorting best scores */

  /* description line */
  char *bline;
  /* information for final alignment */
  struct a_struct aln_d;	/* these values are used by -m9 */

  int a_res_cnt;
  struct a_res_str a_res;	/* need only a_res, not a_res[2], because different frames
				   for the same sequence are stored separately */
  int have_ares;
  float percent, gpercent;
};

struct stat_str {
  int score;
  int n1;
  double comp;
  double H;
  double escore;
  int segnum;
  int seglen;
};


