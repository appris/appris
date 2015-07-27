/* Concurrent read version */

/* $Id: p_mw.h 28 2008-06-30 16:31:45Z pearson $ */
/* $Revision$  */

#ifndef FSEEK_T_DEF
#ifndef USE_FSEEKO
typedef long fseek_t;
#else
typedef off_t fseek_t;
#endif
#endif

struct beststr {
  struct seq_record *seq;		/* sequence number */
  struct rstruct rst;		/* results info */
  int n1;

  int frame;
  double zscore;

  int rscore;	/* score from shuffled sequence */
  int sw_score;	/* optimal score from alignment */

  double r_escore;
  struct a_struct aln_d;
  float percent, gpercent;

  int a_res_cnt;
  char *aln_code;
  int aln_code_n;
  char *ann_code;
  int ann_code_n;
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

/* this structure passes library sequences to the worker threads
   and returns scores */

#include "w_mw.h"

/*
struct pbuf_head {
  int buf_cnt;
  unsigned char *start;
  struct sqs2 *buf;
};
*/
