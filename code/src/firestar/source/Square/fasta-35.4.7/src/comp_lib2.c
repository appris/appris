/* copyright (c) 1996, 1997, 1998, 1999, 2002  William R. Pearson and the
   U. of Virginia */

/*  $Id: comp_lib2.c 27 2008-06-30 16:27:31Z pearson $ */
/*  $Revision: 116 $  */

/*
 * Jan 17, 2007 - remove #ifdef PRSS - begin better statistics in place
 * for small libraries, related libraries
 *
 * Concurrent read version
 *
 *	Feb 20, 1998 modifications for prss3
 *
 *	December, 1998 - DNA searches are now down with forward and reverse
 *			 strands
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include <limits.h>
#include <float.h>
#include <math.h>

#ifdef UNIX
#include <unistd.h>
#include <sys/types.h>
#include <signal.h>
#endif

#include "defs.h"
#include "mm_file.h"

#include "mw.h"			/* defines beststr */
#include "structs.h"		/* mngmsg, libstruct */
#include "param.h"		/* pstruct, thr_str, buf_head, rstruct */
#include "thr_bufs.h"
#include "drop_func.h"

#define XTERNAL
#include "uascii.h"

char *mp_verstr="35.04";

/********************************/
/* global variable declarations */
/********************************/

extern int fa_max_workers;

#ifdef SUPERFAMNUM
int nsfnum;
int sfnum[10];
extern int sfn_cmp(int *q, int *s);
int nsfnum_n;
int sfnum_n[10];
#endif

/********************************/
/* extern variable declarations */
/********************************/
extern char *prog_func;		/* function label */
extern char *verstr, *iprompt0, *iprompt1, *iprompt2, *refstr;

/********************************/
/*extern function declarations  */
/********************************/
struct lmf_str *openlib(char *, int, int *, int, struct lmf_str *);

int re_getlib(unsigned char *, unsigned char *, 
	      int, int, int, int, int, long *, long *, 
	      struct lmf_str *m_fptr);

void closelib(struct lmf_str *m_fptr);

void irand(int);
int nrand(int);

extern int ann_scan(unsigned char *, int, unsigned char **, int);
extern int scanseq(unsigned char *seq, int n, char *str);
extern void re_ascii(int *qascii, int *sascii);
extern int recode(unsigned char *seq, int n, int *qascii, int nsq);
extern void revcomp(unsigned char *seq, int n, int *c_nt);

extern void init_ascii(int is_ext, int *sascii, int nsq, int is_dna);
extern void qshuffle(unsigned char *aa0, int n0, int nm0);
extern void free_pam2p(int **);

/* initialize environment (doinit.c) */
extern void initenv (int argc, char **argv, struct mngmsg *m_msg,
		     struct pstruct *ppst, unsigned char **aa0);

/* print timing information */
extern void ptime (FILE *, time_t);

#ifdef COMP_MLIB 
#define QGETLIB (q_file_p->getlib)
#endif

#define GETLIB (m_file_p->getlib)

void
save_best(struct buf_head *lib_buf, const struct mngmsg *, struct pstruct *ppst,
	  struct db_str *, FILE *fdata, int *, struct hist_str *, void **);

void
save_shuf(struct buf_head *lib_buf, int nitt, int shuff_max);

void preserve_seq(struct buf2_str *, struct seq_record *, struct beststr *);
void
buf_do_work(unsigned char **aa0, int n0, struct buf_head *lib_bhead_p,
	    struct pstruct *ppst, void **f_str);
void
buf_qshuf_work(unsigned char *aa0s, int n0, struct buf_head *lib_bhead_p,
	       struct pstruct *ppst, void *qf_str, int score_ix);
void
buf_shuf_work(unsigned char **aa0, int n0, unsigned char *aa1s, int zs_win,
	      struct buf_head *lib_bhead_p, struct pstruct *ppst,
	      void **f_str, int score_ix);

/* statistics functions */
extern int
process_hist(struct stat_str *sptr, int nstat, 
	     const struct mngmsg *m_msg,
	     struct pstruct *ppst,
	     struct hist_str *hist, void **, int);
extern void addhistz(double, struct hist_str *); /* scaleswn.c */
void selectbestz(struct beststr **, int, int );
extern double (*find_zp)(int score, double escore, int length, double comp,void *);

void last_stats(const unsigned char *, int, 
		struct stat_str *sptr, int nstats,
		struct beststr **bestp_arr, int nbest,
		const struct mngmsg *m_msg, struct pstruct *ppst, 
		struct hist_str *histp, void *);

int last_calc( unsigned char **aa0, unsigned char *aa1, int maxn,
	       struct beststr **bestp_arr, int nbest,
	       const struct mngmsg *m_msg, struct pstruct *ppst, 
	       void **f_str, void *rs_str);

void scale_scores(struct beststr **bestp_arr, int nbest,
		  struct db_str,struct pstruct *ppst, void *);

extern int shuffle(unsigned char *, unsigned char *, int);
extern int wshuffle(unsigned char *, unsigned char *, int, int);

extern void set_db_size(int, struct db_str *, struct hist_str *);

/* display functions */
extern void
showbest (FILE *fp, unsigned char **aa0, unsigned char *aa1, int maxn,
	  struct beststr **bestp_arr, int nbest,
	  int qlib, struct mngmsg *m_msg,struct pstruct *ppst,
	  struct db_str db, char **gstring2p, void **f_str);

extern void
showalign (FILE *fp, unsigned char **aa0, unsigned char *aa1, int maxn,
	   struct beststr **bestp_arr, int nbest, int qlib, 
	   const struct mngmsg *m_msg, const struct pstruct *ppst,
	   char **gstring2p, void **f_str);

/* misc functions */
void h_init(struct pstruct *, struct mngmsg *, char *);		/* doinit.c */
void last_init(struct mngmsg *, struct pstruct *); /* initfa/sw.c */
void last_params(unsigned char *, int, struct mngmsg *, struct pstruct *);
int validate_params(const unsigned char *, int, const struct mngmsg *,
		    const struct pstruct *,
		    const int *lascii, const int *pascii);

void s_abort(char *, char *);		/* compacc.c */

/* initfa/sw.c */
void resetp(struct mngmsg *, struct pstruct *); 

void gettitle(char *, char *, int);	/* nxgetaa.c */
void libchoice(char *lname, int, struct mngmsg *); /* lib_sel.c */
void libselect(char *lname, struct mngmsg *);	/* lib_sel.c */
void query_parm(struct mngmsg *, struct pstruct *); /* initfa/sw.c */
void selectbestz(struct beststr **, int, int);

/* compacc.c */
void prhist(FILE *, const struct mngmsg *, struct pstruct *, 
	    struct hist_str hist, int nstats, struct db_str, char *, char **, char **);
void printsum(FILE *, struct db_str db);
int reset_maxn(struct mngmsg *, int);	/* set m_msg.maxt, maxn from maxl */

FILE *outfd;			/* Output file */

/* this information is global for fsigint() */
extern time_t s_time();			/* fetches time */
time_t tstart, tscan, tprev, tdone;	/* Timing */
#ifdef COMP_MLIB
time_t ttscan, ttdisp;
#endif
time_t tdstart, tddone;

static struct db_str qtt = {0l, 0l, 0};

#ifdef COMP_THR
/***************************************/
/* thread global variable declarations */
/***************************************/

/* functions for getting/sending buffers to threads (thr_sub.c) */
extern void init_thr(int , struct thr_str *, const struct mngmsg *, struct pstruct *,
		     unsigned char *, int);
extern void start_thr(void);
extern void get_rbuf(struct buf_head **lib_buf, int max_wor_buf);
extern void put_rbuf(struct buf_head *lib_buf, int max_work_buf);
extern void wait_rbuf(int max_work_buf);
extern void rbuf_done(int nthreads);
extern void put_rbuf_done(int nthreads, struct buf_head *lib_buf, 
			  int max_work_buf);
#undef XTERNAL
#include "thr.h"
#endif

struct buf_head *lib_buf2_list;

/* these variables must be global for comp_thr.c so that savebest()
   can use them */
static struct beststr **bestp_arr;	/* array of pointers */
static int nbest;	/* number of best scores */

static struct stat_str *stats; /* array of scores for statistics from real
			     (or shuffled) sequences*/
static struct stat_str *qstats;	/* array of scores for shuffled query stats */
  struct stat_str *rstats;	/* array of scores from shuffled library */

  /* these variables are global so they can be set both by the main()
     program and savebest() in threaded mode.
  */
static int nstats, nqstats, nrstats, kstats;
static double zbestcut;		/* cut off for best z-score */
static int bestfull;		/* index for selectbest() */
static int stats_done=0;	/* flag for z-value processing */

void fsigint();

int
main (int argc, char *argv[]) 
{
  unsigned char *aa0[6], *aa0s, *aa1, *aa1ptr, *aa1a = NULL;
  unsigned char *aa1save;	/* aa1shuff and aa1save must be distinct */
  unsigned char *aa1shuff, *aa1shuff_b=NULL;	/* for new unthreaded version */
  int n1;
  int j;
  int shuff_mult, n1lib_req;
  struct beststr *bbp;

  int *n1tot_ptr=NULL, *n1tot_cur;
  int n1tot_cnt=0;
  int n1tot_v, aa1_loff;

  long loffset, l_off;	/* loffset is the coordinate of first residue
			   when lcont > 0; l_off is not used in the
			   main loop, only in showbest and showalign */
  struct a_res_str *next_ares, *old_ares; /* used to free-up old a_res */

  /* status/parameter information */
  char info_lib_range[MAX_FN];
  char *info_lib_range_p;
  char info_pgm_abbr[MAX_SSTR];
  char info_qlabel[MAX_FN];
  char *info_gstring2p[2];
  char info_gstring3[MAX_STR];
  char *info_hstring_p[2];
  char l_bline[MAX_SSTR];
#ifdef COMP_MLIB
  char q_bline[MAX_STR];
  fseek_t qseek;
  int qlib;
  struct lmf_str *q_file_p;
  int sstart, sstop, is;
#endif
  long next_q_offset;
  int lstart, lstop;
  int id;
  struct lmf_str *m_file_p;

  int t_best, t_rbest, t_qrbest;	/* best score of two/six frames */
  double t_escore, t_rescore, t_qrescore; /* best evalues of two/six frames */
  double db_tt;
  int i_score;

  struct pstruct pst;
  void *f_str[6], *qf_str;	/* different f_str[]'s for forward,reverse */
  int have_f_str=0;

  /* move buf_head, buf2_str outside COMP_THR to be used consistently for everything */
  struct buf_head *lib_bhead_p, *t_lib_bhead_p;

  struct buf2_str *lib_buf2_p;
  struct seq_record *current_seq_p;
  int buf2_cnt;

  long ntbuff;			/* current length (in residues/bytes of buffer */
  int max_buf2_cnt;		/* number of entries per buffer for threads */
  int max_do_cnt;		/* number of entries to use for shuffles/opts-aligns */
  int ave_seq_len, buf_siz;
  int max_work_buf;
  int filled_worker_bufs;
#ifdef COMP_THR
  int empty_reader_bufs;
  int t_reader_buf_readp;
  struct thr_str *work_info;
#endif

  struct mngmsg m_msg;		/* Message from host to manager */
  int iln, itt;			/* index into library names */
  char rline[MAX_FN];
  char argv_line[MAX_STR];
  int t_quiet;

  struct rstruct  rst;		/* results structure */
  struct rstruct  rrst;		/* results structure for shuffle*/
  int i;

  FILE *fdata=NULL;		/* file for full results */
  fseek_t lseek;		/* seek into library of current sequence */
  char libstr[MAX_UID];		/* required because getlib() does not return
				   lcont > 0 */
  struct beststr *best;		/* array of best scores */

  /* save sequence meta info for sequences that are not currently available */
  struct seq_record *best_seqs;

  int n_libstr;			/* length of libstr */
  int jstats;
  int leng;			/* leng is length of the descriptive line */
  int maxn;			/* size of the library sequence examined */
  int maxl;			/* size of library buffer */
  int qlcont;			/* continued query sequence */
  int lcont, ocont, maxt;	/* continued sequence */
  int igncnt=0;			/* count for ignoring sequences warning */
  double zscore;		/* tmp value */
  char *bp;			/* general purpose string ptr */
  
  /* this is necessary because of an SGI Irix 64 issue */
  info_gstring2p[0] = calloc(MAX_STR,sizeof(char));
  info_gstring2p[1] = calloc(MAX_STR,sizeof(char));
  info_hstring_p[0] = calloc(MAX_STR,sizeof(char));
  info_hstring_p[1] = calloc(MAX_STR,sizeof(char));

  /* Initialization */

#if defined(UNIX)
  m_msg.quiet= !isatty(1);
#else
  m_msg.quiet = 0;
#endif

#ifdef PGM_DOC
  argv_line[0]='\0';
  for (i=0; i<argc; i++) {
    strncat(argv_line," ",sizeof(argv_line)-strlen(argv_line)-1);
    if (strchr(argv[i],' ')) {
      strncat(argv_line,"\"",sizeof(argv_line)-strlen(argv_line)-1);
      strncat(argv_line,argv[i],sizeof(argv_line)-strlen(argv_line)-1);
      strncat(argv_line,"\"",sizeof(argv_line)-strlen(argv_line)-1);
    }
    else {
      strncat(argv_line,argv[i],sizeof(argv_line)-strlen(argv_line)-1);
    }
  }
  argv_line[sizeof(argv_line)-1]='\0';

  fprintf(stdout, "#%s\n",argv_line);
#endif

  /* first initialization routine - nothing is known */
  h_init(&pst, &m_msg, info_pgm_abbr);
  
  m_msg.db.length = m_msg.ldb.length = qtt.length = 0l;
  m_msg.db.entries = m_msg.db.carry = 
    m_msg.ldb.entries = m_msg.ldb.carry = qtt.entries = qtt.carry = 0;
  m_msg.pstat_void = NULL;
  m_msg.hist.entries = 0;

  for (iln=0; iln<MAX_LF; iln++) m_msg.lb_mfd[iln]=NULL;

  f_str[5] = f_str[4] = f_str[3] = f_str[2] = f_str[1] = f_str[0] = NULL;

  aa0[0] = NULL;
  /* second initialiation - get commmand line arguments */
  initenv (argc, argv, &m_msg, &pst,&aa0[0]);

  if (m_msg.markx & MX_M11OUT) {
    fprintf(stdout, "#:lav\n\nd {\n   \"%s\"\n}",argv_line+1);
  }

#ifdef COMP_THR
  if ((work_info=
       (struct thr_str *)calloc(fa_max_workers,sizeof(struct thr_str)))==NULL) {
    fprintf(stderr, " cannot allocate work_info[%d]\n",fa_max_workers);
    exit(1);
  }
#else
  fa_max_workers = 1;
#endif

  tstart = tscan = s_time();
  tdstart = time(NULL);

  /* Allocate space for the query and library sequences */
  /* pad aa0[] with an extra SEQ_PAD chars for ALTIVEC padding */
  if (aa0[0]==NULL) {
    if ((aa0[0] = (unsigned char *)malloc((m_msg.max_tot+1+SEQ_PAD)*sizeof(unsigned char)))
	== NULL)
      s_abort ("Unable to allocate query sequence", "");
    *aa0[0]=0;
    aa0[0]++;
  }
  aa0[5]=aa0[4]=aa0[3]=aa0[2]=aa0[1]=aa0[0];

  if ((aa1save = (unsigned char *)malloc((m_msg.max_tot+1)*sizeof (char))) == NULL) {
    s_abort ("Unable to allocate library overlap", "");
  }
  *aa1save=0;
  aa1save++;

  if (m_msg.markx & MX_HTML) {
#ifdef HTML_HEAD    
    fprintf(stdout,"<html>\n<head>\n<title>%s Results</title>\n</head>\n<body>\n",prog_func);
#endif
    fprintf(stdout,"<pre>\n");
  }

  if (m_msg.std_output) {
    fprintf(stdout,"%s\n",iprompt0);
    fprintf(stdout," %s%s\n",verstr,refstr);
  }
  if (m_msg.markx & MX_HTML) fputs("</pre>\n",stdout);

  /* get query library name if not in argv[1] */
  if (m_msg.tname[0] == '\0') {
      if (m_msg.quiet == 1)
	s_abort("Query sequence undefined","");
    l1:	fputs (iprompt1, stdout);
      fflush  (stdout);
      if (fgets (m_msg.tname, MAX_FN, stdin) == NULL)
	s_abort ("Unable to read query library name","");
      m_msg.tname[MAX_FN-1]='\0';
      if ((bp=strchr(m_msg.tname,'\n'))!=NULL) *bp='\0';
      if (m_msg.tname[0] == '\0') goto l1;
  }

  /* Open query library */
  if ((q_file_p= openlib(m_msg.tname, m_msg.qdnaseq,qascii,!m_msg.quiet,NULL))==NULL) {
    s_abort(" cannot open library ",m_msg.tname);
  }

  /* Fetch first sequence */
  qlib = 0;
  m_msg.q_offset = next_q_offset = 0l;
  qlcont = 0;
  m_msg.n0 = 
    QGETLIB (aa0[0], MAXTST, m_msg.qtitle, sizeof(m_msg.qtitle),
	     &qseek, &qlcont,q_file_p,&m_msg.q_off);
  if ((bp=strchr(m_msg.qtitle,' '))!=NULL) *bp='\0';
  strncpy(info_qlabel,m_msg.qtitle,sizeof(info_qlabel));
  if (bp != NULL) *bp = ' ';
  info_qlabel[sizeof(info_qlabel)-1]='\0';

  /* if annotations are included in sequence, remove them */
  if (m_msg.ann_flg) {
    m_msg.n0 = ann_scan(aa0[0],m_msg.n0,&m_msg.aa0a,m_msg.qdnaseq);
  }

  /* if protein and term_code set, add '*' if not there */
  if (m_msg.term_code && !(m_msg.qdnaseq==SEQT_DNA || m_msg.qdnaseq==SEQT_RNA) &&
      aa0[0][m_msg.n0-1]!='*') {
    aa0[0][m_msg.n0++]='*';
    aa0[0][m_msg.n0]=0;
  }

  /* check for subset */
  if (q_file_p->opt_text[0]!='\0') {
    if (q_file_p->opt_text[0]=='-') {
      sstart=0; sscanf(&q_file_p->opt_text[1],"%d",&sstop);
    }
    else {
      sscanf(&q_file_p->opt_text[0],"%d-%d",&sstart,&sstop);
      sstart--;
      if (sstop <= 0 ) sstop = BIGNUM;
    }

    for (id=0,is=sstart; is<min(m_msg.n0,sstop); ) {
      aa0[0][id++]=aa0[0][is++];
    }
    aa0[0][id]=0;
    m_msg.n0 = min(m_msg.n0,sstop)-sstart;
    m_msg.q_off += sstart;
  }

  /* do this all the time */
  /* #if defined(SW_ALTIVEC) || defined(SW_SSE2) */
  /* for ALTIVEC, must pad with SEQ_PAD NULL's */
  for (id=0; id<SEQ_PAD; id++) {aa0[0][m_msg.n0+id]=0;}
  /* #endif */

  /* check to see if query has been segmented */
  if (qlcont) {
    next_q_offset = m_msg.q_offset + m_msg.n0 - m_msg.q_overlap;
  }
  else {
    next_q_offset = 0l;
  }

  /* this probably cannot happen any more */
  if (m_msg.n0 > MAXTST) {
    fprintf(stderr," sequence truncated to %d\n %s\n",MAXTST,m_msg.sqnam);
    if (m_msg.std_output) 
      fprintf(stdout," sequence truncated to %d\n %s\n",MAXTST,m_msg.sqnam);
    aa0[0][MAXTST]='\0';
    m_msg.n0=MAXTST;
  }

  /* check for protein/DNA alphabet type */
  if (m_msg.qdnaseq == SEQT_UNK) {
  /* cannot change the alphabet mapping if a matrix has been set */
  /* do automatic sequence recognition,but only for sequences > 20 residues */
    if ( !pst.pam_set && m_msg.n0 > 20 &&
	(float)scanseq(aa0[0],m_msg.n0,"ACGTUNacgtun")/(float)m_msg.n0 >0.85) {
      pascii = nascii;
      m_msg.qdnaseq = SEQT_DNA;
    }
    else {	/* its protein */
      pascii = aascii;
      m_msg.qdnaseq = SEQT_PROT;
    }
    /* modify qascii to use encoded version 
       cannot use memcpy() because it loses annotations 
    */
    re_ascii(qascii,pascii);
    init_ascii(pst.ext_sq_set,qascii,pst.nsq,m_msg.qdnaseq);
    m_msg.n0 = recode(aa0[0],m_msg.n0,qascii, pst.nsqx);
  }

  /* check sequence length -- cannot do before now because query alphabet may change */
  if (m_msg.n0 <= 0)
    s_abort ("Query sequence length <= 0: ", m_msg.tname);

#ifdef SUPERFAMNUM
  m_msg.nqsfnum = nsfnum;
  for (i=0; i <= nsfnum & i<10; i++) m_msg.qsfnum[i] = sfnum[i];
  m_msg.nqsfnum_n = nsfnum_n;
  for (i=0; i <= nsfnum_n & i<10; i++) m_msg.qsfnum_n[i] = sfnum_n[i];
#endif

  /* reset algorithm parameters for alphabet */
  resetp (&m_msg, &pst);

#ifndef COMP_MLIB
  gettitle(m_msg.tname,m_msg.qtitle,sizeof(m_msg.qtitle));
  if (m_msg.tname[0]=='-' || m_msg.tname[0]=='@') {
    strncmp(m_msg.tname,m_msg.qtitle,sizeof(m_msg.tname));
    if ((bp=strchr(m_msg.tname,' '))!=NULL) *bp='\0';
  }
#endif

  /* get library file names from argv[2] or by prompting */
  if (strlen (m_msg.lname) == 0) {
    if (m_msg.quiet == 1) s_abort("Library name undefined","");
    libchoice(m_msg.lname,sizeof(m_msg.lname),&m_msg);
  }
  
  /* map library abbreviation to file name(s) */
  libselect(m_msg.lname, &m_msg);

  /* Get additional parameters here */
  if (!m_msg.quiet) query_parm (&m_msg, &pst);
  
  last_init(&m_msg, &pst);

  /* Allocate space for saved scores */
  if ((best = 
       (struct beststr *)calloc((MAX_BEST+1),sizeof(struct beststr)))==NULL)
    s_abort("Cannot allocate best struct","");
  if ((bestp_arr = 
       (struct beststr **)malloc((MAX_BEST+1)*sizeof(struct beststr *)))==NULL)
    s_abort("Cannot allocate bestp_arr","");

  /* initialize high score boundary */
  bestp_arr[0] = &best[0];
  best[0].rst.score[0]=best[0].rst.score[1]=best[0].rst.score[2]= INT_MAX;
  best[0].rst.escore=FLT_MIN;	/* for E()-values, lower is best */
  best[0].zscore=FLT_MAX;	/* for Z-scores, bigger is best */
  best++; bestp_arr++;

  /* save best score sequence info */
  if ((best_seqs =
       (struct seq_record *)calloc((MAX_BEST+1),sizeof(struct seq_record)))==NULL)
    s_abort("Cannot allocate best_seqs","");

  /* allocate space for sampled scores */
  if ((stats =
       (struct stat_str *)calloc(MAX_STATS,sizeof(struct stat_str)))==NULL)
    s_abort ("Cannot allocate stats struct","");

  /* allocate space for shuffled library scores */
  if ((rstats =
       (struct stat_str *)calloc(m_msg.shuff_max,sizeof(struct stat_str)))==NULL)
    s_abort ("Cannot allocate rstats struct","");

#ifdef UNIX
  /* set up signals now that input is done */
  signal(SIGHUP,SIG_IGN);
#endif

  /* **************************************************************** */
  /* begin setting things up for threads */
  /* **************************************************************** */
  /* 
     This section defines max_buf2_cnt, the average number of entries per buffer, 
     and max_work_buf, the total number of buffers

     Use a 2 Mbyte (DEF_WORKER_BUF) buffer for each worker.  For
     proteins, that means 5,000 sequences of length 400 (average).
     For DNA, that means 2,000 sequences of length 1000.

     To accommodate larger libraries in memory, use more buffers, not
     bigger buffers.
  */

  if (m_msg.ldnaseq== SEQT_DNA) {
    ave_seq_len = AVE_NT_LEN;
    max_buf2_cnt = DEF_WORKER_BUF/AVE_NT_LEN;
  }
  else {
    ave_seq_len = AVE_AA_LEN;
    max_buf2_cnt = DEF_WORKER_BUF/AVE_AA_LEN;
  }

  /* however - buffer sizes should be a function of the number of
     workers so that all the workers are kept busy.  Assuming a 10,000
     entry library is the smallest we want to schedule, then

  if (max_buf2_cnt > 10000/fa_max_workers) 
    max_buf2_cnt = 10000/(2*fa_max_workers);

  */

  max_buf2_cnt /= m_msg.thr_fact;

  /* finally, max_work_buf should be mod 6 for tfasta/s/f */
  max_buf2_cnt -= (max_buf2_cnt % 6);

  /* max_work_buf is the number of buffers - if the worker buffers are
     small, then make lots more buffers */

  max_work_buf = (DEF_WORKER_BUF * 2 * fa_max_workers)/(ave_seq_len * max_buf2_cnt);
  if (max_work_buf < 2*fa_max_workers) max_work_buf = 2*fa_max_workers;

  max_work_buf -= (max_work_buf%fa_max_workers);

#ifndef COMP_THR
  /* if not threaded, only one (larger) buffer */
  max_buf2_cnt *= max_work_buf;
  max_work_buf = 1;
#endif

  /* first allocate space for buffer headers */
  if ((lib_buf2_list = (struct buf_head *)calloc((size_t)(max_work_buf),sizeof(struct buf_head))) == NULL) {
    fprintf(stderr," cannot allocate lib_buf2_list[%d]\n",max_work_buf);
    exit(1);
  }

#ifdef COMP_THR
  if ((worker_buf = (struct buf_head **)calloc((size_t)(max_work_buf),sizeof(struct buf_head *))) == NULL) {
    fprintf(stderr," cannot allocate **worker_buf[%d]\n",max_work_buf);
    exit(1);
  }

  if ((reader_buf = (struct buf_head **)calloc((size_t)(max_work_buf),sizeof(struct buf_head *))) == NULL) {
    fprintf(stderr," cannot allocate **reader_buf[%d]\n",max_work_buf);
    exit(1);
  }
#endif

  /* allocate space for library buffers and results */

  /* there are four structures/buffers used to keep track of
     sequences/results:

     (1) lib_buf2_list[] is a bhead_str array, which simply stores
         whether the results are ready and the number of results
         available.

     (2) lib_buf2_list[].buf2 is a buf2_str, which points to a
         sequence, and records frame and scores
  
     (3) lib_buf2_list[].buf2[0]->seq points to a seq_record, which
         contains information about a sequence, including an entry
         into the sequence buffer described next.

     (4) lib_buf2_list[].start points to the start of the sequence
         buffer.

  */

  buf_siz = max(max_buf2_cnt*ave_seq_len, m_msg.max_tot * 4);
  if (buf_siz < m_msg.max_tot) buf_siz = m_msg.max_tot;

  for (i=0; i<max_work_buf; i++) {
    /* allocate max_buf2_cnt buf2_str's into each buf2 */
    if ((lib_buf2_list[i].buf2 =(struct buf2_str *)calloc((size_t)(max_buf2_cnt+1),
						   sizeof(struct buf2_str)))
        ==NULL) {
      fprintf(stderr," cannot allocate buffer struct %d %d\n",i,max_buf2_cnt+1);
      exit(1);
    }

    if ((lib_buf2_list[i].buf2[0].seq_b = 
	 lib_buf2_list[i].buf2[0].seq =
	 (struct seq_record *)calloc((size_t)(max_buf2_cnt+1),
				     sizeof(struct seq_record)))
        ==NULL) {
      fprintf(stderr," cannot allocate seq_record buffer %d %d\n",i,max_buf2_cnt+1);
      exit(1);
    }

    if ((lib_buf2_list[i].start =
         (unsigned char *)calloc((size_t)(buf_siz),sizeof(unsigned char)))
        ==NULL) {
      fprintf(stderr," cannot allocate buffer %d: %d\n",i, buf_siz);
      exit(1);
    }

    /* make certain there is a '\0' at the beginning */
    lib_buf2_list[i].start++;
    lib_buf2_list[i].buf2[0].seq_b->aa1b = 
      lib_buf2_list[i].buf2[0].seq->aa1b = lib_buf2_list[i].start;

    /* make certain there is a '\0' at the beginning */

    lib_buf2_list[i].have_results=0;

#ifdef COMP_THR
    reader_buf[i] = &lib_buf2_list[i];
#endif
  }

  /* initialization of global variables for threads/buffers */

#if defined(COMP_THR)
#ifdef DEBUG
  fprintf(stderr," max_work_buf: %d\n", max_work_buf);
#endif
  num_reader_bufs = max_work_buf;
  filled_worker_bufs = empty_reader_bufs = 0;
  num_worker_bufs = 0;
  reader_done = 0;
  worker_buf_workp = 0;
  worker_buf_readp = 0;
  reader_buf_workp = 0;
  reader_buf_readp = 0;

  start_thread = 1;	/* keeps threads from starting */
#endif

  /* Label the output */
  if ((bp = (char *) strchr (m_msg.lname, ' ')) != NULL) *bp = '\0';
  if (m_msg.ltitle[0] == '\0') {
    strncpy(m_msg.ltitle,m_msg.lname,sizeof(m_msg.ltitle));
    m_msg.ltitle[sizeof(m_msg.ltitle)-1]='\0';
  }

  if (m_msg.dfile[0]) fdata=fopen(m_msg.dfile,"w");

#ifdef COMP_MLIB
  if (m_msg.std_output) {
    printf("Query: %s\n", m_msg.tname);
  /*   if (m_msg.nln > 0) printf("searching %s library\n\n",m_msg.lbnames[0]); */
  }
  
#endif

  /* main loop for doing a search, getting the next query */
  while(1) {

    /* Initialize bestp_arr */
    for (nbest = 0; nbest < MAX_BEST; nbest++)
      bestp_arr[nbest] = &best[nbest];
    nbest = 0;

    m_msg.db.length = 0l;
    m_msg.db.entries = m_msg.db.carry = 0;
    qlib++;
    stats_done = 0;

    maxl = m_msg.max_tot - m_msg.n0 -2;	/* maxn = max library sequence space */

    maxn = reset_maxn(&m_msg,maxl);
    pst.maxlen = maxn;

    outfd = stdout;  
    zbestcut = -FLT_MAX;
    nstats = 0;
    nrstats = 0;

    /* get the last parameters */
    last_params(aa0[0],m_msg.n0, &m_msg, &pst);

    if (!validate_params(aa0[0],m_msg.n0, &m_msg, &pst,
			 lascii, pascii)) {
      fprintf(stderr," *** ERROR *** validate_params() failed:\n -- %s\n", argv_line);
      exit(1);
    }

    /*
      if our function returns approximate E()-scores, we do not need to
      work with raw scores and later calculate z-scores.  When
      approx. E()-scores are calculated, we still need various
      statistics structures, but we can get them immediately.  In this
      case, find_zp() must produce a z_score (large positive is good)
      from an e_score.
    */

    if (m_msg.escore_flg) {
      pst.zsflag_f = process_hist(stats,nstats,&m_msg,&pst,
				  &m_msg.hist,&m_msg.pstat_void,0);
      stats_done=1;
    }

#ifndef COMP_THR
    if (m_msg.qshuffle) {
      if ((aa0s=(unsigned char *)calloc(m_msg.n0+2+SEQ_PAD,sizeof(char)))==NULL) {
	fprintf(stderr,"cannot allocate aa0s[%d]\n",m_msg.n0+2);
	exit(1);
      }
      *aa0s='\0';
      aa0s++;
      memcpy(aa0s,aa0[0],m_msg.n0);
      qshuffle(aa0s,m_msg.n0,m_msg.nm0);

#if defined(SW_ALTIVEC) || defined(SW_SSE2)
      /* for ALTIVEC, must pad with 16 NULL's */
      for (id=0; id<SEQ_PAD; id++) {aa0s[m_msg.n0+id]=0;}
#endif
    }

    /* storage for shuffled library sequences */
    if ((aa1shuff_b = aa1shuff = (unsigned char *)malloc((maxn+2)*sizeof (char))) == NULL) {
      s_abort ("Unable to allocate shuffled library sequence", "");
    }
    *aa1shuff=0;
    aa1shuff++;
    irand(0);

    /* previous versions of FASTA have stored the reverse complement in
       the same array as the forward query sequence.  This version
       changes that, by allocating separate space for the reverse complement,
       and thus reducing the demand for a large MAXLIB/MAXTRN for long queries
    */
    if (m_msg.qframe == 2) {
      if ((aa0[1]=(unsigned char *)calloc(m_msg.n0+2+SEQ_PAD,sizeof(char)))==NULL) {
	fprintf(stderr,"cannot allocate aa0[1][%d]\n",m_msg.n0+2);
	exit(1);
      }
      *aa0[1] = '\0';
      aa0[1]++;
      memcpy(aa0[1],aa0[0],m_msg.n0+1);
#if defined(SW_ALTIVEC) || defined(SW_SSE2)
      /* for ALTIVEC, must pad with 16 NULL's */
      for (id=0; id<SEQ_PAD; id++) {aa0[1][m_msg.n0+id]=0;}
#endif
      revcomp(aa0[1],m_msg.n0,&pst.c_nt[0]);
    }
    /* set aa1 for serial - threaded points aa1 to buffer */

    aa1 = aa0[0] + m_msg.n0+1;	/* modified now that aa0[1] is done separately */
    *aa1++='\0';
#else
    init_thr(fa_max_workers, work_info, &m_msg, &pst, aa0[0], max_work_buf);
#endif

    /* allocate space for shuffled query scores (if needed */
    if (m_msg.qshuffle && qstats==NULL) {
      if ((qstats =
	   (struct stat_str *)calloc(m_msg.shuff_max+1,sizeof(struct stat_str)))==NULL)
	s_abort ("Cannot allocate qstats struct","");
    }
    nqstats = 0;

    if (m_msg.markx & MX_HTML) fputs("<pre>\n",stdout);

    /* rline[] is a tmp string */
    if (m_msg.qdnaseq == SEQT_DNA || m_msg.qdnaseq == SEQT_RNA) {
      strncpy(rline,(m_msg.qframe==1)? " (forward-only)" : "\0",sizeof(rline));
      rline[sizeof(rline)-1]='\0';
    }
    else rline[0]='\0';

    leng = (int)strlen(m_msg.qtitle);
    if (leng > 50) leng -= 10;

    sprintf (&m_msg.qtitle[leng], " - %d %s", m_msg.n0, m_msg.sqnam);
    m_msg.seqnm = 0;


    if (m_msg.std_output) {
#ifdef COMP_MLIB
      printf("%3d>>>%s%s\n", qlib,
	     m_msg.qtitle,
	     (m_msg.revcomp ? " (reverse complement)" : rline));
#else
      printf("%.50s: %d %s%s\n %s\n",
	     info_qlabel, m_msg.n0, m_msg.sqnam,
	     (m_msg.revcomp ? " (reverse complement)" : rline));
#endif
      /* check for annotation */
      if (m_msg.ann_flg && m_msg.aa0a != NULL) {
	printf("Annotation: ");
	for (j=0; j<m_msg.n0; j++) {
	  if (m_msg.aa0a[j] && m_msg.ann_arr[m_msg.aa0a[j]] != ' ' ) {
	    printf("|%ld:%c%c",
		   j+m_msg.q_off,m_msg.ann_arr[m_msg.aa0a[j]],pst.sq[aa0[0][j]]);
	  }
	}
	printf("\n");
      }
      printf("Library: %.60s", m_msg.ltitle);
      fflush(stdout);
    }

    n_libstr=MAX_UID;

    fflush(outfd);

    tprev = s_time();
  
    if (fdata) fprintf(fdata,">>>%ld %3d\t%-50s\n",qtt.entries,m_msg.n0,m_msg.qtitle);

    qtt.length += m_msg.n0;
    qtt.entries++;

#ifndef COMP_THR
    /* initialize the comparison function, returning f_str */
    init_work (aa0[0], m_msg.n0, &pst, &f_str[0]);
    have_f_str=1;

    f_str[5] = f_str[4] = f_str[3] = f_str[2] = f_str[1] =  f_str[0];

    /* label library size limits */
    if (pst.n1_low > 0 && pst.n1_high < BIGNUM) 
      sprintf(info_lib_range," (range: %d-%d)",pst.n1_low,pst.n1_high);
    else if (pst.n1_low > 0) 
      sprintf(info_lib_range," (range: >%d)",pst.n1_low);
    else if (pst.n1_high < BIGNUM)
      sprintf(info_lib_range," (range: <%d)",pst.n1_high);
    else
      info_lib_range[0]='\0';
    info_lib_range[sizeof(info_lib_range)-1]='\0';
    info_lib_range_p = info_lib_range;

    if (m_msg.qframe == 2) {
      init_work ( aa0[1], m_msg.n0, &pst, &f_str[1]);
    }
    if (m_msg.qshuffle) {
      init_work ( aa0s, m_msg.n0, &pst, &qf_str);
    }

    lib_bhead_p = lib_buf2_list;	/* equivalent to un-threaded get_rbuf() */
#else	/* COMP_THR */
    start_thr();

    /* now open the library and start reading */
    /* get a buffer and fill it up */
    get_rbuf(&lib_bhead_p,max_work_buf);
    empty_reader_bufs++;
#endif

    lib_bhead_p->buf2_cnt = 0;
    lib_bhead_p->have_results = 0;
    lib_bhead_p->buf2_type=BUF2_DOWORK;
    if (pst.zsflag > 10) {
      lib_bhead_p->buf2_type |= BUF2_DOSHUF;
    }
    lib_buf2_p = lib_bhead_p->buf2;
    ntbuff = 0;

    /* open the library - start the search */

    for (iln = 0; iln < m_msg.nln; iln++) {
      if ((m_msg.lb_mfd[iln] = m_file_p=
	   openlib(m_msg.lbnames[iln], m_msg.ldnaseq, lascii, !m_msg.quiet, m_msg.lb_mfd[iln]))
	  ==NULL) {
	fprintf(stderr," cannot open library %s\n",m_msg.lbnames[iln]);
	continue;
      }

      loffset = 0l;
      lcont = 0;
      ocont = 0;
      n1tot_v = n1tot_cnt = 0;
      n1tot_cur = n1tot_ptr = NULL;

      /* get next buffer to read into */
      maxt = maxn;

      /* read sequence directly into buffer */

      current_seq_p = lib_buf2_p->seq;
      aa1ptr = aa1 = current_seq_p->aa1b;

      while ((n1=GETLIB(aa1ptr, maxt, libstr, n_libstr,
			&lseek,	&lcont, m_file_p,
			&(current_seq_p->l_off)))>=0) {

	/* check for subset on library*/
	if (m_file_p->opt_text[0]!='\0') {
	  if (m_file_p->opt_text[0]=='-') {
	    lstart=0; sscanf(&m_file_p->opt_text[1],"%d",&lstop);
	  }
	  else {
	    sscanf(&m_file_p->opt_text[0],"%d-%d",&lstart,&lstop);
	    lstart--;
	    if (lstop <= 0 ) lstop = BIGNUM;
	  }
	  for (id=0,is=lstart; is<min(n1,lstop); ) aa1[id++]=aa1[is++];
	  aa1[id]=0;
	  n1 = min(n1,lstop)-lstart;
	  current_seq_p->l_off += lstart;
	}

#ifdef DEBUG
	/* check for out of range sequence */
	for (id=0; id<n1; id++) {
	  if (aa1[id] > pst.nsq_e) {
	    fprintf(stderr," *** ERROR *** %s[%d] = %d > %d out of range\n",libstr, id, aa1[id], pst.nsq_e);
	    aa1[id] = 1;
	  }
	}
#endif

	current_seq_p->n1 = n1;
	current_seq_p->m_file_p = (void *)m_file_p;
	current_seq_p->cont = ocont+1;
	current_seq_p->l_offset = loffset;
	current_seq_p->lseek = lseek;

#ifdef SUPERFAMNUM
	current_seq_p->nsfnum = nsfnum;
	if ((current_seq_p->sfnum[0]=sfnum[0])>0 &&
	    (current_seq_p->sfnum[1]=sfnum[1])>0 &&
	    (current_seq_p->sfnum[2]=sfnum[2])>0 &&
	    (current_seq_p->sfnum[3]=sfnum[3])>0 &&
	    (current_seq_p->sfnum[4]=sfnum[4])>0 &&
	    (current_seq_p->sfnum[5]=sfnum[5])>0 &&
	    (current_seq_p->sfnum[6]=sfnum[6])>0 &&
	    (current_seq_p->sfnum[7]=sfnum[7])>0 &&
	    (current_seq_p->sfnum[8]=sfnum[8])>0 &&
	    (current_seq_p->sfnum[9]=sfnum[9])>0) ;
#endif

	if ((bp=strchr(libstr,' '))!=NULL) *bp='\0';
	strncpy(current_seq_p->libstr,libstr,MAX_UID);	/* get old libstr for lcont>0 */

	if (m_msg.term_code && !lcont &&
	    m_msg.ldnaseq==SEQT_PROT && aa1ptr[n1-1]!=m_msg.term_code) {
	  aa1ptr[n1++]=m_msg.term_code;
	  aa1ptr[n1]=0;
	}


#ifdef DEBUG
	if (n_libstr <= MAX_UID) {
	  if ((bp=strchr(current_seq_p->libstr,' '))!=NULL) *bp='\0';
	}

	if (aa1[-1]!='\0' || aa1ptr[n1]!='\0') {
	  fprintf(stderr,"%s: aa1[%d] at %ld:%lld  missing NULL boundaries: %d %d\n",
		  current_seq_p->libstr,n1,
		  m_msg.db.entries+1,current_seq_p->lseek,
		  aa1[-1],aa1ptr[n1]);
	}
#endif

	/* check for a continued sequence and provide a pointer to 
	   the n1_tot array if lcont || ocont */
	n1tot_v += n1;
	if (lcont && !ocont) {	/* get a new pointer */
	  if (n1tot_cnt <= 0) {
	    if ((n1tot_ptr=calloc(1000,sizeof(int)))==NULL) {
	      fprintf(stderr," cannot allocate n1tot_ptr\n");
	      exit(1);
	    }
	    else {n1tot_cnt=1000;}
	  }
	  n1tot_cnt--;
	  n1tot_cur = n1tot_ptr++;
	}
	current_seq_p->n1tot_p = n1tot_cur;

	/* skip based on size range */
	/*
	if (n1tot_v < m_msg.n1_low || n1tot_v > m_msg.n1_high) {
	  goto loop2;
	}
	*/

	m_msg.db.entries++;
	m_msg.db.length += n1;
	if (m_msg.db.length > LONG_MAX) {
	  m_msg.db.length -= LONG_MAX; m_msg.db.carry++;
	}

	/* don't count long sequences more than once */
	if (aa1!=aa1ptr) {	/* this is a continuation */
	  current_seq_p->n1 = n1 += m_msg.l_overlap;	/* corrected 28-June-2008 */
	  m_msg.db.entries--;
	}

#ifdef DEBUG
	/* This finds most reasons for core dumps */
	if (pst.debug_lib)
	  for (i=0; i<n1; i++) {
	    if (aa1[i]>pst.nsq || aa1[i] <= 0) {
	      fprintf(stderr,
		      "%s residue[%d/%d] %d range (%d) lcont/ocont: %d/%d\n%s\n",
		      current_seq_p->libstr,i,current_seq_p->n1,aa1[i],pst.nsq,lcont,ocont,aa1ptr+i);
	      aa1[i]=0;
	      n1=i-1;
	      break;
	    }
	  }
#endif

#ifdef PROGRESS
	if (!m_msg.quiet) 
	  if (m_msg.db.entries % 200 == 199) {
	    fputc('.',stderr);
	    if (m_msg.db.entries % 10000 == 9999) fputc('\n',stderr);
	    else if (m_msg.db.entries % 1000 == 999) fputc(' ',stderr);

	  }
#endif

	if (n1<=1) {
	  /*	if (igncnt++ <10)
		fprintf(stderr,"Ignoring: %s\n",current_seq_p->libstr);
	  */
	  goto loop2;
	}

	ntbuff += n1+1;
	for (itt=m_msg.revcomp; itt<=m_msg.nitt1; itt++) {
	  lib_buf2_p->frame = itt;
	  lib_buf2_p++;		/* point to next buf2 */
	  lib_bhead_p->buf2_cnt++;

	  /* point to the current sequence */
	  lib_buf2_p->seq = current_seq_p;
	} /* for (itt .. */

	if (lcont) {
	  memcpy(aa1save,&aa1[n1-m_msg.l_overlap],m_msg.l_overlap);
	}

	/* if the buffer is filled */
	if (lib_bhead_p->buf2_cnt >= max_buf2_cnt || ntbuff >= buf_siz - maxn) {

	  /* fprintf(stderr," new empty buffer at: %lld\n", current_seq_p->lseek); */

#ifdef COMP_THR		/* if COMP_THR - fill and empty buffers */
	  /* provide filled buffer to workers */
	  lib_bhead_p->have_data = 1;
	  put_rbuf(lib_bhead_p,max_work_buf);
	  filled_worker_bufs++;
	  empty_reader_bufs--;

	  /* get an empty buffer to fill */
	  get_rbuf(&lib_bhead_p,max_work_buf);
	  empty_reader_bufs++;
#else			/* just do the searches */
	  if (lib_bhead_p->buf2_type & BUF2_DOWORK) {
	    buf_do_work(aa0, m_msg.n0, lib_bhead_p, &pst, f_str);
	    if (m_msg.qshuffle)
	      buf_qshuf_work(aa0s,m_msg.n0, lib_bhead_p, &pst, qf_str, pst.score_ix);
	  }
	  if (lib_bhead_p->buf2_type & BUF2_DOSHUF) {
	    buf_shuf_work(aa0,m_msg.n0, aa1shuff, pst.zs_win,
			  lib_bhead_p, &pst, f_str, pst.score_ix);
	  }
#endif

	  /* "empty" buffers have results that must be processed */
	  if (lib_bhead_p->buf2_cnt && lib_bhead_p->have_results) {
#ifdef COMP_THR
	    filled_worker_bufs--;
#endif
	    save_best(lib_bhead_p,&m_msg, &pst, &m_msg.ldb, fdata,m_msg.qsfnum,
		      &m_msg.hist, &m_msg.pstat_void );

	    /* this section of code is only used for re-cycled buffers */
	    if (lib_bhead_p->have_best_save) {
	      lib_buf2_p = lib_bhead_p->buf2;
	      while (lib_bhead_p->buf2_cnt--) {
		if (lib_buf2_p->best_save != NULL) preserve_seq(lib_buf2_p, best_seqs, best);
		lib_buf2_p->best_save = NULL;
		lib_buf2_p++;
	      }
	      lib_bhead_p->have_best_save = 0;
	    }
	  }

	  /* now the buffer is truly empty, fill it up */
	  lib_bhead_p->buf2_cnt = 0;
	  lib_bhead_p->have_results = 0;
	  lib_bhead_p->buf2_type = BUF2_DOWORK;
	  if (pst.zsflag > 10) 	  lib_bhead_p->buf2_type |= BUF2_DOSHUF;

	  lib_buf2_p = lib_bhead_p->buf2;
	  current_seq_p = lib_buf2_p->seq = lib_bhead_p->buf2[0].seq_b;
	  aa1 = current_seq_p->aa1b;
	  ntbuff = 0;
	}
	else {	/* room left in current buffer, increment ptrs */
	  current_seq_p++;
	  lib_buf2_p->seq = current_seq_p;
	  aa1 = current_seq_p->aa1b = current_seq_p[-1].aa1b+n1+1;
	}

      loop2: 
	if (lcont) {
	  maxt = m_msg.maxt3;
	  memcpy(aa1,aa1save,m_msg.l_overlap);
	  aa1ptr= &aa1[m_msg.l_overlap];	/* aa1ptr is where the next GETLIB sequence goes */
					/* aa1 is the beginning of the sequence for do_work() */
	  loffset += n1 - m_msg.l_overlap;	/* this must be n1, which is the old value, not current_seq_p->n1 */ 
	  ocont = lcont;
	}
	else {
	  maxt = maxn;
	  aa1ptr=aa1;
	  if (ocont) *n1tot_cur = n1tot_v;
	  ocont = 0;
	  loffset = 0l;
	  n1tot_v = 0;
	  n1tot_cur = NULL;
	}
      } /* end while((n1=getlib())) */
    } /* end iln=1..nln */

      /* all done */

#ifdef COMP_THR
    /* check last buffers for any results */
    if (lib_bhead_p->buf2_cnt > 0) {
      lib_bhead_p->have_data = 1;
      put_rbuf(lib_bhead_p,max_work_buf);
      empty_reader_bufs--;
    }

    info_lib_range_p = work_info[0].info_lib_range;

    /* wait for the threads to finish */
    wait_rbuf(max_work_buf - empty_reader_bufs);

    for (i=0, t_reader_buf_readp = reader_buf_readp-1;
	 i < num_reader_bufs - empty_reader_bufs ; i++, t_reader_buf_readp--) {
      if (t_reader_buf_readp < 0) t_reader_buf_readp = max_work_buf - 1;
      t_lib_bhead_p = reader_buf[t_reader_buf_readp];
      /*
      fprintf(stderr," buf2_cnt[%d]: %d [data/res:%d/%d]\n",
	      t_reader_buf_readp,t_lib_bhead_p->buf2_cnt,
      	      t_lib_bhead_p->have_data, t_lib_bhead_p->have_results);
      */
      if (t_lib_bhead_p->buf2_cnt > 0 && t_lib_bhead_p->have_results) {
	save_best(t_lib_bhead_p,&m_msg, &pst, &m_msg.ldb, fdata,m_msg.qsfnum,
		  &m_msg.hist, &m_msg.pstat_void);
	t_lib_bhead_p->buf2_cnt = t_lib_bhead_p->have_results = 0;
      }
    }
    /*
    fprintf(stderr," reader_buf_readp: %d nbest: %d\n",
	    reader_buf_readp,nbest);
    */
#else
    if (lib_bhead_p->buf2_type & BUF2_DOWORK) {
      buf_do_work(aa0, m_msg.n0, lib_bhead_p, &pst, f_str);
      if (m_msg.qshuffle)
	buf_qshuf_work(aa0s,m_msg.n0, lib_bhead_p, &pst, qf_str, pst.score_ix);
    }

    if (lib_bhead_p->buf2_type & BUF2_DOSHUF) {
      buf_shuf_work(aa0, m_msg.n0, aa1shuff, pst.zs_win,
		    lib_bhead_p, &pst, f_str, pst.score_ix);
    }

    if (lib_bhead_p->buf2_cnt > 0) {
      save_best(lib_bhead_p,&m_msg, &pst, &m_msg.ldb, fdata,m_msg.qsfnum,
		&m_msg.hist, &m_msg.pstat_void);
    }
    lib_bhead_p->buf2_cnt = lib_bhead_p->have_results = 0;
#endif


#ifdef PROGRESS
    if (!m_msg.quiet)
      if (m_msg.db.entries >= 200) {fprintf(stderr," Done!\n");}
#endif

    m_msg.nbr_seq = m_msg.db.entries;
    get_param(&pst, info_gstring2p,info_gstring3);

    /* *************************** */
    /* analyze the last results    */
    /* *************************** */
    
#ifndef SAMP_STATS
    if (!stats_done && nstats > 0) {
#endif
      /* we ALWAYS do this if SAMP_STATS, because the statistics may have changed */

      pst.zsflag_f = process_hist(stats,nstats,&m_msg, &pst,&m_msg.hist,
				  &m_msg.pstat_void,stats_done);

      if (m_msg.pstat_void != NULL) {
	stats_done = 1;
	for (i = 0; i < nbest; i++) {
	  bestp_arr[i]->zscore =
	    (*find_zp)(bestp_arr[i]->rst.score[pst.score_ix],
		       bestp_arr[i]->rst.escore, bestp_arr[i]->seq->n1, 
		       bestp_arr[i]->rst.comp, m_msg.pstat_void);
	}
#ifndef SAMP_STATS
      }
      else pst.zsflag = -1;
#endif
    }

    /* if there are not many scores, produce better statistics by shuffling */
    /* but only if statistics are enabled (7-July-2008) */
    if (pst.zsflag > -1 && nbest > 0 && nbest < m_msg.shuff_max) {

      /* (1) get the sequences into a buffer - the sequence
             information is currently in the bestp_arr - find out how
             many we have, and how many we will need - the number to shuffle */

      /* figure out how much space we need */
      n1lib_req = 0;
      for (i = 0; i < nbest; i++) {
	if (bestp_arr[i]->seq->aa1b == NULL) {
	  n1lib_req += bestp_arr[i]->seq->n1+ 2;
	  /* 	  n1lib_req += m_msg.max_tot + 2; */
	  /* 	  fprintf(stderr, " %d needs %d %d\n",i,bestp_arr[i]->seq->n1,n1lib_req); */
	}
      }

#ifndef COMP_THR
      if (n1lib_req >= maxn) { /* we need new space, aa1shuff is too small */
	if ((aa1shuff = aa1shuff_b = 
	     (unsigned char *)realloc(aa1shuff_b, n1lib_req*sizeof(char)))==NULL) {
	  fprintf(stderr," *** cannot realloc aa1shuff[%d]\n",n1lib_req);
	  exit(1);
	}
	*aa1shuff = '\0';
	aa1shuff++;
      }
#else
      if ((aa1shuff = aa1shuff_b = 
	   (unsigned char *)calloc(n1lib_req,sizeof(char)))==NULL) {
	fprintf(stderr," *** cannot realloc aa1shuff[%d]\n",n1lib_req);
	exit(1);
      }
      *aa1shuff = '\0';
      aa1shuff++;
#endif

      shuff_mult = (m_msg.shuff_max / nbest) + 1;
      
#ifdef COMP_THR
      max_do_cnt = min(max_buf2_cnt,m_msg.shuff_max / (2 * fa_max_workers));
      /* we don't have a left over one, so we need one */
      if (empty_reader_bufs == 0) {
	get_rbuf(&lib_bhead_p,max_work_buf);
	empty_reader_bufs++;
      }
#else
      max_do_cnt = max_buf2_cnt;
      lib_bhead_p = lib_buf2_list;	/* equivalent to un-threaded get_rbuf() */
#endif
      lib_bhead_p->buf2_cnt = 0;
      lib_bhead_p->have_results = 0;
      lib_bhead_p->buf2_type=BUF2_DOSHUF;
      lib_buf2_p = lib_bhead_p->buf2;

      for (i = 0; i < nbest; i++) {
	bbp = bestp_arr[i];
	if (bbp->seq->aa1b == NULL) {
	  /* get the sequence */
	  (bbp->seq->m_file_p->ranlib)(l_bline, sizeof(l_bline),
				       bbp->seq->lseek,bbp->seq->libstr,bbp->seq->m_file_p);
	  n1 = re_getlib(aa1save,NULL, maxn,m_msg.maxt3,
			 m_msg.l_overlap,bbp->seq->cont,m_msg.term_code,
			 &loffset,&l_off,bbp->seq->m_file_p);

	  /* 	  fprintf(stderr, " %d gets %d %d\n",i,bestp_arr[i]->seq->n1,n1); */

	  memcpy(aa1shuff, aa1save, n1+1);
	  bbp->seq->aa1b = aa1shuff;
	  aa1shuff += n1 + 1;
	}

	for (j = 0; j < shuff_mult; j++ ) {
	  for (itt = m_msg.revcomp; itt <= m_msg.nitt1; itt++) {
	    /* this invalidates lib_buf2_p->seq */
	    lib_buf2_p->seq = bbp->seq;
	    lib_buf2_p->frame = itt;
	    lib_buf2_p++;		/* point to next buf2 */
	    lib_bhead_p->buf2_cnt++;
	  } /* for (itt .. */

	  if (lib_bhead_p->buf2_cnt >= max_do_cnt) {
	    /* (2) send sequences for shuffling */
#ifdef COMP_THR		/* if COMP_THR - fill and empty buffers */
	    /* provide empty buffer to workers */
	    lib_bhead_p->have_data = 1;
	    put_rbuf(lib_bhead_p,max_work_buf);
	    empty_reader_bufs--;

	    /* get an empty buffer to fill */
	    get_rbuf(&lib_bhead_p,max_work_buf);
	    empty_reader_bufs++;
#else			/* just do the searches */
	    if (lib_bhead_p->buf2_type & BUF2_DOSHUF) {
	      buf_shuf_work(aa0,m_msg.n0, aa1save, pst.zs_win,
			    lib_bhead_p, &pst, f_str, pst.score_ix);
	    }
#endif
	    /* (3) save results in the rstats structure */
	    if (lib_bhead_p->buf2_cnt > 0 && lib_bhead_p->have_results) {
	      save_shuf(lib_bhead_p,m_msg.nitt1,m_msg.shuff_max);
	    }

	    lib_bhead_p->buf2_cnt = 0;
	    lib_bhead_p->have_results = 0;
	    lib_bhead_p->buf2_type=BUF2_DOSHUF;
	    lib_buf2_p = lib_bhead_p->buf2;
	  }
	}
      }			/* done with bestp_arr[] */

#ifdef COMP_THR		/* if COMP_THR - fill and empty buffers */
      /* check last buffers for any results */
      if (lib_bhead_p->buf2_cnt > 0) {
	put_rbuf(lib_bhead_p,max_work_buf);
	empty_reader_bufs--;
      }
    
      /* wait for the threads to finish */

      wait_rbuf(max_work_buf - empty_reader_bufs);
      /*
      fprintf(stderr, " num_reader[%d]-empty[%d]: %d\tnrstats: %d\n",
	      num_reader_bufs,empty_reader_bufs,
	      num_reader_bufs-empty_reader_bufs, nrstats);
      */

    for (i=0, t_reader_buf_readp = reader_buf_readp-1;
	 i < num_reader_bufs - empty_reader_bufs ; i++, t_reader_buf_readp--) {
      if (t_reader_buf_readp < 0) t_reader_buf_readp = max_work_buf - 1;
      t_lib_bhead_p = reader_buf[t_reader_buf_readp];

      /*
      fprintf(stderr," buf2_cnt[%d]: %d [data/res:%d/%d]\n",
	      t_reader_buf_readp,t_lib_bhead_p->buf2_cnt,
      	      t_lib_bhead_p->have_data, t_lib_bhead_p->have_results);
      */

      if (t_lib_bhead_p->buf2_cnt > 0 && t_lib_bhead_p->have_results) {
	  save_shuf(t_lib_bhead_p,m_msg.nitt1, m_msg.shuff_max);
	  t_lib_bhead_p->buf2_cnt = t_lib_bhead_p->have_results = 0;
	}

      /*
	fprintf(stderr," nrstats: %d\n",nrstats);
      */
      }
#else	/* just do the searches */
      /* aa1save is used for shuffles, not aa1shuf, because aa1shuf
	 has library sequences */
      buf_shuf_work(aa0,m_msg.n0, aa1save, pst.zs_win,
		    lib_bhead_p, &pst, f_str, pst.score_ix);

      save_shuf(lib_bhead_p,m_msg.nitt1,m_msg.shuff_max);
      lib_bhead_p->buf2_cnt = lib_bhead_p->have_results = 0;
#endif
      /* (4) analyze rstats */
      if (pst.zsflag < 10) pst.zsflag += 10;
      pst.zsflag_f = process_hist(rstats,nrstats,&m_msg, &pst,&m_msg.hist,
				  &m_msg.pstat_void,0);
    }

    if (!pst.zdb_size_set) pst.zdb_size = m_msg.ldb.entries;

#ifdef COMP_THR
    /* before I call last_calc/showbest/showalign, I need init_work() to
       get an f_str. This duplicates some code above, which is used in
       the non-threaded version
    */

    if (!have_f_str) {
      init_work(aa0[0],m_msg.n0,&pst,&f_str[0]);
      have_f_str = 1;
      f_str[5] = f_str[4] = f_str[3] = f_str[2] = f_str[1] =  f_str[0];

      if (m_msg.qframe == 2) {
	if ((aa0[1]=(unsigned char *)calloc((size_t)m_msg.n0+2,
					    sizeof(unsigned char)))==NULL) {
	  fprintf(stderr," cannot allocate aa0[1][%d] for alignments\n",
		  m_msg.n0+2);
	}
	*aa0[1]='\0';
	aa0[1]++;
	memcpy(aa0[1],aa0[0],m_msg.n0+1);
	revcomp(aa0[1],m_msg.n0,&pst.c_nt[0]);
	init_work(aa0[1],m_msg.n0,&pst,&f_str[1]);
      }
    }
#endif

    /* now we have one set of scaled scores for in bestp_arr  -
       for FASTS/F, we need to do some additional processing */

    if (!m_msg.qshuffle) {
      last_stats(aa0[0], m_msg.n0, stats,nstats, bestp_arr,nbest,
		 &m_msg, &pst, &m_msg.hist, &m_msg.pstat_void);
    }
    else {
      last_stats(aa0[0], m_msg.n0,
		 qstats,nqstats, bestp_arr,nbest, &m_msg, &pst, 
		 &m_msg.hist, &m_msg.pstat_void);
    }

    /* here is a contradiction: if pst.zsflag < 0, then m_msg.pstat_void
       should be NULL; if it is not, then process_hist() has been called */
    if (pst.zsflag < 0 && m_msg.pstat_void != NULL) pst.zsflag = 1;

    if (m_msg.last_calc_flg) {
      /* last_calc may need coefficients from last_stats() */
      nbest = last_calc(aa0, aa1save, maxn, bestp_arr, nbest, &m_msg, &pst,
			f_str, m_msg.pstat_void);
    }

    /* in addition to scaling scores, this sorts bestp_arr[nbest] */
    scale_scores(bestp_arr,nbest,m_msg.db, &pst,m_msg.pstat_void);

    get_param(&pst, info_gstring2p,info_gstring3);

    /* **************************************************************** */
    /* label Library: output                                            */
    /* **************************************************************** */

    if (m_msg.std_output) {
      if (m_msg.db.carry==0) {
	fprintf(stdout, " %7ld residues in %5ld sequences\n", m_msg.db.length, m_msg.db.entries);
      }
      else {
	db_tt = (double)m_msg.db.carry*(double)LONG_MAX + (double)m_msg.db.length;
	fprintf(stdout, " %.0f residues in %5ld library sequences\n", db_tt, m_msg.db.entries);
      }

      prhist (stdout, &m_msg, &pst, m_msg.hist, nstats, m_msg.ldb, info_lib_range_p,
	      info_gstring2p, info_hstring_p);

      tscan = s_time();
      printf (" Scan time: ");
      ptime(stdout,tscan-tprev);
      printf ("\n");
    }
#ifdef COMP_MLIB
    ttscan += tscan-tprev;
#endif

  l3:
    if (!m_msg.quiet) {
      printf("Enter filename for results [%s]: ", m_msg.outfile);
      fflush(stdout);
    }

    rline[0]='\0';
    if (!m_msg.quiet && fgets(rline,sizeof(rline),stdin)==NULL) goto end_l;
    if ((bp=strchr(rline,'\n'))!=NULL) *bp = '\0';
    if (rline[0]!='\0') strncpy(m_msg.outfile,rline,sizeof(m_msg.outfile));
    if (m_msg.outfile[0]!='\0') {
      if (!m_msg.outfd) {
	if ((outfd=fopen(m_msg.outfile,"w"))==NULL) {
	  fprintf(stderr," could not open %s\n",m_msg.outfile);
	  if (!m_msg.quiet) goto l3;
	  else goto l4;
	}
	m_msg.outfd = outfd;
      }
      else {
	outfd = m_msg.outfd;
      }

#ifdef PGM_DOC
      fputs(argv_line,outfd);
      fputc('\n',outfd);
#endif  
      fputs(iprompt0,outfd);
      fprintf(outfd," %s%s\n",verstr,refstr);

      fprintf(outfd,"Query: %s%s, %d %s\n",
	      info_qlabel, (m_msg.revcomp ? "-" : "\0"), m_msg.n0,
	      m_msg.sqnam);

#ifdef COMP_MLIB
      fprintf(outfd,"%3d>>>%s - %d %s%s\n", qlib,
	   m_msg.qtitle, m_msg.n0, m_msg.sqnam,
	   (m_msg.revcomp ? " (reverse complement)" : rline));
#else
      fprintf(outfd,"%.50s: %d %s%s\n %s\n",
	   info_qlabel, m_msg.n0, m_msg.sqnam,
	   (m_msg.revcomp ? " (reverse complement)" : rline));
#endif
      fprintf(outfd, "Library: %.60s%s", m_msg.ltitle,info_lib_range_p);
      if (m_msg.db.carry==0) {
	fprintf(outfd, " %7ld residues in %5ld sequences\n", m_msg.db.length, m_msg.db.entries);
      }
      else {
	db_tt = (double)m_msg.db.carry*(double)LONG_MAX + (double)m_msg.db.length;
	fprintf(outfd, " %.0f residues in %5ld library sequences\n", db_tt, m_msg.db.entries);
      }

      prhist(outfd, &m_msg, &pst,m_msg.hist, nstats, m_msg.db, info_lib_range_p,
	     info_gstring2p, info_hstring_p);
    }

  l4:   
    if (m_msg.markx & MX_HTML) {
      fputs("</pre>\n<p>\n<hr>\n<p>\n",outfd);
    }

    /* code from p2_complib.c to pre-calculate -m 9 alignment info -
       requires -q with -m 9 */

    if (m_msg.quiet || m_msg.markx & MX_M9SUMM) {

      /* to determine how many sequences to re-align (either for
	 do_opt() or calc_id() we need to modify m_msg.mshow to get
	 the correct number of alignments */

      if (m_msg.mshow_flg != 1 && pst.zsflag >= 0) {
	for (i=0; i<nbest && bestp_arr[i]->rst.escore< m_msg.e_cut; i++) {}
	m_msg.mshow = i;
      }

      if (m_msg.mshow <= 0) { /* no results to display */
	fprintf(outfd,"!! No sequences with E() < %f\n",m_msg.e_cut);
	m_msg.nshow = 0;
	goto end_l;
      }
    }

    if (m_msg.do_showbest) {
      showbest (stdout, aa0, aa1save, maxn, bestp_arr, nbest, qtt.entries, &m_msg, &pst,
		m_msg.db, info_gstring2p, f_str);

      if (outfd != stdout) {
	t_quiet = m_msg.quiet;
	m_msg.quiet = -1;	/* should guarantee 1..nbest shown */
	showbest (outfd, aa0, aa1save, maxn, bestp_arr, nbest, qtt.entries, &m_msg, &pst,
		  m_msg.db, info_gstring2p, f_str);
	m_msg.quiet = t_quiet;
      }
    }

    if (m_msg.nshow > 0) {
      rline[0]='N';
      if (!m_msg.quiet){
	printf(" Display alignments also? (y/n) [n] "); fflush(stdout);
	if (fgets(rline,sizeof(rline),stdin)==NULL) goto end_l;
      }
      else rline[0]='Y';

      if (toupper((int)rline[0])=='Y') {
	if (!m_msg.quiet && m_msg.do_showbest) {
	  printf(" number of alignments [%d]? ",m_msg.nshow);
	  fflush(stdout);
	  if (fgets(rline,sizeof(rline),stdin)==NULL) goto end_l;
	  if (rline[0]!=0) sscanf(rline,"%d",&m_msg.nshow);
	  m_msg.ashow=m_msg.nshow;
	}

	if (m_msg.markx & (MX_AMAP+ MX_HTML + MX_M9SUMM) && !(m_msg.markx & MX_M10FORM)) {
	  fprintf(outfd,"\n>>>%s%s, %d %s vs %s library\n",
		  info_qlabel,(m_msg.revcomp ? "_rev":"\0"), m_msg.n0,
		  m_msg.sqnam,m_msg.lname);
	}

	if (m_msg.markx & MX_M10FORM) {
	  fprintf(outfd,"\n>>>%s%s, %d %s vs %s library\n",
		  info_qlabel,(m_msg.revcomp ? "-":"\0"), m_msg.n0, m_msg.sqnam,
		  m_msg.lname);
	  fprintf(outfd,"; pg_name: %s\n",argv[0]);
	  fprintf(outfd,"; pg_ver: %s\n",mp_verstr);
	  fprintf(outfd,"; pg_argv:");
	  for (i=0; i<argc; i++)
	    fprintf(outfd," %s",argv[i]);
	  fputc('\n',outfd);
	  fputs(info_gstring3,outfd);
	  fputs(info_hstring_p[0],outfd);
	  fputs(info_hstring_p[1],outfd);
	}

	showalign (outfd, aa0, aa1save, maxn, bestp_arr, nbest, qtt.entries,
		   &m_msg, &pst, info_gstring2p, f_str);

	fflush(outfd);
      }
    }

  end_l:

    if (fdata) {
      fprintf(fdata,"/** Algorithm : %s  **/\n",info_gstring2p[0]);
      fprintf(fdata,"/** Parameters : %s  **/\n",info_gstring2p[1]);
      fprintf(fdata,"%3ld%-50s\n",qtt.entries-1,m_msg.qtitle);
      fflush(fdata);
    }
    
#if defined(COMP_THR)
    rbuf_done(fa_max_workers);

    for (i=0; i<max_work_buf; i++) {
      lib_buf2_list[i].buf2_cnt=
	lib_buf2_list[i].have_results=
	lib_buf2_list[i].have_best_save = 0;
    }

    num_worker_bufs = 0;
    num_reader_bufs = max_work_buf;
    reader_done = 0;
    reader_wait = 1;
    worker_buf_workp = 0;
    worker_buf_readp = 0;
    reader_buf_workp = 0;
    reader_buf_readp = 0;

    start_thread = 1;	/* stop thread from starting again */
    
#endif
    /* clean up best_seqs */
    memset(best_seqs,0,(MAX_BEST+1)*sizeof(struct seq_record));

    /* re-initialize lib_buf2_list buffers */
    for (lib_bhead_p = lib_buf2_list; 
	 lib_bhead_p < lib_buf2_list+max_work_buf; ) {

      /* memset(lib_bhead_p->start-1,0,(size_t)buf_siz); */
      current_seq_p = lib_bhead_p->buf2[0].seq_b;
      /* this wipes out lib_bhead_p->buf2[0].seq_b */
      memset(lib_bhead_p->buf2,0,(size_t)(max_buf2_cnt+1)*sizeof(struct buf2_str));
      /* replace it */
      lib_bhead_p->buf2[0].seq_b = lib_bhead_p->buf2[0].seq = current_seq_p;
      lib_bhead_p->buf2[0].seq->aa1b = lib_bhead_p->start;
      lib_bhead_p++->have_results = 0;
    }

    /* re-initialize library counts */
    m_msg.ldb.length = 0l;
    m_msg.ldb.entries = m_msg.ldb.carry = 0;

    /* clean up alignment encodings */
    for (i=0; i < m_msg.nshow; i++) {
      if (bestp_arr[i]->have_ares & 0x2) {
	next_ares = bestp_arr[i]->a_res.next;
	free(bestp_arr[i]->a_res.res);
	bestp_arr[i]->a_res.nres = 0;
	bestp_arr[i]->a_res.res = NULL;
	while (next_ares != NULL) {
	  old_ares = next_ares;
	  next_ares = next_ares->next;
	  free(old_ares->res);
	  free(old_ares);
	}
      }
      bestp_arr[i]->have_ares = 0;
    }

#ifdef DEBUG
    /* check to see if there are ANY un-reset have_ares */
    for (i=0; i< nbest; i++) {
      if (bestp_arr[i]->have_ares) {
	fprintf(stderr," Un-reset have_ares[%d]: %d\n",i,bestp_arr[i]->have_ares);
	bestp_arr[i]->have_ares = 0;
      }
    }
#endif

    if (m_msg.qframe == 2) free(aa0[1]-1);

    if (have_f_str) {
      if (f_str[1]!=f_str[0]) {
	close_work (aa0[1], m_msg.n0, &pst, &f_str[1]);
      }
      close_work (aa0[0], m_msg.n0, &pst, &f_str[0]);
      have_f_str = 0;
#ifndef COMP_THR
      if (m_msg.qshuffle) close_work (aa0s, m_msg.n0, &pst, &qf_str);
#endif
      if (pst.pam_pssm) {
	free_pam2p(pst.pam2p[0]);
	free_pam2p(pst.pam2p[1]);
      }
    }

    if (aa1shuff_b != NULL) {
      free(aa1shuff_b);
      aa1shuff_b = NULL;
    }

    for (iln=0; iln < m_msg.nln; iln++) {
      if (m_msg.lb_mfd[iln]!=NULL) {
	closelib(m_msg.lb_mfd[iln]);
      }
    }

    tddone = time(NULL);
    tdone = s_time();
    fflush(outfd);

    ttdisp += tdone-tscan;

    /* reset pst parameters to original */
    pst.zsflag = m_msg.zsflag;
    pst.n1_low = m_msg.n1_low;
    pst.n1_high = m_msg.n1_high;

    /*     maxn = m_msg.max_tot; */

    m_msg.q_offset = next_q_offset;

    m_msg.n0 = 
      QGETLIB (aa0[0], MAXTST, m_msg.qtitle, sizeof(m_msg.qtitle),
	       &qseek, &qlcont,q_file_p,&m_msg.q_off);
    if (m_msg.n0 <= 0) break;
    if ((bp=strchr(m_msg.qtitle,' '))!=NULL) *bp='\0';
    strncpy(info_qlabel, m_msg.qtitle,sizeof(info_qlabel));
    if (bp != NULL) *bp=' ';
    info_qlabel[sizeof(info_qlabel)-1]='\0';

#ifdef SUPERFAMNUM
    m_msg.nqsfnum = nsfnum;
    for (i=0; i <= nsfnum & i<10; i++) m_msg.qsfnum[i] = sfnum[i];
    m_msg.nqsfnum_n = nsfnum_n;
    for (i=0; i <= nsfnum_n & i<10; i++) m_msg.qsfnum_n[i] = sfnum_n[i];
#endif

    if (m_msg.ann_flg) {
      m_msg.n0 = ann_scan(aa0[0],m_msg.n0,&m_msg.aa0a,m_msg.qdnaseq);
    }

    if (m_msg.term_code && m_msg.qdnaseq==SEQT_PROT &&
	aa0[0][m_msg.n0-1]!=m_msg.term_code) {
      aa0[0][m_msg.n0++]=m_msg.term_code;
      aa0[0][m_msg.n0]=0;
    }

#if defined(SW_ALTIVEC) || defined(SW_SSE2)
    /* for ALTIVEC, must pad with 16 NULL's */
    for (id=0; id<SEQ_PAD; id++) {aa0[0][m_msg.n0+id]=0;}
#endif

    if (m_msg.outfd) {fputc('\n',stdout);}
  }	/* end of while(1) for multiple queries */

  if (m_msg.markx & MX_M10FORM)
    fprintf(outfd,">>><<<\n");

  tdone = s_time();
  if ( m_msg.markx & MX_HTML) fputs("<p><pre>\n",outfd); 
  if (m_msg.std_output) {
    printsum(outfd, m_msg.db);
  }
  if ( m_msg.markx & MX_HTML) fputs("</pre>\n",outfd);
#ifdef HTML_HEAD
  if (m_msg.markx & MX_HTML) fprintf(outfd,"</body>\n</html>\n");
#endif
  if (m_msg.std_output && outfd!=stdout) printsum(stdout,m_msg.db);

  exit(0);
}   /* End of main program */

void
printsum(FILE *fd, struct db_str ntt)
{
  double db_tt;
  char tstr1[26], tstr2[26];

  strncpy(tstr1,ctime(&tdstart),sizeof(tstr1));
  strncpy(tstr2,ctime(&tddone),sizeof(tstr1));
  tstr1[24]=tstr2[24]='\0';

  /* Print timing to output file as well */
  fprintf(fd, "\n\n%ld residues in %ld query   sequences\n", qtt.length, qtt.entries);
  if (ntt.carry == 0) 
    fprintf(fd, "%ld residues in %ld library sequences\n", ntt.length, ntt.entries);
  else {
    db_tt = (double)ntt.carry*(double)LONG_MAX + (double)ntt.length;
    fprintf(fd, "%.0f residues in %ld library sequences\n", db_tt, ntt.entries);
  }

#ifndef COMP_THR
  fprintf(fd," Scomplib [%s]\n start: %s done: %s\n",mp_verstr,tstr1,tstr2);
#else
  fprintf(fd," Tcomplib [%s] (%d proc)\n start: %s done: %s\n", mp_verstr,
    fa_max_workers,tstr1,tstr2);
#endif
#ifndef COMP_MLIB
  fprintf(fd," Scan time: ");
  ptime(fd, tscan - tprev);
  fprintf (fd," Display time: ");
  ptime (fd, tdone - tscan);
#else
  fprintf(fd," Total Scan time: ");
  ptime(fd, ttscan);
  fprintf (fd," Total Display time: ");
  ptime (fd, ttdisp);
#endif
  fprintf (fd,"\n");
  fprintf (fd, "\nFunction used was %s [%s]\n", prog_func,verstr);
}

void fsigint()
{
  struct db_str db;

  db.entries = db.length = db.carry = 0;
  tdone = s_time();
  tddone = time(NULL);

  printf(" /*** interrupted ***/\n");
  if (outfd!=stdout) fprintf(outfd,"/*** interrupted ***/\n");
  fprintf(stderr,"/*** interrupted ***/\n");

  printsum(stdout,db);
  if (outfd!=stdout) printsum(outfd,db);

  exit(1);
}

/* save the sequence meta-data for a sequence if we need to re-use the
   rbuf/wbuf buffer */

void preserve_seq(struct buf2_str *lib_buf2_p,
		  struct seq_record *best_seqs,
		  struct beststr *best) {
  struct seq_record *dest_seq_p, *saved_seq_p;
  struct beststr *next_bbp;

  saved_seq_p = lib_buf2_p->best_save->seq;
  dest_seq_p = &best_seqs[lib_buf2_p->best_save - best];
  lib_buf2_p->best_save->seq = dest_seq_p;

  for (next_bbp = lib_buf2_p->best_save->bbp_link;
       (next_bbp != NULL) && (next_bbp->seq == saved_seq_p)
	 && (next_bbp->n1 == saved_seq_p->n1);
       next_bbp = next_bbp->bbp_link) {
    next_bbp->seq = dest_seq_p;
  }

  memcpy(dest_seq_p,lib_buf2_p->seq,sizeof(struct seq_record));
  dest_seq_p->aa1b = NULL;
}

/* in the current version (fasta_35_01) save_best is used by both
   threaded and unthreaded versions */

#define COPY_RST_P(d,s) 		\
{ d->rst.score[0] = s->rst.score[0];	\
  d->rst.score[1] = s->rst.score[1];	\
  d->rst.score[2] = s->rst.score[2];	\
  d->rst.comp = s->rst.comp;		\
  d->rst.H = s->rst.H;			\
  d->rst.escore = s->rst.escore;	\
  d->rst.segnum = s->rst.segnum;	\
  d->rst.seglen = s->rst.seglen;	\
}

void
save_best(struct buf_head *lib_bhead_p, const struct mngmsg *m_msp, struct pstruct *ppst, 
	  struct db_str *ldb, FILE *fdata, int *qsfnum, struct hist_str *histp,
	  void **pstat_voidp
	  )
{
  double zscore;
  int i_score;
  struct beststr *bbp;
  struct buf2_str *rbuf_p, *lib_buf2_p;
  int i, t_best, t_rbest, t_qrbest, tm_best, t_n1, sc_ix;
  double e_score, tm_escore, t_rescore, t_qrescore;
  int buf2_cnt;
  int jstats;

  sc_ix = ppst->score_ix;

  lib_buf2_p = lib_bhead_p->buf2;
  buf2_cnt = lib_bhead_p->buf2_cnt;
  
  t_best = t_rbest = t_qrbest = -BIGNUM;
  tm_escore = t_rescore = t_qrescore = FLT_MAX;

  while (buf2_cnt--) { /* count down the number of results */
    rbuf_p = lib_buf2_p++;	/* step through the results buffer */

    if (rbuf_p->rst.score[0] == -BIGNUM) continue;

    i_score = rbuf_p->rst.score[sc_ix];
    e_score = rbuf_p->rst.escore;

    /* need to look for frame 0 if TFASTA, then save stats at frame 6 */
    if (fdata) {
      fprintf(fdata,
	      "%-12s %5d %6d %d %.5f %.5f %4d %4d %4d %g %d %d %8lld\n",
	      rbuf_p->seq->libstr,
	      0,
	      rbuf_p->seq->n1,rbuf_p->frame,rbuf_p->rst.comp,rbuf_p->rst.H,
	      rbuf_p->rst.score[0],rbuf_p->rst.score[1],rbuf_p->rst.score[2],
	      rbuf_p->rst.escore, rbuf_p->rst.segnum, rbuf_p->rst.seglen,
	      rbuf_p->seq->lseek);
    }

    t_n1 = rbuf_p->seq->n1;
    if (i_score > t_best) tm_best = t_best = i_score;
    if (e_score < tm_escore) tm_escore = e_score;

    if (m_msp->qshuffle) {
      if (rbuf_p->qr_score > t_qrbest)
	t_qrbest = rbuf_p->qr_score;
      if (rbuf_p->qr_escore < t_qrescore)
	t_qrescore = rbuf_p->qr_escore;
      
      if (rbuf_p->frame == m_msp->nitt1 && nqstats < m_msp->shuff_max) {
	qstats[nqstats].n1 = rbuf_p->seq->n1;	/* save the best score */
	qstats[nqstats].comp =  rbuf_p->rst.comp;
	qstats[nqstats].H = rbuf_p->rst.H;
	qstats[nqstats].escore = t_qrescore;
	qstats[nqstats++].score = t_qrbest;
	t_qrbest = -BIGNUM;	/* reset t_qrbest, t_qrescore */
	t_qrescore = FLT_MAX;
      }
    }

    if (ppst->zsflag >= 10 && rbuf_p->r_score > t_rbest) {
      t_rbest = rbuf_p->r_score;
      t_rescore = rbuf_p->r_escore;
    }

    /* statistics done for best score of set */

    if (rbuf_p->frame == m_msp->nitt1) {
      ldb->entries++;
      ldb->length += t_n1;
      if (ldb->length > LONG_MAX) {
	ldb->length -= LONG_MAX; ldb->carry++;
      }

      if (nstats < MAX_STATS ) {
	stats[nstats].n1 = t_n1;
	stats[nstats].comp = rbuf_p->rst.comp;
	stats[nstats].H = rbuf_p->rst.H;
	if (ppst->zsflag >= 10) {
	  tm_best = t_rbest;
	  tm_escore = t_rescore;
	  t_rbest = -BIGNUM;
	  t_rescore = FLT_MAX;
	}
	stats[nstats].escore  = tm_escore;
	stats[nstats++].score = tm_best;
	t_best = -BIGNUM;
	tm_escore = FLT_MAX;
      }
      else if (ppst->zsflag > 0) {
	if (!stats_done) {
	  ppst->zsflag_f = process_hist(stats,nstats,m_msp, ppst,
				      histp, pstat_voidp,0);
	  kstats = nstats;
	  stats_done = 1;
	  for (i=0; i<MAX_BEST; i++) {
	    bestp_arr[i]->zscore = 
	      (*find_zp)(bestp_arr[i]->rst.score[ppst->score_ix],
			 bestp_arr[i]->rst.escore, bestp_arr[i]->seq->n1,
			 bestp_arr[i]->rst.comp, *pstat_voidp);
	  }
	}
#ifdef SAMP_STATS
	else {
	  if (!m_msp->escore_flg) {
	    jstats = nrand(++kstats);
	    if (jstats < MAX_STATS) {
	      stats[jstats].n1 = t_n1;
	      stats[jstats].comp = rbuf_p->rst.comp;
	      stats[jstats].H = rbuf_p->rst.H;
	      if (ppst->zsflag >= 10) {
		tm_best = t_rbest;
	      }
	      stats[jstats].score = tm_best;
	    }
	  }
	}
#endif
      }
    }

    /* calculate a z-score if possible, for histogram, and save threshold */
    if (stats_done) {

      zscore=(*find_zp)(i_score, e_score, rbuf_p->seq->n1,(double)rbuf_p->rst.comp,
			*pstat_voidp);

      if (rbuf_p->frame == m_msp->nitt1) {
	addhistz((*find_zp)(t_best, tm_escore, rbuf_p->seq->n1, (double) rbuf_p->rst.comp,
			    *pstat_voidp), histp);
	t_best = t_rbest = -BIGNUM;
	tm_escore = t_rescore = FLT_MAX;
      }
    }
    else zscore = (double) i_score;

    if (zscore > zbestcut) {
      if (nbest >= MAX_BEST) {
	bestfull = nbest-MAX_BEST/4;
	selectbestz(bestp_arr,bestfull-1,nbest);
	zbestcut = bestp_arr[bestfull-1]->zscore;
	nbest = bestfull;
      }
      bbp = bestp_arr[nbest++];

      COPY_RST_P(bbp, rbuf_p);

      bbp->seq = rbuf_p->seq;
      bbp->n1 = rbuf_p->seq->n1;
      if (rbuf_p->best_save) {	/* have we saved this seq before ? */
	/* yes - link back to the previous save, if it points to the
	   same seq_p */
	if (rbuf_p->best_save->seq == rbuf_p->seq) {
	  bbp->bbp_link = rbuf_p->best_save;
	}
	else {	/* break the link */
	  bbp->bbp_link = NULL;
	}
      }
      rbuf_p->best_save = bbp;	/* note where this seq is saved */	

      lib_bhead_p->have_best_save = 1;

      bbp->zscore = zscore;
      bbp->frame = rbuf_p->frame;
    }
  }
}

void
save_shuf(struct buf_head *lib_bhead_p, int nitt1, int shuff_max)
{
  struct buf2_str *rbuf_p, *lib_buf2_p;
  int t_rbest;
  double t_rescore;
  int buf2_cnt, jstats;
  static int kstats=0;

  lib_buf2_p = lib_bhead_p->buf2;
  buf2_cnt = lib_bhead_p->buf2_cnt;
  
  t_rbest = -BIGNUM;

  while (buf2_cnt--) { /* count down the number of results */
    rbuf_p = lib_buf2_p++;	/* step through the results buffer */

    if (rbuf_p->r_score > t_rbest) {
      t_rbest = rbuf_p->r_score;
      t_rescore = rbuf_p->r_escore;
    }

    /* statistics done for best score of set */

    if (rbuf_p->frame == nitt1) {
      if (nrstats < shuff_max ) { kstats = jstats = nrstats++; }
      else {	/* randomly replace */
	jstats = nrand(++kstats);
	if (jstats >= shuff_max) goto done;
      }

      rstats[jstats].n1 = rbuf_p->seq->n1;
      rstats[jstats].comp = rbuf_p->rst.comp;
      rstats[jstats].H = rbuf_p->rst.H;
      rstats[jstats].escore  = t_rescore;
      rstats[jstats].score = t_rbest;
done:
      t_rbest = -BIGNUM;
    }
  }
}


void
buf_do_work(unsigned char **aa0,  int n0,
	    struct buf_head *lib_bhead_p, 
	    struct pstruct *ppst, void **f_str) {
  
  int buf2_cnt;
  struct buf2_str *lib_buf2_p;

  buf2_cnt = lib_bhead_p->buf2_cnt;
  lib_buf2_p = lib_bhead_p->buf2;

  while (buf2_cnt--) {

    lib_buf2_p->rst.score[0] =
      lib_buf2_p->rst.score[1] =
      lib_buf2_p->rst.score[2] = -BIGNUM;

    if (lib_buf2_p->seq->n1 < ppst->n1_low ||
	lib_buf2_p->seq->n1 > ppst->n1_high ) {
      lib_buf2_p++;
      continue;
    }

    /*
    if (lib_buf2_p->seq->aa1b[0] == '\0') {
      fprintf(stderr,"invalid aa1 NULL at: %s frame: %d\n",
	      lib_buf2_p->seq->libstr, lib_buf2_p->frame);
    }
    */

    do_work (aa0[lib_buf2_p->frame], n0,
	     lib_buf2_p->seq->aa1b, lib_buf2_p->seq->n1,
	     lib_buf2_p->frame, ppst, f_str[lib_buf2_p->frame], 0,
	     &(lib_buf2_p->rst));

    lib_buf2_p++;
  }
  lib_bhead_p->have_results = 1;
}

void
buf_qshuf_work(unsigned char *aa0s,  int n0,
	       struct buf_head *lib_bhead_p, 
	       struct pstruct *ppst, void *qf_str,
	       int ix_score)
{
  int buf2_cnt;
  struct buf2_str *lib_buf2_p;
  struct rstruct rrst;

  buf2_cnt = lib_bhead_p->buf2_cnt;
  lib_buf2_p = lib_bhead_p->buf2;
  while (buf2_cnt--) {
    rrst.score[0] = rrst.score[1] = rrst.score[2] = -BIGNUM;

    if (lib_buf2_p->seq->n1 < ppst->n1_low ||
	lib_buf2_p->seq->n1 > ppst->n1_high ) {
      lib_buf2_p++;
      continue;
    }

    do_work (aa0s, n0,
	     lib_buf2_p->seq->aa1b, lib_buf2_p->seq->n1,
	     lib_buf2_p->frame, ppst, qf_str, 1,
	     &rrst);
    lib_buf2_p->qr_score = rrst.score[ix_score];
    lib_buf2_p->qr_escore = rrst.escore;

    lib_buf2_p++;
  }
}

void
buf_shuf_work(unsigned char **aa0,  int n0,
	      unsigned char *aa1s, int zs_win, 
	      struct buf_head *lib_bhead_p, 
	      struct pstruct *ppst, void **f_str,
	      int ix_score)
{
  int buf2_cnt;
  struct buf2_str *lib_buf2_p;
  struct rstruct rrst;

  buf2_cnt = lib_bhead_p->buf2_cnt;
  lib_buf2_p = lib_bhead_p->buf2;
  while (buf2_cnt--) {
    rrst.score[0] = rrst.score[1] = rrst.score[2] = -BIGNUM;

    if (lib_buf2_p->seq->n1 < ppst->n1_low ||
	lib_buf2_p->seq->n1 > ppst->n1_high ) {
      lib_buf2_p++;
      continue;
    }

    if (zs_win > 0) 
      wshuffle(lib_buf2_p->seq->aa1b,aa1s,lib_buf2_p->seq->n1,zs_win);
    else
      shuffle(lib_buf2_p->seq->aa1b,aa1s,lib_buf2_p->seq->n1);

    do_work (aa0[lib_buf2_p->frame], n0,
	     aa1s, lib_buf2_p->seq->n1,
	     lib_buf2_p->frame, ppst, f_str[lib_buf2_p->frame], 0,
	     &rrst);
    lib_buf2_p->rst.comp = rrst.comp;
    lib_buf2_p->rst.H = rrst.H;
    lib_buf2_p->r_score = rrst.score[ix_score];
    lib_buf2_p->r_escore = rrst.escore;

    lib_buf2_p++;
  }
  lib_bhead_p->have_results = 1;
}
