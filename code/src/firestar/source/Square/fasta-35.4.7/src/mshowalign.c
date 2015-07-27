
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/*  $Id: mshowalign.c 27 2008-06-30 16:27:31Z pearson $ */
/* $Revision: 50 $  */

/* mshowalign.c - show sequence alignments in pvcomplib */

/* 
   this is a merged version of showalign.c that works properly with
   both the comp_lib (serial, threaded) and PCOMPLIB parallel versions
   of the programs.

   In the serial and current threaded versions of the programs,
   showalign gets a list of high scoring sequences and must
   re_getlib() the sequence, do_walign(), and then calculate the
   alignment.

   In the PCOMPLIB parallel versions, the worker programs do the
   aligning, so showalign() must send them the appropriate messages to
   have the alignment done, and then collect the alignment results

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "msg.h"
#include "structs.h"
#include "param.h"

#ifdef PCOMPLIB
#ifdef PVM_SRC
#include "pvm3.h"
extern int pinums[];
#endif
#ifdef MPI_SRC
#include "mpi.h"
#endif
#include "p_mw.h"
#else
#include "mm_file.h"
#include "mw.h"
#endif

#ifndef PCOMPLIB

/* used to position the library sequence for re_getlib - also gets
   description */
#define RANLIB (m_fptr->ranlib)

extern struct lmf_str *
re_openlib(struct lmf_str *, int outtty);

int
re_getlib(unsigned char *aa1, unsigned char **aa1a,
	  int maxn, int maxt,
	  int loff, int cont, int term_code,
	  long *l_offset, long *l_off, 
	  struct lmf_str *m_fptr);

#include "drop_func.h"

#endif


extern void calc_astruct(struct a_struct *aln_p, struct a_res_str *a_res_p);

extern void cal_coord(int n0, int n1, long qoffset, long loffset,
		      struct a_struct *aln);

void initseq(char **, char **, char **, int);
void initseq_ann(char **, char **, int);

void freeseq(char **, char **, char **);
void freeseq_ann(char **, char **);

void do_show(FILE *fp, int n0, int n1, int score,
	     char *name0, char *name1, int nml,
	     const struct mngmsg *m_msp, const struct pstruct *ppst,
	     char *seqc0, char *seqc0a,  char *seqc1, char *seqc1a,
	     char *seqca, int nc,
	     float percent, float gpercent, int lc,
	     struct a_struct *aln);

void
do_lav(FILE *fp, struct a_struct *aln, char *seqc, float percent, int is_mirror);

extern void discons(FILE *fd, const struct mngmsg *m_msg,
		    char *seqc0, char *seqc0a,
		    char *seqc1, char *seqc1a,
		    char *seqca, int nc, 
		    int n0, int n1, char *name0, char *name1, int nml,
		    struct a_struct *aln);

extern void disgraph(FILE *fd, int n0, int n1,
		     float percent, int score,
		     int min0, int min1, int max0, int max1, long sq0off,
		     char *name0, char *name1, int nml, int llen, int markx);

extern double (*find_zp)(int score, double escore, int length, double comp,void *);
extern double zs_to_bit(double, int, int);
extern double s_to_bit(int score, int n0, int  n1, void *pu);
extern double bit_to_E (double bit, int n0, int n1, long db_size, void *pu);
extern double zs_to_E(double zs, int n1, int dnaseq, long db_size, struct db_str db);

extern void
do_url1(FILE *, const struct mngmsg *, const struct pstruct *, char *, int,
	struct a_struct , long);

#ifndef A_MARK
#define A_MARK ">>"
#endif

static char l_name[200];	/* link name */

#ifdef PCOMPLIB
#define BBP_INFO(info) bbp->seq->info
#else
#define BBP_INFO(info) bbp->seq->info
#endif

/* this version does not check for m_msg->e_cut because nshow/nbest has
   already been set to limit on e_cut */

void showalign (FILE *fp,
#ifndef PCOMPLIB
		unsigned char **aa0, unsigned char *aa1save, int maxn,
#endif
		struct beststr **bptr, int nbest, int qlib, 
		const struct mngmsg *m_msp, struct pstruct *ppst, 
		char *info_gstring2
#ifndef PCOMPLIB
		, void **f_str
#endif
)
{
  unsigned char *aa1, *aa1a;
  char tmp_str[20];
  char info_str[200];
  char bline[2048], *bl_ptr, *bp, fmt[40];
  int tmp_len, l_llen;
  int t_have_ares;
  char name0[80], name0s[80], name1[200];
  int istart = 0, istop, i = 0, ib, nml;
  int l_ashow;
  int  n1tot;
  struct beststr *bbp;
  struct a_res_str *cur_ares_p;
  int nc, lc, maxc, maxca, seqc_max, annc_max;
  double zscore, bits;
  struct a_struct l_aln, *l_aln_p;
  float percent, gpercent;
  /* strings, lengths for conventional alignment */
  char *seqc0, *seqc0a, *seqc1, *seqc1a, *seqca, *seqc;
  int seqc_len, seqca_len;
  /* strings, lengths, for encoded alignment for MX10 */
  char *seq_code=NULL, *ann_code=NULL;
  int seq_code_len=0, ann_code_len=0;
  long loffset, l_off;
  int lsw_score;
  char html_pre_E[120], html_post_E[120];

#ifdef PCOMPLIB
  struct stage2_str liblist;
  struct qmng_str qm_msg;
#ifdef MPI_SRC
  int int_msg_b[10];
  MPI_Status mpi_status;
#endif
#else
  int n1;
  struct lmf_str *m_fptr;
  int ngap;
#endif

#ifdef PCOMPLIB
  /* this function has its own copy of qm_msg, so we must fill it
     appropriately */
  qm_msg.n0 = m_msp->n0;
  strncpy(qm_msg.libstr,m_msp->qtitle,sizeof(qm_msg.libstr));
#endif

  memcpy(&l_aln, &(m_msp->aln),sizeof(l_aln));
  l_aln_p = &l_aln;  /*   aln_p = &m_msp->aln; */

  /* set the name0,1 label length */
  if (m_msp->markx & MX_M10FORM) nml = 12;
  else nml = m_msp->nmlen;

  if (strlen(m_msp->qtitle) > 0) {
    if (m_msp->qtitle[0]=='>') strncpy(name0s,&m_msp->qtitle[1],sizeof(name0s));
    else strncpy(name0s,m_msp->qtitle,sizeof(name0s));
  }
  else {
    strncpy(name0s,m_msp->tname,sizeof(name0s));
  }
  name0s[sizeof(name0s)-1]='\0';

  if ((bp=strchr(name0s,' '))!=NULL) *bp='\0';

  if (m_msp->revcomp) name0[nml-1]='-';

  if (m_msp->markx & MX_HTML) {
    strncpy(html_pre_E,"<font color=\"darkred\">",sizeof(html_pre_E));
    strncpy(html_post_E,"</font>",sizeof(html_post_E));

  }
  else {
    html_pre_E[0] = html_post_E[0] = '\0';
  }

  l_llen = m_msp->aln.llen;
  if ((m_msp->markx & MX_M9SUMM) && m_msp->show_code != SHOW_CODE_ID) {
    l_llen += 40;
    if (l_llen > 200) l_llen=200;
  }

  sprintf(fmt,"%s%%-%ds (%%d %s)\n",A_MARK,l_llen-5,m_msp->sqnam);

  if (!(m_msp->markx&MX_M10FORM)) fprintf (fp,"\n");

  l_ashow = m_msp->ashow;
  if (l_ashow < 0) l_ashow = m_msp->nshow;
  istart = 0; istop = min(min(nbest,l_ashow),m_msp->nshow);

  for (ib=istart; ib<istop; ib++) {
    bbp = bptr[ib];

#ifdef SHOWUN
    if (BBP_INFO(nsfnum) > 0 && sfn_cmp(m_msp->qsfnum,BBP_INFO(sfnum))) {
      istop = min(istop+1,nbest);
      continue;
    }
#endif
    if (m_msp->quiet==1 && ppst->zsflag>=0 
	&& bbp->rst.escore < m_msp->e_low) continue;

#ifndef PCOMPLIB
    /* get the alignment and score by re-aligning */

    if ((m_fptr=re_openlib(bbp->seq->m_file_p,!m_msp->quiet))==NULL)
      exit(1);

    /* get the description - do not "edit" it yet */

    if (!(m_msp->markx & MX_M10FORM)){
      if (m_msp->long_info) {tmp_len = sizeof(bline)-1;}
      else {tmp_len = l_llen-5;}
      RANLIB(bline,tmp_len,bbp->seq->lseek,bbp->seq->libstr,bbp->seq->m_file_p);
      bline[tmp_len]='\0';
    }
    else {
      RANLIB(bline,sizeof(bline),bbp->seq->lseek,bbp->seq->libstr,bbp->seq->m_file_p);
      bline[sizeof(bline)-1]='\0';
    }

    if (bbp->seq->aa1b == NULL || (m_msp->ann_flg && bbp->seq->aa1_ann==NULL)) {

      n1 = re_getlib(aa1save, m_msp->ann_flg ? &(bbp->seq->aa1_ann) : NULL, maxn,m_msp->maxt3,
		     m_msp->l_overlap,bbp->seq->cont,m_msp->term_code,
		     &loffset,&l_off,bbp->seq->m_file_p);
      aa1 = aa1save;
      aa1a = bbp->seq->aa1_ann;
    }
    else {
      n1 = bbp->seq->n1;
      aa1 = bbp->seq->aa1b;
      aa1a = bbp->seq->aa1_ann;
      l_off = bbp->seq->l_off;
      loffset = bbp->seq->l_offset;
    }

#ifdef DEBUG
    if (n1 != BBP_INFO(n1)) {
      fprintf(stderr," library sequence: %s lengths differ: %d != %d\n",
	      bline,BBP_INFO(n1), n1);
      fprintf(stderr, "offset is: %lld\n",bbp->seq->lseek);
    }
#endif

    if (!bbp->have_ares) {
      /* here we do not bother to make a copy of a_res, because it will only be used once */
      do_walign(aa0[bbp->frame],m_msp->n0, aa1, n1, bbp->frame, ppst,
		f_str[bbp->frame], &bbp->a_res, &t_have_ares);
      bbp->have_ares = t_have_ares;
    }
    else {
      pre_cons(aa1,n1,bbp->frame,f_str[bbp->frame]);
    }

    cur_ares_p = &bbp->a_res;

    aln_func_vals(bbp->frame, l_aln_p);

#else	/* PCOMPLIB - get the alignment information from a worker */

    /* we have a sequence that we need an alignment for -
       send a message to the appropriate worker to produce an alignment 
       qm_msp->slist == 1  -> one alignment
       qm_msp->s_func == DO_ALIGN_FLG -> use the alignment function
       send mngmsg (MSEQTYPE)
       then send number of sequence to be aligned
    */

    qm_msg.slist = 1;
    qm_msg.s_func = DO_ALIGN_FLG;

    liblist.seqnm = bbp->seq->seqnm;
    liblist.frame = bbp->frame;
#ifdef PVM_SRC
    pvm_initsend(PvmDataRaw);
    pvm_pkbyte((char *)&qm_msg,sizeof(struct qmng_str),1);
    pvm_send(pinums[bbp->seq->wrkr],MSEQTYPE);

    pvm_initsend(PvmDataRaw);
    pvm_pkbyte((char *)&liblist,sizeof(struct stage2_str),1);
    pvm_send(pinums[bbp->seq->wrkr],LISTTYPE);
#endif
#ifdef MPI_SRC
    MPI_Send(&qm_msg,sizeof(struct qmng_str),MPI_BYTE,bbp->seq->wrkr,
	     MSEQTYPE,MPI_COMM_WORLD);
    MPI_Send(&liblist,sizeof(struct stage2_str),MPI_BYTE,bbp->seq->wrkr,
	     LISTTYPE,MPI_COMM_WORLD);
#endif
    /* information should be sent */
    /* pick up description */
    strncpy(bline,bbp->seq->bline,l_llen-5);
    bline[l_llen-5]='\0';
#endif	/* PCOMPLIB */

    if (strlen(bline)==0) {
      bline[0]='>';
      strncpy(&bline[1],m_msp->lname,l_llen-5);
      bline[l_llen-5]='\0';
    }

    /* re-format bline */
    while ((bp=strchr(bline,'\n'))!=NULL) *bp=' ';
    if (m_msp->long_info) {
      tmp_len = strlen(bline);
      bl_ptr = bline;
      if (!(m_msp->markx & MX_M10FORM)) while (tmp_len > l_llen) {
	for (i=l_llen; i>10; i--)
	  if (bl_ptr[i]==' ') {
	    bl_ptr[i]='\n';
	    break;
	  }
	if (i <= 10) break;
	tmp_len -= i;
	bl_ptr += i;
      }
      bline[sizeof(bline)-1]='\0';
    }

    n1tot = (BBP_INFO(n1tot_p)) ? *BBP_INFO(n1tot_p) : BBP_INFO(n1);

    strncpy(name1,bline,sizeof(name1));
    name1[sizeof(name1)-1] = '\0';

    if ((!m_msp->markx & MX_M10FORM)) name1[nml]='\0';
    if ((bp = strchr(name1,' '))!=NULL) *bp = '\0';

  /* l_name is used to build an HTML link from the bestscore line to
     the alignment.  It can also be used to discriminate multiple hits
     from the same long sequence.  Text must match that in showbest.c */

    strncpy(name1,bline,sizeof(name1));
    name1[sizeof(name1)-1]='\0';
    if ((bp = strchr(name1,' '))!=NULL) *bp = '\0';
    strncpy(l_name,name1,sizeof(l_name));
    l_name[sizeof(l_name)-1]='\0';
    if ((bp=strchr(&l_name[3],'|'))!=NULL) *bp='\0';
    if (m_msp->nframe > 2) sprintf(&l_name[strlen(l_name)],"_%d",bbp->frame+1);
    else if (m_msp->qframe >= 0 && bbp->frame == 1)
      strncat(l_name,"_r",sizeof(l_name));
    if (bbp->seq->cont-1 > 0) {
      sprintf(tmp_str,":%d",bbp->seq->cont-1);
      strncat(l_name,tmp_str,sizeof(l_name)-strlen(l_name));
    }

    if (!(m_msp->markx & MX_M10FORM)) name1[nml]='\0';

    /* print out score information; */

    if (m_msp->markx & MX_HTML ) {
      fprintf (fp,"<A name=%s>\n<tt><pre>\n",l_name);
    }
    strncpy(name0,name0s,nml);
    name0[nml]='\0';

    if (ppst->zsflag%10 == 6) {
      sprintf(info_str," comp: %.5f H: %.5f",bbp->rst.comp,bbp->rst.H);
    }
    else info_str[0]='\0';

    if (m_msp->markx & MX_M11OUT) {
      fprintf (fp, "s {\n   \"%s\" %ld %ld \n   \"%s\" %ld %ld\n}\n",
	       name0, m_msp->q_off,  m_msp->q_off+m_msp->n0-1, 
	       name1, bbp->seq->l_off, bbp->seq->l_off + bbp->seq->n1 - 1);
      fprintf (fp, "h {\n   \"%s\"\n   \"%s\"\n}\n", m_msp->qtitle, bline);
    }
    else if ((m_msp->markx & MX_ATYPE)!=7 && !(m_msp->markx & MX_M10FORM)) {
      fprintf (fp, fmt,bp=bline,n1tot);	 /* provides >>id  description (length) line */
    }
    else if (m_msp->markx & MX_M10FORM) {
      fprintf (fp,">>%s\n",bline);
    }

#ifdef PCOMPLIB
    /*  get the sw_score, alignment  information,  get seqc0, seqc1 */

#ifdef PVM_SRC
    /* get alignment lengths, percents */
    pvm_recv(pinums[bbp->seq->wrkr],ALN1TYPE);
    pvm_upkint(&nc,1,1);
    pvm_upkint(&lc,1,1);
    pvm_upkint(&maxc,1,1);

    pvm_upkfloat(&percent,1,1);
    pvm_upkfloat(&gpercent,1,1);

    pvm_upkint(&bbp->sw_score,1,1);
    pvm_upkbyte((char *)l_aln_p,sizeof(struct a_struct),1);

    initseq(&seqc0, &seqc1, &seqca, maxc);
    if (m_msp->ann_flg && m_msp->have_ann) {
      initseq_ann(&seqc0a, &seqc1a, maxc);
    }
    else { seqc0a = seqc1a = NULL;}

    pvm_recv(pinums[bbp->seq->wrkr],ALN2TYPE);
    pvm_upkbyte(seqc0,maxc,1);
    pvm_upkbyte(seqc1,maxc,1);
    pvm_upkbyte(seqca,maxc,1);
    if (m_msp->ann_flg && m_msp->have_ann) {
      pvm_upkbyte(seqc0a,maxc,1);
      pvm_upkbyte(seqc1a,maxc,1);
    }
#endif	/* PVM_SRC */
#ifdef MPI_SRC
    MPI_Recv(int_msg_b,4,MPI_INT,bbp->seq->wrkr,ALN1TYPE,MPI_COMM_WORLD,
	     &mpi_status);
    nc = int_msg_b[0];
    lc = int_msg_b[1];
    maxc = int_msg_b[2];
    bbp->sw_score = int_msg_b[3];
    MPI_Recv(&percent,1,MPI_FLOAT,bbp->seq->wrkr,ALN2TYPE,MPI_COMM_WORLD,
	     &mpi_status);
    MPI_Recv(&gpercent,1,MPI_FLOAT,bbp->seq->wrkr,ALN2TYPE,MPI_COMM_WORLD,
	     &mpi_status);
    MPI_Recv(l_aln_p,sizeof(struct a_struct),MPI_BYTE,
	     bbp->seq->wrkr,ALN3TYPE,MPI_COMM_WORLD,&mpi_status);

    initseq(&seqc0, &seqc1, &seqca, maxc);
    if (m_msp->ann_flg) {
      initseq_ann(&seqc0a, &seqc1a, maxc);
    }
    else { seqc0a = seqc1a = NULL;}


    MPI_Recv(seqc0,maxc,MPI_BYTE,bbp->seq->wrkr,ALN2TYPE,MPI_COMM_WORLD,&mpi_status);
    MPI_Recv(seqc1,maxc,MPI_BYTE,bbp->seq->wrkr,ALN3TYPE,MPI_COMM_WORLD,&mpi_status);
    MPI_Recv(seqca,maxc,MPI_BYTE,bbp->seq->wrkr,ALN3TYPE,MPI_COMM_WORLD,&mpi_status);
    if (m_msp->ann_flg) {
      MPI_Recv(seqc0a,maxc,MPI_BYTE,bbp->seq->wrkr,ALN2TYPE,MPI_COMM_WORLD,&mpi_status);
      MPI_Recv(seqc1a,maxc,MPI_BYTE,bbp->seq->wrkr,ALN2TYPE,MPI_COMM_WORLD,&mpi_status);
    }
#endif	/* MPI_SRC */

    loffset = bbp->seq->l_offset-l_off;
    lsw_score = bbp->sw_score;
#else	/* not PCOMPLIB */

    while (cur_ares_p != NULL && cur_ares_p->nres > 0) {

      /* estimate space for alignment consensus */
      if (m_msp->aln.showall==1) {
	maxc = cur_ares_p->nres + max(cur_ares_p->min0,cur_ares_p->min1)+
	  max((m_msp->n0-cur_ares_p->max0),(n1-cur_ares_p->max1))+4;
      }
      else {
	maxc = cur_ares_p->nres + 4*m_msp->aln.llen+4;
      }

      /* get space to put the sequence alignment consensus */
      initseq(&seqc0, &seqc1, &seqca, maxc);
      if (m_msp->ann_flg && (m_msp->aa0a != NULL || aa1a!=NULL)) {
	initseq_ann(&seqc0a, &seqc1a, maxc);
      }
      else { seqc0a = seqc1a = NULL;}

#ifdef LALIGN
      if ((m_msp->markx & MX_M11OUT) == MX_M11OUT) {	/* lav output - skip lots o stuff */
	lsw_score = cur_ares_p->sw_score;
	zscore = (*find_zp)(lsw_score, 0.0, BBP_INFO(n1), 0.0, m_msp->pstat_void);
	bits = zs_to_bit(zscore, m_msp->n0, BBP_INFO(n1));

	calc_astruct(l_aln_p, cur_ares_p);

	cal_coord(m_msp->n0,BBP_INFO(n1),
		  m_msp->q_offset+(m_msp->q_off-1)+(m_msp->sq0off-1),
		  loffset+(l_off-1)+(m_msp->sq1off-1),
		  l_aln_p);

	lc=calc_code(aa0[bbp->frame],m_msp->n0,
		     aa1,n1, 
		     l_aln_p, cur_ares_p,
		     ppst, seqc0,maxc,  m_msp->ann_arr,
		     m_msp->aa0a,bbp->seq->aa1_ann,
		     seqc0a, maxc,
		     f_str[bbp->frame]);

	seqc_len = strlen(seqc0);
	if (seqc0a) seqca_len = strlen(seqc0a);

	if (lc > 0) percent = (100.0*(float)l_aln_p->nident)/(float)lc;
	else percent = -1.00;

	fprintf (fp, "a {\n  s %d %.1f\n", lsw_score, bits);
	do_lav(fp, l_aln_p, seqc0, percent, 0);

	if (ppst->nseq == 1) {
	  fprintf (fp, "a {\n  s %d %.1f\n", lsw_score, bits);
	  do_lav(fp, l_aln_p, seqc0, percent, 1);
	}

	cur_ares_p = cur_ares_p->next;
	continue;
      }
#endif		/* ifdef LALIGN */

      /* build consensus from res, nres (done by workers if PCOMPLIB) */
      nc=calcons_a(aa0[bbp->frame],m_msp->n0,aa1,n1,
		   &lc,l_aln_p, cur_ares_p, ppst, 
		   seqc0, seqc1, seqca, 
		   m_msp->ann_arr,
		   m_msp->aa0a, seqc0a, aa1a, seqc1a,
		   f_str[bbp->frame]);

      /* PCOMPLIB workers return percent, gpercent, so calculate it here */
      if (lc > 0) percent = (100.0*(float)l_aln_p->nident)/(float)lc;
      else percent = -1.00;
      ngap = l_aln_p->ngap_q + l_aln_p->ngap_l;
#ifndef SHOWSIM
      if (lc-ngap> 0) gpercent =(100.0*(float)l_aln_p->nident)/(float)(lc-ngap);
#else
      if (lc > 0) gpercent =(100.0*(float)l_aln_p->nsim)/(float)lc;
#endif
      else gpercent = -1.00;

      lsw_score = cur_ares_p->sw_score;
#endif	/* not PCOMPLIB */

      if (max(strlen(seqc0),strlen(seqc1)) > nc) {
	fprintf(stderr," mshowalign: nc/maxc: %d/%d seqc0/1: %lu/%lu\n",
		nc,maxc,strlen(seqc0),strlen(seqc1));
      }

      /* here PCOMPLIB/comp_lib logic is the same */

#ifdef DEBUG
      /*
	if (lsw_score < bbp->rst.score[ppst->score_ix]) {
	fprintf(stderr," *** warning - SW score=%d < opt score=%d ***\n",
	lsw_score, bbp->rst.score[ppst->score_ix]);
	}
      */
#endif
      cal_coord(m_msp->n0,BBP_INFO(n1),
		m_msp->q_offset+(m_msp->q_off-1)+(m_msp->sq0off-1),
		bbp->seq->l_offset+(bbp->seq->l_off-1)+(m_msp->sq1off-1),
		l_aln_p);

      zscore = (*find_zp)(lsw_score, 0.0, BBP_INFO(n1), 0.0, m_msp->pstat_void);
      /*      bits = s_to_bit(lsw_score, m_msp->n0, BBP_INFO(n1), m_msp->pstat_void); */
      bits = zs_to_bit(zscore, m_msp->n0, BBP_INFO(n1));

      if ((m_msp->markx & MX_ATYPE)!=7 && !(m_msp->markx & MX_M10FORM)) {

	/* this code makes sense for library searches, but not for
	   multiple non-intersecting alignments */

#ifndef LALIGN
	if (m_msp->nframe > 2) 
	  fprintf (fp, "Frame: %d",bbp->frame+1);
	else if (m_msp->nframe > 1) 
	  fprintf (fp, "Frame: %c",(bbp->frame? 'r': 'f'));
	else if (m_msp->qframe >= 0 && bbp->frame > 0 ) {
	  fputs("rev-comp",fp);
	  name0[nml-1]='\0';
	  strcat(name0,"-");
	}

	if (m_msp->arelv > 0)
	  fprintf (fp, " %s: %3d", m_msp->alab[0],bbp->rst.score[0]);
	if (m_msp->arelv > 1)
	  fprintf (fp, " %s: %3d", m_msp->alab[1],bbp->rst.score[1]);
	if (m_msp->arelv > 2)
	  fprintf (fp, " %s: %3d", m_msp->alab[2],bbp->rst.score[2]);
	fprintf (fp,"%s",info_str);
	if (ppst->zsflag>=0) 
	  fprintf (fp, "  Z-score: %4.1f  bits: %3.1f %sE(): %4.2g%s", 
		   bbp->zscore,zs_to_bit(bbp->zscore,m_msp->n0,BBP_INFO(n1)),
		   html_pre_E,bbp->rst.escore,html_post_E);
	fprintf (fp, "\n");
#else /* LALIGN */
	if ((m_msp->markx & MX_M11OUT) == 0) {
	  fprintf (fp, " %s score: %d; ", m_msp->alabel, lsw_score);
	  fprintf (fp," %3.1f bits; E(%ld) <  %.2g\n", bits, ppst->zdb_size,
		   zs_to_E(zscore, BBP_INFO(n1), ppst->dnaseq, ppst->zdb_size, m_msp->db));
	}
#endif
      }
      else if (m_msp->markx & MX_M10FORM) {
#ifndef LALIGN
	if (m_msp->qframe > -1) {
	  if (m_msp->nframe > 2) {
	    fprintf(fp,"; %s_frame: %d\n",m_msp->f_id0,bbp->frame+1);
	  }
	  else {
	    fprintf(fp,"; %s_frame: %c\n",m_msp->f_id0,(bbp->frame > 0? 'r':'f'));
	  }
	}
	fprintf (fp, "; %s_%s: %3d\n", m_msp->f_id0,m_msp->alab[0],bbp->rst.score[0]);
	if (m_msp->arelv > 1)
	  fprintf (fp,"; %s_%s: %3d\n", m_msp->f_id0,m_msp->alab[1],bbp->rst.score[1]);
	if (m_msp->arelv > 2)
	  fprintf (fp,"; %s_%s: %3d\n", m_msp->f_id0,m_msp->alab[2],bbp->rst.score[2]);
	if (info_str[0]) fprintf(fp,"; %s_info: %s\n",m_msp->f_id0,info_str);
	if (ppst->zsflag>=0) 
	  fprintf (fp,"; %s_z-score: %4.1f\n; %s_bits: %3.1f\n; %s_expect: %6.2g\n",
		   m_msp->f_id0,bbp->zscore,
		   m_msp->f_id0,zs_to_bit(bbp->zscore, m_msp->n0, BBP_INFO(n1)),
		   m_msp->f_id0,bbp->rst.escore);
#else
	if ((m_msp->markx & MX_M11OUT) == 0) {
	  fprintf (fp,"; %s_%s score: %d\n", m_msp->f_id0, m_msp->alab[0], lsw_score);
	  fprintf (fp,"; %s_z-score: %4.1f\n; %s_bits %3.1f\n; %s_expect: %6.2g\n",
		   m_msp->f_id0, zscore, m_msp->f_id0, bits, m_msp->f_id0, 
		   zs_to_E(zscore, BBP_INFO(n1), ppst->dnaseq, ppst->zdb_size, m_msp->db));
	}
#endif	
      }

      do_show(fp, m_msp->n0, BBP_INFO(n1), lsw_score, name0, name1, nml,
	      m_msp, ppst, seqc0, seqc0a, seqc1, seqc1a, seqca,
	      nc, percent, gpercent, lc, l_aln_p);

      /* display the encoded alignment left over from showbest()*/

      if ((m_msp->markx & MX_M10FORM) &&
	  (m_msp->markx & MX_M9SUMM) && 
	  (m_msp->show_code == SHOW_CODE_ALIGN)) {

#ifdef PCOMPLIB
	seq_code = bbp->aln_code;
	seq_code_len = bbp->aln_code_n;
	ann_code = bbp->ann_code;
	ann_code_len = bbp->ann_code_n;
#else
	seq_code = cur_ares_p->aln_code;
	seq_code_len = cur_ares_p->aln_code_n;
	ann_code = cur_ares_p->ann_code;
	ann_code_len = cur_ares_p->ann_code_n;
#endif

	if (seq_code_len > 0 && seq_code != NULL) {
	  fprintf(fp,"; al_code: %s\n",seq_code);
	  free(seq_code);
	  if (ann_code_len > 0 && ann_code != NULL) {
	    fprintf(fp,"; al_code_ann: %s\n",ann_code);
	    free(ann_code);
	  }
	}
      }

      if (m_msp->markx & MX_HTML) fprintf(fp,"</pre></tt>\n<hr>\n");
      fflush(fp);

      freeseq(&seqc0,&seqc1, &seqca);
      freeseq_ann(&seqc0a, &seqc1a);

#ifndef PCOMPLIB

      cur_ares_p = cur_ares_p->next;

      if (cur_ares_p != NULL) {
	if (m_msp->markx & MX_HTML) {
	  fprintf (fp,"<tt><pre>\n");
	}
	else {
	  fprintf(fp,">--\n");
	}
      }		/* done finishing up  */
    }		/* while (cur_ares_p) */
    /* we are done displaying the alignment - be sure to free a_res
       memory */
#endif
  }
  if (fp!=stdout) fprintf(fp,"\n");
}

void do_show(FILE *fp, int n0,int n1, int score,
	     char *name0, char *name1, int nml,
	     const struct mngmsg *m_msp, const struct pstruct *ppst,
	     char *seqc0, char *seqc0a,  char *seqc1, char *seqc1a,
	     char *seqca, int nc,
	     float percent, float gpercent, int lc,
	     struct a_struct *aln)
{
  int tmp;

  if (m_msp->markx & MX_AMAP && (m_msp->markx & MX_ATYPE)==7)
    disgraph(fp, n0, n1, percent, score,
	     aln->amin0, aln->amin1, aln->amax0, aln->amax1, m_msp->sq0off,
	     name0, name1, nml, aln->llen, m_msp->markx);
  else if (m_msp->markx & MX_M10FORM) {
    if (ppst->sw_flag && m_msp->arelv>0)
      fprintf(fp,"; %s_score: %d\n",m_msp->f_id1,score);
    fprintf(fp,"; %s_ident: %5.3f\n",m_msp->f_id1,percent/100.0);
#ifndef SHOWSIM
    fprintf(fp,"; %s_gident: %5.3f\n",m_msp->f_id1,gpercent/100.0);
#else
    fprintf(fp,"; %s_sim: %5.3f\n",m_msp->f_id1,gpercent/100.0);
#endif

    fprintf(fp,"; %s_overlap: %d\n",m_msp->f_id1,lc);
    discons(fp, m_msp,
	    seqc0, seqc0a,
	    seqc1, seqc1a, seqca, nc,
	    n0, n1, name0, name1, nml, aln);
  }
  else {
#ifndef LALIGN
    fprintf(fp,"%s score: %d; ",m_msp->alabel, score);
#endif
#ifndef SHOWSIM
    fprintf(fp,"%4.1f%% identity (%4.1f%% ungapped) in %d %s overlap (%ld-%ld:%ld-%ld)\n",
	    percent,gpercent,lc,m_msp->sqnam,aln->d_start0,aln->d_stop0,
	    aln->d_start1,aln->d_stop1);
#else
    fprintf(fp,"%4.1f%% identity (%4.1f%% similar) in %d %s overlap (%ld-%ld:%ld-%ld)\n",
	    percent,gpercent,lc,m_msp->sqnam,aln->d_start0,aln->d_stop0,
	    aln->d_start1,aln->d_stop1);
#endif

    if (m_msp->markx & MX_HTML) {
      do_url1(fp, m_msp, ppst, l_name,n1,*aln,aln->l_offset);
    }

    if (m_msp->markx & MX_AMAP && (m_msp->markx & MX_ATYPE)!=7) {
      fputc('\n',fp);
      tmp = n0;

      if (m_msp->qdnaseq == SEQT_DNA && m_msp->ldnaseq== SEQT_PROT)
	tmp /= 3;

      disgraph(fp, tmp, n1, percent, score,
	       aln->amin0, aln->amin1,
	       aln->amax0, aln->amax1,
	       m_msp->sq0off,
	       name0, name1, nml, aln->llen,m_msp->markx);
    }

    discons(fp, m_msp,
	    seqc0, seqc0a, seqc1, seqc1a, seqca, nc,
	    n0, n1, name0, name1, nml, aln);

    fputc('\n',fp);

  }
}

void 
do_lav(FILE *fp, struct a_struct *aln_p, char *seqc,
       float percent, int is_mirror) {
  int cur_b0, cur_b1, cur_e0, cur_e1;
  int ipercent;
  long len;
  char *seqc_p, *num_e;

  ipercent = (int)(percent+0.5);

  cur_b0 = aln_p->d_start0;
  cur_b1 = aln_p->d_start1;
  cur_e0 = aln_p->d_stop0;
  cur_e1 = aln_p->d_stop1;

  if (!is_mirror) {
    fprintf (fp, "  b %d %d\n  e %d %d\n", 
	     cur_b0, cur_b1, cur_e0, cur_e1);
  }
  else  {
    fprintf (fp, "  b %d %d\n  e %d %d\n", 
	     cur_b1, cur_b0, cur_e1, cur_e0);
  }

  seqc_p = seqc;
  
  while (*seqc_p) {
    if (*seqc_p == '=') {
      len = strtol(seqc_p+1, &num_e, 10);
      cur_e0 = cur_b0 + len - 1;
      cur_e1 = cur_b1 + len - 1;
      if (!is_mirror) {
	fprintf(fp, "  l %d %d %d %d %d\n",
		cur_b0, cur_b1, cur_e0, cur_e1,
		ipercent);
      }
      else  {
	fprintf(fp, "  l %d %d %d %d %d\n",
		cur_b1, cur_b0, cur_e1, cur_e0,
		ipercent);
      }
      cur_b0 = cur_e0 + 1;
      cur_b1 = cur_e1 + 1;
    }
    else if (*seqc_p == '+') {
      len = strtol(seqc_p+1, &num_e, 10);
      cur_b0 += len;
    }
    else {
      len = strtol(seqc_p+1, &num_e, 10);
      cur_b1 += len;
    }
    seqc_p = num_e;
  }

  fprintf (fp, "}\n");
}


#ifndef MPI_SRC
void	/* initialize consensus arrays */
initseq(char **seqc0, char **seqc1, char **seqca, int seqsiz)
{
  *seqc0=(char *)calloc((size_t)seqsiz*3,sizeof(char));
  if (*seqc0==NULL)
    {fprintf(stderr,"cannot allocate consensus arrays %d\n",seqsiz);
     exit(1);}
  *seqc1=*seqc0 + seqsiz;
  *seqca=*seqc1 + seqsiz;
}

void freeseq(char **seqc0, char **seqc1, char **seqca)
{
  free(*seqc0);
  *seqc0 = *seqc1 = *seqca = NULL;
}

void	/* initialize consensus annotation arrays */
initseq_ann(char **seqc0a, char **seqc1a, int seqsiz)
{
  *seqc0a=(char *)calloc((size_t)seqsiz*5,sizeof(char));
  if (*seqc0a==NULL)
    {fprintf(stderr,"cannot allocate consensus arrays %d\n",seqsiz);
     exit(1);}
  *seqc1a=*seqc0a + seqsiz;
}

void freeseq_ann(char **seqc0a, char **seqc1a)
{
  if (*seqc0a != NULL) {
    free(*seqc0a);
    *seqc0a = *seqc1a = NULL;
  }
}
#endif
