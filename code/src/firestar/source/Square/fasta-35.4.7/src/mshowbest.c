
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/*  $Id: mshowbest.c 27 2008-06-30 16:27:31Z pearson $ */
/* $Revision: 47 $  */

/*   29-Oct-2003 - changes so that bbp->seq->cont < 0 => aa1 sequence is
     already in aa1, no re_openlib or re_getlib required
*/

/*   14-May-2003 Changes to use a more consistent coordinate numbering
     system for displays.  aln->d_start[01] is now consistently used
     to report the start of the alignment in all functions, and
     mshowbest.c has been modified to use d_start[01] instead of
     d_start[01]-1.  aln->min[01] now starts at 0 for all functions;
     instead of 1 for some functions (dropnfa.c, dropgsw.c, dropfs2.c
     earlier).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "structs.h"
#include "param.h"

#ifndef PCOMPLIB
#include "mm_file.h"
#include "mw.h"
#else
#include "p_mw.h"
#endif


#define MAX_BLINE 256

#ifndef PCOMPLIB
/* function calls necessary to re_getlib() the sequence and, do
   alignments, if necessary
*/

#define RANLIB (m_fptr->ranlib)

int
re_getlib(unsigned char *, unsigned char **, 
	  int, int, int, int, int, long *, long *, 
	  struct lmf_str *m_fptr);

#include "drop_func.h"

struct lmf_str *re_openlib(struct lmf_str *, int outtty);
#endif

extern void cal_coord(int n0, int n1, long qoffset, long loffset,
		      struct a_struct *aln);

extern void calc_astruct(struct a_struct *aln_p, struct a_res_str *a_res_p);

void header_aux(FILE *);
void show_aux(FILE *, struct beststr *);
void w_abort (char *p, char *p1);

/* BBP_INFO get stuff directly from beststr or from beststr->desptr */
#ifdef PCOMPLIB
#define BBP_INFO(info) bbp->seq->info
#else
#define BBP_INFO(info) bbp->seq->info
#endif

extern double zs_to_bit(double, int, int);

/* showbest() shows a list of high scoring sequence descriptions, and
   their rst.scores.  If -m 9, then an additional complete set of
   alignment information is provided.

   If PCOMPLIB or m_msg.quiet then the number of high scores to be
   shown is pre-determined by m_msg.mshow before showbest is called.

   The comp_lib2.c version re_getlib()'s the sequence for its
   discription, and then does another alignment for -m 9 (Thus, it
   needs an f_str.  The PCOMPLIB version has everything available in
   beststr before showbest() is called.
*/

void showbest (FILE *fp, 
#ifndef PCOMPLIB
	       unsigned char **aa0, unsigned char *aa1save, int maxn,
#endif
	       struct beststr **bptr,int nbest,
	       int qlib, struct mngmsg *m_msp,
	       struct pstruct *ppst, struct db_str db,
	       char **info_gstring2
#ifndef PCOMPLIB
	       ,void **f_str
#endif
)
{
  unsigned char *aa1, *aa1a;
  int ntmp = 0;
  char bline[MAX_BLINE], fmt[40], pad[MAX_BLINE], rline[40];
  char l_name[128];
  int istart = 0, istop, ib;
  int nshow;
  int quiet;
  int r_margin;
  struct beststr *bbp;
  int n1tot;
  char *bp;
  char rel_label[12];
  char tmp_str[20], *seq_code, *ann_code;
  int seq_code_len, ann_code_len;
  long loffset;		/* loffset is offset from beginning of real sequence */
  long l_off;		/* l_off is the the virtual coordinate of residue 1 */
  int n0, n1;
  struct rstruct rst;
  int lc, seqc_max, annc_max, nident, ngap;
  float percent, gpercent;
  struct a_struct *aln_p;
  int *tres;
  int gi_num;
  char html_pre_E[120], html_post_E[120];

#ifndef PCOMPLIB
  struct lmf_str *m_fptr;
#endif

  strncpy(rel_label,"\0",2);
#ifdef SHOWREL
  strncpy(rel_label," related",sizeof(rel_label));
#endif
#ifdef SHOWUN
  strncpy(rel_label," unrelated",sizeof(rel_label));
#endif
  rel_label[sizeof(rel_label)-1]='\0';

#ifdef PCOMPLIB
  quiet = 1;
#else
  quiet = m_msp->quiet;
#endif

  n0 = m_msp->n0;

  if (m_msp->aln.llen > MAX_BLINE) m_msp->aln.llen = MAX_BLINE;

  if (ppst->zsflag < 0) r_margin = 10;
  else if (ppst->zsflag>=0  && m_msp->srelv > 1 ) r_margin = 19;
  else r_margin = 10;

  if (m_msp->markx & MX_M9SUMM && m_msp->show_code == SHOW_CODE_ID) {
#ifdef SHOWSIM
    r_margin += 15;
#else
    r_margin += 10;
#endif
  }

  if (m_msp->markx & MX_HTML) {
    strncpy(html_pre_E,"<font color=\"darkred\">",sizeof(html_pre_E));
    strncpy(html_post_E,"</font>",sizeof(html_post_E));

  }
  else {
    html_pre_E[0] = html_post_E[0] = '\0';
  }


  if (m_msp->nframe < 0) {
#ifndef SUPERFAMNUM
    sprintf(fmt,"%%-%ds (%%4d)",m_msp->aln.llen-r_margin);
#else
    sprintf(fmt,"%%-%ds [%%4d](%%4d)",m_msp->aln.llen-(r_margin+4));
#endif
  }
  else { sprintf(fmt,"%%-%ds (%%4d)",m_msp->aln.llen-(r_margin+4)); }

  memset(pad,' ',m_msp->aln.llen-(r_margin+6));
  pad[m_msp->aln.llen-(r_margin+12)]='\0';

  if (quiet != -1) {	/* quiet is set to -1 in comp_lib2.c to force
			   all significant hits to be shown */
    nshow = 20;
    if (m_msp->mshow == -1) {nshow = nbest;}		/* show all */
    /* show specified number */
    else if (m_msp->mshow_flg) {
      nshow = min (m_msp->mshow, nshow);
    }
  }
  else nshow = m_msp->nshow;

  if (quiet==0) istop = 20;
  else istop = nshow;

  if (quiet==0) {
    printf(" How many scores would you like to see? [%d] ",m_msp->nshow);
    fflush(stdout);
    if (fgets(rline,20,stdin)==NULL) exit(0);
    nshow = m_msp->nshow;
    if (rline[0]!='\n' && rline[0]!=0) sscanf(rline,"%d",&nshow);
    if (nshow<=0) nshow = min(20,nbest);
  }

  if ((bp = strchr (m_msp->qtitle, '\n')) != NULL) *bp = '\0';
/*   fprintf (fp, "%3d %s\n", qlib,m_msp->qtitle); */

  if (m_msp->markx & MX_HTML) fprintf(fp,"<p><tt><pre>\n");

  if (ppst->zsflag >= 0) {
    if (bptr[0]->rst.escore < m_msp->e_cut) {
      if (m_msp->z_bits==1) {/* show bit score */
	fprintf(fp,"\nThe best%s scores are:%s%s bits %sE(%ld)%s",
		rel_label,pad,m_msp->label,html_pre_E,ppst->zdb_size,html_post_E);
      }
      else {/* show z-score */
	fprintf(fp,"\nThe best%s scores are:%s%s z-sc %sE(%ld)%s",
		rel_label,pad,m_msp->label,html_pre_E,ppst->zdb_size,html_post_E);
      }
      header_aux(fp);
      if (m_msp->markx & MX_M9SUMM) {
	if (m_msp->show_code == SHOW_CODE_ID) {
#ifdef SHOWSIM
	  fprintf(fp," %%_id  %%_sim  alen");
#else
	  fprintf(fp," %%_id  alen");
#endif
	}
	else {
	if (m_msp->markx & MX_HTML && m_msp->show_code !=1) { fprintf(fp,"<!-- ");}
#ifndef SHOWSIM
	  fprintf(fp,"\t%%_id  %%_gid %4s  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs ",m_msp->f_id1);
#else
	  fprintf(fp,"\t%%_id  %%_sim %4s  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs ",m_msp->f_id1);
#endif
	}
	if (m_msp->show_code == SHOW_CODE_ALIGN) {	fprintf(fp," aln_code"); }
	if (m_msp->markx & MX_HTML && m_msp->show_code!=1) { fprintf(fp," -->");}
      }
      fprintf(fp,"\n");
    }
    else {
      fprintf(fp,"!! No library sequences with E() < %.2g\n",m_msp->e_cut);
      m_msp->nshow = 0;
      if (m_msp->markx & MX_HTML) fprintf(fp,"<p></tt></pre>\n");
      return;
    }
  }
  else {
    fprintf(fp,"\nThe best%s scores are:%s%s",rel_label,pad,m_msp->label);
    header_aux(fp);
    if (m_msp->markx & MX_M9SUMM) {
      if (m_msp->show_code == SHOW_CODE_ID) {
#ifdef SHOWSIM
	fprintf(fp," %%_id  %%_sm  alen");
#else
	fprintf(fp," %%_id  alen");
#endif
      }
      else {
#ifndef SHOWSIM
	fprintf(fp,"\t%%_id  %%_gid %4s  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs ",m_msp->f_id1);
#else
	fprintf(fp,"\t%%_id  %%_sim %4s  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs ",m_msp->f_id1);
#endif	/* SHOWSIM */
      }
    }
    if (m_msp->show_code == SHOW_CODE_ALIGN) {	fprintf(fp," aln_code"); }
    fprintf(fp,"\n");
  }

  istart = 0;
l1:
  istop = min(nbest,nshow);
  for (ib=istart; ib<istop; ib++) {
    bbp = bptr[ib];

#ifndef PCOMPLIB
#ifdef DEBUG
    if (bbp->seq->n1 != bbp->n1 ) {
      fprintf(stderr, " *** lib len error [%d!=%d] *** %s score %d\n",
	      bbp->seq->n1,bbp->n1, bbp->seq->libstr, bbp->rst.score[0]);
    }
#endif
#endif

#ifdef SUPERFAMNUM
    if (BBP_INFO(nsfnum) > 0 && sfn_cmp(m_msp->qsfnum_n,BBP_INFO(sfnum))) continue;
#ifdef SHOWUN
    if (BBP_INFO(nsfnum) > 0 && sfn_cmp(m_msp->qsfnum,BBP_INFO(sfnum))) {
      istop = min(istop+1,nbest);
    /*
      fprintf(stderr,"skipping %d: %d==%d\n",ib,m_msp->qsfnum,BBP_INFO(sfnum));
      */
      continue;
    }
#endif	/* SHOWUN */
#ifdef SHOWREL
    if (BBP_INFO(nsfnum) == 0 || (BBP_INFO(nsfnum) > 0 && !sfn_cmp(m_msp->qsfnum,BBP_INFO(sfnum)))) {
      istop = min(istop+1,nbest);
      continue;
    }
#endif	/* SHOWREL */
#endif	/* SUPERFAMNUM */
    if (quiet==1 && ppst->zsflag>=0) {
      if (bbp->rst.escore > m_msp->e_cut) {
	nshow = ib;
	goto done;
      }
      else if (bbp->rst.escore < m_msp->e_low) continue;
    }

#ifndef PCOMPLIB
    if ((m_fptr=re_openlib(bbp->seq->m_file_p,!m_msp->quiet))==NULL) {
      fprintf(stderr,"*** cannot re-open %s\n",bbp->seq->m_file_p->lb_name);
      exit(1);
    }
    RANLIB(bline,m_msp->aln.llen,bbp->seq->lseek,bbp->seq->libstr,m_fptr);
#else	/* PCOMPLIB */
  strncpy(bline,BBP_INFO(bline),m_msp->aln.llen-r_margin);
  bline[m_msp->aln.llen]='\0';
#endif

  /* l_name is used to build an HTML link from the bestscore line to
     the alignment.  It can also be used to discriminate multiple hits
     from the same long sequence.  This requires that fast_pan use -m 6. */

  strncpy(l_name,bline,sizeof(l_name)); /* get rid of text after second "|" */
  l_name[sizeof(l_name)-1]='\0';
  if ((bp=strchr(l_name,' '))!=NULL) *bp=0;
  if ((bp=strchr(&l_name[3],'|'))!=NULL) *bp='\0';
  if (m_msp->nframe > 2) sprintf(&l_name[strlen(l_name)],"_%d",bbp->frame+1);
  else if (m_msp->nframe > 0 && bbp->frame == 1)
    strncat(l_name,"_r",sizeof(l_name));
  if (bbp->seq->cont-1 > 0) {
    sprintf(tmp_str,":%d",bbp->seq->cont-1);
    strncat(l_name,tmp_str,sizeof(l_name)-strlen(l_name));
  }


#ifndef PCOMPLIB
  aln_p = &(m_msp->aln);

  if (m_msp->stages>1 || m_msp->markx & MX_M9SUMM) {
    if (bbp->seq->aa1b == NULL || (m_msp->ann_flg && bbp->seq->aa1_ann==NULL)) {

      /* get the sequence but don't save it */
      n1 = re_getlib(aa1save,
		     m_msp->ann_flg ? &(bbp->seq->aa1_ann) : NULL, 
		     maxn,m_msp->maxt3,
		     m_msp->l_overlap,bbp->seq->cont,m_msp->term_code,
		     &loffset,&l_off,bbp->seq->m_file_p);
      aa1 = aa1save;
      aa1a = bbp->seq->aa1_ann;
    }
    else {
      n1 = bbp->seq->n1;
      aa1 = bbp->seq->aa1b;
      aa1a = bbp->seq->aa1_ann;
      loffset = bbp->seq->l_offset;
      l_off = bbp->seq->l_off;
    }

    if (! m_msp->markx & MX_M9SUMM) {
      do_opt (aa0[bbp->frame], m_msp->n0, aa1, n1, bbp->frame, ppst, f_str[bbp->frame], &rst);
      bbp->rst.score[2]=rst.score[2];
    }
    else {
      if (!bbp->have_ares) {	/* showbest() can be called more than once */

	do_walign(aa0[bbp->frame],m_msp->n0, aa1, n1, bbp->frame, 
		  ppst, f_str[bbp->frame],
		  &bbp->a_res, &bbp->have_ares);
      
	/* if do_walign does not provide a fresh a_res,
	   then copy the re-used a_res to a new location */
	if (bbp->have_ares && !(bbp->have_ares & 0x2)) {
	  if ((tres = calloc(bbp->a_res.nres+1,sizeof(int)))!=NULL) {
	    memcpy(tres,bbp->a_res.res,sizeof(int)*bbp->a_res.nres);
	    bbp->a_res.res = tres;
	    bbp->have_ares |= 0x2;	/* set 0x2 if has local copy */
	  }
	  else {
	    bbp->have_ares = 0;		/* could not allocate memory */
	  }
	}
      }
      else {
	pre_cons(aa1,n1,bbp->frame,f_str[bbp->frame]);
      }

      aln_func_vals(bbp->frame, aln_p);

      seqc_max = bbp->a_res.nres + 4*m_msp->aln.llen+4;
      seq_code = NULL;
      seq_code_len = 0;
      if (m_msp->show_code == SHOW_CODE_ALIGN) {
	seq_code=(char *)calloc(seqc_max,sizeof(char));
	/* if we have an annotation string, allocate space for the encoded annotation */
	if (m_msp->ann_arr[0] != '\0') {
	  /* the annotation encoding can be considerably longer than the alignment encoding */
	  annc_max = 4*seqc_max;
	  ann_code=(char *)calloc(annc_max,sizeof(char));
	}
	else {
	  ann_code = NULL;
	  annc_max = 0;
	}
	if (seq_code != NULL) {
	  calc_astruct(aln_p, &bbp->a_res);

      /* we need this for offset information for calc_code, but it is
	 incomplete so we must do it again */
	  cal_coord(m_msp->n0,BBP_INFO(n1),
		    m_msp->q_offset + (m_msp->q_off-1) + (m_msp->sq0off-1),
		    loffset + (l_off-1) + (m_msp->sq1off-1),
		    aln_p);

	  lc=calc_code(aa0[bbp->frame],m_msp->n0,
		       aa1,n1, 
		       &m_msp->aln,&bbp->a_res,
		       ppst,seq_code,seqc_max,
		       m_msp->ann_arr,
		       m_msp->aa0a, aa1a,
		       ann_code, annc_max,
		       f_str[bbp->frame]);
	  seq_code_len = strlen(seq_code);
	  if (ann_code != NULL) ann_code_len = strlen(ann_code);
	  else ann_code_len = 0;
	}
      }
      else {
	lc=calc_id(aa0[bbp->frame],m_msp->n0,aa1,n1,
		   &m_msp->aln, &bbp->a_res,
		   ppst,f_str[bbp->frame]);
      }
      m_msp->aln.a_len = lc;

      nident = m_msp->aln.nident;
      if (lc > 0) percent = (100.0*(float)nident)/(float)lc;
      else percent = -1.00;

      ngap = m_msp->aln.ngap_q + m_msp->aln.ngap_l;
#ifndef SHOWSIM
      if (lc-ngap > 0) gpercent = (100.0*(float)nident)/(float)(lc-ngap);
      else gpercent = -1.00;
#else
      if (lc-ngap > 0) gpercent = (100.0*(float)m_msp->aln.nsim)/(float)(lc);
      else gpercent = -1.00;
#endif	/* SHOWSIM */

    }
  }
#endif	/* PCOMPLIB */

  n1tot = (BBP_INFO(n1tot_p)) ? *BBP_INFO(n1tot_p) : BBP_INFO(n1);

  bp = bline;
  if ((m_msp->markx & MX_HTML) && !strncmp(bline,"gi|",3)) {
    bp = strchr(bline+4,'|')+1;
    *(bp-1) = 0;
    gi_num = atoi(bline+3);
  }

#ifndef SUPERFAMNUM
  bp[m_msp->aln.llen-r_margin]='\0';
#else
  bp[m_msp->aln.llen-r_margin-5]='\0';
#endif

  if (m_msp->nframe == -1) bp[m_msp->aln.llen-r_margin]='\0';
  else bp[m_msp->aln.llen-(r_margin+4)]='\0';

#ifndef SUPERFAMNUM
  fprintf (fp, fmt,bp,n1tot);
#else
  if (m_msp->nframe == -1) {
    fprintf (fp, fmt,bp,BBP_INFO(sfnum[0]),n1tot);
  }
  else {fprintf (fp, fmt,bp,n1tot);}
#endif

  if (m_msp->nframe > 2) fprintf (fp, " [%d]", bbp->frame+1);
  else if (m_msp->nframe >= 0) fprintf(fp," [%c]",(bbp->frame > 0 ?'r':'f'));

  if (m_msp->srelv == 1) fprintf (fp, " %4d", bbp->rst.score[ppst->score_ix]);
  else {
    if (m_msp->srelv-1 > 0) fprintf (fp, " %4d", bbp->rst.score[0]);
    if (m_msp->srelv-1 > 1 || m_msp->stages>1)
      fprintf (fp, " %4d", bbp->rst.score[1]);
    fprintf (fp, " %4d", bbp->rst.score[ppst->score_ix]);
  }

  if (ppst->zsflag>=0) { 
    if (m_msp->z_bits==1) {
      fprintf (fp, " %.1f %s%7.2g%s",zs_to_bit(bbp->zscore,m_msp->n0,BBP_INFO(n1)),html_pre_E,bbp->rst.escore,html_post_E);
    }
    else fprintf (fp, " %.1f %s%7.2g%s",bbp->zscore,html_pre_E,bbp->rst.escore,html_post_E);
  }
  show_aux(fp,bbp);

#ifdef PCOMPLIB
  n1 = BBP_INFO(n1);
  percent = bbp->percent;
  gpercent = bbp->gpercent;
  aln_p = &(bbp->aln_d);
  seq_code = bbp->aln_code;
  seq_code_len = bbp->aln_code_n;
  ann_code = bbp->ann_code;
  ann_code_len = bbp->ann_code_n;
  loffset = bbp->seq->l_offset;
  l_off = 0;
#endif

  if (m_msp->markx & MX_M9SUMM) {
    if (m_msp->show_code != SHOW_CODE_ID) {
      /* we need the coordinates for annotated SHOW_CODE_ALIGN */
      cal_coord(m_msp->n0,BBP_INFO(n1),
		m_msp->q_offset + (m_msp->q_off-1) + (m_msp->sq0off-1),
		loffset + (l_off-1) + (m_msp->sq1off-1),
		aln_p);

      if (m_msp->markx & MX_HTML) fprintf(fp,"<!-- ");
      /*            %_id  %_sim s-w alen an0  ax0  pn0  px0  an1  ax1  pn1  px1 gapq gapl fs  */
      /*                    alignment    min  max            min  max */
      /*                    sequence coordinate    min  max            min  max */
      fprintf(fp,"\t%5.3f %5.3f %4d %4d %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %3d %3d %3d",
	      percent/100.0,gpercent/100.0, 
#ifndef PCOMPLIB
	      bbp->a_res.sw_score,
#else
	      bbp->sw_score,
#endif
	      aln_p->a_len,
	      aln_p->d_start0,aln_p->d_stop0,
	      aln_p->q_offset+1, aln_p->q_offset+m_msp->n0,
	      aln_p->d_start1,aln_p->d_stop1,
	      aln_p->l_offset+1, aln_p->l_offset+BBP_INFO(n1),
	      aln_p->ngap_q,aln_p->ngap_l,aln_p->nfs);
      if (m_msp->show_code == SHOW_CODE_ALIGN
	  && seq_code_len > 0 && seq_code != NULL) {
	fprintf(fp,"\t%s",seq_code);
	if (ann_code_len > 0 && ann_code != NULL) {
	  fprintf(fp,"\t%s",ann_code);
	}
	/*      fprintf(fp," [%2d:%d]",bbp->wrkr,bbp->seqnm); */


	/* if we are doing MX_M10FORM and -m 9c, then we want to keep
	   the alignment code string for the alignment output - otherwise, we
	   can free() it

	   If PCOMPLIB, then it is stored in bbp->ann_code*, otherwise, it's in a_res.
	*/

#ifndef PCOMPLIB
	if (m_msp->markx & MX_M10FORM)  {
	  /* save encoded alignments in a_res */
	  if ((bbp->a_res.aln_code = (char *)calloc(seq_code_len+1,sizeof(char)))!=NULL) {
	    strncpy(bbp->a_res.aln_code,seq_code,seq_code_len+1);
	    bbp->a_res.aln_code[seq_code_len] = '\0';
	    bbp->a_res.aln_code_n = seq_code_len;
	  }
	}
	/* always free the originally allocated encoding */
	free(seq_code);
	seq_code = NULL;
	seq_code_len = 0;
#else
	/* only free it if not to be used */
	if (!(m_msp->markx & MX_M10FORM))  {
	  free(bbp->aln_code);
	  bbp->aln_code_n = 0;
	}
#endif

	/* also clean up ann_code(_n) */

	if (ann_code_len > 0 && ann_code != NULL) {
#ifndef PCOMPLIB
	  if (m_msp->markx & MX_M10FORM)  {
	    /* save encoded annotations in a_res */
	    if ((bbp->a_res.ann_code = (char *)calloc(ann_code_len+1,sizeof(char)))!=NULL) {
	      strncpy(bbp->a_res.ann_code,ann_code,ann_code_len+1);
	      bbp->a_res.ann_code[ann_code_len] = '\0';
	      bbp->a_res.ann_code_n = ann_code_len;
	    }
	    else {
	      bbp->a_res.ann_code = NULL;
	      bbp->a_res.ann_code_n = 0;
	    }
	  }
	  free(ann_code);
	  ann_code = NULL;
	  ann_code_len = 0;
#else
	  if (!(m_msp->markx & MX_M10FORM))  {
	    free(bbp->ann_code);
	    bbp->ann_code_n = 0;
	  }
#endif
	}
      }
      if (m_msp->markx & MX_HTML) fprintf(fp," -->");
    }
    else {
#ifdef SHOWSIM
      fprintf(fp," %5.3f %5.3f %4d", percent/100.0,(float)aln_p->nsim/(float)aln_p->a_len,aln_p->a_len);
#else
      fprintf(fp," %5.3f %4d", percent/100.0,aln_p->a_len);
#endif
    }
  }
  if (m_msp->markx & MX_HTML) fprintf(fp," <A HREF=\"#%s\">align</A>",l_name);
  fprintf (fp, "\n");
  fflush(fp);
  }

  if (quiet==0) {
    printf(" More scores? [0] ");
    fflush(stdout);
    if (fgets(rline,20,stdin)==NULL) exit(0);
    ntmp = 0;
    if (rline[0]!='\n' && rline[0]!=0) sscanf(rline,"%d",&ntmp);
    if (ntmp<=0) ntmp = 0;
    if (ntmp>0) {
      istart = istop;
      nshow += ntmp;
      goto l1;
    }
  }
  else if (quiet == 1)
    if (ib < nbest && (ppst->zsflag>=0 && bbp->rst.escore < m_msp->e_cut)) {
      if (m_msp->mshow_flg && istop >= m_msp->mshow) goto done;
      istart=istop;
      nshow += 10;
      goto l1;
    }

 done:
  m_msp->nshow = nshow;

  if (m_msp->markx & MX_HTML) fprintf(fp,"</pre></tt><p><hr><p>\n");
  if (fp!=stdout) fprintf(fp,"\n");
}

/*
  q[] has one set of sfnums, 0 terminated
  s[] has second
  return first match or 0
*/
