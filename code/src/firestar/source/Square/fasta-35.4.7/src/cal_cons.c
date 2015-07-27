/* cal_cons.c - routines for printing translated alignments for
   fasta, ssearch, ggsearch, glsearch  */

/* copyright (c) 1998, 1999, 2007 by William R. Pearson and the University of Virginia */

/*  $Id: cal_cons.c 27 2008-06-30 16:27:31Z pearson $ */
/* $Revision: 97 $  */

/* removed from dropgsw2.c, dropnfa.c April, 2007 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "defs.h"
#include "param.h"

#ifdef FASTA
#include "dropnfa.h"
#endif

#if defined(SSEARCH) || defined(OSEARCH)
#include "dropgsw2.h"
#endif

#ifdef LALIGN
#include "dropgsw2.h"
#endif

#include "a_mark.h"

static void update_code(char *al_str, int al_str_max, int op, int op_cnt);
extern void aancpy(char *to, char *from, int count, struct pstruct *ppst);

int calcons_a(const unsigned char *aa0, int n0,
	      const unsigned char *aa1, int n1,
	      int *nc,
	      struct a_struct *aln,
	      struct a_res_str *a_res,
	      struct pstruct *ppst,
	      char *seqc0, char *seqc1, char *seqca,
	      const char *ann_arr,
	      const unsigned char *aa0a, char *seqc0a,
	      const unsigned char *aa1a, char *seqc1a,
	      struct f_struct *f_str)
{
  int i0, i1, nn1;
  int op, lenc, nd, ns, itmp;
  const unsigned char *aa1p;
  char *sp0, *sp0a, *sp1, *sp1a, *spa, *sq;
  int *rp;
  int smins, mins, ntmp;
  int have_ann = 0;

  /* seqc0a/seqc1a are always allocated (or not) together, so just check one */
  have_ann = (seqc0a != NULL);

  if (ppst->ext_sq_set) {
    sq = ppst->sqx;
  }
  else {
    sq = ppst->sq;
  }

#ifndef TFASTA
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  aln->amin0 = a_res->min0;
  aln->amax0 = a_res->max0;
  aln->amin1 = a_res->min1;
  aln->amax1 = a_res->max1;

#ifndef LCAL_CONS
  /* will we show all the start ?*/
  if (min(a_res->min0,a_res->min1)<aln->llen || aln->showall==1)
    if (a_res->min0>=a_res->min1) {              /* aa0 extends more to left */
      smins=0;
      if (aln->showall==1) mins = a_res->min0;
      else mins = min(a_res->min0,aln->llcntx);
      aancpy(seqc0,(char *)aa0+a_res->min0-mins,mins,ppst);
      aln->smin0 = a_res->min0-mins;
      if ((mins-a_res->min1)>0) {
	memset(seqc1,' ',mins-a_res->min1);
	aancpy(seqc1+mins-a_res->min1,(char *)aa1p,a_res->min1,ppst);
	aln->smin1 = 0;
      }
      else {
	aancpy(seqc1,(char *)aa1p+a_res->min1-mins,mins,ppst);
	aln->smin1 = a_res->min1-mins;
      }
    }
    else {
      smins=0;
      if (aln->showall == 1) mins=a_res->min1;
      else mins = min(a_res->min1,aln->llcntx);
      aancpy(seqc1,(char *)(aa1p+a_res->min1-mins),mins,ppst);
      aln->smin1 = a_res->min1-mins;
      if ((mins-a_res->min0)>0) {
	memset(seqc0,' ',mins-a_res->min0);
	aancpy(seqc0+mins-a_res->min0,(char *)aa0,a_res->min0,ppst);
	aln->smin0 = 0;
      }
      else {
	aancpy(seqc0,(char *)aa0+a_res->min0-mins,mins,ppst);
	aln->smin0 = a_res->min0-mins;
      }
    }
  else {
    mins= min(aln->llcntx,min(a_res->min0,a_res->min1));
    smins=mins;
    aln->smin0=a_res->min0 - smins;
    aln->smin1=a_res->min1 - smins;
    aancpy(seqc0,(char *)aa0+a_res->min0-mins,mins,ppst);
    aancpy(seqc1,(char *)aa1p+a_res->min1-mins,mins,ppst);
  }
  /* set the alignment code to zero for context */
  memset(seqca,0,mins);
  if (have_ann) {
    memset(seqc0a,' ',mins);
    memset(seqc1a,' ',mins);
  }
#else		/* no flanking context */
  smins = mins = 0;
  aln->smin0=a_res->min0;
  aln->smin1=a_res->min1;
#endif

/* now get the middle */

  spa = seqca+mins;
  sp0 = seqc0+mins;
  sp1 = seqc1+mins;

  if (have_ann) {
    sp0a = seqc0a+mins;
    sp1a = seqc1a+mins;
  }

  rp = a_res->res;
  lenc = aln->nident = aln->nsim = aln->ngap_q = aln->ngap_l = aln->nfs = op = 0;
  i0 = a_res->min0;
  i1 = a_res->min1;
  
  while (i0 < a_res->max0 || i1 < a_res->max1) {
    if (op == 0 && *rp == 0) {	/* match/mismatch */
      op = *rp++;
      lenc++;

      if ((itmp=ppst->pam2[0][aa0[i0]][aa1p[i1]])<0) { *spa = M_NEG; }
      else if (itmp == 0) { *spa = M_ZERO;}
      else {*spa = M_POS;}

      /* nsim is incremented only for PAM scores >=0; thus, nsim can
	 be < nident because of identities (X:X and N:N that have
	 scores of -1 */

      if (*spa == M_POS || *spa==M_ZERO) {
	aln->nsim++;
      }

      if (have_ann) {
	if (aa0a) *sp0a++ = ann_arr[aa0a[i0]];
	else *sp0a++ = ' ';
	if (aa1a) *sp1a++ = ann_arr[aa1a[i1]];
	else *sp1a++ = ' ';
      }

      *sp0 = sq[aa0[i0++]];
      *sp1 = sq[aa1p[i1++]];

      /* nident counts all identities, regardless of pam score */
      if (toupper(*sp0) == toupper(*sp1)) {
	aln->nident++;
	*spa = M_IDENT;
      }
      else {
	if (ppst->nt_align) {
	  if ((toupper(*sp0) == 'T' && toupper(*sp1) == 'U') ||
	      (toupper(*sp0)=='U' && toupper(*sp1)=='T')) {
	    aln->nident++;
	    *spa = M_IDENT;
	  }
	  else if (toupper(*sp0) == 'N') aln->ngap_q++;
	  else if (toupper(*sp1) == 'N') aln->ngap_l++;
	}
	/*  additional code for 'C'=='U', 'K'=='O'
	else {
	  if ((toupper(*sp0) == 'C' && toupper(*sp1) == 'U') ||
	      (toupper(*sp0)=='U' && toupper(*sp1)=='C')  ||
	      (toupper(*sp0) == 'K' && toupper(*sp1) == 'O') ||
	      (toupper(*sp0)=='O' && toupper(*sp1)=='K')) {
	    aln->nident++;
	    *spa = M_IDENT;
	  }
	}
	*/
      }
      sp0++; sp1++; spa++;
    }
    else {	
      if (op==0) op = *rp++;
      if (op>0) {	/* insertion in aa0 */
	*sp0++ = '-';
	*sp1++ = sq[aa1p[i1++]];
	*spa++ = M_DEL;
	if (have_ann) {
	  *sp0a++ = ' ';
	  if (aa1a) *sp1a++ = ann_arr[aa1a[i1]];
	}
	op--;
	lenc++;
	aln->ngap_q++;
      }
      else {		/* insertion in aa1 */
	if (have_ann) {
	  if (aa0a) *sp0a++ = ann_arr[aa0a[i0]];
	  *sp1a++ = ' ';
	}

	*sp0++ = sq[aa0[i0++]];
	*sp1++ = '-';
	*spa++ = M_DEL;
	op++;
	lenc++;
	aln->ngap_l++;
      }
    }
  }

  *nc = lenc;
  if (have_ann) {
    *sp0a = *sp1a = '\0';
  }
  *spa = '\0';

#ifndef LCAL_CONS	/* have context around alignment */
/*	now we have the middle, get the right end */
  if (!aln->llcntx_flg) {
    ns = mins + lenc + aln->llen;	/* show an extra line? */
    ns -= (itmp = ns %aln->llen);	/* itmp = left over on last line */
    if (itmp>aln->llen/2) ns += aln->llen;  /* more than 1/2 , use another*/
    nd = ns - (mins+lenc);		/* this much extra */
  }
  else nd = aln->llcntx;

  if (nd > max(n0-a_res->max0,nn1-a_res->max1)) 
    nd = max(n0-a_res->max0,nn1-a_res->max1);
  
  if (aln->showall==1) {
    nd = max(n0-a_res->max0,nn1-a_res->max1);	/* reset for showall=1 */
    /* get right end */
    aancpy(seqc0+mins+lenc,(char *)aa0+a_res->max0,n0-a_res->max0,ppst);
    aancpy(seqc1+mins+lenc,(char *)aa1p+a_res->max1,nn1-a_res->max1,ppst);
    /* fill with blanks - this is required to use one 'nc' */
    memset(seqc0+mins+lenc+n0-a_res->max0,' ',nd-(n0-a_res->max0));
    memset(seqc1+mins+lenc+nn1-a_res->max1,' ',nd-(nn1-a_res->max1));
  }
  else {
    if ((nd-(n0-a_res->max0))>0) {
      aancpy(seqc0+mins+lenc,(char *)aa0+a_res->max0,(n0-a_res->max0),ppst);
      memset(seqc0+mins+lenc+n0-a_res->max0,' ',nd-(n0-a_res->max0));
    }
    else {
      aancpy(seqc0+mins+lenc,(char *)aa0+a_res->max0,nd,ppst);
    }

    if ((nd-(nn1-a_res->max1))>0) {
      aancpy(seqc1+mins+lenc,(char *)aa1p+a_res->max1,nn1-a_res->max1,ppst);
      memset(seqc1+mins+lenc+nn1-a_res->max1,' ',nd-(nn1-a_res->max1));
    }
    else {
      aancpy(seqc1+mins+lenc,(char *)aa1p+a_res->max1,nd,ppst);
    }
  }
  if (nd > 0 && have_ann) {
    memset(seqc0a+mins+lenc,' ',nd);
    memset(seqc1a+mins+lenc,' ',nd);
    /*
    ntmp = nd-(n0-a_res->max0);
    if (ntmp > 0) memset(seqc0a+mins+lenc+n0-a_res->max0,' ',ntmp);
    ntmp = nd-(nn1-a_res->max1);
    if (ntmp > 0) memset(seqc1a+mins+lenc+nn1-a_res->max1,' ',ntmp);
    */
  }
#else
  nd = 0;
#endif

  /*  fprintf(stderr,"%d\n",mins+lenc+nd); */

  lenc = mins + lenc + nd;

  return lenc;
}

void
calc_astruct(struct a_struct *aln_p, struct a_res_str *a_res_p) {

  aln_p->amin0 = a_res_p->min0;
  aln_p->amax0 = a_res_p->max0;
  aln_p->amin1 = a_res_p->min1;
  aln_p->amax1 = a_res_p->max1;
}

static void
update_code(char *al_str, int al_str_max, int op, int op_cnt) {

  char op_char[5]={"=-+*"};
  char tmp_cnt[20];

  sprintf(tmp_cnt,"%c%d",op_char[op],op_cnt);
  strncat(al_str,tmp_cnt,al_str_max);
}


/* build an array of match/ins/del - length strings */
int calc_code(const unsigned char *aa0, int n0,
	      const unsigned char *aa1, int n1,
	      struct a_struct *aln,
	      struct a_res_str *a_res,
	      struct pstruct *ppst,
	      char *al_str, int al_str_n, 
	      const char *ann_arr, 
	      const unsigned char *aa0a,
	      const unsigned char *aa1a,
	      char *ann_str, int ann_str_n,
	      struct f_struct *f_str)
{
  int i0, i1, nn1;
  int op, lenc;
  int p_op, op_cnt;
  const unsigned char *aa1p;
  char sp0, sp1, sp0a, *sq;
  int *rp;
  int have_ann=0;
  char ann_ch0, ann_ch1;
  char tmp_astr[MAX_STR];
  int sim_code;
  char *sim_sym= aln_map_sym[MX_ACC];

  if (aa0a != NULL && aa1a != NULL) { have_ann = 2;}
  else if (aa0a != NULL || aa1a != NULL) { have_ann = 1;}
  else {have_ann = 0;}

  if (ppst->ext_sq_set) {
    sq = ppst->sqx;
  }
  else {
    sq = ppst->sq;
  }

#ifndef TFASTA
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  /* this is already done by calc_astruct */
  /*
  aln->amin0 = a_res->min0;
  aln->amax0 = a_res->max0;
  aln->amin1 = a_res->min1;
  aln->amax1 = a_res->max1;
  */

  rp = a_res->res;
  lenc = aln->nident = aln->nsim = aln->ngap_q = aln->ngap_l = aln->nfs = op = p_op = 0;
  op_cnt = 0;

  i0 = a_res->min0;
  i1 = a_res->min1;
  
  while (i0 < a_res->max0 || i1 < a_res->max1) {
    if (op == 0 && *rp == 0) {

      sim_code = M_NEG;
      if (ppst->pam2[0][aa0[i0]][aa1p[i1]]>=0) { 
	aln->nsim++;
	sim_code = M_POS;
      }

      sp0 = sq[aa0[i0]];
      sp1 = sq[aa1p[i1]];

      if (p_op == 0 || p_op==3) {
	if (sp0 != '*' && sp1 != '*') {
	  if (p_op == 3) {
	    update_code(al_str,al_str_n-strlen(al_str),p_op,op_cnt);
	    op_cnt = 1; p_op = 0;
	  }
	  else {op_cnt++;}
	}
	else {
	  update_code(al_str,al_str_n-strlen(al_str),p_op,op_cnt);
	  op_cnt = 1; p_op = 3;
	}
      }
      else {
	update_code(al_str,al_str_n-strlen(al_str),p_op,op_cnt);
	op_cnt = 1; p_op = 0;
      }

      op = *rp++;
      lenc++;

      if (toupper(sp0) == toupper(sp1)) {
	aln->nident++;
	sim_code = M_IDENT;
      }
      else if (ppst->nt_align) {
	if ((toupper(sp0) == 'T' && toupper(sp1) == 'U') ||
	    (toupper(sp0)=='U' && toupper(sp1)=='T')) {
	  aln->nident++; sim_code = M_IDENT;
	}
	else if (toupper(sp0) == 'N') aln->ngap_q++;
	else if (toupper(sp1) == 'N') aln->ngap_l++;
      }

      /* check for an annotation */
      if (have_ann) {
	ann_ch0 = ann_ch1 = '\0';
	if (have_ann == 2 && (ann_arr[aa0a[i0]] != ' ' || ann_arr[aa1a[i1]] != ' ')) {
	  ann_ch0 = ann_arr[aa0a[i0]];
	  if (ann_ch0 == ' ') ann_ch0 = 'X';
	  ann_ch1 = ann_arr[aa1a[i1]];
	  if (ann_ch1 == ' ') ann_ch1 = 'X';
	}
	else if (aa0a != NULL && ann_arr[aa0a[i0]]!=' ') {
	  ann_ch0 = ann_arr[aa0a[i0]];
	  ann_ch1 = 'X';
	}
	else if (aa1a != NULL && ann_arr[aa1a[i1]]!=' ') {
	  ann_ch0 = 'X';
	  ann_ch1 = ann_arr[aa1a[i1]];
	}
	if (ann_ch0) {
	  sprintf(tmp_astr, "|%ld:%ld:%c%c:%c%c%c",
		  aln->q_offset+i0+1,aln->l_offset+i1+1,
		  ann_ch0,ann_ch1,sim_sym[sim_code],sp0,sp1);
	  strncat(ann_str, tmp_astr, ann_str_n - strlen(ann_str) - 1);
	  ann_str[ann_str_n-1] = '\0';
	}
      }
      i0++; i1++;
    }
    else {
      if (op==0) op = *rp++;
      if (op>0) {
	if (p_op == 1) { op_cnt++;}
	else {
	  update_code(al_str,al_str_n - strlen(al_str),p_op,op_cnt);
	  op_cnt = 1; p_op = 1;
	}
	op--; lenc++; i1++; aln->ngap_q++;
      }
      else {
	if (p_op == 2) { op_cnt++;}
	else {
	  update_code(al_str,al_str_n - strlen(al_str),p_op,op_cnt);
	  op_cnt = 1; p_op = 2;
	}
	op++; lenc++; i0++; aln->ngap_l++;
      }
    }
  }
  update_code(al_str,al_str_n - strlen(al_str),p_op,op_cnt);

  return lenc;
}

int calc_id(const unsigned char *aa0, int n0,
	    const unsigned char *aa1, int n1,
	    struct a_struct *aln, 
	    struct a_res_str *a_res,
	    struct pstruct *ppst,
	    struct f_struct *f_str)
{
  int i0, i1, nn1;
  int op, lenc;
  int sp0, sp1;
  const unsigned char *aa1p;
  int *rp;
  char *sq;
  
  if (ppst->ext_sq_set) { sq = ppst->sqx; }
  else { sq = ppst->sq; }

#ifndef TFASTA
  aa1p = aa1;
  nn1 = n1;
#else
  aa1p = f_str->aa1x;
  nn1 = f_str->n10;
#endif

  aln->amin0 = a_res->min0;
  aln->amax0 = a_res->max0;
  aln->amin1 = a_res->min1;
  aln->amax1 = a_res->max1;

  rp = a_res->res;
  lenc = aln->nident = aln->nsim = aln->ngap_q = aln->ngap_l = aln->nfs = op = 0;
  i0 = a_res->min0;
  i1 = a_res->min1;

  while (i0 < a_res->max0 || i1 < a_res->max1) {
    if (op == 0 && *rp == 0) {
      op = *rp++;
      lenc++;

      if (ppst->pam2[0][aa0[i0]][aa1p[i1]]>=0) { aln->nsim++;}

      sp0 = sq[aa0[i0++]];
      sp1 = sq[aa1p[i1++]];
      if (toupper(sp0) == toupper(sp1)) {aln->nident++;}
      else if (ppst->nt_align) {
	if ((toupper(sp0)=='T' && toupper(sp1)== 'U')||
	  (toupper(sp0)=='U' && toupper(sp1)=='T')) {aln->nident++;}
	else if (toupper(sp0) == 'N') aln->ngap_q++;
	else if (toupper(sp1) == 'N') aln->ngap_l++;
      }
    }
    else {
      if (op==0) op = *rp++;
      if (op>0) {op--; lenc++; i1++; aln->ngap_q++;}
      else {op++; lenc++; i0++; aln->ngap_l++; }
    }
  }
  return lenc;
}

