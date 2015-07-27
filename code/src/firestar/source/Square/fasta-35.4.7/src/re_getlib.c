/* re_getlib.c - re-acquire a sequence given lseek, lcont */

/* $Id: re_getlib.c 27 2008-06-30 16:27:31Z pearson $  */
/* $Revision: 37 $  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "mm_file.h"
#define XTERNAL
#include "uascii.h"

#define GETLIB (m_fptr->getlib)

/* modified Feb, 2008 to provide aa1a - annotation string */
extern int ann_scan(unsigned char *aa0, int n0, unsigned char **aa0a_p, int seqtype);

int
re_getlib(unsigned char *aa1,
	  unsigned char **aa1a,
	  int maxn,	/* longest aa1 */
	  int maxt3,	/* alternate maxn */
	  int loff,	/* overlap */
	  int lcont,
	  int term_code,
	  long *loffset,	/* offset from real start of sequence */
	  long *l_off_p,	/* coordinate of sequence start */
	  struct lmf_str *m_fptr) {

  unsigned char *aa1ptr;
  int *sascii_save;
  int dummy_flg;
  int icont, maxt, ccont, n1;
  char libstr[20];
  fseek_t lmark; 
  
  aa1ptr = aa1;
  icont=0;

  /* no longer do selection */
  m_fptr->sel_acc_p = NULL;

  *loffset = 0l;
  maxt = maxn;
  n1 = -1;

  /* to process sequences in pieces properly, if lcont > 0, then we
     must read all but the last sequence using the scanning sascii,
     and then read the last piece using the ann_ascii */

  if (lcont > 1) {
    for (ccont=0; ccont<lcont-1; ccont++) {

      n1= GETLIB(aa1ptr,maxt,libstr,sizeof(libstr),&lmark,&icont,m_fptr,l_off_p);

      if (term_code && m_fptr->lib_aa && aa1ptr[n1-1]!=term_code) {
	aa1ptr[n1++]=term_code;
	aa1ptr[n1]=0;
      }

      if (aa1ptr!=aa1) n1 += loff;

      if (icont) {
	maxt = maxt3;
	memcpy(aa1,&aa1[n1-loff],loff);
	aa1ptr= &aa1[loff];
	*loffset += n1 - loff;
      }
      else {
	maxt = maxn;
	aa1ptr=aa1;
      }
    }
  }

  /* for the last one, replace m_fptr->sascii with ann_ascii[], and
     read the sequence */

  /* change sascii matrix only if there are annotations - otherwise
     l_ann_ascii is not initialized */
  if (aa1a != NULL) {
    sascii_save = m_fptr->sascii;
    m_fptr->sascii = l_ann_ascii;

    n1= GETLIB(aa1ptr,maxt,libstr,sizeof(libstr),&lmark,&icont,m_fptr,l_off_p);
    m_fptr->sascii = sascii_save;

    /* now capture the annotation string */

    n1 = ann_scan(aa1ptr,n1,aa1a,0);
  }
  else {
    n1= GETLIB(aa1ptr,maxt,libstr,sizeof(libstr),&lmark,&icont,m_fptr,l_off_p);
  }    

  if (term_code && m_fptr->lib_aa && aa1ptr[n1-1]!=term_code) {
    aa1ptr[n1++]=term_code;
    aa1ptr[n1]=0;
  }

  if (aa1ptr!=aa1) n1 += loff;

  return n1;
}
