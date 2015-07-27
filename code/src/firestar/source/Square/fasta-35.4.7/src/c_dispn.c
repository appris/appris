/*	dispn.c	associated subroutines for matching sequences */

/*  $Id: c_dispn.c 27 2008-06-30 16:27:31Z pearson $ */
/* $Revision: 28 $  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "defs.h"
#include "structs.h"
#include "param.h"

#define XTERNAL

#define YES 1
#define NO 0

#define MAXOUT 201

/* the seqca[] array has the following codes:
   0 - no alignment symbol
   1 - align; pam < 0
   2 - align; pam == 0
   3 - align; pam > 0
   4 - align; ident
   5 - align; del

   the map_sym arrays determine the value to be displayed with each
   type of aligned residue
*/

#include "a_mark.h"

void
discons(FILE *fd, const struct mngmsg *m_msg, 
	char *seqc0, char *seqc0a,
	char *seqc1, char *seqc1a, char *seqca, int nc,
	int n0, int n1, char *name0, char *name1, int nml,
	struct a_struct *aln)
{
  char line[3][MAXOUT];	/* alignment lines [0,2], similarity code [1] */
  char cline[2][MAXOUT+10], *clinep[2];	/* coordinate line */
  int il, i, lend, loff, id;
  int del0, del1, ic, ll0, ll1, ll01, cl0, cl1, rl0, rl1;
  int ic_save;
  char *map_sym_p;
  int l_llen;
  int ioff0, ioff00, ioff1, ioff10;
  long qqoff, lloff;
  int llsgn, llfact, qlsgn, qlfact, qfx0, qfxn, lfx0, lfxn;
  int have_res;
  char *name01;
  char blank[MAX_UID], afmt[32];

  memset(blank,' ',sizeof(blank)-1);
  blank[sizeof(blank)-1]='\0';

  if (nml > 6) {
    blank[nml-6]='\0';
    sprintf(afmt,"%%-%ds %%s\n",nml);
  }
  else {
    blank[0]='\0';
    strncpy(afmt,"%-6s %s\n",sizeof(afmt));
  }

  clinep[0]=cline[0]+1;
  clinep[1]=cline[1]+1;

  if (aln->qlfact == 0) {qlfact = 1;}
  else qlfact = aln->qlfact;
  if (aln->qlrev == 1) {
    qlsgn = -1;
    qfx0 = 0;
    qfxn = 1;
  }
  else {
    qlsgn = 1;
    qfx0 = 1;
    qfxn = 0;
  }

  if (aln->llfact == 0) {llfact = 1;}
  else llfact = aln->llfact;

  if (aln->llrev == 1) {
    llsgn = -1;
    lfx0 = 0;
    lfxn = 1;
  }
  else {
    llsgn = 1;
    lfx0 = 1;
    lfxn = 0;
  }

  l_llen = aln->llen;
  if ((m_msg->markx & MX_M9SUMM) && m_msg->show_code != 1) { l_llen += 40; }

  if ((m_msg->markx & MX_ATYPE)==2) name01=name1;
  else name01 = "\0";

  ioff0=aln->smin0;
  ioff00 = ioff0;
  ioff1=aln->smin1;
  ioff10 = ioff1;
  
  if (m_msg->markx& MX_AMAP && (m_msg->markx & MX_ATYPE)==7) return;

  /* set *map_sym_p to correct match symbol */
  map_sym_p = aln_map_sym[MX_A0];
  if ((m_msg->markx&MX_ATYPE)==1) {map_sym_p = aln_map_sym[MX_A1];}
  else if ((m_msg->markx&MX_ATYPE)==2) {map_sym_p = aln_map_sym[MX_A2];}
  else if (m_msg->markx&MX_M10FORM) { map_sym_p = aln_map_sym[MX_A10]; }

  if (m_msg->markx & MX_ASEP) {
    fprintf(fd,">%s ..\n",name0);
    for (i=0; i<nc && seqc0[i]; i++) {
   /* if (seqc0[i]=='-') fputc('.',fd); else */
      fputc(seqc0[i],fd);
      if (i%50 == 49) fputc('\n',fd);
    }
    if ((i-1)%50 != 49) fputc('\n',fd);
    fprintf(fd,">%s ..\n",name1);
    for (i=0; i<nc && seqc1[i]; i++) {
    /* if (seqc1[i]=='-') fputc('.',fd); else */
      fputc(seqc1[i],fd);
      if (i%50 == 49) fputc('\n',fd);
    }
    if ((i-1)%50 != 49) fputc('\n',fd);
    return;
  }

  if (m_msg->markx & MX_M10FORM) {
    fprintf(fd,">%s ..\n",name0);
    fprintf(fd,"; sq_len: %d\n",n0);
    fprintf(fd,"; sq_offset: %ld\n",aln->q_offset+1);
    fprintf(fd,"; sq_type: %c\n",m_msg->sqtype[0]);
    fprintf(fd,"; al_start: %ld\n",aln->d_start0);
    fprintf(fd,"; al_stop: %ld\n",aln->d_stop0);
    /* in the past, this al_display_start does not include sq0off */
    fprintf(fd,"; al_display_start: %ld\n",
	    aln->q_offset+qlsgn*ioff0*aln->llmult+qfx0);

    have_res = 0;
    for (i=0; i<nc && seqc0[i]; i++) {
      if (!have_res && seqc0[i]==' ') fputc('-',fd);
      else if (seqc0[i]==' ') break;
      else {
	have_res = 1;
	fputc(seqc0[i],fd);
      }
      if (i%50 == 49) fputc('\n',fd);
    }
    if ((i-1)%50!=49 || seqc0[i-1]==' ') fputc('\n',fd);
    fprintf(fd,">%s ..\n",name1);
    fprintf(fd,"; sq_len: %d\n",n1);
    fprintf(fd,"; sq_offset: %ld\n",aln->l_offset+1);
    fprintf(fd,"; sq_type: %c\n",m_msg->sqtype[0]);
    fprintf(fd,"; al_start: %ld\n",aln->d_start1);
    fprintf(fd,"; al_stop: %ld\n",aln->d_stop1);
    /* in the past, this al_display_start does not include sq1off */
    fprintf(fd,"; al_display_start: %ld\n",aln->l_offset+llsgn*ioff1+lfx0);

    have_res = 0;
    for (i=0; i<nc && seqc1[i]; i++) {
      if (!have_res && seqc1[i]==' ') fputc('-',fd);
      else if (seqc1[i]==' ') break;
      else {
	have_res = 1;
	fputc(seqc1[i],fd);
      }
      if (i%50 == 49) fputc('\n',fd);
    }
    if ((i-1)%50!=49 || seqc1[i-1]==' ') fputc('\n',fd);
#ifdef M10_CONS
    fprintf(fd,"; al_cons:\n");
    for (i=0,del0=0,id=ioff0; id-del0<aln->amax0 && i < nc; i++,id++) {
      if (seqc0[i] == '\0' || seqc1[i] == '\0') break;
      if (seqc0[i]=='-' || seqc0[i]==' ' || seqc0[i]=='\\') del0++;
      else if (seqc0[i]=='/') del0++;
      if (id-del0<aln->amin0) fputc(' ',fd);
      else if (seqc0[i]=='-'||seqc1[i]=='-') fputc('-',fd);
      else fputc(map_sym_p[seqca[i]],fd);

      if (i%50 == 49) fputc('\n',fd);
    }
    if ((i-1)%50!=49 || seqc1[i-1]==' ') fputc('\n',fd);
#endif
    return;
  }

  memset(line[0],' ',MAXOUT);
  memset(line[1],' ',MAXOUT);
  memset(line[2],' ',MAXOUT);

  /* cl0 indicates whether a coordinate should be printed over the first
     sequence; cl1 indicates a coordinate for the second;
  */

  ic = 0; del0=del1=0;
  for (il=0; il<(nc+l_llen-1)/l_llen; il++) {
    loff=il*l_llen;
    lend=min(l_llen,nc-loff);

    ll0 = NO; ll1 = NO;

    memset(cline[0],' ',MAXOUT+1);
    memset(cline[1],' ',MAXOUT+1);

    ic_save = ic;
    for (i=0; i<lend; i++, ic++,ioff0++,ioff1++) {
      cl0 =  cl1 = rl0 = rl1 = YES;
      if ((line[0][i]=seqc0[ic])=='-' || seqc0[ic]=='\\') {
	del0++; cl0=rl0=NO;
      }
      else if (seqc0[ic]=='/') {
	del0++; cl0=rl0=NO;
      }
      if ((line[2][i]=seqc1[ic])=='-' || seqc1[ic]=='\\') {
	del1++; cl1=rl1=NO;
      }
      else if (seqc1[ic]=='/') {
	del1++; cl1=rl1=NO;
      }

      if (seqc0[ic]==' ') {del0++; cl0=rl0=NO;}
      else ll0 = YES;
      if (seqc1[ic]==' ') {del1++; cl1=rl1=NO;}
      else ll1 = YES;

      /* the old version used qoffset, this version uses q_offset+q_off */
      qqoff = aln->q_offset + (long)qlsgn*ioff00 +
	(long)qlsgn*qlfact*(ioff0-del0-ioff00);
      if (cl0 && qqoff%10 == 9)  {
	sprintf(&clinep[0][i-qfxn],"%8ld",qqoff+1l);
	clinep[0][i+8-qfxn]=' ';
	rl0 = NO;
      }
      else if (cl0 && qqoff== -1) {
	sprintf(&clinep[0][i-qfxn],"%8ld",0l);
	clinep[0][i+8-qfxn]=' ';
	rl0 = NO;
      }
      else if (rl0 && (qqoff+1)%10 == 0) {
	sprintf(&clinep[0][i-qfxn],"%8ld",qqoff+1);
	clinep[0][i+8-qfxn]=' ';
      }
      
      /* the lloff coordinate of a residue is the sum of:
	 m_msg->sq1off-1	 - the user defined coordinate
	 aln->l_offset	- the offset into the library sequence
	 llsgn*ioff10	- the offset into the beginning of the alignment
	 		  (given in the "natural" coordinate system,
			   except for tfasta3 which provides context)
	 llsgn*llfact*(ioff1-del1-ioff10)
			- the position in the consensus aligment, -gaps
      */

      /* it seems like this should be done in cal_coord() */

      lloff = aln->l_offset + aln->frame +
	(long)llsgn*aln->llmult*ioff10 +
	(long)llsgn*llfact*(ioff1-del1-ioff10);

      if (cl1 && lloff%10 == 9)  {
	sprintf(&clinep[1][i-lfxn],"%8ld",lloff+1l);
	clinep[1][i+8-lfxn]=' ';
	rl1 = NO;
      }
      else if (cl1 && lloff== -1) {
	sprintf(&clinep[1][i],"%8ld",0l);
	clinep[1][i+8-lfxn]=' ';
	rl1 = NO;
      }
      else if (rl1 && (lloff+1)%10 == 0) {
	sprintf(&clinep[1][i-lfxn],"%8ld",lloff+1);
	clinep[1][i+8-lfxn]=' ';
      }

      line[1][i] = ' ';
      if (ioff0-del0 >= aln->amin0 && ioff0-del0 <= aln->amax0) {
	if (seqca[ic]==4) {line[1][i]=map_sym_p[4];}
	else if ((m_msg->markx&MX_ATYPE)==2) line[1][i]=line[2][i];
	else line[1][i] = map_sym_p[seqca[ic]];
      }
      else if ((m_msg->markx&MX_ATYPE)==2) line[1][i]=line[2][i];
    }

    if (m_msg->ann_flg) {
      for (ic=ic_save,i=0; i<lend; ic++,i++) {
	if (seqc0a && seqc0a[ic]!= ' ') clinep[0][i+7-qfxn] = seqc0a[ic];
	if (seqc1a && seqc1a[ic]!= ' ') clinep[1][i+7-qfxn] = seqc1a[ic];
      }
    }

    line[0][lend]=line[1][lend]=line[2][lend]=0;
    clinep[0][lend+7]=clinep[1][lend+7]=0;
    
    ll01 = ll0&&ll1;
    if ((m_msg->markx&MX_ATYPE)==2 && (!aln->showall || ll0)) ll1=0;
    fprintf(fd,"\n");
    if (ll0) fprintf(fd,"%s%s\n",blank,clinep[0]);
    if (ll0) fprintf(fd,afmt,name0,line[0]);
    if (ll01) fprintf(fd,afmt,name01,line[1]);
    if (ll1) fprintf(fd,afmt,name1,line[2]);
    if (ll1) fprintf(fd,"%s%s\n",blank,clinep[1]);
  }
}

static float gscale= -1.0;

void
disgraph(FILE *fd, int n0,int n1, float percent, int score,
	 int min0, int min1, int max0, int max1, long sq0off,
	 char *name0, char *name1, int nml,
	 int mlen, int markx)
{
  int i, gstart, gstop, gend;
  int llen;
  char line[MAXOUT+1];
  char afmt[16], afmtf[64];

  if (nml > 6) {
    sprintf(afmt,"%%-%ds",nml);
  }
  else {
    strncpy(afmt,"%-6s",sizeof(afmt));
  }
  strncpy(afmtf,afmt,sizeof(afmtf));
  strncat(afmtf," %4ld-%4ld:     %5.1f%%:%s:\n",sizeof(afmtf));

  llen = mlen - 10;
  memset(line,' ',llen);

  line[llen-1]='\0';
  if (gscale < 0.0) {
    gscale = (float)llen/(float)n0;
    if ((markx&MX_ATYPE) == 7 ) 
      fprintf(fd,afmtf,name0,sq0off,sq0off+n0-1,100.0,line);
  }

  gstart = (int)(gscale*(float)min0+0.5);
  gstop = (int)(gscale*(float)max0+0.5);
  gend = gstop+(int)(gscale*(float)(n1-max1));

  if (gstop >= llen) gstop = llen-1;
  if (gend >= llen) gend = llen-1;
  for (i=0; i<gstart; i++) line[i]=' ';
  for (; i<gstop; i++) line[i]='-';
  for (; i<llen; i++) line[i]=' ';

  line[gend]=':';
  line[llen]='\0';

  if (markx & MX_AMAP) {
    if ((markx & MX_ATYPE)==7) {	/* markx==4 - no alignment */
      strncpy(afmtf,afmt,sizeof(afmtf));
      strncat(afmtf," %4ld-%4ld:%4d %5.1f%%:%s\n",sizeof(afmtf));
      fprintf(fd,afmtf,name1,min0+sq0off,max0+sq0off-1,score,percent,line);
    }
    else {
      afmtf[0]='>';
      strncpy(&afmtf[1],afmt,sizeof(afmtf)-1);
      strncat(afmtf," %4ld-%4ld:%s\n",sizeof(afmtf));
      fprintf(fd,afmtf, name1,min0+sq0off,max0+sq0off-1,line);
    }
  }
}

void
aancpy(char *to, char *from, int count, struct pstruct *ppst)
{
  char *tp, *sq;
  int nsq;

  if (ppst->ext_sq_set) {
    nsq = ppst->nsqx;
    sq = ppst->sqx;
  }
  else {
    nsq = ppst->nsq;
    sq = ppst->sq;
  }

  tp=to;
  while (count-- && *from) {
    if (*from <= nsq) *tp++ = sq[*(from++)];
    else *tp++ = *from++;
  }
  *tp='\0';
}

void
r_memcpy(dest,src,cnt)
     char *dest, *src;
     int cnt;
{
  while (cnt--) *dest++ = *src++;
}

void
l_memcpy(dest,src,cnt)
     char *dest, *src;
     int cnt;
{
  dest = dest+cnt;
  src = src+cnt;
  while (cnt--) *--dest = *--src;
}

/* this routine now indexs from 1 (rather than 0) because sq starts
   with a 0 */

void
cal_coord(int n0, int n1,
	  long q_offset,
	  long l_offset,
	  struct a_struct *aln)
{
  int l_lsgn, q_lsgn, q_fx0, q_fxn, l_fx0, l_fxn;

  if (aln->qlrev == 1) {
    q_offset += n0;
    q_lsgn = -1;
    q_fx0 = 0;
    q_fxn = 1;
  }
  else {
    q_lsgn = 1;
    q_fx0 = 1;
    q_fxn = 0;
  }

  if (aln->llrev == 1) {
    l_offset += n1;
    l_lsgn = -1;
    l_fx0 = 0;
    l_fxn = 1;
  }
  else {
    l_lsgn = 1;
    l_fx0 = 1;
    l_fxn = 0;
  }
  aln->q_offset = q_offset;
  aln->l_offset = l_offset;
  aln->d_start0 = q_offset+q_lsgn*aln->amin0+q_fx0;
  aln->d_stop0  = q_offset+q_lsgn*aln->amax0+q_fxn;
  aln->d_start1 = l_offset+l_lsgn*aln->amin1*aln->llmult+l_fx0+aln->frame;
  aln->d_stop1  = l_offset+l_lsgn*aln->amax1*aln->llmult+l_fxn+aln->frame;
}
