/* lav2svg.c - produce postscript from lav output */

/* $Id: lav2svg.c 27 2008-06-30 16:27:31Z pearson $ */
/* $Revision: 64 $  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "defs.h"

void openplt(long, long, int, int, char *, char *);
void closeplt();
extern void drawdiag(long n0, long n1);
extern void closepl();
extern void move(int, int);
extern void cont(int, int);
extern void clsline();

int have_bits=0, have_zdb = 0;
long zdb_size = 1;
int g_n0, g_n1;

long pminx, pmaxx, pminy, pmaxy;
int max_x=540, max_y=540;

double fxscal, fyscal, fxoff, fyoff;
#define SX(x) (int)((double)(x)*fxscal+fxoff+6)
#define SY(y) (max_y + 24 - (int)((double)(y)*fyscal+fyoff))

void xaxis(long, int, char *);
void yaxis(long, int, char *);
void legend();
void linetype(int);
void opnline(int s, double bits);
void newline();
void clsline();
void move(int, int);
void draw(int, int);

double bit_to_E(double bits);

void get_str(FILE *, char *str, size_t len);
void get_str2(FILE *, char *, size_t, char *, size_t );
void get_seq_info(FILE *file, 
		  char *str0, size_t len0, int *n0_begin, int *n0_end,
		  char *str1, size_t len1, int *n1_begin, int *n1_end);
void del1(char *);
void do_alignment(FILE *, int, int);

/* for getopt() */
#ifdef UNIX
#include <unistd.h>
#else
extern int optind;
extern char *optarg;
#endif

int 
main(int argc, char **argv) {

  char line[MAX_STR];
  char pgm_desc[MAX_STR];
  char s_name0[MAX_STR], s_name1[MAX_STR];
  char s_desc0[MAX_STR], s_desc1[MAX_STR];
  int p0_beg, p1_beg, p0_end, p1_end;
  int open_plt = 0;
  int copt;

  /* check options */
  while ((copt = getopt(argc, argv, "BZ:")) != -1) {
    switch (copt) {
    case 'B':
      have_bits = 1;
      break;
    case 'Z':
      sscanf(optarg, "%ld", &zdb_size);
      have_zdb = 1;
      break;
    case '?':
    default:
      fprintf(stderr," usage -  ps_lav -B -Z db_size\n");
    }
  }

  while (fgets(line,sizeof(line), stdin)!=NULL) {
    if (line[0] == '#' || line[0] == '\n') continue;
    else {
      switch(line[0]) {
      case 'd':
	get_str(stdin, pgm_desc, sizeof(pgm_desc));
	break;
      case 'h':
	get_str2(stdin, 
		 s_desc0, sizeof(s_desc0),
		 s_desc1, sizeof(s_desc1));
	break;
      case 's':
	get_seq_info(stdin, 
		     s_name0, sizeof(s_name0), &p0_beg, &p0_end,
		     s_name1, sizeof(s_name1), &p1_beg, &p1_end);
	g_n0 = p0_end - p0_beg + 1;
	g_n1 = p1_end - p1_beg + 1;
	break;
      case 'a':
	if (!open_plt) {
	  openplt(g_n0, g_n1, p0_beg, p1_beg,  s_desc0, s_desc1);
	  if ((g_n0 == g_n1) && (p0_beg == p1_beg) && (p0_end == p0_end) &&
	      strcmp(s_name0, s_name1) == 0) {
	    drawdiag(p0_end-p0_beg + 1, p1_end - p1_beg + 1);
	  }
	  open_plt = 1;
	}

	do_alignment(stdin, p0_beg, p1_beg);
	break;
      }
    }
  }
  if (!open_plt) {
    openplt(g_n0, g_n1, p0_beg, p1_beg,  s_desc0, s_desc1);
    if ((g_n0 == g_n1) && (p0_beg == p1_beg) && (p0_end == p0_end) &&
	strcmp(s_name0, s_name1) == 0) {
      drawdiag(p0_end-p0_beg + 1, p1_end - p1_beg + 1);
    }
    open_plt = 1;
  }
  closeplt();
  exit(0);
}

void
get_str(FILE *file, char *str, size_t len) {

  char line[MAX_STR], *bp, *bp1;
  int have_quote = 0;

  str[0] = '\0';

  while (fgets(line, sizeof(line), file)!=NULL) {
    if (line[0] == '}') return;
    if ((bp = strchr(line,'\n'))!=NULL) *bp = '\0';
    if (have_quote == 0 && (bp=strchr(line,'\"'))!=NULL) {
      have_quote == 1;
      if ((bp1 = strchr(bp+1, '\"'))!=NULL) {
	*bp1 = '\0';
	have_quote = 2;
      }
      strncat(str, bp+1, len-1);
      len =- strlen(bp+1);
    }
    else if (have_quote == 1) {
      if ((bp = strchr(line, '\"'))!=NULL) *bp = '\0';
      strncat(str, line, len-1);
      len =- strlen(line);
    }
  }
}

void
get_str2(FILE *file, 
	 char *str0, size_t len0,
	 char *str1, size_t len1) {

  char line[MAX_STR], *bp, *bp1;
  int have_quote0 = 0;
  int have_quote1 = 0;

  str0[0] = str1[0] = '\0';

  while (fgets(line, sizeof(line), file)!=NULL) {
    if (line[0] == '}') return;
    if ((bp = strchr(line,'\n'))!=NULL) *bp = '\0';
    if (have_quote0 == 0 && (bp=strchr(line,'\"'))!=NULL) {
      have_quote0 == 1;
      if ((bp1 = strchr(bp+1, '\"'))!=NULL) {
	*bp1 = '\0';
	have_quote0 = 2;
      }
      strncat(str0, bp+1, len0-1);
      len0 =- strlen(bp+1);
    }
    else if (have_quote0 == 1) {
      if ((bp = strchr(line, '\"'))!=NULL) *bp = '\0';
      strncat(str0, line, len0-1);
      len0 =- strlen(line);
    }
    else if (have_quote1 == 0 && (bp=strchr(line,'\"'))!=NULL) {
      have_quote1 == 1;
      if ((bp1 = strchr(bp+1, '\"'))!=NULL) {
	*bp1 = '\0';
	have_quote1 = 2;
      }
      strncat(str1, bp+1, len1-1);
      len1 =- strlen(bp+1);
    }
    else if (have_quote1 == 1) {
      if ((bp = strchr(line, '\"'))!=NULL) *bp = '\0';
      strncat(str1, line, len1-1);
      len1 =- strlen(line);
    }
  }
}

void
get_seq_info(FILE *file, 
	     char *str0, size_t len0, int *n0_begin, int *n0_end,
	     char *str1, size_t len1, int *n1_begin, int *n1_end) {

  char line[MAX_STR], *bp;
  int have_quote0 = 0;
  int have_quote1 = 0;

  str0[0] = str1[0] = '\0';

  fgets(line, sizeof(line), file);
  if (line[0] == '}') return;
  sscanf(line, "%s %d %d", str0, n0_begin, n0_end);
  fgets(line, sizeof(line), file);
  if (line[0] == '}') return;
  sscanf(line, "%s %d %d", str1, n1_begin, n1_end);

  if ((bp = strchr(str0+1,'\"'))!=NULL) *bp = 0;
  if ((bp = strchr(str1+1,'\"'))!=NULL) *bp = 0;
  if (str0[0] == '\"') {del1(str0);}
  if (str1[0] == '\"') {del1(str1);}

  fgets(line, sizeof(line), file);	/* get the last } */
}

void
del1(char *str) {
  char *str_p;

  str_p = str++;

  while (*str++) {*str_p++ = *str++;}

  *str_p = '\0';
}

void
do_alignment(FILE *file, int p0_beg, int p1_beg) {

  char line[MAX_STR], *bp;
  int score, s0_beg, s0_end, s1_beg, s1_end, percent;
  int have_line = 0;
  double bits;

  while (fgets(line, sizeof(line), file)!=NULL) {
    if ( strchr(line,'}') != NULL) {
      clsline();
      return;
    }
    if ((bp=strchr(line,'s')) != NULL)  sscanf(bp+1, "%d %lf", &score, &bits);
    else if ((bp=strchr(line,'b')) != NULL) sscanf(bp+1, "%d %d", &s0_beg, &s1_beg);
    else if ((bp=strchr(line, 'e')) != NULL) sscanf(bp+1, "%d %d", &s0_end, &s1_end);
    else if ((bp=strchr(line, 'l')) != NULL) {
      sscanf(bp+1, "%d %d %d %d %d",
	     &s0_beg, &s1_beg, &s0_end, &s1_end, &percent);
      if (have_line) {
	draw(SX(s0_beg-p0_beg+1), SY(s1_beg-p1_beg+1));
	draw(SX(s0_end-p0_beg+1), SY(s1_end-p1_beg+1));
      }
      else {
	opnline(score, bits);
	move(SX(s0_beg - p0_beg + 1), SY(s1_beg - p1_beg + 1));
	draw(SX(s0_end - p0_beg + 1), SY(s1_end - p1_beg + 1));
	have_line = 1;
      }
    }
  }
}

#include <math.h>

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

/* produce e_val from bit score */

double
bit_to_E (double bit)
{
  double a_n0, a_n1, p_val;

  a_n0 = (double)g_n0;
  a_n1 = (double)g_n1;

  p_val = a_n0 * a_n1 / pow(2.0, bit);
  if (p_val > 0.01) p_val = 1.0 - exp(-p_val);

  return (double)zdb_size * p_val;
}

/* functions/variables for postscript plots */

/* black blue cyan green lt_green */
char *line_color[]={"black","blue","brown","green","red"};

int *linarr;
int nlinarr=5;

char lvstr[MAX_STR];

double elinval[4]={1e-4,1e-2,1.0,100.0};
double blinval[4]={40.0,30.0,20.0,10.0};
int ilinval[4]={200,100,50,25};

#define DIAG 1
#define INS0 2
#define INS1 4

void
openplt(long n0, long n1, int sq0off, int sq1off, 
	char *xtitle, char *ytitle)
{
  char *getenv(), *sptr;
  time_t tt;
  int tick = 6;

  tt = time(NULL);

  if (strlen(lvstr)>0) {
    sscanf(lvstr,"%lg %lg %lg",&elinval[0],&elinval[1],&elinval[2]);
  }
  else if ((sptr=getenv("LINEVAL"))!=NULL && strlen(sptr)>0) {
    sscanf(sptr,"%lg %lg %lg",&elinval[0],&elinval[1],&elinval[2]);
  }
	
  printf("<?xml version=\"1.0\" standalone=\"no\"?>\n");
  printf("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n");
  printf("\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\n");

  printf("<svg width=\"564\" height=\"612\" version=\"1.1\"\n");
  printf("xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n\n");

  pmaxx = n0;
  pmaxy = n1;

  fxscal = (double)(max_x-1)/(double)(n1);
  fyscal = (double)(max_y-1)/(double)(n0);

  if (fxscal > fyscal) fxscal = fyscal;
  else fyscal = fxscal;

  if (fyscal * n0 < (double)max_y/5.0) 
    fyscal = (double)(max_y-1)/((double)(n0)*5.0);

  fxscal *= 0.9; fxoff = (double)(max_x-1)/11.0;
  fyscal *= 0.9; fyoff = (double)(max_y-1)/11.0;

  /*   printf("currentlinewidth 1.5 mul setlinewidth\n"); */

  printf("<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\"\n",
	 SX(0),SY(n1+1), SX(n0+1)-SX(0), SY(0) - SY(n1+1));
  printf("stroke=\"black\" stroke-width=\"2.0\" fill=\"none\" />\n");

  xaxis(n0,sq1off, xtitle);
  yaxis(n1,sq0off, ytitle);
  legend();
}
	
void
drawdiag(n0,n1)
	long n0, n1;
{
  
	/* 	printf("currentlinewidth 1.5 mul setlinewidth\n"); */
	newline("stroke=\"black\" stroke-width=\"1.5\"");
	move(SX(0),SY(0));
	draw(SX(n0+1),SY(n1+1));
	clsline(n0,n1,10000);
}

/* tick array - values */
int tarr[] = {10,20,50,100,200,500,1000,2000,5000};
int ntarr = sizeof(tarr);

void
xaxis(long n, int offset, char *title)
{
  int i, jm, tick;
  long js, jo, jl;
  char numstr[20],*bp;

  tick = 6;

  /* search for the correct increment for the tick array */
  for (i=0; i<ntarr; i++) {
    /* seek to divide into 20 or fewer divisions */
    if ((jm = n/tarr[i])<21) goto found;
  }
  jm = n/5000l;
  i=ntarr-1;
 found:
  /* js is the start of the value - modify to accomodate offset */
  js = tarr[i];

  /* jo is the offset */
  jo = offset%tarr[i];	/* figure out offset in tarr[i] increments */

  /* jl is the label */
  jl = offset/tarr[i];	/* figure out offset in tarr[i] increments */
  jl *= tarr[i];

  newline("stroke=\"black\" stroke-width=\"1.5\"");
  for (i=1; i<=jm; i++) {
    move(SX((long)i*js - jo),SY(0));
    draw(SX((long)i*js - jo),SY(0)+tick);
  }
  clsline(n,n,10000);

  sprintf(numstr,"%ld",js + jl );
  printf("<text x=\"%d\" y=\"%d\" text-anchor=\"middle\">%s</text>\n",SX((long)js-jo),SY(0)+tick+16,numstr);

  sprintf(numstr,"%ld",jm*js+jl);
  printf("<text x=\"%d\" y=\"%d\" text-anchor=\"middle\">%s</text>\n",SX((long)jm*js-jo),SY(0)+tick+16,numstr);

  printf("<text x=\"%d\" y=\"%d\" text-anchor=\"middle\">%s</text>\n",SX(n/2),SY(0)-tick+42, title);
}
		
void
yaxis(long n, int offset, char *title)
{
  int i, jm, tick;
  long js, jo, jl;
  char numstr[20],*bp;
	
  tick = 6;

  for (i=0; i<ntarr; i++) {
    if ((jm = n/tarr[i])<21) goto found;
  }
  jm = n/5000l;
  i=ntarr-1;

 found:
  js = (long)tarr[i];

  /* jo is the offset */
  jo = offset%tarr[i];	/* figure out offset in tarr[i] increments */
  /* jl is the label */
  jl = offset/tarr[i];	/* figure out offset in tarr[i] increments */
  jl *= tarr[i];

  newline("stroke=\"black\" stroke-width=\"1.5\"");
  for (i=1; i<=jm; i++) {
    move(SX(0),SY((long)i*js-jo));
    draw(SX(0)-tick,SY((long)i*js-jo));
  }
  clsline(n,n,10000);

  sprintf(numstr,"%ld",js+jl);
  printf("<text x=\"%d\" y=\"%d\" text-anchor=\"end\">%s</text>\n",SX(0)-tick-4,SY(js-jo)+4,numstr);

  sprintf(numstr,"%ld",(long)jm*js+jl);
  printf("<text x=\"%d\" y=\"%d\" text-anchor=\"end\">%s</text>\n",SX(0)-tick-4,SY((long)jm*js-jo)+4,numstr);

  /* make a path for the label */

  printf("<g transform=\"rotate(-90, 142, 142)\">\n");
  printf("<text x=\"0\" y=\"38\" text-anchor=\"middle\">%s</text>\n",title);
  printf("</g>\n");
}

void
legend()
{
  int i, last, del;
  int ixp, iyp;
  char numstr[10];
  char optstr[128];
  int xpos[]={144,144,288,288,432};
  int ypos[]={36,18,36,18,27};

  int tick = 6;

  if (have_zdb || have_bits) last = 5;
  else last = 4;

  if (have_zdb)  printf("<text x=\"%d\" y=\"%d\">E(): </text>",54,max_y + 66 - 27);
  else if (have_bits) printf("<text x=\"%d\" y=\"%d\">bits: </text>",54,max_y + 66 - 27);

  del = 10;
  for (i=0; i<last ; i++) {
    
    sprintf(optstr,"stroke-width=\"1.5\" stroke=\"%s\"",line_color[i]);
    newline(optstr);
    /*    linetype(i);*/
    move(xpos[i]-12,max_y + 66 - ypos[i]);
    draw(xpos[i]+48,max_y + 66 - ypos[i]);
    clsline(1000,1000,10000);

    if (have_zdb) {
      if (i==4) sprintf(numstr,"&gt;%.1lg",elinval[3]);
      else sprintf(numstr,"&lt;%.1lg",elinval[i]);
    }
    else if (have_bits) {
      if (i==4) sprintf(numstr,"&lt;%.1lf",blinval[3]);
      else sprintf(numstr,"&gt;=%.1lf",blinval[i]);
    }
    else {
      if (i==3) sprintf(numstr,"&lt;%d",ilinval[3]);
      else sprintf(numstr,"&gt;%d",ilinval[i]);
    }

    printf("<text align=\"center\" x=\"%d\" y=\"%d\">%s</text>\n",xpos[i]+54, max_y + 66 - ypos[i],numstr);
  }
}

void
linetype(type)
     int type;
{
  printf(" stroke=\"%s\"",line_color[type]);
}

void
closeplt()
{
  printf("</svg>\n");
}

void
opnline(int s, double bits)
{
  double e_val;



  if (have_zdb) {
    e_val = bit_to_E(bits);
    printf("<!-- score: %d; bits: %.1g; E(): %.1g -->\n",s,bits,e_val);
    printf("<path ");
    if (e_val < elinval[0]) linetype(0);
    else if (e_val < elinval[1]) linetype(1);
    else if (e_val < elinval[2]) linetype(2);
    else if (e_val < elinval[3]) linetype(3);
    else linetype(4);
  }
  else if (have_bits) {
    printf("<!-- score: %d; bits: %.1g -->\n",s,bits);
    printf("<path ");
    if (bits >= blinval[0]) linetype(0);
    else if (bits >= blinval[1]) linetype(1);
    else if (bits >= blinval[2]) linetype(2);
    else if (bits >= blinval[3]) linetype(3);
    else linetype(4);
  }
  else {
    printf("<!-- score: %d -->\n",s);
    printf("<path ");
    if (s > ilinval[0]) linetype(0);
    else if (s> ilinval[1]) linetype(1);
    else if (s> ilinval[2]) linetype(2);
    else linetype(3);
  }

  printf(" d=\"");
}

void
newline(char *options)
{
  if (options != NULL)  {
      printf("<path %s d=\"",options);
  }
  else printf("<path stroke=\"black\" d=\"");
}

void
clsline(x,y,s)
     long x, y;
     int s;
{
  printf("\" fill=\"none\" />\n");
}

void
move(x,y)
     int x, y;
{
  printf(" M %d %d",x,y);
}

void
draw(x,y)
	int x, y;
{
  printf(" L %d %d",x,y);
}

void
draw_str(str)
     char *str;
{
  char *bp;

  for (bp = strchr(str,'('); (bp!=NULL); bp = strchr(bp+1,'(')) *bp=' ';
  for (bp = strchr(str,')'); (bp!=NULL); bp = strchr(bp+1,')')) *bp=' ';

  printf("(%s) show\n",str);
}

void
draw_sstr(str)
     char *str;
{
  char *bp;

  printf("(%s) show\n",str);
}

void cal_coord(int n0, int n1, 
	       long *a_start0, long *a_stop0, long *a_start1, long *a_stop1 )
{}
