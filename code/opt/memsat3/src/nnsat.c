/*********************************************************
 *    MEMSAT - MEMbrane protein Structure And Topology   *
 * Integral Membrane Protein Topology Prediction Program *
 *           Copyright (C) 2006 David T. Jones           *
 *********************************************************/
/*
  This program is copyright and may not be distributed without
  permission of the author unless specifically permitted under
  the terms of the license agreement.

  THIS SOFTWARE MAY ONLY BE USED FOR NON-COMMERCIAL PURPOSES. PLEASE CONTACT
  THE AUTHOR IF YOU REQUIRE A LICENSE FOR COMMERCIAL USE.

*/

/* Version 3.5 */

/* Program created: August 6th 1992 */
/* Last edit: January 30th 2008 */

/*
 * Description: This program uses dynamic programming to
 * predict the secondary structure and topology of integral
 * membrane sequences based on multiple sequence profiles.
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

#define FALSE 0
#define TRUE 1

#define BIG (100000000)
#define VBIG (1e32F)

#define MAXSEQLEN 2000

#define noDIAGRAM

#define noTRAPEZOID

/* Maximum number of predicted helices */
int MAXNHEL = 20;

/* Minimum length of 'loop' */
int MINLLEN = 2;

/* Minimum length of helix */
int MINHLEN = 19;

/* Maximum length of helix */
int MAXHLEN = 25;

/* Maximum length of "topogenic" loop */
int LIMITLOOP = 60;

/* Minimum score for helix */
int MINHSC = 4000;


int quietflg;

#define SQR(x) ((x)*(x))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))
#define CH malloc_verify(), printf("Heap OK at line : %d.\n",__LINE__);

const float ZERO = 0.0, ONE = 1.0, TWO = 2.0, THREE = 3.0, FOUR = 4.0;

char seq[MAXSEQLEN], sstruc[MAXSEQLEN];
int tpltsc[MAXSEQLEN][5];

const char *rnames[] =
{
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "???"
};

enum aacodes
{
    ALA, ARG, ASN, ASP, CYS,
    GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
    UNK, GAP
};

struct helix
{
    int length, start;
} structure[50], bestst[50];

int mat[50][MAXSEQLEN];
short length[50][MAXSEQLEN], path[50][MAXSEQLEN];

int spscmat[20][33], spfilt = FALSE;

int seqlen, endsig;
float sigscore[50];

const char *rescodes = "ARNDCQEGHILKMFPSTWYVX";

int s_inside[20], s_outside[20], s_insideh[20], s_outsideh[20], s_xhelix[20];
int m_inside[20], m_outside[20], m_insideh[20], m_outsideh[20], m_xhelix[20];

FILE *dfp;

/* Dump a rude message to standard error and exit */
void
  fail(char *errstr)
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
}

/* Return score for given predicted helix */
int
  helscore(int start, int length, int inflg, int nhelix)
{
    int aa, i, lby4;
    int outer=0, tot=0;

    lby4 = length / 4;

    for (i = 0; i < lby4; i++)
    {
	if (inflg)
	{
	    tot += tpltsc[start + i][3];
	    tot += tpltsc[start + length - i - 1][2];
	}
	else
	{
	    tot += tpltsc[start + i][2];
	    tot += tpltsc[start + length - i - 1][3];
	}
    }

    for (i = lby4; i <= length - lby4 - 1; i++)
	tot += tpltsc[i + start][4];

    return tot;
}

/* Print string representation of predicted structure */
char *
  str_model(int inflg, int nhelix, struct helix structure[50])
{
    int i, j, lby4;

    for (i = 0; i < structure[0].start; i++)
	sstruc[i] = inflg ? '+' : '-';

    if (endsig >6 && endsig < structure[0].start && !inflg)
	for (i=0; i<endsig; i++)
	    sstruc[i]= 'S';

    for (i = 0; i < nhelix; i++)
    {
	if (i)
	    for (j = structure[i - 1].start + structure[i - 1].length; j < structure[i].start; j++)
		sstruc[j] = inflg ? '+' : '-';
	for (j = 0; j < structure[i].length; j++)
	    sstruc[j + structure[i].start] = 'X';
	lby4 = structure[i].length / 4;
	for (j = 0; j < lby4; j++)
	{
	    sstruc[structure[i].start + j] = inflg ? 'I' : 'O';
	    sstruc[structure[i].start + structure[i].length - j - 1] = inflg ? 'O' : 'I';
	}
	for (j = lby4; j <= structure[i].length - lby4 - 1; j++)
	    sstruc[structure[i].start + j] = 'X';
	inflg = !inflg;
    }
    for (i = structure[nhelix - 1].start + structure[nhelix - 1].length; i < seqlen; i++)
	sstruc[i] = inflg ? '+' : '-';

    sstruc[seqlen] = '\0';

    return sstruc;
}

/* Calculate score/path matrix by dynamic programming */
void
  calcmat(int inflg, int nhelix)
{
    int aa, i, j, h, k, l, maxj, maxl;
    int maxsc, hsc, lsc, intog;

    intog = inflg ^ (nhelix & 1);
    for (h = nhelix - 1; h >= 0; h--, intog = !intog)
    {
	for (i = 0; i < seqlen; i++)
	    mat[h][i] = -BIG;
	for (i = seqlen - MINHLEN - MINLLEN; i >= MINLLEN; i--)
	{
	    maxsc = -BIG;
	    maxj = maxl = 0;
	    for (l = MINHLEN; l <= MAXHLEN; l++)
		if (i + l + MINLLEN <= seqlen)
		{
		    hsc = helscore(i, l, !intog, nhelix);
		    if (h == nhelix - 1)
		    {
			lsc = 0;
			for (k = i + l; k < seqlen; k++)
			{
			    if (intog)
				lsc += tpltsc[k][0];
			    else
				lsc += tpltsc[k][1];
			}

			if (k - i - l + 1 > LIMITLOOP)
			    lsc = 0;

			if (hsc + lsc > maxsc)
			{
			    maxsc = hsc + lsc;
			    maxl = l;
			}
		    }
		    else
		    {
			/* Calculate initial loop score */
			lsc = 0;
			for (k = i + l; k < i + l + MINLLEN - 1; k++)
			    lsc += intog ? tpltsc[k][0] : tpltsc[k][1];
			for (j = i + l + MINLLEN; j < seqlen - MINHLEN - MINLLEN; j++)
			{
			    /* Add extra loop residue score */
			    lsc += intog ? tpltsc[j - 1][0] : tpltsc[j - 1][1];

			    if (j - i - l > LIMITLOOP)
				lsc = 0;

			    if (hsc + lsc + mat[h + 1][j] > maxsc)
			    {
				maxsc = hsc + lsc + mat[h + 1][j];
				maxl = l;
				maxj = j;
			    }
			}
		    }
		}
	    mat[h][i] = maxsc;
	    length[h][i] = maxl;
	    path[h][i] = maxj;
	}
    }
}

/* Handle first loop region */
void firstlp(int starthel, int inflg, int nhelix)
{
    int aa, i;
    int lsc;

    lsc = 0;
    for (i = 0; i < seqlen - (nhelix - starthel) * (MINHLEN + MINLLEN); i++)
    {
	lsc += inflg ? tpltsc[i][0] : tpltsc[i][1];

    	if (i+1 > LIMITLOOP)
    	    lsc = 0;

    	mat[starthel][i+1] += lsc;
    }
}

/* Trace back highest scoring path through matrix */
int
  trace_back(int starthel, int inflg, int nhelix)
{
    int i, h, res, intog = inflg, maxsc = -BIG, hsc, weak_h=FALSE;

    for (i = MINLLEN; i < seqlen - (nhelix - starthel) * (MINHLEN + MINLLEN); i++)
	if (mat[starthel][i] > maxsc)
	{
	    maxsc = mat[starthel][i];
	    res = i;
	}

    for (h = 0; h < nhelix-starthel; h++, intog = !intog)
    {
	if (quietflg)
	{
	    hsc = helscore(res, length[h+starthel][res], intog, nhelix);
	}
	else
	{
	    if (intog)
		printf("Helix %d from %3d (in) to %3d (out) : %ld\n", h + 1, res+1, res + length[h+starthel][res], hsc = helscore(res, length[h+starthel][res], intog, nhelix));
	    else
		printf("Helix %d from %3d (out) to %3d (in) : %ld\n", h + 1, res+1, res + length[h+starthel][res], hsc = helscore(res, length[h+starthel][res], intog, nhelix));
	}

	if (dfp)
	    fprintf(dfp, "%d %d\n", res, res + length[h+starthel][res]);
	if (hsc < MINHSC)
	    weak_h = TRUE;
	structure[h].start = res;
	structure[h].length = length[h+starthel][res];
	res = path[h+starthel][res];
    }

    if (!quietflg)
	puts(str_model(inflg, nhelix-starthel, structure));
    else
	str_model(inflg, nhelix-starthel, structure);

    return (nhelix == 1) ? maxsc : (weak_h ? -BIG : maxsc);
}

/* Convert AA letter to numeric code (0-22) */
int
  aanum(int ch)
{
    static int aacvs[] =
    {
	999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 22, 11, 10, 12, 2,
	22, 14, 5, 1, 15, 16, 22, 19, 17, 22, 18, 21
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 22);
}

/* Read PSI AA frequency data */
int             getnnvec(char *seq, FILE * lfil)
{
    int             i, j, naa, transtab[20];
    float in, out, hmid, sig;
    char            buf[256], *p;

    naa = 0;
    while (!feof(lfil))
    {
	if (!fgets(buf, 256, lfil))
	    break;
	seq[naa] = aanum(buf[0]);
	if (sscanf(buf + 3, "%f%f%f%f", &in, &out, &hmid, &sig) != 4)
	    fail("Bad NN format!");

	hmid -= 0.5;
	sig -= 0.15;
	
	if (naa < 30)
	{
	    tpltsc[naa][0] = ((in - out) - hmid) * 500 - 1000 * sig;
	    tpltsc[naa][1] = ((out - in) - hmid) * 500 + 1000 * sig;
	    tpltsc[naa][2] = (hmid + (out - in)) * 2000;
	    tpltsc[naa][3] = (hmid + (in - out)) * 2000;
	    tpltsc[naa][4] = hmid * 2000;
	}
	else
	{
	    tpltsc[naa][0] = ((in - out) - hmid) * 500;
	    tpltsc[naa][1] = ((out - in) - hmid) * 500;
	    tpltsc[naa][2] = (hmid + (out - in)) * 2000;
	    tpltsc[naa][3] = (hmid + (in - out)) * 2000;
	    tpltsc[naa][4] = hmid * 2000;
	}

	if (naa < 50)
	    sigscore[naa] = sig;

	naa++;
    }

    return naa;
}

void
usage(char *cmdname)
{
    puts("+-------------------------------------------------------------+");
    puts("| MEMSAT - MEMbrane Structure And Topology prediction program |");
    puts("| Written by David T. Jones                                   |");
    puts("| Department of Computer Science                              |");
    puts("| University College London                                   |");
    puts("| Gower Street                                                |");
    puts("| London                                                      |");
    puts("| WC1E 6BT                                                    |");
    puts("+-------------------------------------------------------------+");
    puts("\nVERSION 3.5 (ACADEMIC VERSION)\n");
    printf("usage: %s {options} profile\n", cmdname);
    puts("\nwhere options are as follows:");
    puts("-hnnn = set maximum number of helices");
    puts("-dnnn = set minimum sequence length per helix");
    puts("-lnnn = set minimum loop length");
    puts("-mnnn = set minimum length of helix");
    puts("-xnnn = set maximum length of helix");
    puts("-snnn = set helix score cutoff");
    exit(1);
}

main(int argc, char **argv)
{
    int aa, n, b, nb, h, i, j, t, scmax, score, nhelix;
    FILE *ifp;
    char name[512], dirname[512], fname[512], *cmdstr, *pred, *cp;
    int finalsc[50][2], inflg, maxsc = -BIG, maxt, maxnh;

    for (cmdstr = *argv++, argc--; argc && **argv == '-'; argv++, argc--)
	switch (*(*argv + 1))
	{
 	case 'q':
	    quietflg = TRUE;
	    break;
 	case 'h':
	    MAXNHEL = atoi(*argv + 2);
	    break;
 	case 'l':
	    MINLLEN = atoi(*argv + 2);
	    break;
 	case 'z':
	    LIMITLOOP = atoi(*argv + 2);
	    break;
 	case 'm':
	    MINHLEN = atoi(*argv + 2);
	    break;
 	case 'x':
	    MAXHLEN = atoi(*argv + 2);
	    break;
 	case 's':
	    MINHSC = atoi(*argv + 2);
	    break;
	default:
	    usage(cmdstr);
	}

    if (argc != 1)
	usage(cmdstr);

    ifp = fopen(argv[0], "r");
    if (!ifp)
	fail("Unable to open neural network input file!");
    seqlen = getnnvec(seq, ifp);

    if (!quietflg)
    {
	puts("MEMSAT Version 3.5 (ACADEMIC VERSION) - Copyright 2006 D.T. Jones\n");
	puts("\nCOMMERCIAL USE OF THIS PROGRAM IS FORBIDDEN\n");

	printf("%s\n%d residues read from file.\n\n", argv[0], seqlen);
    }

    if (seqlen < MINHLEN+2*MINLLEN)
	fail("Target sequence too short!");

#ifdef DIAGRAM
    dfp = fopen("memdiag.out", "w");
#endif

    endsig = 1;

    for (i=1; i<MIN(seqlen, 50); i++)
    {
	sigscore[i] += sigscore[i-1];
/*	printf("%d %f\n", i, sigscore[i]); */
	if (sigscore[i] > sigscore[endsig])
	    endsig = i;
    }
    
    if (endsig > 6)
	printf("Potential signal peptide at N-terminus: length %d, score %f.\n\n", endsig+1, sigscore[endsig]);

    nhelix = 1;
    calcmat(inflg = TRUE, nhelix);
    firstlp(0, inflg, nhelix);
    score = trace_back(0, inflg, nhelix);
    if (score > maxsc)
    {
    	memcpy(bestst, structure, sizeof(structure));
    	maxsc = score;
    	maxt = inflg;
    	maxnh = 1;
    }
    finalsc[0][inflg] = score;
    if (!quietflg)
	printf("Score = %f\n", score/1000.0);
    calcmat(inflg = FALSE, nhelix);
    firstlp(0, inflg, nhelix);
    score = trace_back(0, inflg, nhelix);
    if (score > maxsc)
    {
    	memcpy(bestst, structure, sizeof(structure));
    	maxsc = score;
    	maxt = inflg;
    	maxnh = 1;
    }
    finalsc[0][inflg] = score;
    if (!quietflg)
	printf("Score = %f\n", score/1000.0);
    nhelix = MAX(1, MIN(MAXNHEL, (seqlen - MINLLEN) / (MINHLEN + MINLLEN)));
    calcmat(inflg = TRUE, nhelix);
    for (h=0; h<nhelix-1; h++, inflg = !inflg)
    {
    	firstlp(h, inflg, nhelix);
    	score = trace_back(h, inflg, nhelix);
	if (score > maxsc)
	{
	    memcpy(bestst, structure, sizeof(structure));
	    maxsc = score;
	    maxt = inflg;
	    maxnh = nhelix-h;
	}
	finalsc[nhelix-h-1][inflg] = score;
	if (!quietflg)
	    printf("Score = %f\n", score/1000.0);
    }
    calcmat(inflg = FALSE, nhelix);
    for (h=0; h<nhelix-1; h++, inflg = !inflg)
    {
    	firstlp(h, inflg, nhelix);
    	score = trace_back(h, inflg, nhelix);
	if (score > maxsc)
	{
	    memcpy(bestst, structure, sizeof(structure));
	    maxsc = score;
	    maxt = inflg;
	    maxnh = nhelix-h;
	}
	finalsc[nhelix-h-1][inflg] = score;
	if (!quietflg)
	    printf("Score = %f\n", score/1000.0);
    }

#if 0
	if (dfp)
	    fprintf(dfp, "%g %d %d 1\n", insc / 1000.0, seqlen, nhelix);
#endif

    printf("Summary of topology analysis:\n");
    for (nhelix = 1; nhelix <= MAX(1, MIN(MAXNHEL, (seqlen - MINLLEN) / (MINHLEN + MINLLEN))); nhelix++)
    {
	if (nhelix == 1)
	    printf("%2d helix   (+) : Score = %g\n", nhelix, finalsc[nhelix-1][1] / 1000.0);
	else
	    printf("%2d helices (+) : Score = %g\n", nhelix, finalsc[nhelix-1][1] / 1000.0);
	if (nhelix == 1)
	    printf("%2d helix   (-) : Score = %g\n", nhelix, finalsc[nhelix-1][0] / 1000.0);
	else
	    printf("%2d helices (-) : Score = %g\n", nhelix, finalsc[nhelix-1][0] / 1000.0);
    }

    puts("\nFINAL PREDICTION");
    puts("================");

    t = maxt;
    printf("1: %s ", t ? "(in)" : "(out)");
    printf("%d-%d\t(%.2f)\n", bestst[0].start+1, bestst[0].start+1 + bestst[0].length - 1, helscore(bestst[0].start, bestst[0].length, t, maxnh)/1000.0);
    t = !t;
    for (i = 1; i < maxnh; i++, t = !t)
	printf("%d: %d-%d\t(%.2f)\n", i + 1, bestst[i].start+1, bestst[i].start+1 + bestst[i].length - 1, helscore(bestst[i].start, bestst[i].length, t, maxnh)/1000.0);

    putchar('\n');
    pred = str_model(maxt, maxnh, bestst);

    nb = seqlen / 60 + 1;
    j = 1;
    for (b = 0; b < nb; b++)
    {
	for (i = 0; i < 58; i++)
	{
	    if (b * 60 + i + 3 > seqlen)
		break;
	    j = b * 60 + i + 3;
	    if (!(j % 10))
	    {
		printf("%3d", j);
		i += 2;
	    }
	    else
		printf(" ");
	}
	putchar('\n');

	for (i = 0; i < 60; i++)
	{
	    if (b * 60 + i >= seqlen)
		break;
	    j = b * 60 + i + 1;
	    putchar(pred[j - 1]);
	}

	putchar('\n');

	for (i = 0; i < 60; i++)
	{
	    if (b * 60 + i >= seqlen)
		break;
	    j = b * 60 + i + 1;
	    putchar(rescodes[seq[j - 1]]);
	}

	putchar('\n');
	putchar('\n');
    }

    if (dfp)
	fclose(dfp);

    return 0;
}
