/* Neural Network Identification of Transmembrane Topology (MEMSAT3) */
/* Copyright (C) 2006 David T. Jones - Last Edited : Oct 2006 */

/*
  THIS SOFTWARE MAY ONLY BE USED FOR NON-COMMERCIAL PURPOSES. PLEASE CONTACT
  THE AUTHOR IF YOU REQUIRE A LICENSE FOR COMMERCIAL USE.
*/

/* Prediction Module */

/*
 * Description: This program provides an implementation of a neural network
 * containing one hidden layer, which uses the generalized backpropagation
 * delta rule for learning.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#include "sigmem_net.h"

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

enum
{
    FALSE, TRUE
};

#define MAXSEQLEN 50000

#define SQR(x) ((x)*(x))

#define REAL float

void           *
                calloc(), *malloc();

char           *wtfnm;

int             profile[MAXSEQLEN][20];

int             nwtsum, fwt_to[TOTAL], lwt_to[TOTAL];
REAL            blrate = ILRATE, lrate = ILRATE, alpha = IALPHA;
REAL            activation[TOTAL], bias[TOTAL], netinput[TOTAL], *weight[TOTAL];
REAL            bias[TOTAL];

int             nhelix, nsheet, ncoil, restot, nblocks;

int             profile[MAXSEQLEN][20];
char            seq[MAXSEQLEN], ssstruc[MAXSEQLEN], ssstrel[MAXSEQLEN];
int             seqlen;

char           *rnames[] =
{
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "UNK"
};

enum aacodes
{
    ALA, ARG, ASN, ASP, CYS,
    GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
    UNK
};

/* Generate a small random weight (+/- RWMAX) */
#define rndwt() ((ran0() - 0.5) * RWMAX)

/* Generate random int 0..X-1 */
#define rndint(x) ((int)((REAL)(x)*ran0()))

/* Generate a small 'noise' value (0-0.05) */
#define noise() (0.05*ran0())

void
                err(s)
     char           *s;
{
    fprintf(stderr, "%s\n", s);
}

void
                fail(s)
     char           *s;
{
    fprintf(stderr, "%s\n", s);
    exit(1);
}

long            modulus = 2147483399, seed = 1;

/* Randomize seed */
void
                randomize()
{
    long            time();

    seed = time(NULL) % (modulus + 1);
}

/* Get random long int (ref. CACM, June 1988, Page 742) */
long
                randl()
{
    register long   k;

    k = seed / 52774;
    seed = 40692 * (seed - k * 52774) - k * 3791;
    if (seed < 0)
	seed += modulus;
    return (seed);
}

/* Generate random REAL 0..1 */
REAL
ran0()
{
    static REAL     y, maxran, v[97];
    static short    iff;
    int             j;

    if (!iff)
    {
	iff = TRUE;
	maxran = modulus + 1;
	for (j = 0; j < 97; j++)
	    randl();
	for (j = 0; j < 97; j++)
	    v[j] = randl();
	y = randl();
    }
    j = 97.0 * y / maxran;

    y = v[j];
    v[j] = randl();
    return (y / maxran);
}

/*
 * Back Propagation Engine - as described by McClelland & Rumelhart /
 * Sejnowski & Rosenberg
 */

/* logistic 'squashing' function (+/- 1.0) */
#define logistic(x) ((REAL)1.0 / ((REAL)1.0 + (REAL)exp(-(x))))

void
                compute_output()
{
    int             i, j;
    REAL            netinp;

    for (i = NUM_IN; i < TOTAL; i++)
    {
	netinp = bias[i];

	for (j = fwt_to[i]; j < lwt_to[i]; j++)
	    netinp += activation[j] * weight[i][j];

	/* Trigger neuron */
	activation[i] = logistic(netinp);
    }
}

/*
 * load weights - load all link weights from a disk file
 */
void
                load_wts(fname)
     char           *fname;
{
    int             i, j;
    double          t;
    FILE           *ifp;

    if (!(ifp = fopen(fname, "r")))
    {
	printf("Creating new file : %s ...\n", fname);
	return;
    }

    /* Load input units to hidden layer weights */
    for (i = NUM_IN; i < NUM_IN + NUM_HID; i++)
	for (j = fwt_to[i]; j < lwt_to[i]; j++)
	{
	    fscanf(ifp, "%lf", &t);
	    weight[i][j] = t;
	}

    /* Load hidden layer to output units weights */
    for (i = NUM_IN + NUM_HID; i < TOTAL; i++)
	for (j = fwt_to[i]; j < lwt_to[i]; j++)
	{
	    fscanf(ifp, "%lf", &t);
	    weight[i][j] = t;
	}

    /* Load bias weights */
    for (j = NUM_IN; j < TOTAL; j++)
    {
	fscanf(ifp, "%lf", &t);
	bias[j] = t;
    }

    fclose(ifp);
}

/* Initialize network - wire up units and make random weights */
void
                init()
{
    int             i, j;

    for (i = NUM_IN; i < TOTAL; i++)
	if (!(weight[i] = calloc(TOTAL - NUM_OUT, sizeof(REAL))))
	    fail("init: Out of Memory!");

    /* Connect input units to hidden layer */
    for (i = NUM_IN; i < NUM_IN + NUM_HID; i++)
    {
	fwt_to[i] = 0;
	lwt_to[i] = NUM_IN;
	for (j = 0; j < NUM_IN; j++)
	    weight[i][j] = rndwt();
    }

    /* Connect hidden units to output layer */
    for (i = NUM_IN + NUM_HID; i < TOTAL; i++)
    {
	fwt_to[i] = NUM_IN;
	lwt_to[i] = NUM_IN + NUM_HID;
	for (j = NUM_IN; j < NUM_IN + NUM_HID; j++)
	    weight[i][j] = rndwt();
    }

    /* Randomize bias weights */
    for (j = NUM_IN; j < TOTAL; j++)
	bias[j] = rndwt();
}

static char    *
                instr(ct, cs)
     char           *ct, *cs;
{
    int             l = strlen(cs);

    for (; *ct; ct++)
	if (!strncmp(ct, cs, l))
	    return (ct);
    return (NULL);
}

char           *
                dupstr(s)
     char           *s;
{
    char           *new;

    if (!(new = malloc(strlen(s) + 1)))
	fail("dupstr : Out of memory!");
    return (strcpy(new, s));
}

/* Skip some lines */
void
                skiplns(nl, ifp)
     int             nl;
     FILE           *ifp;
{
    while (nl--)
	while (!feof(ifp) && fgetc(ifp) != '\n');
}

/* Convert AA letter to numeric code (0-20) */
int
                aanum(ch)
     int             ch;
{
    static int      aacvs[] =
    {
	999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2,
	20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 20
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 20);
}

const char     *rescodes = "ARNDCQEGHILKMFPSTWYVXXX";

/* Make prediction */
void
                predict()
{
    int             aa, i, j, k, n, winpos, changes;
    char            pred, predsst[5000], bestpred[5000], *che = "CHE";
    float           score_c[5000], score_h[5000], score_e[5000], bestsc,
                    score;

    for (winpos = 0; winpos < seqlen; winpos++)
    {
	for (j = 0; j < NUM_IN; j++)
	    activation[j] = 0.0;
	for (j = WINL; j <= WINR; j++)
	{
	    if (j + winpos >= 0 && j + winpos < seqlen)
		for (aa=0; aa<20; aa++)
		    activation[(j - WINL) * 21 + aa] = logistic(profile[j+winpos][aa]/100.0);
	    else
		activation[(j - WINL) * 21 + 20] = 1.0;
	}
	
	compute_output();
	
	if (activation[TOTAL - NUM_OUT + 3] > 0.5)
	    pred = 'S';
	else if (activation[TOTAL - NUM_OUT + 2] > 0.4)
	{
	    if (activation[TOTAL - NUM_OUT] > activation[TOTAL - NUM_OUT + 1])
		pred = 'I';
	    else
		pred = 'O';
	}
	else
	{
	    if (activation[TOTAL - NUM_OUT] > activation[TOTAL - NUM_OUT + 1])
		pred = '+';
	    else
		pred = '-';
	}
	
	printf("%c %c %5.3f %5.3f %5.3f %5.3f\n", seq[winpos], pred, activation[TOTAL - NUM_OUT], activation[TOTAL - NUM_OUT + 1], activation[TOTAL - NUM_OUT + 2], activation[TOTAL - NUM_OUT + 3]);
    }
}

#define CH malloc_verify(), printf("Heap OK at line : %d.\n",__LINE__);


/* Read PSI AA frequency data */
int             getmtx(FILE *lfil)
{
    int             i, j, naa;
    char            buf[256], *p;

    if (fscanf(lfil, "%d", &naa) != 1)
      fail("Bad mtx file - no sequence length!");

    if (naa > MAXSEQLEN)
      fail("Input sequence too long!");

    if (fscanf(lfil, "%s", seq) != 1)
      fail("Bad mtx file - no sequence!");

    while (!feof(lfil))
      {
	if (!fgets(buf, 65536, lfil))
	  fail("Bad mtx file!");
	if (!strncmp(buf, "-32768 ", 7))
	  {
	    for (j=0; j<naa; j++)
	      {
		if (sscanf(buf, "%*d%d%*d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%*d%d", &profile[j][ALA],  &profile[j][CYS], &profile[j][ASP],  &profile[j][GLU],  &profile[j][PHE],  &profile[j][GLY],  &profile[j][HIS],  &profile[j][ILE],  &profile[j][LYS],  &profile[j][LEU],  &profile[j][MET],  &profile[j][ASN],  &profile[j][PRO],  &profile[j][GLN],  &profile[j][ARG],  &profile[j][SER],  &profile[j][THR],  &profile[j][VAL],  &profile[j][TRP],  &profile[j][TYR]) != 20)
		  fail("Bad mtx format!");
		if (!fgets(buf, 65536, lfil))
		  break;
	      }
	  }
      }

    return naa;
}

main(argc, argv)
     int             argc;
     char          **argv;

{
    int             i, niters;
    char name[160];
    
    FILE           *ifp;

    /* malloc_debug(3); */
    if (argc != 3)
	fail("usage : mem_pred weight-file mtx-file");
    randomize();
    init();
    load_wts(wtfnm = argv[1]);

    ifp = fopen(argv[2], "r");

    if (!ifp)
	fail("Cannot open mtx file!");
    seqlen = getmtx(ifp);

    fclose(ifp);

    predict();

    return 0;
}
