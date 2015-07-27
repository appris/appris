
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/*  $Id: nrandom.c 27 2008-06-30 16:27:31Z pearson $ */
/* $Revision: 28 $  */

#include <stdlib.h>
#include <time.h>

void 
irand(n)	/* initialize random number generator */
     int n;
{
  if (n == 0) {
    n = time(NULL);
    n = n % 16381;
    if ((n % 2)==0) n++;
  }
  srandom(n);
}

int
nrand(n)	/* returns a random number between 0 and n-1 where n < 2^24) */
     int n;
{
  int rn;

  rn = random();
  rn = (rn % n);
  return rn;
}

