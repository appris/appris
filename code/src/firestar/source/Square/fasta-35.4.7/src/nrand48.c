
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/*  $Id: nrand48.c 27 2008-06-30 16:27:31Z pearson $ */
/* $Revision: 28 $  */

#include <stdlib.h>
#include <time.h>

void 
irand(int n)	/* initialize random number generator */
{
  if (n == 0) {
    n = time(NULL);
    n = n % 16381;
    if ((n % 2)==0) n++;
  }
  srand48(n);
}

int
nrand(int n)	/* returns a random number between 0 and n-1 where n < 64K) */
{
  int rn;

  rn = lrand48();
  rn = rn >> 16;
  rn = (rn % n);
  return rn;
}

