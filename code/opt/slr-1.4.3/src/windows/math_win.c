/*
 *  Copyright 2003-2007 Tim Massingham (tim.massingham@ebi.ac.uk)
 *  Funded by EMBL - European Bioinformatics Institute
 */
/*
 *  This file is part of SLR ("Sitewise Likelihood Ratio")
 *
 *  SLR is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SLR is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with SLR.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <assert.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "math_win.h"

unsigned long __ul_nan[2] = {0xffffffff,0x7fffffff};



/* Calculates log(1+x) accurately when x close to 0.
 */
double log1p ( const double x){
  double res;

  assert(x>-1.);
  
  if ( fabs(x)>sqrt(DBL_EPSILON)){
    res = log(1.+x);
  } else {
    res = x * ( 1. - x/2.);
  }

  return res;   
} 


/*  Calculates exp(x)-1 accurately when x is close to zero
 */
double expm1 ( const double x){
  double res;

  if ( fabs(x)>sqrt(DBL_EPSILON)){
    res = exp(x)-1.;
  } else {
    res = x * ( 1. + x/2.);
  }

  return res;
}


int finite ( const double x){
  if ( -HUGE_VAL==x ||HUGE_VAL == x || NAN==x){
    return 0;
  }
  
  return 1;
}

double tgamma ( double x) {
	static double c[6] = {76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
	double cum,tmp,y; 
	int j;
	
	tmp = x+5.5;
	tmp -= (x+0.5)*log(tmp);
	y=x+1.;
	cum = 1.000000000190015;
	for ( j=0 ; j<6 ; j++){
		cum += c[j]/y;
		y++;
	}
	
	return exp(-tmp+log(2.5066282746310005*cum/x));
}
