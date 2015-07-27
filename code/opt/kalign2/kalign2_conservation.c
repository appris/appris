/*
	kalign2_conservation.c
	
	Released under GPL - see the 'COPYING' file   
	
	Copyright (C) 2006 Timo Lassmann <timolassmann@gmail.com>
	
	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
    
	Please send bug reports, comments etc. to:
	timolassmann@gmail.com
*/

#include <math.h>

#include "kalign2.h"

void entrophy(int* prof,int len)
{
	int i,j;
	float shannon = 0.0;
	float log_two = log(2);
	fprintf(stderr,"%f\n",log_two);
	for ( i = 0; i < len; i++){
		shannon = 0.0;
		//prof[3] = 10;
		//prof[23] += 10;
		for ( j = 0; j < 23;j++){
			if(prof[j]){
			shannon += (float)prof[j]* log((float)prof[j]/(float)prof[23])/log_two;
		
		//	fprintf(stderr,"%f += %d/%d * %f\n",shannon,prof[j],prof[23],log((float)prof[j]/(float)prof[23])/log_two);
			}
		}
		fprintf(stderr,"%f	",shannon);
		
		if (prof[23] < 23){
			shannon = -shannon / (log((float)prof[23])/log_two);
		}else{
			shannon = -shannon / (log((float)23)/log_two);
		}
		fprintf(stderr,"%f\n",shannon);
		prof+=64;
	}
}
