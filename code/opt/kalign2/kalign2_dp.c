/*
	kalign2_dp.c 
	
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

#include "kalign2.h"

int* f_only_pp_dyn(int* path, struct dp_matrix *dp,const float* fprof1,const float* fprof2,const int len_a,const int len_b,int fdim,int stride)
{
//	unsigned int freq[26];
	
	struct states* s = 0;
	char** trace = 0;
	char* tracep = 0;
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	register int i = 0;
	register int j = 0;
	register int c = 0;
	register int f = 0;

	s = dp->s;
	
	trace = dp->tb;

	trace[len_a][len_b] = 32;

	fprof1 += len_a * stride;
	

	s[len_b].a = 0.0;
	s[len_b].ga = -FLOATINFTY;
	s[len_b].gb = -FLOATINFTY;
	//init of first row;
	tracep = trace[len_a];
	
	j = len_b;
	while(--j){
		s[j].a = -FLOATINFTY;
		//s[j].ga = 0;	

		s[j].ga = s[j+1].a;//+prof2[29];
		if (s[j+1].ga > s[j].ga){
			s[j].ga = s[j+1].ga ;
		}
		s[j].gb = -FLOATINFTY;
		tracep[j] = 8;
	}
	
	s[0].a = -FLOATINFTY;
	s[0].ga = -FLOATINFTY;
	s[0].gb = -FLOATINFTY;
	i = len_a;
	while(--i){
		
		fprof1 -= stride;
		
		tracep = trace[i];
  		pa = s[len_b].a;
		pga = s[len_b].ga;
		pgb = s[len_b].gb;
		s[len_b].a = -FLOATINFTY;
		s[len_b].ga = -FLOATINFTY;
		//s[len_b].gb = 0;

		s[len_b].gb = pa;//+prof1[29];
		if(pgb > s[len_b].gb){
			s[len_b].gb = pgb;
		}
	
		tracep[len_b] = 16;
		
		j = len_b;
		
		fprof2 += len_b * stride;
		
		while(--j){
			
			fprof2 -= stride;
			ca = s[j].a;

			c = 1;
			if((pga) > pa){
				pa = pga;
				c = 2;
			}
			if((pgb) > pa){
				pa = pgb;
				c = 4;
			}
					
			for (f = 0; f < fdim;f++){
//				fprintf(stderr,"%d	%d:	%d\n",i,j,fprof1[pga] * fprof2[pga+fdim]);
				pa += fprof1[f] * fprof2[f+fdim];
			}

			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = s[j+1].a;
			if (s[j+1].ga > s[j].ga){
				s[j].ga = s[j+1].ga;
				c |= 8;
			}
			
			pgb = s[j].gb;
			
			s[j].gb = ca;
			if(pgb > s[j].gb){
				s[j].gb = pgb;
				c |= 16;
			}
			tracep[j] = c;
			pa = ca;

		}
			
		fprof2 -= stride;
		//LAST CELL (0)
		ca = s[0].a;

		c = 1;
		if((pga) > pa){
			pa = pga;
			c = 2;
		}
		if((pgb) > pa){
			pa = pgb;
			c = 4;
		}
				
		for (f = 0; f < fdim;f++){
//			fprintf(stderr,"%d	%d:	%d\n",i,j,fprof1[pga] * fprof2[pga+fdim]);
			pa += fprof1[f] * fprof2[f+fdim];
		}
		
		s[0].a = pa;
		
		s[0].ga = -FLOATINFTY;
		
		pgb = s[0].gb;
		s[0].gb = ca;
 		if(pgb> s[0].gb){
			s[0].gb = pgb;
			c |= 16;
		}
		tracep[0] = c;
		
	}
	
	fprof1 -= stride;

	tracep = trace[0];
	j = len_b;

	fprof2 += len_b * stride;

	pa = s[j].a;
	pga = s[j].ga;
	pgb = s[j].gb;
	s[j].a = -FLOATINFTY;
	s[j].ga = -FLOATINFTY;
	//s[j].gb = -INFTY;
	s[len_b].gb = pa;//+prof1[29];
	if(pgb > s[len_b].gb){
		s[len_b].gb = pgb;
	}



	while(--j){

		fprof2 -= stride;

		ca = s[j].a;

		c = 1;

		if((pga) > pa){
			pa = pga;
			c = 2;
		}

		if((pgb) > pa){
			pa = pgb;
			c = 4;
		}

		for (f = 0; f < fdim;f++){
			pa += fprof1[f] * fprof2[f+fdim];
//			fprintf(stderr,"%d	%d:	%d\n",i,j,fprof1[pga] * fprof2[pga+fdim]);
		}


		s[j].a = pa;
		pga = s[j].ga;
		s[j].ga = s[j+1].a;
		if (s[j+1].ga > s[j].ga){
			s[j].ga = s[j+1].ga;
			c |= 8;
		}
		pgb = s[j].gb;
		s[j].gb = -FLOATINFTY;

		tracep[j] = c;
		pa = ca;
	}

	fprof2 -= stride;

	ca = s[0].a;

	c = 1;

	if((pga) > pa){
		pa = pga;
		c = 2;
	}
	if((pgb) > pa){
		pa = pgb;
		c = 4;
	}

	for (f = 0; f < fdim;f++){
		pa += fprof1[f] * fprof2[f+fdim];
//		fprintf(stderr,"%d	%d:	%d\n",i,j,fprof1[pga] * fprof2[pga+fdim]);
	}

	s[0].a = pa;

	s[0].ga = s[1].a;
	if (s[1].ga > s[0].ga){
		s[0].ga = s[1].ga;
		c |= 8;
	}

	pgb = s[0].gb;
	s[0].gb = ca;
	if(pgb> s[0].gb){
		s[0].gb = pgb;
		c |= 16;
	}
	tracep[0] = c;

	pgb = s[0].gb;
	c = 2;
	if(s[0].ga > pgb){
		pgb = s[0].ga;
		c = 1;
	}
	if(s[0].a >= pgb){
		pgb = s[0].a;
		c = 0;
	}

	//fprintf(stderr,"SCORE:%d\n",ca);
	ca = c;
	
	i = 0;
	j = 0;
	f = 1;
	while(trace[i][j] < 32){
	//	fprintf(stderr,"%d->%d	%d:%d	%d:%d\n",c,trace[i][j],i,j,len_a,len_b);
		switch(f){
			case 0:
				if (trace[i][j] & 2){
					f = 1;
					if(i+1!= len_a){
						path[c+1] |= 16;
	//					fprintf(stderr,"GAP_CLOSE\n");
					}else{
						path[c+1] |= 32+16;
					}
				}else if (trace[i][j] & 4){
					f = 2;
					if(j+1!= len_b){
						path[c+1] |= 16;
	//					fprintf(stderr,"GAP_CLOSE\n");
					}else{
						path[c+1] |= 32+16;
					}
				}

				//path[c] = 0;
				i++;
				j++;
			break;
			case 1:
				if(trace[i][j] & 8){
					f = 1;
					if(i!=0 && i!= len_a){
	//				/	fprintf(stderr,"GAP_EXT\n");
						if(!(path[c]&16)){
							path[c] |= 8;
						}
					}else{
						if(!(path[c]&16)){
							path[c] |= 32+8;
						}
					}
				}else{
					f = 0;
					if(i!=0 && i!= len_a){
	//					fprintf(stderr,"GAP_OPEN\n");
						path[c] |= 4;
					}else{
						path[c] |= 32+4;
					}
				}
				path[c] |= 1;
				j++;
			break;
			case  2:
				if(trace[i][j] & 16){
					f = 2;
					if(j !=0 && j != len_b){
	//					fprintf(stderr,"GAP_EXT\n");
						if(!(path[c]&16)){
							path[c] |= 8;
						}
					}else{
						if(!(path[c]&16)){
							path[c] |= 32+8;
						}
					}
				}else{
					f = 0;
					if(j!=0 && j != len_b){
	//					fprintf(stderr,"GAP_OPEN\n");
						path[c] |= 4;
					}else{
						path[c] |= 32+4;
					}
					
				}
				path[c] |= 2;
				i++;
			break;
		}
		c++;
	}
	path[0] = c-1;
	path[c] = 3;
	path[c+1] = pgb;
	return path;
}



int* fpp_dyn(int* path, struct dp_matrix *dp,const float* prof1,const float* prof2,const float* fprof1,const float* fprof2,const int len_a,const int len_b,int fdim,int stride)
{
	unsigned int freq[26];
	
	struct states* s = 0;
	char** trace = 0;
	char* tracep = 0;
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	register int i = 0;
	register int j = 0;
	register int c = 0;
	register int f = 0;
	

	s = dp->s;
	
	trace = dp->tb;

	trace[len_a][len_b] = 32;

	prof1 +=  len_a << 6;
	
	fprof1 += len_a * stride;
	

	s[len_b].a = 0;
	s[len_b].ga = -FLOATINFTY;
	s[len_b].gb = -FLOATINFTY;
	//init of first row;
	tracep = trace[len_a];
	
	j = len_b;
	while(--j){
		s[j].a = -FLOATINFTY;
		//s[j].ga = 0;	
		
		s[j].ga = s[j+1].a+prof2[29];//+prof2[29];
		if (s[j+1].ga+prof2[29] > s[j].ga){
			s[j].ga = s[j+1].ga+prof2[29];
		}
		s[j].gb = -FLOATINFTY;
		tracep[j] = 8;
	}
	
	s[0].a = -FLOATINFTY;
	s[0].ga = -FLOATINFTY;
	s[0].gb = -FLOATINFTY;
	i = len_a;
	while(--i){
		prof1 -= 64;
		
		fprof1 -= stride;

		c = 1;
		for (j = 26; j--;){
			if(prof1[j]){
				freq[c] = j;
				c++;	
			}
		}
		freq[0] = c;
		
		tracep = trace[i];
		pa = s[len_b].a;
		pga = s[len_b].ga;
		pgb = s[len_b].gb;
		s[len_b].a = -FLOATINFTY;
		s[len_b].ga = -FLOATINFTY;
		//s[len_b].gb = 0;
		
		s[len_b].gb = pa+prof1[29];//+prof1[29];
		if(pgb+prof1[29] > s[len_b].gb){
			s[len_b].gb = pgb+prof1[29];
		}
	
		tracep[len_b] = 16;
		
		j = len_b;
		prof2 += len_b << 6;
		
		fprof2 += len_b * stride;
		
		while(--j){
			prof2 -= 64;
			
			fprof2 -= stride;
			ca = s[j].a;

			c = 1;
			if((pga += prof2[91]) > pa){
				pa = pga;
				c = 2;
			}
			if((pgb += prof1[91]) > pa){
				pa = pgb;
				c = 4;
			}
			
			prof2 += 32;
			for (f = freq[0];--f;){
				pa += prof1[freq[f]]*prof2[freq[f]];
			}
			prof2 -= 32;
			
			for (f = 0; f < fdim;f++){
			//	fprintf(stderr,"%d	%d:	%d\n",i,j,fprof1[pga] * fprof2[pga+fdim]);
				pa += fprof1[f] * fprof2[f+fdim];
			}

			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = s[j+1].a+prof2[27];
			if (s[j+1].ga+prof2[28] > s[j].ga){
				s[j].ga = s[j+1].ga+prof2[28];
				c |= 8;
			}
			
			pgb = s[j].gb;
			
			s[j].gb = ca+prof1[27];
			if(pgb+prof1[28] > s[j].gb){
				s[j].gb = pgb+prof1[28];
				c |= 16;
			}
			tracep[j] = c;
			pa = ca;

		}
	
		prof2 -= 64;
		
		fprof2 -= stride;
		//LAST CELL (0)
		ca = s[0].a;

		c = 1;
		if((pga+=prof2[91]) > pa){
			pa = pga;
			c = 2;
		}
		if((pgb+=prof1[91]) > pa){
			pa = pgb;
			c = 4;
		}
		
		prof2 += 32;
		for (f = freq[0];--f;){
			pa += prof1[freq[f]]*prof2[freq[f]];
		}
		prof2 -= 32;
		
		for (f = 0; f < fdim;f++){
		//	fprintf(stderr,"%d	%d:	%d\n",i,j,fprof1[pga] * fprof2[pga+fdim]);
			pa += fprof1[f] * fprof2[f+fdim];
		}
		
		s[0].a = pa;
		
		s[0].ga = -FLOATINFTY;
		
		pgb = s[0].gb;
		s[0].gb = ca+prof1[27]+prof1[29];
 		if(pgb+prof1[29] > s[0].gb){
			s[0].gb = pgb+prof1[29];
			c |= 16;
		}
		tracep[0] = c;	
		
	}
	prof1 -= 64;
	
	fprof1 -= stride;
	
	c = 1;
	for (j = 26; j--;){
		if(prof1[j]){
			freq[c] = j;
			c++;	
		}
	}
	freq[0] = c;
	
	tracep = trace[0];
	j = len_b;
	prof2 += len_b << 6;
	
	fprof2 += len_b * stride;
	
	pa = s[j].a;
	pga = s[j].ga;
	pgb = s[j].gb;
	s[j].a = -FLOATINFTY;
	s[j].ga = -FLOATINFTY;
	//s[j].gb = -INFTY;
	s[len_b].gb = pa+prof1[29];//+prof1[29];
	if(pgb+prof1[29] > s[len_b].gb){
		s[len_b].gb = pgb+prof1[29];
	}


	
	while(--j){
		prof2 -= 64;
		
		fprof2 -= stride;
		
		ca = s[j].a;

		c = 1;

		if((pga+=prof2[91]) > pa){
			pa = pga;
			c = 2;
		}

		if((pgb+=prof1[91]) > pa){
			pa = pgb;
			c = 4;
		}
		
		prof2+=32;
		for (f = freq[0];--f;){
			pa += prof1[freq[f]]*prof2[freq[f]];
		}
		prof2-=32;
		
		
		for (f = 0; f < fdim;f++){
			pa += fprof1[f] * fprof2[f+fdim];
//			fprintf(stderr,"%d	%d:	%d\n",i,j,fprof1[pga] * fprof2[pga+fdim]);
		}
		
		
		s[j].a = pa;
		pga = s[j].ga;
		s[j].ga = s[j+1].a+prof2[27]+prof2[29];
		if (s[j+1].ga+prof2[29] > s[j].ga){
			s[j].ga = s[j+1].ga+prof2[29];
			c |= 8;
		}	
		pgb = s[j].gb;
		s[j].gb = -FLOATINFTY;	
		
		tracep[j] = c;
		pa = ca;
	}
	prof2 -= 64;
	
	fprof2 -= stride;

	ca = s[0].a;
	
	c = 1;
	
	if((pga+=prof2[91]) > pa){
		pa = pga;
		c = 2;
	}
	if((pgb+=prof1[91]) > pa){
		pa = pgb;
		c = 4;
	}
	prof2+=32;
	for (f = freq[0];--f;){
		pa += prof1[freq[f]]*prof2[freq[f]];
	}
	prof2-=32;
			
	for (f = 0; f < fdim;f++){
		pa += fprof1[f] * fprof2[f+fdim];
//		fprintf(stderr,"%d	%d:	%d\n",i,j,fprof1[pga] * fprof2[pga+fdim]);
	}
	
	s[0].a = pa;
	
	s[0].ga = s[1].a+prof2[27]+prof2[29];
	if (s[1].ga+prof2[29] > s[0].ga){
		s[0].ga = s[1].ga+prof2[29];
		c |= 8;
	}
	
	pgb = s[0].gb;
	s[0].gb = ca+prof1[27]+prof1[29];
	if(pgb +prof1[29]> s[0].gb){
		s[0].gb = pgb+prof1[29];
		c |= 16;
	}	
	tracep[0] = c;

	pgb = s[0].gb;
	c = 2;
	if(s[0].ga > pgb){
		pgb = s[0].ga;
		c = 1;
	}
	if(s[0].a >= pgb){
		pgb = s[0].a;
		c = 0;
	}
	
	//fprintf(stderr,"SCORE:%d\n",ca);
	f = c;
	
	i = 0;
	j = 0;
	c = 1;
	while(trace[i][j] < 32){
	//	fprintf(stderr,"%d->%d	%d:%d	%d:%d\n",c,trace[i][j],i,j,len_a,len_b);
		switch(f){
			case 0:
				if (trace[i][j] & 2){
					f = 1;
					if(i+1!= len_a){
						path[c+1] |= 16;
	//					fprintf(stderr,"GAP_CLOSE\n");
					}else{
						path[c+1] |= 32+16;
					}
				}else if (trace[i][j] & 4){
					f = 2;
					if(j+1!= len_b){
						path[c+1] |= 16;
	//					fprintf(stderr,"GAP_CLOSE\n");
					}else{
						path[c+1] |= 32+16;
					}
				}

				//path[c] = 0;
				i++;
				j++;
			break;
			case 1:
				if(trace[i][j] & 8){
					f = 1;
					if(i!=0 && i!= len_a){
	//				/	fprintf(stderr,"GAP_EXT\n");
						if(!(path[c]&16)){
							path[c] |= 8;
						}
					}else{
						if(!(path[c]&16)){
							path[c] |= 32+8;
						}
					}
				}else{
					f = 0;
					if(i!=0 && i!= len_a){
	//					fprintf(stderr,"GAP_OPEN\n");
						path[c] |= 4;
					}else{
						path[c] |= 32+4;
					}
				}
				path[c] |= 1;
				j++;
			break;
			case  2:
				if(trace[i][j] & 16){
					f = 2;
					if(j !=0 && j != len_b){
	//					fprintf(stderr,"GAP_EXT\n");
						if(!(path[c]&16)){
							path[c] |= 8;
						}
					}else{
						if(!(path[c]&16)){
							path[c] |= 32+8;
						}
					}
				}else{
					f = 0;
					if(j!=0 && j != len_b){
	//					fprintf(stderr,"GAP_OPEN\n");
						path[c] |= 4;
					}else{
						path[c] |= 32+4;
					}
					
				}
				path[c] |= 2;
				i++;
			break;
		}
		c++;
	}
	path[0] = c-1;
	path[c] = 3;
	path[c+1] = pgb;
	return path;
}
/*
int* dna_pp_dyn(int* path, struct dp_matrix *dp,const int* prof1,const int* prof2,const int len_a,const int len_b)
{

	struct states* s = 0;
	char** trace = 0;
	char* tracep = 0;
	register int pa = 0;
	register int pga = 0;
	register int pgb = 0;
	register int ca = 0;
	register int i = 0;
	register int j = 0;
	register int c = 0;

	s = dp->s;
	
	trace = dp->tb;

	trace[len_a][len_b] = 32;

	prof1 +=  len_a * 22;

	s[len_b].a = 0;
	s[len_b].ga = -INFTY;
	s[len_b].gb = -INFTY;
	//init of first row;
	tracep = trace[len_a];
	
	j = len_b;
	while(--j){
		s[j].a = -INFTY;
		//s[j].ga = 0;	
		
		s[j].ga = s[j+1].a+prof2[10];//+prof2[29];
		if (s[j+1].ga+prof2[10] > s[j].ga){
			s[j].ga = s[j+1].ga+prof2[10];
		}
		s[j].gb = -INFTY;
		tracep[j] = 8;
	}
	
	s[0].a = -INFTY;
	s[0].ga = -INFTY;
	s[0].gb = -INFTY;
	i = len_a;
	while(--i){
		prof1 -= 22;

		tracep = trace[i];
		pa = s[len_b].a;
		pga = s[len_b].ga;
		pgb = s[len_b].gb;
		s[len_b].a = -INFTY;
		s[len_b].ga = -INFTY;
		//s[len_b].gb = 0;
		
		s[len_b].gb = pa+prof1[10];//+prof1[29];
		if(pgb+prof1[10] > s[len_b].gb){
			s[len_b].gb = pgb+prof1[10];
		}
	
		tracep[len_b] = 16;
		
		j = len_b;
		prof2 += len_b *22;
		while(--j){
			prof2 -= 22;
			ca = s[j].a;

			c = 1;
			if((pga += prof2[30]) > pa){
				pa = pga;
				c = 2;
			}
			if((pgb += prof1[30]) > pa){
				pa = pgb;
				c = 4;
			}
			
			prof2 += 11;
			for (pga = 8;pga--;){
				pa += prof1[pga]*prof2[pga];
			}
			prof2 -= 11;

			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = s[j+1].a+prof2[8];
			if (s[j+1].ga+prof2[9] > s[j].ga){
				s[j].ga = s[j+1].ga+prof2[9];
				c |= 8;
			}
			
			pgb = s[j].gb;
			
			s[j].gb = ca+prof1[8];
			if(pgb+prof1[9] > s[j].gb){
				s[j].gb = pgb+prof1[9];
				c |= 16;
			}
			tracep[j] = c;
			pa = ca;

		}
	
		prof2 -= 22;
		//LAST CELL (0)
		ca = s[0].a;

		c = 1;
		if((pga+=prof2[30]) > pa){
			pa = pga;
			c = 2;
		}
		if((pgb+=prof1[30]) > pa){
			pa = pgb;
			c = 4;
		}
		
		prof2 += 11;
		for (pga = 8;pga--;){
			pa += prof1[pga]*prof2[pga];
		}
		prof2 -= 11;
		
		s[0].a = pa;
		
		s[0].ga = -INFTY;
		
		pgb = s[0].gb;
		s[0].gb = ca+prof1[8]+prof1[10];
 		if(pgb+prof1[10] > s[0].gb){
			s[0].gb = pgb+prof1[10];
			c |= 16;
		}
		tracep[0] = c;	
		
	}
	prof1 -= 22;
	
	tracep = trace[0];
	j = len_b;
	prof2 += len_b *22;
	pa = s[j].a;
	pga = s[j].ga;
	pgb = s[j].gb;
	s[j].a = -INFTY;
	s[j].ga = -INFTY;
	//s[j].gb = -INFTY;
	s[len_b].gb = pa+prof1[10];//+prof1[29];
	if(pgb+prof1[10] > s[len_b].gb){
		s[len_b].gb = pgb+prof1[10];
	}


	
	while(--j){
		prof2 -= 22;
		ca = s[j].a;

		c = 1;

		if((pga+=prof2[30]) > pa){
			pa = pga;
			c = 2;
		}

		if((pgb+=prof1[30]) > pa){
			pa = pgb;
			c = 4;
		}
		
		prof2+=11;
		
		for (pga = 8;pga--;){
			pa += prof1[pga]*prof2[pga];
		}
		prof2-=11;
		
		s[j].a = pa;
		pga = s[j].ga;
		s[j].ga = s[j+1].a+prof2[2]+prof2[10];
		if (s[j+1].ga+prof2[10] > s[j].ga){
			s[j].ga = s[j+1].ga+prof2[10];
			c |= 8;
		}	
		pgb = s[j].gb;
		s[j].gb = -INFTY;	
		
		tracep[j] = c;
		pa = ca;
	}
	prof2 -= 22;

	ca = s[0].a;
	
	c = 1;
	
	if((pga+=prof2[30]) > pa){
		pa = pga;
		c = 2;
	}
	if((pgb+=prof1[30]) > pa){
		pa = pgb;
		c = 4;
	}
	prof2+=11;
	for (pga = 8;pga--;){
		pa += prof1[pga]*prof2[pga];
	}
	prof2-=11;
	
	s[0].a = pa;
	
	s[0].ga = s[1].a+prof2[8]+prof2[10];
	if (s[1].ga+prof2[10] > s[0].ga){
		s[0].ga = s[1].ga+prof2[10];
		c |= 8;
	}
	
	pgb = s[0].gb;
	s[0].gb = ca+prof1[8]+prof1[10];
	if(pgb +prof1[10]> s[0].gb){
		s[0].gb = pgb+prof1[10];
		c |= 16;
	}	
	tracep[0] = c;

	pgb = s[0].gb;
	c = 2;
	if(s[0].ga > pgb){
		pgb = s[0].ga;
		c = 1;
	}
	if(s[0].a >= pgb){
		pgb = s[0].a;
		c = 0;
	}
	
	//fprintf(stderr,"SCORE:%d\n",ca);
	ca = c;
	
	i = 0;
	j = 0;
	c = 1;
	while(trace[i][j] < 32){
	//	fprintf(stderr,"%d->%d	%d:%d	%d:%d\n",c,trace[i][j],i,j,len_a,len_b);
		switch(ca){
			case 0:
				if (trace[i][j] & 2){
					ca = 1;
					if(i+1!= len_a){
						path[c+1] |= 16;
	//					fprintf(stderr,"GAP_CLOSE\n");
					}else{
						path[c+1] |= 32+16;
					}
				}else if (trace[i][j] & 4){
					ca = 2;
					if(j+1!= len_b){
						path[c+1] |= 16;
	//					fprintf(stderr,"GAP_CLOSE\n");
					}else{
						path[c+1] |= 32+16;
					}
				}

				//path[c] = 0;
				i++;
				j++;
			break;
			case 1:
				if(trace[i][j] & 8){
					ca = 1;
					if(i!=0 && i!= len_a){
	//				/	fprintf(stderr,"GAP_EXT\n");
						if(!(path[c]&16)){
							path[c] |= 8;
						}
					}else{
						if(!(path[c]&16)){
							path[c] |= 32+8;
						}
					}
				}else{
					ca = 0;
					if(i!=0 && i!= len_a){
	//					fprintf(stderr,"GAP_OPEN\n");
						path[c] |= 4;
					}else{
						path[c] |= 32+4;
					}
				}
				path[c] |= 1;
				j++;
			break;
			case  2:
				if(trace[i][j] & 16){
					ca = 2;
					if(j !=0 && j != len_b){
	//					fprintf(stderr,"GAP_EXT\n");
						if(!(path[c]&16)){
							path[c] |= 8;
						}
					}else{
						if(!(path[c]&16)){
							path[c] |= 32+8;
						}
					}
				}else{
					ca = 0;
					if(j!=0 && j != len_b){
	//					fprintf(stderr,"GAP_OPEN\n");
						path[c] |= 4;
					}else{
						path[c] |= 32+4;
					}
					
				}
				path[c] |= 2;
				i++;
			break;
		}
		c++;
	}
	path[0] = c-1;
	path[c] = 3;
	path[c+1] = pgb;
	return path;
}


int* pp_dyn2(int* path, struct dp_matrix *dp,const int* prof1,const int* prof2,const int len_a,const int len_b)
{
	unsigned int freq[26];
	
	struct states* s = 0;
	char** trace = 0;
	char* tracep = 0;
	register int pa = 0;
	register int pga = 0;
	register int pgb = 0;
	register int ca = 0;
	register int i = 0;
	register int j = 0;
	register int c = 0;

	s = dp->s;
	
	trace = dp->tb;

	trace[len_a][len_b] = 32;

	prof1 +=  len_a << 6;

	s[len_b].a = 0;
	s[len_b].ga = -INFTY;
	s[len_b].gb = -INFTY;
	//init of first row;
	tracep = trace[len_a];
	
	j = len_b;
	while(--j){
		s[j].a = -INFTY;
	
		s[j].ga = s[j+1].a+prof2[29];//+prof2[29];
		if (s[j+1].ga+prof2[29] > s[j].ga){
			s[j].ga = s[j+1].ga+prof2[29];
		}
		s[j].gb = -INFTY;
		tracep[j] = 8;
	}
	
	s[0].a = -INFTY;
	s[0].ga = -INFTY;
	s[0].gb = -INFTY;
	i = len_a;
	while(--i){
		prof1 -= 64;

		c = 1;
		for (j = 23; j--;){
			if(prof1[j]){
				freq[c] = j;
				c++;	
			}
		}
		freq[0] = c;
		
		tracep = trace[i];
		pa = s[len_b].a;
		pga = s[len_b].ga;
		pgb = s[len_b].gb;
		s[len_b].a = -INFTY;
		s[len_b].ga = -INFTY;
		
		s[len_b].gb = pa+prof1[29];
		if(pgb+prof1[29] > s[len_b].gb){
			s[len_b].gb = pgb+prof1[29];
		}
	
		tracep[len_b] = 16;
		
		j = len_b;
		prof2 += len_b << 6;
		while(--j){
			prof2 -= 64;
			ca = s[j].a;

			c = 1;
			if((pga += prof2[91]) > pa){
				pa = pga;
				c = 2;
			}
			if((pgb += prof1[91]) > pa){
				pa = pgb;
				c = 4;
			}
			
			prof2 += 32;
			for (pga = freq[0];--pga;){
				pgb = freq[pga];
				pa += prof1[pgb]*prof2[pgb];
			}
			prof2 -= 32;

			s[j].a = pa;

			pga = s[j].ga;
			
			s[j].ga = s[j+1].a+prof2[27];
			if (s[j+1].ga+prof2[28] > s[j].ga){
				s[j].ga = s[j+1].ga+prof2[28];
				c |= 8;
			}
			
			pgb = s[j].gb;
			
			s[j].gb = ca+prof1[27];
			if(pgb+prof1[28] > s[j].gb){
				s[j].gb = pgb+prof1[28];
				c |= 16;
			}
			tracep[j] = c;
			pa = ca;

		}
	
		prof2 -= 64;
		//LAST CELL (0)
		ca = s[0].a;

		c = 1;
		if((pga+=prof2[91]) > pa){
			pa = pga;
			c = 2;
		}
		if((pgb+=prof1[91]) > pa){
			pa = pgb;
			c = 4;
		}
		
		prof2 += 32;
		for (pga = freq[0];--pga;){
			pgb = freq[pga];
			pa += prof1[pgb]*prof2[pgb];
		}
		prof2 -= 32;
		
		s[0].a = pa;
		
		s[0].ga = -INFTY;
		
		pgb = s[0].gb;
		s[0].gb = ca+prof1[27]+prof1[29];
 		if(pgb+prof1[29] > s[0].gb){
			s[0].gb = pgb+prof1[29];
			c |= 16;
		}
		tracep[0] = c;	
		
	}
	prof1 -= 64;
	
	c = 1;
	for (j = 23; j--;){
		if(prof1[j]){
			freq[c] = j;
			c++;	
		}
	}
	freq[0] = c;
	
	tracep = trace[0];
	j = len_b;
	prof2 += len_b << 6;
	pa = s[j].a;
	pga = s[j].ga;
	pgb = s[j].gb;
	s[j].a = -INFTY;
	s[j].ga = -INFTY;
	
	s[len_b].gb = pa+prof1[29];
	if(pgb+prof1[29] > s[len_b].gb){
		s[len_b].gb = pgb+prof1[29];
	}


	
	while(--j){
		prof2 -= 64;
		ca = s[j].a;

		c = 1;

		if((pga+=prof2[91]) > pa){
			pa = pga;
			c = 2;
		}

		if((pgb+=prof1[91]) > pa){
			pa = pgb;
			c = 4;
		}
		
		prof2+=32;
		
		for (pga = freq[0];--pga;){
			pgb = freq[pga];
			pa += prof1[pgb]*prof2[pgb];
		}
		prof2-=32;
		
		s[j].a = pa;
		pga = s[j].ga;
		s[j].ga = s[j+1].a+prof2[27]+prof2[29];
		if (s[j+1].ga+prof2[29] > s[j].ga){
			s[j].ga = s[j+1].ga+prof2[29];
			c |= 8;
		}	
		pgb = s[j].gb;
		s[j].gb = -INFTY;	
		
		tracep[j] = c;
		pa = ca;
	}
	prof2 -= 64;

	ca = s[0].a;
	
	c = 1;
	
	if((pga+=prof2[91]) > pa){
		pa = pga;
		c = 2;
	}
	if((pgb+=prof1[91]) > pa){
		pa = pgb;
		c = 4;
	}
	prof2+=32;
	for (pga = freq[0];--pga;){	
		pgb = freq[pga];
		pa += prof1[pgb]*prof2[pgb];
	}
	prof2-=32;
	
	s[0].a = pa;
	
	s[0].ga = s[1].a+prof2[27]+prof2[29];
	if (s[1].ga+prof2[29] > s[0].ga){
		s[0].ga = s[1].ga+prof2[29];
		c |= 8;
	}
	
	pgb = s[0].gb;
	s[0].gb = ca+prof1[27]+prof1[29];
	if(pgb +prof1[29]> s[0].gb){
		s[0].gb = pgb+prof1[29];
		c |= 16;
	}	
	tracep[0] = c;

	pgb = s[0].gb;
	c = 2;
	if(s[0].ga > pgb){
		pgb = s[0].ga;
		c = 1;
	}
	if(s[0].a >= pgb){
		pgb = s[0].a;
		c = 0;
	}

	ca = c;
	
	int ga = 1; 
	int gb = 1;
	i = 0;
	j = 0;
	c = 1;
	while(trace[i][j] < 32){
		if(i ==0 || j == 0){
			path[c] |= 128;
		}
		if(i ==len_a || j == len_b){
			path[c] |= 64;
		}
	
		switch(ca){
			case 0:
				if (trace[i][j] & 2){
					ca = 1;
				}else if (trace[i][j] & 4){
					ca = 2;
				}
				path[c] = 0;
				i++;
				j++;
			break;
			case 1:
				if(trace[i][j] & 8){
					ca = 1;
				}else{
					path[c-(gb-1)] |= gb << 16;
					gb = 0;
					ca = 0;
				}
				path[c] |= 1;
				j++;
				gb++;
			break;
			case  2:
				if(trace[i][j] & 16){
					ca = 2;
				}else{
					path[c-(ga-1)] |= ga << 16;
					ga = 0;
					ca = 0;
				}
				path[c] |= 2;
				i++;
				ga++;
			break;
		}
		c++;
	}
	if (ca == 1){
		path[c-(gb-1)] |= (gb-1) << 16;	
	}
	if(ca == 2){
		path[c-(ga-1)] |= (ga-1) << 16;	
	}
	path[0] = c-1;
	path[c] = 3;
	path[c+1] = pgb;
	return path;
}

int* ps_dyn2(int* path, struct dp_matrix *dp,const int* prof1,const int* seq2,const int len_a,const int len_b,int sip)
{
	
	struct states* s = 0;
	char** trace = 0;
	char* tracep = 0;
	register int pa = 0;
	register int pga = 0;
	register int pgb = 0;
	register int ca = 0;
	register int i = 0;
	register int j = 0;
	register int c = 0;
	
	const int open = gpo * sip;
	const int ext = gpe *sip; 

	

	s = dp->s;
	
	trace = dp->tb;

	trace[len_a][len_b] = 32;

	prof1 +=  len_a << 6;

	s[len_b].a = 0;
	s[len_b].ga = -INFTY;
	s[len_b].gb = -INFTY;
	tracep = trace[len_a];
	j = len_b;
	

	while(--j){
		s[j].a = -INFTY;
		s[j].ga = s[j+1].a-tgpe;
		if (s[j+1].ga-tgpe > s[j].ga){
			s[j].ga = s[j+1].ga-tgpe;
		}
		
		
		
		s[j].gb = -INFTY;
		tracep[j] = 8;
	}
	
	s[0].a = -INFTY;
	s[0].ga = -INFTY;
	s[0].gb = -INFTY;
	i = len_a;
	while(--i){
		prof1 -= 64;
		
		tracep = trace[i];
		pa = s[len_b].a;
		pga = s[len_b].ga;
		pgb = s[len_b].gb;
		s[len_b].a = -INFTY;
		s[len_b].ga = -INFTY;

		s[len_b].gb = pa+prof1[29];
		if(pgb+prof1[29] > s[len_b].gb){
			s[len_b].gb = pgb+prof1[29];
		}
		
		
		
		
		tracep[len_b] = 16;
		
		j = len_b;
		
		while(--j){

			ca = s[j].a;

			c = 1;
			if((pga -= open) > pa){
				pa = pga;
				c = 2;
			}
			if((pgb += prof1[91]) > pa){
				pa = pgb;
				c = 4;
			}
			
			pa += prof1[32 + seq2[j]];

			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = s[j+1].a-open;
			if (s[j+1].ga-ext > s[j].ga){
				s[j].ga = s[j+1].ga-ext;
				c |= 8;
			}
			
			pgb = s[j].gb;
			
			s[j].gb = ca+prof1[27];
			if(pgb+prof1[28] > s[j].gb){
				s[j].gb = pgb+prof1[28];
				c |= 16;
			}
			tracep[j] = c;
			pa = ca;

		}
	
		ca = s[0].a;

		c = 1;
		if((pga-=open) > pa){
			pa = pga;
			c = 2;
		}
		if((pgb+=prof1[91]) > pa){
			pa = pgb;
			c = 4;
		}
		pa += prof1[32+seq2[0]];
		s[0].a = pa;
		
		s[0].ga = -INFTY;
		
		pgb = s[0].gb;
		s[0].gb = ca+prof1[27]+prof1[29];
 		if(pgb+prof1[29] > s[0].gb){
			s[0].gb = pgb+prof1[29];
			c |= 16;
		}
		tracep[0] = c;	
		
	}
	prof1 -= 64;
	

	
	tracep = trace[0];
	j = len_b;
	pa = s[j].a;
	pga = s[j].ga;
	pgb = s[j].gb;
	s[j].a = -INFTY;
	s[j].ga = -INFTY;

	s[len_b].gb = pa+prof1[29];
	if(pgb+prof1[29] > s[len_b].gb){
		s[len_b].gb = pgb+prof1[29];
	}

	
	while(--j){

		ca = s[j].a;

		c = 1;

		if((pga-=open) > pa){
			pa = pga;
			c = 2;
		}

		if((pgb+=prof1[91]) > pa){
			pa = pgb;
			c = 4;
		}
		pa += prof1[32+seq2[j]];		
		s[j].a = pa;
		pga = s[j].ga;
		s[j].ga = s[j+1].a-(open+tgpe);
		if (s[j+1].ga-tgpe > s[j].ga){
			s[j].ga = s[j+1].ga-tgpe;
			c |= 8;
		}	
		pgb = s[j].gb;
		s[j].gb = -INFTY;	
		
		tracep[j] = c;
		pa = ca;
	}


	ca = s[0].a;
	
	c = 1;
	
	if((pga-=open) > pa){
		pa = pga;
		c = 2;
	}
	if((pgb+=prof1[91]) > pa){
		pa = pgb;
		c = 4;
	}
	pa += prof1[32+seq2[0]];	
	s[0].a = pa;
	
	s[0].ga = s[1].a-(open+tgpe);
	if (s[1].ga-tgpe > s[0].ga){
		s[0].ga = s[1].ga-tgpe;
		c |= 8;
	}
	
	pgb = s[0].gb;
	s[0].gb = ca+prof1[27]+prof1[29];
	if(pgb+prof1[29] > s[0].gb){
		s[0].gb = pgb+prof1[29];
		c |= 16;
	}	
	tracep[0] = c;


	pgb = s[0].gb;
	c = 2;
	if(s[0].ga > pgb){
		pgb = s[0].ga;
		c = 1;
	}
	if(s[0].a >= pgb){
		pgb = s[0].a;
		c = 0;
	}
	
	ca = c;
	int ga = 1; 
	int gb = 1;
	i = 0;
	j = 0;
	c = 1;
	while(trace[i][j] < 32){
		if(i ==0 || j == 0){
			path[c] |= 128;
		}
		if(i ==len_a || j == len_b){
			path[c] |= 64;
		}
	
		switch(ca){
			case 0:
				if (trace[i][j] & 2){
					ca = 1;
				}else if (trace[i][j] & 4){
					ca = 2;
				}
				path[c] = 0;
				i++;
				j++;
			break;
			case 1:
				if(trace[i][j] & 8){
					ca = 1;
				}else{
					path[c-(gb-1)] |= gb << 16;
					gb = 0;
					ca = 0;
				}
				path[c] |= 1;
				j++;
				gb++;
			break;
			case  2:
				if(trace[i][j] & 16){
					ca = 2;
				}else{
					path[c-(ga-1)] |= ga << 16;
					ga = 0;
					ca = 0;
				}
				path[c] |= 2;
				i++;
				ga++;
			break;
		}
		c++;
	}
	if (ca == 1){
		path[c-(gb-1)] |= (gb-1) << 16;	
	}
	if(ca == 2){
		path[c-(ga-1)] |= (ga-1) << 16;	
	}
	path[0] = c-1;
	path[c] = 3;
	path[c+1] = pgb;
	return path;
}

int* ss_dyn2(int**subm,int* path, struct dp_matrix *dp,const int* seq1,const int* seq2,const int len_a,const int len_b)
{
	struct states* s = 0;
	int *subp = 0;
	char** trace = 0;
	char* tracep = 0;
	register int pa = 0;
	register int pga = 0;
	register int pgb = 0;
	register int ca = 0;
	register int i = 0;
	register int j = 0;
	register int c = 0;
	
	s = dp->s;
	
	trace = dp->tb;

	trace[len_a][len_b] = 32;

	s[len_b].a = 0;
	s[len_b].ga = -INFTY;
	s[len_b].gb = -INFTY;
	

	tracep = trace[len_a];
	j = len_b;
	

	while(--j){
		s[j].a = -INFTY;

		s[j].ga = s[j+1].a-tgpe;
		if (s[j+1].ga-tgpe > s[j].ga){
			s[j].ga = s[j+1].ga-tgpe;
		}
				
		s[j].gb = -INFTY;
		tracep[j] = 8;
	}
	
	s[0].a = -INFTY;
	s[0].ga = -INFTY;
	s[0].gb = -INFTY;
	
	i = len_a;
	while(--i){
		
		tracep = trace[i];
		pa = s[len_b].a;
		pga = s[len_b].ga;
		pgb = s[len_b].gb;
		
		s[len_b].a = -INFTY;
		s[len_b].ga = -INFTY;
		
		s[len_b].gb = pa-tgpe;
		if(pgb-tgpe > s[len_b].gb){
			s[len_b].gb = pgb-tgpe;
		}
		
		
		tracep[len_b] = 16;
		j = len_b;
		subp = subm[seq1[i]];
		while(--j){
			ca = s[j].a;
			
			c = 1;
			if((pga -= gpo) > pa){
				pa = pga;
				c = 2;
			}
			if((pgb -= gpo) > pa){
				pa = pgb;
				c = 4;
			}

			pa += subp[seq2[j]];

			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = s[j+1].a-gpo;
			if (s[j+1].ga-gpe > s[j].ga){
				s[j].ga = s[j+1].ga-gpe;
				c |= 8;
			}
			
			pgb = s[j].gb;
			
			s[j].gb = ca-gpo;
			if(pgb-gpe > s[j].gb){
				s[j].gb = pgb-gpe;
				c |= 16;
			}
			tracep[j] = c;
			pa = ca;

		}
	
		ca = s[0].a;

		c = 1;
		if((pga-=gpo) > pa){
			pa = pga;
			c = 2;
		}
		if((pgb-=gpo) > pa){
			pa = pgb;
			c = 4;
		}
		
		pa += subp[seq2[0]];
		
		s[0].a = pa;
		
		s[0].ga = -INFTY;
		
		pgb = s[0].gb;
		s[0].gb = ca-(gpo+tgpe);
 		if(pgb-tgpe > s[0].gb){
			s[0].gb = pgb-tgpe;
			c |= 16;
		}
		tracep[0] = c;			
	}

	subp = subm[seq1[0]];
	tracep = trace[0];
	j = len_b;
	pa = s[j].a;
	pga = s[j].ga;
	pgb = s[j].gb;
	s[j].a = -INFTY;
	s[j].ga = -INFTY;
	
	s[j].gb = pa-tgpe;
	if(pgb-tgpe > s[j].gb){
		s[j].gb = pgb-tgpe;
	}

	while(--j){

		ca = s[j].a;

		c = 1;

		if((pga-=gpo) > pa){
			pa = pga;
			c = 2;
		}

		if((pgb-=gpo) > pa){
			pa = pgb;
			c = 4;
		}
		
		pa += subp[seq2[j]];
		
		s[j].a = pa;
		
		pga = s[j].ga;
		s[j].ga = s[j+1].a-(gpo+tgpe);
		if (s[j+1].ga-tgpe > s[j].ga){
			s[j].ga = s[j+1].ga-tgpe;
			c |= 8;
		}	
		pgb = s[j].gb;
		s[j].gb = -INFTY;	
		tracep[j] = c;
		pa = ca;
	}
	
	ca = s[0].a;
	
	c = 1;
	
	if((pga-=gpo) > pa){
		pa = pga;
		c = 2;
	}
	if((pgb-=gpo) > pa){
		pa = pgb;
		c = 4;
	}

	pa += subp[seq2[0]];
	
	s[0].a = pa;
	
	
	s[0].ga = s[1].a-(gpo+tgpe);
	if (s[1].ga-tgpe > s[0].ga){
		s[0].ga = s[1].ga-tgpe;
		c |= 8;
	}
	
	pgb = s[0].gb;
	s[0].gb = ca-(gpo+tgpe);
	if(pgb-tgpe > s[0].gb){
		s[0].gb = pgb-tgpe;
		c |= 16;
	}	
	tracep[0] = c;


	pgb = s[0].gb;
	c = 2;
	if(s[0].ga > pgb){
		pgb = s[0].ga;
		c = 1;
	}
	if(s[0].a >= pgb){
		pgb = s[0].a;
		c = 0;
	}
	
	ca = c;
	
	int ga = 1; 
	int gb = 1;
	i = 0;
	j = 0;
	c = 1;
	while(trace[i][j] < 32){
		if(i ==0 || j == 0){
			path[c] |= 128;
		}
		if(i ==len_a || j == len_b){
			path[c] |= 64;
		}
	
		switch(ca){
			case 0:
				if (trace[i][j] & 2){
					ca = 1;
				}else if (trace[i][j] & 4){
					ca = 2;
				}
				path[c] = 0;
				i++;
				j++;
			break;
			case 1:
				if(trace[i][j] & 8){
					ca = 1;
				}else{
					path[c-(gb-1)] |= gb << 16;
					gb = 0;
					ca = 0;
				}
				path[c] |= 1;
				j++;
				gb++;
			break;
			case  2:
				if(trace[i][j] & 16){
					ca = 2;
				}else{
					path[c-(ga-1)] |= ga << 16;
					ga = 0;
					ca = 0;
				}
				path[c] |= 2;
				i++;
				ga++;
			break;
		}
		c++;
	}
	if (ca == 1){
		path[c-(gb-1)] |= (gb-1) << 16;	
	}
	if(ca == 2){
		path[c-(ga-1)] |= (ga-1) << 16;	
	}
	path[0] = c-1;
	path[c] = 3;
	path[c+1] = pgb;
	return path;
}



int* aapp_dyn(int* path, struct dp_matrix *dp,const int* prof1,const int* prof2,const int len_a,const int len_b,const int mmbonus)
{
	unsigned int freq[26];
	
	struct states* s = 0;
	char** trace = 0;
	char* tracep = 0;
	register int pa = 0;
	register int pga = 0;
	register int pgb = 0;
	register int ca = 0;
	register int i = 0;
	register int j = 0;
	register int c = 0;

	s = dp->s;
	
	trace = dp->tb;

	trace[len_a][len_b] = 32;

	prof1 +=  len_a << 6;

	s[len_b].a = 0;
	s[len_b].ga = -INFTY;
	s[len_b].gb = -INFTY;
	//init of first row;
	tracep = trace[len_a];
	
	j = len_b;
	while(--j){
		s[j].a = -INFTY;
		
		s[j].ga = s[j+1].a+prof2[29];
		if (s[j+1].ga+prof2[29] > s[j].ga){
			s[j].ga = s[j+1].ga+prof2[29];
		}
		s[j].gb = -INFTY;
		tracep[j] = 8;
	}
	
	s[0].a = -INFTY;
	s[0].ga = -INFTY;
	s[0].gb = -INFTY;
	i = len_a;
	while(--i){
		prof1 -= 64;

		c = 1;
		for (j = 26; j--;){
			if(prof1[j]){
				freq[c] = j;
				c++;	
			}
		}
		freq[0] = c;
		
		tracep = trace[i];
		pa = s[len_b].a + mmbonus;
		pga = s[len_b].ga;
		pgb = s[len_b].gb;
		s[len_b].a = -INFTY;
		s[len_b].ga = -INFTY;
		
		s[len_b].gb = pa+prof1[29];
		if(pgb+prof1[29] > s[len_b].gb){
			s[len_b].gb = pgb+prof1[29];
		}
	
		tracep[len_b] = 16;
		
		j = len_b;
		prof2 += len_b << 6;
		while(--j){
			prof2 -= 64;
			ca = s[j].a;

			c = 1;
			if((pga += prof2[91]) > pa){
				pa = pga;
				c = 2;
			}
			if((pgb += prof1[91]) > pa){
				pa = pgb;
				c = 4;
			}

			prof2 += 32;
			for (pga = freq[0];--pga;){
				pgb = freq[pga];
				pa += prof1[pgb]*prof2[pgb];
			}
			prof2 -= 32;

			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = s[j+1].a+prof2[27];
			if (s[j+1].ga+prof2[28] > s[j].ga){
				s[j].ga = s[j+1].ga+prof2[28];
				c |= 8;
			}
			
			pgb = s[j].gb;
			
			s[j].gb = ca+prof1[27];
			if(pgb+prof1[28] > s[j].gb){
				s[j].gb = pgb+prof1[28];
				c |= 16;
			}
			tracep[j] = c;
			pa = ca+ mmbonus;

		}
	
		prof2 -= 64;
		//LAST CELL (0)
		ca = s[0].a;

		c = 1;
		if((pga+=prof2[91]) > pa){
			pa = pga;
			c = 2;
		}
		if((pgb+=prof1[91]) > pa){
			pa = pgb;
			c = 4;
		}
		
		prof2 += 32;
		for (pga = freq[0];--pga;){
			pgb = freq[pga];
			pa += prof1[pgb]*prof2[pgb];
		}
		prof2 -= 32;
		
		s[0].a = pa;
		
		s[0].ga = -INFTY;
		
		pgb = s[0].gb;
		s[0].gb = ca+prof1[27]+prof1[29];
 		if(pgb+prof1[29] > s[0].gb){
			s[0].gb = pgb+prof1[29];
			c |= 16;
		}
		tracep[0] = c;	
		
	}
	prof1 -= 64;
	
	c = 1;
	for (j = 26; j--;){
		if(prof1[j]){
			freq[c] = j;
			c++;	
		}
	}
	freq[0] = c;
	
	tracep = trace[0];
	j = len_b;
	prof2 += len_b << 6;
	pa = s[j].a+ mmbonus;
	pga = s[j].ga;
	pgb = s[j].gb;
	s[j].a = -INFTY;
	s[j].ga = -INFTY;

	s[len_b].gb = pa+prof1[29];
	if(pgb+prof1[29] > s[len_b].gb){
		s[len_b].gb = pgb+prof1[29];
	}


	
	while(--j){
		prof2 -= 64;
		ca = s[j].a;

		c = 1;

		if((pga+=prof2[91]) > pa){
			pa = pga;
			c = 2;
		}

		if((pgb+=prof1[91]) > pa){
			pa = pgb;
			c = 4;
		}	
		prof2+=32;
		
		for (pga = freq[0];--pga;){
			pgb = freq[pga];
			pa += prof1[pgb]*prof2[pgb];
		}
		prof2-=32;
		
		s[j].a = pa;
		pga = s[j].ga;
		s[j].ga = s[j+1].a+prof2[27]+prof2[29];
		if (s[j+1].ga+prof2[29] > s[j].ga){
			s[j].ga = s[j+1].ga+prof2[29];
			c |= 8;
		}	
		pgb = s[j].gb;
		s[j].gb = -INFTY;	
		
		tracep[j] = c;
		pa = ca+ mmbonus;
	}
	prof2 -= 64;

	ca = s[0].a;
	
	c = 1;
	
	if((pga+=prof2[91]) > pa){
		pa = pga;
		c = 2;
	}
	if((pgb+=prof1[91]) > pa){
		pa = pgb;
		c = 4;
	}
		
	prof2+=32;
	for (pga = freq[0];--pga;){	
		pgb = freq[pga];
		pa += prof1[pgb]*prof2[pgb];
	}
	prof2-=32;
	
	s[0].a = pa;
	
	s[0].ga = s[1].a+prof2[27]+prof2[29];
	if (s[1].ga+prof2[29] > s[0].ga){
		s[0].ga = s[1].ga+prof2[29];
		c |= 8;
	}
	
	pgb = s[0].gb;
	s[0].gb = ca+prof1[27]+prof1[29];
	if(pgb +prof1[29]> s[0].gb){
		s[0].gb = pgb+prof1[29];
		c |= 16;
	}	
	tracep[0] = c;

	pgb = s[0].gb;
	c = 2;
	if(s[0].ga > pgb){
		pgb = s[0].ga;
		c = 1;
	}
	if(s[0].a >= pgb){
		pgb = s[0].a;
		c = 0;
	}
	
	//fprintf(stderr,"SCORE:%d\n",ca);
	ca = c;
	
	i = 0;
	j = 0;
	c = 1;
	while(trace[i][j] < 32){
	//	fprintf(stderr,"%d->%d	%d:%d	%d:%d\n",c,trace[i][j],i,j,len_a,len_b);
		switch(ca){
			case 0:
				if (trace[i][j] & 2){
					ca = 1;
					if(i+1!= len_a){
						path[c+1] |= 16;
	//					fprintf(stderr,"GAP_CLOSE\n");
					}else{
						path[c+1] |= 32+16;
					}
				}else if (trace[i][j] & 4){
					ca = 2;
					if(j+1!= len_b){
						path[c+1] |= 16;
	//					fprintf(stderr,"GAP_CLOSE\n");
					}else{
						path[c+1] |= 32+16;
					}
				}

				//path[c] = 0;
				i++;
				j++;
			break;
			case 1:
				if(trace[i][j] & 8){
					ca = 1;
					if(i!=0 && i!= len_a){
	//				/	fprintf(stderr,"GAP_EXT\n");
						if(!(path[c]&16)){
							path[c] |= 8;
						}
					}else{
						if(!(path[c]&16)){
							path[c] |= 32+8;
						}
					}
				}else{
					ca = 0;
					if(i!=0 && i!= len_a){
	//					fprintf(stderr,"GAP_OPEN\n");
						path[c] |= 4;
					}else{
						path[c] |= 32+4;
					}
				}
				path[c] |= 1;
				j++;
			break;
			case  2:
				if(trace[i][j] & 16){
					ca = 2;
					if(j !=0 && j != len_b){
	//					fprintf(stderr,"GAP_EXT\n");
						if(!(path[c]&16)){
							path[c] |= 8;
						}
					}else{
						if(!(path[c]&16)){
							path[c] |= 32+8;
						}
					}
				}else{
					ca = 0;
					if(j!=0 && j != len_b){
	//					fprintf(stderr,"GAP_OPEN\n");
						path[c] |= 4;
					}else{
						path[c] |= 32+4;
					}
					
				}
				path[c] |= 2;
				i++;
			break;
		}
		c++;
	}
	path[0] = c-1;
	path[c] = 3;
	path[c+1] = pgb;
	
	return path;
}*/





int* pp_dyn(int* path, struct dp_matrix *dp,const float* prof1,const float* prof2,const int len_a,const int len_b)
{
	unsigned int freq[26];
	
	struct states* s = 0;
	char** trace = 0;
	char* tracep = 0;
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	register int i = 0;
	register int j = 0;
	register int c = 0;
	register int f = 0;

	s = dp->s;
	
	trace = dp->tb;

	trace[len_a][len_b] = 32;

	prof1 +=  len_a << 6;

	s[len_b].a = 0.0;
	s[len_b].ga = -FLOATINFTY;
	s[len_b].gb = -FLOATINFTY;
	//init of first row;
	tracep = trace[len_a];
	
	j = len_b;
	while(--j){
		s[j].a = -FLOATINFTY;
		
		s[j].ga = s[j+1].a+prof2[29];
		if (s[j+1].ga+prof2[29] > s[j].ga){
			s[j].ga = s[j+1].ga+prof2[29];
		}
		s[j].gb = -INFTY;
		tracep[j] = 8;
	}
	
	s[0].a = -FLOATINFTY;
	s[0].ga = -FLOATINFTY;
	s[0].gb = -FLOATINFTY;
	i = len_a;
	while(--i){
		prof1 -= 64;

		c = 1;
		for (j = 26; j--;){
			if(prof1[j]){
				freq[c] = j;
				c++;	
			}
		}
		freq[0] = c;
		
		tracep = trace[i];
		pa = s[len_b].a;
		pga = s[len_b].ga;
		pgb = s[len_b].gb;
		s[len_b].a = -FLOATINFTY;
		s[len_b].ga = -FLOATINFTY;
		
		s[len_b].gb = pa+prof1[29];
		if(pgb+prof1[29] > s[len_b].gb){
			s[len_b].gb = pgb+prof1[29];
		}
	
		tracep[len_b] = 16;
		
		j = len_b;
		prof2 += len_b << 6;
		while(--j){
			prof2 -= 64;
			ca = s[j].a;

			c = 1;
			if((pga += prof2[91]) > pa){
				pa = pga;
				c = 2;
			}
			if((pgb += prof1[91]) > pa){
				pa = pgb;
				c = 4;
			}
			
			prof2 += 32;
			for (f = freq[0];--f;){
				pa += prof1[freq[f]]*prof2[freq[f]];
			}
			prof2 -= 32;

			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = s[j+1].a+prof2[27];
			if (s[j+1].ga+prof2[28] > s[j].ga){
				s[j].ga = s[j+1].ga+prof2[28];
				c |= 8;
			}
			
			pgb = s[j].gb;
			
			s[j].gb = ca+prof1[27];
			if(pgb+prof1[28] > s[j].gb){
				s[j].gb = pgb+prof1[28];
				c |= 16;
			}
			tracep[j] = c;
			pa = ca;

		}
	
		prof2 -= 64;
		//LAST CELL (0)
		ca = s[0].a;

		c = 1;
		if((pga+=prof2[91]) > pa){
			pa = pga;
			c = 2;
		}
		if((pgb+=prof1[91]) > pa){
			pa = pgb;
			c = 4;
		}
		
		prof2 += 32;
		for (f = freq[0];--f;){
			pa += prof1[freq[f]]*prof2[freq[f]];
		}
		prof2 -= 32;
		
		s[0].a = pa;
		
		s[0].ga = -FLOATINFTY;
		
		pgb = s[0].gb;
		s[0].gb = ca+prof1[27]+prof1[29];
 		if(pgb+prof1[29] > s[0].gb){
			s[0].gb = pgb+prof1[29];
			c |= 16;
		}
		tracep[0] = c;	
		
	}
	prof1 -= 64;
	
	c = 1;
	for (j = 26; j--;){
		if(prof1[j]){
			freq[c] = j;
			c++;	
		}
	}
	freq[0] = c;
	
	tracep = trace[0];
	j = len_b;
	prof2 += len_b << 6;
	pa = s[j].a;
	pga = s[j].ga;
	pgb = s[j].gb;
	s[j].a = -FLOATINFTY;
	s[j].ga = -FLOATINFTY;

	s[len_b].gb = pa+prof1[29];
	if(pgb+prof1[29] > s[len_b].gb){
		s[len_b].gb = pgb+prof1[29];
	}


	
	while(--j){
		prof2 -= 64;
		ca = s[j].a;

		c = 1;

		if((pga+=prof2[91]) > pa){
			pa = pga;
			c = 2;
		}

		if((pgb+=prof1[91]) > pa){
			pa = pgb;
			c = 4;
		}
		
		prof2+=32;
		
		for (f = freq[0];--f;){
			pa += prof1[freq[f]]*prof2[freq[f]];
		}
		prof2-=32;
		
		s[j].a = pa;
		pga = s[j].ga;
		s[j].ga = s[j+1].a+prof2[27]+prof2[29];
		if (s[j+1].ga+prof2[29] > s[j].ga){
			s[j].ga = s[j+1].ga+prof2[29];
			c |= 8;
		}	
		pgb = s[j].gb;
		s[j].gb = -FLOATINFTY;	
		
		tracep[j] = c;
		pa = ca;
	}
	prof2 -= 64;

	ca = s[0].a;
	
	c = 1;
	
	if((pga+=prof2[91]) > pa){
		pa = pga;
		c = 2;
	}
	if((pgb+=prof1[91]) > pa){
		pa = pgb;
		c = 4;
	}
	prof2+=32;
	for (f = freq[0];--f;){
		pa += prof1[freq[f]]*prof2[freq[f]];
	}
	prof2-=32;
	
	s[0].a = pa;
	
	s[0].ga = s[1].a+prof2[27]+prof2[29];
	if (s[1].ga+prof2[29] > s[0].ga){
		s[0].ga = s[1].ga+prof2[29];
		c |= 8;
	}
	
	pgb = s[0].gb;
	s[0].gb = ca+prof1[27]+prof1[29];
	if(pgb +prof1[29]> s[0].gb){
		s[0].gb = pgb+prof1[29];
		c |= 16;
	}	
	tracep[0] = c;

	pgb = s[0].gb;
	c = 2;
	if(s[0].ga > pgb){
		pgb = s[0].ga;
		c = 1;
	}
	if(s[0].a >= pgb){
		pgb = s[0].a;
		c = 0;
	}
	
	//fprintf(stderr,"SCORE:%d\n",ca);
	f = c;
	
	i = 0;
	j = 0;
	c = 1;
	while(trace[i][j] < 32){
	//	fprintf(stderr,"%d->%d	%d:%d	%d:%d\n",c,trace[i][j],i,j,len_a,len_b);
		switch(f){
			case 0:
				if (trace[i][j] & 2){
					f = 1;
					if(i+1!= len_a){
						path[c+1] |= 16;
	//					fprintf(stderr,"GAP_CLOSE\n");
					}else{
						path[c+1] |= 32+16;
					}
				}else if (trace[i][j] & 4){
					f = 2;
					if(j+1!= len_b){
						path[c+1] |= 16;
	//					fprintf(stderr,"GAP_CLOSE\n");
					}else{
						path[c+1] |= 32+16;
					}
				}

				//path[c] = 0;
				i++;
				j++;
			break;
			case 1:
				if(trace[i][j] & 8){
					f = 1;
					if(i!=0 && i!= len_a){
	//				/	fprintf(stderr,"GAP_EXT\n");
						if(!(path[c]&16)){
							path[c] |= 8;
						}
					}else{
						if(!(path[c]&16)){
							path[c] |= 32+8;
						}
					}
				}else{
					f = 0;
					if(i!=0 && i!= len_a){
	//					fprintf(stderr,"GAP_OPEN\n");
						path[c] |= 4;
					}else{
						path[c] |= 32+4;
					}
				}
				path[c] |= 1;
				j++;
			break;
			case  2:
				if(trace[i][j] & 16){
					f = 2;
					if(j !=0 && j != len_b){
	//					fprintf(stderr,"GAP_EXT\n");
						if(!(path[c]&16)){
							path[c] |= 8;
						}
					}else{
						if(!(path[c]&16)){
							path[c] |= 32+8;
						}
					}
				}else{
					f = 0;
					if(j!=0 && j != len_b){
	//					fprintf(stderr,"GAP_OPEN\n");
						path[c] |= 4;
					}else{
						path[c] |= 32+4;
					}
					
				}
				path[c] |= 2;
				i++;
			break;
		}
		c++;
	}
	path[0] = c-1;
	path[c] = 3;
	path[c+1] = pgb;
	return path;
}


int* ps_dyn(int* path, struct dp_matrix *dp,const float* prof1,const int* seq2,const int len_a,const int len_b,int sip)
{
	struct states* s = 0;
	char** trace = 0;
	char* tracep = 0;
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	register int i = 0;
	register int j = 0;
	register int c = 0;
	register int f = 0;
	
	const float open = gpo * sip;
	const float ext = gpe *sip;

	s = dp->s;
	
	trace = dp->tb;

	trace[len_a][len_b] = 32;

	prof1 +=  len_a << 6;

	s[len_b].a = 0.0;
	s[len_b].ga = -FLOATINFTY;
	s[len_b].gb = -FLOATINFTY;
	//init of first row;
	tracep = trace[len_a];
	j = len_b;
	

	while(--j){
		s[j].a = -FLOATINFTY;
		//s[j].ga = 0;	
		
		s[j].ga = s[j+1].a-tgpe;//-topen;
		if (s[j+1].ga-tgpe > s[j].ga){
			s[j].ga = s[j+1].ga-tgpe;
		}
		
		s[j].gb = -FLOATINFTY;
		tracep[j] = 8;
	}
	
	s[0].a = -FLOATINFTY;
	s[0].ga = -FLOATINFTY;
	s[0].gb = -FLOATINFTY;
	i = len_a;
	while(--i){
		prof1 -= 64;
		
		tracep = trace[i];
		pa = s[len_b].a;
		pga = s[len_b].ga;
		pgb = s[len_b].gb;
		s[len_b].a = -FLOATINFTY;
		s[len_b].ga = -FLOATINFTY;
		//s[len_b].gb = 0;
		s[len_b].gb = pa+prof1[29];//+prof1[29];
		if(pgb+prof1[29] > s[len_b].gb){
			s[len_b].gb = pgb+prof1[29];
		}

		tracep[len_b] = 16;
		
		j = len_b;
		
		while(--j){

			ca = s[j].a;

			c = 1;
			if((pga -= open) > pa){
				pa = pga;
				c = 2;
			}
			if((pgb += prof1[91]) > pa){
				pa = pgb;
				c = 4;
			}
			
			pa += prof1[32 + seq2[j]];

			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = s[j+1].a-open;
			if (s[j+1].ga-ext > s[j].ga){
				s[j].ga = s[j+1].ga-ext;
				c |= 8;
			}
			
			pgb = s[j].gb;
			
			s[j].gb = ca+prof1[27];
			if(pgb+prof1[28] > s[j].gb){
				s[j].gb = pgb+prof1[28];
				c |= 16;
			}
			tracep[j] = c;
			pa = ca;

		}
	
		//LAST CELL (0)
		ca = s[0].a;

		c = 1;
		if((pga-=open) > pa){
			pa = pga;
			c = 2;
		}
		if((pgb+=prof1[91]) > pa){
			pa = pgb;
			c = 4;
		}
		pa += prof1[32+seq2[0]];
		s[0].a = pa;
		
		s[0].ga = -FLOATINFTY;
		
		pgb = s[0].gb;
		s[0].gb = ca+prof1[27]+prof1[29];
 		if(pgb+prof1[29] > s[0].gb){
			s[0].gb = pgb+prof1[29];
			c |= 16;
		}
		tracep[0] = c;	
		
	}
	prof1 -= 64;
	

	
	tracep = trace[0];
	j = len_b;
	pa = s[j].a;
	pga = s[j].ga;
	pgb = s[j].gb;
	s[j].a = -FLOATINFTY;
	s[j].ga = -FLOATINFTY;
	//s[j].gb = -INFTY;
	s[len_b].gb = pa+prof1[29];//+prof1[29];
	if(pgb+prof1[29] > s[len_b].gb){
		s[len_b].gb = pgb+prof1[29];
	}

	
	while(--j){

		ca = s[j].a;

		c = 1;

		if((pga-=open) > pa){
			pa = pga;
			c = 2;
		}

		if((pgb+=prof1[91]) > pa){
			pa = pgb;
			c = 4;
		}
		pa += prof1[32+seq2[j]];		
		s[j].a = pa;
		pga = s[j].ga;
		s[j].ga = s[j+1].a-(open+tgpe);
		if (s[j+1].ga-tgpe > s[j].ga){
			s[j].ga = s[j+1].ga-tgpe;
			c |= 8;
		}	
		pgb = s[j].gb;
		s[j].gb = -INFTY;	
		
		tracep[j] = c;
		pa = ca;
	}


	ca = s[0].a;
	
	c = 1;
	
	if((pga-=open) > pa){
		pa = pga;
		c = 2;
	}
	if((pgb+=prof1[91]) > pa){
		pa = pgb;
		c = 4;
	}
	pa += prof1[32+seq2[0]];	
	s[0].a = pa;
	
	s[0].ga = s[1].a-(open+tgpe);
	if (s[1].ga-tgpe > s[0].ga){
		s[0].ga = s[1].ga-tgpe;
		c |= 8;
	}
	
	pgb = s[0].gb;
	s[0].gb = ca+prof1[27]+prof1[29];
	if(pgb+prof1[29] > s[0].gb){
		s[0].gb = pgb+prof1[29];
		c |= 16;
	}	
	tracep[0] = c;


	pgb = s[0].gb;
	c = 2;
	if(s[0].ga > pgb){
		pgb = s[0].ga;
		c = 1;
	}
	if(s[0].a >= pgb){
		pgb = s[0].a;
		c = 0;
	}
	
	//fprintf(stderr,"SCORE:%d\n",ca);
	f = c;
	
	i = 0;
	j = 0;
	c = 1;
	while(trace[i][j] < 32){
	//	fprintf(stderr,"%d->%d	%d:%d	%d:%d\n",c,trace[i][j],i,j,len_a,len_b);
		switch(f){
			case 0:
				if (trace[i][j] & 2){
					f = 1;
					if(i+1!= len_a){
						path[c+1] |= 16;
	//					fprintf(stderr,"GAP_CLOSE\n");
					}else{
						path[c+1] |= 32+16;
					}
				}else if (trace[i][j] & 4){
					f = 2;
					if(j+1!= len_b){
						path[c+1] |= 16;
	//					fprintf(stderr,"GAP_CLOSE\n");
					}else{
						path[c+1] |= 32+16;
					}
				}

				//path[c] = 0;
				i++;
				j++;
			break;
			case 1:
				if(trace[i][j] & 8){
					f = 1;
					if(i!=0 && i!= len_a){
	//				/	fprintf(stderr,"GAP_EXT\n");
						if(!(path[c]&16)){
							path[c] |= 8;
						}
					}else{
						if(!(path[c]&16)){
							path[c] |= 32+8;
						}
					}
				}else{
					f = 0;
					if(i!=0 && i!= len_a){
	//					fprintf(stderr,"GAP_OPEN\n");
						path[c] |= 4;
					}else{
						path[c] |= 32+4;
					}
				}
				path[c] |= 1;
				j++;
			break;
			case  2:
				if(trace[i][j] & 16){
					f = 2;
					if(j !=0 && j != len_b){
	//					fprintf(stderr,"GAP_EXT\n");
						if(!(path[c]&16)){
							path[c] |= 8;
						}
					}else{
						if(!(path[c]&16)){
							path[c] |= 32+8;
						}
					}
				}else{
					f = 0;
					if(j!=0 && j != len_b){
	//					fprintf(stderr,"GAP_OPEN\n");
						path[c] |= 4;
					}else{
						path[c] |= 32+4;
					}
					
				}
				path[c] |= 2;
				i++;
			break;
		}
		c++;
	}
	path[0] = c-1;
	path[c] = 3;
	path[c+1] = pgb;
	return path;
}

int* ss_dyn(float**subm,int* path, struct dp_matrix *dp,const int* seq1,const int* seq2,const int len_a,const int len_b)
{
	struct states* s = 0;
	const float *subp = 0;
	char** trace = 0;
	char* tracep = 0;
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	register int i = 0;
	register int j = 0;
	register int c = 0;
	register int f = 0;
	
	s = dp->s;
	
	trace = dp->tb;

	trace[len_a][len_b] = 32;

	s[len_b].a = 0.0;
	s[len_b].ga = -FLOATINFTY;
	s[len_b].gb = -FLOATINFTY;
	
	//init of first row;
	tracep = trace[len_a];
	j = len_b;
	

	while(--j){
		s[j].a = -FLOATINFTY;
		//s[j].ga = 0;	
		s[j].ga = s[j+1].a-tgpe;//-gpo;
		if (s[j+1].ga-tgpe > s[j].ga){
			s[j].ga = s[j+1].ga-tgpe;
		}
		
		
		
		s[j].gb = -FLOATINFTY;
		tracep[j] = 8;
	}
	
	s[0].a = -FLOATINFTY;
	s[0].ga = -FLOATINFTY;
	s[0].gb = -FLOATINFTY;
	
	i = len_a;
	while(--i){
		
		tracep = trace[i];
		pa = s[len_b].a;
		pga = s[len_b].ga;
		pgb = s[len_b].gb;
		
		s[len_b].a = -FLOATINFTY;
		s[len_b].ga = -FLOATINFTY;
		//s[len_b].gb = 0;
		
		s[len_b].gb = pa-tgpe;//-gpo;
		if(pgb-tgpe > s[len_b].gb){
			s[len_b].gb = pgb-tgpe;
		}
		
		
		tracep[len_b] = 16;
		j = len_b;
		subp = subm[seq1[i]];
		while(--j){
			ca = s[j].a;
			
			c = 1;
			if((pga -= gpo) > pa){
				pa = pga;
				c = 2;
			}
			if((pgb -= gpo) > pa){
				pa = pgb;
				c = 4;
			}

			pa += subp[seq2[j]];

			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = s[j+1].a-gpo;
			if (s[j+1].ga-gpe > s[j].ga){
				s[j].ga = s[j+1].ga-gpe;
				c |= 8;
			}
			
			pgb = s[j].gb;
			
			s[j].gb = ca-gpo;
			if(pgb-gpe > s[j].gb){
				s[j].gb = pgb-gpe;
				c |= 16;
			}
			tracep[j] = c;
			pa = ca;

		}
	
		//LAST CELL (0)
		ca = s[0].a;

		c = 1;
		if((pga-=gpo) > pa){
			pa = pga;
			c = 2;
		}
		if((pgb-=gpo) > pa){
			pa = pgb;
			c = 4;
		}
		
		pa += subp[seq2[0]];
		
		s[0].a = pa;
		
		s[0].ga = -FLOATINFTY;
		
		pgb = s[0].gb;
		s[0].gb = ca-(gpo+tgpe);
 		if(pgb-tgpe > s[0].gb){
			s[0].gb = pgb-tgpe;
			c |= 16;
		}
		tracep[0] = c;			
	}

	subp = subm[seq1[0]];
	tracep = trace[0];
	j = len_b;
	pa = s[j].a;
	pga = s[j].ga;
	pgb = s[j].gb;
	s[j].a = -FLOATINFTY;
	s[j].ga = -FLOATINFTY;
	
	s[j].gb = pa-tgpe;//-gpo;
	if(pgb-tgpe > s[j].gb){
		s[j].gb = pgb-tgpe;
	}
	
	//s[j].gb = -INFTY;
	while(--j){

		ca = s[j].a;

		c = 1;

		if((pga-=gpo) > pa){
			pa = pga;
			c = 2;
		}

		if((pgb-=gpo) > pa){
			pa = pgb;
			c = 4;
		}
		
		pa += subp[seq2[j]];
		
		s[j].a = pa;
		
		pga = s[j].ga;
		s[j].ga = s[j+1].a-(gpo+tgpe);
		if (s[j+1].ga-tgpe > s[j].ga){
			s[j].ga = s[j+1].ga-tgpe;
			c |= 8;
		}	
		pgb = s[j].gb;
		s[j].gb = -FLOATINFTY;	
		tracep[j] = c;
		pa = ca;
	}
	
	ca = s[0].a;
	
	c = 1;
	
	if((pga-=gpo) > pa){
		pa = pga;
		c = 2;
	}
	if((pgb-=gpo) > pa){
		pa = pgb;
		c = 4;
	}

	pa += subp[seq2[0]];
	
	s[0].a = pa;
	
	
	s[0].ga = s[1].a-(gpo+tgpe);
	if (s[1].ga-tgpe > s[0].ga){
		s[0].ga = s[1].ga-tgpe;
		c |= 8;
	}
	
	pgb = s[0].gb;
	s[0].gb = ca-(gpo+tgpe);
	if(pgb-tgpe > s[0].gb){
		s[0].gb = pgb-tgpe;
		c |= 16;
	}	
	tracep[0] = c;


	pgb = s[0].gb;
	c = 2;
	if(s[0].ga > pgb){
		pgb = s[0].ga;
		c = 1;
	}
	if(s[0].a >= pgb){
		pgb = s[0].a;
		c = 0;
	}
	
	f = c;
	
	i = 0;
	j = 0;
	c = 1;
	while(trace[i][j] < 32){
	//	fprintf(stderr,"%d->%d	%d:%d	%d:%d\n",c,trace[i][j],i,j,len_a,len_b);
		switch(f){
			case 0:
				if (trace[i][j] & 2){
					f = 1;
					if(i+1!= len_a){
						path[c+1] |= 16;
	//					fprintf(stderr,"GAP_CLOSE\n");
					}else{
						path[c+1] |= 32+16;
					}
				}else if (trace[i][j] & 4){
					f = 2;
					if(j+1!= len_b){
						path[c+1] |= 16;
	//					fprintf(stderr,"GAP_CLOSE\n");
					}else{
						path[c+1] |= 32+16;
					}
				}

				//path[c] = 0;
				i++;
				j++;
			break;
			case 1:
				if(trace[i][j] & 8){
					f = 1;
					if(i!=0 && i!= len_a){
	//				/	fprintf(stderr,"GAP_EXT\n");
						if(!(path[c]&16)){
							path[c] |= 8;
						}
					}else{
						if(!(path[c]&16)){
							path[c] |= 32+8;
						}
					}
				}else{
					f = 0;
					if(i!=0 && i!= len_a){
	//					fprintf(stderr,"GAP_OPEN\n");
						path[c] |= 4;
					}else{
						path[c] |= 32+4;
					}
				}
				path[c] |= 1;
				j++;
			break;
			case  2:
				if(trace[i][j] & 16){
					f = 2;
					if(j !=0 && j != len_b){
	//					fprintf(stderr,"GAP_EXT\n");
						if(!(path[c]&16)){
							path[c] |= 8;
						}
					}else{
						if(!(path[c]&16)){
							path[c] |= 32+8;
						}
					}
				}else{
					f = 0;
					if(j!=0 && j != len_b){
	//					fprintf(stderr,"GAP_OPEN\n");
						path[c] |= 4;
					}else{
						path[c] |= 32+4;
					}
					
				}
				path[c] |= 2;
				i++;
			break;
		}
		c++;
	}
	path[0] = c-1;
	path[c] = 3;
	path[c+1] = pgb;
	return path;
}
