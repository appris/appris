/*
	kalign2_hirschberg_large.c
	
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
#include "kalign2_hirschberg_large.h"
#define MAX(a, b) (a > b ? a : b)
#define MAX3(a,b,c) MAX(MAX(a,b),c)
//#include <emmintrin.h>

float local_gpo;
float local_gpe;
float local_tgpe;

int** hirschberg_large_alignment(struct alignment* aln,int* tree,float**submatrix, int** map,int window,float strength)
{
	struct hirsch_large_mem* hm = 0;
	int i,j,g,a,b,c;
	int len_a;
	int len_b;
	

	
	float** profile = 0;
	
	float** subm = 0;
	subm = malloc(sizeof(float*)*32);
	for(i = 0; i < 32;i++){
		subm[i] = malloc(sizeof(float)*32);
		for (j = 0; j < 32;j++){
			subm[i][j] = (float)submatrix[i][j];
		}
	}
	local_gpo = (float)gpo;
	local_gpe = (float)gpe;
	local_tgpe = (float)tgpe;
	

	profile = malloc(sizeof(float*)*numprofiles);
	for ( i = 0;i< numprofiles;i++){
		profile[i] = 0;
	}

	map = malloc(sizeof(int*)*numprofiles);
	for ( i = 0;i < numprofiles;i++){
		map[i] = 0;
	}
	
	hm = hirsch_large_mem_alloc(hm,1024);

	fprintf(stderr,"\nAlignment:\n");

	for (i = 0; i < (numseq-1);i++){
		a = tree[i*3];
		b = tree[i*3+1];
		c = tree[i*3+2];
		fprintf(stderr,"\r%8.0f percent done",(float)(i) /(float)numseq * 100);
		//fprintf(stderr,"Aligning:%d %d->%d	done:%f\n",a,b,c,((float)(i+1)/(float)numseq)*100);
		len_a = aln->sl[a];
		len_b = aln->sl[b];

		
		g = (len_a > len_b)? len_a:len_b;
		map[c] = malloc(sizeof(int) * (g+2));
		if(g > hm->size){
			hm = hirsch_large_mem_realloc(hm,g);
		}

		for (j = 0; j < (g+2);j++){
			map[c][j] = -1;
		}

		if (a < numseq){
			profile[a] = make_large_profile(profile[a],aln->s[a],len_a,subm);
		}else{
			set_large_gap_penalties(profile[a],len_a,aln->nsip[b]);
		}
		if (b < numseq){
			profile[b] = make_large_profile(profile[b],aln->s[b],len_b,subm);
		}else{		
			set_large_gap_penalties(profile[b],len_b,aln->nsip[a]);
		}
		
		hm->starta = 0;
		hm->startb = 0;
		hm->enda = len_a;
		hm->endb = len_b;
		hm->len_a = len_a;
		hm->len_b = len_b;
		
		hm->f[0].a = 0.0;
		hm->f[0].ga =  -FLOATINFTY;
		hm->f[0].gb = -FLOATINFTY;
		hm->b[0].a = 0.0;
		hm->b[0].ga =  -FLOATINFTY;
		hm->b[0].gb =  -FLOATINFTY;
	//	fprintf(stderr,"LENA:%d	LENB:%d	numseq:%d\n",len_a,len_b,numseq);
		if(a < numseq){
			if(b < numseq){
				map[c] = hirsch_large_ss_dyn(subm,aln->s[a],aln->s[b],hm,map[c]);
			}else{
				hm->enda = len_b;
				hm->endb = len_a;
				hm->len_a = len_b;
				hm->len_b = len_a;
				map[c] = hirsch_large_ps_dyn(profile[b],aln->s[a],hm,map[c],aln->nsip[b]);
				map[c] = mirror_hirsch_path(map[c],len_a,len_b);
			}
		}else{
			if(b < numseq){
				map[c] = hirsch_large_ps_dyn(profile[a],aln->s[b],hm,map[c],aln->nsip[a]);
			}else{
				if(len_a < len_b){
					map[c] = hirsch_large_pp_dyn(profile[a],profile[b],hm,map[c]);
				}else{
					hm->enda = len_b;
					hm->endb = len_a;
					hm->len_a = len_b;
					hm->len_b = len_a;
					map[c] = hirsch_large_pp_dyn(profile[b],profile[a],hm,map[c]);
					map[c] = mirror_hirsch_path(map[c],len_a,len_b);
				}
			}
		}
		
		map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);

		if(i != numseq-2){
			profile[c] = malloc(sizeof(float)*64*(map[c][0]+2));
			profile[c] = large_update(profile[a],profile[b],profile[c],map[c],aln->nsip[a],aln->nsip[b]);
		}
			
		aln->sl[c] = map[c][0];
	
		aln->nsip[c] = aln->nsip[a] + aln->nsip[b];
		aln->sip[c] = malloc(sizeof(int)*(aln->nsip[a] + aln->nsip[b]));
		g =0;
		for (j = aln->nsip[a];j--;){
			aln->sip[c][g] = aln->sip[a][j];
			g++;
		}
		for (j = aln->nsip[b];j--;){
			aln->sip[c][g] = aln->sip[b][j];
			g++;
		}

		free(profile[a]);
		free(profile[b]);
	}
	fprintf(stderr,"\r%8.0f percent done\n",100.0);
	free(profile);
	hirsch_large_mem_free(hm);
	for (i = 32;i--;){
		free(subm[i]);
		free(submatrix[i]);
	}
	free(subm);
	free(submatrix);
	return map;
}


int* hirsch_large_ss_dyn(float**subm, const int* seq1,const int* seq2,struct hirsch_large_mem* hm, int* hirsch_path)
{
	int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
	float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
	int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

	if(hm->starta  >= hm->enda){
		return hirsch_path;
	}
	if(hm->startb  >= hm->endb){
		return hirsch_path;
	}


	hm->enda = mid;

	//fprintf(stderr,"Forward:%d-%d	%d-%d\n",hm->starta,hm->enda,hm->startb,hm->endb);
	hm->f = foward_large_hirsch_ss_dyn(subm,seq1,seq2,hm);

	hm->starta = mid;
	hm->enda = old_cor[1];
	//fprintf(stderr,"Backward:%d-%d	%d-%d\n",hm->starta,hm->enda,hm->startb,hm->endb);
	hm->b = backward_large_hirsch_ss_dyn(subm,seq1,seq2,hm);


	hirsch_path = hirsch_large_align_two_ss_vector(subm,seq1,seq2,hm,hirsch_path,input_states,old_cor);
	return hirsch_path;
}

int* hirsch_large_align_two_ss_vector(float**subm,const int* seq1,const int* seq2,struct hirsch_large_mem* hm,int* hirsch_path,float input_states[],int old_cor[])
{
	struct large_states* f = hm->f;
 	struct large_states* b = hm->b;
	int i,j,c;
	int transition = -1;
	
	
	//code:
	// a -> a = 1
	// a -> ga = 2
	// a -> gb = 3
	// ga ->ga = 4
	// ga -> a = 5
	//gb->gb = 6;
	//gb->a = 7;
	
	//int max = -FLOATINFTY;
	float max = -FLOATINFTY;	
	//float middle =  (hm->endb - hm->startb)/2 + hm->startb;
	float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
	float sub = 0.0;
	
	//i = hm->startb;
	i = old_cor[2];
	c = -1;
	//for(i = hm->startb; i < hm->endb;i++){
	for(i = old_cor[2]; i < old_cor[3];i++){
	
		sub = abs(middle -i);
		sub /= 1000; 
	//	fprintf(stderr,"%d-%d	%f\n",hm->startb,hm->endb,sub);
		if(f[i].a+b[i].a-sub > max){
			max = f[i].a+b[i].a-sub;
	//		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
			transition = 1;
			c = i;
		}
		if(f[i].a+b[i].ga-local_gpo-sub > max){
			max = f[i].a+b[i].ga-local_gpo-sub;
	//		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
			transition = 2;
			c = i;
		}
		if(f[i].a+b[i].gb -local_gpo -sub > max){
			max = f[i].a+b[i].gb - local_gpo-sub;
	//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
			transition = 3;
			c = i;
		}
		if(f[i].ga+b[i].a - local_gpo-sub > max){
			max = f[i].ga+b[i].a - local_gpo-sub;
	//		fprintf(stderr,"gap_a->aligned:%d + %d + %d(local_gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
			transition = 5;
			c = i;
		}


		if(hm->startb == 0){
			if(f[i].gb+b[i].gb - local_tgpe-sub > max){
				max = f[i].gb+b[i].gb -local_tgpe-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}else{
			if(f[i].gb+b[i].gb - local_gpe -sub> max){
				max = f[i].gb+b[i].gb - local_gpe-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}
		if(f[i].gb+b[i].a - local_gpo-sub > max){
			max = f[i].gb+b[i].a - local_gpo-sub;
	//		fprintf(stderr,"gap_b->aligned:%d + %d + %d(local_gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
			transition = 7;
			c = i;
		}
	}
	//i = hm->endb;
	i = old_cor[3];
	sub = abs(middle -i);
	sub /= 1000; 
	
	if(f[i].a+b[i].gb-local_gpo-sub > max){
		max = f[i].a+b[i].gb - local_gpo-sub;
	//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
		transition = 3;
		c = i;
	}
	if(hm->endb == hm->len_b){
		if(f[i].gb+b[i].gb -local_tgpe-sub > max){
			max = f[i].gb+b[i].gb - local_tgpe-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}	
	}else{
		if(f[i].gb+b[i].gb - gpe-sub > max){
			max = f[i].gb+b[i].gb - gpe-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(local_gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}
	}
	
	
	//fprintf(stderr,"Transition:%d	at:%d\n",transition,c);
	
	j = hirsch_path[0];
	switch(transition){
		case 1: //a -> a = 1
			
			hirsch_path[old_cor[4]] = c;
			hirsch_path[old_cor[4]+1] = c+1;
			
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
	//		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			break;
		case 2:// a -> ga = 2
			
			hirsch_path[old_cor[4]] = c;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			

			//backward:
			hm->starta = old_cor[4];
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = 0;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			break;
		case 3:// a -> gb = 3
			
			hirsch_path[old_cor[4]] = c;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_large_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			break;
		case 5://ga -> a = 5
			hirsch_path[old_cor[4]+1] = c+1;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = 0;
			hm->b[0].gb = -FLOATINFTY;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4];
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_large_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			break;
		case 6://gb->gb = 6;
			
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_large_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			break;
		case 7://gb->a = 7;
			
			hirsch_path[old_cor[4]+1] = c+1;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_large_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			break;
	}
		
	return hirsch_path;
}



struct large_states* foward_large_hirsch_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_large_mem* hm)
{
	struct large_states* s = hm->f;
	float *subp = 0;
	const int starta = hm->starta;
	const int enda = hm->enda;
	const int startb = hm->startb;
	const int endb = hm->endb;
	
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	register int i = 0;
	register int j = 0;
	

	s[startb].a = s[0].a;
	s[startb].ga = s[0].ga;
	s[startb].gb = s[0].gb;
	if(startb == 0){
		for (j = startb+1; j < endb;j++){

			s[j].a = -FLOATINFTY;
			
			//s[j].ga = s[j-1].a-local_tgpe;
			//if (s[j-1].ga-local_tgpe > s[j].ga){
			//	s[j].ga = s[j-1].ga-local_tgpe;
			//}
			//if(s[j-1].ga > s[j-1].a){
			//	s[j].ga = s[j-1].ga-local_tgpe;
			//}else{
			//	s[j].ga = s[j-1].a-local_tgpe;
			//}
			s[j].ga = MAX(s[j-1].ga,s[j-1].a)-local_tgpe;
			
			s[j].gb = -FLOATINFTY;
		}		
	}else{

		for (j = startb+1; j < endb;j++){
			s[j].a = -FLOATINFTY;
			
			//s[j].ga = s[j-1].a-local_gpo;
			//if (s[j-1].ga - gpe > s[j].ga){
			//	s[j].ga = s[j-1].ga-gpe;
			//}
			//if(s[j-1].ga - gpe >s[j-1].a-local_gpo){
			//	s[j].ga = s[j-1].ga-gpe;
			//}else{
			//	s[j].ga = s[j-1].a-local_gpo;
			//}
			s[j].ga = MAX(s[j-1].ga - local_gpe,s[j-1].a-local_gpo);
			
			s[j].gb = -FLOATINFTY;
		}
	}
	s[endb].a = -FLOATINFTY;
	s[endb].ga = -FLOATINFTY;
	s[endb].gb = -FLOATINFTY;

	seq2--;
	for (i = starta;i < enda;i++){
		subp = subm[seq1[i]];

		pa = s[startb].a;
		pga = s[startb].ga;
		pgb = s[startb].gb;
		if(startb == 0){
			s[startb].a = -FLOATINFTY;
			s[startb].ga = -FLOATINFTY;
		
			//s[startb].gb = pa-local_tgpe;
			//if(pgb - local_tgpe > s[startb].gb){
			//	s[startb].gb = pgb-local_tgpe;
			//}
			//if(pgb > pa){
			//	s[startb].gb = pgb-local_tgpe;
			//}else{
			//	s[startb].gb = pa-local_tgpe;
			//}
			s[startb].gb = MAX(pgb,pa) - local_tgpe;
		}else{
			s[startb].a = -FLOATINFTY;
			s[startb].ga = -FLOATINFTY;
		
			//s[startb].gb = pa-local_gpo;
			//if(pgb -gpe > s[startb].gb){
			//	s[startb].gb = pgb -gpe;
			//}
			//if(pgb - gpe > pa - local_gpo){
			//	s[startb].gb = pgb - gpe;
			//}else{
			//	s[startb].gb = pa - local_gpo;
			//}
			s[startb].gb = MAX(pgb - local_gpe,pa - local_gpo);
			
		}
		for (j = startb+1; j <= endb;j++){
			ca = s[j].a;
			//if((pga -= local_gpo) > pa){
			//	pa = pga;
			//}
			//if((pgb -= local_gpo) > pa){
			//	pa = pgb;
			//}
			pa = MAX3(pa,pga-local_gpo,pgb-local_gpo);
			
			pa += subp[seq2[j]];
			
			s[j].a = pa;
			
			pga = s[j].ga;
			
			//s[j].ga = s[j-1].a-local_gpo;
			//if (s[j-1].ga-gpe > s[j].ga){
			//	s[j].ga = s[j-1].ga-gpe;
			//}
			
			//if(s[j-1].ga-gpe >s[j-1].a-local_gpo){
			//	s[j].ga = s[j-1].ga-gpe;
			//}else{
			//	s[j].ga = s[j-1].a-local_gpo;
			//}
			s[j].ga = MAX(s[j-1].ga-local_gpe,s[j-1].a-local_gpo);
			
			pgb = s[j].gb;
			
			//s[j].gb = ca-local_gpo;
			//if(pgb-gpe > s[j].gb){
			//	s[j].gb = pgb-gpe;
			//}
			
			//if(pgb-gpe >  ca-local_gpo){
			//	s[j].gb = pgb-gpe;
			//}else{
			//	s[j].gb = ca-local_gpo;
			//}
			s[j].gb = MAX(pgb-local_gpe ,ca-local_gpo);
			
			pa = ca;
		}
	}
	return s;
}

struct large_states* backward_large_hirsch_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_large_mem* hm)
{

	struct large_states* s = hm->b;
	float *subp = 0;
	const int starta = hm->starta;
	const int enda = hm->enda;
	const int startb = hm->startb;
	const int endb = hm->endb;
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	register int i = 0;
	register int j = 0;

	s[endb].a = s[0].a ;
	s[endb].ga = s[0].ga;
	s[endb].gb = s[0].gb;
	
	
	//init of first row;
	
	//j = endb-startb;
	if(endb == hm->len_b){
		for(j = endb-1;j > startb;j--){
			s[j].a = -FLOATINFTY;
			
			//s[j].ga = s[j+1].a-local_tgpe;
			//if (s[j+1].ga-local_tgpe > s[j].ga){
			//	s[j].ga = s[j+1].ga-local_tgpe;
			//}
			//if(s[j+1].ga > s[j+1].a){
			//	s[j].ga = s[j+1].ga-local_tgpe;
			//}else{
			//	s[j].ga = s[j+1].a-local_tgpe;
			//}
			s[j].ga = MAX(s[j+1].ga,s[j+1].a)-local_tgpe;
			
			s[j].gb = -FLOATINFTY;
		}
	}else{
		for(j = endb-1;j > startb;j--){
			s[j].a = -FLOATINFTY;
			
			//s[j].ga = s[j+1].a-local_gpo;
			//if (s[j+1].ga-gpe > s[j].ga){
			//	s[j].ga = s[j+1].ga-gpe;
			//}
			//if(s[j+1].ga-gpe > s[j+1].a-local_gpo){
			//	s[j].ga = s[j+1].ga-gpe;
			//}else{
			//	s[j].ga = s[j+1].a-local_gpo;
			//}
			s[j].ga = MAX(s[j+1].ga-local_gpe,s[j+1].a-local_gpo);
			
			s[j].gb = -FLOATINFTY;
		}
	}

	
	s[startb].a = -FLOATINFTY;
	s[startb].ga = -FLOATINFTY;
	s[startb].gb = -FLOATINFTY;

	i = enda-starta;
	seq1+= starta;
	while(i--){
		subp = subm[seq1[i]];
		pa = s[endb].a;
		pga = s[endb].ga;
		pgb = s[endb].gb;
		s[endb].a = -FLOATINFTY;
		s[endb].ga = -FLOATINFTY;

		if(endb == hm->len_b){
			//s[endb].gb = pa-local_tgpe;
			//if(pgb-local_tgpe > s[endb].gb){
			//	s[endb].gb = pgb-local_tgpe;
			//}
			//if(pgb > pa){
			//	s[endb].gb = pgb-local_tgpe;
			//}else{
			//	s[endb].gb = pa-local_tgpe;
			//}
			s[endb].gb = MAX(pgb,pa)-local_tgpe;
			
			
		}else{
			//s[endb].gb = pa-local_gpo;
			//if(pgb-gpe > s[endb].gb){
			//	s[endb].gb = pgb-gpe;
			//}
			//if(pgb-gpe  > pa-local_gpo){
			//	s[endb].gb = pgb-gpe;
			//}else{
			//	s[endb].gb = pa-local_gpo;
			//}
			s[endb].gb = MAX(pgb-local_gpe,pa-local_gpo);
			
		}

		for(j = endb-1;j >= startb;j--){

			ca = s[j].a;
			//if((pga -= local_gpo) > pa){
			//	pa = pga;
			//}
			//if((pgb -= local_gpo) > pa){
			//	pa = pgb;
			//}
			pa = MAX3(pa,pga - local_gpo,pgb-local_gpo);
			
			pa += subp[seq2[j]];

			s[j].a = pa;
			
			pga = s[j].ga;
			
			//s[j].ga = s[j+1].a-local_gpo;
			//if (s[j+1].ga-gpe > s[j].ga){
			//	s[j].ga = s[j+1].ga-gpe;
			//}
			
			
			//if(s[j+1].ga-gpe >s[j+1].a-local_gpo){
			//	s[j].ga = s[j+1].ga-gpe;
			//}else{
			//	s[j].ga = s[j+1].a-local_gpo;
			//}
			s[j].ga = MAX(s[j+1].ga-local_gpe,s[j+1].a-local_gpo);
			
			pgb = s[j].gb;
			
			//s[j].gb = ca-local_gpo;
			//if(pgb-gpe > s[j].gb){
			//	s[j].gb = pgb-gpe;
			//}
			//if(pgb-gpe >  ca-local_gpo){
			//	s[j].gb = pgb-gpe;
			//}else{
			//	s[j].gb = ca-local_gpo;
			//}
			s[j].gb = MAX(pgb-local_gpe,ca-local_gpo);
			
			pa = ca;
		}
	}		
	return s;
}


int* hirsch_large_ps_dyn(const float* prof1,const int* seq2,struct hirsch_large_mem* hm, int* hirsch_path,int sip)
{
	int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
	float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
	int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};


	if(hm->starta  >= hm->enda){
		return hirsch_path;
	}
	if(hm->startb  >= hm->endb){
		return hirsch_path;
	}

	hm->enda = mid;
	hm->f = foward_large_hirsch_ps_dyn(prof1,seq2,hm,sip);
	
	/*int i;
	fprintf(stderr,"FOWARD\n");
	for (i = hm->startb; i <= hm->endb;i++){
		fprintf(stderr,"%d	%d	%d\n",hm->f[i].a,hm->f[i].ga,hm->f[i].gb);
	}*/

	hm->starta = mid;
	hm->enda = old_cor[1];
	hm->b = backward_large_hirsch_ps_dyn(prof1,seq2,hm,sip);
	
	/*fprintf(stderr,"BaCKWARD\n");
	for (i = hm->startb; i <= hm->endb;i++){
		fprintf(stderr,"%d	%d	%d\n",hm->b[i].a,hm->b[i].ga,hm->b[i].gb);
	}*/

	hirsch_path = hirsch_large_align_two_ps_vector(prof1,seq2,hm,hirsch_path,input_states,old_cor,sip);
	return hirsch_path;
}



int* hirsch_large_align_two_ps_vector(const float* prof1,const int* seq2,struct hirsch_large_mem* hm,int* hirsch_path,float input_states[],int old_cor[],int sip)
{
	struct large_states* f = hm->f;
 	struct large_states* b = hm->b;
	int i,j,c;
	int transition = -1;
	
	const float open = local_gpo * sip;
	
	
	//code:
	// a -> a = 1
	// a -> ga = 2
	// a -> gb = 3
	// ga ->ga = 4
	// ga -> a = 5
	//gb->gb = 6;
	//gb->a = 7;
	
	//int max = -FLOATINFTY;
	float max = -FLOATINFTY;	
	//float middle =  (hm->endb - hm->startb)/2 + hm->startb;
	float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
	float sub = 0.0;
	
	
	prof1+= ((old_cor[4]+1)<<6);
	
	//i = hm->startb;
	i = old_cor[2];
	c = -1;
	//for(i = hm->startb; i < hm->endb;i++){
	for(i = old_cor[2]; i < old_cor[3];i++){
		sub = abs(middle -i);
		sub /= 1000;		
		if(f[i].a+b[i].a-sub> max){
			max = f[i].a+b[i].a-sub;
	//		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
			transition = 1;
			c = i;
		}
		if(f[i].a+b[i].ga-open-sub > max){
			max = f[i].a+b[i].ga-open-sub;
	//		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
			transition = 2;
			c = i;
		}
		if(f[i].a+b[i].gb+prof1[27]-sub > max){
			max = f[i].a+b[i].gb+prof1[27]-sub;
	//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
			transition = 3;
			c = i;
		}
		if(f[i].ga+b[i].a-open-sub > max){
			max = f[i].ga+b[i].a-open-sub;
	//		fprintf(stderr,"gap_a->aligned:%d + %d + %d(local_gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
			transition = 5;
			c = i;
		}


		if(hm->startb == 0){
			if(f[i].gb+b[i].gb+prof1[29]-sub > max){
				max = f[i].gb+b[i].gb+prof1[29]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}else{
			if(f[i].gb+b[i].gb+prof1[28]-sub > max){
				max = f[i].gb+b[i].gb+prof1[28]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}
		if(f[i].gb+b[i].a+prof1[27]-sub > max){
			max = f[i].gb+b[i].a+prof1[27]-sub;
	//		fprintf(stderr,"gap_b->aligned:%d + %d + %d(local_gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
			transition = 7;
			c = i;
		}
	}
	//i = hm->endb;
	i = old_cor[3];
	
	sub = abs(middle -i);
	sub /= 1000; 
	if(f[i].a+b[i].gb+prof1[27]-sub > max){
		max = f[i].a+b[i].gb+prof1[27]-sub;
	//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
		transition = 3;
		c = i;
	}
	if(hm->endb == hm->len_b){
		if(f[i].gb+b[i].gb+prof1[29]-sub > max){
			max = f[i].gb+b[i].gb+prof1[29]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}	
	}else{
		if(f[i].gb+b[i].gb+prof1[28]-sub > max){
			max = f[i].gb+b[i].gb+prof1[28]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}
	}
	
	
	
	prof1-= ((old_cor[4]+1)<<6);
	
	//fprintf(stderr,"Transition:%d	at:%d\n",transition,c);
	
	j = hirsch_path[0];
	switch(transition){
		case 1: //a -> a = 1
			
			hirsch_path[old_cor[4]] = c;
			hirsch_path[old_cor[4]+1] = c+1;
			
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
	//		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
			break;
		case 2:// a -> ga = 2
			
			hirsch_path[old_cor[4]] = c;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

			//backward:
			hm->starta = old_cor[4];
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = 0;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
			break;
		case 3:// a -> gb = 3
			
			hirsch_path[old_cor[4]] = c;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_large_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
			break;
		case 5://ga -> a = 5
			hirsch_path[old_cor[4]+1] = c+1;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = 0;
			hm->b[0].gb = -FLOATINFTY;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4];
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_large_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
			break;
		case 6://gb->gb = 6;
			
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_ps_dyn(prof1,seq2,hm,hirsch_path,sip);			


			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_large_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
			break;
		case 7://gb->a = 7;
			
			hirsch_path[old_cor[4]+1] = c+1;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_large_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
			break;
	}
		
	return hirsch_path;
}

struct large_states* foward_large_hirsch_ps_dyn(const float* prof1,const int* seq2,struct hirsch_large_mem* hm,int sip)
{
//	unsigned int freq[26];
	struct large_states* s = hm->f;
	const int starta = hm->starta;
	const int enda = hm->enda;
	const int startb = hm->startb;
	const int endb = hm->endb;
	
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	register int i = 0;
	register int j = 0;
	
	const float open = local_gpo * sip;
	const float ext = local_gpe *sip; 
	const float text = local_tgpe * sip;
	
	
	
	prof1 += (starta)<< 6;
	s[startb].a = s[0].a;
	s[startb].ga = s[0].ga;
	s[startb].gb = s[0].gb;
	if(startb == 0){
		for (j = startb+1; j < endb;j++){
			s[j].a = -FLOATINFTY;
			
			//s[j].ga = s[j-1].a-text;
			//if (s[j-1].ga-text > s[j].ga){
			//	s[j].ga = s[j-1].ga-text;
			//}
			//if(s[j-1].ga > s[j-1].a){
			//	s[j].ga = s[j-1].ga-text;
			//}else{
			//	s[j].ga = s[j-1].a-text;
			//}
			s[j].ga = MAX(s[j-1].ga,s[j-1].a) - text;
			s[j].gb = -FLOATINFTY;
		}	
	}else{

		for (j = startb+1; j < endb;j++){
			s[j].a = -FLOATINFTY;
			
			//s[j].ga = s[j-1].a-open;
			//if (s[j-1].ga-ext > s[j].ga){
			//	s[j].ga = s[j-1].ga-ext;
			//}
			//if(s[j-1].ga-ext > s[j-1].a-open){
			//	s[j].ga = s[j-1].ga-ext;
			//}else{
			//	s[j].ga = s[j-1].a-open;
			//}	
			s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);
			
			s[j].gb = -FLOATINFTY;
		}
	}
	
	
	s[endb].a = -FLOATINFTY;
	s[endb].ga = -FLOATINFTY;
	s[endb].gb = -FLOATINFTY;
	seq2--;

	for (i = starta;i < enda;i++){
		prof1 += 64;
		//pa = 1;
		//for (j = 26; j--;){
		//	if(prof1[j]){
		//		freq[pa] = j;
		//		pa++;	
		//	}
		//}
		//freq[0] = pa;
		pa = s[startb].a;
		pga = s[startb].ga;
		pgb = s[startb].gb;
		if(startb == 0){
			s[startb].a = -FLOATINFTY;
			s[startb].ga = -FLOATINFTY;
		
			//s[startb].gb = pa+prof1[29];
			//if(pgb+prof1[29] > s[startb].gb){
			//	s[startb].gb = pgb+prof1[29];
			//}
			//if(pgb > pa){
			//	s[startb].gb = pgb+prof1[29];
			//}else{
			//	s[startb].gb = pa+prof1[29];
			//}
			s[startb].gb = MAX(pgb,pa)+prof1[29];
		}else{
			s[startb].a = -FLOATINFTY;
			s[startb].ga = -FLOATINFTY;
		
			//s[startb].gb = pa+prof1[27];
			//if(pgb+prof1[28] > s[startb].gb){
			//	s[startb].gb = pgb+prof1[28];
			//}
			//if(pgb+prof1[28] > pa+prof1[27]){
			//	s[startb].gb = pgb+prof1[28];
			//}else{
			//	s[startb].gb = pa+prof1[27];
			//}
			s[startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
		}
		for (j = startb+1; j <= endb;j++){
			ca = s[j].a;
			
			//if((pga -= open) > pa){
			//	pa = pga;
			//}

			//if((pgb += prof1[-37]) > pa){
			//	pa = pgb;
			//}
			pa = MAX3(pa,pga -open,pgb + prof1[-37]);
			
			pa += prof1[32 + seq2[j]];


			s[j].a = pa;
			
			pga = s[j].ga;
			
			//s[j].ga = s[j-1].a-open;
			//if (s[j-1].ga-ext > s[j].ga){
			//	s[j].ga = s[j-1].ga-ext;
			//}
			
			//if (s[j-1].ga-ext  > s[j-1].a-open){
			//	s[j].ga = s[j-1].ga-ext;
			//}else{
			//	s[j].ga = s[j-1].a-open;
			//}
			s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);
			
			
			pgb = s[j].gb;
			
			//s[j].gb = ca+prof1[27];
			//if(pgb+prof1[28] > s[j].gb){
			//	s[j].gb = pgb+prof1[28];
			//}
			//if(pgb+prof1[28] > ca+prof1[27]){
			//	s[j].gb = pgb+prof1[28];
			//}else{
			//	s[j].gb = ca+prof1[27];
			//}
			s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);
			
			
			pa = ca;
		}	
	}
	prof1 -= enda << 6;
	return s;
}

struct large_states* backward_large_hirsch_ps_dyn(const float* prof1,const int* seq2,struct hirsch_large_mem* hm,int sip)
{
	//unsigned int freq[26];
	struct large_states* s = hm->b;
	const int starta = hm->starta;
	const int enda = hm->enda;
	const int startb = hm->startb;
	const int endb = hm->endb;
	
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	register int i = 0;
	register int j = 0;
	
	const float open = local_gpo * sip;
	const float ext = local_gpe *sip; 
	const float text = local_tgpe * sip;
	

	prof1 += (enda+1) << 6;

	s[endb].a = s[0].a;
	s[endb].ga = s[0].ga;
	s[endb].gb = s[0].gb;
	
	
	//init of first row;
	//j = endb-startb;
	if(endb == hm->len_b){
		for(j = endb-1;j > startb;j--){
			s[j].a = -FLOATINFTY;
			//s[j].ga = s[j+1].a-text;
			//if (s[j+1].ga-text > s[j].ga){
			//	s[j].ga = s[j+1].ga-text;
			//}
			//if(s[j+1].ga > s[j+1].a){
			//	s[j].ga = s[j+1].ga-text;
			//}else{
			//	s[j].ga = s[j+1].a-text;
			//}
			s[j].ga = MAX(s[j+1].ga,s[j+1].a)-text;
			
			s[j].gb = -FLOATINFTY;
		}
	}else{
		for(j = endb-1;j > startb;j--){
			s[j].a = -FLOATINFTY;
			
			//s[j].ga = s[j+1].a-open;
			//if (s[j+1].ga-ext > s[j].ga){
			//	s[j].ga = s[j+1].ga-ext;
			//}
			//if(s[j+1].ga-ext > s[j+1].a-open){
			//	s[j].ga = s[j+1].ga-ext;
			//}else{
			//	s[j].ga = s[j+1].a-open;
			//}
			s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
			s[j].gb = -FLOATINFTY;
		}
	}
	
	s[startb].a = -FLOATINFTY;
	s[startb].ga = -FLOATINFTY;
	s[startb].gb = -FLOATINFTY;

	i = enda-starta;
	while(i--){
		prof1 -= 64;

		//pa = 1;
		//for (j = 26; j--;){
		//	if(prof1[j]){
		//		freq[pa] = j;
		//		pa++;	
		//	}
		//}
		//freq[0] = pa;
		
		pa = s[endb].a;
		pga = s[endb].ga;
		pgb = s[endb].gb;
		s[endb].a = -FLOATINFTY;
		s[endb].ga = -FLOATINFTY;

		if(endb == hm->len_b){
			//s[endb].gb = pa+prof1[29];
			//if(pgb+prof1[29] > s[endb].gb){
			//	s[endb].gb = pgb+prof1[29];
			//}
			//if(pgb > pa){
			//	s[endb].gb = pgb+prof1[29];
			//}else{
			//	s[endb].gb = pa+prof1[29];
			//}
			s[endb].gb = MAX(pgb,pa) +prof1[29];
		}else{
			//s[endb].gb = pa+prof1[27];
			//if(pgb+prof1[28] > s[endb].gb){
			//	s[endb].gb = pgb+prof1[28];
			//}
			//if(pgb+prof1[28] > pa+prof1[27]){
			//	s[endb].gb = pgb+prof1[28];
			//}else{
			//	s[endb].gb = pa+prof1[27];
			//}
			s[endb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
			
			
		}

		for(j = endb-1;j >= startb;j--){
			ca = s[j].a;
			//if((pga -= open) > pa){
			//	pa = pga;
			//}
			//if((pgb += prof1[91]) > pa){
			//	pa = pgb;
			//}

			pa = MAX3(pa,pga - open,pgb +prof1[91]);
			pa += prof1[32 + seq2[j]];

			s[j].a = pa;
			
			pga = s[j].ga;
			
			//s[j].ga = s[j+1].a-open;
			//if (s[j+1].ga-ext > s[j].ga){
			//	s[j].ga = s[j+1].ga-ext;
			//}
			//if (s[j+1].ga-ext  > s[j+1].a-open){
			//	s[j].ga = s[j+1].ga-ext;
			//}else{
			//	s[j].ga = s[j+1].a-open;
			//}
			s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
			
			
			pgb = s[j].gb;
			
			//s[j].gb = ca+prof1[27];
			//if(pgb+prof1[28] > s[j].gb){
			//	s[j].gb = pgb+prof1[28];
			//}
			//if(pgb+prof1[28] > ca+prof1[27]){
			//	s[j].gb = pgb+prof1[28];
			//}else{
			//	s[j].gb = ca+prof1[27];
			//}
			s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);
			
			pa = ca;
		}
	}		
	return s;
}




int* hirsch_large_pp_dyn(const float* prof1,const float* prof2,struct hirsch_large_mem* hm, int* hirsch_path)
{
	int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
	float input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
	int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

	
	//fprintf(stderr,"starta:%d enda:%d startb:%d endb:%d mid:%d\n",hm->starta,hm->enda,hm->startb,hm->endb,mid);
	
	
	if(hm->starta  >= hm->enda){
		return hirsch_path;
	}
	if(hm->startb  >= hm->endb){
		return hirsch_path;
	}

	hm->enda = mid;
	hm->f = foward_large_hirsch_pp_dyn(prof1,prof2,hm);
	/*int i;
	fprintf(stderr,"FOWARD\n");
	for (i = hm->startb; i <= hm->endb;i++){
		fprintf(stderr,"%d	%d	%d\n",hm->f[i].a,hm->f[i].ga,hm->f[i].gb);
	}*/

	hm->starta = mid;
	hm->enda = old_cor[1];
	hm->b = backward_large_hirsch_pp_dyn(prof1,prof2,hm);
	/*fprintf(stderr,"BaCKWARD\n");

	for (i = hm->startb; i <= hm->endb;i++){
		fprintf(stderr,"%d	%d	%d\n",hm->b[i].a,hm->b[i].ga,hm->b[i].gb);
	}*/

	hirsch_path = hirsch_large_align_two_pp_vector(prof1,prof2,hm,hirsch_path,input_states,old_cor);
	return hirsch_path;
}



int* hirsch_large_align_two_pp_vector(const float* prof1,const float* prof2,struct hirsch_large_mem* hm,int* hirsch_path,float input_states[],int old_cor[])
{
	struct large_states* f = hm->f;
 	struct large_states* b = hm->b;
	int i,j,c;
	int transition = -1;
	
	
	//code:
	// a -> a = 1
	// a -> ga = 2
	// a -> gb = 3
	// ga ->ga = 4
	// ga -> a = 5
	//gb->gb = 6;
	//gb->a = 7;
	
	//int max = -FLOATINFTY;
	float max = -FLOATINFTY;	
	//float middle =  (hm->endb - hm->startb)/2 + hm->startb;
	float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
	float sub = 0.0;
	

	prof1+= ((old_cor[4]+1) << 6);
	//prof2 += 64 * (hm->startb);
	//i = hm->startb;
	prof2 += old_cor[2] << 6;
	i = old_cor[2];
	c = -1;
	//for(i = hm->startb; i < hm->endb;i++){
	for(i = old_cor[2]; i < old_cor[3];i++){
		sub = abs(middle -i);
		sub /= 1000; 
		prof2 += 64;
		//fprintf(stderr,"%d	%d	%d \n",f[i].a,b[i].a,max);
		if(f[i].a+b[i].a-sub > max){
			max = f[i].a+b[i].a-sub;
	//		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
			transition = 1;
			c = i;
		}
		if(f[i].a+b[i].ga+prof2[27]-sub > max){
			max = f[i].a+b[i].ga+prof2[27]-sub;
	//		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
			transition = 2;
			c = i;
		}
		if(f[i].a+b[i].gb+prof1[27] -sub> max){
			max = f[i].a+b[i].gb+prof1[27]-sub;
	//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
			transition = 3;
			c = i;
		}
		if(f[i].ga+b[i].a+prof2[27]-sub > max){
			max = f[i].ga+b[i].a+prof2[27]-sub;
	//		fprintf(stderr,"gap_a->aligned:%d + %d + %d(local_gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
			transition = 5;
			c = i;
		}


		if(hm->startb == 0){
			if(f[i].gb+b[i].gb+prof1[29]-sub > max){
				max = f[i].gb+b[i].gb+prof1[29]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}else{
			if(f[i].gb+b[i].gb+prof1[28]-sub > max){
				max = f[i].gb+b[i].gb+prof1[28]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}
		if(f[i].gb+b[i].a+prof1[27]-sub > max){
			max = f[i].gb+b[i].a+prof1[27]-sub;
	//		fprintf(stderr,"gap_b->aligned:%d + %d + %d(local_gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
			transition = 7;
			c = i;
		}
	}
	//i = hm->endb;
	i = old_cor[3];
	sub = abs(middle -i);
	sub /= 1000; 
	if(f[i].a+b[i].gb+prof1[27]-sub > max){
		max = f[i].a+b[i].gb+prof1[27]-sub;
	//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
		transition = 3;
		c = i;
	}
	if(hm->endb == hm->len_b){
		if(f[i].gb+b[i].gb+prof1[29]-sub > max){
			max = f[i].gb+b[i].gb+prof1[29]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}	
	}else{
		if(f[i].gb+b[i].gb+prof1[28]-sub > max){
			max = f[i].gb+b[i].gb+prof1[28]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}
	}
	
	
	
	prof1-= (old_cor[4]+1)<<6;
	//prof2 -= hm->endb << 6;
	prof2 -= old_cor[3] << 6;
	
	//fprintf(stderr,"Transition:%d	at:%d\n",transition,c);
	//if(transition == -1){
	//	exit(0);
	//}
	
	j = hirsch_path[0];
	switch(transition){
		case 1: //a -> a = 1
			
			hirsch_path[old_cor[4]] = c;
			hirsch_path[old_cor[4]+1] = c+1;
			
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			//fprintf(stderr,"Using this for start:%ld	%ld	%ld\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
		case 2:// a -> ga = 2
			
			hirsch_path[old_cor[4]] = c;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4];
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = 0;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
		case 3:// a -> gb = 3
			
			hirsch_path[old_cor[4]] = c;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_large_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
		case 5://ga -> a = 5
			hirsch_path[old_cor[4]+1] = c+1;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = 0;
			hm->b[0].gb = -FLOATINFTY;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4];
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_large_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
		case 6://gb->gb = 6;
			
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_large_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
		case 7://gb->a = 7;
			
			hirsch_path[old_cor[4]+1] = c+1;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_large_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_large_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
	}
		
	return hirsch_path;
}

struct large_states* foward_large_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_large_mem* hm)
{
	int freq[26];

	/*union print{
		__m128i m;
		int tmp[4];
	} output;
	
	__m128i xmm1;
	__m128i xmm2;*/

	
	struct large_states* s = hm->f;
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	register int i = 0;
	register int j = 0;
	register int c = 0;
	
	prof1 += (hm->starta) << 6;
	prof2 +=  (hm->startb) << 6;
	s[hm->startb].a = s[0].a;
	s[hm->startb].ga = s[0].ga;
	s[hm->startb].gb = s[0].gb;
	/*if(s[hm->startb].ga == -FLOATINFTY && s[hm->startb].a == -FLOATINFTY){
		for (j = hm->startb+1; j < hm->endb;j++){
			prof2+=64;
			s[j].a = -FLOATINFTY;
			s[j].ga = -FLOATINFTY;
			s[j].gb = -FLOATINFTY;
		}	
		prof2+=64;	
	}else{
	*/
	if(hm->startb == 0){
		for (j = hm->startb+1; j < hm->endb;j++){
			prof2+=64;
			s[j].a = -FLOATINFTY;
			//if (s[j-1].ga > s[j-1].a){
			//	s[j].ga = s[j-1].ga+prof2[29];
			//}else{
			//	s[j].ga = s[j-1].a+prof2[29];
			//}
			s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof2[29];
			s[j].gb = -FLOATINFTY;
		}	
		prof2+=64;	
	}else{

		for (j = hm->startb+1; j < hm->endb;j++){
			prof2+=64;
			s[j].a = -FLOATINFTY;
			//if(s[j-1].ga+prof2[28] > s[j-1].a+prof2[27]){
			//	s[j].ga = s[j-1].ga+prof2[28];
			//}else{
			//	s[j].ga = s[j-1].a+prof2[27];
			//}
			s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);
			
			
			s[j].gb = -FLOATINFTY;
		}
		prof2+=64;
	}
	//}
	
	prof2 -= (hm->endb-hm->startb) << 6;
	
	s[hm->endb].a = -FLOATINFTY;
	s[hm->endb].ga = -FLOATINFTY;
	s[hm->endb].gb = -FLOATINFTY;


	for (i = hm->starta;i < hm->enda;i++){
		prof1 += 64;
		c = 1;
		for (j = 0;j < 26; j++){
			if(prof1[j]){
				freq[c] = j;
				c++;	
			}
		}
		freq[0] = c;
			
		pa = s[hm->startb].a;
		pga = s[hm->startb].ga;
		pgb = s[hm->startb].gb;
		s[hm->startb].a = -FLOATINFTY;
		s[hm->startb].ga = -FLOATINFTY;
		
		//if(pgb == -FLOATINFTY && pa == -FLOATINFTY){
		//	s[hm->startb].gb = -FLOATINFTY;
		//}else{
		if(hm->startb == 0){
			//if(pgb > pa ){
			//	s[hm->startb].gb = pgb+prof1[29];
			//}else{
			//	s[hm->startb].gb = pa+prof1[29];
			//}
			s[hm->startb].gb = MAX(pgb,pa)+ prof1[29];
		}else{
			//if(pgb+prof1[28] >  pa+prof1[27]){
			//	s[hm->startb].gb = pgb+prof1[28];
			//}else{
			//	s[hm->startb].gb = pa+prof1[27];
			//}
			s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
		}
		//}
		for (j = hm->startb+1; j <= hm->endb;j++){
			prof2 += 64;
			ca = s[j].a;
			

			//if((pga += prof2[-37]) > pa){
			//	pa = pga;
			//}
		

			//if((pgb += prof1[-37]) > pa){
			//	pa = pgb;
			//}
			//pa = MAX(pgb + prof1[-37],pa);
			//pa = MAX(pga + prof2[-37],pa);

			pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);
			
			
			prof2 += 32;
			for (c = 1;c < freq[0];c++){
				pa += prof1[freq[c]]*prof2[freq[c]];
			}
			prof2 -= 32;

			s[j].a = pa;
			
			pga = s[j].ga;
			
			//if (s[j-1].ga+prof2[28] > s[j-1].a+prof2[27]){
			//	s[j].ga = s[j-1].ga+prof2[28];
			//}else{
			//	s[j].ga = s[j-1].a+prof2[27];
			//}
			
			
			/*xmm1 =  _mm_set_epi32 (ca,s[j-1].a,pgb,s[j-1].ga);
			xmm2 = _mm_set_epi32 (prof1[27],prof2[27],prof1[28],prof2[28]);
			xmm1 = _mm_add_epi32 (xmm1,xmm2);
			xmm2 =  _mm_srli_si128(xmm1, 8);
			output.m = _mm_cmpgt_epi32(xmm1,xmm2);
			output.m = _mm_or_si128( _mm_andnot_si128(output.m,xmm2),_mm_and_si128(output.m,xmm1));
			
			s[j].ga =output.tmp[0];
			s[j].gb = output.tmp[1];*/
			//output.m = _mm_add_epi32 (xmm1,xmm2);
			//_mm_store_si128(dst_ptr, xmm3);
			//fprintf(stderr,"%d	%d	%d	%d	%d	%d	%d	%d\n",output.tmp[0],output.tmp[1],output.tmp[2],output.tmp[3],s[j-1].ga+prof2[28],s[j-1].a+prof2[27],pgb+prof1[28] ,ca+prof1[27]);
			
			
			s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);
			
			
			
			pgb = s[j].gb;
			

			//if(pgb+prof1[28] > ca+prof1[27]){
			//	s[j].gb = pgb+prof1[28];
			//}else{
			//	s[j].gb = ca+prof1[27];
			//}
			s[j].gb = MAX(pgb+prof1[28] ,ca+prof1[27]);
			//fprintf(stderr,"%d	%d	%d	%d\n",output.tmp[0],output.tmp[1],s[j].ga,s[j].gb );
			pa = ca;
		}
		prof2 -= (hm->endb-hm->startb) << 6;
		
	}
	prof1 -=  (hm->enda) << 6;
	return s;
}

struct large_states* backward_large_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_large_mem* hm)
{
	int freq[26];
	struct large_states* s = hm->b;
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	register int i = 0;
	register int j = 0;
	register int c = 0;
	
	prof1 += (hm->enda+1) << 6;
	prof2 += (hm->endb+1) << 6;
	s[hm->endb].a = s[0].a;
	s[hm->endb].ga = s[0].ga;
	s[hm->endb].gb = s[0].gb;
	
	
	//init of first row;
	//j = endb-startb;

	/*if(s[hm->endb].ga == -FLOATINFTY && s[hm->endb].a == -FLOATINFTY){
		for(j = hm->endb-1;j > hm->startb;j--){
			prof2 -= 64;
			s[j].a = -FLOATINFTY;
			s[j].ga = -FLOATINFTY;
			s[j].gb = -FLOATINFTY;
		}
		prof2 -= 64;
	}else{*/
	
	if(hm->endb == hm->len_b){
		
		for(j = hm->endb-1;j > hm->startb;j--){
			prof2 -= 64;
			s[j].a = -FLOATINFTY;
			//if(s[j+1].ga > s[j+1].a){
			//	s[j].ga = s[j+1].ga+prof2[29];
			//}else{
			//	s[j].ga = s[j+1].a+prof2[29];
			//}
			s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof2[29];
			s[j].gb = -FLOATINFTY;
		}
		prof2 -= 64;
	}else{
		for(j = hm->endb-1;j > hm->startb;j--){
			prof2 -= 64;
			s[j].a = -FLOATINFTY;
			//if(s[j+1].ga+prof2[28] > s[j+1].a+prof2[27]){
			//	s[j].ga = s[j+1].ga+prof2[28];
			//}else{
			//	s[j].ga = s[j+1].a+prof2[27];
			//}
			s[j].ga = MAX(s[j+1].ga+prof2[28],s[j+1].a+prof2[27]);
			s[j].gb = -FLOATINFTY;
		}
		prof2 -= 64;
	}
	//}
	
	s[hm->startb].a = -FLOATINFTY;
	s[hm->startb].ga = -FLOATINFTY;
	s[hm->startb].gb = -FLOATINFTY;
	
	i = hm->enda-hm->starta;
	while(i--){
		prof1 -= 64;

		c = 1;
		for (j = 0;j < 26; j++){
			if(prof1[j]){
				freq[c] = j;
				c++;	
			}
		}
		freq[0] = c;
		
		pa = s[hm->endb].a;
		pga = s[hm->endb].ga;
		pgb = s[hm->endb].gb;
		s[hm->endb].a = -FLOATINFTY;
		s[hm->endb].ga = -FLOATINFTY;
		//if(pgb == -FLOATINFTY && pa == -FLOATINFTY){
		//	s[hm->endb].gb = -FLOATINFTY;
		//}else{
		if(hm->endb == hm->len_b){
			//if(pgb > pa){
			//	s[hm->endb].gb = pgb+prof1[29];
			//}else{
			//	s[hm->endb].gb = pa+prof1[29];
			//}	
			s[hm->endb].gb = MAX(pgb,pa)+prof1[29];
		}else{
			//if(pgb+prof1[28] > pa+prof1[27]){
			//	s[hm->endb].gb = pgb+prof1[28];
			//}else{
			//	s[hm->endb].gb = pa+prof1[27];
			//}
			s[hm->endb].gb = MAX(pgb+prof1[28] ,pa+prof1[27]);
		}
		//}
		//j = endb-startb;
		prof2 += (hm->endb-hm->startb) << 6;
		//while(j--){
		for(j = hm->endb-1;j >= hm->startb;j--){
			prof2 -= 64;
			ca = s[j].a;
			
			//pa = MAX(pga + prof2[91],pa);
			//pa = MAX(pgb + prof1[91],pa);
			pa = MAX3(pa,pga + prof2[91],pgb + prof1[91]);
			
			//if((pga += prof2[91]) > pa){
			//	pa = pga;
			//}
			//if((pgb += prof1[91]) > pa){
			//	pa = pgb;
			//}

			prof2 += 32;
			for (c = 1;c < freq[0];c++){
				pa += prof1[freq[c]]*prof2[freq[c]];
			}
			prof2 -= 32;

			s[j].a = pa;
			
			pga = s[j].ga;

			//if (s[j+1].ga+prof2[28] > s[j+1].a+prof2[27]){
			//	s[j].ga = s[j+1].ga+prof2[28];
			//}else{
			//	s[j].ga = s[j+1].a+prof2[27];
			//}
			s[j].ga = MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);

			pgb = s[j].gb;

			s[j].gb = MAX(pgb+prof1[28], ca+prof1[27]);
			//if(pgb+prof1[28] > ca+prof1[27]){
			//	s[j].gb = pgb+prof1[28];
			//}else{
			//	s[j].gb = ca+prof1[27];
			//}

			pa = ca;
		}
	}		
	return s;
}

float* make_large_profile(float* prof, int* seq,int len,float** subm)
{
	int i,j,c;	
	prof = malloc(sizeof(float)*(len+2)*64);
	prof +=  (64 *(len+1));

	for (i = 0;i < 64;i++){
		prof[i] = 0.0;
	}
	prof[23+32] = -local_gpo;
	prof[24+32] = -local_gpe;
	prof[25+32] = -local_tgpe;

	
	i = len;
	while(i--){
		prof -= 64;

		for (j = 0;j < 64;j++){
			prof[j] = 0.0;
		}
		c = seq[i];
		
		prof[c] += 1.0;
		
		prof += 32;
		
		for(j = 23;j--;){
			prof[j] = subm[c][j];
		}
		prof[23] = -local_gpo;
		prof[24] = -local_gpe;
		prof[25] = -local_tgpe;
		
		prof -= 32;
	}
	prof -= 64;
	for (i = 0;i < 64;i++){
		prof[i] = 0.0;
	}
	prof[23+32] = -local_gpo;
	prof[24+32] = -local_gpe;
	prof[25+32] = -local_tgpe;	
	return prof;
}

void set_large_gap_penalties(float* prof,int len,int nsip)
{
	int i;
	
	prof +=  (64 *(len+1));
	prof[27] = prof[55]*nsip;//gap open or close
	prof[28] = prof[56]*nsip;//gap extention
		
	prof[29] = prof[57]*nsip;//gap open or close
	i = len+1;
	while(i--){
		prof -= 64;
		prof[27] = prof[55]*nsip;//gap open or close
		prof[28] = prof[56]*nsip;//gap extention
		
		prof[29] = prof[57]*nsip;//gap open or close
	}
}


float* large_update(float* profa,float* profb,float* newp,int* path,int sipa,int sipb)
{
	int i,j,c;
	for (i = 64; i--;){
		newp[i] = profa[i] + profb[i];
	}
	
	profa += 64;
	profb += 64;
	newp += 64;

	c = 1;
	
	while(path[c] != 3){
		//Idea: limit the 'virtual' number of residues of one type to x.
		// i.e. only allow a maximum of 10 alanines to be registered in each column
		// the penalty for aligning a 'G' to this column will stay stable even when many (>10) alanines are present.
		// the difference in score between the 'correct' (all alanine) and incorrect (alanines + glycine) will not increase 
		// with the number of sequences. -> see Durbin pp 140
		
		if (!path[c]){
			//fprintf(stderr,"Align	%d\n",c);
			for (i = 64; i--;){
				newp[i] = profa[i] + profb[i];
			}
				
			
			profa += 64;
			profb += 64;
		}
		
		if (path[c] & 1){
			//fprintf(stderr,"Gap_A:%d\n",c);
			//printf("open:%d	ext:%d	%d	%d\n",si->nsip[a] * local_gpo,si->nsip[a] * gpe,si->nsip[a] * profb[41],si->nsip[a] * profb[46]);
			for (i = 64; i--;){
				newp[i] = profb[i];
			}
			profb += 64;
			if(!(path[c] & 20)){
				if(path[c] & 32){
					newp[25] += sipa;//1;
					i = local_tgpe*sipa;
				}else{
					newp[24] += sipa;//1;
					i = local_gpe*sipa;
				}
				
				for (j = 32; j < 55;j++){
					newp[j] -=i;
				}
			}else{
			if (path[c] & 16){ 
	//			fprintf(stderr,"close_open");
				if(path[c] & 32){
					newp[25] += sipa;//1;
					i = local_tgpe*sipa;
					newp[23] += sipa;//1;
					i += local_gpo*sipa;
				}else{
					newp[23] += sipa;//1;
					i = local_gpo*sipa;
				}
								
				for (j = 32; j < 55;j++){
					newp[j] -=i;
				}
			}
			if (path[c] & 4){ 
	//			fprintf(stderr,"Gap_open");
				if(path[c] & 32){
					newp[25] += sipa;//1;
					i = local_tgpe*sipa;
					newp[23] += sipa;//1;
					i += local_gpo*sipa;
				}else{
					newp[23] += sipa;//1;
					i = local_gpo*sipa;
				}
				for (j = 32; j < 55;j++){
					newp[j] -=i;
				}
			}
			}	
		}
		if (path[c] & 2){
			//fprintf(stderr,"Gap_B:%d\n",c);
			//printf("open:%d	ext:%d	%d	%d\n",si->nsip[b] * local_gpo,si->nsip[b] * gpe,profa[26],profa[27]);
			for (i = 64; i--;){
				newp[i] = profa[i];
			}
			profa+=64;
			if(!(path[c] & 20)){
				if(path[c] & 32){
					newp[25] += sipb;//1;
					i = local_tgpe*sipb;
				}else{
					newp[24] += sipb;//1;
					i = local_gpe*sipb;
				}
				for (j = 32; j < 55;j++){
					newp[j] -=i;
				}
			}else{
			if (path[c] & 16){
	//			fprintf(stderr,"close_open");
				if(path[c] & 32){
					newp[25] += sipb;//1;
					i =  local_tgpe*sipb;
					newp[23] += sipb;//1;
					i +=  local_gpo*sipb;
				}else{
					newp[23] += sipb;//1;
					i =  local_gpo*sipb;
				}
				for (j = 32; j < 55;j++){
					newp[j] -=i;
				}
			}
			if (path[c] & 4){
	//			fprintf(stderr,"Gap_open");
				if(path[c] & 32){
					newp[25] += sipb;//1;
					i = local_tgpe*sipb;
					newp[23] += sipb;//1;
					i += local_gpo*sipb;
				}else{
					newp[23] += sipb;//1;
					i = local_gpo*sipb;
				}
				
				for (j = 32; j < 55;j++){
					newp[j] -=i;
				}
			}
			}		
		}
		newp += 64;
		c++;
	}
	for (i = 64; i--;){
		newp[i] =  profa[i] + profb[i];
	}	
	newp -= (path[0]+1) *64;
	return newp;
}


struct hirsch_large_mem* hirsch_large_mem_alloc(struct hirsch_large_mem* hm,int x)
{

	// a=((typeof(a))(((int)(((void *)malloc(c+15))+15))&-16)). 
	hm = (struct hirsch_large_mem *) malloc(sizeof(struct hirsch_large_mem));
	hm->starta = 0;
	hm->startb = 0;
	hm->enda = 0;
	hm->endb = 0;
	hm->size = x;
	hm->len_a = 0;
	hm->len_b = 0;
	hm->f = malloc(sizeof(struct large_states)* (x+1));
	hm->b = malloc(sizeof(struct large_states)* (x+1));
	return hm;
}

struct hirsch_large_mem* hirsch_large_mem_realloc(struct hirsch_large_mem* hm,int x)
{
	hm->starta = 0;
	hm->startb = 0;
	hm->enda = 0;
	hm->endb = 0;
	hm->len_a = 0;
	hm->len_b = 0;
	hm->size = x;
	hm->f = realloc(hm->f,sizeof(struct large_states)* (x+1));
	hm->b = realloc(hm->b,sizeof(struct large_states)* (x+1));
	return hm;
}

void hirsch_large_mem_free(struct hirsch_large_mem* hm)
{
	free(hm->f);
	free(hm->b);
	free(hm);
}


