/*
	kalign2_hirschberg.c
	
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
#include "kalign2_hirschberg.h"
#define MAX(a, b) (a > b ? a : b)
#define MAX3(a,b,c) MAX(MAX(a,b),c)
//#include <emmintrin.h>

int** hirschberg_alignment(struct alignment* aln,int* tree,float**submatrix, int** map,int window,float strength)
{
	struct hirsch_mem* hm = 0;
	int i,j,g,a,b,c;
	int len_a;
	int len_b;
	float** profile = 0;

	profile = malloc(sizeof(float*)*numprofiles);
	for ( i = 0;i< numprofiles;i++){
		profile[i] = 0;
	}

	map = malloc(sizeof(int*)*numprofiles);
	for ( i = 0;i < numprofiles;i++){
		map[i] = 0;
	}
	
	hm = hirsch_mem_alloc(hm,1024);

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
			hm = hirsch_mem_realloc(hm,g);
		}

		for (j = 0; j < (g+2);j++){
			map[c][j] = -1;
		}

		if (a < numseq){
			profile[a] = make_profile(profile[a],aln->s[a],len_a,submatrix);
		}else{
			set_gap_penalties(profile[a],len_a,aln->nsip[b],strength,aln->nsip[a]);
			//smooth_gaps(profile[a],len_a,window,strength);
			
			//increase_gaps(profile[a],len_a,window,strength);
		}
		if (b < numseq){
			profile[b] = make_profile(profile[b],aln->s[b],len_b,submatrix);
		}else{		
			set_gap_penalties(profile[b],len_b,aln->nsip[a],strength,aln->nsip[b]);
			//smooth_gaps(profile[b],len_b,window,strength);
			//increase_gaps(profile[b],len_b,window,strength);
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
				map[c] = hirsch_ss_dyn(submatrix,aln->s[a],aln->s[b],hm,map[c]);
			}else{
				hm->enda = len_b;
				hm->endb = len_a;
				hm->len_a = len_b;
				hm->len_b = len_a;
				map[c] = hirsch_ps_dyn(profile[b],aln->s[a],hm,map[c],aln->nsip[b]);
				map[c] = mirror_hirsch_path(map[c],len_a,len_b);
			}
		}else{
			if(b < numseq){
				map[c] = hirsch_ps_dyn(profile[a],aln->s[b],hm,map[c],aln->nsip[a]);
			}else{
				if(len_a < len_b){
					map[c] = hirsch_pp_dyn(profile[a],profile[b],hm,map[c]);
				}else{
					hm->enda = len_b;
					hm->endb = len_a;
					hm->len_a = len_b;
					hm->len_b = len_a;
					map[c] = hirsch_pp_dyn(profile[b],profile[a],hm,map[c]);
					map[c] = mirror_hirsch_path(map[c],len_a,len_b);
				}
			}
		}
		
		map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);

		if(i != numseq-2){
			profile[c] = malloc(sizeof(float)*64*(map[c][0]+2));
			profile[c] = update(profile[a],profile[b],profile[c],map[c],aln->nsip[a],aln->nsip[b]);
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
	hirsch_mem_free(hm);
	for (i = 32;i--;){
		free(submatrix[i]);
	}
	free(submatrix);
	return map;
}


int** hirschberg_alignment_against_a(struct alignment* aln,int* tree,float**submatrix, int** map,int window,float strength)
{
	struct hirsch_mem* hm = 0;
	int i,j,g,a,b,c;
	int len_a;
	int len_b;
	float** profile = 0;

	profile = malloc(sizeof(float*)*numprofiles);
	for ( i = 0;i< numprofiles;i++){
		profile[i] = 0;
	}

	map = malloc(sizeof(int*)*numprofiles);
	for ( i = 0;i < numprofiles;i++){
		map[i] = 0;
	}
	
	hm = hirsch_mem_alloc(hm,1024);

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
			hm = hirsch_mem_realloc(hm,g);
		}

		for (j = 0; j < (g+2);j++){
			map[c][j] = -1;
		}

		if (a < numseq){
			profile[a] = make_profile(profile[a],aln->s[a],len_a,submatrix);
		}else{
			set_gap_penalties(profile[a],len_a,aln->nsip[b],0,aln->nsip[a]);
			//smooth_gaps(profile[a],len_a,window,strength);
			
			//increase_gaps(profile[a],len_a,window,strength);
		}
		if (b < numseq){
			profile[b] = make_profile(profile[b],aln->s[b],len_b,submatrix);
		}else{		
			set_gap_penalties(profile[b],len_b,aln->nsip[a],0,aln->nsip[b]);
			//smooth_gaps(profile[b],len_b,window,strength);
			//increase_gaps(profile[b],len_b,window,strength);
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
				map[c] = hirsch_ss_dyn(submatrix,aln->s[a],aln->s[b],hm,map[c]);
			}else{
				hm->enda = len_b;
				hm->endb = len_a;
				hm->len_a = len_b;
				hm->len_b = len_a;
				map[c] = hirsch_ps_dyn(profile[b],aln->s[a],hm,map[c],aln->nsip[b]);
				map[c] = mirror_hirsch_path(map[c],len_a,len_b);
			}
		}else{
			if(b < numseq){
				map[c] = hirsch_ps_dyn(profile[a],aln->s[b],hm,map[c],aln->nsip[a]);
			}else{
				if(len_a < len_b){
					map[c] = hirsch_pp_dyn(profile[a],profile[b],hm,map[c]);
				}else{
					hm->enda = len_b;
					hm->endb = len_a;
					hm->len_a = len_b;
					hm->len_b = len_a;
					map[c] = hirsch_pp_dyn(profile[b],profile[a],hm,map[c]);
					map[c] = mirror_hirsch_path(map[c],len_a,len_b);
				}
			}
		}
		
		map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);

		if(i != numseq-2){
			profile[c] = malloc(sizeof(float)*64*(map[c][0]+2));
			profile[c] = update_only_a(profile[a],profile[b],profile[c],map[c],aln->nsip[a],aln->nsip[b]);
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
	hirsch_mem_free(hm);
	for (i = 32;i--;){
		free(submatrix[i]);
	}
	free(submatrix);
	return map;
}


int* hirsch_ss_dyn(float**subm, const int* seq1,const int* seq2,struct hirsch_mem* hm, int* hirsch_path)
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
	hm->f = foward_hirsch_ss_dyn(subm,seq1,seq2,hm);

	hm->starta = mid;
	hm->enda = old_cor[1];
	//fprintf(stderr,"Backward:%d-%d	%d-%d\n",hm->starta,hm->enda,hm->startb,hm->endb);
	hm->b = backward_hirsch_ss_dyn(subm,seq1,seq2,hm);


	hirsch_path = hirsch_align_two_ss_vector(subm,seq1,seq2,hm,hirsch_path,input_states,old_cor);
	return hirsch_path;
}

int* hirsch_align_two_ss_vector(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[])
{
	struct states* f = hm->f;
 	struct states* b = hm->b;
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
	
	//int max = -INFTY;
	float max = -INFTY;
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
		if(f[i].a+b[i].ga-gpo-sub > max){
			max = f[i].a+b[i].ga-gpo-sub;
	//		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
			transition = 2;
			c = i;
		}
		if(f[i].a+b[i].gb -gpo-sub > max){
			max = f[i].a+b[i].gb - gpo-sub;
	//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
			transition = 3;
			c = i;
		}
		if(f[i].ga+b[i].a - gpo-sub > max){
			max = f[i].ga+b[i].a - gpo-sub;
	//		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
			transition = 5;
			c = i;
		}


		if(hm->startb == 0){
			if(f[i].gb+b[i].gb - tgpe-sub > max){
				max = f[i].gb+b[i].gb -tgpe-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}else{
			if(f[i].gb+b[i].gb - gpe -sub> max){
				max = f[i].gb+b[i].gb - gpe-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}
		if(f[i].gb+b[i].a - gpo-sub > max){
			max = f[i].gb+b[i].a - gpo-sub;
	//		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
			transition = 7;
			c = i;
		}
	}
	//i = hm->endb;
	i = old_cor[3];
	sub = abs(middle -i);
	sub /= 1000; 
	
	if(f[i].a+b[i].gb-gpo-sub > max){
		max = f[i].a+b[i].gb - gpo-sub;
	//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
		transition = 3;
		c = i;
	}
	if(hm->endb == hm->len_b){
		if(f[i].gb+b[i].gb -tgpe-sub > max){
			max = f[i].gb+b[i].gb - tgpe-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}	
	}else{
		if(f[i].gb+b[i].gb - gpe-sub > max){
			max = f[i].gb+b[i].gb - gpe-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
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
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
	//		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			break;
		case 2:// a -> ga = 2
			
			hirsch_path[old_cor[4]] = c;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			

			//backward:
			hm->starta = old_cor[4];
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = 0.0;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			break;
		case 3:// a -> gb = 3
			
			hirsch_path[old_cor[4]] = c;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0.0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			break;
		case 5://ga -> a = 5
			hirsch_path[old_cor[4]+1] = c+1;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = 0.0;
			hm->b[0].gb = -FLOATINFTY;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4];
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			break;
		case 6://gb->gb = 6;
			
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0.0;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0.0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
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
			hm->b[0].gb = 0.0;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			break;
	}
		
	return hirsch_path;
}



struct states* foward_hirsch_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm)
{
	struct states* s = hm->f;
	float *subp = 0;
	const int starta = hm->starta;
	const int enda = hm->enda;
	const int startb =hm->startb;
	const int endb = hm->endb;
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	register float xa = 0;
	register float xga = 0;
	register int i = 0;
	register int j = 0;
	

	s[startb].a = s[0].a;
	s[startb].ga = s[0].ga;
	s[startb].gb = s[0].gb;
	if(startb){
		for (j = startb+1; j < endb;j++){
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j-1].ga - gpe,s[j-1].a-gpo);
			s[j].gb = -FLOATINFTY;
		}
	}else{
		for (j = startb+1; j < endb;j++){
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j-1].ga,s[j-1].a)-tgpe;
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
		s[startb].a = -FLOATINFTY;
		s[startb].ga = -FLOATINFTY;
		
		xa = s[startb].a;
		xga = s[startb].ga;
		
		if(startb){
			s[startb].gb = MAX(pgb - gpe,pa - gpo);
		}else{
			s[startb].gb = MAX(pgb,pa) - tgpe;
		}
		for (j = startb+1; j < endb;j++){
			ca = s[j].a;
			pa = MAX3(pa,pga-gpo,pgb-gpo);
			pa += subp[seq2[j]];
			
			s[j].a = pa;
			
			pga = s[j].ga;
			//s[j].ga = MAX(s[j-1].ga-gpe,s[j-1].a-gpo);
			s[j].ga = MAX(xga-gpe,xa-gpo);
			
			pgb = s[j].gb;
			s[j].gb = MAX(pgb-gpe ,ca-gpo);
			
			pa = ca;
			
			xa = s[j].a;
			xga = s[j].ga;
			
		}
		ca = s[j].a;
		pa = MAX3(pa,pga-gpo,pgb-gpo);
		pa += subp[seq2[j]];
			
		s[j].a = pa;
			
		s[j].ga = -FLOATINFTY;//MAX(s[j-1].ga-gpe,s[j-1].a-gpo);
		if (endb != hm->len_b){
			s[j].gb = MAX(s[j].gb-gpe ,ca-gpo);
		}else{
			s[j].gb = MAX(s[j].gb,ca)-tgpe;
		}

	}
	return s;
}

struct states* backward_hirsch_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm)
{

	struct states* s = hm->b;
	float *subp = 0;
	const int starta = hm->starta;
	const int enda = hm->enda;
	const int startb =hm->startb;
	const int endb = hm->endb;
	
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	
	register float xa = 0;
	register float xga = 0;
	
	register int i = 0;
	register int j = 0;

	s[endb].a = s[0].a ;
	s[endb].ga = s[0].ga;
	s[endb].gb = s[0].gb;
	
	
	//init of first row;
	
	//j = endb-startb;
	if(endb != hm->len_b){
		for(j = endb-1;j > startb;j--){
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j+1].ga-gpe,s[j+1].a-gpo);	
			s[j].gb = -FLOATINFTY;
		}
	}else{
		for(j = endb-1;j > startb;j--){
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j+1].ga,s[j+1].a)-tgpe;
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

		xa = s[endb].a;
		xga = s[endb].ga;
		
		if(endb != hm->len_b){
			s[endb].gb = MAX(pgb-gpe,pa-gpo);
		}else{
			s[endb].gb = MAX(pgb,pa)-tgpe;
		}

		for(j = endb-1;j > startb;j--){

			ca = s[j].a;

			pa = MAX3(pa,pga - gpo,pgb-gpo);
			
			pa += subp[seq2[j]];

			s[j].a = pa;
			
			pga = s[j].ga;

			//s[j].ga = MAX(s[j+1].ga-gpe,s[j+1].a-gpo);
			
			s[j].ga = MAX(xga-gpe,xa-gpo);

			pgb = s[j].gb;
			s[j].gb = MAX(pgb-gpe,ca-gpo);
			
			pa = ca;
			xa = s[j].a;
			xga = s[j].ga;
		}
		ca = s[j].a;

		pa = MAX3(pa,pga - gpo,pgb-gpo);
			
		pa += subp[seq2[j]];

		s[j].a = pa;
		
		s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga-gpe,s[j+1].a-gpo);
		
		if(startb){
			s[j].gb = MAX(s[j].gb-gpe,ca-gpo);
		}else{
			s[j].gb = MAX(s[j].gb,ca)-tgpe;
		}

		
	}		
	return s;
}


int* hirsch_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm, int* hirsch_path,int sip)
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
	hm->f = foward_hirsch_ps_dyn(prof1,seq2,hm,sip);
	
	/*int i;
	fprintf(stderr,"FOWARD\n");
	for (i = hm->startb; i <= hm->endb;i++){
		fprintf(stderr,"%d	%d	%d\n",hm->f[i].a,hm->f[i].ga,hm->f[i].gb);
	}*/

	hm->starta = mid;
	hm->enda = old_cor[1];
	hm->b = backward_hirsch_ps_dyn(prof1,seq2,hm,sip);
	
	/*fprintf(stderr,"BaCKWARD\n");
	for (i = hm->startb; i <= hm->endb;i++){
		fprintf(stderr,"%d	%d	%d\n",hm->b[i].a,hm->b[i].ga,hm->b[i].gb);
	}*/

	hirsch_path = hirsch_align_two_ps_vector(prof1,seq2,hm,hirsch_path,input_states,old_cor,sip);
	return hirsch_path;
}



int* hirsch_align_two_ps_vector(const float* prof1,const int* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],int sip)
{
	struct states* f = hm->f;
 	struct states* b = hm->b;
	int i,j,c;
	int transition = -1;
	
	const float open = gpo * sip;
	
	
	//code:
	// a -> a = 1
	// a -> ga = 2
	// a -> gb = 3
	// ga ->ga = 4
	// ga -> a = 5
	//gb->gb = 6;
	//gb->a = 7;
	
	//int max = -INFTY;
	float max = -INFTY;
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
	//		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
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
		if(f[i].gb+b[i].a+prof1[-37]-sub > max){
			max = f[i].gb+b[i].a+prof1[-37]-sub;
	//		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
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
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
	//		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
			break;
		case 2:// a -> ga = 2
			
			hirsch_path[old_cor[4]] = c;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

			//backward:
			hm->starta = old_cor[4];
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = 0.0;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
			break;
		case 3:// a -> gb = 3
			
			hirsch_path[old_cor[4]] = c;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0.0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
			break;
		case 5://ga -> a = 5
			hirsch_path[old_cor[4]+1] = c+1;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = 0.0;
			hm->b[0].gb = -FLOATINFTY;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4];
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
			break;
		case 6://gb->gb = 6;
			
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0.0;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);			


			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0.0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
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
			hm->b[0].gb = 0.0;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
			break;
	}
		
	return hirsch_path;
}

struct states* foward_hirsch_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm,int sip)
{
	struct states* s = hm->f;

	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	
	register float xa = 0;
	register float xga = 0;
	
	register int i = 0;
	register int j = 0;
	
	const float open = gpo * sip;
	const float ext = gpe *sip;
	const float text = tgpe * sip;
	
	
	
	prof1 += (hm->starta)<< 6;
	s[hm->startb].a = s[0].a;
	s[hm->startb].ga = s[0].ga;
	s[hm->startb].gb = s[0].gb;
	if(hm->startb){
		for (j = hm->startb+1; j < hm->endb;j++){
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);	
			s[j].gb = -FLOATINFTY;
		}
	}else{
		for (j = hm->startb+1; j < hm->endb;j++){
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j-1].ga,s[j-1].a) - text;
			s[j].gb = -FLOATINFTY;
		}

	}
	
	
	s[hm->endb].a = -FLOATINFTY;
	s[hm->endb].ga = -FLOATINFTY;
	s[hm->endb].gb = -FLOATINFTY;
	seq2--;

	for (i = hm->starta;i < hm->enda;i++){
		prof1 += 64;

		pa = s[hm->startb].a;
		pga = s[hm->startb].ga;
		pgb = s[hm->startb].gb;
		s[hm->startb].a = -FLOATINFTY;
		s[hm->startb].ga = -FLOATINFTY;
		
		xa = s[hm->startb].a;
		xga = s[hm->startb].ga;
		
		
		if(hm->startb){
			s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
		}else{
			s[hm->startb].gb = MAX(pgb,pa)+prof1[29];
		}
		for (j = hm->startb+1; j < hm->endb;j++){
			ca = s[j].a;

			pa = MAX3(pa,pga -open,pgb + prof1[-37]);
			
			pa += prof1[32 + seq2[j]];


			s[j].a = pa;
			
			pga = s[j].ga;
			
			//s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);
			s[j].ga = MAX(xga-ext,xa-open);
			
				
			pgb = s[j].gb;
			
			s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);
				
			pa = ca;
			xa = s[j].a;
			xga = s[j].ga;
			
		}
		ca = s[j].a;

		pa = MAX3(pa,pga -open,pgb + prof1[-37]);
			
		pa += prof1[32 + seq2[j]];


		s[j].a = pa;

		s[j].ga = -FLOATINFTY;//MAX(s[j-1].ga-ext,s[j-1].a-open);
				
		if (hm->endb != hm->len_b){
			s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);
		}else{
			s[j].gb = MAX(s[j].gb,ca)+ prof1[29];
		}
		
	}
	prof1 -= hm->enda << 6;
	return s;
}

struct states* backward_hirsch_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm,int sip)
{
	struct states* s = hm->b;
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	
	register float xa = 0;
	register float xga = 0;
	
	register int i = 0;
	register int j = 0;
	
	const float open = gpo * sip;
	const float ext = gpe *sip;
	const float text = tgpe * sip;
	

	prof1 += (hm->enda+1) << 6;

	s[hm->endb].a = s[0].a;
	s[hm->endb].ga = s[0].ga;
	s[hm->endb].gb = s[0].gb;
	
	if(hm->endb != hm->len_b){
		for(j = hm->endb-1;j > hm->startb;j--){
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
			s[j].gb = -FLOATINFTY;
		}
	}else{
		for(j = hm->endb-1;j > hm->startb;j--){
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j+1].ga,s[j+1].a)-text;
			s[j].gb = -FLOATINFTY;
		}
	}
	
	s[hm->startb].a = -FLOATINFTY;
	s[hm->startb].ga = -FLOATINFTY;
	s[hm->startb].gb = -FLOATINFTY;

	i = hm->enda-hm->starta;
	while(i--){
		prof1 -= 64;
		pa = s[hm->endb].a;
		pga = s[hm->endb].ga;
		pgb = s[hm->endb].gb;
		s[hm->endb].a = -FLOATINFTY;
		s[hm->endb].ga = -FLOATINFTY;
		
		xa = s[hm->endb].a;
		xga = s[hm->endb].ga;
		

		if(hm->endb != hm->len_b){
			s[hm->endb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
		}else{
			s[hm->endb].gb = MAX(pgb,pa) +prof1[29];
		}

		for(j = hm->endb-1;j > hm->startb;j--){
			ca = s[j].a;

			pa = MAX3(pa,pga - open,pgb +prof1[91]);
			pa += prof1[32 + seq2[j]];

			s[j].a = pa;
			
			pga = s[j].ga;

			//s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
			s[j].ga = MAX(xga-ext,xa-open);
				
			pgb = s[j].gb;

			s[j].gb = MAX(pgb+prof1[28],ca+prof1[27]);
			
			pa = ca;
			xa = s[j].a;
			xga = s[j].ga;
			
			
		}
		ca = s[j].a;

		pa = MAX3(pa,pga - open,pgb +prof1[91]);
		pa += prof1[32 + seq2[j]];

		s[j].a = pa;
			

		s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga-ext,s[j+1].a-open);
		if(hm->startb){
			s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);
		}else{
			s[j].gb = MAX(s[j].gb,ca)+prof1[29];
		}

	}		
	return s;
}




int* hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm, int* hirsch_path)
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
	hm->f = foward_hirsch_pp_dyn(prof1,prof2,hm);
	/*int i;
	fprintf(stderr,"FOWARD\n");
	for (i = hm->startb; i <= hm->endb;i++){
		fprintf(stderr,"%d	%d	%d\n",hm->f[i].a,hm->f[i].ga,hm->f[i].gb);
	}*/

	hm->starta = mid;
	hm->enda = old_cor[1];
	hm->b = backward_hirsch_pp_dyn(prof1,prof2,hm);
	/*fprintf(stderr,"BaCKWARD\n");

	for (i = hm->startb; i <= hm->endb;i++){
		fprintf(stderr,"%d	%d	%d\n",hm->b[i].a,hm->b[i].ga,hm->b[i].gb);
	}*/

	hirsch_path = hirsch_align_two_pp_vector(prof1,prof2,hm,hirsch_path,input_states,old_cor);
	return hirsch_path;
}



int* hirsch_align_two_pp_vector(const float* prof1,const float* prof2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[])
{
	struct states* f = hm->f;
 	struct states* b = hm->b;
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
	
	//int max = -INFTY;
	float max = -INFTY;	
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
		if(f[i].ga+b[i].a+prof2[-37]-sub > max){
			max = f[i].ga+b[i].a+prof2[-37]-sub;
	//		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
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
		if(f[i].gb+b[i].a+prof1[-37]-sub > max){
			max = f[i].gb+b[i].a+prof1[-37]-sub;
	//		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
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
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			//fprintf(stderr,"Using this for start:%ld	%ld	%ld\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
		case 2:// a -> ga = 2
			
			hirsch_path[old_cor[4]] = c;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4];
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = 0.0;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
		case 3:// a -> gb = 3
			
			hirsch_path[old_cor[4]] = c;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0.0;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = -FLOATINFTY;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0.0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
		case 5://ga -> a = 5
			hirsch_path[old_cor[4]+1] = c+1;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = 0.0;
			hm->b[0].gb = -FLOATINFTY;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4];
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
		case 6://gb->gb = 6;
			
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -FLOATINFTY;
			hm->b[0].ga = -FLOATINFTY;
			hm->b[0].gb = 0.0;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -FLOATINFTY;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = 0.0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
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
			hm->b[0].gb = 0.0;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0.0;
			hm->f[0].ga = -FLOATINFTY;
			hm->f[0].gb = -FLOATINFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
	}
		
	return hirsch_path;
}

struct states* foward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm)
{
	unsigned int freq[26];

	struct states* s = hm->f;
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	
	register float xa = 0;
	register float xga = 0;
	
	register int i = 0;
	register int j = 0;
	register int c = 0;
	
	prof1 += (hm->starta) << 6;
	prof2 +=  (hm->startb) << 6;
	s[hm->startb].a = s[0].a;
	s[hm->startb].ga = s[0].ga;
	s[hm->startb].gb = s[0].gb;
	if(hm->startb){
		for (j = hm->startb+1; j < hm->endb;j++){
			prof2+=64;
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);
			s[j].gb = -FLOATINFTY;
		}
		prof2+=64;
	}else{
		for (j = hm->startb+1; j < hm->endb;j++){
			prof2+=64;
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof2[29];
			s[j].gb = -FLOATINFTY;
		}	
		prof2+=64;	
	}

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
		
		xa = s[hm->startb].a;
		xga = s[hm->startb].ga;
		

		if(hm->startb){
			s[hm->startb].gb = MAX(pgb+prof1[28],pa+prof1[27]);
		}else{
			s[hm->startb].gb = MAX(pgb,pa)+ prof1[29];
		}
		for (j = hm->startb+1; j < hm->endb;j++){
			prof2 += 64;
			ca = s[j].a;
			
			pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);
			
			prof2 += 32;
			for (c = 1;c < freq[0];c++){
				pa += prof1[freq[c]]*prof2[freq[c]];
			}
			prof2 -= 32;

			s[j].a = pa;
			
			pga = s[j].ga;
			
			//s[j].ga = MAX(s[j-1].ga+prof2[28],s[j-1].a+prof2[27]);
			s[j].ga = MAX(xga+prof2[28],xa+prof2[27]);
				
			pgb = s[j].gb;

			s[j].gb = MAX(pgb+prof1[28] ,ca+prof1[27]);

			pa = ca;
			
			
			xa = s[j].a;
			xga = s[j].ga;
		}
		prof2 += 64;
		ca = s[j].a;
			
		pa = MAX3(pa,pga + prof2[-37],pgb + prof1[-37]);

		prof2 += 32;
		for (c = 1;c < freq[0];c++){
			pa += prof1[freq[c]]*prof2[freq[c]];
		}
		prof2 -= 32;

		s[j].a = pa;

		s[j].ga = -FLOATINFTY;

		if (hm->endb != hm->len_b){
			s[j].gb = MAX(s[j].gb+prof1[28] ,ca+prof1[27]);
		}else{
			s[j].gb = MAX(s[j].gb,ca)+ prof1[29];
		}
		prof2 -= (hm->endb-hm->startb) << 6;
		
	}
	prof1 -=  (hm->enda) << 6;
	return s;
}

struct states* backward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm)
{
	unsigned int freq[26];
	struct states* s = hm->b;
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	
	register float xa = 0;
	register float xga = 0;
	
	register int i = 0;
	register int j = 0;
	register int c = 0;

	prof1 += (hm->enda+1) << 6;
	prof2 += (hm->endb+1) << 6;
	s[hm->endb].a = s[0].a;
	s[hm->endb].ga = s[0].ga;
	s[hm->endb].gb = s[0].gb;
	if(hm->endb != hm->len_b){
		for(j = hm->endb-1;j > hm->startb;j--){
			prof2 -= 64;
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j+1].ga+prof2[28],s[j+1].a+prof2[27]);
			s[j].gb = -FLOATINFTY;
		}
		prof2 -= 64;
	}else{
		for(j = hm->endb-1;j > hm->startb;j--){
			prof2 -= 64;
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof2[29];
			s[j].gb = -FLOATINFTY;
		}
		prof2 -= 64;
	}

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
		
		xa = s[hm->endb].a;
		xga = s[hm->endb].ga;
		
		if(hm->endb != hm->len_b){
			s[hm->endb].gb = MAX(pgb+prof1[28] ,pa+prof1[27]);
		}else{
			s[hm->endb].gb = MAX(pgb,pa)+prof1[29];
		}

		prof2 += (hm->endb-hm->startb) << 6;
		for(j = hm->endb-1;j > hm->startb;j--){
			prof2 -= 64;
			ca = s[j].a;

			pa = MAX3(pa,pga + prof2[91],pgb + prof1[91]);

			prof2 += 32;
			for (c = 1;c < freq[0];c++){
				pa += prof1[freq[c]]*prof2[freq[c]];
			}
			prof2 -= 32;

			s[j].a = pa;
			
			pga = s[j].ga;

			//s[j].ga = MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);
			s[j].ga = MAX(xga+prof2[28], xa+prof2[27]);

			pgb = s[j].gb;

			s[j].gb = MAX(pgb+prof1[28], ca+prof1[27]);

			pa = ca;
			xa = s[j].a;
			xga = s[j].ga;
		}
		prof2 -= 64;
		ca = s[j].a;

		pa = MAX3(pa,pga + prof2[91],pgb + prof1[91]);
		prof2 += 32;
		for (c = 1;c < freq[0];c++){
			pa += prof1[freq[c]]*prof2[freq[c]];
		}
		prof2 -= 32;
		s[j].a = pa;
		
		//pga = s[j].ga;
		s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);

		//pgb = s[j].gb;
		if(hm->startb){
			s[j].gb = MAX(s[j].gb+prof1[28], ca+prof1[27]);
		}else{
			s[j].gb = MAX(s[j].gb,ca)+prof1[29];
		}

		//pa = ca;
	}		
	return s;
}


int* mirror_hirsch_path(int* hirsch_path,int len_a,int len_b)
{
	int* np = 0;
	
	int i;
	np =malloc(sizeof(int)*(len_a+2));
	for(i =0; i < len_a+2;i++){
		np[i] = -1;
	}

	for(i = 1; i <= len_b;i++){
		if(hirsch_path[i] != -1){
			np[hirsch_path[i]] = i;
		}
	}

	free(hirsch_path);
	return np;
}

int* add_gap_info_to_hirsch_path(int* hirsch_path,int len_a,int len_b)
{
	int i,j;
	int a = 0;
	int b = 0;

	int* np = 0;
	np =malloc(sizeof(int)*(len_a+len_b+2));
	for(i =0; i < len_a+len_b+2;i++){
		np[i] = 0;
	}

	j = 1;
	b = -1;
	if(hirsch_path[1] == -1){
		np[j] = 2;
		j++;
	}else{
		if(hirsch_path[1] != 1){
			for ( a = 0;a < hirsch_path[1] -1;a++){
				np[j] = 1;
				j++;
			}
			np[j] = 0;
			j++;
		}else{
			np[j] = 0;
			j++;
		}
	}
	b = hirsch_path[1];
	
	/*for ( i= 0;i <= len_a;i++){
		fprintf(stderr,"%d,",hirsch_path[i]);
	} 
	fprintf(stderr,"\n");*/
	
	for(i = 2; i <= len_a;i++){
	
		if(hirsch_path[i] == -1){
			np[j] = 2;
			j++;
		}else{
			if(hirsch_path[i]-1 != b && b != -1){
				for ( a = 0;a < hirsch_path[i] - b-1;a++){
					np[j] = 1;
					j++;
				}
				np[j] = 0;
				j++;
			}else{
				np[j] = 0;
				j++;
			}
		}
		b = hirsch_path[i];
	}
	
	
	
	
	
	if(hirsch_path[len_a] < len_b && hirsch_path[len_a] != -1){
	//	fprintf(stderr,"WARNING:%d	%d\n",hirsch_path[len_a],len_b);
		for ( a = 0;a < len_b - hirsch_path[len_a];a++){
			np[j] = 1;
			j++;
		}
	} 
	np[0] = j-1;
	np[j] = 3;
	np = realloc(np,sizeof(int)* (np[0]+2));
	//for ( i= 0;i <= np[0];i++){
	//	fprintf(stderr,"%d,",np[i]);
	//} 
	//fprintf(stderr,"\n");

	free(hirsch_path);

	//add gap info..
	i = 2;
	while(np[i] != 3){
		if ((np[i-1] &3) && !(np[i] & 3)){
			if(np[i-1] & 8){
				np[i-1] += 8;
			}else{
				np[i-1] |= 16;
			}
		}else if (!(np[i-1] & 3) &&(np[i] &3)){
			np[i] |= 4;
		}else if ((np[i-1] & 1) && (np[i] & 1)){
			np[i] |= 8;
		}else if ((np[i-1] & 2) && (np[i] & 2)){
			np[i] |= 8;
		}
		i++;
	}
	//add terminal gap...
	i = 1;
	while(np[i] != 0){
		np[i] |= 32;
		i++;
	}
	j = i;
	i = np[0];
	while(np[i] != 0){
		np[i] |= 32;
		i--;
	}
	//for ( i= 0;i <= np[0];i++){
	//	fprintf(stderr,"%d,",np[i]);
	//} 
	//fprintf(stderr,"\n");
	return np;
}

/*
int* foward_pp_dyn(int* path, struct dp_matrix *dp,const int* prof1,const int* prof2,const int len_a,const int len_b)
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

	trace[0][0] = 32;


	s[0].a = 0;
	s[0].ga = -INFTY;
	s[0].gb = -INFTY;
	//init of first row;
	tracep = trace[0];

	for (j = 1; j < len_b;j++){
		s[j].a = -INFTY;
		
		s[j].ga = s[j-1].a+prof2[29];
		if (s[j-1].ga+prof2[29] > s[j].ga){
			s[j].ga = s[j-1].ga+prof2[29];
		}
		s[j].gb = -INFTY;
		tracep[j] = 8;
	}
	
	s[len_b].a = -INFTY;
	s[len_b].ga = -INFTY;
	s[len_b].gb = -INFTY;
	
	for ( i = 1;i <len_a;i++){
		prof1 += 64;

		c = 1;
		for (j = 26; j--;){
			if(prof1[j]){
				freq[c] = j;
				c++;	
			}
		}
		freq[0] = c;
		
		tracep = trace[i];
		pa = s[0].a;
		pga = s[0].ga;
		pgb = s[0].gb;
		s[0].a = -INFTY;
		s[0].ga = -INFTY;
		
		s[0].gb = pa+prof1[29];
		if(pgb+prof1[29] > s[0].gb){
			s[0].gb = pgb+prof1[29];
		}
	
		tracep[0] = 16;

		for (j = 1; j < len_b;j++){
			prof2 += 64;
			ca = s[j].a;

			c = 1;
			if((pga += prof2[-37]) > pa){
				pa = pga;
				c = 2;
			}
			if((pgb += prof1[-37]) > pa){
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
			
			s[j].ga = s[j-1].a+prof2[27];
			if (s[j-1].ga+prof2[28] > s[j].ga){
				s[j].ga = s[j-1].ga+prof2[28];
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
		
	
		prof2 += 64;
		//LAST CELL (0)
		ca = s[len_b].a;

		c = 1;
		if((pga+=prof2[-37]) > pa){
			pa = pga;
			c = 2;
		}
		if((pgb+=prof1[-37]) > pa){
			pa = pgb;
			c = 4;
		}
		
		prof2 += 32;
		for (pga = freq[0];--pga;){
			pgb = freq[pga];
			pa += prof1[pgb]*prof2[pgb];
		}
		prof2 -= 32;
		
		s[len_b].a = pa;
		
		s[len_b].ga = -INFTY;
		
		pgb = s[len_b].gb;
		s[len_b].gb = ca+prof1[27]+prof1[29];
 		if(pgb+prof1[29] > s[len_b].gb){
			s[len_b].gb = pgb+prof1[29];
			c |= 16;
		}
		tracep[len_b] = c;	
		prof2 -= len_b << 6;
		
	}
	prof1 += 64;
	
	c = 1;
	for (j = 26; j--;){
		if(prof1[j]){
			freq[c] = j;
			c++;	
		}
	}
	freq[0] = c;
	
	tracep = trace[len_a];

	pa = s[0].a;
	pga = s[0].ga;
	pgb = s[0].gb;
	s[0].a = -INFTY;
	s[0].ga = -INFTY;

	s[0].gb = pa+prof1[29];
	if(pgb+prof1[29] > s[0].gb){
		s[0].gb = pgb+prof1[29];
	}
	tracep[0] = 16;

	for (j = 1;j< len_b;j++){	

		prof2 += 64;
		ca = s[j].a;

		c = 1;

		if((pga+=prof2[-37]) > pa){
			pa = pga;
			c = 2;
		}

		if((pgb+=prof1[-37]) > pa){
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
		s[j].ga = s[j-1].a+prof2[27]+prof2[29];
		if (s[j-1].ga+prof2[29] > s[j].ga){
			s[j].ga = s[j-1].ga+prof2[29];
			c |= 8;
		}	
		pgb = s[j].gb;
		s[j].gb = -INFTY;	
		
		tracep[j] = c;
		pa = ca;
	}
	prof2 += 64;

	ca = s[len_b].a;
	
	c = 1;
	
	if((pga+=prof2[-37]) > pa){
		pa = pga;
		c = 2;
	}
	if((pgb+=prof1[-37]) > pa){
		pa = pgb;
		c = 4;
	}
	prof2+=32;
	for (pga = freq[0];--pga;){	
		pgb = freq[pga];
		pa += prof1[pgb]*prof2[pgb];
	}
	prof2-=32;
	
	s[len_b].a = pa;
	
	s[len_b].ga = s[len_b-1].a+prof2[27]+prof2[29];
	if (s[len_b-1].ga+prof2[29] > s[len_b].ga){
		s[len_b].ga = s[len_b-1].ga+prof2[29];
		c |= 8;
	}
	
	pgb = s[len_b].gb;
	s[len_b].gb = ca+prof1[27]+prof1[29];
	if(pgb +prof1[29]> s[len_b].gb){
		s[len_b].gb = pgb+prof1[29];
		c |= 16;
	}	
	tracep[len_b] = c;

	pgb = s[len_b].gb;
	c = 2;
	if(s[len_b].ga > pgb){
		pgb = s[len_b].ga;
		c = 1;
	}
	if(s[len_b].a >= pgb){
		pgb = s[len_b].a;
		c = 0;
	}
	
	ca = c;
	
	i = len_a;
	j = len_b;
	c = 1;
	while(trace[i][j] < 32){
	//	fprintf(stderr,"%d->%d	%d:%d	%d:%d\n",c,trace[i][j],i,j,len_a,len_b);
		switch(ca){
			case 0:
				if (trace[i][j] & 2){
					ca = 1;
					if(i-1!= 0){
						path[c+1] |= 16;
	//					fprintf(stderr,"GAP_CLOSE\n");
					}else{
						path[c+1] |= 32+16;
					}
				}else if (trace[i][j] & 4){
					ca = 2;
					if(j-1!= 0){
						path[c+1] |= 16;
	//					fprintf(stderr,"GAP_CLOSE\n");
					}else{
						path[c+1] |= 32+16;
					}
				}

				//path[c] = 0;
				i--;
				j--;
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
				j--;
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
					if(j !=0 && j != len_b){
	//					fprintf(stderr,"GAP_OPEN\n");
						path[c] |= 4;
					}else{
						path[c] |= 32+4;
					}
					
				}
				path[c] |= 2;
				i--;
			break;
		}
		c++;
	}
	
	
	
	path[0] = c-1;
	path[c] = 3;
	path[c+1] = pgb;
	
	j = path[0];
	for(i =0 ;i < path[0]/2;i++){
		c = path[i+1];
		path[i+1] = path[j-i];
		path[j -i] = c;
	}
	return path;
}

int* backward_pp_dyn(int* path, struct dp_matrix *dp,const int* prof1,const int* prof2,const int len_a,const int len_b)
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

	prof1+= 64;
	prof2 += 64;


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

*/
