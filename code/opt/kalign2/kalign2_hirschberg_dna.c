/*
	kalign2_hirschberg_dna.c
	
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
#include "kalign2_hirschberg_dna.h"
#define MAX(a, b) (a > b ? a : b)
#define MAX3(a,b,c) MAX(MAX(a,b),c)



int** dna_alignment(struct alignment* aln,int* tree,float**submatrix, int** map,float strength)
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
		//fprintf(stderr,"Aligning:%d %d->%d	done:%0.2f\n",a,b,c,((float)(i+1)/(float)numseq)*100);
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
			profile[a] = dna_make_profile(profile[a],aln->s[a],len_a,submatrix);
		}
		if (b < numseq){
			profile[b] = dna_make_profile(profile[b],aln->s[b],len_b,submatrix);
		}
		fprintf(stderr,"Saving mem...\n");
	
		dna_set_gap_penalties(profile[a],len_a,aln->nsip[b],strength,aln->nsip[a]);
		dna_set_gap_penalties(profile[b],len_b,aln->nsip[a],strength,aln->nsip[b]);

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
				map[c] = hirsch_dna_ss_dyn(submatrix,aln->s[a],aln->s[b],hm,map[c]);
			}else{
				hm->enda = len_b;
				hm->endb = len_a;
				hm->len_a = len_b;
				hm->len_b = len_a;
				map[c] = hirsch_dna_ps_dyn(profile[b],aln->s[a],hm,map[c],aln->nsip[b]);
				map[c] = mirror_hirsch_path(map[c],len_a,len_b);
			}
		}else{
			if(b < numseq){
				map[c] = hirsch_dna_ps_dyn(profile[a],aln->s[b],hm,map[c],aln->nsip[a]);
			}else{
				if(len_a < len_b){
					map[c] = hirsch_dna_pp_dyn(profile[a],profile[b],hm,map[c]);
				}else{
					hm->enda = len_b;
					hm->endb = len_a;
					hm->len_a = len_b;
					hm->len_b = len_a;
					map[c] = hirsch_dna_pp_dyn(profile[b],profile[a],hm,map[c]);
					map[c] = mirror_hirsch_path(map[c],len_a,len_b);
				}
			}
		}
		map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);

		if(i != numseq-2){
			profile[c] = malloc(sizeof(float)*22*(map[c][0]+2));
			profile[c] = dna_update(profile[a],profile[b],profile[c],map[c],aln->nsip[a],aln->nsip[b]);
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
	//free(profile[numprofiles-1]);
	free(profile);
	hirsch_mem_free(hm);
	for (i = 32;i--;){
		free(submatrix[i]);
	}
	free(submatrix);
	return map;
}

int** dna_alignment_against_a(struct alignment* aln,int* tree,float**submatrix, int** map,float strength)
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
		//fprintf(stderr,"Aligning:%d %d->%d	done:%0.2f\n",a,b,c,((float)(i+1)/(float)numseq)*100);
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
			profile[a] = dna_make_profile(profile[a],aln->s[a],len_a,submatrix);
		}
		if (b < numseq){
			profile[b] = dna_make_profile(profile[b],aln->s[b],len_b,submatrix);
		}
		
	
		dna_set_gap_penalties(profile[a],len_a,1,strength,1);//aln->nsip[b]);
		dna_set_gap_penalties(profile[b],len_b,1,strength,1);//aln->nsip[a]);

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
				map[c] = hirsch_dna_ss_dyn(submatrix,aln->s[a],aln->s[b],hm,map[c]);
			}else{
				hm->enda = len_b;
				hm->endb = len_a;
				hm->len_a = len_b;
				hm->len_b = len_a;
				map[c] = hirsch_dna_ps_dyn(profile[b],aln->s[a],hm,map[c],aln->nsip[b]);
				map[c] = mirror_hirsch_path(map[c],len_a,len_b);
			}
		}else{
			if(b < numseq){
				map[c] = hirsch_dna_ps_dyn(profile[a],aln->s[b],hm,map[c],1);//aln->nsip[a]);
			}else{
				if(len_a < len_b){
					map[c] = hirsch_dna_pp_dyn(profile[a],profile[b],hm,map[c]);
				}else{
					hm->enda = len_b;
					hm->endb = len_a;
					hm->len_a = len_b;
					hm->len_b = len_a;
					map[c] = hirsch_dna_pp_dyn(profile[b],profile[a],hm,map[c]);
					map[c] = mirror_hirsch_path(map[c],len_a,len_b);
				}
			}
		}

		
		map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);

		if(i != numseq-2){
			profile[c] = malloc(sizeof(float)*22*(map[c][0]+2));
			profile[c] = dna_update_only_a(profile[a],profile[b],profile[c],map[c],aln->nsip[a],aln->nsip[b]);
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
	//free(profile[numprofiles-1]);
	free(profile);
	hirsch_mem_free(hm);
	for (i = 32;i--;){
		free(submatrix[i]);
	}
	free(submatrix);
	return map;
}



int* hirsch_dna_ss_dyn(float**subm, const int* seq1,const int* seq2,struct hirsch_mem* hm, int* hirsch_path)
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
	hm->f = foward_hirsch_dna_ss_dyn(subm,seq1,seq2,hm);

	hm->starta = mid;
	hm->enda = old_cor[1];
	//fprintf(stderr,"Backward:%d-%d	%d-%d\n",hm->starta,hm->enda,hm->startb,hm->endb);
	hm->b = backward_hirsch_dna_ss_dyn(subm,seq1,seq2,hm);


	hirsch_path = hirsch_align_two_dna_ss_vector(subm,seq1,seq2,hm,hirsch_path,input_states,old_cor);
	return hirsch_path;
}

int* hirsch_align_two_dna_ss_vector(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[])
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
	float middle =  (hm->endb - hm->startb)/2 + hm->startb;
	float sub = 0.0;
	
	i = hm->startb;
	c = -1;
	for(i = hm->startb; i < hm->endb;i++){
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
	i = hm->endb;
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
		if(f[i].gb+b[i].gb - tgpe-sub > max){
			max = f[i].gb+b[i].gb - tgpe-sub;
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
			hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

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
			hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
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
			hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			

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
			hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
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
			hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
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
			hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
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
			hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

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
			hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
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
			hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

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
			hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
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
			hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);

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
			hirsch_path = hirsch_dna_ss_dyn(subm,seq1,seq2,hm,hirsch_path);
			break;
	}
		
	return hirsch_path;
}



struct states* foward_hirsch_dna_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm)
{
	struct states* s = hm->f;
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
			s[j].ga = MAX(s[j-1].ga,s[j-1].a)-tgpe;
			
			s[j].gb = -FLOATINFTY;
		}		
	}else{

		for (j = startb+1; j < endb;j++){
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j-1].ga - gpe,s[j-1].a-gpo);
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
		if(startb == 0){
			s[startb].gb = MAX(pgb,pa) - tgpe;
		}else{
			s[startb].gb = MAX(pgb - gpe,pa - gpo);
		}
		for (j = startb+1; j < endb;j++){
			ca = s[j].a;
			pa = MAX3(pa,pga-gpo,pgb-gpo);
			pa += subp[seq2[j]];
			
			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = MAX(s[j-1].ga-gpe,s[j-1].a-gpo);
			
			pgb = s[j].gb;
			
			s[j].gb = MAX(pgb-gpe ,ca-gpo);
			
			pa = ca;
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

struct states* backward_hirsch_dna_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm)
{

	struct states* s = hm->b;
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
			s[j].ga = MAX(s[j+1].ga,s[j+1].a)-tgpe;
			s[j].gb = -FLOATINFTY;
		}
	}else{
		for(j = endb-1;j > startb;j--){
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j+1].ga-gpe,s[j+1].a-gpo);	
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
			s[endb].gb = MAX(pgb,pa)-tgpe;
		}else{
			s[endb].gb = MAX(pgb-gpe,pa-gpo);
		}

		for(j = endb-1;j > startb;j--){

			ca = s[j].a;
			pa = MAX3(pa,pga - gpo,pgb-gpo);
			
			pa += subp[seq2[j]];

			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = MAX(s[j+1].ga-gpe,s[j+1].a-gpo);
			
			pgb = s[j].gb;

			s[j].gb = MAX(pgb-gpe,ca-gpo);
			
			pa = ca;
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


int* hirsch_dna_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm, int* hirsch_path,int sip)
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
	hm->f = foward_hirsch_dna_ps_dyn(prof1,seq2,hm,sip);
	
	/*int i;
	fprintf(stderr,"FOWARD\n");
	for (i = hm->startb; i <= hm->endb;i++){
		fprintf(stderr,"%d	%d	%d\n",hm->f[i].a,hm->f[i].ga,hm->f[i].gb);
	}*/

	hm->starta = mid;
	hm->enda = old_cor[1];
	hm->b = backward_hirsch_dna_ps_dyn(prof1,seq2,hm,sip);
	
	/*fprintf(stderr,"BaCKWARD\n");
	for (i = hm->startb; i <= hm->endb;i++){
		fprintf(stderr,"%d	%d	%d\n",hm->b[i].a,hm->b[i].ga,hm->b[i].gb);
	}*/

	hirsch_path = hirsch_align_two_dna_ps_vector(prof1,seq2,hm,hirsch_path,input_states,old_cor,sip);
	return hirsch_path;
}



int* hirsch_align_two_dna_ps_vector(const float* prof1,const int* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],int sip)
{
	struct states* f = hm->f;
 	struct states* b = hm->b;
	int i,j,c;
	int transition = -1;
	
	const int open = gpo * sip;
	
	
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
	float middle =  (hm->endb - hm->startb)/2 + hm->startb;
	float sub = 0.0;
	
	
	prof1+= (22 * (old_cor[4]+1));
	
	i = hm->startb;
	c = -1;
	for(i = hm->startb; i < hm->endb;i++){
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
		if(f[i].a+b[i].gb+prof1[8]-sub > max){
			max = f[i].a+b[i].gb+prof1[8]-sub;
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
			if(f[i].gb+b[i].gb+prof1[10]-sub > max){
				max = f[i].gb+b[i].gb+prof1[10]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}else{
			if(f[i].gb+b[i].gb+prof1[9]-sub > max){
				max = f[i].gb+b[i].gb+prof1[9]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}
		if(f[i].gb+b[i].a+prof1[8-22]-sub > max){
			max = f[i].gb+b[i].a+prof1[8-22]-sub;
	//		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
			transition = 7;
			c = i;
		}
	}
	i = hm->endb;
	sub = abs(middle -i);
	sub /= 1000; 
	if(f[i].a+b[i].gb+prof1[8]-sub > max){
		max = f[i].a+b[i].gb+prof1[8]-sub;
	//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
		transition = 3;
		c = i;
	}
	if(hm->endb == hm->len_b){
		if(f[i].gb+b[i].gb+prof1[10]-sub > max){
			max = f[i].gb+b[i].gb+prof1[10]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}	
	}else{
		if(f[i].gb+b[i].gb+prof1[9]-sub > max){
			max = f[i].gb+b[i].gb+prof1[0]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}
	}
	
	
	
	prof1-= (22 * (old_cor[4]+1));
	
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
			hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

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
			hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
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
			hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

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
			hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
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
			hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

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
			hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
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
			hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

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
			hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
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
			hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);			


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
			hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
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
			hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);

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
			hirsch_path = hirsch_dna_ps_dyn(prof1,seq2,hm,hirsch_path,sip);
			break;
	}
		
	return hirsch_path;
}

struct states* foward_hirsch_dna_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm,int sip)
{
	//unsigned int freq[26];
	struct states* s = hm->f;
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
	
	const float open = gpo * sip;
	const float ext = gpe *sip;
	const float text = tgpe * sip;
	
	
	
	prof1 += (starta) * 22;
	s[startb].a = s[0].a;
	s[startb].ga = s[0].ga;
	s[startb].gb = s[0].gb;
	if(startb == 0){
		for (j = startb+1; j < endb;j++){
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j-1].ga,s[j-1].a) - text;
			s[j].gb = -FLOATINFTY;
		}	
	}else{
		for (j = startb+1; j < endb;j++){
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);	
			s[j].gb = -FLOATINFTY;
		}
	}
	
	
	s[endb].a = -FLOATINFTY;
	s[endb].ga = -FLOATINFTY;
	s[endb].gb = -FLOATINFTY;
	seq2--;

	for (i = starta;i < enda;i++){
		prof1 += 22;
		pa = s[startb].a;
		pga = s[startb].ga;
		pgb = s[startb].gb;
		s[startb].a = -FLOATINFTY;
		s[startb].ga = -FLOATINFTY;
		if(startb == 0){
			s[startb].gb = MAX(pgb,pa)+prof1[10];
		}else{
			s[startb].gb = MAX(pgb+prof1[9],pa+prof1[8]);
		}
		for (j = startb+1; j < endb;j++){
			ca = s[j].a;
			pa = MAX3(pa,pga -open,pgb + prof1[-14]);
			pa += prof1[11 + seq2[j]];


			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = MAX(s[j-1].ga-ext,s[j-1].a-open);
			
			pgb = s[j].gb;

			s[j].gb = MAX(pgb+prof1[9],ca+prof1[8]);

			pa = ca;
		}	
		ca = s[j].a;

		pa = MAX3(pa,pga -open,pgb + prof1[-14]);
			
		pa += prof1[11 + seq2[j]];


		s[j].a = pa;

		s[j].ga = -FLOATINFTY;//MAX(s[j-1].ga-ext,s[j-1].a-open);
				
		if (hm->endb != hm->len_b){
			s[j].gb = MAX(s[j].gb+prof1[9] ,ca+prof1[8]);
		}else{
			s[j].gb = MAX(s[j].gb,ca)+ prof1[10];
		}
	}
	prof1 -= 22 * (enda);
	return s;
}

struct states* backward_hirsch_dna_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm,int sip)
{
	//unsigned int freq[26];
	struct states* s = hm->b;
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
	
	const float open = gpo * sip;
	const float ext = gpe *sip;
	const float text = tgpe * sip;
	

	prof1 += (enda+1) * 22;

	s[endb].a = s[0].a;
	s[endb].ga = s[0].ga;
	s[endb].gb = s[0].gb;
	
	
	//init of first row;
	//j = endb-startb;
	if(endb == hm->len_b){
		for(j = endb-1;j > startb;j--){
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j+1].ga,s[j+1].a)-text;	
			s[j].gb = -FLOATINFTY;
		}
	}else{
		for(j = endb-1;j > startb;j--){
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);
			s[j].gb = -FLOATINFTY;
		}
	}
	
	s[startb].a = -FLOATINFTY;
	s[startb].ga = -FLOATINFTY;
	s[startb].gb = -FLOATINFTY;

	i = enda-starta;
	while(i--){
		prof1 -= 22;

		pa = s[endb].a;
		pga = s[endb].ga;
		pgb = s[endb].gb;
		s[endb].a = -FLOATINFTY;
		s[endb].ga = -FLOATINFTY;

		if(endb == hm->len_b){
			s[endb].gb = MAX(pgb,pa) +prof1[10];
		}else{
			s[endb].gb = MAX(pgb+prof1[9],pa+prof1[8]);
		}

		for(j = endb-1;j > startb;j--){
			ca = s[j].a;
			pa = MAX3(pa,pga - open,pgb +prof1[30]);
			pa += prof1[11 + seq2[j]];

			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = MAX(s[j+1].ga-ext,s[j+1].a-open);

			pgb = s[j].gb;

			s[j].gb = MAX(pgb+prof1[9],ca+prof1[8]);
			
			pa = ca;
		}
		ca = s[j].a;

		pa = MAX3(pa,pga - open,pgb +prof1[30]);
		pa += prof1[11 + seq2[j]];

		s[j].a = pa;
			

		s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga-ext,s[j+1].a-open);
		if(hm->startb){
			s[j].gb = MAX(s[j].gb+prof1[9], ca+prof1[8]);
		}else{
			s[j].gb = MAX(s[j].gb,ca)+prof1[10];
		}
	}		
	return s;
}




int* hirsch_dna_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm, int* hirsch_path)
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
	hm->f = foward_hirsch_dna_pp_dyn(prof1,prof2,hm);
	/*int i;
	fprintf(stderr,"FOWARD\n");
	for (i = hm->startb; i <= hm->endb;i++){
		fprintf(stderr,"%d	%d	%d\n",hm->f[i].a,hm->f[i].ga,hm->f[i].gb);
	}*/

	hm->starta = mid;
	hm->enda = old_cor[1];
	hm->b = backward_hirsch_dna_pp_dyn(prof1,prof2,hm);
	/*fprintf(stderr,"BaCKWARD\n");

	for (i = hm->startb; i <= hm->endb;i++){
		fprintf(stderr,"%d	%d	%d\n",hm->b[i].a,hm->b[i].ga,hm->b[i].gb);
	}*/

	hirsch_path = hirsch_align_two_dna_pp_vector(prof1,prof2,hm,hirsch_path,input_states,old_cor);
	return hirsch_path;
}



int* hirsch_align_two_dna_pp_vector(const float* prof1,const float* prof2,struct hirsch_mem* hm,int* hirsch_path, float input_states[],int old_cor[])
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
	float middle =  (hm->endb - hm->startb)/2 + hm->startb;
	float sub = 0.0;
	
	
	prof1+= (22 * (old_cor[4]+1));
	prof2 += (22 * (hm->startb));
	
	i = hm->startb;
	c = -1;
	for(i = hm->startb; i < hm->endb;i++){
		sub = abs(middle -i);
		sub /= 1000; 
		prof2 += 22;
		if(f[i].a+b[i].a-sub > max){
			max = f[i].a+b[i].a-sub;
	//		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
			transition = 1;
			c = i;
		}
		if(f[i].a+b[i].ga+prof2[8]-sub > max){
			max = f[i].a+b[i].ga+prof2[8]-sub;
	//		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
			transition = 2;
			c = i;
		}
		if(f[i].a+b[i].gb+prof1[8]-sub > max){
			max = f[i].a+b[i].gb+prof1[8]-sub;
	//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
			transition = 3;
			c = i;
		}
		if(f[i].ga+b[i].a+prof2[-14]-sub > max){
			max = f[i].ga+b[i].a+prof2[-14]-sub;
	//		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
			transition = 5;
			c = i;
		}


		if(hm->startb == 0){
			if(f[i].gb+b[i].gb+prof1[10]-sub > max){
				max = f[i].gb+b[i].gb+prof1[10]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}else{
			if(f[i].gb+b[i].gb+prof1[9]-sub > max){
				max = f[i].gb+b[i].gb+prof1[9]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}
		if(f[i].gb+b[i].a+prof1[-14]-sub > max){
			max = f[i].gb+b[i].a+prof1[-14]-sub;
	//		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
			transition = 7;
			c = i;
		}
	}
	i = hm->endb;
	sub = abs(middle -i);
	sub /= 1000; 
	if(f[i].a+b[i].gb+prof1[8]-sub > max){
		max = f[i].a+b[i].gb+prof1[8]-sub;
	//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
		transition = 3;
		c = i;
	}
	if(hm->endb == hm->len_b){
		if(f[i].gb+b[i].gb+prof1[10]-sub > max){
			max = f[i].gb+b[i].gb+prof1[10]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}	
	}else{
		if(f[i].gb+b[i].gb+prof1[9]-sub > max){
			max = f[i].gb+b[i].gb+prof1[9]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}
	}
	
	
	
	prof1-= (22 * (old_cor[4]+1));
	prof2 -= (hm->endb *22);
	
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
			hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);

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
			hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);
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
			hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);

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
			hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);
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
			hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);

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
			hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);
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
			hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);

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
			hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);
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
			hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);

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
			hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);
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
			hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);

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
			hirsch_path = hirsch_dna_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
	}
		
	return hirsch_path;
}

struct states* foward_hirsch_dna_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm)
{
	struct states* s = hm->f;

	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	register int i = 0;
	register int j = 0;
	
	
	
	prof1 += (hm->starta) * 22;
	prof2 +=  (hm->startb) * 22;
	s[hm->startb].a = s[0].a;
	s[hm->startb].ga = s[0].ga;
	s[hm->startb].gb = s[0].gb;
	if(hm->startb == 0){
		for (j = hm->startb+1; j < hm->endb;j++){
			prof2+=22;
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j-1].ga,s[j-1].a)+prof2[10];
			s[j].gb = -FLOATINFTY;
		}	
		prof2 += 22;	
	}else{

		for (j = hm->startb+1; j < hm->endb;j++){
			prof2+=22;
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j-1].ga+prof2[9],s[j-1].a+prof2[8]);
			s[j].gb = -FLOATINFTY;
		}
		prof2 += 22;
	}
	
	prof2 -= (hm->endb-hm->startb) * 22;
	
	s[hm->endb].a = -FLOATINFTY;
	s[hm->endb].ga = -FLOATINFTY;
	s[hm->endb].gb = -FLOATINFTY;


	for (i = hm->starta;i < hm->enda;i++){
		prof1 += 22;

		pa = s[hm->startb].a;
		pga = s[hm->startb].ga;
		pgb = s[hm->startb].gb;
		s[hm->startb].a = -FLOATINFTY;
		s[hm->startb].ga = -FLOATINFTY;
		if(hm->startb == 0){
			s[hm->startb].gb = MAX(pgb,pa)+ prof1[10];
		}else{
			s[hm->startb].gb = MAX(pgb+prof1[9],pa+prof1[8]);
		}
		for (j = hm->startb+1; j < hm->endb;j++){
			prof2 += 22;
			ca = s[j].a;
			pa = MAX3(pa,pga + prof2[-14],pgb + prof1[-14]);

			prof2 += 11;

			pa += prof1[0]*prof2[0];
			pa += prof1[1]*prof2[1];
			pa += prof1[2]*prof2[2];
			pa += prof1[3]*prof2[3];
			pa += prof1[4]*prof2[4];
			pa += prof1[5]*prof2[5];
			pa += prof1[6]*prof2[6];
			pa += prof1[7]*prof2[7];
			
			
			prof2 -= 11;	

			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = MAX(s[j-1].ga+prof2[9],s[j-1].a+prof2[8]);
			
			pgb = s[j].gb;

			s[j].gb = MAX(pgb+prof1[9] ,ca+prof1[8]);

			pa = ca;
		}
		prof2 += 22;
		ca = s[j].a;
			
		pa = MAX3(pa,pga + prof2[-14],pgb + prof1[-14]);
		prof2 += 11;

		pa += prof1[0]*prof2[0];
		pa += prof1[1]*prof2[1];
		pa += prof1[2]*prof2[2];
		pa += prof1[3]*prof2[3];
		pa += prof1[4]*prof2[4];
		pa += prof1[5]*prof2[5];
		pa += prof1[6]*prof2[6];
		pa += prof1[7]*prof2[7];
		
		prof2 -= 11;	

		s[j].a = pa;

		s[j].ga = -FLOATINFTY;

		if (hm->endb != hm->len_b){
			s[j].gb = MAX(s[j].gb+prof1[9] ,ca+prof1[8]);
		}else{
			s[j].gb = MAX(s[j].gb,ca)+ prof1[10];
		}
		
		
		prof2 -= (hm->endb-hm->startb) * 22;
		
	}
	prof1 -= 22 * (hm->enda);
	return s;
}

struct states* backward_hirsch_dna_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm)
{
	struct states* s = hm->b;
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	register int i = 0;
	register int j = 0;

	prof1 += (hm->enda+1) * 22;
	prof2 += (hm->endb+1) * 22;
	s[hm->endb].a = s[0].a;
	s[hm->endb].ga = s[0].ga;
	s[hm->endb].gb = s[0].gb;
	
	
	//init of first row;
	//j = endb-startb;
	if(hm->endb == hm->len_b){
		
		for(j = hm->endb-1;j > hm->startb;j--){
			prof2 -= 22;
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j+1].ga,s[j+1].a)+prof2[10];
			s[j].gb = -FLOATINFTY;
		}
		prof2 -= 22;
	}else{
		for(j = hm->endb-1;j > hm->startb;j--){
			prof2 -= 22;
			s[j].a = -FLOATINFTY;
			s[j].ga = MAX(s[j+1].ga+prof2[9],s[j+1].a+prof2[8]);
			s[j].gb = -FLOATINFTY;
		}
		prof2 -= 22;
	}
	
	s[hm->startb].a = -FLOATINFTY;
	s[hm->startb].ga = -FLOATINFTY;
	s[hm->startb].gb = -FLOATINFTY;

	i = hm->enda-hm->starta;
	while(i--){
		prof1 -= 22;

		pa = s[hm->endb].a;
		pga = s[hm->endb].ga;
		pgb = s[hm->endb].gb;
		s[hm->endb].a = -FLOATINFTY;
		s[hm->endb].ga = -FLOATINFTY;

		if(hm->endb == hm->len_b){
			s[hm->endb].gb = MAX(pgb,pa)+prof1[10];
		}else{
			s[hm->endb].gb = MAX(pgb+prof1[9] ,pa+prof1[8]);
		}
		//j = endb-startb;
		prof2 += (hm->endb-hm->startb) *22;
		//while(j--){
		for(j = hm->endb-1;j > hm->startb;j--){
			prof2 -= 22;
			ca = s[j].a;

			pa = MAX3(pa,pga + prof2[30],pgb + prof1[30]);

			prof2 += 11;
			pa += prof1[0]*prof2[0];
			pa += prof1[1]*prof2[1];
			pa += prof1[2]*prof2[2];
			pa += prof1[3]*prof2[3];
			pa += prof1[4]*prof2[4];
			pa += prof1[5]*prof2[5];
			pa += prof1[6]*prof2[6];
			pa += prof1[7]*prof2[7];
			prof2 -= 11;

			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = MAX(s[j+1].ga+prof2[9], s[j+1].a+prof2[8]);
			
			pgb = s[j].gb;

			s[j].gb = MAX(pgb+prof1[9], ca+prof1[8]);

			pa = ca;
		}
		prof2 -= 22;
		ca = s[j].a;

		pa = MAX3(pa,pga + prof2[30],pgb + prof1[30]);
		
		prof2 += 11;
		pa += prof1[0]*prof2[0];
		pa += prof1[1]*prof2[1];
		pa += prof1[2]*prof2[2];
		pa += prof1[3]*prof2[3];
		pa += prof1[4]*prof2[4];
		pa += prof1[5]*prof2[5];
		pa += prof1[6]*prof2[6];
		pa += prof1[7]*prof2[7];
		prof2 -= 11;

		s[j].a = pa;
		
		//pga = s[j].ga;
		s[j].ga = -FLOATINFTY;//MAX(s[j+1].ga+prof2[28], s[j+1].a+prof2[27]);

		//pgb = s[j].gb;
		if(hm->startb){
			s[j].gb = MAX(s[j].gb+prof1[9], ca+prof1[8]);
		}else{
			s[j].gb = MAX(s[j].gb,ca)+prof1[10];
		}
	}		
	return s;
}
