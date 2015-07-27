/*
	kalign2_feature.c 
	
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
#include "kalign2_feature.h"

static int stride;
static int dim;
static int gpo_pos;
static int gpe_pos;
static int tgpe_pos;

int** feature_hirschberg_alignment(struct alignment* aln,int* tree,float**submatrix, int** map,struct feature_matrix* fm)
{
	struct hirsch_mem* hm = 0;
	int i,j,g,a,b,c;
	int len_a;
	int len_b;
	float** profile = 0;
	
	stride = (26+fm->mdim)*2 + 3; 
	dim = 26+fm->mdim;
	gpo_pos = (dim << 1) + 0;
	gpe_pos = (dim << 1) + 1;
	tgpe_pos = (dim << 1) + 2;

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
			profile[a] = make_unified_profile(profile[a],aln,a,submatrix,fm);
		}
		set_unified_gap_penalties(profile[a],len_a,aln->nsip[b]);
		
		if (b < numseq){
			profile[b] = make_unified_profile(profile[b],aln,b,submatrix,fm);
		}
		set_unified_gap_penalties(profile[b],len_b,aln->nsip[a]);
		
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
		//dim = 26;
	//	fprintf(stderr,"LENA:%d	LENB:%d	numseq:%d\n",len_a,len_b,numseq);
		if(len_a < len_b){
		//	fprintf(stderr,"normal\n");
			map[c] = feature_hirsch_pp_dyn(profile[a],profile[b],hm,map[c]);
		}else{
		//	fprintf(stderr,"goofy\n");
			hm->enda = len_b;
			hm->endb = len_a;
			hm->len_a = len_b;
			hm->len_b = len_a;
			map[c] = feature_hirsch_pp_dyn(profile[b],profile[a],hm,map[c]);
			map[c] = mirror_hirsch_path(map[c],len_a,len_b);
		}
		map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);

		if(i != numseq-2){
			profile[c] = malloc(sizeof(float)*stride*(map[c][0]+2));
			profile[c] = feature_hirschberg_update(profile[a],profile[b],profile[c],map[c],aln->nsip[a],aln->nsip[b]);
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
	free_feature_matrix(fm);
	return map;
}


float* feature_hirschberg_update(const float* profa,const float* profb,float* newp,int* path,int sipa,int sipb)
{
	int i,j,c;
	for (i = stride; i--;){
		newp[i] = profa[i] + profb[i];
	}
	
	profa += stride;
	profb += stride;
	newp += stride;

	c = 1;
	
	while(path[c] != 3){
		//Idea: limit the 'virtual' number of residues of one type to x.
		// i.e. only allow a maximum of 10 alanines to be registered in each column
		// the penalty for aligning a 'G' to this column will stay stable even when many (>10) alanines are present.
		// the difference in score between the 'correct' (all alanine) and incorrect (alanines + glycine) will not increase 
		// with the number of sequences. -> see Durbin pp 140
		
		if (!path[c]){
			//fprintf(stderr,"Align	%d\n",c);
			for (i = stride; i--;){
				newp[i] = profa[i] + profb[i];
			}
				
			
			profa += stride;
			profb += stride;
		}
		
		if (path[c] & 1){
			//fprintf(stderr,"Gap_A:%d\n",c);
			//printf("open:%d	ext:%d	%d	%d\n",si->nsip[a] * gpo,si->nsip[a] * gpe,si->nsip[a] * profb[41],si->nsip[a] * profb[46]);
			for (i = stride; i--;){
				newp[i] = profb[i];
			}
			profb += stride;
			if(!(path[c] & 20)){
				if(path[c] & 32){
					newp[25] += sipa;//1;
					i = tgpe*sipa;
				}else{
					newp[24] += sipa;//1;
					i = gpe*sipa;
				}
				
				for (j = dim; j < dim+23;j++){
					newp[j] -=i;
				}
			}else{
			if (path[c] & 16){ 
	//			fprintf(stderr,"close_open");
				if(path[c] & 32){
					newp[25] += sipa;//1;
					i = tgpe*sipa;
					newp[23] += sipa;//1;
					i += gpo*sipa;
				}else{
					newp[23] += sipa;//1;
					i = gpo*sipa;
				}
								
				for (j = dim; j < dim+23;j++){
					newp[j] -=i;
				}
			}
			if (path[c] & 4){ 
	//			fprintf(stderr,"Gap_open");
				if(path[c] & 32){
					newp[25] += sipa;//1;
					i = tgpe*sipa;
					newp[23] += sipa;//1;
					i += gpo*sipa;
				}else{
					newp[23] += sipa;//1;
					i = gpo*sipa;
				}
				for (j = dim; j < dim+23;j++){
					newp[j] -=i;
				}
			}
			}		
			
		}
		if (path[c] & 2){
			//fprintf(stderr,"Gap_B:%d\n",c);
			//printf("open:%d	ext:%d	%d	%d\n",si->nsip[b] * gpo,si->nsip[b] * gpe,profa[26],profa[27]);
			for (i = stride; i--;){
				newp[i] = profa[i];
			}
			profa+=stride;
			if(!(path[c] & 20)){
				if(path[c] & 32){
					newp[25] += sipb;//1;
					i = tgpe*sipb;
				}else{
					newp[24] += sipb;//1;
					i = gpe*sipb;
				}
				for (j = dim; j < dim+23;j++){
					newp[j] -=i;
				}
			}else{
			if (path[c] & 16){
	//			fprintf(stderr,"close_open");
				if(path[c] & 32){
					newp[25] += sipb;//1;
					i =  tgpe*sipb;
					newp[23] += sipb;//1;
					i +=  gpo*sipb;
				}else{
					newp[23] += sipb;//1;
					i =  gpo*sipb;
				}
				for (j = dim; j < dim+23;j++){
					newp[j] -=i;
				}
			}
			if (path[c] & 4){
	//			fprintf(stderr,"Gap_open");
				if(path[c] & 32){
					newp[25] += sipb;//1;
					i = tgpe*sipb;
					newp[23] += sipb;//1;
					i += gpo*sipb;
				}else{
					newp[23] += sipb;//1;
					i = gpo*sipb;
				}
				
				for (j = dim; j < dim+23;j++){
					newp[j] -=i;
				}
			}
			}		
		}
		newp += stride;
		c++;
	}
	for (i = stride; i--;){
		newp[i] =  profa[i] + profb[i];
	}	
	newp -= (path[0]+1) * stride;
	return newp;
}



float* make_unified_profile(float* prof,struct alignment* aln, int num,float** subm,struct feature_matrix* fm)
{
	struct feature* f = aln->ft[num];
	int i,j,c;
	int* seq = aln->s[num];
	//detemine minimim width of profile... 
	//stride = (26+fm->mdim)*2 + 3; 
	
	int len = aln->sl[num];
	prof = malloc(sizeof(float)*(len+2)*stride);
	prof +=  (stride *(len+1));
	for (i = 0;i < stride;i++){
		prof[i] = 0;
	}
	prof[23+dim] = -gpo;
	prof[24+dim] = -gpe;
	prof[25+dim] = -tgpe;

	
	i = len;
	while(i--){
		prof -= stride;

		for (j = 0;j < stride;j++){
			prof[j] = 0;
		}
		c = seq[i];
		
		prof[c] += 1;
		
		prof += dim;
		
		for(j = 0; j < 23;j++){
			prof[j] = subm[c][j];
		}
		prof[23] = -gpo;
		prof[24] = -gpe;
		prof[25] = -tgpe;
		prof -= dim;
	}
	prof -= stride;
	for (i = 0;i < stride;i++){
		prof[i] = 0;
	}
	prof[23+dim] = -gpo;
	prof[24+dim] = -gpe;
	prof[25+dim] = -tgpe;	

	while(f){
		if(f->color != -1){
			if(f->start < len && f->end < len){
				for (i = f->start;i <= f->end;i++){
					prof[i*stride+26 + f->color] += 1;
					//prof[i*stride+dim+26 + f->color] += 75;
					//fprintf(stderr,"FOUND on  %d : %s	%s\n",num,f->type,f->note);
					for ( j = 0 ; j < fm->mdim ;j++){
						prof[i*stride+dim+26+j] += fm->m[f->color][j];
					}
				}
			}
		}
		f = f->next;
	}	

	//exit(0);
 	return prof;
}


void set_unified_gap_penalties(float* prof,int len,int nsip)
{
	int i;
	
	prof +=  (stride *(len+1));
	prof[gpo_pos] = prof[dim+23]*nsip;
	prof[gpe_pos] = prof[dim+24]*nsip;
	prof[tgpe_pos] = prof[dim+25]*nsip;
	i = len+1;
	while(i--){
		prof -= stride;
		prof[gpo_pos] = prof[dim+23]*nsip;
		prof[gpe_pos] = prof[dim+24]*nsip;
		prof[tgpe_pos] = prof[dim+25]*nsip;
	}
}



int* feature_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm, int* hirsch_path)
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
	hm->f = feature_foward_hirsch_pp_dyn(prof1,prof2,hm);
	/*int i;
	fprintf(stderr,"FOWARD\n");
	for (i = hm->startb; i <= hm->endb;i++){
		fprintf(stderr,"%d	%d	%d\n",hm->f[i].a,hm->f[i].ga,hm->f[i].gb);
	}*/

	hm->starta = mid;
	hm->enda = old_cor[1];
	hm->b = feature_backward_hirsch_pp_dyn(prof1,prof2,hm);
	/*fprintf(stderr,"BaCKWARD\n");

	for (i = hm->startb; i <= hm->endb;i++){
		fprintf(stderr,"%d	%d	%d\n",hm->b[i].a,hm->b[i].ga,hm->b[i].gb);
	}*/

	hirsch_path = feature_hirsch_align_two_pp_vector(prof1,prof2,hm,hirsch_path,input_states,old_cor);
	return hirsch_path;
}



int* feature_hirsch_align_two_pp_vector(const float* prof1,const float* prof2,struct hirsch_mem* hm,int* hirsch_path, float input_states[],int old_cor[])
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
	float max = -FLOATINFTY;	
	//float middle =  (hm->endb - hm->startb)/2 + hm->startb;
	float middle =  (old_cor[3] - old_cor[2])/2 + old_cor[2];
	float sub = 0.0;
	

	prof1+= (stride * (old_cor[4]+1));
	//prof2 += stride * (hm->startb);
	//i = hm->startb;
	prof2 += stride * (old_cor[2]);
	i = old_cor[2];
	c = -1;
	//for(i = hm->startb; i < hm->endb;i++){
	for(i = old_cor[2]; i < old_cor[3];i++){
		sub = abs(middle -i);
		sub /= 1000; 
		prof2 += stride;
		//fprintf(stderr,"%d	%d	%d \n",f[i].a,b[i].a,max);
		if(f[i].a+b[i].a-sub > max){
			max = f[i].a+b[i].a-sub;
	//		fprintf(stderr,"aligned->aligned:%d + %d = %d\n",f[i].a,b[i].a,f[i].a+b[i].a);
			transition = 1;
			c = i;
		}
		if(f[i].a+b[i].ga+prof2[gpo_pos]-sub > max){
			max = f[i].a+b[i].ga+prof2[gpo_pos]-sub;
	//		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
			transition = 2;
			c = i;
		}
		if(f[i].a+b[i].gb+prof1[gpo_pos] -sub> max){
			max = f[i].a+b[i].gb+prof1[gpo_pos]-sub;
	//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
			transition = 3;
			c = i;
		}
		if(f[i].ga+b[i].a+prof2[gpo_pos]-sub > max){
			max = f[i].ga+b[i].a+prof2[gpo_pos]-sub;
	//		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
			transition = 5;
			c = i;
		}


		if(hm->startb == 0){
			if(f[i].gb+b[i].gb+prof1[tgpe_pos]-sub > max){
				max = f[i].gb+b[i].gb+prof1[tgpe_pos]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}else{
			if(f[i].gb+b[i].gb+prof1[gpe_pos]-sub > max){
				max = f[i].gb+b[i].gb+prof1[gpe_pos]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}
		if(f[i].gb+b[i].a+prof1[gpo_pos]-sub > max){
			max = f[i].gb+b[i].a+prof1[gpo_pos]-sub;
	//		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
			transition = 7;
			c = i;
		}
	}
	//i = hm->endb;
	i = old_cor[3];
	sub = abs(middle -i);
	sub /= 1000; 
	if(f[i].a+b[i].gb+prof1[gpo_pos]-sub > max){
		max = f[i].a+b[i].gb+prof1[gpo_pos]-sub;
	//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
		transition = 3;
		c = i;
	}
	if(hm->endb == hm->len_b){
		if(f[i].gb+b[i].gb+prof1[tgpe_pos]-sub > max){
			max = f[i].gb+b[i].gb+prof1[tgpe_pos]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}	
	}else{
		if(f[i].gb+b[i].gb+prof1[gpe_pos]-sub > max){
			max = f[i].gb+b[i].gb+prof1[gpe_pos]-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}
	}
	
	
	
	prof1-= (stride * (old_cor[4]+1));
	//prof2 -= hm->endb << 6;
	prof2 -= old_cor[3] * stride;
	
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
	//		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = feature_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

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
			hirsch_path = feature_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
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
			hirsch_path = feature_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

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
			hirsch_path = feature_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
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
			hirsch_path = feature_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

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
			hirsch_path = feature_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
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
			hirsch_path = feature_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

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
			hirsch_path = feature_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
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
			hirsch_path = feature_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

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
			hirsch_path = feature_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
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
			hirsch_path = feature_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

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
			hirsch_path = feature_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
	}
		
	return hirsch_path;
}

struct states* feature_foward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm)
{
	unsigned int freq[dim];
	struct states* s = hm->f;
	//const int starta = hm->starta;
	//const int enda = hm->enda;
	//const int startb = hm->startb;
	//const int endb = hm->endb;
	
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	register int i = 0;
	register int j = 0;
	register int c = 0;
	
	
	
	prof1 += (hm->starta) * stride;
	prof2 +=  (hm->startb) * stride;
	s[hm->startb].a = s[0].a;
	s[hm->startb].ga = s[0].ga;
	s[hm->startb].gb = s[0].gb;
	if(hm->startb == 0){
		for (j = hm->startb+1; j < hm->endb;j++){
			prof2+=stride;
			s[j].a = -FLOATINFTY;
			
			//s[j].ga = s[j-1].a+prof2[tgpe_pos];
			//if (s[j-1].ga+prof2[tgpe_pos] > s[j].ga){
			//	s[j].ga = s[j-1].ga+prof2[tgpe_pos];
			//}
			if(s[j-1].ga > s[j-1].a){
				s[j].ga = s[j-1].ga+prof2[tgpe_pos];
			}else{
				s[j].ga = s[j-1].a+prof2[tgpe_pos];
			}

			
			s[j].gb = -FLOATINFTY;
		}	
		prof2+=stride;	
	}else{

		for (j = hm->startb+1; j < hm->endb;j++){
			prof2+=stride;
			s[j].a = -FLOATINFTY;
			
			//s[j].ga = s[j-1].a+prof2[gpo_pos];
			//if (s[j-1].ga+prof2[gpe_pos] > s[j].ga){
			//	s[j].ga = s[j-1].ga+prof2[gpe_pos];
			//}

			if(s[j-1].ga+prof2[gpe_pos] > s[j-1].a+prof2[gpo_pos]){
				s[j].ga = s[j-1].ga+prof2[gpe_pos];
			}else{
				s[j].ga = s[j-1].a+prof2[gpo_pos];
			}
			
			s[j].gb = -FLOATINFTY;
		//	prof2+=64;
		}
		prof2+=stride;
	}
	
	prof2 -= (hm->endb-hm->startb) * stride;
	
	s[hm->endb].a = -FLOATINFTY;
	s[hm->endb].ga = -FLOATINFTY;
	s[hm->endb].gb = -FLOATINFTY;


	for (i = hm->starta;i < hm->enda;i++){
		prof1 += stride;
		c = 1;
		for (j = 0;j < dim; j++){
			if(prof1[j]){
				freq[c] = j;
				c++;
			}
		}
		freq[0] = c;
			
		pa = s[hm->startb].a;
		pga = s[hm->startb].ga;
		pgb = s[hm->startb].gb;
		if(hm->startb == 0){
			s[hm->startb].a = -FLOATINFTY;
			s[hm->startb].ga = -FLOATINFTY;
		
			//s[hm->startb].gb = pa+prof1[tgpe_pos];
			//if(pgb+prof1[tgpe_pos] > s[hm->startb].gb){
			//	s[hm->startb].gb = pgb+prof1[tgpe_pos];
			//}
			if(pgb > pa){
				s[hm->startb].gb = pgb+prof1[tgpe_pos];
			}else{
				s[hm->startb].gb = pa+prof1[tgpe_pos];
			}
			
		}else{
			s[hm->startb].a = -FLOATINFTY;
			s[hm->startb].ga = -FLOATINFTY;
		
			//s[hm->startb].gb = pa+prof1[gpo_pos];
			//if(pgb+prof1[gpe_pos] > s[hm->startb].gb){
			//	s[hm->startb].gb = pgb+prof1[gpe_pos];
			//}
			if(pgb+prof1[gpe_pos] > pa+prof1[gpo_pos]){
				s[hm->startb].gb = pgb+prof1[gpe_pos];
			}else{
				s[hm->startb].gb = pa+prof1[gpo_pos];
			}
			
		}
		for (j = hm->startb+1; j <= hm->endb;j++){
			prof2 += stride;
			ca = s[j].a;
			
			/*pga += prof2[-37];
			pga = pa - pga;
			pa = pa -((pga>>31)&pga);
			
			pgb += prof1[-37];
			pa = pa -(((pa - pgb)>>31)&(pa -pgb));*/
			//fprintf(stderr,"%d	%d	%d	%p	%d\n",i,j,gpo_pos-stride,prof2);
			if((pga += prof2[gpo_pos-stride]) > pa){
				pa = pga;
			}

			if((pgb += prof1[gpo_pos-stride]) > pa){
				pa = pgb;
			}
			
			prof2 += dim;
			for (c = 1;c < freq[0];c++){
				pa += prof1[freq[c]]*prof2[freq[c]];
			}
			prof2 -= dim;	

			s[j].a = pa;
			
			pga = s[j].ga;
			
			//s[j].ga = s[j-1].a+prof2[gpo_pos];
			//if (s[j-1].ga+prof2[gpe_pos] > s[j].ga){
			//	s[j].ga = s[j-1].ga+prof2[gpe_pos];
			//}
			if(s[j-1].ga+prof2[gpe_pos] > s[j-1].a+prof2[gpo_pos]){
				s[j].ga = s[j-1].ga+prof2[gpe_pos];
			}else{
				s[j].ga = s[j-1].a+prof2[gpo_pos];
			}
			
			
			pgb = s[j].gb;
			
			//s[j].gb = ca+prof1[gpo_pos];
			//if(pgb+prof1[gpe_pos] > s[j].gb){
			//	s[j].gb = pgb+prof1[gpe_pos];
			//}
			if(pgb+prof1[gpe_pos] > ca+prof1[gpo_pos]){
				s[j].gb = pgb+prof1[gpe_pos];
			}else{
				s[j].gb = ca+prof1[gpo_pos];
			}
			
			pa = ca;
		}
		prof2 -= (hm->endb-hm->startb) * stride;
		
	}
	prof1 -= stride * (hm->enda);
	return s;
}

struct states* feature_backward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm)
{
	unsigned int freq[dim];
	struct states* s = hm->b;
	//const int starta = hm->starta;
	//const int enda = hm->enda;
	//const int startb = hm->startb;
	//const int endb = hm->endb;
	
	register float pa = 0;
	register float pga = 0;
	register float pgb = 0;
	register float ca = 0;
	register int i = 0;
	register int j = 0;
	register int c = 0;

	prof1 += (hm->enda+1) * stride;
	prof2 += (hm->endb+1) * stride;
	s[hm->endb].a = s[0].a;
	s[hm->endb].ga = s[0].ga;
	s[hm->endb].gb = s[0].gb;
	
	
	//init of first row;
	//j = endb-startb;
	if(hm->endb == hm->len_b){
		
		for(j = hm->endb-1;j > hm->startb;j--){
			prof2 -= stride;
			s[j].a = -FLOATINFTY;
			
			//s[j].ga = s[j+1].a+prof2[tgpe_pos];
			//if (s[j+1].ga+prof2[tgpe_pos] > s[j].ga){
			//	s[j].ga = s[j+1].ga+prof2[tgpe_pos];
			//}
			if(s[j+1].ga > s[j+1].a){
				s[j].ga = s[j+1].ga+prof2[tgpe_pos];
			}else{
				s[j].ga = s[j+1].a+prof2[tgpe_pos];
			}
			
			s[j].gb = -FLOATINFTY;
		}
		prof2 -= stride;
	}else{
		for(j = hm->endb-1;j > hm->startb;j--){
			prof2 -= stride;
			s[j].a = -FLOATINFTY;
			
			//s[j].ga = s[j+1].a+prof2[gpo_pos];
			//if (s[j+1].ga+prof2[gpe_pos] > s[j].ga){
			//	s[j].ga = s[j+1].ga+prof2[gpe_pos];
			//}
			if(s[j+1].ga+prof2[gpe_pos] > s[j+1].a+prof2[gpo_pos]){
				s[j].ga = s[j+1].ga+prof2[gpe_pos];
			}else{
				s[j].ga = s[j+1].a+prof2[gpo_pos];
			}
			s[j].gb = -FLOATINFTY;
		//	prof2 -= 64;
		}
		prof2 -= stride;
	}
	
	s[hm->startb].a = -FLOATINFTY;
	s[hm->startb].ga = -FLOATINFTY;
	s[hm->startb].gb = -FLOATINFTY;
//	prof2 -= (endb -startb) << 6;

	i = hm->enda-hm->starta;
	while(i--){
		prof1 -= stride;

		c = 1;
		for (j = 0;j < dim; j++){
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

		if(hm->endb == hm->len_b){
			//s[hm->endb].gb = pa+prof1[tgpe_pos];
			//if(pgb+prof1[tgpe_pos] > s[hm->endb].gb){
			//	s[hm->endb].gb = pgb+prof1[tgpe_pos];
			//}
			if(pgb > pa){
				s[hm->endb].gb = pgb+prof1[tgpe_pos];
			}else{
				s[hm->endb].gb = pa+prof1[tgpe_pos];
			}
		}else{
			//s[hm->endb].gb = pa+prof1[gpo_pos];
			//if(pgb+prof1[gpe_pos] > s[hm->endb].gb){
			//	s[hm->endb].gb = pgb+prof1[gpe_pos];
			//}
			if(pgb+prof1[gpe_pos] > pa+prof1[gpo_pos]){
				s[hm->endb].gb = pgb+prof1[gpe_pos];
			}else{
				s[hm->endb].gb = pa+prof1[gpo_pos];
			}
		}
		//j = endb-startb;
		prof2 += (hm->endb-hm->startb) * stride;
		//while(j--){
		for(j = hm->endb-1;j >= hm->startb;j--){
			prof2 -= stride;
			ca = s[j].a;
			if((pga += prof2[stride+ gpo_pos]) > pa){
				pa = pga;
			}
			if((pgb += prof1[stride+gpo_pos]) > pa){
				pa = pgb;
			}
			
			prof2 += dim;
			for (c = 1;c < freq[0];c++){
				pa += prof1[freq[c]]*prof2[freq[c]];
			}
			prof2 -= dim;

			s[j].a = pa;
			
			pga = s[j].ga;
			
			//s[j].ga = s[j+1].a+prof2[gpo_pos];
			//if (s[j+1].ga+prof2[gpe_pos] > s[j].ga){
			//	s[j].ga = s[j+1].ga+prof2[gpe_pos];
			//}
			if(s[j+1].ga+prof2[gpe_pos] > s[j+1].a+prof2[gpo_pos]){
				s[j].ga = s[j+1].ga+prof2[gpe_pos];
			}else{
				s[j].ga = s[j+1].a+prof2[gpo_pos];
			}
			
			pgb = s[j].gb;
			
			//s[j].gb = ca+prof1[gpo_pos];
			//if(pgb+prof1[gpe_pos] > s[j].gb){
			//	s[j].gb = pgb+prof1[gpe_pos];
			//}
			if(pgb+prof1[gpe_pos] > ca+prof1[gpo_pos]){
				s[j].gb = pgb+prof1[gpe_pos];
			}else{
				s[j].gb = ca+prof1[gpo_pos];
			}
			pa = ca;
		}
	}		
	return s;
}




struct feature_matrix* get_feature_matrix(struct feature_matrix* fm, struct alignment* aln,struct parameters* param)
{
	struct utype_ufeat* utf = 0;
	struct feature* n = 0;
	struct feature* p = 0;
	
	int i = 0;
	int j = 0;
	char* requested_feature = param->feature_type;

	
	utf = get_unique_features(aln,utf);
	
	fm = malloc(sizeof(struct feature_matrix));
	
	
	
	
	if (byg_start(requested_feature,"allALL")!= -1){
		n = utf->f;
		i = 0;
		while(n){
			n->color = i;
			i++;
			n = n->next;
		}
	}else if(byg_start(requested_feature,"maxplpMAXPLP")!= -1){
		n = utf->f;
		i = 0;
		while(n){
			if(byg_start("SIGNAL PEPTIDE",n->note)!= -1){
				n->color = 0;
			}
			if(byg_start("TRANSMEMBRANE",n->note)!= -1){
				n->color = 1;
			}
			if(byg_start("TRANSLOCATED LOOP",n->note)!= -1){
				n->color = 2;
			}
			if(byg_start("CYTOPLASMIC LOOP",n->note)!= -1){
				n->color = 3;
			}
			n = n->next;
		}
		i = 4;
	}else{
		n = utf->f;
		i = 0;
		while(n){
			if(check_identity(requested_feature,n->type)!= -1){
				//fprintf(stderr,"%s	%s\n",requested_feature,n->type);
				n->color = i;
				i++;
			}else{
				n->color = -1;
			}
			n = n->next;
		}
	}
	/*if(!i){
		fprintf(stderr,"WARNING: no feature of type '%s' was found in the input file.\n",requested_feature);
		fprintf(stderr,"	\nAvailable features are:\n\n");
		fprintf(stderr,"	Type		Feature\n");
		fprintf(stderr,"	-----------------------------\n");
		n = utf->f;
		while(n){
			fprintf(stderr,"	%s		%s\n",n->type,n->note);
			n = n->next;
		}
		
		free_utf(utf);
		free_aln(aln);
		free(fm);
		return 0;
	}*/
	
	if(byg_start(requested_feature,"maxplp")!= -1){
		fm->mdim = 4;
		fm->stride = fm->mdim << 1;
		fm->m = malloc(sizeof(float*)*fm->mdim);
		for (i = 0;i < fm->mdim;i++){
			fm->m[i] = malloc(sizeof(float)*fm->mdim);
		}
	/*1: 0.60 0.20 0.15 0.05
	2: 0.20 0.60 0.10 0.10
	3: 0.15 0.10 0.50 0.25
	4: 0.05 0.10 0.25 0.60*/

		fm->m[0][0] = 60;
		fm->m[0][1] = 20;
		fm->m[0][2] = 15;
		fm->m[0][3] = 5;
		fm->m[1][0] = 20;
		fm->m[1][1] = 60;
		fm->m[1][2] = 10;
		fm->m[1][3] = 10;
		fm->m[2][0] = 15;
		fm->m[2][1] = 10;
		fm->m[2][2] = 50;
		fm->m[2][3] = 25;
		fm->m[3][0] = 5;
		fm->m[3][1] = 10;
		fm->m[3][2] = 25;
		fm->m[3][3] = 60;
		
	}else if(byg_start(requested_feature,"wumanber")!= -1){
		fm->mdim = i;
		fm->stride = fm->mdim << 1;
		fm->m = malloc(sizeof(float*)*fm->mdim);
		for (i = 0;i < fm->mdim;i++){
			fm->m[i] = malloc(sizeof(float)*fm->mdim);
			for (j = 0;j < fm->mdim;j++){
				fm->m[i][j] = 0;
			}
		}
		for (i = 0;i < fm->mdim;i++){
			fm->m[i][i] = 100;
		}	
	//	fprintf(stderr,"WU	%d	\n",fm->mdim);
	}else{
		fm->mdim = i;
		fm->stride = fm->mdim << 1;
		fm->m = malloc(sizeof(float*)*fm->mdim);
		for (i = 0;i < fm->mdim;i++){
			fm->m[i] = malloc(sizeof(float)*fm->mdim);
			for (j = 0;j < fm->mdim;j++){
				fm->m[i][j] = param->diff_feature_score;
			}
		}
		for (i = 0;i < fm->mdim;i++){
			fm->m[i][i] = param->same_feature_score;
		}	
		/*for (i = 0;i < fm->mdim;i++){
			for (j = 0;j < fm->mdim;j++){
				fprintf(stderr,"%f ",fm->m[i][j]);
			}
			fprintf(stderr,"\n");
		}*/

	}
	//float fr = 0.0;
	for (i = numseq;i--;){
                n = aln->ft[i];
	//	fprintf(stderr,"SEQUENCE %d\n",i);
                while(n){
                        p = utf->f;
                        while(p){
                        	if(check_identity(requested_feature,n->type)!= -1){
                                if(check_identity(n->note,p->note)!= -1){
					n->color = p->color;
             //                           fr += n->end - n->start+1;
            //                           fprintf(stderr,"SEQ:%d	FEATURE FOUND:%s	%s	%d-%d	 color:%d	\n",i,n->note,p->note,n->start,n->end,n->color);
                                        break;
                                }
                                }
                                p = p->next;
                        }
                        n = n->next;
                }
        }
//	fprintf(stderr,"%f\n",fr);
	
	//float res = 0.0;
	
	//for (i = 0; i < numseq;i++){
	//	res += aln->sl[i];
	//}
	//fprintf(stdout,"%f	%f	%f\n",fr,res,fr/res);
	//exit(0);
	
	/*
	n = utf->t;
	fprintf(stderr,"TYPES:	we use:%d\n",i);
	while(n){
		fprintf(stderr,"%s\n",n->type);
		n = n->next;
	}
	
	n = utf->f;
	fprintf(stderr,"Features:\n");
	i = 0;
	while(n){
		fprintf(stderr,"%d:	%s:%s	col:%d\n",i,n->type,n->note,n->color);
		i++;
		n = n->next;
	}
	fprintf(stderr,"REQUESTED FEATURE:%s\n",requested_feature);

	for (i = 0;i < fm->mdim;i++){
		for (j = 0;j < fm->mdim;j++){
			fprintf(stderr,"%d ",fm->m[i][j]);
		}
		fprintf(stderr,"\n");
	}	
	fprintf(stderr,"\n");*/
	
	
	free_utf(utf);
	return fm; 
}

struct utype_ufeat* get_unique_features(struct alignment* aln,struct utype_ufeat* utf)
{
	int i;
	utf = malloc(sizeof(struct utype_ufeat)*1);
	utf->t = 0;
	utf->f = 0;
	for (i =0; i < numseq;i++){
		utf = traverse_ft(utf,aln->ft[i]);
	}
	return utf;
}


struct utype_ufeat* traverse_ft(struct utype_ufeat* utf,struct feature* n)
{
	if (n != NULL){
		utf->t = add_unique_type(utf->t,n);
		utf->f = add_unique_feature(utf->f,n);
		traverse_ft(utf,n->next);
	}
	return utf;
}


struct feature* add_unique_feature(struct feature *n, struct feature *toadd)
{
	int i;
		
        if (n == NULL){
		n = (struct feature*) malloc(sizeof(struct feature));
		n->type = malloc(sizeof(char)* (strlen(toadd->type)+1));
		for ( i= 0;i < strlen(toadd->type);i++){
			n->type[i] = toadd->type[i];
		}
		n->type[i] = 0;
		
		n->note = malloc(sizeof(char)* (strlen(toadd->note)+1));
		for ( i= 0;i < strlen(toadd->note);i++){
			n->note[i] = toadd->note[i];
		}
		n->note[i] = 0;
		
		n->start = toadd->end - toadd->start;
		n->end = 0;
		
		n->next = 0;
	}else{
	
		if((check_identity(toadd->note,n->note)== -1)){
			n->next = add_unique_feature(n->next,toadd);
		}else{
			n->start += toadd->end - toadd->start;
		}
	}
        return n;
}

struct feature* add_unique_type(struct feature *n, struct feature *toadd)
{
	int i;
		
        if (n == NULL){
		n = (struct feature*) malloc(sizeof(struct feature));
		n->type = malloc(sizeof(char)* (strlen(toadd->type)+1));
		for ( i= 0;i < strlen(toadd->type);i++){
			n->type[i] = toadd->type[i];
		}
		n->type[i] = 0;
		
		n->note = malloc(sizeof(char)* (strlen(toadd->note)+1));
		for ( i= 0;i < strlen(toadd->note);i++){
			n->note[i] = toadd->note[i];
		}
		n->note[i] = 0;
		
		n->start = 0;
		n->end = 0;
		
		n->next = 0;
	}else{
		if((check_identity(toadd->type,n->type)== -1)){
			n->next = add_unique_type(n->next,toadd);
		}
	}
        return n;
}






