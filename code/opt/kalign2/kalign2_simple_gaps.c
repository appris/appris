/*
	kalign2_simple_gaps.c
	
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

int** simple_hirschberg_alignment(struct alignment* aln,int* tree,int**submatrix, int** map)
{
	struct hirsch_mem* hm = 0;
	int i,j,g,a,b,c;
	int len_a;
	int len_b;
	int** profile = 0;

	profile = malloc(sizeof(int*)*numprofiles);
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
	//	fprintf(stderr,"Aligning:%d %d->%d	done:%f\n",a,b,c,((float)(i+1)/(float)numseq)*100);
		len_a = aln->sl[a];
		len_b = aln->sl[b];

		
		g = (len_a > len_b)? len_a:len_b;
		map[c] = malloc(sizeof(int) * (g+2));
		if(g > hm->size){
			hm = hirsch_mem_realloc(hm,g);
		}

		for (j = 0; j < (g+2);j++){
		//	hirsch_path[j] = -1;
			map[c][j] = -1;
		//	map[c][j] = 0;
		}
	//	map[c][0] = len_a;
		//map[c][len_a+len_b+1] = 3;

		if (a < numseq){
			profile[a] = simple_make_profile(profile[a],aln->s[a],len_a,submatrix);
		}
		if (b < numseq){
			profile[b] = simple_make_profile(profile[b],aln->s[b],len_b,submatrix);
		}

		hm->starta = 0;
		hm->startb = 0;
		hm->enda = len_a;
		hm->endb = len_b;
		hm->len_a = len_a;
		hm->len_b = len_b;
		
		hm->f[0].a = 0;
		hm->f[0].ga =  -INFTY;
		hm->f[0].gb = -INFTY;
		hm->b[0].a = 0;
		hm->b[0].ga =  -INFTY;
		hm->b[0].gb =  -INFTY;
	//	fprintf(stderr,"LENA:%d	LENB:%d	numseq:%d\n",len_a,len_b,numseq);
		if(len_a < len_b){
			map[c] = simple_hirsch_pp_dyn(profile[a],profile[b],hm,map[c]);
		}else{
			hm->enda = len_b;
			hm->endb = len_a;
			hm->len_a = len_b;
			hm->len_b = len_a;
			map[c] = simple_hirsch_pp_dyn(profile[b],profile[a],hm,map[c]);
			map[c] = mirror_hirsch_path(map[c],len_a,len_b);
		}
					
		map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);

		if(i != numseq-2){
			profile[c] = malloc(sizeof(int)*64*(map[c][0]+2));
			profile[c] = simple_update(profile[a],profile[b],profile[c],map[c]);
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


int* simple_make_profile(int* prof, int* seq,int len,int** subm)
{
	int i,j,c;	
	prof = malloc(sizeof(int)*(len+2)*64);
	prof +=  (64 *(len+1));

	for (i = 0;i < 64;i++){
		prof[i] = 0;
	}
	prof[23+32] = -gpo;
	prof[26] = 1;

	
	i = len;
	while(i--){
		prof -= 64;

		for (j = 0;j < 64;j++){
			prof[j] = 0;
		}
		prof[26] = 1;//number of residues // both additive
		c = seq[i];
		
		prof[c] += 1;
		
		prof += 32;
		
		for(j = 23;j--;){
			prof[j] = subm[c][j];
		}
		prof[23] = -gpo;
		
		prof -= 32;
	}
	prof -= 64;
	for (i = 0;i < 64;i++){
		prof[i] = 0;
	}

	prof[23+32] = -gpo;
	prof[26] = 1;
	return prof;
}




int* simple_update(int* profa,int* profb,int* newp,int* path)
{
	int i,c;
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
			for (i = 64; i--;){
				newp[i] = profa[i] + profb[i];
			}
				
			
			profa += 64;
			profb += 64;
		}
		if (path[c] & 1){
			for (i = 64; i--;){
				newp[i] = profb[i];
			}
			profb += 64;			
		}
		if (path[c] & 2){
			for (i = 64; i--;){
				newp[i] = profa[i];
			}
			profa+=64;			
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



int* simple_hirsch_pp_dyn(const int* prof1,const int* prof2,struct hirsch_mem* hm, int* hirsch_path)
{
	int mid = ((hm->enda - hm->starta) / 2)+ hm->starta;
	int input_states[6] = {hm->f[0].a,hm->f[0].ga,hm->f[0].gb,hm->b[0].a,hm->b[0].ga,hm->b[0].gb};
	int old_cor[5] = {hm->starta,hm->enda,hm->startb,hm->endb,mid};

	
	//fprintf(stderr,"starta:%d enda:%d startb:%d endb:%d mid:%d\n",hm->starta,hm->enda,hm->startb,hm->endb,mid);
	
	
	if(hm->starta  >= hm->enda){
		return hirsch_path;
	}
	if(hm->startb  >= hm->endb){
		return hirsch_path;
	}

	hm->enda = mid;
	hm->f = simple_foward_hirsch_pp_dyn(prof1,prof2,hm);
	/*int i;
	fprintf(stderr,"FOWARD\n");
	for (i = hm->startb; i <= hm->endb;i++){
		fprintf(stderr,"%d	%d	%d\n",hm->f[i].a,hm->f[i].ga,hm->f[i].gb);
	}*/

	hm->starta = mid;
	hm->enda = old_cor[1];
	hm->b = simple_backward_hirsch_pp_dyn(prof1,prof2,hm);
	/*fprintf(stderr,"BaCKWARD\n");

	for (i = hm->startb; i <= hm->endb;i++){
		fprintf(stderr,"%d	%d	%d\n",hm->b[i].a,hm->b[i].ga,hm->b[i].gb);
	}*/

	hirsch_path = simple_hirsch_align_two_pp_vector(prof1,prof2,hm,hirsch_path,input_states,old_cor);
	return hirsch_path;
}



int* simple_hirsch_align_two_pp_vector(const int* prof1,const int* prof2,struct hirsch_mem* hm,int* hirsch_path,int input_states[],int old_cor[])
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
	

	prof1+= (64 * (old_cor[4]+1));
	prof2 += 64 * (hm->startb);
	i = hm->startb;
	c = -1;
	for(i = hm->startb; i < hm->endb;i++){
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
		if(f[i].a+b[i].ga+prof2[23]*prof1[26]-sub > max){
			max = f[i].a+b[i].ga+prof2[23]*prof1[26]-sub;
	//		fprintf(stderr,"aligned->gap_a:%d + %d +%d = %d\n",f[i].a,b[i].ga,prof1[27],f[i].a+b[i].ga+prof2[27]);
			transition = 2;
			c = i;
		}
		if(f[i].a+b[i].gb+prof1[23]*prof2[26] -sub> max){
			max = f[i].a+b[i].gb+prof1[23]*prof2[26]-sub;
	//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
			transition = 3;
			c = i;
		}
		if(f[i].ga+b[i].a+prof2[23]*prof1[26]-sub > max){
			max = f[i].ga+b[i].a+prof2[23]*prof1[26]-sub;
	//		fprintf(stderr,"gap_a->aligned:%d + %d + %d(gpo) = %d\n",f[i].ga,b[i].a,prof2[27],f[i].ga+b[i].a+prof2[27]);
			transition = 5;
			c = i;
		}


		if(hm->startb == 0){
			if(f[i].gb+b[i].gb-sub > max){
				max = f[i].gb+b[i].gb-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}else{
			if(f[i].gb+b[i].gb-sub > max){
				max = f[i].gb+b[i].gb-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
				transition = 6;
				c = i;
			}
		}
		if(f[i].gb+b[i].a+prof1[23]*prof2[26]-sub > max){
			max = f[i].gb+b[i].a+prof1[23]*prof2[26]-sub;
	//		fprintf(stderr,"gap_b->aligned:%d + %d + %d(gpo) = %d\n",f[i].gb,b[i].a,prof1[27],f[i].gb+b[i].a+prof1[27]);
			transition = 7;
			c = i;
		}
	}
	i = hm->endb;
	sub = abs(middle -i);
	sub /= 1000; 
	if(f[i].a+b[i].gb+prof1[23]*prof2[26]-sub > max){
		max = f[i].a+b[i].gb+prof1[23]*prof2[26]-sub;
	//		fprintf(stderr,"aligned->gap_b:%d + %d +%d = %d\n",f[i].a,b[i].gb,prof1[27],f[i].a+b[i].gb+prof1[27]);
		transition = 3;
		c = i;
	}
	if(hm->endb == hm->len_b){
		if(f[i].gb+b[i].gb-sub > max){
			max = f[i].gb+b[i].gb-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}	
	}else{
		if(f[i].gb+b[i].gb-sub > max){
			max = f[i].gb+b[i].gb-sub;
	//			fprintf(stderr,"gap_b->gap_b:%d + %d +%d(gpe) =%d \n",f[i].gb, b[i].gb, prof1[28],f[i].gb+b[i].gb+prof1[28]);
			transition = 6;
			c = i;
		}
	}
	
	
	
	prof1-= (64 * (old_cor[4]+1));
	prof2 -= hm->endb << 6;
	
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
			hm->b[0].ga = -INFTY;
			hm->b[0].gb = -INFTY;
	//		fprintf(stderr,"Using this for start:%d	%d	%d\n",hm->f[0].a,hm->f[0].ga,hm->f[0].gb);
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = simple_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0;
			hm->f[0].ga = -INFTY;
			hm->f[0].gb = -INFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = simple_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
		case 2:// a -> ga = 2
			
			hirsch_path[old_cor[4]] = c;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0;
			hm->b[0].ga = -INFTY;
			hm->b[0].gb = -INFTY;
			
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = simple_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4];
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = -INFTY;
			hm->f[0].ga = 0;
			hm->f[0].gb = -INFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d  what:%d-%d	%d-%d\n",c+1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = simple_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
		case 3:// a -> gb = 3
			
			hirsch_path[old_cor[4]] = c;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4],c);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = 0;
			hm->b[0].ga = -INFTY;
			hm->b[0].gb = -INFTY;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = simple_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -INFTY;
			hm->f[0].ga = -INFTY;
			hm->f[0].gb = 0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = simple_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
		case 5://ga -> a = 5
			hirsch_path[old_cor[4]+1] = c+1;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);

			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -INFTY;
			hm->b[0].ga = 0;
			hm->b[0].gb = -INFTY;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4];
			
			hm->startb = old_cor[2];
			hm->endb = c-1;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = simple_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0;
			hm->f[0].ga = -INFTY;
			hm->f[0].gb = -INFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = simple_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
		case 6://gb->gb = 6;
			
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -INFTY;
			hm->b[0].ga = -INFTY;
			hm->b[0].gb = 0;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = simple_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c;
			hm->endb = old_cor[3];
			hm->f[0].a = -INFTY;
			hm->f[0].ga = -INFTY;
			hm->f[0].gb = 0;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = simple_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
		case 7://gb->a = 7;
			
			hirsch_path[old_cor[4]+1] = c+1;
	//		fprintf(stderr,"Aligning:%d-%d\n",old_cor[4]+1,c+1);
			//foward:
			hm->f[0].a = input_states[0];
			hm->f[0].ga = input_states[1];
			hm->f[0].gb = input_states[2];
			hm->b[0].a = -INFTY;
			hm->b[0].ga = -INFTY;
			hm->b[0].gb = 0;
			
			hm->starta = old_cor[0];
			hm->enda = old_cor[4]-1;
			hm->startb = old_cor[2];
			hm->endb = c;
			//fprintf(stderr,"Following first: %d  what:%d-%d	%d-%d\n",c-1,hm->starta,hm->enda,hm->startb,hm->endb);
			hirsch_path = simple_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);

			//backward:
			hm->starta = old_cor[4]+1;
			hm->enda = old_cor[1];
			hm->startb = c+1;
			hm->endb = old_cor[3];
			hm->f[0].a = 0;
			hm->f[0].ga = -INFTY;
			hm->f[0].gb = -INFTY;
			hm->b[0].a = input_states[3];
			hm->b[0].ga = input_states[4];
			hm->b[0].gb = input_states[5];
	
			//fprintf(stderr,"Following last: %d\n",c+1);
			hirsch_path = simple_hirsch_pp_dyn(prof1,prof2,hm,hirsch_path);
			break;
	}
		
	return hirsch_path;
}



struct states* simple_foward_hirsch_pp_dyn(const int* prof1,const int* prof2,struct hirsch_mem* hm)
{
	unsigned int freq[23];
	struct states* s = hm->f;

	register int pa = 0;
	register int pga = 0;
	register int pgb = 0;
	register int ca = 0;
	register int i = 0;
	register int j = 0;
	
	
	
	prof1 += (hm->starta) << 6;
	prof2 +=  (hm->startb) << 6;
	s[hm->startb].a = s[0].a;
	s[hm->startb].ga = s[0].ga;
	s[hm->startb].gb = s[0].gb;
	if(hm->startb == 0){
		for (j = hm->startb+1; j < hm->endb;j++){
			prof2+=64;
			s[j].a = -INFTY;
			
			s[j].ga = s[j-1].a;
			if (s[j-1].ga > s[j].ga){
				s[j].ga = s[j-1].ga;
			}
			s[j].gb = -INFTY;
		}	
		prof2+=64;	
	}else{

		for (j = hm->startb+1; j < hm->endb;j++){
			prof2+=64;
			s[j].a = -INFTY;
			
			s[j].ga = s[j-1].a+prof2[23]*prof1[26];
			if (s[j-1].ga > s[j].ga){
				s[j].ga = s[j-1].ga;
			}
			s[j].gb = -INFTY;
		//	prof2+=64;
		}
		prof2+=64;
	}
	
	prof2 -= (hm->endb-hm->startb) << 6;
	
	s[hm->endb].a = -INFTY;
	s[hm->endb].ga = -INFTY;
	s[hm->endb].gb = -INFTY;


	for (i = hm->starta;i < hm->enda;i++){
		prof1 += 64;
		pa = 1;
 		for (j = 23; j--;){
			if(prof1[j]){
				freq[pa] = j;
				pa++;	
			}
		}
		freq[0] = pa;
			
		pa = s[hm->startb].a;
		pga = s[hm->startb].ga;
		pgb = s[hm->startb].gb;
		if(hm->startb == 0){
			s[hm->startb].a = -INFTY;
			s[hm->startb].ga = -INFTY;
		
			s[hm->startb].gb = pa;
				if(pgb> s[hm->startb].gb){
				s[hm->startb].gb = pgb;
			}
		}else{
			s[hm->startb].a = -INFTY;
			s[hm->startb].ga = -INFTY;
		
			s[hm->startb].gb = pa+prof1[23]*prof2[26];
			if(pgb > s[hm->startb].gb){
				s[hm->startb].gb = pgb;
			}
		}
		for (j = hm->startb+1; j <= hm->endb;j++){
			prof2 += 64;
			ca = s[j].a;
			
			if((pga += prof2[23-64]*prof1[26-64]) > pa){ 
				pa = pga;
			}

			if((pgb += prof1[23-64]*prof2[26-64]) > pa){
				pa = pgb;
			}
			
			prof2 += 32;
			for (pga = freq[0];--pga;){
				pgb = freq[pga];
				pa += prof1[pgb]*prof2[pgb];
			}
			prof2 -= 32;	

			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = s[j-1].a+prof2[23]*prof1[26];
			if (s[j-1].ga> s[j].ga){
				s[j].ga = s[j-1].ga;
			}
			
			pgb = s[j].gb;
			
			s[j].gb = ca+prof1[23]*prof2[26];
			if(pgb > s[j].gb){
				s[j].gb = pgb;
			}
			pa = ca;
		}
		prof2 -= (hm->endb-hm->startb) << 6;
		
	}
	prof1 -= 64 * (hm->enda);
	return s;
}

struct states* simple_backward_hirsch_pp_dyn(const int* prof1,const int* prof2,struct hirsch_mem* hm)
{
	unsigned int freq[23];
	struct states* s = hm->b;
	
	register int pa = 0;
	register int pga = 0;
	register int pgb = 0;
	register int ca = 0;
	register int i = 0;
	register int j = 0;

	prof1 += (hm->enda+1) << 6;
	prof2 += (hm->endb+1) << 6;
	s[hm->endb].a = s[0].a;
	s[hm->endb].ga = s[0].ga;
	s[hm->endb].gb = s[0].gb;
	
	
	//init of first row;
	//j = endb-startb;
	if(hm->endb == hm->len_b){
		
		for(j = hm->endb-1;j > hm->startb;j--){
			prof2 -= 64;
			s[j].a = -INFTY;
			
			s[j].ga = s[j+1].a;
			if (s[j+1].ga> s[j].ga){
				s[j].ga = s[j+1].ga;
			}
			s[j].gb = -INFTY;
		}
		prof2 -= 64;
	}else{
		for(j = hm->endb-1;j > hm->startb;j--){
			prof2 -= 64;
			s[j].a = -INFTY;
			
			s[j].ga = s[j+1].a+prof2[23]*prof1[26];
			if (s[j+1].ga > s[j].ga){
				s[j].ga = s[j+1].ga;
			}
			s[j].gb = -INFTY;
		//	prof2 -= 64;
		}
		prof2 -= 64;
	}
	
	s[hm->startb].a = -INFTY;
	s[hm->startb].ga = -INFTY;
	s[hm->startb].gb = -INFTY;
//	prof2 -= (endb -startb) << 6;

	i = hm->enda-hm->starta;
	while(i--){
		prof1 -= 64;

		pa = 1;
		for (j = 23; j--;){
			if(prof1[j]){
				freq[pa] = j;
				pa++;	
			}
		}
		freq[0] = pa;
		
		pa = s[hm->endb].a;
		pga = s[hm->endb].ga;
		pgb = s[hm->endb].gb;
		s[hm->endb].a = -INFTY;
		s[hm->endb].ga = -INFTY;

		if(hm->endb == hm->len_b){
			s[hm->endb].gb = pa;
			if(pgb > s[hm->endb].gb){
				s[hm->endb].gb = pgb;
			}
		}else{
			s[hm->endb].gb = pa+prof1[23]*prof2[26];
			if(pgb> s[hm->endb].gb){
				s[hm->endb].gb = pgb;
			}
		}
		//j = endb-startb;
		prof2 += (hm->endb-hm->startb) << 6;
		//while(j--){
		for(j = hm->endb-1;j >= hm->startb;j--){
			prof2 -= 64;
			ca = s[j].a;
			if((pga += prof2[64+23]*prof1[26]) > pa){
				pa = pga;
			}
			if((pgb += prof1[64+23]*prof2[26]) > pa){
				pa = pgb;
			}
			
			prof2 += 32;
			for (pga = freq[0];--pga;){
				pgb = freq[pga];
				pa += prof1[pgb]*prof2[pgb];
			}
			prof2 -= 32;

			s[j].a = pa;
			
			pga = s[j].ga;
			
			s[j].ga = s[j+1].a+prof2[23]*prof1[26];
			if (s[j+1].ga > s[j].ga){
				s[j].ga = s[j+1].ga;
			}
			
			pgb = s[j].gb;
			
			s[j].gb = ca+prof1[23]*prof2[26];
			if(pgb > s[j].gb){
				s[j].gb = pgb;
			}
			pa = ca;
		}
	}		
	return s;
}


