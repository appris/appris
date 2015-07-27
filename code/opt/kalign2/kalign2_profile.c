/*
	kalign2_profile.c 
	
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


/*
void add_feature_information_from_alignment(int* path,int* fprof1,int* fprof2,int weight)
{
 	int i = 0;
 	int j = 0;
 	int c = 1;
	while(path[c] != 3){		
		if (!path[c]){
			fprof1[i] +=1;
			fprof1[i+1] +=weight;
			fprof2[j] +=1;
			fprof2[j+1] +=weight;
			i+=2;
			j+=2;	
		}
		if (path[c] & 1){
			j+=2;
		}
		if (path[c] & 2){
			i+=2;		
		}
		c++;
	}
	free(path);
}*/


float* update2(const float* profa, const float* profb,float* newp,int* path,int sipa,int sipb,float internal_gap_weight)
{
	int i,c;
	int* gap_len = 0;
	int gap_cost = 0;
	
	gap_len = malloc(sizeof(int)* (path[0]+1));
	gap_len[0] = 0;
	
	//fprintf(stderr,"%d len,,,,\n",path[0]);
	for(i = 1; i <= path[0];i++){
	//	fprintf(stderr,"%d,%d	",i,path[i]);
		gap_len[i] = (path[i] >> 16);
		path[i]  = path[i] & 0x0000ffff;
	//	fprintf(stderr,"%d	%d\n",path[i],gap_len[i]);
	}
	//gap_len[path[0]] = 0;
//	int len = 0;
	c = 1;
	
	/*while(path[c] != 3){	
		fprintf(stderr,"%d %d	%d\n",c,path[c],gap_len[c]);
				
		c++;
	}
	exit(0);*/


	while(path[c] != 3){		
		gap_cost = 0;
		if (!path[c]){
			while(!path[c] && path[c] != 3){
		//		fprintf(stderr,"Align	%d	%d\n",c,path[c]);
				for (i = 64; i--;){
					newp[i] = profa[i] + profb[i];
				}
				profa += 64;
				profb += 64;
				newp += 64; 
				c++;
			}
		}else if (path[c] & 1){
			//fprintf(stderr,"%d\n",gap_len[c]);
			if(path[c]  & 128){//N terminal gap !!!!!!1
				for (i = 0; i < gap_len[c]-1;i++){
					gap_cost += profb[29+64*i];
			//		fprintf(stderr,"i:%d	%d\n",i,gap_cost);
				}				
				gap_cost += profb[27+64*i];
			//	fprintf(stderr,"i:%d	%d\n",i,gap_cost);
			}else if(path[c]  & 64){//c terminal gap !!!!!!1
			//	fprintf(stderr,"c terminal gap\n");
				gap_cost += profb[27+64];
			//	fprintf(stderr,"i:%d	%d\n",0,gap_cost);
				for (i = 1; i < gap_len[c];i++){
					gap_cost += profb[29+64*i];
			//		fprintf(stderr,"i:%d	%d\n",i,gap_cost);
				}
			}else{
			//	fprintf(stderr,"middle gap\n");
				gap_cost += profb[27+64];
			//	fprintf(stderr,"i:%d	%d\n",0,gap_cost);
				for (i = 1; i < gap_len[c]-1;i++){
					gap_cost += profb[28+64*i];
			//		fprintf(stderr,"i:%d	%d\n",i,gap_cost);
				}
				gap_cost += profb[27+64*i];
			//	fprintf(stderr,"i:%d	%d\n",i,gap_cost);
			}
			//fprintf(stderr,"gap_A	%d	%d	length:%d	cost:%d\n",c,path[c],gap_len[c],gap_cost);
			gap_cost /=  gap_len[c];
			gap_cost *= internal_gap_weight;
			
			while(path[c] & 1 && path[c] != 3){
		//		fprintf(stderr,"gap_A	%d	%d	cost:%d\n",c,path[c],gap_cost);
				for (i = 64; i--;){
					newp[i] = profb[i];
				}
				newp[23] += gap_cost;
				for (i = 32; i < 55;i++){
					newp[i] += gap_cost;
				}
				profb +=64;
				newp += 64;
				c++;
			}		
		}else if (path[c] & 2){
			//fprintf(stderr,"%d\n",gap_len[c]);
			if(path[c]  & 128){//N terminal gap !!!!!!1
				for (i = 0; i < gap_len[c]-1;i++){
					gap_cost += profa[29+64*i];
			//		fprintf(stderr,"i:%d	%d\n",i,gap_cost);
				}
				gap_cost += profa[27+64*i];
			//	fprintf(stderr,"i:%d	%d\n",i,gap_cost);
			}else if(path[c]  & 64){//c terminal gap !!!!!!1
			//	fprintf(stderr,"c terminal gap\n");
				gap_cost += profa[27+64];
			//	fprintf(stderr,"i:%d	%d\n",c-1,gap_cost);
				for (i = 1; i < gap_len[c];i++){
					gap_cost += profa[29+64*i];
			//		fprintf(stderr,"i:%d	%d\n",i,gap_cost);
				}
			}else{
			//	fprintf(stderr,"middle gap\n");
				gap_cost += profa[27+64];
			//	fprintf(stderr,"i:%d	%d\n",c-1,gap_cost);
				for (i = 1; i < gap_len[c]-1;i++){
					gap_cost += profa[28+64*i];
			//		fprintf(stderr,"i:%d	%d\n",i,gap_cost);
				}
				gap_cost += profa[27+64*i];
			//	fprintf(stderr,"i:%d	%d\n",i,gap_cost);
			}
			
			gap_cost /=  gap_len[c];
			
			gap_cost *= internal_gap_weight;
			
			while(path[c] & 2 && path[c] != 3){
		//		fprintf(stderr,"gap_b	%d	%d cost:%d\n",c,path[c],gap_cost);
				for (i = 64; i--;){
					newp[i] = profa[i];
				}
				newp[23] += gap_cost;
				for (i = 32;i < 55;i++){
					newp[i] += gap_cost;
				}
				profa +=64;
				newp += 64;
				c++;
			}			
		}
	}
	for (i = 64; i--;){
		newp[i] =  profa[i] + profb[i];
	}	
	newp -= path[0] *64; 
	
	free(gap_len);
	
	return newp;
}


void smooth_gaps(float* prof,int len,int window,float strength)
{
	float tmp_gpo;
	float tmp_gpe;
	float tmp_tgpe;
	int i,j;
	if(!(window &1)){
		window--;
	}
	for ( i = (window/2); i < len - (window/2);i++){
		tmp_gpo = 0.0;
		tmp_gpe = 0.0;
		tmp_tgpe = 0.0;
		for (j = -(window/2); j < (window/2);j++){
			tmp_gpo += (float)prof[27+((i+j)*64)]*strength;
			tmp_gpe += (float) prof[28+((i+j)*64)]*strength;
			tmp_tgpe += (float) prof[29+((i+j)*64)]*strength;
		}
		tmp_gpo /= window;
		tmp_gpe /= window;
		tmp_tgpe /= window;
		prof[27+(i*64)] =  prof[27+(i*64)]*(1.0-strength) + tmp_gpo;
		prof[28+(i*64)] =  prof[28+(i*64)]*(1.0-strength) + tmp_gpe;
		prof[29+(i*64)] =  prof[29+(i*64)]*(1.0-strength) + tmp_tgpe;
	}
}


void increase_gaps(float* prof,int len,int window,float strength)
{
	float* mod = 0;
	int i,j,c;
	int start_pos = 0;
	int end_pos = 0; 
	
	mod = malloc(sizeof(float)*window);
	for ( i = 0; i < window;i++){
		mod[i] = (strength - i*(float)strength / (float) window) - (0.5*strength);
	}
	//only gpo first....
	for ( i = 0; i < len;i++){
// 	//	fprintf(stderr,"(%0.2f:%0.2f) ",prof[26],prof[23]);
		prof[26] = 0.0;
		prof+= 64;
	}
	prof -= len << 6;
	
	
	
	for ( i = 0; i < len;i++){
		
		if(prof[23]!= 0){
			
			start_pos = i-window;
			if(start_pos < 0){
				c = start_pos + window;	
			}else{
				c = window;
			}
			
			for ( j = c;j--;){
				prof[26 - (64*(j+1))] +=  mod[j];
			}
			end_pos = i+window;
			if(end_pos > len){
				c = len - i;	
			}else{
				c = window;
			}
			//fprintf(stderr,"%d %d\n",i,c);
			for (j = 0;j < c;j++){
				prof[26 +(64*(j+1))] +=  mod[j];
			}
		}
		prof+= 64;
	}
	prof -= len << 6;

	for ( i = 0; i < len;i++){
// 	//	fprintf(stderr,"(%0.2f:%0.2f) ",prof[26],prof[23]);
		prof[27] = prof[27] * (prof[26]+1.0);
		prof[28] = prof[28] * (prof[26]+1.0);
		prof[29] = prof[29] * (prof[26]+1.0);
		prof+= 64;
	}
	prof -= len << 6;
	
	free(mod);
}


void set_gap_penalties2(float* prof,int len,int nsip,int window,float strength)
{
	int i,j;
	float tmp_gpo;
	float tmp_gpe;
	float tmp_tgpe;
	
	prof +=  (64 *(len));
	
	prof[27] = prof[55]*nsip*-gpo;
	prof[28] = prof[55]*nsip*-gpe;
	prof[29] = prof[55]*nsip*-tgpe;

	i = len;
	while(i--){
		prof -= 64;
		prof[27] = prof[55]*nsip*-gpo;
		prof[28] = prof[55]*nsip*-gpe;
		
		prof[29] = prof[55]*nsip*-tgpe;
	}
	if(!(window &1)){
		window--;
	}
	
	
	for ( i = (window/2); i < len - (window/2);i++){
		tmp_gpo = 0.0;
		tmp_gpe = 0.0;
		tmp_tgpe = 0.0;
		for (j = -(window/2); j < (window/2);j++){
			tmp_gpo += (float)prof[27+((i+j)*64)]*strength;
			tmp_gpe += (float) prof[28+((i+j)*64)]*strength;
			tmp_tgpe += (float) prof[29+((i+j)*64)]*strength;
		}
		tmp_gpo /= window;
		tmp_gpe /= window;
		tmp_tgpe /= window;
		prof[27+(i*64)] =  prof[27+(i*64)]*(1-strength) + tmp_gpo;
		prof[28+(i*64)] =  prof[28+(i*64)]*(1-strength) + tmp_gpe;
		prof[29+(i*64)] =  prof[29+(i*64)]*(1-strength) + tmp_tgpe;
	}
	
	/*for ( i = 2; i < len-2;i++){
                prof[27+(i*64)]  = (prof[27+((i-2)*64)] +prof[27+((i-1)*64)] + prof[27+(i*64)] + prof[27+((i+1)*64)] +prof[27+((i+2)*64)])/ 5;
        }*/
       /* for ( i = 2; i < len-2;i++){
                prof[28+(i*64)]  = (prof[28+((i-2)*64)] + prof[28+((i-1)*64)] + prof[28+(i*64)] + prof[28+((i+1)*64)] +prof[28+((i+2)*64)])/ 5;
        }
        for ( i = 2; i < len-2;i++){
                prof[29+(i*64)]  = (prof[29+((i-2)*64)] + prof[29+((i-1)*64)] + prof[29+(i*64)] + prof[29+((i+1)*64)] +prof[29+((i+2)*64)])/ 5;
        }*/

}



float* make_profile2(float* prof, int* seq,int len,float** subm)
{
	int i,j,c;	
	prof = malloc(sizeof(float)*(len+1)*64);
	prof +=  (64 *len);
	
	for (i = 0;i < 64;i++){
		prof[i] = 0;
	}
	prof[55] = 1;

	i = len;
	while(i--){
		prof -= 64;

		for (j = 0;j < 64;j++){
			prof[j] = 0;
		}
		c = seq[i];
		
		prof[c] += 1;
					
		prof += 32;
		for(j = 23;j--;){
			prof[j] = subm[c][j];
			
		}
		prof[23] = 1;
		prof -= 32;

	}
	return prof;
}


float*  feature_update(const float* profa, const float* profb,float* newp,int* path,int stride)
{
	int i,c;
	c = 1;
	while(path[c] != 3){
		if (!path[c]){
			for (i = stride; i--;){
				newp[i] = profa[i] + profb[i];
			}
			profa += stride;
			profb += stride;
		}
		if (path[c] & 1){
			for (i = stride; i--;){
				newp[i] = profb[i];
			}
			profb += stride;	
			
		}
		if (path[c] & 2){
			for (i = stride; i--;){
				newp[i] = profa[i];
			}
			profa+=stride;			
		}
		newp += stride;
		c++;
	}
	for (i = stride; i--;){
		newp[i] =  profa[i] + profb[i];
	}	
	newp -= path[0] *stride;
	return newp;
}

float* make_wu_profile(float* prof,float* wu,int len)
{
	int i;	
	
	
	prof = malloc(sizeof(float)*(len+1)*2);
	
	
	for (i = 0;i < (len+1)*2;i++){
		prof[i] = 0;
	}
	for (i = 0; i < len;i++){
		if(!wu[i]){
			prof[i<<1] = 1;
			prof[(i<<1)+1] = 1;
			
		}else{
			prof[i<<1] = wu[i]+1;
			prof[(i<<1)+1] = wu[i]+1;
		}
	}
	return prof;
}


float* make_feature_profile(float* prof,struct feature* f,int len,struct feature_matrix* fm)
{
	int i,j;	
	
	
	prof = malloc(sizeof(int)*(len+1)*fm->stride);
	
	
	for (i = 0;i < (len+1)*fm->stride;i++){
		prof[i] = 0;
	}
	
	while(f){
		if(f->color != -1){
			if(f->start < len && f->end < len){
				for (i = f->start-1;i < f->end;i++){
					prof[i*fm->stride + f->color] += 1;
					for ( j =fm->mdim ;j < fm->stride;j++){
						prof[i*fm->stride+j] += fm->m[f->color][j-fm->mdim];
					}
				}
			}
		}
		f = f->next;
	}	
	return prof;
}



float* make_profile(float* prof, int* seq,int len, float** subm)
{
	int i,j,c;	
	prof = malloc(sizeof(float)*(len+2)*64);
	prof +=  (64 *(len+1));

	for (i = 0;i < 64;i++){
		prof[i] = 0;
	}
	prof[23+32] = -gpo;
	prof[24+32] = -gpe;
	prof[25+32] = -tgpe;

	
	i = len;
	while(i--){
		prof -= 64;

		for (j = 0;j < 64;j++){
			prof[j] = 0;
		}
		c = seq[i];
		
		prof[c] += 1;
		
		prof += 32;
		
		for(j = 23;j--;){
			prof[j] = subm[c][j];
		}
		prof[23] = -gpo;
		prof[24] = -gpe;
		prof[25] = -tgpe;
		
		prof -= 32;
	}
	prof -= 64;
	for (i = 0;i < 64;i++){
		prof[i] = 0;
	}
	prof[23+32] = -gpo;
	prof[24+32] = -gpe;
	prof[25+32] = -tgpe;	
	return prof;
}

float* dna_make_profile(float* prof, int* seq,int len,float** subm)
//int* make_profile(int* prof, int* seq,int len)
{
	int i,j,c;	
	prof = malloc(sizeof(float)*(len+2)*22);
	prof +=  (22 *(len+1));
	//fprintf(stderr,"Len:%d	%d\n",len,64*len);
	//for (i = 64;i--;){
	for (i = 0;i < 22;i++){
		prof[i] = 0;
	}
	prof[5+11] = -gpo;
	prof[6+11] = -gpe;
	prof[7+11] = -tgpe;

	
	i = len;
	while(i--){
		prof -= 22;
		//fprintf(stderr,"-64\n");
		//for (j = 64; j--;){
		for (j = 0;j < 22;j++){
			prof[j] = 0;
		}
		c = seq[i];
		
		prof[c] += 1;
		
		//n = feature[i];
		//prof[n+23] = 1;
		
		prof += 11;
		for(j = 5;j--;){
			prof[j] = subm[c][j];
		}
		prof[5] = -gpo;
		prof[6] = -gpe;
		prof[7] = -tgpe;
		prof -= 11;
	}
	prof -= 22;
	for (i = 0;i < 22;i++){
		prof[i] = 0;
	}
	prof[5+11] = -gpo;
	prof[6+11] = -gpe;
	prof[7+11] = -tgpe;
	
	return prof;
}




float* update(const float* profa, const float* profb,float* newp,int* path,int sipa,int sipb)
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
			//printf("open:%d	ext:%d	%d	%d\n",si->nsip[a] * gpo,si->nsip[a] * gpe,si->nsip[a] * profb[41],si->nsip[a] * profb[46]);
			for (i = 64; i--;){
				newp[i] = profb[i];
			}
			profb += 64;
			#ifndef SIMPLE 
			if(!(path[c] & 20)){
				if(path[c] & 32){
					newp[25] += sipa;//1;
					i = tgpe*sipa;
				}else{
					newp[24] += sipa;//1;
					i = gpe*sipa;
				}
				
				for (j = 32; j < 55;j++){
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
								
				for (j = 32; j < 55;j++){
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
				for (j = 32; j < 55;j++){
					newp[j] -=i;
				}
			}
			}
			#endif
			
			
		}
		if (path[c] & 2){
			//fprintf(stderr,"Gap_B:%d\n",c);
			//printf("open:%d	ext:%d	%d	%d\n",si->nsip[b] * gpo,si->nsip[b] * gpe,profa[26],profa[27]);
			for (i = 64; i--;){
				newp[i] = profa[i];
			}
			profa+=64;
			#ifndef SIMPLE 
			if(!(path[c] & 20)){
				if(path[c] & 32){
					newp[25] += sipb;//1;
					i = tgpe*sipb;
				}else{
					newp[24] += sipb;//1;
					i = gpe*sipb;
				}
				for (j = 32; j < 55;j++){
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
				for (j = 32; j < 55;j++){
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
				
				for (j = 32; j < 55;j++){
					newp[j] -=i;
				}
			}
			}
			#endif
			
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


float* update_only_a(const float* profa, const float* profb,float* newp,int* path,int sipa,int sipb)
{
	int i,c;
	for (i = 64; i--;){
		newp[i] = profa[i];// + profb[i];
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
				newp[i] = profa[i];// + profb[i];
			}
				
			
			profa += 64;
			profb += 64;
		}
		
		if (path[c] & 1){
			//fprintf(stderr,"Gap_A:%d\n",c);
			//printf("open:%d	ext:%d	%d	%d\n",si->nsip[a] * gpo,si->nsip[a] * gpe,si->nsip[a] * profb[41],si->nsip[a] * profb[46]);
			for (i = 64; i--;){
				newp[i] = 0.0;//profb[i];
			}
			profb += 64;
		}
		if (path[c] & 2){
			//fprintf(stderr,"Gap_B:%d\n",c);
			//printf("open:%d	ext:%d	%d	%d\n",si->nsip[b] * gpo,si->nsip[b] * gpe,profa[26],profa[27]);
			for (i = 64; i--;){
				newp[i] = profa[i];
			}
			profa+=64;
			
		}
		newp += 64;
		c++;
	}
	for (i = 64; i--;){
		newp[i] =  profa[i];// + profb[i];
	}	
	newp -= (path[0]+1) *64;
	return newp;
}



float* dna_update(const float* profa, const float* profb, float* newp,int* path,int sipa,int sipb)
{
	int i,j,c;
	
	for (i = 22; i--;){
		newp[i] = profa[i] + profb[i];
	}
	
	profa += 22;
	profb += 22;
	newp += 22;
	
	
	c = 1;
	while(path[c] != 3){
		//Idea: limit the 'virtual' number of residues of one type to x.
		// i.e. only allow a maximum of 10 alanines to be registered in each column
		// the penalty for aligning a 'G' to this column will stay stable even when many (>10) alanines are present.
		// the difference in score between the 'correct' (all alanine) and incorrect (alanines + glycine) will not increase 
		// with the number of sequences. -> see Durbin pp 140
		
		if (!path[c]){
			//fprintf(stderr,"Align	%d\n",c);
			for (i = 22; i--;){
				newp[i] = profa[i] + profb[i];
			}
				
			
			profa += 22;
			profb += 22;
		}
		if (path[c] & 1){
			//fprintf(stderr,"Gap_A:%d\n",c);
			//printf("open:%d	ext:%d	%d	%d\n",si->nsip[a] * gpo,si->nsip[a] * gpe,si->nsip[a] * profb[41],si->nsip[a] * profb[46]);
			for (i = 22; i--;){
				newp[i] = profb[i];
			}
			profb += 22;
			if(!(path[c]&20)){
				if(path[c]&32){
					newp[7] += sipa;//1;
					i = tgpe*sipa;
				}else{
					newp[6] += sipa;//1;
					i = gpe*sipa;
				}
				
				for (j = 11; j < 16;j++){
					newp[j] -=i;
				}
			}else{
			if (path[c] & 16){ 
	//			fprintf(stderr,"close_open");
				if(path[c]&32){
					newp[7] += sipa;//1;
					i = tgpe*sipa;
					newp[5] += sipa;//1;
					i += gpo*sipa;
				}else{
					newp[5] += sipa;//1;
					i = gpo*sipa;
				}
								
				for (j = 11; j < 16;j++){
					newp[j] -=i;
				}
			}
			if (path[c] & 4){ 
	//			fprintf(stderr,"Gap_open");
				if(path[c]&32){
					newp[7] += sipa;//1;
					i = tgpe*sipa;
					newp[5] += sipa;//1;
					i += gpo*sipa;
				}else{
					newp[5] += sipa;//1;
					i = gpo*sipa;
				}
				for (j = 11; j < 16; j++){
					newp[j] -=i;
				}
			}
			}
			
			
		}
		if (path[c] & 2){
			//fprintf(stderr,"Gap_B:%d\n",c);
			//printf("open:%d	ext:%d	%d	%d\n",si->nsip[b] * gpo,si->nsip[b] * gpe,profa[26],profa[27]);
			for (i = 22; i--;){
				newp[i] = profa[i];
			}
			profa+=22;
			if(!(path[c]&20)){
				if(path[c]&32){
					newp[7] += sipb;//1;
					i = tgpe*sipb;
				}else{
					newp[6] += sipb;//1;
					i = gpe*sipb;
				}
				for (j = 11; j < 16;j++){
					newp[j] -=i;
				}
			}else{
			if (path[c] & 16){
	//			fprintf(stderr,"close_open");
				if(path[c]&32){
					newp[7] += sipb;//1;
					i =  tgpe*sipb;
					newp[5] += sipb;//1;
					i +=  gpo*sipb;
				}else{
					newp[5] += sipb;//1;
					i =  gpo*sipb;
				}
				for (j = 11; j < 16;j++){
					newp[j] -=i;
				}
			}
			if (path[c] & 4){
	//			fprintf(stderr,"Gap_open");
				if(path[c]&32){
					newp[7] += sipb;//1;
					i = tgpe*sipb;
					newp[5] += sipb;//1;
					i += gpo*sipb;
				}else{
					newp[5] += sipb;//1;
					i = gpo*sipb;
				}
				
				for (j = 11; j < 16;j++){
					newp[j] -=i;
				}
			}
			}
			
		}
		newp += 22;
		c++;
	}
	for (i = 22; i--;){
		newp[i] =  profa[i] + profb[i];
	}	
	newp -= (path[0]+1) *22; 
	return newp;
}


float* dna_update_only_a(const float* profa, const float* profb, float* newp,int* path,int sipa,int sipb)
{
	int i,c;
	
	for (i = 22; i--;){
		newp[i] = profa[i];// + profb[i];
	}
	
	profa += 22;
	profb += 22;
	newp += 22;
	
	
	c = 1;
	while(path[c] != 3){
		//Idea: limit the 'virtual' number of residues of one type to x.
		// i.e. only allow a maximum of 10 alanines to be registered in each column
		// the penalty for aligning a 'G' to this column will stay stable even when many (>10) alanines are present.
		// the difference in score between the 'correct' (all alanine) and incorrect (alanines + glycine) will not increase 
		// with the number of sequences. -> see Durbin pp 140
		
		if (!path[c]){
			//fprintf(stderr,"Align	%d\n",c);
			for (i = 22; i--;){
				newp[i] = profa[i];//+ profb[i];
			}
				
			
			profa += 22;
			profb += 22;
		}
		if (path[c] & 1){
			//fprintf(stderr,"Gap_A:%d\n",c);
			//printf("open:%d	ext:%d	%d	%d\n",si->nsip[a] * gpo,si->nsip[a] * gpe,si->nsip[a] * profb[41],si->nsip[a] * profb[46]);
			for (i = 22; i--;){
				newp[i] = 0.0f;//profb[i];
			}
			profb += 22;
			newp[5] = 1000000;
			newp[6] = 1000000;
			newp[7] = 1000000;
		}
		if (path[c] & 2){
			//fprintf(stderr,"Gap_B:%d\n",c);
			//printf("open:%d	ext:%d	%d	%d\n",si->nsip[b] * gpo,si->nsip[b] * gpe,profa[26],profa[27]);
			for (i = 22; i--;){
				newp[i] = profa[i];
			}
			profa+=22;
		}
		newp += 22;
		c++;
	}
	for (i = 22; i--;){
		newp[i] =  profa[i];// + profb[i];
	}	
	newp -= (path[0]+1) *22; 
	return newp;
}


void dna_set_gap_penalties(float* prof,int len,int nsip,float strength,int nsip_c)
{
	int i;
	int j;
	float res = (float)nsip_c;
	float local_res = 0;
	float w = 0.0;
	prof +=  (22 *(len+1));
	
	local_res = 0;
	for(j = 0; j < 5;j++){
		local_res+=  prof[j];
	}
		
	w = 1.0 + (((local_res - 1.0 )/ res) * strength);
	
	prof[8] = prof[16]*nsip*w;//gap open or close
	prof[9] = prof[17]*nsip*w;//gap extention
	prof[10] = prof[18]*nsip*w;//gap open or close

	//prof[30] = prof[58]*nsip;//gap extention
	
	
	i = len+1;
	while(i--){
		prof -= 22;
		local_res = 0;
		for(j = 0; j < 5;j++){
			local_res+=  prof[j];
		}
	
	
		w = 1.0 + (((local_res - 1.0 )/ res) * strength);
		
		prof[8] = prof[16]*nsip*w;//gap open or close
		prof[9] = prof[17]*nsip*w;//gap extention
		
		prof[10] = prof[18]*nsip*w;//gap open or close
	//	prof[30] = prof[58]*nsip;//gap extention
	}
}

void set_gap_penalties(float* prof,int len,int nsip,float strength,int nsip_c)
{
	int i;
	int j;
	float res = (float)nsip_c;
	float local_res = 0;
	float w = 0.0;
	
	prof +=  (64 *(len+1));
	
	
	local_res = 0;
	for(j = 0; j < 23;j++){
		local_res+=  prof[j];
	}
	
	
	w = 1.0 + (((local_res - 1.0 )/ res) * strength);
	
	prof[27] = prof[55]*nsip*w;//gap open or close  23
	prof[28] = prof[56]*nsip*w;//gap extention 24 
		
	prof[29] = prof[57]*nsip*w;//gap open or close 25 
	i = len+1;
	while(i--){
		prof -= 64;

		local_res = 0;
		for(j = 0; j < 23;j++){
			local_res+=  prof[j];
		}
	
	
		w = 1.0 + (((local_res - 1.0 )/ res) * strength);
			
		prof[27] = prof[55]*nsip*w;//gap open or close
		prof[28] = prof[56]*nsip*w;//gap extention
		
		prof[29] = prof[57]*nsip*w;//gap open or close
	}
}

