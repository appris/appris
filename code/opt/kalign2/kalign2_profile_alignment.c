/*
	kalign2_profile_alignment.c
	
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
#include "kalign2_profile_alignment.h"
#include "kalign2_hirschberg.h"

void profile_alignment_main(struct alignment* aln,struct parameters* param,float** submatrix)
{
	float** dm = 0;
	int* tree = 0;
	struct aln_tree_node* tree2 = 0;
	int i,j;
	int tmp_numseq;
	int tmp_numprofiles;
	
	local_numseq = 0;
	local_numprofiles = 0;
	
	//determine number of profiles that were inputted....
	
	while(aln->sl[local_numseq+numseq]){
		local_numseq++;
	}
	
	local_numprofiles = (local_numseq << 1) - 1;
	//fprintf(stderr,"%d	%d\n",local_numseq,local_numprofiles);
	
	for (i = 0;i < numseq;i++){
	//	fprintf(stderr,"%d	%d	%d\n",i,aln->s[i][0],aln->s[i][1]);
		aln->s[i] = assign_gap_codes(aln->s[i],aln->sl[i]); 
	}
	
	if(param->dna == 1){
		if(byg_start(param->tree,"njNJ") != -1){
			dm =  dna_profile_distance(aln,dm,param,1);
		}else{
			dm =  dna_profile_distance(aln,dm,param,0);
		}
	}else{
		if(byg_start(param->tree,"njNJ") != -1){
			dm =  protein_profile_wu_distance(aln,dm,param,1);
		}else{
			dm =  protein_profile_wu_distance(aln,dm,param,0);
		}
	}
	/*for ( i=0; i < local_numseq;i++){
		for (j = 0;j < local_numseq;j++){
			fprintf(stderr,"%f ",dm[i][j]);
		}
		fprintf(stderr,"\n");
	}*/
	
	tmp_numseq = numseq;
	tmp_numprofiles = numprofiles;
	
	numseq = local_numseq;
 	numprofiles = local_numprofiles;
	
	if(byg_start(param->tree,"njNJ") != -1){
		tree2 = real_nj(dm,param->ntree);
	}else{
		tree2 = real_upgma(dm,param->ntree);
	}
	
	
	
	//WAs here need too add tree2 -> treee..... 
	
	
	tree = malloc(sizeof(int)*(numseq*3+1));
	for ( i = 1; i < (numseq*3)+1;i++){
		tree[i] = 0;
	}
	tree[0] = 1;
	tree = readtree(tree2,tree);
	for (i = 0; i < (numseq*3);i++){
		tree[i] = tree[i+1]+ tmp_numseq;
	}
	//exit(0);
	
	numseq = tmp_numseq;
	numprofiles = tmp_numprofiles;
	
	int** map = 0;
	
	map =  hirschberg_profile_alignment(aln,tree,submatrix, map);
	//clear up sequence array to be reused as gap array....
	int *p = 0;
	for (i = 0; i < numseq;i++){
		p = aln->s[i];
		for (j = 0; j < aln->sl[i];j++){
			p[j] = 0;
		}
	}
	//clear up
	int a,b,c;
	for (i = 0; i < (local_numseq-1)*3;i +=3){
		a = tree[i];
		b = tree[i+1];
		c =  tree[i+2];
		aln = make_seq(aln,a,b,map[c]);
	}

	for (i = 0; i < numseq;i++){
		aln->nsip[i] = 0;
	}
	aln =  sort_sequences(aln,tree,param->sort);

	//for (i = 0; i < numseq;i++){
	//	fprintf(stderr,"%d	%d	%d\n",i,aln->nsip[i],aln->sip[i][0]);
	//}
	
	
	output(aln,param);
	
	
	free(tree2->links);
	free(tree2->internal_lables);
	free(tree2);
	

	free(map);
	free(tree);
	exit(0);
}


int** hirschberg_profile_alignment(struct alignment* aln,int* tree,float**submatrix, int** map)
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

	for (i = 0; i < (local_numseq-1);i++){
		a = tree[i*3];
		b = tree[i*3+1];
		c = tree[i*3+2];
		fprintf(stderr,"\r%8.0f percent done",(float)(i) /(float)local_numseq * 100);
		//fprintf(stderr,"Aligning:%d %d->%d	done:%f\n",a,b,c,((float)(i)/(float)local_numseq)*100);
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
		
		if (a < numseq+local_numseq){
			profile[a] = make_profile_from_alignment(profile[a],a,aln,submatrix);
		}
		if (b < numseq+local_numseq){
			profile[b] = make_profile_from_alignment(profile[b],b,aln,submatrix);
		}
		
		
		set_gap_penalties(profile[b],len_b,aln->nsip[a],0,aln->nsip[b]);
		set_gap_penalties(profile[a],len_a,aln->nsip[b],0,aln->nsip[a]);
		
		hm->starta = 0;
		hm->startb = 0;
		hm->enda = len_a;
		hm->endb = len_b;
		hm->len_a = len_a;
		hm->len_b = len_b;
		
		hm->f[0].a = 0;
		hm->f[0].ga =  -FLOATINFTY;
		hm->f[0].gb = -FLOATINFTY;
		hm->b[0].a = 0;
		hm->b[0].ga =  -FLOATINFTY;
		hm->b[0].gb =  -FLOATINFTY;
	//	fprintf(stderr,"LENA:%d	LENB:%d	numseq:%d\n",len_a,len_b,numseq);
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
		
		map[c] = add_gap_info_to_hirsch_path(map[c],len_a,len_b);

		if(i != local_numseq-2){
			//fprintf(stderr,"updating....\n");
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


int* assign_gap_codes(int* seq,int len)
{
	int i;
	if(seq[0] < 0 && seq[1] < 0){
		seq[0] = -2;
	}
	
	for(i = 1; i < len;i++){
		if(seq[i-1] < 0 && seq[i] < 0){
			seq[i] = -2;
		}
		if(seq[i-1] < 0 && seq[i] >= 0){
			seq[i-1] = -1;
		}
	}
	i = 0;
	while(seq[i] < 0){
		if(seq[i] == -2){
			seq[i] = -3;
		}
		i++;
	}
	i = len-1;
	while(seq[i] < 0){
		if(seq[i] == -2){
			seq[i] = -3;
		}
		i--;
	}

	
	return seq;
}



float* make_profile_from_alignment(float* prof, int num,struct alignment* aln,float** subm)
{
	int i,j,c;	
	int a;
	int len = aln->sl[num];
	int* seq = 0;
	prof = malloc(sizeof(float)*(len+2)*64);
	for ( i = 0; i < (len+2)*64;i++){
		prof[i] = 0;
	}
	
	for ( a = 0; a < aln->nsip[num];a++){
	//	fprintf(stderr,"SEQ:%d\n",a);
		seq = aln->s[aln->sip[num][a]];
		prof +=  (64 *(len+1));

		prof[23+32] -= gpo;
		prof[24+32] -= gpe;
		prof[25+32] -= tgpe;
	
		
		i = len;
		while(i--){
			prof -= 64;
	
			c = seq[i];
			if(c >= 0){
	//			if(i == 0){
	//				fprintf(stderr,"%d	\n",c);
	//			}
				prof[c] += 1;
				prof += 32;
				for(j = 23;j--;){
					prof[j] += subm[c][j];
				}
				prof[23] -= gpo;
				prof[24] -= gpe;
				prof[25] -= tgpe;
				prof -= 32;
			}else if(c == -1){
				prof[23] += 1;
				for (j = 32; j < 55;j++){
					prof[j] -= gpo;
				}
			}else if(c == -2){
				prof[24] += 1;
				for (j = 32; j < 55;j++){
					prof[j] -= gpe;
				}
			}else if(c == -3){
				prof[25] += 1;
				for (j = 32; j < 55;j++){
					prof[j] -= tgpe;
				}
			}
		}
		prof -= 64;
		prof[23+32] -= gpo;
		prof[24+32] -= gpe;
		prof[25+32] -= tgpe;
	}
	return prof;
}


float** protein_profile_wu_distance(struct alignment* aln,float** dm,struct parameters* param, int nj)
{
	struct bignode* hash[1024];
	int*p =0;
	int i,j,m,n,a,b;
	unsigned int hv;
		
	int** local_seq = 0;
	int* local_sl = 0;
	
	local_seq = malloc(sizeof(int*)*numseq);
	local_sl = malloc(sizeof(int)*numseq);
	for(i = 0; i< numseq;i++){
		local_seq[i] = malloc(sizeof(int)*aln->sl[i]);
		a = 0;
		p = aln->s[i];
		for (j = 0;j < aln->sl[i];j++){
			if(p[j] >= 0){
				local_seq[i][a] = p[j];
				a++;
			}
		}
		local_sl[i] = a;
	}
	//determine number of profiles that were inputted....

	
	for (i = 0;i < 1024;i++){
		hash[i] = 0;
	}	

	if (nj){
		dm = malloc (sizeof(float*)*local_numprofiles);
		for (i = local_numprofiles;i--;){
			dm[i] = malloc (sizeof (float)*(local_numprofiles));
			for (j = local_numprofiles;j--;){
				dm[i][j] = 0.0f;
			}
		}
	}else{
		dm = malloc (sizeof(float*)*local_numseq);
		for (i = local_numseq;i--;){
			dm[i] = malloc (sizeof (float)*(local_numseq));
			for (j = local_numseq;j--;){
				dm[i][j] = 0.0f;
			}
		}
	}
	fprintf(stderr,"Distance Calculation:\n");
	b = (numseq*(numseq-1))/2;
	a = 1;	
	
	//fprintf(stderr,"%d	%d	%8.0f\n",a,b,(float)a /(float)b * 100);
		
	for (i = 0; i < numseq-1;i++){
		m = is_member(aln,i);
		p = local_seq[i];
		//p = malloc(sizeof(int)
				
		for (j = local_sl[i]-2;j--;){
	//		hv = (p[j+1] << 5) + p[j+2];
	//		hash[hv] = insert_hash(hash[hv],j);
			hv = (p[j] << 5) + p[j+1];
			hash[hv] = big_insert_hash(hash[hv],j);
			
			hv = (p[j] << 5) + p[j+2];
			hash[hv] = big_insert_hash(hash[hv],j);
		//	hv = (si->s[i][j+1] << 5) +t2;
		//	hash[hv] = insert_hash(hash[hv],j);
		}
		for (j = i+1; j < numseq;j++){
			n = is_member(aln,j);
			if(n != m){
				//fprintf(stderr,"%d	%d\n",n,m);
				p = local_seq[j];
				dm[m][n] += protein_wu_distance_calculation(hash,p,local_sl[j],local_sl[j]+local_sl[i],param->zlevel);
				//fprintf(stderr,"%d-%d(%d	%d):%f\n",m,n,i,j,dm[m][n]);
				//exit(0);
				//dm[i][j] /= min;
				dm[n][m] = dm[m][n];
				
			}
			fprintf(stderr,"\r%8.0f percent done",(float)a /(float)b * 100);
			a++;
			
		}
		
		
		for (j = 1024;j--;){
			if (hash[j]){
				big_remove_nodes(hash[j]);
				hash[j] = 0;
			}
		}	
	}
	
	for(i = 0; i< numseq;i++){
		free(local_seq[i]);
	}
	free(local_seq);
	free(local_sl);
	return dm;
}

float** dna_profile_distance(struct alignment* aln,float** dm,struct parameters* param, int nj)
{
	struct bignode* hash[1024];
	
	int *p = 0;
	int i,j,a,b,m,n;
	unsigned int hv;
	int** local_seq = 0;
	int* local_sl = 0;
	
	local_seq = malloc(sizeof(int*)*numseq);
	local_sl = malloc(sizeof(int)*numseq);
	for(i = 0; i< numseq;i++){
		local_seq[i] = malloc(sizeof(int)*aln->sl[i]);
		a = 0;
		p = aln->s[i];
		for (j = 0;j < aln->sl[i];j++){
			if(p[j] >= 0){
				local_seq[i][a] = p[j];
				a++;
			}
		}
		local_sl[i] = a;
	}
	
	fprintf(stderr,"Distance Calculation:\n");
	
	
	for (i = 0;i < 1024;i++){
		hash[i] = 0;
	}	

	if (nj){
		dm = malloc (sizeof(float*)*local_numprofiles);
		for (i = local_numprofiles;i--;){
			dm[i] = malloc (sizeof (float)*(local_numprofiles));
			for (j = local_numprofiles;j--;){
				dm[i][j] = 0.0f;
			}
		}
	}else{
		dm = malloc (sizeof(float*)*local_numseq);
		for (i = local_numseq;i--;){
			dm[i] = malloc (sizeof (float)*(local_numseq));
			for (j = local_numseq;j--;){
				dm[i][j] = 0.0f;
			}
		}
	}

	b = (numseq*(numseq-1))/2;
	a = 1;	
	
	for (i = 0; i < numseq-1;i++){
		m = is_member(aln,i);
		p = local_seq[i];
		for (j = local_sl[i]-5;j--;){
			hv = ((p[j]&3)<<8) + ((p[j+1]&3)<<6) + ((p[j+2]&3)<<4)  + ((p[j+3]&3)<<2) + (p[j+4]&3);//ABCDE
			hash[hv] = big_insert_hash(hash[hv],j);
			hv = ((p[j]&3)<<8) + ((p[j+1]&3)<<6) + ((p[j+2]&3)<<4)  + ((p[j+3]&3)<<2) + (p[j+5]&3);//ABCDF
			hash[hv] = big_insert_hash(hash[hv],j);
			hv = ((p[j]&3)<<8) + ((p[j+1]&3)<<6) + ((p[j+2]&3)<<4)  + ((p[j+4]&3)<<2) + (p[j+5]&3);//ABCEF
			hash[hv] = big_insert_hash(hash[hv],j);
			hv = ((p[j]&3)<<8) + ((p[j+1]&3)<<6) + ((p[j+3]&3)<<4)  + ((p[j+4]&3)<<2) + (p[j+5]&3);//ABDEF
			hash[hv] = big_insert_hash(hash[hv],j);
			hv = ((p[j]&3)<<8) + ((p[j+2]&3)<<6) + ((p[j+3]&3)<<4) + ((p[j+4]&3)<<2) + (p[j+5]&3);//ACDEF
			hash[hv] = big_insert_hash(hash[hv],j);
		}
		for (j = i+1; j < numseq;j++){
			n = is_member(aln,j);
			if(n != m){
				p = local_seq[j];
				//min =  (si->sl[i] > si->sl[j]) ?si->sl[j] :si->sl[i];
				dm[m][n] += dna_distance_calculation(hash,p,local_sl[j],local_sl[j]+local_sl[i],param->zlevel);
				//dm[i][j] /= (aln->sl[i] > aln->sl[j]) ? aln->sl[j] : aln->sl[i];
				dm[n][m] = dm[m][n];
			}
			fprintf(stderr,"\r%8.0f percent done",(float)a /(float)b * 100);
			a++;
		}
		
		for (j = 1024;j--;){
			if (hash[j]){
				big_remove_nodes(hash[j]);
				hash[j] = 0;
			}
		}
	}
	
	for(i = 0; i< numseq;i++){
		free(local_seq[i]);
	}
	free(local_seq);
	free(local_sl);
	return dm;
}

int is_member(struct alignment* aln,int test)
{
	int i,j;
	for (i = numseq;i<numseq+local_numseq;i++){
		for(j = 0;j < aln->nsip[i];j++){
			if(aln->sip[i][j] == test){
				return i-numseq;
			}
		}
	}
	return -1;
}

