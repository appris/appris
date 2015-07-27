/*
	kalign2_alignment_types.c 
	
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

int** default_alignment(struct alignment* aln,int* tree,float**submatrix, int** map)
{
	struct dp_matrix *dp = 0;
	int i,j,g,a,b,c;
	int len_a;
	int len_b;
	float** profile = 0;
	float* profa = 0;
	float* profb = 0;
	
	profile = malloc(sizeof(float*)*numprofiles);
	for ( i = 0;i< numprofiles;i++){
		profile[i] = 0;
	}

	map = malloc(sizeof(int*)*numprofiles);
	for ( i = 0;i < numprofiles;i++){
		map[i] = 0;
	}
	
	
	dp = dp_matrix_alloc(dp,511,511);
	
	fprintf(stderr,"\nAlignment:\n");
	
	//c = numseq;
	for (i = 0; i < (numseq-1);i++){
		a = tree[i*3];
		b = tree[i*3+1];
		c = tree[i*3+2];
		fprintf(stderr,"\r%8.0f percent done",(float)(i) /(float)numseq * 100);
		//fprintf(stderr,"Aligning:%d %d->%d	%d	%d\n",a,b,c,numseq,i);
		len_a = aln->sl[a];
		len_b = aln->sl[b];
		dp = dp_matrix_realloc(dp,len_a,len_b);
	
		map[c] = malloc(sizeof(int) * (len_a+len_b+2));
		for (j = len_a+len_b+2;j--;){
			map[c][j] = 0;
		}
		if (a < numseq){
			profile[a] = make_profile(profile[a],aln->s[a],len_a,submatrix);
		}
		if (b < numseq){
			profile[b] = make_profile(profile[b],aln->s[b],len_b,submatrix);
		}
		profa = profile[a]+64;
		profb = profile[b]+64;
	
		set_gap_penalties(profile[a],len_a,aln->nsip[b],0,aln->nsip[a]);
		set_gap_penalties(profile[b],len_b,aln->nsip[a],0,aln->nsip[b]);
		if(aln->nsip[a] == 1){
			if(aln->nsip[b] == 1){
				map[c] = ss_dyn(submatrix,map[c],dp,aln->s[a],aln->s[b],len_a,len_b);
			}else{
				map[c] = ps_dyn(map[c],dp,profb,aln->s[a],len_b,len_a,aln->nsip[b]);
				map[c] = mirror_path(map[c]);
			}
		}else{
			if(aln->nsip[b] == 1){
				map[c] = ps_dyn(map[c],dp,profa,aln->s[b],len_a,len_b,aln->nsip[a]);
			}else{
				if (len_a > len_b){			
					map[c] = pp_dyn(map[c],dp,profa,profb,len_a,len_b);
				}else{
					map[c] = pp_dyn(map[c],dp,profb,profa,len_b,len_a);
					map[c] = mirror_path(map[c]);
				}
			}
		}
			
		profile[c] = malloc(sizeof(float)*64*(len_a+len_b+2));

		profile[c] = update(profile[a],profile[b],profile[c],map[c],aln->nsip[a],aln->nsip[b]);

	
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
	free(profile[numprofiles-1]);
	free(profile);
	
	dp_matrix_free(dp);
	for (i = 32;i--;){
		free(submatrix[i]);
	}
	free(submatrix);
	return map;
}

/*
int** aa_alignment(struct alignment* aln,int* tree,int**submatrix, int** map,int mmbonus)
{
	struct dp_matrix *dp = 0;
	int i,j,g,a,b,c;
	int len_a;
	int len_b;
	int** profile = 0;
	int* profa = 0;
	int* profb = 0;
	
	int pbonus = 0;
	
	profile = malloc(sizeof(int*)*numprofiles);
	for ( i = 0;i< numprofiles;i++){
		profile[i] = 0;
	}

	map = malloc(sizeof(int*)*numprofiles);
	for ( i = 0;i < numprofiles;i++){
		map[i] = 0;
	}
	
	
	dp = dp_matrix_alloc(dp,511,511);
	c = numseq;
	for (i = 0; i < (numseq-1);i++){
		a = tree[i*3];
		b = tree[i*3+1];
		c = tree[i*3+2];
		fprintf(stderr,"Aligning:%d %d->%d\n",a,b,c);
		len_a = aln->sl[a];
		len_b = aln->sl[b];
		dp = dp_matrix_realloc(dp,len_a,len_b);
	
		map[c] = malloc(sizeof(int) * (len_a+len_b+2));
		for (j = len_a+len_b+2;j--;){
			map[c][j] = 0;
		}
		if (a < numseq){
			profile[a] = make_profile(profile[a],aln->s[a],len_a,submatrix);
		}
		if (b < numseq){
			profile[b] = make_profile(profile[b],aln->s[b],len_b,submatrix);
		}
		profa = profile[a];
		profb = profile[b];
	
		set_gap_penalties(profa,len_a,aln->nsip[b]);
		set_gap_penalties(profb,len_b,aln->nsip[a]);
		
		pbonus = mmbonus * aln->nsip[a] * aln->nsip[b];
		
		if (len_a > len_b){
			map[c] = aapp_dyn(map[c],dp,profa,profb,len_a,len_b,pbonus);
		}else{
			map[c] = aapp_dyn(map[c],dp,profb,profa,len_b,len_a,pbonus);
			map[c] = mirror_path(map[c]);
		}
			
		profile[c] = malloc(sizeof(int)*64*(len_a+len_b+2));
		profile[c] = update(profa,profb,profile[c],map[c],aln->nsip[a],aln->nsip[b]);
	
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
		free(profa);
		free(profb);
	}
	
	free(profile[numprofiles-1]);
	free(profile);
	
	dp_matrix_free(dp);
	for (i = 32;i--;){
		free(submatrix[i]);
	}
	free(submatrix);
	return map;
}*/

/*
int** alter_gaps_alignment(struct alignment* aln,int* tree,int**submatrix, int** map,int n,float range,int weight)
{
	struct dp_matrix *dp = 0;
	int i,j,g,a,b,c;
	
	int org_gpo = gpo;
	int org_gpe = gpe;
	int org_tgpe = tgpe;
	
	float gpo_step = 0;
	float gpe_step = 0;
	float tgpe_step = 0;
	
	int len_a;
	int len_b;
	int** profile = 0;
	int* profa = 0;
	int* profb = 0;
	int* path = 0;
	
	int* fprofa = 0;
	int* fprofb = 0;
	
	if(!(n &1)){
		n--;
	}
	
	float per = 0.0;
	
	per =(float) range*2/(n+1);

	gpo_step = (float)gpo * per;
	gpe_step = (float)gpe * per;
	tgpe_step = (float)tgpe * per;
	
	
	profile = malloc(sizeof(int*)*numprofiles);
	for ( i = 0;i< numprofiles;i++){
		profile[i] = 0;
	}
	
	map = malloc(sizeof(int*)*numprofiles);
	for ( i = 0;i < numprofiles;i++){
		map[i] = 0;
	}
	
	dp = dp_matrix_alloc(dp,511,511);
	c = numseq;
	
	for (i = 0; i < (numseq-1);i++){
		a = tree[i*3];
		b = tree[i*3+1];
		c = tree[i*3+2];
		fprintf(stderr,"Aligning:%d %d->%d\n",a,b,c);
		len_a = aln->sl[a];
		len_b = aln->sl[b];
		dp = dp_matrix_realloc(dp,len_a,len_b);
	
		map[c] = malloc(sizeof(int) * (len_a+len_b+2));
		for (j = len_a+len_b+2;j--;){
			map[c][j] = 0;
		}
		if (a < numseq){
			profile[a] = make_profile(profile[a],aln->s[a],len_a,submatrix);
		}
		if (b < numseq){
			profile[b] = make_profile(profile[b],aln->s[b],len_b,submatrix);
		}
		profa = profile[a];
		profb = profile[b];
		
		fprofa = malloc(sizeof(int)*(len_a+1)*2);
		for (j = 0;j < (len_a+1)*2;j++){
			fprofa[j] = 0;
		}
		fprofb = malloc(sizeof(int)*(len_b+1)*2);
		for (j = 0;j < (len_b+1)*2;j++){
			fprofb[j] = 0;
		}
		
		gpo = org_gpo - ((int)gpo_step* (n/2));
		gpe = org_gpe - ((int)gpe_step* (n/2));
		tgpe = org_tgpe - ((int)tgpe_step* (n/2));
		
		for (j = 0; j < n;j++){
			set_gap_penalties(profa,len_a,aln->nsip[b]);
			set_gap_penalties(profb,len_b,aln->nsip[a]);
			
			path = malloc(sizeof(int) * (len_a+len_b+2));
			for (g = len_a+len_b+2;g--;){
				path[g] = 0;
			}
			
			if(aln->nsip[a] == 1){
				if(aln->nsip[b] == 1){
					path = ss_dyn(submatrix,path,dp,aln->s[a],aln->s[b],len_a,len_b);
				}else{
					path = ps_dyn(path,dp,profb,aln->s[a],len_b,len_a,aln->nsip[b]);
					path = mirror_path(path);
				}
			}else{
				if(aln->nsip[b] == 1){
					path = ps_dyn(path,dp,profa,aln->s[b],len_a,len_b,aln->nsip[a]);
				}else{
					if (len_a > len_b){
						path = pp_dyn(path,dp,profa,profb,len_a,len_b);
					}else{
						path = pp_dyn(path,dp,profb,profa,len_b,len_a);
						path = mirror_path(path);
					}
				}
			}
			fprintf(stderr,"Test alignment with gpo:%d	gpe:%d	tgpe:%d\n",gpo,gpe,tgpe);
			
			
			add_feature_information_from_alignment(path,fprofa,fprofb,weight/n);
			

			gpo += 	(int)gpo_step;
			gpe +=  (int)gpe_step;
			tgpe += (int)tgpe_step;
		}
		gpo = org_gpo;
		gpe = org_gpe;
		tgpe = org_tgpe;
	
		set_gap_penalties(profa,len_a,aln->nsip[b]);
		set_gap_penalties(profb,len_b,aln->nsip[a]);
		

		
		if (len_a > len_b){
		//	map[c] = f_only_pp_dyn(map[c],dp,fprofa,fprofb,len_a,len_b,1,2);
			map[c] = fpp_dyn(map[c],dp,profa,profb,fprofa,fprofb,len_a,len_b,1,2);
		}else{
		//	map[c] = f_only_pp_dyn(map[c],dp,fprofb,fprofa,len_b,len_a,1,2);
			map[c] = fpp_dyn(map[c],dp,profb,profa,fprofb,fprofa,len_b,len_a,1,2);
			map[c] = mirror_path(map[c]);
		}
					
		profile[c] = malloc(sizeof(int)*64*(len_a+len_b+2));
		profile[c] = update(profa,profb,profile[c],map[c],aln->nsip[a],aln->nsip[b]);
		
		
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
		free(profa);
		free(profb);
		
		free(fprofa);
		free(fprofb);
	}
	
	free(profile[numprofiles-1]);
	free(profile);
		
	dp_matrix_free(dp);
	for (i = 32;i--;){
		free(submatrix[i]);
	}
	free(submatrix);
	return map;
}*/

/*
int** test_alignment(struct alignment* aln,int* tree,float **submatrix, int** map,float internal_gap_weight,int window,float strength)
{
	struct dp_matrix *dp = 0;
	int i,j,g,a,b,c;
	int len_a;
	int len_b;
	float** profile = 0;
	float* profa = 0;
	float* profb = 0;
	
	profile = malloc(sizeof(float*)*numprofiles);
	for ( i = 0;i< numprofiles;i++){
		profile[i] = 0;
	}

	map = malloc(sizeof(int*)*numprofiles);
	for ( i = 0;i < numprofiles;i++){
		map[i] = 0;
	}
	
	
	dp = dp_matrix_alloc(dp,511,511);
	c = numseq;
	for (i = 0; i < (numseq-1);i++){
		a = tree[i*3];
		b = tree[i*3+1];
		c = tree[i*3+2];
		fprintf(stderr,"Aligning:%d %d->%d\n",a,b,c);
		len_a = aln->sl[a];
		len_b = aln->sl[b];
		dp = dp_matrix_realloc(dp,len_a,len_b);
	
		map[c] = malloc(sizeof(int) * (len_a+len_b+2));
		for (j = len_a+len_b+2;j--;){
			map[c][j] = 0;
		}
		if (a < numseq){
			profile[a] = make_profile2(profile[a],aln->s[a],len_a,submatrix);
		}
		if (b < numseq){
			profile[b] = make_profile2(profile[b],aln->s[b],len_b,submatrix);
		}
		profa = profile[a];
		profb = profile[b];
	
		set_gap_penalties2(profa,len_a,aln->nsip[b],window,strength);
		set_gap_penalties2(profb,len_b,aln->nsip[a],window,strength);

		if(aln->nsip[a] == 1){
			if(aln->nsip[b] == 1){
				map[c] = ss_dyn2(submatrix,map[c],dp,aln->s[a],aln->s[b],len_a,len_b);
			}else{
			//	map[c] = ps_dyn2(map[c],dp,profb,aln->s[a],len_b,len_a,aln->nsip[b]);
				
				map[c] = pp_dyn2(map[c],dp,profb,profa,len_b,len_a);
				map[c] = mirror_path(map[c]);
			}
		}else{
			if(aln->nsip[b] == 1){
			//	map[c] = ps_dyn2(map[c],dp,profa,aln->s[b],len_a,len_b,aln->nsip[a]);
				map[c] = pp_dyn2(map[c],dp,profa,profb,len_a,len_b);
			}else{
				if (len_a > len_b){
					map[c] = pp_dyn2(map[c],dp,profa,profb,len_a,len_b);
				}else{
					map[c] = pp_dyn2(map[c],dp,profb,profa,len_b,len_a);
					map[c] = mirror_path(map[c]);
				}
			}
		}
					
		profile[c] = malloc(sizeof(float)*64*(len_a+len_b+2));
		profile[c] = update2(profa,profb,profile[c],map[c],aln->nsip[a],aln->nsip[b],internal_gap_weight);
		
		
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
		free(profa);
		free(profb);
	}
	
	free(profile[numprofiles-1]);
	free(profile);
	
	dp_matrix_free(dp);
	for (i = 32;i--;){
		free(submatrix[i]);
	}
	free(submatrix);
	return map;
}*/
/*
int** feature_alignment(struct alignment* aln,int* tree,int**submatrix, int** map,struct feature_matrix* fm)
{
	struct dp_matrix *dp = 0;
	int i,j,g,a,b,c;
	int len_a;
	int len_b;
	int** profile = 0;
	int* profa = 0;
	int* profb = 0;
	
	int** fprofile = 0;
	int* fprofa = 0;
	int* fprofb = 0;
	
	profile = malloc(sizeof(int*)*numprofiles);
	for ( i = 0;i< numprofiles;i++){
		profile[i] = 0;
	}
	
	fprofile = malloc(sizeof(int*)*numprofiles);
	for ( i = 0;i< numprofiles;i++){
		fprofile[i] = 0;
	}

	map = malloc(sizeof(int*)*numprofiles);
	for ( i = 0;i < numprofiles;i++){
		map[i] = 0;
	}
	
	dp = dp_matrix_alloc(dp,511,511);
	c = numseq;
	//if(!param->dna){
	for (i = 0; i < (numseq-1);i++){
		a = tree[i*3];
		b = tree[i*3+1];
		c = tree[i*3+2];
		fprintf(stderr,"Aligning:%d %d->%d\n",a,b,c);
		len_a = aln->sl[a];
		len_b = aln->sl[b];
		dp = dp_matrix_realloc(dp,len_a,len_b);
	
		map[c] = malloc(sizeof(int) * (len_a+len_b+2));
		for (j = len_a+len_b+2;j--;){
			map[c][j] = 0;
		}
		if (a < numseq){
			profile[a] = make_profile(profile[a],aln->s[a],len_a,submatrix);
		//	fprintf(stderr,"Making feature profile for %d	(%s)\n",a,aln->sn[a]);
			fprofile[a] = make_feature_profile(fprofile[a],aln->ft[a],len_a,fm);		
		}
		if (b < numseq){
			profile[b] = make_profile(profile[b],aln->s[b],len_b,submatrix);
		//	fprintf(stderr,"Making feature profile for %d	(%s)\n",b,aln->sn[b]);
			fprofile[b] = make_feature_profile(fprofile[b],aln->ft[b],len_b,fm);
		}
		//profa = profile[a];
		//profb = profile[b];
		profa = profile[a]+64;
		profb = profile[b]+64;
	
		fprofa = fprofile[a];
		fprofb = fprofile[b];
	
		set_gap_penalties(profile[a],len_a,aln->nsip[b]);
		set_gap_penalties(profile[b],len_b,aln->nsip[a]);

		if (len_a > len_b){
			map[c] = fpp_dyn(map[c],dp,profa,profb,fprofa,fprofb,len_a,len_b,fm->mdim,fm->stride);
		}else{
			map[c] = fpp_dyn(map[c],dp,profb,profa,fprofb,fprofa,len_b,len_a,fm->mdim,fm->stride);
			map[c] = mirror_path(map[c]);
		}
					
		profile[c] = malloc(sizeof(int)*64*(len_a+len_b+2));
		profile[c] = update(profile[a],profile[b],profile[c],map[c],aln->nsip[a],aln->nsip[b]);
		
		fprofile[c] = malloc(sizeof(int)*fm->stride*(len_a+len_b+2));
		fprofile[c] = feature_update(fprofa,fprofb,fprofile[c],map[c],fm->stride); 
		
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
		
		free(fprofa);
		free(fprofb);
		
	}
	
	free(profile[numprofiles-1]);
	free(profile);
	
	free(fprofile[numprofiles-1]);
	free(fprofile );
	
	dp_matrix_free(dp);
	for (i = 32;i--;){
		free(submatrix[i]);
	}
	free(submatrix);
	free_feature_matrix(fm);
	return map;
}*/

struct ntree_data* ntree_sub_alignment(struct ntree_data* ntree_data,int* tree,int num)
{
	struct dp_matrix *dp = 0;
	struct alignment* aln = 0;
	int i,j,g,a,b,c;
	int len_a;
	int len_b;
	float** local_profile = 0;
	float* profa = 0;
	float* profb = 0;

	int** local_map = 0;
	int* local_sl = 0;
	int* local_nsip = 0;
	int** local_sip = 0;
	
	int* which_to_alloc = 0;
	
	aln = ntree_data->aln;
	
	which_to_alloc = malloc(sizeof(int*)*numprofiles);
	for ( i = 0;i< numprofiles;i++){
		which_to_alloc[i] = 0;
	}
	
	local_profile = malloc(sizeof(float*)*numprofiles);
	local_sl = malloc(sizeof(int)*numprofiles);
	local_nsip = malloc(sizeof(int)*numprofiles);
	local_sip = malloc(sizeof(int*)*numprofiles);
	
	
	for (i = 0; i < num-1;i++){
		a = tree[i*3+1];
		if(!which_to_alloc[a]){
			which_to_alloc[a] = 1;
		}
		b = tree[i*3+2];
		if(!which_to_alloc[b]){
			which_to_alloc[b] = 1;
		}
		c = tree[i*3+3];
		if(!which_to_alloc[c]){
			which_to_alloc[c] = 2;
		}
	}
	//for ( i = 0;i< numprofiles;i++){
	//	fprintf(stderr,"alloc?:%d	%d\n",i,which_to_alloc[i]);
	//}
	
//	exit(0);
	for ( i = 0;i< numprofiles;i++){
		if(which_to_alloc[i] == 1){
			local_profile[i] = ntree_data->profile[i];
			local_sl[i] = aln->sl[i];
			local_nsip[i] = aln->nsip[i];
			local_sip[i] =  malloc(sizeof(int*)*aln->nsip[i]);
			for(j = 0;j < aln->nsip[i];j++){
				local_sip[i][j] = aln->sip[i][j];
			}
		}else{
			local_profile[i] = 0;
			local_sl[i] = 0;
			local_nsip[i] = 0;
			local_sip[i] = 0;
		}
	}
	/*
	for ( i = 0;i< numprofiles;i++){
		local_profile[i] = ntree_data->profile[i];
		local_sl[i] = aln->sl[i];
		local_nsip[i] = aln->nsip[i];
		if(aln->sip[i]){
			fprintf(stderr,"Allocing..:%d\n",aln->nsip[i]);
			local_sip[i] =  malloc(sizeof(int*)*aln->nsip[i]);
			for(j = 0;j < aln->nsip[i];j++){
				local_sip[i][j] = aln->sip[i][j];
			}
		}else{
			local_sip[i] = 0;
		}
	}*/

	local_map = malloc(sizeof(int*)*numprofiles);
	for ( i = 0;i < numprofiles;i++){
		local_map[i] = 0;
	}
	
	
	dp = dp_matrix_alloc(dp,511,511);
	c = numseq;
	for (i = 0; i < num-1;i++){
		a = tree[i*3+1];
		b = tree[i*3+2];
		c = tree[i*3+3];
	//	fprintf(stderr,"Aligning:%d %d->%d\n",a,b,c);
		len_a = local_sl[a];
		len_b = local_sl[b];
		dp = dp_matrix_realloc(dp,len_a,len_b);
	
		local_map[c] = malloc(sizeof(int) * (len_a+len_b+2));
		for (j = len_a+len_b+2;j--;){
			local_map[c][j] = 0;
		}
		if (a < numseq){
			local_profile[a] = make_profile(local_profile[a],aln->s[a],len_a,ntree_data->submatrix);
		}
		if (b < numseq){
			local_profile[b] = make_profile(local_profile[b],aln->s[b],len_b,ntree_data->submatrix);
		}
		profa = local_profile[a];
		profb = local_profile[b];
	
		set_gap_penalties(profa,len_a,local_nsip[b],0,local_nsip[a]);
		set_gap_penalties(profb,len_b,local_nsip[a],0,local_nsip[b]);

		if(local_nsip[a] == 1){
			if(local_nsip[b] == 1){
				local_map[c] = ss_dyn(ntree_data->submatrix,local_map[c],dp,aln->s[a],aln->s[b],len_a,len_b);
			}else{
				local_map[c] = ps_dyn(local_map[c],dp,profb,aln->s[a],len_b,len_a,local_nsip[b]);
				local_map[c] = mirror_path(local_map[c]);
			}
		}else{
			if(local_nsip[b] == 1){
				local_map[c] = ps_dyn(local_map[c],dp,profa,aln->s[b],len_a,len_b,local_nsip[a]);
			}else{
				if (len_a > len_b){
					local_map[c] = pp_dyn(local_map[c],dp,profa,profb,len_a,len_b);
				}else{
					local_map[c] = pp_dyn(local_map[c],dp,profb,profa,len_b,len_a);
					local_map[c] = mirror_path(local_map[c]);
				}
			}
		}
			
		local_profile[c] = malloc(sizeof(float)*64*(len_a+len_b+2));
		local_profile[c] = update(profa,profb,local_profile[c],local_map[c],local_nsip[a],local_nsip[b]);
	
		local_sl[c] = local_map[c][0];
	
		local_nsip[c] = local_nsip[a] + local_nsip[b];
		local_sip[c] = malloc(sizeof(int)*(local_nsip[a] + local_nsip[b]));
		g =0;
		for (j = local_nsip[a];j--;){
			local_sip[c][g] = local_sip[a][j];
			g++;
		}
		for (j = local_nsip[b];j--;){
			local_sip[c][g] = local_sip[b][j];
			g++;
		}
	//	free(profa);
	//	free(profb);
	}
	
	if(ntree_data->profile[c]){
		if(ntree_data->map[c][ntree_data->map[c][0]+2]  < local_map[c][local_map[c][0]+2]){
			fprintf(stderr,"%d\n",local_map[c][local_map[c][0]+2]);
			//remove old map,profile,etc..
			for (i = 0; i < num-1;i++){
				c = tree[i*3+3];
				free(ntree_data->map[c]);
				free(ntree_data->profile[c]);
				free(aln->sip[c]);
				ntree_data->map[c] = malloc(sizeof(int)*(local_map[c][0]+3));
				for (j = 0; j < local_map[c][0]+3;j++){
					ntree_data->map[c][j] = local_map[c][j];
				}
				aln->sip[c] = malloc(sizeof(int)*local_nsip[c]);
				aln->nsip[c] = local_nsip[c];
				for (j = 0; j < local_nsip[c];j++){
					aln->sip[c][j] = local_sip[c][j];
				}
				aln->sl[c] = local_sl[c];
				
			}
			ntree_data->profile[c] = malloc(sizeof(int)*64*(aln->sl[c]+1));
			for (i = 0; i < (64*(aln->sl[c]+1));i++){
				ntree_data->profile[c][i] = local_profile[c][i];
			}
			ntree_data->tree[0] -= (tree[0]-1);
			for (j = 1; j < tree[0];j++){
				ntree_data->tree[ntree_data->tree[0]+j-1] = tree[j];
			}
			ntree_data->tree[0] += (tree[0]-1);

		}else{
			fprintf(stderr,"no improvement\n");
		}
	}else{
		fprintf(stderr,"%d\n",local_map[c][local_map[c][0]+2]);
		for (i = 0; i < num-1;i++){
			c = tree[i*3+3];
			ntree_data->map[c] = malloc(sizeof(int)*(local_map[c][0]+3));
			for (j = 0; j < local_map[c][0]+3;j++){
				ntree_data->map[c][j] = local_map[c][j];
			}

			aln->sip[c] = malloc(sizeof(int)*local_nsip[c]);
			aln->nsip[c] = local_nsip[c];
			for (j = 0; j < local_nsip[c];j++){
				aln->sip[c][j] = local_sip[c][j];
			}
			aln->sl[c] = local_sl[c];
		}
		ntree_data->profile[c] = malloc(sizeof(int)*64*(aln->sl[c]+1));
		for (i = 0; i < (64*(aln->sl[c]+1));i++){
			ntree_data->profile[c][i] = local_profile[c][i];
		}
		for (j = 1; j < tree[0];j++){
			ntree_data->tree[ntree_data->tree[0]+j-1] = tree[j];
		}
		ntree_data->tree[0] += tree[0]-1;
	}

	for ( i = 0;i< numprofiles;i++){
		if(which_to_alloc[i] == 1){
			free(local_sip[i]);
			if(i < numseq){
				free(local_profile[i]);	
			}
		}
		if(which_to_alloc[i] == 2){
			free(local_profile[i]);
			free(local_map[i]);
			free(local_sip[i]);
		}
		
	}

	free(which_to_alloc);
	free(local_map);
	free(local_sip);
	free(local_nsip);
	free(local_profile);
	free(local_sl);
	
	dp_matrix_free(dp);
	return ntree_data;
}

struct ntree_data* ntree_alignment(struct ntree_data* ntree_data)
{
	int i;
	ntree_data->profile = malloc(sizeof(float*)*numprofiles);
	for ( i = 0;i< numprofiles;i++){
		ntree_data->profile[i] = 0;
	}
	
	ntree_data->map = malloc(sizeof(int*)*numprofiles);
	for ( i = 0;i < numprofiles;i++){
		ntree_data->map[i] = 0;
	}

	ntree_data =  alignntree(ntree_data,ntree_data->realtree);
	
	for ( i = 0;i< numprofiles;i++){
		if(ntree_data->profile[i]){
			free(ntree_data->profile[i]);
		}
	}
	free(ntree_data->profile);
	
	for (i = 32;i--;){
		free(ntree_data->submatrix[i]);
	}
	free(ntree_data->submatrix);
	free_real_tree(ntree_data->realtree);
	return ntree_data;
}

