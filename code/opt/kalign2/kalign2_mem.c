/*
	kalign2_mem.c 
	
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

void freesimpletree(struct tree_node* p)
{
	if(p->left){
		freesimpletree(p->left);
	}
	if(p->right){
		freesimpletree(p->right);
	}
	free(p);
}

void free_real_tree(struct aln_tree_node* p)
{
	int i = 0;
	while(p->links[i]){
		free_real_tree(p->links[i]);
		i++;	
	}
	free(p->internal_lables);
	free(p->links);
	free(p);
}


void free_feature_matrix(struct feature_matrix* fm)
{
	int i;
	for (i = 0;i < fm->mdim;i++){
		free(fm->m[i]);
	}
	free(fm->m);
	free(fm);
}


void free_utf(struct utype_ufeat* utf)
{
	free_ft(utf->t);
	free_ft(utf->f);
	free(utf);
}

/*#ifndef MEMORY 
void* malloc(int size)
{
	void* p;
	p = (void*)malloc(size);
	
	if (!p){
		fprintf(stderr,"Out of memory!\n");
		exit(0);
	}
	return p;
}
#endif*/
struct names* names_alloc(struct names* n)
{
	int i;
	n = malloc(sizeof(struct names));
	n->start = malloc(sizeof(int)*numseq);
	n->end = malloc(sizeof(int)*numseq);
	n->len = malloc(sizeof(int)*numseq);
	
	for (i = 0; i < numseq;i++){
		n->start[i] = 0;
 		n->end[i] = 0;//aln->lsn[i];
 		n->len[i] = 0;
	}
	return n;
}

void names_free(struct names* n)
{
	free(n->start);
	free(n->end);
	free(n->len);
	free(n);
}


struct alignment* aln_alloc(struct alignment* aln)
{
	int i;
	aln = (struct alignment*) malloc(sizeof(struct alignment));
	aln->s = malloc(sizeof(int*) * (numseq ));
	aln->seq = malloc(sizeof(char*) * (numseq ));
	aln->ft =  malloc(sizeof(struct feature* ) * (numseq));
	aln->si  =  malloc(sizeof(struct sequence_information* ) * (numseq));
	aln->sl = malloc(sizeof(unsigned int) * (numprofiles));
	aln->sip = malloc(sizeof(unsigned int*)* numprofiles);
	
	aln->nsip = malloc(sizeof(unsigned int)* numprofiles);
	aln->sn = malloc(sizeof(char*) * numseq);
	aln->lsn = malloc(sizeof(unsigned int) * numseq);
	for (i =0;i < numprofiles;i++){
		aln->sip[i] = 0;
		aln->nsip[i] = 0;
		aln->sl[i] = 0;
	}
	
	for(i =0;i < numseq;i++){
		aln->lsn[i] = 0;
		aln->ft[i] = 0;
		aln->si[i] = 0;
		aln->sip[i] = malloc(sizeof(int)*1);
		aln->nsip[i] = 1;
		aln->sip[i][0] = i;
	}
	return aln;
}


void free_aln(struct alignment* aln)
{
	int i;
	for (i = numseq;i--;){
		free(aln->s[i]);
		free(aln->seq[i]);
		free(aln->sn[i]);
	}

	if(aln->ft){
		for(i = numseq;i--;){
			free_ft(aln->ft[i]);
		}
		free(aln->ft);
	}
	if(aln->si){
		free(aln->si);
	}

	for (i = numprofiles;i--;){
		if(aln->sip[i]){
			free(aln->sip[i]);
		}
	}
	free(aln->seq);
	free(aln->s);
	free(aln->sn);
	free(aln->sl);
	free(aln->lsn);
	free(aln->sip);
	free(aln->nsip);
	free(aln);
}


void free_param(struct parameters* param)
{
	free(param->infile);
	free(param);
}

void free_ft(struct feature* n)
{
	struct feature* old_n = 0;
	 if (n != NULL){
	 	old_n = n;
	 	n= n ->next;
	 	free(old_n->type);
	 	free(old_n->note);
 	 	free(old_n);
		free_ft(n);
	}
}

struct hirsch_mem* hirsch_mem_alloc(struct hirsch_mem* hm,int x)
{

	// a=((typeof(a))(((int)(((void *)malloc(c+15))+15))&-16)). 
	hm = (struct hirsch_mem *) malloc(sizeof(struct hirsch_mem));
	hm->starta = 0;
	hm->startb = 0;
	hm->enda = 0;
	hm->endb = 0;
	hm->size = x;
	hm->len_a = 0;
	hm->len_b = 0;
	hm->f = malloc(sizeof(struct states)* (x+1));
	hm->b = malloc(sizeof(struct states)* (x+1));
	return hm;
}

struct hirsch_mem* hirsch_mem_realloc(struct hirsch_mem* hm,int x)
{
	hm->starta = 0;
	hm->startb = 0;
	hm->enda = 0;
	hm->endb = 0;
	hm->len_a = 0;
	hm->len_b = 0;
	hm->size = x;
	hm->f = realloc(hm->f,sizeof(struct states)* (x+1));
	hm->b = realloc(hm->b,sizeof(struct states)* (x+1));
	return hm;
}

void hirsch_mem_free(struct hirsch_mem* hm)
{
	free(hm->f);
	free(hm->b);
	free(hm);
}


struct dp_matrix* dp_matrix_realloc(struct dp_matrix *dp,int x,int y)
{
	int i;
	if (x > y){
		y = x;
	}else{
		x =y;
	}
	if ( x > dp->x || y > dp->y){
		//printf("REALLOCING:%d-%d	%d-%d\n",x,y,dp->x,dp->y);
		i = 1;
		while (i <= y){
			i <<= 1;
		//	printf("i:%d	y:%d\n",i,y);
		}
		y = i-1;
		i = 1;
		while (i <= x){
			i <<= 1;
			//printf("i:%d	y:%d\n",i,y);
		}
		x = i-1;
		//printf("NEWX:%d	NEWY:%d\n",x,y);
		dp->s = realloc(dp->s,sizeof(struct states)* (y+1));
		dp->tb  = (char**) realloc (dp->tb,sizeof (char*)*(x+1));
		dp->tb_mem = (void*) realloc(dp->tb_mem,sizeof(char) * (x+1) * (y+1));
		dp->tb[0] = (char*) dp->tb_mem;
  		for (i = 1; i <= x; i++){
			dp->tb[i] = dp->tb[0] +(i*(y+1));
		}
		dp->x = x;
		dp->y = y;
	}
	return dp;
}

struct dp_matrix* dp_matrix_alloc(struct dp_matrix *dp,int x,int y)
{
	int i;
	dp = (struct dp_matrix *) malloc(sizeof(struct dp_matrix));
	dp->x = x;
	dp->y = y;
	dp->s = malloc(sizeof(struct states)* (y+1));
	dp->tb = (char**) malloc(sizeof(char*) * (x+1));
	dp->tb_mem = (void *) malloc(sizeof(char) * (x+1) * (y+1));
	dp->tb[0] = (char*) dp->tb_mem;
	for ( i = 1; i <= x;i++){
		dp->tb[i] = dp->tb[0] +(i*(y+1));
	}
	return dp;
}

void dp_matrix_free(struct dp_matrix *dp)
{
	free(dp->s);
	free(dp->tb);
	free(dp->tb_mem);
	free(dp);
}
