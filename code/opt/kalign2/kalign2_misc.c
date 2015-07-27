/*
	kalign2_misc.c 
	
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



#include <string.h>
#include <ctype.h>
#include "kalign2.h"

void print_tree(struct aln_tree_node* p,struct alignment* aln,char* outfile)
{
	FILE *fout = NULL;
	if ((fout = fopen(outfile, "w")) == NULL){
		fprintf(stderr,"can't open output\n");
		exit(0);
	}
	//fprintf(stderr,"\n\n%s\n",outfile);
	/*if(byg_start("xml",outfile) != -1){
		fprintf(fout,"<?xml version=\"1.0\" encoding=\"UTF-8\"?> <phyloxml xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"http://www.phyloxml.org/schema/phyloxml.xsd\"><phylogeny>\n");
		
		print_phyloxml_tree(p,aln,fout);
		
		fprintf(fout,"</phylogeny></phyloxml>\n");
		
	}else{*/
		print_newick_tree(p,aln,fout);
		fprintf(fout,";");
	//}
	fclose(fout);
}



void print_newick_tree(struct aln_tree_node* p,struct alignment* aln,FILE *fout)
{
	int j;
		
	if(p->links[0]){
	
		fprintf(fout,"(");
		print_newick_tree(p->links[0],aln,fout);
	}
	if(p->num < numseq){
		//If you want to print the actual names of the sequences 
		for (j = 0; j < aln->lsn[p->num];j++){
			if(isspace((int)aln->sn[p->num][j])){
				fprintf(fout,"_");
			}else{
				fprintf(fout,"%c",aln->sn[p->num][j]);
			}
		}
		//If you want to print the just the number of the sequence
		//fprintf(stdout,"%d",p->num);
		
	}else{
		fprintf(fout,",");
	}
	if(p->links[1]){
		print_newick_tree(p->links[1],aln,fout);
		fprintf(fout,")");
	}
}

void print_phyloxml_tree(struct aln_tree_node* p,struct alignment* aln,FILE *fout)
{
	int j;
		
	if(p->links[0]){
	
		fprintf(fout,"<clade>\n");
		print_phyloxml_tree(p->links[0],aln,fout);
	}
	if(p->num < numseq){
		//If you want to print the actual names of the sequences 
		fprintf(fout,"<clade>\n<name>");
		for (j = 0; j < aln->lsn[p->num];j++){
			fprintf(fout,"%c",aln->sn[p->num][j]);
		}
		fprintf(fout,"</name>\n</clade>\n");
		//If you want to print the just the number of the sequence
		//fprintf(stdout,"%d",p->num);
		
	}else{
		//fprintf(fout,",");
	}
	if(p->links[1]){
		print_phyloxml_tree(p->links[1],aln,fout);
		fprintf(fout,"</clade>\n");
	}
}

struct alignment* sort_sequences(struct alignment* aln,int* tree,char* sort)
{
	int i,j,a,b,c;
	int choice = 0;
	
	if(sort){
		if (byg_start("input",sort) != -1){
			choice = 0;
		}else if (byg_start("tree",sort) != -1){
			choice = 1;
		}else if (byg_start("gaps",sort) != -1){
			choice = 2;
		}else{
			choice = 3;
		}
	}
	//fprintf(stderr,"COICE:%d\n",choice);
	switch(choice){
		case 0:
			for (i = 0; i < numseq;i++){
				aln->nsip[i] = i;
			}
			break;
		case 1:
			c = 0;
			for (i = 0; i < (numseq-1)*3;i +=3){
				//fprintf(stderr,"TREE %d	%d	%d\n",tree[i],tree[i+1],tree[i+2]);
				if(tree[i]  < numseq){
					aln->nsip[c] = tree[i];
					c++;
				}
				if(tree[i+1]  < numseq){
					aln->nsip[c] = tree[i+1];
					c++;
				}
			}
			break;
		case 2:
			for (i = 0; i < numseq;i++){
				a = 1000000;
				b = -1;
				for (j =0; j<numseq;j++){
					if(aln->nsip[j] < a){
						a = aln->nsip[j];
						b = j;
					}
				}	
				tree[i] = b;
				aln->nsip[b] = 1000000;
			}
			for (i = 0; i < numseq;i++){
				aln->nsip[i] = tree[i];
			}
			break;
		case 3:
			aln = sort_in_relation(aln,sort);
			break;
		default:
			for (i = 0; i < numseq;i++){
				aln->nsip[i] = i;
			}
			break;
	}
	
	/*for (i = 0; i < numseq;i++){
		fprintf(stderr,"%d\n",aln->nsip[i]);
	}*/
	
	return aln;
}

struct alignment* sort_in_relation(struct alignment* aln,char* sort)
{
	int i,j,c;
	int target = -1;
	int id = 0;
	int positions = 0;
	int posa = 0;
	int posb = 0;
	for (i = 0; i < numseq;i++){
		if (byg_start(sort,aln->sn[i]) != -1){
			target = i;
			aln->sip[i][0] = 1000;
			break;
		}
	}
	if(target == -1){
		target = 0;
		aln->sip[0][0] = 1000;
	}
	for (i = 0; i < numseq;i++){
		if(i != target){
			posa = 0;
			posb =0;
			c = 0;
			id = 0;
			positions = 0;
			for (j = 0; j < aln->sl[i];j++){
				posa += aln->s[i][j]+1;
				while(posa > posb){
					posb += aln->s[target][c]+1;
					c++;
				}
				if(posa == posb){
					if((int) aln->seq[i][j] == (int) aln->seq[target][c-1]){
						id += 1000;
					}
					positions += 1;
				}
			}
			if(positions){
				aln->sip[i][0] = id/positions;
			}else{
				aln->sip[i][0] = 0;
			}
		}
	}
	for (i = 0; i < numseq;i++){
		aln->nsip[i] = i;
	}
	quickSort(aln, numseq);
	return aln;
}

void quickSort(struct alignment* aln, int array_size)
{
	q_sort(aln, 0, array_size - 1);
}

void q_sort(struct alignment* aln, int left, int right)
{
	int pivot, l_hold, r_hold;
	int pivot2;
	l_hold = left;
	r_hold = right;
	pivot2 = aln->nsip[left];
	pivot = aln->sip[left][0];// numbers[left];
	while (left < right){
		while ((aln->sip[right][0] <= pivot) && (left < right)){
			right--;
		}
		if (left != right){
			aln->sip[left][0] = aln->sip[right][0];
			aln->nsip[left] = aln->nsip[right];
			left++;
		}
		while ((aln->sip[left][0] >= pivot) && (left < right)){
			left++;
		}
		if (left != right){
			aln->sip[right][0] = aln->sip[left][0];
			aln->nsip[right] = aln->nsip[left];
			right--;
		}
	}
	aln->sip[left][0] = pivot;
	aln->nsip[left] = pivot2;
	pivot = left;
	left = l_hold;
	right = r_hold;
	if (left < pivot){
		q_sort(aln, left, pivot-1);
	}
	if (right > pivot){
		q_sort(aln, pivot+1, right);
	}
}

int* readtree(struct aln_tree_node* p,int* tree)
{
	if(p->links[0]){
		tree = readtree(p->links[0],tree);
	}
	if(p->links[1]){
		tree = readtree(p->links[1],tree);
	}
	
	if(p->links[0]){
		if(p->links[1]){
			tree[tree[0]] = p->links[0]->num;
			tree[tree[0]+1] = p->links[1]->num;
			tree[tree[0]+2] = p->num;
			tree[0] +=3;
			free(p->links[0]->internal_lables);
			free(p->links[0]->links);
			free(p->links[0]);
			free(p->links[1]->internal_lables);
			free(p->links[1]->links);
			free(p->links[1]);
		}
	}
	return tree;
}

struct alignment* make_dna(struct alignment* aln)
{

	//int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
	int i,j;
	int* p;
	
	for(i = 0;i < numseq;i++){
		p = aln->s[i];
		for (j = 0; j < aln->sl[i];j++){
			switch(p[j]){
				case 2: //C
					p[j] = 1;
					break;
				case 6: //G
					p[j] = 2;
					break;
				case 17: //T  or U 
					p[j] = 3;
					break;
				case 12: // N
					p[j] = 4;
					break;
				case 20: // X
					p[j] = 4;
					break;
				case 23://O whatever that is...
					p[j] = 4;
					break;
			}
		//	printf("%d\n",p[j]);
		}
	}
	return aln;
}

float** read_matrix(float** subm,struct parameters* param)
{
	int i,j;
	int m_pos = 0;
	short *matrix_pointer = 0;
	short blosum50mt[]={
  5,
 -2,  5,
 -1, -3, 13,
 -2,  5, -4,  8,
 -1,  1, -3,  2,  6,
 -3, -4, -2, -5, -3,  8,
  0, -1, -3, -1, -3, -4,  8,
 -2,  0, -3, -1,  0, -1, -2, 10,
 -1, -4, -2, -4, -4,  0, -4, -4,  5,
 -1,  0, -3, -1,  1, -4, -2,  0, -3,  6,
 -2, -4, -2, -4, -3,  1, -4, -3,  2, -3,  5,
 -1, -3, -2, -4, -2,  0, -3, -1,  2, -2,  3,  7,
 -1,  4, -2,  2,  0, -4,  0,  1, -3,  0, -4, -2,  7,
 -1, -2, -4, -1, -1, -4, -2, -2, -3, -1, -4, -3, -2, 10,
 -1,  0, -3,  0,  2, -4, -2,  1, -3,  2, -2,  0,  0, -1,  7,
 -2, -1, -4, -2,  0, -3, -3,  0, -4,  3, -3, -2, -1, -3,  1,  7,
  1,  0, -1,  0, -1, -3,  0, -1, -3,  0, -3, -2,  1, -1,  0, -1,  5,
  0,  0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  2,  5,
  0, -4, -1, -4, -3, -1, -4, -4,  4, -3,  1,  1, -3, -3, -3, -3, -2,  0,  5,
 -3, -5, -5, -5, -3,  1, -3, -3, -3, -3, -2, -1, -4, -4, -1, -3, -4, -3, -3, 15,
 -1, -1, -2, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1, -2, -1, -1, -1,  0, -1, -3, -1,
 -2, -3, -3, -3, -2,  4, -3,  2, -1, -2, -1,  0, -2, -3, -1, -1, -2, -2, -1,  2, -1,  8,
 -1,  2, -3,  1,  5, -4, -2,  0, -3,  1, -3, -1,  0, -1,  4,  0,  0, -1, -3, -2, -1, -2,  5};
 
	short blosum62mt[]={
  40,
 -20,  40,
  0, -30,  90,
 -20,  40, -30,  60,
 -10,  10, -40,  20,  50,
 -20, -30, -20, -30, -30,  60,
  0, -10, -30, -10, -20, -30,  60,
 -20,  0, -30, -10,  0, -10, -20,  80,
 -10, -30, -10, -30, -30,  0, -40, -30,  40,
 -10,  0, -30, -10,  10, -30, -20, -10, -30,  50,
 -10, -40, -10, -40, -30,  0, -40, -30,  20, -20,  40,
 -10, -30, -10, -30, -20,  0, -30, -20,  10, -10,  20,  50,
 -20,  30, -30,  10,  0, -30,  0,  10, -30,  0, -30, -20,  60,
 -10, -20, -30, -10, -10, -40, -20, -20, -30, -10, -30, -20, -20,  70,
 -10,  0, -30,  0,  20, -30, -20,  0, -30,  10, -20,  0,  0, -10,  50,
 -10, -10, -30, -20,  0, -30, -20,  0, -30,  20, -20, -10,  0, -20,  10,  50,
  10,  0, -10,  0,  0, -20,  0, -10, -20,  0, -20, -10,  10, -10,  0, -10,  40,
  0, -10, -10, -10, -10, -20, -20, -20, -10, -10, -10, -10,  0, -10, -10, -10,  10,  50,
  0, -30, -10, -30, -20, -10, -30, -30,  30, -20,  10,  10, -30, -20, -20, -30, -20,  0,  40,
 -30, -40, -20, -40, -30,  10, -20, -20, -30, -30, -20, -10, -40, -40, -20, -30, -30, -20, -30, 110,
  0, -10, -20, -10, -10, -10, -10, -10, -10, -10, -10, -10, -10, -20, -10, -10,  0,  0, -10, -20, -10,
 -20, -30, -20, -30, -20,  30, -30,  20, -10, -20, -10, -10, -20, -30, -10, -20, -20, -20, -10,  20, -10,  70,
 -10,  10, -30,  10,  40, -30, -20,  0, -30,  10, -30, -10,  0, -10,  30,  0,  0, -10, -20, -30, -10, -20,  40};
	
	short gon250mt[]={
	24,
	0,   0,
	5,   0, 115,
	-3,   0, -32,  47,
	 0,   0, -30,  27,  36,
	-23,   0,  -8, -45, -39,  70,
	5,   0, -20,   1,  -8, -52,  66,
	-8,   0, -13,   4,   4,  -1, -14,  60,
	-8,   0, -11, -38, -27,  10, -45, -22,  40,
	-4,   0, -28,   5,  12, -33, -11,   6, -21,  32,
	-12,   0, -15, -40, -28,  20, -44, -19,  28, -21,  40,
	-7,   0,  -9, -30, -20,  16, -35, -13,  25, -14,  28,  43,
	-3,   0, -18,  22,   9, -31,   4,  12, -28,   8, -30, -22,  38,
	 3,   0, -31,  -7,  -5, -38, -16, -11, -26,  -6, -23, -24,  -9,  76,
	-2,   0, -24,   9,  17, -26, -10,  12, -19,  15, -16, -10,   7,  -2,  27,
	-6,   0, -22,  -3,   4, -32, -10,   6, -24,  27, -22, -17,   3,  -9,  15,  47,
	11,   0,   1,   5,   2, -28,   4,  -2, -18,   1, -21, -14,   9,   4,   2,  -2,  22,
	 6,   0,  -5,   0,  -1, -22, -11,  -3,  -6,   1, -13,  -6,   5,   1,   0,  -2,  15,  25,
	 1,   0,   0, -29, -19,   1, -33, -20,  31, -17,  18,  16, -22, -18, -15, -20, -10,   0,  34,
	-36,   0, -10, -52, -43,  36, -40,  -8, -18, -35,  -7, -10, -36, -50, -27, -16, -33, -35, -26, 142,
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
	-22,   0,  -5, -28, -27,  51, -40,  22,  -7, -21,   0,  -2, -14, -31, -17, -18, -19, -19, -11,  41,   0,  78,
	 0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0};
	if(param->sub_matrix){
		if(byg_start(param->sub_matrix,"blosum62BLOSUM62") != -1){
			matrix_pointer = blosum62mt;
			//m_pos = 0;
			//for (i = 0;i < 23;i++){
			//	for (j = 0;j <= i;j++){
			//		matrix_pointer[m_pos] = matrix_pointer[m_pos] * 10;
			//		m_pos++;
			//	}
			//}
			gpo = 55;
			gpe = 8;
			tgpe = 1;
		}
		if(byg_start(param->sub_matrix,"blosum50BLOSUM50") != -1){
			matrix_pointer = blosum50mt;
			m_pos = 0;
			for (i = 0;i < 23;i++){
				for (j = 0;j <= i;j++){
					matrix_pointer[m_pos] = matrix_pointer[m_pos] * 10;
					m_pos++;
				}
			}
			gpo = 55;
			gpe = 8;
			tgpe = 1;
		}
		//vogt....
		
	}else{
		if(!param->dna){
			// gpo:5.494941        gpe:0.852492        tgpe:0.442410       bonus: 3.408872     z-cutoff: 58.823309 -> 0.829257 accuracy on bb3 
			gpo = 54.94941;
			gpe = 8.52492;
			tgpe = 4.42410;
			
			//gpo = 54;
			//gpe = 8;
			//tgpe = 4;
			//-gpo 10.9898        -gpe 0.852492      -tgpe  0.442410    -bonus    0.2   -zcutoff     58.823309  
		//	param->secret = 0.2;
			matrix_pointer = gon250mt;
		}else{
			//gpo = 400;
		//	gpe =  30;
			//tgpe = 30;
			
			//param->gpo = 43.4;
			//param->gpe = 3.94;
			//param->tgpe = 29.26;
			
			//gpo = 43.4 *5;
			gpo = 217;
			gpe = 39.4;
			tgpe =  292.6;
			//param->secret = 28.3;
			param->zlevel = 61.08;
			param->internal_gap_weight = 49.14;
			
		}
	}
	if(param->gpo!= -1){
		//param->gpo *= 5;
		gpo = param->gpo;
	}
	if(param->gpe != -1){
		//param->gpe *= 10;
		gpe = param->gpe;
	}
	if(param->tgpe != -1){
		//param->tgpe *= 10;
		tgpe = param->tgpe;
	}

//	if(param->secret != -1){
//		//param->secret *= 10;
//	}else{
	if(param->secret == -1){
		if(!param->dna){
			param->secret = 0.2;
		}else{
			param->secret = 283.0;
		}
	}

	
	//fprintf(stderr,"%d	%d	%d	%d\n",gpo,gpe,tgpe,	 (int)param->secret);
	subm = malloc(sizeof (float*) * 32);
	for (i = 32;i--;){
		subm[i] = malloc(sizeof(float) * 32);
		for (j = 32;j--;){
			subm[i][j] = param->secret;//0;//gpe << 1;//-5;// better on Balibase
		}
	}
	if(param->dna){
		/*subm[0][0] += 10;
		subm[0][1] += 6;
		subm[1][0] += 6;
		subm[1][1] += 10;
		subm[2][2] += 10;
		subm[2][3] += 6;
		subm[3][2] += 6;
		subm[3][3] += 10;*/
//		     A    C    G    T    .    N
//	A   91 -114  -31 -123    0  -43
		subm[0][0] += 91;
		subm[0][1] += -114;
		subm[0][2] += -31;
		subm[0][3] += -123;

//	C -114  100 -125  -31    0  -43
		subm[1][0] += -114;
		subm[1][1] += 100;
		subm[1][2] += -125;
		subm[1][3] += -31;

//	G  -31 -125  100 -114    0  -43
		subm[2][0] += -31;
		subm[2][1] += -125;
		subm[2][2] += 100;
		subm[2][3] += -114;

//	T -123  -31 -114   91    0  -43
		subm[3][0] += -123;
		subm[3][1] += -31;
		subm[3][2] += -114;
		subm[3][3] += 91;

//	.    0    0    0    0    0    0
//	N  -43  -43  -43  -43    0  -43


		/*for (i = 0; i < 4;i++){
			for (j = 0;j < 4;j++){
				if(i == j){
					subm[i][j] += 1;
				}else{
					subm[i][j] -= 3;
				}
			}
		}*/
		
	}else{

		m_pos = 0;
		for (i = 0;i < 23;i++){
			for (j = 0;j <= i;j++){
				if (i == j){
				//	subm[i][j] += blosum62mt[m_pos]*10;
					subm[i][j] += matrix_pointer[m_pos];
				}else{
				//	subm[i][j] += blosum62mt[m_pos]*10;
				//	subm[j][i] += blosum62mt[m_pos]*10;
					subm[i][j] += matrix_pointer[m_pos];
					subm[j][i] += matrix_pointer[m_pos];
				}
				m_pos++;
			}
		}
		/*for (i = 0; i < 23;i++){
			for (j = 0; j < 23;j++){
				fprintf(stderr,"%d ",subm[i][j]);
			} 
			fprintf(stderr,"\n");
		}
		fprintf(stderr,"\n");*/
	}
	return subm;
}

struct alignment* make_seq(struct alignment* aln,int a,int b,int* path)
{
	int c;
	int i;
	int posa = 0;
	int posb = 0;

	int* gap_a = 0;
	int* gap_b = 0;

	gap_a = malloc ((path[0]+1)*sizeof(int));
	gap_b = malloc ((path[0]+1)*sizeof(int));

	for (i = path[0]+1;i--;){
		gap_a[i] = 0;
		gap_b[i] = 0;
	}
	c = 1;
	while(path[c] != 3){
		if (!path[c]){
			posa++;
			posb++;
		}
		if (path[c] & 1){
			gap_a[posa] += 1;
			posb++;
		}
		if (path[c] & 2){
			gap_b[posb] += 1;
			posa++;
		}
		c++;
	}	
	for (i = aln->nsip[a];i--;){
		update_gaps(aln->sl[aln->sip[a][i]],aln->s[aln->sip[a][i]],path[0],gap_a);
	}
	for (i = aln->nsip[b];i--;){
		update_gaps(aln->sl[aln->sip[b][i]],aln->s[aln->sip[b][i]],path[0],gap_b);
	}
	free(gap_a);
	free(gap_b);
	free(path);
	return aln;
}


void update_gaps(int old_len,int*gis,int new_len,int *newgaps)
{
	unsigned int i,j;
	int add = 0;
	int rel_pos = 0;
	for (i = 0; i <= old_len;i++){
		add = 0;
		for (j = rel_pos;j <= rel_pos + gis[i];j++){
			if (newgaps[j] != 0){
				add += newgaps[j];
			}
		}
		rel_pos += gis[i]+1;
		gis[i] += add;
	}
}

int* mirror_path(int* path)
{
	int c = 1;
	while(path[c] != 3){
		if (path[c] & 1){
			path[c] += 1;
		}else if (path[c] & 2){
			path[c] -= 1;
		}
		c++;
	}
	return path;
}

struct node* insert(struct node *n, int pos)
{
        if (n == NULL){
		n = (struct node*) malloc(sizeof(struct node));
		n->next = 0;
		n->pos = pos;
	}else{
		n->next = insert(n->next,pos);
	}
        return n;
}

struct bignode* big_insert_hash(struct bignode *n,const unsigned int pos)
{
	struct bignode* p = 0;
	if(n){
		if(n->num < NODESIZE){
			n->pos[n->num] = pos;
			n->num++;
			return n;
		}else{
			p = (struct bignode*) malloc(sizeof(struct bignode));
			p->pos[0] = pos;
			p->num = 1;
			p->next = n;
		}
	}else{
		p = (struct bignode*) malloc(sizeof(struct bignode));
		p->pos[0] = pos;
		p->num = 1;
		p->next = n;
	}
	return p;
}

void big_remove_nodes(struct bignode *n)
{
	struct bignode* p;
	while(n){
		p = n;
		n = n->next;
		free(p);
	}
}

void big_print_nodes(struct bignode *n)
{
	int i;
	while(n){
		for (i = 0; i < n->num;i++){
			fprintf(stderr,"%d ",n->pos[i]);
		}
		n = n->next;
	}
}

struct node* insert_hash(struct node *n, int pos)
{
	struct node* p;
	p = (struct node*) malloc(sizeof(struct node));
	p->pos = pos;
	p->next = n;
	return p;
}

void remove_nodes(struct node *n)
{
	struct node* p;
	while(n){
		p = n;
		n = n->next;
		free(p);
	}
}
