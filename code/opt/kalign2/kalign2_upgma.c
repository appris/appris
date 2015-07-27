/*
	kalign2_upgma.c
	
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

struct aln_tree_node* real_upgma(float **dm)
{
	int i,j;
	int *as = 0;
	float max;
	int node_a = 0;
	int node_b = 0;
	int cnode = numseq;
	
	struct aln_tree_node** tree = 0;
	struct aln_tree_node* tmp = 0;
	
	as = malloc(sizeof(int)*numseq);
	for (i = numseq; i--;){
		as[i] = i+1;
	}
	
	tree = malloc(sizeof(struct aln_tree_node*)*numseq);
	for (i=0;i < numseq;i++){
		tree[i] = malloc(sizeof(struct aln_tree_node));
		tree[i]->done = 1;
		tree[i]->num = i;
		tree[i]->path = 0;
		tree[i]->profile = 0;
		tree[i]->seq = 0;//seq[i];
		tree[i]->len = 0;//len[i]; 
		/*
		Needs to be +2 because:
		at n = 3 is is possible to get a perfectly balanced binary tree with 4 sequences at intermediate nodes
		*/
		tree[i]->links = malloc(sizeof(struct aln_tree_node*)*3);
		
		for ( j =0;j < 3;j++){
			tree[i]->links[j] = 0;
		}
	}
	
	while (cnode != numprofiles){
		max = -INFTY;
		for (i = 0;i < numseq-1; i++){
			if (as[i]){
			for ( j = i + 1;j < numseq;j++){
				if (as[j]){
				if (dm[i][j] > max){
					max = dm[i][j];
					node_a = i;
					node_b = j;
				}
				}
			}
			}
		}
		tmp = malloc(sizeof(struct aln_tree_node));
		tmp->done = 0;
		tmp->path = 0;
		tmp->profile = 0;
		tmp->num = cnode;
		tmp->seq = 0;
		tmp->len = 0;
		tmp->links = malloc(sizeof(struct aln_tree_node*)*(3));
		tmp->links[0] = tree[node_a];
		tmp->links[1] = tree[node_b];
		tmp->links[2] =0;

		tree[node_a] = tmp;
		tree[node_b] = 0;
				
		/*deactivate  sequences to be joined*/
		as[node_a] = cnode+1;
		as[node_b] = 0;
		cnode++;    
		
		/*calculate new distances*/
		for (j = numseq;j--;){
			if (j != node_b){
				dm[node_a][j] = (dm[node_a][j] + dm[node_b][j])/2;
			}
		}
		dm[node_a][node_a] = 0.0f;
		for (j = numseq;j--;){
			dm[j][node_a] = dm[node_a][j];
			dm[j][node_b] = 0.0f;
			dm[node_b][j] = 0.0f;
		}		
	}
	tmp = tree[node_a];
	
	for (i = numprofiles;i--;){
		free(dm[i]);
	}
	free(dm);
	
	
	free(tree);
	free(as);
	return tmp;
}

int* nj(float **dm,int* tree)
{
	int i,j;
	//float **dm = 0;
	float *r = 0;
	float *r_div = 0;
	int *active = 0;
	int node = 0;
	float min = 0;
	int join_a = 0;
	int join_b = 0;
	int leaves = 0;
	int c =0;

	leaves = numseq;
	
	r = malloc ((numseq*2-1) *sizeof(float));
	r_div = malloc ((numseq*2-1) *sizeof(float));
	active = malloc((numseq*2-1)*sizeof(int));
	for ( i = 0;i < numseq*2-1;i++){
		active[i] = 0;
	}
	for ( i = 0;i < numseq;i++){
		active[i] = 1;
	}
	
	node = numseq;
	while (node != numseq*2 -1){
		for (i = 0;i<numseq*2-1;i++){
			if (active[i]){
				r[i] = 0;
				for (j = 0;j < numseq*2-1;j++){
					if (active[j]){
						r[i] += (i<j) ?dm[i][j]:dm[j][i];
					}
				}
				r_div[i] = r[i] / (leaves-2);
			}
		}
		for ( j = 0;j < numseq*2-1;j++){
			if (active[j]){
			for ( i = j+1;i < numseq*2-1;i++){
				if (active[i]){
				dm[i][j] = dm[j][i] - (r[i] + r[j])/2;
				}
			}
			}
		}
		min = -INFTY;
		for ( j = 0;j < numseq*2-1;j++){
			if (active[j]){
			for ( i = j+1;i < numseq*2-1;i++){
				if (active[i]){
					if (dm[i][j] > min){
						min = dm[i][j];
						join_a = j;
						join_b = i;
					}
				}
			}
			}
		}
		//join_a always smaller than join_b && both smaller than node
		dm[join_a][node] =  dm[join_a][join_b]/2 + (r_div[join_a] - r_div[join_b])/2;
		dm[join_b][node] =  dm[join_a][join_b] - dm[join_a][node];

		active[join_a] = 0;
		active[join_b] = 0;
		tree[c] = join_a;
		tree[c+1] = join_b;
		tree[c+2] = node;

		for (i = 0;i<numseq*2-1;i++){
			if (active[i]){
				dm[i][node] = (i>join_a) ? dm[join_a][i]: dm[i][join_a];
				dm[i][node] -= dm[join_a][node];
				dm[i][node] += (i > join_b) ? dm[join_b][i] : dm[i][join_b] ;
				dm[i][node] -= dm[join_b][node];
				dm[i][node] /= 2;
			}
		}
		active[node] = 1;
		c += 3;
		node++;

	}


	for (i = numprofiles;i--;){
		free(dm[i]);
	}
	free(dm);

	free(r);
	free(r_div);
	free(active);

	return tree;
}

int* upgma(float **dm,int* tree)
{
	int i,j,t;
	int *as = 0;
	float max;
	int node_a = 0;
	int node_b = 0;
	int cnode = numseq;

	as = malloc(sizeof(int)*numseq);
	for (i = numseq; i--;){
		as[i] = i+1;
	}

	
	t = 0;
	while (cnode != numprofiles){
		max = -INFTY;
		for (i = 0;i < numseq-1; i++){
			if (as[i]){
			for ( j = i + 1;j < numseq;j++){
				if (as[j]){
				if (dm[i][j] > max){
					max = dm[i][j];
					node_a = i;
					node_b = j;
				}
				}
			}
			}
		}
		
		tree[t] = as[node_a]-1;
		tree[t+1] = as[node_b]-1;
		tree[t+2] = cnode;
		t += 3;	
		
		/*deactivate  sequences to be joined*/
		as[node_a] = cnode+1;
		as[node_b] = 0;
		cnode++;    
		
		/*calculate new distances*/
		for (j = numseq;j--;){
			if (j != node_b){
				dm[node_a][j] = (dm[node_a][j] + dm[node_b][j])/2;
			}
		}
		dm[node_a][node_a] = 0.0f;
		for (j = numseq;j--;){
			dm[j][node_a] = dm[node_a][j];
			dm[j][node_b] = 0.0f;
			dm[node_b][j] = 0.0f;
		}		
	}
	free(as);
	return tree;
}
