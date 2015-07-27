/*
	kalign2_main.c 
	
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



unsigned int numseq = 0;
unsigned int numprofiles = 0;
float gpo = 0;
float gpe = 0;
float tgpe = 0;

int main(int argc,char **argv)
{
	int i;
	int* tree = 0;
	int a, b, c;
	
	struct alignment* aln = 0;
	struct parameters* param = 0;
	struct aln_tree_node* tree2 = 0;
	
	param = malloc(sizeof(struct parameters));
	
	param =  interface(param,argc,argv);
	
	aln = detect_and_read_sequences(aln,param);
	
	if(param->ntree > numseq){
		param->ntree = numseq;
	}

	//DETECT DNA
	if(param->dna == -1){
		for (i = 0; i < numseq;i++){
			param->dna = byg_detect(aln->s[i],aln->sl[i]);
			if(param->dna){
				break;
			}
		}
	}
	//param->dna = 0;
	//fprintf(stderr,"DNA:%d\n",param->dna);
	//exit(0);
	
	if(param->dna == 1){
		//brief sanity check...
		for (i = 0; i < numseq;i++){
			if(aln->sl[i] < 6){
				fprintf(stderr,"Dna/Rna alignments are only supported for sequences longer than 6.");
				free(param);
				free_aln(aln);
				exit(0);
			}
		}
		aln =  make_dna(aln);
	}

	int j;
	
	if(param->reformat){
		for (i = 0 ;i < numseq;i++){
			aln->nsip[i] = i;
			for (j = 0; j < aln->sl[i];j++){
				aln->s[i][j] = 0;
			}
		}
		param->format = "fasta";//param->reformat;
		output(aln,param);
		exit(1);
	}
	
	
	
	//fast distance calculation;
	float** submatrix = 0;
	submatrix = read_matrix(submatrix,param); // sets gap penalties as well.....
	
	if(!param->quiet){
		parameter_message(param);
	}
	
	if(byg_start(param->alignment_type,"profPROFprofilePROFILE") != -1){
		profile_alignment_main(aln,param,submatrix);
	}
	
	float** dm = 0;
	if(param->ntree > 1){
		if(byg_start(param->distance,"pairclustalPAIRCLUSTAL") != -1){
			if(byg_start(param->tree,"njNJ") != -1){
				dm = protein_pairwise_alignment_distance(aln,dm,param,submatrix,1);
			}else{
				dm = protein_pairwise_alignment_distance(aln,dm,param,submatrix,0);
			}
		}else if(byg_start("wu",param->alignment_type) != -1){
			dm =  protein_wu_distance2(aln,dm,param);
		//	param->feature_type = "wumanber";
		}else if(param->dna == 1){
			if(byg_start(param->tree,"njNJ") != -1){
				dm =  dna_distance(aln,dm,param,1);
			}else{
				dm =  dna_distance(aln,dm,param,0);
			}
		}else{
			if(byg_start(param->tree,"njNJ") != -1){
				dm =  protein_wu_distance(aln,dm,param,1);
			}else{
				dm =  protein_wu_distance(aln,dm,param,0);
			}
		}
		/*int j; 
		for (i = 0; i< numseq;i++){
			for (j = 0; j< numseq;j++){
				fprintf(stderr,"%f	",dm[i][j]);
			}
			fprintf(stderr,"\n");
		}*/

		if(byg_start(param->tree,"njNJ") != -1){
			tree2 = real_nj(dm,param->ntree);
		}else{
			tree2 = real_upgma(dm,param->ntree);
		}
		if(param->print_tree){
			print_tree(tree2,aln,param->print_tree);
		}
	}

	tree = malloc(sizeof(int)*(numseq*3+1));
	for ( i = 1; i < (numseq*3)+1;i++){
		tree[i] = 0;
	}
	tree[0] = 1; 
	
	if(param->ntree < 2){
		tree[0] = 0;
		tree[1] = 1;
		
		c = numseq;
		tree[2] = c;
		a = 2;
		for ( i = 3; i < (numseq-1)*3;i+=3){
			tree[i] = c;
			tree[i+1] = a;
			c++;
			tree[i+2] = c;
			a++;
		}
	}else if(param->ntree > 2){
		ntreeify(tree2,param->ntree);
	}else{
		tree = readtree(tree2,tree);
		for (i = 0; i < (numseq*3);i++){
			tree[i] = tree[i+1];
		}
		free(tree2->links);
		free(tree2->internal_lables);
		free(tree2);
	}

	


	//get matrices... 
	struct feature_matrix* fm = 0;
	
	struct ntree_data* ntree_data = 0;
	
	int** map = 0;
	if(param->ntree > 2){
		ntree_data = malloc(sizeof(struct ntree_data));
		ntree_data->realtree = tree2;
		ntree_data->aln = aln;
		ntree_data->profile = 0;
		ntree_data->map = 0;
		ntree_data->ntree = param->ntree;
		ntree_data->submatrix = submatrix;
		ntree_data->tree = tree; 
		
		ntree_data = ntree_alignment(ntree_data);
		map = ntree_data->map;
		tree = ntree_data->tree;
		for (i = 0; i < (numseq*3);i++){
			tree[i] = tree[i+1];
		}
		free(ntree_data);
	}else if (param->feature_type){
		fm = get_feature_matrix(fm,aln,param);
		if(!fm){
			for (i = 32;i--;){
				free(submatrix[i]);
			}
			free(submatrix);
			free_param(param);
			free(map);
			free(tree);
			exit(0);
		}
		
		map = feature_hirschberg_alignment(aln,tree,submatrix,map,fm);
		//exit(0);
		//map =  feature_alignment(aln,tree,submatrix, map,fm);

	}else if (byg_start("pairwise",param->alignment_type) != -1){
		if(param->dna == 1){
			map = dna_alignment_against_a(aln,tree,submatrix, map,param->gap_inc);
		}else{
			map = hirschberg_alignment_against_a(aln,tree,submatrix, map,param->smooth_window,param->gap_inc);
		}
		//map =  default_alignment(aln,tree,submatrix, map);
	}else if (byg_start("fast",param->alignment_type) != -1){
		map =  default_alignment(aln,tree,submatrix, map);
	}else if(param->dna == 1){
		map =  dna_alignment(aln,tree,submatrix, map,param->gap_inc);
	/*}else if (byg_start("test",param->alignment_type) != -1){
		map =  test_alignment(aln,tree,submatrix, map,param->internal_gap_weight,param->smooth_window,param->gap_inc);
	}else if (param->aa){
		map =  aa_alignment(aln,tree,submatrix, map,param->aa);
	}else if (param->alter_gaps){
		map = alter_gaps_alignment(aln,tree,submatrix,map,param->alter_gaps,param->alter_range,param->alter_weight);
	}else if (byg_start("altergaps",param->alignment_type) != -1){
		map = alter_gaps_alignment(aln,tree,submatrix,map,param->alter_gaps,param->alter_range,param->alter_weight);
	}else if(byg_start("simple",param->alignment_type) != -1){
		map =  simple_hirschberg_alignment(aln,tree,submatrix, map);*/
	}else if(byg_start("advanced",param->alignment_type) != -1){
		map =  advanced_hirschberg_alignment(aln,tree,submatrix, map,param->smooth_window,param->gap_inc,param->internal_gap_weight);
	}else{
		map =  hirschberg_alignment(aln,tree,submatrix, map,param->smooth_window,param->gap_inc);
	}
	
	
	//clear up sequence array to be reused as gap array....
	int *p = 0;
	for (i = 0; i < numseq;i++){
		p = aln->s[i];
		for (a = 0; a < aln->sl[i];a++){
			p[a] = 0;
		}
	}
	//clear up
	
	for (i = 0; i < (numseq-1)*3;i +=3){
		a = tree[i];
		b = tree[i+1];
		aln = make_seq(aln,a,b,map[tree[i+2]]);
	}
	
	//for (i = 0; i < numseq;i++){
	//	fprintf(stderr,"%s	%d\n",aln->sn[i],aln->nsip[i]);
	//}

	
	for (i = 0; i < numseq;i++){
		aln->nsip[i] = 0;
	}
	
	
	
	aln =  sort_sequences(aln,tree,param->sort);

	//for (i = 0; i < numseq;i++){
	//	fprintf(stderr,"%d	%d	%d\n",i,aln->nsip[i],aln->sip[i][0]);
	//}
	
	
	output(aln,param);
/*	if(!param->format){
		fasta_output(aln,param->outfile);
	}else{
		if (byg_start("msf",param->format) != -1){
			msf_output(aln,param->outfile);
		}else if (byg_start("clustal",param->format) != -1){
			clustal_output(aln,param->outfile);
		}else if (byg_start("macsim",param->format) != -1){
			macsim_output(aln,param->outfile,param->infile[0]);
		}
	}
	free_param(param);*/
	
	free(map);
	free(tree);
	return 0;
}



