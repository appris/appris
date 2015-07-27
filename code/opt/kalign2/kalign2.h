/*
	kalign2.h
	
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

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#define INFTY FLT_MAX
#define FLOATINFTY FLT_MAX

#define NODESIZE 16


#ifdef MEMORY 
#define tmalloc malloc
#endif


extern unsigned int numseq;
extern unsigned int numprofiles;
extern float gpo;
extern float gpe;
extern float tgpe;

struct feature_matrix{
	float** m;
	int mdim;
	int stride;
};

struct utype_ufeat{
	struct feature *t;
	struct feature *f;
};

struct parameters{
	char **infile;
	char *input;
	char *outfile;
	char* format;
	//int reformat;
	char* feature_type;
	char* alignment_type;
	char* feature_mode;
	char* distance;
	char* tree;
	char* sort;
	char* sub_matrix;
	char* print_tree;
	char* print_svg_tree;
	float gpo;
	float gpe;
	float tgpe;
	float secret;
	float zlevel;
	float same_feature_score;
	float diff_feature_score;
	
	int reformat;
	int id;
	int aa;
	int alter_gaps;
	int ntree;
	int help_flag;
	int quiet;
	
	int dna;
	float alter_range;
	int alter_weight;
	float internal_gap_weight;
	int smooth_window;
	float gap_inc;
};

struct node{
	struct node *next;
        int pos;
};

struct names{
	int* start;
	int* end;
	int* len; 
};

struct bignode{
	struct bignode *next;
	unsigned int pos[NODESIZE];
	unsigned int num;
};

struct bignode* big_insert_hash(struct bignode *n,const unsigned int pos);
void big_remove_nodes(struct bignode *n);
void big_print_nodes(struct bignode *n);



struct alignment{
	//struct node** seq;
	struct feature** ft;
	struct sequence_info** si;
	unsigned int** sip;
	unsigned int* nsip;
	unsigned int* sl;
	unsigned int* lsn;
	int** s;
	char**seq;
	char** sn;
};

struct sequence_info{
	struct sequence_info* next;
	char* name;
	char* value;
};

struct feature{
	struct feature *next;
	char* type;
	char* note;
	int start;
	int end;
	int color;
};

struct hirsch_mem{
	struct states* f;
	struct states* b;
	int starta;
	int startb;
	int enda;
	int endb;
	int size;
	int len_a;
	int len_b;
};

struct dp_matrix{
	struct states* s;
	void* tb_mem;
	char** tb;
	int x;
	int y;
};

struct states{
	float a;
	float ga;
	float gb;
	float x;
};

struct aln_tree_node{
	struct aln_tree_node** links;
	int* internal_lables;
	int* path;
	int* profile;
	int* seq;
	int len;
	int done;
	int num;
};

struct tree_node{
	struct tree_node* left;
	struct tree_node*right;
 	int label;
 	int edge;
};

struct ntree_data{
	struct aln_tree_node* realtree;
	struct alignment* aln;
	float** profile;
	int** map;
	float**submatrix;
	int* tree;
	int ntree;
};

struct alignment* sort_sequences(struct alignment* aln,int* tree,char* sort);

struct aln_tree_node* real_upgma(float **dm,int ntree);

int* readtree(struct aln_tree_node* p,int* tree);

struct parameters* interface(struct parameters* param,int argc,char **argv);
void parameter_message(struct parameters* param);

struct dp_matrix* dp_matrix_alloc(struct dp_matrix *dp,int x,int y);
struct dp_matrix* dp_matrix_realloc(struct dp_matrix *dp,int x,int y);
void dp_matrix_free(struct dp_matrix *dp);

struct alignment* detect_and_read_sequences(struct alignment* aln,struct parameters* param);
void output(struct alignment* aln,struct parameters* param);

int* upgma(float **dm,int* tree);
int* nj(float **dm,int* tree);
void print_simple_phylip_tree(struct aln_tree_node* p);


struct alignment* make_dna(struct alignment* aln);

float** read_matrix(float** subm,struct parameters* param);

int* f_only_pp_dyn(int* path, struct dp_matrix *dp,const float* fprof1,const float* fprof2,const int len_a,const int len_b,int fdim,int stride);

int* fpp_dyn(int* path, struct dp_matrix *dp,const float* prof1,const float* prof2,const float* fprof1,const float* fprof2,const int len_a,const int len_b,int fdim,int stride);

int* dna_pp_dyn(int* path, struct dp_matrix *dp,const int* prof1,const int* prof2,const int len_a,const int len_b);

int* pp_dyn(int* path, struct dp_matrix *dp,const float* prof1,const float* prof2,const int len_a,const int len_b);
int* ps_dyn(int* path, struct dp_matrix *dp,const float* prof1,const int* seq2,const int len_a,const int len_b,int sip);
int* ss_dyn(float**subm,int* path, struct dp_matrix *dp,const int* seq1,const int* seq2,const int len_a,const int len_b);

int* mirror_path(int* path);

float* make_profile(float* prof,int* seq,int len, float** subm);
float* dna_make_profile(float* prof,int* seq,int len, float** subm);

float* update(const float*profa, const float* profb,float* newp,int* path,int sipa,int sipb);
float* update_only_a(const float* profa, const float* profb,float* newp,int* path,int sipa,int sipb);
float* dna_update(const float*profa,const float* profb,float* newp,int* path,int sipa,int sipb);
float* dna_update_only_a(const float* profa, const float* profb, float* newp,int* path,int sipa,int sipb);


void set_gap_penalties(float* prof,int len,int nsip,float strength,int nsip_c);
void dna_set_gap_penalties(float* prof,int len,int nsip,float strength,int nsip_c);

float** protein_pairwise_alignment_distance(struct alignment* aln,float** dm,struct parameters* param,float**subm, int nj);
float get_distance_from_pairwise_alignment(int* path,int* seq1,int* seq2);

float** protein_wu_distance2(struct alignment* si,float** dm,struct parameters* param);
float protein_wu_distance_calculation2(struct node* hash[],int* seq,int seqlen,int diagonals,int mode);

float** protein_wu_distance(struct alignment* si,float** dm,struct parameters* param, int nj);
//float protein_wu_distance_calculation(struct node* hash[],int* seq,int seqlen,int diagonals,int mode);

float protein_wu_distance_calculation(struct bignode* hash[], const int* seq, const int seqlen,const int diagonals, const float mode);

float** dna_distance(struct alignment* si,float** dm,struct parameters* param,int nj);
float dna_distance_calculation(struct bignode* hash[],int* p,int seqlen,int diagonals,float mode);


int byg_detect(int* text,int n);

int check_identity(char* n,char*m);

int byg_count(char* pattern,char*text);
int byg_start(char* pattern,char*text);
int byg_end(char* pattern,char*text);

struct node* insert(struct node *n, int pos);
struct node* insert_hash(struct node *n, int pos);

void remove_nodes(struct node *n);

#ifndef MEMORY 
void* tmalloc(int size);
#endif


struct alignment* aln_alloc(struct alignment* aln);
void free_aln(struct alignment* aln);
void free_param(struct parameters* param);
void free_ft(struct feature* n);

int* pp_dyn2(int* path, struct dp_matrix *dp,const int* prof1,const int* prof2,const int len_a,const int len_b);
int* ps_dyn2(int* path, struct dp_matrix *dp,const int* prof1,const int* seq2,const int len_a,const int len_b,int sip);
int* ss_dyn2(int**subm,int* path, struct dp_matrix *dp,const int* seq1,const int* seq2,const int len_a,const int len_b);


float* make_profile2(float* prof, int* seq,int len, float** subm);
void set_gap_penalties2(float* prof,int len,int nsip,int window,float strength);
float* update2(const float* profa,const float* profb,float* newp,int* path,int sipa,int sipb,float internal_gap_weight);

struct feature_matrix* get_feature_matrix(struct feature_matrix* fm, struct alignment* aln,struct parameters*param);
void free_utf(struct utype_ufeat* utf);
void free_feature_matrix(struct feature_matrix* fm);

struct utype_ufeat* get_unique_features(struct alignment* aln,struct utype_ufeat* utf);
struct utype_ufeat* traverse_ft(struct utype_ufeat* utf,struct feature* n);
struct feature* add_unique_feature(struct feature *n, struct feature *toadd);
struct feature* add_unique_type(struct feature *n, struct feature *toadd);

int** default_alignment(struct alignment* aln,int* tree, float**submatrix, int** map);
int** feature_alignment(struct alignment* aln,int* tree,float**submatrix, int** map,struct feature_matrix* fm);
int** test_alignment(struct alignment* aln,int* tree,float**submatrix, int** map,float internal_gap_weight,int window,float strength);



struct ntree_data* ntree_alignment(struct ntree_data* ntree_data);
struct ntree_data* ntree_sub_alignment(struct ntree_data* ntree_data,int* tree,int num);

float* make_feature_profile(float* prof,struct feature* f,int len,struct feature_matrix* fm);
float*  feature_update(const float* profa, const float* profb,float* newp,int* path,int stride);

void printtree(struct aln_tree_node* p);

struct ntree_data* alignntree(struct ntree_data* ntree_data,struct aln_tree_node* p);

//int** alignntree(struct alignment* aln,int** submatrix, struct aln_tree_node* p,int** map,int ntree);
void ntreeify(struct aln_tree_node* p,int ntree);

struct tree_node* simpleinsert(struct tree_node* p,int target, int new_edge,int leaf_label);
void printsimpleTree(struct tree_node* p);
int* ticker(int* milometer,int elements);
int* readsimpletree(struct tree_node* p,int* tree);
int add_label_simpletree(struct tree_node* p,int* nodes,int i);
//int** find_best_topology(struct alignment* aln,int**submatrix,int** map,int* leaves,int* nodes,int ntree);
void free_real_tree(struct aln_tree_node* p);
struct ntree_data* find_best_topology(struct ntree_data* ntree_data,int* leaves,int* nodes);
void freesimpletree(struct tree_node* p);

struct aln_tree_node* real_nj(float **dm,int ntree);

//int** alter_gaps_alignment(struct alignment* aln,int* tree,int**submatrix, int** map,int n,float range,int weight);
//void add_feature_information_from_alignment(int* path,int* fprof1,int* fprof2,int weight);

struct alignment* protein_wu_sw(struct node* hash[],struct alignment* aln,int a,int b);
float protein_wu_distance_calculation3(struct node* hash[],int* seq,int seqlen,int diagonals,int mode);
float* make_wu_profile(float* prof,float* wu,int len);

//int** aa_alignment(struct alignment* aln,int* tree,int**submatrix, int** map,int mmbonus);
//int* aapp_dyn(int* path, struct dp_matrix *dp,const int* prof1,const int* prof2,const int len_a,const int len_b,const int mmbonus);




int** hirschberg_alignment(struct alignment* aln,int* tree,float**submatrix, int** map,int window,float strength);
int** hirschberg_alignment_against_a(struct alignment* aln,int* tree,float**submatrix, int** map,int window,float strength);
//int* foward_pp_dyn(int* path, struct dp_matrix *dp,const float* prof1,const float* prof2,const int len_a,const int len_b);
//int* backward_pp_dyn(int* path, struct dp_matrix *dp,const float* prof1,const float* prof2,const int len_a,const int len_b);




struct hirsch_mem* hirsch_mem_alloc(struct hirsch_mem* hm,int x);
struct hirsch_mem* hirsch_mem_realloc(struct hirsch_mem* hm,int x);
void hirsch_mem_free(struct hirsch_mem* hm);

int* mirror_hirsch_path(int* hirsch_path,int len_a,int len_b);
int* add_gap_info_to_hirsch_path(int* hirsch_path,int len_a,int len_b);

//DNA alignment via hirsch/Myer Miller

int** dna_alignment(struct alignment* aln,int* tree,float**submatrix, int** map,float strength);
int** dna_alignment_against_a(struct alignment* aln,int* tree,float**submatrix, int** map,float strength);

struct alignment* make_seq(struct alignment* aln,int a,int b,int* path);
void update_gaps(int old_len,int*gis,int new_len,int *newgaps);
//void print_alignment(struct alignment* aln);

struct alignment* sort_in_relation(struct alignment* aln,char* sort);
void quickSort(struct alignment* aln, int array_size);
void q_sort(struct alignment* aln, int left, int right);

void smooth_gaps(float* prof,int len,int window,float strength);



int** advanced_hirschberg_alignment(struct alignment* aln,int* tree,float**submatrix, int** map,int window,float strength,float internal_gap_weight);

int** simple_hirschberg_alignment(struct alignment* aln,int* tree,float**submatrix, int** map);

float* simple_make_profile(float* prof, int* seq,int len, float** subm);
float* simple_update(float* profa,float* profb, float* newp,int* path);

int* simple_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm, int* hirsch_path);
int* simple_hirsch_align_two_pp_vector(const float* prof1,const float* prof2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);
struct states* simple_foward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm);
struct states* simple_backward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm);

int** feature_hirschberg_alignment(struct alignment* aln,int* tree,float**submatrix, int** map,struct feature_matrix* fm);

void profile_alignment_main(struct alignment* aln,struct parameters* param,float** submatrix);


void increase_gaps(float* prof,int len,int window,float strength);


struct names* names_alloc(struct names* n);
void names_free(struct names* n);


void print_tree(struct aln_tree_node* p,struct alignment* aln,char* outfile);
void print_newick_tree(struct aln_tree_node* p,struct alignment* aln, FILE *fout);
void print_phyloxml_tree(struct aln_tree_node* p,struct alignment* aln,FILE *fout);


struct alignment* phylo (struct alignment* aln,char* outfile);




