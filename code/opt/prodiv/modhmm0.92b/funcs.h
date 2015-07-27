


/* function declarations */

/* readhmm */
int readhmm(FILE*, struct hmm_multi_s*);

/* readhmm_multialpha */
int readhmm_multialpha(FILE*, struct hmm_multi_s*);
void transform_singlehmmfile_to_multi(FILE *hmmfile, FILE *outfile);
int readhmm_check(FILE *hmmfile);
void copy_hmm_struct(struct hmm_multi_s *hmm, struct hmm_multi_s *retrain_hmm);

/* readseqs */
void get_sequences_std(FILE*, struct sequences_s*, struct hmm_s*);
void get_labeled_sequences_std(FILE*, struct sequences_s*, struct hmm_s*);
void get_sequences_fasta(FILE*, struct sequences_s*);
void get_sequences_msa_std(FILE*, FILE*, struct msa_sequences_s*, struct hmm_s*, int, struct replacement_letter_s*);
void get_sequences_msa_prf(FILE *seqfile, FILE *priorfile, struct msa_sequences_s *msa_seq_infop,
			   struct hmm_s *hmmp, int lead_seq);


/* readseqs_multi */
int seqfile_has_labels(FILE *seqfile);
void get_sequence_fasta_multi(FILE *seqfile, struct sequences_multi_s *seq_infop, int seq_nr);
void get_sequences_std_multi(FILE *seqfile, struct sequences_multi_s *seq_infop, struct hmm_multi_s *hmmp, int seq_nr);
void get_sequences_msa_std_multi(FILE*, FILE*, struct msa_sequences_multi_s*, struct hmm_multi_s*,
				 int, struct replacement_letter_multi_s*);
void get_sequences_msa_prf_multi(FILE *seqfile, FILE *priorfile, struct msa_sequences_multi_s *msa_seq_infop,
				 struct hmm_multi_s *hmmp);

/* savehmm */
int savehmm(FILE*, struct hmm_multi_s*);
int savehmm_multialpha(FILE*, struct hmm_multi_s*);


/* core_algorithms */
int forward(struct hmm_s*, struct letter_s*, struct forward_s**, double**, int);
int backward(struct hmm_s*, struct letter_s*, struct backward_s**, double*, int);
int viterbi(struct hmm_s*, struct letter_s*, struct viterbi_s**, int);
int one_best(struct hmm_s*, struct letter_s*, struct one_best_s**, double**, int, char*);
int msa_forward(struct hmm_s*, struct msa_sequences_s*, int,
		int, int, struct forward_s**, double**, int, int, int, double*);
int msa_backward(struct hmm_s*, struct msa_sequences_s*, int,
		int, struct backward_s**, double*, int, int, int, double*);
int msa_viterbi(struct hmm_s*, struct msa_sequences_s*, int,
		int, int, struct viterbi_s**, int, int, int, double*);
int msa_one_best(struct hmm_s*, struct msa_sequences_s*, int,
		int, int, struct one_best_s**, double**, int, char*, int, int, double*);


/* core_algorithms_multialpha */
int forward_multi(struct hmm_multi_s*, struct letter_s*,  struct letter_s*,  struct letter_s*,  struct letter_s*,
	    struct forward_s**, double**, int, int);
int backward_multi(struct hmm_multi_s*, struct letter_s*,  struct letter_s*,  struct letter_s*,  struct letter_s*,
	     struct backward_s**, double*, int, int);
int viterbi_multi(struct hmm_multi_s*, struct letter_s*,  struct letter_s*,  struct letter_s*,  struct letter_s*,
	    struct viterbi_s**, int, int);
int one_best_multi(struct hmm_multi_s*, struct letter_s*,  struct letter_s*,  struct letter_s*,  struct letter_s*,
	     struct one_best_s**, double**, int, char*, int);
int msa_forward_multi(struct hmm_multi_s*, struct msa_sequences_multi_s*, int,
		int, int, struct forward_s**, double**, int, int, int, int, double*, double*, double*, double*);
int msa_backward_multi(struct hmm_multi_s*, struct msa_sequences_multi_s*, int,
		       int, struct backward_s**, double*, int, int, int, int, double*, double*, double*, double*);
int msa_viterbi_multi(struct hmm_multi_s*, struct msa_sequences_multi_s*, int,
		int, int, struct viterbi_s**, int, int, int, int, double*, double*, double*, double*);
int msa_one_best_multi(struct hmm_multi_s*, struct msa_sequences_multi_s*, int,
		int, int, struct one_best_s**, double**, int, char*, int, int, int, double*, double*, double*, double*);

/* tm_core_algorithms */
int tm_viterbi(struct hmm_s*, struct letter_s*, struct viterbi_s**, struct aa_distrib_mtx_s*, int);


/* training_algorithms */
void baum_welch_std(struct hmm_s*, struct sequence_s*, int, int, int);
void baum_welch_dirichlet(struct hmm_s*, struct sequence_s*, int, int, int, int, int);
void extended_baum_welch_dirichlet(struct hmm_s*, struct sequence_s*, int, int, int, int, int);
void msa_baum_welch_dirichlet(struct hmm_s*, struct msa_sequences_s*, int, int, int, int, int, int, int, int, int, int, double*);
void extended_msa_baum_welch_dirichlet(struct hmm_s*, struct msa_sequences_s*, int, int, int, int, int, int, int, int, int, int,
				       double*);

/* training_algorithms */
void baum_welch_std_multi(struct hmm_multi_s *hmmp, struct sequence_multi_s *seqsp, int nr_seqs, int annealing, int use_labels,
			  int multi_scoring_method, int use_prior);
void baum_welch_dirichlet_multi(struct hmm_multi_s *hmmp, struct sequence_multi_s *seqsp, int nr_seqs, int annealing, int use_labels,
				int use_transition_pseudo_counts, int use_emission_pseudo_counts, int multi_scoring_method,
				int use_prior);
void msa_baum_welch_dirichlet_multi(struct hmm_multi_s *hmmp, struct msa_sequences_multi_s *msa_seq_infop, int nr_seqs,
				    int annealing,
				    int use_gap_shares, int use_lead_columns, int use_labels, int use_transition_pseudo_counts,
				    int use_emission_pseudo_counts, int normalize, int scoring_method, int use_nr_occ,
				    int multi_scoring_method, double *aa_freqs, double *aa_freqs_2, double *aa_freqs_3,
				    double *aa_freqs_4, int use_prior);
void extended_msa_baum_welch_dirichlet_multi(struct hmm_multi_s *hmmp, struct msa_sequences_multi_s *msa_seq_infop,
					     int nr_seqs, int annealing,
					     int use_gap_shares, int use_lead_columns, int use_labels,
					     int use_transition_pseudo_counts,
					     int use_emission_pseudo_counts, int normalize, int scoring_method, int use_nr_occ,
					     int multi_scoring_method, double *aa_freqs, double *aa_freqs_2, double *aa_freqs_3,
					     double *aa_freqs_4, int use_prior);



/* std_funcs */
void* malloc_or_die(int);
void init_float_mtx(double*, double, int);
void init_viterbi_s_mtx(struct viterbi_s*, double, int);
void printhelp_modhmms();
void printhelp_modhmms_msa();
void printhelp_hmmtrain();
void printhelp_hmmtrain_msa();
void printhelp_modhmms_multialpha();
void printhelp_modhmms_msa_multialpha();
void printhelp_hmmtrain_multialpha();
void printhelp_hmmtrain_msa_multialpha();
void printhelp_modhmms_tm_multialpha();
void printhelp_modhmms_tm_msa_multialpha();
void printhelp_hmmtrain_tm_multialpha();
void printhelp_hmmtrain_tm_msa_multialpha();
void printhelp_modhmms_tm();
void printhelp_modhmms_tm_msa();
void printhelp_hmmtrain_tm();
void printhelp_hmmtrain_tm_msa();

void printhelp_chmmtrain();
void printhelp_chmmtrain_msa();
void printhelp_chmmtrain_multialpha();
void printhelp_chmmtrain_msa_multialpha();
void printhelp_add_alphabet();
void printhelp_add2profilehmm();
void printhelp_cal();
void printhelp_opt();

int get_mtx_index(int,int,int);
int get_alphabet_index(struct letter_s*, char*, int);
int get_alphabet_index_msa_query(char*, char*, int);
int get_replacement_letter_index(struct letter_s*, struct replacement_letter_s*);
int get_replacement_letter_index_multi(struct letter_s *c, struct replacement_letter_multi_s *replacement_letters, int alphabet);
int get_alphabet_index_single(char*, char, int);
int get_replacement_letter_index_single(char*, struct replacement_letter_s*);
int get_seq_length(struct letter_s*);
int path_length(int, int, struct hmm_s*, int);
int path_length_multi(int, int, struct hmm_multi_s*, int);
void print_seq(struct letter_s*, FILE*, int, char*, int);
struct path_element* get_end_path_start(int l, struct hmm_s *hmmp);
struct path_element* get_end_path_start_multi(int l, struct hmm_multi_s *hmmp);
char* get_profile_vertex_type(int, int*);
void get_replacement_letters(FILE*, struct replacement_letter_s*);
void get_aa_distrib_mtx(FILE *distribmtxfile, struct aa_distrib_mtx_s *aa_distrib_mtxp);
void get_replacement_letters_multi(FILE *replfile, struct replacement_letter_multi_s *replacement_lettersp);
char* letter_as_string(struct letter_s*);
char* sequence_as_string(struct letter_s*);
void get_viterbi_label_path(struct viterbi_s *cur, struct hmm_s *hmmp,
			    struct viterbi_s *viterbi_mtxp, int row, int row_size, char *labels, int *ip);
void get_viterbi_label_path_multi(struct viterbi_s *cur, struct hmm_multi_s *hmmp,
				  struct viterbi_s *viterbi_mtxp, int row, int row_size, char *labels, int *ip);
void get_viterbi_path(struct viterbi_s *cur, struct hmm_s *hmmp,
		      struct viterbi_s *viterbi_mtxp, int row, int row_size, int *path, int *ip);
void get_viterbi_path_multi(struct viterbi_s *cur, struct hmm_multi_s *hmmp,
			    struct viterbi_s *viterbi_mtxp, int row, int row_size, int *path, int *ip);
void itoa(char* s, int nr);
void ftoa(char* s, double nr, int prec);
int read_subst_matrix(double **mtx, FILE *substmtxfile);
int read_subst_matrix_multi(double **mtxpp, double **mtxpp_2, double **mtxpp_3, double **mtxpp_4, FILE *substmtxfile);
int read_prior_file(struct emission_dirichlet_s *em_di, struct hmm_s *hmmp, FILE *priorfile);
int read_frequencies(FILE *freqfile, double **aa_freqs);
int read_frequencies_multi(FILE *freqfile, double **aa_freqsp, double **aa_freqsp_2, double **aa_freqsp_3, double **aa_freqsp_4);
int read_prior_file_multi(struct emission_dirichlet_s *em_di, struct hmm_multi_s *hmmp, FILE *priorfile, int alphabet);
int locked_state(struct hmm_s *hmmp, int v);
int locked_state_multi(struct hmm_multi_s *hmmp, int v);
int get_best_reliability_score(double reliability_score_1, double reliability_score_2, double reliability_score_3);
void hmm_garbage_collection(FILE *hmmfile, struct hmm_s *hmmp);
void hmm_garbage_collection_multi(FILE *hmmfile, struct hmm_multi_s *hmmp);
void hmm_garbage_collection_multi_no_dirichlet(FILE *hmmfile, struct hmm_multi_s *hmmp);
void msa_seq_garbage_collection_multi(struct msa_sequences_multi_s *msa_seq_info, int nr_alphabets);
void seq_garbage_collection_multi(struct sequences_multi_s *seq_info, int nr_alphabets);
void get_msa_labels(FILE *labelfile, struct msa_sequences_s *msa_seq_infop, struct hmm_s *hmmp);
void get_msa_labels_all_columns(FILE *labelfile, struct msa_sequences_s *msa_seq_infop, struct hmm_s *hmmp);
int update_shares_prior(struct emission_dirichlet_s *em_di, struct hmm_s *hmmp,
			struct msa_sequences_s *msa_seq_infop, int l);
int replacement_letter(struct letter_s *cur_letterp, struct replacement_letter_s *replacement_letters, 
		       struct msa_sequences_s *msa_seq_infop, struct hmm_s *hmmp, int seq_pos);
void get_labels_multi(FILE *labelfile, struct sequences_multi_s *seq_infop, struct hmm_multi_s *hmmp, int seq_nr);
void get_msa_labels_multi(FILE *labelfile, struct msa_sequences_multi_s *msa_seq_infop, struct hmm_multi_s *hmmp);
void get_msa_labels_all_columns_multi(FILE *labelfile, struct msa_sequences_multi_s *msa_seq_infop, struct hmm_multi_s *hmmp);
int update_shares_prior_multi(struct emission_dirichlet_s *em_di, struct hmm_multi_s *hmmp,
			struct msa_sequences_multi_s *msa_seq_infop, int l, int alphabet);
int replacement_letter_multi(struct letter_s *cur_letterp, struct replacement_letter_multi_s *replacement_letters, 
		       struct msa_sequences_multi_s *msa_seq_infop, struct hmm_multi_s *hmmp, int seq_pos, int alphabet);
int get_nr_alphabets(FILE *hmmfile);
void get_set_of_labels(struct hmm_s *hmmp);
void get_set_of_labels_multi(struct hmm_multi_s *hmmp);
void get_reverse_msa_seq_multi(struct msa_sequences_multi_s *msa_seq_infop, struct msa_sequences_multi_s *reverse_msa_seq_infop,
			       struct hmm_multi_s *hmmp);
void get_reverse_seq_multi(struct sequence_multi_s *seqs, struct letter_s **reverse_seq_1,
			   struct letter_s **reverse_seq_2, struct letter_s **reverse_seq_3,
			   struct letter_s **reverse_seq_4, struct hmm_multi_s *hmmp, int seq_len);

/* std calculation funcs */
double get_single_gaussian_statescore(double mu, double sigma_square, double letter);
double get_dp_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
			 int p, double *emissions,  int vertex, int normalize, double *gap_shares);
double get_dp_picasso_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
				 int p, double *emissions,  int vertex, int normalize, double *gap_shares, double *aa_freqs);
double get_sjolander_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
			       int p, double *emissions, int vertex, int normalize, double *gap_shares);
double get_sjolander_reversed_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
					 int p, double *emissions,  int vertex, int normalize, double *gap_shares);
double get_picasso_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
			 int p, double *emissions,  int vertex, int normalize, double *gap_shares, double *aa_freqs);
double get_picasso_sym_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
			 int p, double *emissions,  int vertex, int normalize, double *gap_shares, double *aa_freqs);
double get_subst_mtx_product_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
					int p, double *emissions, int vertex, double *subst_mtx);
double get_subst_mtx_dot_product_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
					    int p, double *emissions,  int vertex, int normalize, double *gap_shares,
					    int query_index, double *subst_mtx);
double get_subst_mtx_dot_product_prior_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
						  int p, double *emissions,  int vertex, int normalize, double *gap_shares,
						  int query_index, double *subst_mtx);

void add_to_E_dot_product(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
			  int k, int a_size, int normalize);
void add_to_E_dot_product_picasso(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
				  int k, int a_size, int normalize);
void add_to_E_picasso(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
		      int k, int a_size, int normalize);
void add_to_E_picasso_sym(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
		      int k, int a_size, int normalize);
void add_to_E_sjolander_score(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
			      int k, int a_size, int normalize);
void add_to_E_sjolander_reversed_score(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
				       int k, int a_size, int normalize);
void add_to_E_subst_mtx_product(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
				int k, int a_size, int normalize, double *subst_mtx);
void add_to_E_subst_mtx_dot_product(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
				    int k, int a_size, int normalize, double *subst_mtx, char *alphabet);
void add_to_E_subst_mtx_dot_product_prior(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
					  int k, int a_size, int normalize, double *subst_mtx, char *alphabet);

void add_to_E_dot_product_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
				int k, int a_size, int normalize);
void add_to_E_dot_product_picasso_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
					 int k, int a_size, int normalize);
void add_to_E_picasso_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
			     int k, int a_size, int normalize);
void add_to_E_picasso_sym_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
				 int k, int a_size, int normalize);
void add_to_E_sjolander_score_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
				int k, int a_size, int normalize);
void add_to_E_sjolander_reversed_score_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
					      int k, int a_size, int normalize);
void add_to_E_subst_mtx_product_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
				       int k, int a_size, int normalize, double *subst_mtx);
void add_to_E_subst_mtx_dot_product_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
					   int k, int a_size, int normalize, double *subst_mtx, char *alphabet);
void add_to_E_subst_mtx_dot_product_prior_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
					  int k, int a_size, int normalize, double *subst_mtx, char *alphabet);

void update_labelings(struct one_best_s *cur_rowp, char *vertex_labels, 
		      int *sorted_v_list, int seq_len, int c, char *labels, int nr_of_labels, int nr_v);
void deallocate_row_labelings(struct one_best_s *prev_rowp, int nr_v);


/* debug_funcs */
void dump_align_matrix(int nr_rows, int nr_cols, struct align_mtx_element_s *mtx);
void dump_trans_matrix(int,int,double*);
void dump_int_trans_matrix(int nr_rows, int nr_cols, double *mtx);
void dump_emiss_matrix(int,int,double*);
void dump_post_prob_matrix(int nr_rows, int nr_cols, double *mtx);
void dump_forward_matrix(int,int,struct forward_s*);
void dump_backward_matrix(int,int,struct backward_s*);
void dump_viterbi_matrix(int nr_rows, int nr_cols, struct viterbi_s *mtx);
void dump_one_best_matrix(int, int, struct one_best_s*);
void dump_scaling_array(int,double*);
void dump_from_trans_array(int,struct path_element**);
void dump_to_trans_array(int,struct path_element**);
void dump_viterbi_path(struct viterbi_s*, struct hmm_s*, struct viterbi_s*, int, int);
void dump_viterbi_label_path(struct viterbi_s*, struct hmm_s*, struct viterbi_s*, int, int);
void dump_T_matrix(int,int,double*);
void dump_E_matrix(int,int,double*);
void dump_modules(struct hmm_s*);
void dump_distrib_groups(int*, int);
void dump_trans_tie_groups(struct transition_s*, int);
void dump_prior_struct(struct emission_dirichlet_s*);
void dump_silent_vertices(struct hmm_s*);
void dump_silent_vertices_multi(struct hmm_multi_s *hmmp);
void dump_locked_vertices(struct hmm_s*);
void dump_seqs(struct sequences_s*);
void dump_seqs_multi(struct sequences_multi_s*);
void dump_msa_seqs(struct msa_sequences_s*, int);
void dump_msa_seqs_multi(struct msa_sequences_s*, struct hmm_multi_s*);
void dump_to_silent_trans_array(int, int**);
void dump_aa_distrib_mtx(struct aa_distrib_mtx_s *aa_distrib_mtxp);
void dump_v_list(int*);
void dump_labeling(char*, int);
void dump_label_tmp_list(int *list);
void check_for_corrupt_values(int nr_rows, int nr_cols, double *mtx, char *name);
void dump_subst_mtx(int a_size, double *mtx);
void dump_multi_modules(struct hmm_multi_s *hmmp);
void dump_weights(double *total_weights, int nr_seqs);
