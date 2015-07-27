#include <time.h>
#include <stdio.h>

/* various constants */
#define END -999 /* used in to_trans_array, from_trans_array, silent_vertices, etc
		  * for marking the end of a transition group */
#define TOT_END -1010 /* used as a more final end marker than END */
#define DEFAULT 99.0 /* used in log_matrices to represent probability 0.0 */
#define NO_PREV -1 /* used in viterbi to indicate no max value is found yet */
#define NOPROB -1 /* used in forward, backward and viterbi for signaling to
		   * caller that the given sequence has 0 prob of being produced
		   * by the hmm */
#define OK 0 /* all went as expected in the function that returns this value */
#define SILENT 888.0 /* if emission probability = SILENT (in emissions and log emissions)
		      * then this state is a silent state */

#define MAX_SEQ_NAME_SIZE 100 /* maximal size of a sequence name, if larger it will
			       * simply be truncated */
#define MAX_GAP_SIZE 20 /* maximum gap size, meaning paths between states longer than this value will
			 * be disregarded */
#define ENDLETTER '\0'

#define NO -1
#define YES 0
#define OH_MY_YES 100

/* hmm file types */
#define SINGLE_HMM 0
#define MULTI_HMM 1

/* core alg scoring methods for msa */
#define DOT_PRODUCT 1
#define SUBST_MTX_DOT_PRODUCT 3
#define SUBST_MTX_PRODUCT 2
#define SUBST_MTX_DOT_PRODUCT_PRIOR 4
#define SJOLANDER 5
#define PICASSO 6
#define PICASSO_SYM 7
#define DOT_PRODUCT_PICASSO 8
#define SJOLANDER_REVERSED 9

/* multi alphabet scoring methods */
#define JOINT_PROB 10
#define AVERAGE 11

/* training algorithms */
#define BW_STD 0
#define CML_STD 1

/* types of modules */
#define SINGLENODE 2
#define CLUSTER 3
#define SINGLELOOP 4
#define FORWARD_STD 5
#define FORWARD_ALT 6
#define PROFILE7 7
#define PROFILE9 9

/* types of vertices */
#define STANDARDV 101
#define SILENTV 102
#define STARTV 100
#define LOCKEDV 104
#define PROFILEV 103
#define ENDV 999

/* sequence formats */
#define STANDARD 0 /* default */
#define FASTA 1
#define LABELED 3
#define PROFILE 4
#define MSA_STANDARD 2

/* score algorithms */
#define FORW_ALG 10
#define VITERBI_ALG 11
#define ONE_BEST_ALG 12

/* max alphabet size */
#define MAX_NR_ALPHABETS 4;

/* alphabet types */
#define DISCRETE 1
#define CONTINUOUS 2

/* model maker options */
#define FAST 1
#define QUERY 2
#define HAND 3


/* weighting schemes */
#define NONE -2
#define HENIKOFF 1

//#define DEBUG
//#define DEBUG2

/* Structure: hmm_s
 * 
 * Declaration of a general HMM
 */

struct hmm_s {
  
  /* general info */
  char name[100];  /* name of the HMM */
  struct time_t *constr_t; /* time of construction */
  char alphabet[1000]; /* the alphabet */
  int alphabet_type; /* discrete or continuous */
  int a_size; /* size of alphabet */
  int nr_m; /* nr of modules */
  int nr_v; /* nr of vertices (states) */
  int nr_t; /* nr of transitions */
  int nr_d; /* nr of distribution groups */
  int nr_dt; /* total nr of distribution ties */
  int nr_ttg; /* nr of transition tie groups */
  int nr_tt; /* total nr of transition ties */
  int nr_ed; /* nr of emission dirichlet structs */
  int nr_tp;
  int startnode; /* nr of startnode */
  int endnode; /* nr of endnode */
  struct replacement_letter_s *replacement_letters; /* pointer to wildcard letters */
  struct module_s **modules; /* pointer to array of pointers to the modules */

  /* data structures */
  int *silent_vertices; /* pointer to array of the silent vertices */
  int *locked_vertices; /* pointer to arrray of locked vertices */
  char *vertex_labels; /* pointer to array of vertex labels */
  char *labels; /* pointer to array which contains the set of vertex labels */
  int nr_labels; /* the number of labels (set count) for this hmm */
  double *vertex_trans_prior_scalers; /* pointer to values that decide how to weight prior distribution contribution
				       * when updating vertex parameters */
  double *vertex_emiss_prior_scalers; /* --------------------------------- " --------------------------------------- */
  double *transitions; /* pointer to transition probabilities matrix,
		       * from_v on vertical, to_v on horizontal */
  double *log_transitions; /* pointer to log of trans probs, for viterbi */
  double *emissions; /* pointer to emission probabilities matrix, letters on horizontal,
		     *  vertices on vertical */
  double *log_emissions; /* pointer to log of emiss probs, for viterbi */
  double *tot_transitions; /* pointer to transition probability matrix storing the total prob of going from a to b, adding all
			    * paths via silent vertices */
  double *max_log_transitions; /* pointer to log of maximal transprob for going from a to b of all possible paths via silent states*/
  struct path_element **from_trans_array; /* pointer to an array that for each vertex stores
					   * a pointer to the paths to the vertices that has
					   * a transition to that vertex (including via one or more
					   * silent states */
  struct path_element **to_trans_array; /* pointer to an array that for each vertex stores
					 * a pointer to the paths to the vertices that it has
					 * a transition to (including via one or more silent
					 * states) */
  int **to_silent_trans_array; /* pointer to array that for each vertex stores a pointer to the silent vertices this vertex
				* has a DIRECT transition to */
  struct path_element **tot_to_trans_array;  /* pointer to an array that for each vertex stores
					      * a pointer to the vertices that it has
					      * a transition to (including via one or more silent
					      * states) */
  struct path_element **tot_from_trans_array; /* pointer to an array that for each vertex stores
					       * a pointer to the vertices that has
					       * a transition to that vertex (including via one or more
					       * silent states */
  int *distrib_groups; /* a pointer to the arrays for the distribution groups,
			* i.e. groups of vertices that have identical emission
			* probabilities */
  struct transition_s *trans_tie_groups; /* a pointer to the arrays for the transition tie groups, 
					   * i.e. groups of transitions that have identical transition probabilities */
  struct emission_dirichlet_s *emission_dirichlets; /* pointer to the different dirichlet
						     * mixtures */
  struct emission_dirichlet_s **ed_ps; /* pointer to array of pointers (one for each state)
					* to the different dirichlet mixtures */
  double *subst_mtx; /* substitution matrix array giving the probability of two alphabet letters  being related */
};

struct hmm_multi_s {
  
  /* general info */
  char name[100];  /* name of the HMM */
  struct time_t *constr_t; /* time of construction */
  int nr_alphabets;
  char alphabet[1000]; /* the alphabet */
  char alphabet_2[1000]; /* the alphabet */
  char alphabet_3[1000]; /* the alphabet */
  char alphabet_4[1000]; /* the alphabet */
  int alphabet_type; /* discrete or continuous */
  int alphabet_type_2; /* discrete or continuous */
  int alphabet_type_3; /* discrete or continuous */
  int alphabet_type_4; /* discrete or continuous */
  int a_size; /* size of alphabet */
  int a_size_2; /* size of alphabet */
  int a_size_3; /* size of alphabet */
  int a_size_4; /* size of alphabet */
  int nr_m; /* nr of modules */
  int nr_v; /* nr of vertices (states) */
  int nr_t; /* nr of transitions */
  int nr_d; /* nr of distribution groups */
  int nr_dt; /* total nr of distribution ties */
  int nr_ttg; /* nr of transition tie groups */
  int nr_tt; /* total nr of transition ties */
  int nr_ed; /* nr of emission dirichlet structs */
  int nr_ed_2; /* nr of emission dirichlet structs */
  int nr_ed_3; /* nr of emission dirichlet structs */
  int nr_ed_4; /* nr of emission dirichlet structs */
  int nr_tp;
  int startnode; /* nr of startnode */
  int endnode; /* nr of endnode */
  struct replacement_letter_multi_s *replacement_letters; /* pointer to wildcard letters */
  
  struct module_multi_s **modules; /* pointer to array of pointers to the modules */

  /* data structures */
  int *silent_vertices; /* pointer to array of the silent vertices */
  int *locked_vertices; /* pointer to arrray of locked vertices */
  char *vertex_labels; /* pointer to array of vertex labels */
  char *labels; /* pointer to array which contains the set of vertex labels */
  int nr_labels; /* the number of labels (set count) for this hmm */
  double *vertex_trans_prior_scalers; /* pointer to values that decide how to weight prior distribution contribution
				       * when updating vertex parameters */
  double *vertex_emiss_prior_scalers; /* --------------------------------- " --------------------------------------- */
  double *vertex_emiss_prior_scalers_2; /* --------------------------------- " --------------------------------------- */
  double *vertex_emiss_prior_scalers_3; /* --------------------------------- " --------------------------------------- */
  double *vertex_emiss_prior_scalers_4; /* --------------------------------- " --------------------------------------- */
  double *transitions; /* pointer to transition probabilities matrix,
		       * from_v on vertical, to_v on horizontal */
  double *log_transitions; /* pointer to log of trans probs, for viterbi */
  double *emissions; /* pointer to emission probabilities matrix, letters on horizontal,
		      *  vertices on vertical */
  double *emissions_2; /* pointer to emission probabilities matrix, letters on horizontal,
			*  vertices on vertical */
  double *emissions_3; /* pointer to emission probabilities matrix, letters on horizontal,
			*  vertices on vertical */
  double *emissions_4; /* pointer to emission probabilities matrix, letters on horizontal,
			*  vertices on vertical */
  double *log_emissions; /* pointer to log of emiss probs, for viterbi */
  double *log_emissions_2; /* pointer to log of emiss probs, for viterbi */
  double *log_emissions_3; /* pointer to log of emiss probs, for viterbi */
  double *log_emissions_4; /* pointer to log of emiss probs, for viterbi */
  double *tot_transitions; /* pointer to transition probability matrix storing the total prob of going from a to b, adding all
			   * paths via silent vertices */
  double *max_log_transitions; /* pointer to log of maximal transprob for going from a to b of all possible paths via silent states*/
  struct path_element **from_trans_array; /* pointer to an array that for each vertex stores
					   * a pointer to the paths to the vertices that has
					   * a transition to that vertex (including via one or more
					   * silent states */
  struct path_element **to_trans_array; /* pointer to an array that for each vertex stores
					 * a pointer to the paths to the vertices that it has
					 * a transition to (including via one or more silent
					 * states) */
  int **to_silent_trans_array; /* pointer to array that for each vertex stores a pointer to the silent vertices this vertex
				* has a DIRECT transition to */
  
  struct path_element **tot_to_trans_array;  /* pointer to an array that for each vertex stores
					      * a pointer to the vertices that it has
					      * a transition to (including via one or more silent
					      * states) */
  struct path_element **tot_from_trans_array; /* pointer to an array that for each vertex stores
					       * a pointer to the vertices that has
					       * a transition to that vertex (including via one or more
					       * silent states */
  int *distrib_groups; /* a pointer to the arrays for the distribution groups,
			* i.e. groups of vertices that have identical emission
			* probabilities */
  struct transition_s *trans_tie_groups; /* a pointer to the arrays for the transition tie groups, 
					   * i.e. groups of transitions that have identical transition probabilities */
  struct emission_dirichlet_s *emission_dirichlets; /* pointer to the different dirichlet
						     * mixtures */
  struct emission_dirichlet_s *emission_dirichlets_2; /* pointer to the different dirichlet
						       * mixtures */
  struct emission_dirichlet_s *emission_dirichlets_3; /* pointer to the different dirichlet
						       * mixtures */
  struct emission_dirichlet_s *emission_dirichlets_4; /* pointer to the different dirichlet
						       * mixtures */
  struct emission_dirichlet_s **ed_ps; /* pointer to array of pointers (one for each state)
					* to the different dirichlet mixtures */
  struct emission_dirichlet_s **ed_ps_2; /* pointer to array of pointers (one for each state)
					  * to the different dirichlet mixtures */
  struct emission_dirichlet_s **ed_ps_3; /* pointer to array of pointers (one for each state)
					  * to the different dirichlet mixtures */
  struct emission_dirichlet_s **ed_ps_4; /* pointer to array of pointers (one for each state)
					  * to the different dirichlet mixtures */
  double *subst_mtx; /* substitution matrix array giving the probability of two alphabet letters  being related */
  double *subst_mtx_2; /* substitution matrix array giving the probability of two alphabet letters  being related */
  double *subst_mtx_3; /* substitution matrix array giving the probability of two alphabet letters  being related */
  double *subst_mtx_4; /* substitution matrix array giving the probability of two alphabet letters  being related */
};


/* Structure: null_model_s
 *
 * Declaration of null model struct
 */
struct null_model_s {
  double trans_prob;
  int a_size;
  char alphabet[1000];
  double *emissions;
};

/* Structure: null_model_multi_s
 *
 * Declaration of null model multi struct
 */
struct null_model_multi_s {
  int nr_alphabets;
  double trans_prob;
  int a_size;
  int a_size_2;
  int a_size_3;
  int a_size_4;
  char alphabet[1000];
  char alphabet_2[1000];
  char alphabet_3[1000];
  char alphabet_4[1000];
  double *emissions;
  double *emissions_2;
  double *emissions_3;
  double *emissions_4;
  
};


/* Structure: module_s
 *
 * Declaration of module
 */
struct module_s {
  char name[50];
  int type;
  int v_type;
  int nr_v;
  int *vertices;
  char priorfile_name[200];
};


/* Structure: module_multi_s
 *
 * Declaration of module
 */
struct module_multi_s {
  char name[50];
  int type;
  int v_type;
  int nr_v;
  int *vertices;
  char priorfile_name[200];
  char priorfile_name_2[200];
  char priorfile_name_3[200];
  char priorfile_name_4[200];
};

/* Structure emission_dirichlet_s
 *
 * Declaration of dirichlet mixture 
 */
struct emission_dirichlet_s {
  char name[200];
  int nr_components;
  int alphabet_size;
  double *q_values; /* each component's "probability" */
  double *alpha_sums; /* sums of the prior values */
  double *logbeta_values; /* precalculated beta values B(alpha) for each alpha */
  double *prior_values; /* matrix with prior values */
};


/* Structure: viterbi_s
 *
 * Declaration of viterbi matrix element
 */

struct viterbi_s {
  double prob; /* is really a log prob */
  int prev;
  struct path_element *prevp;
};

/* Structure: forward_s
 *
 * Declaration of forward matrix element
 */

struct forward_s {
  double prob;
  //int distance_to_next;
};

/* Structure: backward_s
 *
 * Declaration of backward matrix element
 */

struct backward_s {
  double prob;
  //int distance_to_next;
};

/* Structure: one_best_s
 *
 * Declaration of one_best matrix element
 */

struct one_best_s {
  double prob;
  int is_updated;
  char *labeling;
};



/* Structure: letter_prob_s
 *
 * Declaration of letter probability distribution element
 */
struct letter_prob_s {
  char letter;
  double share;
};


/* Structure: path_element
 *
 * Declaration of trans_array elements
 */
struct path_element {
  int vertex;
  struct path_element *next;
};

/* Structure: letter_s
 *
 * Declaration of alphabet symbol
 */
struct letter_s {
  char letter[5];
  char label;
  double cont_letter;
};

/* Structure: sequence_s
 *
 * Declaration of sequence info holder
 */
struct sequence_s {
  char name[MAX_SEQ_NAME_SIZE];
  int length;
  double weight;
  struct letter_s *seq;
};

/* Structure: sequence_multi_s
 *
 * Declaration of sequence info holder for multiple alphabet sequences
 */
struct sequence_multi_s {
  char name[MAX_SEQ_NAME_SIZE];
  int length;
  double weight;
  struct letter_s *seq_1;
  struct letter_s *seq_2;
  struct letter_s *seq_3;
  struct letter_s *seq_4;
};


/* Structure: sequences_s
 *
 * Declaration of struct for info about the sequences
 */

struct sequences_s {
  int nr_seqs;
  int longest_seq;
  int shortest_seq;
  int avg_seq_len;
  struct sequence_s *seqs;
};


/* Structure: sequences_multi_s
 *
 * Declaration of struct for info about the sequences
 */
struct sequences_multi_s {
  int nr_alphabets;
  int nr_seqs;
  int longest_seq;
  int shortest_seq;
  int avg_seq_len;
  struct sequence_multi_s *seqs;
};


/* Structure: MSA_letter_s
 *
 * Declaration of struct for MSA_letter
 */
struct msa_letter_s {
  double nr_occurences; /* non integer occurences exist */
  double share;
  double prior_share;
  char label;
  char query_letter[5];
  double cont_letter;
};

/* Structure: msa_sequences_s
 *
 * Declaration of struct for msa sequence info
 */
struct msa_sequences_s {
  int nr_seqs;
  int msa_seq_length;
  int nr_lead_columns;
  struct msa_letter_s *msa_seq;
  int **gaps;
  int *lead_columns_start;
  int *lead_columns_end;
  double *gap_shares;
};

/* Structure: msa_sequences_multi_s
 *
 * Declaration of struct for msa sequence info
 */
struct msa_sequences_multi_s {
  int nr_alphabets;
  int nr_seqs;
  int msa_seq_length;
  int nr_lead_columns;
  struct msa_letter_s *msa_seq_1;
  struct msa_letter_s *msa_seq_2;
  struct msa_letter_s *msa_seq_3;
  struct msa_letter_s *msa_seq_4;
  int **gaps;
  int *lead_columns_start;
  int *lead_columns_end;
  double *gap_shares;
};


/* Structure: replacement_letter_s
 *
 * Declaration of struct for replacement letter info
 */
struct replacement_letter_s {
  int nr_rl;
  struct letter_s *letters;
  double *probs;
};

/* Structure: replacement_letter_multi_s
 *
 * Declaration of struct for replacement letter info
 */
struct replacement_letter_multi_s {
   int nr_alphabets;
   int nr_rl_1;
   int nr_rl_2;
   int nr_rl_3;
   int nr_rl_4;
   struct letter_s *letters_1;
   double *probs_1;
   struct letter_s *letters_2;
   double *probs_2;
   struct letter_s *letters_3;
   double *probs_3;
   struct letter_s *letters_4;
   double *probs_4;
 };

/* Structure: aa_distrib_mtx_s
 *
 * Declaration of struct for amino acid distribution matrix
 */
struct aa_distrib_mtx_s {
  int a_size;
  double *inside_values;
  double *outside_values;
  double *membrane_values;
};


/* Structure: v_list_element_s
 *
 * Declaration of struct for n_best v_list elements
 */
struct v_list_element_s {
  int vertex; /* vertex nr */
  int address; /* pointer address of this vertex's labeling */
};

/* Structure: transition_s
 *
 * Declaration of struct for transition (used for transition ties)
 */
struct transition_s {
  int from_v;
  int to_v;
};

struct align_mtx_element_s {
  int score;
  char last;
};

struct alignment_s {
  int target_pos;
  int template_pos;
  char target_letter[5];
  char template_letter[5];
};
