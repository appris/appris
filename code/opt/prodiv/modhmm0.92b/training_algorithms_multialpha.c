#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
//#include <double.h>


#include "structs.h"
#include "funcs.h"

#define INNER_BW_THRESHOLD 0.1
#define OUTER_BW_THRESHOLD 0.1
#define CML_THRESHOLD 0.2
#define TRUE 1
#define FALSE -1
#define STARTRANDOM 0

#define REST_LETTER_INDEX 0.5

/* for simulated annealing */
#define INIT_TEMP 1.0 
#define INIT_COOL 0.8
#define ANNEAL_THRESHOLD 0.1
#define DONE -10
#define ACTIVE 25

/* for transition matrix pseudo count */
#define TRANSITION_PSEUDO_VALUE 0.1
#define EMISSION_PSEUDO_VALUE 1.0

//#define DEBUG_BW
//#define DEBUG_BW_TRANS
//#define DEBUG_BW2
//#define DEBUG_EXTBW
//#define DEBUG_PRIORS
//#define DEBUG_Tkl

extern int verbose;

void update_emiss_mtx_std_multi(struct hmm_multi_s*, double*, int, int);
void update_emiss_mtx_std_continuous_multi(struct hmm_multi_s*, double*, int, int);
void update_emiss_mtx_pseudocount_multi(struct hmm_multi_s*, double*, int, int);
void update_emiss_mtx_prior_multi(struct hmm_multi_s*, double*, int, struct emission_dirichlet_s*, int);
void update_trans_mtx_std_multi(struct hmm_multi_s*, double*, int);
void update_trans_mtx_pseudocount_multi(struct hmm_multi_s*, double*, int);
void update_tot_trans_mtx_multi(struct hmm_multi_s*);
void recalculate_emiss_expectations_multi(struct hmm_multi_s*, double*, int);
void recalculate_trans_expectations_multi(struct hmm_multi_s *hmmp, double *T);
double add_Eka_contribution_multi(struct hmm_multi_s*, struct letter_s*, struct forward_s*,
				  struct backward_s*, int, int, int);
double add_Eka_contribution_continuous_multi(struct hmm_multi_s*, struct letter_s*, struct forward_s*,
					     struct backward_s*, int, int, int, double*, int);
double add_Eka_contribution_msa_multi(struct hmm_multi_s*, struct msa_sequences_multi_s*, struct forward_s*,
				      struct backward_s*, int, int, int, int);
void add_Tkl_contribution_multi(struct hmm_multi_s*, struct letter_s*, struct letter_s*, struct letter_s*,
				struct letter_s*, struct forward_s*,
				struct backward_s*, double*, int, int,
				struct path_element*, int, int, int, int, double*, int use_labels, int multi_scoring_method);
void add_Tkl_contribution_msa_multi(struct hmm_multi_s*, struct msa_sequences_multi_s*, struct forward_s*,
				    struct backward_s*, double*,
				    int, int, struct path_element*, double*, int use_gap_shares, int use_lead_columns, int i,
				    int use_labels, int scoring_method, int normalize, int multi_scoring_method,
				    double *aa_freqs, double *aa_freqs_2, double *aa_freqs_3, double *aa_freqs_4);
void random_walk_multi(struct hmm_multi_s*, double*, double*, char**, int, int, int);
int silent_state_multi(int, struct hmm_multi_s*);
void anneal_E_matrix_multi(double temperature, double *E, struct hmm_multi_s *hmmp, int alphabet);
void anneal_T_matrix_multi(double temperature, double *T, struct hmm_multi_s *hmmp);
void calculate_TE_contributions_multi(double *T, double *E, double *E_2, double *E_3, double *E_4,
				      double *T_lab, double *E_lab, double *E_lab_2, double *E_lab_3, double *E_lab_4,
				      double *T_ulab, double *E_ulab, double *E_ulab_2, double *E_ulab_3, double *E_ulab_4,
				      double *emissions, double *emissions_2, double *emissions_3, double *emissions_4,
				      double *transitions, int nr_v, int a_size, int a_size_2, int a_size_3, int a_size_4,
				      double *emiss_prior_scalers, double *emiss_prior_scalers_2, double *emiss_prior_scalers_3,
				      double *emiss_prior_scalers_4, int rd, int nr_alphabets);
void add_to_E_multi(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p, int k, int a_size, int normalize,
		    double *subst_mtx, int alphabet, int scoring_method, int use_nr_occ, int alphabet_type, double *emissions);



/************** baum-welch training algorithm ************************/




/* implementation of the baum-welch training algorithm using dirichlet prior mixture to
 * calculate update of emission (and transition) matrices */
void baum_welch_dirichlet_multi(struct hmm_multi_s *hmmp, struct sequence_multi_s *seqsp, int nr_seqs, int annealing, int use_labels,
				int use_transition_pseudo_counts, int use_emission_pseudo_counts,
				int multi_scoring_method, int use_prior)
{
  double *T, *E, *E_2, *E_3, *E_4; /* matrices for the estimated number of times
		 * each transition (T) and emission (E) is used */
  struct forward_s *forw_mtx; /* forward matrix */
  struct backward_s *backw_mtx; /* backward matrix */
  double *forw_scale; /* scaling array */
  int s,p,k,l,a,d; /* loop counters, s loops over the sequences, p over the
			* positions in the sequence, k and l over states, a over the alphabet
			* and d over the distribution groups */
  struct path_element *lp;
  double t_res, t_res_1, t_res_2, t_res_3; /* for temporary results */
  double t_res_4, t_res_5, t_res_6; /* for temporary results */
  double e_res, e_res_1, e_res_2, e_res_3; /* for temporary results */

  int seq_len; /* length of the seqences */
  int a_index, a_index_2, a_index_3, a_index_4; /* holds current letters index in the alphabet */
  struct letter_s *seq, *seq_2, *seq_3, *seq_4; /* pointer to current sequence */
  double old_log_likelihood, new_log_likelihood; /* to calculate when to stop */
  double likelihood; /* temporary variable for calculating likelihood of a sequence */
  int max_nr_iterations, iteration;

  /* dirichlet prior variables */
  struct emission_dirichlet_s *priorp;
  struct emission_dirichlet_s *priorp_2;
  struct emission_dirichlet_s *priorp_3;
  struct emission_dirichlet_s *priorp_4;

  /* simulated annealing variables */
  double temperature;
  double cooling_factor;
  int annealing_status;
 

  /* some initialization */
  old_log_likelihood = 9999.0;
  new_log_likelihood = 9999.0;
  max_nr_iterations = 20;
  iteration = 1;
  if(annealing == YES) {
    temperature = INIT_TEMP;
    cooling_factor = INIT_COOL;
    annealing_status = ACTIVE;
  }
  else {
    annealing_status = DONE;
  }
  

  do {  
    T = (double*)(malloc_or_die(hmmp->nr_v * hmmp->nr_v * 
			       sizeof(double)));
    if(hmmp->alphabet_type == DISCRETE) {
      E = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size * 
				  sizeof(double)));
    }
    else {
      E = (double*)(malloc_or_die(hmmp->nr_v * (hmmp->a_size + 1) * 
				  sizeof(double)));
    }
    if(hmmp->nr_alphabets > 1) {
      if(hmmp->alphabet_type_2 == DISCRETE) {
	E_2 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_2 * 
				      sizeof(double)));
      }
      else {
	E_2 = (double*)(malloc_or_die(hmmp->nr_v * (hmmp->a_size_2 + 1) * 
				      sizeof(double)));
      }
    }
    if(hmmp->nr_alphabets > 2) {
      if(hmmp->alphabet_type_3 == DISCRETE) {
	E_3 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_3 * 
				      sizeof(double)));
      }
      else {
	E_3 = (double*)(malloc_or_die(hmmp->nr_v * (hmmp->a_size_3 + 1) * 
				      sizeof(double)));
      }
    }
    if(hmmp->nr_alphabets > 3) {
      if(hmmp->alphabet_type_4 == DISCRETE) {
	E_4 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_4 * 
				      sizeof(double)));
      }
      else {
	E_4 = (double*)(malloc_or_die(hmmp->nr_v * (hmmp->a_size_4 + 1) * 
				      sizeof(double)));
      }
    }

    old_log_likelihood = new_log_likelihood;
    new_log_likelihood = 0.0;
    for(s = 0; s < nr_seqs; s++) {
      /* Convert sequence to 1...L for easier indexing */
      seq_len = (seqsp + s)->length;
      seq = (struct letter_s*) (malloc_or_die((seq_len + 2) * sizeof(struct letter_s)));
      memcpy(seq+1, (seqsp + s)->seq_1, seq_len * sizeof(struct letter_s));
      if(hmmp->nr_alphabets > 1) {
	seq_2 = (struct letter_s*) (malloc_or_die((seq_len + 2) * sizeof(struct letter_s)));
	memcpy(seq_2+1, (seqsp + s)->seq_2, seq_len * sizeof(struct letter_s));
      }
      if(hmmp->nr_alphabets > 2) {
	seq_3 = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
	memcpy(seq_3+1, (seqsp + s)->seq_3, seq_len * sizeof(struct letter_s));
      }
      if(hmmp->nr_alphabets > 3) {
	seq_4 = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
	memcpy(seq_4+1, (seqsp + s)->seq_4, seq_len * sizeof(struct letter_s));
      }

      
      /* calculate forward and backward matrices */
      forward_multi(hmmp, (seqsp + s)->seq_1,(seqsp + s)->seq_2, (seqsp + s)->seq_3, (seqsp + s)->seq_4,
		    &forw_mtx, &forw_scale, use_labels, multi_scoring_method);
      backward_multi(hmmp, (seqsp + s)->seq_1, (seqsp + s)->seq_2, (seqsp + s)->seq_3, (seqsp + s)->seq_4,
		     &backw_mtx, forw_scale, use_labels, multi_scoring_method);
      /* memory for forw_mtx, scale_mtx and
       * backw_mtx is allocated in the functions */
      
      /* update new_log_likelihood */
      likelihood = log10((forw_mtx +
			  get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob);
      for(k = 0; k <= seq_len; k++) {
	likelihood = likelihood + log10(*(forw_scale + k));
      }
#ifdef DEBUG_BW
      dump_scaling_array(k-1,forw_scale);
      printf("likelihood = %f\n", likelihood);
#endif      
      new_log_likelihood += likelihood;
      
      for(k = 0; k < hmmp->nr_v-1; k++) /* k = from vertex */ {
      lp = *(hmmp->to_trans_array + k);
	while(lp->vertex != END) /* l = to-vertex */ {
	  for(p = 1; p <= seq_len; p++) {
	    
	     /* get alphabet index for c (add replacement letter stuff here) */
	    if(hmmp->alphabet_type == DISCRETE) {
	      a_index = get_alphabet_index(&seq[p], hmmp->alphabet, hmmp->a_size);
	    }	    
	    if(hmmp->nr_alphabets > 1 && hmmp->alphabet_type_2 == DISCRETE) {
	      a_index_2 = get_alphabet_index(&seq_2[p], hmmp->alphabet_2, hmmp->a_size_2);
	    }
	    if(hmmp->nr_alphabets > 2 && hmmp->alphabet_type_3 == DISCRETE) {
	      a_index_3 = get_alphabet_index(&seq_3[p], hmmp->alphabet_3, hmmp->a_size_3);
	    }
	    if(hmmp->nr_alphabets > 3 && hmmp->alphabet_type_4 == DISCRETE) {
	      a_index_4 = get_alphabet_index(&seq_4[p], hmmp->alphabet_4, hmmp->a_size_4);
	    } 

	    /* add T[k][l] contribution for this sequence */
	    add_Tkl_contribution_multi(hmmp, seq+1, seq_2+1, seq_3+1, seq_4+1, forw_mtx, backw_mtx,
				       forw_scale, p, k, lp, a_index, a_index_2, a_index_3, a_index_4, T, use_labels,
				       multi_scoring_method);
	    
	    /* continuous? */
	    
	  }
	  /* move on to next path */
	  while(lp->next != NULL) {
	    lp++;
	  }
	  lp++;
	}
	/* calculate E[k][a] contribution from this sequence */
	if(silent_state_multi(k, hmmp) != 0) {
	  for(p = 1; p <= seq_len; p++) {
	    if(hmmp->alphabet_type == DISCRETE) {
	      a_index = get_alphabet_index(&seq[p], hmmp->alphabet, hmmp->a_size);
	      *(E + get_mtx_index(k, a_index, hmmp->a_size)) +=
		add_Eka_contribution_multi(hmmp, seq+1, forw_mtx, backw_mtx, p, k, multi_scoring_method);
	    }
	    else {
	      add_Eka_contribution_continuous_multi(hmmp, seq+1, forw_mtx, backw_mtx, p, k, multi_scoring_method, E, 1);
	    }
	    if(hmmp->nr_alphabets > 1 && hmmp->alphabet_type_2 == DISCRETE) {
	      a_index_2 = get_alphabet_index(&seq_2[p], hmmp->alphabet_2, hmmp->a_size_2);
	      *(E_2 + get_mtx_index(k, a_index_2, hmmp->a_size_2)) +=
		add_Eka_contribution_multi(hmmp, seq_2+1, forw_mtx, backw_mtx, p, k, multi_scoring_method);
	    }
	    else if(hmmp->nr_alphabets > 1) {
	      add_Eka_contribution_continuous_multi(hmmp, seq_2+1, forw_mtx, backw_mtx, p, k, multi_scoring_method, E_2, 2);
	    }
	    if(hmmp->nr_alphabets > 2 && hmmp->alphabet_type_3 == DISCRETE) {
	      a_index_3 = get_alphabet_index(&seq_3[p], hmmp->alphabet_3, hmmp->a_size_3);
	      *(E_3 + get_mtx_index(k, a_index_3, hmmp->a_size_3)) +=
		add_Eka_contribution_multi(hmmp, seq_3+1, forw_mtx, backw_mtx, p, k, multi_scoring_method);
	    }
	    else  if(hmmp->nr_alphabets > 2) {
	      add_Eka_contribution_continuous_multi(hmmp, seq_3+1, forw_mtx, backw_mtx, p, k, multi_scoring_method, E_3, 3);
	    }
	    if(hmmp->nr_alphabets > 3 && hmmp->alphabet_type_4 == DISCRETE) {
	      a_index_4 = get_alphabet_index(&seq_4[p], hmmp->alphabet_4, hmmp->a_size_4);
	      *(E_4 + get_mtx_index(k, a_index_4, hmmp->a_size_4)) +=
		add_Eka_contribution_multi(hmmp, seq_4+1, forw_mtx, backw_mtx, p, k, multi_scoring_method);
	    }
	    else if(hmmp->nr_alphabets > 3) {
	      add_Eka_contribution_continuous_multi(hmmp, seq_4+1, forw_mtx, backw_mtx, p, k, multi_scoring_method, E_4, 4);
	    }
	  }
	}
      }
      /* some garbage collection */
      free(seq);
      if(hmmp->nr_alphabets > 1) {
	free(seq_2);
      }
      if(hmmp->nr_alphabets > 2) {
	free(seq_3);
      }
      if(hmmp->nr_alphabets > 3) {
	free(seq_4);
      }
      free(forw_mtx);
      free(backw_mtx);
      free(forw_scale); 
    }
    if(verbose == YES) {
      printf("log likelihood rd %d: %f\n", iteration, new_log_likelihood);
    }
    
#ifdef DEBUG_BW2
    dump_T_matrix(hmmp->nr_v, hmmp->nr_v, T);
    dump_E_matrix(hmmp->nr_v, hmmp->a_size, E);
    //dump_E_matrix(hmmp->nr_v, hmmp->a_size_2 + 1, E_2);
#endif
    
    /* check if likelihood change is small enough, then we are done */
    if(fabs(new_log_likelihood - old_log_likelihood) < INNER_BW_THRESHOLD && annealing_status == DONE) {
	break;
    }
    
    /* if simulated annealing is used, scramble results in E and T matrices */
    if(annealing == YES && temperature > ANNEAL_THRESHOLD) {
      anneal_E_matrix_multi(temperature, E, hmmp, 1);
      if(hmmp->nr_alphabets > 1) {
	anneal_E_matrix_multi(temperature, E_2, hmmp, 2);
      }
      if(hmmp->nr_alphabets > 2) {
	anneal_E_matrix_multi(temperature, E_3, hmmp, 3);
      }
      if(hmmp->nr_alphabets > 3) {
	anneal_E_matrix_multi(temperature, E_4, hmmp, 4);
      }
      anneal_T_matrix_multi(temperature, T, hmmp);
      temperature = temperature * cooling_factor;
    }
    
    if(temperature < ANNEAL_THRESHOLD) {
      annealing_status = DONE;
    }
    
    /* recalculate emission expectations according to distribution groups 
     * by simply taking the mean of the expected emissions within this group
     * for each letter in the alphabet and replacing each expectation for the
     * letter with this value for every member of the distribution group */
    recalculate_emiss_expectations_multi(hmmp, E, 1);
    if(hmmp->nr_alphabets > 1) {
      recalculate_emiss_expectations_multi(hmmp, E_2, 2);
    }
    if(hmmp->nr_alphabets > 2) {
      recalculate_emiss_expectations_multi(hmmp, E_3, 3);
    }
    if(hmmp->nr_alphabets > 3) {
      recalculate_emiss_expectations_multi(hmmp, E_4, 4);
    }
    
    /* recalculate transition expectations for tied transitions according
     * to the same scheme as for emission distribution groups */
    recalculate_trans_expectations_multi(hmmp, T);
    
    for(k = 0; k < hmmp->nr_v-1; k++) /* k = from-vertex */ {
      /* update transition matrix */
      if(use_transition_pseudo_counts == YES) {
	update_trans_mtx_pseudocount_multi(hmmp, T, k);
      }
      else {
	update_trans_mtx_std_multi(hmmp, T, k);
      }
      
      
#ifdef DEBUG_PRIORS
      printf("Starting emission matrix update\n");
#endif
      
      /* update emission matrix using Dirichlet prior files if they exist*/
      priorp = *(hmmp->ed_ps + k);
      if(priorp != NULL && use_prior == YES && hmmp->alphabet_type == DISCRETE) {
#ifdef DEBUG_PRIORS
	printf("k = %d\n", k);
	printf("value = %x\n", priorp);
#endif
	update_emiss_mtx_prior_multi(hmmp, E, k, priorp, 1);
      }
      else if(use_emission_pseudo_counts == YES && hmmp->alphabet_type == DISCRETE)
	/* update emissions matrix "normally" when dirichlet file is missing */ {
	update_emiss_mtx_pseudocount_multi(hmmp, E, k, 1);
      }
      else if(hmmp->alphabet_type == DISCRETE) {
	update_emiss_mtx_std_multi(hmmp, E, k, 1);
      }
      else {
	update_emiss_mtx_std_continuous_multi(hmmp, E, k, 1);
      }
      
      if(hmmp->nr_alphabets > 1) {
	priorp = *(hmmp->ed_ps_2 + k);
	if(priorp != NULL && use_prior == YES && hmmp->alphabet_type_2 == DISCRETE) {
	  update_emiss_mtx_prior_multi(hmmp, E_2, k, priorp, 2);
	}
	else if(use_emission_pseudo_counts == YES && hmmp->alphabet_type_2 == DISCRETE)
	  /* update emissions matrix "normally" when dirichlet file is missing */ {
	  update_emiss_mtx_pseudocount_multi(hmmp, E_2, k, 2);
	}
	else if(hmmp->alphabet_type_2 == DISCRETE) {
	  update_emiss_mtx_std_multi(hmmp, E_2, k, 2);
	}
	else {
	  update_emiss_mtx_std_continuous_multi(hmmp, E_2, k, 2);
	}
      }


      if(hmmp->nr_alphabets > 2) {
	priorp = *(hmmp->ed_ps_3 + k);
	if(priorp != NULL && use_prior == YES && hmmp->alphabet_type_3 == DISCRETE) {
	  update_emiss_mtx_prior_multi(hmmp, E_3, k, priorp, 3);
	}
	else if(use_emission_pseudo_counts == YES && hmmp->alphabet_type_3 == DISCRETE)
	  /* update emissions matrix "normally" when dirichlet file is missing */ {
	  update_emiss_mtx_pseudocount_multi(hmmp, E_3, k, 3);
	}
	else if(hmmp->alphabet_type_3 == DISCRETE) {
	  update_emiss_mtx_std_multi(hmmp, E_3, k, 3);
	}
	else {
	  update_emiss_mtx_std_continuous_multi(hmmp, E_3, k, 3);
	}
      }

      if(hmmp->nr_alphabets > 3) {
	priorp = *(hmmp->ed_ps_4 + k);
	if(priorp != NULL && use_prior == YES && hmmp->alphabet_type_4 == DISCRETE) {
	  update_emiss_mtx_prior_multi(hmmp, E_4, k, priorp, 4);
	}
	else if(use_emission_pseudo_counts == YES && hmmp->alphabet_type_4 == DISCRETE) 
	  /* update emissions matrix "normally" when dirichlet file is missing */ {
	  update_emiss_mtx_pseudocount_multi(hmmp, E_4, k, 4);
	}
	else if(hmmp->alphabet_type_4 == DISCRETE) {
	  update_emiss_mtx_std_multi(hmmp, E_4, k, 4);
	}
	else {
	  update_emiss_mtx_std_continuous_multi(hmmp, E_4, k, 4);
	}
      }
    }
   
#ifdef DEBUG_BW
    dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
    dump_emiss_matrix(hmmp->nr_v, hmmp->a_size, hmmp->emissions);
#endif    
    
    /* some garbage collection */
    free(E);
    if(hmmp->nr_alphabets > 1) {
      free(E_2);
    }
    if(hmmp->nr_alphabets > 2) {
      free(E_3);
    }
    if(hmmp->nr_alphabets > 3) {
      free(E_4);
    }
    free(T);
    max_nr_iterations--;
    iteration++;
  }
  while(max_nr_iterations > 0); /* break condition is also when log_likelihood_difference is
				 * smaller than THRESHOLD, checked inside the loop for
				 * better efficiency */
  
  
#ifdef DEBUG_BW
  dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
  dump_emiss_matrix(hmmp->nr_v, hmmp->a_size, hmmp->emissions);
#endif    
  
}

/* implementation of the baum-welch training algorithm using dirichlet prior mixture to
 * calculate update of emission (and transition) matrices and using a multiple sequence
 * alignment as the training sequence */
void msa_baum_welch_dirichlet_multi(struct hmm_multi_s *hmmp, struct msa_sequences_multi_s *msa_seq_infop, int nr_seqs,
				    int annealing, int use_gap_shares, int use_lead_columns, int use_labels,
				    int use_transition_pseudo_counts, int use_emission_pseudo_counts, int normalize,
				    int scoring_method, int use_nr_occ, int multi_scoring_method, double *aa_freqs,
				    double *aa_freqs_2, double *aa_freqs_3, double *aa_freqs_4, int use_prior)
{
  struct msa_sequences_multi_s *msa_seq_infop_start;
  double *T, *E, *E_2, *E_3, *E_4; /* matrices for the estimated number of times
				    * each transition (T) and emission (E) is used */
  struct forward_s *forw_mtx; /* forward matrix */
  struct backward_s *backw_mtx; /* backward matrix */
  double *forw_scale; /* scaling array */
  int s,p,k,l,a,d,i; /* loop counters, s loops over the sequences, p over the
		      * positions in the sequence, k and l over states, a over the alphabet,
		      * d over the distribution groups and i is a slush variable  */
  struct path_element *lp;
  double t_res, t_res_1, t_res_2, t_res_3; /* for temporary results */
  double t_res_4, t_res_5, t_res_6; /* for temporary results */
  double e_res, e_res_1, e_res_2, e_res_3; /* for temporary results */

  int seq_len; /* length of the sequences */
  int a_index; /* holds current letters index in the alphabet */
  double old_log_likelihood, new_log_likelihood; /* to calculate when to stop */
  double likelihood; /* temporary variable for calculating likelihood of a sequence */
  int max_nr_iterations, iteration;
  double Eka_base;
  int query_index; /* index of query seq */

  /* dirichlet prior variables */
  struct emission_dirichlet_s *priorp;
  struct emission_dirichlet_s *priorp_2;
  struct emission_dirichlet_s *priorp_3;
  struct emission_dirichlet_s *priorp_4;

  /* simulated annealing varialbles */
  double temperature;
  double cooling_factor;
  int annealing_status;

  /* help variables for add_to_E */
  int alphabet_nr;
  int alphabet;
  int a_size;
  double *E_cur;
  double *subst_mtx;
  struct msa_letter_s *msa_seq;
  double *tmp_emissions;
  int alphabet_type;

  /* remember start of sequences */
  msa_seq_infop_start = msa_seq_infop;
 
  old_log_likelihood = 9999.0;
  new_log_likelihood = 9999.0;
  max_nr_iterations = 20;
  iteration = 1;
  if(annealing == YES) {
    temperature = INIT_TEMP;
    cooling_factor = INIT_COOL;
    annealing_status = ACTIVE;
  }
  else {
    annealing_status = DONE;
  }
  
#ifdef DEBUG_BW2
  check_for_corrupt_values(hmmp->nr_v, hmmp->a_size, hmmp->emissions , "emiss");
  check_for_corrupt_values(hmmp->nr_v, hmmp->nr_v, hmmp->transitions , "trans");
#endif
  
  do {
#ifdef DEBUG_BW2
    printf("starting baum-welch loop\n");
#endif
    /* initialize matrices */
    T = (double*)(malloc_or_die(hmmp->nr_v * hmmp->nr_v * 
			       sizeof(double)));

    if(hmmp->alphabet_type == DISCRETE) {
      E = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size * 
				  sizeof(double)));
    }
    else {
      E = (double*)(malloc_or_die(hmmp->nr_v * (hmmp->a_size + 1) * 
				  sizeof(double)));
    }
    if(hmmp->nr_alphabets > 1) {
      if(hmmp->alphabet_type_2 == DISCRETE) {
	E_2 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_2 * 
				      sizeof(double)));
      }
      else {
	E_2 = (double*)(malloc_or_die(hmmp->nr_v * (hmmp->a_size_2 + 1) * 
				      sizeof(double)));
      }
    }
    if(hmmp->nr_alphabets > 2) {
      if(hmmp->alphabet_type_3 == DISCRETE) {
	E_3 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_3 * 
				      sizeof(double)));
      }
      else {
	E_3 = (double*)(malloc_or_die(hmmp->nr_v * (hmmp->a_size_3 + 1) * 
				      sizeof(double)));
      }
    }
    if(hmmp->nr_alphabets > 3) {
      if(hmmp->alphabet_type_4 == DISCRETE) {
	E_4 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_4 * 
				      sizeof(double)));
      }
      else {
	E_4 = (double*)(malloc_or_die(hmmp->nr_v * (hmmp->a_size_4 + 1) * 
				      sizeof(double)));
      }
    }

    /* reset sequence pointer */
    msa_seq_infop = msa_seq_infop_start;
    

    old_log_likelihood = new_log_likelihood;
    new_log_likelihood = 0.0;
    
    for(s = 0; s < nr_seqs; s++) {
      if(use_lead_columns == YES) {
	seq_len = msa_seq_infop->nr_lead_columns;
      }
      else {
	seq_len = msa_seq_infop->msa_seq_length;
      }
      
      /* calculate forward and backward matrices
       * memory for forw_mtx, scale_mtx and
       * backw_mtx is allocated in the functions */
#ifdef DEBUG_BW2
      printf("running forward for seq %d\n", s + 1);
#endif
      msa_forward_multi(hmmp, msa_seq_infop, use_lead_columns, use_gap_shares, NO, &forw_mtx, &forw_scale, use_labels, normalize,
			scoring_method, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
#ifdef DEBUG_BW2
      dump_forward_matrix(seq_len + 2, hmmp->nr_v, forw_mtx);
      printf("running backward\n");
#endif
      msa_backward_multi(hmmp, msa_seq_infop, use_lead_columns, use_gap_shares, &backw_mtx, forw_scale, use_labels, normalize,
			 scoring_method, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
#ifdef DEBUG_BW2
      check_for_corrupt_values(seq_len + 2, hmmp->nr_v, forw_mtx , "F");
      check_for_corrupt_values(seq_len + 2, hmmp->nr_v, backw_mtx , "B");
      printf("done with backward\n");
#endif
      /* update new_log_likelihood */
      likelihood = log10((forw_mtx +
			  get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob);
      for(k = 0; k <= seq_len; k++) {
	likelihood = likelihood + log10(*(forw_scale + k));
      }
#ifdef DEBUG_BW2
      dump_scaling_array(k-1,forw_scale);
      printf("likelihood = %f\n", likelihood);
#endif      
      new_log_likelihood += likelihood;
      
      for(k = 0; k < hmmp->nr_v-1; k++) /* k = from vertex */ {
	lp = *(hmmp->to_trans_array + k);
	while(lp->vertex != END) /* l = to-vertex */ {
	  i = 0;
	  while(1) {
	    if(use_lead_columns == NO) {
	      p = i;
	    }
	    else {
	      p = *(msa_seq_infop->lead_columns_start + i);
	    }
	    /* add T[k][l] contribution for the msa-sequence */
	    add_Tkl_contribution_msa_multi(hmmp, msa_seq_infop, forw_mtx, backw_mtx,
					   forw_scale, p, k, lp, T, use_gap_shares, use_lead_columns, i, use_labels, scoring_method,
					   normalize, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
	    i++;
	    if(use_lead_columns == NO) {
	      if(i >= seq_len) {
		break;
	      }
	    }
	    else {
	      if(*(msa_seq_infop->lead_columns_start + i) == END) {
		break;
	      }
	    }
	  }
	  /* move on to next path */
	  while(lp->next != NULL) {
	    lp++;
	  }
	  lp++;
	}
	
	/* calculate E[k][a] contribution from this sequence */
	if(silent_state_multi(k, hmmp) != 0) {
	  i = 0;
	  while(1) {
	    
	    /* get correct index for this letter column */
	    if(use_lead_columns == NO) {
	      p = i;
	    }
	    else {
	      p = *(msa_seq_infop->lead_columns_start + i);
	    }

	    /* get basic scoring result */
	    Eka_base = add_Eka_contribution_msa_multi(hmmp, msa_seq_infop, forw_mtx, backw_mtx, p, k,
						      i, use_lead_columns);
	    
	    /* loop over the alphabets */
	    for(alphabet_nr = 1; alphabet_nr <= hmmp->nr_alphabets; alphabet_nr++) {
	      if(alphabet_nr == 1) {
		alphabet = hmmp->alphabet;
		subst_mtx = hmmp->subst_mtx;
		a_size = hmmp->a_size;
		E_cur = E;
		msa_seq = msa_seq_infop->msa_seq_1;
		alphabet_type = hmmp->alphabet_type;
		tmp_emissions = hmmp->emissions;
	      }
	      else if(alphabet_nr == 2) {
		alphabet = hmmp->alphabet_2;
		subst_mtx = hmmp->subst_mtx_2;
		a_size = hmmp->a_size_2;
		E_cur = E_2;
		msa_seq = msa_seq_infop->msa_seq_2;
		alphabet_type = hmmp->alphabet_type_2;
		tmp_emissions = hmmp->emissions_2;
	      }
	      else if(alphabet_nr == 3) {
		alphabet = hmmp->alphabet_3;
		subst_mtx = hmmp->subst_mtx_3;
		a_size = hmmp->a_size_3;
		E_cur = E_3;
		msa_seq = msa_seq_infop->msa_seq_3;
		alphabet_type = hmmp->alphabet_type_3;
		tmp_emissions = hmmp->emissions_3;
	      }
	      else if(alphabet_nr == 4) {
		alphabet = hmmp->alphabet_4;
		subst_mtx = hmmp->subst_mtx_4;
		a_size = hmmp->a_size_4;
		E_cur = E_4;
		msa_seq = msa_seq_infop->msa_seq_4;
		alphabet_type = hmmp->alphabet_type_4;
		tmp_emissions = hmmp->emissions_4;
	      }
	      
	      /* get result and add to matrix according to scoring method */
	      add_to_E_multi(E_cur, Eka_base, msa_seq, p, k, a_size, normalize, subst_mtx,
			     alphabet, scoring_method, use_nr_occ, alphabet_type, tmp_emissions);
	    }
	    /* update loop index, check if we are done */
	    i++;
	    if(use_lead_columns == NO) {
	      if(i >= seq_len) {
		break;
	      }
	    }
	    else {
	      if(*(msa_seq_infop->lead_columns_start + i) == END) {
		break;
	      }
	    }
	  }
	}
      }
      
      /* some garbage collection */
      free(forw_mtx);
      free(backw_mtx);
      free(forw_scale);

      msa_seq_infop++;
    }
	
    if(verbose == YES) {
      printf("log likelihood rd %d: %f\n", iteration, new_log_likelihood);
    }
   
     
 
#ifdef DEBUG_BW2
    dump_T_matrix(hmmp->nr_v, hmmp->nr_v, T);
    dump_E_matrix(hmmp->nr_v, hmmp->a_size, E);
    if(hmmp->nr_alphabets > 1) {
      dump_E_matrix(hmmp->nr_v, hmmp->a_size_2 + 1, E_2);
    }
#endif
    /* check if likelihood change is small enough, then we are done */
    if(fabs(new_log_likelihood - old_log_likelihood) < INNER_BW_THRESHOLD && annealing_status == DONE) {
      break;
    }
    
    /* if simulated annealing is used, scramble results in E and T matrices */
    if(annealing == YES && temperature > ANNEAL_THRESHOLD) {
      anneal_E_matrix_multi(temperature, E, hmmp, 1);
      if(hmmp->nr_alphabets > 1) {
	anneal_E_matrix_multi(temperature, E_2, hmmp, 2);
      }
      if(hmmp->nr_alphabets > 2) {
	anneal_E_matrix_multi(temperature, E_3, hmmp, 3);
      }
      if(hmmp->nr_alphabets > 3) {
	anneal_E_matrix_multi(temperature, E_4, hmmp, 4);
      }
      anneal_T_matrix_multi(temperature, T, hmmp);
      temperature = temperature * cooling_factor;
    }

    if(temperature < ANNEAL_THRESHOLD) {
      annealing_status = DONE;
    }

    /* recalculate emission expectations according to distribution groups 
     * by simply taking the mean of the expected emissions within this group
     * for each letter in the alphabet and replacing each expectation for the
     * letter with this value for every member of the distribution group */
    recalculate_emiss_expectations_multi(hmmp, E, 1);
    if(hmmp->nr_alphabets > 1) {
      recalculate_emiss_expectations_multi(hmmp, E_2, 2);
    }
    if(hmmp->nr_alphabets > 2) {
      recalculate_emiss_expectations_multi(hmmp, E_3, 3);
    }
    if(hmmp->nr_alphabets > 3) {
      recalculate_emiss_expectations_multi(hmmp, E_4, 4);
    }
    
    /* recalculate transition expectations for tied transitions according
     * to the same scheme as for emission distribution groups */
    recalculate_trans_expectations_multi(hmmp, T);


    for(k = 0; k < hmmp->nr_v-1; k++) /* k = from-vertex */ {
      /* update transition matrix */
      if(use_transition_pseudo_counts == YES) {
	update_trans_mtx_pseudocount_multi(hmmp, T, k);
      }
      else {
	update_trans_mtx_std_multi(hmmp, T, k);
      }
      
      
#ifdef DEBUG_PRIORS
      printf("Starting emission matrix update\n");
#endif
      

      /* update emission matrix using Dirichlet prior files if they exist*/
      priorp = *(hmmp->ed_ps + k);
      if(priorp != NULL && use_prior == YES && hmmp->alphabet_type == DISCRETE) {
#ifdef DEBUG_PRIORS	
	printf("k = %d\n", k);
	printf("value = %x\n", priorp);
#endif
	update_emiss_mtx_prior_multi(hmmp, E, k, priorp, 1);
      }
      else if(use_emission_pseudo_counts == YES && hmmp->alphabet_type == DISCRETE)
	/* update emissions matrix "normally" when dirichlet file is missing */ {
	update_emiss_mtx_pseudocount_multi(hmmp, E, k, 1);
      }
      else if(hmmp->alphabet_type == DISCRETE) {
	update_emiss_mtx_std_multi(hmmp, E, k, 1);
      }
      else {
	update_emiss_mtx_std_continuous_multi(hmmp, E, k, 1);
      }
      

      if(hmmp->nr_alphabets > 1) {
	priorp = *(hmmp->ed_ps_2 + k);
	if(priorp != NULL && use_prior == YES && hmmp->alphabet_type_2 == DISCRETE) {
	  update_emiss_mtx_prior_multi(hmmp, E_2, k, priorp, 2);
	}
	else if(use_emission_pseudo_counts == YES && hmmp->alphabet_type_2 == DISCRETE)
	  /* update emissions matrix "normally" when dirichlet file is missing */ {
	  update_emiss_mtx_pseudocount_multi(hmmp, E_2, k, 2);
	}
	else if(hmmp->alphabet_type_2 == DISCRETE) {
	  update_emiss_mtx_std_multi(hmmp, E_2, k, 2);
	}
	else {
	  update_emiss_mtx_std_continuous_multi(hmmp, E_2, k, 2);
	}
      }

      if(hmmp->nr_alphabets > 2) {
	priorp = *(hmmp->ed_ps_3 + k);
	if(priorp != NULL && use_prior == YES && hmmp->alphabet_type_3 == DISCRETE) {
	  update_emiss_mtx_prior_multi(hmmp, E_3, k, priorp, 3);
	}
	else if(use_emission_pseudo_counts == YES && hmmp->alphabet_type_3 == DISCRETE)
	  /* update emissions matrix "normally" when dirichlet file is missing */ {
	  update_emiss_mtx_pseudocount_multi(hmmp, E_3, k, 3);
	}
	else if(hmmp->alphabet_type_3 == DISCRETE) {
	  update_emiss_mtx_std_multi(hmmp, E_3, k, 3);
	}
	else {
	  update_emiss_mtx_std_continuous_multi(hmmp, E_3, k, 3);
	}
      }

      if(hmmp->nr_alphabets > 3) {
	priorp = *(hmmp->ed_ps_4 + k);
	if(priorp != NULL && use_prior == YES && hmmp->alphabet_type_4 == DISCRETE) {
	  update_emiss_mtx_prior_multi(hmmp, E_4, k, priorp, 4);
	}
	else if(use_emission_pseudo_counts == YES && hmmp->alphabet_type_4 == DISCRETE) 
	  /* update emissions matrix "normally" when dirichlet file is missing */ {
	  update_emiss_mtx_pseudocount_multi(hmmp, E_4, k, 4);
	}
	else if(hmmp->alphabet_type_4 == DISCRETE) {
	  update_emiss_mtx_std_multi(hmmp, E_4, k, 4);
	}
	else {
	  update_emiss_mtx_std_continuous_multi(hmmp, E_4, k, 4);
	}
      }
    }
    
#ifdef DEBUG_BW2
    dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
    dump_emiss_matrix(hmmp->nr_v, hmmp->a_size, hmmp->emissions);
    if(hmmp->nr_alphabets > 1) {
      dump_emiss_matrix(hmmp->nr_v, hmmp->a_size_2, hmmp->emissions_2);
    }
#endif    
    
    /* some garbage collection */
    free(E);
    if(hmmp->nr_alphabets > 1) {
      free(E_2);
    }
    if(hmmp->nr_alphabets > 2) {
      free(E_3);
    }
    if(hmmp->nr_alphabets > 3) {
      free(E_4);
    }
    free(T);
    max_nr_iterations--;
    iteration++;
#ifdef DEBIG_BW2
    printf("end of baum-welch-loop\n");
#endif
  }
  while(max_nr_iterations > 0); /* break condition is also when log_likelihood_difference is
				 * smaller than THRESHOLD, checked inside the loop for
				 * better efficiency */
  
  
#ifdef DEBUG_BW
  dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
  dump_emiss_matrix(hmmp->nr_v, hmmp->a_size, hmmp->emissions);
#endif    
}


/* implementation of the conditional maximum likelihood version of the
 * baum-welch training algorithm using dirichlet prior mixture to
 * calculate update of emission (and transition) matrices */
void extended_baum_welch_dirichlet_multi(struct hmm_multi_s *hmmp, struct sequence_multi_s *seqsp, int nr_seqs,
					 int annealing, int use_labels,
					 int use_transition_pseudo_counts, int use_emission_pseudo_counts,
					 int multi_scoring_method, int use_prior)
{
  double *T, *E, *E_2, *E_3, *E_4; /* matrices for the estimated number of times
		 * each transition (T) and emission (E) is used */
  double *T_lab, *E_lab, *E_lab_2, *E_lab_3, *E_lab_4, *T_ulab, *E_ulab, *E_ulab_2, *E_ulab_3, *E_ulab_4;
  struct forward_s *forw_mtx; /* forward matrix */
  struct backward_s *backw_mtx; /* backward matrix */
  double *forw_scale; /* scaling array */
  int s,p,k,l,a,d; /* loop counters, s loops over the sequences, p over the
			* positions in the sequence, k and l over states, a over the alphabet
			* and d over the distribution groups */
  struct path_element *lp;
  double t_res, t_res_1, t_res_2, t_res_3; /* for temporary results */
  double t_res_4, t_res_5, t_res_6; /* for temporary results */
  double e_res, e_res_1, e_res_2, e_res_3; /* for temporary results */
  double t_res_ulab;

  int seq_len; /* length of the seqences */
  int a_index, a_index_2, a_index_3, a_index_4; /* holds current letters index in the alphabet */
  struct letter_s *seq, *seq_2, *seq_3, *seq_4; /* pointer to current sequence */
  double old_log_likelihood_lab, new_log_likelihood_lab;
  double old_log_likelihood_ulab, new_log_likelihood_ulab; /* to calculate when to stop */
  double likelihood; /* temporary variable for calculating likelihood of a sequence */
  int max_nr_iterations, iteration;

  /* dirichlet prior variables */
  struct emission_dirichlet_s *priorp;

  /* simulated annealing variables */
  double temperature;
  double cooling_factor;
  int annealing_status;
  
  /* some initialization */
  old_log_likelihood_lab = 9999.0;
  new_log_likelihood_lab = 9999.0;
  old_log_likelihood_ulab = 9999.0;
  new_log_likelihood_ulab = 9999.0;
  max_nr_iterations = 70;
  iteration = 1;
  if(annealing == YES) {
    temperature = INIT_TEMP;
    cooling_factor = INIT_COOL;
    annealing_status = ACTIVE;
  }
  else {
    annealing_status = DONE;
  }
  

  do {
    /* allocate per iteration matrices */
    T = (double*)(malloc_or_die(hmmp->nr_v * hmmp->nr_v * sizeof(double)));
    E = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size * sizeof(double)));
    T_lab = (double*)(malloc_or_die(hmmp->nr_v * hmmp->nr_v * sizeof(double)));
    E_lab = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size * sizeof(double)));
    T_ulab = (double*)(malloc_or_die(hmmp->nr_v * hmmp->nr_v * sizeof(double)));
    E_ulab = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size * sizeof(double))); 

    if(hmmp->nr_alphabets > 1) {
      E_2 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_2 * sizeof(double)));
      E_lab_2 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_2 * sizeof(double)));
      E_ulab_2 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_2 * sizeof(double)));
    }
    if(hmmp->nr_alphabets > 2) {
      E_3 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_3 * sizeof(double)));
      E_lab_3 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_3 * sizeof(double)));
      E_ulab_3 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_3 * sizeof(double)));
    }
    if(hmmp->nr_alphabets > 3) {
      E_4 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_4 * sizeof(double)));
      E_lab_4 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_4 * sizeof(double)));
      E_ulab_4 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_4 * sizeof(double)));
    }

    old_log_likelihood_ulab = new_log_likelihood_ulab;
    new_log_likelihood_ulab = 0.0;
    old_log_likelihood_lab = new_log_likelihood_lab;
    new_log_likelihood_lab = 0.0;
    for(s = 0; s < nr_seqs; s++) {
      /* Convert sequence to 1...L for easier indexing */
      seq_len = (seqsp + s)->length;
      seq = (struct letter_s*) (malloc_or_die((seq_len + 2) * sizeof(struct letter_s)));
      memcpy(seq+1, (seqsp + s)->seq_1, seq_len * sizeof(struct letter_s));
      if(hmmp->nr_alphabets > 1) {
	seq_2 = (struct letter_s*) (malloc_or_die((seq_len + 2) * sizeof(struct letter_s)));
	memcpy(seq_2+1, (seqsp + s)->seq_2, seq_len * sizeof(struct letter_s));
      }
      if(hmmp->nr_alphabets > 2) {
	seq_3 = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
	memcpy(seq_3+1, (seqsp + s)->seq_3, seq_len * sizeof(struct letter_s));
      }
      if(hmmp->nr_alphabets > 3) {
	seq_4 = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
	memcpy(seq_4+1, (seqsp + s)->seq_4, seq_len * sizeof(struct letter_s));
      }

      
      /* calculate forward and backward matrices */
      forward_multi(hmmp, (seqsp + s)->seq_1,(seqsp + s)->seq_2, (seqsp + s)->seq_3, (seqsp + s)->seq_4,
		    &forw_mtx, &forw_scale, NO, multi_scoring_method);
      backward_multi(hmmp, (seqsp + s)->seq_1, (seqsp + s)->seq_2, (seqsp + s)->seq_3, (seqsp + s)->seq_4,
		     &backw_mtx, forw_scale, NO, multi_scoring_method);
      
      /* memory for forw_mtx, scale_mtx and
       * backw_mtx is allocated in the functions */
      
      /* update new_log_likelihood */
      likelihood = log10((forw_mtx +
			  get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob);
      for(k = 0; k <= seq_len; k++) {
	likelihood = likelihood + log10(*(forw_scale + k));
      }
#ifdef DEBUG_BW
      dump_scaling_array(k-1,forw_scale);
      printf("likelihood = %f\n", likelihood);
#endif      
      new_log_likelihood_ulab += likelihood;
      
      for(k = 0; k < hmmp->nr_v-1; k++) /* k = from vertex */ {
	lp = *(hmmp->to_trans_array + k);
	while(lp->vertex != END) /* l = to-vertex */ {
	  for(p = 1; p <= seq_len; p++) {
	    
	    /* get alphabet index for c*/
	    a_index = get_alphabet_index(&seq[p], hmmp->alphabet, hmmp->a_size);
	    if(hmmp->nr_alphabets > 1) {
	      a_index_2 = get_alphabet_index(&seq_2[p], hmmp->alphabet_2, hmmp->a_size_2);
	    }
	    if(hmmp->nr_alphabets > 2) {
	      a_index_3 = get_alphabet_index(&seq_3[p], hmmp->alphabet_3, hmmp->a_size_3);
	    }
	    if(hmmp->nr_alphabets > 3) {
	      a_index_4 = get_alphabet_index(&seq_4[p], hmmp->alphabet_4, hmmp->a_size_4);
	    } 
	    
	    /* add T[k][l] contribution for this sequence */
	    add_Tkl_contribution_multi(hmmp, seq+1, seq_2+1, seq_3+1, seq_4+1, forw_mtx, backw_mtx,
				       forw_scale, p, k, lp, a_index, a_index_2, a_index_3, a_index_4, T_ulab, NO,
				       multi_scoring_method);
	  }
	  /* move on to next path */
	  while(lp->next != NULL) {
	    lp++;
	  }
	  lp++;
	}

	/* calculate E[k][a] contribution from this sequence */
	if(silent_state_multi(k, hmmp) != 0) {
	  for(p = 1; p <= seq_len; p++) {
	    a_index = get_alphabet_index(&seq[p], hmmp->alphabet, hmmp->a_size);
	    if(hmmp->nr_alphabets > 1) {
	      a_index_2 = get_alphabet_index(&seq_2[p], hmmp->alphabet_2, hmmp->a_size_2);
	    }
	    if(hmmp->nr_alphabets > 2) {
	      a_index_3 = get_alphabet_index(&seq_3[p], hmmp->alphabet_3, hmmp->a_size_3);
	    }
	    if(hmmp->nr_alphabets > 3) {
	      a_index_4 = get_alphabet_index(&seq_4[p], hmmp->alphabet_4, hmmp->a_size_4);
	    }
	    /* get result and add to matrix */
	    *(E_ulab + get_mtx_index(k, a_index, hmmp->a_size)) +=
	      add_Eka_contribution_multi(hmmp, seq+1, forw_mtx, backw_mtx, p, k, multi_scoring_method);
	    if(hmmp->nr_alphabets > 1) {
	      *(E_ulab_2 + get_mtx_index(k, a_index_2, hmmp->a_size_2)) +=
		add_Eka_contribution_multi(hmmp, seq_2+1, forw_mtx, backw_mtx, p, k, multi_scoring_method);
	    }
	    if(hmmp->nr_alphabets > 2) {
	      *(E_ulab_3 + get_mtx_index(k, a_index_3, hmmp->a_size_3)) +=
		add_Eka_contribution_multi(hmmp, seq_3+1, forw_mtx, backw_mtx, p, k, multi_scoring_method);
	    }
	    if(hmmp->nr_alphabets > 3) {
	      *(E_ulab_4 + get_mtx_index(k, a_index_4, hmmp->a_size_4)) +=
		add_Eka_contribution_multi(hmmp, seq_4+1, forw_mtx, backw_mtx, p, k, multi_scoring_method);
	    }
	  }
	}
      }

      t_res_ulab = (forw_mtx + get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob;
      /* some garbage collection */
      free(forw_mtx);
      free(backw_mtx);
      free(forw_scale);

     
      /********* calculations using labels *************/
      
      /* calculate forward and backward matrices */
      forward_multi(hmmp, (seqsp + s)->seq_1,(seqsp + s)->seq_2, (seqsp + s)->seq_3, (seqsp + s)->seq_4,
		    &forw_mtx, &forw_scale, YES, multi_scoring_method);
      backward_multi(hmmp, (seqsp + s)->seq_1, (seqsp + s)->seq_2, (seqsp + s)->seq_3, (seqsp + s)->seq_4,
		     &backw_mtx, forw_scale, YES, multi_scoring_method);
      /* memory for forw_mtx, scale_mtx and
       * backw_mtx is allocated in the functions */
      
      /* update new_log_likelihood */
      likelihood = log10((forw_mtx +
			  get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob);
      for(k = 0; k <= seq_len; k++) {
	likelihood = likelihood + log10(*(forw_scale + k));
      }
#ifdef DEBUG_BW
      dump_scaling_array(k-1,forw_scale);
      printf("likelihood = %f\n", likelihood);
#endif      
      new_log_likelihood_lab += likelihood;
      
      for(k = 0; k < hmmp->nr_v-1; k++) /* k = from vertex */ {
	lp = *(hmmp->to_trans_array + k);
	while(lp->vertex != END) /* l = to-vertex */ {
	  for(p = 1; p <= seq_len; p++) {
	    
	     /* get alphabet index for c*/
	    a_index = get_alphabet_index(&seq[p], hmmp->alphabet, hmmp->a_size);
	    if(hmmp->nr_alphabets > 1) {
	      a_index_2 = get_alphabet_index(&seq_2[p], hmmp->alphabet_2, hmmp->a_size_2);
	    }
	    if(hmmp->nr_alphabets > 2) {
	      a_index_3 = get_alphabet_index(&seq_3[p], hmmp->alphabet_3, hmmp->a_size_3);
	    }
	    if(hmmp->nr_alphabets > 3) {
	      a_index_4 = get_alphabet_index(&seq_4[p], hmmp->alphabet_4, hmmp->a_size_4);
	    } 
	    
	    /* add T[k][l] contribution for this sequence */
	    add_Tkl_contribution_multi(hmmp, seq+1, seq_2+1, seq_3+1, seq_4+1, forw_mtx, backw_mtx,
				       forw_scale, p, k, lp, a_index, a_index_2, a_index_3, a_index_4, T_lab, YES,
				       multi_scoring_method);
	  }
	  /* move on to next path */
	  while(lp->next != NULL) {
	    lp++;
	  }
	  lp++;
	}
	/* calculate E[k][a] contribution from this sequence */
	if(silent_state_multi(k, hmmp) != 0) {
	  for(p = 1; p <= seq_len; p++) {
	    a_index = get_alphabet_index(&seq[p], hmmp->alphabet, hmmp->a_size);
	    if(hmmp->nr_alphabets > 1) {
	      a_index_2 = get_alphabet_index(&seq_2[p], hmmp->alphabet_2, hmmp->a_size_2);
	    }
	    if(hmmp->nr_alphabets > 2) {
	      a_index_3 = get_alphabet_index(&seq_3[p], hmmp->alphabet_3, hmmp->a_size_3);
	    }
	    if(hmmp->nr_alphabets > 3) {
	      a_index_4 = get_alphabet_index(&seq_4[p], hmmp->alphabet_4, hmmp->a_size_4);
	    } 
	    /* get result and add to matrix */
	    *(E_lab + get_mtx_index(k, a_index, hmmp->a_size)) +=
	      add_Eka_contribution_multi(hmmp, seq+1, forw_mtx, backw_mtx, p, k, multi_scoring_method);
	    if(hmmp->nr_alphabets > 1) {
	      *(E_lab_2 + get_mtx_index(k, a_index_2, hmmp->a_size_2)) +=
		add_Eka_contribution_multi(hmmp, seq_2+1, forw_mtx, backw_mtx, p, k, multi_scoring_method);
	    }
	    if(hmmp->nr_alphabets > 2) {
	      *(E_lab_3 + get_mtx_index(k, a_index_3, hmmp->a_size_3)) +=
		add_Eka_contribution_multi(hmmp, seq_3+1, forw_mtx, backw_mtx, p, k, multi_scoring_method);
	    }
	    if(hmmp->nr_alphabets > 3) {
	      *(E_lab_4 + get_mtx_index(k, a_index_4, hmmp->a_size_4)) +=
		add_Eka_contribution_multi(hmmp, seq_4+1, forw_mtx, backw_mtx, p, k, multi_scoring_method);
	    }
	  }
	}
      }

      /* some garbage collection */
      free(seq);
      if(hmmp->nr_alphabets > 1) {
	free(seq_2);
      }
      if(hmmp->nr_alphabets > 2) {
	free(seq_3);
      }
      if(hmmp->nr_alphabets > 3) {
	free(seq_4);
      }
      free(forw_mtx);
      free(backw_mtx);
      free(forw_scale);

      
    }
    if(verbose == YES) { 
      printf("log likelihood diff rd %d: %f\n", iteration, new_log_likelihood_ulab - new_log_likelihood_lab);
    }
    
#ifdef DEBUG_BW
    dump_T_matrix(hmmp->nr_v, hmmp->nr_v, T);
    dump_E_matrix(hmmp->nr_v, hmmp->a_size, E);
#endif
    
    /* recalculate emission expectations according to distribution groups 
     * by simply taking the mean of the expected emissions within this group
     * for each letter in the alphabet and replacing each expectation for the
     * letter with this value for every member of the distribution group */
    recalculate_emiss_expectations_multi(hmmp, E_lab, 1);
    if(hmmp->nr_alphabets > 1) {
      recalculate_emiss_expectations_multi(hmmp, E_lab_2, 2);
    }
    if(hmmp->nr_alphabets > 2) {
      recalculate_emiss_expectations_multi(hmmp, E_lab_3, 3);
    }
    if(hmmp->nr_alphabets > 3) {
      recalculate_emiss_expectations_multi(hmmp, E_lab_4, 4);
    }

    recalculate_emiss_expectations_multi(hmmp, E_ulab, 1);
    if(hmmp->nr_alphabets > 1) {
      recalculate_emiss_expectations_multi(hmmp, E_ulab_2, 2);
    }
    if(hmmp->nr_alphabets > 2) {
      recalculate_emiss_expectations_multi(hmmp, E_ulab_3, 3);
    }
    if(hmmp->nr_alphabets > 3) {
      recalculate_emiss_expectations_multi(hmmp, E_ulab_4, 4);
    }
    
    /* recalculate transition expectations for tied transitions according
     * to the same scheme as for emission distribution groups */
    recalculate_trans_expectations_multi(hmmp, T_lab);
    recalculate_trans_expectations_multi(hmmp, T_ulab);
    
    
    /* update real T end E matrices */
    calculate_TE_contributions_multi(T, E, E_2, E_3, E_4, T_lab, E_lab, E_lab_2, E_lab_3, E_lab_4, T_ulab, E_ulab,
				     E_ulab_2, E_ulab_3, E_ulab_4, hmmp->emissions, hmmp->emissions_2, hmmp->emissions_3,
				     hmmp->emissions_4, hmmp->transitions, hmmp->nr_v, hmmp->a_size,
				     hmmp->a_size_2, hmmp->a_size_3, hmmp->a_size_4, hmmp->vertex_emiss_prior_scalers,
				     hmmp->vertex_emiss_prior_scalers_2, hmmp->vertex_emiss_prior_scalers_3,
				     hmmp->vertex_emiss_prior_scalers_4, iteration, hmmp->nr_alphabets);
    
    /* check if likelihood change is small enough, then we are done */
    if(fabs((new_log_likelihood_ulab - new_log_likelihood_lab) - (old_log_likelihood_ulab - old_log_likelihood_lab))
       < CML_THRESHOLD && annealing_status == DONE) {
      break;
    }
    
    /* if simulated annealing is used, scramble results in E and T matrices */
    if(annealing == YES && temperature > ANNEAL_THRESHOLD) {
      anneal_E_matrix_multi(temperature, E, hmmp, 1);
      if(hmmp->nr_alphabets > 1) {
	anneal_E_matrix_multi(temperature, E_2, hmmp, 2);
      }
      if(hmmp->nr_alphabets > 2) {
	anneal_E_matrix_multi(temperature, E_3, hmmp, 3);
      }
      if(hmmp->nr_alphabets > 3) {
	anneal_E_matrix_multi(temperature, E_4, hmmp, 4);
      }
      anneal_T_matrix_multi(temperature, T, hmmp);
      temperature = temperature * cooling_factor;
    }

    if(temperature < ANNEAL_THRESHOLD) {
      annealing_status = DONE;
    }

    for(k = 0; k < hmmp->nr_v-1; k++) /* k = from-vertex */ {
      /* update transition matrix */
      if(use_transition_pseudo_counts == YES) {
	update_trans_mtx_pseudocount_multi(hmmp, T, k);
      }
      else {
	update_trans_mtx_std_multi(hmmp, T, k);
      }
      
      
#ifdef DEBUG_PRIORS
      printf("Starting emission matrix update\n");
#endif
      
      /* update emission matrix using Dirichlet prior files if they exist*/
      priorp = *(hmmp->ed_ps + k);
      if(priorp != NULL && use_prior == YES) {
#ifdef DEBUG_PRIORS	
	printf("k = %d\n", k);
	printf("value = %x\n", priorp);
#endif
	update_emiss_mtx_prior_multi(hmmp, E, k, priorp, 1);
      }
      else if(use_emission_pseudo_counts == YES) /* update emissions matrix "normally" when dirichlet file is missing */ {
	update_emiss_mtx_pseudocount_multi(hmmp, E, k, 1);
      }
      else {
	update_emiss_mtx_std_multi(hmmp, E, k, 1);
      }
      

      if(hmmp->nr_alphabets > 1) {
	priorp = *(hmmp->ed_ps_2 + k);
	if(priorp != NULL && use_prior == YES) {
	  update_emiss_mtx_prior_multi(hmmp, E_2, k, priorp, 2);
	}
	else if(use_emission_pseudo_counts == YES) /* update emissions matrix "normally" when dirichlet file is missing */ {
	  update_emiss_mtx_pseudocount_multi(hmmp, E_2, k, 2);
	}
	else {
	  update_emiss_mtx_std_multi(hmmp, E_2, k, 2);
	}
      }

      if(hmmp->nr_alphabets > 2) {
	priorp = *(hmmp->ed_ps_3 + k);
	if(priorp != NULL && use_prior == YES) {
	  update_emiss_mtx_prior_multi(hmmp, E_3, k, priorp, 3);
	}
	else if(use_emission_pseudo_counts == YES) /* update emissions matrix "normally" when dirichlet file is missing */ {
	  update_emiss_mtx_pseudocount_multi(hmmp, E_3, k, 3);
	}
	else {
	  update_emiss_mtx_std_multi(hmmp, E_3, k, 3);
	}
      }

      if(hmmp->nr_alphabets > 3) {
	priorp = *(hmmp->ed_ps_4 + k);
	if(priorp != NULL && use_prior == YES) {
	  update_emiss_mtx_prior_multi(hmmp, E_4, k, priorp, 4);
	}
	else if(use_emission_pseudo_counts == YES) /* update emissions matrix "normally" when dirichlet file is missing */ {
	  update_emiss_mtx_pseudocount_multi(hmmp, E_4, k, 4);
	}
	else {
	  update_emiss_mtx_std_multi(hmmp, E_4, k, 4);
	}
      }
    }
    
#ifdef DEBUG_BW
    dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
    dump_emiss_matrix(hmmp->nr_v, hmmp->a_size, hmmp->emissions);
#endif    
    
    /* some garbage collection */
    free(E);
    if(hmmp->nr_alphabets > 1) {
      free(E_2);
    }
    if(hmmp->nr_alphabets > 2) {
      free(E_3);
    }
    if(hmmp->nr_alphabets > 3) {
      free(E_4);
    }
    free(T);
    free(T_lab);
    free(E_lab);
    if(hmmp->nr_alphabets > 1) {
      free(E_lab_2);
    }
    if(hmmp->nr_alphabets > 2) {
      free(E_lab_3);
    }
    if(hmmp->nr_alphabets > 3) {
      free(E_lab_4);
    }
    free(T_ulab);
    free(E_ulab);
    if(hmmp->nr_alphabets > 1) {
      free(E_ulab_2);
    }
    if(hmmp->nr_alphabets > 2) {
      free(E_ulab_3);
    }
    if(hmmp->nr_alphabets > 3) {
      free(E_ulab_4);
    }
    max_nr_iterations--;
    iteration++;
  }
  while(max_nr_iterations > 0); /* break condition is also when log_likelihood_difference is
				 * smaller than THRESHOLD, checked inside the loop for
				 * better efficiency */
#ifdef DEBUG_BW2
  printf("exiting\n");
#endif
#ifdef DEBUG_BW
  dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
  dump_emiss_matrix(hmmp->nr_v, hmmp->a_size, hmmp->emissions);
#endif

}



/* implementation of the baum-welch training algorithm using dirichlet prior mixture to
 * calculate update of emission (and transition) matrices and using a multiple sequence
 * alignment as the training sequence */
void extended_msa_baum_welch_dirichlet_multi(struct hmm_multi_s *hmmp, struct msa_sequences_multi_s *msa_seq_infop,
					     int nr_seqs, int annealing,
					     int use_gap_shares, int use_lead_columns, int use_labels, int use_transition_pseudo_counts,
					     int use_emission_pseudo_counts, int normalize, int scoring_method, int use_nr_occ,
					     int multi_scoring_method, double *aa_freqs,
					     double *aa_freqs_2, double *aa_freqs_3, double *aa_freqs_4, int use_prior)
{
  struct msa_sequences_multi_s *msa_seq_infop_start;
  double *T, *E, *E_2, *E_3, *E_4; /* matrices for the estimated number of times
		 * each transition (T) and emission (E) is used */
  double *T_lab, *E_lab, *E_lab_2, *E_lab_3, *E_lab_4, *T_ulab, *E_ulab, *E_ulab_2, *E_ulab_3, *E_ulab_4;
  struct forward_s *forw_mtx; /* forward matrix */
  struct backward_s *backw_mtx; /* backward matrix */
  double *forw_scale; /* scaling array */
  int s,p,k,l,a,d,i; /* loop counters, s loops over the sequences, p over the
		    * positions in the sequence, k and l over states, a over the alphabet,
		    * d over the distribution groups and i is a slush variable  */
  struct path_element *lp;
  double t_res, t_res_1, t_res_2, t_res_3; /* for temporary results */
  double t_res_4, t_res_5, t_res_6; /* for temporary results */
  double e_res, e_res_1, e_res_2, e_res_3; /* for temporary results */

  int seq_len; /* length of the seqences */
  int a_index, a_index_2, a_index_3, a_index_4; /* holds current letters index in the alphabet */
  struct letter_s *seq; /* pointer to current sequence */
  double old_log_likelihood_lab, new_log_likelihood_lab;
  double old_log_likelihood_ulab, new_log_likelihood_ulab; /* to calculate when to stop */
  double likelihood; /* temporary variable for calculating likelihood of a sequence */
  int max_nr_iterations, iteration;
  double Eka_base;

  /* dirichlet prior variables */
  struct emission_dirichlet_s *priorp;

  /* simulated annealing varialbles */
  double temperature;
  double cooling_factor;
  int annealing_status;

   /* help variables for add_to_E */
  int alphabet_nr;
  int alphabet;
  int a_size;
  double *E_cur;
  double *subst_mtx;
  struct msa_letter_s *msa_seq;

  /* remember start of sequence pointer */
  msa_seq_infop_start = msa_seq_infop;

  old_log_likelihood_lab = 9999.0;
  new_log_likelihood_lab = 9999.0;
  old_log_likelihood_ulab = 9999.0;
  new_log_likelihood_ulab = 9999.0;
  max_nr_iterations = 70;
  iteration = 1;
  if(annealing == YES) {
    temperature = INIT_TEMP;
    cooling_factor = INIT_COOL;
    annealing_status = ACTIVE;
  }
  else {
    annealing_status = DONE;
  }
  
  do {
#ifdef DEBUG_BW2
    printf("starting baum-welch loop\n");
#endif
    /* allocate per iteration matrices */
    T = (double*)(malloc_or_die(hmmp->nr_v * hmmp->nr_v * sizeof(double)));
    E = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size * sizeof(double)));
    T_lab = (double*)(malloc_or_die(hmmp->nr_v * hmmp->nr_v * sizeof(double)));
    E_lab = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size * sizeof(double)));
    T_ulab = (double*)(malloc_or_die(hmmp->nr_v * hmmp->nr_v * sizeof(double)));
    E_ulab = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size * sizeof(double))); 

    if(hmmp->nr_alphabets > 1) {
      E_2 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_2 * sizeof(double)));
      E_lab_2 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_2 * sizeof(double)));
      E_ulab_2 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_2 * sizeof(double)));
    }
    if(hmmp->nr_alphabets > 2) {
      E_3 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_3 * sizeof(double)));
      E_lab_3 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_3 * sizeof(double)));
      E_ulab_3 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_3 * sizeof(double)));
    }
    if(hmmp->nr_alphabets > 3) {
      E_4 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_4 * sizeof(double)));
      E_lab_4 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_4 * sizeof(double)));
      E_ulab_4 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_4 * sizeof(double)));
    }
    
    old_log_likelihood_ulab = new_log_likelihood_ulab;
    new_log_likelihood_ulab = 0.0;
    old_log_likelihood_lab = new_log_likelihood_lab;
    new_log_likelihood_lab = 0.0;
    
    /* reset sequence pointer */
    msa_seq_infop = msa_seq_infop_start;

    for(s = 0; s < nr_seqs; s++) {
      if(use_lead_columns == YES) {
	seq_len = msa_seq_infop->nr_lead_columns;
      }
      else {
	seq_len = msa_seq_infop->msa_seq_length;
      }

      /* calculate for unlabeled sequences */

      /* calculate forward and backward matrices
       * memory for forw_mtx, scale_mtx and
       * backw_mtx is allocated in the functions */
#ifdef DEBUG_BW2      
      printf("running forward unlabeled\n");
#endif
      msa_forward_multi(hmmp, msa_seq_infop, use_lead_columns, use_gap_shares, NO, &forw_mtx, &forw_scale, NO, normalize,
			scoring_method, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
#ifdef DEBUG_BW2     
      printf("running backward unlabeled\n");
#endif
      msa_backward_multi(hmmp, msa_seq_infop, use_lead_columns, use_gap_shares, &backw_mtx, forw_scale, NO, normalize,
			 scoring_method, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
#ifdef  DEBUG_BW2
      printf("done with backward unlabeled\n");
#endif
      /* update new_log_likelihood */
      likelihood = log10((forw_mtx +
			  get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob);
      for(k = 0; k <= seq_len; k++) {
	likelihood = likelihood + log10(*(forw_scale + k));
      }
#ifdef DEBUG_BW
      dump_scaling_array(k-1,forw_scale);
      printf("likelihood = %f\n", likelihood);
#endif      
      new_log_likelihood_ulab += likelihood;
      for(k = 0; k < hmmp->nr_v-1; k++) /* k = from vertex */ {
	lp = *(hmmp->to_trans_array + k);
	while(lp->vertex != END) /* l = to-vertex */ {
	  i = 0;
	  while(1) {
	    if(use_lead_columns == NO) {
	      p = i;
	    }
	    else {
	      p = *(msa_seq_infop->lead_columns_start + i);
	    }
	    add_Tkl_contribution_msa_multi(hmmp, msa_seq_infop, forw_mtx, backw_mtx,
					   forw_scale, p, k, lp, T_ulab, use_gap_shares, use_lead_columns, i, NO, scoring_method,
					   normalize, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
	    i++;
	    if(use_lead_columns == NO) {
	      if(i >= seq_len) {
		break;
	      }
	    }
	    else {
	      if(*(msa_seq_infop->lead_columns_start + i) == END) {
		break;
	      }
	    }
	  }
	  /* move on to next path */
	  while(lp->next != NULL) {
	    lp++;
	  }
	  lp++;
	}
	
	
	/* calculate E[k][a] contribution from this sequence */
	if(silent_state_multi(k, hmmp) != 0) {
	  i = 0;
	  while(1) {
	    
	    /* get correct incex for this letter column */
	    if(use_lead_columns == NO) {
	      p = i;
	    }
	    else {
	      p = *(msa_seq_infop->lead_columns_start + i);
	    }

	    /* get basic scoring result */
	      Eka_base = add_Eka_contribution_msa_multi(hmmp, msa_seq_infop, forw_mtx, backw_mtx, p, k,
							i, use_lead_columns);

	    /* loop over the alphabets */
	    for(alphabet_nr = 1; alphabet_nr <= hmmp->nr_alphabets; alphabet_nr++) {
	      if(alphabet_nr == 1) {
		alphabet = hmmp->alphabet;
		subst_mtx = hmmp->subst_mtx;
		a_size = hmmp->a_size;
		E_cur = E_ulab;
		msa_seq = msa_seq_infop->msa_seq_1;
	      }
	      else if(alphabet_nr == 2) {
		alphabet = hmmp->alphabet_2;
		subst_mtx = hmmp->subst_mtx_2;
		a_size = hmmp->a_size_2;
		E_cur = E_ulab_2;
		msa_seq = msa_seq_infop->msa_seq_2;
	      }
	      else if(alphabet_nr == 3) {
		alphabet = hmmp->alphabet_3;
		subst_mtx = hmmp->subst_mtx_3;
		a_size = hmmp->a_size_3;
		E_cur = E_ulab_3;
		msa_seq = msa_seq_infop->msa_seq_3;
	      }
	      else if(alphabet_nr == 4) {
		alphabet = hmmp->alphabet_4;
		subst_mtx = hmmp->subst_mtx_4;
		a_size = hmmp->a_size_4;
		E_cur = E_ulab_4;
		msa_seq = msa_seq_infop->msa_seq_4;
	      }
	      
	      /* get result and add to matrix according to scoring method */
	      add_to_E_multi(E_cur, Eka_base, msa_seq, p, k, a_size, normalize, subst_mtx,
			     alphabet, scoring_method, use_nr_occ, DISCRETE, NULL);
	    }
	    /* update loop index, check if we are done */
	    i++;
	    if(use_lead_columns == NO) {
	      if(i >= seq_len) {
		break;
	      }
	    }
	    else {
	      if(*(msa_seq_infop->lead_columns_start + i) == END) {
		break;
	      }
	    }
	  }
	}
      }
      
#ifdef DEBUG_BW
      dump_T_matrix(hmmp->nr_v, hmmp->nr_v, T_ulab);
      dump_E_matrix(hmmp->nr_v, hmmp->a_size, E_ulab);
#endif


      /* some garbage collection */
      free(forw_mtx);
      free(backw_mtx);
      free(forw_scale);
     
      
      /* calculate for labeled seqs */
      
      /* calculate forward and backward matrices
       * memory for forw_mtx, scale_mtx and
       * backw_mtx is allocated in the functions */
#ifdef DEBUG_BW2
      printf("running forward labeled\n");
#endif
      msa_forward_multi(hmmp, msa_seq_infop, use_lead_columns, use_gap_shares, NO, &forw_mtx,
			&forw_scale, YES, normalize, scoring_method, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3,
			aa_freqs_4);
#ifdef DEBUG_BW2
      printf("running backward labeled\n");
#endif
      msa_backward_multi(hmmp, msa_seq_infop, use_lead_columns, use_gap_shares, &backw_mtx, forw_scale,
			 YES, normalize, scoring_method, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
#ifdef DEBUG_BW2     
      printf("done with backward labeled\n");
#endif
      /* update new_log_likelihood */
      likelihood = log10((forw_mtx +
			  get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob);
      for(k = 0; k <= seq_len; k++) {
	likelihood = likelihood + log10(*(forw_scale + k));
      }
#ifdef DEBUG_BW
      dump_scaling_array(k-1,forw_scale);
      printf("likelihood = %f\n", likelihood);
#endif      
      new_log_likelihood_lab += likelihood;
      
      for(k = 0; k < hmmp->nr_v-1; k++) /* k = from vertex */ {
	lp = *(hmmp->to_trans_array + k);
	while(lp->vertex != END) /* l = to-vertex */ {
	  i = 0;
	  while(1) {
	    if(use_lead_columns == NO) {
	      p = i;
	    }
	    else {
	      p = *(msa_seq_infop->lead_columns_start + i);
	    }
	    /* add T[k][l] contribution for the msa-sequence */
	    add_Tkl_contribution_msa_multi(hmmp, msa_seq_infop, forw_mtx, backw_mtx,
					   forw_scale, p, k, lp, T_lab, use_gap_shares, use_lead_columns, i, YES, scoring_method,
					   normalize, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
	    i++;
	    if(use_lead_columns == NO) {
	      if(i >= seq_len) {
		break;
	      }
	    }
	    else {
	      if(*(msa_seq_infop->lead_columns_start + i) == END) {
		break;
	      }
	    }
	  }
	  /* move on to next path */
	  while(lp->next != NULL) {
	    lp++;
	  }
	  lp++;
	}
	
	/* calculate E[k][a] contribution from this sequence */
	if(silent_state_multi(k, hmmp) != 0) {
	  i = 0;
	  while(1) {
	    
	    /* get correct incex for this letter column */
	    if(use_lead_columns == NO) {
	      p = i;
	    }
	    else {
	      p = *(msa_seq_infop->lead_columns_start + i);
	    }

	    /* get basic scoring result */
	      Eka_base = add_Eka_contribution_msa_multi(hmmp, msa_seq_infop, forw_mtx, backw_mtx, p, k,
							i, use_lead_columns);

	    /* loop over the alphabets */
	    for(alphabet_nr = 1; alphabet_nr <= hmmp->nr_alphabets; alphabet_nr++) {
	      if(alphabet_nr == 1) {
		alphabet = hmmp->alphabet;
		subst_mtx = hmmp->subst_mtx;
		a_size = hmmp->a_size;
		E_cur = E_lab;
		msa_seq = msa_seq_infop->msa_seq_1;
	      }
	      else if(alphabet_nr == 2) {
		alphabet = hmmp->alphabet_2;
		subst_mtx = hmmp->subst_mtx_2;
		a_size = hmmp->a_size_2;
		E_cur = E_lab_2;
		msa_seq = msa_seq_infop->msa_seq_2;
	      }
	      else if(alphabet_nr == 3) {
		alphabet = hmmp->alphabet_3;
		subst_mtx = hmmp->subst_mtx_3;
		a_size = hmmp->a_size_3;
		E_cur = E_lab_3;
		msa_seq = msa_seq_infop->msa_seq_3;
	      }
	      else if(alphabet_nr == 4) {
		alphabet = hmmp->alphabet_4;
		subst_mtx = hmmp->subst_mtx_4;
		a_size = hmmp->a_size_4;
		E_cur = E_lab_4;
		msa_seq = msa_seq_infop->msa_seq_4;
	      }
	      
	      /* get result and add to matrix according to scoring method */
	      add_to_E_multi(E_cur, Eka_base, msa_seq, p, k, a_size, normalize, subst_mtx,
			     alphabet, scoring_method, use_nr_occ, DISCRETE, NULL);
	    }
	    /* update loop index, check if we are done */
	    i++;
	    if(use_lead_columns == NO) {
	      if(i >= seq_len) {
		break;
	      }
	    }
	    else {
	      if(*(msa_seq_infop->lead_columns_start + i) == END) {
		break;
	      }
	    }
	  }
	}
      }
      
#ifdef DEBUG_BW
      dump_T_matrix(hmmp->nr_v, hmmp->nr_v, T_lab);
      dump_E_matrix(hmmp->nr_v, hmmp->a_size, E_lab);
#endif

      /* some garbage collection */
      free(forw_mtx);
      free(backw_mtx);
      free(forw_scale);
   
      msa_seq_infop++;
    }
    
    if(verbose == YES) {
      printf("log likelihood diff rd %d: %f\n", iteration, new_log_likelihood_ulab - new_log_likelihood_lab);
    }
    
    /* recalculate emission expectations according to distribution groups 
     * by simply taking the mean of the expected emissions within this group
     * for each letter in the alphabet and replacing each expectation for the
     * letter with this value for every member of the distribution group */
    recalculate_emiss_expectations_multi(hmmp, E_lab, 1);
    if(hmmp->nr_alphabets > 1) {
      recalculate_emiss_expectations_multi(hmmp, E_lab_2, 2);
    }
    if(hmmp->nr_alphabets > 2) {
      recalculate_emiss_expectations_multi(hmmp, E_lab_3, 3);
    }
    if(hmmp->nr_alphabets > 3) {
      recalculate_emiss_expectations_multi(hmmp, E_lab_4, 4);
    }
    
    recalculate_emiss_expectations_multi(hmmp, E_ulab, 1);
    if(hmmp->nr_alphabets > 1) {
      recalculate_emiss_expectations_multi(hmmp, E_ulab_2, 2);
    }
    if(hmmp->nr_alphabets > 2) {
      recalculate_emiss_expectations_multi(hmmp, E_ulab_3, 3);
    }
    if(hmmp->nr_alphabets > 3) {
      recalculate_emiss_expectations_multi(hmmp, E_ulab_4, 4);
    }
    
    /* recalculate transition expectations for tied transitions according
     * to the same scheme as for emission distribution groups */
    recalculate_trans_expectations_multi(hmmp, T_lab);
    recalculate_trans_expectations_multi(hmmp, T_ulab);
    
#ifdef DEBUG_EXTBW
    dump_T_matrix(hmmp->nr_v, hmmp->nr_v, T_lab);
    dump_T_matrix(hmmp->nr_v, hmmp->nr_v, T_ulab);
    dump_E_matrix(hmmp->nr_v, hmmp->a_size, E_lab);
#endif

    /* update real T end E matrices */
    calculate_TE_contributions_multi(T, E, E_2, E_3, E_4, T_lab, E_lab, E_lab_2, E_lab_3, E_lab_4, T_ulab, E_ulab,
				     E_ulab_2, E_ulab_3, E_ulab_4, hmmp->emissions, hmmp->emissions_2, hmmp->emissions_3,
				     hmmp->emissions_4, hmmp->transitions, hmmp->nr_v, hmmp->a_size,
				     hmmp->a_size_2, hmmp->a_size_3, hmmp->a_size_4, hmmp->vertex_emiss_prior_scalers,
				     hmmp->vertex_emiss_prior_scalers_2, hmmp->vertex_emiss_prior_scalers_3,
				     hmmp->vertex_emiss_prior_scalers_4, iteration, hmmp->nr_alphabets);
    
    /* check if likelihood change is small enough, then we are done */
    if(fabs((new_log_likelihood_ulab - new_log_likelihood_lab) - (old_log_likelihood_ulab - old_log_likelihood_lab)) <
       CML_THRESHOLD && annealing_status == DONE) {
      break;
    }
    
    /* if simulated annealing is used, scramble results in E and T matrices */
    if(annealing == YES && temperature > ANNEAL_THRESHOLD) {
      anneal_E_matrix_multi(temperature, E, hmmp, 1);
      if(hmmp->nr_alphabets > 1) {
	anneal_E_matrix_multi(temperature, E_2, hmmp, 2);
      }
      if(hmmp->nr_alphabets > 2) {
	anneal_E_matrix_multi(temperature, E_3, hmmp, 3);
      }
      if(hmmp->nr_alphabets > 3) {
	anneal_E_matrix_multi(temperature, E_4, hmmp, 4);
      }
      anneal_T_matrix_multi(temperature, T, hmmp);
      temperature = temperature * cooling_factor;
    }

    if(temperature < ANNEAL_THRESHOLD) {
      annealing_status = DONE;
    }
    
    for(k = 0; k < hmmp->nr_v-1; k++) /* k = from-vertex */ {
      /* update transition matrix */
      if(use_transition_pseudo_counts == YES) {
	update_trans_mtx_pseudocount_multi(hmmp, T, k);
      }
      else {
	update_trans_mtx_std_multi(hmmp, T, k);
      }
      
#ifdef DEBUG_PRIORS
      printf("Starting emission matrix update\n");
#endif
      
      /* update emission matrix using Dirichlet prior files if they exist*/
      priorp = *(hmmp->ed_ps + k);
      if(priorp != NULL && use_prior == YES) {
#ifdef DEBUG_PRIORS	
	printf("k = %d\n", k);
	printf("value = %x\n", priorp);
#endif
	update_emiss_mtx_prior_multi(hmmp, E, k, priorp, 1);
      }
      else if(use_emission_pseudo_counts == YES) /* update emissions matrix "normally" when dirichlet file is missing */ {
	update_emiss_mtx_pseudocount_multi(hmmp, E, k, 1);
      }
      else {
	update_emiss_mtx_std_multi(hmmp, E, k, 1);
      }
      

      if(hmmp->nr_alphabets > 1) {
	priorp = *(hmmp->ed_ps_2 + k);
	if(priorp != NULL && use_prior == YES) {
	  update_emiss_mtx_prior_multi(hmmp, E_2, k, priorp, 2);
	}
	else if(use_emission_pseudo_counts == YES) /* update emissions matrix "normally" when dirichlet file is missing */ {
	  update_emiss_mtx_pseudocount_multi(hmmp, E_2, k, 2);
	}
	else {
	  update_emiss_mtx_std_multi(hmmp, E_2, k, 2);
	}
      }

      if(hmmp->nr_alphabets > 2) {
	priorp = *(hmmp->ed_ps_3 + k);
	if(priorp != NULL && use_prior == YES) {
	  update_emiss_mtx_prior_multi(hmmp, E_3, k, priorp, 3);
	}
	else if(use_emission_pseudo_counts == YES) /* update emissions matrix "normally" when dirichlet file is missing */ {
	  update_emiss_mtx_pseudocount_multi(hmmp, E_3, k, 3);
	}
	else {
	  update_emiss_mtx_std_multi(hmmp, E_3, k, 3);
	}
      }

      if(hmmp->nr_alphabets > 3) {
	priorp = *(hmmp->ed_ps_4 + k);
	if(priorp != NULL && use_prior == YES) {
	  update_emiss_mtx_prior_multi(hmmp, E_4, k, priorp, 4);
	}
	else if(use_emission_pseudo_counts == YES) /* update emissions matrix "normally" when dirichlet file is missing */ {
	  update_emiss_mtx_pseudocount_multi(hmmp, E_4, k, 4);
	}
	else {
	  update_emiss_mtx_std_multi(hmmp, E_4, k, 4);
	}
      }
    }
    
#ifdef DEBUG_BW
    dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
    dump_emiss_matrix(hmmp->nr_v, hmmp->a_size, hmmp->emissions);
#endif    
    
    /* some garbage collection */
    free(E);
    if(hmmp->nr_alphabets > 1) {
      free(E_2);
    }
    if(hmmp->nr_alphabets > 2) {
      free(E_3);
    }
    if(hmmp->nr_alphabets > 3) {
      free(E_4);
    }
    free(T);
    free(T_lab);
    free(E_lab);
    if(hmmp->nr_alphabets > 1) {
      free(E_lab_2);
    }
    if(hmmp->nr_alphabets > 2) {
      free(E_lab_3);
    }
    if(hmmp->nr_alphabets > 3) {
      free(E_lab_4);
    }
    free(T_ulab);
    free(E_ulab);
    if(hmmp->nr_alphabets > 1) {
      free(E_ulab_2);
    }
    if(hmmp->nr_alphabets > 2) {
      free(E_ulab_3);
    }
    if(hmmp->nr_alphabets > 3) {
      free(E_ulab_4);
    }
    
    max_nr_iterations--;
    iteration++;
#ifdef DEBUG_BW2    
    printf("end of baum-welch-loop\n");
#endif
  }
  while(max_nr_iterations > 0); /* break condition is also when log_likelihood_difference is
				 * smaller than THRESHOLD, checked inside the loop for
				 * better efficiency */
  
  
#ifdef DEBUG_BW
  dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
  dump_emiss_matrix(hmmp->nr_v, hmmp->a_size, hmmp->emissions);
#endif    
}




/*************************************************************/
/************** help functions *******************************/
/*************************************************************/

/* calculates T[k][l] contribution */
void add_Tkl_contribution_multi(struct hmm_multi_s *hmmp, struct letter_s *seq, struct letter_s *seq_2,
				struct letter_s *seq_3, struct letter_s *seq_4, struct forward_s *forw_mtx,
				struct backward_s *backw_mtx, double *forw_scale, int p, int k,
				struct path_element *lp_const, int a_index, int a_index_2, int a_index_3, int a_index_4,
				double *T, int use_labels, int multi_scoring_method)
{
  
  double t_res, t_res_1, t_res_2, t_res_3, t_res_3_1, t_res_3_2, t_res_3_4, t_res_4, t_res_5, t_res_6, t_res_3_temp;
  struct path_element *lp_shadow, *lp, *lp_end;
  int from_v_end;
  int j;
  
  lp = lp_const;
  lp_shadow = lp;
  
  /* calculate T[k][l] contribution using scaled values*/
  t_res_1 = (forw_mtx + get_mtx_index(p-1, k, hmmp->nr_v))->prob;
  t_res_2 = *(hmmp->transitions + get_mtx_index(k, lp->vertex, hmmp->nr_v));
  while(lp->next != NULL) {
    t_res_2 = t_res_2 * *(hmmp->transitions + get_mtx_index(lp->vertex, (lp+1)->vertex, hmmp->nr_v));
    lp++;
  }
  if(use_labels == YES && *(hmmp->vertex_labels + lp->vertex) != seq[p-1].label && seq[p-1].label != '.') {
    t_res_3 = 0.0;
  }
  else {
    if(multi_scoring_method == JOINT_PROB) {
      if(hmmp->alphabet_type == DISCRETE) {
	if(a_index < 0) {
	  /* letter is not in alphabet => it is a replacement letter.
	   * This is not implemented yet, simple solution is to ignore this value, i.e. set letter to 'X' */
	  t_res_3 = 0.0;
	}
	else {
	  t_res_3 = *(hmmp->emissions + get_mtx_index(lp->vertex, a_index, hmmp->a_size));
	}
      }
      else {
	t_res_3 = 0.0;
	for(j = 0; j < hmmp->a_size / 3; j++) {
	  t_res_3 += get_single_gaussian_statescore(*(hmmp->emissions + get_mtx_index(lp->vertex, (j * 3), hmmp->a_size)),
						    *(hmmp->emissions + get_mtx_index(lp->vertex, (j * 3 + 1), hmmp->a_size)),
						    seq[p-1].cont_letter) *
	    *((hmmp->emissions) + (lp->vertex * (hmmp->a_size)) + (j * 3 + 2));
	}
      }
      if(hmmp->nr_alphabets > 1) {
	if(hmmp->alphabet_type_2 == DISCRETE) {
	  if(a_index_2 < 0) {
	    /* letter is not in alphabet => it is a replacement letter.
	     * This is not implemented yet, simple solution is to ignore this value, i.e. set letter to 'X' */
	    t_res_3 = 0.0;
	  }
	  else {
	    t_res_3 *= *(hmmp->emissions_2 + get_mtx_index(lp->vertex, a_index_2, hmmp->a_size_2));
	  }	
	}
	else {
	  t_res_3_temp = 0.0;
	  for(j = 0; j < hmmp->a_size_2 / 3; j++) {
	    t_res_3_temp += get_single_gaussian_statescore(*(hmmp->emissions_2 + get_mtx_index(lp->vertex, (j * 3), hmmp->a_size_2)),
						      *(hmmp->emissions_2 + get_mtx_index(lp->vertex, (j * 3 + 1), hmmp->a_size_2)),
						      seq_2[p-1].cont_letter) *
	      *((hmmp->emissions_2) + (lp->vertex * (hmmp->a_size_2)) + (j * 3 + 2));
	  }
	  t_res_3 *= t_res_3_temp;
	}
      }
      if(hmmp->nr_alphabets > 2) {
	if(hmmp->alphabet_type_2 == DISCRETE) {
	  if(a_index_3 < 0) {
	    /* letter is not in alphabet => it is a replacement letter.
	     * This is not implemented yet, simple solution is to ignore this value, i.e. set letter to 'X' */
	    t_res_3 = 0.0;
	  }
	  else {
	    t_res_3 *= *(hmmp->emissions_3 + get_mtx_index(lp->vertex, a_index_3, hmmp->a_size_3));
	  }
	}
	else {
	  t_res_3_temp = 0.0;
	  for(j = 0; j < hmmp->a_size_3 / 3; j++) {
	    t_res_3_temp += get_single_gaussian_statescore(*(hmmp->emissions_3 + get_mtx_index(lp->vertex, (j * 3), hmmp->a_size_3)),
						      *(hmmp->emissions_3 + get_mtx_index(lp->vertex, (j * 3 + 1), hmmp->a_size_3)),
						      seq_3[p-1].cont_letter) *
	      *((hmmp->emissions_3) + (lp->vertex * (hmmp->a_size_3)) + (j * 3 + 2));
	  }
	  t_res_3 *= t_res_3_temp;
	}
      }
      if(hmmp->nr_alphabets > 3) {
	if(hmmp->alphabet_type_3 == DISCRETE) {
	  if(a_index_4 < 0) {
	    /* letter is not in alphabet => it is a replacement letter.
	     * This is not implemented yet, simple solution is to ignore this value, i.e. set letter to 'X' */
	    t_res_3 = 0.0;
	  }
	  else {
	    t_res_3 *= *(hmmp->emissions_4 + get_mtx_index(lp->vertex, a_index_4, hmmp->a_size_4));
	  }
	}
	else {
	  t_res_3_temp = 0.0;
	  for(j = 0; j < hmmp->a_size_4 / 3; j++) {
	    t_res_3_temp += get_single_gaussian_statescore(*(hmmp->emissions_4 + get_mtx_index(lp->vertex, (j * 3), hmmp->a_size_4)),
						      *(hmmp->emissions_4 + get_mtx_index(lp->vertex, (j * 3 + 1), hmmp->a_size_4)),
						      seq_4[p-1].cont_letter) *
	      *((hmmp->emissions_4) + (lp->vertex * (hmmp->a_size_4)) + (j * 3 + 2));
	  }
	  t_res_3 *= t_res_3_temp;
	}
      }
    }
  }
  t_res_4 =  (backw_mtx + get_mtx_index(p, lp->vertex, hmmp->nr_v))->prob;
  t_res = t_res_1 * t_res_2 * t_res_3 * t_res_4;
  if(t_res == 0) {
    return ; /* no reason to update with zero value */
  }

#ifdef DEBUG_BW_TRANS
  printf("for T[%d][%d]\n", k, lp->vertex);
  sequence_as_string(seq);
  printf("forw_mtx value = %f\n", t_res_1);
  printf("transitions value = %f\n", t_res_2);
  printf("emissions value = %f\n", t_res_3);
  printf("backw_mtx value = %f\n", t_res_4);
#endif  
  
  /* divide by scaled probability for sequence s */
  t_res_5 = (forw_mtx + get_mtx_index(get_seq_length(seq)+1, hmmp->nr_v-1, hmmp->nr_v))->prob;
#ifdef DEBUG_BW_TRANS 
  printf("seq_length = %d\n", get_seq_length(seq));
  printf("t_res_5 = %f\n", t_res_5);
#endif
  
  t_res = t_res / t_res_5;
  
  /* divide result with scaling factor for position p,
   * since it is not included in the contribution, but is
   * included in the probability for the sequence */
  t_res_6 = *(forw_scale + p);
  t_res = t_res / t_res_6;
  
#ifdef DEBUG_BW_TRANS
  printf("tot prob for sequence = %f\n", t_res_5);
  printf("scale p = %f\n", t_res_6);
  printf("res = %f\n", t_res);
#endif
  
  /* if last letter of sequence and a path from the emitting state for this letter to end state exists
     update this path as well with the same value */
  if(p == get_seq_length(seq) && path_length_multi(lp->vertex, hmmp->nr_v-1, hmmp, 0) > 0) {
    from_v_end = lp->vertex;
    lp_end = get_end_path_start_multi(lp->vertex, hmmp);
    *(T + get_mtx_index(from_v_end, lp_end->vertex, hmmp->nr_v)) += t_res;
    from_v_end = lp_end->vertex;
    lp_end++;
    while(lp_end->next != NULL) {
      *(T + get_mtx_index(from_v_end, lp_end->vertex, hmmp->nr_v)) += t_res;
      from_v_end = lp_end->vertex;
      lp_end++;
    }
  }

  /* update T-matrix for transition indices that correspond to current path */
  lp = lp_shadow; 
  *(T + get_mtx_index(k, lp->vertex, hmmp->nr_v)) += t_res;
  lp++;
  while(lp_shadow->next != NULL) {
    *(T + get_mtx_index(lp_shadow->vertex, lp->vertex, hmmp->nr_v)) += t_res;
    lp_shadow = lp;
    lp++;   
  }
 
#ifdef DEBUG_BW_TRANS
  printf("adding result to T-mtx index: %d\n", get_mtx_index(k,lp_shadow->vertex,hmmp->nr_v));
#endif
}

/* calculates T[k][l] contribution */
void add_Tkl_contribution_msa_multi(struct hmm_multi_s *hmmp, struct msa_sequences_multi_s *msa_seq_infop,
				    struct forward_s *forw_mtx, struct backward_s *backw_mtx,
				    double *forw_scale, int p, int k,
				    struct path_element *lp_const, double *T, int use_gap_shares,
				    int use_lead_columns, int i, int use_labels, int scoring_method, int normalize,
				    int multi_scoring_method, double *aa_freqs_1, double *aa_freqs_2,
				    double *aa_freqs_3, double *aa_freqs_4)
{
  double t_res, t_res_1, t_res_2, t_res_3, t_res_4, t_res_5, t_res_6, temp_res, t_res_3_tot;
  struct path_element *lp_shadow, *lp, *lp_end;
  int a_index, a_index2;
  int query_index;
  double default_share, rest_share;
  double seq_normalizer;
  double state_normalizer;
  double subst_mtx_normalizer;
  int from_v_end;
  
  int alphabet;

  int a_size, a_size_1;
  struct msa_letter_s *msa_seq;
  double *emissions;
  double *subst_mtx;
  double *aa_freqs;
  int alphabet_type;
  int j;

  lp = lp_const;
  lp_shadow = lp;
  
  
  /* calculate T[k][l] contribution using scaled values*/
  
  /* get f_i value */
  t_res_1 = (forw_mtx + get_mtx_index(i, k, hmmp->nr_v))->prob;

  /* get a_kl value*/
  t_res_2 = *(hmmp->transitions + get_mtx_index(k, lp->vertex, hmmp->nr_v));
  while(lp->next != NULL) {
    t_res_2 = t_res_2 * *(hmmp->transitions + get_mtx_index(lp->vertex, (lp+1)->vertex, hmmp->nr_v));
    lp++;
  }

  /* calculate e_i(x) */
  t_res_3_tot = 1.0;
  a_size_1 = hmmp->a_size;
  for(alphabet = 1; alphabet <= hmmp->nr_alphabets; alphabet++) {
    seq_normalizer = 0.0;
    state_normalizer = 0.0;
    subst_mtx_normalizer = 0.0;
    if(alphabet == 1) {
      if(hmmp->alphabet_type == DISCRETE) {
	query_index = get_alphabet_index_msa_query((msa_seq_infop->msa_seq_1 + (p * (hmmp->a_size+1)))->query_letter,
						   hmmp->alphabet, hmmp->a_size);
	if(query_index < 0) {
	  query_index = hmmp->a_size; /* if letter is wild card, use default column in subst matrix */
	}
      }
      a_size = hmmp->a_size;
      msa_seq = msa_seq_infop->msa_seq_1;
      emissions = hmmp->emissions;
      subst_mtx = hmmp->subst_mtx;
      alphabet_type = hmmp->alphabet_type;
      aa_freqs = aa_freqs_1;
    }
    if(alphabet == 2) {
      if(hmmp->alphabet_type_2 == DISCRETE) {
	query_index = get_alphabet_index_msa_query((msa_seq_infop->msa_seq_2 + (p * (hmmp->a_size_2+1)))->query_letter,
						   hmmp->alphabet_2, hmmp->a_size_2);
	if(query_index < 0) {
	  query_index = hmmp->a_size_2; /* if letter is wild card, use default column in subst matrix */
	}
      }
      a_size = hmmp->a_size_2;
      msa_seq = msa_seq_infop->msa_seq_2;
      emissions = hmmp->emissions_2;
      subst_mtx = hmmp->subst_mtx_2;
      alphabet_type = hmmp->alphabet_type_2;
      aa_freqs = aa_freqs_2;
    }
    if(alphabet == 3) {
      if(hmmp->alphabet_type_3 == DISCRETE) {
	query_index = get_alphabet_index_msa_query((msa_seq_infop->msa_seq_3 + (p * (hmmp->a_size_3+1)))->query_letter,
						   hmmp->alphabet_3, hmmp->a_size_3);
	if(query_index < 0) {
	  query_index = hmmp->a_size_3; /* if letter is wild card, use default column in subst matrix */
	}
      }
      a_size = hmmp->a_size_3;
      msa_seq = msa_seq_infop->msa_seq_3;
      emissions = hmmp->emissions_3;
      subst_mtx = hmmp->subst_mtx_3;
      alphabet_type = hmmp->alphabet_type_3;
      aa_freqs = aa_freqs_3;
    }
    if(alphabet == 4) {
      if(hmmp->alphabet_type_4 == DISCRETE) {
	query_index = get_alphabet_index_msa_query((msa_seq_infop->msa_seq_4 + (p * (hmmp->a_size_4+1)))->query_letter,
						   hmmp->alphabet_4, hmmp->a_size_4);
	if(query_index < 0) {
	  query_index = hmmp->a_size_4; /* if letter is wild card, use default column in subst matrix */
	}
      }
      a_size = hmmp->a_size_4;
      msa_seq = msa_seq_infop->msa_seq_4;
      emissions = hmmp->emissions_4;
      subst_mtx = hmmp->subst_mtx_4;
      alphabet_type = hmmp->alphabet_type_4;
      aa_freqs = aa_freqs_4;
    }
    
    default_share = 1.0 / (double)(a_size);
    t_res_3 = 0.0;

    /* use first alphabet here since the labels are placed in the first alphabet */
    if(use_labels == YES && *(hmmp->vertex_labels + lp->vertex) !=
       (msa_seq_infop->msa_seq_1 + get_mtx_index(p,0, a_size_1+1))->label &&
       (msa_seq_infop->msa_seq_1 + get_mtx_index(p,0, a_size_1+1))->label != '.') {
      t_res_3 = 0.0;
    }
    else if(alphabet_type == CONTINUOUS) {
      if((msa_seq + get_mtx_index(p, 0, a_size + 1))->nr_occurences > 0.0) {
	t_res_3 = 0.0;
	for(j = 0; j < a_size / 3; j++) {
	  t_res_3 += get_single_gaussian_statescore(*(emissions + get_mtx_index(lp->vertex, (j * 3), a_size)),
						    *(emissions + get_mtx_index(lp->vertex, (j * 3 + 1), a_size)),
						    (msa_seq + get_mtx_index(p, 0, a_size + 1))->share) *
	    *((emissions) + (lp->vertex * (a_size)) + (j * 3 + 2));
	}
      }
      else {
	t_res_3 = 0.0;
      }
    }
    else if(scoring_method == DOT_PRODUCT) {
      t_res_3 = get_dp_statescore(a_size, use_gap_shares, NO, msa_seq, p, emissions,
				  lp->vertex, normalize, msa_seq_infop->gap_shares);
    }
    else if(scoring_method == DOT_PRODUCT_PICASSO) {
      t_res_3 = get_dp_picasso_statescore(a_size, use_gap_shares, NO, msa_seq, p, emissions,
					 lp->vertex ,normalize, msa_seq_infop->gap_shares, aa_freqs);
    }
    else if(scoring_method == PICASSO) {
      t_res_3 = get_picasso_statescore(a_size, use_gap_shares, NO, msa_seq, p, emissions,
				       lp->vertex ,normalize, msa_seq_infop->gap_shares, aa_freqs);
    }
    else if(scoring_method == PICASSO_SYM) {
      t_res_3 = get_picasso_sym_statescore(a_size, use_gap_shares, NO, msa_seq, p, emissions,
					   lp->vertex ,normalize, msa_seq_infop->gap_shares, aa_freqs);
    }
    else if(scoring_method == SJOLANDER) {
      t_res_3 = get_sjolander_statescore(a_size, use_gap_shares, NO, msa_seq, p, emissions,
					 lp->vertex, normalize, msa_seq_infop->gap_shares);
    }
    else if(scoring_method == SJOLANDER_REVERSED) {
      t_res_3 = get_sjolander_statescore(a_size, use_gap_shares, NO, msa_seq, p, emissions,
					 lp->vertex, normalize, msa_seq_infop->gap_shares);
    }
    else if(scoring_method == SUBST_MTX_PRODUCT) {
      t_res_3 = get_subst_mtx_product_statescore(a_size, use_gap_shares, NO, msa_seq, p, emissions,
						 lp->vertex, subst_mtx);
    }
    else if(scoring_method == SUBST_MTX_DOT_PRODUCT) {
      t_res_3 = get_subst_mtx_dot_product_statescore(a_size, use_gap_shares, NO, msa_seq, p, emissions,
						     lp->vertex, normalize, msa_seq_infop->gap_shares,
						     query_index, subst_mtx);
    }
    else if(scoring_method == SUBST_MTX_DOT_PRODUCT_PRIOR) {
      t_res_3 = get_subst_mtx_dot_product_prior_statescore(a_size, use_gap_shares, NO, msa_seq, p,
							   emissions, lp->vertex, normalize, msa_seq_infop->gap_shares,
							   query_index, subst_mtx);
    }
    
    
    if(multi_scoring_method == JOINT_PROB) {
      t_res_3_tot *= t_res_3;
    }
    else {
      printf("Error: only joint prob is multiscoring is implemented\n");
    }
  }
  
  
  /* get b_i+1 value */
  t_res_4 = (backw_mtx + get_mtx_index(i+1, lp->vertex, hmmp->nr_v))->prob;
  
  t_res = t_res_1 * t_res_2 * t_res_3_tot * t_res_4;
  
#ifdef DEBUG_Tkl
  printf("for T[%d][%d]\n", k, lp->vertex);
  printf("forw_mtx value = %f\n", t_res_1);
  printf("transitions value = %f\n", t_res_2);
  printf("emissions value = %f\n", t_res_3_tot);
  printf("backw_mtx value = %f\n", t_res_4);
  printf("t_res = %f\n", t_res);
#endif  
  if(t_res == 0) {
    return ; /* no reason to update with zero value */
  } 
   
  /* divide by scaled probability for sequence s */
  if(use_lead_columns == NO) {
    t_res_5 = (forw_mtx + get_mtx_index((msa_seq_infop->msa_seq_length)+1,
					hmmp->nr_v-1, hmmp->nr_v))->prob;
  }
  else {
    t_res_5 = (forw_mtx + get_mtx_index((msa_seq_infop->nr_lead_columns)+1,
					hmmp->nr_v-1, hmmp->nr_v))->prob;
  }
  t_res = t_res / t_res_5;
  
  /* divide result with scaling factor for position p+1,
   * since it is not included in the contribution, but is
   * included in the probability for the sequence */
  t_res_6 = *(forw_scale + i+1);
  t_res = t_res / t_res_6;
  
#ifdef DEBUG_Tkl
  printf("tot prob for sequence = %f\n", t_res_5);
  printf("scale p = %f\n", t_res_6);
  printf("res = %f\n", t_res);
  printf("p = %d\n", p);
  printf("seq_length = %d\n", msa_seq_infop->msa_seq_length);
#endif
 

  /* if last letter of sequence and a path from the emitting state for this letter to end state exists
     update this path as well with the same value */
  if(p == msa_seq_infop->msa_seq_length - 1 && path_length_multi(lp->vertex, hmmp->nr_v-1, hmmp, 0) > 0) {
    from_v_end = lp->vertex;
    lp_end = get_end_path_start_multi(lp->vertex, hmmp);
    *(T + get_mtx_index(from_v_end, lp_end->vertex, hmmp->nr_v)) += t_res;
    from_v_end = lp_end->vertex;
    lp_end++;
    while(lp_end->next != NULL) {
      *(T + get_mtx_index(from_v_end, lp_end->vertex, hmmp->nr_v)) += t_res;
      from_v_end = lp_end->vertex;
      lp_end++;
    }
  }
  
 

  /* update T-matrix for transition indices that corresponds to current path */
  lp = lp_shadow; 
  *(T + get_mtx_index(k, lp->vertex, hmmp->nr_v)) += t_res;
  lp++;
  while(lp_shadow->next != NULL) {
    *(T + get_mtx_index(lp_shadow->vertex, lp->vertex, hmmp->nr_v)) += t_res;
    lp_shadow = lp;
    lp++;   
  }
 
#ifdef DEBUG_Tkl
  printf("adding result to T-mtx index: %d\n\n", get_mtx_index(k,lp_shadow->vertex,hmmp->nr_v));
#endif
}


double add_Eka_contribution_multi(struct hmm_multi_s *hmmp, struct letter_s *seq, struct forward_s *forw_mtx,
				  struct backward_s *backw_mtx, int p, int k, int multi_scoring_method)
{
  double e_res, e_res_1, e_res_2, e_res_3;
  
  /* get contribution from this position in the sequence */
  e_res_1 = (forw_mtx + get_mtx_index(p, k, hmmp->nr_v))->prob;
  e_res_2 = (backw_mtx + get_mtx_index(p, k, hmmp->nr_v))->prob;
  e_res = e_res_1 * e_res_2;
  if(e_res == 0) {
    return 0.0; /* no use updating with a zero value */
  }

#ifdef DEBUG_BW
  printf("forw_mtx = %f\n", e_res_1);
  printf("backw_mtx = %f\n", e_res_2);
  printf("res = %f\n", e_res);
#endif	    
  
  /* divide with total probability of current sequence */
  e_res_3 = (forw_mtx + get_mtx_index(get_seq_length(seq)+1, hmmp->nr_v-1, hmmp->nr_v))->prob;
  e_res = e_res / e_res_3;
  
#ifdef DEBUG_BW
  printf("total prob = %f\n", e_res_3);
  printf("e_res = %f\n", e_res);
#endif
  
  if(multi_scoring_method == JOINT_PROB) {
    /* nothing more needs to be done, this updating procedure will update with the same number for all alphabets */
  }

  return e_res;
}

double add_Eka_contribution_continuous_multi(struct hmm_multi_s *hmmp, struct letter_s *seq, struct forward_s *forw_mtx,
					     struct backward_s *backw_mtx, int p, int k, int multi_scoring_method, double *E,
					     int alphabet)
{
  double e_res, e_res_1, e_res_2, e_res_3;
  double mean_value, varians;
  double continuous_score_all, continuous_score_j, gamma_p_j;
  int j;
  double *emissions;
  int a_size;
  

  if(alphabet == 1) {
    emissions = hmmp->emissions;
    a_size = hmmp->a_size;
  }
  else if(alphabet == 2) {
    emissions = hmmp->emissions_2;
    a_size = hmmp->a_size_2;
  }
  else if(alphabet == 3) {
    emissions = hmmp->emissions_3;
    a_size = hmmp->a_size_3;
  }
  else if(alphabet == 4) {
    emissions = hmmp->emissions_4;
    a_size = hmmp->a_size_4;
  }
  else {
    printf("strange alphabet nr: %d\n", alphabet);
    exit(0);
  }

  /* get contribution from this position in the sequence */
  e_res_1 = (forw_mtx + get_mtx_index(p, k, hmmp->nr_v))->prob;
  e_res_2 = (backw_mtx + get_mtx_index(p, k, hmmp->nr_v))->prob;
  e_res = e_res_1 * e_res_2;
  if(e_res == 0) {
    return 0.0; /* no use updating with a zero value */
  }

#ifdef DEBUG_BW
  printf("forw_mtx = %f\n", e_res_1);
  printf("backw_mtx = %f\n", e_res_2);
  printf("res = %f\n", e_res);
#endif	    
  

/* divide with total probability of current sequence */
  e_res_3 = (forw_mtx + get_mtx_index(get_seq_length(seq)+1, hmmp->nr_v-1, hmmp->nr_v))->prob;
  e_res = e_res / e_res_3;
  
  mean_value = (seq + p - 1)->cont_letter;
 
  
#ifdef DEBUG_BW
  printf("mean = %f\n", mean_value);
  printf("varians = %f\n", varians);
  printf("e_res = %f\n", e_res);
  printf("e_res_3 = %f\n", e_res_3);
#endif
  
  if(multi_scoring_method == JOINT_PROB) {
    /* nothing more needs to be done, this updating procedure will update with the same number for all alphabets */
  }
  
  continuous_score_all = 0.0;
  for(j = 0; j < a_size / 3; j++) {
    continuous_score_all += 
      get_single_gaussian_statescore(*(emissions + get_mtx_index(k, (j * 3), a_size)),
				     *(emissions + get_mtx_index(k, (j * 3 + 1), a_size)),
				     mean_value) *
      *((emissions) + (k * (a_size)) + (j * 3 + 2));
  }
  
  for(j = 0; j < a_size / 3; j++) {
    continuous_score_j = 
      get_single_gaussian_statescore(*(emissions + get_mtx_index(k, (j * 3), a_size)),
				     *(emissions + get_mtx_index(k, (j * 3 + 1), a_size)),
				     mean_value) *
      *((emissions) + (k * (a_size)) + (j * 3 + 2));
    varians = pow((seq + p - 1)->cont_letter - *(emissions + get_mtx_index(k, j * 3, a_size)), 2);
    if(continuous_score_all > 0.0) {
      gamma_p_j = e_res * continuous_score_j / continuous_score_all;
    }
    else {
      gamma_p_j = 0.0;
    }
    *(E + get_mtx_index(k, j * 3, a_size + 1)) += mean_value * gamma_p_j;
    *(E + get_mtx_index(k, j * 3 + 1, a_size + 1)) += varians * gamma_p_j;
    *(E + get_mtx_index(k, j * 3 + 2, a_size + 1)) += gamma_p_j;
  }
  
  *(E + get_mtx_index(k, j * 3, a_size + 1)) += e_res;

}


double add_Eka_contribution_msa_multi(struct hmm_multi_s *hmmp, struct msa_sequences_multi_s *msa_seq_infop,
				      struct forward_s *forw_mtx, struct backward_s *backw_mtx, int p,
				      int k, int i, int use_lead_columns)
{
  double e_res, e_res_1, e_res_2, e_res_3;
  
  /* get contribution from this position in the sequence */
  e_res_1 = (forw_mtx + get_mtx_index(i+1, k, hmmp->nr_v))->prob;
  e_res_2 = (backw_mtx + get_mtx_index(i+1, k, hmmp->nr_v))->prob;
  e_res = e_res_1 * e_res_2;
  if(e_res == 0) {
    return 0.0; /* no use updating with a zero value */
  }
  
#ifdef DEBUG_BW
  printf("forw_mtx = %f\n", e_res_1);
  printf("backw_mtx = %f\n", e_res_2);
  printf("res = %f\n", e_res);
#endif	    
  
  /* divide with total probability of current sequence */
  if(use_lead_columns == NO) {
    e_res_3 = (forw_mtx + get_mtx_index((msa_seq_infop->msa_seq_length)+1,
					hmmp->nr_v-1, hmmp->nr_v))->prob;
  }
  else {
    e_res_3 = (forw_mtx + get_mtx_index((msa_seq_infop->nr_lead_columns)+1,
				      hmmp->nr_v-1, hmmp->nr_v))->prob;
  }
  e_res = e_res / e_res_3;
  
#ifdef DEBUG_BW
  printf("total prob = %f\n", e_res_3);
  printf("e_res = %f\n", e_res);
#endif
  
  return e_res;
}


void recalculate_emiss_expectations_multi(struct hmm_multi_s *hmmp, double *E, int alphabet)
{
  int a, d, k, l;
  double e_res;
  int a_size;

  if(alphabet == 1) {
    a_size = hmmp->a_size;
  }
  if(alphabet == 2) {
    a_size = hmmp->a_size_2;
  }
  if(alphabet == 3) {
    a_size = hmmp->a_size_3;
  }
 if(alphabet == 4) {
    a_size = hmmp->a_size_4;
  } 
  

  for(a = 0; a < a_size;a++) {
    k = 0;
    for(d = 0; d < hmmp->nr_d; d++) {
      e_res = 0;
      l = 0;
      while(*(hmmp->distrib_groups + k) != END) {
	e_res += *(E + get_mtx_index(*(hmmp->distrib_groups + k),
				     a, a_size));
	k++;
	l++;
      }
#ifdef DEBUG_BW
      printf("e_res = %f, l = %d\n", e_res, l);
#endif
      e_res = e_res/(double)l;
      k = k - l;
      while(*(hmmp->distrib_groups + k) != END) {
	*(E + get_mtx_index(*(hmmp->distrib_groups + k),
			    a, a_size)) = e_res;
	k++;
      }
      k++;
    }
  }
}

void recalculate_trans_expectations_multi(struct hmm_multi_s *hmmp, double *T)
{
  int t,k,l;
  double t_res;
  k = 0;
  for(t = 0; t < hmmp->nr_ttg; t++) {
    t_res = 0;
    l = 0;
    while((hmmp->trans_tie_groups + k)->from_v != END) {
      t_res += *(T + get_mtx_index((hmmp->trans_tie_groups + k)->from_v, (hmmp->trans_tie_groups + k)->to_v, hmmp->nr_v));
      k++;
      l++;
    }
    t_res = t_res/(double)l;
    k = k - l;
    while((hmmp->trans_tie_groups + k)->from_v != END) {
      *(T + get_mtx_index((hmmp->trans_tie_groups + k)->from_v, (hmmp->trans_tie_groups + k)->to_v, hmmp->nr_v)) = t_res;
      k++;
    }
    k++;
  }
}


void update_trans_mtx_std_multi(struct hmm_multi_s *hmmp, double *T, int k)
{
  int l;
  double t_res_1, t_res_2;
  int i;

#ifdef DEBUG_BW_TRANS
  dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
  dump_T_matrix(hmmp->nr_v, hmmp->nr_v, T);
#endif

  t_res_1 = 0;
  for(l = 0; l < hmmp->nr_v-1; l++) {
    t_res_1 += *(T + get_mtx_index(k,l,hmmp->nr_v));
  }
  if(t_res_1 == 0.0) {
    /* no dividing with 0 */
  }
  else {
    for(l = 0; l < hmmp->nr_v-1; l++) /* l = to-vertex, k = from-vertex */  {
      t_res_2 = *(T + get_mtx_index(k, l, hmmp->nr_v));
      *(hmmp->transitions + get_mtx_index(k,l,hmmp->nr_v)) = t_res_2/t_res_1;
      if(t_res_2 != 0.0 ) {
	*(hmmp->log_transitions + get_mtx_index(k,l,hmmp->nr_v)) = 
	  log10(t_res_2/t_res_1);
      }
      else {
	*(hmmp->log_transitions + get_mtx_index(k,l,hmmp->nr_v)) = DEFAULT;
      }
    }
  }
  update_tot_trans_mtx_multi(hmmp);
}

void update_trans_mtx_pseudocount_multi(struct hmm_multi_s *hmmp, double *T, int k)
{
  int l;
  double t_res_1, t_res_2;
  int i;
  double pseudo_value;

  pseudo_value = TRANSITION_PSEUDO_VALUE;
  t_res_1 = 0.0;
  for(l = 0; l < hmmp->nr_v-1; l++) {
    t_res_1 += *(T + get_mtx_index(k,l,hmmp->nr_v));
    if(*(T + get_mtx_index(k,l,hmmp->nr_v)) != 0.0) {
      t_res_1 += pseudo_value;
    }
  }
  
 
  if(t_res_1 == 0.0) {
    /* no dividing with 0 */
  }

  else {
    for(l = 0; l < hmmp->nr_v-1; l++) /* l = to-vertex */  {
      t_res_2 = *(T + get_mtx_index(k, l, hmmp->nr_v));
     
      if(t_res_2 != 0.0) {
	t_res_2 += pseudo_value;
      }
      *(hmmp->transitions + get_mtx_index(k,l,hmmp->nr_v)) = t_res_2/t_res_1;
      if(t_res_2 != 0.0 ) {
	*(hmmp->log_transitions + get_mtx_index(k,l,hmmp->nr_v)) = 
	  log10(t_res_2/t_res_1);
      }
      else {
	*(hmmp->log_transitions + get_mtx_index(k,l,hmmp->nr_v)) = DEFAULT;
      }
    }
  }
  update_tot_trans_mtx_multi(hmmp);
}

void update_emiss_mtx_std_multi(struct hmm_multi_s *hmmp, double *E, int k, int alphabet)
{
  /* NOTE: k = current vertex */
  double e_res_1, e_res_2;
  int a_index; 
  int a_size;
  double *emissions;
  double *log_emissions;
  
  if(alphabet == 1) {
    a_size = hmmp->a_size;
    emissions = hmmp->emissions;
    log_emissions = hmmp->log_emissions;
  }
  if(alphabet == 2) {
    a_size = hmmp->a_size_2;
    emissions = hmmp->emissions_2;
    log_emissions = hmmp->log_emissions_2;
  }
  if(alphabet == 3) {
    a_size = hmmp->a_size_3;
    emissions = hmmp->emissions_3;
    log_emissions = hmmp->log_emissions_3;
  }
  if(alphabet == 4) {
    a_size = hmmp->a_size_4;
    emissions = hmmp->emissions_4;
    log_emissions = hmmp->log_emissions_4;
  }

  if(silent_state_multi(k, hmmp) == YES) {
    return; /* this is a silent state, no updating should be done */
  }
  else if(locked_state_multi(hmmp, k) == YES) {
    return; /* this state's parameters are locked, don't update */
  }


  e_res_1 = 0;
  for(a_index = 0; a_index < a_size; a_index++) {
    e_res_1 += *(E + get_mtx_index(k, a_index, a_size));
  }
  if(e_res_1 == 0.0) {
  }
  else {
    for(a_index = 0; a_index < a_size; a_index++) {
      e_res_2 = *(E + get_mtx_index(k, a_index, a_size));
      *(emissions + get_mtx_index(k, a_index, a_size)) = e_res_2/e_res_1;
      if(e_res_2 != 0) {
	*(log_emissions + get_mtx_index(k, a_index, a_size)) =
	  log10(e_res_2/e_res_1);
      }
      else {
	*(log_emissions + get_mtx_index(k, a_index, a_size)) = DEFAULT;
      }
    }
  }
}

void update_emiss_mtx_std_continuous_multi(struct hmm_multi_s *hmmp, double *E, int k, int alphabet)
{
  /* NOTE: k = current vertex */
  double e_res_1, e_res_2;
  int a_index; 
  int a_size;
  double *emissions;
  double *log_emissions;
  double gamma_all, gamma_j, mean_j, var_j;
  

  if(alphabet == 1) {
    a_size = hmmp->a_size;
    emissions = hmmp->emissions;
    log_emissions = hmmp->log_emissions;
  }
  if(alphabet == 2) {
    a_size = hmmp->a_size_2;
    emissions = hmmp->emissions_2;
    log_emissions = hmmp->log_emissions_2;
  }
  if(alphabet == 3) {
    a_size = hmmp->a_size_3;
    emissions = hmmp->emissions_3;
    log_emissions = hmmp->log_emissions_3;
  }
  if(alphabet == 4) {
    a_size = hmmp->a_size_4;
    emissions = hmmp->emissions_4;
    log_emissions = hmmp->log_emissions_4;
  }
  if(silent_state_multi(k, hmmp) == YES) {
    return; /* this is a silent state, no updating should be done */
  }
  else if(locked_state_multi(hmmp, k) == YES) {
    return; /* this state's parameters are locked, don't update */
  }
  else {
    gamma_all = *(E + get_mtx_index(k, a_size, a_size + 1));
    for(a_index = 0; a_index < a_size; a_index += 3) {
      gamma_j = *(E + get_mtx_index(k, a_index + 2, a_size + 1));
      mean_j = *(E + get_mtx_index(k, a_index, a_size + 1));
      var_j = *(E + get_mtx_index(k, a_index + 1, a_size + 1));
      if(gamma_all == 0.0) {
	*(emissions + get_mtx_index(k, a_index + 2, a_size)) = 0.0;
	*(log_emissions + get_mtx_index(k, a_index + 2, a_size)) = DEFAULT;
      }
      else {
	*(emissions + get_mtx_index(k, a_index + 2, a_size)) = gamma_j / gamma_all;
	*(log_emissions + get_mtx_index(k, a_index + 2, a_size)) = log10(gamma_j / gamma_all);
      }

      if(gamma_j == 0.0) {
	*(emissions + get_mtx_index(k, a_index, a_size)) = 0.0;
	*(log_emissions + get_mtx_index(k, a_index, a_size)) = DEFAULT;
	*(emissions + get_mtx_index(k, a_index + 1, a_size)) = 0.0;
	*(log_emissions + get_mtx_index(k, a_index + 1, a_size)) = DEFAULT;
      }
      else {
	*(emissions + get_mtx_index(k, a_index, a_size)) = mean_j / gamma_j;
	*(log_emissions + get_mtx_index(k, a_index, a_size)) = log10(mean_j / gamma_j);
	*(emissions + get_mtx_index(k, a_index + 1, a_size)) = var_j / gamma_j;
	*(log_emissions + get_mtx_index(k, a_index + 1, a_size)) = log10(var_j / gamma_j);
      }
    }
  }
}


void update_emiss_mtx_pseudocount_multi(struct hmm_multi_s *hmmp, double *E, int k, int alphabet)
{
  /* NOTE: k = current vertex */
  double e_res_1, e_res_2;
  int a_index;
  int pseudo_value;
  int a_size;
  double *emissions;
  double *log_emissions;

  if(alphabet == 1) {
    a_size = hmmp->a_size;
    emissions = hmmp->emissions;
    log_emissions = hmmp->log_emissions;
  }
  if(alphabet == 2) {
    a_size = hmmp->a_size_2;
    emissions = hmmp->emissions_2;
    log_emissions = hmmp->log_emissions_2;
  }
  if(alphabet == 3) {
    a_size = hmmp->a_size_3;
    emissions = hmmp->emissions_3;
    log_emissions = hmmp->log_emissions_3;
  }
  if(alphabet == 4) {
    a_size = hmmp->a_size_4;
    emissions = hmmp->emissions_4;
    log_emissions = hmmp->log_emissions_4;
  }
  
  pseudo_value = EMISSION_PSEUDO_VALUE;
  if(silent_state_multi(k,hmmp) == YES) {
    return; /* this is a silent state, no updating should be done */
  }
  else if(locked_state_multi(hmmp, k) == YES) {
    return; /* this state's parameters are locked, don't update */
  }
  
  e_res_1 = 0;
  for(a_index = 0; a_index < a_size; a_index++) {
    e_res_1 += *(E + get_mtx_index(k, a_index, a_size));
    if(*(E + get_mtx_index(k, a_index, a_size)) != 0.0) {
      e_res_1 += pseudo_value;
    }
  }
  if(e_res_1 == 0.0) {
  }
  else {
    for(a_index = 0; a_index < a_size; a_index++) {
      e_res_2 = *(E + get_mtx_index(k, a_index, a_size));
      if(e_res_2 != 0.0) {
	e_res_2 += pseudo_value;
      }
      *(emissions + get_mtx_index(k, a_index, a_size)) = e_res_2/e_res_1;
      if(e_res_2 != 0) {
	*(log_emissions + get_mtx_index(k, a_index, a_size)) =
	  log10(e_res_2/e_res_1);
      }
      else {
	*(log_emissions + get_mtx_index(k, a_index, a_size)) = DEFAULT;
      }
    }
  }
}


void update_emiss_mtx_prior_multi(struct hmm_multi_s *hmmp, double *E, int k, struct emission_dirichlet_s *priorp, int alphabet)
{

  /* NOTE: k = current vertex */
  int nr_components, comps, a_index; 
  double scaling_factor, X_sum, *X_values, ed_res1, E_sums, *logbeta_an_values;
  double q_value, exponent, prior_prob, tot_prior_prob;
  double prior_scaler;
  int a_size;
  double *emissions;
  double *log_emissions;


  if(alphabet == 1) {
    a_size = hmmp->a_size;
    emissions = hmmp->emissions;
    log_emissions = hmmp->log_emissions;
  }
  if(alphabet == 2) {
    a_size = hmmp->a_size_2;
    emissions = hmmp->emissions_2;
    log_emissions = hmmp->log_emissions_2;
  }
  if(alphabet == 3) {
    a_size = hmmp->a_size_3;
    emissions = hmmp->emissions_3;
    log_emissions = hmmp->log_emissions_3;
  }
  if(alphabet == 4) {
    a_size = hmmp->a_size_4;
    emissions = hmmp->emissions_4;
    log_emissions = hmmp->log_emissions_4;
  }

  /* Note that this function expects that a prior struct is loaded and ready for this particular alphabet */

  prior_scaler = *(hmmp->vertex_emiss_prior_scalers + k);
  if(*(emissions + get_mtx_index(k, 0, a_size)) == SILENT) {
    return; /* this is a silent state, no updating should be done */
  }
  else if(locked_state_multi(hmmp, k) == YES) {
    return; /* this state's parameters are locked, don't update */
  }
  
  nr_components = priorp->nr_components;
  logbeta_an_values = malloc_or_die(nr_components * sizeof(double));
  scaling_factor = 0.0 - FLT_MAX;
  X_sum = 0.0;
  X_values = malloc_or_die(a_size * sizeof(double));
  
  /* calculate logB(alpha + n) for all components +
   * calculate scaling factor for logB(alpha + n) - logB(alpha) */
  for(comps = 0; comps < nr_components; comps++) {
    ed_res1 = 0;
    E_sums = 0;
    for(a_index = 0; a_index < a_size; a_index++) {
      ed_res1 += lgamma(*(priorp->prior_values +
			  get_mtx_index(comps, a_index, a_size)) +
			*(E + get_mtx_index(k,a_index,a_size)));
      E_sums += *(E + get_mtx_index(k,a_index, a_size));
    }
    ed_res1 = ed_res1 - lgamma(*(priorp->alpha_sums + comps) + E_sums);
    *(logbeta_an_values + comps) = ed_res1;
    if((ed_res1 = ed_res1 - *(priorp->logbeta_values + comps)) > scaling_factor) {
      scaling_factor = ed_res1;
    }
  }
#ifdef DEBUG_PRIORS
  printf("ed_res1(top) = %f\n", ed_res1);
  printf("scaling_factor = %f\n", scaling_factor);
#endif

  /* calculate all the Xi's */
  for(a_index = 0; a_index < a_size; a_index++) {
    *(X_values + a_index) = 0;
    for(comps = 0; comps < nr_components; comps++) {
      q_value = *(priorp->q_values + comps);
      exponent = (*(logbeta_an_values + comps) - *(priorp->logbeta_values + comps) - 
		  scaling_factor);
      prior_prob = (*(priorp->prior_values + get_mtx_index(comps,a_index, a_size)) * prior_scaler +
		    *(E + get_mtx_index(k,a_index,a_size)));
      tot_prior_prob = (*(priorp->alpha_sums + comps) + E_sums);
      *(X_values + a_index) += q_value * exp(exponent) * prior_prob / tot_prior_prob;
#ifdef DEBUG_PRIORS
    printf("q_value = %f\n", q_value);
    printf("exponent = %f\n", exponent);
    printf("prior_prob = %f\n", prior_prob);
    printf("tot_prior_prob = %f\n", tot_prior_prob);
    printf("X_values[%d] = %f\n", a_index, *(X_values + a_index));
#endif
    }
    X_sum += *(X_values + a_index);
  }
  
  /* update emission matrix */
  for(a_index = 0; a_index < a_size; a_index++) {
    ed_res1 = *(X_values + a_index) / X_sum;
#ifdef DEBUG_PRIORS
    printf("ed_res1 = %f\n", ed_res1);
    printf("X_sum = %f\n", X_sum);
#endif
    if(ed_res1 != 0.0) {
      *(emissions + get_mtx_index(k, a_index, a_size)) = ed_res1;
      *(log_emissions + get_mtx_index(k, a_index, a_size)) = log10(ed_res1);
    }
    else {
      *(emissions + get_mtx_index(k, a_index, a_size)) = ed_res1;
      *(log_emissions + get_mtx_index(k, a_index, a_size)) = DEFAULT;
    }
  }

  free(logbeta_an_values);
  free(X_values);
}


void anneal_E_matrix_multi(double temperature, double *E, struct hmm_multi_s *hmmp, int alphabet)
{
  int i,j;
  double rand_nr;
  int a_size;
  
  if(alphabet == 1) {
    if(hmmp->alphabet_type == DISCRETE) {
      a_size = hmmp->a_size;
    }
    else {
      a_size = hmmp->a_size + 1;
    }
  }
  if(alphabet == 2) {
    if(hmmp->alphabet_type_2 == DISCRETE) {
      a_size = hmmp->a_size_2;
    }
    else {
      a_size = hmmp->a_size_2 + 1;
    }
  }
  if(alphabet == 3) {
    if(hmmp->alphabet_type_3 == DISCRETE) {
      a_size = hmmp->a_size_3;
    }
    else {
      a_size = hmmp->a_size_3 + 1;
    }
  }
  if(alphabet == 4) {
    if(hmmp->alphabet_type_4 == DISCRETE) {
      a_size = hmmp->a_size_4;
    }
    else {
      a_size = hmmp->a_size_4 + 1;
    }
  }
  

  srand(time(0));
  for(i = 1; i < hmmp->nr_v - 1; i++) {
    for(j = 0; j < a_size; j++) {
      rand_nr = (double)rand()/RAND_MAX;
      rand_nr = rand_nr * temperature;      
      *(E + get_mtx_index(i, j, a_size)) +=  *(E + get_mtx_index(i, j, a_size)) * rand_nr;
    }
  }
}

void anneal_T_matrix_multi(double temperature, double *T, struct hmm_multi_s *hmmp)
{
  int i,j;
  double rand_nr;

  srand(time(0));
  for(i = 1; i < hmmp->nr_v - 1; i++) {
    for(j = 0; j < hmmp->nr_v; j++) {
      rand_nr = (double)rand()/RAND_MAX;
      rand_nr = rand_nr * temperature;
      *(T + get_mtx_index(i, j, hmmp->nr_v)) +=  *(T + get_mtx_index(i, j, hmmp->nr_v)) * rand_nr;
    }
  }
}

void calculate_TE_contributions_multi(double *T, double *E, double *E_2, double *E_3, double *E_4,
				      double *T_lab, double *E_lab, double *E_lab_2, double *E_lab_3, double *E_lab_4,
				      double *T_ulab, double *E_ulab, double *E_ulab_2, double *E_ulab_3, double *E_ulab_4,
				      double *emissions, double *emissions_2, double *emissions_3, double *emissions_4,
				      double *transitions, int nr_v, int a_size, int a_size_2, int a_size_3, int a_size_4,
				      double *emiss_prior_scalers, double *emiss_prior_scalers_2, double *emiss_prior_scalers_3,
				      double *emiss_prior_scalers_4, int rd, int nr_alphabets) {
  int v,w;
  int a;
  int x,y, y_2, y_3, y_4;
  double rowsum;
  double T_divider, E_divider, E_divider_2, E_divider_3, E_divider_4, max_T_ulab;
  double max_E_ulab, max_E_ulab_2, max_E_ulab_3, max_E_ulab_4;
  double E_limiter, E_limiter_2, E_limiter_3, E_limiter_4, T_limiter;
  double DIVIDER_SCALER = 7.0;


#ifdef DEBUG_EXTBW
  printf("matrices before update emissions mtx\n");
  dump_E_matrix(nr_v, a_size, emissions);
#endif

#ifdef DEBUG_EXTBW
  printf("matrices before update ulab\n");
  dump_E_matrix(nr_v, a_size, E_ulab);
#endif

#ifdef DEBUG_EXTBW
  printf("matrices before update lab\n");
  dump_E_matrix(nr_v, a_size, E_lab);
#endif


  /* copy current emission and transition matrix values to E and T matrices + scale values */
  for(v = 0; v < nr_v-1; v++) {
    for(w = 1; w < nr_v-1; w++) {
      x = get_mtx_index(v,w,nr_v);
      *(T + x) = *(transitions + x);
    }
  }
  for(v = 1; v < nr_v-1; v++) {
    for(a = 0; a < a_size; a++) {
      y = get_mtx_index(v,a,a_size);
      *(E + y) = *(emissions + y);
    }
  }
  if(nr_alphabets > 1) {
    for(v = 1; v < nr_v-1; v++) {
      for(a = 0; a < a_size_2; a++) {
	y = get_mtx_index(v,a,a_size_2);
	*(E_2 + y) = *(emissions_2 + y);
      }
    }
  }
  if(nr_alphabets > 2) {
    for(v = 1; v < nr_v-1; v++) {
      for(a = 0; a < a_size_3; a++) {
	y = get_mtx_index(v,a,a_size_3);
	*(E_3 + y) = *(emissions_3 + y);
      }
    }
  }
  if(nr_alphabets > 3) {
    for(v = 1; v < nr_v-1; v++) {
      for(a = 0; a < a_size_4; a++) {
	y = get_mtx_index(v,a,a_size_4);
	*(E_4 + y) = *(emissions_4 + y);
      }
    }
  }
  
#ifdef DEBUG_EXTBW
  printf("matrices after emission copy update E matrices\n");
  dump_E_matrix(nr_v, a_size, E);
  dump_T_matrix(nr_v, nr_v, T);
#endif

  /* check matrices for potential negative values + compensate */
  T_divider = 1.0;
  E_divider = 1.0;
  E_divider_2 = 1.0;
  E_divider_3 = 1.0;
  E_divider_4 = 1.0;
  for(v = 0; v < nr_v-1; v++) {
    for(w = 1; w < nr_v-1; w++) {
      x = get_mtx_index(v,w,nr_v);
      T_limiter = *(T + x);
      if(*(T + x) > (1.0 - *(T + x))) {
	T_limiter = 1.0 - *(T + x);
      }
      if(*(T + x) != 0.0 && (*(T_ulab + x) - *(T_lab + x)) / T_limiter > T_divider) {
	T_divider = (*(T_ulab + x) - *(T_lab + x)) / T_limiter;
#ifdef DEBUG_EXTBW
	printf("needed T_divider\n");
#endif
      }
    }
  }
  for(v = 1; v < nr_v-1; v++) {
    for(a = 0; a < a_size; a++) {
      y = get_mtx_index(v,a,a_size);
      E_limiter = *(E + y);
      if(*(E + y) > (1.0 - *(E + y))) {
	E_limiter = 1.0 - *(E + y);
      }
      if(*(E + y) != 0.0 && (*(E_ulab + y) - *(E_lab + y)) / E_limiter > E_divider) {
	E_divider = (*(E_ulab + y) - *(E_lab + y)) / E_limiter;
#ifdef DEBUG_EXTBW	
	printf("needed E_divider\n");
	printf("E_divider = %f\n", E_divider);
#endif	
      }
    }
  }
  if(nr_alphabets > 1) {
    for(v = 1; v < nr_v-1; v++) {
      for(a = 0; a < a_size_2; a++) {
	y = get_mtx_index(v,a,a_size_2);
	E_limiter_2 = *(E_2 + y);
	if(*(E_2 + y) > (1.0 - *(E_2 + y))) {
	  E_limiter_2 = 1.0 - *(E_2 + y);
	}
	if(*(E_2 + y) != 0.0 && (*(E_ulab_2 + y) - *(E_lab_2 + y)) / E_limiter > E_divider_2) {
	  E_divider_2 = (*(E_ulab_2 + y) - *(E_lab_2 + y)) / *(E_2 + y);
#ifdef DEBUG_EXTBW	
	  printf("needed E_divider_2\n");
	  printf("E_divider_2 = %f\n", E_divider_2);
#endif	
	}
      }
    }
  }
  if(nr_alphabets > 2) {
    for(v = 1; v < nr_v-1; v++) {
      for(a = 0; a < a_size_3; a++) {
	y = get_mtx_index(v,a,a_size_3);
	E_limiter_3 = *(E_3 + y);
	if(*(E_3 + y) > (1.0 - *(E_3 + y))) {
	  E_limiter_3 = 1.0 - *(E_3 + y);
	}
	if(*(E_3 + y) != 0.0 && (*(E_ulab_3 + y) - *(E_lab_3 + y)) / E_limiter_3 > E_divider_3) {
	  E_divider_3 = (*(E_ulab_3 + y) - *(E_lab_3 + y)) / *(E_3 + y);
#ifdef DEBUG_EXTBW	
	  printf("needed E_divider\n");
#endif	
	}
      }
    }
  }
  if(nr_alphabets > 3) {
    for(v = 1; v < nr_v-1; v++) {
      for(a = 0; a < a_size_4; a++) {
	y = get_mtx_index(v,a,a_size_4);
	E_limiter_4 = *(E_4 + y);
	if(*(E_4 + y) > (1.0 - *(E_4 + y))) {
	E_limiter_4 = 1.0 - *(E_4 + y);
	}
	if(*(E_4 + y) != 0.0 && (*(E_ulab_4 + y) - *(E_lab_4 + y)) / E_limiter_4 > E_divider_4) {
	  E_divider_4 = (*(E_ulab_4 + y) - *(E_lab_4 + y)) / *(E_4 + y);
#ifdef DEBUG_EXTBW	
	  printf("needed E_divider\n");
#endif	
	}
      }
    }
  }

  T_divider = T_divider * DIVIDER_SCALER;;
  E_divider = E_divider * DIVIDER_SCALER;
  E_divider_2 = E_divider_2 * DIVIDER_SCALER;
  E_divider_3 = E_divider_3 * DIVIDER_SCALER;
  E_divider_4 = E_divider_4 * DIVIDER_SCALER;
  
  /* add T_lab - T_ulab and E_lab - E_ulab values to T and E matrices */
  for(v = 0; v < nr_v-1; v++) {
    for(w = 1; w < nr_v-1; w++) {
      x = get_mtx_index(v,w,nr_v);
      *(T + x) += (*(T_lab + x) - *(T_ulab + x)) / T_divider;
      if(*(T + x) < 0.0) {
	*(T + x) = 0.0;
      }
    }
  }
  for(v = 1; v < nr_v-1; v++) {
    for(a = 0; a < a_size; a++) {
      y = get_mtx_index(v,a,a_size);
      *(E + y) += (*(E_lab + y) - *(E_ulab + y)) / (E_divider * *(emiss_prior_scalers + v));
      if(*(E + y) < 0.0) {
	*(E + y) = 0.0;
      }
    }
  }
  if(nr_alphabets > 1) {
    for(v = 1; v < nr_v-1; v++) {
      for(a = 0; a < a_size_2; a++) {
	y = get_mtx_index(v,a,a_size_2);
	*(E_2 + y) += (*(E_lab_2 + y) - *(E_ulab_2 + y)) / (E_divider_2 * *(emiss_prior_scalers_2 + v));
	if(*(E_2 + y) < 0.0) {
	  *(E_2 + y) = 0.0;
	}
      }
    }
  }
  if(nr_alphabets > 2) {
    for(v = 1; v < nr_v-1; v++) {
      for(a = 0; a < a_size_3; a++) {
	y = get_mtx_index(v,a,a_size_3);
	*(E_3 + y) += (*(E_lab_3 + y) - *(E_ulab_3 + y)) / (E_divider_3 * *(emiss_prior_scalers_3 + v));
	if(*(E_3 + y) < 0.0) {
	  *(E_3 + y) = 0.0;
	}
      }
    }
  }
  if(nr_alphabets > 3) {
    for(v = 1; v < nr_v-1; v++) {
      for(a = 0; a < a_size_4; a++) {
	y = get_mtx_index(v,a,a_size_4);
	*(E_4 + y) += (*(E_lab_4 + y) - *(E_ulab_4 + y)) / (E_divider_4 * *(emiss_prior_scalers_4 + v));
	if(*(E_4 + y) < 0.0) {
	  *(E_4 + y) = 0.0;
	}
      }
    }
  }
  

#ifdef DEBUG_EXTBW
  printf("matrices after lab/ulab update E matrices\n");
  dump_E_matrix(nr_v, a_size, E);

  printf("T matrices after update\n");
  dump_T_matrix(nr_v, nr_v, T);
#endif
  
  /* weight matrices with prior distribution */
  /* not implemented, not obvious if it will work */


  /* normalize matrices */
  for(v = 0; v < nr_v-1; v++) {
    rowsum = 0.0;
    for(w = 1; w < nr_v-1; w++) {
      x = get_mtx_index(v,w,nr_v);
      rowsum += *(transitions + x);
    }
    if(rowsum != 0.0) {
      for(w = 1; w < nr_v-1; w++) {
	x = get_mtx_index(v,w,nr_v);
	*(transitions + x) = *(transitions + x) / rowsum;
      }
    }
  }
  for(v = 1; v < nr_v-1; v++) {
    rowsum = 0.0;
    for(a = 0; a < a_size; a++) {
      y = get_mtx_index(v,a,a_size);
      rowsum += *(emissions + y);
    }
    if(rowsum != 0.0) {
      for(a = 0; a < a_size; a++) {
	y = get_mtx_index(v,a,a_size);
	*(emissions + y) = *(emissions + y) / rowsum;
      }
    }
  }
  if(nr_alphabets > 1) {
    for(v = 1; v < nr_v-1; v++) {
      rowsum = 0.0;
      for(a = 0; a < a_size_2; a++) {
	y = get_mtx_index(v,a,a_size_2);
	rowsum += *(emissions_2 + y);
      }
      if(rowsum != 0.0) {
	for(a = 0; a < a_size_2; a++) {
	  y = get_mtx_index(v,a,a_size_2);
	  *(emissions_2 + y) = *(emissions_2 + y) / rowsum;
	}
      }
    }
  }
  if(nr_alphabets > 2) {
    for(v = 1; v < nr_v-1; v++) {
      rowsum = 0.0;
      for(a = 0; a < a_size_3; a++) {
	y = get_mtx_index(v,a,a_size_3);
	rowsum += *(emissions_3 + y);
      }
      if(rowsum != 0.0) {
	for(a = 0; a < a_size_3; a++) {
	  y = get_mtx_index(v,a,a_size_3);
	  *(emissions_3 + y) = *(emissions_3 + y) / rowsum;
	}
      }
    }
  }
  if(nr_alphabets > 3) {
    for(v = 1; v < nr_v-1; v++) {
      rowsum = 0.0;
      for(a = 0; a < a_size_4; a++) {
	y = get_mtx_index(v,a,a_size_4);
	rowsum += *(emissions_4 + y);
      }
      if(rowsum != 0.0) {
	for(a = 0; a < a_size_4; a++) {
	  y = get_mtx_index(v,a,a_size_4);
	  *(emissions_4 + y) = *(emissions_4 + y) / rowsum;
	}
      }
    }
  }
#ifdef DEBUG_EXTBW
  printf("E matrices after update\n");
  dump_E_matrix(nr_v, a_size, E);

  printf("T matrices after update\n");
  dump_T_matrix(nr_v, nr_v, T);
  printf("***************************************\n");
#endif

}

/****************************some utility methods***********************************************/
int silent_state_multi(int k, struct hmm_multi_s *hmmp)
{ 
  if(*(hmmp->emissions + get_mtx_index(k,0,hmmp->a_size)) == SILENT) {
    return YES;
  }
  else {
    return NO;
  }
}

void update_tot_trans_mtx_multi(struct hmm_multi_s *hmmp)
{
  int v,w;
  struct path_element *wp;
  double t_res;
  
#ifdef DEBUG_BW
  //hmmp->tot_transitions = (double*)malloc_or_die(hmmp->nr_v * hmmp->nr_v * sizeof(double)); 
#endif
  
  /***************** changed to loop over trans to end state as well, may not always work *********************/

  for(v = 0; v < hmmp->nr_v;v++) {
    wp = *(hmmp->from_trans_array + v);
    while(wp->vertex != END) /* w = from-vertex */ {
      t_res = 1.0;
      w = wp->vertex;
      while(wp->next != NULL) {
	t_res = t_res * *((hmmp->transitions) +
			    get_mtx_index(wp->vertex, (wp + 1)->vertex, hmmp->nr_v));
	/* probability of transition from w to v via silent states */
	wp++;
      }
      t_res = t_res * *((hmmp->transitions) +
			  get_mtx_index(wp->vertex, v, hmmp->nr_v));
      *(hmmp->tot_transitions + get_mtx_index(w,v,hmmp->nr_v)) = t_res;
      wp++;
    }
  }
#ifdef DEBUG_BW
  dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->tot_transitions);
  dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
#endif
}

/* the general method for adding values to E matrix */
void add_to_E_multi(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p, int k, int a_size, int normalize,
		    double *subst_mtx, int alphabet, int scoring_method, int use_nr_occ, int alphabet_type, double *emissions)
{

  if(alphabet_type == CONTINUOUS) {
    add_to_E_continuous(E, Eka_base, msa_seq, p, k, a_size, emissions);
  }
  else if(scoring_method == DOT_PRODUCT && use_nr_occ == YES) {
    add_to_E_dot_product_nr_occ(E, Eka_base, msa_seq, p, k, a_size, normalize);
  }
  else if(scoring_method == DOT_PRODUCT_PICASSO && use_nr_occ == YES) {
    add_to_E_dot_product_picasso_nr_occ(E, Eka_base, msa_seq, p, k, a_size, normalize);
  }
  else if(scoring_method == PICASSO && use_nr_occ == YES) {
    add_to_E_picasso_nr_occ(E, Eka_base, msa_seq, p, k, a_size, normalize);
  }
  else if(scoring_method == PICASSO_SYM && use_nr_occ == YES) {
    add_to_E_picasso_sym_nr_occ(E, Eka_base, msa_seq, p, k, a_size, normalize);
  }
  else if(scoring_method == SJOLANDER && use_nr_occ == YES) {
    add_to_E_sjolander_score_nr_occ(E, Eka_base, msa_seq, p, k, a_size, normalize);
  }
  else if(scoring_method == SJOLANDER_REVERSED && use_nr_occ == YES) {
    add_to_E_sjolander_reversed_score_nr_occ(E, Eka_base, msa_seq, p, k, a_size, normalize);
  }
  else if(scoring_method == SUBST_MTX_PRODUCT && use_nr_occ == YES) {
    add_to_E_subst_mtx_product_nr_occ(E, Eka_base, msa_seq, p, k, a_size, normalize, subst_mtx);
  }
  else if(scoring_method == SUBST_MTX_DOT_PRODUCT && use_nr_occ == YES) {
    add_to_E_subst_mtx_dot_product_nr_occ(E, Eka_base, msa_seq, p, k, a_size, normalize, subst_mtx,
					  alphabet);
  }
  else if(scoring_method == SUBST_MTX_DOT_PRODUCT_PRIOR && use_nr_occ == YES) {
    add_to_E_subst_mtx_dot_product_prior_nr_occ(E, Eka_base, msa_seq, p, k, a_size, normalize, subst_mtx,
						alphabet);
  }
  else if(scoring_method == DOT_PRODUCT) {
    add_to_E_dot_product(E, Eka_base, msa_seq, p, k, a_size, normalize);
  }
  else if(scoring_method == DOT_PRODUCT_PICASSO) {
    add_to_E_dot_product_picasso(E, Eka_base, msa_seq, p, k, a_size, normalize);
  }
  else if(scoring_method == PICASSO) {
    add_to_E_picasso(E, Eka_base, msa_seq, p, k, a_size, normalize);
  }
  else if(scoring_method == PICASSO_SYM) {
    add_to_E_picasso_sym(E, Eka_base, msa_seq, p, k, a_size, normalize);
  }
  else if(scoring_method == SJOLANDER) {
    add_to_E_sjolander_score(E, Eka_base, msa_seq, p, k, a_size, normalize);
  }
  else if(scoring_method == SJOLANDER_REVERSED) {
    add_to_E_sjolander_reversed_score(E, Eka_base, msa_seq, p, k, a_size, normalize);
  }
  else if(scoring_method == SUBST_MTX_PRODUCT) {
    add_to_E_subst_mtx_product(E, Eka_base, msa_seq, p, k, a_size, normalize, subst_mtx);
  }
  else if(scoring_method == SUBST_MTX_DOT_PRODUCT) {
    add_to_E_subst_mtx_dot_product(E, Eka_base, msa_seq, p, k, a_size, normalize, subst_mtx,
				   alphabet);
  }
  else if(scoring_method == SUBST_MTX_DOT_PRODUCT_PRIOR) {
    add_to_E_subst_mtx_dot_product_prior(E, Eka_base, msa_seq, p, k, a_size, normalize, subst_mtx,
					 alphabet);
  }
  else {
    printf("Error: Unrecognized scoring method\n");
    exit(0);
  }
}

