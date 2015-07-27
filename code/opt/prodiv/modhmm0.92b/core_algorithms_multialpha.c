#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "structs.h"
#include "funcs.h"

#define V_LIST_END  -99
#define V_LIST_NEXT  -9

//#define DEBUG_LABELING_UPDATE
//#define DEBUG_DEALLOCATE_LABELINGS
//#define DEBUG_FW
//#define DEBUG_BW
//#define DEBUG_VI
//#define PATH

extern int verbose;

double dot_product_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
			 struct msa_sequences_multi_s *msa_seq_infop, int c, int v, int normalize, int multi_scoring_method);
double dot_prouct_picasso_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
				struct msa_sequences_multi_s *msa_seq_infop, int c, int v, int normalize, int multi_scoring_method,
				double *aa_freqs, double *aa_freqs_2, double *aa_freqs_3, double *aa_freqs_4);
double picasso_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
			 struct msa_sequences_multi_s *msa_seq_infop, int c, int v, int normalize, int multi_scoring_method,
			 double *aa_freqs, double *aa_freqs_2, double *aa_freqs_3, double *aa_freqs_4);
double picasso_sym_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
			 struct msa_sequences_multi_s *msa_seq_infop, int c, int v, int normalize, int multi_scoring_method,
			 double *aa_freqs, double *aa_freqs_2, double *aa_freqs_3, double *aa_freqs_4);
double sjolander_score_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
			     struct msa_sequences_multi_s *msa_seq_infop, int c, int v, int normalize, int multi_scoring_method);
double sjolander_reversed_score_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
				      struct msa_sequences_multi_s *msa_seq_infop, int c, int v, int normalize,
				      int multi_scoring_method);
double subst_mtx_product_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
			       struct msa_sequences_multi_s *msa_seq_infop, int c, int v, int normalize, int multi_scoring_method);
double subst_mtx_dot_product_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
				   struct msa_sequences_multi_s *msa_seq_infop, int c, int v, int a_index, int a_index_2,
				   int a_index_3, int a_index_4,
				   int normalize, int multi_scoring_method);
double subst_mtx_dot_product_prior_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
					 struct msa_sequences_multi_s *msa_seq_infop, int c, int v, int a_index, int a_index_2,
					 int a_index_3, int a_index_4,
					 int normalize, int multi_scoring_method);
double get_msa_emission_score_multi(struct msa_sequences_multi_s *msa_seq_infop, struct hmm_multi_s *hmmp, int c, int v,
				    int use_labels, int normalize, int a_index, int a_index_2,
				    int a_index_3, int a_index_4, int scoring_method, int use_gap_shares,
				    int use_prior_shares, int multi_scoring_method, double *aa_freqs, double *aa_freqs_2,
				    double *aa_freqs_3, double *aa_freqs_4);
double get_single_emission_score_multi(struct hmm_multi_s *hmmp, struct letter_s *seq, struct letter_s *seq_2,
				       struct letter_s *seq_3, struct letter_s *seq_4, int c, int v, int replacement_letter_c,
				       int replacement_letter_c_2, int replacement_letter_c_3,
				       int replacement_letter_c_4, int use_labels, int a_index, int a_index_2,
				       int a_index_3, int a_index_4, int multi_scoring_method);


/************************* the forward algorithm **********************************/
int forward_multi(struct hmm_multi_s *hmmp, struct letter_s *s, struct letter_s *s_2, struct letter_s *s_3,
		  struct letter_s *s_4, struct forward_s **ret_forw_mtxpp,
		  double **ret_scale_fspp, int use_labels, int multi_scoring_method)
{
  struct forward_s *forw_mtxp, *cur_rowp, *prev_rowp; /* pointers to forward matrix */
  double *scale_fsp; /* pointer to array of scaling factors */
  struct letter_s *seq, *seq_2, *seq_3, *seq_4; /* pointer to the sequence */
  int i,j; /* loop variables */
  int seq_len; /* length of the sequence */
  int c, v, w; /* looping indices, c loops over the sequence,
		  v and w over the vertices in the HMM */
  double row_sum, res, t_res1, t_res2, t_res3; /* temporary variables to calculate probabilities */
  double t_res3_1, t_res3_2, t_res3_3, t_res3_4;
  int nr_v;
  struct path_element *wp; /* for traversing the paths from v to w */
  int a_index, a_index_2, a_index_3, a_index_4; /* holds the current letter's position in the alphabet */
  int replacement_letter_c, replacement_letter_c_2, replacement_letter_c_3, replacement_letter_c_4;
  
#ifdef DEBUG_FW
  printf("running forward\n");
#endif
  /* Allocate memory for matrix and scaling factors + some initial setup:
   * Note 1: forward probability matrix has the sequence indices vertically
   * and the vertex indices horizontally meaning it will be filled row by row
   * Note 2: *ret_forw_mtxpp and *ret_scale_fspp are allocated here, but must be
   * freed by caller */
  nr_v = hmmp->nr_v;
  seq_len = get_seq_length(s);
  *ret_forw_mtxpp = (struct forward_s*)(malloc_or_die((seq_len+2) *
						      hmmp->nr_v *
						      sizeof(struct forward_s)));
  forw_mtxp = *ret_forw_mtxpp;
  *ret_scale_fspp = (double*)(malloc_or_die((seq_len+2) * sizeof(double)));
  scale_fsp = *ret_scale_fspp;
  
  /* Convert sequence to 1...L for easier indexing */
  seq = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
  memcpy(seq+1, s, seq_len * sizeof(struct letter_s));
  if(hmmp->nr_alphabets > 1) {
    seq_2 = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
    memcpy(seq_2+1, s_2, seq_len * sizeof(struct letter_s));
  }
  if(hmmp->nr_alphabets > 2) {
    seq_3 = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
    memcpy(seq_3+1, s_3, seq_len * sizeof(struct letter_s));
  }
  if(hmmp->nr_alphabets > 3) {
    seq_4 = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
    memcpy(seq_4+1, s_4, seq_len * sizeof(struct letter_s));
  }

  /* Initialize the first row of the matrix */
  forw_mtxp->prob = 1.0; /* sets index (0,0) to 1.0,
			    the rest are already 0.0 as they should be */
  *scale_fsp = 1.0;

  /* Fill in middle rows */
  prev_rowp = forw_mtxp;
  cur_rowp = forw_mtxp + hmmp->nr_v;
  for(c = 1; c <= seq_len; c++) {
    
    /* get alphabet index for c*/
    if(hmmp->alphabet_type == DISCRETE) {
      replacement_letter_c = NO;
      a_index = get_alphabet_index(&seq[c], hmmp->alphabet, hmmp->a_size);
      if(a_index < 0) {
	a_index = get_replacement_letter_index_multi(&seq[c], hmmp->replacement_letters, 1);
	if(a_index < 0) {
	  printf("Letter '%s' is not in alphabet\n", (&seq[c])->letter);
	  return NOPROB;
	}
	replacement_letter_c = YES;
      }
    }
    if(hmmp->nr_alphabets > 1) {
      if(hmmp->alphabet_type_2 == DISCRETE) {
	replacement_letter_c_2 = NO;
	a_index_2 = get_alphabet_index(&seq_2[c], hmmp->alphabet_2, hmmp->a_size_2);
	if(a_index_2 < 0) {
	  a_index_2 = get_replacement_letter_index_multi(&seq_2[c], hmmp->replacement_letters, 2);
	  if(a_index_2 < 0) {
	    printf("Letter '%s' is not in alphabet\n", (&seq_2[c])->letter);
	    return NOPROB;
	  }
	  replacement_letter_c_2 = YES;
	}
      }
    }
    if(hmmp->nr_alphabets > 2) {
      if(hmmp->alphabet_type_3 == DISCRETE) {
	replacement_letter_c_3 = NO;
	a_index_3 = get_alphabet_index(&seq_3[c], hmmp->alphabet_3, hmmp->a_size_3);
	if(a_index_3 < 0) {
	  a_index_3 = get_replacement_letter_index_multi(&seq_3[c], hmmp->replacement_letters, 3);
	  if(a_index_3 < 0) {
	    printf("Letter '%s' is not in alphabet\n", (&seq_3[c])->letter);
	    return NOPROB;
	  }
	  replacement_letter_c_3 = YES;
	}
      }
    }
    if(hmmp->nr_alphabets > 3) {
      if(hmmp->alphabet_type_4 == DISCRETE) {
	replacement_letter_c_4 = NO;
	a_index_4 = get_alphabet_index(&seq_4[c], hmmp->alphabet_4, hmmp->a_size_4);
	if(a_index_4 < 0) {
	  a_index_4 = get_replacement_letter_index_multi(&seq_4[c], hmmp->replacement_letters, 4);
	  if(a_index_4 < 0) {
	    printf("Letter '%s' is not in alphabet\n", (&seq_4[c])->letter);
	    return NOPROB;
	  }
	  replacement_letter_c_4 = YES;
	}
      }
    }
#ifdef DEBUG_FW  
    printf("seq[c] = %s\n", &seq[c]);
    if(hmmp->alphabet_type == DISCRETE) {
      printf("a_index = %d\n", a_index);
    }
#endif
    /* calculate sum of probabilities */
    row_sum = 0;
    for(v = 1; v < hmmp->nr_v - 1; v++) /* v = to-vertex */ {
#ifdef DEBUG_FW
      printf("prob to vertex %d\n", v);
#endif
      res = 0.0;
      wp = *(hmmp->tot_from_trans_array + v);
      while((w = wp->vertex) != END) /* w = from-vertex */ {
	/* calculate intermediate results */
	res += (prev_rowp + w)->prob * *(hmmp->tot_transitions + (w * nr_v + v));
	if(*(hmmp->tot_transitions + (w * nr_v + v)) < 0) {
	  printf("found model transition prob from %d to %d < 0.0\n", w, v);
	  exit(0);
	}
	wp++;
#ifdef DEBUG_FW
	printf("prev = %f: ", (prev_rowp + w)->prob);
	printf("trans = %f\n",  *(hmmp->tot_transitions + (w * nr_v + v)));
#endif
      }
     
      
      /* calculate prob of producing letter l in v */
      t_res3 = get_single_emission_score_multi(hmmp, seq, seq_2, seq_3, seq_4, c, v, replacement_letter_c, replacement_letter_c_2,
					       replacement_letter_c_3, replacement_letter_c_4, use_labels, a_index,
					       a_index_2, a_index_3, a_index_4, multi_scoring_method);
      
      res = res * t_res3;
      row_sum += res;
	
      /* save result in matrices */
      (cur_rowp + v)->prob = res;
      
      
#ifdef DEBUG_FW
      printf("res = %f\n", res);
#endif
    }
    
    /* scale the results, row_sum = the total probability of 
     * having produced the sequence up to and including character c */
    scale_fsp[c] = row_sum;
#ifdef DEBUG_FW
    printf("rowsum = %f\n",row_sum);
    printf("scaling set to: %f\n", scale_fsp[c]);
#endif
    if(row_sum == 0.0) {
      printf("Sequence cannot be produced by this hmm\n");
      sequence_as_string(s);
      printf("pos = %d\n",c);
      return NOPROB;
    }
    for(v = 0; v < hmmp->nr_v; v++) {
      if((cur_rowp + v)->prob != 0){
	(cur_rowp + v)->prob = ((cur_rowp + v)->prob)/row_sum; /* scaling */
      }
    }
    
    /* move row pointers one row forward */
    prev_rowp = cur_rowp;
    cur_rowp = cur_rowp + hmmp->nr_v;
  }

 
  /* Fill in transition to end state */
  res = 0;
  wp = *(hmmp->tot_from_trans_array + hmmp->nr_v - 1);
  while(wp->vertex != END) {
    t_res1 = (prev_rowp + wp->vertex)->prob;
    t_res2 = *((hmmp->tot_transitions) + get_mtx_index(wp->vertex, hmmp->nr_v-1, hmmp->nr_v));
    if(t_res2 > 1.0) {
      t_res2 = 1.0;
    }
    res += t_res1 * t_res2;
    wp++;
  }


#ifdef DEBUG_FW
  printf("res = %f\n", res);
#endif
  (cur_rowp + hmmp->nr_v - 1)->prob = res; /* obs: no scaling performed here */

#ifdef DEBUG_FW
  dump_forward_matrix(seq_len + 2, hmmp->nr_v, forw_mtxp);
  dump_scaling_array(seq_len + 1, scale_fsp);
#endif
  
  /* Garbage collection and return */
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
  return OK;
}



/**************************** the backward algorithm **********************************/
int backward_multi(struct hmm_multi_s *hmmp, struct letter_s *s, struct letter_s *s_2, struct letter_s *s_3,
		   struct letter_s *s_4, struct backward_s **ret_backw_mtxpp,
		   double *scale_fsp, int use_labels, int multi_scoring_method)
{
  struct backward_s *backw_mtxp, *cur_rowp, *prev_rowp; /* pointers to backward matrix */
  struct letter_s *seq, *seq_2, *seq_3, *seq_4; /* pointer to the sequence */
  int i; /* loop index */
  int seq_len; /* length of the sequence */
  int c, v, w; /* looping indices, c loops over the sequence,
		  v and w over the vertices in the HMM */
  int a_index, a_index_2, a_index_3, a_index_4; /* holds the current letter's position in the alphabet */
  double row_sum, res,t_res1,t_res2, t_res3; /* temporary variables to calculate results */
  double t_res3_1, t_res3_2, t_res3_3, t_res3_4;
  struct path_element *wp;
  int replacement_letter_c, replacement_letter_c_2, replacement_letter_c_3, replacement_letter_c_4;
  
  /* Allocate memory for matrix + some initial setup:
   * Note 1: probability matrix has the sequence indices vertically
   * and the vertex indices horizontally meaning it will be filled row by row
   * Note 2: *ret_backw_mtxpp is allocated here, but must be
   * freed by caller
   * Note 3: *scale_fspp is the scaling array produced by forward() meaning backward()
   * must be called after forward() */
  seq_len = get_seq_length(s);
  *ret_backw_mtxpp = (struct backward_s*)(malloc_or_die((seq_len+2) *
							hmmp->nr_v * 
							sizeof(struct backward_s)));
  backw_mtxp = *ret_backw_mtxpp;
  
  /* Convert sequence to 1...L for easier indexing */
  seq = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
  memcpy(seq+1, s, seq_len * sizeof(struct letter_s));
  if(hmmp->nr_alphabets > 1) {
    seq_2 = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
    memcpy(seq_2+1, s_2, seq_len * sizeof(struct letter_s));
  }
  if(hmmp->nr_alphabets > 2) {
    seq_3 = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
    memcpy(seq_3+1, s_3, seq_len * sizeof(struct letter_s));
  }
  if(hmmp->nr_alphabets > 3) {
    seq_4 = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
    memcpy(seq_4+1, s_4, seq_len * sizeof(struct letter_s));
  }


  /* Initialize the last row of the matrix */
  (backw_mtxp + get_mtx_index(seq_len + 1, hmmp->nr_v - 1,
			      hmmp->nr_v))->prob = 1.0; /* sets index
							 * (seq_len+1,nr_v-1) i.e.
							 * the lower right
							 * corner of the matrix to 
							 * 1.0,the rest are already
							 * 0.0 as they should be*/
  
  /* Fill in next to last row in matrix (i.e. add prob for path to end state for all states
    * that have a transition path to the end state) */
  prev_rowp = backw_mtxp + get_mtx_index(seq_len + 1 , 0, hmmp->nr_v);
  cur_rowp = prev_rowp - hmmp->nr_v;
  wp = *(hmmp->from_trans_array + hmmp->nr_v - 1);
  
  while(wp->vertex != END) {
    w = wp->vertex;
    t_res1 = 1.0;
    while(wp->next != NULL) {
      t_res1 = t_res1 * *(hmmp->transitions +
			  get_mtx_index(wp->vertex, (wp + 1)->vertex, hmmp->nr_v));
      wp++;
    }
    (cur_rowp + w)->prob += t_res1;
    if((cur_rowp + w)->prob > 1.0) {
      (cur_rowp + w)->prob = 1.0;
    }
    wp++;
  }
  prev_rowp = cur_rowp;
  cur_rowp = cur_rowp - hmmp->nr_v;

  
  /* Fill in first rows moving upwards in matrix */
  for(c = seq_len; c >= 1; c--) {
    
    /* get alphabet index for c */
#ifdef DEBUG_BW
    printf("c = %d\n", c);
#endif

    /* get alphabet index for c*/
    if(hmmp->alphabet_type == DISCRETE) {
      replacement_letter_c = NO;
      a_index = get_alphabet_index(&seq[c], hmmp->alphabet, hmmp->a_size);
      if(a_index < 0) {
	a_index = get_replacement_letter_index_multi(&seq[c], hmmp->replacement_letters, 1);
	if(a_index < 0) {
	  printf("Letter '%s' is not in alphabet\n", (&seq[c])->letter);
	  return NOPROB;
	}
	replacement_letter_c = YES;
      }
    }
    if(hmmp->nr_alphabets > 1) {
      if(hmmp->alphabet_type_2 == DISCRETE) {
	replacement_letter_c_2 = NO;
	a_index_2 = get_alphabet_index(&seq_2[c], hmmp->alphabet_2, hmmp->a_size_2);
	if(a_index_2 < 0) {
	  a_index_2 = get_replacement_letter_index_multi(&seq_2[c], hmmp->replacement_letters, 2);
	  if(a_index_2 < 0) {
	    printf("Letter '%s' is not in alphabet\n", (&seq_2[c])->letter);
	    return NOPROB;
	  }
	  replacement_letter_c_2 = YES;
	}
      }
    }
    if(hmmp->nr_alphabets > 2) {
      if(hmmp->alphabet_type_3 == DISCRETE) {
	replacement_letter_c_3 = NO;
	a_index_3 = get_alphabet_index(&seq_3[c], hmmp->alphabet_3, hmmp->a_size_3);
	if(a_index_3 < 0) {
	  a_index_3 = get_replacement_letter_index_multi(&seq_3[c], hmmp->replacement_letters, 3);
	  if(a_index_3 < 0) {
	    printf("Letter '%s' is not in alphabet\n", (&seq_3[c])->letter);
	    return NOPROB;
	  }
	  replacement_letter_c_3 = YES;
	}
      }
    }
    if(hmmp->nr_alphabets > 3) {
      if(hmmp->alphabet_type_4 == DISCRETE) {
	replacement_letter_c_4 = NO;
	a_index_4 = get_alphabet_index(&seq_4[c], hmmp->alphabet_4, hmmp->a_size_4);
	if(a_index_4 < 0) {
	  a_index_4 = get_replacement_letter_index_multi(&seq_4[c], hmmp->replacement_letters, 4);
	  if(a_index_4 < 0) {
	  printf("Letter '%s' is not in alphabet\n", (&seq_4[c])->letter);
	  return NOPROB;
	  }
	  replacement_letter_c_4 = YES;
	}
      }
    }
    
    
    /* calculate sum of probabilities */
    for(v = 0; v < hmmp->nr_v - 1; v++) /* v = from-vertex */{
#ifdef DEBUG_BW
      printf("prob passing through vertex %d:\n", v);
#endif
      
      res = 0;
      wp = *(hmmp->tot_to_trans_array + v);
      while(wp->vertex != END) 	/* w = to-vertex */ {
	/* total probability of transiting from v to w on all possible path (via silent states) */
	t_res2 = *((hmmp->tot_transitions) + get_mtx_index(v, wp->vertex, hmmp->nr_v));
	if(t_res2 < 0) {
	  printf("found model transition prob from %d to %d < 0.0\n", v, wp->vertex);
	  exit(0);
	}
	t_res1 = (prev_rowp + wp->vertex)->prob; /* probability of having produced the
						  * sequence after c passing through
						  * vertex w */
	

	/* calculate prob of producing letter l in v */
	t_res3 = get_single_emission_score_multi(hmmp, seq, seq_2, seq_3, seq_4, c, wp->vertex,
						 replacement_letter_c, replacement_letter_c_2,
						 replacement_letter_c_3, replacement_letter_c_4, use_labels, a_index,
						 a_index_2, a_index_3, a_index_4, multi_scoring_method);
	

	res += t_res1 * t_res2 * t_res3;
	wp++;
      }
#ifdef DEBUG_BW
      printf("prev = %f: ", t_res1);
      printf("trans = %f\n", t_res2);
#endif
      
      /* save result in matrices */
      (cur_rowp + v)->prob = res;
#ifdef DEBUG_BW
      printf("res = %f\n", res);
#endif
    }
    
    /* scale the results using the scaling factors produced by forward() */
#ifdef DEBUG_BW
    printf("scaling set to: %f\n", scale_fsp[c]);
    printf("and c = %d\n", c);
    dump_scaling_array(seq_len + 1, scale_fsp);
#endif
    for(v = 0; v < hmmp->nr_v; v++) {
      if((cur_rowp + v)->prob != 0){
	if(scale_fsp[c] != 0.0) { 
	  (cur_rowp + v)->prob = ((cur_rowp + v)->prob)/scale_fsp[c]; /* scaling */
	}
	else {
	  printf("Sequence cannot be produced by this hmm\n"); 
	  sequence_as_string(s);
	  return NOPROB;
	}
      }
    }
    
    /* move row pointers one row backward */
    prev_rowp = cur_rowp;
    cur_rowp = cur_rowp - hmmp->nr_v;
  }

#ifdef DEBUG_BW
  dump_backward_matrix(seq_len + 2, hmmp->nr_v, backw_mtxp);
  dump_scaling_array(seq_len + 1, scale_fsp);
#endif

  /* Garbage collection and return */
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
  return OK;
}


/************************* the viterbi algorithm **********************************/
int viterbi_multi(struct hmm_multi_s *hmmp, struct letter_s *s, struct letter_s *s_2, struct letter_s *s_3,
		  struct letter_s *s_4, struct viterbi_s **ret_viterbi_mtxpp, int use_labels, int multi_scoring_method)
{
  struct viterbi_s *viterbi_mtxp, *cur_rowp;
  struct viterbi_s *prev_rowp; /* pointers to viterbi matrix */
  struct letter_s *seq, *seq_2, *seq_3, *seq_4; /* pointer to the sequence */
  int i; /* loop index */
  int prev; /* nr of previous vertex in viterbi path */
  int seq_len; /* length of the sequence */
  int c, v, w; /* looping indices, c loops over the sequence,
		  v and w over the vertices in the HMM */
  int a_index, a_index_2, a_index_3, a_index_4; /* holds the current letter's position in the alphabet */
  double max, res, t_res1, t_res2, t_res3; /* temporary variables to calculate probabilities */
  double t_res3_1, t_res3_2, t_res3_3, t_res3_4;
  struct path_element *wp, *from_vp, *prevp; /* pointers to elements in to_trans_array and
				     * from_trans_array are used for retrieving the path
				     * from state a to state b */
  int replacement_letter_c, replacement_letter_c_2, replacement_letter_c_3, replacement_letter_c_4;
  int MARKOV_CHAIN; /* set to yes means all emission probs = 1.0, independent of the letter and the state distribution */


  MARKOV_CHAIN = NO;
  /* Allocate memory for matrix + some initial setup:
   * Note 1: viterbi probability matrix has the vertex indices horizontally
   * and the sequence indices vertically meaning it will be filled row by row
   * Note 2: *ret_viterbi_mtxpp is allocated here, but must be
   * freed by caller
   * Note 3: The viterbi algorithm uses the to_trans_array to remember the path
   * taken from a state to another (by letting the prevp point to the path in
   * the to_trans_array that was taken to reach this state) */
  seq_len = get_seq_length(s);
  *ret_viterbi_mtxpp = (struct viterbi_s*)
    (malloc_or_die((seq_len+2) * hmmp->nr_v * sizeof(struct viterbi_s)));
  init_viterbi_s_mtx(*ret_viterbi_mtxpp, DEFAULT, (seq_len+2) * hmmp->nr_v);
  viterbi_mtxp = *ret_viterbi_mtxpp;
  
  /* Convert sequence to 1...L for easier indexing */
  seq = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
  memcpy(seq+1, s, seq_len * sizeof(struct letter_s));
  if(hmmp->nr_alphabets > 1) {
    seq_2 = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
    memcpy(seq_2+1, s_2, seq_len * sizeof(struct letter_s));
  }
  if(hmmp->nr_alphabets > 2) {
    seq_3 = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
    memcpy(seq_3+1, s_3, seq_len * sizeof(struct letter_s));
  }
  if(hmmp->nr_alphabets > 3) {
    seq_4 = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
    memcpy(seq_4+1, s_4, seq_len * sizeof(struct letter_s));
  }

  /* Initialize the first row of the matrix */
  viterbi_mtxp->prob = 0.0; /* sets index (0,0) to 0.0 (i.e. prob 1.0),
			    the rest are already DEFAULT as they should be */

  /* Fill in middle rows */
  prev_rowp = viterbi_mtxp;
  cur_rowp = viterbi_mtxp + hmmp->nr_v;
  
  for(c = 1; c <= seq_len; c++) {
    
    /* get alphabet index for c*/
    if(hmmp->alphabet_type == DISCRETE) {
      replacement_letter_c = NO;
      a_index = get_alphabet_index(&seq[c], hmmp->alphabet, hmmp->a_size);
      if(a_index < 0) {
	a_index = get_replacement_letter_index_multi(&seq[c], hmmp->replacement_letters, 1);
	if(a_index < 0) {
	  printf("Letter '%s' is not in alphabet\n", (&seq[c])->letter);
	  return NOPROB;
	}
	replacement_letter_c = YES;
      }
    }
    if(hmmp->nr_alphabets > 1) {
      if(hmmp->alphabet_type_2 == DISCRETE) {
	replacement_letter_c_2 = NO;
	a_index_2 = get_alphabet_index(&seq_2[c], hmmp->alphabet_2, hmmp->a_size_2);
	if(a_index_2 < 0) {
	  a_index_2 = get_replacement_letter_index_multi(&seq_2[c], hmmp->replacement_letters, 2);
	  if(a_index_2 < 0) {
	    printf("Letter '%s' is not in alphabet\n", (&seq_2[c])->letter);
	    return NOPROB;
	  }
	  replacement_letter_c_2 = YES;
	}
      }
    }
    if(hmmp->nr_alphabets > 2) {
      if(hmmp->alphabet_type_3 == DISCRETE) {
	replacement_letter_c_3 = NO;
	a_index_3 = get_alphabet_index(&seq_3[c], hmmp->alphabet_3, hmmp->a_size_3);
	if(a_index_3 < 0) {
	  a_index_3 = get_replacement_letter_index_multi(&seq_3[c], hmmp->replacement_letters, 3);
	  if(a_index_3 < 0) {
	    printf("Letter '%s' is not in alphabet\n", (&seq_3[c])->letter);
	    return NOPROB;
	  }
	  replacement_letter_c_3 = YES;
	}
      }
    }
    if(hmmp->nr_alphabets > 3) {
      if(hmmp->alphabet_type_4 == DISCRETE) {
	replacement_letter_c_4 = NO;
	a_index_4 = get_alphabet_index(&seq_4[c], hmmp->alphabet_4, hmmp->a_size_4);
	if(a_index_4 < 0) {
	  a_index_4 = get_replacement_letter_index_multi(&seq_4[c], hmmp->replacement_letters, 4);
	  if(a_index_4 < 0) {
	    printf("Letter '%s' is not in alphabet\n", (&seq_4[c])->letter);
	    return NOPROB;
	  }
	  replacement_letter_c_4 = YES;
	}
      }
    }
    
    
    /* calculate sum of probabilities */
    for(v = 1; v < hmmp->nr_v - 1; v++) /* v = to-vertex */ {
#ifdef DEBUG_VI
      printf("prob to vertex %d:\n", v);
#endif
      prev = NO_PREV;
      max = DEFAULT;
      wp = *(hmmp->tot_from_trans_array + v);
      res = 0;
      while((w = wp->vertex) != END) /* w = from-vertex */ {
	/* calculate intermediate results */
	from_vp = wp;
	t_res1 = (prev_rowp + w)->prob; /* probability of having produced the
					 * sequence up to the last letter ending
					 * in vertex w */
	/* calculate intermediate results */
	t_res2 = *(hmmp->max_log_transitions + (w * hmmp->nr_v + v));
	if(t_res1 != DEFAULT && t_res2 != DEFAULT) {
	  res = t_res1 + t_res2;
	}
	else {
	  res = DEFAULT;
	}

	if(prev == NO_PREV && res != DEFAULT &&
	   (use_labels == NO || seq[c].label == *(hmmp->vertex_labels + v) || seq[c].label == '.')) {
	  max = res;
	  prev = from_vp->vertex;
	  prevp = from_vp;
	  continue;
	}
#ifdef DEBUG_VI
	printf("seqlabel = %c, vertexlable = %c\n", seq[c].label,*(hmmp->vertex_labels + v)); 
#endif
	if(res > max && res != DEFAULT &&
	   (use_labels == NO || seq[c].label == *(hmmp->vertex_labels + v) || seq[c].label == '.')) {
	  max = res;
	  prev = from_vp->vertex;
	  prevp = from_vp;
	  }
	wp++;
#ifdef DEBUG_VI
	printf("prev in pos(%d) = %f: ",from_vp->vertex, t_res1);
	printf("trans from %d to %d = %f\n", from_vp->vertex, v, t_res2);
#endif
      }
      
#ifdef DEBUG_VI
      printf("max before setting: %f\n", max);
#endif

      /* add the logprob of reaching state v to the logprob of producing letter l in v*/
      if(MARKOV_CHAIN == YES) {
	if((*((hmmp->log_emissions) + (v * (hmmp->a_size)) + a_index)) != DEFAULT &&
	   max != DEFAULT) {
	  max = max + 0.0;
	}
	else {
	  max = DEFAULT;
	}
      }
      else {
	/* calculate prob of producing letter l in v */
	t_res3 = get_single_emission_score_multi(hmmp, seq, seq_2, seq_3, seq_4, c, v, replacement_letter_c, replacement_letter_c_2,
						 replacement_letter_c_3, replacement_letter_c_4, use_labels, a_index,
						 a_index_2, a_index_3, a_index_4, multi_scoring_method);
      }
      
      if(t_res3 == 0.0) {
	max = DEFAULT;
      }
      else {
	t_res3 = log10(t_res3);
	if(max != DEFAULT) {
	  max = max + t_res3;
	}
      }
     
      /* save result in matrices */
      (cur_rowp + v)->prob = max;
      (cur_rowp + v)->prev = prev;
      (cur_rowp + v)->prevp = prevp;
#ifdef DEBUG_VI
      printf("max after setting: %f\n", max);
#endif
    }
    
    /* move row pointers one row forward */
    prev_rowp = cur_rowp;
    cur_rowp = cur_rowp + hmmp->nr_v;
  }

  /* Fill in transition to end state */
#ifdef DEBUG_VI
  printf("\ntransition to end state:\n");
#endif
  max = DEFAULT;
  prev = NO_PREV;
  prevp = NULL;
  wp = *(hmmp->from_trans_array + hmmp->nr_v-1);
  while(wp->vertex != END) /* w = from-vertex */ {
    from_vp = wp;
    t_res1 = (prev_rowp + wp->vertex)->prob; /* probability of having produced the
					      * sequence up to the last letter ending
					      * in vertex w */
    t_res2 = 0.0;
    while(wp->next != NULL) {
      t_res2 = t_res2 + *((hmmp->log_transitions) +
			  get_mtx_index(wp->vertex, (wp+1)->vertex, hmmp->nr_v));
      wp++;
    }
    t_res2 = t_res2 + *((hmmp->log_transitions) +
			get_mtx_index(wp->vertex, hmmp->nr_v - 1, hmmp->nr_v));
    
#ifdef DEBUG_VI
    printf("prev = %f: ", t_res1);
    printf("trans = %f\n", t_res2);
#endif
    if(t_res1 != DEFAULT && t_res2 != DEFAULT) {
      res = t_res1 + t_res2;
    }
    else {
      res = DEFAULT;
    }
    if(prev == NO_PREV && res != DEFAULT) {
      prev = from_vp->vertex;
      prevp = from_vp;
      max = res;
      continue;
    }
    if(res > max && res != DEFAULT) {
      prev = from_vp->vertex;
      prevp = from_vp;
      max = res;
    }
    wp++;
  }
#ifdef DEBUG_VI
  printf("res = %f\n", max);
#endif
  if(max != DEFAULT) {
    (cur_rowp + hmmp->nr_v-1)->prob = max;
    (cur_rowp + hmmp->nr_v-1)->prev = prev;
    (cur_rowp + hmmp->nr_v-1)->prevp = prevp;
  }
  else {
#ifdef DEBUG_VI
    dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
    dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->log_transitions);
#endif    
    printf("Sequence cannot be produced by this hmm\n"); 
    sequence_as_string(s);
    return NOPROB;
  }

#ifdef PATH
  //dump_viterbi_matrix(seq_len + 2, hmmp->nr_v, viterbi_mtxp);
  //dump_viterbi_path((cur_rowp + w)->prevp);
  printf("normalized log likelihood for most probable path = %f\n",
	 0.0 - (((cur_rowp + hmmp->nr_v - 1)->prob) / seq_len));
  printf("and most probable path is: ");
  //dump_viterbi_path((cur_rowp + hmmp->nr_v - 1), hmmp, viterbi_mtxp, seq_len + 1, hmmp->nr_v);
  printf("%d\n",hmmp->nr_v - 1);
  printf("log prob = %f\n", (cur_rowp + hmmp->nr_v-1)->prob);
  printf("real prob = %f\n", pow(10,(cur_rowp + hmmp->nr_v-1)->prob));
  dump_viterbi_label_path((cur_rowp + hmmp->nr_v - 1), hmmp, viterbi_mtxp, seq_len + 1, hmmp->nr_v);
  printf("\n");
#endif

   /* Garbage collection and return */
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
  return OK;
}


/************************* the one_best algorithm **********************************/
int one_best_multi(struct hmm_multi_s *hmmp, struct letter_s *s, struct letter_s *s_2, struct letter_s *s_3,
		   struct letter_s *s_4, struct one_best_s **ret_one_best_mtxpp,
		   double **ret_scale_fspp, int use_labels, char *best_labeling, int multi_scoring_method)
{
  struct one_best_s *one_best_mtxp, *cur_rowp, *prev_rowp; /* pointers to forward matrix */
  double *scale_fsp; /* pointer to array of scaling factors */
  struct letter_s *seq, *seq_2, *seq_3, *seq_4; /* pointer to the sequence */
  int *sorted_v_list; /* final list for the association of same-labeling-elements */
  int v_list_index;
  int i,j; /* loop variables */
  int seq_len; /* length of the sequence */
  int c, v, w; /* looping indices, c loops over the sequence,
		  v and w over the vertices in the HMM */
  double row_sum, res, t_res1, t_res2, t_res3; /* temporary variables to calculate probabilities */
  double t_res3_1, t_res3_2, t_res3_3, t_res3_4;
  int nr_v;
  struct path_element *wp; /* for traversing the paths from v to w */
  int a_index, a_index_2, a_index_3, a_index_4; /* holds the current letter's position in the alphabet */
  int replacement_letter_c, replacement_letter_c_2, replacement_letter_c_3, replacement_letter_c_4;

  double scaled_result;
  
  /* Allocate memory for matrix and scaling factors + some initial setup:
   * Note 1: one_best probability matrix has the sequence indices vertically
   * and the vertex indices horizontally meaning it will be filled row by row
   * Note 2: *ret_one_best_mtxpp and *ret_scale_fspp are allocated here, but must be
   * freed by caller */
  nr_v = hmmp->nr_v;
  seq_len = get_seq_length(s);
  *ret_one_best_mtxpp = (struct one_best_s*)(malloc_or_die((seq_len+2) *
						       hmmp->nr_v *
						       sizeof(struct one_best_s)));
  one_best_mtxp = *ret_one_best_mtxpp;
  *ret_scale_fspp = (double*)(malloc_or_die((seq_len+2) * sizeof(double)));
  scale_fsp = *ret_scale_fspp;
  sorted_v_list = (int*)(malloc_or_die((hmmp->nr_v * 2 + 1) * sizeof(int)));

  /* Convert sequence to 1...L for easier indexing */
  seq = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
  memcpy(seq+1, s, seq_len * sizeof(struct letter_s));
  if(hmmp->nr_alphabets > 1) {
    seq_2 = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
    memcpy(seq_2+1, s_2, seq_len * sizeof(struct letter_s));
  }
  if(hmmp->nr_alphabets > 2) {
    seq_3 = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
    memcpy(seq_3+1, s_3, seq_len * sizeof(struct letter_s));
  }
  if(hmmp->nr_alphabets > 3) {
    seq_4 = (struct letter_s*) malloc_or_die(((seq_len + 2) * sizeof(struct letter_s)));
    memcpy(seq_4+1, s_4, seq_len * sizeof(struct letter_s));
  }
  
  /* Initialize the first row of the matrix */
  one_best_mtxp->prob = 1.0; /* sets index (0,0) to 1.0,
				the rest are already 0.0 as they should be */ 
  one_best_mtxp->labeling = (char*)(malloc_or_die(1 * sizeof(char)));
  for(i = 1; i < hmmp->nr_v; i++) {
    (one_best_mtxp+i)->labeling = NULL;
  }
  *scale_fsp = 1.0;
  
  /* create initial sorted v-list*/
  *(sorted_v_list) = 0; /* 0 is always the number of the start state */
  *(sorted_v_list + 1) = V_LIST_NEXT;
  *(sorted_v_list + 2) = V_LIST_END;
  

  /* Fill in middle rows */
  prev_rowp = one_best_mtxp;
  cur_rowp = one_best_mtxp + hmmp->nr_v;
  for(c = 1; c <= seq_len; c++) {
    
    /* get alphabet index for c*/
    if(hmmp->alphabet_type == DISCRETE) {
      replacement_letter_c = NO;
      a_index = get_alphabet_index(&seq[c], hmmp->alphabet, hmmp->a_size);
      if(a_index < 0) {
	a_index = get_replacement_letter_index_multi(&seq[c], hmmp->replacement_letters, 1);
	if(a_index < 0) {
	  printf("Letter '%s' is not in alphabet\n", (&seq[c])->letter);
	  return NOPROB;
	}
	replacement_letter_c = YES;
      }
    }
    if(hmmp->nr_alphabets > 1) {
      if(hmmp->alphabet_type_2 == DISCRETE) {
	replacement_letter_c_2 = NO;
	a_index_2 = get_alphabet_index(&seq_2[c], hmmp->alphabet_2, hmmp->a_size_2);
	if(a_index_2 < 0) {
	  a_index_2 = get_replacement_letter_index_multi(&seq_2[c], hmmp->replacement_letters, 2);
	  if(a_index_2 < 0) {
	    printf("Letter '%s' is not in alphabet\n", (&seq_2[c])->letter);
	    return NOPROB;
	  }
	  replacement_letter_c_2 = YES;
	}
      }
    }
    if(hmmp->nr_alphabets > 2) {
      if(hmmp->alphabet_type_3 == DISCRETE) {
	replacement_letter_c_3 = NO;
	a_index_3 = get_alphabet_index(&seq_3[c], hmmp->alphabet_3, hmmp->a_size_3);
	if(a_index_3 < 0) {
	  a_index_3 = get_replacement_letter_index_multi(&seq_3[c], hmmp->replacement_letters, 3);
	  if(a_index_3 < 0) {
	    printf("Letter '%s' is not in alphabet\n", (&seq_3[c])->letter);
	    return NOPROB;
	  }
	  replacement_letter_c_3 = YES;
	}
      }
    }
    if(hmmp->nr_alphabets > 3) {
      if(hmmp->alphabet_type_4 == DISCRETE) {
	replacement_letter_c_4 = NO;
	a_index_4 = get_alphabet_index(&seq_4[c], hmmp->alphabet_4, hmmp->a_size_4);
	if(a_index_4 < 0) {
	  a_index_4 = get_replacement_letter_index_multi(&seq_4[c], hmmp->replacement_letters, 4);
	  if(a_index_4 < 0) {
	    printf("Letter '%s' is not in alphabet\n", (&seq_4[c])->letter);
	    return NOPROB;
	  }
	  replacement_letter_c_4 = YES;
	}
      }
    }
    
#ifdef DEBUG_FW  
    printf("seq[c] = %s\n", &seq[c]);
    printf("a_index = %d\n", a_index);
    printf("v_list dump out ");
    dump_v_list(sorted_v_list);
   
#endif

    /* calculate probabilities */
    row_sum = 0.0;
    for(v = 0; v < hmmp->nr_v - 1; v++) /* v = to-vertex */ {
#ifdef DEBUG_FW
      printf("prob to vertex %d\n", v);
#endif
      v_list_index = 0;
      (cur_rowp + v)->labeling = NULL;
      
      while(*(sorted_v_list + v_list_index) != V_LIST_END) {
	res = 0.0;
	while(*(sorted_v_list + v_list_index) != V_LIST_NEXT) {
	  w = *(sorted_v_list + v_list_index); /* w = from-vertex */
	  /* calculate intermediate results */
	  res += (prev_rowp + w)->prob * *(hmmp->tot_transitions + (w * nr_v + v));
	  if(*(hmmp->tot_transitions + (w * nr_v + v)) < 0) {
	    printf("Error: found model transition prob from %d to %d < 0.0\n", w, v);
	    exit(0);
	  }
#ifdef DEBUG_FW
	  printf("prev = %f: ", (prev_rowp + w)->prob);
	  printf("trans = %f\n", *(hmmp->tot_transitions + (w * nr_v + v)));
#endif
	  v_list_index++;
	}
	
	if(res == 0.0) {
	  v_list_index++;
	  continue;
	}

	/* calculate prob of producing letter l in v */
	t_res3 = get_single_emission_score_multi(hmmp, seq, seq_2, seq_3, seq_4, c, v, replacement_letter_c, replacement_letter_c_2,
						 replacement_letter_c_3, replacement_letter_c_4, use_labels, a_index,
						 a_index_2, a_index_3, a_index_4, multi_scoring_method);
	
	
	res = res * t_res3;
	
	
	/* check if this score is best so far */
#ifdef DEBUG_FW
	printf("best score = %f\n",(cur_rowp + v)->prob);
#endif
	if(res > (cur_rowp + v)->prob) {
#ifdef DEBUG_FW	  
	  printf("updating best score to %f\n", res);  
#endif
	  /* save result in matrices */
	  (cur_rowp + v)->prob = res;
	  /* set pointer to point to current labeling */
	  (cur_rowp + v)->labeling = (prev_rowp + w)->labeling;
	  if((cur_rowp + v)->labeling == NULL) {
	    printf("Error: NULL labeling when updating best score\n");
	    exit(0);
	  }
	}
	
#ifdef DEBUG_FW
	printf("res = %f\n", res);
#endif
	v_list_index++;
      }
    }
    
    /* update labeling pointers */
    update_labelings(cur_rowp, hmmp->vertex_labels, sorted_v_list, seq_len, c, hmmp->labels, hmmp->nr_labels, hmmp->nr_v);
    deallocate_row_labelings(prev_rowp, hmmp->nr_v);
    
    /* scale the results, row_sum = the total probability of 
     * having produced the labelings up to and including character c */
    row_sum = 0.0;
    for(v = 0; v < hmmp->nr_v; v++) {
      row_sum = row_sum + (cur_rowp + v)->prob;
#ifdef DEBUG_FW
      dump_labeling((cur_rowp + v)->labeling, c);
      printf("c = %d\n", c);
#endif
    }
    scale_fsp[c] = row_sum;
#ifdef DEBUG_FW
    printf("rowsum = %f\n",row_sum);
    printf("scaling set to: %f\n", scale_fsp[c]);
#endif
    if(row_sum == 0.0) {
      printf("Sequence cannot be produced by this hmm\n"); 
      sequence_as_string(s);
      exit(0);
      return NOPROB;
    }
    for(v = 0; v < hmmp->nr_v; v++) {
      if((cur_rowp + v)->prob != 0.0){
	(cur_rowp + v)->prob = ((cur_rowp + v)->prob)/row_sum; /* scaling */
      }
    }
    
    /* move row pointers one row forward */
    prev_rowp = cur_rowp;
    cur_rowp = cur_rowp + hmmp->nr_v;
  }

  /* Fill in transition to end state */
  v_list_index = 0;
  (cur_rowp + hmmp->nr_v - 1)->labeling = NULL;
  while(*(sorted_v_list + v_list_index) != V_LIST_END) {
    res = 0.0;
    while(*(sorted_v_list + v_list_index) != V_LIST_NEXT) {
      w = *(sorted_v_list + v_list_index); /* w = from-vertex */
      t_res1 = (prev_rowp + w)->prob;
      t_res2 = *((hmmp->tot_transitions) + get_mtx_index(w, hmmp->nr_v-1, hmmp->nr_v));
      if(t_res2 > 1.0) {
	t_res2 = 1.0;
      }
      res += t_res1 * t_res2;
      v_list_index++;
    }
    
    /* check if this score is best so far */
    if(res > (cur_rowp + hmmp->nr_v - 1)->prob) {
      /* save result in matrices */
      (cur_rowp + hmmp->nr_v - 1)->prob = res;
      /* set pointer to point to current labeling */
      (cur_rowp + hmmp->nr_v - 1)->labeling = (prev_rowp + w)->labeling;
    }
    
    v_list_index++;
  }

#ifdef DEBUG_FW
  dump_one_best_matrix(seq_len + 2, hmmp->nr_v, one_best_mtxp);
  dump_scaling_array(seq_len + 1, scale_fsp);
#endif
  
  /* store results */
  scaled_result = (cur_rowp + hmmp->nr_v - 1)->prob;
  memcpy(best_labeling, ((cur_rowp + hmmp->nr_v - 1)->labeling) + 1, seq_len * sizeof(char));
  best_labeling[seq_len] = '\0';
#ifdef DEBUG_PATH
  printf("seq_len = %d\n", seq_len);
  printf("best labeling = %s\n", best_labeling);
#endif
  
  /* Garbage collection and return */
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
  free(sorted_v_list);
  /* FREE labelings, except result labeling */
  deallocate_row_labelings(prev_rowp, hmmp->nr_v);
  
  return OK;
}





/*************************************************************************************
 ************************** the forward, backward and viterbi algorithms *************
 ************************** for scoring an msa against the hmm ***********************
 *************************************************************************************/


/************************* the msa_forward algorithm **********************************/
int msa_forward_multi(struct hmm_multi_s *hmmp, struct msa_sequences_multi_s *msa_seq_infop, int use_lead_columns,
		      int use_gap_shares, int use_prior_shares, struct forward_s **ret_forw_mtxpp, double **ret_scale_fspp,
		      int use_labels, int normalize, int scoring_method, int multi_scoring_method, double *aa_freqs,
		      double *aa_freqs_2, double *aa_freqs_3, double *aa_freqs_4)
{
  struct forward_s *forw_mtxp, *cur_rowp, *prev_rowp; /* pointers to forward matrix */
  double *scale_fsp; /* pointer to array of scaling factors */
  int seq_len; /* length of the sequence */
  int c, v, w; /* looping indices, c loops over the sequence,
		  v and w over the vertices in the HMM */
  int i,j; /* general loop indices */
  int a_index, a_index_2, a_index_3, a_index_4; /* holds the current letter's position in the alphabet */
  double row_sum, res, t_res1, t_res2, t_res3; /* temporary variables to calculate probabilities */
  struct path_element *wp; /* for traversing the paths from v to w */
  int nr_v;

#ifdef DEBUG_FW
  printf("entering msa_forward\n");
#endif

  /* Allocate memory for matrix and scaling factors + some initial setup:
   * Note 1: forward probability matrix has the sequence indices vertically
   * and the vertex indices horizontally meaning it will be filled row by row
   * Note 2: *ret_forw_mtxpp and *ret_scale_fspp are allocated here, but must be
   * freed by caller */
  nr_v = hmmp->nr_v;
  if(use_lead_columns == YES) {
    seq_len = msa_seq_infop->nr_lead_columns;
  }
  else {
    seq_len = msa_seq_infop->msa_seq_length;
  }
  *ret_forw_mtxpp = (struct forward_s*)(malloc_or_die((seq_len+2) *
						      hmmp->nr_v *
						      sizeof(struct forward_s)));
  forw_mtxp = *ret_forw_mtxpp;
  *ret_scale_fspp = (double*)(malloc_or_die((seq_len+2) * sizeof(double)));
  scale_fsp = *ret_scale_fspp;

  /* Initialize the first row of the matrix */
  forw_mtxp->prob = 1.0; /* sets index (0,0) to 1.0,
			    the rest are already 0.0 as they should be */
  *scale_fsp = 1.0;
  
  /* Fill in middle rows */
  prev_rowp = forw_mtxp;
  cur_rowp = forw_mtxp + hmmp->nr_v;
  j = 0;
  if(use_lead_columns == YES) {
    c = *(msa_seq_infop->lead_columns_start + j);
  }
  else {
    c = 0;
  }
  while(c != END && c < msa_seq_infop->msa_seq_length)  {
    a_index_2 = -1;
    a_index_3 = -1;
    a_index_4 = -1;
    
    if(hmmp->alphabet_type == DISCRETE) {
      
      a_index = get_alphabet_index((msa_seq_infop->msa_seq_1 + (c * (hmmp->a_size+1)))->query_letter, hmmp->alphabet, hmmp->a_size);
      if(a_index < 0) {
	a_index = hmmp->a_size; /* if letter is wild card, use default column in subst matrix */
      }
    }
    if(hmmp->nr_alphabets > 1) {
      if(hmmp->alphabet_type_2 == DISCRETE) {
	a_index_2 = get_alphabet_index((msa_seq_infop->msa_seq_2 + (c * (hmmp->a_size_2+1)))->query_letter,
				       hmmp->alphabet_2, hmmp->a_size_2);
	if(a_index_2 < 0) {
	  a_index_2 = hmmp->a_size_2; /* if letter is wild card, use default column in subst matrix */
	}
      }
    }
    if(hmmp->nr_alphabets > 2) {
      if(hmmp->alphabet_type_3 == DISCRETE) {
	a_index_3 = get_alphabet_index((msa_seq_infop->msa_seq_3 + (c * (hmmp->a_size_3+1)))->query_letter,
				       hmmp->alphabet_3, hmmp->a_size_3);
	if(a_index_3 < 0) {
	  a_index_3 = hmmp->a_size_3; /* if letter is wild card, use default column in subst matrix */
	}
      }
    }
    if(hmmp->nr_alphabets > 3) {
      if(hmmp->alphabet_type_4 == DISCRETE) {
	a_index_4 = get_alphabet_index((msa_seq_infop->msa_seq_4 + (c * (hmmp->a_size_4+1)))->query_letter,
				       hmmp->alphabet_4, hmmp->a_size_4);
	if(a_index_4 < 0) {
	  a_index_4 = hmmp->a_size_4; /* if letter is wild card, use default column in subst matrix */
	}
      }
    }
    /* calculate sum of probabilities */
    row_sum = 0;
    for(v = 1; v < hmmp->nr_v - 1; v++) /* v = to-vertex */ {
#ifdef DEBUG_FW
      printf("prob to vertex %d:\n", v);
#endif
      res = 0;
      wp = *(hmmp->tot_from_trans_array + v);
      while((w = wp->vertex) != END) /* w = from-vertex */ {
	/* calculate intermediate results */
	res += (prev_rowp + w)->prob * *((hmmp->tot_transitions) + (w * nr_v + v));
	if(*((hmmp->tot_transitions) + (w * nr_v + v)) < 0) {
	  printf("found model transition prob from %d to %d < 0.0\n", w, v);
	  exit(0);
	}
	wp++;
#ifdef DEBUG_FW
	printf("prev = %f: ",(prev_rowp + w)->prob );
	printf("trans = %f\n",*((hmmp->tot_transitions) + (w * nr_v + v)) );
#endif
      }
      
      /* calculate the prob of producing letters l in v*/
      t_res3 = get_msa_emission_score_multi(msa_seq_infop, hmmp, c, v, use_labels, normalize, a_index, a_index_2,
					    a_index_3, a_index_4, scoring_method, use_gap_shares,
					    use_prior_shares, multi_scoring_method, aa_freqs, aa_freqs_2,
					    aa_freqs_3, aa_freqs_4);
      
      res = res * t_res3;
      row_sum += res;

      /* save result in matrices */
      (cur_rowp + v)->prob = res;
#ifdef DEBUG_FW
      printf("letter = %f\n", t_res3);
      printf("res = %f\n", res);
#endif
    }
    
    /* scale the results, row_sum = the total probability of 
     * having produced the sequence up to and including character c */
     if(use_lead_columns == YES) {
      scale_fsp[j+1] = row_sum;
    }
    else {
      scale_fsp[c+1] = row_sum;
    }
    
#ifdef DEBUG_FW
    printf("rowsum for row %d = %f\n", c, row_sum);
    printf("scaling set to: %f\n", scale_fsp[c+1]);
#endif
    if(row_sum == 0.0) {
      printf("Probability for this msa = 0.0, pos = %d\n", c);
      printf("in forward\n");
      exit(0);
      return NOPROB;
    }
    for(v = 0; v < hmmp->nr_v; v++) {
      if((cur_rowp + v)->prob != 0){
	(cur_rowp + v)->prob = ((cur_rowp + v)->prob)/row_sum; /* scaling */
      }
    }
    
    /* move row pointers one row forward */
    prev_rowp = cur_rowp;
    cur_rowp = cur_rowp + hmmp->nr_v;

    /* update current column */
    if(use_lead_columns == YES) {
      j++;
      c = *(msa_seq_infop->lead_columns_start + j);
    }
    else {
      c++;
    }
  }


  /* Fill in transition to end state */
#ifdef DEBUG_FW
  printf("\ntransition to end state:\n");
#endif
  res = 0;
  wp = *(hmmp->tot_from_trans_array + hmmp->nr_v - 1);
  while(wp->vertex != END) {
    t_res1 = (prev_rowp + wp->vertex)->prob;
    t_res2 = *((hmmp->tot_transitions) + get_mtx_index(wp->vertex, hmmp->nr_v - 1, hmmp->nr_v));
    if(t_res2 > 1.0) {
      t_res2 = 1.0;
    }
    res += t_res1 * t_res2;
   
    wp++;
  }
#ifdef DEBUG_FW
  printf("res = %f\n", res);
#endif
  (cur_rowp + hmmp->nr_v - 1)->prob = res; /* obs: no scaling performed here */

#ifdef DEBUG_FW
  dump_forward_matrix(seq_len + 2, hmmp->nr_v, forw_mtxp);
  dump_scaling_array(seq_len + 1, scale_fsp);
  printf("exiting msa_forward\n\n\n");
#endif

  /* Garbage collection and return */
  return OK;
}




/**************************** the msa_backward algorithm **********************************/
int msa_backward_multi(struct hmm_multi_s *hmmp, struct msa_sequences_multi_s *msa_seq_infop, int use_lead_columns,
		       int use_gap_shares, struct backward_s **ret_backw_mtxpp, double *scale_fsp, int use_labels, int normalize,
		       int scoring_method, int multi_scoring_method, double *aa_freqs, double *aa_freqs_2, double *aa_freqs_3,
		       double *aa_freqs_4)
{
  struct backward_s *backw_mtxp, *cur_rowp, *prev_rowp; /* pointers to backward matrix */
  int seq_len; /* length of the sequence */
  int c, v, w; /* looping indices, c loops over the sequence,
		  v and w over the vertices in the HMM */
  int i,j; /* general loop indices */
  int all_columns; /* boolean for checking whether all columns have been looped over */
  int a_index, a_index_2, a_index_3, a_index_4; /* holds the current letter's position in the alphabet */
  double row_sum, res,t_res1,t_res2, t_res3; /* temporary variables to calculate results */
  struct path_element *wp;
  
  /* Allocate memory for matrix + some initial setup:
   * Note 1: probability matrix has the sequence indices vertically
   * and the vertex indices horizontally meaning it will be filled row by row
   * Note 2: *ret_backw_mtxpp is allocated here, but must be
   * freed by caller
   * Note 3: *scale_fspp is the scaling array produced by forward() meaning backward()
   * must be called after forward() */

  if(use_lead_columns == YES) {
    seq_len = msa_seq_infop->nr_lead_columns;
  }
  else {
    seq_len = msa_seq_infop->msa_seq_length;
  }
  *ret_backw_mtxpp = (struct backward_s*)(malloc_or_die((seq_len+2) *
							hmmp->nr_v * 
							sizeof(struct backward_s)));
  backw_mtxp = *ret_backw_mtxpp;
  
  /* Initialize the last row of the matrix */
  (backw_mtxp + get_mtx_index(seq_len + 1, hmmp->nr_v - 1,
			      hmmp->nr_v))->prob = 1.0; /* sets index
							 * (seq_len+1,nr_v-1) i.e.
							 * the lower right
							 * corner of the matrix to 
							 * 1.0,the rest are already
							 * 0.0 as they should be*/
 
  /* Fill in next to last row in matrix (i.e. add prob for path to end state for all states
   * that have a transition path to the end state) */
  prev_rowp = backw_mtxp + get_mtx_index(seq_len + 1, 0, hmmp->nr_v);
  cur_rowp = prev_rowp - hmmp->nr_v;
  wp = *(hmmp->from_trans_array + hmmp->nr_v - 1);
  while(wp->vertex != END) {
    w = wp->vertex;
    t_res1 = 1.0;
    while(wp->next != NULL) {
      t_res1 = t_res1 * *(hmmp->transitions +
			  get_mtx_index(wp->vertex, (wp + 1)->vertex, hmmp->nr_v));
      wp++;
    }
    (cur_rowp + w)->prob += t_res1;
    if((cur_rowp + w)->prob > 1.0) {
      (cur_rowp + w)->prob = 1.0;
    }
    wp++;
  }
  prev_rowp = cur_rowp;
  cur_rowp = cur_rowp - hmmp->nr_v;


  /* Fill in first rows moving upwards in matrix */
  all_columns = NO;
  j = 1;
  if(use_lead_columns == YES) {
    c = *(msa_seq_infop->lead_columns_end - j);
  }
  else {
    c = seq_len - 1;
  }

  while(c >= 0 && all_columns == NO) {
    a_index_2 = -1;
    a_index_3 = -1;
    a_index_4 = -1;
    if(hmmp->alphabet_type == DISCRETE) {
      a_index = get_alphabet_index((msa_seq_infop->msa_seq_1 + (c * (hmmp->a_size+1)))->query_letter, hmmp->alphabet, hmmp->a_size);
      if(a_index < 0) {
	a_index = hmmp->a_size; /* if letter is wild card, use default column in subst matrix */
      }
    }
    if(hmmp->nr_alphabets > 1) {
      if(hmmp->alphabet_type_2 == DISCRETE) {
	a_index_2 = get_alphabet_index((msa_seq_infop->msa_seq_2 + (c * (hmmp->a_size_2+1)))->query_letter,
				       hmmp->alphabet_2, hmmp->a_size_2);
	if(a_index_2 < 0) {
	  a_index_2 = hmmp->a_size_2; /* if letter is wild card, use default column in subst matrix */
	}
      }
    }
    if(hmmp->nr_alphabets > 2) {
      if(hmmp->alphabet_type_3 == DISCRETE) {
	a_index_3 = get_alphabet_index((msa_seq_infop->msa_seq_3 + (c * (hmmp->a_size_3+1)))->query_letter,
				       hmmp->alphabet_3, hmmp->a_size_3);
	if(a_index_3 < 0) {
	  a_index_3 = hmmp->a_size_3; /* if letter is wild card, use default column in subst matrix */
	}
      }
    }
    if(hmmp->nr_alphabets > 3) {
      if(hmmp->alphabet_type_4 == DISCRETE) {
	a_index_4 = get_alphabet_index((msa_seq_infop->msa_seq_4 + (c * (hmmp->a_size_4+1)))->query_letter,
				       hmmp->alphabet_4, hmmp->a_size_4);
	if(a_index_4 < 0) {
	  a_index_4 = hmmp->a_size_4; /* if letter is wild card, use default column in subst matrix */
	}
      }
    }
    /* calculate sum of probabilities */
    for(v = 0; v < hmmp->nr_v - 1; v++) /* v = from-vertex */{
#ifdef DEBUG_BW
      printf("prob passing through vertex %d:\n", v);
#endif
      
      res = 0;
      wp = *(hmmp->tot_to_trans_array + v);
      while(wp->vertex != END) 	/* w = to-vertex */ {
	/* probability of transiting from v to w on this particular path */
	t_res2 = *((hmmp->tot_transitions) + get_mtx_index(v, wp->vertex, hmmp->nr_v));
	if(t_res2 < 0) {
	  printf("found model transition prob from %d to %d < 0.0\n", v, wp->vertex);
	  exit(0);
	}
	t_res1 = (prev_rowp + wp->vertex)->prob; /* probability of having produced the
						  * sequence after c passing through
						  * vertex w */
	
	/* calculate the prob of producing letters l in v*/
	t_res3 = get_msa_emission_score_multi(msa_seq_infop, hmmp, c, wp->vertex, use_labels, normalize, a_index, a_index_2,
					      a_index_3, a_index_4, scoring_method, use_gap_shares,
					      NO, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);	
	
	res += t_res1 * t_res2 * t_res3; /* multiply the probabilities and
					  * add to previous sum */
	wp++;
      }
#ifdef DEBUG_BW
      printf("prev = %f: ", t_res1);
      printf("trans = %f\n", t_res2);
#endif
     
      

      /* save result in matrices */
      (cur_rowp + v)->prob = res;
#ifdef DEBUG_BW
      printf("res = %f\n", res);
#endif
    }
    
    /* scale the results using the scaling factors produced by forward() */
#ifdef DEBUG_BW
    printf("c = %d\n", c);
    dump_scaling_array(seq_len + 1, scale_fsp);
#endif
    for(v = 0; v < hmmp->nr_v; v++) {
      if((cur_rowp + v)->prob != 0){
	if(use_lead_columns == YES) {
	  if(scale_fsp[msa_seq_infop->nr_lead_columns + 1 - j] != 0.0) { 
	    (cur_rowp + v)->prob =
	      ((cur_rowp + v)->prob)/scale_fsp[msa_seq_infop->nr_lead_columns + 1 - j]; /* scaling */
	  }
	  else {
	    printf("This msa cannot be produced by hmm\n");
	    return NOPROB;
	  }
	}
	else {
	  if(scale_fsp[c+1] != 0.0) { 
	    (cur_rowp + v)->prob = ((cur_rowp + v)->prob)/scale_fsp[c+1]; /* scaling */
	  }
	  else {
	    printf("This msa cannot be produced by hmm\n");
	    return NOPROB;
	  }
	}
      }
    }
    
    /* move row pointers one row backward */
    prev_rowp = cur_rowp;
    cur_rowp = cur_rowp - hmmp->nr_v;


    /* update current column */
    if(use_lead_columns == YES) {
      if(c == *(msa_seq_infop->lead_columns_start)) {
	all_columns = YES;
      }
      j++;
      c = *(msa_seq_infop->lead_columns_end - j);
    }
    else {
      c = c - 1;
    }
  }

#ifdef DEBUG_BW
  dump_backward_matrix(seq_len + 2, hmmp->nr_v, backw_mtxp);
  dump_scaling_array(seq_len + 1, scale_fsp);
#endif

  /* Garbage collection and return */
  return OK;
}


/************************* the msa_viterbi algorithm **********************************/
int msa_viterbi_multi(struct hmm_multi_s *hmmp, struct msa_sequences_multi_s *msa_seq_infop, int use_lead_columns,
		      int use_gap_shares, int use_prior_shares, struct viterbi_s **ret_viterbi_mtxpp, int use_labels, int normalize,
		      int scoring_method, int multi_scoring_method, double *aa_freqs, double *aa_freqs_2, double *aa_freqs_3,
		      double *aa_freqs_4)
{
  struct viterbi_s *viterbi_mtxp, *cur_rowp, *prev_rowp; /* pointers to viterbi matrix */
  int seq_len; /* length of the sequence */
  int c, v, w; /* looping indices, c loops over the sequence,
		  v and w over the vertices in the HMM */
  int i,j; /* general loop indices */
  int a_index, a_index_2, a_index_3, a_index_4; /* holds the current letter's position in the alphabet */
  double max, res, t_res1, t_res2, t_res3; /* temporary variables to calculate probabilities */
  struct path_element *wp, *from_vp, *prevp; /* for traversing the paths from v to w */
  int nr_v;
  int prev;
  double seq_normalizer, state_normalizer;
  
  /* Allocate memory for matrix  + some initial setup:
   * Note 1: viterbi probability matrix has the sequence indices vertically
   * and the vertex indices horizontally meaning it will be filled row by row
   * Note 2: *ret_viterbi_mtxpp and *ret_scale_fspp are allocated here, but must be
   * freed by caller */
  nr_v = hmmp->nr_v;
  if(use_lead_columns == YES) {
    seq_len = msa_seq_infop->nr_lead_columns;
  }
  else {
    seq_len = msa_seq_infop->msa_seq_length;
  }
  *ret_viterbi_mtxpp = (struct viterbi_s*)(malloc_or_die((seq_len+2) *
							 hmmp->nr_v *
							 sizeof(struct viterbi_s)));
  init_viterbi_s_mtx(*ret_viterbi_mtxpp, DEFAULT, (seq_len+2) * hmmp->nr_v);
  viterbi_mtxp = *ret_viterbi_mtxpp;
  
  /* Initialize the first row of the matrix */
  viterbi_mtxp->prob = 0.0; /* sets index (0,0) to 0.0 (i.e. prob 1.0),
			       the rest are already DEFAULT as they should be */
  
  /* Fill in middle rows */
  prev_rowp = viterbi_mtxp;
  cur_rowp = viterbi_mtxp + hmmp->nr_v;
  j = 0;
  if(use_lead_columns == YES) {
    c = *(msa_seq_infop->lead_columns_start + j);
  }
  else {
    c = 0;
  }

 

  while(c != END && c < msa_seq_infop->msa_seq_length)  {
    a_index_2 = -1;
    a_index_3 = -1;
    a_index_4 = -1;
    if(hmmp->alphabet_type == DISCRETE) {
      a_index = get_alphabet_index((msa_seq_infop->msa_seq_1 + (c * (hmmp->a_size+1)))->query_letter, hmmp->alphabet, hmmp->a_size);
      if(a_index < 0) {
	a_index = hmmp->a_size; /* if letter is wild card, use default column in subst matrix */
      }
    }
    if(hmmp->nr_alphabets > 1) {
      if(hmmp->alphabet_type_2 == DISCRETE) {
	a_index_2 = get_alphabet_index((msa_seq_infop->msa_seq_2 + (c * (hmmp->a_size_2+1)))->query_letter,
				       hmmp->alphabet_2, hmmp->a_size_2);
	if(a_index_2 < 0) {
	  a_index_2 = hmmp->a_size_2; /* if letter is wild card, use default column in subst matrix */
	}
      }
    }
    if(hmmp->nr_alphabets > 2) {
      if(hmmp->alphabet_type_3 == DISCRETE) {
	a_index_3 = get_alphabet_index((msa_seq_infop->msa_seq_3 + (c * (hmmp->a_size_3+1)))->query_letter,
				       hmmp->alphabet_3, hmmp->a_size_3);
	if(a_index_3 < 0) {
	  a_index_3 = hmmp->a_size_3; /* if letter is wild card, use default column in subst matrix */
	}
      }
    }
    if(hmmp->nr_alphabets > 3) {
      if(hmmp->alphabet_type_4 == DISCRETE) {
	a_index_4 = get_alphabet_index((msa_seq_infop->msa_seq_4 + (c * (hmmp->a_size_4+1)))->query_letter,
				       hmmp->alphabet_4, hmmp->a_size_4);
	if(a_index_4 < 0) {
	  a_index_4 = hmmp->a_size_4; /* if letter is wild card, use default column in subst matrix */
	}
      }
    }
    
    /* calculate sum of probabilities */
#ifdef DEBUG_VI
    printf("label nr %d for pos %d = %c\n", j, c, (msa_seq_infop->msa_seq_1 + (c * (hmmp->a_size+1)))->label);
#endif    
    for(v = 1; v < hmmp->nr_v - 1; v++) /* v = to-vertex */ {
#ifdef DEBUG_VI
      printf("prob to vertex %d:\n", v);
#endif
      prev = NO_PREV;
      max = DEFAULT;
      wp = *(hmmp->tot_from_trans_array + v);
      res = 0;
      while((w = wp->vertex) != END) /* w = from-vertex */ {
	/* calculate intermediate results */
	from_vp = wp;
	t_res1 = (prev_rowp + w)->prob; /* probability of having produced the
					 * sequence up to the last letter ending
					 * in vertex w */
	/* calculate intermediate results */
	t_res2 = *(hmmp->max_log_transitions + (w * hmmp->nr_v + v));
	if(t_res1 != DEFAULT && t_res2 != DEFAULT) {
	  res = t_res1 + t_res2;
	}
	else {
	  res = DEFAULT;
	}
	

	if(prev == NO_PREV && res != DEFAULT &&
	   (use_labels == NO || (msa_seq_infop->msa_seq_1 + (c * (hmmp->a_size+1)))->label == *(hmmp->vertex_labels + v) ||
	    (msa_seq_infop->msa_seq_1 + (c * (hmmp->a_size+1)))->label == '.')) {
	  max = res;
	  prev = from_vp->vertex;
	  prevp = from_vp;
	  continue;
	}
	  
	if(res > max && res != DEFAULT &&
	   (use_labels == NO || (msa_seq_infop->msa_seq_1 + (c * (hmmp->a_size+1)))->label == *(hmmp->vertex_labels + v) ||
	    (msa_seq_infop->msa_seq_1 + (c * (hmmp->a_size+1)))->label == '.')) {
	  max = res;
	  prev = from_vp->vertex;
	  prevp = from_vp;
	  }
	wp++;
#ifdef DEBUG_VI
	printf("prev in pos(%d) = %f: ",from_vp->vertex, t_res1);
	printf("trans from %d to %d = %f\n", from_vp->vertex, v, t_res2);
	printf("tot score = %f\n", res);
#endif
      }
      
#ifdef DEBUG_VI
      printf("max before setting: %f\n", max);
#endif
      
      /* calculate the prob of producing letters l in v*/
      t_res3 = get_msa_emission_score_multi(msa_seq_infop, hmmp, c, v, use_labels, normalize, a_index, a_index_2,
					    a_index_3, a_index_4, scoring_method, use_gap_shares,
					    use_prior_shares, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3,
					    aa_freqs_4);

      if(t_res3 == 0.0) {
	max = DEFAULT;
      }
      else {
	t_res3 = log10(t_res3);
	if(max != DEFAULT) {
	  max = max + t_res3;
	}
      }
    
      /* save result in matrices */
      (cur_rowp + v)->prob = max;
      (cur_rowp + v)->prev = prev;
      (cur_rowp + v)->prevp = prevp;
#ifdef DEBUG_VI
      printf("max after setting: %f\n", max);
#endif
    } 
    
    /* move row pointers one row viterbi */
    prev_rowp = cur_rowp;
    cur_rowp = cur_rowp + hmmp->nr_v;

    /* update current column */
    if(use_lead_columns == YES) {
      j++;
      c = *(msa_seq_infop->lead_columns_start + j);
    }
    else {
      c++;
    }
  }
  

  /* Fill in transition to end state */
#ifdef DEBUG_VI
  printf("\ntransition to end state:\n");
#endif
  max = DEFAULT;
  prev = NO_PREV;
  prevp = NULL;
  wp = *(hmmp->from_trans_array + hmmp->nr_v-1);
  while(wp->vertex != END) /* w = from-vertex */ {
    from_vp = wp;
    t_res1 = (prev_rowp + wp->vertex)->prob; /* probability of having produced the
					      * sequence up to the last letter ending
					      * in vertex w */
    t_res2 = 0.0;
    while(wp->next != NULL) {
      t_res2 = t_res2 + *((hmmp->log_transitions) +
			  get_mtx_index(wp->vertex, (wp+1)->vertex, hmmp->nr_v));
      wp++;
    }
    t_res2 = t_res2 + *((hmmp->log_transitions) +
			get_mtx_index(wp->vertex, hmmp->nr_v - 1, hmmp->nr_v));
    
#ifdef DEBUG_VI
    printf("prev = %f: ", t_res1);
    printf("trans = %f\n", t_res2);
#endif
    if(t_res1 != DEFAULT && t_res2 != DEFAULT) {
      res = t_res1 + t_res2;
    }
    else {
      res = DEFAULT;
    }
    if(prev == NO_PREV && res != DEFAULT) {
      prev = from_vp->vertex;
      prevp = from_vp;
      max = res;
      continue;
    }
    if(res > max && res != DEFAULT) {
      prev = from_vp->vertex;
      prevp = from_vp;
      max = res;
    }
    wp++;
  }
#ifdef DEBUG_VI
  printf("res = %f\n", max);
#endif
  if(max != DEFAULT) {
    (cur_rowp + hmmp->nr_v-1)->prob = max;
    (cur_rowp + hmmp->nr_v-1)->prev = prev;
    (cur_rowp + hmmp->nr_v-1)->prevp = prevp;
  }
  else {
#ifdef DEBUG_VI
    dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
    dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->log_transitions);
#endif    
    printf("This msa cannot be produced by hmm\n");
    printf("inne\n");
    return NOPROB;
  }
  
#ifdef PATH
  //dump_viterbi_matrix(seq_len + 2, hmmp->nr_v, viterbi_mtxp);
  //dump_viterbi_path((cur_rowp + w)->prevp);
  printf("normalized log likelihood for most probable path = %f\n",
	 0.0 - (((cur_rowp + hmmp->nr_v - 1)->prob) / seq_len));
  printf("and most probable path is: ");
  dump_viterbi_path((cur_rowp + hmmp->nr_v - 1), hmmp, viterbi_mtxp, seq_len + 1, hmmp->nr_v);
  printf("%d\n",hmmp->nr_v - 1);
  printf("log prob = %f\n", (cur_rowp + hmmp->nr_v-1)->prob);
  printf("real prob = %f\n", pow(10,(cur_rowp + hmmp->nr_v-1)->prob));
  dump_viterbi_label_path((cur_rowp + hmmp->nr_v - 1), hmmp, viterbi_mtxp, seq_len + 1, hmmp->nr_v);
  printf("\n");
#endif

  /* Garbage collection and return */
  return OK;
}


/************************* the msa_one_best algorithm **********************************/
int msa_one_best_multi(struct hmm_multi_s *hmmp, struct msa_sequences_multi_s *msa_seq_infop, int use_lead_columns,
		       int use_gap_shares, int use_prior_shares, struct one_best_s **ret_one_best_mtxpp, double **ret_scale_fspp,
		       int use_labels, char *best_labeling, int normalize, int scoring_method, int multi_scoring_method,
		       double *aa_freqs, double *aa_freqs_2, double *aa_freqs_3, double *aa_freqs_4)
{
  struct one_best_s *one_best_mtxp, *cur_rowp, *prev_rowp; /* pointers to one_best matrix */
  double *scale_fsp; /* pointer to array of scaling factors */
  int seq_len; /* length of the sequence */
  int c, v, w; /* looping indices, c loops over the sequence,
		  v and w over the vertices in the HMM */
  int i,j; /* general loop indices */
  int a_index, a_index_2, a_index_3, a_index_4; /* holds the current letter's position in the alphabet */
  double row_sum, res, t_res1, t_res2, t_res3; /* temporary variables to calculate probabilities */
  struct path_element *wp; /* for traversing the paths from v to w */
  int nr_v;
  double scaled_result;

  int *sorted_v_list;
  int v_list_index;
  
  /* Allocate memory for matrix and scaling factors + some initial setup:
   * Note 1: one_best probability matrix has the sequence indices vertically
   * and the vertex indices horizontally meaning it will be filled row by row
   * Note 2: *ret_one_best_mtxpp and *ret_scale_fspp are allocated here, but must be
   * freed by caller */
  nr_v = hmmp->nr_v;
  if(use_lead_columns == YES) {
    seq_len = msa_seq_infop->nr_lead_columns;
  }
  else {
    seq_len = msa_seq_infop->msa_seq_length;
  }
  *ret_one_best_mtxpp = (struct one_best_s*)(malloc_or_die((seq_len+2) *
						      hmmp->nr_v *
						      sizeof(struct one_best_s)));
  one_best_mtxp = *ret_one_best_mtxpp;
  *ret_scale_fspp = (double*)(malloc_or_die((seq_len+2) * sizeof(double)));
  scale_fsp = *ret_scale_fspp;
  sorted_v_list = (int*)(malloc_or_die((hmmp->nr_v * 2 + 1) * sizeof(int)));
  
  /* Initialize the first row of the matrix */
  one_best_mtxp->prob = 1.0; /* sets index (0,0) to 1.0,
			    the rest are already 0.0 as they should be */
  one_best_mtxp->labeling = (char*)(malloc_or_die(1 * sizeof(char)));
  for(i = 1; i < hmmp->nr_v; i++) {
    (one_best_mtxp+i)->labeling = NULL;
  }
  *scale_fsp = 1.0;
  
  /* create initial sorted v-list*/
  *(sorted_v_list) = 0; /* 0 is always the number of the start state */
  *(sorted_v_list + 1) = V_LIST_NEXT;
  *(sorted_v_list + 2) = V_LIST_END;
  



  /* Fill in middle rows */
  prev_rowp = one_best_mtxp;
  cur_rowp = one_best_mtxp + hmmp->nr_v;
  j = 0;
  if(use_lead_columns == YES) {
    c = *(msa_seq_infop->lead_columns_start + j);
  }
  else {
    c = 0;
  }
  while(c != END && c < msa_seq_infop->msa_seq_length)  {
    a_index_2 = -1;
    a_index_3 = -1;
    a_index_4 = -1;
    if(hmmp->alphabet_type == DISCRETE) {
      a_index = get_alphabet_index((msa_seq_infop->msa_seq_1 + (c * (hmmp->a_size+1)))->query_letter, hmmp->alphabet, hmmp->a_size);
      if(a_index < 0) {
	a_index = hmmp->a_size; /* if letter is wild card, use default column in subst matrix */
      }
    }
    if(hmmp->nr_alphabets > 1) {
      if(hmmp->alphabet_type_2 == DISCRETE) {
	a_index_2 = get_alphabet_index((msa_seq_infop->msa_seq_2 + (c * (hmmp->a_size_2+1)))->query_letter,
				       hmmp->alphabet_2, hmmp->a_size_2);
	if(a_index_2 < 0) {
	a_index_2 = hmmp->a_size_2; /* if letter is wild card, use default column in subst matrix */
	}
      }
    }
    if(hmmp->nr_alphabets > 2) {
      if(hmmp->alphabet_type_3 == DISCRETE) {
	a_index_3 = get_alphabet_index((msa_seq_infop->msa_seq_3 + (c * (hmmp->a_size_3+1)))->query_letter,
				       hmmp->alphabet_3, hmmp->a_size_3);
	if(a_index_3 < 0) {
	  a_index_3 = hmmp->a_size_3; /* if letter is wild card, use default column in subst matrix */
	}
      }
    }
    if(hmmp->nr_alphabets > 3) {
      if(hmmp->alphabet_type_4 == DISCRETE) {
	a_index_4 = get_alphabet_index((msa_seq_infop->msa_seq_4 + (c * (hmmp->a_size_4+1)))->query_letter,
				       hmmp->alphabet_4, hmmp->a_size_4);
	if(a_index_4 < 0) {
	  a_index_4 = hmmp->a_size_4; /* if letter is wild card, use default column in subst matrix */
	}
      }
    }
    
    /* calculate sum of probabilities */
#ifdef DEBUG_FW
    printf("label for pos %d = %c\n", c, (msa_seq_infop->msa_seq_1 + (c * (hmmp->a_size+1)))->label);
#endif    
    for(v = 1; v < hmmp->nr_v - 1; v++) /* v = to-vertex */ {
#ifdef DEBUG_FW
      printf("prob to vertex %d:\n", v);
#endif
      
      v_list_index = 0;
      (cur_rowp + v)->labeling = NULL;
      while(*(sorted_v_list + v_list_index) != V_LIST_END) {
	res = 0.0;
	while(*(sorted_v_list + v_list_index) != V_LIST_NEXT) {
	  w = *(sorted_v_list + v_list_index); /* w = from-vertex */
	  /* calculate intermediate results */
	  res += (prev_rowp + w)->prob * *(hmmp->tot_transitions + (w * nr_v + v));
	  if(*(hmmp->tot_transitions + (w * nr_v + v)) < 0) {
	    printf("found model transition prob from %d to %d < 0.0\n", w, v);
	    exit(0);
	  }
#ifdef DEBUG_FW
	  printf("prev = %f: ", (prev_rowp + w)->prob);
	  printf("trans = %f\n", *(hmmp->tot_transitions + (w * nr_v + v)));
#endif
	  v_list_index++;
	}
	

	/* calculate the prob of producing letters l in v*/
	t_res3 = get_msa_emission_score_multi(msa_seq_infop, hmmp, c, v, use_labels, normalize, a_index, a_index_2,
					      a_index_3, a_index_4, scoring_method, use_gap_shares,
					      use_prior_shares, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3,
					      aa_freqs_4);
       
	res = res * t_res3;
	
	/* check if this score is best so far */
#ifdef DEBUG_FW
	printf("best score = %f\n",(cur_rowp + v)->prob);
#endif	
	if(res > (cur_rowp + v)->prob) {
#ifdef DEBUG_FW
	  printf("updating best score to %f\n", res);
#endif	  
	  /* save result in matrices */
	  (cur_rowp + v)->prob = res;
	  /* set pointer to point to current labeling */
	  (cur_rowp + v)->labeling = (prev_rowp + w)->labeling;
	  if((cur_rowp + v)->labeling == NULL) {
	    printf("Error: NULL labeling when updating best score\n");
	    exit(0);
	  }
	}
	
#ifdef DEBUG_FW
	printf("res = %f\n", res);
#endif
	v_list_index++;
      }
    }
    
    /* update labeling pointers */
    if(use_lead_columns == YES) {
      update_labelings(cur_rowp, hmmp->vertex_labels, sorted_v_list, seq_len, j+1, hmmp->labels, hmmp->nr_labels, hmmp->nr_v);
    }
    else {
      update_labelings(cur_rowp, hmmp->vertex_labels, sorted_v_list, seq_len, c+1, hmmp->labels, hmmp->nr_labels, hmmp->nr_v);
    }
    deallocate_row_labelings(prev_rowp, hmmp->nr_v);

    /* scale the results, row_sum = the total probability of 
     * having produced the labelings up to and including character c */
    row_sum = 0.0;
    for(v = 0; v < hmmp->nr_v; v++) {
      row_sum = row_sum + (cur_rowp + v)->prob;
    }
   
    if(use_lead_columns == YES) {
      scale_fsp[j+1] = row_sum;
    }
    else {
      scale_fsp[c+1] = row_sum;
    }
    
    
#ifdef DEBUG_FW
    printf("rowsum for row %d= %f\n", c, row_sum);
    printf("scaling set to: %f\n", scale_fsp[c+1]);
#endif
    if(row_sum == 0.0) {
      printf("Probability for this msa = 0.0\n");
      return NOPROB;
    }
    for(v = 0; v < hmmp->nr_v; v++) {
      if((cur_rowp + v)->prob != 0){
	(cur_rowp + v)->prob = ((cur_rowp + v)->prob)/row_sum; /* scaling */
      }
    }
    
    /* move row pointers one row one_best */
    prev_rowp = cur_rowp;
    cur_rowp = cur_rowp + hmmp->nr_v;

    /* update current column */
    if(use_lead_columns == YES) {
      j++;
      c = *(msa_seq_infop->lead_columns_start + j);
    }
    else {
      c++;
    }
  }


  /* Fill in transition to end state */
  v_list_index = 0;
  (cur_rowp + hmmp->nr_v - 1)->labeling = NULL;
  while(*(sorted_v_list + v_list_index) != V_LIST_END) {
    res = 0.0;
    while(*(sorted_v_list + v_list_index) != V_LIST_NEXT) {
      w = *(sorted_v_list + v_list_index); /* w = from-vertex */
      t_res1 = (prev_rowp + w)->prob;
      t_res2 = *((hmmp->tot_transitions) + get_mtx_index(w, hmmp->nr_v-1, hmmp->nr_v));
      if(t_res2 > 1.0) {
	t_res2 = 1.0;
      }
      res += t_res1 * t_res2;
      v_list_index++;
    }
    
    /* check if this score is best so far */
    if(res > (cur_rowp + hmmp->nr_v - 1)->prob) {
      /* save result in matrices */
      (cur_rowp + hmmp->nr_v - 1)->prob = res;
      /* set pointer to point to current labeling */
      (cur_rowp + hmmp->nr_v - 1)->labeling = (prev_rowp + w)->labeling;  
    }
    
    v_list_index++;
  }
#ifdef DEBUG_FW
  dump_one_best_matrix(seq_len + 2, hmmp->nr_v, one_best_mtxp);
  dump_scaling_array(seq_len + 1, scale_fsp);
#endif
  
  /* store results */
  scaled_result = (cur_rowp + hmmp->nr_v - 1)->prob;
  memcpy(best_labeling, ((cur_rowp + hmmp->nr_v - 1)->labeling) + 1, seq_len * sizeof(char));
  best_labeling[seq_len] = '\0';
#ifdef DEBUG_PATH
  printf("seq_len = %d\n", seq_len);
  printf("best labeling = %s\n", best_labeling);
#endif
  
  /* Garbage collection and return */
  free(sorted_v_list);
  /* FREE labelings, except result labeling */
  deallocate_row_labelings(prev_rowp, hmmp->nr_v);

  return OK;
}




/************************* help functions ***********************************/



double dot_product_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
			 struct msa_sequences_multi_s *msa_seq_infop, 
			 int c, int v, int normalize, int multi_scoring_method)
{
  int i,j;
  double t_res3_1, t_res3_2, t_res3_3, t_res3_4, t_res3_tot;
  double seq_normalizer;
  double state_normalizer;
  
  t_res3_tot = 0.0;
  
  t_res3_1 = 0.0;
  t_res3_2 = 0.0;
  t_res3_3 = 0.0;
  t_res3_4 = 0.0;
  
  if(hmmp->alphabet_type == DISCRETE) {
    t_res3_1 = get_dp_statescore(hmmp->a_size, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_1, c, hmmp->emissions,
				 v, normalize, msa_seq_infop->gap_shares);
  }
  else {
    if((msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->nr_occurences > 0.0) {
      t_res3_1 = 0.0;
      for(j = 0; j < hmmp->a_size / 3; j++) {
	t_res3_1 += get_single_gaussian_statescore((*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3))),
						   (*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 1))),
						   (msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->share) *
	  *((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 2));
      }
    }
    else {
      t_res3_1 = 1.0;
    }
  }
  if(hmmp->nr_alphabets > 1) {
    if(hmmp->alphabet_type_2 == DISCRETE) {
      t_res3_2 = get_dp_statescore(hmmp->a_size_2, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_2, c, hmmp->emissions_2,
				   v, normalize, msa_seq_infop->gap_shares);
    }
    else {
      if((msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->nr_occurences > 0.0) {
	t_res3_2 = 0.0;
	for(j = 0; j < hmmp->a_size_2 / 3; j++) {
	  t_res3_2 += get_single_gaussian_statescore((*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3))),
						     (*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->share) *
	    *((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_2 = 1.0;
      }
    }
  }

  if(hmmp->nr_alphabets > 2) {
    if(hmmp->alphabet_type_3 == DISCRETE) {
      t_res3_3 = get_dp_statescore(hmmp->a_size_3, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_3, c, hmmp->emissions_3,
				   v, normalize, msa_seq_infop->gap_shares);
    }
    else {
      if((msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->nr_occurences > 0.0) {
	t_res3_3 = 0.0;
	for(j = 0; j < hmmp->a_size_3 / 3; j++) {
	  t_res3_3 += get_single_gaussian_statescore((*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3))),
						     (*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->share) *
	    *((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_3 = 1.0;
      }
    }
  }

  if(hmmp->nr_alphabets > 3) {
    if(hmmp->alphabet_type_4 == DISCRETE) {
      t_res3_4 = get_dp_statescore(hmmp->a_size_4, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_4, c, hmmp->emissions_4,
				   v, normalize, msa_seq_infop->gap_shares);
    }
    else {
      if((msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->nr_occurences > 0.0) {
	t_res3_4 = 0.0;
	for(j = 0; j < hmmp->a_size_4 / 3; j++) {
	  t_res3_4 += get_single_gaussian_statescore((*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3))),
						     (*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->share) *
	    *((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_4 = 1.0;
      }
    }
  }

  
  if(multi_scoring_method == JOINT_PROB) {
    t_res3_tot = t_res3_1;
    if(hmmp->nr_alphabets > 1) {
      t_res3_tot *= t_res3_2;
    }
    if(hmmp->nr_alphabets > 2) {
      t_res3_tot *= t_res3_3;
    }
    if(hmmp->nr_alphabets > 3) {
      t_res3_tot *= t_res3_4;
    }
    return t_res3_tot;
  }
  else {
    /* use other multialpha scoring method, not implemented yet */
    printf("Error: only JOINT_PROB scoring is implemented\n");
    exit(0);
  }
}

double dot_product_picasso_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
				 struct msa_sequences_multi_s *msa_seq_infop, 
				 int c, int v, int normalize, int multi_scoring_method, double *aa_freqs,
				 double *aa_freqs_2, double *aa_freqs_3, double *aa_freqs_4)
{
  int i,j;
  double t_res3_1, t_res3_2, t_res3_3, t_res3_4, t_res3_tot;
  double seq_normalizer;
  double state_normalizer;
  
  t_res3_tot = 0.0;
  
  t_res3_1 = 0.0;
  t_res3_2 = 0.0;
  t_res3_3 = 0.0;
  t_res3_4 = 0.0;
  
  
  if(hmmp->alphabet_type == DISCRETE) {
    t_res3_1 = get_dp_picasso_statescore(hmmp->a_size, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_1,
					 c, hmmp->emissions,
					 v, normalize, msa_seq_infop->gap_shares, aa_freqs);
  }
  else {
    if((msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->nr_occurences > 0.0) {
      t_res3_1 = 0.0;
      for(j = 0; j < hmmp->a_size / 3; j++) {
	t_res3_1 += get_single_gaussian_statescore((*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3))),
						   (*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 1))),
						   (msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->share) *
	  *((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 2));
      }
    }
    else {
      t_res3_1 = 1.0;
    }
  }
  if(hmmp->nr_alphabets > 1) {
    if(hmmp->alphabet_type_2 == DISCRETE) {
      t_res3_2 = get_dp_picasso_statescore(hmmp->a_size_2, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_2,
					   c, hmmp->emissions_2,
					   v, normalize, msa_seq_infop->gap_shares, aa_freqs_2);
    }
    else {
      if((msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->nr_occurences > 0.0) {
	t_res3_2 = 0.0;
	for(j = 0; j < hmmp->a_size_2 / 3; j++) {
	  t_res3_2 += get_single_gaussian_statescore((*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3))),
						     (*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->share) *
	    *((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_2 = 1.0;
      }
    }
  }
  
  if(hmmp->nr_alphabets > 2) {
    if(hmmp->alphabet_type_3 == DISCRETE) {
      t_res3_3 = get_dp_picasso_statescore(hmmp->a_size_3, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_3,
					   c, hmmp->emissions_3,
					   v, normalize, msa_seq_infop->gap_shares, aa_freqs_3);
    }
    else {
      if((msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->nr_occurences > 0.0) {
	t_res3_3 = 0.0;
	for(j = 0; j < hmmp->a_size_3 / 3; j++) {
	  t_res3_3 += get_single_gaussian_statescore((*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3))),
						     (*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->share) *
	    *((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_3 = 1.0;
      }
    }
  }
  
  if(hmmp->nr_alphabets > 3) {
    if(hmmp->alphabet_type_4 == DISCRETE) {
      t_res3_4 = get_dp_picasso_statescore(hmmp->a_size_4, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_4,
					   c, hmmp->emissions_4,
					   v, normalize, msa_seq_infop->gap_shares, aa_freqs_4);
    }
    else {
      if((msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->nr_occurences > 0.0) {
	t_res3_4 = 0.0;
	for(j = 0; j < hmmp->a_size_4 / 3; j++) {
	  t_res3_4 += get_single_gaussian_statescore((*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3))),
						     (*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->share) *
	    *((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_4 = 1.0;
      }
    }
  }
  
  
  if(multi_scoring_method == JOINT_PROB) {
    t_res3_tot = t_res3_1;
    if(hmmp->nr_alphabets > 1) {
      t_res3_tot *= t_res3_2;
    }
    if(hmmp->nr_alphabets > 2) {
      t_res3_tot *= t_res3_3;
    }
    if(hmmp->nr_alphabets > 3) {
      t_res3_tot *= t_res3_4;
    }
    return t_res3_tot;
  }
  else {
    /* use other multialpha scoring method, not implemented yet */
    printf("Error: only JOINT_PROB scoring is implemented\n");
    exit(0);
  }
}


double picasso_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
		     struct msa_sequences_multi_s *msa_seq_infop, 
		     int c, int v, int normalize, int multi_scoring_method, double *aa_freqs,
		     double *aa_freqs_2, double *aa_freqs_3, double *aa_freqs_4)
{
  int i,j;
  double t_res3_1, t_res3_2, t_res3_3, t_res3_4, t_res3_tot;
  double seq_normalizer;
  double state_normalizer;
  
  t_res3_tot = 0.0;
  
  t_res3_1 = 0.0;
  t_res3_2 = 0.0;
  t_res3_3 = 0.0;
  t_res3_4 = 0.0;
  
  
  if(hmmp->alphabet_type == DISCRETE) {
    t_res3_1 = get_picasso_statescore(hmmp->a_size, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_1, c, hmmp->emissions,
				      v, normalize, msa_seq_infop->gap_shares, aa_freqs);
    
  }
  else {
    if((msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->nr_occurences > 0.0) {
      t_res3_1 = 0.0;
      for(j = 0; j < hmmp->a_size / 3; j++) {
	t_res3_1 += get_single_gaussian_statescore((*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3))),
						   (*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 1))),
						   (msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->share) *
	  *((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 2));
      }
    }
    else {
      t_res3_1 = 1.0;
    }
  }
  if(hmmp->nr_alphabets > 1) {
    if(hmmp->alphabet_type_2 == DISCRETE) {
      t_res3_2 = get_picasso_statescore(hmmp->a_size_2, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_2, c,
					hmmp->emissions_2, v, normalize, msa_seq_infop->gap_shares, aa_freqs_2);
    }
    else {
      if((msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->nr_occurences > 0.0) {
	t_res3_2 = 0.0;
	for(j = 0; j < hmmp->a_size_2 / 3; j++) {
	  t_res3_2 += get_single_gaussian_statescore((*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3))),
						     (*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->share) *
	    *((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_2 = 1.0;
      }
    }
  }
  
  if(hmmp->nr_alphabets > 2) {
    if(hmmp->alphabet_type_3 == DISCRETE) {
      t_res3_3 = get_picasso_statescore(hmmp->a_size_3, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_3, c,
					hmmp->emissions_3, v, normalize, msa_seq_infop->gap_shares, aa_freqs_3);
    }
    else {
      if((msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->nr_occurences > 0.0) {
	t_res3_3 = 0.0;
	for(j = 0; j < hmmp->a_size_3 / 3; j++) {
	  t_res3_3 += get_single_gaussian_statescore((*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3))),
						     (*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->share) *
	    *((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_3 = 1.0;
      }
    }
  }
  
  if(hmmp->nr_alphabets > 3) {
    if(hmmp->alphabet_type_4 == DISCRETE) {
      t_res3_4 = get_picasso_statescore(hmmp->a_size_4, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_4, c,
					hmmp->emissions_4, v, normalize, msa_seq_infop->gap_shares, aa_freqs_4);
    }
    else {
       if((msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->nr_occurences > 0.0) {
	 t_res3_4 = 0.0;
	for(j = 0; j < hmmp->a_size_4 / 3; j++) {
	  t_res3_4 += get_single_gaussian_statescore((*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3))),
						     (*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->share) *
	    *((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 2));
	}
       }
       else {
	 t_res3_4 = 1.0;
       }
    }
  }
  
  
  if(multi_scoring_method == JOINT_PROB) {
    t_res3_tot = t_res3_1;
    if(hmmp->nr_alphabets > 1) {
      t_res3_tot *= t_res3_2;
    }
    if(hmmp->nr_alphabets > 2) {
      t_res3_tot *= t_res3_3;
    }
    if(hmmp->nr_alphabets > 3) {
      t_res3_tot *= t_res3_4;
    }
    return t_res3_tot;
  }
  else {
    /* use other multialpha scoring method, not implemented yet */
    printf("Error: only JOINT_PROB scoring is implemented\n");
    exit(0);
  }
}

double picasso_sym_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
			 struct msa_sequences_multi_s *msa_seq_infop, 
			 int c, int v, int normalize, int multi_scoring_method, double *aa_freqs,
			 double *aa_freqs_2, double *aa_freqs_3, double *aa_freqs_4)
{
  int i,j;
  double t_res3_1, t_res3_2, t_res3_3, t_res3_4, t_res3_tot;
  double seq_normalizer;
  double state_normalizer;
  
  t_res3_tot = 0.0;
  
  t_res3_1 = 0.0;
  t_res3_2 = 0.0;
  t_res3_3 = 0.0;
  t_res3_4 = 0.0;
  
  
  
  if(hmmp->alphabet_type == DISCRETE) {
    t_res3_1 = get_picasso_sym_statescore(hmmp->a_size, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_1, c,
					  hmmp->emissions,
					  v, normalize, msa_seq_infop->gap_shares, aa_freqs);
  }
  else {
    if((msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->nr_occurences > 0.0) {
      t_res3_1 = 0.0;
      for(j = 0; j < hmmp->a_size / 3; j++) {
	t_res3_1 += get_single_gaussian_statescore((*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3))),
						   (*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 1))),
						   (msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->share) *
	  *((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 2));
      }
    }
    else {
      t_res3_1 = 1.0;
    }
  }
  if(hmmp->nr_alphabets > 1) {
    if(hmmp->alphabet_type_2 == DISCRETE) {
      t_res3_2 = get_picasso_sym_statescore(hmmp->a_size_2, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_2, c,
					  hmmp->emissions_2,
					    v, normalize, msa_seq_infop->gap_shares, aa_freqs_2);
    }
    else {
      if((msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->nr_occurences > 0.0) {
	t_res3_2 = 0.0;
	for(j = 0; j < hmmp->a_size_2 / 3; j++) {
	  t_res3_2 += get_single_gaussian_statescore((*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3))),
						     (*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->share) *
	    *((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_2 = 1.0;
      }
    }
  }
  if(hmmp->nr_alphabets > 2) {
    if(hmmp->alphabet_type_3 == DISCRETE) {
      t_res3_3 = get_picasso_sym_statescore(hmmp->a_size_3, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_3, c,
					    hmmp->emissions_3,
					    v, normalize, msa_seq_infop->gap_shares, aa_freqs_3);
    }
    else {
      if((msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->nr_occurences > 0.0) {
	t_res3_3 = 0.0;
	for(j = 0; j < hmmp->a_size_3 / 3; j++) {
	  t_res3_3 += get_single_gaussian_statescore((*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3))),
						     (*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->share) *
	    *((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_3 = 1.0;
      }
    }
  }

  if(hmmp->nr_alphabets > 3) {
    if(hmmp->alphabet_type_4 == DISCRETE) {
      t_res3_4 = get_picasso_sym_statescore(hmmp->a_size_4, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_4, c,
					    hmmp->emissions_4,
					    v, normalize, msa_seq_infop->gap_shares, aa_freqs_4);
    }
    else {
      if((msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->nr_occurences > 0.0) {
	t_res3_4 = 0.0;
	for(j = 0; j < hmmp->a_size_4 / 3; j++) {
	  t_res3_4 += get_single_gaussian_statescore((*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3))),
						     (*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->share) *
	    *((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_4 = 1.0;
      }
    }

  }
  
  if(multi_scoring_method == JOINT_PROB) {
    t_res3_tot = t_res3_1;
    if(hmmp->nr_alphabets > 1) {
      t_res3_tot *= t_res3_2;
    }
    if(hmmp->nr_alphabets > 2) {
      t_res3_tot *= t_res3_3;
    }
    if(hmmp->nr_alphabets > 3) {
      t_res3_tot *= t_res3_4;
    }
    return t_res3_tot;
  }
  else {
    /* use other multialpha scoring method, not implemented yet */
    printf("Error: only JOINT_PROB scoring is implemented\n");
    exit(0);
  }
}



double sjolander_score_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
			     struct msa_sequences_multi_s *msa_seq_infop, 
			     int c, int v, int normalize, int multi_scoring_method)
{
  int i,j;
  double t_res3_tot, t_res3_1, t_res3_2, t_res3_3, t_res3_4;
  double seq_normalizer;
  double state_normalizer;
  
  
  t_res3_tot = 0.0;
  t_res3_1 = 1.0;
  t_res3_2 = 1.0;
  t_res3_3 = 1.0;
  t_res3_4 = 1.0;
  

  if(hmmp->alphabet_type == DISCRETE) {
    t_res3_1 = get_sjolander_statescore(hmmp->a_size, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_1, c, hmmp->emissions,
					v, normalize, msa_seq_infop->gap_shares);
  }
  else {
    if((msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->nr_occurences > 0.0) {
      t_res3_1 = 0.0;
      for(j = 0; j < hmmp->a_size / 3; j++) {
	t_res3_1 += get_single_gaussian_statescore((*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3))),
						   (*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 1))),
						   (msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->share) *
	  *((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 2));
      }
    }
    else {
      t_res3_1 = 1.0;
    }
  }
  if(hmmp->nr_alphabets > 1) {
    if(hmmp->alphabet_type_2 == DISCRETE) {
      t_res3_2 = get_sjolander_statescore(hmmp->a_size_2, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_2,
					  c, hmmp->emissions_2,
					  v, normalize, msa_seq_infop->gap_shares);
    }
    else {
      if((msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->nr_occurences > 0.0) {
	t_res3_2 = 0.0;
	for(j = 0; j < hmmp->a_size_2 / 3; j++) {
	  t_res3_2 += get_single_gaussian_statescore((*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3))),
						     (*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->share) *
	    *((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_2 = 1.0;
      }
    }
  }

  if(hmmp->nr_alphabets > 2) {
    if(hmmp->alphabet_type_3 == DISCRETE){
      t_res3_3 = get_sjolander_statescore(hmmp->a_size_3, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_3,
					  c, hmmp->emissions_3,
					  v, normalize, msa_seq_infop->gap_shares);
    }
    else {
      if((msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->nr_occurences > 0.0) {
	t_res3_3 = 0.0;
	for(j = 0; j < hmmp->a_size_3 / 3; j++) {
	  t_res3_3 += get_single_gaussian_statescore((*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3))),
						     (*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->share) *
	    *((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_3 = 1.0;
      }
    }
  }
  if(hmmp->nr_alphabets > 3) {
    if(hmmp->alphabet_type_4 == DISCRETE) {
      t_res3_4 = get_sjolander_statescore(hmmp->a_size_4, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_4,
					  c, hmmp->emissions_4,
					  v, normalize, msa_seq_infop->gap_shares);
    }
    else {
      if((msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->nr_occurences > 0.0) {
	t_res3_4 = 0.0;
	for(j = 0; j < hmmp->a_size_4 / 3; j++) {
	  t_res3_4 += get_single_gaussian_statescore((*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3))),
						     (*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->share) *
	    *((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_4 = 1.0;
      }
    }
  }
  
  if(multi_scoring_method == JOINT_PROB) {
    t_res3_tot = t_res3_1;
    if(hmmp->nr_alphabets > 1) {
      t_res3_tot *= t_res3_2;
    }
    if(hmmp->nr_alphabets > 2) {
      t_res3_tot *= t_res3_3;
    }
    if(hmmp->nr_alphabets > 3) {
      t_res3_tot *= t_res3_4;
    }
    return t_res3_tot;
  }
  else {
    /* use other multialpha scoring method, not implemented yet */
    printf("Error: only JOINT_PROB scoring is implemented\n");
    exit(0);
  }
}

double sjolander_reversed_score_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
				      struct msa_sequences_multi_s *msa_seq_infop, 
				      int c, int v, int normalize, int multi_scoring_method)
{
  int i,j;
  double t_res3_tot, t_res3_1, t_res3_2, t_res3_3, t_res3_4;
  double seq_normalizer;
  double state_normalizer;
  
  
  t_res3_tot = 0.0;
  t_res3_1 = 1.0;
  t_res3_2 = 1.0;
  t_res3_3 = 1.0;
  t_res3_4 = 1.0;
  
  if(hmmp->alphabet_type == DISCRETE) {
    t_res3_1 = get_sjolander_reversed_statescore(hmmp->a_size, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_1,
						 c, hmmp->emissions,
						 v, normalize, msa_seq_infop->gap_shares);
  }
  else {
    if((msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->nr_occurences > 0.0) {
      t_res3_1 = 0.0;
      for(j = 0; j < hmmp->a_size / 3; j++) {
	t_res3_1 += get_single_gaussian_statescore((*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3))),
						   (*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 1))),
						   (msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->share) *
	  *((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 2));
      }
    }
    else {
      t_res3_1 = 1.0;
    }
  }
  if(hmmp->nr_alphabets > 1) {
    if(hmmp->alphabet_type_2 == DISCRETE) {
      t_res3_2 = get_sjolander_reversed_statescore(hmmp->a_size_2, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_2,
						   c, hmmp->emissions_2,
						   v, normalize, msa_seq_infop->gap_shares);
    }
    else {
      if((msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->nr_occurences > 0.0) {
	t_res3_2 = 0.0;
	for(j = 0; j < hmmp->a_size_2 / 3; j++) {
	  t_res3_2 += get_single_gaussian_statescore((*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3))),
						     (*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->share) *
	    *((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_2 = 1.0;
      }
    }
  }
  if(hmmp->nr_alphabets > 2) {
    if(hmmp->alphabet_type_3 == DISCRETE) {
      t_res3_3 = get_sjolander_reversed_statescore(hmmp->a_size_3, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_3,
						   c, hmmp->emissions_3,
						   v, normalize, msa_seq_infop->gap_shares);
    }
    else {
      if((msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->nr_occurences > 0.0) {
	t_res3_3 = 0.0;
	for(j = 0; j < hmmp->a_size_3 / 3; j++) {
	  t_res3_3 += get_single_gaussian_statescore((*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3))),
						     (*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->share) *
	    *((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_3 = 1.0;
      }
    }
  }
  if(hmmp->nr_alphabets > 3) {
    if(hmmp->alphabet_type_4 == DISCRETE) {
      t_res3_4 = get_sjolander_reversed_statescore(hmmp->a_size_4, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_4,
						   c, hmmp->emissions_4,
						   v, normalize, msa_seq_infop->gap_shares);
    }
    else {
      if((msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->nr_occurences > 0.0) {
	t_res3_4 = 0.0;
	for(j = 0; j < hmmp->a_size_4 / 3; j++) {
	  t_res3_4 += get_single_gaussian_statescore((*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3))),
						     (*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->share) *
	    *((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_4 = 1.0;
      }
    }
  }
  
  if(multi_scoring_method == JOINT_PROB) {
    t_res3_tot = t_res3_1;
    if(hmmp->nr_alphabets > 1) {
      t_res3_tot *= t_res3_2;
    }
    if(hmmp->nr_alphabets > 2) {
      t_res3_tot *= t_res3_3;
    }
    if(hmmp->nr_alphabets > 3) {
      t_res3_tot *= t_res3_4;
    }
    return t_res3_tot;
  }
  else {
    /* use other multialpha scoring method, not implemented yet */
    printf("Error: only JOINT_PROB scoring is implemented\n");
    exit(0);
  }
}

double subst_mtx_product_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
			       struct msa_sequences_multi_s *msa_seq_infop, int c, int v, int normalize, int multi_scoring_method)
{
  int i,j;
  double t_res3_1, t_res3_2, t_res3_3, t_res3_4, t_res3_tot;
  /* Note: normalization not used for this scoring method */
  
  t_res3_1 = 0.0;
  t_res3_2 = 0.0;
  t_res3_3 = 0.0;
  t_res3_4 = 0.0;
  t_res3_tot = 0.0;
  
  if(hmmp->alphabet_type == DISCRETE) {
    t_res3_1 = get_subst_mtx_product_statescore(hmmp->a_size, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_1,
						c, hmmp->emissions,
						v, hmmp->subst_mtx);
  }
  else {
    if((msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->nr_occurences > 0.0) {
      t_res3_1 = 0.0;
      for(j = 0; j < hmmp->a_size / 3; j++) {
	t_res3_1 += get_single_gaussian_statescore((*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3))),
						   (*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 1))),
						   (msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->share) *
	  *((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 2));
      }
    }
    else {
      t_res3_1 = 1.0;
    }
  }
  if(hmmp->nr_alphabets > 1) {
    if(hmmp->alphabet_type_2 == DISCRETE) {
      t_res3_2 = get_subst_mtx_product_statescore(hmmp->a_size_2, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_2,
						  c, hmmp->emissions_2,
						  v, hmmp->subst_mtx_2);
    }
    else {
      if((msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->nr_occurences > 0.0) {
	t_res3_2 = 0.0;
	for(j = 0; j < hmmp->a_size_2 / 3; j++) {
	  t_res3_2 += get_single_gaussian_statescore((*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3))),
						     (*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->share) *
	    *((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_2 = 1.0;
      }
    }
  }
  
  if(hmmp->nr_alphabets > 2) {
    if(hmmp->alphabet_type_3 == DISCRETE) {
      t_res3_3 = get_subst_mtx_product_statescore(hmmp->a_size_3, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_3,
						  c, hmmp->emissions_3,
						  v, hmmp->subst_mtx_3);
    }
    else {
      if((msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->nr_occurences > 0.0) {
	t_res3_3 = 0.0;
	for(j = 0; j < hmmp->a_size_3 / 3; j++) {
	  t_res3_3 += get_single_gaussian_statescore((*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3))),
						     (*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->share) *
	    *((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_3 = 1.0;
      }
    }
  }
  
  if(hmmp->nr_alphabets > 3) {
    if(hmmp->alphabet_type_4 == DISCRETE) {
      t_res3_4 = get_subst_mtx_product_statescore(hmmp->a_size_4, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_4,
						  c, hmmp->emissions_4,
						  v, hmmp->subst_mtx_4);
    }
    else {
      if((msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->nr_occurences > 0.0) {
	t_res3_4 = 0.0;
	for(j = 0; j < hmmp->a_size_4 / 3; j++) {
	  t_res3_4 += get_single_gaussian_statescore((*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3))),
						     (*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->share) *
	    *((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_4 = 1.0;
      }
    }
  }
  

  if(multi_scoring_method == JOINT_PROB) {
    t_res3_tot = t_res3_1;
    if(hmmp->nr_alphabets > 1) {
      t_res3_tot *= t_res3_2;
    }
    if(hmmp->nr_alphabets > 2) {
      t_res3_tot *= t_res3_3;
    }
    if(hmmp->nr_alphabets > 3) {
      t_res3_tot *= t_res3_4;
    }
    return t_res3_tot;
  }
  else {
    /* use other multialpha scoring method, not implemented yet */
    printf("Error: only JOINT_PROB scoring is implemented\n");
    exit(0);
  }
}

double subst_mtx_dot_product_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
				   struct msa_sequences_multi_s *msa_seq_infop, int c, int v, int a_index, int a_index_2, 
				   int a_index_3, int a_index_4, int normalize,
				   int multi_scoring_method)
{
  int i,j;
  double t_res3_1, t_res3_2, t_res3_3, t_res3_4, t_res3_tot;
  double seq_normalizer;
  double state_normalizer;
  double subst_mtx_normalizer;

  
  

  t_res3_1 = 0.0;
  t_res3_2 = 0.0;
  t_res3_3 = 0.0;
  t_res3_4 = 0.0;
  t_res3_tot = 0.0;

  if(hmmp->alphabet_type == DISCRETE) {
    t_res3_1 = get_subst_mtx_dot_product_statescore(hmmp->a_size, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_1,
						    c, hmmp->emissions,
						    v, normalize, msa_seq_infop->gap_shares, a_index, hmmp->subst_mtx);
  }
  else {
    if((msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->nr_occurences > 0.0) {
      t_res3_1 = 0.0;
      for(j = 0; j < hmmp->a_size / 3; j++) {
	t_res3_1 += get_single_gaussian_statescore((*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3))),
						   (*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 1))),
						   (msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->share) *
	  *((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 2));
      }
    }
    else {
      t_res3_1 = 1.0;
    }
  }
  if(hmmp->nr_alphabets > 1) {
    if(hmmp->alphabet_type_2 == DISCRETE) {
      t_res3_2 = get_subst_mtx_dot_product_statescore(hmmp->a_size_2, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_2,
						      c, hmmp->emissions_2,
						      v, normalize, msa_seq_infop->gap_shares, a_index_2, hmmp->subst_mtx_2);
    }
    else {
      if((msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->nr_occurences > 0.0) {
	t_res3_2 = 0.0;
	for(j = 0; j < hmmp->a_size_2 / 3; j++) {
	  t_res3_2 += get_single_gaussian_statescore((*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3))),
						     (*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->share) *
	    *((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_2 = 1.0;
      }
    }
  }
  
  if(hmmp->nr_alphabets > 2) {
    if(hmmp->alphabet_type_3 == DISCRETE) {
      t_res3_3 = get_subst_mtx_dot_product_statescore(hmmp->a_size_3, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_3,
						      c, hmmp->emissions_3,
						      v, normalize, msa_seq_infop->gap_shares, a_index_3, hmmp->subst_mtx_3);
    }
    else {
      if((msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->nr_occurences > 0.0) {
	t_res3_3 = 0.0;
	for(j = 0; j < hmmp->a_size_3 / 3; j++) {
	  t_res3_3 += get_single_gaussian_statescore((*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3))),
						     (*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->share) *
	    *((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_3 = 1.0;
      }
    }
  }
  if(hmmp->nr_alphabets > 3) {
    if(hmmp->alphabet_type_4 == DISCRETE) {
      t_res3_4 = get_subst_mtx_dot_product_statescore(hmmp->a_size_4, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_4,
						      c, hmmp->emissions_4,
						      v, normalize, msa_seq_infop->gap_shares, a_index_4, hmmp->subst_mtx_4);
    }
    else {
      if((msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->nr_occurences > 0.0) {
	t_res3_4 = 0.0;
	for(j = 0; j < hmmp->a_size_4 / 3; j++) {
	  t_res3_4 += get_single_gaussian_statescore((*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3))),
						     (*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->share) *
	    *((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_4 = 1.0;
      }
    }
  }
  
  
  if(multi_scoring_method == JOINT_PROB) {
    t_res3_tot = t_res3_1;
    if(hmmp->nr_alphabets > 1) {
      t_res3_tot *= t_res3_2;
    }
    if(hmmp->nr_alphabets > 2) {
      t_res3_tot *= t_res3_3;
    }
    if(hmmp->nr_alphabets > 3) {
      t_res3_tot *= t_res3_4;
    }
    return t_res3_tot;
  }
  else {
    /* use other multialpha scoring method, not implemented yet */
    printf("Error: only JOINT_PROB scoring is implemented\n");
    exit(0);
  }
}

double subst_mtx_dot_product_prior_multi(struct hmm_multi_s *hmmp, int use_prior_shares, int use_gap_shares,
					 struct msa_sequences_multi_s *msa_seq_infop, 
					 int c, int v, int a_index, int a_index_2, int a_index_3, int a_index_4,
					 int normalize, int multi_scoring_method)
{
  int i,j;
  double t_res3_1, t_res3_2, t_res3_3, t_res3_4, t_res3_tot;
  double rest_share;
  double default_share;
  double seq_normalizer;
  double state_normalizer;
  double subst_mtx_normalizer;

  t_res3_1 = 0.0;
  t_res3_2 = 0.0;
  t_res3_3 = 0.0;
  t_res3_4 = 0.0;
  t_res3_tot = 0.0;

  
  if(hmmp->alphabet_type == DISCRETE) {
    t_res3_1 = get_subst_mtx_dot_product_prior_statescore(hmmp->a_size, use_gap_shares, use_prior_shares, msa_seq_infop->msa_seq_1,
							  c, hmmp->emissions,
							  v, normalize, msa_seq_infop->gap_shares, a_index, hmmp->subst_mtx); 
  }
  else {
    if((msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->nr_occurences > 0.0) {
      t_res3_1 = 0.0;
      for(j = 0; j < hmmp->a_size / 3; j++) {
	t_res3_1 += get_single_gaussian_statescore((*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3))),
						   (*((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 1))),
						   (msa_seq_infop->msa_seq_1 + get_mtx_index(c, 0, hmmp->a_size+1))->share) *
	  *((hmmp->emissions) + (v * (hmmp->a_size)) + (j * 3 + 2));
      }
    }
    else {
      t_res3_1 = 1.0;
    }
  }
  if(hmmp->nr_alphabets > 1) {
    if(hmmp->alphabet_type_2 == DISCRETE) {
      t_res3_2 = get_subst_mtx_dot_product_prior_statescore(hmmp->a_size_2, use_gap_shares, use_prior_shares,
							    msa_seq_infop->msa_seq_2,
							    c, hmmp->emissions_2,
							    v, normalize, msa_seq_infop->gap_shares, a_index_2, hmmp->subst_mtx_2);
    }
    else {
      if((msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->nr_occurences > 0.0) {
	t_res3_2 = 0.0;
	for(j = 0; j < hmmp->a_size_2 / 3; j++) {
	  t_res3_2 += get_single_gaussian_statescore((*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3))),
						     (*((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_2 + get_mtx_index(c, 0, hmmp->a_size_2+1))->share) *
	    *((hmmp->emissions_2) + (v * (hmmp->a_size_2)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_2 = 1.0;
      }
    }
  }

  if(hmmp->nr_alphabets > 2) {
    if(hmmp->alphabet_type_3 == DISCRETE) {
      t_res3_3 = get_subst_mtx_dot_product_prior_statescore(hmmp->a_size_3, use_gap_shares, use_prior_shares,
							    msa_seq_infop->msa_seq_3,
							    c, hmmp->emissions_3,
							    v, normalize, msa_seq_infop->gap_shares, a_index_3, hmmp->subst_mtx_3);
    }
    else {
      if((msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->nr_occurences > 0.0) {
	t_res3_3 = 0.0;
	for(j = 0; j < hmmp->a_size_3 / 3; j++) {
	  t_res3_3 += get_single_gaussian_statescore((*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3))),
						     (*((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_3 + get_mtx_index(c, 0, hmmp->a_size_3+1))->share) *
	    *((hmmp->emissions_3) + (v * (hmmp->a_size_3)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_3 = 1.0;
      }
    }
  }

  if(hmmp->nr_alphabets > 3) {
    if(hmmp->alphabet_type_4 == DISCRETE) {
      t_res3_4 = get_subst_mtx_dot_product_prior_statescore(hmmp->a_size_4, use_gap_shares, use_prior_shares,
							    msa_seq_infop->msa_seq_4,
							    c, hmmp->emissions_4,
							    v, normalize, msa_seq_infop->gap_shares, a_index_4, hmmp->subst_mtx_4);
    }
    else {
      if((msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->nr_occurences > 0.0) {
	t_res3_4 = 0.0;
	for(j = 0; j < hmmp->a_size_4 / 3; j++) {
	  t_res3_4 += get_single_gaussian_statescore((*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3))),
						     (*((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 1))),
						     (msa_seq_infop->msa_seq_4 + get_mtx_index(c, 0, hmmp->a_size_4+1))->share) *
	    *((hmmp->emissions_4) + (v * (hmmp->a_size_4)) + (j * 3 + 2));
	}
      }
      else {
	t_res3_4 = 1.0;
      }
    }
  }


  if(multi_scoring_method == JOINT_PROB) {
    t_res3_tot = t_res3_1;
    if(hmmp->nr_alphabets > 1) {
      t_res3_tot *= t_res3_2;
    }
    if(hmmp->nr_alphabets > 2) {
      t_res3_tot *= t_res3_3;
    }
    if(hmmp->nr_alphabets > 3) {
      t_res3_tot *= t_res3_4;
    }
    return t_res3_tot;
  }
  else {
    /* use other multialpha scoring method, not implemented yet */
    printf("Error: only JOINT_PROB scoring is implemented\n");
    exit(0);
  }
}


double get_msa_emission_score_multi(struct msa_sequences_multi_s *msa_seq_infop, struct hmm_multi_s *hmmp, int c, int v,
				    int use_labels, int normalize, int a_index, int a_index_2,
				    int a_index_3, int a_index_4, int scoring_method, int use_gap_shares,
				    int use_prior_shares, int multi_scoring_method, double *aa_freqs,
				    double *aa_freqs_2, double *aa_freqs_3, double *aa_freqs_4)
{

  double t_res3;


  /* multiply the prob of reaching state v with the prob of producing letters l in v*/
  t_res3 = 0.0;
  if(use_labels == YES && (msa_seq_infop->msa_seq_1 + (c * (hmmp->a_size+1)))->label != *(hmmp->vertex_labels + v) &&
     (msa_seq_infop->msa_seq_1 + (c * (hmmp->a_size+1)))->label != '.') {
  }
  /* calculate the simple dot product of the hmm-state vector and the msa vector */
  else if(scoring_method == DOT_PRODUCT) {
    t_res3 += dot_product_multi(hmmp, use_prior_shares, use_gap_shares, msa_seq_infop, c, v, normalize, multi_scoring_method);
  }
  else if(scoring_method == DOT_PRODUCT_PICASSO) {
    t_res3 += dot_product_picasso_multi(hmmp, use_prior_shares, use_gap_shares, msa_seq_infop, c, v, normalize,
					multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
  }
  else if(scoring_method == PICASSO) {
    t_res3 += picasso_multi(hmmp, use_prior_shares, use_gap_shares, msa_seq_infop, c, v, normalize, multi_scoring_method,
			    aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
 }
  
  else if(scoring_method == PICASSO_SYM) {
    t_res3 += picasso_sym_multi(hmmp, use_prior_shares, use_gap_shares, msa_seq_infop, c, v, normalize,
				multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
  }
  /* calculate sjolander score for hmm-state vector and msa vector */
  else if(scoring_method == SJOLANDER) {
    t_res3 += sjolander_score_multi(hmmp, use_prior_shares, use_gap_shares, msa_seq_infop, c, v, normalize, multi_scoring_method);
  }
  /* calculate sjolander score for hmm-state vector and msa vector */
  else if(scoring_method == SJOLANDER_REVERSED) {
    t_res3 += sjolander_reversed_score_multi(hmmp, use_prior_shares, use_gap_shares, msa_seq_infop, c, v, normalize,
					     multi_scoring_method);
  }

  /* calculate the joint sum of the emissions-vector and the msa-vector 
   * multiplying each result with the probability of the two amino acids being related, taken from 
   * a given substitution matrix */
  else if(scoring_method == SUBST_MTX_PRODUCT) {
    t_res3 += subst_mtx_product_multi(hmmp, use_prior_shares, use_gap_shares, msa_seq_infop, c, v, normalize,
				      multi_scoring_method);
  }
  /* a simpler and faster form of substitution mtx product, which multiplies the values of amino acids i
   * in the two columns and multiplies this value with a scaling factor taken from the substitition matrix, which
   * depends on the original aa in the query sequence and aa i */
  else if(scoring_method == SUBST_MTX_DOT_PRODUCT) {
    t_res3 += subst_mtx_dot_product_multi(hmmp, use_prior_shares, use_gap_shares, msa_seq_infop, c, v, a_index,
					  a_index_2, a_index_3, a_index_4, normalize,
					  multi_scoring_method);
  }
  /* like subst mtx dot product, but with a prior prob for the part of the seq-column not used in the
   * multiplication */
  else if(scoring_method == SUBST_MTX_DOT_PRODUCT_PRIOR) {
    t_res3 += subst_mtx_dot_product_prior_multi(hmmp, use_prior_shares, use_gap_shares, msa_seq_infop, c, v, a_index,
						a_index_2, a_index_3, a_index_4, normalize,
						multi_scoring_method);
  }

  return t_res3;
}


double get_single_emission_score_multi(struct hmm_multi_s *hmmp, struct letter_s *seq, struct letter_s *seq_2,
				       struct letter_s *seq_3, struct letter_s *seq_4, int c, int v, int replacement_letter_c,
				       int replacement_letter_c_2, int replacement_letter_c_3,
				       int replacement_letter_c_4, int use_labels, int a_index, int a_index_2,
				       int a_index_3, int a_index_4, int multi_scoring_method)
{
  double t_res3;
  double t_res3_1, t_res3_2, t_res3_3, t_res3_4;
  double cur_t_res3;
  int i,j;
  int cur_alphabet;
  double *cur_probs, *cur_emissions;
  int cur_a_size;
  int cur_a_index;
  int cur_replacement_letter;
  int cur_alphabet_type;
  struct letter_s *cur_seq;
  

  for(cur_alphabet = 1; cur_alphabet <= hmmp->nr_alphabets; cur_alphabet++) {
    if(cur_alphabet == 1) {
      cur_a_size = hmmp->a_size;
      cur_a_index = a_index;
      cur_probs = hmmp->replacement_letters->probs_1;
      cur_emissions = hmmp->emissions;
      cur_replacement_letter = replacement_letter_c;
      cur_alphabet_type = hmmp->alphabet_type;
      cur_seq = seq;
    }
    else if(cur_alphabet == 2) {
      cur_a_size = hmmp->a_size_2;
      cur_a_index = a_index_2;
      cur_probs = hmmp->replacement_letters->probs_2;
      cur_emissions = hmmp->emissions_2;
      cur_replacement_letter = replacement_letter_c_2;
      cur_alphabet_type = hmmp->alphabet_type_2;
      cur_seq = seq_2;
    }
    else if(cur_alphabet == 3) {
      cur_a_size = hmmp->a_size_3;
      cur_a_index = a_index_3;
      cur_probs = hmmp->replacement_letters->probs_3;
      cur_emissions = hmmp->emissions_3;
      cur_replacement_letter = replacement_letter_c_3;
      cur_alphabet_type = hmmp->alphabet_type_3;
      cur_seq = seq_3;
    }
    else if(cur_alphabet == 4) {
      cur_a_size = hmmp->a_size_4;
      cur_a_index = a_index_4;
      cur_probs = hmmp->replacement_letters->probs_4;
      cur_emissions = hmmp->emissions_4;
      cur_replacement_letter = replacement_letter_c_4;
      cur_alphabet_type = hmmp->alphabet_type_4;
      cur_seq = seq_4;
    }
    
    if(cur_alphabet_type == DISCRETE) {
      if(cur_replacement_letter == YES) {
	/* count emission prob with dot-product method */
	cur_t_res3 = 0.0;
	if(use_labels == YES && seq[c].label != *(hmmp->vertex_labels + v) && seq[c].label != '.') {
	  cur_t_res3 = 0.0;
	}
	else {
	  for(i = 0; i < cur_a_size; i++) {
	    cur_t_res3 += *(cur_probs + get_mtx_index(cur_a_index, i, cur_a_size)) *
	      *(cur_emissions + get_mtx_index(v,i,cur_a_size));
	  }
	}
      }
      else {
	if(use_labels == YES && seq[c].label != *(hmmp->vertex_labels + v) && seq[c].label != '.') {
	  cur_t_res3 = 0.0;
	}
	else {
	  cur_t_res3 = (*((cur_emissions) + (v * (cur_a_size)) + cur_a_index));
	}      
      }
    }
    else {
      if(use_labels == YES && seq[c].label != *(hmmp->vertex_labels + v) && seq[c].label != '.') {
	cur_t_res3 = 0.0;
      }
      else {
	t_res3 = 0.0;
	for(j = 0; j < cur_a_size / 3; j++) {
	  t_res3 += get_single_gaussian_statescore((*((cur_emissions) + (v * (cur_a_size)) + (j * 3))),
						   (*((cur_emissions) + (v * (cur_a_size)) + (j * 3 + 1))),
						   cur_seq[c].cont_letter) *
	    *((cur_emissions) + (v * (cur_a_size)) + (j * 3 + 2));
	}
      }      
    }
    
    if(cur_alphabet == 1) {
      t_res3_1 = cur_t_res3;
    }
    else if(cur_alphabet == 2) {
      t_res3_2 = cur_t_res3;
    }
    else if(cur_alphabet == 3) {
      t_res3_3 = cur_t_res3;
    }
    else if(cur_alphabet == 4) {
      t_res3_4 = cur_t_res3;
    }
  }
  
  /* calculate total res for this to-vertex (v) */
  if(multi_scoring_method == JOINT_PROB) {
    t_res3 = t_res3_1;
    if(hmmp->nr_alphabets > 1) {
      t_res3 = t_res3 * t_res3_2;
    }
    if(hmmp->nr_alphabets > 2) {
      t_res3 = t_res3 * t_res3_3;
    }
    if(hmmp->nr_alphabets > 3) {
      t_res3 = t_res3 * t_res3_4;
    }
  }

  return t_res3;
}
