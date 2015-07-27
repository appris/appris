#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>




#include "structs.h"
#include "funcs.h"

//#define DEBUG_LABELING_UPDATE
//#define DEBUG_DEALLOCATE_LABELINGS


#define REST_LETTER_INDEX 0.5

#define V_LIST_END  -99
#define V_LIST_NEXT  -9


double get_single_gaussian_statescore(double mu, double sigma_square, double letter)
{

  double res;

  if(sigma_square <= 0.0) {
    return 0.0;
  }
  else {
    res = exp(0.0 - (pow((letter-mu),2) / (2.0 * sigma_square))) / sqrt(sigma_square * 2.0 * 3.141592655);
    //printf("mu = %f, letter = %f, sigma_square = %f,  res = %f\n", mu, letter, sigma_square, res);
    return res;
  }
}

double get_dp_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
			 int p, double *emissions,  int vertex, int normalize, double *gap_shares)
{
  
  double seq_normalizer;
  double state_normalizer;
  double subst_mtx_normalizer;
  int a_index;
  double t_res_3;
  
  seq_normalizer = 0.0;
  state_normalizer = 0.0;
  subst_mtx_normalizer = 0.0;
  
  t_res_3 = 0.0;
  /* scoring using dot-product method */
  for(a_index = 0; a_index < a_size; a_index++) {

    if(use_prior_shares == YES) {
      t_res_3 += *(emissions + (vertex * a_size + a_index)) *
	(msa_seq + (p * (a_size+1) + a_index))->prior_share;
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->prior_share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
    else if(use_gap_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_size, a_size+1))->share == 1.0) {
	printf("Error: all gap column in sequence\n");
	exit(0);
      }
      t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share /
	(1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share);
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share /
			      (1.0 - *(gap_shares + p)), 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
    else {
      t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share;
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
  }
  if(normalize == YES) {
    seq_normalizer = sqrt(seq_normalizer);
    state_normalizer = sqrt(state_normalizer);
    if(t_res_3 != 0.0) {
      t_res_3 = t_res_3 / (seq_normalizer * state_normalizer);
    }
  }

  if(t_res_3 < 0.0) {
    printf("t_res_3 = %f\n", t_res_3);
    printf("Error: got strange dot product state result value\n");
    exit(0);
  }
  
  return t_res_3;
}


double get_dp_picasso_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
				 int p, double *emissions,  int vertex, int normalize, double *gap_shares, double *aa_freqs)
{
  
  double seq_normalizer;
  double state_normalizer;
  double subst_mtx_normalizer;
  int a_index;
  double t_res_3;
  
  seq_normalizer = 0.0;
  state_normalizer = 0.0;
  subst_mtx_normalizer = 0.0;
  
  t_res_3 = 0.0;
  /* scoring using dot-product method */
  for(a_index = 0; a_index < a_size; a_index++) {
    if(use_prior_shares == YES) {
      t_res_3 += *(emissions + (vertex * a_size + a_index)) *
	(msa_seq + (p * (a_size+1) + a_index))->prior_share /  *(aa_freqs + a_index);
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->prior_share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
    else if(use_gap_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_size, a_size+1))->share == 1.0) {
	printf("Error: all gap column in sequence\n");
	exit(0);
      }
      t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share /
	((1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share) *  *(aa_freqs + a_index));
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share /
			      (1.0 - *(gap_shares + p)), 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
    else {
      t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share /  *(aa_freqs + a_index);
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
  }
  if(normalize == YES) {
    seq_normalizer = sqrt(seq_normalizer);
    state_normalizer = sqrt(state_normalizer);
#ifdef DEBUG_BW
    printf("state_normalizer = %f\n", state_normalizer);
    printf("seq_normalizer = %f\n", seq_normalizer);
#endif     
    if(t_res_3 != 0.0) {
      t_res_3 = t_res_3 / (seq_normalizer * state_normalizer);
    }
  }

  if(t_res_3 < 0.0) {
    printf("t_res_3 = %f\n", t_res_3);
    printf("Error: got strange dot product state result value\n");
    exit(0);
  }
  
  return t_res_3;
}


double get_sjolander_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
				int p, double *emissions,  int vertex, int normalize, double *gap_shares)
{
  double seq_normalizer;
  double state_normalizer;
  double subst_mtx_normalizer;
  int a_index;
  double t_res_3;
  
  seq_normalizer = 0.0;
  state_normalizer = 0.0;
  subst_mtx_normalizer = 0.0;

  t_res_3 = 1.0;
  /* scoring using sjolander score method */

  for(a_index = 0; a_index < a_size; a_index++) {
    if(use_prior_shares == YES) {
      t_res_3 *= pow(*(emissions + (vertex * a_size + a_index)), 
		    (msa_seq + (p * (a_size+1) + a_index))->prior_share);
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->prior_share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
     
    }
    else if(use_gap_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_size, a_size+1))->share == 1.0) {
	printf("Error: all gap column in sequence\n");
	exit(0);
      }

      t_res_3 *= pow(*(emissions + get_mtx_index(vertex, a_index, a_size)),
		     (msa_seq + get_mtx_index(p,a_index, a_size+1))->share /
		     (1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share));
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share /
			      (1.0 - *(gap_shares + p)), 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
    else {
      t_res_3 *= pow(*(emissions + get_mtx_index(vertex, a_index, a_size)),
		     (msa_seq + get_mtx_index(p,a_index, a_size+1))->share);

      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
  }

  if(normalize == YES) {
    seq_normalizer = sqrt(seq_normalizer);
    state_normalizer = sqrt(state_normalizer);
#ifdef DEBUG_BW
    printf("state_normalizer = %f\n", state_normalizer);
    printf("seq_normalizer = %f\n", seq_normalizer);
#endif     
    if(t_res_3 != 0.0) {
      t_res_3 = t_res_3 / (seq_normalizer * state_normalizer);
    }
  }
  if(t_res_3 < 0.0) {
    printf("Error: got strange geometric mean state result value\n");
    printf("t_res_3 = %f\n", t_res_3);
    exit(0);
  }
  return t_res_3;
}



double get_sjolander_reversed_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
					 int p, double *emissions,  int vertex, int normalize, double *gap_shares)
{
  double seq_normalizer;
  double state_normalizer;
  double subst_mtx_normalizer;
  int a_index;
  double t_res_3;
  
  seq_normalizer = 0.0;
  state_normalizer = 0.0;
  subst_mtx_normalizer = 0.0;
  

  t_res_3 = 1.0;
  /* scoring using sjolander score method */
  for(a_index = 0; a_index < a_size; a_index++) {
    if(use_prior_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->prior_share != 0.0) {
      	t_res_3 *= pow((msa_seq + get_mtx_index(p,a_index, a_size+1))->prior_share,
      		       *(emissions + get_mtx_index(vertex, a_index, a_size)));
      }
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->prior_share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
    else if(use_gap_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_size, a_size+1))->share == 1.0) {
	printf("Error: all gap column in sequence\n");
	exit(0);
      }
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->share != 0.0) {
	t_res_3 *= pow((msa_seq + get_mtx_index(p,a_index, a_size+1))->share /
		       (1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share),
		       *(emissions + get_mtx_index(vertex, a_index, a_size)));
      }
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share /
			      (1.0 - *(gap_shares + p)), 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
    else {
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->share != 0.0) {
	t_res_3 *= pow((msa_seq + get_mtx_index(p,a_index, a_size+1))->share,
		       *(emissions + get_mtx_index(vertex, a_index, a_size)));
      }
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
      }
    }
  }
  if(normalize == YES) {
    seq_normalizer = sqrt(seq_normalizer);
    state_normalizer = sqrt(state_normalizer);
#ifdef DEBUG_BW
    printf("state_normalizer = %f\n", state_normalizer);
    printf("seq_normalizer = %f\n", seq_normalizer);
#endif     
    if(t_res_3 != 0.0) {
      t_res_3 = t_res_3 / (seq_normalizer * state_normalizer);
    }
  }
  if(t_res_3 < 0.0) {
    printf("Error: got strange geometric mean state result value\n");
    printf("t_res_3 = %f\n", t_res_3);
    exit(0);
  }
  
  return t_res_3;
}


double get_picasso_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
			      int p, double *emissions,  int vertex, int normalize, double *gap_shares, double *aa_freqs)
{
  
  double seq_normalizer;
  double state_normalizer;
  double subst_mtx_normalizer;
  int a_index;
  double t_res_3;
  
  seq_normalizer = 0.0;
  state_normalizer = 0.0;
  subst_mtx_normalizer = 0.0;
  

  t_res_3 = 1.0;
  /* scoring using picasso-product method */

  for(a_index = 0; a_index < a_size; a_index++) {
    if(use_prior_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->prior_share != 0.0  &&
	 *(emissions + get_mtx_index(vertex, a_index, a_size)) != SILENT) {
	t_res_3 *= pow((msa_seq + get_mtx_index(p,a_index, a_size+1))->prior_share / *(aa_freqs + a_index),
		       *(emissions + (vertex * a_size + a_index)));
      }
    }
    else if(use_gap_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_size, a_size+1))->share == 1.0) {
	printf("Error: all gap column in sequence\n");
	exit(0);
      }
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->share != 0.0  &&
	 *(emissions + get_mtx_index(vertex, a_index, a_size)) != SILENT) {
	t_res_3 *= pow(((msa_seq + get_mtx_index(p,a_index, a_size+1))->share /
			(1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share)) / *(aa_freqs + a_index),
		       *(emissions + get_mtx_index(vertex, a_index, a_size)));
      }
    }
    else {
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->share != 0.0 &&
	 *(emissions + get_mtx_index(vertex, a_index, a_size)) != SILENT) {
	t_res_3 *= pow((msa_seq + get_mtx_index(p,a_index, a_size+1))->share / *(aa_freqs + a_index),
		       *(emissions + get_mtx_index(vertex, a_index, a_size)));
      }
    }
  }
  
  if(t_res_3 < 0.0 || t_res_3 > 1000000000000.0) {
    printf("Error: got strange picasso product state result value\n");
    exit(0);
  }
  
  return t_res_3;
}


double get_picasso_sym_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
				  int p, double *emissions,  int vertex, int normalize, double *gap_shares, double *aa_freqs)
{
  
  double seq_normalizer;
  double state_normalizer;
  double subst_mtx_normalizer;
  int a_index;
  double t_res_3;
  
  seq_normalizer = 0.0;
  state_normalizer = 0.0;
  subst_mtx_normalizer = 0.0;
  

  t_res_3 = 1.0;
  /* scoring using picasso-product method */

  for(a_index = 0; a_index < a_size; a_index++) {
    if(use_prior_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->prior_share != 0.0  &&
	 *(emissions + get_mtx_index(vertex, a_index, a_size)) != SILENT) {
	t_res_3 *= pow((msa_seq + get_mtx_index(p,a_index, a_size+1))->prior_share / *(aa_freqs + a_index),
		       *(emissions + (vertex * a_size + a_index))) *
	  pow(*(emissions + (vertex * a_size + a_index)) / *(aa_freqs + a_index),
	      (msa_seq + get_mtx_index(p,a_index, a_size+1))->prior_share);
      }
    }
    else if(use_gap_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_size, a_size+1))->share == 1.0) {
	printf("Error: all gap column in sequence\n");
	exit(0);
      }
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->share != 0.0  &&
	 *(emissions + get_mtx_index(vertex, a_index, a_size)) != SILENT) {
	t_res_3 *= pow(((msa_seq + get_mtx_index(p,a_index, a_size+1))->share /
			(1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share)) / *(aa_freqs + a_index),
		       *(emissions + get_mtx_index(vertex, a_index, a_size))) *
	  pow(*(emissions + (vertex * a_size + a_index)) / *(aa_freqs + a_index),
	      (msa_seq + get_mtx_index(p,a_index, a_size+1))->share /
	      (1.0 - (msa_seq + get_mtx_index(p,a_size, a_size+1))->share));
      }
    }
    else {
      //printf(" *(aa_freqs + a_index) = %f\n",  *(aa_freqs + a_index));
      //printf("(msa_seq + (p * (a_size+1) + a_index))->share = %f\n", (msa_seq + get_mtx_index(p,a_index, a_size+1))->share);
      //printf("*(emissions + get_mtx_index(vertex, a_index, a_size)) = %f\n",*(emissions + get_mtx_index(vertex, a_index, a_size))); 
      if((msa_seq + get_mtx_index(p,a_index, a_size+1))->share != 0.0 &&
	 *(emissions + get_mtx_index(vertex, a_index, a_size)) != SILENT) {
	t_res_3 *= pow((msa_seq + get_mtx_index(p,a_index, a_size+1))->share / *(aa_freqs + a_index),
		       *(emissions + get_mtx_index(vertex, a_index, a_size))) *
	  pow(*(emissions + (vertex * a_size + a_index)) / *(aa_freqs + a_index),
	      (msa_seq + get_mtx_index(p,a_index, a_size+1))->share);
      }
    }
  }

  if(t_res_3 < 0.0 || t_res_3 > 1000000000000.0) {
    printf("Error: got strange picasso product state result value\n");
    exit(0);
  }
  
  return t_res_3;
}


double get_subst_mtx_product_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
					int p, double *emissions, int vertex, double *subst_mtx)
{
  int a_index, a_index2;
  double t_res_3;
  
  t_res_3 = 0.0;
  for(a_index = 0; a_index < a_size; a_index++) {
    for(a_index2 = 0; a_index2 < a_size; a_index2++) {
      if(use_gap_shares == YES) {
	t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	  (msa_seq + get_mtx_index(p,a_index2, a_size+1))->share /
	  (1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share) *
	  *(subst_mtx + (a_index * a_size + a_index2));
      }
      else {
	t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	  (msa_seq + get_mtx_index(p,a_index2, a_size+1))->share /
	  *(subst_mtx + (a_index * a_size + a_index2));
      }
    }
  }
  
  if(t_res_3 < 0.0) {
    printf("Error: got strange subst mtx product state result value\n");
    exit(0);
  }
  
  return t_res_3;
}


double get_subst_mtx_dot_product_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
					    int p, double *emissions,  int vertex, int normalize, double *gap_shares,
					    int query_index, double *subst_mtx)
{
  double seq_normalizer;
  double state_normalizer;
  double subst_mtx_normalizer;
  int a_index;
  double t_res_3;
  
  seq_normalizer = 0.0;
  state_normalizer = 0.0;
  subst_mtx_normalizer = 0.0;
  
  t_res_3 = 0.0;
  /* scoring using subst_mtx_dot-product method */
  for(a_index = 0; a_index < a_size; a_index++) {
    if(use_prior_shares == YES) {
      t_res_3 += *(emissions + (vertex * a_size + a_index)) *
	(msa_seq + (p * (a_size+1) + a_index))->prior_share *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size));
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->prior_share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
	subst_mtx_normalizer +=  pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)), 2);
      }
    }
    else if(use_gap_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_size, a_size+1))->share == 1.0) {
	printf("Error: all gap column in sequence\n");
	  exit(0);
      }
      t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size)) /
	(1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share);
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share /
			      (1.0 - *(gap_shares + p)), 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
	subst_mtx_normalizer +=  pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)), 2);
      }
    }
    else {
      t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size));
      if(normalize == YES) {
	seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share, 2);
	state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
	subst_mtx_normalizer +=  pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)), 2);
      }
    }
  }
  
  if(normalize == YES) {
    seq_normalizer = sqrt(seq_normalizer);
    state_normalizer = sqrt(state_normalizer);
    subst_mtx_normalizer = sqrt(subst_mtx_normalizer);
#ifdef DEBUG
    printf("state_normalizer = %f\n", state_normalizer);
      printf("seq_normalizer = %f\n", seq_normalizer);
#endif      
      if(t_res_3 != 0.0) {
	t_res_3 = t_res_3 / (seq_normalizer * state_normalizer * subst_mtx_normalizer);
      }
  }

  if(t_res_3 < 0.0) {
    printf("Error: got strange subst mtx dot product state result value\n");
    exit(0);
  }
  
  return t_res_3;
}


double get_subst_mtx_dot_product_prior_statescore(int a_size, int use_gap_shares, int use_prior_shares, struct msa_letter_s *msa_seq,
						  int p, double *emissions,  int vertex, int normalize, double *gap_shares,
						  int query_index, double *subst_mtx)
{
  double seq_normalizer;
  double state_normalizer;
  double subst_mtx_normalizer;
  int a_index;
  double t_res_3;
  double rest_share, default_share;
  
  seq_normalizer = 0.0;
  state_normalizer = 0.0;
  subst_mtx_normalizer = 0.0;
  default_share = 1.0 / (double)(a_size);

  t_res_3 = 0.0;
  
  /* scoring using dot-product method */
  rest_share = 1.0;
  for(a_index = 0; a_index < a_size; a_index++) {
    if(use_prior_shares == YES) {
      t_res_3 += *(emissions + (vertex * a_size + a_index)) *
	(msa_seq + (p * (a_size+1) + a_index))->prior_share *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size));
      if(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)) != 0.0) {
	rest_share = rest_share - (msa_seq + (p * (a_size+1) + a_index))->prior_share;
	if(normalize == YES) {
	  seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->prior_share, 2);
	  state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
	  subst_mtx_normalizer +=  pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)), 2);
	}
      }
    }
    else if(use_gap_shares == YES) {
      if((msa_seq + get_mtx_index(p,a_size, a_size+1))->share == 1.0) {
	printf("Error: all gap column in sequence\n");
	exit(0);
      }
      t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size)) /
	(1.0 -(msa_seq + get_mtx_index(p,a_size, a_size+1))->share);
      if(*(subst_mtx + get_mtx_index(a_index, a_index, a_size)) != 0.0) {
	rest_share = rest_share - (msa_seq + (p * (a_size+1) + a_index))->share / 
	  (1.0 - *(gap_shares + p));
	if(normalize == YES) {
	  seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share /
				(1.0 - *(gap_shares + p)), 2);
	  state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
	  subst_mtx_normalizer +=  pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)), 2);
	}
      }
      
    }
    else {
      t_res_3 += *(emissions + get_mtx_index(vertex, a_index, a_size)) *
	(msa_seq + get_mtx_index(p,a_index, a_size+1))->share *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size));
      if(*(subst_mtx + get_mtx_index(a_index, a_index, a_size)) != 0.0) {
	rest_share = rest_share - (msa_seq + (p * (a_size+1) + a_index))->share;
	if(normalize == YES) {
	  seq_normalizer += pow((msa_seq + (p * (a_size+1) + a_index))->share, 2);
	  state_normalizer += pow(*(emissions + (vertex * a_size + a_index)), 2);
	  subst_mtx_normalizer +=  pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)), 2);
	}
      }
    }
  }
  if(rest_share < 0.0) {
    rest_share = 0.0;
  }
  t_res_3 += default_share * rest_share;
  seq_normalizer += pow(rest_share, 2);
  state_normalizer += pow(default_share, 2);
  
  if(normalize == YES) {
    seq_normalizer = sqrt(seq_normalizer);
    state_normalizer = sqrt(state_normalizer);
    subst_mtx_normalizer = sqrt(subst_mtx_normalizer);
#ifdef DEBUG_BW
    printf("state_normalizer = %f\n", state_normalizer);
    printf("seq_normalizer = %f\n", seq_normalizer);
#endif
    if(t_res_3 != 0.0) {
      t_res_3 = t_res_3 / (seq_normalizer * state_normalizer * subst_mtx_normalizer);
    }
  }

  if(t_res_3 < 0.0) {
    printf("Error: got strange subst mtx dot product prior state result value\n");
    exit(0);
  }
  
  return t_res_3;
}



/************************************* add to E methods *********************************************/
void add_to_E_continuous(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
			 int k, int a_size, double *emissions)
{
  double mean_value, varians;
  int i,j;
  double continuous_score_all, continuous_score_j, gamma_p_j;

  mean_value = (msa_seq + get_mtx_index(p,0,a_size+1))->share;
  
  
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
    varians = pow((msa_seq + get_mtx_index(p,0,a_size+1))->share - *(emissions + get_mtx_index(k, j * 3, a_size)), 2);
    if(continuous_score_all > 0.0) {
      gamma_p_j = Eka_base * continuous_score_j / continuous_score_all;
    }
    else {
      gamma_p_j = 0.0;
    }
    *(E + get_mtx_index(k, j * 3, a_size + 1)) += mean_value * gamma_p_j;
    *(E + get_mtx_index(k, j * 3 + 1, a_size + 1)) += varians * gamma_p_j;
    *(E + get_mtx_index(k, j * 3 + 2, a_size + 1)) += gamma_p_j;
  }
  
  *(E + get_mtx_index(k, j * 3, a_size + 1)) += Eka_base;
}



void add_to_E_dot_product(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
			  int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  double prf_column_length;

  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) / prf_column_length;
    }
  }
}

void add_to_E_dot_product_picasso(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
				  int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  double prf_column_length;
  
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) / prf_column_length;
    }
  }
}

void add_to_E_picasso(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
		      int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  double prf_column_length;

  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share);
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) / prf_column_length;
    }
  }
}

void add_to_E_picasso_sym(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
		      int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  double prf_column_length;

  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share);
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) / prf_column_length;
    }
  }
}

void add_to_E_sjolander_score(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
			      int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  double prf_column_length;


  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share);
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) / prf_column_length;
    }
  }
}

void add_to_E_sjolander_reversed_score(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
				       int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  double prf_column_length;
  
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) / prf_column_length;
    }
  }
}

void add_to_E_dot_product_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
				 int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  double prf_column_length;

  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) / prf_column_length;
    }
  }
}

void add_to_E_dot_product_picasso_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
					 int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  double prf_column_length;

  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) / prf_column_length;
    }
  }
}

void add_to_E_picasso_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
				 int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  double prf_column_length;

  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) / prf_column_length;
    }
  }
}

void add_to_E_picasso_sym_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
				 int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  double prf_column_length;

  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) / prf_column_length;
    }
  }
}

void add_to_E_sjolander_score_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
				     int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  double prf_column_length;
  
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) / prf_column_length;
    }
  }
}

void add_to_E_sjolander_reversed_score_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
					      int k, int a_size, int normalize)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  double prf_column_length;
  
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences);
      
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) / prf_column_length;
    }
  }
}

void add_to_E_subst_mtx_product(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
			  int k, int a_size, int normalize, double *subst_mtx)
{
  int a_index, a_index2;
  double prf_column_length;
  double subst_mtx_row_length;
  
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      for(a_index2 = 0; a_index2 < a_size; a_index2++) {
	*(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	  (double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) *
	  *(subst_mtx + get_mtx_index(a_index, a_index2, a_size));
      }
    }
  }
  else {
    printf("Error: no normalizing in subst_mtx_product, yet...\n");
    exit(0);
  }

}
void add_to_E_subst_mtx_product_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
				       int k, int a_size, int normalize, double *subst_mtx)
{
  int a_index, a_index2;
  double prf_column_length;
  double subst_mtx_row_length;
  
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      for(a_index2 = 0; a_index2 < a_size; a_index2++) {
	*(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	  (double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) *
	  *(subst_mtx + get_mtx_index(a_index, a_index2, a_size));
      }
    }
  }
  else {
    printf("Error: no normalizing in subst_mtx_product, yet...\n");
    exit(0);
  }

}

void add_to_E_subst_mtx_dot_product(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
			  int k, int a_size, int normalize, double *subst_mtx, char *alphabet)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  double prf_column_length, subst_mtx_row_length;
  int query_index;
  
  query_index = get_alphabet_index_msa_query((msa_seq + (p * (a_size+1)))->query_letter, alphabet, a_size);
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size));
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    subst_mtx_row_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      subst_mtx_row_length += pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)),2);
    }
    subst_mtx_row_length = sqrt(subst_mtx_row_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size))/ 
	(prf_column_length * subst_mtx_row_length);
    }
  }
}

void add_to_E_subst_mtx_dot_product_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
					   int k, int a_size, int normalize, double *subst_mtx, char *alphabet)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  double prf_column_length, subst_mtx_row_length;
  int query_index;
  
  query_index = get_alphabet_index_msa_query((msa_seq + (p * (a_size+1)))->query_letter, alphabet, a_size);
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size));
    }
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    subst_mtx_row_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      subst_mtx_row_length += pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)),2);
    }
    subst_mtx_row_length = sqrt(subst_mtx_row_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size))/ 
	(prf_column_length * subst_mtx_row_length);
    }
  }
}

void add_to_E_subst_mtx_dot_product_prior(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
					  int k, int a_size, int normalize, double *subst_mtx, char *alphabet)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  double prf_column_length, subst_mtx_row_length;
  int query_index;
  double rli;
  
  rli = REST_LETTER_INDEX;
  
  query_index = get_alphabet_index_msa_query((msa_seq + (p * (a_size+1)))->query_letter, alphabet, a_size);
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size));

      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) * rli;
    }
    
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    subst_mtx_row_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      subst_mtx_row_length += pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)),2);
    }
    subst_mtx_row_length = sqrt(subst_mtx_row_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size))/ 
	(prf_column_length * subst_mtx_row_length);

      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->share) * rli / 
	(prf_column_length * subst_mtx_row_length);
    }
  }
}


void add_to_E_subst_mtx_dot_product_prior_nr_occ(double *E, double Eka_base, struct msa_letter_s *msa_seq, int p,
						 int k, int a_size, int normalize, double *subst_mtx, char *alphabet)
{
  /* k is a state index, p is a sequence position index */
  
  int a_index;
  double prf_column_length, subst_mtx_row_length;
  int query_index;
  double rli;
  
  rli = REST_LETTER_INDEX;
  
  query_index = get_alphabet_index_msa_query((msa_seq + (p * (a_size+1)))->query_letter, alphabet, a_size);
  if(normalize == NO) {
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size));

      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) * rli;
    }
    
  }
  else {
    prf_column_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      prf_column_length += pow(((msa_seq + get_mtx_index(p, a_index, a_size+1))->share),2);
    }
    prf_column_length = sqrt(prf_column_length);
    subst_mtx_row_length = 0.0;
    for(a_index = 0; a_index < a_size; a_index++) {
      subst_mtx_row_length += pow(*(subst_mtx + get_mtx_index(query_index, a_index, a_size)),2);
    }
    subst_mtx_row_length = sqrt(subst_mtx_row_length);
    for(a_index = 0; a_index < a_size; a_index++) {
      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) *
	*(subst_mtx + get_mtx_index(query_index, a_index, a_size))/ 
	(prf_column_length * subst_mtx_row_length);

      *(E + get_mtx_index(k, a_index, a_size)) += Eka_base *
	(double)((msa_seq + get_mtx_index(p, a_index, a_size+1))->nr_occurences) * rli / 
	(prf_column_length * subst_mtx_row_length);
    }
  }
}



/* General versions of the functions needed for keeping track of the labeleings in the one-best algorithm */

void update_labelings(struct one_best_s *cur_rowp, char *vertex_labels, 
		      int *sorted_v_list, int seq_len, int c, char *labels, int nr_of_labels, int nr_v)
{
  int v,w;
  int v_list_index;
  int cur_address;
  int first;
  char *tmp_labeling;
  char cur_label;
  int **same_labeling_lists;
  int *same_labeling_list_indices;
  int i;
  

#ifdef DEBUG_LABELING_UPDATE
  dump_v_list(sorted_v_list);
  printf("nr of labels = %d\n", nr_of_labels);
  printf("dump of nr three\n");
#endif

  v_list_index = 0;

  same_labeling_list_indices = (int*)(malloc_or_die(nr_of_labels * sizeof(int)));
  
  same_labeling_lists = (int**)(malloc_or_die(nr_of_labels * sizeof(int*)));


  for(i = 0; i < nr_of_labels; i++) {
    same_labeling_lists[i] = (int*)(malloc_or_die((nr_v * 2 + 1) * sizeof(int)));
  }
 

  for(v = 0; v < nr_v; v++) {
    (cur_rowp + v)->is_updated = NO;
  }
  
  for(i = 0; i < nr_of_labels; i++) {
    same_labeling_list_indices[i] = 0;
  }


  for(v = 1; v < nr_v-1; v++) {
    /* find all states with same labeling as this state up to previous row */
    
    if((cur_rowp + v)->is_updated == NO && (cur_rowp + v)->labeling != NULL) {
      cur_address = (int)((cur_rowp+v)->labeling);
#ifdef DEBUG_LABELING_UPDATE
      printf("searching vertex %d\n", v);
#endif
      for(i = 0; i < nr_of_labels; i++) {
	if(*(vertex_labels + v) == *(labels + i)) {
	  *(*(same_labeling_lists + i) + same_labeling_list_indices[i]) = v;
	  same_labeling_list_indices[i] += 1;
	  break;
	}
      }
      
      (cur_rowp+v)->is_updated = YES;
      for(w = v+1; w < nr_v-1; w++) {
	if((int)((cur_rowp+w)->labeling) == cur_address && (cur_rowp + w)->is_updated == NO) {
#ifdef DEBUG_LABELING_UPDATE	  
	  printf("found same address, vertex nr = %d\n", w);
#endif
	  for(i = 0; i < nr_of_labels; i++) {
	    if(*(vertex_labels + w) == *(labels + i)) {
	      *(*(same_labeling_lists + i) + same_labeling_list_indices[i]) = w;
	      same_labeling_list_indices[i] += 1;
	      break;
	    }
	  }

	  (cur_rowp+w)->is_updated = YES;
	}
      }

      for(i = 0; i < nr_of_labels; i++) {
	*(*(same_labeling_lists + i) + same_labeling_list_indices[i]) = END;
	same_labeling_list_indices[i] += 1;
      }
    }
  }
  for(i = 0; i < nr_of_labels; i++) {
    *(*(same_labeling_lists + i) + same_labeling_list_indices[i]) = TOT_END;
  }


#ifdef DEBUG_LABELING_UPDATE
  for(i = 0; i < nr_of_labels; i++) {
    printf("same_labeling_lists, label: %c\n", labels[i]);
    dump_label_tmp_list(*(same_labeling_lists + i));
  }
  //exit(0);
#endif
  
  for(i = 0; i < nr_of_labels; i++) {
    same_labeling_list_indices[i] = 0;
    while(*(*(same_labeling_lists + i) + same_labeling_list_indices[i]) != TOT_END) {
      first = YES;
      while(*(*(same_labeling_lists + i) + same_labeling_list_indices[i]) != END) {
	/* update sorted_v_list */
	*(sorted_v_list + v_list_index) = *(*(same_labeling_lists + i) + same_labeling_list_indices[i]);
	v_list_index++;

	/* update pointers and label paths */
	if(first == YES) {
	  tmp_labeling = (cur_rowp + *(*(same_labeling_lists + i) + same_labeling_list_indices[i]))->labeling;
	  (cur_rowp + *(*(same_labeling_lists + i) + same_labeling_list_indices[i]))->labeling =
	    (char*)malloc_or_die((c+1) * sizeof(char));
	  memcpy((cur_rowp + *(*(same_labeling_lists + i) + same_labeling_list_indices[i]))->labeling,
		 tmp_labeling, (c) * sizeof(char));
	  ((cur_rowp + *(*(same_labeling_lists + i) + same_labeling_list_indices[i]))->labeling)[c] = labels[i];
#ifdef DEBUG_LABELING_UPDATE	 
	  printf("added label;labels[%d] = %c\n", i, labels[i]);
#endif
	  first = NO;
	  tmp_labeling = (cur_rowp + *(*(same_labeling_lists + i) + same_labeling_list_indices[i]))->labeling;
	}
	else {
	  (cur_rowp + *(*(same_labeling_lists + i) + same_labeling_list_indices[i]))->labeling = tmp_labeling;
	}
#ifdef DEBUG_LABELING_UPDATE
      printf("label length c = %d\n", c);
      dump_labeling((cur_rowp + *(*(same_labeling_lists + i) + same_labeling_list_indices[i]))->labeling, c);
#endif


	same_labeling_list_indices[i] += 1;
      }
      if(first == NO) {
	*(sorted_v_list + v_list_index) = V_LIST_NEXT;
	v_list_index++;
      }
      same_labeling_list_indices[i] += 1;
    }
  }
  
  *(sorted_v_list + v_list_index) = V_LIST_END;
  
  for(v = 1; v < nr_v; v++) {
    (cur_rowp + v)->is_updated = NO;
  }
  

  /* garbage collection */
  for(i = 0; i < nr_of_labels; i++) {
    free(same_labeling_lists[i]);
  }
  free(same_labeling_list_indices);
  free(same_labeling_lists);
}

void deallocate_row_labelings(struct one_best_s *prev_rowp, int nr_v)
{
  int dealloc_index;
  int cur_address;
  int v,w;
  int *dealloc_list;
  
#ifdef DEBUG_DEALLOCATE_LABELINGS
  printf("starting dealloc\n");
  printf("nr_v = %d\n", nr_v);
#endif

  dealloc_list = (int*)(malloc_or_die((nr_v+1) * sizeof(int)));

  for(v = 0; v < nr_v; v++) {
    (prev_rowp + v)->is_updated = NO;
  }
  
  
  /* deallocate last row's labelings */
  dealloc_index = 0;
  for(v = 0; v < nr_v; v++) {
    /* find all states with same labeling as this state up to previous row */
    if((prev_rowp + v)->is_updated == NO && (prev_rowp + v)->labeling != NULL) {
      cur_address = (int)((prev_rowp+v)->labeling);
      dealloc_list[dealloc_index] = v;
      dealloc_index++;
      (prev_rowp+v)->is_updated = YES;
      for(w = v+1; w < nr_v; w++) {
	if((int)((prev_rowp+w)->labeling) == cur_address) {
#ifdef DEBUG_DEALLOCATE_LABELINGS
	  printf("found same address, vertices %d and %d: %x\n", v, w, (prev_rowp + w)->labeling);
#endif	
	  (prev_rowp+w)->is_updated = YES;
	  (prev_rowp+w)->labeling = NULL;
	}
      }
    }
  }
  dealloc_list[dealloc_index] = END;
  
  for(dealloc_index = 0; dealloc_list[dealloc_index] != END; dealloc_index++) {
#ifdef DEBUG_DEALLOCATE_LABELINGS    
    printf("dealloc_index = %d\n", dealloc_index);
    printf("freeing labeling of vertex %d\n", dealloc_list[dealloc_index]);
#endif
    free((prev_rowp + dealloc_list[dealloc_index])->labeling);
#ifdef DEBUG_DEALLOCATE_LABELINGS    
    printf("done\n");
#endif
  } 
  
  free(dealloc_list);
#ifdef DEBUG_DEALLOCATE_LABELINGS  
  printf("exiting dealloc\n");
#endif
}


