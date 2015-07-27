#include <stdlib.h>
#include <stdio.h>



#include "structs.h"


void dump_align_matrix(int nr_rows, int nr_cols, struct align_mtx_element_s *mtx)
{
  int i,j;
  struct align_mtx_element_s *mtx_tmp;
  mtx_tmp = mtx;
  
  printf("alignment matrix dump:\n");
  printf("nr rows: %d\n", nr_rows);
  printf("nr columns: %d\n", nr_cols);
  printf("scores:\n");
  for(i = 0; i < nr_rows; i++) {
    for(j = 0; j < nr_cols; j++) {
      printf("%d  ", mtx_tmp->score);
      mtx_tmp++;
    }
    printf("\n");
  }
  printf("\n\n");

  mtx_tmp = mtx;

  printf("lasts:\n");
  for(i = 0; i < nr_rows; i++) {
    for(j = 0; j < nr_cols; j++) {
      printf("%c  ", mtx->last);
      mtx++;
    }
    printf("\n");
  }
  printf("\n");
}

void dump_trans_matrix(int nr_rows, int nr_cols, double *mtx)
{
  int i,j;
  printf("transition matrix dump:\n");
  printf("nr rows: %d\n", nr_rows);
  printf("nr columns: %d\n", nr_cols);
  for(i = 0; i < nr_rows; i++) {
    for(j = 0; j < nr_cols; j++) {
      printf("%f  ", *mtx);
      mtx++;
    }
    printf("\n");
  }
  printf("\n");
}

void dump_int_trans_matrix(int nr_rows, int nr_cols, int *mtx)
{
  int i,j;
  printf("transition matrix dump:\n");
  printf("nr rows: %d\n", nr_rows);
  printf("nr columns: %d\n", nr_cols);
  for(i = 0; i < nr_rows; i++) {
    for(j = 0; j < nr_cols; j++) {
      printf("%d  ", *mtx);
      mtx++;
    }
    printf("\n");
  }
  printf("\n");
}

void dump_emiss_matrix(int nr_rows, int nr_cols, double *mtx)
{
  int i,j;
  printf("emission matrix dump:\n");
  for(i = 0; i < nr_rows; i++) {
    for(j = 0; j < nr_cols; j++) {
      printf("%f  ", *mtx);
      mtx++;
    }
    printf("\n");
  }
  printf("\n");
}

void dump_post_prob_matrix(int nr_rows, int nr_cols, double *mtx)
{
  int i,j;
  printf("posterior probability matrix dump:\n");
  for(i = 0; i < nr_rows; i++) {
    for(j = 0; j < nr_cols; j++) {
      printf("%f  ", *mtx);
      mtx++;
    }
    printf("\n");
  }
  printf("\n");
}

void dump_forward_matrix(int nr_rows, int nr_cols, struct forward_s *mtx)
{
  int i,j;
  printf("forward matrix dump:\n");
  for(i = 0; i < nr_rows; i++) {
    for(j = 0; j < nr_cols; j++) {
      printf("%f  ", mtx->prob);
      mtx++;
    }
    printf("\n");
  }
  printf("\n");
}

void dump_viterbi_matrix(int nr_rows, int nr_cols, struct viterbi_s *mtx)
{
  int i,j;
  struct viterbi_s *mtx_2;
  
  mtx_2 = mtx;
  printf("viterbi matrix dump probs:\n");
  for(i = 0; i < nr_rows; i++) {
    for(j = 0; j < nr_cols; j++) {
      //printf("%d  ", mtx->prev);
       printf("%f  ", mtx->prob);
      mtx++;
    }
    printf("\n");
  }
  printf("\n");

  printf("viterbi matrix dump prevs:\n");
  for(i = 0; i < nr_rows; i++) {
    printf("row: %d\n", i);
    for(j = 0; j < nr_cols; j++) {
      printf("col: %d = %d  \n", j, mtx_2->prev);
      //printf("%f  ", mtx->prob);
      mtx_2++;
    }
    printf("nr_cols = %d\n", j);
    printf("\n");
  }
  printf("\n");
}

void dump_one_best_matrix(int nr_rows, int nr_cols, struct one_best_s *mtx)
{
  int i,j;
  printf("one best matrix dump:\n");
  printf("%4d  ", 0);
  for(j = 0; j < nr_cols; j++) {
      printf("%8d  ", j);
  }
  printf("\n");
  for(i = 0; i < nr_rows; i++) {
    printf("%4d  ", i);
    for(j = 0; j < nr_cols; j++) {
      printf("%f  ", mtx->prob);
      mtx++;
    }
    printf("\n");
  }
  printf("\n");
}

void dump_backward_matrix(int nr_rows, int nr_cols, struct backward_s *mtx)
{
  int i,j;
  printf("backward matrix dump:\n");
  for(i = 0; i < nr_rows; i++) {
    for(j = 0; j < nr_cols; j++) {
      printf("%f  ", mtx->prob);
      mtx++;
    }
    printf("\n");
  }
  printf("\n");
}

void dump_scaling_array(int len, double *mtx)
{
  int i;
  printf("scaling array dump\n");
  for(i = 0; i < len; i++) {
    printf("%d: %f ", i, *mtx);
    mtx++;
  }
  printf("\n");
}

void dump_from_trans_array(int nr_v, struct path_element **array)
{
  int i;
  struct path_element *cur;
  

  printf("\nfrom_trans_array_dump:\n");
  printf("nr_v = %d\n", nr_v);
  
  
  for(i = 0; i < nr_v; i++) {
    if((*(array + i))->vertex == END) {
      printf("no transitions to vertex %d\n", i);
    }
    else {
      printf("paths to vertex %d:\n", i);
      cur = *(array + i);
      while(cur->vertex != END) {
	printf("%d ", cur->vertex);
	while(cur->next != NULL) {
	  cur = cur->next;
	  printf("%d ", cur->vertex);
	}
	printf("%d ", i);
	printf("\n");
	cur++;
      }
    }
  }
  printf("\n");
}

void dump_subst_mtx(double *mtx, int a_size)
{
  int i,j;
  int nr_rows = a_size + 1;
  int nr_cols = a_size;
  printf("substitution matrix dump:\n");
  printf("nr rows: %d\n", nr_rows);
  printf("nr columns: %d\n", nr_cols);
  for(i = 0; i < nr_rows; i++) {
    for(j = 0; j < nr_cols; j++) {
      printf("%f  ", *mtx);
      mtx++;
    }
    printf("\n");
  }
  printf("\n");
}

void dump_to_trans_array(int nr_v, struct path_element **array)
{
  int i;
  struct path_element *cur;
  

  printf("\nto_trans_array_dump:\n");
  printf("nr_v = %d\n", nr_v);
  
  for(i = 0; i < nr_v; i++) {
    if((*(array + i))->vertex == END) {
      printf("no transitions from vertex %d\n", i);
    }
    else {
      printf("paths from vertex %d:\n", i);
      cur = *(array + i);
      while(cur->vertex != END) {
	printf("%d ", i);
	printf("%d ", cur->vertex);
	while(cur->next != NULL) {
	  cur = cur->next;
	  printf("%d ", cur->vertex);
	}
	printf("\n");
	cur++;
      }
    }
  }
  printf("\n");
}


void dump_T_matrix(int nr_rows, int nr_cols, double *mtx)
{
  int i,j;
  printf("T matrix dump:\n");
  printf("nr rows: %d\n", nr_rows);
  printf("nr columns: %d\n", nr_cols);
  for(i = 0; i < nr_rows; i++) {
    for(j = 0; j < nr_cols; j++) {
      printf("%f  ", *mtx);
      mtx++;
    }
    printf("\n");
  }
  printf("\n");
}

void dump_E_matrix(int nr_rows, int nr_cols, double *mtx)
{
  int i,j;
  printf("E matrix dump:\n");
  printf("nr rows: %d\n", nr_rows);
  printf("nr columns: %d\n", nr_cols);
  for(i = 0; i < nr_rows; i++) {
    for(j = 0; j < nr_cols; j++) {
      printf("%f  ", *mtx);
      mtx++;
    }
    printf("\n");
  }
  printf("\n");
}

/* dumps most probable state path for given string (misses end state) */
void dump_viterbi_path(struct viterbi_s *cur, struct hmm_s *hmmp,
		       struct viterbi_s *viterbi_mtxp, int row, int row_size)
{
  struct path_element *p_el;

  if(cur->prev == 0) {
    p_el = cur->prevp;
    printf("0 ");
    while(p_el->next != NULL) {
      printf("%d ", p_el->vertex);
      p_el++;
    }
  }
  else {
    dump_viterbi_path(viterbi_mtxp + get_mtx_index(row-1, cur->prev, row_size), hmmp,
		      viterbi_mtxp, row-1, row_size);
    p_el = cur->prevp;
    printf("%d ", cur->prev);
    while(p_el->next != NULL) {
      p_el++;
      printf("%d ", p_el->vertex);
    }
  }
}

void dump_viterbi_label_path(struct viterbi_s *cur, struct hmm_s *hmmp,
			     struct viterbi_s *viterbi_mtxp, int row, int row_size)
{
  struct path_element *p_el;

  if(cur->prev == 0) {
    p_el = cur->prevp;
    printf("0 ");
    while(p_el->next != NULL) {
      printf("%c ", *(hmmp->vertex_labels + p_el->vertex));
      p_el++;
    }
  }
  else {
    dump_viterbi_label_path(viterbi_mtxp + get_mtx_index(row-1, cur->prev, row_size), hmmp,
		      viterbi_mtxp, row-1, row_size);
    p_el = cur->prevp;
    printf("%c ", *(hmmp->vertex_labels + (int)(cur->prev)));
    //printf("%d ", (int)(cur->prev));
    while(p_el->next != NULL) {
      p_el++;
      printf("%c ", *(hmmp->vertex_labels + p_el->vertex));
    }
  }
}

void dump_modules(struct hmm_s *hmmp)
{
  int i,j;

  printf("\nModule dump:\n");
  for(i = 0 ; i < hmmp->nr_m; i++) {
    printf("module: %s", (*(hmmp->modules + i))->name);
    printf("nr_v: %d\n", (*(hmmp->modules + i))->nr_v);
    printf("vertices: ");
    for(j = 0; j < (*(hmmp->modules + i))->nr_v;j++) {
      printf("%d ",*(((*(hmmp->modules + i))->vertices) + j));
    }
    printf("\n");
  }

  
}

void dump_multi_modules(struct hmm_multi_s *hmmp)
{
  int i,j;

  printf("\nModule dump:\n");
  for(i = 0 ; i < hmmp->nr_m; i++) {
    printf("module: %s", (*(hmmp->modules + i))->name);
    printf("nr_v: %d\n", (*(hmmp->modules + i))->nr_v);
    printf("vertices: ");
    for(j = 0; j < (*(hmmp->modules + i))->nr_v;j++) {
      printf("%d ",*(((*(hmmp->modules + i))->vertices) + j));
    }
    printf("\n");
  }  
}

void dump_distrib_groups(int *distrib_groups, int nr_d)
{
  int i,j;

  printf("\ndistribution groups dump:\n");
  for(i = 0 ; i < nr_d; i++) {
    printf("Group %d: ", i+1);
    while(1) {
      if(*distrib_groups == END) {
	distrib_groups++;
	printf("\n");
	break;
      }
      else {
	printf("%d ", *distrib_groups);
	distrib_groups++;
      }
    }
    printf("\n");
  }
  printf("\n");
}

void dump_trans_tie_groups(struct transition_s *trans_tie_groups, int nr_ttg)
{
  int i,j;

  printf("\ntrans tie groups dump:\n");
  for(i = 0 ; i < nr_ttg; i++) {
    printf("Tie %d: ", i+1);
    while(1) {
      if(trans_tie_groups->from_v == END) {
	trans_tie_groups++;
	printf("\n");
	break;
      }
      else {
	printf("%d->%d ", trans_tie_groups->from_v, trans_tie_groups->to_v);
	trans_tie_groups++;
      }
    }
    printf("\n");
  }
  printf("\n");
}

void dump_prior_struct(struct emission_dirichlet_s *emission_priorsp)
{
  int i,j;
  double *mtx = emission_priorsp->prior_values;
  
  printf("prifile: %s\n", emission_priorsp->name);
  printf("nr components: %d\n", emission_priorsp->nr_components);
  printf("q_values: ");
  for(i = 0; i < emission_priorsp->nr_components; i++) {
    printf("%f ", *((emission_priorsp->q_values) + i));
  }
  printf("\n");
  printf("alpha_sums: ");
  for(i = 0; i < emission_priorsp->nr_components; i++) {
    printf("%f ", *((emission_priorsp->alpha_sums) + i));
  }
  printf("\n");
  printf("alpha_values:\n");
  for(i = 0; i < emission_priorsp->nr_components; i++) {
    for(j = 0; j < emission_priorsp->alphabet_size; j++) {
      printf("%f  ", *mtx);
      mtx++;
    }
    printf("\n");
  }
  printf("\n");
  
  printf("logbeta_values: ");
  for(i = 0; i < emission_priorsp->nr_components; i++) {
    printf("%f ", *((emission_priorsp->logbeta_values) + i));
  }
  printf("\n");
}

void dump_silent_vertices(struct hmm_s *hmmp)
{
  int i;

  printf("silent vertices dump: ");
  for(i = 0;;i++) {
    if(*(hmmp->silent_vertices + i) == END) {
      break;
    }
    else {
      printf("%d ", *(hmmp->silent_vertices + i));
    }
  }
  printf("\n");
}

void dump_silent_vertices_multi(struct hmm_multi_s *hmmp)
{
  int i;

  printf("silent vertices dump: ");
  for(i = 0;;i++) {
    if(*(hmmp->silent_vertices + i) == END) {
      break;
    }
    else {
      printf("%d ", *(hmmp->silent_vertices + i));
    }
  }
  printf("\n");
}

void dump_locked_vertices(struct hmm_s *hmmp)
{
  int i;

  printf("locked vertices dump: ");
  for(i = 0;;i++) {
    if(*(hmmp->locked_vertices + i) == END) {
      break;
    }
    else {
      printf("%d ", *(hmmp->locked_vertices + i));
    }
  }
  printf("\n");
}


void dump_seq(struct letter_s *seq)
{
  int i;
  
  while(seq->letter[0] != '\0') {
    i = 0;
    while(seq->letter[i] != '\0') {
      printf("%c", seq->letter[i]);
      i++;
    }
    printf(";");
    seq++;
  }
  printf("\n");
}

void dump_seqs(struct sequences_s *seq_infop)
{
  int i,j;
  struct letter_s *seqsp;


  printf("seqs dump\n");
  for(i = 0; i < seq_infop->nr_seqs; i++) {
    seqsp = (seq_infop->seqs + i)->seq;
    printf("seq_name = %s\n", (seq_infop->seqs + i)->name);
    printf("seq_length = %d\n", (seq_infop->seqs + i)->length);
    while(seqsp->letter[0] != '\0') {
      j = 0;
      while(seqsp->letter[j] != '\0') {
	printf("%c", seqsp->letter[j]);
	j++;
      }
      printf(";");
      seqsp++;
    }
    printf("\n");
  }
}

void dump_seqs_multi(struct sequences_multi_s *seq_infop)
{
  int i,j,k;
  struct letter_s *seqsp;
  
  printf("multi alphabet seqs dump\n");
  printf("nr_seqs = %d\n", seq_infop->nr_seqs);
  printf("longest_seq = %d\n", seq_infop->longest_seq);
  printf("shortest_seq = %d\n", seq_infop->shortest_seq);
  printf("avg_seq_len = %d\n", seq_infop->avg_seq_len);
  for(i = 0; i < seq_infop->nr_seqs; i++) {
    for(k = 0; k < seq_infop->nr_alphabets; k++) {
      if(k == 0) {
	seqsp = (seq_infop->seqs + i)->seq_1;
      }
      if(k == 1) {
	seqsp = (seq_infop->seqs + i)->seq_2;
      }
      if(k == 2) {
	seqsp = (seq_infop->seqs + i)->seq_3;
      }
      if(k == 3) {
	seqsp = (seq_infop->seqs + i)->seq_4;
      }
      printf("seq_name = %s\n", (seq_infop->seqs + i)->name);
      printf("seq_length = %d\n", (seq_infop->seqs + i)->length);
      while(seqsp->letter[0] != '\0') {
	j = 0;
	while(seqsp->letter[j] != '\0') {
	  printf("%c", seqsp->letter[j]);
	  j++;
	}
	printf(";");
	seqsp++;
      }
      printf("\n");
    }
  }
}

void dump_labeled_seqs_multi(struct sequences_multi_s *seq_infop)
{
  int i,j,k;
  struct letter_s *seqsp;

  printf("seqs dump\n");
  for(i = 0; i < seq_infop->nr_seqs; i++) {
    for(k = 0; k < seq_infop->nr_alphabets; k++) {
      if(k == 0) {
	seqsp = (seq_infop->seqs + i)->seq_1;
      }
      if(k == 1) {
	seqsp = (seq_infop->seqs + i)->seq_2;
      }
      if(k == 2) {
	seqsp = (seq_infop->seqs + i)->seq_3;
      }
      if(k == 3) {
	seqsp = (seq_infop->seqs + i)->seq_4;
      }
      printf("seq_name = %s\n", (seq_infop->seqs + i)->name);
      printf("seq_length = %d\n", (seq_infop->seqs + i)->length);
      while(seqsp->letter[0] != '\0') {
	j = 0;
	while(seqsp->letter[j] != '\0') {
	  printf("%c", seqsp->letter[j]);
	  j++;
	}
	if(k == 0) {
	  printf("/%c",seqsp->label);
	}
	printf(";");
	seqsp++;
      }
      printf("\n");
    }
  }
}

void dump_labeled_seqs(struct sequences_s *seq_infop)
{
  int i,j;
  struct letter_s *seqsp;

  printf("seqs dump\n");
  for(i = 0; i < seq_infop->nr_seqs; i++) {
    seqsp = (seq_infop->seqs + i)->seq;
    printf("seq_name = %s\n", (seq_infop->seqs + i)->name);
    printf("seq_length = %d\n", (seq_infop->seqs + i)->length);
    while(seqsp->letter[0] != '\0') {
      j = 0;
      while(seqsp->letter[j] != '\0') {
	printf("%c", seqsp->letter[j]);
	j++;
      }
      printf("/%c",seqsp->label);
      printf(";");
      seqsp++;
    }
    printf("\n");
  }
}


void dump_msa_seqs(struct msa_sequences_s *msa_seq_infop, int a_size)
{
  int i,j;
  int *cur_pos;

  cur_pos = (int*)(msa_seq_infop->gaps + msa_seq_infop->msa_seq_length);
  printf("msa seqs dump:\n");
  printf("nr_seqs = %d\n", msa_seq_infop->nr_seqs);
  printf("msa_seq_length = %d\n", msa_seq_infop->msa_seq_length);
  printf("nr_lead_columns = %d\n", msa_seq_infop->nr_lead_columns);
  for(i = 0; i < msa_seq_infop->msa_seq_length; i++) {
    printf("pos %d: ", i+1);
    printf("\n");
    printf("label: %c\n", (msa_seq_infop->msa_seq + get_mtx_index(i,0,a_size+1))->label);
    for(j = 0; j < a_size+1; j++) {
      printf("%f ",(msa_seq_infop->msa_seq + get_mtx_index(i,j,a_size+1))->nr_occurences); 
    }
    printf("\n");
    for(j = 0; j < a_size+1; j++) {
      printf("%f ",(msa_seq_infop->msa_seq + get_mtx_index(i,j,a_size+1))->share); 
    }
    printf("\n");
    for(j = 0; j < a_size+1; j++) {
      printf("%f ",(msa_seq_infop->msa_seq + get_mtx_index(i,j,a_size+1))->prior_share); 
    }
    printf("\n");
  }
  printf("\n");
  printf("gaps:\n");
  for(i = 0; i < msa_seq_infop->msa_seq_length; i++) {
    printf("pos %d: ",i+1);
    while(*cur_pos != END) {
      printf("%d ", *cur_pos);
      cur_pos++;
    }
    printf("\n");
    cur_pos++;
  }

  printf("lead columns: ");
  i = 0;
  while(*(msa_seq_infop->lead_columns_start + i) != END) {
    printf("%d ", *(msa_seq_infop->lead_columns_start + i));
    i++;
  }
  printf("\n");

  printf("gap_shares: ");
  for(i = 0; i < msa_seq_infop->msa_seq_length; i++) {
    printf("%f ", *(msa_seq_infop->gap_shares + i));
  }
  printf("\n\n");
}

void dump_msa_seqs_multi(struct msa_sequences_multi_s *msa_seq_infop, struct hmm_multi_s *hmmp)
{
  int i,j;
  int *cur_pos;

  cur_pos = (int*)(msa_seq_infop->gaps + msa_seq_infop->msa_seq_length);
  printf("msa seqs dump:\n");

  printf("nr_seqs = %d\n", msa_seq_infop->nr_seqs);
  printf("msa_seq_length = %d\n", msa_seq_infop->msa_seq_length);
  printf("nr_lead_columns = %d\n\n", msa_seq_infop->nr_lead_columns);
  
  printf("Alphabet 1\n");
  for(i = 0; i < msa_seq_infop->msa_seq_length; i++) {
    printf("pos %d: ", i+1);
    printf("\n");
    printf("label: %c\n", (msa_seq_infop->msa_seq_1 + get_mtx_index(i,0,hmmp->a_size+1))->label);
    for(j = 0; j < hmmp->a_size+1; j++) {
      printf("%f ",(msa_seq_infop->msa_seq_1 + get_mtx_index(i,j,hmmp->a_size+1))->nr_occurences); 
    }
    printf("\n");
    for(j = 0; j < hmmp->a_size+1; j++) {
      printf("%f ",(msa_seq_infop->msa_seq_1 + get_mtx_index(i,j,hmmp->a_size+1))->share); 
    }
    printf("\n");
    for(j = 0; j < hmmp->a_size+1; j++) {
      printf("%f ",(msa_seq_infop->msa_seq_1 + get_mtx_index(i,j,hmmp->a_size+1))->prior_share); 
    }
    printf("\n");
  }
  printf("\n");
  
  if(hmmp->nr_alphabets > 1) {
    printf("Alphabet 2\n");
    for(i = 0; i < msa_seq_infop->msa_seq_length; i++) {
      printf("pos %d: ", i+1);
      printf("\n");
      printf("label: %c\n", (msa_seq_infop->msa_seq_2 + get_mtx_index(i,0,hmmp->a_size_2+1))->label);
      for(j = 0; j < hmmp->a_size_2+1; j++) {
	printf("%f ",(msa_seq_infop->msa_seq_2 + get_mtx_index(i,j,hmmp->a_size_2+1))->nr_occurences); 
      }
      printf("\n");
      for(j = 0; j < hmmp->a_size_2+1; j++) {
	printf("%f ",(msa_seq_infop->msa_seq_2 + get_mtx_index(i,j,hmmp->a_size_2+1))->share); 
      }
      printf("\n");
      for(j = 0; j < hmmp->a_size_2+1; j++) {
	printf("%f ",(msa_seq_infop->msa_seq_2 + get_mtx_index(i,j,hmmp->a_size_2+1))->prior_share); 
      }
      printf("\n");
    }
    printf("\n");
  }

  if(hmmp->nr_alphabets > 2) {
    printf("Alphabet 3\n");
    for(i = 0; i < msa_seq_infop->msa_seq_length; i++) {
      printf("pos %d: ", i+1);
      printf("\n");
      printf("label: %c\n", (msa_seq_infop->msa_seq_3 + get_mtx_index(i,0,hmmp->a_size_3+1))->label);
      for(j = 0; j < hmmp->a_size_3+1; j++) {
	printf("%f ",(msa_seq_infop->msa_seq_3 + get_mtx_index(i,j,hmmp->a_size_3+1))->nr_occurences); 
      }
      printf("\n");
      for(j = 0; j < hmmp->a_size_3+1; j++) {
	printf("%f ",(msa_seq_infop->msa_seq_3 + get_mtx_index(i,j,hmmp->a_size_3+1))->share); 
      }
      printf("\n");
      for(j = 0; j < hmmp->a_size_3+1; j++) {
	printf("%f ",(msa_seq_infop->msa_seq_3 + get_mtx_index(i,j,hmmp->a_size_3+1))->prior_share); 
      }
      printf("\n");
    }
    printf("\n");
  }

  if(hmmp->nr_alphabets > 3) {
    printf("Alphabet 4\n");
    for(i = 0; i < msa_seq_infop->msa_seq_length; i++) {
      printf("pos %d: ", i+1);
      printf("\n");
      printf("label: %c\n", (msa_seq_infop->msa_seq_4 + get_mtx_index(i,0,hmmp->a_size_4+1))->label);
      for(j = 0; j < hmmp->a_size_4+1; j++) {
	printf("%f ",(msa_seq_infop->msa_seq_4 + get_mtx_index(i,j,hmmp->a_size_4+1))->nr_occurences); 
      }
      printf("\n");
      for(j = 0; j < hmmp->a_size_4+1; j++) {
	printf("%f ",(msa_seq_infop->msa_seq_4 + get_mtx_index(i,j,hmmp->a_size_4+1))->share); 
      }
      printf("\n");
      for(j = 0; j < hmmp->a_size_4+1; j++) {
	printf("%f ",(msa_seq_infop->msa_seq_4 + get_mtx_index(i,j,hmmp->a_size_4+1))->prior_share); 
      }
      printf("\n");
    }
    printf("\n");
  }

  
  printf("gaps:\n");
  
  for(i = 0; i < msa_seq_infop->msa_seq_length; i++) {
    printf("pos %d: ",i+1);
    while(*cur_pos != END) {
      printf("%d ", *cur_pos);
      cur_pos++;
    }
    printf("\n");
    cur_pos++;
  }

  
  printf("lead columns: ");
  i = 0;
  while(*(msa_seq_infop->lead_columns_start + i) != END) {
    printf("%d ", *(msa_seq_infop->lead_columns_start + i));
    i++;
  }
  printf("\n");

  printf("gap_shares: ");
  for(i = 0; i < msa_seq_infop->msa_seq_length; i++) {
    printf("%f ", *(msa_seq_infop->gap_shares + i));
  }
  printf("\n\n");
}


void dump_replacement_letters(struct replacement_letter_s *replacement_letters, int a_size)
{
  int i,j;
  int nr_rows, nr_cols;
  double *prob_mtx;

  prob_mtx = replacement_letters->probs;
  nr_rows = replacement_letters->nr_rl;
  nr_cols = a_size;
  
  
  printf("replacement letter dump:\n");
  printf("nr_letters = %d\n", replacement_letters->nr_rl);
  printf("letters: ");
  for(i = 0; i < nr_rows;i++) {
    printf("%c ", (*(replacement_letters->letters + i)));
  }
  printf("\nprobs:\n");
  for(i = 0; i < nr_rows; i++) {
    for(j = 0; j < nr_cols; j++) {
      printf("%f  ", *prob_mtx);
      prob_mtx++;;
    }
    printf("\n");
  }
  printf("\n");
}

void dump_replacement_letters_multi(struct replacement_letter_multi_s *replacement_letters, int alphabet, int a_size)
{
  int i,j;
  int nr_rows, nr_cols;
  double *prob_mtx;

  if(alphabet == 1) {
    prob_mtx = replacement_letters->probs_1;
    nr_rows = replacement_letters->nr_rl_1;
  }
  if(alphabet == 2) {
    prob_mtx = replacement_letters->probs_2;
    nr_rows = replacement_letters->nr_rl_2;
  }
  if(alphabet == 3) {
    prob_mtx = replacement_letters->probs_3;
    nr_rows = replacement_letters->nr_rl_3;
  }
  if(alphabet == 4) {
    prob_mtx = replacement_letters->probs_4;
    nr_rows = replacement_letters->nr_rl_4;
  }
  nr_cols = a_size;
  
  
  printf("replacement letter dump:\n");
  if(alphabet == 1) {
    printf("nr_letters = %d\n", replacement_letters->nr_rl_1);
  }
  if(alphabet == 2) {
    printf("nr_letters = %d\n", replacement_letters->nr_rl_2);
  }
  if(alphabet == 3) {
    printf("nr_letters = %d\n", replacement_letters->nr_rl_3);
  }
  if(alphabet == 4) {
    printf("nr_letters = %d\n", replacement_letters->nr_rl_4);
  }
  printf("letters: ");
  for(i = 0; i < nr_rows;i++) {
    if(alphabet == 1) {
      printf("%c ", (*(replacement_letters->letters_1 + i)));
    }
    if(alphabet == 2) {
      printf("%c ", (*(replacement_letters->letters_2 + i)));
    }
    if(alphabet == 3) {
      printf("%c ", (*(replacement_letters->letters_3 + i)));
    }
    if(alphabet == 4) {
      printf("%c ", (*(replacement_letters->letters_4 + i)));
    }
  }
  printf("\nprobs:\n");
  for(i = 0; i < nr_rows; i++) {
    for(j = 0; j < nr_cols; j++) {
      printf("%f  ", *prob_mtx);
      prob_mtx++;;
    }
    printf("\n");
  }
  printf("\n");
}

void dump_calibration_matrix(int nr_rows, int nr_cols, int *mtx)
{
  int i,j;
  printf("calibration matrix dump:\n");
  printf("nr rows: %d\n", nr_rows);
  printf("nr columns: %d\n", nr_cols);
  for(i = 0; i < nr_rows; i++) {
    for(j = 0; j < nr_cols; j++) {
      printf("%d  ", *mtx);
      mtx++;
    }
    printf("\n");
  }
  printf("\n");
}

void dump_to_silent_trans_array(int nr_v, int **array)
{
  int i;
  int pos;

  printf("To_silent_trans_array_dump:\n");
  for(i = 0; i < nr_v; i++) {
    printf("Transitions from vertex %d: ", i);
    pos = 0;
    while(*(*(array + i) + pos) != END) {
      printf("%d ", *(*(array + i) + pos));
      pos++;
    }
    printf("\n");
  }
  printf("\n");
}

void dump_aa_distrib_mtx(struct aa_distrib_mtx_s *aa_distrib_mtxp)
{
  int i;

  printf("aa_distrib_mtx_dump\n");
  printf("a_size = %d\n", aa_distrib_mtxp->a_size);
  for(i = 0; i < aa_distrib_mtxp->a_size;i++) {
    printf("%f %f %f\n", *(aa_distrib_mtxp->inside_values + i), *(aa_distrib_mtxp->outside_values + i),
	   *(aa_distrib_mtxp->membrane_values + i));
  }
}


void dump_v_list(int *sorted_v_list)
{
  int i = 0;
  printf("v-list dump:\n");
  while(*(sorted_v_list + i) != -99) {
    if(*(sorted_v_list + i) == -9) {
      printf(" | ");
    }
    else {
      printf("%d ", *(sorted_v_list + i));
    }
    i++;
  }
  printf("\n");
}

void dump_labeling(char *labeling, int c)
{
  int i;
  
  printf("labeling dump:\n");
  if(labeling == NULL) {
    printf("NULL\n");
    return;
  }
  for(i = 1; i <= c; i++) {
    printf("%c",labeling[i]);
  }
  printf("\n");
}

void dump_label_tmp_list(int *list)
{
  int i;

  printf("label_list_dump: \n");
  
  for(i = 0; list[i] != TOT_END;i++) {
    printf("%d ",list[i]);
  }
  printf("\n");
}

void dump_labels(char *labels, int nr_labels)
{
  int i;

  printf("label_dump: \n");
  for(i = 0; i < nr_labels; i++) {
    printf("%c ", labels[i]);
  }
  printf("\n");
}

void check_for_corrupt_values(int nr_rows, int nr_cols, double *mtx, char *name)
{
  int v,w;
  int x;

  /* test if this matrix contains strange values */
  for(v = 0; v < nr_rows; v++) {
    for(w = 0; w < nr_cols; w++) {
      x = get_mtx_index(v,w,nr_cols);
      if(*(mtx + x) < 0.0) {
	printf("%s[%d][%d] = %f\n", name, v,w,*(mtx + x));
	if(*(mtx + x) == 0.0) {
	}
      }
      else if(!(*(mtx + x) <= 1000000000000.0)) {
	printf("%s[%d][%d] = %f\n", name, v,w,*(mtx + x));
      }
      else {
	//printf("%s[%d][%d] = %f\n", name, v,w,*(mtx + x));
      }
    }
  } 
}

void dump_weights(double *weights, int nr_seqs)
{
  int i;

  printf("weight dump:\n");
  for(i = 0; i < nr_seqs; i++) {
    printf("%f\n", *(weights + i));
  }
}
