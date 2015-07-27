#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
//#include <double.h>
#include <limits.h>

#include "structs.h"
#include "funcs.h"

//#define DEBUG_ALPHA
//#define DEBUG_REVERSE_SEQ
//#define DEBUG_PRI

#define POS 0


int verbose = NO;



void *malloc_or_die(int mem_size)
{
  void *mem;
  int i;
  
  if(mem_size > 6000 * 6000 * 8) {
    printf("Trying to allocate to much memory: %d bytes\n", mem_size);
    exit(0);
  }

  if((mem = malloc(mem_size)) != NULL) {
    memset(mem, 0, mem_size);
    return mem;
  }
  else {
    perror("Memory trouble");
    exit(-1);
  }
}

void init_float_mtx(double *mtx, double init_nr, int mtx_size)
{
  int i;
  for(i = 0; i < mtx_size; i++) {
    *mtx = init_nr;
    mtx++;
  }
}

void init_viterbi_s_mtx(struct viterbi_s *mtx, double init_nr, int mtx_size)
{
  int i;
  for(i = 0; i < mtx_size; i++) {
    mtx->prob = init_nr;
    mtx++;
  }
}

int get_mtx_index(int row, int col, int row_size)
{
  return row * row_size + col;
}

int get_alphabet_index(struct letter_s *c, char *alphabet, int a_size)
{
  int a_index;
  int found_index;
  int i;

  for(a_index = 0; a_index < a_size; a_index++) {
    
    found_index = NO;

#ifdef DEBUG_ALPHA
    printf("a_index = %d\n", a_index);
    printf("c = %s\n", c->letter);
    printf("Alphabet = %s\n", alphabet);   
#endif

    for(i = 0; c->letter[i] != '\0'; i++) {
      if(*alphabet == ';') {
	found_index = NO;
	break;
      }
      else if(c->letter[i] == *alphabet) {
	found_index = YES;
      }
      else {
	found_index = NO;
	break;
      }
      alphabet++;
    }
    if(found_index == YES) {
      break;
    }
    while(*alphabet != ';') {
      alphabet++;
    }
    alphabet++;
  }
  
  if(a_index >= a_size) {
    return -1;
  }
  return a_index;
}

int get_alphabet_index_msa_query(char *c, char *alphabet, int a_size)
{
  int a_index;
  int found_index;
  int i;

  for(a_index = 0; a_index < a_size; a_index++) {
    
    found_index = NO;

#ifdef DEBUG_ALPHA
    printf("a_index = %d\n", a_index);
    printf("c = %s\n", c);
    printf("Alphabet = %s\n", alphabet);   
#endif

    for(i = 0; c[i] != '\0'; i++) {
      if(*alphabet == ';') {
	found_index = NO;
	break;
      }
      else if(c[i] == *alphabet) {
	found_index = YES;
      }
      else {
	found_index = NO;
	break;
      }
      alphabet++;
    }
    if(found_index == YES) {
      break;
    }
    while(*alphabet != ';') {
      alphabet++;
    }
    alphabet++;
  }
  
  if(a_index >= a_size) {
    return -1;
  }
  return a_index;
}

/* if letter is known to be only one character */
int get_alphabet_index_single(char *alphabet, char letter, int a_size)
{
  int a_index;
  //printf("alphabet = %s\n", alphabet);
  //printf("letter = %c\n", letter);
  for(a_index = 0; a_index < a_size; a_index++) {
    if(letter == *alphabet) {
      return a_index;
    }
    alphabet++;
    alphabet++;
  }
  
  return -1;

}

int get_seq_length(struct letter_s *s)
{
  int i, len;
  
  len = 0;
  while((s->letter)[0] != '\0') {
    len++;
    s++;
  }
  return len;
}

/* checks if there is a path in the hmm from vertex 'v' to vertex 'w'
 * either directly or via silent states */
int path_length(int v, int w, struct hmm_s *hmmp, int length)
{
  int *xp;
  int tot_length;
  int temp_length;
  
  if(length > MAX_GAP_SIZE) {
    return 0;
  }
  
  tot_length = 0;
  if(*(hmmp->transitions + get_mtx_index(v, w, hmmp->nr_v)) != 0.0) {
    tot_length = length + 1; /* direct path */
  }
  for(xp = hmmp->silent_vertices; *xp != END; xp++) {
    if(*(hmmp->transitions + get_mtx_index(v, *xp, hmmp->nr_v)) != 0.0) {
      if((temp_length = path_length(*xp, w, hmmp, length + 1)) > 0) {
	tot_length += length + temp_length; /* path via a silent state */
      }
    }
  }
  return tot_length; /* if 0, then there is no path */
}

/* checks if there is a path in the hmm from vertex 'v' to vertex 'w'
 * either directly or via silent states */
int path_length_multi(int v, int w, struct hmm_multi_s *hmmp, int length)
{
  int *xp;
  int tot_length;
  int temp_length;
  

  if(length > MAX_GAP_SIZE) {
    return 0;
  }
  
  tot_length = 0;
  if(*(hmmp->transitions + get_mtx_index(v, w, hmmp->nr_v)) != 0.0) {
    tot_length = length + 1; /* direct path */
  }
  for(xp = hmmp->silent_vertices; *xp != END; xp++) {
    if(*(hmmp->transitions + get_mtx_index(v, *xp, hmmp->nr_v)) != 0.0) {
      if((temp_length = path_length_multi(*xp, w, hmmp, length + 1)) > 0) {
	tot_length += length + temp_length; /* path via a silent state */
      }
    }
  }
  return tot_length; /* if 0, then there is no path */
}


struct path_element* get_end_path_start(int l, struct hmm_s *hmmp)
{
  struct path_element *pep, *pep_const;
  
  pep_const = *(hmmp->to_trans_array + l);
  pep = pep_const;
  while(pep->vertex != END) {
    while(pep->next != NULL) {
      pep++;
    }
    if(pep->vertex == hmmp->nr_v-1) {
      return pep_const;
    }
    else {
      pep++;
      pep_const = pep;
    }
  }
  return NULL;
}

struct path_element* get_end_path_start_multi(int l, struct hmm_multi_s *hmmp)
{
  struct path_element *pep, *pep_const;
  
  pep_const = *(hmmp->to_trans_array + l);
  pep = pep_const;
  while(pep->vertex != END) {
    while(pep->next != NULL) {
      pep++;
    }
    if(pep->vertex == hmmp->nr_v-1) {
      return pep_const;
    }
    else {
      pep++;
      pep_const = pep;
    }
  }
  return NULL;
}

void print_seq(struct letter_s *seq, FILE *outfile, int seq_nr, char *name, int seq_length)
{
  int i,j;

  fprintf(outfile, "NR %d: %s\n", seq_nr, name);
  fprintf(outfile, ">");
  for(i = 0; i < seq_length; i++) {
    j = 0;
    while(*((seq+i)->letter + j) != '\0') {
      fputc((int)*((seq+i)->letter + j), outfile);
      fputc((int)';', outfile);
      j++;
    }
  }
  fputc('\n', outfile);
}


char* get_profile_vertex_type(int v, int *silent_vertices)
{
  while(*silent_vertices != END) {
    if(v == *silent_vertices) {
      return "silent\n";
    }
    silent_vertices++;
  }

  return "standard\n";
}

void get_aa_distrib_mtx(FILE *distribmtxfile, struct aa_distrib_mtx_s *aa_distrib_mtxp)
{
  int MAX_LINE = 500;
  char s[500];
  char *temp;
  int i;

  if(distribmtxfile == NULL) {
    aa_distrib_mtxp->a_size = -1;
    return;
  }

  i = 0;
  while(fgets(s, MAX_LINE, distribmtxfile) != NULL) {
    if(s[0] == '\n' || s[0] == '#') {
    }
    else {
      aa_distrib_mtxp->a_size = atoi(s);
      aa_distrib_mtxp->inside_values = (double*)malloc_or_die(aa_distrib_mtxp->a_size * sizeof(double));
      aa_distrib_mtxp->outside_values = (double*)malloc_or_die(aa_distrib_mtxp->a_size * sizeof(double));
      aa_distrib_mtxp->membrane_values = (double*)malloc_or_die(aa_distrib_mtxp->a_size * sizeof(double));
      break;
    }
  }
  while(fgets(s, MAX_LINE, distribmtxfile) != NULL){
    if(s[0] == '\n' || s[0] == '#') {
      
    }
    else {
      temp = s+2;
      *(aa_distrib_mtxp->inside_values + i) = strtod(temp, &temp);
      *(aa_distrib_mtxp->outside_values + i) = strtod(temp, &temp);
      *(aa_distrib_mtxp->membrane_values + i) = strtod(temp, &temp);
      *(aa_distrib_mtxp->membrane_values + i) += strtod(temp, &temp);
      i++;
    }
  }
  //dump_aa_distrib_mtx(aa_distrib_mtxp);
}


void get_replacement_letters(FILE *replfile, struct replacement_letter_s *replacement_lettersp)
{
  int MAX_LINE = 1000;
  char s[1000];
  int i,j,k;
  int nr_repl_letters;
  int cur_letter;
  struct letter_s le;
  double prob;
  int a_index;
  int done;
  int a_size;
  char *alphabet;
 

  if(replfile == NULL) {
    replacement_lettersp->nr_rl = 0;
    return;
  }
  
  while(fgets(s, MAX_LINE, replfile) != NULL){
    if(s[0] == '\n' || s[0] == '#') {
    }
    else {
      a_size = atoi(s);
      alphabet = (char*)malloc_or_die(a_size * sizeof(char) * 5);
      break;
    }
  }
  i = 0;
  while(fgets(s, MAX_LINE, replfile) != NULL){
    if(s[0] == '\n' || s[0] == '#') {
      
    }
    else {
      strcpy(alphabet, s);
      break;
    }
  }

  
  while(fgets(s, MAX_LINE, replfile) != NULL){
    if(s[0] == '\n' || s[0] == '#') {
      
    }
    else {
      nr_repl_letters = atoi(s);
      replacement_lettersp->nr_rl = nr_repl_letters;
      replacement_lettersp->letters = (struct letter_s*)(malloc_or_die(nr_repl_letters * sizeof(struct letter_s)));
      replacement_lettersp->probs = (double*)(malloc_or_die(nr_repl_letters * a_size * sizeof(double)));
      break;
    }
  }

  cur_letter = 0;
  while (fgets(s, MAX_LINE, replfile) != NULL) {
    if(s[0] == '#' || s[0] == '\n') {
      
    }
    else {
      i = 0;
      while(s[i] != ' ') {
	*((replacement_lettersp->letters + cur_letter)->letter + i) = s[i];
	i++;
      }
      *((replacement_lettersp->letters + cur_letter)->letter + i) = '\0';
      while(s[i] == ' ' || s[i] == '=') {
	i++;
      }
      done = NO;
      while(done == NO) {
	j = 0;
	while(s[i] != ':') {
	  *(le.letter + j) = s[i];
	  i++;
	  j++;
	}
	*(le.letter + j) = '\0';
	i++;
	
	prob = atof(&s[i]);
	a_index = get_alphabet_index(&le, alphabet, a_size);
	*(replacement_lettersp->probs + get_mtx_index(cur_letter, a_index, a_size)) = prob;
	while(s[i] != ' ' && s[i] != '\n') {
	  i++;
	}
	if(s[i] == '\n') {
	  done = YES;
	}
	i++;
      }
      cur_letter++;
    }
  }
  free(alphabet);
  //dump_replacement_letters(replacement_lettersp, a_size);
}

void get_replacement_letters_multi(FILE *replfile, struct replacement_letter_multi_s *replacement_lettersp)
{
  int MAX_LINE = 1000;
  char s[1000];
  int i,j,k,l,m;
  int nr_repl_letters;
  int cur_letter;
  struct letter_s le;
  double prob;
  int a_index;
  int done;
  int a_size;
  char *alphabet;
  int nr_alphabets;
  int letterindex;
  

  if(replfile == NULL) {
    replacement_lettersp->nr_alphabets = 0;
    return;
  }    
  replacement_lettersp->nr_rl_1 = 0;
  replacement_lettersp->nr_rl_2 = 0;
  replacement_lettersp->nr_rl_3 = 0;
  replacement_lettersp->nr_rl_4 = 0;
 
  while(fgets(s, MAX_LINE, replfile) != NULL){
    if(s[0] == '\n' || s[0] == '#') {
    }
    else {
      nr_alphabets = atoi(s);
      replacement_lettersp->nr_alphabets = nr_alphabets;
      break;
    }
  }  

  for(l = 0; l < nr_alphabets; l++) {
    while(fgets(s, MAX_LINE, replfile) != NULL){
      if(s[0] == '\n' || s[0] == '#') {
      }
      else {
	a_size = atoi(s);
	alphabet = (char*)malloc_or_die(a_size * sizeof(char) * 5);
	break;
      }
    }
    

    i = 0;
    while(fgets(s, MAX_LINE, replfile) != NULL){
      if(s[0] == '\n' || s[0] == '#') {
	
      }
      else {
	strcpy(alphabet, s);
	break;
      }
    }
    
    while(fgets(s, MAX_LINE, replfile) != NULL){
      if(s[0] == '\n' || s[0] == '#') {
      }
      else {
	if(l == 0) {
	  nr_repl_letters = atoi(s);
	  replacement_lettersp->nr_rl_1 = nr_repl_letters;
	  replacement_lettersp->letters_1 = (struct letter_s*)(malloc_or_die(nr_repl_letters * sizeof(struct letter_s)));
	  replacement_lettersp->probs_1 = (double*)(malloc_or_die(nr_repl_letters * a_size * sizeof(double)));
	  break;
	}
	if(l == 1) {
	  nr_repl_letters = atoi(s);
	  replacement_lettersp->nr_rl_2 = nr_repl_letters;
	  replacement_lettersp->letters_2 = (struct letter_s*)(malloc_or_die(nr_repl_letters * sizeof(struct letter_s)));
	  replacement_lettersp->probs_2 = (double*)(malloc_or_die(nr_repl_letters * a_size * sizeof(double)));
	  break;
	}
	if(l == 2) {
	  nr_repl_letters = atoi(s);
	  replacement_lettersp->nr_rl_3 = nr_repl_letters;
	  replacement_lettersp->letters_3 = (struct letter_s*)(malloc_or_die(nr_repl_letters * sizeof(struct letter_s)));
	  replacement_lettersp->probs_3 = (double*)(malloc_or_die(nr_repl_letters * a_size * sizeof(double)));
	  break;
	}
	if(l == 3) {
	  nr_repl_letters = atoi(s);
	  replacement_lettersp->nr_rl_4 = nr_repl_letters;
	  replacement_lettersp->letters_4 = (struct letter_s*)(malloc_or_die(nr_repl_letters * sizeof(struct letter_s)));
	  replacement_lettersp->probs_4 = (double*)(malloc_or_die(nr_repl_letters * a_size * sizeof(double)));
	  break;
	}
      }
    }

    cur_letter = 0;
    letterindex = 0;
    while (fgets(s, MAX_LINE, replfile) != NULL) {
      if(letterindex >= nr_repl_letters) {
	break;
      }
      if(s[0] == '#' || s[0] == '\n') {
	
      }
      else {
	i = 0;
	while(s[i] != ' ') {
	  if(l == 0) {
	    *((replacement_lettersp->letters_1 + cur_letter)->letter + i) = s[i]; 
	  }
	  if(l == 1) {
	    *((replacement_lettersp->letters_2 + cur_letter)->letter + i) = s[i];
	  }
	  if(l == 2) {
	    *((replacement_lettersp->letters_3 + cur_letter)->letter + i) = s[i];
	    }
	  if(l == 3) {
	    *((replacement_lettersp->letters_4 + cur_letter)->letter + i) = s[i];
	  }
	  i++;
	}
	
	if(l == 0) {
	  *((replacement_lettersp->letters_1 + cur_letter)->letter + i) = '\0';
	}
	if(l == 1) {
	  *((replacement_lettersp->letters_2 + cur_letter)->letter + i) = '\0';
	}
	if(l == 2) {
	  *((replacement_lettersp->letters_3 + cur_letter)->letter + i) = '\0';
	}
	if(l == 3) {
	  *((replacement_lettersp->letters_4 + cur_letter)->letter + i) = '\0';
	}
	while(s[i] == ' ' || s[i] == '=') {
	  i++;
	}
	done = NO;
	while(done == NO) {
	  j = 0;
	  while(s[i] != ':') {
	    *(le.letter + j) = s[i];
	    i++;
	    j++;
	  }
	  *(le.letter + j) = '\0';
	  i++;
	  
	  prob = atof(&s[i]);
	  /* get alphabet index for correct alphabet (must hardcode this to the same order as in the hmm */
	  a_index = get_alphabet_index(&le, alphabet, a_size);

	  if(l == 0) {
	    *(replacement_lettersp->probs_1 + get_mtx_index(cur_letter, a_index, a_size)) = prob;
	  }
	  if(l == 1) {
	    *(replacement_lettersp->probs_2 + get_mtx_index(cur_letter, a_index, a_size)) = prob;
	  }
	  if(l == 2) {
	    *(replacement_lettersp->probs_3 + get_mtx_index(cur_letter, a_index, a_size)) = prob;
	  }
	  if(l == 3) {
	    *(replacement_lettersp->probs_4 + get_mtx_index(cur_letter, a_index, a_size)) = prob;
	  }
	  while(s[i] != ' ' && s[i] != '\n') {
	    i++;
	  }
	  if(s[i] == '\n') {
	    done = YES;
	  }
	  i++;
	}
	cur_letter++;
	letterindex++;
      }
      
    }
    free(alphabet);
    //dump_replacement_letters_multi(replacement_lettersp, l+1, a_size);
  }
}

int get_replacement_letter_index(struct letter_s *c, struct replacement_letter_s *replacement_letters)
{
  int a_index;
  for(a_index = 0; a_index < replacement_letters->nr_rl; a_index++) {
    if(strcmp((replacement_letters->letters + a_index)->letter, c->letter) == 0) {
      return a_index;
    }  
  }
  
  return -1;
}

int get_replacement_letter_index_multi(struct letter_s *c, struct replacement_letter_multi_s *replacement_letters, int alphabet)
{
  int a_index;
  if(alphabet == 1) {
    for(a_index = 0; a_index < replacement_letters->nr_rl_1; a_index++) {
      if(strcmp((replacement_letters->letters_1 + a_index)->letter, c->letter) == 0) {
	return a_index;
      }  
    }
  }
  if(alphabet == 2) {
    for(a_index = 0; a_index < replacement_letters->nr_rl_2; a_index++) {
      if(strcmp((replacement_letters->letters_2 + a_index)->letter, c->letter) == 0) {
	return a_index;
      }  
    }
  }
  if(alphabet == 3) {
    for(a_index = 0; a_index < replacement_letters->nr_rl_3; a_index++) {
      if(strcmp((replacement_letters->letters_3 + a_index)->letter, c->letter) == 0) {
	return a_index;
      }  
    }
  }
  if(alphabet == 4) {
    for(a_index = 0; a_index < replacement_letters->nr_rl_4; a_index++) {
      if(strcmp((replacement_letters->letters_4 + a_index)->letter, c->letter) == 0) {
	return a_index;
      }  
    }
  }
  
  return -1;
}

int get_replacement_letter_index_single(char *c, struct replacement_letter_s *replacement_letters)
{
  int a_index;
  for(a_index = 0; a_index < replacement_letters->nr_rl; a_index++) {
    if(((replacement_letters->letters + a_index)->letter)[0] == c) {
      return a_index;
    }  
  }
  
  return -1;
}

char* sequence_as_string(struct letter_s *sequence)
{
  struct letter_s *c;
  char *s;
  int s_length;

  s_length = 0;
  c = sequence;
  while(c->letter[0] != '\0') {
    s_length++;
    c++;
  }
  
  s = (char*)(malloc_or_die(s_length * 5 * sizeof(char)));

  c = sequence;
  while(c->letter[0] != '\0') {
    strcat(s,c->letter);
    c++;
  }
  
  printf("%s\n", s);
  free(s);
  return NULL;
}

void get_viterbi_label_path(struct viterbi_s *cur, struct hmm_s *hmmp,
			    struct viterbi_s *viterbi_mtxp, int row, int row_size, char *labels, int *ip)
{
  struct path_element *p_el;
  
  if(cur->prev == 0) {
    p_el = cur->prevp;
    while(p_el->next != NULL) {
      labels[*ip] = *(hmmp->vertex_labels + p_el->vertex);
      *ip = (*ip) + 1;
      p_el++;
    }
  }
  else {
    get_viterbi_label_path(viterbi_mtxp + get_mtx_index(row-1, cur->prev, row_size), hmmp,
		      viterbi_mtxp, row-1, row_size, labels, ip);
    p_el = cur->prevp;
    labels[*ip] = *(hmmp->vertex_labels + (int)(cur->prev));
    *ip = (*ip) + 1;
    while(p_el->next != NULL) {
      p_el++;
      labels[*ip] = *(hmmp->vertex_labels + p_el->vertex);
      *ip = (*ip) + 1;
    }
  }
}

void get_viterbi_label_path_multi(struct viterbi_s *cur, struct hmm_multi_s *hmmp,
			    struct viterbi_s *viterbi_mtxp, int row, int row_size, char *labels, int *ip)
{
  struct path_element *p_el;
  
  if(cur->prev == 0) {
    p_el = cur->prevp;
    while(p_el->next != NULL) {
      labels[*ip] = *(hmmp->vertex_labels + p_el->vertex);
      *ip = (*ip) + 1;
      p_el++;
    }
  }
  else {
    get_viterbi_label_path_multi(viterbi_mtxp + get_mtx_index(row-1, cur->prev, row_size), hmmp,
				 viterbi_mtxp, row-1, row_size, labels, ip);
    p_el = cur->prevp;
    labels[*ip] = *(hmmp->vertex_labels + (int)(cur->prev));
    *ip = (*ip) + 1;
    while(p_el->next != NULL) {
      p_el++;
      labels[*ip] = *(hmmp->vertex_labels + p_el->vertex);
      *ip = (*ip) + 1;
    }
  }
}


void get_viterbi_path(struct viterbi_s *cur, struct hmm_s *hmmp,
		      struct viterbi_s *viterbi_mtxp, int row, int row_size, int *path, int *ip)
{
  struct path_element *p_el;
  
  //printf("cur->prev = %d\n", cur->prev);
  if(cur->prev == 0) {
    p_el = cur->prevp;
    while(p_el->next != NULL) {
      path[*ip] = p_el->vertex;
      *ip = (*ip) + 1;
      p_el++;
    }
  }
  else {
    get_viterbi_path(viterbi_mtxp + get_mtx_index(row-1, cur->prev, row_size), hmmp,
			   viterbi_mtxp, row-1, row_size, path, ip);
    p_el = cur->prevp;
    path[*ip] = p_el->vertex;
    //printf("%d ", path[*ip]);
    *ip = (*ip) + 1;
    while(p_el->next != NULL) {
      p_el++;
      path[*ip] = p_el->vertex;
      *ip = (*ip) + 1;
    }
  }
}

void get_viterbi_path_multi(struct viterbi_s *cur, struct hmm_multi_s *hmmp,
			    struct viterbi_s *viterbi_mtxp, int row, int row_size, int *path, int *ip)
{
  struct path_element *p_el;
  
  //printf("cur->prev = %d\n", cur->prev);
  if(cur->prev == 0) {
    p_el = cur->prevp;
    while(p_el->next != NULL) {
      path[*ip] = p_el->vertex;
      *ip = (*ip) + 1;
      p_el++;
    }
  }
  else {
    get_viterbi_path_multi(viterbi_mtxp + get_mtx_index(row-1, cur->prev, row_size), hmmp,
			   viterbi_mtxp, row-1, row_size, path, ip);
    p_el = cur->prevp;
    path[*ip] = p_el->vertex;
    //printf("%d ", path[*ip]);
    *ip = (*ip) + 1;
    while(p_el->next != NULL) {
      p_el++;
      path[*ip] = p_el->vertex;
      *ip = (*ip) + 1;
    }
  }
}


void loosen_labels(char *labels, char *loose_labels, int label_looseness, int seq_len)
{
  int *locked_labels;
  int reg_len;
  int i,j;
  char cur;
  
  /* initial memory and copy stuff */
  locked_labels = (int*)(malloc_or_die(seq_len * sizeof(int)));
  memcpy(loose_labels, labels, seq_len * sizeof(char));
  

  /* lock middle labels */
  reg_len = 1;
  cur = labels[0];
  locked_labels[0] = 1;
  for(i = 1; i < seq_len; i++) {
    if(labels[i] != cur) {
      cur = labels[i];
      if(reg_len == i) {
	/* first reg shift, do nothing */
      }
      else {
	if(reg_len % 2 == 1) {
	  locked_labels[i - (reg_len/2 + 1)] = 1;
	}
	else {
	  locked_labels[i - reg_len/2] = 1;
	  locked_labels[i - (reg_len/2 + 1)] = 1;
	}
      }
      reg_len = 1;
    }
    else {
      reg_len++;
    }
  }

  /* loosen labels */
  reg_len = 1;
  cur = labels[0];
  for(i = 1; i < seq_len; i++) {
    if(labels[i] != cur) {
      for(j = 0; j < label_looseness; j++) {
	if(locked_labels[i + j] == 0 && locked_labels[i-1-j] == 0) {
	  //printf("i = %d, j = %d\n", i , j);
	  loose_labels[i + j] = '.';
	  loose_labels[i-1-j] = '.';
	}
	else {
	  break;
	}
      }

      cur = labels[i];
      reg_len = 1;
    }
    else {
      reg_len++;
    }

  }
  free(locked_labels);

  //dump_labeling(labels, seq_len);
  //dump_labeling(loose_labels, seq_len);
}

int read_subst_matrix(double **mtxpp, FILE *substmtxfile)
{
  int MAX_LINE = 1000;
  char s[1000];
  int i,j,k;
  int nr_rows;
  int row_le_index, col_le_index;
  struct letter_s row_le, col_le;
  double prob;
  int row_a_index, col_a_index;
  int done;
  int a_size;
  char *alphabet;
  double *mtxp;

  if(substmtxfile == NULL) {
    return NO;
  }
  
  while(fgets(s, MAX_LINE, substmtxfile) != NULL){
    if(s[0] == '\n' || s[0] == '#') {
    }
    else {
      a_size = atoi(s);
      alphabet = (char*)malloc_or_die(a_size * sizeof(char) * 5);
      *mtxpp = (double*)(malloc_or_die(a_size * (a_size + 1) * sizeof(double)));
      mtxp = *mtxpp;
      break;
    }
  }
  i = 0;
  while(fgets(s, MAX_LINE, substmtxfile) != NULL){
    if(s[0] == '\n' || s[0] == '#') {
      
    }
    else {
      strcpy(alphabet, s);
      break;
    }
  }

  while (fgets(s, MAX_LINE, substmtxfile) != NULL) {
    if(s[0] == '#' || s[0] == '\n') {
      
    }
    else {
      i = 0;
      while(s[i] != ' ') {
	*(row_le.letter + i) = s[i];
	i++;
      }
      *(row_le.letter + i) = '\0';
      while(s[i] == ' ' || s[i] == '=') {
	i++;
      }
      done = NO;
      while(done == NO) {
	j = 0;
	while(s[i] != ':') {
	  *(col_le.letter + j) = s[i];
	  i++;
	  j++;
	}
	*(col_le.letter + j) = '\0';
	i++;
	
	prob = atof(&s[i]);
	row_a_index = get_alphabet_index(&row_le, alphabet, a_size);
	if(row_a_index < 0) {
	  row_a_index = a_size;
	}
	col_a_index = get_alphabet_index(&col_le, alphabet, a_size);
	*(mtxp + get_mtx_index(row_a_index, col_a_index, a_size)) = prob;
	while(s[i] != ' ' && s[i] != '\n') {
	  i++;
	}
	if(s[i] == '\n') {
	  done = YES;
	}
	i++;
      }
    }
  }
  free(alphabet);
  //dump_subst_mtx(mtxp, a_size);
  //exit(0);  
  return YES;
}


int read_subst_matrix_multi(double **mtxpp, double **mtxpp_2, double **mtxpp_3, double **mtxpp_4, FILE *substmtxfile)
{
  int MAX_LINE = 1000;
  char s[1000];
  int i,j,k, l;
  int nr_rows;
  int row_le_index, col_le_index;
  struct letter_s row_le, col_le;
  double prob;
  int row_a_index, col_a_index;
  int done;
  int a_size;
  char *alphabet;
  double *mtxp;
  int nr_alphabets;
  
  if(substmtxfile == NULL) {
    return NO;
  }
  
  while(fgets(s, MAX_LINE, substmtxfile) != NULL){
    if(s[0] == '\n' || s[0] == '#') {
    }
    else {
      nr_alphabets = atoi(s);
      break;
    }
  }
  
  for(l = 0; l < nr_alphabets; l++) {
    while(fgets(s, MAX_LINE, substmtxfile) != NULL){
      if(s[0] == '\n' || s[0] == '#') {
      }
      else {
	a_size = atoi(s);
	alphabet = (char*)malloc_or_die(a_size * sizeof(char) * 5);
	if(l == 0) {
	  *mtxpp = (double*)(malloc_or_die(a_size * (a_size + 1) * sizeof(double)));
	  mtxp = *mtxpp;
	}
	if(l == 1) {
	  *mtxpp_2 = (double*)(malloc_or_die(a_size * (a_size + 1) * sizeof(double)));
	  mtxp = *mtxpp_2;
	}
	if(l == 2) {
	  *mtxpp_3 = (double*)(malloc_or_die(a_size * (a_size + 1) * sizeof(double)));
	  mtxp = *mtxpp_3;
	}
	if(l == 3) {
	  *mtxpp_4 = (double*)(malloc_or_die(a_size * (a_size + 1) * sizeof(double)));
	  mtxp = *mtxpp_4;
	}
	break;
      }
    }
    i = 0;
    while(fgets(s, MAX_LINE, substmtxfile) != NULL){
      if(s[0] == '\n' || s[0] == '#') {
	
      }
      else {
	strcpy(alphabet, s);
	break;
      }
    }
    
    while (fgets(s, MAX_LINE, substmtxfile) != NULL) {
      if(s[0] == '#' || s[0] == '\n') {
	
      }
      else if(strncmp(s, "END", 3) == 0) {
	break;
      }
      else {
	i = 0;
	while(s[i] != ' ') {
	  *(row_le.letter + i) = s[i];
	  i++;
	}
	*(row_le.letter + i) = '\0';
	while(s[i] == ' ' || s[i] == '=') {
	  i++;
	}
	done = NO;
	while(done == NO) {
	  j = 0;
	  while(s[i] != ':') {
	    *(col_le.letter + j) = s[i];
	    i++;
	    j++;
	  }
	  *(col_le.letter + j) = '\0';
	  i++;
	  
	  prob = atof(&s[i]);
	  row_a_index = get_alphabet_index(&row_le, alphabet, a_size);
	  if(row_a_index < 0) {
	    row_a_index = a_size;
	  }
	  col_a_index = get_alphabet_index(&col_le, alphabet, a_size);
	  *(mtxp + get_mtx_index(row_a_index, col_a_index, a_size)) = prob;
	  while(s[i] != ' ' && s[i] != '\n') {
	    i++;
	  }
	  if(s[i] == '\n') {
	    done = YES;
	  }
	  i++;
	}
      }
    }
    free(alphabet);
    //dump_subst_mtx(mtxp, a_size);
    //exit(0);  
  }
  return YES;
}


int read_prior_file(struct emission_dirichlet_s *em_di, struct hmm_s *hmmp, FILE *priorfile)
{
  int j,k;
  double q_value, alpha_value, alpha_sum, logbeta;
  char s[2048];
  char ps[2048];
  char *file_name;
  char *pri;
  
  rewind(priorfile);

  /* put default name */
  strcpy(em_di->name, "default");

  /* put nr of components in struct */
  if(fgets(ps, 2048, priorfile) != NULL) {
  }
  else {
    return -1;
  }
  while(*ps == '#' || *ps == '\n') {
    if(fgets(ps, 2048, priorfile) != NULL) {
    }
    else {
      return -1;
    }
  }
  em_di->nr_components = atoi(&ps[0]);
  
  /* allocate memory for arrays and matrix to this prior struct */
  em_di->q_values = malloc_or_die(em_di->nr_components *
				 sizeof(double));
  em_di->alpha_sums = malloc_or_die(em_di->nr_components *
				   sizeof(double));
  em_di->logbeta_values =
    malloc_or_die(em_di->nr_components * sizeof(double));
  em_di->prior_values = malloc_or_die(em_di->nr_components *
				     hmmp->a_size * sizeof(double));
  
  for(j = 0; j < em_di->nr_components; j++) {
    /* put q-value in array  */
    if(fgets(ps, 2048, priorfile) != NULL) {
    }
    else {
      return -1;
    }
    while(*ps == '#' || *ps == '\n') {
      if(fgets(ps, 2048, priorfile) != NULL) {
      }
      else {
	return -1;
      }
    }
    q_value = atof(&ps[0]);
    *(em_di->q_values + j) = q_value; 
#ifdef DEBUG_PRI
    printf("q_value = %f\n", *(em_di->q_values + j));
#endif
    
    /* put alpha-values of this component in matrix */
    alpha_sum = 0.0;
    k = 0;
    if(fgets(ps, 2048, priorfile) != NULL) {
    }
    else {
      return -1;
    }
    while(*ps == '#' || *ps == '\n') {
      if(fgets(ps, 2048, priorfile) != NULL) {
      }
      else {
	return -1;
      }
    }
    pri = &ps[0];
    for(k = 0; k < hmmp->a_size; k++) {
      alpha_value = strtod(pri, &pri);
      alpha_sum += alpha_value;
      *((em_di->prior_values) +
	get_mtx_index(j, k, hmmp->a_size)) = alpha_value;
    }
    
    /* put sum of alphavalues in array */
    *((em_di->alpha_sums) + j) = alpha_sum; 
    
    /* calculate logB(alpha) for this compoment, store in array*/
    logbeta = 0;
    for(k = 0; k < hmmp->a_size; k++) {
      logbeta += lgamma(*(em_di->prior_values +
			  get_mtx_index(j, k, hmmp->a_size)));
      
#ifdef DEBUG_PRI
      printf("prior_value = %f\n", *((em_di->prior_values) +
				     get_mtx_index(j, k, hmmp->a_size)));
      printf("lgamma_value = %f\n", lgamma(*((em_di->prior_values) +
					     get_mtx_index(j, k, hmmp->a_size))));
#endif
    }
    logbeta = logbeta - lgamma(*(em_di->alpha_sums + j));
    *(em_di->logbeta_values + j) = logbeta;
  }
  
#ifdef DEBUG_PRI
  dump_prior_struct(em_di);
  exit(0);
#endif
}

int read_frequencies(FILE *freqfile, double **aa_freqsp)
{
  char ps[2048];
  int cur;
  int a_size;
  double *aa_freqs;
  
  /* read frequencies */
  if(fgets(ps, 2048, freqfile) != NULL) {
  }
  else {
    return -1;
  }
  while(*ps == '#' || *ps == '\n') {
    if(fgets(ps, 2048, freqfile) != NULL) {
    }
    else {
      return -1;
    }
  }
  
  

  a_size = atoi(ps);

  aa_freqs = (double*)malloc_or_die(a_size * sizeof(double));
  *aa_freqsp = aa_freqs;
  if(freqfile == NULL) {
    printf("Could not read prior frequencies\n");
    exit(0);
  }

  /* read frequencies */
  for(cur = 0; cur < a_size; cur++) {
    if(fgets(ps, 2048, freqfile) != NULL) {
    }
    else {
      return -1;
    }
    while(*ps == '#' || *ps == '\n') {
      if(fgets(ps, 2048, freqfile) != NULL) {
      }
      else {
	return -1;
      }
    }
    *(aa_freqs + cur) = atof(ps);
    //printf("aa_freqs[%d] = %f\n", cur, *(aa_freqs + cur));
  }
  
  return 1;
  
}


int read_frequencies_multi(FILE *freqfile, double **aa_freqsp, double **aa_freqsp_2, double **aa_freqsp_3, double **aa_freqsp_4)
{
  char ps[2048];
  int cur;
  int a_size;
  double *aa_freqs, *aa_freqs_2, *aa_freqs_3, *aa_freqs_4, *aa_freqs_temp;
  int nr_alphabets;
  int i;

  /* read frequencies */
  if(fgets(ps, 2048, freqfile) != NULL) {
  }
  else {
    return -1;
  }
  while(*ps == '#' || *ps == '\n') {
    if(fgets(ps, 2048, freqfile) != NULL) {
    }
    else {
      return -1;
    }
  }
  
  nr_alphabets = atoi(ps);

  for(i = 0; i <= nr_alphabets; i++) {
    if(fgets(ps, 2048, freqfile) != NULL) {
    }
    else {
      return -1;
    }
    while(*ps == '#' || *ps == '\n') {
      if(fgets(ps, 2048, freqfile) != NULL) {
      }
      else {
	return -1;
      }
    }
    a_size = atoi(ps);
    if(i == 0) {
      aa_freqs = (double*)malloc_or_die(a_size * sizeof(double));
      *aa_freqsp = aa_freqs;
      aa_freqs_temp = aa_freqs;
    }
    if(i == 1) {
      aa_freqs_2 = (double*)malloc_or_die(a_size * sizeof(double));
      *aa_freqsp_2 = aa_freqs_2;
      aa_freqs_temp == aa_freqs_2;
    }
    if(i == 2) {
      aa_freqs_3 = (double*)malloc_or_die(a_size * sizeof(double));
      *aa_freqsp_3 = aa_freqs_3;
      aa_freqs_temp == aa_freqs_3;
    }
    if(i == 3) {
      aa_freqs_4 = (double*)malloc_or_die(a_size * sizeof(double));
      *aa_freqsp_4 = aa_freqs_4;
      aa_freqs_temp == aa_freqs_4;
    }
    
    
    /* read frequencies */
    for(cur = 0; cur < a_size; cur++) {
      if(fgets(ps, 2048, freqfile) != NULL) {
      }
      else {
	return -1;
      }
      while(*ps == '#' || *ps == '\n') {
	if(fgets(ps, 2048, freqfile) != NULL) {
	}
	else {
	  return -1;
	}
      }
      *(aa_freqs_temp + cur) = atof(ps);
      //printf("aa_freqs_temp[%d] = %f\n", cur, *(aa_freqs + cur));
    }
  }
  
  return 1;
  
}




int read_prior_file_multi(struct emission_dirichlet_s *em_di, struct hmm_multi_s *hmmp, FILE *priorfile, int alphabet)
{
  int j,k;
  double q_value, alpha_value, alpha_sum, logbeta;
  char s[2048];
  char ps[2048];
  char *file_name;
  char *pri;
  int a_size;
  
  if(alphabet == 1) {
    a_size = hmmp->a_size;
  }
  if(alphabet == 2) {
    a_size = hmmp->a_size_2;
    if(hmmp->nr_alphabets < 2) {
      printf("Trying to read priorfile for alphabet 2, but hmm only has one alphabet\n");
      exit(0);
    }
  }
  if(alphabet == 3) {
    a_size = hmmp->a_size_3;
    if(hmmp->nr_alphabets < 3) {
      printf("Trying to read priorfile for alphabet 3, but hmm only has two alphabets\n");
      exit(0);
    }
  }
  if(alphabet == 4) {
    a_size = hmmp->a_size_4;
    if(hmmp->nr_alphabets < 4) {
      printf("Trying to read priorfile for alphabet 4, but hmm only has three alphabets\n");
      exit(0);
    }
  }

  if(priorfile == NULL) {
    em_di->nr_components = 0;
    return 0;
  }
  rewind(priorfile);
  
  /* put default name */
  strcpy(em_di->name, "default");

  /* put nr of components in struct */
  if(fgets(ps, 2048, priorfile) != NULL) {
  }
  else {
    return -1;
  }
  while(*ps == '#') {
    if(fgets(ps, 2048, priorfile) != NULL) {
    }
    else {
      return -1;
    }
  }
  em_di->nr_components = atoi(&ps[0]);
  
  /* allocate memory for arrays and matrix to this prior struct */
  em_di->q_values = malloc_or_die(em_di->nr_components *
				  sizeof(double));
  em_di->alpha_sums = malloc_or_die(em_di->nr_components *
				    sizeof(double));
  em_di->logbeta_values =
    malloc_or_die(em_di->nr_components * sizeof(double));
  em_di->prior_values = malloc_or_die(em_di->nr_components *
				      a_size * sizeof(double));
  
  for(j = 0; j < em_di->nr_components; j++) {
    /* put q-value in array  */
    if(fgets(ps, 2048, priorfile) != NULL) {
    }
    else {
      return -1;
    }
    while(*ps == '#' || *ps == '\n') {
      if(fgets(ps, 2048, priorfile) != NULL) {
      }
      else {
	return -1;
      }
    }
    q_value = atof(&ps[0]);
    *(em_di->q_values + j) = q_value; 
#ifdef DEBUG_PRI
    printf("q_value = %f\n", *(em_di->q_values + j));
#endif
    
    /* put alpha-values of this component in matrix */
    alpha_sum = 0.0;
    k = 0;
    if(fgets(ps, 2048, priorfile) != NULL) {
    }
    else {
      return -1;
    }
    while(*ps == '#' || *ps == '\n') {
      if(fgets(ps, 2048, priorfile) != NULL) {
      }
      else {
	return -1;
      }
    }
    pri = &ps[0];
    for(k = 0; k < a_size; k++) {
      alpha_value = strtod(pri, &pri);
      alpha_sum += alpha_value;
      *((em_di->prior_values) +
	get_mtx_index(j, k, a_size)) = alpha_value;
    }
    
    /* put sum of alphavalues in array */
    *((em_di->alpha_sums) + j) = alpha_sum; 
    
    /* calculate logB(alpha) for this compoment, store in array*/
    logbeta = 0;
    for(k = 0; k < a_size; k++) {
      logbeta += lgamma(*(em_di->prior_values +
			  get_mtx_index(j, k, a_size)));
      
#ifdef DEBUG_PRI
      printf("prior_value = %f\n", *((em_di->prior_values) +
				     get_mtx_index(j, k, a_size)));
      printf("lgamma_value = %f\n", lgamma(*((em_di->prior_values) +
					     get_mtx_index(j, k, a_size))));
#endif
    }
    logbeta = logbeta - lgamma(*(em_di->alpha_sums + j));
    *(em_di->logbeta_values + j) = logbeta;
  }
  
#ifdef DEBUG_PRI
  dump_prior_struct(em_di);
  exit(0);
#endif
 

}

int read_multi_prior_file_multi(struct emission_dirichlet_s *em_di, struct hmm_multi_s *hmmp, FILE *priorfile, int alphabet)
{

  /* returns negative value if an error in the priorfile is detected, 0 if there is no prior information for the alphabet
   * and a postive value if prior components were read successfully */

  int j,k;
  double q_value, alpha_value, alpha_sum, logbeta;
  char s[2048];
  char ps[2048];
  char *file_name;
  char *pri;
  int a_size;
  int cur_alphabet;

  cur_alphabet = -1;

  if(alphabet == 1) {
    a_size = hmmp->a_size;
  }
  if(alphabet == 2) {
    a_size = hmmp->a_size_2;
    if(hmmp->nr_alphabets < 2) {
      printf("Trying to read priorfile for alphabet 2, but hmm only has one alphabet\n");
      exit(0);
    }
  }
  if(alphabet == 3) {
    a_size = hmmp->a_size_3;
    if(hmmp->nr_alphabets < 3) {
      printf("Trying to read priorfile for alphabet 3, but hmm only has two alphabets\n");
      exit(0);
    }
  }
  if(alphabet == 4) {
    a_size = hmmp->a_size_4;
    if(hmmp->nr_alphabets < 4) {
      printf("Trying to read priorfile for alphabet 4, but hmm only has three alphabets\n");
      exit(0);
    }
  }
  
  if(priorfile == NULL) {
    em_di->nr_components = 0;
    return 0;
  }

  rewind(priorfile);

  /* put default name */
  strcpy(em_di->name, "default");

  /* find correct alphabet */
  
  while(fgets(ps, 2048, priorfile) != NULL) {
    if(strncmp(ps, "ALPHABET:", 9) == 0) {
      cur_alphabet = atoi(&ps[9]);
      if(cur_alphabet == alphabet) {
	break;
      }
    }
  }
  if(cur_alphabet != alphabet) {
    return 0;
  }
  while(fgets(ps, 2048, priorfile) != NULL) {
    if(strncmp(ps, "START", 5) == 0) {
      break;
    }
  }

  /* put nr of components in struct */
  if(fgets(ps, 2048, priorfile) != NULL) {
  }
  else {
    return -1;
  }
  while(*ps == '#' || *ps == '\n') {
    if(fgets(ps, 2048, priorfile) != NULL) {
    }
    else {
      return -1;
    }
  }
  em_di->nr_components = atoi(&ps[0]);
  em_di->alphabet_size = a_size;
  
  /* allocate memory for arrays and matrix to this prior struct */
  em_di->q_values = malloc_or_die(em_di->nr_components *
				  sizeof(double));
  em_di->alpha_sums = malloc_or_die(em_di->nr_components *
				    sizeof(double));
  em_di->logbeta_values =
    malloc_or_die(em_di->nr_components * sizeof(double));
  em_di->prior_values = malloc_or_die(em_di->nr_components *
				      a_size * sizeof(double));
  
  for(j = 0; j < em_di->nr_components; j++) {
    /* put q-value in array  */
    if(fgets(ps, 2048, priorfile) != NULL) {
    }
    else {
      return -1;
    }
    while(*ps == '#' || *ps == '\n') {
      if(fgets(ps, 2048, priorfile) != NULL) {
      }
      else {
	return -1;
      }
    }
    q_value = atof(&ps[0]);
    *(em_di->q_values + j) = q_value; 
#ifdef DEBUG_PRI
    printf("q_value = %f\n", *(em_di->q_values + j));
#endif
    
    /* put alpha-values of this component in matrix */
    alpha_sum = 0.0;
    k = 0;
    if(fgets(ps, 2048, priorfile) != NULL) {
    }
    else {
      return -1;
    }
    while(*ps == '#' || *ps == '\n') {
      if(fgets(ps, 2048, priorfile) != NULL) {
      }
      else {
	return -1;
      }
    }
    pri = &ps[0];
    for(k = 0; k < a_size; k++) {
      alpha_value = strtod(pri, &pri);
      alpha_sum += alpha_value;
      *((em_di->prior_values) +
	get_mtx_index(j, k, a_size)) = alpha_value;
    }
    
    /* put sum of alphavalues in array */
    *((em_di->alpha_sums) + j) = alpha_sum; 
    
    /* calculate logB(alpha) for this compoment, store in array*/
    logbeta = 0;
    for(k = 0; k < a_size; k++) {
      logbeta += lgamma(*(em_di->prior_values +
			  get_mtx_index(j, k, a_size)));
      
#ifdef DEBUG_PRI
      printf("prior_value = %f\n", *((em_di->prior_values) +
				     get_mtx_index(j, k, a_size)));
      printf("lgamma_value = %f\n", lgamma(*((em_di->prior_values) +
					     get_mtx_index(j, k, a_size))));
#endif
    }
    logbeta = logbeta - lgamma(*(em_di->alpha_sums + j));
    *(em_di->logbeta_values + j) = logbeta;
  }
  
#ifdef DEBUG_PRI
  dump_prior_struct(em_di);
  exit(0);
#endif
  
  return 1;

}




int locked_state(struct hmm_s *hmmp, int v)
{
  if(*(hmmp->locked_vertices + v) == YES) {
    return YES;
  }
  else {
    return NO;
  }
}

int locked_state_multi(struct hmm_multi_s *hmmp, int v)
{
  if(*(hmmp->locked_vertices + v) == YES) {
    return YES;
  }
  else {
    return NO;
  }
}

int get_best_reliability_score(double reliability_score_1, double reliability_score_2, double reliability_score_3)
{
  int max = 2;
  if(reliability_score_1 > reliability_score_2 && reliability_score_1 > reliability_score_3) {
    max = 1;
  }
  else if(reliability_score_3 > reliability_score_2) {
    max = 3;
  }
  return max;
}





void itoa(char* s, int nr) {
  int dec, sign;
  char temp[30];
  
  strcpy(temp, fcvt(nr, 0, &dec, &sign));
  if(sign == POS) {
    strncpy(s, temp, dec);
    s[dec] = '\0';
  }
  else {
    s[0] = '-';
    strncpy(s+1, temp, dec);
    s[dec+1] = '\0';
  }


} 

void ftoa(char* s, double nr, int prec) {
  int dec, sign;
  char temp[30];
  int i, pos;

  strcpy(temp, fcvt(nr, prec, &dec, &sign));
  if(sign == POS) {
    if(dec <= 0) {
      s[0] = '0';
      s[1] = '.';
      pos = 2;
      for(i = 0; i > dec; i--) {
	s[pos] = '0';
	pos++;
      } 
      strncpy(s+pos, temp, prec);
      s[pos+prec] = '\0';
    }
    else {
      strncpy(s, temp, dec);
      s[dec] = '.';
      strncpy(s+dec+1, temp+dec, prec);
      s[dec+1+prec] = '\0';
    }
  }
  else {
     if(dec <= 0) {
      s[0] = '-';
      s[1] = '0';
      s[2] = '.';
      pos = 3;
      for(i = 0; i > dec; i--) {
	s[pos] = '0';
	pos++;
      }
      strncpy(s+pos, temp, prec);
      s[pos+prec] = '\0';
    }
    else {
      s[0] = '-';
      strncpy(s+1, temp, dec);
      s[dec+1] = '.';
      strncpy(s+dec+2, temp+dec, prec);
      s[dec+2+prec] = '\0';
    }
  }
}


void hmm_garbage_collection(FILE *hmmfile, struct hmm_s *hmmp)
{
  int i;

  free(hmmp->transitions);
  free(hmmp->log_transitions);
  free(hmmp->emissions);
  free(hmmp->log_emissions);
  for(i = 0; i < hmmp->nr_m; i++) {
    free((*(hmmp->modules + i))->vertices);
  }
  free(hmmp->silent_vertices);
  free(hmmp->vertex_labels);
  free(hmmp->labels);
  free(hmmp->vertex_trans_prior_scalers);
  free(hmmp->vertex_emiss_prior_scalers);
  free(hmmp->modules);
  free(hmmp->to_trans_array);
  free(hmmp->from_trans_array);
  free(hmmp->to_silent_trans_array);
  free(hmmp->tot_transitions);
  free(hmmp->max_log_transitions);
  free(hmmp->tot_from_trans_array);
  free(hmmp->tot_to_trans_array);
  free(hmmp->distrib_groups);
  free(hmmp->trans_tie_groups);
  for(i = 0; i < hmmp->nr_ed; i++) {
    free(hmmp->emission_dirichlets->q_values);
    free(hmmp->emission_dirichlets->alpha_sums);
    free(hmmp->emission_dirichlets->logbeta_values);
    free(hmmp->emission_dirichlets->prior_values);
  }
  free(hmmp->emission_dirichlets);
  free(hmmp->ed_ps);
  fclose(hmmfile);
}

void hmm_garbage_collection_multi(FILE *hmmfile, struct hmm_multi_s *hmmp)
{
  int i;

  free(hmmp->transitions);
  free(hmmp->log_transitions);
  free(hmmp->emissions);
  free(hmmp->log_emissions);
  if(hmmp->nr_alphabets > 1) {
    free(hmmp->emissions_2);
    free(hmmp->log_emissions_2);
  }
  if(hmmp->nr_alphabets > 2) {
    free(hmmp->emissions_3);
    free(hmmp->log_emissions_3);
  }
  if(hmmp->nr_alphabets > 3) {
    free(hmmp->emissions_4);
    free(hmmp->log_emissions_4);
  }
  for(i = 0; i < hmmp->nr_m; i++) {
    free((*(hmmp->modules + i))->vertices);
  }
  free(hmmp->silent_vertices);
  free(hmmp->vertex_labels);
  free(hmmp->labels);
  free(hmmp->vertex_trans_prior_scalers);
  free(hmmp->vertex_emiss_prior_scalers);
  if(hmmp->nr_alphabets > 1) {
    free(hmmp->vertex_emiss_prior_scalers_2);
  }
  if(hmmp->nr_alphabets > 2) {
    free(hmmp->vertex_emiss_prior_scalers_3);
  }
  if(hmmp->nr_alphabets > 3) {
    free(hmmp->vertex_emiss_prior_scalers_4);
  }
  free(hmmp->modules);
  free(hmmp->to_trans_array);
  free(hmmp->from_trans_array);
  free(hmmp->to_silent_trans_array);
  free(hmmp->tot_transitions);
  free(hmmp->max_log_transitions);
  free(hmmp->tot_from_trans_array);
  free(hmmp->tot_to_trans_array);
  free(hmmp->distrib_groups);
  free(hmmp->trans_tie_groups);
  for(i = 0; i < hmmp->nr_ed; i++) {
    free(hmmp->emission_dirichlets->q_values);
    free(hmmp->emission_dirichlets->alpha_sums);
    free(hmmp->emission_dirichlets->logbeta_values);
    free(hmmp->emission_dirichlets->prior_values);
  }
  free(hmmp->emission_dirichlets);
  free(hmmp->ed_ps);
  if(hmmp->nr_alphabets > 1) {
    for(i = 0; i < hmmp->nr_ed_2; i++) {
      free(hmmp->emission_dirichlets_2->q_values);
      free(hmmp->emission_dirichlets_2->alpha_sums);
      free(hmmp->emission_dirichlets_2->logbeta_values);
      free(hmmp->emission_dirichlets_2->prior_values);
    }
    free(hmmp->emission_dirichlets_2);
    free(hmmp->ed_ps_2);
  }
  if(hmmp->nr_alphabets > 2) {
    for(i = 0; i < hmmp->nr_ed_3; i++) {
      free(hmmp->emission_dirichlets_3->q_values);
      free(hmmp->emission_dirichlets_3->alpha_sums);
      free(hmmp->emission_dirichlets_3->logbeta_values);
      free(hmmp->emission_dirichlets_3->prior_values);
    }
    free(hmmp->emission_dirichlets_3);
    free(hmmp->ed_ps_3);
  }
  if(hmmp->nr_alphabets > 3) {
    for(i = 0; i < hmmp->nr_ed_4; i++) {
      free(hmmp->emission_dirichlets_4->q_values);
      free(hmmp->emission_dirichlets_4->alpha_sums);
      free(hmmp->emission_dirichlets_4->logbeta_values);
      free(hmmp->emission_dirichlets_4->prior_values);
    }
    free(hmmp->emission_dirichlets_4);
    free(hmmp->ed_ps_4);
  }
  if(hmmfile != NULL) {
    fclose(hmmfile);
  }
}

void hmm_garbage_collection_multi_no_dirichlet(FILE *hmmfile, struct hmm_multi_s *hmmp)
{
  int i;

  //printf("transitions\n");
  free(hmmp->transitions);
  free(hmmp->log_transitions);
  free(hmmp->emissions);
  //printf("log_emissions\n");
  free(hmmp->log_emissions);
  if(hmmp->nr_alphabets > 1) {
    free(hmmp->emissions_2);
    free(hmmp->log_emissions_2);
  }
  if(hmmp->nr_alphabets > 2) {
    free(hmmp->emissions_3);
    free(hmmp->log_emissions_3);
  }
  if(hmmp->nr_alphabets > 3) {
    free(hmmp->emissions_4);
    free(hmmp->log_emissions_4);
  }
  for(i = 0; i < hmmp->nr_m; i++) {
    //free((*(hmmp->modules + i))->vertices);
  }
  //printf("silent_vertices\n");
  free(hmmp->silent_vertices);
  //printf("vertex_labels\n");
  free(hmmp->vertex_labels);
  //printf("labels\n");
  free(hmmp->labels);
  //printf("trans_prior_scalers\n");
  free(hmmp->vertex_trans_prior_scalers);
  //printf("emiss_prior_scalers\n");
  free(hmmp->vertex_emiss_prior_scalers);
  if(hmmp->nr_alphabets > 1) {
    free(hmmp->vertex_emiss_prior_scalers_2);
  }
  if(hmmp->nr_alphabets > 2) {
    free(hmmp->vertex_emiss_prior_scalers_3);
  }
  if(hmmp->nr_alphabets > 3) {
    free(hmmp->vertex_emiss_prior_scalers_4);
  }
  //printf("modules\n");
  free(hmmp->modules);
  free(hmmp->to_trans_array);
  free(hmmp->from_trans_array);
  free(hmmp->to_silent_trans_array);
  //printf("tot_transitionss\n");
  free(hmmp->tot_transitions);
  free(hmmp->max_log_transitions);
  free(hmmp->tot_from_trans_array);
  free(hmmp->tot_to_trans_array);
  free(hmmp->distrib_groups);
  //printf("trans_tie_groups\n");
  free(hmmp->trans_tie_groups);
  
  if(hmmfile != NULL) {
    fclose(hmmfile);
  }
}

void msa_seq_garbage_collection_multi(struct msa_sequences_multi_s *msa_seq_info, int nr_alphabets)
{
  
}

void seq_garbage_collection_multi(struct sequences_multi_s *seq_info, int nr_alphabets)
{
  
}

struct letter_s* get_reverse_seq(struct letter_s *seq_const, int seq_len)
{
  struct letter_s *reverse_seq, *seq;
  int i;

  seq = seq_const;
  reverse_seq = (struct letter_s*)(malloc_or_die((seq_len+1)* sizeof(struct letter_s)));
  for(i = seq_len - 1; i >= 0; i--) {
    memcpy(reverse_seq + i, seq, sizeof(struct letter_s));
    seq++;
  }
  
  memcpy(reverse_seq + seq_len, seq, sizeof(struct letter_s));
  

#ifdef DEBUG_REVERSE_SEQ
  dump_seq(reverse_seq);
  dump_seq(seq_const);
  printf("reverse_seq = %x\n", reverse_seq);
#endif
  
  return reverse_seq;
}

void get_reverse_seq_multi(struct sequence_multi_s *seqs, struct letter_s **reverse_seq_1,
			   struct letter_s **reverse_seq_2, struct letter_s **reverse_seq_3,
			   struct letter_s **reverse_seq_4, struct hmm_multi_s *hmmp, int seq_len)
{
  /* letter_s pointers allocated here must be freed by caller */

  struct letter_s *reverse_seq, *seq, *seq_const;
  int i,j;
  *reverse_seq_1 = NULL;
  *reverse_seq_2 = NULL;
  *reverse_seq_3 = NULL;
  *reverse_seq_4 = NULL;

  for(j = 0; j < hmmp->nr_alphabets; j++) {
    reverse_seq = (struct letter_s*)(malloc_or_die((seq_len+1)* sizeof(struct letter_s)));
    if(j == 0) {
      seq = seqs->seq_1;
      *reverse_seq_1 = reverse_seq;
    }
    if(j == 1) {
      seq = seqs->seq_2;
      *reverse_seq_2 = reverse_seq;
    }
    if(j == 1) {
      seq = seqs->seq_3;
      *reverse_seq_3 = reverse_seq;
    }
    if(j == 1) {
      seq = seqs->seq_4;
      *reverse_seq_4 = reverse_seq;
    }
    seq_const = seq;
    for(i = seq_len - 1; i >= 0; i--) {
      memcpy(reverse_seq + i, seq, sizeof(struct letter_s));
      seq++;
    }
    
    memcpy(reverse_seq + seq_len, seq, sizeof(struct letter_s));
    
#ifdef DEBUG_REVERSE_SEQ
    printf("seq_len = %d\n", seq_len);
    dump_seq(reverse_seq);
    dump_seq(*reverse_seq_1);
    dump_seq(seq_const);
    printf("reverse_seq = %x\n", reverse_seq);
#endif
  }
}


void get_reverse_msa_seq(struct msa_sequences_s *msa_seq_infop, struct msa_sequences_s *reverse_msa_seq_infop, int a_size)
{

  /* note that this function does not implement the gaps-pointing, everything points to END */
  
  int i,j;

  reverse_msa_seq_infop->nr_seqs = msa_seq_infop->nr_seqs;
  reverse_msa_seq_infop->msa_seq_length = msa_seq_infop->msa_seq_length;
  reverse_msa_seq_infop->nr_lead_columns = msa_seq_infop->nr_lead_columns;
  reverse_msa_seq_infop->msa_seq = (struct msa_letter_s*)(malloc_or_die(msa_seq_infop->msa_seq_length * (a_size + 1) *
									sizeof(struct msa_letter_s)));
  reverse_msa_seq_infop->gaps = (int**)(malloc_or_die(msa_seq_infop->msa_seq_length * sizeof(int*) +
						      msa_seq_infop->msa_seq_length * sizeof(int)));
  reverse_msa_seq_infop->lead_columns_start = (int*)(malloc_or_die((msa_seq_infop->nr_lead_columns +1)* sizeof(int)));
  reverse_msa_seq_infop->lead_columns_end = reverse_msa_seq_infop->lead_columns_start + (msa_seq_infop->nr_lead_columns);
  reverse_msa_seq_infop->gap_shares = (double*)(malloc_or_die(msa_seq_infop->msa_seq_length * sizeof(double)));
  
  /* get sequence data */
  j = 0;
  for(i = msa_seq_infop->msa_seq_length - 1; i >= 0; i--) {
    memcpy(reverse_msa_seq_infop->msa_seq + (i * (a_size + 1)),
					     msa_seq_infop->msa_seq + (j * (a_size + 1)),
								       sizeof(struct msa_letter_s) * (a_size+1));
    j++;
  }
  
  /* get gap data (not implemented) */
  for(i = 0; i < msa_seq_infop->msa_seq_length; i++) {
    *(reverse_msa_seq_infop->gaps + (msa_seq_infop->msa_seq_length + i)) = END;
  }
  for(i = 0; i < msa_seq_infop->msa_seq_length; i++) {
    *(reverse_msa_seq_infop->gaps + i) = (int*)(reverse_msa_seq_infop->gaps + (msa_seq_infop->msa_seq_length + i));
  }

  /* get gap shares */
  j = 0;
  for(i = msa_seq_infop->msa_seq_length - 1; i >= 0; i--) {
    memcpy(reverse_msa_seq_infop->gap_shares + i, msa_seq_infop->gap_shares + j, 1 * sizeof(double));
    j++;
  }

  /* get lead_columns */
  j = 0;
  for(i = msa_seq_infop->nr_lead_columns - 1; i >= 0; i--) {
    *(reverse_msa_seq_infop->lead_columns_start + i) = (msa_seq_infop->msa_seq_length - 1) - *(msa_seq_infop->lead_columns_start + j);
    j++;
  }
  *(reverse_msa_seq_infop->lead_columns_start + msa_seq_infop->nr_lead_columns) = END;

#ifdef DEBUG_REVERSE_SEQ
  dump_msa_seqs(msa_seq_infop, a_size);
  dump_msa_seqs(reverse_msa_seq_infop, a_size);
#endif
}

void get_reverse_msa_seq_multi(struct msa_sequences_multi_s *msa_seq_infop, struct msa_sequences_multi_s *reverse_msa_seq_infop,
			       struct hmm_multi_s *hmmp)
{

  /* note that this function does not implement the gaps-pointing, everything points to END */
  
  int i,j;
  int nr_alphabets;
  int a_size, a_size_2, a_size_3, a_size_4;

  nr_alphabets = hmmp->nr_alphabets;
  a_size = hmmp->a_size;
  a_size_2 = hmmp->a_size_2;
  a_size_3 = hmmp->a_size_3;
  a_size_4 = hmmp->a_size_4;
  


  reverse_msa_seq_infop->nr_seqs = msa_seq_infop->nr_seqs;
  reverse_msa_seq_infop->msa_seq_length = msa_seq_infop->msa_seq_length;
  reverse_msa_seq_infop->nr_lead_columns = msa_seq_infop->nr_lead_columns;
  reverse_msa_seq_infop->msa_seq_1 = (struct msa_letter_s*)(malloc_or_die(msa_seq_infop->msa_seq_length * (a_size + 1) *
									sizeof(struct msa_letter_s)));
  if(nr_alphabets > 1) {
    reverse_msa_seq_infop->msa_seq_2 = (struct msa_letter_s*)(malloc_or_die(msa_seq_infop->msa_seq_length * (a_size_2 + 1) *
									    sizeof(struct msa_letter_s)));
  }
  if(nr_alphabets > 2) {
    reverse_msa_seq_infop->msa_seq_3 = (struct msa_letter_s*)(malloc_or_die(msa_seq_infop->msa_seq_length * (a_size_3 + 1) *
									    sizeof(struct msa_letter_s)));
  }
  if(nr_alphabets > 3) {
    reverse_msa_seq_infop->msa_seq_4 = (struct msa_letter_s*)(malloc_or_die(msa_seq_infop->msa_seq_length * (a_size_4 + 1) *
									    sizeof(struct msa_letter_s)));
  }
  reverse_msa_seq_infop->gaps = (int**)(malloc_or_die(msa_seq_infop->msa_seq_length * sizeof(int*) +
						      msa_seq_infop->msa_seq_length * sizeof(int)));
  reverse_msa_seq_infop->lead_columns_start = (int*)(malloc_or_die((msa_seq_infop->nr_lead_columns +1)* sizeof(int)));
  reverse_msa_seq_infop->lead_columns_end = reverse_msa_seq_infop->lead_columns_start + (msa_seq_infop->nr_lead_columns);
  reverse_msa_seq_infop->gap_shares = (double*)(malloc_or_die(msa_seq_infop->msa_seq_length * sizeof(double)));
  
  /* get sequence data */
  j = 0;
  for(i = msa_seq_infop->msa_seq_length - 1; i >= 0; i--) {
    memcpy(reverse_msa_seq_infop->msa_seq_1 + (i * (a_size + 1)),
	   msa_seq_infop->msa_seq_1 + (j * (a_size + 1)),
	   sizeof(struct msa_letter_s) * (a_size+1));
    j++;
  }
  if(nr_alphabets > 1) {
    j = 0;
    for(i = msa_seq_infop->msa_seq_length - 1; i >= 0; i--) {
      memcpy(reverse_msa_seq_infop->msa_seq_2 + (i * (a_size_2 + 1)),
	     msa_seq_infop->msa_seq_2 + (j * (a_size_2 + 1)),
	     sizeof(struct msa_letter_s) * (a_size_2+1));
      j++;
    }
  }
  if(nr_alphabets > 2) {
    j = 0;
    for(i = msa_seq_infop->msa_seq_length - 1; i >= 0; i--) {
      memcpy(reverse_msa_seq_infop->msa_seq_3 + (i * (a_size_3 + 1)),
	     msa_seq_infop->msa_seq_3 + (j * (a_size_3 + 1)),
	     sizeof(struct msa_letter_s) * (a_size_3+1));
      j++;
    }
  }
  if(nr_alphabets > 3) {
    j = 0;
    for(i = msa_seq_infop->msa_seq_length - 1; i >= 0; i--) {
      memcpy(reverse_msa_seq_infop->msa_seq_4 + (i * (a_size_4 + 1)),
	     msa_seq_infop->msa_seq_4 + (j * (a_size_4 + 1)),
	     sizeof(struct msa_letter_s) * (a_size_4+1));
    j++;
    }
  }



  /* get gap data (not implemented) */
  for(i = 0; i < msa_seq_infop->msa_seq_length; i++) {
    *(reverse_msa_seq_infop->gaps + (msa_seq_infop->msa_seq_length + i)) = END;
  }
  for(i = 0; i < msa_seq_infop->msa_seq_length; i++) {
    *(reverse_msa_seq_infop->gaps + i) = (int*)(reverse_msa_seq_infop->gaps + (msa_seq_infop->msa_seq_length + i));
  }

  /* get gap shares */
  j = 0;
  for(i = msa_seq_infop->msa_seq_length - 1; i >= 0; i--) {
    memcpy(reverse_msa_seq_infop->gap_shares + i, msa_seq_infop->gap_shares + j, 1 * sizeof(double));
    j++;
  }

  /* get lead_columns */
  j = 0;
  for(i = msa_seq_infop->nr_lead_columns - 1; i >= 0; i--) {
    *(reverse_msa_seq_infop->lead_columns_start + i) = (msa_seq_infop->msa_seq_length - 1) - *(msa_seq_infop->lead_columns_start + j);
    j++;
  }
  *(reverse_msa_seq_infop->lead_columns_start + msa_seq_infop->nr_lead_columns) = END;

#ifdef DEBUG_REVERSE_SEQ
  dump_msa_seqs(msa_seq_infop, a_size);
  dump_msa_seqs(reverse_msa_seq_infop, a_size);
#endif
}



void get_msa_labels(FILE *labelfile, struct msa_sequences_s *msa_seq_infop, struct hmm_s *hmmp) {
  char row[30000];
  int i,j;

  rewind(labelfile);
  while(1) {
    if(fgets(row,30000,labelfile) != NULL) {
      if(row[0] == '/') {
	for(i = 1; row[i] != '/';i++) {
	  (msa_seq_infop->msa_seq + (*(msa_seq_infop->lead_columns_start + i-1)) * (hmmp->a_size+1))->label = row[i];
	}
      }
    }
    else {
      break;
    }
  }
}

void get_msa_labels_all_columns(FILE *labelfile, struct msa_sequences_s *msa_seq_infop, struct hmm_s *hmmp) {
  char row[30000];
  int i,j;

  rewind(labelfile);
  while(1) {
    if(fgets(row,30000,labelfile) != NULL) {
      if(row[0] == '/') {
	for(i = 1; row[i] != '/';i++) {
	  (msa_seq_infop->msa_seq + (i-1) * (hmmp->a_size+1))->label = row[i];
	}
      }
    }
    else {
      break;
    }
  }
}

int update_shares_prior(struct emission_dirichlet_s *em_di, struct hmm_s *hmmp,
			struct msa_sequences_s *msa_seq_infop, int l)
{
  int nr_components, comps, a_index; 
  double occ_sums;
  double q_value, scaling_factor, X_sum, *X_values, ed_res1, *logbeta_an_values;
  double exponent, prior_prob, tot_prior_prob;
  
 
  /************the update part ********************/
  nr_components = em_di->nr_components;
  logbeta_an_values = malloc_or_die(nr_components * sizeof(double));
  scaling_factor = -FLT_MAX;
  X_sum = 0.0;
  X_values = malloc_or_die(hmmp->a_size * sizeof(double));
 
 
  /* calculate logB(alpha + n) for all components +
   * calculate scaling factor for logB(alpha + n) - logB(alpha) */
  for(comps = 0; comps < nr_components; comps++) {
    ed_res1 = 0;
    occ_sums = 0;
    for(a_index = 0; a_index < hmmp->a_size; a_index++) {
      ed_res1 += lgamma(*(em_di->prior_values +
			  get_mtx_index(comps, a_index, hmmp->a_size)) +
			(double)((msa_seq_infop->msa_seq +
				 get_mtx_index(l,a_index,hmmp->a_size+1))->nr_occurences));
      occ_sums += (msa_seq_infop->msa_seq + get_mtx_index(l,a_index, hmmp->a_size+1))->nr_occurences;
    }
    ed_res1 = ed_res1 - lgamma(*(em_di->alpha_sums + comps) + (double)(occ_sums));
    *(logbeta_an_values + comps) = ed_res1;
    if((ed_res1 = ed_res1 - *(em_di->logbeta_values + comps)) > scaling_factor) {
      scaling_factor = ed_res1;
    }
  }
  
  /* calculate all the Xi's */
  for(a_index = 0; a_index < hmmp->a_size; a_index++) {
    *(X_values + a_index) = 0;
    for(comps = 0; comps < nr_components; comps++) {
      q_value = *(em_di->q_values + comps);
      exponent = (*(logbeta_an_values + comps) - *(em_di->logbeta_values + comps) - 
		  scaling_factor);
      prior_prob = (*(em_di->prior_values + get_mtx_index(comps,a_index, hmmp->a_size)) +
		    (double)((msa_seq_infop->msa_seq +
			     get_mtx_index(l,a_index,hmmp->a_size+1))->nr_occurences));
      tot_prior_prob = (*(em_di->alpha_sums + comps) + (double)(occ_sums));
      *(X_values + a_index) += q_value * exp(exponent) * prior_prob / tot_prior_prob;
#ifdef DEBUG_PRI
      printf("\nscaling factor = %f\n", scaling_factor);
      printf("a_index = %d\n", a_index);
      printf("q_value = %f\n", q_value);
      printf("exponent = %f\n", exponent);
      printf("prior_prob = %f\n", prior_prob);
      printf("tot_prior_prob = %f\n\n", tot_prior_prob);
#endif
    }
    X_sum += *(X_values + a_index);
  }
  
  /* update share values */
  for(a_index = 0; a_index < hmmp->a_size; a_index++) {
    ed_res1 = *(X_values + a_index) / X_sum;
    if(ed_res1 != 0.0) {
      (msa_seq_infop->msa_seq + get_mtx_index(l, a_index, hmmp->a_size+1))->prior_share = ed_res1;
    }
    else {
      (msa_seq_infop->msa_seq + get_mtx_index(l, a_index, hmmp->a_size+1))->prior_share = ed_res1;
    }
  }  

  /* cleanup */
  free(logbeta_an_values);
  free(X_values);
}


int replacement_letter(struct letter_s *cur_letterp, struct replacement_letter_s *replacement_letters, 
		       struct msa_sequences_s *msa_seq_infop, struct hmm_s *hmmp, int seq_pos)
{
  int i,j,k;
  struct letter_s *repl_letter;
  int same_letter;

  /* find out if letter in cur_letterp is a replacement_letter */
  for(i = 0; i < replacement_letters->nr_rl; i++) {
    repl_letter = replacement_letters->letters + i;
    same_letter = YES;
    j = 0;
    while(*(repl_letter->letter + j) != '\0' && *(cur_letterp->letter + j) != '\0') {
      if(*(repl_letter->letter + j) == *(cur_letterp->letter + j)) {
      }
      else {
	same_letter = NO;
      }
      j++;
    }
    if(*(repl_letter->letter + j) != '\0' || *(cur_letterp->letter + j) != '\0') {
      same_letter = NO;
    }
    else if(same_letter == YES) {
	break;
    }
  }
  if(same_letter == NO) {
    return NO;
  }
  else { /* k represents the regular letter, i represents which repl_letter this is */ 
    for(k = 0; k < hmmp->a_size; k++) {
      (msa_seq_infop->msa_seq + get_mtx_index(seq_pos,k, hmmp->a_size+1))->nr_occurences += 
	*(replacement_letters->probs + get_mtx_index(i,k,hmmp->a_size));
    }
    return YES;
  }
  
}


void get_labels_multi(FILE *labelfile, struct sequences_multi_s *seq_infop, struct hmm_multi_s *hmmp, int seq_nr) {
  char row[30000];
  int i,j;

  rewind(labelfile);
  while(1) {
    if(fgets(row,30000,labelfile) != NULL) {
      if(row[0] == '/') {
	for(i = 1; row[i] != '/';i++) {
	  ((seq_infop->seqs + seq_nr)->seq_1 + (i-1))->label = row[i];
	  if(hmmp->nr_alphabets > 1) {
	    ((seq_infop->seqs + seq_nr)->seq_2 + (i-1))->label = row[i];
	  }
	  if(hmmp->nr_alphabets > 2) {
	    ((seq_infop->seqs + seq_nr)->seq_3 + (i-1))->label = row[i];
	  }
	  if(hmmp->nr_alphabets > 3) {
	    ((seq_infop->seqs + seq_nr)->seq_4 + (i-1))->label = row[i];
	  }
	}
      }
    }
    else {
      break;
    }
  }
  rewind(labelfile);
}

void get_msa_labels_multi(FILE *labelfile, struct msa_sequences_multi_s *msa_seq_infop, struct hmm_multi_s *hmmp) {
  char row[30000];
  int i,j;

  rewind(labelfile);
  while(1) {
    if(fgets(row,30000,labelfile) != NULL) {
      if(row[0] == '/') {
	for(i = 1; row[i] != '/';i++) {
	  (msa_seq_infop->msa_seq_1 + (*(msa_seq_infop->lead_columns_start + i-1)) * (hmmp->a_size+1))->label = row[i];
	  if(hmmp->nr_alphabets > 1) {
	    (msa_seq_infop->msa_seq_2 + (*(msa_seq_infop->lead_columns_start + i-1)) * (hmmp->a_size_2+1))->label = row[i];
	  }
	  if(hmmp->nr_alphabets > 2) {
	    (msa_seq_infop->msa_seq_3 + (*(msa_seq_infop->lead_columns_start + i-1)) * (hmmp->a_size_3+1))->label = row[i];
	  }
	  if(hmmp->nr_alphabets > 3) {
	    (msa_seq_infop->msa_seq_4 + (*(msa_seq_infop->lead_columns_start + i-1)) * (hmmp->a_size_4+1))->label = row[i];
	  }
	}
      }
    }
    else {
      break;
    }
  }
}



void get_msa_labels_all_columns_multi(FILE *labelfile, struct msa_sequences_multi_s *msa_seq_infop, struct hmm_multi_s *hmmp) {
  char row[30000];
  int i,j;

  rewind(labelfile);
  while(1) {
    if(fgets(row,30000,labelfile) != NULL) {
      if(row[0] == '/') {
	for(i = 1; row[i] != '/';i++) {
	  (msa_seq_infop->msa_seq_1 + (i-1) * (hmmp->a_size+1))->label = row[i];
	  if(hmmp->nr_alphabets > 1) {
	    (msa_seq_infop->msa_seq_2 + (i-1) * (hmmp->a_size_2+1))->label = row[i];
	  }
	  if(hmmp->nr_alphabets > 2) {
	    (msa_seq_infop->msa_seq_3 + (i-1) * (hmmp->a_size_3+1))->label = row[i];
	  }
	  if(hmmp->nr_alphabets > 3) {
	    (msa_seq_infop->msa_seq_4 + (i-1) * (hmmp->a_size_4+1))->label = row[i];
	  }
	}
      }
    }
    else {
      break;
    }
  }
}


int update_shares_prior_multi(struct emission_dirichlet_s *em_di, struct hmm_multi_s *hmmp,
			struct msa_sequences_multi_s *msa_seq_infop, int l, int alphabet)
{
  int nr_components, comps, a_index; 
  double occ_sums;
  double q_value, scaling_factor, X_sum, *X_values, ed_res1, *logbeta_an_values;
  double exponent, prior_prob, tot_prior_prob;
  int a_size;
  struct msa_letter_s *msa_seq;

  /* check if this alphabet has a prior */
  if(em_di->nr_components <= 0) {
    //printf("em_di nr comps = %d\n",em_di->nr_components);
    //printf("alphabet = %d\n", alphabet);
    return;
  }
  
  /* set a_size and msa_seq according to alphabet */
  if(alphabet == 1) {
    a_size = hmmp->a_size;
    msa_seq = msa_seq_infop->msa_seq_1;
  }
  if(alphabet == 2) {
    a_size = hmmp->a_size_2;
    msa_seq = msa_seq_infop->msa_seq_2;
  }
  if(alphabet == 3) {
    a_size = hmmp->a_size_3;
    msa_seq = msa_seq_infop->msa_seq_3; 
  }
  if(alphabet == 4) {
    a_size = hmmp->a_size_4;
    msa_seq = msa_seq_infop->msa_seq_4;
  }
 
  /************the update part ********************/
  nr_components = em_di->nr_components;
  logbeta_an_values = malloc_or_die(nr_components * sizeof(double));
  scaling_factor = -FLT_MAX;
  X_sum = 0.0;
  X_values = malloc_or_die(a_size * sizeof(double));
 
 
  /* calculate logB(alpha + n) for all components +
   * calculate scaling factor for logB(alpha + n) - logB(alpha) */
  for(comps = 0; comps < nr_components; comps++) {
    ed_res1 = 0;
    occ_sums = 0;
    for(a_index = 0; a_index < a_size; a_index++) {
      ed_res1 += lgamma(*(em_di->prior_values +
			  get_mtx_index(comps, a_index, a_size)) +
			(double)((msa_seq +
				 get_mtx_index(l,a_index,a_size+1))->nr_occurences));
      occ_sums += (msa_seq + get_mtx_index(l,a_index, a_size+1))->nr_occurences;
    }
    ed_res1 = ed_res1 - lgamma(*(em_di->alpha_sums + comps) + (double)(occ_sums));
    *(logbeta_an_values + comps) = ed_res1;
    if((ed_res1 = ed_res1 - *(em_di->logbeta_values + comps)) > scaling_factor) {
      scaling_factor = ed_res1;
    }
  }
  
  /* calculate all the Xi's */
  for(a_index = 0; a_index < a_size; a_index++) {
    *(X_values + a_index) = 0;
    for(comps = 0; comps < nr_components; comps++) {
      q_value = *(em_di->q_values + comps);
      exponent = (*(logbeta_an_values + comps) - *(em_di->logbeta_values + comps) - 
		  scaling_factor);
      prior_prob = (*(em_di->prior_values + get_mtx_index(comps,a_index, a_size)) +
		    (double)((msa_seq +
			     get_mtx_index(l,a_index,a_size+1))->nr_occurences));
      tot_prior_prob = (*(em_di->alpha_sums + comps) + (double)(occ_sums));
      *(X_values + a_index) += q_value * exp(exponent) * prior_prob / tot_prior_prob;
#ifdef DEBUG_PRI
      printf("\nscaling factor = %f\n", scaling_factor);
      printf("a_index = %d\n", a_index);
      printf("q_value = %f\n", q_value);
      printf("exponent = %f\n", exponent);
      printf("prior_prob = %f\n", prior_prob);
      printf("tot_prior_prob = %f\n\n", tot_prior_prob);
#endif
    }
    X_sum += *(X_values + a_index);
  }
  
  /* update share values */
  //printf("a_size = %d\n",a_size);
  for(a_index = 0; a_index < a_size; a_index++) {
    ed_res1 = *(X_values + a_index) / X_sum;
    if(ed_res1 != 0.0) {
      (msa_seq + get_mtx_index(l, a_index, a_size+1))->prior_share = ed_res1;
    }
    else {
      (msa_seq + get_mtx_index(l, a_index, a_size+1))->prior_share = ed_res1;
    }
    //printf("msa_seq = %f\n",(msa_seq + get_mtx_index(l, a_index, a_size+1))->prior_share);
  }  

  //printf("returning\n");
  /* cleanup */
  free(logbeta_an_values);
  free(X_values);
}


int replacement_letter_multi(struct letter_s *cur_letterp, struct replacement_letter_multi_s *replacement_letters, 
		       struct msa_sequences_multi_s *msa_seq_infop, struct hmm_multi_s *hmmp, int seq_pos, int alphabet)
{
  int i,j,k;
  struct letter_s *repl_letter;
  int same_letter;

  

  /* find out if letter in cur_letterp is a replacement_letter */
  if(alphabet == 1) {
    if(replacement_letters->nr_rl_1 <= 0) {
      return NO;
    }

    for(i = 0; i < replacement_letters->nr_rl_1; i++) {
      repl_letter = replacement_letters->letters_1 + i;
      same_letter = YES;
      j = 0;
      while(*(repl_letter->letter + j) != '\0' && *(cur_letterp->letter + j) != '\0') {
	if(*(repl_letter->letter + j) == *(cur_letterp->letter + j)) {
	}
	else {
	  same_letter = NO;
	}
	j++;
      }
      if(*(repl_letter->letter + j) != '\0' || *(cur_letterp->letter + j) != '\0') {
	same_letter = NO;
      }
      else if(same_letter == YES) {
	break;
      }
    }
    if(same_letter == NO) {
      return NO;
    }
    else { /* k represents the regular letter, i represents which repl_letter this is */ 
      for(k = 0; k < hmmp->a_size; k++) {
	(msa_seq_infop->msa_seq_1 + get_mtx_index(seq_pos,k, hmmp->a_size+1))->nr_occurences += 
	  *(replacement_letters->probs_1 + get_mtx_index(i,k,hmmp->a_size));
      }
      return YES;
    }
  }

  /* find out if letter in cur_letterp is a replacement_letter */
  if(alphabet == 2) {
    if(replacement_letters->nr_rl_2 <= 0) {
      return NO;
    }
    for(i = 0; i < replacement_letters->nr_rl_2; i++) {
      repl_letter = replacement_letters->letters_2 + i;
      same_letter = YES;
      j = 0;
      while(*(repl_letter->letter + j) != '\0' && *(cur_letterp->letter + j) != '\0') {
	if(*(repl_letter->letter + j) == *(cur_letterp->letter + j)) {
	}
	else {
	  same_letter = NO;
	}
	j++;
      }
      if(*(repl_letter->letter + j) != '\0' || *(cur_letterp->letter + j) != '\0') {
	same_letter = NO;
      }
      else if(same_letter == YES) {
	break;
      }
    }
    if(same_letter == NO) {
      return NO;
    }
    else { /* k represents the regular letter, i represents which repl_letter this is */ 
      for(k = 0; k < hmmp->a_size_2; k++) {
	(msa_seq_infop->msa_seq_2 + get_mtx_index(seq_pos,k, hmmp->a_size_2+1))->nr_occurences += 
	  *(replacement_letters->probs_2 + get_mtx_index(i,k,hmmp->a_size_2));
      }
      return YES;
    }
  }

   /* find out if letter in cur_letterp is a replacement_letter */
  if(alphabet == 3) {
    if(replacement_letters->nr_rl_3 <= 0) {
      return NO;
    }
    for(i = 0; i < replacement_letters->nr_rl_3; i++) {
      repl_letter = replacement_letters->letters_3 + i;
      same_letter = YES;
      j = 0;
      while(*(repl_letter->letter + j) != '\0' && *(cur_letterp->letter + j) != '\0') {
	if(*(repl_letter->letter + j) == *(cur_letterp->letter + j)) {
	}
	else {
	  same_letter = NO;
	}
	j++;
      }
      if(*(repl_letter->letter + j) != '\0' || *(cur_letterp->letter + j) != '\0') {
	same_letter = NO;
      }
      else if(same_letter == YES) {
	break;
      }
    }
    if(same_letter == NO) {
      return NO;
    }
    else { /* k represents the regular letter, i represents which repl_letter this is */ 
      for(k = 0; k < hmmp->a_size_3; k++) {
	(msa_seq_infop->msa_seq_3 + get_mtx_index(seq_pos,k, hmmp->a_size_3+1))->nr_occurences += 
	  *(replacement_letters->probs_3 + get_mtx_index(i,k,hmmp->a_size_3));
      }
      return YES;
    }
  }

   /* find out if letter in cur_letterp is a replacement_letter */
  if(alphabet == 4) {
    if(replacement_letters->nr_rl_4 <= 0) {
      return NO;
    }
    for(i = 0; i < replacement_letters->nr_rl_4; i++) {
      repl_letter = replacement_letters->letters_4 + i;
      same_letter = YES;
      j = 0;
      while(*(repl_letter->letter + j) != '\0' && *(cur_letterp->letter + j) != '\0') {
	if(*(repl_letter->letter + j) == *(cur_letterp->letter + j)) {
	}
	else {
	  same_letter = NO;
	}
	j++;
      }
      if(*(repl_letter->letter + j) != '\0' || *(cur_letterp->letter + j) != '\0') {
	same_letter = NO;
      }
      else if(same_letter == YES) {
	break;
      }
    }
    if(same_letter == NO) {
      return NO;
    }
    else { /* k represents the regular letter, i represents which repl_letter this is */ 
      for(k = 0; k < hmmp->a_size_4; k++) {
	(msa_seq_infop->msa_seq_4 + get_mtx_index(seq_pos,k, hmmp->a_size_4+1))->nr_occurences += 
	  *(replacement_letters->probs_4 + get_mtx_index(i,k,hmmp->a_size_4));
      }
      return YES;
    }
  }
  
  return NO;
}

int get_nr_alphabets(FILE *hmmfile)
{
  int MAX_LINE = 4000;
  int i;
  char s[1000];

  if(hmmfile == NULL) {
    return 0;
  }
  else {
    for(i = 0; i < 3; i++) {
      if(fgets(s, MAX_LINE, hmmfile) != NULL) {
	
      }
      else {
	return 0;
      }
    }
    if(fgets(s, MAX_LINE, hmmfile) != NULL) {
      if(strncmp(s, "NR OF", 5) == 0) {
	  return atoi(&s[17]);
      }
      else {
	return 11;
      }
    }
    else {
	return 0;
    } 
  }
}

void get_set_of_labels(struct hmm_s *hmmp)
{
  int i,j;
  int nr_l;
  char *labels;
  int is_listed;

  labels = (char*)(malloc_or_die(hmmp->nr_v * sizeof(char)));

  nr_l = 0;
  for(i = 1; i < hmmp->nr_v-1; i++) {
    is_listed = NO;
    for(j = 0; j < nr_l; j++) {
      if(hmmp->vertex_labels[i] == *(labels + j)) {
	is_listed = YES;
	break;
      }
    }
    
    if(is_listed == NO) {
      *(labels + nr_l) = hmmp->vertex_labels[i];
      nr_l ++;
    }
  }

  hmmp->labels = labels;
  hmmp->nr_labels = nr_l;
}

void get_set_of_labels_multi(struct hmm_multi_s *hmmp)
{
  int i,j;
  int nr_l;
  char *labels;
  int is_listed;

  labels = (char*)(malloc_or_die(hmmp->nr_v * sizeof(char)));

  nr_l = 0;
  for(i = 1; i < hmmp->nr_v-1; i++) {
    is_listed = NO;
    for(j = 0; j < nr_l; j++) {
      if(hmmp->vertex_labels[i] == *(labels + j)) {
	is_listed = YES;
	break;
      }
    }
    
    if(is_listed == NO) {
      *(labels + nr_l) = hmmp->vertex_labels[i];
      nr_l ++;
    }
  }

  hmmp->labels = labels;
  hmmp->nr_labels = nr_l;
}
