#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
//#include <double.h>
#include <limits.h>

#include "structs.h"
#include "funcs.h"

#define MAX_LINE 60000

#define DIRICHLET 12
#define NONE -2

#define LONGEST_SEQ -1 /* Note: this constant is also defined in hmmsearch_msa.c */

//#define DEBUG_SEQ_STD
//#define DEBUG_SEQ_FASTA
//#define DEBUG_MSA_SEQ_STD
//#define DEBUG_MSA_SEQ_PRF
//#define DEBUG_PRI

extern int verbose;


int seqfile_has_labels(FILE *seqfile)
{
  char row[30000];
  
  rewind(seqfile);
  while(1) {
    if(fgets(row,30000,seqfile) != NULL) {
      if(row[0] == '/') {
	rewind(seqfile);
	return YES;
      }
    }
    else {
      rewind(seqfile);
      return NO;
    }
  }
  
}
     

void get_sequence_fasta_multi(FILE *seqfile, struct sequences_multi_s *seq_infop, int seq_nr)
{
  int i,j,a;
  int nr_letters;
  int longest_seq, shortest_seq, avg_seq_len;
  int inside_seq;
  char c;
  struct letter_s *seqsp;


  nr_letters = 0;
  longest_seq = seq_infop->longest_seq;
  shortest_seq = seq_infop->shortest_seq;
  inside_seq = NO;
  
  while(1) {
    i = fgetc(seqfile);
    if(i == '>') {
      break;
    }
    else if(i == '\s' && i == '\n') {
      
    }
    else {
      printf("Sequence file does not seem to be in correct format\n");
      exit(0);
    }
  }
  rewind(seqfile);

  /* Find out how much space to allocate for the sequence = count total number of letters*/
  while((i = fgetc(seqfile)) != EOF) {
    c = (char)i;
    if(((i >= 65 && i <= 90) || (i >= 97 && i <= 122)) && inside_seq == YES) {
      nr_letters++;
    }
    else if(c == '>') {
      while((i = fgetc(seqfile)) != '\n') {
	if(i == EOF) {
	  printf("Sequence file is not in correct FASTA format\n");
	  exit(0);
	}
      }
      inside_seq = YES;
    }
    else if(c == '/') {
      inside_seq = NO;
    }
  }
#ifdef DEBUG_SEQ_FASTA
  printf("nr_letters = %d\n", nr_letters);
#endif
  /* Allocate memory, must be freed by caller*/
  if(nr_letters > seq_infop->longest_seq) {
    seq_infop->longest_seq = nr_letters;
  }
  if(nr_letters < seq_infop->shortest_seq) {
    seq_infop->shortest_seq = nr_letters;
  }
  seq_infop->avg_seq_len += nr_letters;
  (seq_infop->seqs + seq_nr)->seq_1 = (struct letter_s*)(malloc_or_die(nr_letters * 2 * sizeof(struct letter_s)));
  (seq_infop->seqs + seq_nr)->length = nr_letters;
  //printf("%x\n", (seq_infop->seqs + seq_nr)->seq_1);
  //printf("%x\n", nr_letters * 2 * sizeof(struct letter_s));

/* Read sequences into memory */
  seqsp = (struct letter_s*)((seq_infop->seqs + seq_nr)->seq_1);
  rewind(seqfile);
  inside_seq = NO;
  while((c = (char)(fgetc(seqfile))) != '>' && c != EOF) {
    
  }
  if(c == '>') {
    inside_seq = YES;
    while((c = (char)(fgetc(seqfile))) == ' ') {
      
    }
    /* read sequence name */
    a = 0;
    (seq_infop->seqs + seq_nr)->name[a] = c;
    a++;
    while((c = (char)(fgetc(seqfile))) != ' ' && c != '\n') {
      if(a < MAX_SEQ_NAME_SIZE-1) {
	(seq_infop->seqs + seq_nr)->name[a] = c;
	a++;
      }
    }
    (seq_infop->seqs + seq_nr)->name[a] = '\0';
    if(c != '\n') {
      while((c = (char)(fgetc(seqfile))) != '\n') {
	/* read rest of row and forget about it */
      }
    }
    while((c = (char)(fgetc(seqfile))) != EOF) {
      if(c == '\n' || c == ' ') {
	
      }
      else if(c == '/') {
	break;
      }
      else if ((c >= 65 && c <= 90) || (c >= 97 && c <= 122)) {
	if(c >= 97) {
	  c = c - 32; //convert to upper case
	}
	(seqsp->letter)[0] = c;
	(seqsp->letter)[1] = '\0';
	seqsp->label = '.';
	seqsp++;
	//printf("%x\n", seqsp);
      }
    }
  }
  else {
    printf("Sequence file is not in correct FASTA format\n");
    exit(0);
  }
}

 
void get_sequence_std_multi(FILE *seqfile, struct sequences_multi_s *seq_infop, struct hmm_multi_s *hmmp, int seq_nr)
{
  int i,j,a,k,last;
  int nr_letters;
  int longest_seq, shortest_seq, avg_seq_len;
  int temp_seq_len;
  int inside_seq;
  int nr_alphabets, nr_alphabets_temp;
  char c;
  struct letter_s *seqsp;
  struct sequence_multi_s temp_seq;

  nr_letters = 0;
  inside_seq = NO;
  nr_alphabets = 0;
  nr_alphabets_temp = 0;
  last = '!';
  char line[MAX_LINE], *cur_line;
  double letter_val;
  
  cur_line = line;
  while(fgets(line, MAX_LINE, seqfile) != NULL) {
    if(line[0] == '<' || line[0] == '#' || line[0] == '\s' || line[0] == '\n' || line[0] == '/') {
      
    }
    else {
      printf("Sequence file does not seem to be in correct format\n");
      exit(0);
    }
    for(i = 0; line[i] != '\0'; i++) {
      if((line[i] == '>' || line[i] == '+') && i > 0 && line[i-1] != ';') {
	printf("Sequence file does not seem to be in correct format\n");
	printf("All letters must be followed by ';'\n");
	exit(0);
      }
    }
  }
  rewind(seqfile);
  
  /* Find out how much space to allocate for the sequences = count total number of letters*/
  while((i = fgetc(seqfile)) != EOF) {
    c = (char)i;
    if(c == ';' && inside_seq == YES) {
      nr_letters++;
      last = '!';
    }
    else if(c == '<') {
      inside_seq = YES;
      nr_letters = 0;
      if(last == '<') {
	nr_alphabets_temp++;
      }
      else {
	nr_alphabets_temp = 1;
	last = '<';
      }
    }
    else if(c == '#') {
      inside_seq = YES;
      if(last == '#') {
	nr_alphabets_temp++;
      }
      else {
	nr_alphabets_temp = 1;
	last = '#';
      }
    }
    else if((c == '>' || c == '+')&& inside_seq == YES) {
      inside_seq = NO;
      if(nr_alphabets_temp > nr_alphabets) {
	nr_alphabets = nr_alphabets_temp;
	nr_alphabets_temp = 0;
      }
      if(c == '>') {
	switch(nr_alphabets) {
	case 1: hmmp->alphabet_type = DISCRETE; break;
	case 2: hmmp->alphabet_type_2 = DISCRETE; break;
	case 3: hmmp->alphabet_type_3 = DISCRETE; break;
	case 4: hmmp->alphabet_type_4 = DISCRETE; break;
	}
      }
      if(c == '+') {
	switch(nr_alphabets) {
	case 1: hmmp->alphabet_type = CONTINUOUS; break;
	case 2: hmmp->alphabet_type_2 = CONTINUOUS; break;
	case 3: hmmp->alphabet_type_3 = CONTINUOUS; break;
	case 4: hmmp->alphabet_type_4 = CONTINUOUS; break;
	}
      }
      last = '!';
    }
    else if(c == '/') {
      break;
    }
    else {
      last = '!';
    }
  }

#ifdef DEBUG_SEQ_STD
  printf("nr_letters = %d\n", nr_letters);
  printf("nr_alphabets = %d\n", nr_alphabets);
#endif

  /* check if nr of alphabets in sequence file corresponds to nr of alphabets in hmm file */
  if(hmmp->nr_alphabets < nr_alphabets) {
    printf("Warning: HMM has %d alphabets,  while sequence file contains %d alphabets\n", hmmp->nr_alphabets, nr_alphabets);
    printf("Only the first %d alphabets will be read and used\n", hmmp->nr_alphabets);
    nr_alphabets = hmmp->nr_alphabets;
  }
  else if(hmmp->nr_alphabets > nr_alphabets) {
    printf("Error: HMM has %d alphabets,  while sequence file contains %d alphabets\n", hmmp->nr_alphabets, nr_alphabets);
    exit(0);
  }

  /* Allocate memory, must be freed by caller*/
  if(nr_letters > seq_infop->longest_seq) {
    seq_infop->longest_seq = nr_letters;
  }
  if(nr_letters < seq_infop->shortest_seq) {
    seq_infop->shortest_seq = nr_letters;
  }
  seq_infop->avg_seq_len += nr_letters;
  (seq_infop->seqs + seq_nr)->seq_1 = (struct letter_s*)(malloc_or_die(nr_letters * 2 * sizeof(struct letter_s)));
  if(nr_alphabets > 1) {
    (seq_infop->seqs + seq_nr)->seq_2 = (struct letter_s*)(malloc_or_die(nr_letters * 2 * sizeof(struct letter_s)));
  }
  if(nr_alphabets > 2) {
    (seq_infop->seqs + seq_nr)->seq_3 = (struct letter_s*)(malloc_or_die(nr_letters * 2 * sizeof(struct letter_s)));
  }
  if(nr_alphabets > 3) {
    (seq_infop->seqs + seq_nr)->seq_4 = (struct letter_s*)(malloc_or_die(nr_letters * 2 * sizeof(struct letter_s)));
  }
  (seq_infop->seqs + seq_nr)->length = nr_letters;
  
  /* Read sequences into memory */ 
  rewind(seqfile);
  
  for(k = 0; k < nr_alphabets;) {
    if(k == 0) {
      seqsp = (seq_infop->seqs + seq_nr)->seq_1;
    }
    if(k == 1) {
      seqsp = (seq_infop->seqs + seq_nr)->seq_2;
    }
    if(k == 2) {
      seqsp = (seq_infop->seqs + seq_nr)->seq_3;
    }
    if(k == 3) {
      seqsp = (seq_infop->seqs + seq_nr)->seq_4;
    }
    a = 0;
    c = (char)(fgetc(seqfile));
    if(c != '<' && c != '#') {
      /* not a sequence on this line, just continue */
      if(c != '\n') {
	while((c = (char)(fgetc(seqfile))) != '\n') {
	}
      }
      continue;
    }
    else if(c == '<') /* this line is a sequence */ {
      while((c = (char)(fgetc(seqfile))) != '>') {
	if(c == '<') {
	  
	}
	else if(c == ';') {
	  seqsp->letter[a] = '\0';
	  seqsp->label = '.';
	  seqsp++;
	  a = 0;
	}
	else {
	  seqsp->letter[a] = c;
	  a++;
	}
      }
      
      seqsp->letter[0] = '\0';
      seqsp++;
      strcpy((seq_infop->seqs + seq_nr)->name, "\0");
      if(k == 0) {
	(seq_infop->seqs + seq_nr)->length = get_seq_length((seq_infop->seqs + seq_nr)->seq_1);
      }
      k++;
    }
    else if(c == '#') {
      fgets(line, MAX_LINE, seqfile);
      cur_line = line;
      while(*cur_line == '#') {
	cur_line++;
      }
      while(*cur_line != '+') {
	letter_val = strtod(cur_line, &cur_line);
	seqsp->letter[0] = 'X';
	seqsp->letter[1] = '\0';
	seqsp->label = '.';
	seqsp->cont_letter = letter_val;
	seqsp++;
	if(*cur_line != ';') {
	  printf("Strange continuous sequence file format\n");
	}
	else {
	  cur_line = cur_line + 1;
	}
      }
      seqsp->letter[0] = '\0';
      seqsp++;
      strcpy((seq_infop->seqs + seq_nr)->name, "\0");
      if(k == 0) {
	(seq_infop->seqs + seq_nr)->length = get_seq_length((seq_infop->seqs + seq_nr)->seq_1);
      }
      k++;
    }
  }
  
#ifdef DEBUG_SEQ_STD
  printf("exiting read_seqs\n");
#endif
  
}



/* Note: msa_seq_infop->seq and msa_seq_infop->gaps will be allocated here but must be
 * freed by caller */
void get_sequences_msa_std_multi(FILE *seqfile, FILE *priorfile_a1, struct msa_sequences_multi_s *msa_seq_infop,
				 struct hmm_multi_s *hmmp, int lead_seq, struct replacement_letter_multi_s *replacement_letters)
{
  int i,j,k,l,m;
  int done, inside_seq, seq_pos, letter_pos, a_index, nr_seqs, cur_pos, cur_seq;
  int msa_seq_length, msa_length, seq_length, longest_seq, longest_seq_length;
  int seq_index, nr_lead_columns;
  int gaps_per_column, tot_nr_gaps, nr_gaps;
  double occurences_per_column;
  int get_letter_columns, get_query_letters;
  int use_priordistribution_1, use_priordistribution_2, use_priordistribution_3, use_priordistribution_4;
  char c;
  struct letter_s cur_letter;
  int **cur_posp;
  struct emission_dirichlet_s em_di_1, em_di_2, em_di_3, em_di_4;
  int is_first;
  int use_replacement_letters;
  int is_empty;
  char last;
  int nr_alphabets, nr_alphabets_temp, nr_continuous_alphabets, nr_continuous_alphabets_temp;
  char line[MAX_LINE], *cur_line;
  double letter_val;
  int read_priorfile;
  

  /* find out length of alignment and allocate memory for probability matrix by
   * reading all sequence rows and remembering the longest */
  tot_nr_gaps = 0;
  done = NO;
  inside_seq = NO;
  msa_seq_length = 0;
  nr_seqs = 0;
  msa_length = 0;
  seq_length = 0;
  seq_index = 0;
  longest_seq = -1;
  longest_seq_length = 0;
  is_first = YES;
  
  while(1) {
    i = fgetc(seqfile);
    if(i == '<') {
      break;
    }
    else if(i == '#') {
      break;
    }
    else if(i == '\s' && i == '\n') {
     
    }
    else {
      printf("Sequence file does not seem to be in correct format\n");
      exit(0);
    }
  }
  rewind(seqfile);

  if(priorfile_a1 == NULL) {
    use_priordistribution_1 = NONE;
    use_priordistribution_2 = NONE;
    use_priordistribution_3 = NONE;
    use_priordistribution_4 = NONE;
  }
  else {
    read_priorfile = read_multi_prior_file_multi(&em_di_1, hmmp, priorfile_a1, 1);
    if(read_priorfile > 0) {
      use_priordistribution_1 = DIRICHLET;
    }
    else if(read_priorfile < 0) {
      printf("Error: Incorrect priorfile format\n");
      exit(0);
    }
    else {
      use_priordistribution_1 = NONE;
    }
    if(hmmp->nr_alphabets > 1) {
      read_priorfile = read_multi_prior_file_multi(&em_di_2, hmmp, priorfile_a1, 2);
      if(read_priorfile > 0) {
	use_priordistribution_2 = DIRICHLET;
      }
      else if(read_priorfile < 0) {
	printf("Error: Incorrect priorfile format\n");
	exit(0);
      }
      else {
	use_priordistribution_2 = NONE;
      }
    }
    if(hmmp->nr_alphabets > 2) {
      read_priorfile = read_multi_prior_file_multi(&em_di_3, hmmp, priorfile_a1, 3);
      if(read_priorfile > 0) {
	use_priordistribution_3 = DIRICHLET;
      }
      else if(read_priorfile < 0) {
	printf("Error: Incorrect priorfile format\n");
	exit(0);
      }
      else {
	use_priordistribution_3 = NONE;
      }
    }
    if(hmmp->nr_alphabets > 3) {
      read_priorfile = read_multi_prior_file_multi(&em_di_4, hmmp, priorfile_a1, 4);
      if(read_priorfile > 0) {
	use_priordistribution_4 = DIRICHLET;
      }
      else if(read_priorfile < 0) {
	printf("Error: Incorrect priorfile format\n");
	exit(0);
      }
      else {
	use_priordistribution_4 = NONE;
      }
    }
  }
   
  if(replacement_letters->nr_alphabets == 0) {
    use_replacement_letters = NO;
  }
  else {
    use_replacement_letters = YES;
  }

  /* check if file is empty */
  is_empty = YES;
  while(done != YES) {
    c = (char)fgetc(seqfile);
    if((int)c == EOF) {
      break;
    }
    else if (c == '<' || c == '#') {
      is_empty = NO;
      break;
    }
  }
  if(is_empty == YES) {
    if(verbose == YES) {
      printf("File is empty\n");
    } 
  }
  else {
  }
  rewind(seqfile);

  last = '!';
  nr_alphabets = 0;
  nr_alphabets_temp = 0;
  nr_continuous_alphabets = 0;
  nr_continuous_alphabets_temp = 0;
  while(done != YES) {
    c = (char)fgetc(seqfile);
    if((int)c == EOF) {
      done = YES;
      continue;
    }
    if(c == '<' && last != '<') {
      inside_seq = YES;
      seq_length = 0;
      if(is_first == YES) {
	msa_seq_length = 0;
      }
      seq_index++;
      nr_seqs++;
      last = '<';
      nr_alphabets_temp = 1;
    }
    else if(c == '#' && last != '#') {
      inside_seq = YES;
      seq_length = 0;
      if(is_first == YES) {
	msa_seq_length = 0;
      }
      seq_index++;
      last = '#';
      nr_alphabets_temp = 1;
    }
    else if(c == '<') {
      nr_alphabets_temp++;
    }
    else if(c == '#') {
      nr_alphabets_temp++;
      nr_continuous_alphabets_temp++;
    }
    else if(c == '>') {
      if(seq_length > longest_seq_length) {
	longest_seq = seq_index;
	longest_seq_length = seq_length;
      }
      inside_seq = NO;
      is_first = NO;
      if(nr_alphabets_temp > nr_alphabets) {
	nr_alphabets = nr_alphabets_temp;
	switch(nr_alphabets) {
	case 1: hmmp->alphabet_type = DISCRETE; break;
	case 2: hmmp->alphabet_type_2 = DISCRETE; break;
	case 3: hmmp->alphabet_type_3 = DISCRETE; break;
	case 4: hmmp->alphabet_type_4 = DISCRETE; break;
	}
      }
      
      nr_alphabets_temp = 0;
      last = '!';
    }
    else if(c == '+') {
      if(seq_length > longest_seq_length) {
	longest_seq = seq_index;
	longest_seq_length = seq_length;
      }
      inside_seq = NO;
      is_first = NO;
      if(nr_alphabets_temp > nr_alphabets) {
	nr_alphabets = nr_alphabets_temp;
	switch(nr_alphabets) {
	case 1: hmmp->alphabet_type = CONTINUOUS; break;
	case 2: hmmp->alphabet_type_2 = CONTINUOUS; break;
	case 3: hmmp->alphabet_type_3 = CONTINUOUS; break;
	case 4: hmmp->alphabet_type_4 = CONTINUOUS; break;
	}	
      }
      if(nr_continuous_alphabets_temp > nr_continuous_alphabets) {
	nr_continuous_alphabets = nr_continuous_alphabets_temp;
      }
      
      
      nr_alphabets_temp = 0;
      nr_continuous_alphabets_temp = 0;
      last = '!';
    }
    if(c == ';' && inside_seq == YES) {
      if(is_first == YES) {
	msa_seq_length++;
      }
      seq_length++;
      last = '!';
    }
    else if(inside_seq == YES && (c == '_' || c == ' ' || c == '-' || c == '.')) {
      seq_length--;
      last = '!';
    }
  }
  nr_seqs = nr_seqs / (nr_alphabets - nr_continuous_alphabets);
  msa_seq_infop->nr_seqs = nr_seqs;

#ifdef DEBUG_MSA_SEQ_STD
  printf("reached checkpoint 1\n");
  printf("msa_seq_length = %d\n", msa_seq_length);
  printf("nr_alphabets = %d\n", nr_alphabets);
  printf("longest seq = %d\n", longest_seq);
  printf("longest_seq_length = %d\n", longest_seq_length);
  printf("nr_seqs = %d\n", nr_seqs);
#endif

  msa_seq_infop->msa_seq_1 = (struct msa_letter_s*)
    malloc_or_die(msa_seq_length * (hmmp->a_size+1) * sizeof(struct msa_letter_s));
  msa_seq_infop->lead_columns_start = (int*)malloc_or_die((msa_seq_length+1) * sizeof(int));
  if(hmmp->nr_alphabets > 1){
    msa_seq_infop->msa_seq_2 = (struct msa_letter_s*)
      malloc_or_die(msa_seq_length * (hmmp->a_size_2+1) * sizeof(struct msa_letter_s));
  }
  if(hmmp->nr_alphabets > 2){
    msa_seq_infop->msa_seq_3 = (struct msa_letter_s*)
      malloc_or_die(msa_seq_length * (hmmp->a_size_3+1) * sizeof(struct msa_letter_s));
  }
  if(hmmp->nr_alphabets > 3){
    msa_seq_infop->msa_seq_4 = (struct msa_letter_s*)
      malloc_or_die(msa_seq_length * (hmmp->a_size_4+1) * sizeof(struct msa_letter_s));
  }
  
  /* read first alphabet, save nr of occurences of each letter in each position,
   * including nr of gaps */
  rewind(seqfile);
  seq_pos = 0;
  inside_seq = NO;
  seq_index = 0;
  get_letter_columns = NO;
  k = 0;
  l = 0;
  nr_lead_columns = 0;
  if(lead_seq == LONGEST_SEQ) {
    lead_seq = longest_seq;
  }
  
  for(m = 0; m < nr_seqs;) {
    i = fgetc(seqfile);
    if(i == EOF) {
      break;
      
    }
    c = (char)i;
    if(c == '<') {
      seq_pos = 0;
      inside_seq = YES;
      seq_index++;
      if(seq_index == lead_seq) {
	get_letter_columns = YES;
	get_query_letters = YES;
      }
    }
    else if(c == '#') {
      get_letter_columns = YES;
      fgets(line, MAX_LINE, seqfile);
      cur_line = line;
      while(*cur_line == '#') {
	cur_line++;
      }
      while(*cur_line != '+') {
	if(*cur_line == '_' || *cur_line == ' ' ||
	   (*cur_line == '-' && *(cur_line + 1) == ';') || *cur_line == '.') /* gap */ {
	  (msa_seq_infop->msa_seq_1 + get_mtx_index(seq_pos, 0, hmmp->a_size+1))->share = 0.0;
	  (msa_seq_infop->msa_seq_1 + get_mtx_index(seq_pos, 0, hmmp->a_size+1))->nr_occurences = 0.0;
	  cur_line++;
	}
	else if(*cur_line == 'X') {
	  (msa_seq_infop->msa_seq_1 + get_mtx_index(seq_pos, 0, hmmp->a_size+1))->share = 0.0;
	  (msa_seq_infop->msa_seq_1 + get_mtx_index(seq_pos, 0, hmmp->a_size+1))->nr_occurences = -1.0;
	  cur_line++;
	}
	else {
	  letter_val = strtod(cur_line, &cur_line);
	  (msa_seq_infop->msa_seq_1 + get_mtx_index(seq_pos,0, hmmp->a_size+1))->share = letter_val;
	  (msa_seq_infop->msa_seq_1 + get_mtx_index(seq_pos,0, hmmp->a_size+1))->nr_occurences = 1.0;
	  if(get_letter_columns == YES) {
	    *(msa_seq_infop->lead_columns_start + k) = seq_pos; /* letter column */
	    k++;
	    nr_lead_columns++;
	  }
	}
	if(*cur_line != ';') {
	  printf("Strange continuous sequence file format\n");
	  exit(0);
	}
	else {
	  cur_line = cur_line + 1;
	  seq_pos++;
	}
      }
      get_letter_columns = NO;
      break;
    }
    else if(c == '>') {
      inside_seq = NO;
      get_letter_columns = NO;
      get_query_letters = NO;
      m++;
    }
    else if(inside_seq == YES) {
      /* reading one letter */
      letter_pos = 0;
      cur_letter.letter[letter_pos] = c;
      letter_pos++;
      while((c = (char)fgetc(seqfile)) != ';') {
	cur_letter.letter[letter_pos] = c;
	letter_pos++;
      }
      if(letter_pos > 4) {
	printf("max letter size = 4 characters\n");
	exit(0);
      }
      cur_letter.letter[letter_pos] = '\0';
      if(seq_pos + 1 > msa_length) {
	msa_length = seq_pos + 1; /* update nr of positions */
      }
      if(cur_letter.letter[0] == '_' || cur_letter.letter[0] == ' ' ||
	 cur_letter.letter[0] == '-' || cur_letter.letter[0] == '.') /* gap */ {
	(msa_seq_infop->msa_seq_1 + get_mtx_index(seq_pos, hmmp->a_size, hmmp->a_size+1))->nr_occurences
	  += 1;
	seq_pos++;
      }
      else /* non gap character */ {
	if(use_replacement_letters == YES &&
	   replacement_letter_multi(&cur_letter, replacement_letters, msa_seq_infop, hmmp, seq_pos,1) == YES){
	  /* adding of nr_occurences for these letters is done on the fly in replacement_letter-function */
	}
	else /* add occurence of this letter */{
	  a_index = get_alphabet_index(&cur_letter, hmmp->alphabet, hmmp->a_size);
	  //printf("a_index = %d\n", a_index);
	  if(a_index < 0) {
	    printf("Could not read msa seq from file: letter '%s' is not in alphabet\n", cur_letter.letter);
	    exit(0);
	  }
	  else {
	    (msa_seq_infop->msa_seq_1 + get_mtx_index(seq_pos,a_index, hmmp->a_size+1))->nr_occurences += 1;
	  }
	}
	if(get_letter_columns == YES) {
	  *(msa_seq_infop->lead_columns_start + k) = seq_pos; /* letter column */
	  k++;
	  nr_lead_columns++;
	}
	seq_pos++;
      }
      if(get_query_letters == YES) {
	memcpy((msa_seq_infop->msa_seq_1 + l * (hmmp->a_size + 1))->query_letter, &(cur_letter.letter),
	       sizeof(char) * 5); /* query letter */
	l++;
      }
    }
    //printf("m = %d\n",m);
  }
  /* read second alphabet, save nr of occurences of each letter in each position,
   * including nr of gaps */
  if(hmmp->nr_alphabets > 1) {
    seq_pos = 0;
    inside_seq = NO;
    seq_index = 0;
    get_letter_columns = NO;
    l = 0;
    last = '!';
    for(m = 0; m < nr_seqs;) {
      i = fgetc(seqfile);
      if(i == EOF) {
	break;
      }
      c = (char)i;
      //printf("alpha2:c = %c\n",c);
      
      if(c == '<' && last != '<') {
	seq_pos = 0;
	inside_seq = YES;
	seq_index++;
	if(seq_index == lead_seq) {
	  get_query_letters = YES;
	}
	last = c;
      }
      else if(c == '<') {
	last = c;
      }
      else if(c == '#') {
	fgets(line, MAX_LINE, seqfile);
	cur_line = line;
	while(*cur_line == '#') {
	  cur_line++;
	}
	while(*cur_line != '+') {
	  if(*cur_line == '_' || *cur_line == ' ' ||
	     (*cur_line == '-' && *(cur_line + 1) == ';') || *cur_line == '.') /* gap */ {
	    (msa_seq_infop->msa_seq_2 + get_mtx_index(seq_pos, 0, hmmp->a_size_2+1))->share = 0.0;
	    (msa_seq_infop->msa_seq_2 + get_mtx_index(seq_pos, 0, hmmp->a_size_2+1))->nr_occurences = 0.0;
	    cur_line++;
	  }
	  else if(*cur_line == 'X') {
	    (msa_seq_infop->msa_seq_2 + get_mtx_index(seq_pos, 0, hmmp->a_size_2+1))->share = 0.0;
	    (msa_seq_infop->msa_seq_2 + get_mtx_index(seq_pos, 0, hmmp->a_size_2+1))->nr_occurences = -1.0;
	    cur_line++;
	  }
	  else {
	    letter_val = strtod(cur_line, &cur_line);
	    (msa_seq_infop->msa_seq_2 + get_mtx_index(seq_pos,0, hmmp->a_size_2+1))->share = letter_val;
	    (msa_seq_infop->msa_seq_2 + get_mtx_index(seq_pos,0, hmmp->a_size_2+1))->nr_occurences = 1.0;
	  }
	  if(*cur_line != ';') {
	    printf("Strange continuous sequence file format\n");
	    exit(0);
	  }
	  else {
	    cur_line = cur_line + 1;
	    seq_pos++;
	  }
	}
	break;
      }
      else if(c == '>') {
	inside_seq = NO;
	get_query_letters = NO;
	m++;
	last = '!';
      }
      else if(inside_seq == YES) {
	/* reading one letter */
	letter_pos = 0;
	cur_letter.letter[letter_pos] = c;
	letter_pos++;
	while((c = (char)fgetc(seqfile)) != ';') {
	  cur_letter.letter[letter_pos] = c;
	  letter_pos++;
	}
	if(letter_pos > 4) {
	  printf("max letter size = 4 characters\n");
	  exit(0);
	}
	cur_letter.letter[letter_pos] = '\0';
	if(seq_pos + 1 > msa_length) {
	  msa_length = seq_pos + 1; /* update nr of positions */
	}
	if(cur_letter.letter[0] == '_' || cur_letter.letter[0] == ' ' ||
	   cur_letter.letter[0] == '-' || cur_letter.letter[0] == '.') /* gap */ {
	  (msa_seq_infop->msa_seq_2 + get_mtx_index(seq_pos, hmmp->a_size_2, hmmp->a_size_2+1))->nr_occurences
	    += 1;
	  seq_pos++;
	}
	else /* non gap character */ {
	  if(use_replacement_letters == YES &&
	     replacement_letter_multi(&cur_letter, replacement_letters, msa_seq_infop, hmmp, seq_pos,2) == YES){
	    /* adding of nr_occurences for these letters is done on the fly in replacement_letter-function */
	  }
	  else /* add occurence of this letter */{
	    a_index = get_alphabet_index(&cur_letter, hmmp->alphabet_2, hmmp->a_size_2);
	    if(a_index < 0) {
	      printf("Could not read msa seq from file: letter '%s' is not in alphabet\n", cur_letter.letter);
	      exit(0);
	    }
	    else {
	      (msa_seq_infop->msa_seq_2 + get_mtx_index(seq_pos,a_index, hmmp->a_size_2+1))->nr_occurences += 1;
	    }
	  }
	  seq_pos++;
	}
	if(get_query_letters == YES) {
	  memcpy((msa_seq_infop->msa_seq_2 + l * (hmmp->a_size_2 + 1))->query_letter, &(cur_letter.letter),
		 sizeof(char) * 5); /* query letter */
	  l++;
	}
      }
    }
  }

  /* read third alphabet, save nr of occurences of each letter in each position,
   * including nr of gaps */
  if(hmmp->nr_alphabets > 2) {
    seq_pos = 0;
    inside_seq = NO;
    seq_index = 0;
    get_letter_columns = NO;
    k = 0;
    l = 0;
    nr_lead_columns = 0;
    last = '!';
    for(m = 0; m < nr_seqs;) {
      i = fgetc(seqfile);
      if(i == EOF) {
	break;
      }
      c = (char)i;
      //printf("alpha3:c = %c\n",c);
      
      if(c == '<' && last != '<') {
	seq_pos = 0;
	inside_seq = YES;
	seq_index++;
	if(seq_index == lead_seq) {
	  get_query_letters = YES;
	}
	last = c;
      }
      else if(c == '<') {
	last = c;
      }
      else if(c == '#') {
	fgets(line, MAX_LINE, seqfile);
	cur_line = line;
	while(*cur_line == '#') {
	  cur_line++;
	}
	while(*cur_line != '+') {
	  if(*cur_line == '_' || *cur_line == ' ' ||
	     (*cur_line == '-' && *(cur_line + 1) == ';') || *cur_line == '.') /* gap */ {
	    (msa_seq_infop->msa_seq_3 + get_mtx_index(seq_pos, 0, hmmp->a_size_3+1))->share = 0.0;
	    (msa_seq_infop->msa_seq_3 + get_mtx_index(seq_pos, 0, hmmp->a_size_3+1))->nr_occurences = 0.0;
	    seq_pos++;
	    cur_line++;
	  }
	  else if(*cur_line == 'X') {
	    (msa_seq_infop->msa_seq_3 + get_mtx_index(seq_pos, 0, hmmp->a_size+3))->share = 0.0;
	    (msa_seq_infop->msa_seq_3 + get_mtx_index(seq_pos, 0, hmmp->a_size+3))->nr_occurences = -1.0;
	    seq_pos++;
	    cur_line++;
	  }
	  else {
	    letter_val = strtod(cur_line, &cur_line);
	    (msa_seq_infop->msa_seq_3 + get_mtx_index(seq_pos,0, hmmp->a_size_3+1))->share = letter_val;
	    (msa_seq_infop->msa_seq_3 + get_mtx_index(seq_pos,0, hmmp->a_size_3+1))->nr_occurences = 1.0;
	  }
	  if(*cur_line != ';') {
	    printf("Strange continuous sequence file format\n");
	    exit(0);
	  }
	  else {
	    cur_line = cur_line + 1;
	    seq_pos++;
	  }
	}
	break;
      }
      else if(c == '>') {
	inside_seq = NO;
	get_query_letters = NO;
	m++;
	last = '!';
      }
      else if(inside_seq == YES) {
	/* reading one letter */
	letter_pos = 0;
	cur_letter.letter[letter_pos] = c;
	letter_pos++;
	while((c = (char)fgetc(seqfile)) != ';') {
	  cur_letter.letter[letter_pos] = c;
	  letter_pos++;
	}
	if(letter_pos > 4) {
	  printf("max letter size = 4 characters\n");
	  exit(0);
	}
	cur_letter.letter[letter_pos] = '\0';
	if(seq_pos + 1 > msa_length) {
	  msa_length = seq_pos + 1; /* update nr of positions */
	}
	if(cur_letter.letter[0] == '_' || cur_letter.letter[0] == ' ' ||
	   cur_letter.letter[0] == '-' || cur_letter.letter[0] == '.') /* gap */ {
	  (msa_seq_infop->msa_seq_3 + get_mtx_index(seq_pos, hmmp->a_size_3, hmmp->a_size_3+1))->nr_occurences
	    += 1;
	  seq_pos++;
	}
	else /* non gap character */ {
	  if(use_replacement_letters == YES &&
	     replacement_letter_multi(&cur_letter, replacement_letters, msa_seq_infop, hmmp, seq_pos,3) == YES){
	    /* adding of nr_occurences for these letters is done on the fly in replacement_letter-function */
	  }
	  else /* add occurence of this letter */{
	    a_index = get_alphabet_index(&cur_letter, hmmp->alphabet_3, hmmp->a_size_3);
	    if(a_index < 0) {
	      printf("Could not read msa seq from file: letter '%s' is not in alphabet\n", cur_letter.letter);
	      exit(0);
	    }
	    else {
	      (msa_seq_infop->msa_seq_3 + get_mtx_index(seq_pos,a_index, hmmp->a_size_3+1))->nr_occurences += 1;
	    }
	  }
	  seq_pos++;
	}
	if(get_query_letters == YES) {
	  memcpy((msa_seq_infop->msa_seq_3 + l * (hmmp->a_size_3 + 1))->query_letter, &(cur_letter.letter),
		 sizeof(char) * 5); /* query letter */
	  l++;
	}
      }
    }
  }

  /* read fourth alphabet, save nr of occurences of each letter in each position,
   * including nr of gaps */
  if(hmmp->nr_alphabets > 3) {
    seq_pos = 0;
    inside_seq = NO;
    seq_index = 0;
    k = 0;
    l = 0;
    nr_lead_columns = 0;
    last = '!';
    for(m = 0; m < nr_seqs;) {
      i = fgetc(seqfile);
      if(i == EOF) {
	break;
      }
      c = (char)i;
      //printf("alpha4:c = %c\n",c);
      
      if(c == '<' && last != '<') {
	seq_pos = 0;
	inside_seq = YES;
	seq_index++;
	if(seq_index == lead_seq) {
	  get_query_letters = YES;
	}
	last = c;
      }
      else if(c == '<') {
	last = c;
      }
       else if(c == '#') {
	fgets(line, MAX_LINE, seqfile);
	cur_line = line;
	while(*cur_line == '#') {
	  cur_line++;
	}
	while(*cur_line != '+') {
	  if(*cur_line == '_' || *cur_line == ' ' ||
	     (*cur_line == '-' && *(cur_line + 1) == ';') || *cur_line == '.') /* gap */ {
	    (msa_seq_infop->msa_seq_4 + get_mtx_index(seq_pos, 0, hmmp->a_size_4+1))->share = 0.0;
	    (msa_seq_infop->msa_seq_4 + get_mtx_index(seq_pos, 0, hmmp->a_size_4+1))->nr_occurences = 0.0;
	    seq_pos++;
	    cur_line++;
	  }
	  else if(*cur_line == 'X') {
	    (msa_seq_infop->msa_seq_4 + get_mtx_index(seq_pos, 0, hmmp->a_size_4+1))->share = 0.0;
	    (msa_seq_infop->msa_seq_4 + get_mtx_index(seq_pos, 0, hmmp->a_size_4+1))->nr_occurences = -1.0;
	    seq_pos++;
	    cur_line++;
	  }
	  else {
	    letter_val = strtod(cur_line, &cur_line);
	    (msa_seq_infop->msa_seq_4 + get_mtx_index(seq_pos,0, hmmp->a_size_4+1))->share = letter_val;
	    (msa_seq_infop->msa_seq_4 + get_mtx_index(seq_pos,0, hmmp->a_size_4+1))->nr_occurences = 1.0;
	  }
	  if(*cur_line != ';') {
	    printf("Strange continuous sequence file format\n");
	    exit(0);
	  }
	  else {
	    cur_line = cur_line + 1;
	    seq_pos++;
	  }
	}
	break;
      }
      else if(c == '>') {
	inside_seq = NO;
	get_query_letters = NO;
	m++;
	last = '!';
      }
      else if(inside_seq == YES) {
	/* reading one letter */
	letter_pos = 0;
	cur_letter.letter[letter_pos] = c;
	letter_pos++;
	while((c = (char)fgetc(seqfile)) != ';') {
	  cur_letter.letter[letter_pos] = c;
	  letter_pos++;
	}
	if(letter_pos > 4) {
	  printf("max letter size = 4 characters\n");
	  exit(0);
	}
	cur_letter.letter[letter_pos] = '\0';
	if(seq_pos + 1 > msa_length) {
	  msa_length = seq_pos + 1; /* update nr of positions */
	}
	if(cur_letter.letter[0] == '_' || cur_letter.letter[0] == ' ' ||
	   cur_letter.letter[0] == '-' || cur_letter.letter[0] == '.') /* gap */ {
	  (msa_seq_infop->msa_seq_4 + get_mtx_index(seq_pos, hmmp->a_size_4, hmmp->a_size_4+1))->nr_occurences
	    += 1;
	  seq_pos++;
	}
	else /* non gap character */ {
	  if(use_replacement_letters == YES &&
	     replacement_letter_multi(&cur_letter, replacement_letters, msa_seq_infop, hmmp, seq_pos,4) == YES){
	    /* adding of nr_occurences for these letters is done on the fly in replacement_letter-function */
	  }
	  else /* add occurence of this letter */{
	    a_index = get_alphabet_index(&cur_letter, hmmp->alphabet_4, hmmp->a_size_4);
	    if(a_index < 0) {
	      printf("Could not read msa seq from file: letter '%s' is not in alphabet\n", cur_letter.letter);
	      exit(0);
	    }
	    else {
	      (msa_seq_infop->msa_seq_4 + get_mtx_index(seq_pos,a_index, hmmp->a_size_4+1))->nr_occurences += 1;
	    }
	  }
	  seq_pos++;
	}
	if(get_query_letters == YES) {
	  memcpy((msa_seq_infop->msa_seq_4 + l * (hmmp->a_size_4 + 1))->query_letter, &(cur_letter.letter),
		 sizeof(char) * 5); /* query letter */
	  l++;
	}
      }
    }
  }
  
#ifdef DEBUG_MSA_SEQ_STD
  printf("reached checkpoint 2\n");
  printf("total_nr_gaps = %d\n", tot_nr_gaps);
#endif

  *(msa_seq_infop->lead_columns_start + k) = END;
  msa_seq_infop->lead_columns_end = msa_seq_infop->lead_columns_start + k;
  msa_seq_infop->nr_lead_columns = nr_lead_columns;
  //msa_seq_infop->nr_seqs = nr_seqs;
  msa_seq_infop->msa_seq_length = msa_length;
  msa_seq_infop->gap_shares = (double*)malloc_or_die(msa_length * sizeof(double));
 
  tot_nr_gaps = 0;
  gaps_per_column = 0;
  /* go through each position and calculate distribution for all letters in alphabet 1 */
  for(i = 0; i < msa_length; i++) {
    occurences_per_column = 0;
    gaps_per_column = 0;
    for(j = 0; j < hmmp->a_size+1; j++) {
      occurences_per_column += (msa_seq_infop->msa_seq_1 +
				      get_mtx_index(i,j, hmmp->a_size+1))->nr_occurences;
      if(j == hmmp->a_size) {
	tot_nr_gaps += (int)(msa_seq_infop->msa_seq_1 +
			     get_mtx_index(i,j, hmmp->a_size+1))->nr_occurences;
	gaps_per_column = (int)(msa_seq_infop->msa_seq_1 +
				get_mtx_index(i,j, hmmp->a_size+1))->nr_occurences;
      }
    }
    if(use_priordistribution_1 == DIRICHLET) /* calculate share update using dirichlet prior mixture */ {
      update_shares_prior_multi(&em_di_1, hmmp, msa_seq_infop, i,1);
      
    }
    /* simple share update just using the standard quotient */ 
    for(j = 0; j < hmmp->a_size+1; j++) {
      if(hmmp->alphabet_type == DISCRETE) {
	(msa_seq_infop->msa_seq_1 + get_mtx_index(i,j, hmmp->a_size+1))->share =
	  (msa_seq_infop->msa_seq_1 + get_mtx_index(i,j, hmmp->a_size+1))->nr_occurences / 
	  occurences_per_column;
	(msa_seq_infop->msa_seq_1 + get_mtx_index(i,j, hmmp->a_size+1))->label = '.';
      }
    }    
    
    /* calculate gap_share */
    if(hmmp->alphabet_type == DISCRETE) {
      *(msa_seq_infop->gap_shares + i) = (double)(gaps_per_column) / occurences_per_column;
    }
  }

  
  /* go through each position and calculate distribution for all letters in alphabet 2 */
  if(hmmp->nr_alphabets > 1) {
    for(i = 0; i < msa_length; i++) {
      occurences_per_column = 0;
      gaps_per_column = 0;
      for(j = 0; j < hmmp->a_size_2+1; j++) {
	occurences_per_column += (msa_seq_infop->msa_seq_2 +
				      get_mtx_index(i,j, hmmp->a_size_2+1))->nr_occurences;
      }
      if(use_priordistribution_2 == DIRICHLET) /* calculate share update using dirichlet prior mixture */ {
	update_shares_prior_multi(&em_di_2, hmmp, msa_seq_infop, i,2);
	//printf("inne if\n");
      }
      /* simple share update just using the standard quotient */ 
      for(j = 0; j < hmmp->a_size_2+1; j++) {
	if(hmmp->alphabet_type_2 == DISCRETE) {
	  (msa_seq_infop->msa_seq_2 + get_mtx_index(i,j, hmmp->a_size_2+1))->share =
	    (msa_seq_infop->msa_seq_2 + get_mtx_index(i,j, hmmp->a_size_2+1))->nr_occurences / 
	    occurences_per_column;
	  (msa_seq_infop->msa_seq_2 + get_mtx_index(i,j, hmmp->a_size_2+1))->label = '.';
	}
      }    
    }
  }
  /* go through each position and calculate distribution for all letters in alphabet 3 */
  if(hmmp->nr_alphabets > 2) {
    for(i = 0; i < msa_length; i++) {
      occurences_per_column = 0;
      gaps_per_column = 0;
      for(j = 0; j < hmmp->a_size_3+1; j++) {
	occurences_per_column += (msa_seq_infop->msa_seq_3 +
				  get_mtx_index(i,j, hmmp->a_size_3+1))->nr_occurences;
      }
      if(use_priordistribution_3 == DIRICHLET) /* calculate share update using dirichlet prior mixture */ {
	update_shares_prior_multi(&em_di_3, hmmp, msa_seq_infop, i,3);
      }
      /* simple share update just using the standard quotient */ 
      for(j = 0; j < hmmp->a_size_3+1; j++) {
	if(hmmp->alphabet_type_3 == DISCRETE) {
	  (msa_seq_infop->msa_seq_3 + get_mtx_index(i,j, hmmp->a_size_3+1))->share =
	    (msa_seq_infop->msa_seq_3 + get_mtx_index(i,j, hmmp->a_size_3+1))->nr_occurences / 
	    occurences_per_column;
	  (msa_seq_infop->msa_seq_3 + get_mtx_index(i,j, hmmp->a_size_3+1))->label = '.';
	}    
      }
    }
  }
  
  /* go through each position and calculate distribution for all letters in alphabet 4 */
  if(hmmp->nr_alphabets > 3) {
    for(i = 0; i < msa_length; i++) {
      occurences_per_column = 0;
      gaps_per_column = 0;
      for(j = 0; j < hmmp->a_size_4+1; j++) {
	occurences_per_column += (msa_seq_infop->msa_seq_4 +
				  get_mtx_index(i,j, hmmp->a_size_4+1))->nr_occurences;
      }
      if(use_priordistribution_4 == DIRICHLET) /* calculate share update using dirichlet prior mixture */ {
	update_shares_prior_multi(&em_di_4, hmmp, msa_seq_infop, i,4);
      }
      /* simple share update just using the standard quotient */ 
      for(j = 0; j < hmmp->a_size_4+1; j++) {
	if(hmmp->alphabet_type_4 == DISCRETE) {
	  (msa_seq_infop->msa_seq_4 + get_mtx_index(i,j, hmmp->a_size_4+1))->share =
	    (msa_seq_infop->msa_seq_4 + get_mtx_index(i,j, hmmp->a_size_4+1))->nr_occurences / 
	    occurences_per_column;
	  (msa_seq_infop->msa_seq_4 + get_mtx_index(i,j, hmmp->a_size_4+1))->label = '.';
	}
      }    
    }
  }
  
#ifdef DEBUG_MSA_SEQ_STD
  printf("reached checkpoint 3\n");
  printf("total_nr_gaps = %d\n", tot_nr_gaps);
 
#endif
/* allocate memory for gaps array and set initial pointers for each position*/
  msa_seq_infop->gaps = (int**)malloc_or_die(msa_length * sizeof(int*) +
					     (tot_nr_gaps + msa_length) * sizeof(int));
  nr_gaps = 0;
  for(i = 0; i < msa_length; i++) {
    *(msa_seq_infop->gaps + i) = (int*)(msa_seq_infop->gaps + msa_length) + nr_gaps;
    nr_gaps += (int)(msa_seq_infop->msa_seq_1 + get_mtx_index(i, hmmp->a_size, hmmp->a_size+1))->nr_occurences;
    nr_gaps += 1;
  }
#ifdef DEBUG_MSA_SEQ_STD
  printf("%x\n", msa_seq_infop->gaps);
  printf("%x\n", *(msa_seq_infop->gaps));
  printf("%x\n", *(msa_seq_infop->gaps + 1));
  printf("%x\n", *(msa_seq_infop->gaps + 2));
  printf("%x\n", *(msa_seq_infop->gaps + 3));
  printf("nr_seqs: %d\n", nr_seqs);
#endif
  /********** verify that this is done correctly even if there are multiple alphabets */ 

  /* go through every sequence and store which sequences have gaps at what
   * positions in gaps array */
  rewind(seqfile);
  inside_seq = NO;
  seq_pos = 0;
  cur_posp = msa_seq_infop->gaps;
  last = '!';
  for(i = 0; i < nr_seqs;) {
    c = (char)fgetc(seqfile);
    if(c == '<' && last == '<') {
      break;
    }
    else if(c == '<') {
      inside_seq = YES;
      cur_posp = msa_seq_infop->gaps;
      last = '<';
    }
    else if(c == '>') {
      inside_seq = NO;
      i++;
      last = '!';
    }
    else if(inside_seq == YES && (c == '_' || c == ' ' || c == '-' || c == '.')) {
      (**cur_posp) = i+1;
      
#ifdef DEBUG_MSA_SEQ_STD
      printf("posp = %x\n", cur_posp);
      printf("*posp = %x\n", *cur_posp);
      printf("**posp = %x\n", **cur_posp);
#endif
      (*cur_posp)++;
      cur_posp++;
      while((c = (char)fgetc(seqfile)) != ';') {
	
      }
      last = '!';
    }
    else if(inside_seq == YES) {
      while((c = (char)fgetc(seqfile)) != ';') {
	
      }
      cur_posp++;
      last = '!';
    }
  }
  cur_posp = msa_seq_infop->gaps;
  for(i = 0; i < msa_seq_infop->msa_seq_length; i++) {
     (**cur_posp) = END;
     cur_posp++;
  }

  nr_gaps = 0;
  for(i = 0; i < msa_length; i++) {
    *(msa_seq_infop->gaps + i) = (int*)(msa_seq_infop->gaps + msa_length) + nr_gaps;
    nr_gaps += (int)(msa_seq_infop->msa_seq_1 + get_mtx_index(i, hmmp->a_size, hmmp->a_size+1))->nr_occurences;
    nr_gaps += 1;
  }


  
#ifdef DEBUG_MSA_SEQ_STD
  dump_msa_seqs_multi(msa_seq_infop, hmmp);
  printf("reached checkpoint 4\n");
  exit(0);
#endif

  /* cleanup and return */
  if(use_priordistribution_1 == DIRICHLET) {
    if(em_di_1.nr_components > 0) {
      free(em_di_1.q_values);
      free(em_di_1.alpha_sums);
      free(em_di_1.logbeta_values);
      free(em_di_1.prior_values);
    }
  }
  if(use_priordistribution_2 == DIRICHLET) {
    if(em_di_2.nr_components > 0) {
      free(em_di_2.q_values);
      free(em_di_2.alpha_sums);
      free(em_di_2.logbeta_values);
      free(em_di_2.prior_values);
    }
  }
  if(use_priordistribution_3 == DIRICHLET) {
    if(em_di_3.nr_components > 0) {
      free(em_di_3.q_values);
      free(em_di_3.alpha_sums);
      free(em_di_3.logbeta_values);
      free(em_di_3.prior_values);
    }
  }
  if(use_priordistribution_4 == DIRICHLET) {
    if(em_di_4.nr_components > 0) {
      free(em_di_4.q_values);
      free(em_di_4.alpha_sums);
      free(em_di_4.logbeta_values);
      free(em_di_4.prior_values);
    }
  }

  return;
}


/* Note: msa_seq_infop->seq and msa_seq_infop->gaps will be allocated here but must be
 * freed by caller */
void get_sequences_msa_prf_multi(FILE *seqfile, FILE *priorfile, struct msa_sequences_multi_s *msa_seq_infop,
				 struct hmm_multi_s *hmmp)
{
  
  int i,j,k,l,m;
  int done, inside_seq, seq_pos, letter_pos, a_index, nr_seqs, cur_pos, cur_seq;
  int msa_length, seq_length, longest_seq, longest_seq_length;
  int seq_index, nr_lead_columns;
  int gaps_per_column, tot_nr_gaps, nr_gaps;
  double occurences_per_column;
  int get_letter_columns, get_query_letters;
  int use_priordistribution,use_priordistribution_2, use_priordistribution_3, use_priordistribution_4; 
  char c;
  char s[MAX_LINE];
  struct letter_s cur_letter;
  int **cur_posp;
  struct emission_dirichlet_s em_di,em_di_2, em_di_3, em_di_4;
  int is_first;
  int use_replacement_letters;
  int is_empty;
  int col_pos;
  char *seq_ptr;
  char *endptr;
  double nr_occurences;
  int label_pos, cur_letter_pos;
  char label, letter;
  int read_priorfile;
  
  
  /* find out length of alignment and allocate memory for probability matrix by
   * reading all sequence rows and remembering the longest */
  done = NO;
  inside_seq = NO;
  msa_length = 0;
  seq_length = 0;
  seq_index = 0;

if(priorfile == NULL) {
    use_priordistribution = NONE;
    use_priordistribution_2 = NONE;
    use_priordistribution_3 = NONE;
    use_priordistribution_4 = NONE;
  }
  else {
    read_priorfile = read_multi_prior_file_multi(&em_di, hmmp, priorfile, 1);
    if(read_priorfile > 0) {
      use_priordistribution = DIRICHLET;
    }
    else if(read_priorfile < 0) {
      printf("Error: Incorrect priorfile format\n");
      exit(0);
    }
    else {
      use_priordistribution = NONE;
    }
    if(hmmp->nr_alphabets > 1) {
      read_priorfile = read_multi_prior_file_multi(&em_di_2, hmmp, priorfile, 2);
      if(read_priorfile > 0) {
	use_priordistribution_2 = DIRICHLET;
      }
      else if(read_priorfile < 0) {
	printf("Error: Incorrect priorfile format\n");
	exit(0);
      }
      else {
	use_priordistribution_2 = NONE;
      }
    }
    if(hmmp->nr_alphabets > 2) {
      read_priorfile = read_multi_prior_file_multi(&em_di_3, hmmp, priorfile, 3);
      if(read_priorfile > 0) {
	use_priordistribution_3 = DIRICHLET;
      }
      else if(read_priorfile < 0) {
	printf("Error: Incorrect priorfile format\n");
	exit(0);
      }
      else {
	use_priordistribution_3 = NONE;
      }
    }
    if(hmmp->nr_alphabets > 3) {
      read_priorfile = read_multi_prior_file_multi(&em_di_4, hmmp, priorfile, 4);
      if(read_priorfile > 0) {
	use_priordistribution_4 = DIRICHLET;
      }
      else if(read_priorfile < 0) {
	printf("Error: Incorrect priorfile format\n");
	exit(0);
      }
      else {
	use_priordistribution_4 = NONE;
      }
    }
  }

  
  
  /* check if file is empty */
  is_empty = YES;
  while(done != YES) {
    c = (char)fgetc(seqfile);
    if((int)c == EOF) {
      break;
    }
    else {
      is_empty = NO;
      break;
    }
  }
  if(is_empty == YES) {
    if(verbose == YES) {
      printf("File is empty\n");
    } 
  }
  else {
  }
  rewind(seqfile);

  while(fgets(s, MAX_LINE, seqfile) != NULL) {
    if(strncmp(s,"COL",3) == 0) {
      msa_length = atoi(&s[4]);
    }
  }
#ifdef DEBUG_MSA_SEQ_PRF
  printf("reached checkpoint 1\n");
  printf("msa_length = %d\n", msa_length);
#endif

  
  msa_seq_infop->lead_columns_start = (int*)malloc_or_die((msa_length+1) * sizeof(int));
  msa_seq_infop->msa_seq_1 = (struct msa_letter_s*)
    malloc_or_die(msa_length * (hmmp->a_size+1) * sizeof(struct msa_letter_s));

#ifdef DEBUG_MSA_SEQ_PRF
  printf("malloced ok\n");
#endif

  /* read alphabet 1, save nr of occurences of each letter in each position,
   * including nr of gaps */
  rewind(seqfile);
  seq_pos = 0;
  nr_lead_columns = 0;
  k = 0;
  l = 0;
  inside_seq = NO;
  while(fgets(s, MAX_LINE, seqfile) != NULL) {
    if(strncmp(s,"END 1",5 ) == 0) {
      break;
    }
    if(strncmp(s,"START 1",7) == 0) {
      inside_seq = YES;
    }
    else if(strncmp(s,"NR of aligned sequences",23) == 0) {
      nr_seqs = strtol(s + 24, NULL, 10);
    }
    else if(strncmp(s,"COL",3) == 0 && inside_seq == YES) {
      /* split into columns depending on a_size, add '-' */
      col_pos = 11;
      seq_ptr = s + col_pos;
      for(m = 0; m < hmmp->a_size + 1; m++) {
	nr_occurences = strtod(seq_ptr, &endptr);
	if(endptr == seq_ptr) {
	  printf("Error reading column: no frequency was read\n");
	  exit(0);
	}
	else {
	  /* add nr of occurences to datastructure */
	  seq_ptr = endptr;
	  (msa_seq_infop->msa_seq_1 + get_mtx_index(seq_pos,m, hmmp->a_size+1))->nr_occurences = nr_occurences;
	}
      }
      /* read space */
      strtod(seq_ptr, &endptr);
      seq_ptr = endptr;

      /* read label */
      label_pos = 0;
      for(label_pos = 0;seq_ptr[label_pos] != '\n';label_pos++) {
	if(seq_ptr[label_pos] == ' ') {
	  
	}
	else {
	  label = seq_ptr[label_pos];
	  seq_ptr = seq_ptr + label_pos + 1;
	  break;
	}
      }
    
      /* read query letter */
      letter_pos = 0;
      cur_letter_pos = 0;
      for(letter_pos = 0;seq_ptr[letter_pos] != '\n';letter_pos++) {
	if(seq_ptr[letter_pos] == ' ') {
	  
	}
	else {
	  letter = seq_ptr[letter_pos];
	  cur_letter.letter[cur_letter_pos] = letter;
	  cur_letter_pos++;
	  seq_ptr = seq_ptr + letter_pos + 1;
	  break;
	}
      }
    
      if(cur_letter_pos >= 5) {
	printf("Maximum of four characters for one letter\n");
	exit(0);
      }
      else {
	cur_letter.letter[cur_letter_pos] = '\0';
      }

      /* store lead seq and query columns + label */
      if(cur_letter.letter[0] != '-') {
	*(msa_seq_infop->lead_columns_start + k) = seq_pos; /* letter column */
	if(label == '-') {
	  label = '.';
	}
	/* only add label to first msa_letter of the first alphabet */
	(msa_seq_infop->msa_seq_1 + (*(msa_seq_infop->lead_columns_start + k)) * (hmmp->a_size+1))->label = label;
	k++;
	nr_lead_columns++;
      }
      if(1 == 1) {
	memcpy((msa_seq_infop->msa_seq_1 + l * (hmmp->a_size + 1))->query_letter, &(cur_letter.letter),
	       sizeof(char) * 5); /* query letter */
	l++;
      }
      seq_pos++;
    }
  }
  
#ifdef DEBUG_MSA_SEQ_PRF
  printf("reached checkpoint 2\n");
  printf("total_nr_gaps = %d\n", tot_nr_gaps);
#endif
  
  
  

  /* read alphabet 2, save nr of occurences of each letter in each position,
   * including nr of gaps */
  if(hmmp->nr_alphabets > 1) {
    msa_seq_infop->msa_seq_2 = (struct msa_letter_s*)
    malloc_or_die(msa_length * (hmmp->a_size_2+1) * sizeof(struct msa_letter_s));
    rewind(seqfile);
    seq_pos = 0;
    k = 0;
    l = 0;
    inside_seq = NO;
    while(fgets(s, MAX_LINE, seqfile) != NULL) {
      if(strncmp(s,"END 2",5 ) == 0) {
	break;
      }
      if(strncmp(s,"START 2",7) == 0) {
	inside_seq = YES;
      }
      else if(strncmp(s,"NR of aligned sequences",23) == 0) {
	nr_seqs = strtol(s + 24, NULL, 10);
      }
      else if(strncmp(s,"COL",3) == 0 && inside_seq == YES) {
	/* split into columns depending on a_size, add '-' */
	col_pos = 11;
	seq_ptr = s + col_pos;
	for(m = 0; m < hmmp->a_size_2 + 1; m++) {
	  nr_occurences = strtod(seq_ptr, &endptr);
	  if(endptr == seq_ptr) {
	    printf("Error reading column: no frequency was read\n");
	    exit(0);
	  }
	  else {
	    /* add nr of occurences to datastructure */
	    seq_ptr = endptr;
	    (msa_seq_infop->msa_seq_2 + get_mtx_index(seq_pos,m, hmmp->a_size_2+1))->nr_occurences = nr_occurences;
	  }
	}
	/* read space */
	strtod(seq_ptr, &endptr);
	seq_ptr = endptr;
	
	/* read label */
	label_pos = 0;
	for(label_pos = 0;seq_ptr[label_pos] != '\n';label_pos++) {
	  if(seq_ptr[label_pos] == ' ') {
	    
	  }
	  else {
	    label = seq_ptr[label_pos];
	    seq_ptr = seq_ptr + label_pos + 1;
	    break;
	  }
	}
	
	/* read query letter */
	letter_pos = 0;
	cur_letter_pos = 0;
	for(letter_pos = 0;seq_ptr[letter_pos] != '\n';letter_pos++) {
	  if(seq_ptr[letter_pos] == ' ') {
	    
	  }
	  else {
	    letter = seq_ptr[letter_pos];
	    cur_letter.letter[cur_letter_pos] = letter;
	    cur_letter_pos++;
	    seq_ptr = seq_ptr + letter_pos + 1;
	    break;
	  }
	}
	
	if(cur_letter_pos >= 5) {
	  printf("Maximum of four characters for one letter\n");
	  exit(0);
	}
	else {
	  cur_letter.letter[cur_letter_pos] = '\0';
	}
	
	/* store lead seq and query columns + label */
	if(cur_letter.letter[0] != '-') {
	  //*(msa_seq_infop->lead_columns_start + k) = seq_pos; /* letter column */
	  //(msa_seq_infop->msa_seq_2 + (*(msa_seq_infop->lead_columns_start + k)) * (hmmp->a_size_2+1))->label = label;
	  k++;
	}
	if(1 == 1) {
	  //memcpy((msa_seq_infop->msa_seq_2 + l * (hmmp->a_size_2 + 1))->query_letter, &(cur_letter.letter),
	  //sizeof(char) * 5); /* query letter */
	  l++;
	}
	seq_pos++;
      }
    }
  
#ifdef DEBUG_MSA_SEQ_PRF
  printf("reached checkpoint 2.1\n");
  printf("total_nr_gaps = %d\n", tot_nr_gaps);
#endif
  }
  
  if(hmmp->nr_alphabets > 2) {
    msa_seq_infop->msa_seq_3 = (struct msa_letter_s*)
      malloc_or_die(msa_length * (hmmp->a_size_3+1) * sizeof(struct msa_letter_s));
    
    /* read alphabet 1, save nr of occurences of each letter in each position,
     * including nr of gaps */
    rewind(seqfile);
    seq_pos = 0;
    k = 0;
    l = 0;
    inside_seq = NO;
    while(fgets(s, MAX_LINE, seqfile) != NULL) {
      if(strncmp(s,"END 3",5 ) == 0) {
	break;
      }
      if(strncmp(s,"START 3",7) == 0) {
	inside_seq = YES;
      }
      else if(strncmp(s,"NR of aligned sequences",23) == 0) {
	nr_seqs = strtol(s + 24, NULL, 10);
      }
      else if(strncmp(s,"COL",3) == 0 && inside_seq == YES) {
	/* split into columns depending on a_size, add '-' */
	col_pos = 11;
	seq_ptr = s + col_pos;
	for(m = 0; m < hmmp->a_size_3 + 1; m++) {
	  nr_occurences = strtod(seq_ptr, &endptr);
	  if(endptr == seq_ptr) {
	    printf("Error reading column: no frequency was read\n");
	    exit(0);
	  }
	  else {
	    /* add nr of occurences to datastructure */
	    seq_ptr = endptr;
	    (msa_seq_infop->msa_seq_3 + get_mtx_index(seq_pos,m, hmmp->a_size_3+1))->nr_occurences = nr_occurences;
	  }
	}
	/* read space */
	strtod(seq_ptr, &endptr);
	seq_ptr = endptr;
	
	/* read label */
	label_pos = 0;
	for(label_pos = 0;seq_ptr[label_pos] != '\n';label_pos++) {
	  if(seq_ptr[label_pos] == ' ') {
	    
	  }
	  else {
	    label = seq_ptr[label_pos];
	    seq_ptr = seq_ptr + label_pos + 1;
	    break;
	  }
	}
	
	/* read query letter */
	letter_pos = 0;
	cur_letter_pos = 0;
	for(letter_pos = 0;seq_ptr[letter_pos] != '\n';letter_pos++) {
	  if(seq_ptr[letter_pos] == ' ') {
	    
	  }
	  else {
	    letter = seq_ptr[letter_pos];
	    cur_letter.letter[cur_letter_pos] = letter;
	    cur_letter_pos++;
	    seq_ptr = seq_ptr + letter_pos + 1;
	    break;
	  }
	}
	
	if(cur_letter_pos >= 5) {
	  printf("Maximum of four characters for one letter\n");
	  exit(0);
	}
	else {
	  cur_letter.letter[cur_letter_pos] = '\0';
	}

	/* store lead seq and query columns + label */
	if(cur_letter.letter[0] != '-') {
	  //*(msa_seq_infop->lead_columns_start + k) = seq_pos; /* letter column */
	  //(msa_seq_infop->msa_seq_3 + (*(msa_seq_infop->lead_columns_start + k)) * (hmmp->a_size_3+1))->label = label;
	  k++;
	}
	if(1 == 1) {
	  //memcpy((msa_seq_infop->msa_seq_3 + l * (hmmp->a_size_3 + 1))->query_letter, &(cur_letter.letter),
	  //	 sizeof(char) * 5); /* query letter */
	  l++;
	}
	seq_pos++;
      }
    }
  
#ifdef DEBUG_MSA_SEQ_PRF
    printf("reached checkpoint 2.2\n");
    printf("total_nr_gaps = %d\n", tot_nr_gaps);
#endif
  }  

  if(hmmp->nr_alphabets > 3) {
    msa_seq_infop->msa_seq_4 = (struct msa_letter_s*)
      malloc_or_die(msa_length * (hmmp->a_size_4+1) * sizeof(struct msa_letter_s));
    
    /* read alphabet 4, save nr of occurences of each letter in each position,
     * including nr of gaps */
    rewind(seqfile);
    seq_pos = 0;
    k = 0;
    l = 0;
    inside_seq = NO;
    while(fgets(s, MAX_LINE, seqfile) != NULL) {
      if(strncmp(s,"END 4",5 ) == 0) {
	break;
      }
      if(strncmp(s,"START 4",7) == 0) {
	inside_seq = YES;
      }
      else if(strncmp(s,"NR of aligned sequences",23) == 0) {
	nr_seqs = strtol(s + 24, NULL, 10);
      }
      else if(strncmp(s,"COL",3) == 0 && inside_seq == YES) {
	/* split into columns depending on a_size, add '-' */
	col_pos = 11;
	seq_ptr = s + col_pos;
	for(m = 0; m < hmmp->a_size_4 + 1; m++) {
	  nr_occurences = strtod(seq_ptr, &endptr);
	  if(endptr == seq_ptr) {
	    printf("Error reading column: no frequency was read\n");
	    exit(0);
	  }
	  else {
	    /* add nr of occurences to datastructure */
	    seq_ptr = endptr;
	    (msa_seq_infop->msa_seq_4 + get_mtx_index(seq_pos,m, hmmp->a_size_4+1))->nr_occurences = nr_occurences;
	  }
	}
	/* read space */
	strtod(seq_ptr, &endptr);
	seq_ptr = endptr;
	
	/* read label */
	label_pos = 0;
	for(label_pos = 0;seq_ptr[label_pos] != '\n';label_pos++) {
	  if(seq_ptr[label_pos] == ' ') {
	    
	  }
	  else {
	    label = seq_ptr[label_pos];
	    seq_ptr = seq_ptr + label_pos + 1;
	    break;
	  }
	}
	
	/* read query letter */
	letter_pos = 0;
	cur_letter_pos = 0;
	for(letter_pos = 0;seq_ptr[letter_pos] != '\n';letter_pos++) {
	  if(seq_ptr[letter_pos] == ' ') {
	    
	  }
	  else {
	    letter = seq_ptr[letter_pos];
	    cur_letter.letter[cur_letter_pos] = letter;
	    cur_letter_pos++;
	    seq_ptr = seq_ptr + letter_pos + 1;
	    break;
	  }
	}
	
	if(cur_letter_pos >= 5) {
	  printf("Maximum of four characters for one letter\n");
	  exit(0);
	}
	else {
	  cur_letter.letter[cur_letter_pos] = '\0';
	}
	
	/* store lead seq and query columns + label */
	if(cur_letter.letter[0] != '-') {
	  //*(msa_seq_infop->lead_columns_start + k) = seq_pos; /* letter column */
	  //(msa_seq_infop->msa_seq_4 + (*(msa_seq_infop->lead_columns_start + k)) * (hmmp->a_size_4+1))->label = label;
	  k++;
	}
	if(1 == 1) {
	  //memcpy((msa_seq_infop->msa_seq_4 + l * (hmmp->a_size_4 + 1))->query_letter, &(cur_letter.letter),
	  //sizeof(char) * 5); /* query letter */
	  l++;
	}
	seq_pos++;
      }
    }
    
#ifdef DEBUG_MSA_SEQ_PRF
    printf("reached checkpoint 2.3\n");
    printf("total_nr_gaps = %d\n", tot_nr_gaps);
#endif
  }
  
  
  /* from everything should be untouched */
  *(msa_seq_infop->lead_columns_start + k) = END;
  msa_seq_infop->lead_columns_end = msa_seq_infop->lead_columns_start + k;
  msa_seq_infop->nr_lead_columns = nr_lead_columns;
  msa_seq_infop->nr_seqs = nr_seqs;
  msa_seq_infop->msa_seq_length = msa_length;
  msa_seq_infop->gap_shares = (double*)malloc_or_die(msa_length * sizeof(double));
  tot_nr_gaps = 0;
  gaps_per_column = 0;
  
  /* go through each position and calculate distribution for all letters in alphabet */
  for(i = 0; i < msa_length; i++) {
    occurences_per_column = 0;
    gaps_per_column = 0;
    for(j = 0; j < hmmp->a_size+1; j++) {
      occurences_per_column += (msa_seq_infop->msa_seq_1 +
				      get_mtx_index(i,j, hmmp->a_size+1))->nr_occurences;
      if(j == hmmp->a_size) {
	tot_nr_gaps += (int)(msa_seq_infop->msa_seq_1 +
			     get_mtx_index(i,j, hmmp->a_size+1))->nr_occurences;
	gaps_per_column = (int)(msa_seq_infop->msa_seq_1 +
				get_mtx_index(i,j, hmmp->a_size+1))->nr_occurences;
      }
    }
    if(use_priordistribution == DIRICHLET) /* calculate share update using dirichlet prior mixture */ {
      update_shares_prior_multi(&em_di, hmmp, msa_seq_infop, i,1);
    }
    /* simple share update just using the standard quotient */ 
    for(j = 0; j < hmmp->a_size+1; j++) {
      (msa_seq_infop->msa_seq_1 + get_mtx_index(i,j, hmmp->a_size+1))->share =
	(msa_seq_infop->msa_seq_1 + get_mtx_index(i,j, hmmp->a_size+1))->nr_occurences / 
	occurences_per_column;
    }    
    
    /* calculate gap_share */
    *(msa_seq_infop->gap_shares + i) = (double)(gaps_per_column) / occurences_per_column;
  }

  /* go through each position and calculate distribution for all letters in alphabet 2 */
  if(hmmp->nr_alphabets > 1) {
    for(i = 0; i < msa_length; i++) {
      occurences_per_column = 0;
      gaps_per_column = 0;
      for(j = 0; j < hmmp->a_size_2+1; j++) {
      occurences_per_column += (msa_seq_infop->msa_seq_2 +
				get_mtx_index(i,j, hmmp->a_size_2+1))->nr_occurences;
      if(j == hmmp->a_size_2) {
	tot_nr_gaps += (int)(msa_seq_infop->msa_seq_2 +
			     get_mtx_index(i,j, hmmp->a_size_2+1))->nr_occurences;
	gaps_per_column = (int)(msa_seq_infop->msa_seq_2 +
				get_mtx_index(i,j, hmmp->a_size_2+1))->nr_occurences;
      }
    }
      if(use_priordistribution_2 == DIRICHLET) /* calculate share update using dirichlet prior mixture */ {
      update_shares_prior_multi(&em_di_2, hmmp, msa_seq_infop, i,2);
      }
    /* simple share update just using the standard quotient */ 
      for(j = 0; j < hmmp->a_size_2+1; j++) {
	(msa_seq_infop->msa_seq_2 + get_mtx_index(i,j, hmmp->a_size_2+1))->share =
	  (msa_seq_infop->msa_seq_2 + get_mtx_index(i,j, hmmp->a_size_2+1))->nr_occurences / 
	  occurences_per_column;
      }
      
      /* calculate gap_share */
      *(msa_seq_infop->gap_shares + i) = (double)(gaps_per_column) / occurences_per_column;
    }
  }

  /* go through each position and calculate distribution for all letters in alphabet 3 */
  if(hmmp->nr_alphabets > 2) {
    for(i = 0; i < msa_length; i++) {
      occurences_per_column = 0;
      gaps_per_column = 0;
      for(j = 0; j < hmmp->a_size_3+1; j++) {
      occurences_per_column += (msa_seq_infop->msa_seq_3 +
				get_mtx_index(i,j, hmmp->a_size_3+1))->nr_occurences;
      if(j == hmmp->a_size_3) {
	tot_nr_gaps += (int)(msa_seq_infop->msa_seq_3 +
			     get_mtx_index(i,j, hmmp->a_size_3+1))->nr_occurences;
	gaps_per_column = (int)(msa_seq_infop->msa_seq_3 +
				get_mtx_index(i,j, hmmp->a_size_3+1))->nr_occurences;
      }
    }
      if(use_priordistribution_3 == DIRICHLET) /* calculate share update using dirichlet prior mixture */ {
      update_shares_prior_multi(&em_di_3, hmmp, msa_seq_infop, i,3);
      }
    /* simple share update just using the standard quotient */ 
      for(j = 0; j < hmmp->a_size_3+1; j++) {
	(msa_seq_infop->msa_seq_3 + get_mtx_index(i,j, hmmp->a_size_3+1))->share =
	  (msa_seq_infop->msa_seq_3 + get_mtx_index(i,j, hmmp->a_size_3+1))->nr_occurences / 
	  occurences_per_column;
      }
      
      /* calculate gap_share */
      *(msa_seq_infop->gap_shares + i) = (double)(gaps_per_column) / occurences_per_column;
    }
  }
  
  /* go through each position and calculate distribution for all letters in alphabet 4 */
  if(hmmp->nr_alphabets > 3) {
    for(i = 0; i < msa_length; i++) {
      occurences_per_column = 0;
      gaps_per_column = 0;
      for(j = 0; j < hmmp->a_size_4+1; j++) {
      occurences_per_column += (msa_seq_infop->msa_seq_4 +
				get_mtx_index(i,j, hmmp->a_size_4+1))->nr_occurences;
      if(j == hmmp->a_size_4) {
	tot_nr_gaps += (int)(msa_seq_infop->msa_seq_4 +
			     get_mtx_index(i,j, hmmp->a_size_4+1))->nr_occurences;
	gaps_per_column = (int)(msa_seq_infop->msa_seq_4 +
				get_mtx_index(i,j, hmmp->a_size_4+1))->nr_occurences;
      }
    }
      if(use_priordistribution_4 == DIRICHLET) /* calculate share update using dirichlet prior mixture */ {
      update_shares_prior_multi(&em_di_4, hmmp, msa_seq_infop, i,4);
      }
    /* simple share update just using the standard quotient */ 
      for(j = 0; j < hmmp->a_size_4+1; j++) {
	(msa_seq_infop->msa_seq_4 + get_mtx_index(i,j, hmmp->a_size_4+1))->share =
	  (msa_seq_infop->msa_seq_4 + get_mtx_index(i,j, hmmp->a_size_4+1))->nr_occurences / 
	  occurences_per_column;
      }
      
      /* calculate gap_share */
      *(msa_seq_infop->gap_shares + i) = (double)(gaps_per_column) / occurences_per_column;
    }
  }
  
#ifdef DEBUG_MSA_SEQ_PRF
  printf("reached checkpoint 3\n");
  printf("total_nr_gaps = %d\n", tot_nr_gaps);
#endif
/* allocate memory for gaps array and set initial pointers for each position*/
  msa_seq_infop->gaps = (int**)malloc_or_die(msa_length * sizeof(int*) +
					     (tot_nr_gaps + msa_length) * sizeof(int));
  nr_gaps = 0;
  for(i = 0; i < msa_length; i++) {
    *(msa_seq_infop->gaps + i) = (int*)(msa_seq_infop->gaps + msa_length) + i;
    **(msa_seq_infop->gaps + i) = END;
  }
#ifdef DEBUG_MSA_SEQ_PRF
  printf("%x\n", msa_seq_infop->gaps);
  printf("%x\n", *(msa_seq_infop->gaps));
  printf("%x\n", *(msa_seq_infop->gaps + 1));
  printf("%x\n", *(msa_seq_infop->gaps + 2));
  printf("%x\n", *(msa_seq_infop->gaps + 3));
  printf("nr_seqs: %d\n", nr_seqs);
#endif

#ifdef DEBUG_MSA_SEQ_PRF
  dump_msa_seqs_multi(msa_seq_infop, hmmp);
  exit(0);
#endif

  /* cleanup and return */
  if(use_priordistribution == DIRICHLET) {
    free(em_di.q_values);
    free(em_di.alpha_sums);
    free(em_di.logbeta_values);
    free(em_di.prior_values);
  }
  if(use_priordistribution == DIRICHLET && hmmp->nr_alphabets > 1) {
    free(em_di_2.q_values);
    free(em_di_2.alpha_sums);
    free(em_di_2.logbeta_values);
    free(em_di_2.prior_values);
  }
  if(use_priordistribution == DIRICHLET && hmmp->nr_alphabets > 2) {
    free(em_di_3.q_values);
    free(em_di_3.alpha_sums);
    free(em_di_3.logbeta_values);
    free(em_di_3.prior_values);
  }
  if(use_priordistribution == DIRICHLET && hmmp->nr_alphabets > 3) {
    free(em_di_4.q_values);
    free(em_di_4.alpha_sums);
    free(em_di_4.logbeta_values);
    free(em_di_4.prior_values);
  }
  return;
}
