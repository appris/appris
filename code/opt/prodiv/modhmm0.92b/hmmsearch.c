#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "structs.h"
#include "funcs.h"
#include "cmdline_hmmsearch.h"


#define NORM_LOG_LIKE 0
#define LOGODDS 1

#define HMM 20
#define SEQ 21

#define LONGEST_SEQ -1 /* Note: this constant is also defined in read_seqs.c */
#define FIRST_SEQ 1

#define LEAD_SEQ 10
#define ALL_SEQS 11

#define MAX_LINE 500

extern int verbose;


static struct hmm_multi_s hmm;
static struct hmm_multi_s retrain_hmm;
static struct msa_sequences_multi_s *msa_seq_infop;
static struct sequences_multi_s seq_info;
static struct replacement_letter_multi_s replacement_letters;
static struct null_model_multi_s null_model;
double *subst_mtxp;
double *subst_mtxp_2;
double *subst_mtxp_3;
double *subst_mtxp_4;
double *aa_freqs;
double *aa_freqs_2;
double *aa_freqs_3;
double *aa_freqs_4;


void get_null_model_multi(FILE *nullfile);
double get_nullmodel_score_multi(struct letter_s *seq, struct letter_s *seq_2, struct letter_s *seq_3,
				 struct letter_s *seq_4, int seq_len, int multi_scoring_method);
double get_nullmodel_score_msa_multi(int seq_len, int prf_mode, int use_prior,
				     int use_gap_shares, int normalize,
				     int scoring_method, int multi_scoring_method, double *aa_freqs,
				     double *aa_freqs_2, double *aa_freqs_3, double *aa_freqs_4);
void get_scores(int run_max_d, int run_forward, int run_backward, int run_viterbi, int run_n_best,
		int seq_format, int prf_mode, int use_gap_shares, int use_prior, int use_labels, int normalize,
		int scoring_method, int multi_scoring_method, int print_labels, int print_post_probs,
		int print_log_like, int print_log_odds, int print_reversi, int print_path, FILE *outfile);
void get_post_prob_matrix(double **ret_post_prob_matrixp, double forw_score, struct forward_s *forw_mtx,
			  struct backward_s *backw_mtx, int seq_len);

int main(int argc, char* argv[])
{
  /* command line options */
  FILE *outfile; /* file to prints results to */
  FILE *hmmnamefile; /* file to read hmm names from */
  FILE *hmmfile; /* file to read hmms from */
  FILE *seqfile; /* file to read sequences from */
  FILE *nullfile; /* file to read null model from */
  FILE *seqnamefile; /* file for reading sequence names */
  FILE *replfile; /* file for reading special replacement letters */
  FILE *priorfile; /* file to read priordistribution from */
  FILE *substmtxfile; /* file to read substitution matrix from */
  FILE *freqfile;
  int seq_format; /* input sequence format */
  int align_alg, label_alg; /* algorithm for calculating match score */
  int query; /* are hmms or seqs the queryside */
  char seq_name[200]; 
  char hmm_name[200];
  char outfilepath_name[200];
  char cur_savefile_name[200];
  char *pure_seq_name_p, pure_seq_name_a[200]; /* sequences name without path name */
  char *pure_hmm_name_p, pure_hmm_name_a[200]; /* sequences name without path name */
  char *last_slash, *last_dot;
  int prf_mode; /* method for doing the msa-hmm search */
  int lead_seq; /* index of the sequence on which the search is performed */
  int use_gap_shares; /* account for/don't account for gaps use when calculating the share of a symbol */
  int use_prior; /* are we working with priors or not */
  int use_labels;
  int normalize;
  int scoring_method;
  int read_subst_mtx;
  int multi_scoring_method;
  int print_log_like;
  int print_log_odds;
  int print_reversi;
  int print_labels;
  int print_post_probs;
  int print_path;
  int run_forward;
  int run_backward;
  int run_viterbi;
  int run_n_best;
  int run_max_d;
  int hmmfiletype;
  int use_lead_columns;
  
  struct gengetopt_args_info args_info;

  /* temporary variables */
  int i,j; /*standard loop indices */
 
  /*init some variables */
  outfile = NULL;
  hmmnamefile = NULL;
  hmmfile = NULL;
  seqfile = NULL;
  nullfile = NULL;
  seqnamefile = NULL;
  replfile = NULL;
  priorfile = NULL;
  substmtxfile = NULL;
  freqfile = NULL;
  seq_format = STANDARD;
  align_alg = FORW_ALG;
  label_alg = ONE_BEST_ALG;
  query = SEQ;
  outfilepath_name[0] = '\0';
  lead_seq = -1;
  prf_mode = ALL_SEQS;
  use_gap_shares = YES;
  use_prior = NO;
  use_labels = NO;
  normalize = NO;
  scoring_method = DOT_PRODUCT;
  read_subst_mtx = NO;
  multi_scoring_method = JOINT_PROB;
  subst_mtxp = NULL;
  subst_mtxp_2 = NULL;
  subst_mtxp_3 = NULL;
  subst_mtxp_4 = NULL;
  aa_freqs = NULL;
  aa_freqs_2 = NULL;
  aa_freqs_3 = NULL;
  aa_freqs_4 = NULL;
  print_log_like = NO;
  print_log_odds = NO;
  print_reversi = NO;
  print_labels = NO;
  print_post_probs = NO;
  print_path = NO;
  run_forward = NO;
  run_backward = NO;
  run_viterbi = NO;
  run_n_best = NO;
  run_max_d = NO;


  /* parse command line */

  if(cmdline_parser(argc, argv, &args_info) != 0) {
    exit(1);
  }

  /* compulsory options */
  if(args_info.hmmnamefile_given) {
    if((hmmnamefile = fopen(args_info.hmmnamefile_arg, "r")) == NULL) {
      perror(args_info.hmmnamefile_arg);
      exit(0);
    }
    else {
      printf("Opened file %s for reading hmm names\n",args_info.hmmnamefile_arg);
    }
  }
  if(args_info.seqnamefile_given) {
    if((seqnamefile = fopen(args_info.seqnamefile_arg, "r")) == NULL) {
      perror(args_info.seqnamefile_arg);
      exit(0);
    }
    else {
      printf("Opened file %s for reading sequence names\n",args_info.seqnamefile_arg);
    }
  }
  if(args_info.outpath_given) {
    strcpy(outfilepath_name, args_info.outpath_arg);
  }
  
  if(args_info.seqformat_given) {
    if(strcmp(args_info.seqformat_arg, "fa") == 0) {
      seq_format = FASTA; 
    }
    else if(strcmp(args_info.seqformat_arg, "s") == 0) {
      seq_format = STANDARD;
    }
    else if(strcmp(args_info.seqformat_arg, "msa") == 0) {
      seq_format = MSA_STANDARD;
    }
    else if(strcmp(args_info.seqformat_arg, "prf") == 0) {
      seq_format = PROFILE;
    }
    else {
      printf("Incorrect sequence format: %s\n", args_info.seqformat_arg);
      exit(0);
    }
  }

  /* non compulsory options */
  if(args_info.smxfile_given) {
    if((substmtxfile = fopen(args_info.smxfile_arg, "r")) == NULL) {
      perror(args_info.smxfile_arg);
      exit(0);
    }
    else {
      read_subst_mtx = YES;
      printf("Opened file %s for reading substitution matrix\n",args_info.smxfile_arg);
    }
  }
  if(args_info.freqfile_given) {
    if((freqfile = fopen(args_info.freqfile_arg, "r")) == NULL) {
      perror(args_info.freqfile_arg);
      exit(0);
    }
    else {
      
      printf("Opened file %s for reading background frequencies\n",args_info.freqfile_arg);
    }
  }
  if(args_info.replfile_given) {
    if((replfile = fopen(args_info.replfile_arg, "r")) == NULL) {
      perror(args_info.replfile_arg);
      exit(0);
    }
    else {
      printf("Opened file %s for reading replacement letters\n",args_info.replfile_arg);
    }
  }
  if(args_info.priorfile_given) {
    if((priorfile = fopen(args_info.priorfile_arg, "r")) == NULL) {
      perror(args_info.priorfile_arg);
      exit(0);
    }
    else {
      printf("Opened file %s for reading priors\n",args_info.priorfile_arg);
    }
  }
  if(args_info.nullfile_given) {
    if((nullfile = fopen(args_info.nullfile_arg, "r")) == NULL) {
      perror(args_info.nullfile_arg);
      exit(0);
    }
    else {
      printf("Opened file %s for reading null model\n",args_info.nullfile_arg);
    }
  }
  
  if(args_info.anchor_given) {
    if(strcmp(args_info.anchor_arg, "hmm") == 0) {
      query = SEQ;
    }
    else if(strcmp(args_info.anchor_arg, "seq") == 0) {
      query = HMM;
    }
    else {
      printf("Eror: Unrecognized anchor option\n");
      exit(0);
    }
  }
 
  if(args_info.labeloutput_given) {
    print_labels = YES;
    print_post_probs = YES;
  }
  
  if(args_info.alignmentoutput_given) {
    print_reversi = YES;
    print_log_odds = YES;
    print_log_like = YES;
  }
  
  if(args_info.viterbi_given) {
    label_alg = VITERBI_ALG;
    align_alg = VITERBI_ALG;
  }

  if(args_info.forward_given) {
    align_alg = FORW_ALG;
  }
  if(args_info.nbest_given) {
    align_alg = ONE_BEST_ALG;
  }
  if(args_info.max_d_given) {
    run_max_d = YES;
  }
  
  if(args_info.path_given) {
    print_path = YES;
  }

  if(args_info.nopostout_given) {
    print_post_probs = NO;
  }
  if(args_info.nolabelout_given) {
    print_labels = NO;
  }
  if(args_info.nollout_given) {
    print_log_like = NO;
  }
  if(args_info.nooddsout_given) {
    print_log_odds = NO;
  }
  if(args_info.norevout_given) {
    print_reversi = NO;
  }

  if(args_info.alignpostout_given) {
    print_post_probs = YES;
  }
  if(args_info.alignlabelout_given) {
    print_labels = YES;
  }
  if(args_info.labelllout_given) {
    print_log_like = YES;
  }
  if(args_info.labeloddsout_given) {
    print_log_odds = YES;
  }
  if(args_info.labelrevout_given) {
    print_reversi = YES;
  }
 
  /* msa scoring options */
  if(args_info.msascoring_given) {
    if(strcmp(args_info.msascoring_arg, "DP") == 0) {
      scoring_method = DOT_PRODUCT;
    }
    else if(strcmp(args_info.msascoring_arg, "DPPI") == 0) {
      scoring_method = DOT_PRODUCT_PICASSO;
    }
    else if(strcmp(args_info.msascoring_arg, "PI") == 0) {
      scoring_method = PICASSO;
    }
    else if(strcmp(args_info.msascoring_arg, "PIS") == 0) {
      scoring_method = PICASSO_SYM;
    }
    else if(strcmp(args_info.msascoring_arg, "GM") == 0) {
      scoring_method = SJOLANDER;
    }
    else if(strcmp(args_info.msascoring_arg, "GMR") == 0) {
      scoring_method = SJOLANDER_REVERSED;
    }
    //else if(strcmp(args_info.msascoring_arg, "SMP") == 0) {
    //  scoring_method = SUBST_MTX_PRODUCT;
    //}
    //else if(strcmp(args_info.msascoring_arg, "SMDP") == 0) {
    //  scoring_method = SUBST_MTX_DOT_PRODUCT;
    //}
    //else if(strcmp(args_info.msascoring_arg, "SMDPP") == 0) {
    //  scoring_method = SUBST_MTX_DOT_PRODUCT_PRIOR;
    //}
    else {
      printf("Incorrect scoring method option: %s\n", args_info.msascoring_arg);
      exit(0);
    }
  }
  if(args_info.usecolumns_given) {
    if(strcmp(args_info.usecolumns_arg, "all") == 0) {
      prf_mode = ALL_SEQS;
    }
    else {
      lead_seq = atoi(args_info.usecolumns_arg);
      prf_mode = LEAD_SEQ;
      if(lead_seq <= 0) {
	printf("Incorrect use-column option: %s\n", args_info.usecolumns_arg);
	exit(0);
      }
    }
  }

  /* flags */
  if(args_info.nolabels_given) {
    /* checked after seqread */
  }
  if(args_info.verbose_given) {
    verbose = YES;
  }
  


  /* check which algorithms to run */
  if(print_post_probs == YES) {
    run_forward = YES;
    run_backward = YES;
  }
  if(print_labels == YES) {
    if(label_alg == VITERBI_ALG) {
      run_viterbi = YES;
    }
    else if(label_alg == ONE_BEST_ALG) {
      run_n_best = YES;
    }
  }
  if(print_log_like == YES || print_log_odds == YES || print_reversi == YES) {
    if(align_alg == VITERBI_ALG) {
      run_viterbi = YES;
    }
    else if(align_alg == FORW_ALG) {
      run_forward = YES;
    }
  }
  
  
  /* read subst mtx */
  if(substmtxfile != NULL) {  
    read_subst_matrix_multi(&subst_mtxp, &subst_mtxp_2, &subst_mtxp_3, &subst_mtxp_4, substmtxfile);
  }
  
  /* get frequency file */
  if(freqfile != NULL) {
    read_frequencies_multi(freqfile, &aa_freqs, &aa_freqs_2, &aa_freqs_3, &aa_freqs_4);
  }
  
  /* get prior file */

  /* get replacement letters */
  if(replfile != NULL) {
    get_replacement_letters_multi(replfile, &replacement_letters);
  }
  else {
    replacement_letters.nr_alphabets = 0;
    replacement_letters.nr_rl_1 = 0;
    replacement_letters.nr_rl_2 = 0;
    replacement_letters.nr_rl_3 = 0;
    replacement_letters.nr_rl_4 = 0;
    replacement_letters.letters_1 = NULL;
    replacement_letters.probs_1 = NULL;
    replacement_letters.letters_2 = NULL;
    replacement_letters.probs_2 = NULL;
    replacement_letters.letters_3 = NULL;
    replacement_letters.probs_3 = NULL;
    replacement_letters.letters_4 = NULL;
    replacement_letters.probs_4 = NULL;
  }

  if((scoring_method == SUBST_MTX_PRODUCT || scoring_method == SUBST_MTX_DOT_PRODUCT ||
      scoring_method == SUBST_MTX_DOT_PRODUCT_PRIOR)
     && read_subst_mtx == NO) {
    printf("Error: No substitution matrix supplied\n");
    exit(0);
  }
  
  /* read null model */
  get_null_model_multi(nullfile);
  

  if(query == SEQ) {
    /* get hmms and score all seqs against them */
    while(fgets(hmm_name, MAX_LINE, hmmnamefile) != NULL) {
      /* get hmm from file */
      hmm_name[strlen(hmm_name)-1] = '\0';
      if((hmmfile = fopen(hmm_name, "r")) == NULL) {
	perror(hmm_name);
	continue;
      }
      hmmfiletype = readhmm_check(hmmfile);
      if(hmmfiletype == SINGLE_HMM) {
	readhmm(hmmfile, &hmm);
      }
      else if(hmmfiletype == MULTI_HMM) {
	readhmm_multialpha(hmmfile, &hmm);
      }
      hmm.subst_mtx = subst_mtxp;
      hmm.subst_mtx_2 = subst_mtxp_2;
      hmm.subst_mtx_3 = subst_mtxp_3;
      hmm.subst_mtx_4 = subst_mtxp_4;
      hmm.replacement_letters = &replacement_letters;
      
      /* open outfile + print init info */
      if((pure_hmm_name_p = strrchr(hmm_name, '/')) != NULL) {
	strcpy(pure_hmm_name_a, pure_hmm_name_p+1);
      }
      else {
	strcpy(pure_hmm_name_a, hmm_name);
      }
      strcpy(cur_savefile_name, outfilepath_name);
      if(strlen(cur_savefile_name) != 0 && cur_savefile_name[strlen(cur_savefile_name) - 1] != '/') {
	strcat(cur_savefile_name, "/");
      }
      strcat(cur_savefile_name, pure_hmm_name_a);
      strcat(cur_savefile_name, ".res");
      
      if((outfile = fopen(cur_savefile_name, "w")) == NULL) {
	perror(cur_savefile_name);
	continue;
      }
      else {
	fprintf(outfile, "# Scores for HMM: '%s'\n\n", pure_hmm_name_a);
      }
      
      /* rewind file with the sequence names */
      rewind(seqnamefile);

      /* loop over the sequences */
      while(fgets(seq_name, MAX_LINE, seqnamefile) != NULL) {
	/* open seqfile */
	seq_name[strlen(seq_name)-1] = '\0';
	if((seqfile = fopen(seq_name, "r")) == NULL) {
	  perror(seq_name);
	  continue;
	}
	
	/* check sequence file for labels + check nolabel flag */
	if(seqfile_has_labels(seqfile) == YES) {
	  use_labels = YES;
	}
	else {
	  use_labels = NO;
	}
	if(args_info.nolabels_given) {
	  use_labels = NO;
	}
	
	/* allocate space for msa_seq_info structs */
	if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
	  msa_seq_infop = (struct msa_sequences_multi_s*)(malloc_or_die(1 * sizeof(struct msa_sequences_multi_s)));
	}
	else if(seq_format == FASTA || seq_format == STANDARD) {
	  seq_info.seqs = malloc_or_die(sizeof(struct sequence_multi_s) * 1);
	  seq_info.nr_alphabets = hmm.nr_alphabets;
	  seq_info.nr_seqs = 1;
	  seq_info.longest_seq = 0;
	  seq_info.shortest_seq = INT_MAX;
	  seq_info.avg_seq_len = 0;
	}
	

	if(seq_format == FASTA) {
	  if(hmm.nr_alphabets > 1) {
	    printf("Error: fasta is a one alphabet format only\n");
	    exit(0);
	  }
	  else {
	    get_sequence_fasta_multi(seqfile, &seq_info, 0);
	    hmm.alphabet_type = DISCRETE;
	    if(use_labels == YES) {
	      get_labels_multi(seqfile, &seq_info, &hmm, 0);
	    }
	  }
	}
	else if(seq_format == STANDARD) {
	  get_sequence_std_multi(seqfile, &seq_info, &hmm, 0);
	  if(use_labels == YES) {
	    get_labels_multi(seqfile, &seq_info, &hmm, 0);
	  }
	}
	else if(seq_format == MSA_STANDARD) {
	  if(prf_mode == LEAD_SEQ) {
	    get_sequences_msa_std_multi(seqfile, priorfile, msa_seq_infop, &hmm,
					lead_seq, &replacement_letters);
	  }
	  else {
	    get_sequences_msa_std_multi(seqfile, priorfile, msa_seq_infop, &hmm, -1, &replacement_letters);
	  }
	  
	  /* get labels */
	  if(use_labels == YES) {
	    if(prf_mode == LEAD_SEQ) {
	      get_msa_labels_multi(seqfile, msa_seq_infop, &hmm);
	    }
	    else if(prf_mode == ALL_SEQS) {
	      get_msa_labels_all_columns_multi(seqfile, msa_seq_infop, &hmm);
	    }
	    else {
	      printf("Internal error: strange prf_mode value\n");
	      exit(0);
	    }
	  }
	}
	else if(seq_format == PROFILE) {
	  get_sequences_msa_prf_multi(seqfile, priorfile,  msa_seq_infop, &hmm);
	}
	if(seq_format == STANDARD || seq_format == FASTA) {
	  seq_info.avg_seq_len = ((int)(seq_info.avg_seq_len / seq_info.nr_seqs));
	  //dump_labeled_seqs_multi(&seq_info);
	}
	else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
	  //dump_msa_seqs_multi(msa_seq_infop, &hmm);
	  //exit(0);
	}
	fclose(seqfile);

	/* print seq info */
	if(verbose == YES) {
	  printf("Running seq %s on %s\n", seq_name, hmm_name);
	}

	
	/* retrain model if run_max_d option is given */
	if(run_max_d == YES) {
	  if(verbose == YES) {
	    printf("retraining ...\n");
	  }
	  if(seq_format == STANDARD || seq_format == FASTA) {
	    copy_hmm_struct(&hmm, &retrain_hmm);
	    baum_welch_dirichlet_multi(&retrain_hmm, seq_info.seqs, 1, NO, use_labels, NO,
				       NO, multi_scoring_method, NO);
	  }
	  else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
	    if(prf_mode == ALL_SEQS) {
	      use_lead_columns = NO;
	    }
	    else if(prf_mode == LEAD_SEQ) {
	      use_lead_columns = YES;
	    }
	    copy_hmm_struct(&hmm, &retrain_hmm);
	    msa_baum_welch_dirichlet_multi(&retrain_hmm, msa_seq_infop, 1, NO, use_gap_shares, use_lead_columns, use_labels,
					   NO, NO, normalize, scoring_method,
					   NO, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4,
					   NO);
	  }
	  if(verbose == YES) {
	    printf("retraining done\n");
	  }
	}
	
	if((pure_seq_name_p = strrchr(seq_name, '/')) != NULL) {
	  strcpy(pure_seq_name_a, pure_seq_name_p+1);
	}
	else {
	  strcpy(pure_seq_name_a, seq_name);
	}
	fprintf(outfile, "%s\n", pure_seq_name_a);
	

	/* get score info for this seq and print it to outfile */
	if(verbose == YES) {
	  printf("scoring ...\n");
	}
	get_scores(run_max_d, run_forward, run_backward, run_viterbi, run_n_best, seq_format, prf_mode,
		   use_gap_shares, use_prior, use_labels, normalize, scoring_method, multi_scoring_method,
		   print_labels, print_post_probs, print_log_like, print_log_odds, print_reversi, print_path,
		   outfile);
	if(verbose == YES) {
	  printf("scoring done\n");
	}
	
	/* deallocate retrain_hmm */
	if(run_max_d == YES) {
	  hmm_garbage_collection_multi_no_dirichlet(NULL, &retrain_hmm);
	}

	/* deallocate seqinfo */
	if(seq_format == STANDARD || seq_format == FASTA) {
	  free(((seq_info.seqs))->seq_1);
	  if(hmm.nr_alphabets > 1) {
	    free((seq_info.seqs)->seq_2);
	  }
	  if(hmm.nr_alphabets > 2) {
	    free((seq_info.seqs)->seq_3);
	  }
	  if(hmm.nr_alphabets > 3) {
	    free((seq_info.seqs)->seq_4);
	  }
	  free(seq_info.seqs);
	}
	if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
	  free((*(msa_seq_infop)).msa_seq_1);
	  if(hmm.nr_alphabets > 1) {
	    free((*(msa_seq_infop)).msa_seq_2);
	  }
	  if(hmm.nr_alphabets > 2) {
	    free((*(msa_seq_infop)).msa_seq_3);
	  }
	  if(hmm.nr_alphabets > 3) {
	    free((*(msa_seq_infop)).msa_seq_4);
	  }
	  free((*(msa_seq_infop)).gaps);
	}
	free(msa_seq_infop);
      }
      
      /* deallocate hmm_info */
      hmm_garbage_collection_multi(hmmfile, &hmm);
      fclose(outfile);
    }
  }

  
  else if(query == HMM) {
    /* for each sequence
      * score it on all hmms and print results
      * done as if query = seq, but results are saved differently */
    
    /* print init info in outfiles */
    rewind(seqnamefile);
    while(fgets(seq_name, MAX_LINE, seqnamefile) != NULL) {
      /* open outfile + print init info */
      seq_name[strlen(seq_name)-1] = '\0';
      if((pure_seq_name_p = strrchr(seq_name, '/')) != NULL) {
	strcpy(pure_seq_name_a, pure_seq_name_p+1);
      }
      else {
	strcpy(pure_seq_name_a, seq_name);
      }
      strcpy(cur_savefile_name, outfilepath_name);
      if(strlen(cur_savefile_name) != 0 && cur_savefile_name[strlen(cur_savefile_name) - 1] != '/') {
	strcat(cur_savefile_name, "/");
      }
      strcat(cur_savefile_name, pure_seq_name_a);
      strcat(cur_savefile_name, ".res");
      
      if((outfile = fopen(cur_savefile_name, "w")) == NULL) {
	perror(cur_savefile_name);
	continue;
      }
      else {
	fprintf(outfile, "# Scores for sequence: '%s'\n\n", pure_seq_name_a);
	fclose(outfile);
      }
    }
    while(fgets(hmm_name, MAX_LINE, hmmnamefile) != NULL) {
      /* get hmm from file */
      hmm_name[strlen(hmm_name)-1] = '\0';
      if((hmmfile = fopen(hmm_name, "r")) == NULL) {
	perror(hmm_name);
	continue;
      }
      hmmfiletype = readhmm_check(hmmfile);
      if(hmmfiletype == SINGLE_HMM) {
	readhmm(hmmfile, &hmm);
      }
      else if(hmmfiletype == MULTI_HMM) {
	readhmm_multialpha(hmmfile, &hmm);
      }
      hmm.subst_mtx = subst_mtxp;
      hmm.subst_mtx_2 = subst_mtxp_2;
      hmm.subst_mtx_3 = subst_mtxp_3;
      hmm.subst_mtx_4 = subst_mtxp_4;
      hmm.replacement_letters = &replacement_letters;
      
      if((pure_hmm_name_p = strrchr(hmm_name, '/')) != NULL) {
	strcpy(pure_hmm_name_a, pure_hmm_name_p+1);
      }
      else {
	strcpy(pure_hmm_name_a, hmm_name);
      }

      /* rewind file with the sequence names */
      rewind(seqnamefile);

      /* loop over the sequences */
      while(fgets(seq_name, MAX_LINE, seqnamefile) != NULL) {
	/* open seqfile */
	seq_name[strlen(seq_name)-1] = '\0';
	if((seqfile = fopen(seq_name, "r")) == NULL) {
	  perror(seq_name);
	  continue;
	}
	
	/* open outfile + print init info */
	if((pure_seq_name_p = strrchr(seq_name, '/')) != NULL) {
	  strcpy(pure_seq_name_a, pure_seq_name_p+1);
	}
	else {
	  strcpy(pure_seq_name_a, seq_name);
	}
	strcpy(cur_savefile_name, outfilepath_name);
	if(strlen(cur_savefile_name) != 0 && cur_savefile_name[strlen(cur_savefile_name) - 1] != '/') {
	  strcat(cur_savefile_name, "/");
	}
	strcat(cur_savefile_name, pure_seq_name_a);
	strcat(cur_savefile_name, ".res");
      
	if((outfile = fopen(cur_savefile_name, "a")) == NULL) {
	  perror(cur_savefile_name);
	  continue;
	}
	else {
	}

	/* check sequence file for labels + check nolabel flag */
	if(seqfile_has_labels(seqfile) == YES) {
	  use_labels = YES;
	}
	else {
	  use_labels = NO;
	}
	if(args_info.nolabels_given) {
	  use_labels = NO;
	}
	
	/* allocate space for msa_seq_info structs */
	if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
	  msa_seq_infop = (struct msa_sequences_multi_s*)(malloc_or_die(1 * sizeof(struct msa_sequences_multi_s)));
	}
	else if(seq_format == FASTA || seq_format == STANDARD) {
	  seq_info.seqs = malloc_or_die(sizeof(struct sequence_multi_s) * 1);
	  seq_info.nr_alphabets = hmm.nr_alphabets;
	  seq_info.nr_seqs = 1;
	  seq_info.longest_seq = 0;
	  seq_info.shortest_seq = INT_MAX;
	  seq_info.avg_seq_len = 0;
	}
	

	if(seq_format == FASTA) {
	  if(hmm.nr_alphabets > 1) {
	    printf("Error: fasta is a one alphabet format only\n");
	    exit(0);
	  }
	  else {
	    get_sequence_fasta_multi(seqfile, &seq_info, 0);
	    hmm.alphabet_type = DISCRETE;
	    if(use_labels == YES) {
	      get_labels_multi(seqfile, &seq_info, &hmm, 0);
	    }
	  }
	}
	else if(seq_format == STANDARD) {
	  get_sequence_std_multi(seqfile, &seq_info, &hmm, 0);
	  if(use_labels == YES) {
	    get_labels_multi(seqfile, &seq_info, &hmm, 0);
	  }
	}
	else if(seq_format == MSA_STANDARD) {
	  if(prf_mode == LEAD_SEQ) {
	    get_sequences_msa_std_multi(seqfile, priorfile, msa_seq_infop, &hmm,
					lead_seq, &replacement_letters);
	  }
	  else {
	    get_sequences_msa_std_multi(seqfile, priorfile, msa_seq_infop, &hmm, -1, &replacement_letters);
	  }
	  
	  /* get labels */
	  if(use_labels == YES) {
	    if(prf_mode == LEAD_SEQ) {
	      get_msa_labels_multi(seqfile, msa_seq_infop, &hmm);
	    }
	    else if(prf_mode == ALL_SEQS) {
	      get_msa_labels_all_columns_multi(seqfile, msa_seq_infop, &hmm);
	    }
	    else {
	      printf("Internal error: strange prf_mode value\n");
	      exit(0);
	    }
	  }
	}
	else if(seq_format == PROFILE) {
	  get_sequences_msa_prf_multi(seqfile, priorfile,  msa_seq_infop, &hmm);
	}
	if(seq_format == STANDARD || seq_format == FASTA) {
	  seq_info.avg_seq_len = ((int)(seq_info.avg_seq_len / seq_info.nr_seqs));
	  //dump_labeled_seqs_multi(&seq_info);
	}
	else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
	  
	}
	fclose(seqfile);


	fprintf(outfile, "%s:\n", pure_hmm_name_a);
	if(verbose == YES) {
	  printf("Running seq %s on %s\n", seq_name, hmm_name);
	}
	
	/* retrain model if run_max_d option is given */
	if(run_max_d == YES) {
	  if(verbose == YES) {
	    printf("retraining ... \n");
	  }
	  if(seq_format == STANDARD || seq_format == FASTA) {
	    copy_hmm_struct(&hmm, &retrain_hmm);
	    baum_welch_dirichlet_multi(&retrain_hmm, seq_info.seqs, 1, NO, use_labels, NO,
				       NO, multi_scoring_method, NO);
	  }
	  else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
	    if(prf_mode == ALL_SEQS) {
	      use_lead_columns = NO;
	    }
	    else if(prf_mode == LEAD_SEQ) {
	      use_lead_columns = YES;
	    }
	    copy_hmm_struct(&hmm, &retrain_hmm);
	    msa_baum_welch_dirichlet_multi(&retrain_hmm, msa_seq_infop, 1, NO, use_gap_shares, use_lead_columns, use_labels,
					   NO, NO, normalize, scoring_method,
					   NO, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4,
					   NO);
	  }
	  if(verbose == YES) {
	    printf("retraining done\n");
	  }
	}
	  
	
	/* get score info for this seq and print it to outfile */
	if(verbose == YES) {
	  printf("scoring ...\n");
	}
	get_scores(run_max_d, run_forward, run_backward, run_viterbi, run_n_best, seq_format, prf_mode,
		   use_gap_shares, use_prior, use_labels, normalize, scoring_method, multi_scoring_method,
		   print_labels, print_post_probs, print_log_like, print_log_odds, print_reversi, print_path,
		   outfile);
	if(verbose == YES) {
	  printf("scoring done\n");
	}
	
	/* deallocate retrain_hmm */
	if(run_max_d == YES) {
	  hmm_garbage_collection_multi_no_dirichlet(NULL, &retrain_hmm);
	}
	
	/* deallocate seqinfo */
	if(seq_format == STANDARD || seq_format == FASTA) {
	  free(((seq_info.seqs))->seq_1);
	  if(hmm.nr_alphabets > 1) {
	    free((seq_info.seqs)->seq_2);
	  }
	  if(hmm.nr_alphabets > 2) {
	    free((seq_info.seqs)->seq_3);
	  }
	  if(hmm.nr_alphabets > 3) {
	    free((seq_info.seqs)->seq_4);
	  }
	  free(seq_info.seqs);
	}
	if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
	  free((*(msa_seq_infop)).msa_seq_1);
	  if(hmm.nr_alphabets > 1) {
	    free((*(msa_seq_infop)).msa_seq_2);
	  }
	  if(hmm.nr_alphabets > 2) {
	    free((*(msa_seq_infop)).msa_seq_3);
	  }
	  if(hmm.nr_alphabets > 3) {
	    free((*(msa_seq_infop)).msa_seq_4);
	  }
	  free((*(msa_seq_infop)).gaps);
	}
	free(msa_seq_infop);
	fclose(outfile);
      }
      
      /* deallocate hmm_info */
      hmm_garbage_collection_multi(hmmfile, &hmm);
    }
  }


  /* deallocate null model */
  if(nullfile != NULL) {
    if(null_model.nr_alphabets >= 1) {
      free(null_model.emissions);
    }
    if(null_model.nr_alphabets >= 2) {
      free(null_model.emissions_2);
    }
    if(null_model.nr_alphabets >= 3) {
      free(null_model.emissions_3);
    }	
    if(null_model.nr_alphabets >= 4) {
      free(null_model.emissions_4);
    }
    fclose(nullfile);
  }
  
  /* deallocate subst mtx */
  if(substmtxfile != NULL) {
    free(subst_mtxp);
    if(hmm.nr_alphabets > 1) {
      free(subst_mtxp_2);
    }
    if(hmm.nr_alphabets > 2) {
      free(subst_mtxp_3);
    }
    if(hmm.nr_alphabets > 3) {
      free(subst_mtxp_4);
    }
    fclose(substmtxfile);
  }

  /* deallocate freqs */
  if(freqfile != NULL) {
    free(aa_freqs);
    if(hmm.nr_alphabets > 1) {
      free(aa_freqs_2);
    }
    if(hmm.nr_alphabets > 2) {
      free(aa_freqs_3);
    }
    if(hmm.nr_alphabets > 3) {
      free(aa_freqs_4);
    }
    fclose(freqfile);
  }

  /* deallocate replacement letters */
  if(replfile != NULL) {
    if(replacement_letters.nr_rl_1 > 0) {
      free(replacement_letters.letters_1);
      free(replacement_letters.probs_1);
    }
    if(replacement_letters.nr_rl_2 > 0) {
      free(replacement_letters.letters_2);
      free(replacement_letters.probs_2);
    }
    if(replacement_letters.nr_rl_3 > 0) {
      free(replacement_letters.letters_3);
      free(replacement_letters.probs_3);
    }
    if(replacement_letters.nr_rl_4 > 0) {
      free(replacement_letters.letters_4);
      free(replacement_letters.probs_4);
    }
    fclose(replfile);
  }
}





/*************************************************************************************************/
void get_null_model_multi(FILE *nullfile)
{
  int i,j;
  char s[MAX_LINE];
  
  if(nullfile != NULL) {
    /* read nullfile */
    while(1) {
      if(fgets(s, MAX_LINE, nullfile) == NULL) {
	printf("Could not read null model file\n");
	exit(0);
      }
      if(s[0] == '#' || s[0] == '\n') {
	continue;
      }
      else {
	null_model.nr_alphabets = atoi(s);
	break;
      }
    }

    for(i = 0; i < null_model.nr_alphabets; i++) {
      while(1) {
	if(fgets(s, MAX_LINE, nullfile) == NULL) {
	  printf("Could not read null model file\n");
	  exit(0);
	}
	if(s[0] == '#' || s[0] == '\n') {
	  continue;
	}
	else {
	  switch(i) {
	  case 0:
	    null_model.a_size = atoi(s);
	    null_model.emissions = (double*)malloc_or_die(null_model.a_size * sizeof(double));
	    break;
	  case 1:
	    null_model.a_size_2 = atoi(s);
	    null_model.emissions_2 = (double*)malloc_or_die(null_model.a_size_2 * sizeof(double));
	    break;
	  case 2:
	    null_model.a_size_3 = atoi(s);
	    null_model.emissions_3 = (double*)malloc_or_die(null_model.a_size_3 * sizeof(double));
	    break;
	  case 3:
	    null_model.a_size_4 = atoi(s);
	    null_model.emissions_4 = (double*)malloc_or_die(null_model.a_size_4 * sizeof(double));
	    break;
	  }
	  break;
	}
      }
      j = 0;
      switch(i) {
      case 0:
	while(j < null_model.a_size) {
	  if (fgets(s, MAX_LINE, nullfile) != NULL) {
	    if(s[0] != '#' && s[0] != '\n') {
	      *(null_model.emissions + j) = atof(s);
	      j++;
	    }
	  }
	  else {
	    printf("Could not read null model file\n");
	    exit(0);
	  }
	}
	break;
      case 1:
	while(j < null_model.a_size_2) {
	  if (fgets(s, MAX_LINE, nullfile) != NULL) {
	    if(s[0] != '#' && s[0] != '\n') {
	      *(null_model.emissions_2 + j) = atof(s);
	      j++;
	    }
	  }
	  else {
	    printf("Could not read null model file\n");
	    exit(0);
	  }
	}
	break;	
      case 2:
	while(j < null_model.a_size_3) {
	  if (fgets(s, MAX_LINE, nullfile) != NULL) {
	    if(s[0] != '#' && s[0] != '\n') {
	      *(null_model.emissions_3 + j) = atof(s);
	      j++;
	    }
	  }
	  else {
	    printf("Could not read null model file\n");
	    exit(0);
	  }
	}
	break;
      case 3:
	while(j < null_model.a_size_4) {
	  if (fgets(s, MAX_LINE, nullfile) != NULL) {
	    if(s[0] != '#' && s[0] != '\n') {
	      *(null_model.emissions_4 + j) = atof(s);
	      j++;
	    }
	  }
	  else {
	    printf("Could not read null model file\n");
	    exit(0);
	  }
	}
	break;
      }  
    }
    while(1) {
      if(fgets(s, MAX_LINE, nullfile) != NULL) {
	if(s[0] != '#' && s[0] != '\n') {
	  null_model.trans_prob = atof(s);
	  break;
	}
      }
      else {
	printf("Could not read null model file\n");
	exit(0);
      }
    }
  }
  else {
    null_model.a_size = -1;
    null_model.nr_alphabets = -1;
  }
}



/*********************score methods******************************************/
void get_scores(int run_max_d, int run_forward, int run_backward, int run_viterbi, int run_n_best,
		int seq_format, int prf_mode,
		int use_gap_shares, int use_prior, int use_labels, int normalize, int scoring_method,
		int multi_scoring_method, int print_labels, int print_post_probs, int print_log_like,
		int print_log_odds, int print_reversi, int print_path, FILE *outfile)
{
  struct forward_s *forw_mtx;
  struct viterbi_s *viterbi_mtx;
  struct one_best_s *one_best_mtx;
  struct backward_s *backw_mtx;
  struct forward_s *rev_forw_mtx;
  struct viterbi_s *rev_viterbi_mtx;
  double *forw_scale, *rev_forw_scale, *scaling_factors;
  char *labels;
  int use_lead_columns;
  struct letter_s *reverse_seq, *reverse_seq_2, *reverse_seq_3, *reverse_seq_4;
  struct msa_sequences_multi_s reverse_msa_seq_info;
  double nullmodel_score;
  double norm_log_like, logodds, reversi;
  int *path, *path_incrementable;

  int seq_len;
  double viterbi_score, forward_score, backward_score, one_best_score, rev_viterbi_score, rev_forward_score, raw_forward_score;
  double *post_prob_matrix, *post_prob_label_matrix;
  int label_nr;

  int a,b;
  int i,j,k;
  
  struct hmm_multi_s *hmmp;

  if(run_max_d == YES) {
    hmmp = &retrain_hmm;
  }
  else {
    hmmp = &hmm;
  }

  if(prf_mode == ALL_SEQS) {
    use_lead_columns = NO;
  }
  else if(prf_mode == LEAD_SEQ) {
    use_lead_columns = YES;
  }
  else {
    printf("Internal error: strange profile mode\n");
    exit(0);
  }
  
  /* get seq_length */
  if(seq_format == STANDARD || seq_format == FASTA) {
    seq_len = get_seq_length(seq_info.seqs->seq_1);
  }
  else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
    if(use_lead_columns == YES) {
      seq_len = msa_seq_infop->nr_lead_columns;
    }
    else {
      seq_len = msa_seq_infop->msa_seq_length;
    }
  }
  fprintf(outfile, "Seq length: %d\n", seq_len);

  
  /* if needed run forward */
  if(run_forward == YES) {
    if(seq_format == STANDARD || seq_format == FASTA) {
      forward_multi(hmmp, seq_info.seqs->seq_1, seq_info.seqs->seq_2, seq_info.seqs->seq_3, seq_info.seqs->seq_4,
		    &forw_mtx, &forw_scale, use_labels, multi_scoring_method);
      raw_forward_score = (forw_mtx + get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob;
      forward_score = log10((forw_mtx + get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob);
      for (j = seq_len; j > 0; j--) {
	forward_score = forward_score + log10(*(forw_scale + j));
      }
    }
    else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
      msa_forward_multi(hmmp, msa_seq_infop, use_lead_columns, use_gap_shares, use_prior, &forw_mtx, &forw_scale,
			use_labels, normalize,scoring_method, multi_scoring_method, aa_freqs, aa_freqs_2,
			aa_freqs_3, aa_freqs_4);
      raw_forward_score = (forw_mtx + get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob;
      forward_score = log10((forw_mtx + get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob);
      for (j = seq_len; j > 0; j--) {
	forward_score = forward_score + log10(*(forw_scale + j));
      }
    }
  }
  
  /* if needed run backward */
  if(run_backward == YES) {
    if(seq_format == STANDARD || seq_format == FASTA) {
      backward_multi(hmmp, seq_info.seqs->seq_1, seq_info.seqs->seq_2, seq_info.seqs->seq_3, seq_info.seqs->seq_4,
		     &backw_mtx, forw_scale, use_labels, multi_scoring_method);
      backward_score = log10((backw_mtx + get_mtx_index(0, 0, hmmp->nr_v))->prob);
      for (j = seq_len; j > 0; j--) {
	backward_score = backward_score + log10(*(forw_scale + j));
      }
    }
    else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
      msa_backward_multi(hmmp, msa_seq_infop, use_lead_columns, use_gap_shares, &backw_mtx, forw_scale, use_labels, normalize,
			 scoring_method, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
      backward_score = log10((backw_mtx + get_mtx_index(0, 0, hmmp->nr_v))->prob);
      for (j = seq_len; j > 0; j--) {
	backward_score = backward_score + log10(*(forw_scale + j));
      }
    }
  }
  
  /* if needed run viterbi */
  if(run_viterbi == YES || print_path == YES) {
    if(seq_format == STANDARD || seq_format == FASTA) {
      viterbi_multi(hmmp, seq_info.seqs->seq_1, seq_info.seqs->seq_2, seq_info.seqs->seq_3, seq_info.seqs->seq_4,
		    &viterbi_mtx, use_labels, multi_scoring_method);
      viterbi_score = (viterbi_mtx + get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob;
      if(viterbi_score == DEFAULT) {
	viterbi_score = log10(0.0);
      }
    }
    else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
      msa_viterbi_multi(hmmp, msa_seq_infop, use_lead_columns, use_gap_shares, use_prior, &viterbi_mtx, use_labels, normalize,
			scoring_method, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
      viterbi_score = (viterbi_mtx + get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob;
      if(viterbi_score == DEFAULT) {
	viterbi_score = log10(0.0);
      }
    }
  }
  
  /* if needed run n-best */
  if(run_n_best == YES) {
    labels = (char*)malloc_or_die((seq_len + 1) * sizeof(char));
    if(seq_format == STANDARD || seq_format == FASTA) {
      one_best_multi(hmmp, seq_info.seqs->seq_1, seq_info.seqs->seq_2, seq_info.seqs->seq_3, seq_info.seqs->seq_4,
		     &one_best_mtx, &scaling_factors, use_labels, labels, multi_scoring_method);
      one_best_score = log10((one_best_mtx + get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob);
      for (j = seq_len; j > 0; j--) {
	one_best_score = one_best_score + log10(*(scaling_factors + j));
      }
    }
    else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
      msa_one_best_multi(hmmp, msa_seq_infop, use_lead_columns, use_gap_shares, use_prior,
			 &one_best_mtx, &scaling_factors, use_labels, labels,
			 normalize, scoring_method, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
      one_best_score = log10((one_best_mtx + get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob);
      for (j = seq_len; j > 0; j--) {
	one_best_score = one_best_score + log10(*(scaling_factors + j));
      }
    }
  }


  /* get log likelihod score + print + deallocate */
  if(print_log_like == YES) {
    if(run_viterbi == YES) {
      norm_log_like = 0.0 - (viterbi_score/(double)(seq_len));
    }
    else if(run_forward == YES) {
      norm_log_like = 0.0 - (forward_score/(double)(seq_len));
    }
    else {
      
    }
    /* print to outfile */
    fprintf(outfile, "Normalized log likelihood: %.4f\n", norm_log_like);
    
  }
  
  /* get log odds score  + print + deallocate */
  if(print_log_odds == YES) {
    if(seq_format == STANDARD || seq_format == FASTA) {
      nullmodel_score = get_nullmodel_score_multi(seq_info.seqs->seq_1, seq_info.seqs->seq_2, seq_info.seqs->seq_3,
						  seq_info.seqs->seq_4, seq_len, multi_scoring_method);
    }
    else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
      nullmodel_score = get_nullmodel_score_msa_multi(seq_len, use_labels, use_prior, use_gap_shares, normalize, scoring_method,
						      multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
    }
    if(run_viterbi == YES) {
      logodds = viterbi_score - nullmodel_score;
    }
    if(run_forward == YES) {
      logodds = forward_score - nullmodel_score;
    }
    /* print to outfile */
    fprintf(outfile, "Logodds: %.4f\n", logodds);
  }
  /* get reversi score + print + deallocate */
  if(print_reversi == YES) {
    if(seq_format == STANDARD || seq_format == FASTA) {
      get_reverse_seq_multi(seq_info.seqs, &reverse_seq, &reverse_seq_2, &reverse_seq_3, &reverse_seq_4, hmmp, seq_len);
      if(run_forward == YES) {
	forward_multi(hmmp, reverse_seq, reverse_seq_2, reverse_seq_3, reverse_seq_4,
		      &rev_forw_mtx, &rev_forw_scale, use_labels, multi_scoring_method);
	rev_forward_score = log10((rev_forw_mtx + get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob);
	for (j = seq_len; j > 0; j--) {
	  rev_forward_score = rev_forward_score + log10(*(rev_forw_scale + j));
	}
      }
      else if(run_viterbi == YES) {
	viterbi_multi(hmmp, reverse_seq, reverse_seq_2, reverse_seq_3, reverse_seq_4,
		      &rev_viterbi_mtx, use_labels, multi_scoring_method);
	rev_viterbi_score = (rev_viterbi_mtx + get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob;
	if(rev_viterbi_score == DEFAULT) {
	  rev_viterbi_score = log10(0.0);
	}
      }
    }
    else if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
      get_reverse_msa_seq_multi(msa_seq_infop, &reverse_msa_seq_info, hmmp);
      if(run_forward == YES) {
	msa_forward_multi(hmmp, &reverse_msa_seq_info, use_lead_columns, use_gap_shares, use_prior, &rev_forw_mtx,
			  &rev_forw_scale, use_labels, normalize, scoring_method, multi_scoring_method, aa_freqs,
			  aa_freqs_2, aa_freqs_3, aa_freqs_4);
	rev_forward_score = log10((rev_forw_mtx + get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob);
	for (j = seq_len; j > 0; j--) {
	  rev_forward_score = rev_forward_score + log10(*(rev_forw_scale + j));
	}
      }
      else if(run_viterbi == YES) {
	msa_viterbi_multi(hmmp, &reverse_msa_seq_info, use_lead_columns, use_gap_shares, use_prior, &rev_viterbi_mtx,
			  use_labels, normalize,
			  scoring_method, multi_scoring_method, aa_freqs, aa_freqs_2, aa_freqs_3, aa_freqs_4);
	rev_viterbi_score = (rev_viterbi_mtx + get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v))->prob;
	if(rev_viterbi_score == DEFAULT) {
	  rev_viterbi_score = log10(0.0);
	}
      }
    }
    if(run_viterbi == YES) {
      reversi = viterbi_score - rev_viterbi_score;
    }
    if(run_forward == YES) {
      reversi = forward_score - rev_forward_score;
    }
    /* print to outfile */
    fprintf(outfile, "Reversi: %.4f\n", reversi);
    
    /* deallocate */
    if(seq_format == STANDARD || seq_format == FASTA) {
      free(reverse_seq);
      if(hmmp->nr_alphabets > 1) {
	free(reverse_seq_2);
      }
      if(hmmp->nr_alphabets > 2) {
	free(reverse_seq_3);
      }
      if(hmmp->nr_alphabets > 3) {
	free(reverse_seq_4);
      }	  
    }
    if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
      free(reverse_msa_seq_info.msa_seq_1);
	if(hmmp->nr_alphabets > 1) {
	  free(reverse_msa_seq_info.msa_seq_2);
	}
	if(hmmp->nr_alphabets > 2) {
	  free(reverse_msa_seq_info.msa_seq_3);
	}
	if(hmmp->nr_alphabets > 3) {
	  free(reverse_msa_seq_info.msa_seq_4);
	}
	free(reverse_msa_seq_info.gaps);
	free(reverse_msa_seq_info.gap_shares);
	free(reverse_msa_seq_info.lead_columns_start);
    }
    if(run_forward == YES) {
      free(rev_forw_mtx);
      free(rev_forw_scale);
    }
    if(run_viterbi == YES) {
      free(rev_viterbi_mtx);
    }
  }
  
  
  /* get labeling + print + deallocate */
  if(print_labels == YES) {
    if(seq_format == STANDARD || seq_format == FASTA) {
      if(run_viterbi == YES) {
	a = 0;
	labels = (char*)malloc_or_die((seq_len + 1) * sizeof(char));
	get_viterbi_label_path_multi(viterbi_mtx + get_mtx_index(seq_len+1,hmmp->nr_v-1,hmmp->nr_v),
				     hmmp, viterbi_mtx, seq_len + 1, hmmp->nr_v,
				     labels, &a);
	labels[seq_len] = '\0';
      }
    }
    if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
      if(run_viterbi == YES) {
	a = 0;
	labels = (char*)malloc_or_die((seq_len + 1) * sizeof(char));
	get_viterbi_label_path_multi(viterbi_mtx + get_mtx_index(seq_len+1, hmmp->nr_v-1, hmmp->nr_v),
				     hmmp, viterbi_mtx, seq_len + 1, hmmp->nr_v, labels, &a);
	labels[seq_len] = '\0';
      }
    }
    /* print to outfile */
    fprintf(outfile, "Labeling:\n");
    j = 0;
    while(1) {
      if(labels[j] == '\0') {
	fprintf(outfile, "\n");
	break;
      }
      else {
	fprintf(outfile, "%c", labels[j]);
	j++;
      }
      if(j%60 == 0) {
	fprintf(outfile, "\n");
      }
    }
    
    /* deallocate */
    free(labels);
  }

  /* get posterior probs + print + deallocate */
  if(print_post_probs == YES) {
    get_post_prob_matrix(&post_prob_matrix, raw_forward_score, forw_mtx,
			 backw_mtx, seq_len);
    
    //dump_post_prob_matrix(seq_len + 2, hmmp->nr_v, post_prob_matrix);
    post_prob_label_matrix = (double*)(malloc_or_die((seq_len+2) * hmmp->nr_labels * sizeof(double)));
    for(i = 0; i < seq_len + 2; i++) {
      for(j = 1; j < hmmp->nr_v - 1; j++) {
	int label_nr = -1;
	for(k = 0; k < hmmp->nr_labels; k++) {
	  if(hmmp->vertex_labels[j] == hmmp->labels[k]) {
	    label_nr = k;
	    break;
	  }
	}
	if(label_nr < 0) {
	  printf("Internal error: strange label\n");
	  exit(0);
	}
	*(post_prob_label_matrix + get_mtx_index(i,label_nr,hmmp->nr_labels)) +=  
	  *(post_prob_matrix + get_mtx_index(i,j,hmmp->nr_v));
      }
    }
    
    
    //dump_post_prob_matrix(seq_len + 2, hmmp->nr_labels, post_prob_label_matrix);

    fprintf(outfile, "Posterior probabilities:\n");
    fprintf(outfile, "pos\t");
    for(i = 0; i < hmmp->nr_labels; i++) {
      fprintf(outfile, "%c\t", hmmp->labels[i]);
    }
    fprintf(outfile, "\n");
    for(i = 1; i < seq_len + 1; i++) {
      fprintf(outfile, "%d\t", i);
      for(j = 0; j < hmmp->nr_labels; j++) {
	fprintf(outfile, "%.3f\t", *(post_prob_label_matrix + get_mtx_index(i,j,hmmp->nr_labels)));
      }
      fprintf(outfile, "\n");
    }
    fprintf(outfile, "\n");
    free(post_prob_label_matrix);
    free(post_prob_matrix);
  }
  

  /* get path + print + deallocate */
  if(print_path == YES) {
    if(seq_format == STANDARD || seq_format == FASTA) {
      if(print_path == YES) {
	path = (int*)malloc_or_die((seq_len+1) * 2 * sizeof(int));
	path_incrementable = path;
	b = 0;
	get_viterbi_path_multi(viterbi_mtx + get_mtx_index(seq_len+1,hmmp->nr_v-1,hmmp->nr_v), hmmp, viterbi_mtx,
			       seq_len + 1, hmmp->nr_v, path_incrementable, &b);
      }
      
    }
    if(seq_format == MSA_STANDARD || seq_format == PROFILE) {
      if(print_path == YES) {
	path = (int*)malloc_or_die((seq_len+1) * 2 * sizeof(int));
	path_incrementable = path;
	b = 0;
	get_viterbi_path_multi(viterbi_mtx + get_mtx_index(seq_len+1,hmmp->nr_v-1,hmmp->nr_v), hmmp, viterbi_mtx,
			       seq_len + 1, hmmp->nr_v, path_incrementable, &b);
      }
    }
    /* print to outfile */
    fprintf(outfile, "Path:\n");
    if(print_path == YES) {
      for(j = 0; j < seq_len; j++) {
	if(j % 60 == 0 && j != 0) {
	  fprintf(outfile, "\n");
	}
	fprintf(outfile, "%d ", path[j]);
      }
      fprintf(outfile, "\n");
    }
    
    /* deallocate */
    free(path);
  }

  fprintf(outfile, "\n");
  

  /* deallocate result matrix info */
  if(run_forward == YES) {
    free(forw_mtx);
    free(forw_scale);
  }
  if(run_backward == YES) {
    free(backw_mtx);
  }
  if(run_viterbi == YES || print_path == YES) {
    free(viterbi_mtx);
  }
  if(run_n_best == YES) {
    free(one_best_mtx);
    free(scaling_factors);
  }
  
}


/*********************help functions***************************************/
double get_nullmodel_score_multi(struct letter_s *seq, struct letter_s *seq_2, struct letter_s *seq_3,
				 struct letter_s *seq_4, int seq_len, int multi_scoring_method)
{
  int avg_seq_len;
  double emiss_prob, emiss_prob_2, emiss_prob_3, emiss_prob_4;
  double e1, e2, e3, e4;
  double trans_prob;
  double null_model_score;
  int letter, letter_2, letter_3, letter_4;
  int l,i;
  double t_res;
  
  /* calculate null model score for seq */
  if(null_model.nr_alphabets < 0) {
    /* use default null model */
    emiss_prob = 1.0 / (double)hmm.a_size;
    if(hmm.nr_alphabets > 1) {
      emiss_prob_2 = 1.0 / (double)hmm.a_size_2;
    }
    if(hmm.nr_alphabets > 2) {
      emiss_prob_3 = 1.0 / (double)hmm.a_size_3;
    }
    if(hmm.nr_alphabets > 3) {
      emiss_prob_4 = 1.0 / (double)hmm.a_size_4;
    }
    trans_prob = (double)(seq_len)/(double)(seq_len + 1);
    if(multi_scoring_method == JOINT_PROB) {
      if(hmm.nr_alphabets == 1) {
	null_model_score = seq_len * (log10(emiss_prob) + log10(trans_prob));
      }
      else if(hmm.nr_alphabets == 2) {
	null_model_score = seq_len * (log10(emiss_prob) + log10(emiss_prob_2) + log10(trans_prob));
      }
      else if(hmm.nr_alphabets == 3) {
	null_model_score = seq_len * (log10(emiss_prob) + log10(emiss_prob_2) + log10(emiss_prob_3) + log10(trans_prob));
      }
      else if(hmm.nr_alphabets == 2) {
	null_model_score = seq_len * (log10(emiss_prob) + log10(emiss_prob_2) + log10(emiss_prob_3) + log10(emiss_prob_4)
				      + log10(trans_prob));
      }
    }
    else {
      /* use other multialpha scoring method, not implemented yet */
      printf("Error: only JOINT_PROB scoring is implemented\n");
      exit(0);
    }
  }
  else {
    /* use specified null model */
    null_model_score = 0.0;
    for(l = 0; l < seq_len; l++) {
      letter = get_alphabet_index(&seq[l], hmm.alphabet, hmm.a_size);
      if(hmm.nr_alphabets > 1) {
	letter_2 = get_alphabet_index(&seq_2[l], hmm.alphabet_2, hmm.a_size_2);
      }
      if(hmm.nr_alphabets > 2) {
	letter_3 = get_alphabet_index(&seq_3[l], hmm.alphabet_3, hmm.a_size_3);
      }
      if(hmm.nr_alphabets > 3) {
	letter_4 = get_alphabet_index(&seq_4[l], hmm.alphabet_4, hmm.a_size_4);
      }

      if(letter >= 0) {
	e1 = null_model.emissions[letter];
	//null_model_score += log10(null_model.emissions[letter]) + log10(null_model.trans_prob);
      }
      else {
	/* need replacement letters */
	letter = get_replacement_letter_index_multi(&seq[l], &replacement_letters, 1);
	if(letter >= 0) {
	  e1 = 0.0;
	  for(i = 0; i < hmm.a_size; i++) {
	    e1 += *(hmm.replacement_letters->probs_1 + get_mtx_index(letter, i, hmm.a_size)) * null_model.emissions[i];
	  }
	  //null_model_score += log10(t_res) + log10(null_model.trans_prob);
	}
	else {
	  printf("Could not find letter %s when scoring null model\n", &seq[l]);
	  return DEFAULT;
	}
      }
      if(hmm.nr_alphabets > 1) {
	if(letter_2 >= 0) {
	  e2 = null_model.emissions_2[letter_2];
	}
	else {
	  /* need replacement letters */
	  letter_2 = get_replacement_letter_index_multi(&seq_2[l], &replacement_letters, 2);
	  if(letter_2 >= 0) {
	    e2 = 0.0;
	    for(i = 0; i < hmm.a_size_2; i++) {
	      e2 += *(hmm.replacement_letters->probs_2 + get_mtx_index(letter_2, i, hmm.a_size_2)) * null_model.emissions_2[i];
	    }
	  }
	  else {
	    printf("Could not find letter %s when scoring null model\n", &seq_2[l]);
	    return DEFAULT;
	  }
	}
      }
      if(hmm.nr_alphabets > 2) {
	if(letter_3 >= 0) {
	  e3 = null_model.emissions_3[letter_3];
	}
	else {
	  /* need replacement letters */
	  letter_3 = get_replacement_letter_index_multi(&seq_3[l], &replacement_letters, 3);
	  if(letter_3 >= 0) {
	    e3 = 0.0;
	    for(i = 0; i < hmm.a_size_3; i++) {
	      e3 += *(hmm.replacement_letters->probs_3 + get_mtx_index(letter_3, i, hmm.a_size_3)) * null_model.emissions_3[i];
	    }
	  }
	  else {
	    printf("Could not find letter %s when scoring null model\n", &seq_3[l]);
	    return DEFAULT;
	  }
	}
      }
      if(hmm.nr_alphabets > 3) {
	if(letter_4 >= 0) {
	  e4 = null_model.emissions_4[letter_4];
	}
	else {
	  /* need replacement letters */
	  letter_4 = get_replacement_letter_index_multi(&seq_4[l], &replacement_letters, 4);
	  if(letter_4 >= 0) {
	    e4 = 0.0;
	    for(i = 0; i < hmm.a_size_4; i++) {
	      e4 += *(hmm.replacement_letters->probs_4 + get_mtx_index(letter_4, i, hmm.a_size_4)) * null_model.emissions_4[i];
	    }
	  }
	  else {
	    printf("Could not find letter %s when scoring null model\n", &seq[l]);
	    return DEFAULT;
	  }
	}
      }

      if(multi_scoring_method == JOINT_PROB) {
	null_model_score += log10(e1) + log10(null_model.trans_prob);
	if(hmm.nr_alphabets > 1) {
	  null_model_score += log10(e2);
	}
	if(hmm.nr_alphabets > 2) {
	  null_model_score += log10(e3);
	}
	if(hmm.nr_alphabets > 3) {
	  null_model_score += log10(e4);
	}
      }
      else {
	/* use other multialpha scoring method, not implemented yet */
	printf("Error: only JOINT_PROB scoring is implemented\n");
	exit(0);
      }
    }
  }

  return null_model_score;
}






double get_nullmodel_score_msa_multi(int seq_len, int prf_mode, int use_prior,
				     int use_gap_shares, int normalize,
				     int scoring_method, int multi_scoring_method, double *aa_freqs, double *aa_freqs_2,
				     double *aa_freqs_3, double *aa_freqs_4)
{
  int avg_seq_len;
  double emiss_prob, emiss_prob_2, emiss_prob_3, emiss_prob_4;
  double trans_prob;
  double null_model_score;
  int letter;
  int k,l,p,m;
  double col_score, col_score_2, col_score_3, col_score_4;
  int using_default_null_model;
  double seq_normalizer, state_normalizer;

  using_default_null_model = NO;
  /* calculate null model score for seq */
  if(null_model.nr_alphabets < 0) {
    /* use default null model */
    emiss_prob = 1.0 / (double)hmm.a_size;
    if(hmm.nr_alphabets > 1) {
      emiss_prob_2 = 1.0 / (double)hmm.a_size_2;
    }
    if(hmm.nr_alphabets > 2) {
      emiss_prob_3 = 1.0 / (double)hmm.a_size_3;
    }
    if(hmm.nr_alphabets > 3) {
      emiss_prob_4 = 1.0 / (double)hmm.a_size_4;
    }
    trans_prob = (double)(seq_len)/(double)(seq_len + 1);
    null_model.a_size = hmm.a_size;
    null_model.emissions = (double*)malloc_or_die(null_model.a_size * sizeof(double));
    if(hmm.nr_alphabets > 1) {
      null_model.a_size_2 = hmm.a_size_2;
      null_model.emissions_2 = (double*)malloc_or_die(null_model.a_size_2 * sizeof(double));
    }
    if(hmm.nr_alphabets > 2) {
      null_model.a_size_3 = hmm.a_size_3;
      null_model.emissions_3 = (double*)malloc_or_die(null_model.a_size_3 * sizeof(double));
    }
    if(hmm.nr_alphabets > 3) {
      null_model.a_size_4 = hmm.a_size_4;
      null_model.emissions_4 = (double*)malloc_or_die(null_model.a_size_4 * sizeof(double));
    }
    for(k = 0; k < hmm.a_size; k++) {
      null_model.emissions[k] = emiss_prob;
      if(hmm.nr_alphabets > 1) {
	null_model.emissions_2[k] = emiss_prob_2;
      }
      if(hmm.nr_alphabets > 2) {
	null_model.emissions_3[k] = emiss_prob_3;
      }
      if(hmm.nr_alphabets > 3) {
	null_model.emissions_4[k] = emiss_prob_4;
      }
    }
    
    null_model.trans_prob = trans_prob;
    using_default_null_model = YES;
  }
  

  /* NOTE: must include other scoring methods here as well (copy from core_algorithms) */
  /* use specified null model */
  null_model_score = 0.0;
  if(scoring_method == DOT_PRODUCT) {
    l = 0;
    while(1) {

      /* set index for seq pointer */
      if(prf_mode == ALL_SEQS) {
	p = l;
      }
      else {
	p = *(msa_seq_infop->lead_columns_start + l);
      }

      /* get col_score for the different alphabets */
      col_score = 0.0;
      col_score_2 = 0.0;
      col_score_3 = 0.0;
      col_score_4 = 0.0;
      col_score = get_dp_statescore(null_model.a_size, use_gap_shares, use_prior, msa_seq_infop->msa_seq_1,
				    p, null_model.emissions,
				    0, normalize, msa_seq_infop->gap_shares);
      if(hmm.nr_alphabets > 1) {
	col_score_2 = get_dp_statescore(null_model.a_size_2, use_gap_shares, use_prior, msa_seq_infop->msa_seq_2,
					p, null_model.emissions_2,
					0, normalize, msa_seq_infop->gap_shares);
      }
      if(hmm.nr_alphabets > 2) {
	col_score_3 = get_dp_statescore(null_model.a_size_3, use_gap_shares, use_prior, msa_seq_infop->msa_seq_3,
					p, null_model.emissions_3,
					0, normalize, msa_seq_infop->gap_shares);
      }
      if(hmm.nr_alphabets > 3) {
	col_score_4 = get_dp_statescore(null_model.a_size_4, use_gap_shares, use_prior, msa_seq_infop->msa_seq_4,
					p, null_model.emissions_4,
					0, normalize, msa_seq_infop->gap_shares);
      }
      
      /* calculate total column score */
      if(multi_scoring_method == JOINT_PROB) {
	null_model_score += log10(col_score) + log10(null_model.trans_prob);
	if(hmm.nr_alphabets > 1) {
	  null_model_score += log10(col_score_2);
	}
	if(hmm.nr_alphabets > 2) {
	    null_model_score += log10(col_score_3);
	}
	if(hmm.nr_alphabets > 3) {
	  null_model_score += log10(col_score_4);
	}
      }
      else {
	/* use other multialpha scoring method, not implemented yet */
	printf("Error: only JOINT_PROB scoring is implemented\n");
	exit(0);
      }

      /* update seq index and check if we are done */ 
      l++;
      if(prf_mode == ALL_SEQS) {
	if(l >= seq_len) {
	  break;
	}
      }
      else {
	if(*(msa_seq_infop->lead_columns_start + l) == END) {
	  break;
	}
      }
    }
  }

  if(scoring_method == DOT_PRODUCT_PICASSO) {
    l = 0;
    while(1) {

      /* set index for seq pointer */
      if(prf_mode == ALL_SEQS) {
	p = l;
      }
      else {
	p = *(msa_seq_infop->lead_columns_start + l);
      }

      /* get col_score for the different alphabets */
      col_score = 0.0;
      col_score_2 = 0.0;
      col_score_3 = 0.0;
      col_score_4 = 0.0;
      col_score = get_dp_picasso_statescore(null_model.a_size, use_gap_shares, use_prior, msa_seq_infop->msa_seq_1,
				    p, null_model.emissions,
				    0, normalize, msa_seq_infop->gap_shares, aa_freqs);
      if(hmm.nr_alphabets > 1) {
	col_score_2 = get_dp_picasso_statescore(null_model.a_size_2, use_gap_shares, use_prior, msa_seq_infop->msa_seq_2,
					p, null_model.emissions_2,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_2);
      }
      if(hmm.nr_alphabets > 2) {
	col_score_3 = get_dp_picasso_statescore(null_model.a_size_3, use_gap_shares, use_prior, msa_seq_infop->msa_seq_3,
					p, null_model.emissions_3,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_3);
      }
      if(hmm.nr_alphabets > 3) {
	col_score_4 = get_dp_picasso_statescore(null_model.a_size_4, use_gap_shares, use_prior, msa_seq_infop->msa_seq_4,
					p, null_model.emissions_4,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_4);
      }
      
      /* calculate total column score */
      if(multi_scoring_method == JOINT_PROB) {
	null_model_score += log10(col_score) + log10(null_model.trans_prob);
	if(hmm.nr_alphabets > 1) {
	  null_model_score += log10(col_score_2);
	}
	if(hmm.nr_alphabets > 2) {
	    null_model_score += log10(col_score_3);
	}
	if(hmm.nr_alphabets > 3) {
	  null_model_score += log10(col_score_4);
	}
      }
      else {
	/* use other multialpha scoring method, not implemented yet */
	printf("Error: only JOINT_PROB scoring is implemented\n");
	exit(0);
      }

      /* update seq index and check if we are done */ 
      l++;
      if(prf_mode == ALL_SEQS) {
	if(l >= seq_len) {
	  break;
	}
      }
      else {
	if(*(msa_seq_infop->lead_columns_start + l) == END) {
	  break;
	}
      }
    }
  }

  if(scoring_method == PICASSO) {
    l = 0;
    while(1) {

      /* set index for seq pointer */
      if(prf_mode == ALL_SEQS) {
	p = l;
      }
      else {
	p = *(msa_seq_infop->lead_columns_start + l);
      }

      /* get col_score for the different alphabets */
      col_score = 0.0;
      col_score_2 = 0.0;
      col_score_3 = 0.0;
      col_score_4 = 0.0;
      col_score = get_picasso_statescore(null_model.a_size, use_gap_shares, use_prior, msa_seq_infop->msa_seq_1,
				    p, null_model.emissions,
				    0, normalize, msa_seq_infop->gap_shares, aa_freqs);
      if(hmm.nr_alphabets > 1) {
	col_score_2 = get_picasso_statescore(null_model.a_size_2, use_gap_shares, use_prior, msa_seq_infop->msa_seq_2,
					p, null_model.emissions_2,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_2);
      }
      if(hmm.nr_alphabets > 2) {
	col_score_3 = get_picasso_statescore(null_model.a_size_3, use_gap_shares, use_prior, msa_seq_infop->msa_seq_3,
					p, null_model.emissions_3,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_3);
      }
      if(hmm.nr_alphabets > 3) {
	col_score_4 = get_picasso_statescore(null_model.a_size_4, use_gap_shares, use_prior, msa_seq_infop->msa_seq_4,
					p, null_model.emissions_4,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_4);
      }
      
      /* calculate total column score */
      if(multi_scoring_method == JOINT_PROB) {
	null_model_score += log10(col_score) + log10(null_model.trans_prob);
	if(hmm.nr_alphabets > 1) {
	  null_model_score += log10(col_score_2);
	}
	if(hmm.nr_alphabets > 2) {
	    null_model_score += log10(col_score_3);
	}
	if(hmm.nr_alphabets > 3) {
	  null_model_score += log10(col_score_4);
	}
      }
      else {
	/* use other multialpha scoring method, not implemented yet */
	printf("Error: only JOINT_PROB scoring is implemented\n");
	exit(0);
      }

      /* update seq index and check if we are done */ 
      l++;
      if(prf_mode == ALL_SEQS) {
	if(l >= seq_len) {
	  break;
	}
      }
      else {
	if(*(msa_seq_infop->lead_columns_start + l) == END) {
	  break;
	}
      }
    }
  }

  if(scoring_method == PICASSO_SYM) {
    l = 0;
    while(1) {

      /* set index for seq pointer */
      if(prf_mode == ALL_SEQS) {
	p = l;
      }
      else {
	p = *(msa_seq_infop->lead_columns_start + l);
      }

      /* get col_score for the different alphabets */
      col_score = 0.0;
      col_score_2 = 0.0;
      col_score_3 = 0.0;
      col_score_4 = 0.0;
      col_score = get_picasso_sym_statescore(null_model.a_size, use_gap_shares, use_prior, msa_seq_infop->msa_seq_1,
				    p, null_model.emissions,
				    0, normalize, msa_seq_infop->gap_shares, aa_freqs);
      if(hmm.nr_alphabets > 1) {
	col_score_2 = get_picasso_sym_statescore(null_model.a_size_2, use_gap_shares, use_prior, msa_seq_infop->msa_seq_2,
					p, null_model.emissions_2,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_2);
      }
      if(hmm.nr_alphabets > 2) {
	col_score_3 = get_picasso_sym_statescore(null_model.a_size_3, use_gap_shares, use_prior, msa_seq_infop->msa_seq_3,
					p, null_model.emissions_3,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_3);
      }
      if(hmm.nr_alphabets > 3) {
	col_score_4 = get_picasso_sym_statescore(null_model.a_size_4, use_gap_shares, use_prior, msa_seq_infop->msa_seq_4,
					p, null_model.emissions_4,
					0, normalize, msa_seq_infop->gap_shares, aa_freqs_4);
      }
      
      /* calculate total column score */
      if(multi_scoring_method == JOINT_PROB) {
	null_model_score += log10(col_score) + log10(null_model.trans_prob);
	if(hmm.nr_alphabets > 1) {
	  null_model_score += log10(col_score_2);
	}
	if(hmm.nr_alphabets > 2) {
	    null_model_score += log10(col_score_3);
	}
	if(hmm.nr_alphabets > 3) {
	  null_model_score += log10(col_score_4);
	}
      }
      else {
	/* use other multialpha scoring method, not implemented yet */
	printf("Error: only JOINT_PROB scoring is implemented\n");
	exit(0);
      }

      /* update seq index and check if we are done */ 
      l++;
      if(prf_mode == ALL_SEQS) {
	if(l >= seq_len) {
	  break;
	}
      }
      else {
	if(*(msa_seq_infop->lead_columns_start + l) == END) {
	  break;
	}
      }
    }
  }

  if(scoring_method == SJOLANDER) {
    l = 0;
    while(1) {
       /* set index for seq pointer */
      if(prf_mode == ALL_SEQS) {
	p = l;
      }
      else {
	p = *(msa_seq_infop->lead_columns_start + l);
      }

      /* get col_score for the different alphabets */
      col_score = 0.0;
      col_score_2 = 0.0;
      col_score_3 = 0.0;
      col_score_4 = 0.0;
      col_score = get_sjolander_statescore(null_model.a_size, use_gap_shares, use_prior, msa_seq_infop->msa_seq_1,
					   p, null_model.emissions,
					   0, normalize, msa_seq_infop->gap_shares);
      if(hmm.nr_alphabets > 1) {
	col_score_2 = get_sjolander_statescore(null_model.a_size_2, use_gap_shares, use_prior, msa_seq_infop->msa_seq_2,
						 p, null_model.emissions_2,
					       0, normalize, msa_seq_infop->gap_shares);
      }
      if(hmm.nr_alphabets > 2) {
	col_score_3 = get_sjolander_statescore(null_model.a_size_3, use_gap_shares, use_prior, msa_seq_infop->msa_seq_3,
					       p, null_model.emissions_3,
						 0, normalize, msa_seq_infop->gap_shares);
      }
      if(hmm.nr_alphabets > 3) {
	  col_score_4 = get_sjolander_statescore(null_model.a_size_4, use_gap_shares, use_prior, msa_seq_infop->msa_seq_4,
						 p, null_model.emissions_4,
						 0, normalize, msa_seq_infop->gap_shares);
	}
      
      /* calculate total column score */
      if(multi_scoring_method == JOINT_PROB) {
	null_model_score += log10(col_score) + log10(null_model.trans_prob);
	if(hmm.nr_alphabets > 1) {
	  null_model_score += log10(col_score_2);
	}
	if(hmm.nr_alphabets > 2) {
	  null_model_score += log10(col_score_3);
	}
	  if(hmm.nr_alphabets > 3) {
	    null_model_score += log10(col_score_4);
	  }
      }
      else {
	/* use other multialpha scoring method, not implemented yet */
	printf("Error: only JOINT_PROB scoring is implemented\n");
	exit(0);
      }

      /* update seq index and check if we are done */ 
      l++;
      if(prf_mode == ALL_SEQS) {
	if(l >= seq_len) {
	  break;
	}
      }
      else {
	if(*(msa_seq_infop->lead_columns_start + l) == END) {
	  break;
	}
      }
    }
  }

  if(scoring_method == SJOLANDER_REVERSED) {
    l = 0;
    while(1) {
       /* set index for seq pointer */
      if(prf_mode == ALL_SEQS) {
	p = l;
      }
      else {
	p = *(msa_seq_infop->lead_columns_start + l);
      }

      /* get col_score for the different alphabets */
      col_score = 0.0;
      col_score_2 = 0.0;
      col_score_3 = 0.0;
      col_score_4 = 0.0;
      col_score = get_sjolander_reversed_statescore(null_model.a_size, use_gap_shares, use_prior, msa_seq_infop->msa_seq_1,
						    p, null_model.emissions,
						    0, normalize, msa_seq_infop->gap_shares);
      if(hmm.nr_alphabets > 1) {
	col_score_2 = get_sjolander_reversed_statescore(null_model.a_size_2, use_gap_shares, use_prior, msa_seq_infop->msa_seq_2,
							p, null_model.emissions_2,
							0, normalize, msa_seq_infop->gap_shares);
      }
      if(hmm.nr_alphabets > 2) {
	col_score_3 = get_sjolander_reversed_statescore(null_model.a_size_3, use_gap_shares, use_prior, msa_seq_infop->msa_seq_3,
							p, null_model.emissions_3,
							0, normalize, msa_seq_infop->gap_shares);
      }
      if(hmm.nr_alphabets > 3) {
	col_score_4 = get_sjolander_reversed_statescore(null_model.a_size_4, use_gap_shares, use_prior, msa_seq_infop->msa_seq_4,
							p, null_model.emissions_4,
							0, normalize, msa_seq_infop->gap_shares);
      }
      
      /* calculate total column score */
      if(multi_scoring_method == JOINT_PROB) {
	null_model_score += log10(col_score) + log10(null_model.trans_prob);
	if(hmm.nr_alphabets > 1) {
	  null_model_score += log10(col_score_2);
	}
	if(hmm.nr_alphabets > 2) {
	  null_model_score += log10(col_score_3);
	}
	  if(hmm.nr_alphabets > 3) {
	    null_model_score += log10(col_score_4);
	  }
      }
      else {
	/* use other multialpha scoring method, not implemented yet */
	printf("Error: only JOINT_PROB scoring is implemented\n");
	exit(0);
      }

      /* update seq index and check if we are done */ 
      l++;
      if(prf_mode == ALL_SEQS) {
	if(l >= seq_len) {
	  break;
	}
      }
      else {
	if(*(msa_seq_infop->lead_columns_start + l) == END) {
	  break;
	}
      }
    }
  }

  else if(scoring_method == SUBST_MTX_PRODUCT) {
    l = 0;
    while(1) {
      /* set index for seq pointer */
      if(prf_mode == ALL_SEQS) {
	p = l;
      }
      else {
	p = *(msa_seq_infop->lead_columns_start + l);
      }

      /* get col_score for the different alphabets */
      col_score = 0.0;
      col_score_2 = 0.0;
      col_score_3 = 0.0;
      col_score_4 = 0.0;
      col_score = get_subst_mtx_product_statescore(null_model.a_size, use_gap_shares, use_prior, msa_seq_infop->msa_seq_1,
						   p, null_model.emissions,
						   0, hmm.subst_mtx);
      if(hmm.nr_alphabets > 1) {
	col_score_2 = get_subst_mtx_product_statescore(null_model.a_size_2, use_gap_shares, use_prior, msa_seq_infop->msa_seq_2,
						       p, null_model.emissions_2,
						       0, hmm.subst_mtx_2);
      }
      if(hmm.nr_alphabets > 2) {
	col_score_3 = get_subst_mtx_product_statescore(null_model.a_size_3, use_gap_shares, use_prior, msa_seq_infop->msa_seq_3,
						       p, null_model.emissions_3,
						       0, hmm.subst_mtx_3);
      }
      if(hmm.nr_alphabets > 3) {
	col_score_4 = get_subst_mtx_product_statescore(null_model.a_size_4, use_gap_shares, use_prior, msa_seq_infop->msa_seq_4,
						       p, null_model.emissions_4,
						       0, hmm.subst_mtx_4);
      }
      
      /* calculate total column score */
      if(multi_scoring_method == JOINT_PROB) {
	null_model_score += log10(col_score) + log10(null_model.trans_prob);
	if(hmm.nr_alphabets > 1) {
	  null_model_score += log10(col_score_2);
	}
	if(hmm.nr_alphabets > 2) {
	  null_model_score += log10(col_score_3);
	}
	if(hmm.nr_alphabets > 3) {
	  null_model_score += log10(col_score_4);
	}
      }
      else {
	/* use other multialpha scoring method, not implemented yet */
	printf("Error: only JOINT_PROB scoring is implemented\n");
	exit(0);
      }
      
      /* update seq index and check if we are done */ 
      l++;
      if(prf_mode == ALL_SEQS) {
	if(l >= seq_len) {
	  break;
	}
      }
      else {
	if(*(msa_seq_infop->lead_columns_start + l) == END) {
	  break;
	}
      }
    }
  }
    
  else if(scoring_method == SUBST_MTX_DOT_PRODUCT) {
    printf("SMDP scoring not implemented yet\n");
    exit(0);
	   
  }
  else if(scoring_method == SUBST_MTX_DOT_PRODUCT_PRIOR) {
    printf("SMDPP scoring not implemented yet\n");
    exit(0);
  }
  
  if(using_default_null_model == YES) {
    free(null_model.emissions);
    null_model.a_size = -1;
  }
  return null_model_score;
}



void get_post_prob_matrix(double **ret_post_prob_matrixp, double forw_score, struct forward_s *forw_mtx,
			  struct backward_s *backw_mtx, int seq_len)
{
  int i,j;
  double *post_prob_matrix;
  double post_prob_score;
  
  /* must be freed by caller */
  post_prob_matrix = (double*)(malloc_or_die((seq_len+2) * hmm.nr_v * sizeof(double)));
  for(i = 0; i < seq_len + 2; i++) {
    for(j = 0; j < hmm.nr_v; j++) {
      //printf("forw_score = %f, forw_mtx = %f, backw_mtx = %f\n", forw_score,(forw_mtx + get_mtx_index(i,j,hmm.nr_v))->prob,
      //     (backw_mtx + get_mtx_index(i,j,hmm.nr_v))->prob);  
      post_prob_score = (forw_mtx + get_mtx_index(i,j,hmm.nr_v))->prob * (backw_mtx + get_mtx_index(i,j,hmm.nr_v))->prob /
	forw_score;
      *(post_prob_matrix + get_mtx_index(i,j,hmm.nr_v)) = post_prob_score;
    }
  }

  *ret_post_prob_matrixp = post_prob_matrix;
}
