#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "structs.h" /* data structures etc */
#include "funcs.h" /* function header */

#define MAX_LINE 500

//#define DEBUG_RD
//#define DEBUG_RDPRI

extern int verbose;

int read_module_multi(char*, FILE*, struct hmm_multi_s*, struct module_multi_s*, int*, int*);
void check_probs_multi(struct hmm_multi_s*);
void create_to_silent_trans_array_multi(struct hmm_multi_s *hmmp);
void create_from_silent_trans_array_multi(struct hmm_multi_s *hmmp);
void create_from_trans_array_multi(struct hmm_multi_s*);
void create_to_trans_array_multi(struct hmm_multi_s*);
void create_tot_transitions_multi(struct hmm_multi_s*);
void create_tot_trans_arrays_multi(struct hmm_multi_s *hmmp);
void add_all_from_paths_multi(int v, int w, struct hmm_multi_s *hmmp,
		   struct path_element **from_transp, struct path_element *temp_pathp, int length);
void add_all_to_paths_multi(int v, int w, struct hmm_multi_s *hmmp,
		   struct path_element **from_transp, struct path_element *temp_pathp, int length);
int read_prior_files_multi(int, struct emission_dirichlet_s*, int, FILE*);
int read_trans_prior_files_multi(int, void*, struct hmm_multi_s*, FILE*);
int silent_vertex_multi(int v, struct hmm_multi_s *hmmp);


int readhmm_check(FILE *hmmfile) {
  char  s[MAX_LINE];
  int filetype;

  filetype = SINGLE_HMM;
  while(fgets(s, MAX_LINE, hmmfile) != NULL) {
    if(strncmp(s, "NR OF ALPHABETS:",16) == 0) {
      filetype = MULTI_HMM;
      break;
    }
    else if(strncmp(s, "ALPHABET LENGTH:", 16) == 0) {
      filetype = SINGLE_HMM;
      break;
    }
  }
  rewind(hmmfile);
  return filetype;
}


void transform_singlehmmfile_to_multi(FILE *hmmfile, FILE *outfile)
{
  char s[MAX_LINE], s_2[MAX_LINE];
  if(outfile == NULL || hmmfile == NULL) {
    printf("Could not transform singlehmm to multi\n");
    exit(0);
  }

  while(1) {
    if(fgets(s, MAX_LINE, hmmfile) != NULL) {
      if(strncmp(s, "ALPHABET:", 9) == 0) {
	if(fputs("NR OF ALPHABETS: 1\n", outfile) == EOF) {
	  perror("");
	}
	if(fputs("ALPHABET 1:", outfile) == EOF) {
	  perror("");
	}
	if(fputs(&s[9], outfile) == EOF) {
	  perror("");
	}
      }
      else if(strncmp(s, "ALPHABET LENGTH:", 16) == 0) {
	if(fputs("ALPHABET LENGTH 1:", outfile) == EOF) {
	  perror("");
	}
	if(fputs(&s[16], outfile) == EOF) {
	  perror("");
	}
      }
      else if(strncmp(s, "NR OF EMISSION PRIORFILES:", 26) == 0) {
	if(fputs("NR OF EMISSION PRIORFILES 1:", outfile) == EOF) {
	  perror("");
	}
	if(fputs(&s[26], outfile) == EOF) {
	  perror("");
	}
      }
      else if(strncmp(s, "EMISSION PRIORFILES:", 20) == 0) {
	if(fputs("EMISSION PRIORFILES_1:", outfile) == EOF) {
	  perror("");
	}
	if(fputs(&s[20], outfile) == EOF) {
	  perror("");
	}
      }
      else if(strncmp(s, "Emission prior file:", 20) == 0) {
	if(fputs("Emission prior file 1:", outfile) == EOF) {
	  perror("");
	}
	if(fputs(&s[20], outfile) == EOF) {
	  perror("");
	}
      }
      else if(strncmp(s, "Emission prior scaler:", 22) == 0) {
	if(fputs("Emission prior scaler 1:", outfile) == EOF) {
	  perror("");
	}
	if(fputs(&s[22], outfile) == EOF) {
	  perror("");
	}
      }
      else if(strncmp(s, "Nr emissions =", 14) == 0) {
	if(fputs("Nr emissions 1 =", outfile) == EOF) {
	  perror("");
	}
	if(fputs(&s[14], outfile) == EOF) {
	  perror("");
	}
      }
      else if(strncmp(s, "Emission probabilities", 22) == 0) {
	if(fputs("Emission probabilities 1\n", outfile) == EOF) {
	  perror("");
	}
      }
      else {
	if(fputs(s, outfile) == EOF) {
	  perror("");
	}
      }
    }
    else {
      break;
    }
  }
}

void copy_hmm_struct(struct hmm_multi_s *hmmp, struct hmm_multi_s *retrain_hmmp)
{

  memcpy(retrain_hmmp->name, hmmp->name, 100);  /* name of the HMM */
  retrain_hmmp->constr_t = NULL; /* time of construction */
  retrain_hmmp->nr_alphabets = hmmp->nr_alphabets;
  memcpy(retrain_hmmp->alphabet, hmmp->alphabet, 1000); /* the alphabet */
  memcpy(retrain_hmmp->alphabet_2, hmmp->alphabet_2, 1000); /* the alphabet */
  memcpy(retrain_hmmp->alphabet_3, hmmp->alphabet_3, 1000); /* the alphabet */
  memcpy(retrain_hmmp->alphabet_4, hmmp->alphabet_4, 1000); /* the alphabet */
  retrain_hmmp->alphabet_type = hmmp->alphabet_type;
  retrain_hmmp->alphabet_type_2 = hmmp->alphabet_type_2;
  retrain_hmmp->alphabet_type_3 = hmmp->alphabet_type_3;
  retrain_hmmp->alphabet_type_4 = hmmp->alphabet_type_4;
  retrain_hmmp->a_size =  hmmp->a_size;
  retrain_hmmp->a_size_2 =  hmmp->a_size_2;
  retrain_hmmp->a_size_3 =  hmmp->a_size_3;
  retrain_hmmp->a_size_4 =  hmmp->a_size_4;
  retrain_hmmp->nr_m =  hmmp->nr_m;
  retrain_hmmp->nr_v =  hmmp->nr_v;
  retrain_hmmp->nr_t =  hmmp->nr_t;
  retrain_hmmp->nr_d =  hmmp->nr_d;
  retrain_hmmp->nr_dt =  hmmp->nr_dt;
  retrain_hmmp->nr_ttg =  hmmp->nr_ttg;
  retrain_hmmp->nr_tt =  hmmp->nr_tt;
  retrain_hmmp->nr_ed =  hmmp->nr_ed;
  retrain_hmmp->nr_ed_2 =  hmmp->nr_ed_2;
  retrain_hmmp->nr_ed_3 =  hmmp->nr_ed_3;
  retrain_hmmp->nr_ed_4 =  hmmp->nr_ed_4;
  retrain_hmmp->nr_tp =  hmmp->nr_tp;
  retrain_hmmp->startnode =  hmmp->startnode;
  retrain_hmmp->endnode =  hmmp->endnode;
 
  retrain_hmmp->replacement_letters =  hmmp->replacement_letters; /* point to same struct */

  
  retrain_hmmp->modules = (struct module_multi_s**)(malloc_or_die(hmmp->nr_m * sizeof(struct module_multi_s*) +
								  hmmp->nr_m * sizeof(struct module_multi_s)));
  memcpy(retrain_hmmp->modules, hmmp->modules, hmmp->nr_m * sizeof(struct module_multi_s*) +
	 hmmp->nr_m * sizeof(struct module_multi_s));
  
  
  retrain_hmmp->silent_vertices = (int*)(malloc_or_die((hmmp->nr_v + 1) * sizeof(int)));
  retrain_hmmp->locked_vertices = (int*)(malloc_or_die((hmmp->nr_v + 1) * sizeof(int)));
  memcpy(retrain_hmmp->silent_vertices, hmmp->silent_vertices, (hmmp->nr_v + 1) * sizeof(int));
  memcpy(retrain_hmmp->locked_vertices, hmmp->locked_vertices, (hmmp->nr_v + 1) * sizeof(int));
  
  retrain_hmmp->vertex_labels = (char*)(malloc_or_die(hmmp->nr_v * sizeof(char)));
  memcpy(retrain_hmmp->vertex_labels, hmmp->vertex_labels, hmmp->nr_v * sizeof(char));
  retrain_hmmp->labels = (char*)(malloc_or_die(hmmp->nr_v * sizeof(char)));
  memcpy(retrain_hmmp->labels, hmmp->labels, hmmp->nr_v * sizeof(char));
  retrain_hmmp->nr_labels = hmmp->nr_labels;
  
 
  
  retrain_hmmp->vertex_trans_prior_scalers = (double*)(malloc_or_die(hmmp->nr_v * sizeof(double)));
  memcpy(retrain_hmmp->vertex_trans_prior_scalers, hmmp->vertex_trans_prior_scalers, hmmp->nr_v * sizeof(double));
   
  retrain_hmmp->vertex_emiss_prior_scalers = (double*)(malloc_or_die(hmmp->nr_v * sizeof(double)));
  memcpy(retrain_hmmp->vertex_emiss_prior_scalers, hmmp->vertex_emiss_prior_scalers, hmmp->nr_v * sizeof(double));
  if(hmmp->nr_alphabets > 1) {
    retrain_hmmp->vertex_emiss_prior_scalers_2 = (double*)(malloc_or_die(hmmp->nr_v * sizeof(double)));
    memcpy(retrain_hmmp->vertex_emiss_prior_scalers_2, hmmp->vertex_emiss_prior_scalers_2, hmmp->nr_v * sizeof(double));
    
  }
  if(hmmp->nr_alphabets > 2) {
    retrain_hmmp->vertex_emiss_prior_scalers_3 = (double*)(malloc_or_die(hmmp->nr_v * sizeof(double)));
    memcpy(retrain_hmmp->vertex_emiss_prior_scalers_3, hmmp->vertex_emiss_prior_scalers_3, hmmp->nr_v * sizeof(double)); 
  }
  if(hmmp->nr_alphabets > 3) {
    retrain_hmmp->vertex_emiss_prior_scalers_4 = (double*)(malloc_or_die(hmmp->nr_v * sizeof(double)));
    memcpy(retrain_hmmp->vertex_emiss_prior_scalers_4, hmmp->vertex_emiss_prior_scalers_4, hmmp->nr_v * sizeof(double));
  }
  
  retrain_hmmp->transitions = (double*)(malloc_or_die(hmmp->nr_v * hmmp->nr_v * sizeof(double)));
  retrain_hmmp->log_transitions = (double*)(malloc_or_die(hmmp->nr_v * hmmp->nr_v * sizeof(double)));
  memcpy(retrain_hmmp->transitions, hmmp->transitions, hmmp->nr_v * hmmp->nr_v * sizeof(double));
  memcpy(retrain_hmmp->log_transitions, hmmp->log_transitions, hmmp->nr_v * hmmp->nr_v * sizeof(double));
  
  retrain_hmmp->emissions = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size * sizeof(double)));
  retrain_hmmp->log_emissions = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size * sizeof(double)));
  memcpy(retrain_hmmp->emissions, hmmp->emissions, hmmp->nr_v * hmmp->a_size * sizeof(double));
  memcpy(retrain_hmmp->log_emissions, hmmp->log_emissions, hmmp->nr_v * hmmp->a_size * sizeof(double));
  if(hmmp->nr_alphabets > 1) {
    retrain_hmmp->emissions_2 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_2 * 
							sizeof(double)));
    retrain_hmmp->log_emissions_2 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_2 * 
							    sizeof(double)));
    memcpy(retrain_hmmp->emissions_2, hmmp->emissions_2, hmmp->nr_v * hmmp->a_size * sizeof(double));
    memcpy(retrain_hmmp->log_emissions_2, hmmp->log_emissions_2, hmmp->nr_v * hmmp->a_size * sizeof(double));
  }
  if(hmmp->nr_alphabets > 2) {
    retrain_hmmp->emissions_3 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_3 * 
							sizeof(double)));
    retrain_hmmp->log_emissions_3 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_3 * 
							    sizeof(double)));
    memcpy(retrain_hmmp->emissions_3, hmmp->emissions_3, hmmp->nr_v * hmmp->a_size * sizeof(double));
    memcpy(retrain_hmmp->log_emissions_3, hmmp->log_emissions_3, hmmp->nr_v * hmmp->a_size * sizeof(double));
  }
  if(hmmp->nr_alphabets > 3) {
    retrain_hmmp->emissions_4 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_4 * 
							sizeof(double)));
    retrain_hmmp->log_emissions_4 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_4 * 
							    sizeof(double)));
    memcpy(retrain_hmmp->emissions_4, hmmp->emissions_4, hmmp->nr_v * hmmp->a_size * sizeof(double));
    memcpy(retrain_hmmp->log_emissions_4, hmmp->log_emissions_4, hmmp->nr_v * hmmp->a_size * sizeof(double));
  }
  
 

   /* create to_silent_trans_array */
  create_to_silent_trans_array_multi(retrain_hmmp);
  
  /* create from_trans_array */
  create_from_trans_array_multi(retrain_hmmp);
 
  /* create to_trans_array */
  create_to_trans_array_multi(retrain_hmmp);
 
  /* create tot_transitions */
  create_tot_transitions_multi(retrain_hmmp);
 
  /* create tot_to_trans_array and tot_from_trans_arrays */
  create_tot_trans_arrays_multi(retrain_hmmp);




  /* data structures */
  
  retrain_hmmp->distrib_groups = (int*)(malloc_or_die((hmmp->nr_d + hmmp->nr_v) * sizeof(int)));
  memcpy(retrain_hmmp->distrib_groups, hmmp->distrib_groups,(hmmp->nr_d + hmmp->nr_v) * sizeof(int)); 
  
  retrain_hmmp->trans_tie_groups = (int*)(malloc_or_die((hmmp->nr_t + hmmp->nr_ttg) * sizeof(struct transition_s)));
  memcpy(retrain_hmmp->trans_tie_groups, hmmp->trans_tie_groups,(hmmp->nr_t + hmmp->nr_ttg) * sizeof(struct transition_s)); 
 
  
  retrain_hmmp->emission_dirichlets =  hmmp->emission_dirichlets;
  retrain_hmmp->ed_ps = hmmp->ed_ps;
  retrain_hmmp->emission_dirichlets_2 =  hmmp->emission_dirichlets_2;
  retrain_hmmp->ed_ps_2 = hmmp->ed_ps_2;
  retrain_hmmp->emission_dirichlets_3 =  hmmp->emission_dirichlets_3;
  retrain_hmmp->ed_ps_3 = hmmp->ed_ps_3;
  retrain_hmmp->emission_dirichlets_4 =  hmmp->emission_dirichlets_4;
  retrain_hmmp->ed_ps_4 = hmmp->ed_ps_4;
    
  retrain_hmmp->subst_mtx = hmmp->subst_mtx;
  retrain_hmmp->subst_mtx_2 = hmmp->subst_mtx_2;
  retrain_hmmp->subst_mtx_3 = hmmp->subst_mtx_3;
  retrain_hmmp->subst_mtx_4 = hmmp->subst_mtx_4;
 
}

int readhmm_multialpha(FILE *file, struct hmm_multi_s *hmmp)
{
  char  s[MAX_LINE], *c;
  int i,j,k;
  int res;
  int **from_trans_array, **to_trans_array;
  int *from_trans, *to_trans, *cur;
  struct module_multi_s *module;
  struct emission_dirichlet_s *emission_priorsp, *emission_priorsp_2, *emission_priorsp_3, *emission_priorsp_4;
  void *transition_priorsp;
  int nr_priorfiles, nr_trans_priorfiles;
  int silent_counter, locked_counter;
  char *nr_trans_tiesp, *nr_distrib_tiesp;
  struct transition_s *trans_ties;
  struct transition_s trans;

  if(verbose == YES) {
    printf("reading hmm ");
  }
  /* read header */
  if(fgets(s, MAX_LINE, file) != NULL) {
  }
  /* name */
  if(fgets(s, MAX_LINE, file) != NULL) {
    strcpy(hmmp->name, &s[6]);
    if(verbose == YES) {
      printf("%s ... ", hmmp->name);
      fflush(stdout);
    }
  }


  /* creation time */
  fgets(s, MAX_LINE, file);
  /* nr of alphabets */
  if(fgets(s, MAX_LINE, file) != NULL) {
    hmmp->nr_alphabets = atoi(&s[17]);
  }
  
  if(hmmp->nr_alphabets < 1 || hmmp->nr_alphabets > 4) {
    printf("Wrong nr of alphabets\n");
    exit(0);
  }
  /* alphabet 1 */
  if(fgets(s, MAX_LINE, file) != NULL) {
    strcpy(hmmp->alphabet, &s[12]);
  }
  
  /* set alphabet to be = DISCRETE as default, this should be reset when reading the sequences */
  hmmp->alphabet_type = DISCRETE;
  
  /* alphabet length 1 */
  if(fgets(s, MAX_LINE, file) != NULL) {
    hmmp->a_size = atoi(&s[19]);
  }
  if(hmmp->nr_alphabets > 1) {
    /* alphabet 2 */
    if(fgets(s, MAX_LINE, file) != NULL) {
      strcpy(hmmp->alphabet_2, &s[12]);
    }

    /* set alphabet to be = DISCRETE as default, this should be reset when reading the sequences */
    hmmp->alphabet_type_2 = DISCRETE;

    /* alphabet length 2 */
    if(fgets(s, MAX_LINE, file) != NULL) {
      hmmp->a_size_2 = atoi(&s[19]);
    }
  }
  if(hmmp->nr_alphabets > 2) {
    /* alphabet 3 */
    if(fgets(s, MAX_LINE, file) != NULL) {
      strcpy(hmmp->alphabet_3, &s[12]);
    }

    /* set alphabet to be = DISCRETE as default, this should be reset when reading the sequences */
    hmmp->alphabet_type_3 = DISCRETE;

    /* alphabet length 3 */
    if(fgets(s, MAX_LINE, file) != NULL) {
      hmmp->a_size_3 = atoi(&s[19]);
    }
  }
  if(hmmp->nr_alphabets > 3) {
    /* alphabet 4 */
    if(fgets(s, MAX_LINE, file) != NULL) {
      strcpy(hmmp->alphabet_4, &s[12]);
    }

    /* set alphabet to be = DISCRETE as default, this should be reset when reading the sequences */
    hmmp->alphabet_type_4 = DISCRETE;

    /* alphabet length 4 */
    if(fgets(s, MAX_LINE, file) != NULL) {
      hmmp->a_size_4 = atoi(&s[19]);
    }
  }
  /* nr of modules */
  if(fgets(s, MAX_LINE, file) != NULL) {
    hmmp->nr_m = atoi(&s[15]);
  }


  /* nr of vertices */
  if(fgets(s, MAX_LINE, file) != NULL) {
    hmmp->nr_v = atoi(&s[16]);
    hmmp->transitions = (double*)(malloc_or_die(hmmp->nr_v * hmmp->nr_v * 
						sizeof(double)));
    hmmp->log_transitions = (double*)(malloc_or_die(hmmp->nr_v * hmmp->nr_v * 
						    sizeof(double)));
    init_float_mtx(hmmp->log_transitions, DEFAULT, hmmp->nr_v * hmmp->nr_v);
    hmmp->emissions = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size * 
					      sizeof(double)));
    hmmp->log_emissions = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size * 
						  sizeof(double)));
    init_float_mtx(hmmp->log_emissions, DEFAULT, hmmp->nr_v * hmmp->a_size);
    if(hmmp->nr_alphabets > 1) {
      hmmp->emissions_2 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_2 * 
						  sizeof(double)));
      hmmp->log_emissions_2 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_2 * 
						      sizeof(double)));
      init_float_mtx(hmmp->log_emissions_2, DEFAULT, hmmp->nr_v * hmmp->a_size_2);
    }
     
    if(hmmp->nr_alphabets > 2) {
      hmmp->emissions_3 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_3 * 
						sizeof(double)));
      hmmp->log_emissions_3 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_3 * 
						    sizeof(double)));
      init_float_mtx(hmmp->log_emissions_3, DEFAULT, hmmp->nr_v * hmmp->a_size_3);
    }
    if(hmmp->nr_alphabets > 3) {
      hmmp->emissions_4 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_4 * 
						sizeof(double)));
      hmmp->log_emissions_4 = (double*)(malloc_or_die(hmmp->nr_v * hmmp->a_size_4 * 
						  sizeof(double)));
      init_float_mtx(hmmp->log_emissions_4, DEFAULT, hmmp->nr_v * hmmp->a_size_4);
    }
    
    hmmp->modules = (struct module_multi_s**)(malloc_or_die(hmmp->nr_m * sizeof(struct module_multi_s*) +
							    hmmp->nr_m * sizeof(struct module_multi_s)));

    module = (struct module_multi_s*)(hmmp->modules + hmmp->nr_m);
    hmmp->silent_vertices = (int*)(malloc_or_die((hmmp->nr_v + 1) * sizeof(int)));
    hmmp->locked_vertices = (int*)(malloc_or_die((hmmp->nr_v + 1) * sizeof(int)));
  
    for(i = 0; i < hmmp->nr_v; i++) {
      *(hmmp->locked_vertices + i) = NO;
    }
    hmmp->vertex_labels = (char*)(malloc_or_die(hmmp->nr_v * sizeof(char)));
    hmmp->vertex_trans_prior_scalers = (double*)(malloc_or_die(hmmp->nr_v * sizeof(double)));
    
    hmmp->vertex_emiss_prior_scalers = (double*)(malloc_or_die(hmmp->nr_v * sizeof(double)));
    if(hmmp->nr_alphabets > 1) {
      hmmp->vertex_emiss_prior_scalers_2 = (double*)(malloc_or_die(hmmp->nr_v * sizeof(double)));
    }
    if(hmmp->nr_alphabets > 2) {
      hmmp->vertex_emiss_prior_scalers_3 = (double*)(malloc_or_die(hmmp->nr_v * sizeof(double)));
    }
    if(hmmp->nr_alphabets > 3) {
      hmmp->vertex_emiss_prior_scalers_4 = (double*)(malloc_or_die(hmmp->nr_v * sizeof(double)));
    }
  }
  /* nr of transitions */
  if(fgets(s, MAX_LINE, file) != NULL) {
    hmmp->nr_t = atoi(&s[19]);
  }
  /* nr of distribution groups */
  if(fgets(s, MAX_LINE, file) != NULL) {
    hmmp->nr_d = atoi(&s[27]);
    hmmp->distrib_groups = (int*)(malloc_or_die((hmmp->nr_d + hmmp->nr_v) * sizeof(int)));
  }

  /* nr of trans tie groups */
  if(fgets(s, MAX_LINE, file) != NULL) {
    hmmp->nr_ttg = atoi(&s[29]);
    hmmp->trans_tie_groups = (int*)(malloc_or_die((hmmp->nr_t + hmmp->nr_ttg) * sizeof(struct transition_s)));
  }

  /* nr of emission priorfiles */
  if(fgets(s, MAX_LINE, file) != NULL) {
    nr_priorfiles = atoi(&s[29]);
    hmmp->nr_ed = nr_priorfiles;
    emission_priorsp = malloc_or_die(nr_priorfiles * sizeof(struct emission_dirichlet_s));
    hmmp->emission_dirichlets = emission_priorsp;
    hmmp->ed_ps = malloc_or_die(hmmp->nr_v * sizeof(struct emission_dirichlet_s*));
  }
  /* read the emission priorfiles */
  if(read_prior_files_multi(nr_priorfiles, emission_priorsp, hmmp->a_size, file) < 0) {
    printf("Could not read emission priorfiles\n");
    exit(-1);
  }

  if(hmmp->nr_alphabets > 1) {
    /* nr of emission priorfiles */
    if(fgets(s, MAX_LINE, file) != NULL) {
      nr_priorfiles = atoi(&s[29]);
      hmmp->nr_ed_2 = nr_priorfiles;
      emission_priorsp_2 = malloc_or_die(nr_priorfiles * sizeof(struct emission_dirichlet_s));
      hmmp->emission_dirichlets_2 = emission_priorsp_2;
      hmmp->ed_ps_2 = malloc_or_die(hmmp->nr_v * sizeof(struct emission_dirichlet_s*));
    }
    /* read the emission priorfiles */
    if(read_prior_files_multi(nr_priorfiles, emission_priorsp_2, hmmp->a_size_2, file) < 0) {
      printf("Could not read emission priorfiles\n");
      exit(-1);
    }
  }

  if(hmmp->nr_alphabets > 2) {
    /* nr of emission priorfiles */
    if(fgets(s, MAX_LINE, file) != NULL) {
      nr_priorfiles = atoi(&s[29]);
      hmmp->nr_ed_3 = nr_priorfiles;
      emission_priorsp_3 = malloc_or_die(nr_priorfiles * sizeof(struct emission_dirichlet_s));
      hmmp->emission_dirichlets_3 = emission_priorsp_3;
      hmmp->ed_ps_3 = malloc_or_die(hmmp->nr_v * sizeof(struct emission_dirichlet_s*));
    }
    /* read the emission priorfiles */
    if(read_prior_files_multi(nr_priorfiles, emission_priorsp_3, hmmp->a_size_3, file) < 0) {
      printf("Could not read emission priorfiles\n");
      exit(-1);
    }
  }

  if(hmmp->nr_alphabets > 3) {
    /* nr of emission priorfiles */
    if(fgets(s, MAX_LINE, file) != NULL) {
      nr_priorfiles = atoi(&s[29]);
      hmmp->nr_ed_4 = nr_priorfiles;
      emission_priorsp_4 = malloc_or_die(nr_priorfiles * sizeof(struct emission_dirichlet_s));
      hmmp->emission_dirichlets_4 = emission_priorsp_4;
      hmmp->ed_ps_4 = malloc_or_die(hmmp->nr_v * sizeof(struct emission_dirichlet_s*));
    }
    /* read the emission priorfiles */
    if(read_prior_files_multi(nr_priorfiles, emission_priorsp_4, hmmp->a_size_4, file) < 0) {
      printf("Could not read emission priorfiles\n");
      exit(-1);
    }
  }

   /* nr of transition priorfiles */
  if(fgets(s, MAX_LINE, file) != NULL) {
    nr_trans_priorfiles = atoi(&s[29]);
    transition_priorsp = NULL;
    /* not implemented yet */
  }
  /* read the transition priorfiles */
  if(read_trans_prior_files_multi(nr_trans_priorfiles, transition_priorsp, hmmp, file) < 0) {
    printf("Could not read transition priorfiles\n");
    exit(-1);
  }

  /* empty row */
  if(fgets(s, MAX_LINE, file) != NULL) {
  }
  /* empty row */
  if(fgets(s, MAX_LINE, file) != NULL) {
  }
  /* reads ****************Modules*****************/
  if(fgets(s, MAX_LINE, file) != NULL) {
  }
  /* read the modules */
  silent_counter = 0;
  locked_counter = 0;
  for(i = 0; i < hmmp->nr_m; i++) {
    *(hmmp->modules + i) = module;
    if((res = read_module_multi(s, file, hmmp, module, &silent_counter, &locked_counter)) < 0) {
      printf("Could not read modules\n");
      exit(-1);
    }
    module++;
  }

  *(hmmp->silent_vertices + silent_counter) = END;
  *(hmmp->locked_vertices + hmmp->nr_v) = END;
#ifdef DEBUG_RD
  //dump_locked_vertices(hmmp);
  //dump_silent_vertices(hmmp);
  //dump_multi_modules(hmmp);
#endif
  
  /* empty row */
  if(fgets(s, MAX_LINE, file) != NULL) {
  }
  
  /* reads ****************Emission distribution groups*****************/
  if(fgets(s, MAX_LINE, file) != NULL) {
  }
  

  /* read the distribution groups */
  cur = hmmp->distrib_groups;
  for(i = 0; i < hmmp->nr_d; i++) {
    j = 0;
    if(fgets(s, MAX_LINE, file) != NULL) {
      while(1) {
	if(s[j] == ':') {
	  break;
	}
	j++;
      }
      j++;
      j++;
      while(1) {
	*cur = atoi(&s[j]);
	cur++;
	while(s[j] != ' ' && s[j] != '\n') {
	  j++;
	}
	while(s[j] == ' ') {
	  j++;
	}
	if(s[j] == '\n') {
	  break;
	}
      }
      *cur = END;
      cur++;
    }
  }
  
  /* empty row */
  if(fgets(s, MAX_LINE, file) != NULL) {
  }

  /* reads ****************Transition tie groups*****************/
  if(fgets(s, MAX_LINE, file) != NULL) {
  }
  /* read the trans tie groups */
  trans_ties = hmmp->trans_tie_groups;
  for(i = 0; i < hmmp->nr_ttg; i++) {
    if(fgets(s, MAX_LINE, file) != NULL && s[0] != '\n') {
      j = 0;
      while(1) {
	if(s[j] == ':') {
	  break;
	}
	j++;
      }
      j++;
      j++;
      while(1) {
	trans.from_v = atoi(&s[j]);
	while(s[j] != '>') {
	  j++;
	}
	j++;
	trans.to_v = atoi(&s[j]);
	memcpy(trans_ties, &trans, sizeof(struct transition_s));
	trans_ties++;
	while(s[j] != ' ' && s[j] != '\n') {
	  j++;
	}
	while(s[j] == ' ') {
	  j++;
	}
	if(s[j] == '\n') {
	  break;
	}
      }
      trans.to_v = END;
      trans.from_v = END;
      memcpy(trans_ties, &trans, sizeof(struct transition_s));
      trans_ties++;
    }
    else {
      hmmp->nr_ttg = i;
      break;
    }
  }
#ifdef DEBUG_RD
  dump_distrib_groups(hmmp->distrib_groups, hmmp->nr_d);
  dump_trans_tie_groups(hmmp->trans_tie_groups, hmmp->nr_ttg);
#endif
  /* create to_silent_trans_array */
  create_to_silent_trans_array_multi(hmmp);
  
  /* create from_trans_array */
  create_from_trans_array_multi(hmmp);
 
  /* create to_trans_array */
  create_to_trans_array_multi(hmmp);
 
  /* create tot_transitions */
  create_tot_transitions_multi(hmmp);
 
  /* create tot_to_trans_array and tot_from_trans_arrays*/
  create_tot_trans_arrays_multi(hmmp);

  /* get the set of labels and the number of labels */
  get_set_of_labels_multi(hmmp);

#ifdef DEBUG_RD 
  dump_emiss_matrix(hmmp->nr_v, hmmp->a_size, hmmp->emissions);
  dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
  printf("hmmp->emission_dirichlets = %x\n", hmmp->emission_dirichlets);
  for(i = 0; i < hmmp->nr_v; i++) {
    printf("hmmp->ed_ps for vertex %d = %x\n", i, *(hmmp->ed_ps + i));
  }
#endif

  /* make sure all probabilities are legal*/
  //check_probs_multi(hmmp);

  if(verbose == YES) {
    printf("done\n");
  }
}


/************************read_module_multi *************************************/
int read_module_multi(char *s, FILE *file, struct hmm_multi_s *hmmp, struct module_multi_s *modulep,
		int *silent_counter, int *locked_counter)
{
  int nr_v, nr_t, nr_e, nr_e2, nr_e3, nr_e4, nr_et;
  int i,j,k;
  int from_v, to_v;
  double prob, log_prob;
  char type[50];
  char prifile_name[500], prifile_name_2[500], prifile_name_3[500], prifile_name_4[500];
  char *p, *probp;
  struct emission_dirichlet_s *priorsp, *priorsp_2, *priorsp_3, *priorsp_4;
  int silent_vertex;

  //i = 0;

  /* module name */
  if(fgets(s, MAX_LINE, file) != NULL) {
    strcpy(modulep->name, &s[8]);
#ifdef DEBUG_RD
    printf("module name %s", s);
#endif 
  }
  /* module type */
  if(fgets(s, MAX_LINE, file) != NULL) {
    strcpy(type, &s[6]);
    if(strncmp(type, "Singlenode", 10) == 0) {
      modulep->type = SINGLENODE;
    }
    else if(strncmp(type, "Cluster", 7) == 0) {
      modulep->type = CLUSTER;
    }
    else if(strncmp(type, "Forward_std", 11) == 0) {
      modulep->type = FORWARD_STD;
    }
    else if(strncmp(type, "Forward_alt", 11) == 0) {
      modulep->type = FORWARD_ALT;
    }
    else if(strncmp(type, "Singleloop", 10) == 0) {
      modulep->type = SINGLELOOP;
    }
    else if(strncmp(type, "Profile7", 8) == 0) {
      modulep->type = PROFILE7;
    }
    else if(strncmp(type, "Profile9", 8) == 0) {
      modulep->type = PROFILE9;
    }
    else {
      printf("Error: module is of unknown type\n");
      exit(-1);
    }
   
  }

  /* nr vertices */
  if(fgets(s, MAX_LINE, file) != NULL) {
    modulep->nr_v = atoi(&s[12]);
    modulep->vertices = (int*)(malloc_or_die(modulep->nr_v * sizeof(int)));
  }
  /* emission prior file */
  if(fgets(s, MAX_LINE, file) != NULL) {
    strcpy(prifile_name, (&s[23]));
    if((p = strstr(prifile_name, "\n")) != NULL) {
      *p = '\0';
    }
    if(strncmp(prifile_name, "null", 4) == 0) {
      strcpy(modulep->priorfile_name, "null");
      priorsp = NULL;
    }
    else { 
      strcpy(modulep->priorfile_name, prifile_name);
      strcat(modulep->priorfile_name, "\0");
      for(i = 0; i < hmmp->nr_ed; i++) {
	priorsp = (hmmp->emission_dirichlets + i);
	if((strncmp(prifile_name, priorsp->name, 200)) == 0) {
	  /* keep this priorsp */
	  break;
	}
	else {
#ifdef DEBUG_RD
	  printf("prifile_name = %s\n", prifile_name);
	  printf("priorsp->name = %s\n", priorsp->name);
#endif
	}
	if(i == hmmp->nr_ed-1) /* no name equals priorfile_name */{
	  printf("Couldn't find emission priorfile '%s'\n", prifile_name);
	  exit(-1);
	}
      }
    }
  }

  if(hmmp->nr_alphabets > 1) {
    /* emission prior file */
    if(fgets(s, MAX_LINE, file) != NULL) {
      strcpy(prifile_name_2, (&s[23]));
      if((p = strstr(prifile_name_2, "\n")) != NULL) {
	*p = '\0';
      }
      if(strncmp(prifile_name_2, "null", 4) == 0) {
	strcpy(modulep->priorfile_name_2, "null");
	priorsp_2 = NULL;
      }
      else { 
	strcpy(modulep->priorfile_name_2, prifile_name_2);
	strcat(modulep->priorfile_name_2, "\0");
	for(i = 0; i < hmmp->nr_ed_2; i++) {
	  priorsp_2 = (hmmp->emission_dirichlets_2 + i);
	  if((strncmp(prifile_name_2, priorsp_2->name, 200)) == 0) {
	    /* keep this priorsp */
	    break;
	  }
	  else {
#ifdef DEBUG_RD
	    printf("prifile_name_2 = %s\n", prifile_name_2);
	    printf("priorsp_2->name = %s\n", priorsp_2->name);
#endif
	  }
	  if(i == hmmp->nr_ed_2 - 1) /* no name equals priorfile_name */{
	    printf("Couldn't find emission priorfile '%s'\n", prifile_name_2);
	    exit(-1);
	  }
	}
      }
    }
  }
  
  if(hmmp->nr_alphabets > 2) {
    /* emission prior file */
    if(fgets(s, MAX_LINE, file) != NULL) {
      strcpy(prifile_name_3, (&s[23]));
      if((p = strstr(prifile_name_3, "\n")) != NULL) {
	*p = '\0';
      }
      if(strncmp(prifile_name_3, "null", 4) == 0) {
	strcpy(modulep->priorfile_name_3, "null");
	priorsp_3 = NULL;
      }
      else { 
	strcpy(modulep->priorfile_name_3, prifile_name_3);
	strcat(modulep->priorfile_name_3, "\0");
	for(i = 0; i < hmmp->nr_ed_3; i++) {
	  priorsp_3 = (hmmp->emission_dirichlets_3 + i);
	  if((strncmp(prifile_name_3, priorsp_3->name, 200)) == 0) {
	    /* keep this priorsp */
	    break;
	  }
	  else {
#ifdef DEBUG_RD
	    printf("prifile_name_3 = %s\n", prifile_name_3);
	    printf("priorsp_3->name = %s\n", priorsp_3->name);
#endif
	  }
	  if(i == hmmp->nr_ed_3 - 1) /* no name equals priorfile_name */{
	    printf("Couldn't find emission priorfile '%s'\n", prifile_name_3);
	    exit(-1);
	  }
	}
      }
    }
  }

  if(hmmp->nr_alphabets > 3) {
    /* emission prior file */
    if(fgets(s, MAX_LINE, file) != NULL) {
      strcpy(prifile_name_4, (&s[23]));
      if((p = strstr(prifile_name_4, "\n")) != NULL) {
	*p = '\0';
      }
      if(strncmp(prifile_name_4, "null", 4) == 0) {
	strcpy(modulep->priorfile_name_4, "null");
	priorsp_4 = NULL;
      }
      else { 
	strcpy(modulep->priorfile_name_4, prifile_name_4);
	strcat(modulep->priorfile_name_4, "\0");
	for(i = 0; i < hmmp->nr_ed_4; i++) {
	  priorsp_4 = (hmmp->emission_dirichlets_4 + i);
	  if((strncmp(prifile_name_4, priorsp_4->name, 200)) == 0) {
	    /* keep this priorsp */
	    break;
	  }
	  else {
#ifdef DEBUG_RD
	    printf("prifile_name_4 = %s\n", prifile_name_4);
	    printf("priorsp_4->name = %s\n", priorsp_4->name);
#endif
	  }
	  if(i == hmmp->nr_ed_4 - 1) /* no name equals priorfile_name */{
	    printf("Couldn't find emission priorfile '%s'\n", prifile_name_4);
	    exit(-1);
	  }
	}
      }
    }
  }

  /* transition prior file */
  if(fgets(s, MAX_LINE, file) != NULL) {
    /* not implemented yet */
  }
  /* empty line */
  if(fgets(s, MAX_LINE, file) != NULL) {
  }

#ifdef DEBUG_RD
  printf("modulep->nr_v = %d\n", modulep->nr_v);
#endif  
/* loop over the vertices */
   for(i = 0; i < modulep->nr_v; i++) {
    /* Vertex nr */
    if(fgets(s, MAX_LINE, file) != NULL) {
#ifdef DEBUG_RD      
      printf("Vertex nr: %s", s);
#endif
      from_v = atoi(&s[7]);
      *(modulep->vertices + i) = from_v;
      /* connect this vertex to its priorfile */
      *(hmmp->ed_ps + from_v) = priorsp;
      if(hmmp->nr_alphabets > 1) {
	*(hmmp->ed_ps_2 + from_v) = priorsp_2;
      }
      if(hmmp->nr_alphabets > 2) {
	*(hmmp->ed_ps_3 + from_v) = priorsp_3;
      }
      if(hmmp->nr_alphabets > 3) {
	*(hmmp->ed_ps_4 + from_v) = priorsp_4;
      }
    }
    else {
      printf("Error in hmm spec\n");
      exit(0);
    }
    /* Vertex type */
    if(fgets(s, MAX_LINE, file) != NULL) {
#ifdef DEBUG_RD
      printf("Vertex type: %s", s);
#endif
      strcpy(type, &s[13]);
      if(modulep->type == PROFILE7 || modulep->type == PROFILE9) {
	modulep->v_type = PROFILEV;
      }
      if(strncmp(type, "standard", 8) == 0) {
	if(modulep->type != PROFILE7 && modulep->type != PROFILE9) {
	  modulep->v_type = STANDARDV;
	  silent_vertex = NO;
	}
      }
      else if(strncmp(type, "silent", 6) == 0) {
	silent_vertex = YES;
	if(modulep->type != PROFILE7 && modulep->type != PROFILE9) {
	  modulep->v_type = SILENTV;
	}
	*(hmmp->silent_vertices + *silent_counter) = from_v;
	*silent_counter = *silent_counter + 1;
      }
      else if(strncmp(type, "locked", 5) == 0) {
	modulep->v_type = LOCKEDV;
	*(hmmp->locked_vertices + from_v) = YES;
	*locked_counter = *locked_counter + 1;
	silent_vertex = NO;
      }
      else if(strncmp(type, "start", 5) == 0) {
	modulep->v_type = STARTV;
	hmmp->startnode = from_v;
      }
      else if(strncmp(type, "end", 3) == 0) {
	modulep->v_type = ENDV;
	hmmp->endnode = from_v;
      }
      else {
	printf("Error: vertex type is undefined\n");
	printf("vertex type = %s\n", type);
	printf("from_v = %d\n", from_v);
	exit(-1);
      }
    }
    /* Vertex label */
    if(fgets(s, MAX_LINE, file) != NULL) {
      *(hmmp->vertex_labels + from_v) = s[14];
    }
    /* transition prior scaler */
    if(fgets(s, MAX_LINE, file) != NULL) {
      *(hmmp->vertex_trans_prior_scalers + from_v) = atof(&(s[25]));
    }
    /* emission prior scaler */
    if(fgets(s, MAX_LINE, file) != NULL) {
      *(hmmp->vertex_emiss_prior_scalers + from_v) = atof(&(s[25]));
    }
    if(hmmp->nr_alphabets > 1) {
      /* emission prior scaler */
      if(fgets(s, MAX_LINE, file) != NULL) {
	*(hmmp->vertex_emiss_prior_scalers_2 + from_v) = atof(&(s[25]));
      }
    }
    if(hmmp->nr_alphabets > 2) {
      /* emission prior scaler */
      if(fgets(s, MAX_LINE, file) != NULL) {
	*(hmmp->vertex_emiss_prior_scalers_3 + from_v) = atof(&(s[25]));
      }
    }
    if(hmmp->nr_alphabets > 3) {
      /* emission prior scaler */
      if(fgets(s, MAX_LINE, file) != NULL) {
	*(hmmp->vertex_emiss_prior_scalers_4 + from_v) = atof(&(s[25]));
      }
    }
    /* Nr transitions */
    if(fgets(s, MAX_LINE, file) != NULL) {
      nr_t = atoi(&s[17]);
    }
    /* Nr end transitions */
    if(fgets(s, MAX_LINE, file) != NULL) {
      nr_et = atoi(&s[21]);
    }
    /* Nr emissions */
    if(fgets(s, MAX_LINE, file) != NULL) {
      nr_e = atoi(&s[17]);
    }
    if(hmmp->nr_alphabets > 1) {
      /* Nr emissions */
      if(fgets(s, MAX_LINE, file) != NULL) {
	nr_e2 = atoi(&s[17]);
      }
    }
    if(hmmp->nr_alphabets > 2) {
      /* Nr emissions */
      if(fgets(s, MAX_LINE, file) != NULL) {
	nr_e3 = atoi(&s[17]);
      }
    }
    if(hmmp->nr_alphabets > 3) {
      /* Nr emissions */
      if(fgets(s, MAX_LINE, file) != NULL) {
	nr_e4 = atoi(&s[17]);
      }
    }
    /* read transition probabilities */
    fgets(s, MAX_LINE, file);
    for(j = 0; j < nr_t; j++) {
      if(fgets(s, MAX_LINE, file) != NULL) {
	to_v = atoi(&s[8]);
	if(to_v < 10 ) {
	  prob = (double)(atof(&s[11]));
	}
	else if(to_v < 100) {
	  prob = (double)(atof(&s[12]));
	 }
	else if(to_v < 1000) {
	  prob = (double)(atof(&s[13]));
	}
	else if(to_v < 10000) {
	  prob = (double)(atof(&s[14]));
	 }
	else {
	  printf("Sorry, reader cannot handle HMMs with more than 10000 states\n");
	   exit(0);
	}
	if(prob != 0.0) {
	  log_prob = log10(prob);
	}
	else {
	  log_prob = DEFAULT;
	}
#ifdef DEBUG_RD
	printf("prob from %d to %d = %f\n", from_v, to_v, prob);
#endif
	*(hmmp->transitions + get_mtx_index(from_v, to_v, hmmp->nr_v)) = prob;
	*(hmmp->log_transitions + get_mtx_index(from_v, to_v, hmmp->nr_v)) = log_prob;
      }
    }
    
    /* read end transition probabilities */
    fgets(s, MAX_LINE, file);
    for(j = 0; j < nr_et; j++) {
      if(fgets(s, MAX_LINE, file) != NULL) {
	to_v = atoi(&s[8]);
	if(to_v < 10 ) {
	  prob = (double)(atof(&s[11]));
	}
	else if(to_v < 100) {
	  prob = (double)(atof(&s[12]));
	}
	else if(to_v < 1000) {
	  prob = (double)(atof(&s[13]));
	}
	else if(to_v < 10000) {
	  prob = (double)(atof(&s[14]));
	}
	else {
	  printf("Sorry, reader cannot handle HMMs with more than 10000 states\n");
	  exit(0);
	}
	if(prob != 0.0) {
	  log_prob = log10(prob);
	}
	else {
	  log_prob = DEFAULT;
	}
#ifdef DEBUG_RD
	printf("prob from %d to %d = %f\n", from_v, to_v, prob);
#endif
	*(hmmp->transitions + get_mtx_index(from_v, to_v, hmmp->nr_v)) = prob;
	*(hmmp->log_transitions + get_mtx_index(from_v, to_v, hmmp->nr_v)) = log_prob;
      }
    }
    /* read emission probabilities */
    fgets(s, MAX_LINE, file);
    for(j = 0; j < nr_e; j++) {
      if(fgets(s, MAX_LINE, file) != NULL) {
#ifdef DEBUG_RD
	printf("%s", s);
#endif
	k = 0;
	while(s[k] != ' ') {
	  k++;
	}
	if(k > 7) {
	  printf("Cannot read hmm file, please check hmm specification\n");
	}
	k++;
	prob = (double)(atof(&s[k]));
	if(prob != 0.0) {
	  log_prob = log10(prob);
	}
	else {
	  log_prob = DEFAULT;
	}
	if(silent_vertex == YES) {
	  prob = SILENT;
	  log_prob = SILENT;
	}
	*(hmmp->emissions + get_mtx_index(from_v, j, hmmp->a_size)) = prob;
	*(hmmp->log_emissions + get_mtx_index(from_v, j, hmmp->a_size)) = log_prob;
      }
    }
    fgets(s, MAX_LINE, file);
    
    if(hmmp->nr_alphabets > 1) {
      fgets(s, MAX_LINE, file);
      for(j = 0; j < nr_e2; j++) {
	if(fgets(s, MAX_LINE, file) != NULL) {
	  k = 0;
	  while(s[k] != ' ') {
	    k++;
	  }
	  if(k > 7) {
	    printf("Cannot read hmm file, please check hmm specification\n");
	  }
	  k++;
	  prob = (double)(atof(&s[k]));
	  if(prob != 0.0) {
	    log_prob = log10(prob);
	  }
	  else {
	    log_prob = DEFAULT;
	  }
	  if(silent_vertex == YES) {
	    prob = SILENT;
	    log_prob = SILENT;
	  }
	  *(hmmp->emissions_2 + get_mtx_index(from_v, j, hmmp->a_size_2)) = prob;
	  *(hmmp->log_emissions_2 + get_mtx_index(from_v, j, hmmp->a_size_2)) = log_prob;
	}
      }
      fgets(s, MAX_LINE, file);
    }
    if(hmmp->nr_alphabets > 2) {
      fgets(s, MAX_LINE, file);
      for(j = 0; j < nr_e3; j++) {
	if(fgets(s, MAX_LINE, file) != NULL) {
	  k = 0;
	  while(s[k] != ' ') {
	    k++;
	  }
	  if(k > 7) {
	    printf("Cannot read hmm file, please check hmm specification\n");
	  }
	  k++;
	  prob = (double)(atof(&s[k]));
	  if(prob != 0.0) {
	    log_prob = log10(prob);
	  }
	  else {
	    log_prob = DEFAULT;
	  }
	  if(silent_vertex == YES) {
	    prob = SILENT;
	    log_prob = SILENT;
	  }
	  *(hmmp->emissions_3 + get_mtx_index(from_v, j, hmmp->a_size_3)) = prob;
	  *(hmmp->log_emissions_3 + get_mtx_index(from_v, j, hmmp->a_size_3)) = log_prob;
	}
      }
      fgets(s, MAX_LINE, file);
    }

    if(hmmp->nr_alphabets > 3) {
      fgets(s, MAX_LINE, file);
      for(j = 0; j < nr_e4; j++) {
	if(fgets(s, MAX_LINE, file) != NULL) {
	  k = 0;
	  while(s[k] != ' ') {
	    k++;
	  }
	  if(k > 7) {
	    printf("Cannot read hmm file, please check hmm specification\n");
	  }
	  k++;
	  prob = (double)(atof(&s[k]));
	  if(prob != 0.0) {
	    log_prob = log10(prob);
	  }
	  else {
	    log_prob = DEFAULT;
	  }
	  if(silent_vertex == YES) {
	    prob = SILENT;
	    log_prob = SILENT;
	  }
	  
	  *(hmmp->emissions_4 + get_mtx_index(from_v, j, hmmp->a_size_4)) = prob;
	  *(hmmp->log_emissions_4 + get_mtx_index(from_v, j, hmmp->a_size_4)) = log_prob;
	}
      }
      fgets(s, MAX_LINE, file);
    }
    
    silent_vertex = NO;
  }
  
  /* read ---------------------------------------- */
  fgets(s, MAX_LINE, file);
  return 0;
  
}


void create_to_silent_trans_array_multi(struct hmm_multi_s *hmmp)
{
  int v,w;
  int malloc_size;
  int *values;

  malloc_size = 0;
  for(v = 0; v < hmmp->nr_v; v++) {
    for(w = 0; w < hmmp->nr_v; w++) {
      if(*(hmmp->transitions + get_mtx_index(v,w,hmmp->nr_v)) != 0 && silent_vertex_multi(w,hmmp) == YES) {
	malloc_size++;
      }
    }
    malloc_size++;
  }
  
  hmmp->to_silent_trans_array = (int**)malloc_or_die(hmmp->nr_v * sizeof(int*) + malloc_size * sizeof(int));
  values = (int*)(hmmp->to_silent_trans_array + hmmp->nr_v);

  for(v = 0; v < hmmp->nr_v; v++) {
    *(hmmp->to_silent_trans_array + v) = values;
    for(w = 0; w < hmmp->nr_v; w++) {
      if(*(hmmp->transitions + get_mtx_index(v,w,hmmp->nr_v)) != 0 && silent_vertex_multi(w,hmmp) == YES) {
	*values = w;
	values++;
      }
    }
    *values = END;
    values++;
  }

#ifdef DEBUG_RD
  dump_to_silent_trans_array(hmmp->nr_v, hmmp->to_silent_trans_array);
#endif
}

/* Go through transmission matrix and get all probabilities that are not 0
 * into from_trans_array */
void create_from_trans_array_multi(struct hmm_multi_s *hmmp)
{
  int v,w,*xp;
  int has_to_trans;
  int array_head_size, array_tail_size;
  struct path_element **from_trans_array, *from_trans, *temp_path;

  array_tail_size = 0;
  array_head_size = hmmp->nr_v;

  /* estimate how much space we need to store transitions */
  array_tail_size = (hmmp->nr_t/hmmp->nr_v + 3 + MAX_GAP_SIZE) * MAX_GAP_SIZE/2 * hmmp->nr_v;
  
#ifdef DEBUG_RD
  printf("array_head_size, array_tail_size = %d, %d\n", array_head_size, array_tail_size);
#endif  
  from_trans_array = (struct path_element**)
    (malloc_or_die(array_head_size * sizeof(struct path_element*) +
		   (array_tail_size + hmmp->nr_v) * sizeof(struct path_element)));
  from_trans = (struct path_element*)(from_trans_array + hmmp->nr_v);
  hmmp->from_trans_array = from_trans_array;
  
  /* find all paths and add them to from_trans_array */
  for(v = 0; v < hmmp->nr_v; v++) /* to-vertex */ {
    *from_trans_array = from_trans;
    if(silent_vertex_multi(v, hmmp) == YES) {
      from_trans->vertex = END;
      from_trans->next = NULL;
      from_trans++;
      from_trans_array++;
      continue;
    }
    for(w = 0; w < hmmp->nr_v; w++) /* from-vertex */ {
      if(silent_vertex_multi(w,hmmp) == YES) {
	continue;
      }
      temp_path = (struct path_element*)(malloc_or_die(1000 * sizeof(struct path_element)));
      add_all_from_paths_multi(w, v, hmmp, &from_trans, temp_path, 0);
      free(temp_path);
    }
    from_trans->vertex = END;
    from_trans->next = NULL;
    from_trans++;
    from_trans_array++;
  }
#ifdef DEBUG_RD
  dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
  dump_from_trans_array(hmmp->nr_v, hmmp->from_trans_array);
#endif
}

void add_all_from_paths_multi(int v, int w, struct hmm_multi_s *hmmp,
			struct path_element **from_transp, struct path_element *temp_pathp,
			int length)
{
  int i,j;
  int *xp;
  struct path_element p_el, cur_p_el;
  
  if(length > MAX_GAP_SIZE) {
    return;
  }

  cur_p_el.vertex = v;
  cur_p_el.next = NULL;
  memcpy(temp_pathp + length, &cur_p_el, sizeof(struct path_element));

  if(*(hmmp->transitions + get_mtx_index(v, w, hmmp->nr_v)) != 0.0) {
#ifdef DEBUG_RD
    printf("adding path: ");
#endif
    /* direct path to w, add total path */
    for(i = 0; i < length; i++) {
      p_el.vertex = (temp_pathp + i)->vertex;
#ifdef DEBUG_RD
      printf("%d ", p_el.vertex);
#endif
      p_el.next = (*from_transp) + 1;
      memcpy(*from_transp, &p_el, sizeof(struct path_element));
      (*from_transp)++;
    }
    memcpy(*from_transp, &cur_p_el, sizeof(struct path_element));
#ifdef DEBUG_RD
    printf("%d %d\n", cur_p_el.vertex, w);
#endif    
    (*from_transp)++;
  }
  xp = *(hmmp->to_silent_trans_array + v);
  while(*xp != END) {
    add_all_from_paths_multi(*xp, w, hmmp, from_transp, temp_pathp, length + 1);
    xp++;
  }

}


/* Go through transmission matrix and get all probabilities that are not 0
 * into to_trans_array */
void create_to_trans_array_multi(struct hmm_multi_s *hmmp)
{
  int v,w,*xp;
  int has_to_trans;
  int array_head_size, array_tail_size;
  struct path_element **to_trans_array, *to_trans, *temp_path;

  array_tail_size = 0;
  array_head_size = hmmp->nr_v;
  
  /* estimate how much space we need to store transitions */
  array_tail_size = (hmmp->nr_t/hmmp->nr_v + 3 + MAX_GAP_SIZE) * MAX_GAP_SIZE/2 * hmmp->nr_v;
#ifdef DEBUG_RD
  printf("array_tail_size = %d\n", array_tail_size);
#endif
  to_trans_array = (struct path_element**)
    (malloc_or_die(array_head_size * sizeof(struct path_element*) +
		   (array_tail_size + hmmp->nr_v) * sizeof(struct path_element)));
  to_trans = (struct path_element*)(to_trans_array + hmmp->nr_v);
  hmmp->to_trans_array = to_trans_array;
  
  /* find all paths and add them to to_trans_array */
  for(v = 0; v < hmmp->nr_v; v++) /* from-vertex */ {
    *to_trans_array = to_trans;
    if(silent_vertex_multi(v, hmmp) == YES) {
      to_trans->vertex = END;
      to_trans->next = NULL;
      to_trans++;
      to_trans_array++;
      continue;
    }
    for(w = 0; w < hmmp->nr_v; w++) /* to-vertex */ {
      if(silent_vertex_multi(w,hmmp) == YES) {
	continue;
      }
      temp_path = (struct path_element*)(malloc_or_die(1000 * sizeof(struct path_element)));
      add_all_to_paths_multi(v, w, hmmp, &to_trans, temp_path, 0);
      free(temp_path);
    }
    to_trans->vertex = END;
    to_trans->next = NULL;
    to_trans++;
    to_trans_array++;
  }
 
 
#ifdef DEBUG_RD
  printf("array_head_size, array_tail_size = %d, %d\n", array_head_size, array_tail_size);
  dump_to_trans_array(hmmp->nr_v, hmmp->to_trans_array);
#endif 
}

void add_all_to_paths_multi(int v, int w, struct hmm_multi_s *hmmp,
		  struct path_element **to_transp, struct path_element *temp_pathp, int length)
{
  int i,j;
  int *xp;
  struct path_element p_el;
  
  if(length > MAX_GAP_SIZE) {
    return;
  }

  if(*(hmmp->transitions + get_mtx_index(v, w, hmmp->nr_v)) != 0.0) {
    /* direct path to w, add total path */
    for(i = 0; i < length; i++) {
      p_el.vertex = (temp_pathp + i)->vertex;
      p_el.next = (*to_transp) + 1;
      memcpy(*to_transp, &p_el, sizeof(struct path_element));
      (*to_transp)++;
    }
    p_el.vertex = w;
    p_el.next = NULL;
    memcpy(*to_transp, &p_el, sizeof(struct path_element));
    (*to_transp)++;
  }
  
  xp = *(hmmp->to_silent_trans_array + v);
  while(*xp != END) {
    (temp_pathp + length)->vertex = *xp;
    (temp_pathp + length)->next = NULL;
    add_all_to_paths_multi(*xp, w, hmmp, to_transp, temp_pathp, length + 1);
    xp++;
  }
  

}


int silent_vertex_multi(int k, struct hmm_multi_s *hmmp)
{
#ifdef DEBUG_RD
  printf("startnode = %d\n",hmmp->startnode);
  printf("endnode = %d\n",hmmp->endnode);
  dump_silent_vertices_multi(hmmp);
  dump_emiss_matrix(hmmp->nr_v, hmmp->a_size, hmmp->emissions);
#endif
  if(*(hmmp->emissions + get_mtx_index(k,0,hmmp->a_size)) == SILENT && k != hmmp->startnode && k != hmmp->endnode) {
    return YES;
  }
  else {
    return NO;
  }
}



/* check all probabilities and abort if some prob > 1.0 or < 0.0 */
void check_probs_multi(struct hmm_multi_s *hmmp)
{
  int i,j;
  double sum;
  double diff;
  double prob;
  /* transition probabilities first */
  for(i = 0; i < hmmp->nr_v; i++) {
    sum = 0;
    for(j = 0; j < hmmp->nr_v; j++) {
      prob = *((hmmp->transitions) + (i*hmmp->nr_v) + j);
      if(prob > 1.0 || prob < 0.0) {
	printf("Illegal probabilities (prob < 0.0 or prob > 1.0)\n");
	exit(-1);
      }
      else {
	sum += prob;
      }
    }
    diff = 1.0 - sum;
    /* maybe something about auto correcting the probabilities
     * will be implemented later */
  }

  /* then emission probabilities */
  for(i = 0; i < hmmp->nr_v; i++) {
    sum = 0;
    for(j = 0; j < hmmp->a_size; j++) {
      prob = *((hmmp->emissions) + (i*hmmp->a_size) + j);
      if((prob > 1.0 || prob < 0.0) && prob != SILENT) {
	printf("Illegal probabilities (prob < 0.0 or prob > 1.0)\n");
	exit(-1);
      }
      else {
	sum += prob;
      }
    }
    diff = 1.0 - sum;
    /* maybe something about auto correcting the probabilities
     * will be implemented later */
  }

#ifdef DEBUG_RD
  //dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->transitions);
  //dump_trans_matrix(hmmp->nr_v, hmmp->nr_v, hmmp->log_transitions);
  //dump_emiss_matrix(hmmp->nr_v, hmmp->a_size, hmmp->emissions);
  //dump_emiss_matrix(hmmp->nr_v, hmmp->a_size, hmmp->log_emissions);
#endif

}

int read_prior_files_multi(int nr_priorfiles, struct emission_dirichlet_s *emission_priorsp,
			   int a_size, FILE *file)
{
  int i,j,k;
  double q_value, alpha_value, alpha_sum, logbeta;
  char s[2048], *p;
  char ps[2048];
  char *file_name;
  char *pri;
  FILE *priorfile;
  
 
  /* find \n sign and remove it */
  if(fgets(s, 2048, file) != NULL) {
    p = s;
    while((p = strstr(p, "\n")) != NULL) {
      strncpy(p, " ", 1);
    }
  }
  
  /* read all before first filename */
  strtok(s," ");

 
  if((file_name = strtok(NULL, " ")) == NULL) {
    /* done */
    return 0;
  }
 
  for(i = 0; i < nr_priorfiles; i++) {
    /* open the priorfile */
    if((file_name = strtok(NULL, " ")) == NULL) {
      /* done */
      return 0;
    }
    
    else {
      if((priorfile = fopen(file_name,"r")) != NULL) {
	printf("Opened priorfile %s\n", file_name);
      }
      else {
	perror(file_name);
	return -1;
      }
    }

    /* put name i struct */
    strcpy(emission_priorsp->name, file_name);

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
    (emission_priorsp + i)->nr_components = atoi(&ps[0]);
    (emission_priorsp + i)->alphabet_size = a_size;
    
    
    /* allocate memory for arrays and matrix to this prior struct */
    (emission_priorsp + i)->q_values = malloc_or_die((emission_priorsp + i)->nr_components *
						     sizeof(double));
    (emission_priorsp + i)->alpha_sums = malloc_or_die((emission_priorsp + i)->nr_components *
						     sizeof(double));
    (emission_priorsp + i)->logbeta_values =
      malloc_or_die((emission_priorsp + i)->nr_components * sizeof(double));
    (emission_priorsp + i)->prior_values = malloc_or_die((emission_priorsp + i)->nr_components *
							 a_size * sizeof(double));

    for(j = 0; j < (emission_priorsp + i)->nr_components; j++) {
      /* put q-value in array, skip empty and comment lines  */
      while(1) {
	if(fgets(ps, 2048, priorfile) != NULL) {
	  if(*ps == '#' || *ps == '\n') {
	    
	  }
	  else {
	    break;
	  }
	}
	else {
	  printf("Prior file has incorrect format\n");
	}
      }
      q_value = atof(&ps[0]);
      *((emission_priorsp + i)->q_values + j) = q_value; 
#ifdef DEBUG_RDPRI
      printf("q_value = %f\n", *(((emission_priorsp + i)->q_values) + j));
#endif
      
      /* put alpha-values of this component in matrix */
      alpha_sum = 0.0;
      k = 0;
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
      
      pri = &ps[0];
      for(k = 0; k < a_size; k++) {
	alpha_value = strtod(pri, &pri);
	alpha_sum += alpha_value;
	*(((emission_priorsp + i)->prior_values) +
	  get_mtx_index(j, k,a_size)) = alpha_value;
      }
      
      /* put sum of alphavalues in array */
      *(((emission_priorsp + i)->alpha_sums) + j) = alpha_sum; 
      
      /* calculate logB(alpha) for this component, store in array*/
      logbeta = 0;
      for(k = 0; k < a_size; k++) {
	logbeta += lgamma(*(((emission_priorsp + i)->prior_values) +
			    get_mtx_index(j, k, a_size)));
	
#ifdef DEBUG_RDPRI
	printf("prior_value = %f\n", *(((emission_priorsp + i)->prior_values) +
				       get_mtx_index(j, k, a_size)));
	printf("lgamma_value = %f\n", lgamma(*(((emission_priorsp + i)->prior_values) +
					       get_mtx_index(j, k, a_size))));
#endif
      }
      logbeta = logbeta - lgamma(*(((emission_priorsp + i)->alpha_sums) + j));
      *(((emission_priorsp + i)->logbeta_values) + j) = logbeta;
    }
#ifdef DEBUG_RDPRI
    dump_prior_struct(emission_priorsp + i);
#endif
    
    /* some cleanup before continuing with next prior file */
    fclose(priorfile);
    emission_priorsp++;

  }
  return 0;
}

int read_trans_prior_files_multi(int nr_priorfiles, void *emission_priorsp,
		     struct hmm_multi_s *hmmp, FILE *file)
{
  char s[2048];
  
  /* not implemented yet */
  if(fgets(s, 2048, file) != NULL) {
  }
  return 0;
}

void create_tot_transitions_multi(struct hmm_multi_s *hmmp)
{
  int v,w;
  struct path_element *wp;
  double t_res;
  double log_t_res, cur_value;

  hmmp->tot_transitions = (double*)malloc_or_die(hmmp->nr_v * hmmp->nr_v * sizeof(double));
  hmmp->max_log_transitions = (double*)malloc_or_die(hmmp->nr_v * hmmp->nr_v * sizeof(double));
  init_float_mtx(hmmp->max_log_transitions, DEFAULT, hmmp->nr_v * hmmp->nr_v);

  for(v = 0; v < hmmp->nr_v; v++) {
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
      /* tot_transitions */
      *(hmmp->tot_transitions + get_mtx_index(w,v,hmmp->nr_v)) += t_res;
      
      /* max_log_transitions */
      if(t_res != 0.0) {
	log_t_res = log10(t_res);
	cur_value = *(hmmp->max_log_transitions + get_mtx_index(w,v,hmmp->nr_v));
	if(cur_value == DEFAULT || log_t_res > cur_value) {
	  *(hmmp->max_log_transitions + get_mtx_index(w,v,hmmp->nr_v)) = log_t_res;
	}
      }
      wp++;
    }
  }
}

void create_tot_trans_arrays_multi(struct hmm_multi_s *hmmp)
{
  int v,w;
  struct path_element *p_elp;
  int malloc_size;

  malloc_size = 0;
  for(v = 0; v < hmmp->nr_v; v++) {
    for(w = 0; w < hmmp->nr_v; w++) {
      if(*(hmmp->tot_transitions + get_mtx_index(v,w,hmmp->nr_v)) != 0.0) {
	malloc_size++;
      }
    }
  }
  
  hmmp->tot_to_trans_array = (struct path_element**)malloc_or_die(hmmp->nr_v * sizeof(struct path_element*) +
								  (malloc_size + hmmp->nr_v) * sizeof(struct path_element));
  
  hmmp->tot_from_trans_array = (struct path_element**)malloc_or_die(hmmp->nr_v * sizeof(struct path_element*) +
								    (malloc_size + hmmp->nr_v) * sizeof(struct path_element));

  /* fill in tot_to_trans_array */
  p_elp = (struct path_element*)(hmmp->tot_to_trans_array + hmmp->nr_v);
  for(v = 0; v < hmmp->nr_v; v++) {
    *(hmmp->tot_to_trans_array + v) = p_elp;
    for(w = 0; w < hmmp->nr_v; w++) {
      if(*(hmmp->tot_transitions + get_mtx_index(v,w,hmmp->nr_v)) != 0.0) {
	p_elp->vertex = w;
	p_elp->next = NULL;
	p_elp++;
      }
    }
    p_elp->vertex = END;
    p_elp->next = NULL;
    p_elp++;
  }
  

  /* fill in tot_from_trans_array */
  p_elp = (struct path_element*)(hmmp->tot_from_trans_array + hmmp->nr_v);
  for(v = 0; v < hmmp->nr_v; v++) {
    *(hmmp->tot_from_trans_array + v) = p_elp;
    for(w = 0; w < hmmp->nr_v; w++) {
      if(*(hmmp->tot_transitions + get_mtx_index(w,v,hmmp->nr_v)) != 0.0) {
	p_elp->vertex = w;
	p_elp->next = NULL;
	p_elp++;
      }
    }
    p_elp->vertex = END;
    p_elp->next = NULL;
    p_elp++;
  }
}
