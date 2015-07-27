/*
	kalign2_input.h 
	
	Released under GPL - see the 'COPYING' file   
	
	Copyright (C) 2006 Timo Lassmann <timolassmann@gmail.com>
	
	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
    
	Please send bug reports, comments etc. to:
	timolassmann@gmail.com
*/


#include <unistd.h>
#include <string.h>
#include <ctype.h>

#define SEEK_START 0
#define SEEK_END 2

struct alignment* read_sequences(struct alignment* aln,char* string);
struct alignment* read_sequences_from_swissprot(struct alignment* aln,char* string);
struct alignment* read_sequences_uniprot_xml(struct alignment* aln,char* string);
struct alignment* read_sequences_macsim_xml(struct alignment* aln,char* string);
struct feature* read_ft(struct feature* ft,char* p);
struct alignment* read_sequences_clustal(struct alignment* aln,char* string);
struct alignment* read_sequences_stockholm(struct alignment* aln,char* string);

struct alignment* read_alignment(struct alignment* aln,char* string);
struct alignment* read_alignment_from_swissprot(struct alignment* aln,char* string);
struct alignment* read_alignment_uniprot_xml(struct alignment* aln,char* string);
struct alignment* read_alignment_macsim_xml(struct alignment* aln,char* string);
struct feature* read_ft(struct feature* ft,char* p);
struct alignment* read_alignment_clustal(struct alignment* aln,char* string);
struct alignment* read_alignment_stockholm(struct alignment* aln,char* string);

char* get_input_into_string(char* string,char* infile);



int count_sequences_macsim(char* string);
int count_sequences_swissprot(char* string);
int count_sequences_uniprot(char* string);
int count_sequences_stockholm(char* string);
int count_sequences_clustalw(char* string);
int count_sequences_fasta(char* string);

static char  usage[] = "\n\
        Usage: kalign2   [INFILE] [OUTFILE] [OPTIONS]\n\
        \n\
	Options:\n\n\
        -s,	-gapopen          Gap open penalty\n\
        	-gap_open\n\
        	-gpo\n\
        	\n\
        -e,	-gapextension     Gap extension penalty\n\
        	-gap_ext\n\
        	-gpe\n\
        	\n\
        -t,	-terminal_gap_extension_penalty	Terminal gap penalties\n\
        	-tgpe\n\
        	\n\
        -m,	-matrix_bonus     A constant added to the substitution matrix.\n\
        	-bonus\n\
        	\n\
        -c,	-sort            The order in which the sequences appear in the output alignment.\n\
		                   <input, tree, gaps.>\n\
		\n\
        -g,	-feature          Selects feature mode and specifies which features are to be used:\n\
		                   e.g. all, maxplp, STRUCT, PFAM-A....\n\
           	-same_feature_score          Score for aligning same features\n\
		-diff_feature_score          Penalty for aligning different features\n\
        	\n\
        -d,	-distance         Distance method.\n\
		                   <wu,pair>\n\
		\n\
        -b,	-guide-tree       Guide tree method.\n\
		-tree             <nj,upgma>\n\
		\n\
	-z,	-zcutoff         Parameter used in the wu-manber based distance calculation\n\
		\n\
        -i,	-input            The input file.\n\
        	-infile\n\
        	-in\n\
        	\n\
        -o,	-output           The output file.\n\
        	-outfile\n\
        	-out\n\
        	\n\
        -a,	-gap_inc           Parameter increases gap penalties depending on the number of existing gaps\n\
        	\n\
        -f,	-format           The output format:\n\
		                   <fasta, msf, aln, clu, macsim>\n\
		\n\
	-q,	-quiet            Print nothing to STDERR.\n\
		                  Read nothing from STDIN\n\
	\n\
	Examples:\n\n\
	Using pipes:\n\
		kalign2 [OPTIONS] < [INFILE]   > [OUTFILE]\n\
		more [INFILE] |  kalign2 [OPTIONS] > [OUTFILE]\n\
         \n\
	Relaxed gap penalties:\n\
		kalign2 -gpo 60 -gpe 9 -tgpe 0 -bonus 0 < [INFILE]   > [OUTFILE]\n\
         \n\
        Feature alignment with pairwise alignment based distance method and NJ guide tree:\n\
        	kalign2 -in test.xml -distance pair -tree nj -sort gaps -feature STRUCT -format macsim -out test.macsim\n\
        ";


