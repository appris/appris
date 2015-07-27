/*
	kalign2_profile_alignment.h
	
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
int local_numseq;
int local_numprofiles;

int* assign_gap_codes(int* seq,int len);


int is_member(struct alignment* aln,int test);

float** dna_profile_distance(struct alignment* aln,float** dm,struct parameters* param, int nj);
float** protein_profile_wu_distance(struct alignment* aln,float** dm,struct parameters* param, int nj);


int** hirschberg_profile_alignment(struct alignment* aln,int* tree,float**submatrix, int** map);
float* make_profile_from_alignment(float* prof, int num,struct alignment* aln,float** subm);


