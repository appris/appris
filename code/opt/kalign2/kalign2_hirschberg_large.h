/*
	kalign2_hirschberg_large.h
	
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

struct large_states{
	float a;
	float ga;
	float gb;
	float x;
};



struct hirsch_large_mem{
	struct large_states* f;
	struct large_states* b;
	int starta;
	int startb;
	int enda;
	int endb;
	int size;
	int len_a;
	int len_b;
};

int* hirsch_large_pp_dyn(const float* prof1,const float* prof2,struct hirsch_large_mem* hm, int* hirsch_path);
struct large_states* foward_large_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_large_mem* hm);
struct large_states* backward_large_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_large_mem* hm);
int* hirsch_large_align_two_pp_vector(const float* prof1,const float* prof2,struct hirsch_large_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);

int* hirsch_large_ps_dyn(const float* prof1,const int* seq2,struct hirsch_large_mem* hm, int* hirsch_path,int sip);
struct large_states* foward_large_hirsch_ps_dyn(const float* prof1,const int* seq2,struct hirsch_large_mem* hm,int sip);
struct large_states* backward_large_hirsch_ps_dyn(const float* prof1,const int* seq2,struct hirsch_large_mem* hm,int sip);
int* hirsch_large_align_two_ps_vector(const float* prof1,const int* seq2,struct hirsch_large_mem* hm,int* hirsch_path,float input_states[],int old_cor[],int sip);


int* hirsch_large_ss_dyn(float**subm, const int* seq1,const int* seq2,struct hirsch_large_mem* hm, int* hirsch_path);
struct large_states* foward_large_hirsch_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_large_mem* hm);
struct large_states* backward_large_hirsch_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_large_mem* hm);
int* hirsch_large_align_two_ss_vector(float**subm,const int* seq1,const int* seq2,struct hirsch_large_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);

float* make_large_profile(float* prof, int* seq,int len,float** subm);
void set_large_gap_penalties(float* prof,int len,int nsip);
float* large_update(float* profa,float* profb,float* newp,int* path,int sipa,int sipb);

struct hirsch_large_mem* hirsch_large_mem_alloc(struct hirsch_large_mem* hm,int x);
struct hirsch_large_mem* hirsch_large_mem_realloc(struct hirsch_large_mem* hm,int x);
void hirsch_large_mem_free(struct hirsch_large_mem* hm);


