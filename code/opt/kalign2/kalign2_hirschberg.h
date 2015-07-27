/*
	kalign2_hirschberg.h
	
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

int* hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm, int* hirsch_path);
struct states* foward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm);
struct states* backward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm);
int* hirsch_align_two_pp_vector(const float* prof1,const float* prof2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);

struct states* backward_hirsch_pp_dynm(const float* prof1,const float* prof2,struct hirsch_mem* hm);

struct states* backward_hirsch_pp_dyn0(const float* prof1,const float* prof2,struct hirsch_mem* hm);

int* hirsch_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm, int* hirsch_path,int sip);
struct states* foward_hirsch_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm,int sip);
struct states* backward_hirsch_ps_dyn(const float* prof1,const int* seq2,struct hirsch_mem* hm,int sip);
int* hirsch_align_two_ps_vector(const float* prof1,const int* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[],int sip);


int* hirsch_ss_dyn(float**subm, const int* seq1,const int* seq2,struct hirsch_mem* hm, int* hirsch_path);
struct states* foward_hirsch_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm);
struct states* backward_hirsch_ss_dyn(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm);
int* hirsch_align_two_ss_vector(float**subm,const int* seq1,const int* seq2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);
