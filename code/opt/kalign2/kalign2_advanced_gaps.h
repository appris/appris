/*
	kalign2_advanced_gaps.h
	
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

void advanced_smooth_gaps(float* prof,int len,int window,float strength);
float* advanced_make_profile(float* prof, int* seq,int len,float** subm);
float* advanced_update(const float* profa,const float* profb,float* newp,int* path,int sipa,int sipb,float internal_gap_weight);


int* advanced_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm, int* hirsch_path);
int* advanced_hirsch_align_two_pp_vector(const float* prof1,const float* prof2,struct hirsch_mem* hm,int* hirsch_path,float input_states[],int old_cor[]);
struct states* advanced_foward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm);
struct states* advanced_backward_hirsch_pp_dyn(const float* prof1,const float* prof2,struct hirsch_mem* hm);
