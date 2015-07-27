/*
	kalign2_output.h
	
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

#include <ctype.h>

void aln_output(struct alignment* aln,struct parameters* param);
void msf_output(struct alignment* aln,char* outfile);
void fasta_output(struct alignment* aln,char* outfile);
void clustal_output(struct alignment* aln,char* outfile);
void macsim_output(struct alignment* aln,char* outfile,char* infile);

struct names* get_meaningful_names(struct alignment* aln,int id);





