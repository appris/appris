/*
	kalign2_output.c
	
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

#include "kalign2.h"
#include "kalign2_output.h"

void output(struct alignment* aln,struct parameters* param)
{
	if(!param->format){
		fasta_output(aln,param->outfile);
	}else{
		if (byg_start(param->format,"alnALNclustalCLUSTALclustalwCLUSTALWclustalWClustalW") != -1){
			aln_output(aln,param);
		}else if (byg_start(param->format,"msfMSFgcgGCGpileupPILEUP") != -1){
			msf_output(aln,param->outfile);
		}else if (byg_start(param->format,"eclu") != -1){
			clustal_output(aln,param->outfile);
		}else if (byg_start("macsim",param->format) != -1){
			macsim_output(aln,param->outfile,param->infile[0]);
		}else{
			fasta_output(aln,param->outfile);
		}
	}
	free_param(param);
}

void macsim_output(struct alignment* aln,char* outfile,char* infile)
{
	int i,j,f;
	int tmp;
	struct feature *fn = 0;
	FILE *fout = NULL;
	if(outfile){
		if ((fout = fopen(outfile, "w")) == NULL){
			fprintf(stderr,"can't open output\n");
			exit(0);
		}
	}else{
		fout = stdout; 
	}
	fprintf(fout,"<?xml version=\"1.0\"?>\n<!DOCTYPE macsim SYSTEM \"http://www-bio3d-igbmc.u-strasbg.fr/macsim.dtd\">\n<macsim>\n<alignment>\n<aln-name>");
	if(infile){
		fprintf(fout,"%s.kalign",infile);
	}else{
		fprintf(fout,"kalign alignment");
	}
	fprintf(fout,"</aln-name>\n");

	for (i =0;i< numseq;i++){
		//c = aln->sl[i];
		f = aln->nsip[i];
		
		fprintf(fout,"<sequence seq-type=\"Protein\">\n");
		fprintf(fout,"<seq-name>");
		for (j =0; j < aln->lsn[f];j++){
			if(!iscntrl((int)aln->sn[f][j])){
				fprintf(fout,"%c",aln->sn[f][j]);
			}
		}
		fprintf(fout,"</seq-name>");
		fprintf(fout,"<seq-info>\n");
		fprintf(fout,"<accession>1aab_</accession>\n");
		fprintf(fout,"<nid>1aab_</nid>\n");
		fprintf(fout,"<ec>0.0.0.0</ec>\n");
		fprintf(fout,"<group>0</group>\n");
		if(aln->ft){
		if(aln->ft[f]){
			
			fprintf(fout,"<ftable>\n");
			fn = aln->ft[f];
			while(fn){
				fprintf(fout,"<fitem><ftype>%s</ftype><fstart>%d</fstart><fstop>%d</fstop><fnote>%s</fnote></fitem>\n",fn->type,fn->start,fn->end,fn->note);
				fn = fn->next;
			}
			fprintf(fout,"</ftable>\n</seq-info>\n");
		}
		}
		fprintf(fout,"<seq-data>\n");

		for (j = 0; j < aln->sl[f];j++){
			tmp = aln->s[f][j];
			while (tmp){
				fprintf(fout,"-");
				tmp--;
			}
			fprintf(fout,"%c",aln->seq[f][j]);
		}
		tmp =aln->s[f][aln->sl[f]];
		while (tmp){
			fprintf(fout,"-");
			tmp--;
		}
		fprintf(fout,"\n");
		fprintf(fout,"</seq-data>\n");
		fprintf(fout,"</sequence>\n");
	}
	fprintf(fout,"</alignment>\n");
	fprintf(fout,"</macsim>\n");
	if(outfile){
		fclose(fout);
	}
	free_aln(aln);
}


void msf_output(struct alignment* aln,char* outfile)
{
	int i,j,c,f,g;
	int max = 0;
	int aln_len = 0;
	int tmp;
	char** linear_seq = 0;
	FILE *fout = NULL;
	
	linear_seq = malloc(sizeof(char*)*numseq);
	
	aln_len = 0;
	for (j = 0; j <= aln->sl[0];j++){
		aln_len+= aln->s[0][j];
	}
	aln_len += aln->sl[0];
	
	for (i =0;i< numseq;i++){
		linear_seq[i] = malloc(sizeof(char)*(aln_len+1));
		
		c = 0;
		for (j = 0; j < aln->sl[i];j++){
			tmp = aln->s[i][j];
			while (tmp){
				linear_seq[i][c] ='-';
				c++;
				tmp--;
			}
			linear_seq[i][c] = aln->seq[i][j];
			c++;
		}
		
		tmp =aln->s[i][aln->sl[i]];
		while (tmp){
			linear_seq[i][c] ='-';
			c++;
			tmp--;
		}		
		linear_seq[i][c] = 0;
	}

	if(outfile){
		if ((fout = fopen(outfile, "w")) == NULL){
			fprintf(stderr,"can't open output\n");
			exit(0);
		}
	}else{
		fout= stdout;
	}
	fprintf(fout,"PileUp\n\n\n\n   MSF:   %d  Type: P    Check:  7038   ..\n\n",aln_len);

	for (j = 0; j< numseq;j++){
		if( aln->lsn[j] > max){
			max = aln->lsn[j];
		}
	}

	for (i = 0; i< numseq;i++){
			f = aln->nsip[i];
			fprintf(fout," Name: ");
			for (c = 0; c < aln->lsn[f];c++){
				if(!iscntrl((int)aln->sn[f][c])){
					fprintf(fout,"%c",aln->sn[f][c]);
				}
			}
			while(c < max+3){
				fprintf(fout," ");
				c++;
			}
			fprintf(fout,"Len:   ");
			fprintf(fout,"%d",aln_len);
			fprintf(fout,"  Check:  2349  Weight:  1.00\n");
			
	}
	fprintf(fout,"\n\n//\n\n");

	for (i = 0; i+60 < aln_len;i +=60){
		for (j = 0; j< numseq;j++){
			f = aln->nsip[j];
			for (c = 0; c < aln->lsn[f];c++){
				if(!iscntrl((int)aln->sn[f][c])){
					fprintf(fout,"%c",aln->sn[f][c]);
				}
			}
			while(c < max+3){
				fprintf(fout," ");
				c++;
			}
			g = 1;
			for (c = 0; c < 60;c++){
				fprintf(fout,"%c",linear_seq[f][c+i]);
				if (g == 10){
					fprintf(fout," ");
					g = 0;
				}
				g++;
			}
			fprintf(fout,"\n");
			
		}
		fprintf(fout,"\n\n");
	}
	for (j = 0; j< numseq;j++){
		f = aln->nsip[j];
		
		for (c = 0; c< aln->lsn[f];c++){
			if(!iscntrl((int)aln->sn[f][c])){
				fprintf(fout,"%c",aln->sn[f][c]);
			}
		}
		
		while(c < max+3){
			fprintf(fout," ");
			c++;
		}
		
		g = 1;
		for (c = i; c< aln_len;c++){
			fprintf(fout,"%c",linear_seq[f][c]);
			if (g == 10){
				fprintf(fout," ");
				g = 0;
			}
			g++;
		}
		fprintf(fout,"\n");
	
	}
	fprintf(fout,"\n\n");
	if(outfile){
		fclose(fout);
	}
	
	for (i =0;i< numseq;i++){
		free(linear_seq[i]);
	}
	
	free(linear_seq);
	free_aln(aln);
}


void clustal_output(struct alignment* aln,char* outfile)
{
	int i,j,c,f;
	int tmp;
	int aln_len = 0;
	char** linear_seq = 0;
	
	FILE* fout = NULL;
	
	linear_seq = malloc(sizeof(char*)*numseq);

	aln_len = 0;
	
	for (j = 0; j <= aln->sl[0];j++){
		aln_len+= aln->s[0][j];
	}
	
	aln_len += aln->sl[0];
	
	for (i =0;i< numseq;i++){
		linear_seq[i] = malloc(sizeof(char)*(aln_len+1));
		
		c = 0;
		for (j = 0; j < aln->sl[i];j++){
			tmp = aln->s[i][j];
			while (tmp){
				linear_seq[i][c] ='-';
				c++;
				tmp--;
			}
			linear_seq[i][c] = aln->seq[i][j];
			c++;
		}
		
		tmp =aln->s[i][aln->sl[i]];
		while (tmp){
			linear_seq[i][c] ='-';
			c++;
			tmp--;
		}		
		linear_seq[i][c] = 0;
	}


	if(outfile){
		if ((fout = fopen(outfile, "w")) == NULL){
			fprintf(stderr,"can't open output\n");
			exit(0);
		}
	}else{
		fout = stdout;
	}

	fprintf(fout,"Kalign (2.0) alignment in ClustalW format\n\n\n");


	for (i = 0; i+60 < aln_len;i +=60){
		for (j = 0; j< numseq;j++){
			f = aln->nsip[j];
			for (c = 0; c < aln->lsn[f];c++){
				if(!iscntrl((int)aln->sn[f][c])){
					fprintf(fout,"%c",aln->sn[f][c]);
				}
			}
			while(c < 18){
				fprintf(fout," ");
				c++;
			}
			
			for (c = 0; c < 60;c++){
				fprintf(fout,"%c",linear_seq[f][c+i]);
			}
			fprintf(fout,"\n");
		}
		fprintf(fout,"\n\n");
	}
	for (j = 0; j< numseq;j++){
		f = aln->nsip[j];
		for (c = 0; c< aln->lsn[f];c++){
			if(!iscntrl((int)aln->sn[f][c])){
				fprintf(fout,"%c",aln->sn[f][c]);
			}
		}
		while(c < 18){
			fprintf(fout," ");
			c++;
		}
	
		for (c = i; c< aln_len;c++){
			fprintf(fout,"%c",linear_seq[f][c]);
		}
		fprintf(fout,"\n");
	}
	fprintf(fout,"\n\n");
	if(outfile){
		fclose(fout);
	}
	for (i =0;i< numseq;i++){
		free(linear_seq[i]);
	}
	free(linear_seq);
	free_aln(aln);
}

void aln_output(struct alignment* aln,struct parameters* param)
{
	char* outfile = param->outfile;
	int i,j,c,f;
	int tmp;
	int aln_len = 0;
	
	//int namestart = 0;
	int max_name_len = 20;
	int tmp_len = 0;
	char** linear_seq = 0;
	
	struct names* n;
	
	n = get_meaningful_names(aln,param->id);
	
	//namestart = get_meaningful_names(aln,param->id);
	
	c = -1;
	for (i = 0; i< numseq;i++){
		if(n->len[i] > c){
			c = n->len[i];
		}
		/*f = 0;
		for (j = namestart;j < aln->lsn[i];j++){
			if(isspace((int)aln->sn[i][j])){
				break;
			}
			f++;
		}
		if(f > c){
			c = f;
		}
		}*/
	}
	
	if(c < max_name_len){
		max_name_len = c;//this is know the maximum length of a unique name isdjgbv skj
	}
		
	FILE* fout = NULL;
	
	linear_seq = malloc(sizeof(char*)*numseq);

	aln_len = 0;
	for (j = 0; j <= aln->sl[0];j++){
		aln_len+= aln->s[0][j];
	}
	aln_len += aln->sl[0];
	
	for (i =0;i< numseq;i++){
		linear_seq[i] = malloc(sizeof(char)*(aln_len+1));
		
		c = 0;
		for (j = 0; j < aln->sl[i];j++){
			tmp = aln->s[i][j];
			while (tmp){
				linear_seq[i][c] ='-';
				c++;
				tmp--;
			}
			linear_seq[i][c] = aln->seq[i][j];
			c++;
		}
		
		tmp =aln->s[i][aln->sl[i]];
		while (tmp){
			linear_seq[i][c] ='-';
			c++;
			tmp--;
		}		
		linear_seq[i][c] = 0;
	}
	
	if(outfile){
		if ((fout = fopen(outfile, "w")) == NULL){
			fprintf(stderr,"can't open output\n");
			exit(0);
		}
	}else{
		fout = stdout;
	}
	
	fprintf(fout,"Kalign (2.0) alignment in ClustalW format\n\n\n");

	for (i = 0; i+60 < aln_len;i +=60){
		for (j = 0; j< numseq;j++){
			f = aln->nsip[j];
			tmp_len = (max_name_len < n->len[f]) ? max_name_len:n->len[f];
			for (c = 0; c < tmp_len;c++){
				if(isspace((int)aln->sn[f][c+n->start[f]])){
					break;
				}
				
				if(!iscntrl((int)aln->sn[f][c+n->start[f]])){
					fprintf(fout,"%c",aln->sn[f][c+n->start[f]]);
				}
			}
			
			while(c < max_name_len+5){
				fprintf(fout," ");
				c++;
			}
			
			for (c = 0; c < 60;c++){
				fprintf(fout,"%c",linear_seq[f][c+i]);
			}
			fprintf(fout,"\n");
		}
		fprintf(fout,"\n\n");
	}
	
	for (j = 0; j< numseq;j++){
		f = aln->nsip[j];
		tmp_len = (max_name_len < n->len[f]) ? max_name_len:n->len[f];
		for (c = 0; c< tmp_len;c++){
			if(isspace((int)aln->sn[f][c+n->start[f]])){
				break;
			}
			
			if(!iscntrl((int)aln->sn[f][c+n->start[f]])){
				fprintf(fout,"%c",aln->sn[f][c+n->start[f]]);
			}
		}
		
		while(c < max_name_len + 5){
			fprintf(fout," ");
			c++;
		}
	
		for (c = i; c < aln_len;c++){
			fprintf(fout,"%c",linear_seq[f][c]);
		}
		fprintf(fout,"\n");
	}
	fprintf(fout,"\n\n");
	if(outfile){
		fclose(fout);
	}
	
	names_free(n);
	
	for (i =0;i< numseq;i++){
		free(linear_seq[i]);
	}
	free(linear_seq);
	free_aln(aln);
}

struct names* get_meaningful_names(struct alignment* aln,int id)
{

	struct names* n = 0;
	int i,j,c;
	int min_len = 0;
	int start = 0;
	int globalstart = 1000000;
	
	n = names_alloc(n);
	for (i = 0; i < numseq;i++){
 		n->end[i] = aln->lsn[i];
	}
	
		
	if (id == -1){
		for(i =0; i < numseq-1;i++){
			for (j = i+1; j < numseq;j++){
				min_len = (aln->lsn[i] < aln->lsn[j])? aln->lsn[i] : aln->lsn[j];
				start = 0;
				for (c = 0; c < min_len;c++){
					if(isalnum((int)aln->sn[i][c]) && isalnum((int)aln->sn[j][c])){
						if( aln->sn[i][c] != aln->sn[j][c]){
							break;
						}
					}else{
						if(aln->sn[i][c] == aln->sn[j][c]){
							if(aln->sn[i][c] != '_' && aln->sn[i][c] != '-'){
								start = c+1;
							}
						}else{
							break;
						}
					}
				}
					
				//fprintf(stderr,"%s\n%s\nstart: %d\n\n",aln->sn[i],aln->sn[j],start);
				
				if (start < globalstart){
					globalstart = start;
				}
			} 
		}
		for (i = 0; i < numseq;i++){
			n->start[i] = globalstart;
			for (j = n->start[i]; j < aln->lsn[i];j++){
				if(!isalnum((int)aln->sn[i][j]) && aln->sn[i][j] != '_' && aln->sn[i][j] != '-'){
					n->end[i] = j;
					break;
				}
			}
		}

	}else{
		for(i =0; i < numseq;i++){
			start = 0;
			min_len = 0;
			for (j = 0; j < aln->lsn[i];j++){
				if((isalnum((int)aln->sn[i][j]) || aln->sn[i][j] == '_' || aln->sn[i][j] == '-')&& start == 0 ){
					n->start[i] = j;
					min_len++;
					start = 1;
				}else if ((!isalnum((int)aln->sn[i][j]) && aln->sn[i][j] != '_' && aln->sn[i][j] != '-')&& start == 1) {
					if(id == min_len){
						n->end[i] = j;
						break;
					}
					start = 0;
			
				}
			}
			if(id > min_len){
				fprintf(stderr,"Warning: sequence %d has no %dth word in the identifier line:\n%s\n",i,id,aln->sn[i]);
				n->start[i] = 0;
			}
		}
	}	
	
	for (i = 0; i < numseq;i++){
		//fprintf(stderr,"%s\n%d-%d\n",aln->sn[i],n->start[i],n->end[i]);
 		n->len[i] = n->end[i] - n->start[i];
	}
	
	return n;
}


void fasta_output(struct alignment* aln,char* outfile)
{
	int i,j,c,f;
	int tmp;
	FILE *fout = NULL;
	if(outfile){
		if ((fout = fopen(outfile, "w")) == NULL){
			fprintf(stderr,"can't open output\n");
			exit(0);
		}
	}else{
		fout = stdout;
	}
	for (i = 0; i < numseq;i++){
		f = aln->nsip[i];
		fprintf(fout,">%s\n",aln->sn[f]);
		c = 0;
		for (j = 0; j < aln->sl[f];j++){
			tmp = aln->s[f][j];
			while (tmp){
				fprintf(fout,"-");
				c++;
				if(c == 60 && j != aln->sl[f]-1){
					fprintf(fout,"\n");
					c = 0;
				}
				tmp--;
			}
			fprintf(fout,"%c",aln->seq[f][j]);
			c++;
			if(c == 60 && j != aln->sl[f]-1){
				fprintf(fout,"\n");
				c = 0;
			}
		}
		tmp = aln->s[f][aln->sl[f]];
		while (tmp){
			fprintf(fout,"-");
			c++;
 			if(c == 60 && j != aln->sl[f]-1){
				fprintf(fout,"\n");
				c = 0;
			}
			tmp--;
		}
		fprintf(fout,"\n");
	}
	if(outfile){
		fclose(fout);
	}
	free_aln(aln);
}

