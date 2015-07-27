/*
	kalign2_string_matching.c
	
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

#include <string.h>

int byg_detect(int* text,int n)
{
	int Tc;
	int i  = 0;
	int s = 0;
	int T[256];
	for (i = 0;i < 256;i++){ 
		T[i] = 0; 
	}
	int mb = 1;
	//char *unique_aa = "EFILPQXZ";//permissiv
	//ABCDEFGHIJKLMNOPQRSTUVWXYZ
	char *unique_aa = "BDEFHIJKLMNOPQRSVWYZ";//restrictive
	int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
	for (i= 0;i < 20;i++){
		T[(int)aacode[unique_aa[i]-65]] |= 1;
	}
	for (i = 0;i < n;i++){
	//	fprintf(stderr,"%d\n",text[i]);
		if(text[i] != -1){
			s <<= 1;
			s |= 1;
			Tc = T[text[i]];
			s &= Tc;
			if(s & mb){
				return 0;
			}
		}
	}
	return 1;
}

int check_identity(char* n,char*m)
{
	int len_n;
	int len_m;
	int i;
	
	len_n = strlen(n);
	len_m = strlen(m);
	if(len_m != len_n){
		return -1;
	}
	for (i = 0; i < len_n;i++){
		if(n[i] != m[i]){
			return -1;
		}
	}
	return 1;
	
}


int byg_count(char* pattern,char*text)
{
	int Tc;
	int count = 0;
	int i  = 0;
	int s = 0;
	int T[256];
	for (i = 0;i < 256;i++){ 
		T[i] = 0; 
	}
	
	int m = strlen(pattern);
	int n = strlen (text);
	int mb = (1 << (m-1));
	
	for (i= 0;i < m;i++){
		T[(int)pattern[i]] |= (1 << i);
	}

	for (i = 0;i < n;i++){
		s <<= 1;
		s |= 1;
		Tc = T[(int)text[i]];
		s &= Tc;
		if(s & mb){
			count++;
		}
	}
	return count;
}

int byg_end(char* pattern,char*text)
{
	int Tc;
	int i  = 0;
	int s = 0;
	int T[256];
	for (i = 0;i < 256;i++){ 
		T[i] = 0; 
	}
	
	int m = strlen(pattern);
	int n = strlen (text);
	int mb = (1 << (m-1));

	for (i= 0;i < m;i++){
		T[(int)pattern[i]] |= (1 << i);
	}

	for (i = 0;i < n;i++){
		s <<= 1;
		s |= 1;
		if(!text[i]){
			return -1;
		}
		Tc = T[(int)text[i]];
		s &= Tc;
		if(s & mb){
			return i+1;
		}
	}
	return -1;
}

int byg_start(char* pattern,char*text)
{
	int Tc;
	int i  = 0;
	int s = 0;
	int T[256];
	for (i = 0;i < 256;i++){ 
		T[i] = 0; 
	}
	
	int m = strlen(pattern);
	int n = strlen(text);
	int mb = (1 << (m-1));
	
	for (i= 0;i < m;i++){
		T[(int)pattern[i]] |= (1 << i);
	}

	for (i = 0;i < n;i++){
		s <<= 1;
		s |= 1;
		Tc = T[(int)text[i]];
		s &= Tc;
		if(s & mb){
			return i-m+1;
		}
	}
	return -1;
}

