/*
 *  Copyright 2003-2007 Tim Massingham (tim.massingham@ebi.ac.uk)
 *  Funded by EMBL - European Bioinformatics Institute
 */
/*
 *  This file is part of SLR ("Sitewise Likelihood Ratio")
 *
 *  SLR is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SLR is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with SLR.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>

void warnx ( const char * fmt, ...){
	va_list args;

	assert(NULL!=fmt);
	va_start(args,fmt);

	fputs("Warning: ",stderr);
	vfprintf(stderr,fmt,args);
	
	va_end(args);
}


void err ( const int errnum, const char * fmt, ...){
        va_list args;

	assert(NULL!=fmt);
        va_start(args,fmt);

        fprintf(stderr, "Error: ");
        vfprintf(stderr,fmt,args);
        
        va_end(args);
	abort();
}


