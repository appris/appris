/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifndef UPGMADATA_H
#define UPGMADATA_H

struct distMatRow
{
    distMatRow(double* p, int n)
    {
        ptrToDistMatRow = p;
        numDists = n;
    }
    double* ptrToDistMatRow;
    int numDists; 
};

const int BLANKDIST = -1;
const int INTERNAL = -1;

#endif
