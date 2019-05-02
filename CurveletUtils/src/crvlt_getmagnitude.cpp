/*
 *  crvlt_getmagnitude.cpp
 *  
 *
 *  Created by Tobias Gebäck on 1.2.2008.
 *  Copyright 2008 ETH Zurich. All rights reserved.
 *
 */

#include "cu_mexaux.h"

#include "CurveletUtils.h"

using namespace std;
using namespace fdct_wrapping_ns;
using namespace CurveletUtils;


void GetFinestSizes(const CurveletData &c, const IntNumVec &levs, int &m, int &n)
{
  int mmax=0, nmax=0;
  for (int j=0; j<levs.m(); j++) {
    int lev = levs(j); 
    for (int k=0; k<c[lev].size(); k++) {
      if (c[lev][k].m() > mmax) mmax = c[lev][k].m();
      if (c[lev][k].n() > nmax) nmax = c[lev][k].n();
    }
  }  
  
  m=mmax;
  n=nmax;
}



void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  if (nrhs < 2)
    mexErrMsgTxt("At least 2 inputs required");
  
  if (!mxIsCell(prhs[0]))
    mexErrMsgTxt("input 1 must be a Curvelet data structure");
  CurveletData c;
  mex2cpp(prhs[0], c);
  
  IntNumVec levs;
  mex2cpp(prhs[1], levs);
  for (int k=0; k<levs.m(); k++)   // go from matlab's 1-based index to 0-based
    levs(k)--;
	
	int csize = 1;
	switch (nrhs) {
		case 3:
			if (!mxIsEmpty(prhs[2]))
				mex2cpp(prhs[2], csize);
			break;
	}
	if (csize < 0) mexErrMsgTxt("Bad value for CurveletSize!");
	
  int m, n; 
  GetFinestSizes(c, levs, m, n);   // get sizes from finest specified level
  
	// Extract fields
	DblNumMat mags;
	GetMagnitude(m, n, c, levs, mags, csize);
	
	cpp2mex(mags, plhs[0]);
	
  return;
}
