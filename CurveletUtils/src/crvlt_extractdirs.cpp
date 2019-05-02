/*
 *  crvlt_extractdirs.cpp
 *
 *  Created by Tobias Gebaeck on 12.9.2007.
 *  Copyright 2007 ETH Zurich. All rights reserved.
 *
 */

#include "mex.h"
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
	if (c.size() == 0)
		mexErrMsgTxt("input 1 must be a valid Curvelet data structre");
  
  IntNumVec levs;
  mex2cpp(prhs[1], levs);
  for (int k=0; k<levs.m(); k++)   // go from matlab's 1-based index to 0-based
    levs(k)--;

	int nrflds = 2;
	DirFieldOptions opts = {1, 1, 3, CF_Abs};  // Default values
	switch (nrhs) {
		case 7:
			if (!mxIsEmpty(prhs[6]))
				mex2cpp(prhs[6], (int&)opts.FilterType);
		case 6:
			if (!mxIsEmpty(prhs[5]))
				mex2cpp(prhs[5], opts.LocMaxGap);
		case 5:
			if (!mxIsEmpty(prhs[4]))
				mex2cpp(prhs[4], opts.CurveletSize);
		case 4:
			if (!mxIsEmpty(prhs[3]))
				mex2cpp(prhs[3], opts.CircSumOffset);
		case 3:
			if (!mxIsEmpty(prhs[2]))
				mex2cpp(prhs[2], nrflds);
			break;
	}
	if (nrflds <= 0) mexErrMsgTxt("Bad value for number of fields!");
	if (opts.LocMaxGap < 0) mexErrMsgTxt("Bad value for LocMaxGap!");
	if (opts.CurveletSize < 0) mexErrMsgTxt("Bad value for CurveletSize!");
	//if (opts.CircSumOffset < 0) mexErrMsgTxt("Bad value for CircSumOffset!");
	
  int m, n; 
  GetFinestSizes(c, levs, m, n);   // get sizes from finest specified level
  CurvUtilError cerr;
	
	// Extract fields
	vector< vector<DblNumMat> > field;
	if (opts.CircSumOffset < 0) {
		// Get magnitudes from all directions, but best direction using 1 nearest
		opts.CircSumOffset = -opts.CircSumOffset;
		ExtractMultipleDirFields(m, n, c, levs, nrflds, field, &opts);
		DblNumMat mags;
		GetMagnitude(m, n, c, levs, mags, opts.CurveletSize);
		
		// Recompute magnitudes for field
		for (int fn=0; fn<nrflds; fn++) {
			for (int xi=0; xi<m; xi++) {
				for (int yi=0; yi<n; yi++) {
					double fact = mags(xi,yi) / sqrt(field[fn][0](xi,yi) * field[fn][0](xi,yi) + field[fn][1](xi,yi) * field[fn][1](xi,yi));
					field[fn][0](xi,yi) *= fact;
					field[fn][1](xi,yi) *= fact;
				}
			}
		}
	}
	else
		cerr = ExtractMultipleDirFields(m, n, c, levs, nrflds, field, &opts);	
	
	if (cerr != CUError_NoError)
		mexErrMsgTxt("Error in ExtractMultipleDirFields!");
	
	cpp2mex(field, plhs[0]);
	
  return;
}
