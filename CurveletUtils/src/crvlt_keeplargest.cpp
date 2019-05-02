/*
 *  crvlt_keeplargest.cpp
 *
 *  Created by Tobias Gebaeck on 17.10.2007.
 *  Copyright 2007 ETH Zurich. All rights reserved.
 *
 */

#include "mex.h"
#include "cu_mexaux.h"

#include "CurveletUtils.h"

using namespace std;
using namespace fdct_wrapping_ns;
using namespace CurveletUtils;


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  if (nrhs != 2)
    mexErrMsgTxt("2 inputs required");
  
	if (nlhs < 1 || nlhs > 2)
		mexErrMsgTxt("1 or 2 outputs required");
	
  if (!mxIsCell(prhs[0]))
    mexErrMsgTxt("input 1 must be a Curvelet data structure");
  CurveletData c;
  mex2cpp(prhs[0], c);
  
  int Ncoeffs;
	double percentage;
  mex2cpp(prhs[1], percentage);  // Get double
	if (percentage >= 1.0)   // interpret as number of coefficients
		Ncoeffs = (int) percentage;
	else {    // interpret as fraction of coefficients
		int nel = countnnz(c);
		Ncoeffs = (int) (nel * percentage);
	}
	
	double thrsh = FindNthCoeffSize(c, Ncoeffs);
	ThresholdCoeffs(c, thrsh);
	
	cpp2mex(c, plhs[0]);
	if (nlhs == 2) 
		cpp2mex(thrsh, plhs[1]);
  
  return;
}


