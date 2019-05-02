/*
 *  cu_thresholding.cpp
 *
 *  Created by Tobias Gebäck on 17.10.2007.
 *  Copyright 2007 ETH Zurich. All rights reserved.
 *
 */

#include "CurveletUtils.h"
#include <algorithm>
using namespace std;

namespace CurveletUtils {
	using namespace fdct_wrapping_ns;
	
template <class T> T Max(const NumVec<T> &v) { T mv = 0; for (int k=0; k<v.m(); k++) mv = (v(k) > mv ? v(k) : mv); return mv; }

// > function passed to nth_element
bool compareFcn(double a, double b) { return (a>b); }

// find the size of the Ncoeffs:th coefficient
double FindNthCoeffSize(const CurveletData &c, int Ncoeffs)
{
	vector<double> vals;
	int nel;
	
	if (Ncoeffs <= 0)
		return HUGE_VAL;  // return a very large threshold
	
	// count nr of elements
	nel = countnel(c);
	if (Ncoeffs >= nel)
		return 0.0;
	
	// copy absolute values to a vector
	vals.resize(nel, 0.0);
	vector<double>::iterator it = vals.begin();
	for (int j=0; j<c.size(); j++) {
		for (int k=0; k<c[j].size(); k++) {
			cpx *data = c[j][k].data();
			for (int idx=0; idx < c[j][k].m() * c[j][k].n(); idx++)
				*it++ = abs(data[idx]);
		}
	}
	
	// find nth largest element using STL call
	nth_element(vals.begin(), vals.begin()+Ncoeffs-1, vals.end(), compareFcn);
	
	return vals[Ncoeffs-1];
}


/* threshold coefficients */
void ThresholdCoeffs(CurveletData &c, double thrsh)
{
	for (int j=0; j<c.size(); j++) {
		for (int k=0; k<c[j].size(); k++) {
			cpx *data = c[j][k].data();
			for (int id=0; id < c[j][k].m() * c[j][k].n(); id++) {
				if (abs(data[id]) < thrsh)
					data[id] = cpx(0.0, 0.0);
			}
		}
	}
}

// Threshold coefficients on selected levels
CurvUtilError ThresholdCoeffs(CurveletData &c, double thrsh, const IntNumVec &levels)
{
	for (int levi=0; levi<levels.m(); levi++) {
		int j = levels(levi);
		if (j<0 || j >= c.size())
			return CUError_BadParams;
		
		for (int k=0; k<c[j].size(); k++) {
			cpx *data = c[j][k].data();
			for (int id=0; id < c[j][k].m() * c[j][k].n(); id++) {
				if (abs(data[id]) < thrsh)
					data[id] = cpx(0.0, 0.0);
			}
		}
	}
	
	return CUError_NoError;
}




} // end namespace
