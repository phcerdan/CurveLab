/*
 *  basic_utils.cpp
 *
 *  Created by Tobias Gebaeck on 12.9.2007.
 *  Copyright 2007 ETH Zurich. All rights reserved.
 *
 */

#include "CurveletUtils.h"
#include "fdct_wrapping.hpp"

using namespace fdct_wrapping_ns;

namespace CurveletUtils {

// Create an all-zero curvelet data structure, from image with dimensions m*n
CurvUtilError zero_fdct(int m, int n, CurveletData &cdata, int allcurvelets/* =0 */, int nscales/* = -1 */, int nangles_coarse/* =-1*/)
{
  if (nscales < 0) nscales = (int) ceil(log2(min(m,n)) - 3);
  if (nangles_coarse < 0) nangles_coarse = 16;

  if (nscales == 0 || nangles_coarse == 0)
    return CUError_BadParams;
  
  CpxNumMat x(m,n);  // all zero matrix
  fdct_wrapping(m, n, nscales, nangles_coarse, allcurvelets, x, cdata);
  
  return CUError_NoError;
}


// Count the number of non-zeros, and total number of elements in curvelet data structure
void countnnz(const CurveletData &cdata, int &nnz, int &nel)
{
	nnz = 0;
	nel = 0;
	for (int j=0; j<cdata.size(); j++) {
		for (int l=0; l<cdata[j].size(); l++) {
			nel += cdata[j][l].m() * cdata[j][l].n();
			for (int idx=0; idx<cdata[j][l].m() * cdata[j][l].n(); idx++) {
				if (cdata[j][l].data()[idx] != cpx(0.0,0.0))
					nnz++;
			}
		}
	}
}

// Count the number of non-zero coefficients, and total number, on selected (zero-based) levels
CurvUtilError countnnz(const CurveletData &cdata, const IntNumVec &levels, int &nnz, int &nel)
{
	nnz = 0;
	nel = 0;
	for (int j=0; j<levels.m(); j++) {
		int lev = levels(j);
		if (lev < 0 && lev >= cdata.size())
			return CUError_BadParams;
		
		for (int l=0; l<cdata[lev].size(); l++) {
			nel += cdata[lev][l].m() * cdata[lev][l].n();
			for (int idx=0; idx<cdata[lev][l].m() * cdata[lev][l].n(); idx++) {
				if (cdata[lev][l].data()[idx] != cpx(0.0,0.0))
					nnz++;
			}
		}
	}	
	
	return CUError_NoError;
}

// Count number of nonzero curvelet coefficients
int countnnz(const CurveletData &cdata)
{
	int nel, nnz;
	countnnz(cdata, nnz, nel);
	return nnz;
}

// Count number of nonzero curvelet coefficients on selected (zero-based) levels
CurvUtilError countnnz(const CurveletData &cdata, const IntNumVec &levels, int &nnz)
{
	int nel;
	return countnnz(cdata, levels, nnz, nel);	
}

// Count number of curvelet coefficients
int countnel(const CurveletData &cdata)
{
	int nel = 0;
	for (int j=0; j<cdata.size(); j++) {
		for (int l=0; l<cdata[j].size(); l++) {
			nel += cdata[j][l].m() * cdata[j][l].n();
		}
	}
	
	return nel;	
}

// Count number of curvelet coefficients on selected (zero-based) levels
CurvUtilError countnel(const CurveletData &cdata, const IntNumVec &levels, int &nel)
{
	nel = 0;
	for (int j=0; j<levels.m(); j++) {
		int lev = levels(j);
		if (lev < 0 && lev >= cdata.size())
			return CUError_BadParams;
		
		for (int l=0; l<cdata[lev].size(); l++) {
			nel += cdata[lev][l].m() * cdata[lev][l].n();
		}
	}	
	
	return CUError_NoError;	
}


// Maps coordinates x,y to indices on level tolev in CurveletData structure cdata
// x and y in [0 1)
// dir is 0 or 1 for North/South resp. East/West
void MapXYToLevel(const CurveletData &cdata, double x, double y, int tolev, int &ti, int &tj, int dir)
{
	int L1,L2;
	switch (dir) {
	case 0:
		L1 = cdata[tolev][0].m();
		L2 = cdata[tolev][0].n();
		ti = int(x * L1 + 0.5) % L1;
		tj = int(y * L2 + 0.5) % L2;
		break;
	case 1:
		L1 = cdata[tolev][cdata[tolev].size()/4].m();
		L2 = cdata[tolev][cdata[tolev].size()/4].n();
		ti = int(x * L1 + 0.5) % L1;
		tj = int(y * L2 + 0.5) % L2;
		break;
	}
}

	
// Maps indices (fi,fj) on level fromlev and dir fromdir (actual direction index), to level tolev and 
// either N/S or E/W grid, (dir=0 for N/S and dir=1 for E/W)
void MapToLevel(const CurveletData &cdata, int fromlev, int fromdir, int fi, int fj, int tolev, int &ti, int &tj, int dir)
{
	MapXYToLevel(cdata, fi / (double) cdata[fromlev][fromdir].m(), fj / (double) cdata[fromlev][fromdir].n(), tolev, ti, tj, dir);
}

	
// Map coordinates (x,y) to grid with size m*n, assuming grid is on [0,1) x [0,1)
void MapXYToGrid(int m, int n, double x, double y, int &ti, int &tj)
{
	ti = (int(x * m + 0.5)) % m;
	tj = (int(y * n + 0.5)) % n;
}

// Map indices (fi,fj) on level fromlev and dir fromdir to a grid with size m*n, assumed to be on [0,1) x [0,1)
void MapToGrid(int m, int n, const CurveletData &cdata, int fromlev, int fromdir, int fi, int fj, int &ti, int &tj)
{
	MapXYToGrid(m, n, fi / (double) cdata[fromlev][fromdir].m(), fj / (double) cdata[fromlev][fromdir].n(), ti, tj);
}

// Get quadrant (0..3 for N,E,S,W) for direction diridx on level level
int GetQuadrant(const CurveletData &cdata, int level, int diridx)
{
	return diridx / (cdata[level].size() / 4);
}


// Compute variance and mean of a vector vec
void VarAndMean(const DblNumVec &vec, double &mean, double &var)
{
	double sum=0.0, sqsum=0.0;
	for (int k=0; k<vec.m(); k++) {
		sum += vec(k);
		sqsum += vec(k)*vec(k);
	}
	mean = sum/vec.m();
	var = sqsum/vec.m() - mean*mean;
}

}
