/*
 *  CurveletUtils.h
 *
 *  Created by Tobias Gebäck on 12.9.2007.
 *  Copyright 2007 ETH Zurich. All rights reserved.
 *
 */

using namespace std;

#include "nummat.hpp"
#include "numvec.hpp"

#ifndef NO3D_CURVELETS
namespace fdct3d_ns {
	#include "numtns.hpp"
}
#endif

#include "fdct_wrapping_inc.hpp"

namespace CurveletUtils {
	
	typedef std::vector< std::vector<fdct_wrapping_ns::CpxNumMat> > CurveletData;
#ifndef NO3D_CURVELETS
	typedef std::vector< std::vector<fdct3d_ns::CpxNumTns> > CurveletData3D;
#endif
	
	enum CurvUtilError { CUError_NoError = 0, CUError_BadParams = 1 };
	
	void VarAndMean(const fdct_wrapping_ns::DblNumVec &vec, double &mean, double &var);
	
	/* Mapping of indices and coordinates to closest index on other level */
	void MapXYToLevel(const CurveletData &cdata, double x, double y, int tolev, int &ti, int &tj, int dir);
	void MapToLevel(const CurveletData &cdata, int fromlev, int fromdir, int fi, int fj, int tolev, int &ti, int &tj, int dir);
	void MapXYToGrid(int m, int n, double x, double y, int &ti, int &tj);
	void MapToGrid(int m, int n, const CurveletData &cdata, int fromlev, int fromdir, int fi, int fj, int &ti, int &tj);
	int GetQuadrant(const CurveletData &cdata, int level, int diridx);
	
	/* Compute an all zero fdct, i.e. create the correct data structure */
	CurvUtilError zero_fdct(int m, int n, CurveletData &cdata, int allcurvelets=0, int nscales=-1, int nangles_coarse=-1);
	
	/* count nr of nonzero elements in CurveletData */
	void countnnz(const CurveletData &cdata, int &nnz, int &nel);
	CurvUtilError countnnz(const CurveletData &cdata, const fdct_wrapping_ns::IntNumVec &levels, int &nnz, int &nel);
	int countnnz(const CurveletData &cdata);
	CurvUtilError countnnz(const CurveletData &cdata, const fdct_wrapping_ns::IntNumVec &levels, int &nnz);
	
	/* count nr of elements only */
	int countnel(const CurveletData &cdata);
	CurvUtilError countnel(const CurveletData &cdata, const fdct_wrapping_ns::IntNumVec &levels, int &nel);
	
	/* combinations of real/imaginary parts of curvelets */
	enum CurvFilterType { CF_Abs = 0, CF_Real, CF_Imag, CF_Sum }; 
	 
	/* Options passed to ExtractMultipleDirFields */
	struct DirFieldOptions {
		int CurveletSize;  // The number of positions away that a curvelet coefficient influences
		int CircSumOffset; // running circular sum over -CSO..CSO
		int LocMaxGap;  // the required gap between maximums (for multiple field extraction)
		CurvFilterType FilterType;
	};
	
	/* Extract a directional field describing preferred direction at each grid point */
	CurvUtilError ExtractMultipleDirFields(int m, int n, const CurveletData &cdata, const fdct_wrapping_ns::IntNumVec &levels, int nrfields, std::vector< std::vector<fdct_wrapping_ns::DblNumMat> > &dirfields, DirFieldOptions *opts = NULL);
	
	/* Extract magnitude and direction (as index 0..7) at each grid point */
	CurvUtilError ExtractDirAndMag(int m, int n, const CurveletData &cdata, const fdct_wrapping_ns::IntNumVec &levels, fdct_wrapping_ns::DblNumMat &mag, fdct_wrapping_ns::IntNumMat &dir);
	
	/* Extract directional field, using only projections onto x- and y-axes */
	CurvUtilError ExtractDirections(int m, int n, const CurveletData &cdata, const fdct_wrapping_ns::IntNumVec &levels, fdct_wrapping_ns::DblNumMat &xcmp, fdct_wrapping_ns::DblNumMat &ycmp);
	
	/* Get magnitudes on specified levels, summed over all directions */
	CurvUtilError GetMagnitude(int m, int n, const CurveletData &cdata, const fdct_wrapping_ns::IntNumVec &levels, fdct_wrapping_ns::DblNumMat &mag, int crvlt_size);
	
	/* For mirror-extended curvelets: Get magnitudes on specified levels, summed over all directions */
	CurvUtilError GetMagnitudeME(int m, int n, const CurveletData &cdata, const fdct_wrapping_ns::IntNumVec &levels, fdct_wrapping_ns::DblNumMat &mag, int crvlt_size);
	
#ifndef NO3D_CURVELETS
	/* Extract directional field (normals to surfaces) at each grid point on (m x n x p) grid */
	CurvUtilError ExtractDirFields3D(int m, int n, int p, const CurveletData3D &cdata3, const fdct_wrapping_ns::IntNumVec &levels, int nrfields,  std::vector< std::vector<fdct3d_ns::DblNumTns> > &dirfields, DirFieldOptions *opts = NULL);
	
	/* Get magnitudes of 3D curvelets, summed over all directions, on (m x n x p) grid */
	CurvUtilError GetMagnitude3D(int m, int n, int p, const CurveletData3D &cdata3, const fdct_wrapping_ns::IntNumVec &levels, fdct3d_ns::DblNumTns &mag, int crvlt_size);
#endif
	
	/* find Nth largest coefficient */
	double FindNthCoeffSize(const CurveletData &c, int Ncoeffs);
	
	/* threshold coefficients */
	void ThresholdCoeffs(CurveletData &c, double thrsh);
	CurvUtilError ThresholdCoeffs(CurveletData &c, double thrsh, const fdct_wrapping_ns::IntNumVec &levels);
	
	void ChooseByCurvature(CurveletData &c, double T1, double T2, int lastklev);
	
}
