/*
 *  extractdirs.cpp
 *
 *  Created by Tobias Gebaeck on 12.9.2007.
 *  Copyright 2007 ETH Zurich. All rights reserved.
 *
 */

#include "CurveletUtils.h"
#include "numvec.hpp"
#include <algorithm>

using namespace std;
using namespace fdct_wrapping_ns;
using namespace fdct3d_ns;

namespace CurveletUtils {

const double PI = 3.141592653589793;

// compute the direction of the curvelet (the direction along the valleys/ridges)
// in the form of cosine and sine of angle
void ComputeAngle(int idx, int nangles, double &cosa, double &sina)
{
  int nangpquad = nangles/4;
  double tana = 1.0 - (1.0 + 2.0*(idx % nangpquad)) / nangpquad;
  double scf = sqrt(1 + tana*tana);
	
  switch (idx / (nangpquad)) { 
		case 0:
		case 2:
			cosa = 1.0 / scf;
			sina = tana / scf;
			break;	
    case 1:
    case 3:
      cosa = tana / scf;
			sina = -1.0 / scf; 
      break;
  }
}

	
// Compute the direction of a 3D curvelet (normal to its plane), based on its index
// Angles are spherical coordinates (theta, phi)
	void ComputeAngle3D(int idx, int totangles, int angperdir, double &theta, double &phi)
	{
		int angperface = angperdir*angperdir;
		int face = idx / angperface;
		int k1 = idx % angperdir;
		int k2 = (idx % angperface) / angperdir;
		double alpha = (k1 - angperdir/2 + 0.5) / (angperdir/2.0);
		double beta = (k2 - angperdir/2 + 0.5) / (angperdir/2.0); 
		
		switch (face) {
			case 0:
				theta = PI/2.0 - atan(alpha); phi = PI/2.0 - atan(beta); break;
			case 1:
				theta = PI/2.0 - atan(beta); phi = atan(alpha); break;
			case 2:
				theta = atan(sqrt(alpha*alpha + beta*beta)); phi = atan2(beta, alpha); break;
			case 3:
				theta = PI/2.0 + atan(alpha); phi = -PI/2.0 - atan(beta); break;
			case 4:
				theta = PI/2.0 + atan(beta); phi = PI + atan(alpha); break;
			case 5:
				theta = PI - atan(sqrt(alpha*alpha + beta*beta)); phi = PI + atan2(beta,alpha); break;
		}
		
		if (phi > PI) phi -= 2.0*PI;
		if (phi <= -PI) phi += 2.0*PI;
	}

// Map coordinates (ix,iy) with max dimensions (m,n) to level given by levmat
// return indices on this level in (levx,levy)
void MapCoordToLevel(const CpxNumMat &levmat, int ix, int iy, int m, int n, int& levx, int& levy)
{
  levx = (int) (ix / double(m) * double(levmat.m()) + 0.5) % levmat.m();
  levy = (int) (iy / double(n) * double(levmat.n()) + 0.5) % levmat.n();
}

	void MapCoordToLevel3D(const CpxNumTns &levtns, int ix, int iy, int iz, int m, int n, int p, int &levx, int &levy, int &levz)
	{
		levx = (int) (ix / double(m) * double(levtns.m()) + 0.5) % levtns.m();
		levy = (int) (iy / double(n) * double(levtns.n()) + 0.5) % levtns.n();
		levz = (int) (iz / double(p) * double(levtns.p()) + 0.5) % levtns.p();
	}
	
// Map coordinates (ix,iy) with max dimensions (m,n) to level of actual size (mp,np)
// return indices on this level in (levx,levy)
void MapCoordToLevelME(int mp, int np, int ix, int iy, int m, int n, int& levx, int& levy)
{
	levx = (int) (ix / double(m) * double(mp) + 0.5);
	levy = (int) (iy / double(n) * double(np) + 0.5);
	if (levx >= mp) levx = mp-1;
	if (levy >= np) levy = np-1;
}

	
	
// Find maximum number of angles on selected levels, and the stepsize in angles for other levels
// i.e. each angle on level levels(lx) maps to angfact(lx) levels on finest selected level
void GetAngleFactors(const CurveletData &cdata, const IntNumVec &levels, IntNumVec &angfact, int &maxang)
{
	maxang = 0;
  angfact.resize(levels.m());
  for (int levx=0; levx<levels.m(); levx++) {
    int lev = levels(levx);
    if (cdata[lev].size() > maxang)
      maxang = cdata[lev].size();
  }
  for (int levx=0; levx<levels.m(); levx++)
    angfact(levx) = maxang / cdata[levels(levx)].size();
	
}


	
// datastructure for sorting, keeping original indices
struct CoSortEntry {
	int idx;
	double val;
};

// order function, passed to std::sort
bool cosortgt(const CoSortEntry &a, const CoSortEntry &b) {
	return (a.val>b.val);
}
	
// perform a running, periodic sum of first len entries
// in vec, going sumofs steps in each direction
void circularsum(DblNumVec &vec, int len, int sumofs) 
{
	DblNumVec tmpvec(vec);
	clear(vec);
	for (int k=0; k<len; k++) {
		for (int m=k-sumofs; m<=k+sumofs; m++) {
			vec(k) += tmpvec((m+len) % len);
		}
	}
}
	
	
// returns true if vec(idx) is a local maximum, regarding
// first len entries of vec as a periodic vector
bool islocalmax(DblNumVec &vec, int len, int idx)
{
	int lidx = (idx - 1 + len) % len;
	int uidx = (idx + 1) % len;
	return (vec(idx) >= vec(lidx) && vec(idx) >= vec(uidx));
}

// computes distance between indices idx1 and idx2,
// if they are indices into a periodic array of size period
int periodicdistance(int idx1, int idx2, int period)
{
	return min(abs(idx1-idx2), abs(abs(idx1 - idx2) - period));
}

	
/*
	ExtractMultipleDirFields()
 
 Extracts major direcional fields from the levels (0-based indices) of curvelet data cdata.
 Extracts nrfields fields, each of size (m,n), into dirfields[0..nrfields-1][0..1](0..m-1,0..n-1)
 (the second index is x- and y-direction)
 
 The function sums magnitudes of curvelet coefficients on selected levels, and at each point
 picks the direction with largest coefficient (and second, third,... largest if nrfields > 1)
 
 Options:
 DirFieldOptions.CurveletSize				determines the number of grid points (on its own grid) away
																		that a curvelet influences
 DirFieldOptions.CircSumOffset      determines number of directions away in both directions that 
																		are included in a running sum over the directions
 DirFieldOptions.LocMaxGap					the required minimal distance in direction numbers between
																		successive maxima saved in fields
 
*/
	
CurvUtilError ExtractMultipleDirFields(int m, int n, const CurveletData &cdata, const IntNumVec &levels, int nrfields, std::vector< std::vector<DblNumMat> > &dirfields, DirFieldOptions *a_opts)
{ 
	// Check inputs
	if (nrfields <= 0 || m<=0 || n<=0)
		return CUError_BadParams;
	for (int k=0; k<levels.m(); k++) {
		if (levels(k) >= cdata.size())
			return CUError_BadParams;
	}
		
	// Get parameters (or defaults)
	DirFieldOptions *opts;
	if (a_opts == NULL) {
		opts = new DirFieldOptions;
		opts->CircSumOffset = 1;
		opts->CurveletSize = 1;
		opts->LocMaxGap = 3;
		opts->FilterType = CF_Abs;
	}
	else {
		opts = a_opts;
		if (opts->CurveletSize < 0 || opts->CircSumOffset < 0 || opts->LocMaxGap < 0)
			return CUError_BadParams;
	}
	
	// initialize selected number of fields
	dirfields.resize(nrfields);   
	for (int k=0; k<nrfields; k++) {
		dirfields[k].resize(2); 
		dirfields[k][0].resize(m,n); 
		dirfields[k][1].resize(m,n);
	}
	
	// find maximum number of angles, and the number of angles for each fine angle at each level
  int maxang;
	IntNumVec angfact;
	GetAngleFactors(cdata, levels, angfact, maxang);
	int mangd2 = max(maxang/2, 1); 
	
	// At each pixel, add coefficient magnitudes and pick maximum
	DblNumVec csum(maxang);
  for (int xi=0; xi<m; xi++) {
    for (int yi=0; yi<n; yi++) {
      clear(csum);
      for (int levx=0; levx<levels.m(); levx++) {
        int lev = levels(levx);
        int nang = max((int) cdata[lev].size()/2, 1);  // loop over half the angles, since other half has same magnitudes (assuming real-valued image)
        for (int angx=0; angx<nang; angx++) {
          int lx, ly;
          MapCoordToLevel(cdata[lev][angx], xi, yi, m, n, lx, ly);
          // Add neighboring curvelets
          double cf = 0.0;
					const CpxNumMat *mat = &(cdata[lev][angx]);
					int lox = max(lx-opts->CurveletSize, 0);
					int hix = min(lx+opts->CurveletSize, cdata[lev][angx].m()-1);
					int loy = max(ly-opts->CurveletSize, 0);
					int hiy = min(ly+opts->CurveletSize, cdata[lev][angx].n()-1);
					
					switch (opts->FilterType) {
						case CF_Abs:
							for (int lxi=lox; lxi <= hix; lxi++)
								for (int lyi=loy; lyi <= hiy; lyi++)
									cf += abs((*mat)(lxi,lyi));
							break;
						case CF_Real:
							for (int lxi=lox; lxi <= hix; lxi++)
								for (int lyi=loy; lyi <= hiy; lyi++)
									cf += abs(real((*mat)(lxi,lyi)));							
							break;
						case CF_Imag:
							for (int lxi=lox; lxi <= hix; lxi++)
								for (int lyi=loy; lyi <= hiy; lyi++)
									cf += abs(imag((*mat)(lxi,lyi)));							
							break;
					  case CF_Sum:
							for (int lxi=lox; lxi <= hix; lxi++)
								for (int lyi=loy; lyi <= hiy; lyi++)
									cf += abs(real((*mat)(lxi,lyi))) + abs(imag((*mat)(lxi,lyi)));							
							break;
					}
					
					// Add to all relevant (sub-)angles
          for (int idx=angx*angfact(levx); idx<(angx+1)*angfact(levx); idx++)
            csum(idx) += cf;
        }
      }
			
      // running (circular) sum, and find maximal directions
			circularsum(csum, mangd2, min(opts->CircSumOffset, mangd2));
			std::vector<CoSortEntry> svec;
			svec.resize(mangd2);
      for (int idx=0; idx<mangd2; idx++) {
				svec[idx].val = csum(idx);
				svec[idx].idx = idx;
      }
			sort(svec.begin(), svec.end(), cosortgt);
      
			// compute angles and field for maximal directions
	    double cosa, sina;
			int diridx = 0;
			for (int idx=0; idx<nrfields; idx++) {
				while (diridx < svec.size() && !islocalmax(csum, mangd2, svec[diridx].idx) ) {
					// only accept second direction if far from previously selected directions and is local maximum
					bool ok = false;
					for (int fi=0; fi<idx && !ok; fi++) {
						if (periodicdistance(svec[diridx].idx, svec[fi].idx, svec.size()) > opts->LocMaxGap)
							ok = true;
					}
					if (ok) break;
					
					diridx++;
				}
				if (diridx == svec.size())
					break;
				
				//compute angles and field for next largest direction
				ComputeAngle(svec[diridx].idx, maxang, cosa, sina);
				dirfields[idx][0](xi,yi) = svec[diridx].val * cosa;
				dirfields[idx][1](xi,yi) = svec[diridx].val * sina;
				
				diridx++;
			}			
    }
  }  
	
	if (a_opts == NULL) 
		delete opts;
	
  return CUError_NoError;
}

	
	
/*
 Sums the curvelet coefficient magnitudes on specified levels, over all directions, onto a grid of size (m,n)
 */
CurvUtilError GetMagnitude(int m, int n, const CurveletData &cdata, const IntNumVec &levels, DblNumMat &mag, int crvlt_size)
{
	// Check inputs
	if (m<=0 || n<=0)
		return CUError_BadParams;
	for (int k=0; k<levels.m(); k++) {
		if (levels(k) >= cdata.size())
			return CUError_BadParams;
	}
	if (crvlt_size < 0)
		return CUError_BadParams;
	
	mag.resize(m,n);

	// At each pixel, add coefficient magnitudes and pick maximum
	double magsum;
  for (int xi=0; xi<m; xi++) {
    for (int yi=0; yi<n; yi++) {
      magsum = 0.0;
      for (int levx=0; levx<levels.m(); levx++) {
        int lev = levels(levx);
        int nang = max((int) cdata[lev].size()/2, 1);  // loop over half the angles, since other half has same magnitudes (assuming real-valued image)
        for (int angx=0; angx<nang; angx++) {
          int lx, ly;
          MapCoordToLevel(cdata[lev][angx], xi, yi, m, n, lx, ly);
          
					// Add neighboring curvelets, to compute average
          double avg = 0.0;
					const CpxNumMat *mat = &(cdata[lev][angx]);
					int nel = (1 + 2*crvlt_size) * (1 + 2*crvlt_size);
					int nely = 1 + 2*crvlt_size;
          for (int lxi=lx-crvlt_size; lxi <= lx+crvlt_size; lxi++) {
						if (lxi < 0)
							nel -= nely;
						else if	(lxi >= cdata[lev][angx].m()) 
							nel -=nely;
						else {
							for (int lyi=ly-crvlt_size; lyi <= ly+crvlt_size; lyi++) {
								if (lyi < 0)
									nel -= 1;
								else if (lyi >= cdata[lev][angx].n())
									nel -= 1;
								else
									avg += abs((*mat)(lxi,lyi));          
							}
						}
					}
					magsum += avg / nel;
        }
      }
			
			mag(xi,yi) = 2.0 * magsum;
		}
	}

  return CUError_NoError;
}
	

	
	/*
	 Sums the _mirror-extended_ curvelet coefficient magnitudes on specified levels, over all directions, onto a grid of size (m,n)
	 */
	CurvUtilError GetMagnitudeME(int m, int n, const CurveletData &cdata, const IntNumVec &levels, DblNumMat &mag, int crvlt_size)
	{
		// Check inputs
		if (m<=0 || n<=0)
			return CUError_BadParams;
		for (int k=0; k<levels.m(); k++) {
			if (levels(k) >= cdata.size())
				return CUError_BadParams;
		}
		if (crvlt_size < 0)
			return CUError_BadParams;
		
		mag.resize(m,n);
		
		// At each pixel, add coefficient magnitudes and pick maximum
		double magsum;
		for (int xi=0; xi<m; xi++) {
			for (int yi=0; yi<n; yi++) {
				magsum = 0.0;
				for (int levx=0; levx<levels.m(); levx++) {
					int lev = levels(levx);
					int nang = cdata[lev].size();
					for (int angx=0; angx<nang; angx++) {
						int lx, ly;
						int mp = (int)(cdata[lev][angx].m() / 2) + 1; // actual sizes for mirror-extended levels
						int np = (int)(cdata[lev][angx].n() / 2) + 1;
						MapCoordToLevelME(mp, np, xi, yi, m, n, lx, ly);
						
						// Add neighboring curvelets, to compute average
						double avg = 0.0;
						const CpxNumMat *mat = &(cdata[lev][angx]);
						int nel = (1 + 2*crvlt_size) * (1 + 2*crvlt_size);
						int nely = 1 + 2*crvlt_size;
						for (int lxi=lx-crvlt_size; lxi <= lx+crvlt_size; lxi++) {
							if (lxi < 0)
								nel -= nely;
							else if	(lxi >= mp) 
								nel -=nely;
							else {
								for (int lyi=ly-crvlt_size; lyi <= ly+crvlt_size; lyi++) {
									if (lyi < 0)
										nel -= 1;
									else if (lyi >= np)
										nel -= 1;
									else {
										avg += abs((*mat)(lxi,lyi));          
										avg += abs((*mat)(lxi, ((*mat).n() - lyi) % (*mat).n()));  // add direction pi - theta
									}
								}
							}
						}
						magsum += avg / nel;
					}
				}
				
				mag(xi,yi) = 2.0 * magsum;
			}
		}
		
		return CUError_NoError;
	}
	
													 
}

