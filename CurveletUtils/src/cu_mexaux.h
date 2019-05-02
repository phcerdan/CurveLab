/*
   Copyright (C) 2004 Caltech
   Written by Lexing Ying
	 Modified 2007 by Tobias Gebaeck
*/

#ifndef _CU_MEXAUX_H_
#define _CU_MEXAUX_H_

#include "mex.h"
#include "matrix.h"

using namespace std;

#include "nummat.hpp"
#include "offmat.hpp"
#include "numvec.hpp"

#ifndef NO3D_CURVELETS
namespace fdct3d_ns {
#include "numtns.hpp"
}
#endif

FDCT_WRAPPING_NS_BEGIN_NAMESPACE

#ifndef NO3D_CURVELETS
using namespace fdct3d_ns;
#endif

inline void mex2cpp(const mxArray*& md, int& cd);
inline void cpp2mex(const int& cd, mxArray*& md);

inline void mex2cpp(const mxArray*& md, double& cd);
inline void cpp2mex(const double& cd, mxArray*& md);

inline void mex2cpp(const mxArray*& md, CpxOffMat& cd);
inline void cpp2mex(const CpxOffMat& cd, mxArray*& md);

inline void mex2cpp(const mxArray*& md, CpxNumMat& cd);
inline void cpp2mex(const CpxNumMat& cd, mxArray*& md);

#ifndef NO3D_CURVELETS
inline void mex2cpp(const mxArray*& md, DblNumTns& cd);
inline void cpp2mex(const DblNumTns& cd, mxArray*& md);

inline void mex2cpp(const mxArray*& md, CpxNumTns& cd);
inline void cpp2mex(const CpxNumTns& cd, mxArray*& md);
#endif

template <class T> inline void mex2cpp(const mxArray*& md, std::vector<T>& cd);
template <class T> inline void cpp2mex(const std::vector<T>& cd, mxArray*& md);

//----------------------int
inline void mex2cpp(const mxArray*& md, int& cd)
{
  cd = int(mxGetScalar(md));
  return;
}
inline void cpp2mex(const int& cd, mxArray*& md)
{
  md = mxCreateDoubleScalar(cd);
  return;
}

//----------------------double
inline void mex2cpp(const mxArray*& md, double& cd)
{
  cd = mxGetScalar(md);
  return;
}
inline void cpp2mex(const double& cd, mxArray*& md)
{
  md = mxCreateDoubleScalar(cd);
  return;
}

//----------------------cpxoffmat
inline void mex2cpp(const mxArray*& md, CpxOffMat& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md);
  double* xr = mxGetPr(md);
  double* xi = mxGetPi(md);
  int s = -m/2;
  int t = -n/2;
  cd.resize(m,n);
  if(xr!=NULL && xi!=NULL) {
	 int cnt = 0;
	 for(int j=t; j<t+n; j++)
		for(int i=s; i<s+m; i++) {
		  cd(i,j) = cpx(xr[cnt], xi[cnt]);
		  cnt++;
		}
  } else if(xr!=NULL && xi==NULL) {
	 int cnt = 0;
	 for(int j=t; j<t+n; j++)
		for(int i=s; i<s+m; i++) {
		  cd(i,j) = cpx(xr[cnt], 0);
		  cnt++;
		}
  } else if(xr==NULL && xi!=NULL) {
	 int cnt = 0;
	 for(int j=t; j<t+n; j++)
		for(int i=s; i<s+m; i++) {
		  cd(i,j) = cpx(0, xi[cnt]);
		  cnt++;
		}
  }
  return;
}
inline void cpp2mex(const CpxOffMat& cd, mxArray*& md)
{
  int m = cd.m();
  int n = cd.n();
  int s = -m/2;
  int t = -n/2;
  md = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
  double* xr = mxGetPr(md);
  double* xi = mxGetPi(md);
  int cnt = 0;
  for(int j=t; j<t+n; j++)
	 for(int i=s; i<s+m; i++) {
		xr[cnt] = real(cd(i,j));
		xi[cnt] = imag(cd(i,j));
		cnt++;
	 }
  return;
}


//----------------------cpxnummat
inline void mex2cpp(const mxArray*& md, CpxNumMat& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md);
  double* xr = mxGetPr(md);
  double* xi = mxGetPi(md);
  cd.resize(m,n);
  if(xr!=NULL && xi!=NULL) {
	 int cnt = 0;
	 for(int j=0; j<n; j++)
		for(int i=0; i<m; i++) {
		  cd(i,j) = cpx(xr[cnt], xi[cnt]);
		  cnt++;
		}
  } else if(xr!=NULL && xi==NULL) {
	 int cnt = 0;
	 for(int j=0; j<n; j++)
		for(int i=0; i<m; i++) {
		  cd(i,j) = cpx(xr[cnt], 0);
		  cnt++;
		}
  } else if(xr==NULL && xi!=NULL) {
	 int cnt = 0;
	 for(int j=0; j<n; j++)
		for(int i=0; i<m; i++) {
		  cd(i,j) = cpx(0, xi[cnt]);
		  cnt++;
		}
  }
  return;
}
inline void cpp2mex(const CpxNumMat& cd, mxArray*& md)
{
  int m = cd.m();
  int n = cd.n();
  md = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
  double* xr = mxGetPr(md);
  double* xi = mxGetPi(md);
  int cnt = 0;
  for(int j=0; j<n; j++)
	 for(int i=0; i<m; i++) {
		xr[cnt] = real(cd(i,j));
		xi[cnt] = imag(cd(i,j));
		cnt++;
	 }
  return;
}


//----------------------dblnummat
inline void mex2cpp(const mxArray*& md, DblNumMat& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md);
  double* xr = mxGetPr(md);
  cd.resize(m,n);
  
  if (xr != NULL) {
    int cnt = 0;
    for(int j=0; j<n; j++)
      for(int i=0; i<m; i++) {
        cd(i,j) = xr[cnt];
        cnt++;
      }
  }
  return;
}

inline void cpp2mex(const DblNumMat& cd, mxArray*& md)
{
  int m = cd.m();
  int n = cd.n();
  md = mxCreateDoubleMatrix(m, n, mxREAL);
  double* xr = mxGetPr(md);
  int cnt = 0;
  for(int j=0; j<n; j++)
    for(int i=0; i<m; i++) {
      xr[cnt] = cd(i,j);
      cnt++;
    }
  return;
}

#ifndef NO3D_CURVELETS

//----------------------cpxnumtns
inline void mex2cpp(const mxArray*& md, CpxNumTns& cd)
{
  const int* dims = mxGetDimensions(md);
  int m = dims[0];
  int n = dims[1];
  int p = dims[2];
  double* xr = mxGetPr(md);
  double* xi = mxGetPi(md);
  cd.resize(m,n,p);
  if(xr!=NULL && xi!=NULL) {
		int cnt = 0;
		for(int k=0; k<p; k++)
			for(int j=0; j<n; j++)
				for(int i=0; i<m; i++) {
					cd(i,j,k) = cpx(xr[cnt], xi[cnt]);
					cnt++;
				}
  } else if(xr!=NULL && xi==NULL) {
		int cnt = 0;
		for(int k=0; k<p; k++)
			for(int j=0; j<n; j++)
				for(int i=0; i<m; i++) {
					cd(i,j,k) = cpx(xr[cnt], 0);
					cnt++;
				}
  } else if(xr==NULL && xi!=NULL) {
		int cnt = 0;
		for(int k=0; k<p; k++)
			for(int j=0; j<n; j++)
				for(int i=0; i<m; i++) {
					cd(i,j,k) = cpx(0, xi[cnt]);
					cnt++;
				}
  }
  return;
}
inline void cpp2mex(const CpxNumTns& cd, mxArray*& md)
{
  int m = cd.m();
  int n = cd.n();
  int p = cd.p();
  int ndim = 3;
  int dims[3];  dims[0] = m;  dims[1] = n;  dims[2] = p;
  md = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxCOMPLEX);
  double* xr = mxGetPr(md);
  double* xi = mxGetPi(md);
  int cnt = 0;
  for(int k=0; k<p; k++)
		for(int j=0; j<n; j++)
			for(int i=0; i<m; i++) {
				xr[cnt] = real(cd(i,j,k));
				xi[cnt] = imag(cd(i,j,k));
				cnt++;
			}
  //cd.resize(0,0,0);
  return;
}

//----------------------dblnumtns
inline void mex2cpp(const mxArray*& md, DblNumTns& cd)
{
  const int* dims = mxGetDimensions(md);
  int m = dims[0];
  int n = dims[1];
  int p = dims[2];
  double* xr = mxGetPr(md);
  cd.resize(m,n,p);
  if(xr!=NULL) {
		int cnt = 0;
		for(int k=0; k<p; k++)
			for(int j=0; j<n; j++)
				for(int i=0; i<m; i++) {
					cd(i,j,k) = xr[cnt];
					cnt++;
				}
  }
  return;
}

inline void cpp2mex(const DblNumTns& cd, mxArray*& md)
{
  int m = cd.m();
  int n = cd.n();
  int p = cd.p();
  int ndim = 3;
  int dims[3];  dims[0] = m;  dims[1] = n;  dims[2] = p;
  md = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
  double* xr = mxGetPr(md);
  int cnt = 0;
  for(int k=0; k<p; k++)
		for(int j=0; j<n; j++)
			for(int i=0; i<m; i++) {
				xr[cnt] = cd(i,j,k);
				cnt++;
			}
  //cd.resize(0,0,0);
  return;
}

#endif

//----------------------IntNumMat
inline void mex2cpp(const mxArray*& md, IntNumMat& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md);
  double* xr = mxGetPr(md);
  cd.resize(m,n);
  
  if (xr != NULL) {
    int cnt = 0;
    for(int j=0; j<n; j++) {
      for(int i=0; i<m; i++) {
        cd(i,j) = (int) xr[cnt++];
      }
		}
  }
	return;
}

inline void cpp2mex(const IntNumMat& cd, mxArray*& md)
{
  int m = cd.m();
  int n = cd.n();
  md = mxCreateDoubleMatrix(m, n, mxREAL);
  double* xr = mxGetPr(md);
  int cnt = 0;
  for(int j=0; j<n; j++) {
    for(int i=0; i<m; i++) {
      xr[cnt++] = (double) cd(i,j);
    }
	}
	return;
}



//----------------------IntNumVec
inline void mex2cpp(const mxArray*& md, IntNumVec& cd)
{
  int m = mxGetM(md);
  int n = mxGetN(md);
  assert(m==1 || n==1);
  double* xr = mxGetPr(md);
  cd.resize(m*n);
  
  if (xr != NULL) {
    for(int j=0, cnt=0; j<m*n; j++)
      cd(j) = (int) xr[cnt++];
  }
  return;
}

inline void cpp2mex(const IntNumVec& cd, mxArray*& md)
{
  int m = cd.m();
  md = mxCreateDoubleMatrix(m, 1, mxREAL);
  double* xr = mxGetPr(md);
  for(int j=0, cnt=0; j<m; j++)
    xr[cnt++] = (double) cd(j);

  return;
}


//----------------------vector<...>
template <class T> inline void mex2cpp(const mxArray*& md, std::vector<T>& cd)
{
  int m = mxGetM(md); 
  int n = mxGetN(md); 
	assert(min(m,n)==1);
	
	int len = max(m,n);
  cd.resize(len);
  for(int ci=0; ci<len; ci++) {
	 const mxArray*tt = mxGetCell(md, ci);
	 mex2cpp(tt, cd[ci]);
  }
  return;
}
template <class T> inline void cpp2mex(const std::vector<T>& cd, mxArray*& md)
{
  int n = cd.size();
  md = mxCreateCellMatrix(1, n);
  for(int ci=0; ci<n; ci++) {
	 mxArray* ss;	 cpp2mex(cd[ci], ss);
	 mxSetCell(md, ci, ss);
  }
  return;
}

FDCT_WRAPPING_NS_END_NAMESPACE

#endif
