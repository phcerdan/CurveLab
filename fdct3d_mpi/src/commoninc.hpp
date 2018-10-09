/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
        Written by Lexing Ying
*/

#ifndef _COMMONINC_HPP_
#define _COMMONINC_HPP_

// STL stuff
#include <fstream>
#include <iostream>
#include <sstream>

#include <cassert>
#include <cfloat>
#include <cmath>
#include <complex>
#include <string>

#include <algorithm>
#include <deque>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>
#include <cstring>

// FFT stuff
#include "fftw.h"

// typedef double double;
typedef std::complex<double> cpx;

// AUX functions
inline int pow2(int l) {
  assert(l >= 0);
  return (1 << l);
}

// petsc stuff
#include "petscsnes.h"
#define iC(fun)     \
  {                 \
    int ierr = fun; \
    CHKERRQ(ierr);  \
  }
#define iA(expr)                                             \
  {                                                          \
    if (!(expr)) SETERRQ(1, "Assertion: " #expr " failed!"); \
  }

#endif
