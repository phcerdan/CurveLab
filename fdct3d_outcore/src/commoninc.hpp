/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
        Written by Lexing Ying
*/

#ifndef _FDCT3DINC_HPP_
#define _FDCT3DINC_HPP_

// STL stuff
#include <string.h>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <algorithm>
#include <deque>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>
using namespace std;

// FFT stuff
#include "fftw.h"

// typedef double double;
typedef complex<double> cpx;

// AUX functions
inline int pow2(int l) {
  assert(l >= 0);
  return (1 << l);
}

#endif
