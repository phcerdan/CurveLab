/*
   Copyright (C) 2004 Caltech
   Written by Lexing Ying
*/

#ifndef _FDCT_WRAPPING_INC_HPP_
#define _FDCT_WRAPPING_INC_HPP_

// STL stuff
#include <fstream>
#include <iostream>
#include <sstream>

#include <cassert>
#include <cfloat>
#include <cmath>
#include <complex>
#include <string>
#include <cstring>

#include <algorithm>
#include <deque>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>
// using namespace std;
// FFT stuff
#include "fftw.h"

#define FDCT_WRAPPING_NS_BEGIN_NAMESPACE namespace fdct_wrapping_ns {
#define FDCT_WRAPPING_NS_END_NAMESPACE }

FDCT_WRAPPING_NS_BEGIN_NAMESPACE

// Complex number
typedef std::complex<double> cpx;

// AUX functions
inline int pow2(int l) {
  assert(l >= 0);
  return (1 << l);
}

FDCT_WRAPPING_NS_END_NAMESPACE

#endif
