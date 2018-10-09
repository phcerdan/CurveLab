/*
  Copyright (C) 2004 Caltech
  Written by Lexing Ying
*/

#ifndef _FDCT_USFFT_INC_HPP_
#define _FDCT_USFFT_INC_HPP_

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
// FFT stuff
#include "fftw.h"

#define FDCT_USFFT_NS_BEGIN_NAMESPACE namespace fdct_usfft_ns {
#define FDCT_USFFT_NS_END_NAMESPACE }

FDCT_USFFT_NS_BEGIN_NAMESPACE

typedef std::complex<double> cpx;

// AUX functions
inline int pow2(int l) {
  assert(l >= 0);
  return (1 << l);
}

FDCT_USFFT_NS_END_NAMESPACE

#endif
