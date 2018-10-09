/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
        Written by Lexing Ying
*/

#ifndef _CPXCRVLETOCR_HPP_
#define _CPXCRVLETOCR_HPP_

#include "numtns.hpp"

class CpxCrvletOcr {
 protected:
  std::vector<std::vector<int> > _nxs, _nys, _nzs;
  char _name[100];
  int _maxnb;
  int _count;

  int _clock;
  std::vector<std::vector<CpxNumTns> > _blocks;
  std::vector<std::vector<int> > _szvec;
  std::vector<std::vector<int> > _tmvec;

 public:
  CpxCrvletOcr(const char* name);  //  CpxCrvletOcr(const CpxCrvletOcr& D);
  ~CpxCrvletOcr();  //  CpxCrvletOcr& operator=(const CpxCrvletOcr& D);
  int setup(std::vector<std::vector<int> > nxs, std::vector<std::vector<int> > nys,
            std::vector<std::vector<int> > nzs, int ma);
  CpxNumTns& block(int s, int w);
  // access
  std::vector<std::vector<int> >& nxs() { return _nxs; }
  std::vector<std::vector<int> >& nys() { return _nys; }
  std::vector<std::vector<int> >& nzs() { return _nzs; }
};

#endif
