/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
        Written by Lexing Ying
*/

#ifndef _CPXCRVLETPRTD_HPP_
#define _CPXCRVLETPRTD_HPP_

#include "numtns.hpp"

//-----------------------------------------
// Complex Curvelet Partitioned
class CpxCrvletPrtd {
 protected:
  std::vector<std::vector<int> > _nx, _ny, _nz;  // size of
  std::vector<std::vector<int> > _owners;
  std::vector<std::vector<int> > _sizes;
  std::vector<std::vector<bool> > _exists;
  std::vector<std::vector<CpxNumTns> > _blocks;

 public:
  CpxCrvletPrtd() { ; }
  CpxCrvletPrtd(const CpxCrvletPrtd& D);
  ~CpxCrvletPrtd() { ; }
  CpxCrvletPrtd& operator=(const CpxCrvletPrtd& D);
  int setup(std::vector<std::vector<int> > nx, std::vector<std::vector<int> > ny,
            std::vector<std::vector<int> > nz, std::vector<std::vector<int> >& owners);
  int expand(std::vector<std::vector<bool> >& newexists);
  int scatter(std::vector<std::vector<bool> >& newexists);
  int shift(std::vector<std::vector<int> >& newowners);
  int discard();
  int combine();
  // access
  std::vector<std::vector<int> >& nx() { return _nx; }
  std::vector<std::vector<int> >& ny() { return _ny; }
  std::vector<std::vector<int> >& nz() { return _nz; }

  std::vector<std::vector<int> >& owners() { return _owners; }
  std::vector<std::vector<int> >& sizes() { return _sizes; }
  std::vector<std::vector<bool> >& exists() { return _exists; }

  CpxNumTns& block(int s, int w) {
    assert(_exists[s][w] == true);
    return _blocks[s][w];
  }
  double globalenergy();
  int check();

  // extra
  int mpirank() const {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
  }
  int mpisize() const {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
  }
};

#endif
