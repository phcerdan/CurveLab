/*
   Copyright (C) 2004 Caltech
   Written by Lexing Ying
*/

#include "fdct_wrapping.hpp"
#include "fdct_wrapping_inline.hpp"

using namespace std;
using namespace fdct_wrapping_ns;

int optionsCreate(const char* optfile, std::map<std::string, std::string>& options) {
  options.clear();
  std::ifstream fin(optfile);
  assert(fin.good());
  std::string name;
  fin >> name;
  while (fin.good()) {
    char cont[100];
    fin.getline(cont, 99);
    options[name] = std::string(cont);
    fin >> name;
  }
  fin.close();
  return 0;
}

int main(int argc, char** argv) {
  clock_t ck0, ck1;

  assert(argc == 2);
  // get options
  std::map<std::string, std::string> opts;
  optionsCreate(argv[1], opts);

  // get input data
  std::map<std::string, std::string>::iterator mi;

  int m;
  mi = opts.find("-m");
  assert(mi != opts.end());
  {
    std::istringstream ss((*mi).second);
    ss >> m;
  }
  int n;
  mi = opts.find("-n");
  assert(mi != opts.end());
  {
    std::istringstream ss((*mi).second);
    ss >> n;
  }

  int nbscales;
  mi = opts.find("-nbscales");
  assert(mi != opts.end());
  {
    std::istringstream ss((*mi).second);
    ss >> nbscales;
  }

  int nbangles_coarse;
  mi = opts.find("-nbangles_coarse");
  assert(mi != opts.end());
  {
    std::istringstream ss((*mi).second);
    ss >> nbangles_coarse;
  }

  int ac;
  mi = opts.find("-ac");
  assert(mi != opts.end());
  {
    std::istringstream ss((*mi).second);
    ss >> ac;
  }

  srand48((long)time(NULL));
  CpxNumMat x(m, n);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++) x(i, j) = cpx(drand48(), drand48());

  ck0 = clock();

  // fdct_wrapping_
  std::vector<std::vector<CpxNumMat> > c;  // std::vector<int> extra;
  fdct_wrapping(m, n, nbscales, nbangles_coarse, ac, x, c);
  ck1 = clock();
  std::cout << "FDCT_WRAPPING_  takes " << double(ck1 - ck0) / CLOCKS_PER_SEC
       << " seconds" << std::endl;
  ck0 = ck1;

  // ifdct_wrapping_
  CpxNumMat y(x);
  clear(y);
  ifdct_wrapping(m, n, nbscales, nbangles_coarse, ac, c, y);
  ck1 = clock();
  std::cout << "IFDCT_WRAPPING_ takes " << double(ck1 - ck0) / CLOCKS_PER_SEC
       << " seconds" << std::endl;
  ck0 = ck1;

  CpxNumMat e(m, n);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++) e(i, j) = x(i, j) - y(i, j);
  std::cerr << "accuracy of inversion " << sqrt(energy(e) / (m * n)) << std::endl;

  std::vector<std::vector<double> > sx, sy;
  std::vector<std::vector<double> > fx, fy;
  std::vector<std::vector<int> > nx, ny;
  fdct_wrapping_param(m, n, nbscales, nbangles_coarse, ac, sx, sy, fx, fy, nx,
                      ny);

  return 0;
}
