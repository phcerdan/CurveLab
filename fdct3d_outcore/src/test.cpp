/* FDCT3D (Fast 3d Curvelet Transform)
   Copyright (C) 2004 Caltech
        Written by Lexing Ying
*/

#include "fdct3d.hpp"
#include "fdct3dinline.hpp"

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
  time_t tm0, tm1;

  assert(argc == 2);
  std::map<std::string, std::string> opts;
  optionsCreate(argv[1], opts);
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
  int p;
  mi = opts.find("-p");
  assert(mi != opts.end());
  {
    std::istringstream ss((*mi).second);
    ss >> p;
  }

  int nbscales;
  mi = opts.find("-nbscales");
  assert(mi != opts.end());
  {
    std::istringstream ss((*mi).second);
    ss >> nbscales;
  }

  int nbdstz_coarse;
  mi = opts.find("-nbdstz_coarse");
  assert(mi != opts.end());
  {
    std::istringstream ss((*mi).second);
    ss >> nbdstz_coarse;
  }

  srand48((long)time(NULL));
  CpxNumTns x(m, n, p);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++) x(i, j, k) = cpx(drand48(), 0);

  tm0 = time(NULL);

  std::vector<std::vector<double> > fxs, fys, fzs;
  std::vector<std::vector<int> > nxs, nys, nzs;
  fdct3d_param(m, n, p, nbscales, nbdstz_coarse, fxs, fys, fzs, nxs, nys, nzs);
  tm1 = time(NULL);
  std::cout << "fdct3d_param " << difftime(tm1, tm0) << " seconds" << std::endl;
  tm0 = tm1;

  // 1. fdct3d
  CpxCrvletOcr c("tmpc");
  CpxNumTns w;
  fdct3d_forward(m, n, p, nbscales, nbdstz_coarse, x, c, w);
  tm1 = time(NULL);
  std::cout << "fdct3d_forward " << difftime(tm1, tm0) << " seconds" << std::endl;
  tm0 = tm1;

  // 2. ifdct3d
  // fdct3d_inverse(m, n, p, nbscales, nbdstz_coarse, c,w, x);
  // tm1 = time(NULL);  std::cout<<"fdct3d_inverse "<<difftime(tm1,tm0)<<"
  // seconds"<<std::endl;  tm0 = tm1;

  CpxNumTns newx(x);
  clear(newx);
  fdct3d_inverse(m, n, p, nbscales, nbdstz_coarse, c, w, newx);
  tm1 = time(NULL);
  std::cout << "fdct3d_inverse " << difftime(tm1, tm0) << " seconds" << std::endl;
  tm0 = tm1;  // std::cerr<<energy(newx)<<std::endl;

  double mv = 0.0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      for (int k = 0; k < p; k++) mv = std::max(mv, abs(newx(i, j, k) - x(i, j, k)));
  std::cerr << "max error " << mv << std::endl;

  return 0;
}
