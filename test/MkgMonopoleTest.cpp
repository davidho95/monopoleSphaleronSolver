#include "LATfield2.hpp"
#include "GradientChecker.hpp"
#include "../src/MkgMonopoleSolver.hpp"

int main() {
  parallel.initialize(1,1);

  srand(time(NULL) + parallel.rank());
  // srand(1);

  int dim = 3;
  int sz = 5;
  int latSize[dim] = {sz, sz, sz};
  int haloSize = 1;
  int numCpts = 5;

  double tol = 1e-5;

  LATfield2::Lattice lat(dim, latSize, haloSize);
  LATfield2::Site site(lat);

  double solverTol = 0;
  int maxIterations = 0;

  double electricCharge = std::rand() % 10 / 10.;
  double vev = std::rand() % 10 / 10.;
  double selfCoupling = std::rand() % 10 / 10.;
  double g = std::rand() % 10 / 10.;
  int southPos = std::rand() % (sz - 1);
  int northPos = std::rand() % sz;
  while (northPos <= southPos)
  {
    northPos = std::rand() % sz;
  }

  monsta::MkgMonopoleSolver solver(lat, numCpts, solverTol, maxIterations,
    electricCharge, vev, selfCoupling, g, southPos, northPos);
  monsta::GradientChecker<double> gradChecker(tol, true);
  gradChecker.checkGradients(solver);
  
  return 0;
}