#include "LATfield2.hpp"
#include "GradientChecker.hpp"
#include "../src/Phi4Solver.hpp"

int main() {
  parallel.initialize(1,1);

  srand(time(NULL) + parallel.rank());

  int dim = 3;
  int sz = 3;
  int latSize[dim] = {sz, sz, sz};
  int haloSize = 1;
  int numCpts = 1;

  double tol = 1e-6;

  LATfield2::Lattice lat(dim, latSize, haloSize);

  monsta::Phi4Solver solver(lat,1,1,1,1,1);
  monsta::GradientChecker<double> gradChecker(tol, false);

  gradChecker.checkGradients(solver);

  return 0;
}