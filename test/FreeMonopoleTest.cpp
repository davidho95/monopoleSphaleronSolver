#include "LATfield2.hpp"
#include "GradientChecker.hpp"
#include "../src/FreeMonopoleSolver.hpp"

int main() {
  parallel.initialize(1,1);

  srand(time(NULL) + parallel.rank());

  int dim = 3;
  int sz = 5;
  int latSize[dim] = {sz, sz, sz};
  int haloSize = 1;
  int numCpts = 3;

  double tol = 1e-6;

  LATfield2::Lattice lat(dim, latSize, haloSize);
  LATfield2::Site site(lat);


  double solverTol = 0;
  int maxIterations = 0;
  double g = 1;
  int southPos = std::rand() % sz;
  int northPos = std::rand() % sz;
  while (northPos < southPos) {
    northPos = std::rand() % sz;
  }

  monsta::FreeMonopoleSolver solver(lat, numCpts, solverTol, maxIterations, g, southPos, northPos);
  monsta::GradientChecker<double> gradChecker(tol, true);

  gradChecker.checkGradients(solver);

  return 0;
}