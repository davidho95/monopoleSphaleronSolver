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
  LATfield2::Field<double> testField(lat, numCpts);
  LATfield2::Site site(lat);

  // Initialise random test field
  for (site.first(); site.test(); site.next())
  {
    for (int iCpt = 0; iCpt < numCpts; iCpt++)
    {
      testField(site,iCpt) = (std::rand() % 100) / 100.0;
    }
  }
  testField.updateHalo();

  monsta::Phi4Solver solver(1,1,1,1);
  monsta::GradientChecker<double> gradChecker(tol);

  gradChecker.checkGradients(solver, testField);

  return 0;
}