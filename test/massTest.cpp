#include "LATfield2.hpp"
#include <complex>
#include "../src/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../src/Matrix.hpp"
#include "../src/TheoryChecker.hpp"
#include "../src/GradDescentSolverBBStep.hpp"
#include "../src/Su2Tools.hpp"
#include "../src/MonopoleFileTools.hpp"
#include "../src/MonopoleFieldTools.hpp"
#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{

  std::string outputPath;
  int sz = 16;

  for (int i=1 ; i < argc ; i++ ){
    if ( argv[i][0] != '-' )
      continue;
    switch(argv[i][1]) {
      case 'p':
        outputPath = argv[++i];
        break;
      case 's':
        sz =  atoi(argv[++i]);
        break;
    }
  }

  parallel.initialize(2,2);

  srand(time(NULL) + parallel.rank());

  int dim = 3;
  int latSize[dim] = {sz, sz, sz};
  int haloSize = 1;
  int numMatrices = 4;
  int numRows = 2;
  int numCols = 2;

  double gaugeCoupling = 1;
  double vev = 1.0/sqrt(2);
  double selfCoupling = 1;

  LATfield2::Lattice lattice(dim, latSize, haloSize);
  LATfield2::Field<complex<double> > field(lattice, numMatrices, numRows, numCols, 0);
  field.updateHalo();

  LATfield2::Site site(lattice);

  for (int jj = 0; jj < 100; jj++)
  {
    selfCoupling = 0.01*(100 - jj);
    monsta::GeorgiGlashowSu2Theory theory(gaugeCoupling, vev, selfCoupling, {-1, 1, 1}, true);

    if (jj == 0)
    {
      for (site.first(); site.test(); site.next())
      {
        int xCoord = site.coord(0);
        int yCoord = site.coord(1);
        int zCoord = site.coord(2);
        double r = sqrt(pow(xCoord - sz/2, 2) + pow(yCoord - sz/2, 2) + pow(zCoord - sz/2, 2));
        for (int ii = 0; ii < numMatrices - 1; ii++)
        {
          monsta::Matrix su2Mat = monsta::vecToSu2({0, 0.1, 0.2});
          field(site, ii, 0, 0) = su2Mat(0, 0);
          field(site, ii, 0, 1) = su2Mat(0, 1);
          field(site, ii, 1, 0) = su2Mat(1, 0);
          field(site, ii, 1, 1) = su2Mat(1, 1);
        }
        field(site, 3, 0, 0) = vev;
        field(site, 3, 0, 1) = 0;
        field(site, 3, 1, 0) = 0;
        field(site, 3, 1, 1) = 0;
      }
      field.updateHalo();
    }

    double initialStep = 0.01;
    double maxStepSize = 0.1 / vev;
    double tol = 1e-4;
    int maxNumSteps = 10000;

    monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize);
    // solver.setVerbosity(false);
    solver.solve(theory, field);

    double E = theory.computeEnergy(field);

    if (parallel.rank() == 0)
    {
      ofstream fileStream;
      fileStream.open(outputPath + "/energies.txt", std::ios_base::app);
      fileStream << E << endl;
      fileStream.close();
    }
  }
}