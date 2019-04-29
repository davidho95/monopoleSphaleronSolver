#include "LATfield2.hpp"
#include <complex>
#include "../src/GeorgiGlashowSu2TheoryUnitaryOptimised.hpp"
#include "../src/GeorgiGlashowSu2EomTheory.hpp"
#include "../src/Matrix.hpp"
#include "../src/GradDescentSolverBBStep.hpp"
#include "../src/Su2Tools.hpp"
#include "../src/MonopoleFileTools.hpp"
#include "../src/MonopoleFieldTools.hpp"
#include "../src/TheoryChecker.hpp"
#include <iostream>
#include <fstream>
#include <ctime>

int main(int argc, char **argv)
{
  clock_t begin = clock();

  std::string outputPath;
  int sz = 16;
  int n = 2;
  int m = 2;

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
      case 'n':
        n = atoi(argv[++i]);
        break;
      case 'm':
        m = atoi(argv[++i]);
        break;
    }
  }

  cout << m << " " << n << endl;

  parallel.initialize(n,m);

  srand(1);//time(NULL) + parallel.rank());

  int dim = 3;
  int latSize[dim] = {sz, sz, sz};
  int haloSize = 1;
  int numMatrices = 4;
  int numRows = 2;
  int numCols = 2;

  double gaugeCoupling = 1;
  double vev = 1;//sqrt(10);
  double selfCoupling = 1;

  LATfield2::Lattice lattice(dim, latSize, haloSize);
  LATfield2::Field<complex<double> > field(lattice, numMatrices, numRows, numCols, 0);
  field.updateHalo();

  LATfield2::Site site(lattice);

  for (site.first(); site.test(); site.next())
  {
    int xCoord = site.coord(0);
    int yCoord = site.coord(1);
    int zCoord = site.coord(2);
    double r = sqrt(pow(xCoord - sz/2, 2) + pow(yCoord - sz/2, 2) + pow(zCoord - sz/2, 2));
    for (int ii = 0; ii < numMatrices - 1; ii++)
    {
      // monsta::Matrix su2Mat = monsta::vecToSu2({0, 0.001, 0.01});
      monsta::Matrix su2Mat = monsta::vecToSu2({rand() % 10, rand() % 10, rand() % 10});
      field(site, ii, 0, 0) = su2Mat(0, 0);
      field(site, ii, 0, 1) = su2Mat(0, 1);
      field(site, ii, 1, 0) = su2Mat(1, 0);
      field(site, ii, 1, 1) = su2Mat(1, 1);
    }
    field(site, 3, 0, 0) = vev / sqrt(2);
    field(site, 3, 0, 1) = 0;
    field(site, 3, 1, 0) = 0;
    field(site, 3, 1, 1) = 0;
  }
  field.updateHalo();

  double initialStep = 0.001;
  double maxStepSize = 0.05;
  double tol = 1e-6;
  int maxNumSteps = 10000;

  monsta::GeorgiGlashowSu2Theory theory(gaugeCoupling, vev, selfCoupling, {-1, -1, -1}, true);
  monsta::GeorgiGlashowSu2EomTheory eomTheory(gaugeCoupling, vev, selfCoupling, {-1, -1, -1}, true);

  monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize);
  // solver.setVerbosity(false);
  // solver.solve(theory, field);

  // site.setCoord(rand() % sz, rand() % sz, rand() % sz);
  // (theory.getLocalGradient(field, site, 1)).print();
  // monsta::Matrix(field, site, 1).print();

  monsta::TheoryChecker checker(tol);
  checker.checkGradients(eomTheory, field);

  // double gradSq = 0;

  // for (site.first(); site.test(); site.next())
  // {
  //   for (int matIdx = 0; matIdx < 4; matIdx++)
  //   {
  //     monsta::Matrix gradMat = theory.getLocalGradient(field, site, matIdx);
  //     monsta::Matrix fieldMat(field, site, matIdx);
  //     if (matIdx < 3)
  //     {
  //       // cout << 2.*sqrt(0.5*monsta::trace(gradMat*conjugateTranspose(gradMat))) + monsta::trace(gradMat*conjugateTranspose(fieldMat)) << endl;
  //       gradSq += real(2.*sqrt(0.5*monsta::trace(gradMat*conjugateTranspose(gradMat))) + monsta::trace(gradMat*conjugateTranspose(fieldMat)));
  //     }
  //     else
  //     {
  //       gradSq += pow(real(gradMat(0, 0)), 2);
  //     }
  //   }
  // }

  // parallel.sum(gradSq);
  // COUT << gradSq << endl;
  double E = eomTheory.computeEnergy(field);
  COUT << E << endl;

  clock_t end = clock();
  COUT << double(end - begin) / CLOCKS_PER_SEC << endl;
}