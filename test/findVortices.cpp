#include "LATfield2.hpp"
#include <complex>
#include "../src/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../src/GeorgiGlashowSu2EomTheory.hpp"
#include "../src/Matrix.hpp"
#include "../src/GradDescentSolverBBStep.hpp"
#include "../src/GradDescentSolverFluxPreserving.hpp"
#include "../src/GradDescentSolverCorePreserving.hpp"
#include "../src/Su2Tools.hpp"
#include "../src/MonopoleFileTools.hpp"
#include "../src/MonopoleFieldTools.hpp"
#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{
  std::string outputPath;
  std::string inputPath;
  int sz = 16;
  int n = 2;
  int m = 2;
  double vev = 1;
  double gaugeCoupling = 1;
  double selfCoupling = 1;
  double fluxQuanta = 1;

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
      case 'i':
        inputPath = argv[++i];
        break;
      case 'v':
        vev = atof(argv[++i]);
        break;
      case 'g':
        gaugeCoupling = atof(argv[++i]);
        break;
      case 'l':
        selfCoupling = atof(argv[++i]);
        break;
      case 'b':
        fluxQuanta = atoi(argv[++i]);
    }
  }

  parallel.initialize(n,m);

  int dim = 3;
  int xSz = 1;
  int ySz = sz;
  int zSz = sz;
  int latSize[dim] = {xSz, ySz, zSz};
  int haloSize = 2;
  int numMatrices = 4;
  int numRows = 2;
  int numCols = 2;

  LATfield2::Lattice lattice(dim, latSize, haloSize);
  LATfield2::Field<complex<double> > field(lattice, numMatrices, numRows, numCols, 0);
  monsta::GeorgiGlashowSu2Theory theory(gaugeCoupling, vev, selfCoupling, {0, 0, 0}, false);

  LATfield2::Site site(lattice);

  double initialStep = 0.01;
  double maxStepSize = 0.1;
  double tol = 1e-6;
  int maxNumSteps = 200000;

  monsta::setVacuumField(field, theory);

  for (site.first(); site.test(); site.next())
  {
    if (site.coord(1) >= sz/2 - 2 && site.coord(1) <= sz/2 + 2 && site.coord(2) >= sz/2 - 2 && site.coord(2) <= sz/2 + 2)
    {
      monsta::Matrix su2Mat = monsta::vecToSu2({0.1,0.1,0});
      for (int ii = 0; ii < 3; ii++)
      {
        field(site, ii, 0, 0) = su2Mat(0, 0);
        field(site, ii, 0, 1) = su2Mat(0, 1);
        field(site, ii, 1, 0) = su2Mat(1, 0);
        field(site, ii, 1, 1) = su2Mat(1, 1);
      }
    }
  }
  // COUT << lattice.size(0) << endl;
  monsta::addConstantMagneticField(field, theory, -fluxQuanta);

  double flux = 0;
  for (site.first(); site.test(); site.next())
  {
    flux += theory.getMagneticField(field, site, 0);
  }
  parallel.sum(flux);

  double pi = 4*atan(1);
  flux = flux * theory.getGaugeCoupling() / (4*pi);
  COUT << flux << endl;

  monsta::GradDescentSolverCorePreserving corePreservingSolver(tol, maxNumSteps, initialStep, maxStepSize);
  monsta::GradDescentSolver solver(tol, 200, initialStep, maxStepSize);
  solver.solve(theory, field);

  flux = 0;
  for (site.first(); site.test(); site.next())
  {
    flux += theory.getMagneticField(field, site, 0);
  }
  parallel.sum(flux);

  flux = flux * theory.getGaugeCoupling() / (4*pi);
  COUT << flux << endl;

  corePreservingSolver.solve(theory, field);

  monsta::writeRawField(field, outputPath + "/rawData");
  monsta::writeCoords(field, outputPath + "/coords");
  monsta::writeHiggsFieldUnitary(field, outputPath + "/higgsData");
  monsta::writeMagneticField(field, outputPath + "/magneticFieldData", theory);
  monsta::writeEnergyDensity(field, outputPath + "/energyData", theory);
  monsta::writeUnitaryGaugeField(field, outputPath + "/gaugeData");


}
