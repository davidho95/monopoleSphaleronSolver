#include "LATfield2.hpp"
#include <complex>
#include "../src/ElectroweakTheory.hpp"
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
  int sz = 16;
  int n = 2;
  int m = 2;
  double xAspect = 1;
  double vev = 1;
  double gaugeCoupling = 1;
  double mixingAngle = 2*atan(1) / 3;
  double selfCoupling = 1;
  std::string inputPath;
  std::string outputPath;

  for (int i=1 ; i < argc ; i++ ){
    if ( argv[i][0] != '-' )
      continue;
    switch(argv[i][1]) {
        case 'p':
        outputPath = argv[++i];
        break;
        case 'i':
        inputPath = argv[++i];
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
      case 'v':
        vev = atof(argv[++i]);
        break;
      case 'g':
        gaugeCoupling = atof(argv[++i]);
        break;
      case 'l':
        selfCoupling = atof(argv[++i]);
        break;
      case 'q':
        mixingAngle = atof(argv[++i]);
        break;
      case 'x':
        xAspect = atof(argv[++i]);
        break;
    }
  }
  parallel.initialize(n,m);

  srand(1);

  int dim = 3;
  int xSz = round(xAspect*sz);
  int ySz = sz;
  int zSz = sz;
  int latSize[dim] = {xSz, ySz, zSz};
  int haloSize = 2;
  int numMatrices = 4;
  int numRows = 2;
  int numCols = 2;


  LATfield2::Lattice lattice(dim, latSize, haloSize);
  LATfield2::Field<complex<double> > field(lattice, numMatrices, numRows, numCols, 0);
  monsta::ElectroweakTheory theory(gaugeCoupling, mixingAngle, vev, selfCoupling);
  monsta::GeorgiGlashowSu2Theory ggTheory(gaugeCoupling, vev, selfCoupling);

  monsta::setVacuumField(field, theory);
  // monsta::readRawField(field, inputPath + "/rawData");
  monsta::addConstantMagneticField(field, theory, 1);
  theory.applyBoundaryConditions(field);

  double tol = 1e-8;
  // monsta::TheoryChecker gradChecker(tol);
  // gradChecker.checkGradients(theory, field);

  double initialStep = 0.01;
  double maxStepSize = 0.05;
  int maxNumSteps = 10000;

  LATfield2::Site site(lattice);

  for(site.first(); site.test(); site.next())
  {
    for (int ii = 0; ii < 3; ii++)
    {
      // std::vector<double> su2Vec = monsta::su2ToVec(monsta::Matrix(field, site, ii));
      // field(site, 3, (ii + 1) % 2, (ii + 1) / 2) = cos(su2Vec[2]) - 1i*sin(su2Vec[2]);
      // field(site, 3, (ii + 1) % 2, (ii + 1) / 2) = -su2Vec[2];
    }
  }
  theory.applyBoundaryConditions(field);

  double E = theory.computeEnergy(field);
  COUT << E << endl;

  monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize);
  // solver.solve(theory, field);

  monsta::writeCoords(field, outputPath + "/coords");
  monsta::writeEnergyDensity(field, outputPath + "/energyData", theory);
  monsta::writeHiggsField(field, outputPath + "/higgsData", theory);
  monsta::writeMagneticField(field, outputPath + "/magneticFieldData", theory);

}