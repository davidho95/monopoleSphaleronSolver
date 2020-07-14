#include "LATfield2.hpp"
#include <complex>
#include "../../include/electroweakSphaleron/ElectroweakTheory.hpp"
#include "../../include/monopoleSphaleron/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../../include/Matrix.hpp"
#include "../../include/GradDescentSolverBBStep.hpp"
#include "../../include/Su2Tools.hpp"
#include "../../include/MonopoleFileTools.hpp"
#include "../../include/electroweakSphaleron/ElectroweakFieldTools.hpp"
#include <iostream>
#include <fstream>
#include <bits/stdc++.h> 
#include <iostream> 
#include <sys/stat.h> 
#include <sys/types.h> 

int main(int argc, char **argv)
{
  srand(time(NULL) + parallel.rank());
  std::string outputPath;
  std::string inputPath;
  int sz = 16;
  int n = 2;
  int m = 2;
  double vev = 1;
  double gaugeCoupling = 1;
  double selfCoupling = 1;
  double tanSqMixingAngle = 0.286;
  double fluxQuanta = 1;
  double zAspect = 1;

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
        break;
      case 'i':
        inputPath = argv[++i];
        break;
      case 'z':
        zAspect = atof(argv[++i]);
        break;
      case 'q':
        tanSqMixingAngle = atof(argv[++i]);
        break;

    }
  }

  parallel.initialize(n,m);

  int dim = 3;
  int xSz = sz;
  int ySz = sz;
  int zSz = 2;
  int latSize[dim] = {xSz, ySz, zSz};
  int haloSize = 1;
  int numMatrices = 4;
  int numRows = 2;
  int numCols = 2;

  // cout << zSz << endl;

  LATfield2::Lattice lattice(dim, latSize, haloSize);
  LATfield2::Field<complex<double> > field(lattice, numMatrices, numRows, numCols, 0);
  monsta::ElectroweakTheory theory(gaugeCoupling, tanSqMixingAngle, vev, selfCoupling);

  LATfield2::Site site(lattice);

  double initialStep = 0.01;
  double maxStepSize = 0.01;
  double tol = 1e-6;
  int minNumSteps = 10000;
  int maxNumSteps = 200000;

  monsta::setVacuumField(field, theory);
  theory.applyBoundaryConditions(field);

  // theory.applyBoundaryConditions(field);
  // double E = theory.computeEnergy(field);
  // cout << E << endl;

  if (inputPath != "")
  {
    minNumSteps = 1000;
    monsta::readRawField(field, inputPath + "/rawData");
    theory.applyBoundaryConditions(field);
  }
  else
  {
    for (site.first(); site.test(); site.next())
    {
      int yCoord = site.coord(1) - ySz / 2;
      int zCoord = site.coord(2) - zSz / 2;

      double r = sqrt(pow(yCoord,2) + pow(zCoord,2));

      for (int ii = 0; ii < 3; ii++)
      {
        // srand(site.coord(0)*xSz+site.coord(1));
        std::vector<double> su2Vec = {1e-6*(rand() % 2000 - 1000), 1e-6*(rand() % 2000 - 1000), 1e-6*(rand() % 2000 - 1000)};
        monsta::Matrix su2Mat = monsta::vecToSu2(su2Vec);
        field(site, ii, 0, 0) = su2Mat(0, 0);
        field(site, ii, 0, 1) = su2Mat(0, 1);
        field(site, ii, 1, 0) = su2Mat(1, 0);
        field(site, ii, 1, 1) = su2Mat(1, 1);
        theory.postProcess(field, site, ii);
      }
    }
  }

  monsta::addConstantMagneticField(field, theory, fluxQuanta, 2);
  theory.applyBoundaryConditions(field);

  monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize, minNumSteps, true);
  solver.solve(theory, field);
  solver = monsta::GradDescentSolver(tol, maxNumSteps, initialStep, maxStepSize, minNumSteps);
  solver.solve(theory, field);
  monsta::writeRawField(field, outputPath + "/rawData");
  monsta::writeCoords(field, outputPath + "/coords");
  monsta::writeHiggsMagnitude(field, outputPath + "/higgsData");
  monsta::writeMagneticField(field, outputPath + "/magneticFieldData", theory);
  monsta::writeEnergyDensity(field, outputPath + "/energyData", theory);
}
