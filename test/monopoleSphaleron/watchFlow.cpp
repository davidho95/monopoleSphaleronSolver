#include "LATfield2.hpp"
#include <complex>
#include "../src/GeorgiGlashowSu2TheoryNonUnitary.hpp"
#include "../src/GeorgiGlashowSu2EomTheory.hpp"
#include "../src/Matrix.hpp"
#include "../src/GradDescentSolverBBStep.hpp"
#include "../src/Su2Tools.hpp"
#include "../src/MonopoleFileTools.hpp"
#include "../src/MonopoleFieldTools.hpp"
#include "../src/TheoryChecker.hpp"
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
  double fluxQuanta = 1;
  int maxNumSteps = 200000;
  int outputIncrement = 1000;
  double xAspect = 1;

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
        break;
      case 'N':
        maxNumSteps = atoi(argv[++i]);
        break;
      case 'D':
        outputIncrement = atoi(argv[++i]);
        break;
      case 'x':
        xAspect = atof(argv[++i]);
    }
  }

  parallel.initialize(n,m);

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
  monsta::GeorgiGlashowSu2Theory theory(gaugeCoupling, vev, selfCoupling);

  LATfield2::Site site(lattice);

  double initialStep = 0.01;
  double maxStepSize = 0.025;
  double tol = 1e-10;
  int minNumSteps = 0;

  if (inputPath == "")
  {
    monsta::setVacuumField(field, theory);
    for (site.first(); site.test(); site.next())
    {
      monsta::Matrix higgsMat = monsta::vecToSu2LieAlg({0, 0, vev/sqrt(2)});//{ 0.01*(rand() % 100), 0.01*(rand() % 100), 0.01*(rand() % 100)});
      field(site, 3, 0, 0) = higgsMat(0, 0);
      field(site, 3, 0, 1) = higgsMat(0, 1);
      field(site, 3, 1, 0) = higgsMat(1, 0);
      field(site, 3, 1, 1) = higgsMat(1, 1);
    }
    theory.applyBoundaryConditions(field);

    for (site.first(); site.test(); site.next())
    {
      int yCoord = site.coord(1) - ySz / 2;
      int zCoord = site.coord(2) - zSz / 2;

      double r = sqrt(pow(yCoord,2) + pow(zCoord,2));

      // std::vector<double> su2Vec = {1e-3*exp(-pow(r,2) / 10), -1e-3*exp(-pow(r,2) / 10), 0};
      // std::vector<double> su2Vec = {0.001*(1e-5*double(rand() % 628318 - 314159)), 0.001*(1e-5*double(rand() % 628318 - 314159)), 0.001*(1e-5*double(rand() % 628318 - 314159))};
      for (int ii = 0; ii < 3; ii++)
      {
        // double randNum = 1e-6*(rand() % 2000 - 1000);
        std::vector<double> su2Vec = {1e-6*(rand() % 2000 - 1000), 1e-6*(rand() % 2000 - 1000), 1e-6*(rand() % 2000 - 1000)};
        monsta::Matrix su2Mat = monsta::vecToSu2(su2Vec);
        // su2Mat.print();
        field(site, ii, 0, 0) = su2Mat(0, 0);
        field(site, ii, 0, 1) = su2Mat(0, 1);
        field(site, ii, 1, 0) = su2Mat(1, 0);
        field(site, ii, 1, 1) = su2Mat(1, 1);

        theory.postProcess(field, site, ii);
      }

    }

    monsta::addConstantMagneticField(field, theory, -fluxQuanta);
  }
  else
  {
    monsta::readRawField(field, inputPath + "/rawData");
    theory.applyBoundaryConditions(field);
  }

  double E = theory.computeEnergy(field);
  COUT << E << endl;

  int numIters = 0;


  std::string currentOutputPath = outputPath + "/evolution0";

  mkdir(currentOutputPath.c_str(), 0777);

  monsta::writeRawField(field, currentOutputPath + "/rawData");
  monsta::writeCoords(field, currentOutputPath + "/coords");
  monsta::writeHiggsFieldMagnitude(field, currentOutputPath + "/higgsData");
  monsta::writeMagneticField(field, currentOutputPath + "/magneticFieldData", theory);
  monsta::writeEnergyDensity(field, currentOutputPath + "/energyData", theory);
  // monsta::writeUnitaryGaugeField(field, currentOutputPath + "/gaugeData");
  // monsta::writeUnitaryWField(field, currentOutputPath + "/wData");
  monsta::GradDescentSolver solver(tol, outputIncrement, initialStep, maxStepSize, outputIncrement);
  while (numIters < maxNumSteps)
  {
    numIters += outputIncrement;
    solver.solve(theory, field);

    currentOutputPath = outputPath + "/evolution" + std::to_string(numIters / outputIncrement);

    mkdir(currentOutputPath.c_str(), 0777);

    monsta::writeRawField(field, currentOutputPath + "/rawData");
    monsta::writeCoords(field, currentOutputPath + "/coords");
    monsta::writeHiggsFieldMagnitude(field, currentOutputPath + "/higgsData", theory);
    monsta::writeMagneticField(field, currentOutputPath + "/magneticFieldData", theory);
    monsta::writeEnergyDensity(field, currentOutputPath + "/energyData", theory);
    // monsta::writeUnitaryGaugeField(field, currentOutputPath + "/gaugeData");
    // monsta::writeUnitaryWField(field, currentOutputPath + "/wData");
  }

  cout << srand << endl;


}
