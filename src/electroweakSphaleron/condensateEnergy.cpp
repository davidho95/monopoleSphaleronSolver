#include "LATfield2.hpp"
#include <complex>
#include "../../include/monopoleSphaleron/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../../include/monopoleSphaleron/GeorgiGlashowSu2EomTheory.hpp"
#include "../../include/Matrix.hpp"
#include "../../include/GradDescentSolverBBStep.hpp"
#include "../../include/electroweakSphaleron/GradDescentSolverChigusaElectroweak.hpp"
#include "../../include/Su2Tools.hpp"
#include "../../include/MonopoleFileTools.hpp"
#include "../../include/monopoleSphaleron/MonopoleFieldTools.hpp"
#include <iostream>
#include <fstream>
#include <ctime>
#include <bits/stdc++.h> 
#include <iostream> 
#include <sys/stat.h> 
#include <sys/types.h> 


int main(int argc, char **argv)
{
  clock_t begin = clock();

  std::string outputPath;
  std::string inputPath;
  int sz = 16;
  int n = 2;
  int m = 2;
  double vev = 1;
  double gaugeCoupling = 1;
  double selfCoupling = 1;
  int sep = sz/2;
  double startVev;
  double endVev;
  double numIncrements;
  double xAspect = 1;
  double tanSqMixingAngle = 0.268;
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
      case 'q':
        tanSqMixingAngle = atof(argv[++i]);
        break;
      case 'b':
        fluxQuanta = atoi(argv[++i]);
        break;
      case 'S':
        startVev = atof(argv[++i]);
        break;
      case 'E':
        endVev = atof(argv[++i]);
        break;
      case 'I':
        numIncrements = atof(argv[++i]);
        break;
      case 'x':
        xAspect = atof(argv[++i]);
        break;
    }
  }

  parallel.initialize(n,m);

  int dim = 3;
  int xSz = 1;
  int ySz = sz;
  int zSz = sz;
  int latSize[dim] = {xSz, ySz, zSz};
  int haloSize = 1;
  int numMatrices = 4;
  int numRows = 2;
  int numCols = 2;


  LATfield2::Lattice lattice(dim, latSize, haloSize);
  LATfield2::Site site(lattice);

  double initialStep = 0.01;
  double maxStepSize = 0.05;
  double tol = 1e-3;
  int maxNumSteps = 200000;
  double abortGrad = 0.5;

  monsta::ElectroweakTheory theory(gaugeCoupling, tanSqMixingAngle, vev, selfCoupling, {0, 0, 0});
  LATfield2::Field<complex<double> > field(lattice, numMatrices, numRows, numCols, 0);
  LATfield2::Field<complex<double> > referenceField(lattice, numMatrices, numRows, numCols, 0);

  monsta::readRawField(field, inputPath + "/rawData");
  theory.applyBoundaryConditions(field);

  std::vector<double> vevs;
  for (int ii = 0; ii < numIncrements; ii++)
  {
    vevs.push_back(endVev - (endVev - startVev)*double(ii)/(numIncrements - 1));
  }

  for (int ii = 0; ii < numIncrements; ii++)
  {
    vev = vevs[ii];
    monsta::GradDescentSolver solver(1e-10, 10000, initialStep, maxStepSize*vev, 3000);

    monsta::ElectroweakTheory theory(gaugeCoupling, tanSqMixingAngle, vev, selfCoupling, {0, 0, 0}, false);

    monsta::setVacuumField(field, theory);
    theory.applyBoundaryConditions(field);

    for (site.first(); site.test(); site.next())
    {
      int yCoord = site.coord(1) - ySz / 2;
      int zCoord = site.coord(2) - zSz / 2;

      double r = sqrt(pow(yCoord,2) + pow(zCoord,2));

      for (int ii = 0; ii < 3; ii++)
      {
        std::vector<double> su2Vec = {1e-6*(rand() % 2000 - 1000), 1e-6*(rand() % 2000 - 1000), 1e-6*(rand() % 2000 - 1000)};
        monsta::Matrix su2Mat = monsta::vecToSu2(su2Vec);
        field(site, ii, 0, 0) = su2Mat(0, 0);
        field(site, ii, 0, 1) = su2Mat(0, 1);
        field(site, ii, 1, 0) = su2Mat(1, 0);
        field(site, ii, 1, 1) = su2Mat(1, 1);

        theory.postProcess(field, site, ii);
      }

    }

    monsta::addConstantMagneticField(field, theory, -fluxQuanta);
    theory.applyBoundaryConditions(field);

    solver.solve(theory, field);

    double E = theory.computeEnergy(field);
    if (parallel.rank() == 1)
    {
      ofstream fileStream;
      fileStream.open(outputPath + "/energies.txt", std::ios_base::app);
      fileStream << E << endl;
      fileStream.close();
    }

    monsta::writeCoords(field, outputPath + "/coords");
    monsta::writeHiggsFieldMagnitude(field, outputPath + "/higgsData");
    monsta::writeMagneticField(field, outputPath + "/magneticFieldData", theory);
    monsta::writeEnergyDensity(field, outputPath + "/energyData", theory);
  }
}
