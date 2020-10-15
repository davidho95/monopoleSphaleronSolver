#include "LATfield2.hpp"
#include <complex>
#include "../../include/monopoleSphaleron/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../../include/GradDescentSolverBBStepNoCheckerboard.hpp"
#include "../../include/Matrix.hpp"
#include "../../include/GradDescentSolverBBStep.hpp"
#include "../../include/Su2Tools.hpp"
#include "../../include/MonopoleFileTools.hpp"
#include "../../include/MonopoleFieldTools.hpp"
#include "../../include/monopoleInstanton/GeorgiGlashowSu2Theory4d.hpp"
#include "../../include/monopoleInstanton/InstantonFileTools.hpp"
#include <iostream>
#include <fstream>
#include <ctime>

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
      case 'v':
        vev = atof(argv[++i]);
        break;
      case 'g':
        gaugeCoupling = atof(argv[++i]);
        break;
      case 'l':
        selfCoupling = atof(argv[++i]);
        break;
      case 'x':
        xAspect = atof(argv[++i]);
        break;
      case 'i':
        inputPath = argv[++i];
        break;
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

  LATfield2::Site site(lattice);

  double initialStep = 0.01;
  double maxStepSize = 0.01;;
  double tol = 1e-3;
  int maxNumSteps = 10000;

  monsta::GeorgiGlashowSu2Theory theory(gaugeCoupling, vev, selfCoupling, {2, 0, 0}, true);
  monsta::GeorgiGlashowSu2Theory4d theory4d(gaugeCoupling, vev, selfCoupling);

  if (inputPath != "")
  {
    COUT << inputPath << endl;
    monsta::readRawField(field, inputPath + "/rawData");
    theory.applyBoundaryConditions(field);
  }
  else
  {
    monsta::setSingleMonopoleInitialConditions(field, theory);
    theory.applyBoundaryConditions(field);
  }


  LATfield2::Field<complex<double> > tempField(lattice, numMatrices, numRows, numCols, 0);

  for (site.first(); site.test(); site.next())
  {
    double pi = 4*atan(1);
    double xCoord = site.coord(0) - xSz/2 - 0.5;
    double yCoord = site.coord(1) - ySz/2 - 0.5;
    double zCoord = site.coord(2) - zSz/2 - 0.5;

    double r = sqrt(pow(xCoord, 2) + pow(yCoord, 2) + pow(zCoord, 2));
    for (int ii = 0; ii < 3; ii++)
    {

      std::vector<double> su2Vec = {0, 0, 0};

      if (ii == 1)
      {
        if (site.coord(0) == xSz/2 && site.coord(1) == ySz/2 && site.coord(2) == zSz/2)
        {
          su2Vec[0] += 0.01;
        }
        if (site.coord(0) == xSz/2+1 && site.coord(1) == ySz/2 && site.coord(2) == zSz/2)
        {
          su2Vec[0] += 0.01;
        }
        if (site.coord(0) == xSz/2 && site.coord(1) == ySz/2 && site.coord(2) == zSz/2+1)
        {
          su2Vec[0] += 0.01;
        }
        if (site.coord(0) == xSz/2+1 && site.coord(1) == ySz/2 && site.coord(2) == zSz/2+1)
        {
          su2Vec[0] += 0.01;
        }
      }
      if (ii == 2)
      {
        if (site.coord(0) == xSz/2 && site.coord(1) == ySz/2 && site.coord(2) == zSz/2)
        {
          su2Vec[1] += 0.01;
        }
        if (site.coord(0) == xSz/2+1 && site.coord(1) == ySz/2 && site.coord(2) == zSz/2)
        {
          su2Vec[1] -= 0.01;
        }
        if (site.coord(0) == xSz/2 && site.coord(1) == ySz/2+1 && site.coord(2) == zSz/2)
        {
          su2Vec[1] += 0.01;
        }
        if (site.coord(0) == xSz/2+1 && site.coord(1) == ySz/2+1 && site.coord(2) == zSz/2)
        {
          su2Vec[1] -= 0.01;
        }
      }

      monsta::Matrix su2Mat = monsta::vecToSu2(su2Vec);

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

  // monsta::circShift2(field, tempField, theory, 8, 1, true);
  // monsta::circShift2(tempField, field, theory, 8, 2, true);

  theory.applyBoundaryConditions(field);
  double E = theory.computeEnergy(field);

  monsta::GradDescentSolverNoCheckerboard solver(tol, maxNumSteps, initialStep, maxStepSize, 1000);

  solver.solve(theory, field);

  monsta::writeRawField(field, outputPath + "/rawData");
  monsta::writeCoords(field, outputPath + "/coords");
  monsta::writeHiggsMagnitude(field, outputPath + "/higgsData");
  monsta::writeMagneticField(field, outputPath + "/magneticFieldData", theory);
  monsta::writeEnergyDensity(field, outputPath + "/energyData", theory);
}