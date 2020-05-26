#include "LATfield2.hpp"
#include <complex>
#include "../src/ElectroweakTheory.hpp"
#include "../src/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../src/Matrix.hpp"
#include "../src/TheoryChecker.hpp"
#include "../src/GradDescentSolverBBStep.hpp"
#include "../src/GradDescentSolverChigusaElectroweak.hpp"
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
        break;; 
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
  monsta::ElectroweakTheory theory(gaugeCoupling, mixingAngle, vev, selfCoupling, {0,0,0});

  monsta::setVacuumField(field, theory);
  LATfield2::Site site(lattice);

  if (inputPath != "")
  {
    monsta::readRawField(field, inputPath + "/rawData");
    theory.applyBoundaryConditions(field);
  }

  monsta::writeCoords(field, outputPath + "/coords");
  monsta::writeEnergyDensity(field, outputPath + "/energyData", theory);
  monsta::writeHiggsField(field, outputPath + "/higgsData", theory);
  monsta::writeMagneticField(field, outputPath + "/magneticFieldData", theory);
}