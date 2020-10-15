#include "LATfield2.hpp"
#include <complex>
#include "../../include/monopoleSphaleron/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../../include/monopoleInstanton/GeorgiGlashowSu2Theory4d.hpp"
#include "../../include/Matrix.hpp"
#include "../../include/GradDescentSolverBBStep.hpp"
#include "../../include/GradDescentSolverBBStepNoCheckerboard.hpp"
#include "../../include/Su2Tools.hpp"
#include "../../include/MonopoleFileTools.hpp"
#include "../../include/monopoleInstanton/InstantonFileTools.hpp"
#include "../../include/MonopoleFieldTools.hpp"
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
  int sep = sz/2;
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
      case 'd':
        sep = atoi(argv[++i]);
        break;
      case 'x':
        xAspect = atof(argv[++i]);
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
  double maxStepSize = 0.02;
  double tol = 1e-5;
  int maxNumSteps = 10000;

  monsta::GeorgiGlashowSu2Theory theory(gaugeCoupling, vev, selfCoupling, {2, 0, 0}, true);
  monsta::GeorgiGlashowSu2Theory4d theory4d(gaugeCoupling, vev, selfCoupling);
  monsta::readRawField(field, inputPath + "/rawData");
  theory.applyBoundaryConditions(field);

  monsta::GradDescentSolverNoCheckerboard solver(tol, maxNumSteps, initialStep, maxStepSize, 2000);

  LATfield2::Field<complex<double> > centredField(lattice, numMatrices, numRows, numCols, 0);
  std::vector<int> monopolePos = monsta::findMonopole(field, theory);
  int shiftNum0 = ((sz/2 - 1 - monopolePos[0]) + sz) % sz;
  int shiftNum1 = ((sz/2 - 1 - monopolePos[1]) + sz) % sz;
  int shiftNum2 = ((sz/2 - 1 - monopolePos[2]) + sz) % sz;
  cout << shiftNum0 << endl;
  monsta::circShift2(field, centredField, theory, shiftNum0, 0, true);

  monsta::GeorgiGlashowSu2Theory periodicTheory(gaugeCoupling, vev, selfCoupling, {0, 0, 0}, false);
  LATfield2::Field<complex<double> > pairField(lattice, numMatrices, numRows, numCols, 0);

  monsta::setPairInitialConds2(field, pairField, periodicTheory, sep);
  // monsta::addConstantMagneticField(pairField, periodicTheory, -1);
  monsta::circShift(pairField, field, periodicTheory, -1, 0, false);
  COUT << shiftNum2 << endl;
  monsta::circShift2(field, pairField, periodicTheory, shiftNum1, 1, false);
  monsta::circShift2(pairField, field, periodicTheory, shiftNum2, 2, false);
  monsta::circShift(field, pairField, periodicTheory, 0, 0, false);


  monsta::scaleVev(pairField, periodicTheory);

  LATfield2::Site site2(pairField.lattice());
  // for(site.first(); site.test(); site.next())
  // {
  //   if (site.coord(0) >= xSz / 2) { continue; }
  //   int xCoord = site.coord(0);
  //   int yCoord = site.coord(1);
  //   int zCoord = site.coord(2);

  //   for (int ii = 0; ii < 4; ii++)
  //   {
  //     monsta::Matrix fieldMat(pairField, site, ii);
  //     if (ii == 0) 
  //     {
  //       site2.setCoord(xSz - xCoord, yCoord, zCoord); 
  //     }
  //     else
  //     {
  //       site2.setCoord(xSz - xCoord - 1, yCoord, zCoord);
  //     }
  //     if (ii != 3 ) { fieldMat = monsta::pauli3*fieldMat*monsta::pauli3; }
  //     pairField(site2, ii, 0, 0) = fieldMat(0, 0);
  //     pairField(site2, ii, 0, 1) = fieldMat(0, 1);
  //     pairField(site2, ii, 1, 0) = fieldMat(1, 0);
  //     pairField(site2, ii, 1, 1) = fieldMat(1, 1);
  //   }
  // }

  solver.solve(periodicTheory, pairField);

  monsta::writeRawField(pairField, outputPath + "/rawData");
  monsta::writeCoords(pairField, outputPath + "/coords");
  monsta::writeHiggsMagnitude(pairField, outputPath + "/higgsData");
  monsta::writeMagneticField(pairField, outputPath + "/magneticFieldData", theory4d);
  monsta::writeEnergyDensity(pairField, outputPath + "/energyData", periodicTheory);

}
