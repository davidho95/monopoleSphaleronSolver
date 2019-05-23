#include "LATfield2.hpp"
#include <complex>
#include "../src/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../src/GeorgiGlashowSu2EomTheory.hpp"
#include "../src/Matrix.hpp"
#include "../src/GradDescentSolverBBStepNoCheckerboard.hpp"
#include "../src/Su2Tools.hpp"
#include "../src/TheoryChecker.hpp"
#include "../src/MonopoleFileTools.hpp"
#include "../src/MonopoleFieldTools.hpp"
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

  srand(time(NULL) + parallel.rank());

  int dim = 3;
  int latSize[dim] = {sz, sz, sz};
  int haloSize = 2;
  int numMatrices = 4;
  int numRows = 2;
  int numCols = 2;

  double gaugeCoupling = 1;
  double vev = 1;
  double selfCoupling = 1;

  LATfield2::Lattice lattice(dim, latSize, haloSize);
  LATfield2::Field<complex<double> > field(lattice, numMatrices, numRows, numCols, 0);
  // field.updateHalo();

  LATfield2::Site site(lattice);

  double initialStep = 0.001;
  double maxStepSize = 0.05/vev;
  double tol = 1e-6;
  int maxNumSteps = 10000;

  monsta::GeorgiGlashowSu2Theory theory(gaugeCoupling, vev, selfCoupling, {2, 0, 0}, true);
  monsta::setSingleMonopoleInitialConditions(field, theory);
  // monsta::readRawField(field, outputPath + "/singleMonopoleData/rawData");
  // theory.applyBoundaryConditions(field);


  // monsta::TheoryChecker checker(tol);
  // checker.checkGradients(eomTheory, field);

  // site.setCoord(0,0,0);
  // cout << eomTheory.getLocalEnergyDensity(field, site) << endl;

  monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize);
  solver.solve(theory, field);

  // cout << eomTheory.computeEnergy(field) << endl;
  // double E = theory.computeEnergy(field);
  // COUT << E << endl;

  LATfield2::Field<complex<double> > centredField(lattice, numMatrices, numRows, numCols, 0);
  std::vector<int> monopolePos = monsta::findMonopoleUnitary(field);
  int shiftNum = ((sz/2 - 1 - monopolePos[0]) + sz) % sz;
  monsta::circShift(field, centredField, theory, shiftNum, 0, true);

  monsta::GeorgiGlashowSu2Theory periodicTheory(gaugeCoupling, vev, selfCoupling, {0, 0, 0}, true);
  monsta::GeorgiGlashowSu2EomTheory eomTheory(gaugeCoupling, vev, selfCoupling, {0, 0, 0}, true);
  LATfield2::Field<complex<double> > pairField(lattice, numMatrices, numRows, numCols, 0);
  // monsta::readRawField(pairField, outputPath + "/rawData");

  // for (int sep = 2; sep < sz - 2; sep++)
  // {
  //   monsta::setPairInitialConds(field, pairField, periodicTheory, sep);
  //   solver.solve(periodicTheory, pairField);
  //   monsta::oneFlux(pairField, periodicTheory);
  //   solver.solve(periodicTheory, pairField);
  //   double E = periodicTheory.computeEnergy(pairField);
  //   if (parallel.rank() == 1)
  //   {
  //     ofstream fileStream;
  //     fileStream.open(outputPath + "/energies.txt", std::ios_base::app);
  //     fileStream << E << endl;
  //     fileStream.close();
  //   }
  // }

  monsta::setPairInitialConds(field, pairField, periodicTheory, 4);
  // monsta::oneFlux(pairField, periodicTheory);
  // monsta::readRawField(pairField, outputPath + "/rawData");
  periodicTheory.applyBoundaryConditions(pairField);

  solver.setParams(tol, maxNumSteps, initialStep, maxStepSize);
  // solver.setParams(1e-4, maxNumSteps, 0.0005, 0.005);
  // monsta::scaleVev(pairField, periodicTheory);
  solver.solve(periodicTheory, pairField);


  double gradSum = eomTheory.computeEnergy(pairField);
  COUT << gradSum << endl;


  // double gradSq = eomTheory.computeEnergy(pairField);
  // double gradSqOld = 1e6;

  // solver.setParams(0.01, 1, initialStep, maxStepSize);

  // while (gradSq < gradSqOld)
  // {
  //   gradSqOld = gradSq;
  //   solver.solve(periodicTheory, pairField);
  //   gradSq = eomTheory.computeEnergy(pairField);
  //   COUT << gradSq << endl;
  // }

  // solver.setParams(1e-6, maxNumSteps, maxStepSize, initialStep);
  // solver.solve(periodicTheory, pairField);

  // eomTheory.applyBoundaryConditions(field);
  // double gradSum = eomTheory.computeEnergy(field); 
  // COUT << gradSum << endl;

  int minSteps = 10;
  solver.setParams(tol, maxNumSteps, initialStep, maxStepSize);

  // periodicTheory.applyBoundaryConditions(pairField);
  // solver.solve(periodicTheory, pairField);
  // double gradSum = eomTheory.computeEnergy(pairField); 
  // COUT << gradSum << endl;

  // for(site.first(); site.test(); site.next())
  // {
  //   for (int ii = 0; ii < 3; ii++)
  //   {
  //     field(site, ii, 0, 0) = 1;
  //     field(site, ii, 1, 1) = 1;
  //   }
  //   field(site, 3, 0, 0) = vev/sqrt(2);
  //   // COUT << real(selfCoupling*pow(2.0*pow(field(site, 3, 0, 0),2) - pow(vev, 2),2)) << endl;
  // }
  // for(site.first(); site.test(); site.next())
  // {
  //   // cout << periodicTheory.getLocalEnergyDensity(field, site) << endl;
  // }

  double pi = 4*std::atan(1);
  double extFieldStrength = 2*pi/pow(sz,1);
  // periodicTheory.applyBoundaryConditions(field);
  // monsta::setConstantMagneticFieldUnitary(pairField, periodicTheory, 0, 0);
  // monsta::addMagneticFieldUnitary(pairField, periodicTheory, extFieldStrength, 0);
  // monsta::oneFlux(pairField, periodicTheory);
  // // double E = eomTheory.computeEnergy(field);
  // // COUT << E << endl;
  // solver.solve(periodicTheory, pairField);
  // solver.solve(periodicTheory, pairField);

  // solver.setParams(1e-4, maxNumSteps, 0.0001, 0.001, minSteps);
  // solver.solve(eomTheory, pairField);

  monsta::writeCoords(pairField, outputPath + "/coords");
  monsta::writeHiggsFieldUnitary(pairField, outputPath + "/higgsData");
  monsta::writeMagneticField(pairField, outputPath + "/magneticFieldData", periodicTheory);
  monsta::writeEnergyDensity(pairField, outputPath + "/energyData", periodicTheory);
  monsta::writeUnitaryGaugeField(pairField, outputPath + "/gaugeData");
  // monsta::writeRawField(pairField, outputPath + "/rawData");

  // double E = periodicTheory.computeEnergy(pairField);
  // COUT << E << endl;

  clock_t end = clock();
  // COUT << double(end - begin) / CLOCKS_PER_SEC << endl;
}