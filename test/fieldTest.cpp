#include "LATfield2.hpp"
#include <complex>
#include "../src/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../src/Matrix.hpp"
#include "../src/TheoryChecker.hpp"
#include "../src/GradDescentSolverBBStep.hpp"
#include "../src/Su2Tools.hpp"
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

  srand(1);//time(NULL) + parallel.rank());

  int dim = 3;
  int latSize[dim] = {sz, sz, sz};
  int haloSize = 1;
  int numMatrices = 4;
  int numRows = 2;
  int numCols = 2;

  double gaugeCoupling = 1/(sqrt(5));
  double vev = 1/sqrt(10);
  double selfCoupling = 0.1;

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
      monsta::Matrix su2Mat = monsta::vecToSu2({0, 0.0001, 0.001});
      // monsta::Matrix su2Mat = monsta::vecToSu2({rand() % 10, rand() % 10, rand() % 10});
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
  // double maxStepSize = 0.1 / vev;
  double tol = 1e-4;
  int maxNumSteps = 10000;

  // solver.setVerbosity(false);
  // solver.solve(theory, field);

  // LATfield2::Field<complex<double> > centredField(lattice, numMatrices, numRows, numCols, 0);
  // std::vector<int> monopolePos = monsta::findMonopoleUnitary(field);
  // int shiftNum = ((sz/2 - 1 - monopolePos[0]) + sz) % sz;
  // monsta::circShift(field, centredField, theory, shiftNum, 0, true);

  std::vector<double> mVals = {};
  for (int ii = 0; ii < 21; ii++)
  {
    mVals.push_back(0.1*(1 - ii / 20.0));
  }
  // for (int ii = 0; ii < 21; ii++)
  // {
  //   mVals.push_back(ii/200.0);
  // }
  
  // std::vector<double> betaVals = {0.01, 0.05, 0.1, 0.5, 1, 2, 3, 4, 5, 6, 7};
  std::vector<double> betaVals = {1};// {7, 6, 5, 4, 3, 2, 1, 0.5, 0.1, 0.05, 0.01};
  for (int ii = 0; ii < mVals.size(); ii++)
  {
    if (ii > 0) { break; }
    // double vev = sqrt(mVals[ii] / (2*selfCoupling));
    // // double selfCoupling = 0.1;//.001;//pow(betaVals[ii],2)/8;
    // COUT << selfCoupling << endl;
    double maxStepSize = selfCoupling >= 1 ? 0.1/sqrt(selfCoupling) : 0.1;
    // if (ii > 0) { break; }
    monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize);
    monsta::GeorgiGlashowSu2Theory theory(gaugeCoupling, vev, selfCoupling, {-1, -1, -1}, true);
    LATfield2::Field<complex<double> > newField(lattice, numMatrices, numRows, numCols, 0);
    monsta::circShift(field, newField, theory, 0, 0, true);

    // monsta::TheoryChecker checker(tol, true);
    // checker.checkGradients(theory, field);

    solver.solve(theory, newField);
    double E = theory.computeEnergy(newField);

    if (parallel.rank() == 1)
    {
      ofstream fileStream;
      fileStream.open(outputPath + "/energies.txt", std::ios_base::app);
      fileStream << E << endl;
      fileStream.close();
    }

    cout << newField(site, 3, 0, 0) << endl;

    monsta::circShift(newField, field, theory, 0, 0, true);

    monsta::writeCoords(field, outputPath + "/coords");
    monsta::writeHiggsFieldUnitary(field, outputPath + "/higgsData");
    monsta::writeMagneticField(field, outputPath + "/magneticFieldData", theory);
    monsta::writeEnergyDensity(field, outputPath + "/energyData", theory);

  }


  clock_t end = clock();
  COUT << double(end - begin) / CLOCKS_PER_SEC << endl;
}