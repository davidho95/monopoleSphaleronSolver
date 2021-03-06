#include "LATfield2.hpp"
#include <complex>
#include "../include/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../include/Matrix.hpp"
#include "../include/GradDescentSolverBBStep.hpp"
#include "../include/Su2Tools.hpp"
#include "../include/MonopoleFileTools.hpp"
#include "../include/MonopoleFieldTools.hpp"
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
  int fluxQuanta = 0;
  int range = 1;
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
        fluxQuanta = atof(argv[++i]);
        break;
      case 'r':
        range = atoi(argv[++i]);
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

  double initialStep = 0.001;
  double maxStepSize = 0.05/vev;
  double tol = 1e-4;
  int maxNumSteps = 10000;

  monsta::GeorgiGlashowSu2Theory theory(gaugeCoupling, vev, selfCoupling, {2, 0, 0}, true);
  monsta::readRawField(field, inputPath + "/rawData");
  theory.applyBoundaryConditions(field);

  monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize);
  // solver.solve(theory, field);

  LATfield2::Field<complex<double> > centredField(lattice, numMatrices, numRows, numCols, 0);
  std::vector<int> monopolePos = monsta::findMonopole(field, theory);
  int shiftNum = ((sz/2 - 1 - monopolePos[0]) + sz) % sz;
  monsta::circShift(field, centredField, theory, shiftNum, 0, true);

  monsta::GeorgiGlashowSu2Theory periodicTheory(gaugeCoupling, vev, selfCoupling, {0, 0, 0}, false);
  LATfield2::Field<complex<double> > pairField(lattice, numMatrices, numRows, numCols, 0);

  double pi = 4*atan(1);
  double predictedMax = round(double(sz) / sqrt(4*pi*(fluxQuanta - 0.5)));

  double oldEnergy = 0;
  bool maxReached = false;
  for (int sep = predictedMax - range; sep <= predictedMax + range; sep++)
  {
    if (sep < 1 || sep > xSz - 1) { continue; }
    monsta::setPairInitialConds2(field, pairField, periodicTheory, sep);
    monsta::addConstantMagneticField(pairField, periodicTheory, -fluxQuanta);
    solver.solve(periodicTheory, pairField);
    double energy = periodicTheory.computeEnergy(pairField);
    if (parallel.rank() == 1)
    {
      ofstream fileStream;
      fileStream.open(outputPath + "/energies.txt", std::ios_base::app);
      fileStream << energy << endl;
      fileStream.close();
    }

    if (energy < oldEnergy) { maxReached = true; }
    oldEnergy = energy;

    if (!maxReached)
    {
      monsta::writeRawField(pairField, outputPath + "/rawData");
      monsta::writeCoords(pairField, outputPath + "/coords");
      monsta::writeHiggsFieldUnitary(pairField, outputPath + "/higgsData");
      monsta::writeMagneticField(pairField, outputPath + "/magneticFieldData", periodicTheory);
      monsta::writeEnergyDensity(pairField, outputPath + "/energyData", periodicTheory);
      monsta::writeUnitaryGaugeField(pairField, outputPath + "/gaugeData");
    }
  }
}