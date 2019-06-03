#include "LATfield2.hpp"
#include <complex>
#include "../src/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../src/Matrix.hpp"
#include "../src/GradDescentSolverBBStepNoCheckerboard.hpp"
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
  std::string inputPath;
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
    }
  }

  parallel.initialize(n,m);

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

  LATfield2::Site site(lattice);

  double initialStep = 0.001;
  double maxStepSize = 0.05/vev;
  double tol = 1e-4;
  int maxNumSteps = 10000;

  monsta::GeorgiGlashowSu2Theory theory(gaugeCoupling, vev, selfCoupling, {2, 0, 0}, true);
  monsta::readRawField(field, inputPath + "/rawData");
  theory.applyBoundaryConditions(field);

  monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize);
  solver.solve(theory, field);

  LATfield2::Field<complex<double> > centredField(lattice, numMatrices, numRows, numCols, 0);
  std::vector<int> monopolePos = monsta::findMonopoleUnitary(field);
  int shiftNum = ((sz/2 - 1 - monopolePos[0]) + sz) % sz;
  monsta::circShift(field, centredField, theory, shiftNum, 0, true);

  monsta::GeorgiGlashowSu2Theory periodicTheory(gaugeCoupling, vev, selfCoupling, {0, 0, 0}, true);
  LATfield2::Field<complex<double> > pairField(lattice, numMatrices, numRows, numCols, 0);

  for (int sep = 0; sep < sz - 1; sep++)
  {
    monsta::setPairInitialConds(field, pairField, periodicTheory, sep);
    // monsta::addConstantMagneticField(pairField, periodicTheory,-1);
    solver.solve(periodicTheory, pairField);
    double E = periodicTheory.computeEnergy(pairField);
    if (parallel.rank() == 1)
    {
      ofstream fileStream;
      fileStream.open(outputPath + "/energies.txt", std::ios_base::app);
      fileStream << E << endl;
      fileStream.close();
    }

    if (sep == sz/2)
    {
      monsta::writeCoords(pairField, outputPath + "/coords");
      monsta::writeHiggsFieldUnitary(pairField, outputPath + "/higgsData");
      monsta::writeMagneticField(pairField, outputPath + "/magneticFieldData", periodicTheory);
      monsta::writeEnergyDensity(pairField, outputPath + "/energyData", periodicTheory);
      monsta::writeUnitaryGaugeField(pairField, outputPath + "/gaugeData");
    }
  }
}