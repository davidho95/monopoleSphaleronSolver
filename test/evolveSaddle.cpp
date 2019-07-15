#include "LATfield2.hpp"
#include <complex>
#include "../src/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../src/GeorgiGlashowSu2EomTheory.hpp"
#include "../src/Matrix.hpp"
#include "../src/GradDescentSolverBBStep.hpp"
#include "../src/GradDescentSolverChigusa.hpp"
#include "../src/GradMinimiser.hpp"
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
  double vev = 1;
  double gaugeCoupling = 1;
  double selfCoupling = 1;
  int sep = sz/2;

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
    }
  }

  parallel.initialize(n,m);

  int dim = 3;
  int latSize[dim] = {sz, sz, sz};
  int haloSize = 2;
  int numMatrices = 4;
  int numRows = 2;
  int numCols = 2;

  LATfield2::Lattice lattice(dim, latSize, haloSize);
  LATfield2::Field<complex<double> > field(lattice, numMatrices, numRows, numCols, 0);

  LATfield2::Site site(lattice);

  double initialStep = 0.001;
  double maxStepSize = 0.05*vev*gaugeCoupling;
  double tol = 5e-4;
  double abortGrad = 0.1;
  int maxNumSteps = 50000;
  double correctionCoeff = 1.3;

  COUT << "Max step: " << maxStepSize << endl;
  COUT << "Correction coefficient" << correctionCoeff << endl;

  monsta::GeorgiGlashowSu2Theory periodicTheory(gaugeCoupling, vev, selfCoupling, {0, 0, 0}, false);
  LATfield2::Field<complex<double> > pairField(lattice, numMatrices, numRows, numCols, 0);
  LATfield2::Field<complex<double> > referenceField(lattice, numMatrices, numRows, numCols, 0);

  monsta::readRawField(pairField, inputPath + "/rawData");
  periodicTheory.applyBoundaryConditions(pairField);

  for (site.first(); site.test(); site.next())
  {
    for (int ii = 0; ii < numMatrices; ii++)
    {
      periodicTheory.postProcess(field, site, ii);
    }
  }

  monsta::readRawField(referenceField, inputPath + "/referenceRawData");

  monsta::GradDescentSolverChigusa chigusaSolver(tol, maxNumSteps, initialStep, maxStepSize, correctionCoeff, abortGrad);
  chigusaSolver.solve(periodicTheory, pairField, referenceField);

  monsta::GeorgiGlashowSu2EomTheory eomTheory(gaugeCoupling, vev, selfCoupling, {0, 0, 0}, false);
  double gradSq = eomTheory.computeEnergy(pairField);
  COUT << gradSq << endl;

  monsta::writeRawField(pairField, outputPath + "/rawData");
  monsta::writeRawField(referenceField, outputPath + "/referenceRawData");
  monsta::writeCoords(pairField, outputPath + "/coords");
  monsta::writeHiggsFieldUnitary(pairField, outputPath + "/higgsData");
  monsta::writeMagneticField(pairField, outputPath + "/magneticFieldData", periodicTheory);
  monsta::writeEnergyDensity(pairField, outputPath + "/energyData", periodicTheory);
  monsta::writeEnergyDensity(pairField, outputPath + "/gradData", eomTheory);
  monsta::writeUnitaryGaugeField(pairField, outputPath + "/gaugeData");

}