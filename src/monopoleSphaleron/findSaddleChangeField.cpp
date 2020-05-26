#include "LATfield2.hpp"
#include <complex>
#include "../include/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../include/GeorgiGlashowSu2EomTheory.hpp"
#include "../include/Matrix.hpp"
#include "../include/GradDescentSolverBBStep.hpp"
#include "../include/GradDescentSolverChigusa.hpp"
#include "../include/GradMinimiser.hpp"
#include "../include/Su2Tools.hpp"
#include "../include/MonopoleFileTools.hpp"
#include "../include/MonopoleFieldTools.hpp"
#include <iostream>
#include <fstream>

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
  double correctionCoeff = 1.5;
  double xAspect = 1;
  double deltaB = 0;

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
        correctionCoeff = atof(argv[++i]);
        break;
      case 'x':
        xAspect = atof(argv[++i]);
        break;
      case 'd':
        deltaB = atoi(argv[++i]);
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
  double maxStepSize = 0.02*gaugeCoupling;
  double tol = 5e-4;
  double abortGrad = 0.1;
  int maxNumSteps = 10000;
  int minGradSteps = 100;

  COUT << "Max step: " << maxStepSize << endl;
  COUT << "Correction coefficient" << correctionCoeff << endl;
  COUT << "Minimum naive gradient descent steps" << minGradSteps << endl;

  monsta::GradMinimiser minimiser(tol, maxNumSteps, initialStep, maxStepSize, minGradSteps);

  monsta::GeorgiGlashowSu2Theory periodicTheory(gaugeCoupling, vev, selfCoupling, {0, 0, 0}, false);
  LATfield2::Field<complex<double> > pairField(lattice, numMatrices, numRows, numCols, 0);
  LATfield2::Field<complex<double> > referenceField(lattice, numMatrices, numRows, numCols, 0);

  monsta::readRawField(pairField, inputPath + "/rawData");
  periodicTheory.applyBoundaryConditions(pairField);

  monsta::scaleVev(pairField, periodicTheory);
  monsta::addConstantMagneticField(pairField, periodicTheory, -1*deltaB);

  minimiser.solve(periodicTheory, pairField);

  for (site.first(); site.test(); site.next())
  {
    for (int ii = 0; ii < 4; ii++)
    {
      monsta::Matrix gradMat = periodicTheory.getLocalGradient(pairField, site, ii);
      monsta::Matrix fieldMat(pairField, site, ii);
      if (ii < 3)
      {
        gradMat = gradMat - 0.5*real(monsta::trace(gradMat*conjugateTranspose(fieldMat)))*fieldMat;
      }
      referenceField(site, ii, 0, 0) = gradMat(0, 0);
      referenceField(site, ii, 0, 1) = gradMat(0, 1);
      referenceField(site, ii, 1, 0) = gradMat(1, 0);
      referenceField(site, ii, 1, 1) = gradMat(1, 1);
    }
  }

  double norm = monsta::innerProduct(referenceField, referenceField);

  for (site.first(); site.test(); site.next())
  {
    for (int ii = 0; ii < referenceField.components(); ii++)
    {
      referenceField(site, ii) = referenceField(site, ii)/ sqrt(norm);
    }
  }


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