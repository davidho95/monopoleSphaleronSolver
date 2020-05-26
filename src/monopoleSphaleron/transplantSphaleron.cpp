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
  double correctionCoeff = 1.5;
  double xAspect = 1;
  int minGradSteps = 50;
  int fluxQuanta = 0;

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
      case 'B':
        fluxQuanta = atof(argv[++i]);
        break;
      case 'x':
        xAspect = atof(argv[++i]);
        break;
      case 'M':
        minGradSteps = atoi(argv[++i]);
        break;
    }
  }

  parallel.initialize(n,m);

  int dim = 3;
  int xSz = round(xAspect*sz);
  int ySz = sz;
  int zSz = sz;
  int biglatSize[dim] = {xSz, ySz, zSz};
  int smallLatSize[dim] = {sz, sz, sz};
  int haloSize = 2;
  int numMatrices = 4;
  int numRows = 2;
  int numCols = 2;

  LATfield2::Lattice bigLattice(dim, biglatSize, haloSize);
  LATfield2::Lattice smallLattice(dim, smallLatSize, haloSize);
  LATfield2::Field<complex<double> > bigField(bigLattice, numMatrices, numRows, numCols, 0);
  LATfield2::Field<complex<double> > smallField(smallLattice, numMatrices, numRows, numCols, 0);

  LATfield2::Site bigSite(bigLattice);
  LATfield2::Site smallSite(smallLattice);

  double initialStep = 0.001;
  double maxStepSize = 0.02*gaugeCoupling;
  double tol = 5e-4;
  double abortGrad = 0.1;
  int maxNumSteps = 10000;

  COUT << "Max step: " << maxStepSize << endl;
  COUT << "Correction coefficient" << correctionCoeff << endl;
  COUT << "Minimum naive gradient descent steps" << minGradSteps << endl;

  monsta::GradMinimiser minimiser(tol, maxNumSteps, initialStep, maxStepSize, minGradSteps);
  // monsta::GradDescentSolver solver(tol, 50, initialStep, maxStepSize, 50);

  monsta::GeorgiGlashowSu2Theory periodicTheory(gaugeCoupling, vev, selfCoupling, {0, 0, 0}, false);
  LATfield2::Field<complex<double> > pairField(bigLattice, numMatrices, numRows, numCols, 0);
  LATfield2::Field<complex<double> > referenceField(bigLattice, numMatrices, numRows, numCols, 0);

  monsta::setVacuumField(bigField, periodicTheory);
  monsta::addConstantMagneticField(bigField, periodicTheory, -fluxQuanta);

  monsta::readRawField(smallField, inputPath + "/rawData");
  periodicTheory.applyBoundaryConditions(smallField);
  monsta::scaleVev(pairField, periodicTheory);

  for(smallSite.first(); smallSite.test(); smallSite.next())
  {
    int smallXCoord = smallSite.coord(0);
    int smallYCoord = smallSite.coord(1);
    int smallZCoord = smallSite.coord(2);

    int bigXCoord = (xSz - sz) / 2 + smallXCoord;
    int bigYCoord = smallYCoord;
    int bigZCoord = smallZCoord;

    bigSite.setCoord(bigXCoord, bigYCoord, bigZCoord);

    for (int ii = 0; ii < bigField.components(); ii++)
    {
      bigField(bigSite, ii) = smallField(smallSite, ii);
    }
  }

  periodicTheory.applyBoundaryConditions(bigField);

  minimiser.solve(periodicTheory, bigField);

  for (bigSite.first(); bigSite.test(); bigSite.next())
  {
    for (int ii = 0; ii < 4; ii++)
    {
      monsta::Matrix gradMat = periodicTheory.getLocalGradient(bigField, bigSite, ii);
      monsta::Matrix fieldMat(bigField, bigSite, ii);
      if (ii < 3)
      {
        gradMat = gradMat - 0.5*real(monsta::trace(gradMat*conjugateTranspose(fieldMat)))*fieldMat;
      }
      referenceField(bigSite, ii, 0, 0) = gradMat(0, 0);
      referenceField(bigSite, ii, 0, 1) = gradMat(0, 1);
      referenceField(bigSite, ii, 1, 0) = gradMat(1, 0);
      referenceField(bigSite, ii, 1, 1) = gradMat(1, 1);
    }
  }

  double norm = monsta::innerProduct(referenceField, referenceField);

  for (bigSite.first(); bigSite.test(); bigSite.next())
  {
    for (int ii = 0; ii < referenceField.components(); ii++)
    {
      referenceField(bigSite, ii) = referenceField(bigSite, ii)/ sqrt(norm);
    }
  }


  monsta::GradDescentSolverChigusa chigusaSolver(tol, maxNumSteps, initialStep, maxStepSize, correctionCoeff, abortGrad);
  chigusaSolver.solve(periodicTheory, bigField, referenceField);

  monsta::GeorgiGlashowSu2EomTheory eomTheory(gaugeCoupling, vev, selfCoupling, {0, 0, 0}, false);
  double gradSq = eomTheory.computeEnergy(bigField);
  COUT << gradSq << endl;

  monsta::writeRawField(bigField, outputPath + "/rawData");
  monsta::writeRawField(referenceField, outputPath + "/referenceRawData");
  monsta::writeCoords(bigField, outputPath + "/coords");
  monsta::writeHiggsFieldUnitary(bigField, outputPath + "/higgsData");
  monsta::writeMagneticField(bigField, outputPath + "/magneticFieldData", periodicTheory);
  monsta::writeEnergyDensity(bigField, outputPath + "/energyData", periodicTheory);
  monsta::writeEnergyDensity(bigField, outputPath + "/gradData", eomTheory);
  monsta::writeUnitaryGaugeField(bigField, outputPath + "/gaugeData");

}