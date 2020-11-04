#include "LATfield2.hpp"
#include <complex>
#include "../../include/Matrix.hpp"
#include "../../include/GradDescentSolverBBStep.hpp"
#include "../../include/GradDescentSolverBBStepNoCheckerboard.hpp"
#include "../../include/monopoleSphaleron/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../../include/monopoleInstanton/GradDescentSolverChigusa4d.hpp"
#include "../../include/Su2Tools.hpp"
#include "../../include/Theory.hpp"
#include "../../include/MonopoleFileTools.hpp"
#include "../../include/monopoleInstanton/GeorgiGlashowSu2Theory4d.hpp"
// #include "../../include/MonopoleFieldTools.hpp"
#include "../../include/monopoleInstanton/InstantonFileTools.hpp"
#include "../../include/monopoleInstanton/InstantonFieldTools.hpp"
#include <iostream>
#include <fstream>
int main(int argc, char **argv)
{
  clock_t begin = clock();

  int sz = 16;
  int n = 2;
  int m = 2;
  double xAspect = 1;
  double vev = 1;
  double gaugeCoupling = 1;
  double selfCoupling = 1;
  int fluxQuanta = 0;
  std::string inputPath;
  std::string outputPath;
  int extraSteps = 100;
  double beta = 1.5;
  double tol = 0.001;
  double maxNumSteps = 500000;
  int contract = 0;

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
      case 'x':
        xAspect = atof(argv[++i]);
        break;
      case 'B':
        fluxQuanta = atoi(argv[++i]);
        break;
      case 'b':
        beta = atof(argv[++i]);
        break;
      case 'e':
        extraSteps = atoi(argv[++i]);
        break;
      case 't':
        tol = atof(argv[++i]);
        break;
      case 'M':
        maxNumSteps = atoi(argv[++i]);
        break;
      case 'c':
        contract = atoi(argv[++i]);
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

  double initialStep = 0.0001;
  double maxStepSize = 0.001;
  double abortGrad = 0.01;

  monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize, 100, true);

  monsta::GeorgiGlashowRadialTheory theory(gaugeCoupling, vev, selfCoupling);
  LATfield2::Field<complex<double> > referenceField(bigLattice, numMatrices, numRows, numCols, 0);

  monsta::setVacuumField(bigField, theory);
  monsta::addConstantMagneticField(bigField, theory, -fluxQuanta);

  monsta::readRawField(smallField, inputPath + "/rawData");
  for(smallSite.first(); smallSite.test(); smallSite.next())
  {
    for (int ii = 0; ii < 4; ii++)
    {
      theory.postProcess(smallField, smallSite, ii);
    }
  }
  theory.applyBoundaryConditions(smallField);

  for(smallSite.first(); smallSite.test(); smallSite.next())
  {
    int smallXCoord = smallSite.coord(0);
    int smallYCoord = smallSite.coord(1);
    int smallZCoord = smallSite.coord(2);

    int bigXCoord = smallXCoord;
    int bigYCoord = smallYCoord;
    int bigZCoord = smallZCoord;

    bigSite.setCoord(bigXCoord, bigYCoord, bigZCoord);

    for (int ii = 0; ii < bigField.components(); ii++)
    {
      bigField(bigSite, ii) = smallField(smallSite, ii);
    }
  }

  theory.applyBoundaryConditions(bigField);

  solver.solve(theory, bigField);

  solver = monsta::GradDescentSolver(tol, extraSteps, initialStep, maxStepSize, extraSteps);
  solver.solve(theory, bigField);

  for (bigSite.first(); bigSite.test(); bigSite.next())
  {
    for (int ii = 0; ii < 4; ii++)
    {
      monsta::Matrix gradMat = theory.getLocalGradient(bigField, bigSite, ii);
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

  bool solved = false;

  monsta::GradDescentSolverChigusa chigusaSolver(tol, maxNumSteps, initialStep, maxStepSize, beta, 0.05);
  solved = chigusaSolver.solve(theory, bigField, referenceField);

  monsta::writeRawField(bigField, outputPath + "/rawData");
  monsta::writeCoords(bigField, outputPath + "/coords");
  monsta::writeEnergyDensity(bigField, outputPath + "/energyData", theory);
  monsta::writeHiggsMagnitude(bigField, outputPath + "/higgsData");
  monsta::writeMagneticField(bigField, outputPath + "/magneticFieldData", theory);
  monsta::writeGradients(bigField, outputPath + "/gradData", theory);
  monsta::writeRawField(referenceField, outputPath + "/referenceRawData");

}