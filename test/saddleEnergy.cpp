#include "LATfield2.hpp"
#include <complex>
#include "../src/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../src/GeorgiGlashowSu2EomTheory.hpp"
#include "../src/Matrix.hpp"
#include "../src/GradDescentSolverBBStep.hpp"
#include "../src/GradDescentSolverChigusa.hpp"
#include "../src/GradDescentSolverCorePreserving.hpp"
#include "../src/GradMinimiser.hpp"
#include "../src/Su2Tools.hpp"
#include "../src/MonopoleFileTools.hpp"
#include "../src/MonopoleFieldTools.hpp"
#include <iostream>
#include <fstream>
#include <ctime>
#include <bits/stdc++.h> 
#include <iostream> 
#include <sys/stat.h> 
#include <sys/types.h> 


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
  double startVev;
  double endVev;
  double numIncrements;
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
      case 'S':
        startVev = atof(argv[++i]);
        break;
      case 'E':
        endVev = atof(argv[++i]);
        break;
      case 'I':
        numIncrements = atof(argv[++i]);
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
  double maxStepSize = 0.05;
  double tol = 5e-4;
  int maxNumSteps = 50000;
  double abortGrad = 1;
  double correctionCoeff = 1.3;

  monsta::GeorgiGlashowSu2Theory periodicTheory(gaugeCoupling, vev, selfCoupling, {0, 0, 0}, false);
  LATfield2::Field<complex<double> > pairField(lattice, numMatrices, numRows, numCols, 0);
  LATfield2::Field<complex<double> > referenceField(lattice, numMatrices, numRows, numCols, 0);

  monsta::readRawField(pairField, inputPath + "/rawData");
  periodicTheory.applyBoundaryConditions(pairField);

  std::vector<double> vevs;
  for (int ii = 0; ii < numIncrements; ii++)
  {
    vevs.push_back(endVev - (endVev - startVev)*double(ii)/(numIncrements - 1));
    // COUT << vevs[ii] << endl;
  }

  for (int ii = 0; ii < numIncrements; ii++)
  {
    vev = vevs[ii];
    gaugeCoupling = 0.5;
    selfCoupling = pow(gaugeCoupling, 2) / 2;
    
    // COUT << selfCoupling << endl;

    monsta::GradMinimiser minimiser(tol, maxNumSteps, initialStep, maxStepSize*vev*gaugeCoupling, 50);
    // monsta::GradDescentSolver solver(tol, 50, initialStep, maxStepSize*vev*gaugeCoupling, 50);

    monsta::GeorgiGlashowSu2Theory theory(gaugeCoupling, vev, selfCoupling, {0, 0, 0}, false);
    monsta::GeorgiGlashowSu2EomTheory eomTheory(gaugeCoupling, vev, selfCoupling, {0, 0, 0}, false);

    monsta::scaleVev(pairField, theory);
    minimiser.solve(theory, pairField);
    // solver.solve(theory, pairField);

    for (site.first(); site.test(); site.next())
    {
      for (int ii = 0; ii < 4; ii++)
      {
        monsta::Matrix gradMat = theory.getLocalGradient(pairField, site, ii);
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

    monsta::GradDescentSolverChigusa chigusaSolver(tol, maxNumSteps, initialStep, maxStepSize*vev*gaugeCoupling, correctionCoeff, abortGrad);
    bool solved = chigusaSolver.solve(theory, pairField, referenceField);

    double E = theory.computeEnergy(pairField);
    if (parallel.rank() == 1)
    {
      ofstream fileStream;
      fileStream.open(outputPath + "/energies1.txt", std::ios_base::app);
      fileStream << E << endl;
      fileStream.close();
    }

    std::string vevString = std::to_string(vev);
    vevString.erase(vevString.find_last_not_of('0') + 1, std::string::npos);
    vevString.replace(1,1,"_");

    std::string currentOutputPath = outputPath + "/saddleData" + vevString;

    mkdir(currentOutputPath.c_str(), 0777);

    monsta::writeRawField(pairField, currentOutputPath + "/rawData");
    monsta::writeCoords(pairField, currentOutputPath + "/coords");
    monsta::writeHiggsFieldUnitary(pairField, currentOutputPath + "/higgsData");
    monsta::writeMagneticField(pairField, currentOutputPath + "/magneticFieldData", theory);
    monsta::writeEnergyDensity(pairField, currentOutputPath + "/energyData", theory);
    monsta::writeEnergyDensity(pairField, currentOutputPath + "/gradData", eomTheory);
    monsta::writeUnitaryGaugeField(pairField, currentOutputPath + "/gaugeData");
  }
}