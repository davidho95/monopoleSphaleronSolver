#include "LATfield2.hpp"
#include <complex>
#include "../../include/monopoleSphaleron/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../../include/monopoleSphaleron/GeorgiGlashowSu2EomTheory.hpp"
#include "../../include/Matrix.hpp"
#include "../../include/GradDescentSolverBBStep.hpp"
#include "../../include/electroweakSphaleron/GradDescentSolverChigusaElectroweak.hpp"
#include "../../include/Su2Tools.hpp"
#include "../../include/MonopoleFileTools.hpp"
#include "../../include/monopoleSphaleron/MonopoleFieldTools.hpp"
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
  int fluxStep = 1;
  int numIncrements = 1;
  int startB = 0;
  double xAspect = 1;
  double correctionCoeff = 1.5;
  double tanSqMixingAngle = 0.268;

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
      case 'q':
        tanSqMixingAngle = atof(argv[++i]);
        break;
      case 'b':
        correctionCoeff = atof(argv[++i]);
        break;
      case 'S':
        startB = atof(argv[++i]);
        break;
      case 'B':
        fluxStep = atoi(argv[++i]);
        break;
      case 'I':
        numIncrements = atoi(argv[++i]);
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
  LATfield2::Site site(lattice);

  double initialStep = 0.001;
  double maxStepSize = 0.005;
  double tol = 5e-4;
  int maxNumSteps = 200000;
  double abortGrad = 0.5;

  monsta::ElectroweakTheory theory(gaugeCoupling, tanSqMixingAngle, vev, selfCoupling, {0, 0, 0});
  LATfield2::Field<complex<double> > field(lattice, numMatrices, numRows, numCols, 0);
  LATfield2::Field<complex<double> > referenceField(lattice, numMatrices, numRows, numCols, 0);

  monsta::readRawField(field, inputPath + "/rawData");
  theory.applyBoundaryConditions(field);

  int currentB = startB;
  for (int ii = 0; ii < numIncrements; ii++)
  {
    currentB += fluxStep;
    monsta::addConstantMagneticField(field, theory, -fluxStep, 2);
    monsta::GradDescentSolver minimiser(0, maxNumSteps, initialStep, maxStepSize*vev, 100, true);
    monsta::scaleVev(field, theory);
    minimiser.solve(theory, field);

    minimiser = monsta::GradDescentSolver(0, 500/(pow(vev,2)), initialStep, maxStepSize*vev, 500/(pow(vev,2)), true);
    minimiser.solve(theory, field);

    for (site.first(); site.test(); site.next())
    {
      for (int ii = 0; ii < 4; ii++)
      {
        monsta::Matrix gradMat = theory.getLocalGradient(field, site, ii);
        monsta::Matrix fieldMat(field, site, ii);
        if (ii < 3)
        {
          gradMat = gradMat - 0.5*real(monsta::trace(gradMat*conjugateTranspose(fieldMat)))*fieldMat;
        }
        else
        {
          for (int jj = 1; jj < 4; jj++)
          {
            gradMat(jj) = gradMat(jj) - real(gradMat(jj)*conj(fieldMat(jj)))*fieldMat(jj);
          }
          // cout << gradMat(0,0) << endl;
        }
        // cout << real(trace(gradMat*gradMat)) << endl;
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

    monsta::GradDescentSolverChigusa chigusaSolver(tol, maxNumSteps, initialStep, maxStepSize*vev, correctionCoeff, abortGrad);
    bool solved = chigusaSolver.solve(theory, field, referenceField);

    double E = theory.computeEnergy(field);
    if (parallel.rank() == 1)
    {
      ofstream fileStream;
      fileStream.open(outputPath + "/energies.txt", std::ios_base::app);
      fileStream << E << endl;
      fileStream.close();
    }

    std::string bString = std::to_string(currentB);
    bString.erase(bString.find_last_not_of('0') + 1, std::string::npos);
    bString.replace(1,1,"_");

    std::string currentOutputPath = outputPath + "/sphaleronData" + bString;

    mkdir(currentOutputPath.c_str(), 0777);

    monsta::writeRawField(field, currentOutputPath + "/rawData");
    monsta::writeCoords(field, currentOutputPath + "/coords");
    monsta::writeHiggsFieldMagnitude(field, currentOutputPath + "/higgsData");
    monsta::writeMagneticField(field, currentOutputPath + "/magneticFieldData", theory);
    monsta::writeEnergyDensity(field, currentOutputPath + "/energyData", theory);
    monsta::writeUnitaryGaugeField(field, currentOutputPath + "/gaugeData");
  }
}
