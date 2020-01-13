#include "LATfield2.hpp"
#include <complex>
#include "../src/ElectroweakTheory.hpp"
#include "../src/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../src/Matrix.hpp"
#include "../src/TheoryChecker.hpp"
#include "../src/GradDescentSolverBBStep.hpp"
#include "../src/GradDescentSolverChigusaElectroweak.hpp"
#include "../src/Su2Tools.hpp"
#include "../src/MonopoleFileTools.hpp"
#include "../src/ElectroweakFieldTools.hpp"
#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{
  int sz = 16;
  int n = 2;
  int m = 2;
  double zAspect = 1;
  double vev = 1;
  double gaugeCoupling = 1;
  double mixingAngle = 2*atan(1) / 3;
  double selfCoupling = 1;
  int fluxQuanta = 0;
  std::string ambjornOlesenInputPath;
  std::string sphaleronInputPath;
  std::string outputPath;
  int extraSteps = 1000;

  for (int i=1 ; i < argc ; i++ ){
    if ( argv[i][0] != '-' )
      continue;
    switch(argv[i][1]) {
      case 'p':
        outputPath = argv[++i];
        break;
      case 'i':
        ambjornOlesenInputPath = argv[++i];
        break;
      case 'j':
        sphaleronInputPath = argv[++i];
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
      case 'q':
        mixingAngle = atof(argv[++i]);
        break;
      case 'z':
        zAspect = atof(argv[++i]);
        break;
      case 'b':
        fluxQuanta = atoi(argv[++i]);
        break; 
      case 'e':
        extraSteps = atoi(argv[++i]);
        break;
    }
  }
  parallel.initialize(n,m);

  srand(1);

  int dim = 3;
  int xSz = sz;
  int ySz = sz;
  int zSz = round(zAspect*sz);
  int latSize[dim] = {xSz, ySz, zSz};
  int haloSize = 1;
  int numMatrices = 4;
  int numRows = 2;
  int numCols = 2;


  LATfield2::Lattice lattice(dim, latSize, haloSize);
  LATfield2::Field<complex<double> > ambjornOlesenField(lattice, numMatrices, numRows, numCols, 0);
  LATfield2::Field<complex<double> > sphaleronField(lattice, numMatrices, numRows, numCols, 0);
  LATfield2::Field<complex<double> > combinedField(lattice, numMatrices, numRows, numCols, 0);
  monsta::ElectroweakTheory theory(gaugeCoupling, mixingAngle, vev, selfCoupling, {0,0,0});
  // monsta::GeorgiGlashowSu2Theory ggTheory(gaugeCoupling, vev, selfCoupling);

  double tol = 1e-6;
  double initialStep = 0.001;
  double maxStepSize = 0.005;
  int maxNumSteps = 100000;

  LATfield2::Site site(lattice);

  monsta::readRawField(ambjornOlesenField, ambjornOlesenInputPath + "/rawData");
  theory.applyBoundaryConditions(ambjornOlesenField);
  monsta::readRawField(sphaleronField, sphaleronInputPath + "/rawData");
  theory.applyBoundaryConditions(sphaleronField);

  // monsta::addConstantMagneticField(sphaleronField, theory, -fluxQuanta, 2);

  // monsta::linearSuperpose(ambjornOlesenField, sphaleronField, combinedField, theory);
  // theory.applyBoundaryConditions(combinedField);

  for (site.first(); site.test(); site.next())
  {
    std::vector<double> coords(3);
    coords[0] = site.coord(0) - xSz/2 + 0.5;
    coords[1] = site.coord(1) - ySz/2 + 0.5;
    coords[2] = site.coord(2) - zSz/2 + 0.5;

    for (int ii = 0; ii < numMatrices; ii++)
    {
      if (abs(coords[0]) < 6 && abs(coords[1]) < 6)
      {
        combinedField(site, ii, 0, 0) = sphaleronField(site, ii, 0, 0);
        combinedField(site, ii, 0, 1) = sphaleronField(site, ii, 0, 1);
        combinedField(site, ii, 1, 0) = sphaleronField(site, ii, 1, 0);
        combinedField(site, ii, 1, 1) = sphaleronField(site, ii, 1, 1);
      }
      else
      {
        combinedField(site, ii, 0, 0) = ambjornOlesenField(site, ii, 0, 0);
        combinedField(site, ii, 0, 1) = ambjornOlesenField(site, ii, 0, 1);
        combinedField(site, ii, 1, 0) = ambjornOlesenField(site, ii, 1, 0);
        combinedField(site, ii, 1, 1) = ambjornOlesenField(site, ii, 1, 1);
      }
    }
  }
  theory.applyBoundaryConditions(combinedField);


  double E = theory.computeEnergy(combinedField);
  COUT << E << endl;

  monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize, 100, true);
  solver.solve(theory, combinedField);

  solver = monsta::GradDescentSolver(tol, extraSteps, initialStep, maxStepSize, extraSteps, true);
  solver.solve(theory, combinedField);

  LATfield2::Field<complex<double> > referenceField(lattice, numMatrices, numRows, numCols, 0);

  for (site.first(); site.test(); site.next())
  {
    for (int ii = 0; ii < 4; ii++)
    {
      monsta::Matrix gradMat = theory.getLocalGradient(combinedField, site, ii);
      monsta::Matrix fieldMat(combinedField, site, ii);
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
  // cout << norm << endl;

  for (site.first(); site.test(); site.next())
  {
    for (int ii = 0; ii < referenceField.components(); ii++)
    {
      referenceField(site, ii) = referenceField(site, ii)/ sqrt(norm);
    }
  }


  monsta::GradDescentSolverChigusa chigusaSolver(5e-4, 333000, initialStep, maxStepSize, 1.2, 0.5);
  chigusaSolver.solve(theory, combinedField, referenceField);

  monsta::writeCoords(combinedField, outputPath + "/coords");
  monsta::writeEnergyDensity(combinedField, outputPath + "/energyData", theory);
  monsta::writeHiggsField(combinedField, outputPath + "/higgsData", theory);
  monsta::writeMagneticField(combinedField, outputPath + "/magneticFieldData", theory);
  monsta::writeRawField(combinedField, outputPath + "/rawData");
  monsta::writeRawField(referenceField, outputPath + "/referenceRawData");

}
