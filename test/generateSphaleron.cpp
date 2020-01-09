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
  std::string inputPath;
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
  LATfield2::Field<complex<double> > field(lattice, numMatrices, numRows, numCols, 0);
  monsta::ElectroweakTheory theory(gaugeCoupling, mixingAngle, vev, selfCoupling, {0,0,0});
  // monsta::GeorgiGlashowSu2Theory ggTheory(gaugeCoupling, vev, selfCoupling);

  monsta::setVacuumField(field, theory);
  double tol = 1e-6;
  double initialStep = 0.001;
  double maxStepSize = 0.005;
  int maxNumSteps = 100000;

  LATfield2::Site site(lattice);

  if (inputPath != "")
  {
    monsta::readRawField(field, inputPath + "/rawData");
    theory.applyBoundaryConditions(field);
  }
  else
  {
    for(site.first(); site.test(); site.next())
    {
      std::vector<double> coords(3);
      coords[0] = site.coord(0) - xSz/2 + 0.5;
      coords[1] = site.coord(1) - ySz/2 + 0.5;
      coords[2] = site.coord(2) - zSz/2 + 0.5;
      double r = sqrt(pow(coords[0], 2) + pow(coords[1], 2) + pow(coords[2], 2) + 0.01);
      std::vector<double> gaugeFns(3);
      gaugeFns[0] = 1/cosh(vev*gaugeCoupling*r/3);
      gaugeFns[1] = 1/cosh(vev*gaugeCoupling*r/3);
      gaugeFns[2] = 1/cosh(vev*gaugeCoupling*r/3);
      double higgsFn = tanh(vev*gaugeCoupling*r/3);
      for (int ii = 0; ii < 3; ii++)
      {
        std::vector<double> su2Vec(3);
        int dir1 = (ii + 1) % 3;
        int dir2 = (ii + 2) % 3;
        su2Vec[dir1] = (1/pow(r, 2) * coords[dir2])*gaugeFns[ii];
        su2Vec[dir2] = (-1/pow(r, 2) * coords[dir1])*gaugeFns[ii];
        monsta::Matrix su2Mat = monsta::vecToSu2(su2Vec);

        field(site, ii, 0, 0) = su2Mat(0, 0);
        field(site, ii, 0, 1) = su2Mat(0, 1);
        field(site, ii, 1, 0) = su2Mat(1, 0);
        field(site, ii, 1, 1) = su2Mat(1, 1);
      }
      field(site, 3, 0, 0) = vev / sqrt(2) * higgsFn;
    }
    theory.applyBoundaryConditions(field);
  }
  monsta::addConstantMagneticField(field, theory, fluxQuanta, 2);
  theory.applyBoundaryConditions(field);

  // monsta::TheoryChecker checker(tol);
  // checker.checkGradients(theory, field);

  double E = theory.computeEnergy(field);
  COUT << E << endl;

  monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize, 100, true);
  solver.solve(theory, field);

  solver = monsta::GradDescentSolver(tol, extraSteps, initialStep, maxStepSize, extraSteps, true);
  solver.solve(theory, field);

  LATfield2::Field<complex<double> > referenceField(lattice, numMatrices, numRows, numCols, 0);

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
  // cout << norm << endl;

  for (site.first(); site.test(); site.next())
  {
    for (int ii = 0; ii < referenceField.components(); ii++)
    {
      referenceField(site, ii) = referenceField(site, ii)/ sqrt(norm);
    }
  }


  monsta::GradDescentSolverChigusa chigusaSolver(5e-4, 50000, initialStep, maxStepSize, 1.2, 0.5);
  chigusaSolver.solve(theory, field, referenceField);

  monsta::writeCoords(field, outputPath + "/coords");
  monsta::writeEnergyDensity(field, outputPath + "/energyData", theory);
  monsta::writeHiggsField(field, outputPath + "/higgsData", theory);
  monsta::writeMagneticField(field, outputPath + "/magneticFieldData", theory);
  monsta::writeRawField(field, outputPath + "/rawData");
  monsta::writeRawField(referenceField, outputPath + "/referenceRawData");

}