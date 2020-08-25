#include "LATfield2.hpp"
#include <complex>
#include "../../include/Matrix.hpp"
#include "../../include/GradDescentSolverBBStep.hpp"
#include "../../include/monopoleSphaleron/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../../include/monopoleInstanton/GradDescentSolverChigusa4d.hpp"
#include "../../include/Su2Tools.hpp"
#include "../../include/Theory.hpp"
#include "../../include/MonopoleFileTools.hpp"
// #include "../../include/MonopoleFieldTools.hpp"
#include "../../include/monopoleInstanton/InstantonFileTools.hpp"
#include "../../include/monopoleInstanton/InstantonFieldTools.hpp"
#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{
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
    }
  }
  parallel.initialize(n,m);

  srand(1);

  int dim = 3;
  int xSz = round(xAspect*sz);
  int ySz = sz;
  int zSz = sz;
  int latSize[dim] = {xSz, ySz, zSz};
  int haloSize = 1;
  int numMatrices = 4;
  int numRows = 2;
  int numCols = 2;


  LATfield2::Lattice lattice(dim, latSize, haloSize);
  LATfield2::Field<complex<double> > field(lattice, numMatrices, numRows, numCols, 0);
  monsta::GeorgiGlashowSu2Theory4d theory(gaugeCoupling, vev, selfCoupling);

  double initialStep = 0.001;
  double maxStepSize = 0.01;
  double tol = 1e-3;
  int maxNumSteps = 20000;

  LATfield2::Site site(lattice);

  monsta::readRawField(field, inputPath + "/rawData");
  for(site.first(); site.test(); site.next())
  {
    for (int ii = 0; ii < 4; ii++)
    {
      theory.postProcess(field, site, ii);
    }
  }
  monsta::scaleVev(field, theory);

  monsta::addConstantMagneticField(field, theory, -fluxQuanta);
  theory.applyBoundaryConditions(field);

  double E = theory.computeEnergy(field);
  COUT << E << endl;

  monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize, 100, true);
  solver.solve(theory, field);

  monsta::writeCoords(field, outputPath + "/coords");
  monsta::writeEnergyDensity(field, outputPath + "/energyData", theory);
  monsta::writeHiggsField(field, outputPath + "/higgsData", theory);
  monsta::writeMagneticField(field, outputPath + "/magneticFieldData", theory);
  monsta::writeRawField(field, outputPath + "/rawData");
  monsta::writeGradients(field, outputPath + "/gradData", theory);

  solver = monsta::GradDescentSolver(tol, extraSteps, initialStep, maxStepSize, extraSteps);
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
      referenceField(site, ii) = referenceField(site, ii)/sqrt(norm);
    }
  }


  monsta::GradDescentSolverChigusa chigusaSolver(tol, maxNumSteps, initialStep, maxStepSize, beta, 0.5);
  chigusaSolver.solve(theory, field, referenceField);

  monsta::writeCoords(field, outputPath + "/coords");
  monsta::writeEnergyDensity(field, outputPath + "/energyData", theory);
  monsta::writeHiggsField(field, outputPath + "/higgsData", theory);
  monsta::writeMagneticField(field, outputPath + "/magneticFieldData", theory);
  monsta::writeRawField(field, outputPath + "/rawData");
  monsta::writeRawField(referenceField, outputPath + "/referenceRawData");

}
