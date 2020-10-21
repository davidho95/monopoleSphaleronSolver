#include "LATfield2.hpp"
#include <complex>
#include "../../include/monopoleSphaleron/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../../include/Matrix.hpp"
#include "../../include/GradDescentSolverBBStep.hpp"
#include "../../include/GradDescentSolverBBStepNoCheckerboard.hpp"
#include "../../include/Su2Tools.hpp"
#include "../../include/MonopoleFileTools.hpp"
#include "../../include/MonopoleFieldTools.hpp"
#include "../../include/monopoleInstanton/GeorgiGlashowSu2Theory4d.hpp"
#include "../../include/monopoleInstanton/GeorgiGlashowRadialTheory.hpp"
#include "../../include/monopoleInstanton/InstantonFieldTools.hpp"
#include "../../include/monopoleInstanton/InstantonFileTools.hpp"
#include "../../include/monopoleInstanton/GradDescentSolverChigusa4d.hpp"
#include <iostream>
#include <fstream>


int main(int argc, char **argv)
{
  std::string outputPath;
  std::string inputPath;
  int sz = 16;
  int n = 2;
  int m = 2;
  double vev = 1;
  double gaugeCoupling = 1;
  double selfCoupling = 1;
  int sep = sz/2;
  double xAspect = 1;
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
      case 'd':
        sep = atoi(argv[++i]);
        break;
      case 'x':
        xAspect = atof(argv[++i]);
        break;
      case 'B':
        fluxQuanta = atoi(argv[++i]);
        break;
    }
  }

  parallel.initialize(n,m);

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

  LATfield2::Site site(lattice);

  double initialStep = 0.001;
  double maxStepSize = 0.002;
  double tol = 1e-5;
  int maxNumSteps = 10000;

  monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize, 100);
  monsta::GeorgiGlashowSu2Theory4d theory(gaugeCoupling, vev, selfCoupling);
  monsta::GeorgiGlashowRadialTheory radialTheory(gaugeCoupling, vev, selfCoupling);
  LATfield2::Field<complex<double> > field(lattice, numMatrices, numRows, numCols, 0);

  monsta::readRawField(field, inputPath + "/rawData");
  for (site.first(); site.test(); site.next())
  {
    for (int ii = 0; ii < 4; ii++)
    {
      theory.postProcess(field, site, ii);
    }
  }
  theory.applyBoundaryConditions(field);

  monsta::scaleVev(field, theory);
  
  int radialLatSize[dim] = {xSz/2, ySz, zSz};

  LATfield2::Lattice radialLattice(dim, radialLatSize, haloSize);
  LATfield2::Field<complex<double> > radialField(radialLattice, numMatrices, numRows, numCols, 0);
  
  LATfield2::Site radialSite(radialLattice);

  for (radialSite.first(); radialSite.test(); radialSite.next())
  {
    int radialCoord = radialSite.coord(0);
    int yCoord = radialSite.coord(1);
    int zCoord = radialSite.coord(2);
    int xCoord = radialCoord + xSz / 2;
    if (site.setCoord(xCoord, yCoord, zCoord))
    {
      for (int ii = 0; ii < 3; ii++)
      {
        monsta::Matrix fieldMat = theory.getSu2Link(field, site, ii);
        radialTheory.setSu2Link(radialField, radialSite, ii, fieldMat);
      }
      double higgsVal = theory.getHiggsMagnitude(field, site);
      radialTheory.setHiggsMagnitude(radialField, radialSite, higgsVal);
    }
  }
  radialTheory.applyBoundaryConditions(radialField);
  monsta::addConstantMagneticField(radialField, radialTheory, -fluxQuanta);

  solver.solve(radialTheory, radialField);

  monsta::writeRawField(radialField, outputPath + "/rawData");
  monsta::writeCoords(radialField, outputPath + "/coords");
  monsta::writeHiggsMagnitude(radialField, outputPath + "/higgsData");
  monsta::writeMagneticField(radialField, outputPath + "/magneticFieldData", radialTheory);
  monsta::writeEnergyDensity(radialField, outputPath + "/energyData", radialTheory);
  monsta::writeGradients(radialField, outputPath + "/gradData", radialTheory);
}