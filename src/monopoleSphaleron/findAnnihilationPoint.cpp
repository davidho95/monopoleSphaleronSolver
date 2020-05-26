#include "LATfield2.hpp"
#include <complex>
#include "../include/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../include/Matrix.hpp"
#include "../include/GradDescentSolverBBStep.hpp"
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
  double vevInit = 1;
  double gaugeCouplingInit = 1;
  double selfCouplingInit = 1;
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
        vevInit = atof(argv[++i]);
        break;
      case 'g':
        gaugeCouplingInit = atof(argv[++i]);
        break;
      case 'l':
        selfCouplingInit = atof(argv[++i]);
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

  double initialStep = 0.01;
  double maxStepSize = 0.05;
  double tol = 1e-6;
  int maxNumSteps = 20000;

  monsta::GeorgiGlashowSu2Theory theory(gaugeCouplingInit, vevInit, selfCouplingInit, {0, 0, 0}, false);
  monsta::readRawField(field, inputPath + "/rawData");
  theory.applyBoundaryConditions(field);

  std::vector<double> vevs;
  double currentVev = vevInit;
  double increment = 0.01;

  monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize);

  while (currentVev > 0)
  { 
    vevs.push_back(currentVev);
    currentVev -= increment;
  }

  double annihilationLimit = 20;

  double vev, gaugeCoupling, selfCoupling;
  for (int ii = 0; ii < vevs.size(); ii++)
  {
    vev = vevs[ii];
    gaugeCoupling = vev;
    selfCoupling = pow(gaugeCoupling, 2) / 2;
    monsta::GeorgiGlashowSu2Theory currentTheory(gaugeCoupling, vev, selfCoupling, {0, 0, 0}, false);

    monsta::scaleVev(field, currentTheory);
    solver.solve(currentTheory, field);

    double E = currentTheory.computeEnergy(field);
    if (E < annihilationLimit) { break; }

    monsta::writeRawField(field, outputPath + "/rawData");
    monsta::writeCoords(field, outputPath + "/coords");
    monsta::writeHiggsFieldUnitary(field, outputPath + "/higgsData");
    monsta::writeMagneticField(field, outputPath + "/magneticFieldData", theory);
    monsta::writeEnergyDensity(field, outputPath + "/energyData", theory);
    monsta::writeUnitaryGaugeField(field, outputPath + "/gaugeData");
  }

  COUT << vev << endl;
}