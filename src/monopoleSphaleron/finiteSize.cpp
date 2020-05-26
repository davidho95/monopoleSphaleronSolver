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
  int fluxQuanta = 0;
  int range = 1;

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
      case 'b':
        fluxQuanta = atof(argv[++i]);
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
  double vevInit = 1;
  double gaugeInit = 1;
  double selfInit = 0.5;

  LATfield2::Lattice lattice(dim, latSize, haloSize);
  LATfield2::Field<complex<double> > field(lattice, numMatrices, numRows, numCols, 0);

  LATfield2::Site site(lattice);

  double initialStep = 0.001;
  double maxStepSize = 0.05;
  double tol = 1e-4;
  int maxNumSteps = 200000;

  monsta::GeorgiGlashowSu2Theory theory(gaugeInit, vevInit, selfInit, {2, 0, 0}, true);
  monsta::setSingleMonopoleInitialConditions(field, theory);
  theory.applyBoundaryConditions(field);

  monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize);
  solver.solve(theory, field);

  LATfield2::Field<complex<double> > centredField(lattice, numMatrices, numRows, numCols, 0);
  std::vector<int> monopolePos = monsta::findMonopole(field, theory);
  int shiftNum = ((sz/2 - 1 - monopolePos[0]) + sz) % sz;
  monsta::circShift(field, centredField, theory, shiftNum, 0, true);

  monsta::GeorgiGlashowSu2Theory periodicTheory(gaugeInit, vevInit, selfInit, {0, 0, 0}, false);
  LATfield2::Field<complex<double> > pairField(lattice, numMatrices, numRows, numCols, 0);

  std::vector<double> vevs = { 1.0, 0.9, 0.8, 0.7 };
  std::vector<double> gVals = { 1.0, 1/sqrt(2), 0.5 };

  for (int sep = 1; sep < sz; sep++)
  {
    monsta::setPairInitialConds2(field, pairField, periodicTheory, sep);
    monsta::addConstantMagneticField(pairField, periodicTheory, -fluxQuanta);
    for (int ii = 0; ii < gVals.size(); ii++)
    {
      double selfCoupling = pow(gVals[ii], 2) / 2;
      monsta::setPairInitialConds2(field, pairField, periodicTheory, sep);
      monsta::addConstantMagneticField(pairField, periodicTheory, -fluxQuanta);
      for (int jj = 0; jj < vevs.size(); jj++)
      {
        monsta::GeorgiGlashowSu2Theory periodicTheory(gVals[ii], vevs[jj], selfCoupling, {0, 0, 0}, false);
        solver.solve(periodicTheory, pairField);
        double energy = periodicTheory.computeEnergy(pairField);

        std::string vevString = std::to_string(vevs[jj]);
        vevString.erase(vevString.find_last_not_of('0') + 1, std::string::npos);
        vevString.replace(1,1,"_");

        std::string gString = std::to_string(gVals[ii]);
        gString.erase(gString.find_last_not_of('0') + 1, std::string::npos);
        gString.replace(1,1,"_");

        if (parallel.rank() == 1)
        {
          ofstream fileStream;
          fileStream.open(outputPath + "/energiesG" + gString + "V" + vevString + ".txt", std::ios_base::app);
          fileStream << energy << endl;
          fileStream.close();
        }
      }
    }
  }
}