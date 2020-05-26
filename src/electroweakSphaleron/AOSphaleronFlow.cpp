#include "LATfield2.hpp"
#include <complex>
#include "../../include/electroweakSphaleron/ElectroweakTheory.hpp"
#include "../../include/monopoleSphaleron/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../../include/Matrix.hpp"
#include "../../include/GradDescentSolverBBStep.hpp"
#include "../../include/electroweakSphaleron/GradDescentSolverChigusaElectroweak.hpp"
#include "../../include/Su2Tools.hpp"
#include "../../include/MonopoleFileTools.hpp"
#include "../../include/electroweakSphaleron/ElectroweakFieldTools.hpp"
#include <iostream>
#include <fstream>
#include <bits/stdc++.h> 
#include <iostream> 
#include <sys/stat.h> 
#include <sys/types.h> 

int main(int argc, char **argv)
{
  srand(time(NULL) + parallel.rank());
  std::string outputPath;
  std::string vacuumInputPath;
  std::string sphaleronInputPath;
  int sz = 16;
  int n = 2;
  int m = 2;
  double vev = 1;
  double gaugeCoupling = 1;
  double selfCoupling = 1;
  double tanSqMixingAngle = 0.286;
  double fluxQuanta = 1;
  double zAspect = 1;
  double numInputFiles = 1;
  int extraSteps = 500;

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
        fluxQuanta = atoi(argv[++i]);
        break;
      case 'i':
        vacuumInputPath = argv[++i];
        break;
      case 'j':
        sphaleronInputPath = argv[++i];
        break;
      case 'z':
        zAspect = atof(argv[++i]);
        break;
      case 'q':
        tanSqMixingAngle = atof(argv[++i]);
        break;
      case 'N':
        numInputFiles = atoi(argv[++i]);
        break;
      case 'e':
        extraSteps = atoi(argv[++i]);
        break;
    }
  }

  parallel.initialize(n,m);

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
  LATfield2::Field<complex<double> > AOVacuumField(lattice, numMatrices, numRows, numCols, 0);
  LATfield2::Field<complex<double> > sphaleronField(lattice, numMatrices, numRows, numCols, 0);
  LATfield2::Field<complex<double> > AOSphaleronField(lattice, numMatrices, numRows, numCols, 0);
  LATfield2::Field<complex<double> > AO(lattice, numMatrices, numRows, numCols, 0);
  monsta::ElectroweakTheory theory(gaugeCoupling, tanSqMixingAngle, vev, selfCoupling);

  LATfield2::Site site(lattice);

  double initialStep = 0.01;
  double maxStepSize = 0.01;
  double tol = 1e-6;
  int minNumSteps = 20000;
  int maxNumSteps = 200000;

  COUT << zAspect << endl;

  if (numInputFiles == 1)
  {
    monsta::readFileWithCoords(AOVacuumField, vacuumInputPath + "/coords.txt", vacuumInputPath + "/rawData.txt");
  }
  else
  {
    for (int ii = 0; ii < numInputFiles; ii++)
    {
      monsta::readFileWithCoords(AOVacuumField, vacuumInputPath + "/coords" + to_string(ii) + ".txt", vacuumInputPath + "/rawData" + to_string(ii) + ".txt");
    }
  }
  theory.applyBoundaryConditions(AOVacuumField);

  for (site.first(); site.test(); site.next())
  {
    for (int ii = 0; ii < numMatrices; ii++)
    {
      theory.postProcess(AOVacuumField, site, ii);
    }
  }

  for (int ii = 0; ii < n; ii++)
  {
    LATfield2::Site copySite(lattice);
    int count = 0;
    for (site.first(); site.test(); site.next())
    {
      if (site.coord(2) < ii * (zSz / n)) { continue; }
      if (site.coord(2) > (ii + 1) * (zSz / n)) { continue; }
      count++;
      // cout << count << endl;
      if (site.coord(2) == 0) { continue; }
      copySite = site - 2;
      for (int ii = 0; ii < AOVacuumField.components(); ii++)
      {
        AOVacuumField(site, ii) = AOVacuumField(copySite, ii);
      }
    }
    theory.applyBoundaryConditions(AOVacuumField);
  }

  double E = theory.computeEnergy(AOVacuumField);
  COUT << E << endl;

  for (int ii = 0; ii < 256; ii++)
  {
    monsta::readFileWithCoords(AOVacuumField, vacuumInputPath + "/coords" + to_string(ii) + ".txt", vacuumInputPath + "/rawData" + to_string(ii) + ".txt");
  }
  theory.applyBoundaryConditions(sphaleronField);
  monsta::addConstantMagneticField(sphaleronField, theory, fluxQuanta, 2);

  E = theory.computeEnergy(sphaleronField);
  COUT << E << endl;

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
        AOSphaleronField(site, ii, 0, 0) = sphaleronField(site, ii, 0, 0);
        AOSphaleronField(site, ii, 0, 1) = sphaleronField(site, ii, 0, 1);
        AOSphaleronField(site, ii, 1, 0) = sphaleronField(site, ii, 1, 0);
        AOSphaleronField(site, ii, 1, 1) = sphaleronField(site, ii, 1, 1);
      }
      else
      {
        AOSphaleronField(site, ii, 0, 0) = AOVacuumField(site, ii, 0, 0);
        AOSphaleronField(site, ii, 0, 1) = AOVacuumField(site, ii, 0, 1);
        AOSphaleronField(site, ii, 1, 0) = AOVacuumField(site, ii, 1, 0);
        AOSphaleronField(site, ii, 1, 1) = AOVacuumField(site, ii, 1, 1);
      }
    }
  }
  theory.applyBoundaryConditions(AOSphaleronField);

  E = theory.computeEnergy(AOSphaleronField);
  COUT << E << endl;

  monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize, 100, true, true);
  solver.solve(theory, AOSphaleronField);

  solver = monsta::GradDescentSolver(tol, maxNumSteps, initialStep, maxStepSize, 100, false, true);
  solver.solve(theory, AOSphaleronField);
}
