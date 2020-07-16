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
  double tol = 5e-5;
  int minNumSteps = 20000;
  int maxNumSteps = 200000;

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

  monsta::readRawField(sphaleronField, sphaleronInputPath + "/rawData");
  theory.applyBoundaryConditions(sphaleronField);
  monsta::addConstantMagneticField(sphaleronField, theory, fluxQuanta, 2);

  for (site.first(); site.test(); site.next())
  {
    std::vector<double> coords(3);
    coords[0] = site.coord(0) - xSz/2 + 0.5;
    coords[1] = site.coord(1) - ySz/2 + 0.5;
    coords[2] = site.coord(2) - zSz/2 + 0.5;

    for (int ii = 0; ii < numMatrices; ii++)
    {
      if (abs(coords[0]) < 6 && abs(coords[1]) < 6 && site.coord(2) < 200)
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

  monsta::GradDescentSolver solver(tol, maxNumSteps, initialStep, maxStepSize, 100, true);
  solver.solve(theory, AOSphaleronField);

  LATfield2::Field<complex<double> > savedField(lattice, numMatrices, numRows, numCols, 0);
  for (site.first(); site.test(); site.next())
  {
    for (int matIdx = 0; matIdx < 4; matIdx++)
    {
      savedField(site, matIdx, 0, 0) = AOSphaleronField(site, matIdx, 0, 0);
      savedField(site, matIdx, 0, 1) = AOSphaleronField(site, matIdx, 0, 1);
      savedField(site, matIdx, 1, 0) = AOSphaleronField(site, matIdx, 1, 0);
      savedField(site, matIdx, 1, 1) = AOSphaleronField(site, matIdx, 1, 1);
    }
  }

  bool solved = false;

  while (!solved)
  {
    for (site.first(); site.test(); site.next())
    {
      for (int matIdx = 0; matIdx < 4; matIdx++)
      {
        AOSphaleronField(site, matIdx, 0, 0) = savedField(site, matIdx, 0, 0);
        AOSphaleronField(site, matIdx, 0, 1) = savedField(site, matIdx, 0, 1);
        AOSphaleronField(site, matIdx, 1, 0) = savedField(site, matIdx, 1, 0);
        AOSphaleronField(site, matIdx, 1, 1) = savedField(site, matIdx, 1, 1);
      }
    }
    solver = monsta::GradDescentSolver(tol, extraSteps, initialStep, maxStepSize, extraSteps);
    solver.solve(theory, AOSphaleronField);

    LATfield2::Field<complex<double> > referenceField(lattice, numMatrices, numRows, numCols, 0);

    for (site.first(); site.test(); site.next())
    {
      for (int ii = 0; ii < 4; ii++)
      {
        monsta::Matrix gradMat = theory.getLocalGradient(AOSphaleronField, site, ii);
        monsta::Matrix fieldMat(AOSphaleronField, site, ii);
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


    monsta::GradDescentSolverChigusa chigusaSolver(5e-4, 300000, initialStep, maxStepSize, 1.2, 0.5);
    solved = chigusaSolver.solve(theory, AOSphaleronField, referenceField);
    extraSteps = extraSteps + 500;
  }

  monsta::writeRawField(AOSphaleronField, outputPath + "/rawData");
  monsta::writeCoords(AOSphaleronField, outputPath + "/coords");
  monsta::writeHiggsMagnitude(AOSphaleronField, outputPath + "/higgsData");
  monsta::writeMagneticField(AOSphaleronField, outputPath + "/magneticFieldData", theory);
  monsta::writeEnergyDensity(AOSphaleronField, outputPath + "/energyData", theory);
}
