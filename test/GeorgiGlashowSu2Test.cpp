#include "LATfield2.hpp"
#include "TheoryChecker.hpp"
#include "../src/GradDescentSolver.hpp"
#include "../src/GeorgiGlashowSu2Theory.hpp"
#include "../src/FreeMonopoleTheory.hpp"
#include "../src/System.hpp"
#include "../src/Matrix.hpp"
#include <iostream>
#include <fstream>

int main() {
  parallel.initialize(1,1);

  srand(time(NULL) + parallel.rank());
  // srand(1);

  int dim = 3;
  int sz = 16;
  int latSize[dim] = {sz, sz, sz};
  int haloSize = 1;

  double tol = 1e-2;

  LATfield2::Lattice lattice(dim, latSize, haloSize);

  double gaugeCoupling = 1;
  double vev = 1;
  double selfCoupling = 8;

  int monopolePos1[dim] = {0, sz/2, sz/2};
  int monopolePos2[dim] = {sz/2, sz/2, sz/2};

  monsta::GeorgiGlashowSu2Theory theory(gaugeCoupling, vev, selfCoupling, true);//, monopolePos1, monopolePos2);
  monsta::System testSystem(theory, lattice);
  testSystem.setConstantField(0);
  int numCpts = testSystem.getNumSpatialCpts();

  int maxIter = 100000;

  monsta::GradDescentSolver solver(tol, maxIter);
  solver.setVerbosity(true);

  // cout << testSystem.getEnergy();
  // testSystem.updateGradientReturnMax();

  monsta::TheoryChecker gradChecker(tol, true);

  LATfield2::Site site(lattice);

  // for (site.first(); site.test(); site.next())
  // {
  //   for (int ii = 0; ii < numCpts; ii++)
  //   {
  //     testSystem.setFieldElement(site, ii, (rand() % (2*628) - 628)*0.01);
  //   }
  // }
  // testSystem.updateFieldHalo();
  // testSystem.updateGradient();

  // for (int ii = 0; ii < 11; ii++) {
  //   cout << testSystem.getFieldElement(site, ii) << " ";
  // }
  // cout << endl;

  // cout << testSystem.getEnergy() << endl;

  // gradChecker.checkGradients(testSystem);

  
  int scalarFieldStart = theory.getScalarFieldStart();

  for (site.first(); site.test(); site.next())
  {
    int xCoord = site.coord(0);
    int yCoord = site.coord(1);
    int zCoord = site.coord(2);
    double norm = sqrt(pow(vev,2)/2.);
    for (int ii = 0; ii < scalarFieldStart; ii++)
    {
      // if (ii % 3 != 1) {
        testSystem.setFieldElement(site, ii, 3.1415926535897932384626 / sqrt(3));
      // }
    }

    double theta = 0.01*(rand() % 314);
    double phi = 0;//3.1415926535897932384626 * (rand() % 2);

    testSystem.setFieldElement(site, 9, norm*sin(theta)*cos(phi));
    testSystem.setFieldElement(site, 10, norm*sin(theta)*sin(phi));
    testSystem.setFieldElement(site, 11, norm*cos(theta));

    for (int ii = scalarFieldStart; ii < numCpts; ii++)
    {
      // testSystem.setFieldElement(site, ii, 0.1);
      // testSystem.setFieldElement(site, ii, (rand() % 200 - 100)*0.01);
      if (xCoord == monopolePos1[0] && yCoord == monopolePos1[1] && zCoord == monopolePos1[2]) {
        // testSystem.setFieldElement(site, ii, 0.0001);
      } else if (xCoord == monopolePos2[0] && yCoord == monopolePos2[1] && zCoord == monopolePos2[2]) {
        // testSystem.setFieldElement(site, ii, 0.0001);
      }
    }
  }
  testSystem.updateFieldHalo();
  testSystem.updateGradient();

  solver.solve(testSystem);

  monsta::GeorgiGlashowSu2Theory theory2(gaugeCoupling, vev, selfCoupling, monopolePos1, monopolePos2);
  monsta::System system2(theory2, lattice);

  for (site.first(); site.test(); site.next())
  { for (int ii = 0; ii < numCpts; ii++)
    {
      int xCoord = site.coord(0);
      int yCoord = site.coord(1);
      int zCoord = site.coord(2);
      system2.setFieldElement(site, ii, testSystem.getFieldElement(site, ii));
      if (xCoord == monopolePos1[0] && yCoord == monopolePos1[1] && zCoord == monopolePos1[2]) {
        system2.setFieldElement(site, ii, 0.00001);
      } else if (xCoord == monopolePos2[0] && yCoord == monopolePos2[1] && zCoord == monopolePos2[2]) {
        system2.setFieldElement(site, ii, 0.00001);
      }
    }
  }
  system2.updateFieldHalo();
  system2.updateGradient();

  solver.solve(system2);

  LATfield2::Field<double> outputField(lattice, numCpts);
  outputField.initialize(lattice, numCpts);
  outputField.alloc();

  for (site.first(); site.test(); site.next())
  {
    for (int ii = 0; ii < numCpts; ii++)
    {
      outputField(site, ii) = system2.getFieldElement(site, ii);
    }
  }
  outputField.updateHalo();



  ofstream fileStream1;
  ofstream fileStream2;
  ofstream fileStream3;
  fileStream1.open("gaugeData.txt");
  fileStream2.open("higgsData.txt");
  fileStream3.open("energyData.txt");

  for (site.first(); site.test(); site.next())
  {
    int xCoord = site.coord(0);
    int yCoord = site.coord(1);
    int zCoord = site.coord(2);
    // fileStream3 << theory.getLocalEnergyDensity(outputField, site, outputField) << endl;
    for (int ii = 0; ii < 3; ii++)
    {
      fileStream1 << theory.getMagneticField(outputField, site, ii) << " ";
    }
    fileStream1 << endl;
    fileStream2 << site.coord(0) << " " << site.coord(1) << " " << site.coord(2) << " ";
    for (int ii = scalarFieldStart; ii < numCpts; ii++)
    {
      fileStream2 << system2.getFieldElement(site, ii) << " ";
    }
    fileStream2 << endl;

    double eGauge = 0;

    int numDimensions = 3;

    int twist = 1;

    // Wilson action
    eGauge += 2/gaugeCoupling*(2 - real(monsta::trace(theory.getPlaquette(outputField, site, 1, 2))));
    eGauge += 2/gaugeCoupling*(2 - real(monsta::trace(theory.getPlaquette(outputField, site, 2, 0))));
    eGauge += 2/gaugeCoupling*(2 - real(monsta::trace(theory.getPlaquette(outputField, site, 0, 1))));
    fileStream3 << eGauge << " ";

    double eKin = 0;
    // Covariant Derivative
    for (int ii = 0; ii < numDimensions; ii++)
    {
      eKin += real(monsta::trace(theory.getKineticTerm(outputField, site, ii)));
    }

    fileStream3 << eKin << " ";

    // Higgs Potential
    double eHiggs = 0;
    double scalarVecNormSq = pow(outputField(site, scalarFieldStart), 2)
                   + pow(outputField(site, scalarFieldStart + 1), 2)
                   + pow(outputField(site, scalarFieldStart + 2), 2);
    eHiggs += selfCoupling*pow((2*scalarVecNormSq - pow(vev,2)),2);
    fileStream3 << eHiggs << " ";
    fileStream3 << endl;
      
    }
  fileStream1.close();
  fileStream2.close();
  fileStream3.close();

  return 0;
}

