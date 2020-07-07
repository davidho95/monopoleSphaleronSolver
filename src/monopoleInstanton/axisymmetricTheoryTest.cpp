#include "LATfield2.hpp"
#include <complex>
#include "../../include/monopoleInstanton/GeorgiGlashowSu2TheoryAxisymmetric2.hpp"
#include "../../include/monopoleInstanton/AxisymmetricFieldTools.hpp"
#include "../../include/TheoryChecker.hpp"
#include "../../include/GradDescentSolverBBStep.hpp"
#include "../../include/MonopoleFileTools.hpp"
#include "../../include/monopoleInstanton/InstantonFileTools.hpp"

int main(int argc, char **argv)
{
  int sz = 16;
  int n = 2;
  int m = 2;
  double gaugeCoupling = 1;
  double selfCoupling = 1;
  double vev = 1;
  double numThetaPoints = 100;
  std::string outputPath;
  std::string inputPath;

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
        sz = atoi(argv[++i]);
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
      case 'N':
        numThetaPoints = atof(argv[++i]);
        break;
      case 'n':
        n = atoi(argv[++i]);
        break;
      case 'm':
        m = atoi(argv[++i]);
        break;
    }
  }

  parallel.initialize(n,m);

  srand(time(NULL) + parallel.rank());

  int dim = 2;
  int zSz = sz;
  int rSz = sz;
  int thetaSz = 1;
  int latSize[dim] = {zSz, rSz};
  int haloSize = 1;
  int numMatrices = 4;
  int numRows = 2;
  int numCols = 2;

  LATfield2::Lattice lattice(dim, latSize, haloSize);
  LATfield2::Field<complex<double> > field(lattice, numMatrices, numRows, numCols, 0);
  monsta::GeorgiGlashowSu2TheoryAxisymmetric theory(gaugeCoupling, vev, selfCoupling, numThetaPoints);

  LATfield2::Site site(lattice);

  if (inputPath != "")
  {
    monsta::readRawField(field, inputPath + "/rawData");

    for (site.first(); site.test(); site.next())
    {
      for (int ii = 0; ii < numMatrices; ii++)
      {
        theory.postProcess(field, site, ii);
      }
    }
    theory.applyBoundaryConditions(field);
    monsta::scaleVev(field, theory);
  }
  else
  {
    monsta::setVacuumField(field, theory);
    // monsta::perturbField(field, theory);
    // monsta::setRandomField(field, theory);

    // for (site.first(); site.test(); site.next())
    // {
    //   std::vector<double> su2Vec0(3);
    //   std::vector<double> su2Vec2(3);
    //   double r = theory.getRFromSite(site);
    //   double zCoord = site.coord(0) - sz/2 + 0.5;
    //   theory.setHiggsMagnitude(field, site, vev/sqrt(2));
    //   // cout << monsta::pi/2*(1 - monsta::sign(r)) << endl;
    //   std::vector<double> su2Vec = {0, 0, monsta::pi/2*(1 + monsta::sign(r))};
    //   theory.setSu2Link(field, site, 2, monsta::vecToSu2(su2Vec));
    //   theory.setSu2Link(field, site, 1, monsta::identity);
    //   theory.setSu2Link(field, site, 0, monsta::identity);
    // } 
    // theory.applyBoundaryConditions(field);
    // monsta::perturbField(field, theory);

    double rMax = sz - sz/2 + 0.5;
    double B = 0.5*monsta::pi / pow(rMax, 2);
    for (site.first(); site.test(); site.next())
    {
      std::vector<double> su2Vec0(3);
      std::vector<double> su2Vec2(3);
      double r = theory.getRFromSite(site);
      double zCoord = site.coord(0) - sz/2 + 0.5;
      theory.setHiggsMagnitude(field, site, vev / sqrt(2));
      // cout << monsta::pi/2*(1 - monsta::sign(r)) << endl;
      std::vector<double> su2Vec = {0, 0, (B*monsta::sign(r)*pow(r,2) + monsta::pi/2)};
      theory.setSu2Link(field, site, 2, monsta::vecToSu2(su2Vec));
      theory.setSu2Link(field, site, 1, monsta::identity);
      theory.setSu2Link(field, site, 0, monsta::identity);
    } 
    theory.applyBoundaryConditions(field);
  }
  monsta::perturbField(field, theory, 1e-4);

  double tol = 1e-3;
  int maxIterations = 10000;
  double initialStepSize = 0.001;
  double maxStepSize = 0.002;

  double E = theory.computeEnergy(field);
  COUT << E << endl;

  // for (site.first(); site.test(); site.next())
  // {
  //   cout << theory.getPlaquetteAveragedMetricFactors(site, 1, 2) << endl;
  // }

  monsta::GradDescentSolver solver(tol, maxIterations, initialStepSize, maxStepSize, 100);
  solver.solve(theory, field);

  // monsta::TheoryChecker checker(tol);
  // checker.checkGradients(theory, field);

  monsta::writeCoords(field, outputPath + "/coords");
  monsta::writeEnergyDensity(field, outputPath + "/energyData", theory);
  monsta::writeHiggsField(field, outputPath + "/higgsData", theory);
  monsta::writeMagneticField(field, outputPath + "/magneticFieldData", theory);
  monsta::writeRawField(field, outputPath + "/rawData");
}