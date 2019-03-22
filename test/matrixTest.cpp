#include "LATfield2.hpp"
#include <complex>
#include "../src/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../src/Matrix.hpp"
#include "../src/TheoryChecker.hpp"
#include "../src/GradDescentSolverBBStep.hpp"
#include "../src/Su2Tools.hpp"
#include "../src/MonopoleFileTools.hpp"
#include "../src/MonopoleFieldTools.hpp"
#include <iostream>
#include <fstream>

int main() {
	parallel.initialize(2,2);

	srand(time(NULL) + parallel.rank());
	// srand(1);

	int dim = 3;
	int sz = 16;
	int latSize[dim] = {sz, sz, sz};
	int haloSize = 1;
	int numMatrices = 4;
	int numRows = 2;
	int numCols = 2;

	double gaugeCoupling = 1;
	double vev = 1;
	double selfCoupling = 1;
	// int monopoleSeparation = 4;

	LATfield2::Lattice lattice(dim, latSize, haloSize);
	LATfield2::Field<complex<double> > field(lattice, numMatrices, numRows, numCols, 0);

	LATfield2::Site site(lattice);

	monsta::GeorgiGlashowSu2Theory theory(gaugeCoupling, vev, selfCoupling, {-1, 1, 1}, true);

	// for (site.first(); site.test(); site.next())
	// {
	// 	int xCoord = site.coord(0);
	// 	int yCoord = site.coord(1);
	// 	int zCoord = site.coord(2);
	// 	double r = sqrt(pow(xCoord - sz/2, 2) + pow(yCoord - sz/2, 2) + pow(zCoord - sz/2, 2));
	// 	for (int ii = 0; ii < numMatrices - 1; ii++)
	// 	{
	// 		// monsta::Matrix su2Mat = monsta::vecToSu2({0, 0, 0});
	// 		monsta::Matrix su2Mat = abs(sz*sqrt(3) - r)/(sz*sqrt(3))*monsta::vecToSu2({(rand() % 20 - 10)*0.001, (rand() % 20 - 10)*0.001, (rand() % 20 - 10)*0.001});
	// 		// monsta::Matrix su2Mat = monsta::vecToSu2({0, zCoord, -yCoord});

	// 		field(site, ii, 0, 0) = su2Mat(0, 0);
	// 		field(site, ii, 0, 1) = su2Mat(0, 1);
	// 		field(site, ii, 1, 0) = su2Mat(1, 0);
	// 		field(site, ii, 1, 1) = su2Mat(1, 1);
	// 	}
	// 	field(site, 3, 0, 0) = vev;// + (rand() % 20 - 10)*0.01;
	// 	field(site, 3, 0, 1) = 0;
	// 	field(site, 3, 1, 0) = 0;
	// 	field(site, 3, 1, 1) = 0;
	// }
	for (site.first(); site.test(); site.next())
	{
		int xCoord = site.coord(0);
		int yCoord = site.coord(1);
		int zCoord = site.coord(2);
		double pi = 3.141592654;
		double r = sqrt(pow(xCoord - sz/2, 2) + pow(yCoord - sz/2, 2) + pow(zCoord - sz/2, 2));
		for (int ii = 0; ii < numMatrices - 1; ii++)
		{
			monsta::Matrix su2Mat = monsta::vecToSu2({0, 0, 0});
			if (ii == 2){
				su2Mat = monsta::vecToSu2({0, 0.1, 0.2});
			}
			// monsta::Matrix su2Mat = monsta::vecToSu2({0, 0, 0});

			// monsta::Matrix su2Mat = monsta::vecToSu2({0, zCoord, -yCoord});

			field(site, ii, 0, 0) = su2Mat(0, 0);
			field(site, ii, 0, 1) = su2Mat(0, 1);
			field(site, ii, 1, 0) = su2Mat(1, 0);
			field(site, ii, 1, 1) = su2Mat(1, 1);
		}
		field(site, 3, 0, 0) = vev;// + (rand() % 20 - 10)*0.01;
		field(site, 3, 0, 1) = 0;
		field(site, 3, 1, 0) = 0;
		field(site, 3, 1, 1) = 0;
	}
	theory.applyBoundaryConditions(field);


	// double bVal = 0.1;
	// int bDir = 0;
	// monsta::setConstantMagneticFieldUnitary(field, theory, bVal, bDir);
	// monsta::addMagneticFieldUnitary(field, theory, 0.001, 0);

	double tol = 1e-4;

	// monsta::TheoryChecker gradChecker(tol);
	// gradChecker.checkGradients(theory, field);

	monsta::GradDescentSolver solver(tol, 1000, 0.01, 0.1 / (vev*selfCoupling));
	solver.solve(theory, field);

	monsta::GeorgiGlashowSu2Theory periodicTheory(gaugeCoupling, vev, selfCoupling, {1, 1, 1}, true);
	LATfield2::Field<complex<double> > centredField(lattice, numMatrices, numRows, numCols, 0);
	std::vector<int> monopolePos = monsta::findMonopoleUnitary(field);
	COUT << monopolePos[0] << " " << monopolePos[1] << " " << monopolePos[2] << endl;
	int shiftNum = ((sz/2 - 1 - monopolePos[0]) + sz) % sz;
	monsta::circShift(field, centredField, theory, shiftNum, 0, true);


	// LATfield2::Field<complex<double> > pairField(lattice, numMatrices, numRows, numCols, 0);
	// LATfield2::Field<complex<double> > finalField(lattice, numMatrices, numRows, numCols, 0);
	// for (int sep = 4; sep <= 4; sep++)
	// {
	// 	monsta::setPairInitialConds(centredField, pairField, periodicTheory, sep);

	// 	solver.setParams(tol, 2000, 0.01, 0.01 / (vev*selfCoupling));
	// 	solver.solve(periodicTheory, pairField);

	// 	monsta::GeorgiGlashowSu2Theory dirichletTheory(gaugeCoupling, vev, selfCoupling, {0, 0, 0}, true);

	// 	monsta::circShift(pairField, finalField, periodicTheory, 0, 0, false); // abusing circshift as a copy
	// 	monsta::addMagneticFieldUnitary(finalField, dirichletTheory, 0.01, 0);

	// 	solver.setParams(tol, 2000, 0.01, 0.01 / (vev*selfCoupling));
	// 	solver.solve(dirichletTheory, finalField);

	// 	double E = dirichletTheory.computeEnergy(finalField);

 //    if (parallel.rank() == 1)
 //    {
 //      ofstream fileStream;
 //      fileStream.open("energies.txt", std::ios_base::app);
 //      fileStream << E << endl;
 //      fileStream.close();
 //    }
	// }

 	monsta::writeCoords(field, "coords");
 	monsta::writeHiggsFieldUnitary(field, "higgsData");
 	monsta::writeMagneticField(field, "magneticFieldData", theory);
 	monsta::writeEnergyDensity(field, "energyData", theory);
 	monsta::writeUnitaryGaugeField(field, "gaugeData");
  monsta::mergeFiles("coords");
  monsta::mergeFiles("higgsData");
  monsta::mergeFiles("magneticFieldData");
  monsta::mergeFiles("gaugeData");
  monsta::mergeFiles("energyData");

	return 0;
}