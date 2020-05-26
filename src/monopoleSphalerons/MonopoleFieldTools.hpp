#ifndef MONOPOLEFIELDTOOLS_HPP
#define MONOPOLEFIELDTOOLS_HPP

#include "LATfield2.hpp"
#include <complex>
#include "Su2Tools.hpp"
#include <cstdio>
#include <iostream>
#include <fstream>
// #include "GeorgiGlashowSu2TheoryUnitary.hpp"
#include "ElectroweakTheory.hpp"


namespace monsta
{

  std::vector<int> findMonopole(LATfield2::Field< std::complex<double> > &field, monsta::Theory &theory)
  {
    LATfield2::Site site(field.lattice());
    LATfield2::Site shiftedSite(field.lattice());


    int monopoleX = -1;
    int monopoleY = -1;
    int monopoleZ = -1;
    for (site.first(); site.test(); site.next())
    {
      double divB = 0;
      for (int ii = 0; ii < 3; ii++)
      {
        shiftedSite = site+ii;
        divB += theory.getMagneticField(field, shiftedSite, ii);
        divB -= theory.getMagneticField(field, site, ii);
      }

      double zeroTol = 0.1;
      if (abs(divB) > zeroTol)
      {
        monopoleX = site.coord(0);
        monopoleY = site.coord(1);
        monopoleZ = site.coord(2);
        break;
      }
    }
    parallel.max(monopoleX);
    parallel.max(monopoleY);
    parallel.max(monopoleZ);

    return {monopoleX, monopoleY, monopoleZ};
  }

  void circShift(LATfield2::Field< std::complex<double> > &inputField,
    LATfield2::Field< std::complex<double> > &outputField, monsta::Theory &theory, int numShifts, int dir, bool chargeConjugate)
  {
    LATfield2::Site siteFrom(inputField.lattice());
    LATfield2::Site siteTo(outputField.lattice());
    int size = inputField.lattice().size(dir);

    for (siteTo.first(); siteTo.test(); siteTo.next())
    {
      std::vector<int> inputCoords = {siteTo.coord(0), siteTo.coord(1), siteTo.coord(2)};
      inputCoords[dir] = (((inputCoords[dir] - numShifts) % size) + size) % size;

      siteFrom.setCoord(inputCoords[0], inputCoords[1], inputCoords[2]);
      int numMatrices = 4;
      for (int ii = 0; ii < 4; ii++)
      {
        monsta::Matrix matrixTo(inputField, siteFrom, ii);
        if (ii < 3 && chargeConjugate && siteTo.coord(dir) < numShifts)
        {
          matrixTo = pauli2*matrixTo*pauli2;
        }
        // matrixTo.print();
        outputField(siteTo, ii, 0, 0) = matrixTo(0, 0);
        outputField(siteTo, ii, 0, 1) = matrixTo(0, 1);
        outputField(siteTo, ii, 1, 0) = matrixTo(1, 0);
        outputField(siteTo, ii, 1, 1) = matrixTo(1, 1);
      }
    }
    // outputField.updateHalo();
    theory.applyBoundaryConditions(outputField);
  };

  void setPairInitialConds(LATfield2::Field< std::complex<double> > &singlePoleField,
    LATfield2::Field< std::complex<double> > &pairField, monsta::Theory &theory, int separation)
  {
    LATfield2::Site siteFrom(singlePoleField.lattice());
    LATfield2::Site siteTo(pairField.lattice());
  for (siteTo.first(); siteTo.test(); siteTo.next())
  {
    int xCoord = siteTo.coord(0);
    int yCoord = siteTo.coord(1);
    int zCoord = siteTo.coord(2);

    int xCoordFrom;
    bool chargeConjugate = false;
    int halfSeparation = separation / 2;
    int xSize = singlePoleField.lattice().size(0);
    if (xCoord < xSize/2)
    {
      if (separation % 2 == 0)
      {
        xCoordFrom = xCoord + halfSeparation;
      }
      else 
      {
        xCoordFrom = xCoord + halfSeparation + 1;
      }
    } else {
      xCoordFrom = xCoord - halfSeparation;
      chargeConjugate = true;
    }

    siteFrom.setCoord(xCoordFrom, yCoord, zCoord);

    int numMatrices = 4;
    {
      for (int ii = 0; ii < numMatrices; ii++)
      {
        monsta::Matrix matrixTo(singlePoleField, siteFrom, ii);
        if (chargeConjugate && ii != 3)
        {
          matrixTo = monsta::pauli2*matrixTo*monsta::pauli2;
        }
        pairField(siteTo, ii, 0, 0) = matrixTo(0, 0);
        pairField(siteTo, ii, 0, 1) = matrixTo(0, 1);
        pairField(siteTo, ii, 1, 0) = matrixTo(1, 0);
        pairField(siteTo, ii, 1, 1) = matrixTo(1, 1);
      }
    }
  }
  pairField.updateHalo();
  theory.applyBoundaryConditions(pairField);
  }

void linearSuperpose(LATfield2::Field< std::complex<double> > &field1,
  LATfield2::Field< std::complex<double> > &field2,
  LATfield2::Field< std::complex<double> > &outputField,
  monsta::Theory &theory)
{
  LATfield2::Site site(outputField.lattice());

  for (site.first(); site.test(); site.next())
  {
    for (int ii = 0; ii < 3; ii++)
    {
      std::vector<double> su2Vec1 = su2ToVec(Matrix(field1, site, ii));
      std::vector<double> su2Vec2 = su2ToVec(Matrix(field2, site, ii));
      std::vector<double> outputVec(3);

      for (int jj = 0; jj < 3; jj++)
      {
        outputVec[jj] = su2Vec1[jj] + su2Vec2[jj];
      }

      Matrix outputMat = vecToSu2(outputVec);

      outputField(site, ii, 0, 0) = outputMat(0, 0);
      outputField(site, ii, 0, 1) = outputMat(0, 1);
      outputField(site, ii, 1, 0) = outputMat(1, 0);
      outputField(site, ii, 1, 1) = outputMat(1, 1);
    }

    outputField(site, 3, 0, 0) = 0.5*(field1(site, 3, 0, 0) + field2(site, 3, 0, 0));
  }

  theory.applyBoundaryConditions(outputField);
}

void interpolateFields(LATfield2::Field< std::complex<double> > &field1,
  LATfield2::Field< std::complex<double> > &field2,
  LATfield2::Field< std::complex<double> > &outputField,
  monsta::Theory &theory, double interpolant)
{
  LATfield2::Site site(outputField.lattice());

  double maxDiff = 0;
  for (site.first(); site.test(); site.next())
  {
    for (int jj = 0; jj < field1.components(); jj++)
    {
      outputField(site, jj) = (1 - interpolant)*field1(site, jj) + interpolant*field2(site, jj);
    }
    for (int ii = 0; ii < 3; ii++)
    {
      theory.postProcess(outputField, site, ii);
    }
  }

  theory.applyBoundaryConditions(outputField);
}

void setPairInitialConds2(LATfield2::Field< std::complex<double> > &singlePoleField,
  LATfield2::Field< std::complex<double> > &pairField, monsta::Theory &theory, int separation)
{
  LATfield2::Site site(singlePoleField.lattice());

  int numRows = singlePoleField.rows();
  int numCols = singlePoleField.cols();
  int numMatrices = singlePoleField.components() / (numRows * numCols);

  LATfield2::Field< std::complex<double> > leftField(singlePoleField.lattice(), numMatrices, numRows, numCols, 0);
  LATfield2::Field< std::complex<double> > rightField(singlePoleField.lattice(), numMatrices, numRows, numCols, 0);

  std::vector<int> monopolePos = monsta::findMonopole(singlePoleField, theory);
  int xSize = singlePoleField.lattice().size(0);
  int desiredLeftPos = xSize/2 - (separation + 1)/2;
  int desiredRightPos = xSize/2 + separation/2;
  int shiftNumLeft = ((desiredLeftPos - monopolePos[0]) + xSize) % xSize;
  int shiftNumRight = ((desiredRightPos - monopolePos[0]) + xSize) % xSize;

  circShift(singlePoleField, leftField, theory, shiftNumLeft, 0, true);
  circShift(singlePoleField, rightField, theory, shiftNumRight, 0, true);


  if (monopolePos[0] < desiredLeftPos || monopolePos[0] > desiredRightPos)
  {
    for (site.first(); site.test(); site.next())
    {

      for (int ii = 0; ii < 3; ii++)
      {
        Matrix conjMat(rightField, site, ii);
        conjMat = pauli2*conjMat*pauli2;

        rightField(site, ii, 0, 0) = conjMat(0, 0);
        rightField(site, ii, 0, 1) = conjMat(0, 1);
        rightField(site, ii, 1, 0) = conjMat(1, 0);
        rightField(site, ii, 1, 1) = conjMat(1, 1);
      }
    }
  }

  linearSuperpose(leftField, rightField, pairField, theory);
}

  void contractPair(LATfield2::Field< std::complex<double> > &pairField,
    LATfield2::Field< std::complex<double> > &contractedField, monsta::Theory &theory)
  {
    LATfield2::Site siteFrom(pairField.lattice());
    LATfield2::Site siteTo(contractedField.lattice());
  for (siteTo.first(); siteTo.test(); siteTo.next())
  {
    int xCoord = siteTo.coord(0);
    int yCoord = siteTo.coord(1);
    int zCoord = siteTo.coord(2);

    int xCoordFrom;
    int xSize = pairField.lattice().size(0);
    if (xCoord < xSize/2)
    {
      xCoordFrom = xCoord;
    }
    else if (xCoord < xSize - 1)
    {
      xCoordFrom = xCoord + 1;
    }
    else
    {
      xCoord = xCoord;
    }

    siteFrom.setCoord(xCoordFrom, yCoord, zCoord);

    int numMatrices = 4;
    {
      for (int ii = 0; ii < numMatrices; ii++)
      {
        monsta::Matrix matrixTo(pairField, siteFrom, ii);
        if (ii != 3)
        {
          matrixTo = monsta::pauli2*matrixTo*monsta::pauli2;
        }
        contractedField(siteTo, ii, 0, 0) = matrixTo(0, 0);
        contractedField(siteTo, ii, 0, 1) = matrixTo(0, 1);
        contractedField(siteTo, ii, 1, 0) = matrixTo(1, 0);
        contractedField(siteTo, ii, 1, 1) = matrixTo(1, 1);
      }
    }
  }
  contractedField.updateHalo();
  theory.applyBoundaryConditions(contractedField);
  }

  void setVacuumField(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2Theory &theory)
  {
    double vev = theory.getVev();

    LATfield2::Site site(field.lattice());

    for (site.first(); site.test(); site.next())
    {
      for (int matIdx = 0; matIdx < 3; matIdx++)
      {
        field(site, matIdx, 0, 0) = 1;
        field(site, matIdx, 0, 1) = 0;
        field(site, matIdx, 1, 0) = 0;
        field(site, matIdx, 1, 1) = 1;
      }
      field(site, 3, 0, 0) = vev / sqrt(2);
      field(site, 3, 0, 1) = 0;
      field(site, 3, 1, 0) = 0;
      field(site, 3, 1, 1) = 0;
    }
    theory.applyBoundaryConditions(field);
  }

  void setRandomField(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2Theory &theory)
  {
    double vev = theory.getVev();
    double pi = 4*atan(1);

    LATfield2::Site site(field.lattice());

    for (site.first(); site.test(); site.next())
    {
      for (int matIdx = 0; matIdx < 3; matIdx++)
      {
        std::vector<double> su2Vec(3);
        for (int ii = 0; ii < 3; ii++)
        {
          su2Vec[ii] = (1e-5*double(rand() % 628318 - 314159));
        }
        Matrix su2Mat = vecToSu2(su2Vec);
        field(site, matIdx, 0, 0) = su2Mat(0, 0);
        field(site, matIdx, 0, 1) = su2Mat(0, 1);
        field(site, matIdx, 1, 0) = su2Mat(1, 0);
        field(site, matIdx, 1, 1) = su2Mat(1, 1);
      }
      field(site, 3, 0, 0) = 0.01*double(rand() % 100);
      field(site, 3, 0, 1) = 0;
      field(site, 3, 1, 0) = 0;
      field(site, 3, 1, 1) = 0;
    }
    theory.applyBoundaryConditions(field);
  }

  void setSingleMonopoleInitialConditions(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2Theory &theory)
  {
    // Assumes APP boundary conditions w/ flipped plaquettes along y-z edge

    LATfield2::Site site(field.lattice());
    double vev = theory.getVev();
    double gaugeCoupling = theory.getGaugeCoupling();
    int xSize = field.lattice().size(0);
    int ySize = field.lattice().size(1);
    int zSize = field.lattice().size(2);

    for (site.first(); site.test(); site.next())
    {
      double pi = 4*atan(1);

      int xCoord = site.coord(0);
      int yCoord = site.coord(1);
      int zCoord = site.coord(2);
      double r = sqrt(pow(xCoord - xSize/2, 2) + pow(yCoord - ySize/2, 2) + pow(zCoord - zSize/2, 2)) + 0.1;

      for (int ii = 0; ii < 3; ii++)
      {
        monsta::Matrix su2Mat = monsta::vecToSu2({0.01*double(rand() % 10), 0.01*double(rand() % 10), 0.01*double(rand() % 10)});

        std::vector<double> su2Vec(3);

        // Continuum gauge field (C-periodic BC's)
        if (ii  == 0) {
          su2Vec = {0, 0.5/gaugeCoupling*(zCoord - zSize/2) / pow(r,2), -0.5*gaugeCoupling*(yCoord - ySize/2) / pow(r,2)};
        }
        if (ii == 1) {
          su2Vec = {-0.5/gaugeCoupling*(zCoord - zSize/2) / pow(r,2), 0, 0.5*gaugeCoupling*(xCoord - xSize/2) / pow(r,2)};
        }
        if (ii == 2) {
          su2Vec = {0.5/gaugeCoupling*(yCoord - ySize/2) / pow(r,2), -0.5*gaugeCoupling*(xCoord - xSize/2) / pow(r,2), 0};
        }

        // Move flux line to centre
        if (yCoord + zCoord == zSize && yCoord < ySize / 2)
        {
          if (ii == 1)
          {
            su2Vec[2] += pi;
          }
          if (ii == 2)
          {
            su2Vec[1] += pi;
          }
        }
        if (yCoord == ySize / 2 && zCoord == zSize / 2)
        {
          if (ii == 2)
          {
            su2Vec[1] += pi;
          }
        }

        su2Mat = monsta::vecToSu2(su2Vec);
        field(site, ii, 0, 0) = su2Mat(0, 0);
        field(site, ii, 0, 1) = su2Mat(0, 1);
        field(site, ii, 1, 0) = su2Mat(1, 0);
        field(site, ii, 1, 1) = su2Mat(1, 1);
      }

      field(site, 3, 0, 0) = vev / sqrt(2);


      field(site, 3, 0, 1) = 0;
      field(site, 3, 1, 0) = 0;
      field(site, 3, 1, 1) = 0;

    }

    theory.applyBoundaryConditions(field);
  }

void addConstantMagneticField(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2Theory &theory,
  int fluxQuanta)
  {

    LATfield2::Site site(field.lattice());
    double vev = theory.getVev();
    double gaugeCoupling = theory.getGaugeCoupling();
    int xSize = field.lattice().size(0);
    int ySize = field.lattice().size(1);
    int zSize = field.lattice().size(2);

    double pi = 4*std::atan(1);
    double flux = fluxQuanta*4*pi;
    int fluxSign = fluxQuanta == 0 ? 0 : fluxQuanta / abs(fluxQuanta);

    for (site.first(); site.test(); site.next())
    {
      int xCoord = site.coord(0);
      int yCoord = site.coord(1);
      int zCoord = site.coord(2);

      for (int ii = 0; ii < 3; ii++)
      {
        std::vector<double> su2Vec = monsta::su2ToVec(Matrix(field, site, ii));
        if (ii == 2)
        {
          su2Vec[2] += 0.5*flux/pow(ySize,2)*(yCoord - 0.5);
          if (yCoord > 0)
          {
            su2Vec[2] -= 0.5*flux/ySize;
          }
        }
        if (ii == 1 && yCoord == 0)
        {
          su2Vec[2] -= 0.5*zCoord*flux/zSize;
        }
        monsta::Matrix su2Mat = monsta::vecToSu2(su2Vec);
        field(site, ii, 0, 0) = su2Mat(0, 0);
        field(site, ii, 0, 1) = su2Mat(0, 1);
        field(site, ii, 1, 0) = su2Mat(1, 0);
        field(site, ii, 1, 1) = su2Mat(1, 1);
      }
    }

    theory.applyBoundaryConditions(field);
  }

  void scaleVev(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2Theory &theory)
  {
    LATfield2::Site site(field.lattice());
    double oldVev = 0;
    for (site.first(); site.test(); site.next())
    {
      if (abs(field(site, 3, 0, 0))*sqrt(2) > oldVev)
      {
        oldVev = abs(field(site, 3, 0, 0))*sqrt(2);
      }
    }
    parallel.max(oldVev);
    
    double higgsRatio = theory.getVev() / oldVev;

    for (site.first(); site.test(); site.next())
    {
      field(site, 3, 0, 0) *= higgsRatio;
    }
    theory.applyBoundaryConditions(field);
  }

  double innerProduct(LATfield2::Field< std::complex<double> > &field1, LATfield2::Field< std::complex<double> > &field2)
  {
    LATfield2::Site site(field1.lattice());
    int numMatrices = 4;

    double output = 0;

    for (site.first(); site.test(); site.next())
    {
      for (int ii = 0; ii < field1.components(); ii++)
      {
        output += real(field1(site, ii)*conj(field2(site, ii)));
      }
      // for (int matIdx = 0; matIdx < 3; matIdx++)
      // {
      //   Matrix gaugeMat1(field1, site, matIdx);
      //   Matrix gaugeMat2(field2, site, matIdx);

      //   output += real(trace(gaugeMat1*gaugeMat2));
      // }
      // Matrix scalarMat1 = field1(site, 3, 0, 0)*pauli3;
      // Matrix scalarMat2 = field2(site, 3, 0, 0)*pauli3;
      // output += real(trace(scalarMat1*scalarMat2));
    }
    parallel.sum(output);
    return output;
  }
}

#endif