#ifndef MONOPOLEFIELDTOOLS_HPP
#define MONOPOLEFIELDTOOLS_HPP

#include "LATfield2.hpp"
#include <complex>
#include "Su2Tools.hpp"
#include <cstdio>
#include <iostream>
#include <fstream>


namespace monsta
{
  std::vector<int> findMonopole(LATfield2::Field< std::complex<double> > &field)
  {
    LATfield2::Site site(field.lattice());
    std::vector<double> higgsVec;
    double minHiggsNorm = 1e10;
    int monopoleX = -1;
    int monopoleY = -1;
    int monopoleZ = -1;

    for (site.first(); site.test(); site.next())
    {
      higgsVec = su2LieAlgToVec(Matrix(field, site, 3));
      double higgsNorm = sqrt(pow(higgsVec[0], 2) + pow(higgsVec[1], 2) + pow(higgsVec[2], 2));
      if (higgsNorm < minHiggsNorm)
      {
        minHiggsNorm = higgsNorm;
      }
    }
    parallel.min(minHiggsNorm);
    for (site.first(); site.test(); site.next())
    {
      higgsVec = su2LieAlgToVec(Matrix(field, site, 3));
      double higgsNorm = sqrt(pow(higgsVec[0], 2) + pow(higgsVec[1], 2) + pow(higgsVec[2], 2));
      if (higgsNorm == minHiggsNorm)
      {
        monopoleX = site.coord(0);
        monopoleY = site.coord(1);
        monopoleZ = site.coord(2);
      }
    }
    parallel.max(monopoleX);
    parallel.max(monopoleY);
    parallel.max(monopoleZ);

    return {monopoleX, monopoleY, monopoleZ};
  }

  std::vector<int> findMonopoleUnitary(LATfield2::Field< std::complex<double> > &field)
  {
    LATfield2::Site site(field.lattice());
    std::vector<double> higgsVec;
    double minHiggsNorm = 1e10;
    int monopoleX = -1;
    int monopoleY = -1;
    int monopoleZ = -1;

    for (site.first(); site.test(); site.next())
    {
      double higgsNorm = real(field(site, 3, 0, 0));
      if (higgsNorm < minHiggsNorm)
      {
        minHiggsNorm = higgsNorm;
      }
    }
    parallel.min(minHiggsNorm);
    for (site.first(); site.test(); site.next())
    {
      double higgsNorm = real(field(site, 3, 0, 0));
      if (higgsNorm == minHiggsNorm)
      {
        monopoleX = site.coord(0);
        monopoleY = site.coord(1);
        monopoleZ = site.coord(2);
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
      // cout << inputCoords[0] << inputCoords[1] << inputCoords[2] << endl;
      inputCoords[dir] = (((inputCoords[dir] - numShifts) % size) + size) % size;
      // cout << inputCoords[0] << inputCoords[1] << inputCoords[2] << endl;

      siteFrom.setCoord(inputCoords[0], inputCoords[1], inputCoords[2]);
      int numMatrices = 4;
      for (int ii = 0; ii < 4; ii++)
      {
        monsta::Matrix matrixTo(inputField, siteFrom, ii);
        if (ii < 3 && chargeConjugate && siteTo.coord(dir) < numShifts)
        {
          matrixTo = monsta::pauli2*matrixTo*monsta::pauli2;
          // if (ii == 3)
          // {
          //   matrixTo = -1.0*matrixTo;
          // }
        }
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

  void transferBoundary(LATfield2::Field< std::complex<double> > &inputField,
    LATfield2::Field< std::complex<double> > &outputField)
  {
    LATfield2::Site site(inputField.lattice());
    int xMax = inputField.lattice().size(0) - 1;
    int yMax = inputField.lattice().size(1) - 1;
    int zMax = inputField.lattice().size(2) - 1;

    for (site.haloFirst(); site.haloTest(); site.haloNext())
    {
      int xCoord = site.coord(0);
      int yCoord = site.coord(1);
      int zCoord = site.coord(2);
      if (xCoord > 0 && yCoord > 0 && zCoord > 0 && xCoord < xMax && yCoord < yMax && zCoord < zMax)
      {
        continue;
      }
      for (int ii = 0; ii < 4; ii++)
      {
        outputField(site, ii, 0, 0) = inputField(site, ii, 0, 0);
        outputField(site, ii, 0, 1) = inputField(site, ii, 0, 1);
        outputField(site, ii, 1, 0) = inputField(site, ii, 1, 0);
        outputField(site, ii, 1, 1) = inputField(site, ii, 1, 1);
      }
    }
  }

  void setConstantMagneticFieldUnitary(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2Theory &theory, double fieldVal, int dir)
  {
    double vev = theory.getVev();

    LATfield2::Site site(field.lattice());
    for (site.first(); site.test(); site.next())
    {
      int coord = site.coord((dir + 1) % 3);
      for (int ii = 0; ii < 3; ii++)
      {
        monsta::Matrix su2Mat = monsta::vecToSu2({0, 0, 0});
        if (ii == (dir + 2) % 3)
        {
          su2Mat = monsta::vecToSu2({0,0,0.5*coord*fieldVal});
        }
        field(site, ii, 0, 0) = su2Mat(0, 0);
        field(site, ii, 0, 1) = su2Mat(0, 1);
        field(site, ii, 1, 0) = su2Mat(1, 0);
        field(site, ii, 1, 1) = su2Mat(1, 1);
      }
      field(site, 3, 0, 0) = vev;
      field(site, 3, 0, 1) = 0;
      field(site, 3, 1, 0) = 0;
      field(site, 3, 1, 1) = 0;
    }
    for (site.haloFirst(); site.haloTest(); site.haloNext())
    {
      for (int ii = 0; ii < 3; ii++)
      {
        int coord = site.coord((dir + 1) % 3);
        monsta::Matrix su2Mat = monsta::vecToSu2({0, 0, 0});
        if (ii == (dir + 2) % 3)
        {
          su2Mat = monsta::vecToSu2({0,0,0.5*coord*fieldVal});
        }
        field(site, ii, 0, 0) = su2Mat(0, 0);
        field(site, ii, 0, 1) = su2Mat(0, 1);
        field(site, ii, 1, 0) = su2Mat(1, 0);
        field(site, ii, 1, 1) = su2Mat(1, 1);
      }
      field(site, 3, 0, 0) = vev;
      field(site, 3, 0, 1) = 0;
      field(site, 3, 1, 0) = 0;
      field(site, 3, 1, 1) = 0;
    }
    theory.applyBoundaryConditions(field);
  }

void addMagneticFieldUnitary(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2Theory &theory, double fieldVal, int dir)
  {
    double vev = theory.getVev();

    LATfield2::Site site(field.lattice());
    for (site.first(); site.test(); site.next())
    {
      int coord = site.coord((dir + 1) % 3);
      for (int ii = 0; ii < 3; ii++)
      {
        monsta::Matrix bFieldMat = monsta::vecToSu2({0, 0, 0});
        if (ii == (dir + 2) % 3)
        {
          bFieldMat = monsta::vecToSu2({0,0,0.5*coord*fieldVal});
        }
        bFieldMat = monsta::Matrix(field, site, ii)*bFieldMat;
        field(site, ii, 0, 0) = bFieldMat(0, 0);
        field(site, ii, 0, 1) = bFieldMat(0, 1);
        field(site, ii, 1, 0) = bFieldMat(1, 0);
        field(site, ii, 1, 1) = bFieldMat(1, 1);
      }
    }
    for (site.haloFirst(); site.haloTest(); site.haloNext())
    {
      for (int ii = 0; ii < 3; ii++)
      {
        int coord = site.coord((dir + 1) % 3);
        monsta::Matrix bFieldMat = monsta::vecToSu2({0, 0, 0});
        if (ii == (dir + 2) % 3)
        {
          bFieldMat = monsta::vecToSu2({0,0,0.5*coord*fieldVal});
        }
        bFieldMat = monsta::Matrix(field, site, ii)*bFieldMat;
        field(site, ii, 0, 0) = bFieldMat(0, 0);
        field(site, ii, 0, 1) = bFieldMat(0, 1);
        field(site, ii, 1, 0) = bFieldMat(1, 0);
        field(site, ii, 1, 1) = bFieldMat(1, 1);
      }
    }
    theory.applyBoundaryConditions(field);
  }
}

#endif