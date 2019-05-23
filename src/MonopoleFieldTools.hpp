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
      double higgsNorm = abs(real(field(site, 3, 0, 0)));
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
      // cout << isLocal << endl;
      int numMatrices = 4;
      for (int ii = 0; ii < 4; ii++)
      {
        monsta::Matrix matrixTo(inputField, siteFrom, ii);
        if (ii < 3 && chargeConjugate && siteTo.coord(dir) < numShifts)
        {
          monsta::Matrix conjMat(2);
          if (dir == 0) {
            conjMat = pauli1;
          }
          if (dir == 1) {
            conjMat = pauli2;
          }
          if (dir == 2) {
            conjMat = pauli3;
          }
          matrixTo = conjMat*matrixTo*conjMat;
          if (ii == 3 && dir == 2)
          {
            matrixTo(0,0) = -1.0*matrixTo(0,0);
          }
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
          matrixTo = monsta::pauli1*matrixTo*monsta::pauli1;
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
    int size1 = field.lattice().size((dir + 1) % 3);
    int size2 = field.lattice().size((dir + 2) % 3);
    for (site.first(); site.test(); site.next())
    {
      int coord1 = site.coord((dir + 1) % 3);
      int coord2 = site.coord((dir + 2) % 3);
      for (int ii = 0; ii < 3; ii++)
      {
        std::vector<double> su2Vec = {0, 0, 0};
        if (ii == (dir + 2) % 3)
        {
          // su2Mat = monsta::vecToSu2({0,0,0.5*(coord1 - size1 / 2)*fieldVal});
        }
        if (ii == (dir + 1) % 3)
        {
          su2Vec[2] -= 0.5*(coord2)*fieldVal;
          // su2Mat = monsta::vecToSu2({0,0,-0.5*(coord2 - size2 / 2)*fieldVal});
        }
        // if (ii == 2 && site.coord(2) == field.lattice().size(2) - 1)
        // {
        //   su2Vec[1] += 3.14159/2;
        // }
        monsta::Matrix su2Mat = monsta::vecToSu2(su2Vec);
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
    // for (site.haloFirst(); site.haloTest(); site.haloNext())
    // {
    //   for (int ii = 0; ii < 3; ii++)
    //   {
    //     int coord1 = site.coord((dir + 1) % 3);
    //     int coord2 = site.coord((dir + 2) % 3);
    //     monsta::Matrix su2Mat = monsta::vecToSu2({0, 0, 0});
    //     if (ii == (dir + 2) % 3)
    //     {
    //       // su2Mat = monsta::vecToSu2({0,0,0.5*(coord1 - size1 / 2)*fieldVal});
    //     }
    //     if (ii == (dir + 1) % 3)
    //     {
    //       su2Mat = monsta::vecToSu2({0,0,-0.5*(coord2 - size2 / 2)*fieldVal});
    //     }
    //     field(site, ii, 0, 0) = su2Mat(0, 0);
    //     field(site, ii, 0, 1) = su2Mat(0, 1);
    //     field(site, ii, 1, 0) = su2Mat(1, 0);
    //     field(site, ii, 1, 1) = su2Mat(1, 1);
    //   }
    //   field(site, 3, 0, 0) = vev / sqrt(2);
    //   field(site, 3, 0, 1) = 0;
    //   field(site, 3, 1, 0) = 0;
    //   field(site, 3, 1, 1) = 0;
    // }
    theory.applyBoundaryConditions(field);
  }

  void scaleUpField(LATfield2::Field< std::complex<double> > &smallField, LATfield2::Field< std::complex<double> > &bigField)
  {
    LATfield2::Site smallSite(smallField.lattice());
    LATfield2::Site bigSite(bigField.lattice());

    int scaleFactor = 2;

    int smallSize = smallField.lattice().size(0);
    for (bigSite.first(); bigSite.test(); bigSite.next())
    {
      int bigXCoord = bigSite.coord(0);
      int bigYCoord = bigSite.coord(1);
      int bigZCoord = bigSite.coord(2);

      if (bigYCoord == smallSize*scaleFactor/2 + 1 && bigZCoord == smallSize*scaleFactor/2 + 1)
      {
        smallSite.setCoord(bigXCoord / scaleFactor + 1, bigYCoord / scaleFactor + 1, bigZCoord / scaleFactor + 1);
      }
      else
      {
        smallSite.setCoord(bigXCoord / scaleFactor, bigYCoord / scaleFactor, bigZCoord / scaleFactor);
      }

      for (int ii = 0; ii < bigField.components(); ii++)
      {
        bigField(bigSite, ii) = smallField(smallSite, ii);
      }
    }
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

void addConstantMagneticField(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2Theory &theory)
  {

    LATfield2::Site site(field.lattice());
    double vev = theory.getVev();
    double gaugeCoupling = theory.getGaugeCoupling();
    int xSize = field.lattice().size(0);
    int ySize = field.lattice().size(1);
    int zSize = field.lattice().size(2);

    double pi = 4*std::atan(1);
    double flux = 4*pi;

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
            su2Vec[2] -= 2*pi/ySize;
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

  void transplantMonopole(LATfield2::Field< std::complex<double> > &smallField, LATfield2::Field< std::complex<double> > &bigField)
  {
    LATfield2::Site smallSite(smallField.lattice());
    LATfield2::Site bigSite(bigField.lattice());

    int scaleFactor = 2;

    int smallSize = smallField.lattice().size(0);

    // for (smallSite.first(); smallSite.test(); smallSite.next())
    // {
    //   int smallXCoord = smallSite.coord(0);
    //   int smallYCoord = smallSite.coord(1);
    //   int smallZCoord = smallSite.coord(2);

    //   int bigXCoord;
    //   if (smallXCoord < smallSize / 2)
    //   {
    //     bigXCoord = smallXCoord;
    //   }
    //   else
    //   {
    //     bigXCoord = smallXCoord + smallSize;
    //   }
    //   int bigYCoord = smallYCoord + smallSize / 2;
    //   int bigZCoord = smallZCoord + smallSize / 2;

    //   bigSite.setCoord(bigXCoord, bigYCoord, bigZCoord);
    //   for (int ii = 0; ii < bigField.components(); ii++)
    //   {
    //     bigField(bigSite, ii) = smallField(smallSite, ii);
    //   }
    // }

    for (bigSite.first(); bigSite.test(); bigSite.next())
    {
      int bigXCoord = bigSite.coord(0);
      int bigYCoord = bigSite.coord(1);
      int bigZCoord = bigSite.coord(2);

      int smallXCoord, smallYCoord, smallZCoord;
      if (bigXCoord < smallSize / 2)
      {
        smallXCoord = bigXCoord; 
      }
      else if (bigXCoord > 3 * smallSize / 2)
      {
        smallXCoord = bigXCoord - (2*smallSize);
      }
      else
      {
        smallXCoord = 0;
      }

      if (bigYCoord - smallSize / 2 < 0 || bigYCoord - smallSize / 2 >= smallSize)
      {
        smallYCoord = bigYCoord - smallSize / 2;
      }
      else
      {
        smallYCoord = 0;
      }

      if (bigZCoord - smallSize / 2 < 0 || bigZCoord - smallSize / 2 >= smallSize)
      {
        smallZCoord = bigZCoord - smallSize / 2;
      }
      else
      {
        smallZCoord = 0;
      }

      smallSite.setCoord(smallXCoord, smallYCoord, smallZCoord);

      for (int ii = 0; ii < bigField.components(); ii++)
      {
        bigField(bigSite, ii) = smallField(smallSite, ii);
      }
    }
  }
}

#endif