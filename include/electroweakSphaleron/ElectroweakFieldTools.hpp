#ifndef ELECTROWEAKFILETOOLS_HPP
#define ELECTROWEAKFILETOOLS_HPP

#include "LATfield2.hpp"
#include <complex>
#include "../Su2Tools.hpp"
#include <cstdio>
#include <iostream>
#include <fstream>
#include "ElectroweakTheory.hpp"


namespace monsta
{

  void setVacuumField(LATfield2::Field< std::complex<double> > &field, monsta::ElectroweakTheory &theory)
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
      field(site, 3, 0, 1) = 1;
      field(site, 3, 1, 0) = 1;
      field(site, 3, 1, 1) = 1;
    }
    theory.applyBoundaryConditions(field);
  }

  void setRandomField(LATfield2::Field< std::complex<double> > &field, monsta::ElectroweakTheory &theory)
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
      for (int ii = 0; ii < 3; ii++)
      {
        double theta = 1e-5*double(rand() % 628318 - 314159);
        field(site, 3, (ii + 1) % 2, (ii + 1) / 2) = cos(theta) + 1i*sin(theta);
      }
      field(site, 3, 1, 1) = 0.01*double(rand() % 100);
    }
    theory.applyBoundaryConditions(field);
  }

  void addConstantMagneticField(LATfield2::Field< std::complex<double> > &field, monsta::ElectroweakTheory &theory,
  int fluxQuanta, int fieldDir = 0)
  {

    LATfield2::Site site(field.lattice());
    double vev = theory.getVev();
    double gaugeCoupling = theory.getGaugeCoupling();
    std::vector<int> latSize = {field.lattice().size(0), field.lattice().size(1), field.lattice().size(2)};
    int xSize = field.lattice().size(0);
    int ySize = field.lattice().size(1);
    int zSize = field.lattice().size(2);

    double pi = 4*std::atan(1);
    double flux = fluxQuanta*4*pi;
    int fluxSign = fluxQuanta == 0 ? 0 : fluxQuanta / abs(fluxQuanta);

    for (site.first(); site.test(); site.next())
    {
      std::vector<int> coords = {site.coord(0), site.coord(1), site.coord(2)};
      int xCoord = site.coord(0);
      int yCoord = site.coord(1);
      int zCoord = site.coord(2);

      int dir1 = (fieldDir + 1) % 3;
      int dir2 = (fieldDir + 2) % 3;

      // SU(2) gauge field
      for (int ii = 0; ii < 3; ii++)
      {
        std::vector<double> su2Vec = monsta::su2ToVec(Matrix(field, site, ii));
        if (ii == dir2)
        {
          su2Vec[2] -= 0.5*flux/pow(latSize[dir1],2)*(coords[dir1] - 0.5);
          if (coords[dir1] > 0)
          {
            su2Vec[2] += 0.5*flux/latSize[dir1];
          }
        }
        if (ii == dir1 && coords[dir1] == 0)
        {
          su2Vec[2] += 0.5*coords[dir2]*flux/latSize[dir2];
        }
        monsta::Matrix su2Mat = monsta::vecToSu2(su2Vec);
        field(site, ii, 0, 0) = su2Mat(0, 0);
        field(site, ii, 0, 1) = su2Mat(0, 1);
        field(site, ii, 1, 0) = su2Mat(1, 0);
        field(site, ii, 1, 1) = su2Mat(1, 1);
      }

      // U(1) gauge field
      for (int ii = 0; ii < 3; ii++)
      {
        double u1Angle = arg(theory.getU1Link(field, site, ii));
        if (ii == dir2)
        {
          u1Angle += 0.5*flux/pow(latSize[dir1],2)*(coords[dir1] - 0.5);
          if (coords[dir1] > 0)
          {
            u1Angle -= 0.5*flux/latSize[dir1];
          }
        }
        if (ii == dir1 && coords[dir1] == 0)
        {
          u1Angle -= 0.5*coords[dir2]*flux/latSize[dir2];
        }
        std::complex<double> u1Field = cos(u1Angle) + 1i*sin(u1Angle);
        field(site, 3, (ii + 1) % 2, (ii + 1) / 2) = u1Field;
      }
    }

    theory.applyBoundaryConditions(field);
  }

  void scaleVev(LATfield2::Field< std::complex<double> > &field, monsta::ElectroweakTheory &theory)
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

  void scaleVev(LATfield2::Field< std::complex<double> > &field, monsta::ElectroweakTheory &theory, double oldVev)
  {
    LATfield2::Site site(field.lattice());
    
    double higgsRatio = theory.getVev() / oldVev;

    for (site.first(); site.test(); site.next())
    {
      field(site, 3, 0, 0) *= higgsRatio;
    }
    theory.applyBoundaryConditions(field);
  }

  void linearSuperpose(LATfield2::Field< std::complex<double> > &field1,
  LATfield2::Field< std::complex<double> > &field2,
  LATfield2::Field< std::complex<double> > &outputField,
  monsta::ElectroweakTheory &theory)
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
      for (int ii = 0; ii < 3; ii++)
      {
        double u1Angle1 = arg(theory.getU1Link(field1, site, ii));
        double u1Angle2 = arg(theory.getU1Link(field2, site, ii));

        double outputAngle = u1Angle1 + u1Angle2;

        std::complex<double> outputU1 = cos(outputAngle) + 1i*sin(outputAngle);
        outputField(site, 3, (ii + 1) % 2, (ii + 1) / 2) = outputU1;
      }

      if (abs(field1(site, 3, 0, 0)) < abs(field2(site, 3, 0, 0)))
      {
        outputField(site, 3, 0, 0) = field1(site, 3, 0, 0);
      }
      else
      {
        outputField(site, 3, 0, 0) = field2(site, 3, 0, 0);
      }
    }

    theory.applyBoundaryConditions(outputField);
  }

  double getCSNum(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, monsta::ElectroweakTheory &theory)
  {
    double csNum = 0;

    // Compute Tr[W \wedge dW]
    for (int spatialIdx1 = 0; spatialIdx1 < 3; spatialIdx1++)
    {
      for (int spatialIdx2 = 0; spatialIdx2 < 3; spatialIdx2++)
      {
        for (int spatialIdx3 = 0; spatialIdx3 < 3; spatialIdx3++)
        {
          if (epsilonTensor(spatialIdx1, spatialIdx2, spatialIdx3) == 0) { continue; }
          Matrix su2Plaquette = theory.getSu2Plaquette(field, site, spatialIdx1, spatialIdx2);
          std::complex<double> determinant = su2Plaquette(0,0)*su2Plaquette(1, 1) - su2Plaquette(0, 1)*su2Plaquette(1, 0);
          for (int rowIdx = 0; rowIdx < 2; rowIdx++)
          {
            for (int colIdx = 0; colIdx < 2; colIdx++)
            {
              su2Plaquette(rowIdx, colIdx) = su2Plaquette(rowIdx, colIdx) / sqrt(determinant);
            }
            su2Plaquette(0, 0) = conj(su2Plaquette(1, 1));
            su2Plaquette(1, 0) = -1.0*conj(su2Plaquette(0, 1));
          }
          std::vector<double> fieldStrengthVec = su2ToVec(su2Plaquette);
          std::vector<double> gaugeFieldVec = su2ToVec(theory.getSu2Link(field, site, spatialIdx3));
          // theory.getSu2Plaquette(field, site, spatialIdx1, spatialIdx2).PRINT();
          // COUT << fieldStrengthVec[0] << fieldStrengthVec[1] << fieldStrengthVec[2] << endl;
          for (int su2Idx = 0; su2Idx < 3; su2Idx++)
          {
            csNum -= epsilonTensor(spatialIdx1, spatialIdx2, spatialIdx3)*fieldStrengthVec[su2Idx]*gaugeFieldVec[su2Idx];
          }
        }
      }
    }

    // cout << csNum << endl;

    // Compute Tr[W \wedge W \wedge W]
    for (int spatialIdx1 = 0; spatialIdx1 < 3; spatialIdx1++)
    {
      for (int spatialIdx2 = 0; spatialIdx2 < 3; spatialIdx2++)
      {
        for (int spatialIdx3 = 0; spatialIdx3 < 3; spatialIdx3++)
        {
          std::vector< std::vector<double> > gaugeFieldVecs;
          gaugeFieldVecs.push_back(su2ToVec(theory.getSu2Link(field, site, spatialIdx1)));
          gaugeFieldVecs.push_back(su2ToVec(theory.getSu2Link(field, site, spatialIdx2)));
          gaugeFieldVecs.push_back(su2ToVec(theory.getSu2Link(field, site, spatialIdx3)));
          for (int su2Idx1 = 0; su2Idx1 < 3; su2Idx1++)
          {
            for (int su2Idx2 = 0; su2Idx2 < 3; su2Idx2++)
            {
              for (int su2Idx3 = 0; su2Idx3 < 3; su2Idx3++)
              {
                csNum -= 2.0/3.0 * epsilonTensor(spatialIdx1, spatialIdx2, spatialIdx3)*epsilonTensor(su2Idx1, su2Idx2, su2Idx3)
                  *gaugeFieldVecs[spatialIdx1][su2Idx1]*gaugeFieldVecs[spatialIdx2][su2Idx2]*gaugeFieldVecs[spatialIdx3][su2Idx3];
              }
            }
          }
        }
      }
    }

    return csNum;
  }

}

#endif