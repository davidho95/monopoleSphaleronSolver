#ifndef AXISYMMETRICFIELDTOOLS_HPP
#define AXISYMMETRICFIELDTOOLS_HPP

#include "LATfield2.hpp"
#include <complex>
#include "../Su2Tools.hpp"
#include "GeorgiGlashowSu2Theory2d.hpp"

namespace monsta
{
  void setVacuumField(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2TheoryAxisymmetric &theory)
  {
    double vev = theory.getVev();

    LATfield2::Site site(field.lattice());

    for (site.first(); site.test(); site.next())
    {
      for (int matIdx = 0; matIdx < 3; matIdx++)
      {
        theory.setSu2Link(field, site, matIdx, identity);
      }
      theory.setHiggsMagnitude(field, site, vev / sqrt(2));
    }
    theory.applyBoundaryConditions(field);
  }

  void setVacuumField(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2TheoryAxisymmetric2d &theory)
  {
    double vev = theory.getVev();

    LATfield2::Site site(field.lattice());

    for (site.first(); site.test(); site.next())
    {
      for (int matIdx = 0; matIdx < 2; matIdx++)
      {
        theory.setSu2Link(field, site, matIdx, identity);
      }
      theory.setSu2Link(field, site, 2, Matrix(2));
      theory.setHiggsMagnitude(field, site, vev / sqrt(2));
    }
    theory.applyBoundaryConditions(field);
  }

  void setRandomField(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2TheoryAxisymmetric &theory)
  {
    double vev = theory.getVev();

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
        theory.setSu2Link(field, site, matIdx, su2Mat);
      }
      theory.setHiggsMagnitude(field, site, 0.01*double(rand() % 100));
    }
    theory.applyBoundaryConditions(field);
  }

  void setRandomField(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2TheoryAxisymmetric2d &theory)
  {
    double vev = theory.getVev();

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
        if (matIdx == 2)
        {
          theory.setSu2Link(field, site, matIdx, vecToSu2LieAlg(su2Vec));
          continue;
        }
        Matrix su2Mat = vecToSu2(su2Vec);
        theory.setSu2Link(field, site, matIdx, su2Mat);
      }
      theory.setHiggsMagnitude(field, site, 0.01*double(rand() % 100));
    }
    theory.applyBoundaryConditions(field);
  }

  void perturbField(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2TheoryAxisymmetric &theory)
  {
    double vev = theory.getVev();

    LATfield2::Site site(field.lattice());

    for (site.first(); site.test(); site.next())
    {
      for (int matIdx = 0; matIdx < 3; matIdx++)
      {
        Matrix linkVal = theory.getSu2Link(field, site, matIdx);
        std::vector<double> linkVec = su2ToVec(linkVal);
        std::vector<double> su2Vec(3);
        for (int ii = 0; ii < 3; ii++)
        {
          su2Vec[ii] = (1e-8*double(rand() % 628318 - 314159)) + linkVec[ii];
        }
        Matrix su2Mat = vecToSu2(su2Vec);
        theory.setSu2Link(field, site, matIdx, su2Mat);
      }
      std::complex<double> higgsVal = theory.getHiggsMagnitude(field, site);
      theory.setHiggsMagnitude(field, site, higgsVal + 0.00001*higgsVal*double(rand() % 100));
    }
    theory.applyBoundaryConditions(field);
  }

  void perturbField(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2TheoryAxisymmetric &theory, double perturbationMagnitude)
  {
    double vev = theory.getVev();

    LATfield2::Site site(field.lattice());

    for (site.first(); site.test(); site.next())
    {
      for (int matIdx = 0; matIdx < 3; matIdx++)
      {
        Matrix linkVal = theory.getSu2Link(field, site, matIdx);
        std::vector<double> linkVec = su2ToVec(linkVal);
        std::vector<double> su2Vec(3);
        for (int ii = 0; ii < 3; ii++)
        {
          su2Vec[ii] = (perturbationMagnitude*1e-3*double(rand() % 628318 - 314159)) + linkVec[ii];
        }
        Matrix su2Mat = vecToSu2(su2Vec);
        theory.setSu2Link(field, site, matIdx, su2Mat);
      }
      std::complex<double> higgsVal = theory.getHiggsMagnitude(field, site);
      theory.setHiggsMagnitude(field, site, higgsVal + perturbationMagnitude*higgsVal*double(rand() % 100));
    }
    theory.applyBoundaryConditions(field);
  }

  void perturbField(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2TheoryAxisymmetric2d &theory, double perturbationMagnitude)
  {
    double vev = theory.getVev();

    LATfield2::Site site(field.lattice());

    for (site.first(); site.test(); site.next())
    {
      for (int matIdx = 0; matIdx < 2; matIdx++)
      {
        Matrix linkVal = theory.getSu2Link(field, site, matIdx);
        std::vector<double> linkVec = su2ToVec(linkVal);
        std::vector<double> su2Vec(3);
        for (int ii = 0; ii < 3; ii++)
        {
          su2Vec[ii] = (perturbationMagnitude*1e-3*double(rand() % 628318 - 314159)) + linkVec[ii];
        }
        Matrix su2Mat = vecToSu2(su2Vec);
        theory.setSu2Link(field, site, matIdx, su2Mat);
      }
      Matrix aThetaVal = theory.getScalarField(field, site, 1);
      std::vector<double> lieAlgVec = su2LieAlgToVec(aThetaVal);
      std::vector<double> perturbedVec(3);
      for (int ii = 0; ii < 3; ii++)
      {
        perturbedVec[ii] = (perturbationMagnitude*1e-3*double(rand() % 628318 - 314159)) + lieAlgVec[ii];
      }
      Matrix perturbedMat = vecToSu2LieAlg(perturbedVec);
      theory.setSu2Link(field, site, 2, perturbedMat);

      std::complex<double> higgsVal = theory.getHiggsMagnitude(field, site);
      theory.setHiggsMagnitude(field, site, higgsVal + perturbationMagnitude*higgsVal*double(rand() % 100));
    }
    theory.applyBoundaryConditions(field);
  }

  void scaleVev(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2TheoryAxisymmetric &theory)
  {
    LATfield2::Site site(field.lattice());
    double oldVev = 0;
    for (site.first(); site.test(); site.next())
    {
      if (abs(theory.getHiggsMagnitude(field, site))*sqrt(2) > oldVev)
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

  // void addConstantMagneticField(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2TheoryAxisymmetric &theory)
  // {
  //   LATfield2::Site site(field.lattice());
  //   double gaugeCoupling = theory.getGaugeCoupling();

  //   int zSize = field.lattice().size(0);
  //   int rSize = field.lattice().size(1);

  //   double flux = fluxQuanta*4*monsta::pi;
  //   int fluxSign = fluxQuanta == 0 ? 0 : fluxQuanta / abs(fluxQuanta);

  //   for (site.first(); site.test(); site.next())
  //   {
  //     int zCoord = site.coord(0);
  //     int rCoord = site.coord(1);

  //     for (int ii = 0; ii < 3; ii++)
  //     {
  //       std::vector<double> su2Vec = monsta::su2ToVec(Matrix(field, site, ii));
  //       if (ii == 2)
  //       {
  //         su2Vec[2] += 0.5*flux/pow(zSize,2)*(yCoord - 0.5);
  //         if (yCoord > 0)
  //         {
  //           su2Vec[2] -= 0.5*flux/ySize;
  //         }
  //       }
  //       if (ii == 1 && yCoord == 0)
  //       {
  //         su2Vec[2] -= 0.5*zCoord*flux/zSize;
  //       }
  //       monsta::Matrix su2Mat = monsta::vecToSu2(su2Vec);
  //       field(site, ii, 0, 0) = su2Mat(0, 0);
  //       field(site, ii, 0, 1) = su2Mat(0, 1);
  //       field(site, ii, 1, 0) = su2Mat(1, 0);
  //       field(site, ii, 1, 1) = su2Mat(1, 1);
  //     }
  //   }

  //   theory.applyBoundaryConditions(field);
  // }
}

#endif