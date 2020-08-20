#ifndef INSTANTONFIELDTOOLS_HPP
#define INSTANTONFIELDTOOLS_HPP

#include "LATfield2.hpp"
#include <complex>
#include "../Su2Tools.hpp"
#include "GeorgiGlashowSu2Theory4d.hpp"


namespace monsta
{

  void setVacuumField(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2Theory4d &theory)
  {
    double vev = theory.getVev();

    LATfield2::Site site(field.lattice());

    for (site.first(); site.test(); site.next())
    {
      for (int matIdx = 0; matIdx < 3; matIdx++)
      {
        theory.setSu2Link(field, site, matIdx, identity);
      }
      theory.setHiggsMagnitude(field, site, vev/sqrt(2));
    }
    theory.applyBoundaryConditions(field);
  }

  void setRandomField(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2Theory4d &theory)
  {
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
        monsta::Matrix su2Mat = monsta::vecToSu2(su2Vec);
        theory.setSu2Link(field, site, matIdx, su2Mat);
      }
      theory.setHiggsMagnitude(field, site, 0.01*double(rand() % 100));
    }
    theory.applyBoundaryConditions(field);
  }

  void perturbField(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2Theory4d &theory, double perturbationMagnitude)
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
      double higgsVal = theory.getHiggsMagnitude(field, site);
      theory.setHiggsMagnitude(field, site, higgsVal + perturbationMagnitude*higgsVal*double(rand() % 100));
    }
    theory.applyBoundaryConditions(field);
  }


  void addConstantMagneticField(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2Theory4d &theory,
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
        theory.setSu2Link(field, site, ii, su2Mat);
      }
    }

    theory.applyBoundaryConditions(field);
  }


  
}

#endif