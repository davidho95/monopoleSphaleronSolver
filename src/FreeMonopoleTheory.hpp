#ifndef FREEMONOPOLETHEORY_HPP
#define FREEMONOPOLETHEORY_HPP

#include "LATfield2.hpp"
#include "Theory.hpp"
#include <cmath>

namespace monsta
{
  class FreeMonopoleTheory: public Theory
  {
  public:
    FreeMonopoleTheory(double magneticCharge, int *southPos, int *northPos);

  private:
    double magneticCharge_;
    int *southPos_;
    int *northPos_;

    double getLocalEnergyDensity(LATfield2::Field<double> &field, LATfield2::Site &site) const;
    double getLocalGradient(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt) const;
  };

  FreeMonopoleTheory::FreeMonopoleTheory(double magneticCharge, int *southPos, int *northPos)
  : Theory(), magneticCharge_(magneticCharge), southPos_(southPos), northPos_(northPos)
  {
    setNumSpatialCpts(3);
    if (southPos_[1] != northPos_[1] || southPos_[2] != northPos_[2])
    {
      throw std::invalid_argument("Monopoles must lie on a line parallel to the x axis.");
    }
  }

  double FreeMonopoleTheory::getLocalEnergyDensity(LATfield2::Field<double> &field, LATfield2::Site &site) const
  {
    int xCoord = site.coord(0);
    int yCoord = site.coord(1);
    int zCoord = site.coord(2);
    double Bx = 0;
    if (yCoord == southPos_[1] && zCoord == southPos_[2])
    {
      if (xCoord >= southPos_[0] && xCoord < northPos_[0]) { Bx = magneticCharge_; }
    }

    double E = 0;
      E += pow(field(site+1, 2) - field(site, 2) - field(site+2, 1) + field(site,1) + Bx, 2);
      E += pow(field(site+2, 0) - field(site, 0) - field(site+0, 2) + field(site,2), 2);
      E += pow(field(site+0, 1) - field(site, 1) - field(site+1, 0) + field(site,0), 2);
    return 0.5*E;
  }

  double FreeMonopoleTheory::getLocalGradient(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt) const
  {
    int xCoord = site.coord(0);
    int yCoord = site.coord(1);
    int zCoord = site.coord(2);
    double Bx = 0;
    double BxShiftY = 0;
    double BxShiftZ = 0;
    if (xCoord >= southPos_[0] && xCoord < northPos_[0])
    {
      if (yCoord == southPos_[1] && zCoord == southPos_[2]) { Bx = magneticCharge_; }
      if (yCoord == southPos_[1] + 1 && zCoord == southPos_[2]) { BxShiftY = magneticCharge_; }
      if (yCoord == southPos_[1] && zCoord == southPos_[2] + 1) { BxShiftZ = magneticCharge_; }
    }

    double grad;
    switch (cpt)
    {
      case 0:
        grad = -field(site-1, 0) - field(site-2, 0) + 4*field(site, 0) - field(site+2, 0) - 
          field(site+1, 0) + field(site-1, 1) - field(site, 1) - field(site+0-1, 1) + 
          field(site+0, 1) + field(site-2, 2) - field(site, 2) - field(site+0-2, 2) + 
          field(site+0, 2);
        break;
      case 1:
        grad = -BxShiftZ + Bx + field(site-0, 0) - field(site-0+1, 0) - field(site, 0) + field(site+1, 0) - 
          field(site-0, 1) - field(site-2, 1) + 4*field(site, 1) - field(site+2, 1) - 
          field(site+0, 1) + field(site-2, 2) - field(site, 2) - field(site+1-2, 2) + 
          field(site+1, 2);
        break;
      case 2:
        grad = BxShiftY - Bx + field(site-0, 0) - field(site-0+2, 0) - field(site, 0) + field(site+2, 0) + 
          field(site-1, 1) - field(site-1+2, 1) - field(site, 1) + field(site+2, 1) - 
          field(site-0, 2) - field(site-1, 2) + 4*field(site, 2) - field(site+1, 2) - 
          field(site+0, 2);
          break;
    }
    return grad;
  }
}

#endif