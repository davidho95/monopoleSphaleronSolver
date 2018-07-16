#ifndef MONOPOLESCALARPLANETHEORY_HPP
#define MONOPOLESCALARPLANETHEORY_HPP

#include "LATfield2.hpp"
#include "Theory.hpp"
#include <cmath>

namespace monsta
{
  class MonopoleScalarPlaneTheory: public Theory
  {
  public:
    MonopoleScalarPlaneTheory(double magneticCharge, int *southPos, int *northPos,
      double electricCharge, double vev, double selfCoupling, int *scalarPlaneZCoords);

  private:
    double magneticCharge_;
    int *southPos_;
    int *northPos_;
    double electricCharge_;
    double vev_;
    double selfCoupling_;
    int *scalarPlaneZCoords_;

    double getLocalEnergyDensity(LATfield2::Field<double> &field, LATfield2::Site &site,
      LATfield2::Field<double> &sinCosField) const;
    double getLocalGradient(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt,
      LATfield2::Field<double> &sinCosField) const;
    void generatePrecomputedValues(LATfield2::Field<double> &field, LATfield2::Field<double> &destinationField) const;
    bool isScalarFieldCoord(int zCoord) const;
  };

  MonopoleScalarPlaneTheory::MonopoleScalarPlaneTheory(
    double magneticCharge, int *southPos, int *northPos,
    double electricCharge, double vev, double selfCoupling,
    int* scalarPlaneZCoords)
  : Theory(), magneticCharge_(magneticCharge), southPos_(southPos), northPos_(northPos),
    electricCharge_(electricCharge), vev_(vev), selfCoupling_(selfCoupling),
    scalarPlaneZCoords_(scalarPlaneZCoords)
  {
    setNumSpatialCpts(5);
    setNumPrecomputedCpts(6);
    if (southPos_[1] != northPos_[1] || southPos_[2] != northPos_[2])
    {
      throw std::invalid_argument("Monopoles must lie on a line parallel to the x axis.");
    }
  }

  void MonopoleScalarPlaneTheory::generatePrecomputedValues(LATfield2::Field<double> &field,
    LATfield2::Field<double> &destinationField) const
  {
    int numGaugeCpts = 3;
    int totalNumCpts = numGaugeCpts * 2;

    LATfield2::Site site(field.lattice());

    for(site.first(); site.test(); site.next())
    {
      int zCoord = site.coord(2);
      if (!isScalarFieldCoord(zCoord)) { continue; }
      for (int iCpt = 0; iCpt < numGaugeCpts; iCpt++)
      {
        destinationField(site, iCpt) = sin(electricCharge_*field(site, iCpt));
        destinationField(site, iCpt+numGaugeCpts) = cos(electricCharge_*field(site, iCpt));
      }
    }
    destinationField.updateHalo();
  }

  double MonopoleScalarPlaneTheory::getLocalEnergyDensity(LATfield2::Field<double> &field, LATfield2::Site &site,
    LATfield2::Field<double> &sinCosField) const
  {
    int xCoord = site.coord(0);
    int yCoord = site.coord(1);
    int zCoord = site.coord(2);
    double Bx = 0;
    if (yCoord == southPos_[1] && zCoord == southPos_[2])
    {
      if (xCoord >= southPos_[0] && xCoord < northPos_[0]) { Bx = magneticCharge_; }
    }

    double sinField[3];
    double cosField[3];

    // Precompute sines and cosines for speed
    for (int iCpt = 0; iCpt < 3; iCpt++)
    {
      sinField[iCpt] = sin(electricCharge_*field(site, iCpt));
      cosField[iCpt] = cos(electricCharge_*field(site, iCpt));
    }

    int numGaugeCpts = 3;
    double E = 0;
    E += (pow(field(site, 0) - field(site+1, 0) - field(site, 1) + field(site+0, 1),2) + 
        pow(field(site, 1) - field(site+2, 1) - field(site, 2) + field(site+1, 2) + Bx,2) + 
        pow(field(site, 0) - field(site+2, 0) - field(site, 2) + field(site+0, 2),2))/2.;
    if (!isScalarFieldCoord(zCoord)) { return E; }
    E += pow(sinCosField(site, 0+numGaugeCpts),2)*pow(field(site, 3),2) + 
        pow(sinCosField(site, 1+numGaugeCpts),2)*pow(field(site, 3),2) + 
        pow(sinCosField(site, 0),2)*pow(field(site, 3),2) + 
        pow(sinCosField(site, 1),2)*pow(field(site, 3),2);
    E += -2*sinCosField(site, 1+numGaugeCpts)*field(site, 3)*field(site+1, 3) + pow(field(site+1, 3),2) - 
        2*sinCosField(site, 0+numGaugeCpts)*field(site, 3)*field(site+0, 3) + pow(field(site+0, 3),2) - 
        2*sinCosField(site, 1)*field(site+1, 3)*field(site, 4) - 
        2*sinCosField(site, 0)*field(site+0, 3)*field(site, 4);
    E += pow(sinCosField(site, 0+numGaugeCpts),2)*pow(field(site, 4),2) + 
        pow(sinCosField(site, 1+numGaugeCpts),2)*pow(field(site, 4),2) +  
        pow(sinCosField(site, 0),2)*pow(field(site, 4),2) + 
        pow(sinCosField(site, 1),2)*pow(field(site, 4),2);
    E += selfCoupling_*pow(-pow(vev_,2) + pow(field(site, 3),2) + pow(field(site, 4),2),2);
    E += 2*sinCosField(site, 1)*field(site, 3)*field(site+1, 4) - 
        2*sinCosField(site, 1+numGaugeCpts)*field(site, 4)*field(site+1, 4) + pow(field(site+1, 4),2) + 
        2*sinCosField(site, 0)*field(site, 3)*field(site+0, 4) - 
        2*sinCosField(site, 0+numGaugeCpts)*field(site, 4)*field(site+0, 4) + pow(field(site+0, 4),2);
    return E;
  }

  double MonopoleScalarPlaneTheory::getLocalGradient(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt,
    LATfield2::Field<double> &sinCosField) const
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

    int numGaugeCpts = 3;
    double grad = 0;
    switch (cpt)
    {
      case 0:
        grad = -field(site-1, 0) - field(site-2, 0) + 4*field(site, 0) - field(site+2,0)
          - field(site+1, 0) + field(site-1, 1) - field(site, 1) - field(site+0-1, 1)
          + field(site+0, 1) + field(site-2, 2) - field(site, 2) - field(site+0-2, 2)
          + field(site+0, 2);
        if (!isScalarFieldCoord(zCoord)) { break; }
        grad += 2*electricCharge_*sinCosField(site, 0)*field(site, 3)*field(site+0, 3)
          - 2*electricCharge_*sinCosField(site, 0+numGaugeCpts)*field(site+0, 3)*field(site, 4)
          + 2*electricCharge_*sinCosField(site, 0+numGaugeCpts)*field(site, 3)*field(site+0, 4)
          + 2*electricCharge_*sinCosField(site, 0)*field(site, 4)*field(site+0, 4);
        break;
      case 1:
        grad = -BxShiftZ + Bx + field(site-0, 0) - field(site-0+1, 0) - field(site, 0) + field(site+1, 0)
          - field(site-0, 1) - field(site-2, 1) + 4*field(site, 1) - field(site+2, 1)
          - field(site+0, 1) + field(site-2, 2) - field(site, 2) - field(site+1-2, 2)
          + field(site+1, 2);
          if (!isScalarFieldCoord(zCoord)) { break; }
        grad += 2*electricCharge_*sinCosField(site, 1)*field(site, 3)*field(site+1, 3)
          - 2*electricCharge_*sinCosField(site, 1+numGaugeCpts)*field(site+1, 3)*field(site, 4)
          + 2*electricCharge_*sinCosField(site, 1+numGaugeCpts)*field(site, 3)*field(site+1, 4)
          + 2*electricCharge_*sinCosField(site, 1)*field(site, 4)*field(site+1, 4);
        break;
      case 2:
        grad = BxShiftY - Bx + field(site-0, 0) - field(site-0+2, 0) - field(site, 0) + field(site+2, 0)
          + field(site-1, 1) - field(site-1+2, 1) - field(site, 1) + field(site+2, 1)
          - field(site-0, 2) - field(site-1, 2) + 4*field(site, 2) - field(site+1, 2)
          - field(site+0, 2);
        break;
      case 3:
        if (!isScalarFieldCoord(zCoord)) { break; }
        grad = -2*(sinCosField(site-0, 0+numGaugeCpts)*field(site-0, 3) + sinCosField(site-1, 1+numGaugeCpts)*field(site-1, 3)
          - 4*field(site,3)
          + 2*selfCoupling_*pow(vev_,2)*field(site, 3) - 2*selfCoupling_*pow(field(site,3),3)
          + sinCosField(site,1+numGaugeCpts)*field(site+1, 3)
          + sinCosField(site, 0+numGaugeCpts)*field(site+0, 3) + sinCosField(site-0, 0)*field(site-0, 4)
          + sinCosField(site-1, 1)*field(site-1, 4)
          - 2*selfCoupling_*field(site,3)*pow(field(site, 4),2)
          - sinCosField(site, 1)*field(site+1, 4) - sinCosField(site, 0)*field(site+0, 4));
        break;
      case 4:
        if (!isScalarFieldCoord(zCoord)) { break; }
        grad = 2*(sinCosField(site-0, 0)*field(site-0, 3) + sinCosField(site-1, 1)*field(site-1, 3)
          - sinCosField(site, 1)*field(site+1, 3) - sinCosField(site, 0)*field(site+0, 3)
          - sinCosField(site-0, 0+numGaugeCpts)*field(site-0, 4) - sinCosField(site-1, 1+numGaugeCpts)*field(site-1, 4)
          + 4*field(site, 4)
          - 2*selfCoupling_*pow(vev_,2)*field(site, 4) + 2*selfCoupling_*pow(field(site, 3),2)*field(site, 4)
          + 2*selfCoupling_*pow(field(site, 4),3)
          - sinCosField(site, 1+numGaugeCpts)*field(site+1, 4) - sinCosField(site, 0+numGaugeCpts)*field(site+0,4));
        break;
    }
    return grad;
  }

  bool MonopoleScalarPlaneTheory::isScalarFieldCoord(int zCoord) const
  {
    return (zCoord == scalarPlaneZCoords_[0] || zCoord == scalarPlaneZCoords_[1]);
  }
}

#endif