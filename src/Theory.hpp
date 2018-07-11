#ifndef THEORY_HPP
#define THEORY_HPP

#include "LATfield2.hpp"
#include <cmath>

namespace monsta
{
  class Theory
  {
  public:
    Theory();

    // virtual bool validateField(LATfield2::Field<double> &field) const = 0;
    double computeEnergy(LATfield2::Field<double> &field) const;
    void updateGradient(LATfield2::Field<double> &field) const;
    double updateGradientReturnMax(LATfield2::Field<double> &field) const;
    double updateGradientReturnStepChange(LATfield2::Field<double> &field) const;
    int getNumSpatialCpts() const { return numSpatialCpts_; }
    void setNumSpatialCpts(int numSpatialCpts) { numSpatialCpts_ = numSpatialCpts; }


  protected:
    int numSpatialCpts_ = 1;

    virtual double getLocalEnergyDensity(LATfield2::Field<double> &field, LATfield2::Site &site) const = 0;
    virtual double getLocalGradient(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt) const = 0;
  };

  Theory::Theory() {}

  double Theory::computeEnergy(LATfield2::Field<double> &field) const
  {
    LATfield2::Site site(field.lattice());

    double E = 0;
    for (site.first(); site.test(); site.next())
    {
      E += getLocalEnergyDensity(field, site);
    }

    return E;
  }

  void Theory::updateGradient(LATfield2::Field<double> &field) const
  {
    LATfield2::Site site(field.lattice());

    for (site.first(); site.test(); site.next())
    {
      for (int iCpt = 0; iCpt < numSpatialCpts_; iCpt++)
      {
        field(site, iCpt + numSpatialCpts_) = getLocalGradient(field, site, iCpt);
      }
    }
  }

  double Theory::updateGradientReturnMax(LATfield2::Field<double> &field) const
  {
    LATfield2::Site site(field.lattice());
    double maxGrad = 0;
    double localGradient;

    for (site.first(); site.test(); site.next())
    {
      for (int iCpt = 0; iCpt < numSpatialCpts_; iCpt++)
      {
        localGradient = getLocalGradient(field, site, iCpt);
        field(site, iCpt + numSpatialCpts_) = localGradient;
        if (abs(localGradient) > abs(maxGrad)) { maxGrad = localGradient; }
      }
    }

    return maxGrad;
  }

  double Theory::updateGradientReturnStepChange(LATfield2::Field<double> &field) const
  {
    LATfield2::Site site(field.lattice());
    double localGradientOld;
    double localGradient;

    double stepChangeNumerator = 0;
    double stepChangeDenominator = 0;

    for (site.first(); site.test(); site.next())
    {
      for (int iCpt = 0; iCpt < numSpatialCpts_; iCpt++)
      {
        localGradientOld = field(site, iCpt + numSpatialCpts_);
        localGradient = getLocalGradient(field, site, iCpt);

        field(site, iCpt + numSpatialCpts_) = localGradient;

        stepChangeNumerator += localGradientOld * (localGradientOld - localGradient);
        stepChangeDenominator += pow((localGradient - localGradientOld),2);
      }
    }

    return stepChangeNumerator / stepChangeDenominator;
  }
}


#endif