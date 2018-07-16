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

    double computeEnergy(LATfield2::Field<double> &field) const;
    void updateGradient(LATfield2::Field<double> &field) const;
    double updateGradientReturnMax(LATfield2::Field<double> &field) const;
    double updateGradientReturnStepChange(LATfield2::Field<double> &field) const;
    int getNumSpatialCpts() const { return numSpatialCpts_; }
    void setNumSpatialCpts(int numSpatialCpts) { numSpatialCpts_ = numSpatialCpts; }
    void setNumPrecomputedCpts(int numPrecomputedCpts) { numPrecomputedCpts_ = numPrecomputedCpts; }


  protected:
    int numSpatialCpts_ = 1;
    int numPrecomputedCpts_ = 1;

    virtual double getLocalEnergyDensity(LATfield2::Field<double> &field, LATfield2::Site &site,
      LATfield2::Field<double> &precomputedValues) const = 0;
    virtual double getLocalGradient(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt,
      LATfield2::Field<double> &precomputedValues) const = 0;
    virtual void generatePrecomputedValues(LATfield2::Field<double> &field,
      LATfield2::Field<double> &destinationField) const;
  };

  Theory::Theory() {}

  double Theory::computeEnergy(LATfield2::Field<double> &field) const
  {
    LATfield2::Lattice &lattice = field.lattice();
    LATfield2::Site site(lattice);

    LATfield2::Field<double> precomputedValues(lattice, numPrecomputedCpts_);
    generatePrecomputedValues(field, precomputedValues);

    double E = 0;
    for (site.first(); site.test(); site.next())
    {
      E += getLocalEnergyDensity(field, site, precomputedValues);
    }

    parallel.sum(E);
    return E;
  }

  void Theory::updateGradient(LATfield2::Field<double> &field) const
  {
    LATfield2::Lattice &lattice = field.lattice();
    LATfield2::Site site(lattice);

    LATfield2::Field<double> precomputedValues(lattice, numPrecomputedCpts_);
    generatePrecomputedValues(field, precomputedValues);

    for (site.first(); site.test(); site.next())
    {
      for (int iCpt = 0; iCpt < numSpatialCpts_; iCpt++)
      {
        field(site, iCpt + numSpatialCpts_) = getLocalGradient(field, site, iCpt, precomputedValues);
      }
    }
  }

  double Theory::updateGradientReturnMax(LATfield2::Field<double> &field) const
  {
    LATfield2::Lattice &lattice = field.lattice();
    LATfield2::Site site(lattice);
    double maxGrad = 0;
    double localGradient;

    LATfield2::Field<double> precomputedValues(lattice, numPrecomputedCpts_);
    generatePrecomputedValues(field, precomputedValues);

    for (site.first(); site.test(); site.next())
    {
      for (int iCpt = 0; iCpt < numSpatialCpts_; iCpt++)
      {
        localGradient = getLocalGradient(field, site, iCpt, precomputedValues);
        field(site, iCpt + numSpatialCpts_) = localGradient;
        if (abs(localGradient) > abs(maxGrad)) { maxGrad = localGradient; }
      }
    }

    parallel.max(maxGrad);
    return maxGrad;
  }

  double Theory::updateGradientReturnStepChange(LATfield2::Field<double> &field) const
  {
    LATfield2::Lattice &lattice = field.lattice();
    LATfield2::Site site(lattice);
    double localGradientOld;
    double localGradient;

    double stepChangeNumerator = 0;
    double stepChangeDenominator = 0;

    LATfield2::Field<double> precomputedValues(lattice, numPrecomputedCpts_);
    generatePrecomputedValues(field, precomputedValues);

    for (site.first(); site.test(); site.next())
    {
      for (int iCpt = 0; iCpt < numSpatialCpts_; iCpt++)
      {
        localGradientOld = field(site, iCpt + numSpatialCpts_);
        localGradient = getLocalGradient(field, site, iCpt, precomputedValues);

        field(site, iCpt + numSpatialCpts_) = localGradient;

        stepChangeNumerator += localGradientOld * (localGradientOld - localGradient);
        stepChangeDenominator += pow((localGradient - localGradientOld),2);
      }
    }

    parallel.sum(stepChangeNumerator);
    parallel.sum(stepChangeDenominator);
    return stepChangeNumerator / stepChangeDenominator;
  }

  void Theory::generatePrecomputedValues(LATfield2::Field<double> &field,
    LATfield2::Field<double> &destinationField) const {};
}


#endif