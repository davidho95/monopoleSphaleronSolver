#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include "LATfield2.hpp"
#include "Theory.hpp"

namespace monsta
{
  class System
  {
  public:
    System(const monsta::Theory &theory, LATfield2::Lattice &lattice);

    void setRandomField();

    void setConstantField(double fieldVal);

    void setFieldElement(LATfield2::Site site, int cpt, double value);

    double getFieldElement(LATfield2::Site site, int cpt);

    double getGradientElement(LATfield2::Site site, int cpt);

    double getEnergy() { return theory_.computeEnergy(field_); }

    void updateFieldHalo() { field_.updateHalo(); }

    void updateGradient() { theory_.updateGradient(field_); }
    
    double updateGradientReturnMax() { return theory_.updateGradientReturnMax(field_); }

    double updateGradientReturnStepChange() { return theory_.updateGradientReturnStepChange(field_); }

    LATfield2::Lattice & getLattice() { return lattice_; }

    const int getNumSpatialCpts() { return numSpatialCpts_; }

  private:
    LATfield2::Lattice &lattice_;
    LATfield2::Field<double> field_;
    const monsta::Theory &theory_;
    const int numSpatialCpts_;
  };

  System::System(const monsta::Theory &theory, LATfield2::Lattice &lattice)
  : theory_(theory), lattice_(lattice), numSpatialCpts_(theory.getNumSpatialCpts())
  {
    int totalNumCpts = 2 * numSpatialCpts_; // Gradients and field values stored together
    field_.initialize(lattice, totalNumCpts);
    field_.alloc();

    setConstantField(0);
  }

  void System::setRandomField()
  {
    LATfield2::Site site(lattice_);

    for (site.first(); site.test(); site.next())
    {
      for (int iCpt = 0; iCpt < numSpatialCpts_; iCpt++)
      {
        field_(site, iCpt) = (std::rand() % 10) / 10.0;
      }
    }
    field_.updateHalo();
    updateGradient();
  }

  void System::setConstantField(double fieldVal)
  {
    LATfield2::Site site(lattice_);

    for (site.first(); site.test(); site.next())
    {
      for (int iCpt = 0; iCpt < numSpatialCpts_; iCpt++)
      {
        field_(site, iCpt) = fieldVal;
      }
    }
    field_.updateHalo();
    updateGradient();
  }

  void System::setFieldElement(LATfield2::Site site, int cpt, double value)
  {
    field_(site, cpt) = value;
  }

  double System::getFieldElement(LATfield2::Site site, int cpt)
  {
    return field_(site, cpt);
  }

  double System::getGradientElement(LATfield2::Site site, int cpt)
  {
    return field_(site, cpt + numSpatialCpts_);
  }
}

#endif