#ifndef THEORY_HPP
#define THEORY_HPP

#include "LATfield2.hpp"
#include "./Matrix.hpp"
#include <cmath>

namespace monsta
{
  class Theory
  {
  public:
    Theory(int numFieldMatrices, int numRows, int numCols);

    double computeEnergy(LATfield2::Field< std::complex<double> > &field) const;
    void updateGradient(LATfield2::Field< std::complex<double> > &field) const;
    std::complex<double> updateGradientReturnMax(LATfield2::Field< std::complex<double> > &field) const;
    double updateGradientReturnStepChange(LATfield2::Field< std::complex<double> > &field) const;
    int getnumFieldMatrices() const { return numFieldMatrices_; }
    int getNumRows() const { return numRows_; }
    virtual void applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const { field.updateHalo(); }
    virtual double getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const = 0;
    virtual std::complex<double> getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const = 0;
    virtual monsta::Matrix getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const = 0;
    virtual void postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const = 0;

  protected:
    int numFieldMatrices_ = 1;
    int numRows_ = 1;
    int numCols_ = 1;

  };

  Theory::Theory(int numFieldMatrices, int numRows, int numCols) : numFieldMatrices_(numFieldMatrices), numRows_(numRows), numCols_(numCols)
  {}

  double Theory::computeEnergy(LATfield2::Field< std::complex<double> > &field) const
  {
    LATfield2::Lattice &lattice = field.lattice();
    LATfield2::Site site(lattice);

    double E = 0;
    for (site.first(); site.test(); site.next())
    {
      E += getLocalEnergyDensity(field, site);
    }

    parallel.sum(E);
    return E;
  }

  void Theory::updateGradient(LATfield2::Field< std::complex<double> > &field) const
  {
    LATfield2::Lattice &lattice = field.lattice();
    LATfield2::Site site(lattice);

    for (site.first(); site.test(); site.next())
    {
      for (int matIdx = 0; matIdx < numFieldMatrices_; matIdx++)
      {
        for (int rowIdx = 0; rowIdx < numRows_; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols_; colIdx++)
          {
            field(site, matIdx + numFieldMatrices_, rowIdx, colIdx) = getLocalGradient(field, site, matIdx, rowIdx, colIdx);
          }
        }
      }
    }
  }

  std::complex<double> Theory::updateGradientReturnMax(LATfield2::Field< std::complex<double> > &field) const
  {
    LATfield2::Lattice &lattice = field.lattice();
    LATfield2::Site site(lattice);
    double maxGrad = 0;
    std::complex<double> localGradient;

    for (site.first(); site.test(); site.next())
    {
      for (int matIdx = 0; matIdx < numFieldMatrices_; matIdx++)
      {
        for (int rowIdx = 0; rowIdx < numRows_; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols_; colIdx++)
          {
            localGradient = getLocalGradient(field, site, matIdx, rowIdx, colIdx);
            field(site, matIdx + numFieldMatrices_, rowIdx, colIdx) = localGradient;
            if (abs(localGradient) > abs(maxGrad)) { maxGrad = abs(localGradient); }
          }
        }
      }
    }

    parallel.max(maxGrad);
    return maxGrad;
  }

  double Theory::updateGradientReturnStepChange(LATfield2::Field< std::complex<double> > &field) const
  {
    LATfield2::Lattice &lattice = field.lattice();
    LATfield2::Site site(lattice);
    std::complex<double> localGradientOld;
    std::complex<double> localGradient;

    double stepChangeNumerator = 0;
    double stepChangeDenominator = 0;

    for (site.first(); site.test(); site.next())
    {
      for (int matIdx = 0; matIdx < numFieldMatrices_; matIdx++)
      {
        for (int rowIdx = 0; rowIdx < numRows_; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols_; colIdx++)
          {
            localGradientOld = field(site, matIdx + numFieldMatrices_, rowIdx, colIdx);
            localGradient = getLocalGradient(field, site, matIdx, rowIdx, colIdx);

            field(site, matIdx + numFieldMatrices_, rowIdx, colIdx) = localGradient;

            stepChangeNumerator += abs(localGradientOld * (localGradientOld - localGradient));
            stepChangeDenominator += pow(abs((localGradient - localGradientOld)),2);
          }
        }
      }
    }

    parallel.sum(stepChangeNumerator);
    parallel.sum(stepChangeDenominator);
    return stepChangeNumerator / stepChangeDenominator;
  }
}


#endif