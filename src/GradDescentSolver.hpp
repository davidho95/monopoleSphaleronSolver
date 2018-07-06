#ifndef GRADDESCENTSOLVER_HPP
#define GRADDESCENTSOLVER_HPP

#include <cstdio>
#include "LATfield2.hpp"

namespace monsta {
  template<class FieldType>
  class GradDescentSolver {
  public:
    bool solnFound_ = false;
    LATfield2::Field<FieldType> field_;
    double energy_ = 0;
    LATfield2::Field<FieldType> gradient_;
    double maxGrad_ = 1e6;

    GradDescentSolver(LATfield2::Lattice &lattice, int numCpts, double tol, int maxIterations)
    : lattice_(lattice), numCpts_(numCpts), tol_(tol), maxIterations_(maxIterations)
    {
      field_.initialize(lattice_, numCpts_);
      field_.alloc();
      gradient_.initialize(lattice_, numCpts_);
      gradient_.alloc();
      prevField_.initialize(lattice_, numCpts_);
      prevField_.alloc();
    }

    void setVerbosity(bool isVerbose) { verbose_ = isVerbose; }

    // Initialise field to random array of doubles
    void setInitialField()
    {
      LATfield2::Site site(lattice_);

      for (site.first(); site.test(); site.next())
      {
        for (int iCpt = 0; iCpt < numCpts_; iCpt++)
        {
          field_(site, iCpt) = (std::rand() % 10) / 10.0;
        }
      }
      field_.updateHalo();
      updateEnergy();
      updateGradient();
      initialized_ = true;
    }

    // Initialise field to array of constant values
    void setInitialField(double initVal)
    {
      LATfield2::Site site(lattice_);

      for (site.first(); site.test(); site.next())
      {
        for (int iCpt = 0; iCpt < numCpts_; iCpt++)
        {
          field_(site, iCpt) = initVal;
        }
      }
      field_.updateHalo();
      updateEnergy();
      updateGradient();
      initialized_ = true;
    }

    // Initialise to copy of supplied field
    void setInitialField(const LATfield2::Field<FieldType> &initField)
    {
      LATfield2::Site site(lattice_);

      for (site.first(); site.test(); site.next())
      {
        for (int iCpt = 0; iCpt < numCpts_; iCpt++)
        {
          field_(site, iCpt) = initField(site, iCpt);
        }
      }
      field_.updateHalo();
      updateEnergy();
      updateGradient();
      initialized_ = true;
    }

    void solve() {
      int numIters = 0;
      while (abs(maxGrad_) > tol_ && numIters < maxIterations_) {
        numIters++;
        iterate();
        if (verbose_)
        {
          cout << maxGrad_ << endl;
        }
      }
      updateEnergy();

      if (numIters < maxIterations_) {
        solnFound_ = true;
        std::cout << "Gradient descent finished in " << numIters << " iterations." << std::endl;
        std::cout << "Minimum energy: " << energy_ << std::endl;
      } else {
        std::cout << "Gradient descent aborted after " << maxIterations_ << " iterations." << std::endl;
        std::cout << "Maximum gradient: " << maxGrad_ << std::endl;
        std::cout << "Energy reached: " << energy_ << std::endl;
      }
    }

  protected:
    LATfield2::Lattice lattice_;
    int numCpts_;
    int maxIterations_;
    double tol_;
    double stepSize_ = 0.01;
    bool initialized_ = false;
    bool verbose_ = true;
    bool storingPrevEnergy_ = true;
    double prevEnergy_;
    LATfield2::Field<FieldType> prevField_;

    virtual double getLocalEnergyDensity(LATfield2::Field<FieldType> &field, LATfield2::Site &site) const = 0;
    virtual double getLocalGradient(LATfield2::Field<FieldType> &field, LATfield2::Site &site, int cpt) const = 0;

    void updateEnergy()
    {
      LATfield2::Site site(lattice_);
      energy_ = 0;
      for (site.first(); site.test(); site.next())
      {
        energy_ += getLocalEnergyDensity(field_, site);
      }
    }

    virtual void updateGradient()
    {
      LATfield2::Site site(lattice_);
      double maxGrad = 0;

      double gradVal;
      double prevGradVal;
      double stepSizeNumerator = 0;
      double stepSizeDenominator = 0;

      for (site.first(); site.test(); site.next())
      {
        for (int iCpt = 0; iCpt < numCpts_; iCpt++)
        {
          if (initialized_) {
            prevGradVal = gradient_(site, iCpt);
          }
          gradVal = getLocalGradient(field_, site, iCpt);
          gradient_(site, iCpt) = gradVal;
          if (abs(gradVal) > maxGrad) { maxGrad = abs(gradVal); }
          if (initialized_) {
            stepSizeNumerator += abs((field_(site, iCpt) - prevField_(site, iCpt)) * (gradient_(site, iCpt) - prevGradVal));
            stepSizeDenominator += pow(gradient_(site, iCpt) - prevGradVal, 2);
          }
        }
      }
      gradient_.updateHalo();

      if (initialized_)
      {
          stepSize_ = stepSizeNumerator / stepSizeDenominator;
      }

      maxGrad_ = maxGrad;
    }

    // Performs a single gradient descent iteration, returning the maximum gradient
    virtual void iterate()
    {
      LATfield2::Site site(lattice_);

      // Update entire field in-place
      double gradVal;
      for (site.first(); site.test(); site.next())
      {
        for (int iCpt = 0; iCpt < numCpts_; iCpt++)
        {
          if (initialized_) {
            prevField_(site, iCpt) = field_(site, iCpt);
          }
          field_(site, iCpt) -= stepSize_ * gradient_(site, iCpt);
        }
      }
      field_.updateHalo();
      prevField_.updateHalo();
      if (storingPrevEnergy_)
      {
        prevEnergy_ = energy_;
        updateEnergy();
      }
      updateGradient();
    }

  };
}

#endif