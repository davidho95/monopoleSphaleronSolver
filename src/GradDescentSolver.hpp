#ifndef GRADDESCENTSOLVER_HPP
#define GRADDESCENTSOLVER_HPP

#include <cstdio>
#include <functional>
#include "LATfield2.hpp"

namespace monsta {
  template<class FieldType>
  class GradDescentSolver {
  public:
    GradDescentSolver(double tol, int maxIterations){
      this->tol_ = tol;
      this->maxIterations_ = maxIterations;
    }

    // Computes energy of a given field
    double getEnergy(LATfield2::Field<FieldType> &field) const
    {
      LATfield2::Lattice &lattice = field.lattice();
      LATfield2::Site site(lattice);

      double energy = 0;
      for (site.first(); site.test(); site.next())
      {
        energy += this->energyDensity(field, site);
      }

      return energy;
    }

    LATfield2::Field<FieldType> getGradient(LATfield2::Field<FieldType> &field) const
    {
      LATfield2::Lattice &lattice = field.lattice();
      LATfield2::Site site(lattice);
      int numCpts = field.components();
      LATfield2::Field<FieldType> grad(lattice, numCpts);

      for (site.first(); site.test(); site.next())
      {
        for (int iCpt = 0; iCpt < numCpts; iCpt++)
        {
          grad(site, iCpt) = this->gradient(field, site, iCpt);
        }
      }

      return grad;
    }

    // Performs a single gradient descent iteration, returning the maximum gradient
    double iterate(LATfield2::Field<FieldType> &field)
    {
      LATfield2::Lattice &lattice = field.lattice();
      LATfield2::Site site(lattice);
      int numCpts = field.components();
      double maxGrad = 0;

      // First pass to get all gradients (not done in-place to preserve gauge invariance)
      LATfield2::Field<FieldType> grad = this->getGradient(field);

      // Update entire field in-place
      double gradVal;
      for (site.first(); site.test(); site.next())
      {
        for (int iCpt = 0; iCpt < numCpts; iCpt++)
        {
          gradVal = grad(site, iCpt);
          if (abs(gradVal) > maxGrad) { maxGrad = abs(gradVal); }
          field(site, iCpt) -= this->stepSize_ * gradVal;
        }
      }
      field.updateHalo();

      return maxGrad;
    }

    void solve(LATfield2::Field<FieldType> &field) {
      double bigNumber = 1e6;

      double energy = this->getEnergy(field);
      double energyOld;
      double energyChange = bigNumber;
      double maxGrad = bigNumber;

      int numIters = 0;

      while (maxGrad > this->tol_ && numIters < this->maxIterations_) {
        numIters++;
        energyOld = energy;
        maxGrad = this->iterate(field);
        energy = this->getEnergy(field);
        // cout << maxGrad << endl;
      }

      if (numIters < this->maxIterations_) {
        this->solnFound_ = true;
        std::cout << "Gradient descent finished in " << numIters << " iterations." << std::endl;
        std::cout << "Minimum energy: " << energy << std::endl;
      } else {
        std::cout << "Gradient descent aborted after " << this->maxIterations_ << " iterations." << std::endl;
        std::cout << "Maximum gradient: " << maxGrad << std::endl;
        std::cout << "Energy reached: " << energy << std::endl;
      }
    }

    bool solnFound_ = false;

  private:
    int maxIterations_;
    double tol_;
    double stepSize_ = 0.01;
    double maxStepSize_ = 1;

    virtual double energyDensity(LATfield2::Field<FieldType> &field, LATfield2::Site &site) const = 0;
    virtual double gradient(LATfield2::Field<FieldType> &field, LATfield2::Site &site, int cpt) const = 0;
  };
}

#endif