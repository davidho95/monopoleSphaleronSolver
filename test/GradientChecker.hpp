#ifndef GRADIENTCHECKER_HPP
#define GRADIENTCHECKER_HPP

#include <cstdio>
#include <cmath>
#include "LATfield2.hpp"
#include "../src/GradDescentSolver.hpp"


namespace monsta {

  template<class FieldType>
  class GradientChecker {
  public:
    GradientChecker(double tol)
    {
      this->tol_ = tol;
    }

    bool checkGradients(const monsta::GradDescentSolver<FieldType> &solver,
      LATfield2::Field<FieldType> &testField) const
    {
      bool checkPassed;
      LATfield2::Lattice &lattice = testField.lattice();
      LATfield2::Site site(lattice);
      int numCpts = testField.components();

      LATfield2::Field<FieldType> suppliedGradients = solver.getGradient(testField);

      double maxErr;
      for (site.first(); site.test(); site.next())
      {
        for (int iCpt = 0; iCpt < numCpts; iCpt++)
        {          
          double suppliedGradVal = suppliedGradients(site, iCpt);
          double fdGradVal = getFdGradient(solver, testField, site, iCpt);
          double relError;
          if (abs(suppliedGradVal - fdGradVal) > this->tol_) {
            relError = (suppliedGradVal - fdGradVal) / fdGradVal;
          } else {
            relError = suppliedGradVal - fdGradVal;
          }
          relError = (suppliedGradVal - fdGradVal) / fdGradVal;
          if (abs(relError) > abs(maxErr)) {
            maxErr = relError;
          }

          checkPassed = (abs(relError) > this->tol_);
        }
      }
      if (checkPassed)
      {
        std::cout << "Gradient check failed: maximum deviation from supplied gradient: " << maxErr << std::endl;
      }
      else
      {
        std::cout << "Gradient check Passed: maximum deviation from supplied gradient: " << maxErr << std::endl;
      }
      return checkPassed;
    }

  private:
    double tol_;
    int testLatticeSize = 5;
    double fdStep_ = 1e-5;

    double getFdGradient(const monsta::GradDescentSolver<FieldType> &solver,
      LATfield2::Field<FieldType> &testField, LATfield2::Site &site, int iCpt) const
    {
      double originalFieldVal = testField(site, iCpt);

      //Compute gradient numerically
      testField(site, iCpt) += this->fdStep_;
      testField.updateHalo();
      double EPlus = solver.getEnergy(testField);
      testField(site, iCpt) -= 2*this->fdStep_;
      testField.updateHalo();
      double EMinus = solver.getEnergy(testField);

      // Restore field to original
      testField(site, iCpt) = originalFieldVal;
      testField.updateHalo();

      return (EPlus - EMinus) / (2*this->fdStep_);
    }
  };
}

#endif