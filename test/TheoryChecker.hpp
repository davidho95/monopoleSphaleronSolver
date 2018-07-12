#ifndef THEORYCHECKER_HPP
#define THEORYCHECKER_HPP

#include <cstdio>
#include <cmath>
#include "LATfield2.hpp"
#include "../src/System.hpp"


namespace monsta {
  class TheoryChecker {
  public:
    TheoryChecker(double tol, bool verbose = false)
    : tol_(tol), verbose_(verbose) {}

    bool checkGradients(monsta::System &system) const
    {
      LATfield2::Site site(system.getLattice());
      bool checkPassed;
      int numCpts = system.getNumSpatialCpts();

      double maxErr = 0;
      int count = 0;
      for (site.first(); site.test(); site.next())
      {
        for (int iCpt = 0; iCpt < numCpts; iCpt++)
        {
          double suppliedGradVal = system.getGradientElement(site, iCpt);
          double fdGradVal = getFdGradient(system, site, iCpt);
          double relError;
          if (abs(suppliedGradVal) > tol_)
          {
            relError = (suppliedGradVal - fdGradVal) / fdGradVal;
          }
          else
          {
            relError = suppliedGradVal - fdGradVal;
          }
          if (abs(relError) > abs(maxErr))
          {
            maxErr = relError;
          }

          if (verbose_ && abs(relError) > tol_)
          {
            std::cout << "(" << iCpt << ";"
                             << site.coord(0) << ","
                             << site.coord(1) << ","
                             << site.coord(2) << ") : ";
            std::cout << suppliedGradVal << " " << fdGradVal << endl;
          }
        }
      }
      checkPassed = (abs(maxErr) < tol_);
      if (checkPassed)
      {
        std::cout << "Gradient check passed: maximum deviation from supplied gradient: " << maxErr << std::endl;
      }
      else
      {
        std::cout << "Gradient check failed: maximum deviation from supplied gradient: " << maxErr << std::endl;
      }
      return checkPassed;
    }

  private:
    double tol_;
    bool verbose_;
    double fdStep_ = 1e-5;

    double getFdGradient(monsta::System &system,
      LATfield2::Site &site, int iCpt) const
    {
      double originalFieldVal = system.getFieldElement(site, iCpt);

      //Compute gradient numerically
      system.setFieldElement(site, iCpt, originalFieldVal + fdStep_);
      system.updateFieldHalo();
      double EPlus = system.getEnergy();
      system.setFieldElement(site, iCpt, originalFieldVal - fdStep_);
      system.updateFieldHalo();
      double EMinus = system.getEnergy();

      // Restore field to original
      system.setFieldElement(site, iCpt, originalFieldVal);

      return (EPlus - EMinus) / (2*fdStep_);
    }
  };
}

#endif