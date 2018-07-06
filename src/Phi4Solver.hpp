#ifndef PHI4SOLVER_HPP
#define PHI4SOLVER_HPP

#include "LATfield2.hpp"
#include "GradDescentSolver.hpp"
#include <cmath>

namespace monsta {
  class Phi4Solver: public GradDescentSolver<double> {
  public:
    Phi4Solver(LATfield2::Lattice &lattice, int numCpts, double tol,
      int maxIterations, double vev, double selfCoupling) 
      : GradDescentSolver(lattice, numCpts, tol, maxIterations)
    {
      this->vev_ = vev;
      this->selfCoupling_ = selfCoupling;
    }

  private:
    double vev_;
    double selfCoupling_;

    double getLocalEnergyDensity(LATfield2::Field<double> &field, LATfield2::Site &site) const
    {
      int numCpts = field.components();
      double absFieldSqr = 0;
      for (int iCpt = 0; iCpt < numCpts; iCpt++) {
        absFieldSqr += pow(field(site,iCpt), 2);
      }
      double E = this->selfCoupling_ * pow(absFieldSqr - pow(this->vev_, 2), 2);
      return E;
    }

    double getLocalGradient(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt) const
    {
      int numCpts = field.components();
      double absFieldSqr = 0;
      for (int iCpt = 0; iCpt < numCpts; iCpt++) {
        absFieldSqr += pow(field(site,iCpt), 2);
      }
      double grad = 4 * this->selfCoupling_ * field(site,cpt)*(absFieldSqr - pow(this->vev_, 2));
      return grad;
    }
  };
}

#endif