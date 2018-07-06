#ifndef FREESPACEMAGNETICMAXWELLSOLVER_HPP
#define FREESPACEMAGNETICMAXWELLSOLVER_HPP

#include "LATfield2.hpp"
#include "GradDescentSolver.hpp"
#include <cmath>

namespace monsta {
  class FreeSpaceMagneticMaxwellSolver: public GradDescentSolver<double>
  {
  public:
    FreeSpaceMagneticMaxwellSolver(LATfield2::Lattice &lattice, int numCpts, double tol,
      int maxIterations) 
      : GradDescentSolver(lattice, numCpts, tol, maxIterations) {}

  private:
    double vev_;
    double selfCoupling_;

    double getLocalEnergyDensity(LATfield2::Field<double> &field, LATfield2::Site &site) const
    {
      double E = 0;
        E += pow(field(site+1, 2) - field(site, 2) - field(site+2, 1) + field(site,1), 2);
        E += pow(field(site+2, 0) - field(site, 0) - field(site+0, 2) + field(site,2), 2);
        E += pow(field(site+0, 1) - field(site, 1) - field(site+1, 0) + field(site,0), 2);
      // cout << E << endl;
      return 0.5*E;
    }

    double getLocalGradient(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt) const
    {
      double grad;
      switch (cpt) {
        case 0:
        grad = -field(site-1, 0) - field(site-2, 0) + 4*field(site, 0) - field(site+2, 0) - 
            field(site+1, 0) + field(site-1, 1) - field(site, 1) - field(site+0-1, 1) + 
            field(site+0, 1) + field(site-2, 2) - field(site, 2) - field(site+0-2, 2) + 
            field(site+0, 2);
          break;
        case 1:
          grad = field(site-0, 0) - field(site-0+1, 0) - field(site, 0) + field(site+1, 0) - 
            field(site-0, 1) - field(site-2, 1) + 4*field(site, 1) - field(site+2, 1) - 
            field(site+0, 1) + field(site-2, 2) - field(site, 2) - field(site+1-2, 2) + 
            field(site+1, 2);
          break;
        case 2:
          grad = field(site-0, 0) - field(site-0+2, 0) - field(site, 0) + field(site+2, 0) + 
            field(site-1, 1) - field(site-1+2, 1) - field(site, 1) + field(site+2, 1) - 
            field(site-0, 2) - field(site-1, 2) + 4*field(site, 2) - field(site+1, 2) - 
            field(site+0, 2);
            break;
      }
      return grad;
    }
  };
}

#endif