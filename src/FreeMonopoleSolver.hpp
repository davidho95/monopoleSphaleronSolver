#ifndef FREEMONOPOLESOLVER_HPP
#define FREEMONOPOLESOLVER_HPP

#include "LATfield2.hpp"
#include "GradDescentSolver.hpp"
#include <cmath>

namespace monsta
{
  class FreeMonopoleSolver: public GradDescentSolver<double>
  {
  public:
    FreeMonopoleSolver(LATfield2::Lattice &lattice, int numCpts, double tol,
      int maxIterations, double magneticCharge, int southPos, int northPos) 
      : GradDescentSolver(lattice, numCpts, tol, maxIterations),
        magneticCharge_(magneticCharge), southPos_(southPos), northPos_(northPos)
      {
        int *latSize = lattice.size();
        yMid_ = latSize[1]/2;
        zMid_ = latSize[2]/2;
      }

  private:
    double magneticCharge_;
    int southPos_;
    int northPos_;
    int yMid_;
    int zMid_;

    double getLocalEnergyDensity(LATfield2::Field<double> &field, LATfield2::Site &site) const
    {
      int xCoord = site.coord(0);
      int yCoord = site.coord(1);
      int zCoord = site.coord(2);
      double Bx = 0;
      if (yCoord == yMid_ && zCoord == zMid_)
      {
        if (xCoord >= southPos_ && xCoord < northPos_) { Bx = magneticCharge_; }
      }

      double E = 0;
        E += pow(field(site+1, 2) - field(site, 2) - field(site+2, 1) + field(site,1), 2);
        E += pow(field(site+2, 0) - field(site, 0) - field(site+0, 2) + field(site,2), 2);
        E += pow(field(site+0, 1) - field(site, 1) - field(site+1, 0) + field(site,0), 2);
      return 0.5*E;
    }

    double getLocalGradient(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt) const
    {
      int xCoord = site.coord(0);
      int yCoord = site.coord(1);
      int zCoord = site.coord(2);
      double Bx = 0;
      double BxShiftY = 0;
      double BxShiftZ = 0;
      if (xCoord >= southPos_ && xCoord < northPos_)
      {
        if (yCoord == yMid_ && zCoord == zMid_) { Bx = magneticCharge_; }
        if (yCoord == yMid_ + 1 && zCoord == zMid_) { BxShiftY = magneticCharge_; }
        if (yCoord == yMid_ && zCoord == zMid_ + 1) { BxShiftZ = magneticCharge_; }
      }

      double grad;
      switch (cpt)
      {
        case 0:
          grad = -field(site-1, 0) - field(site-2, 0) + 4*field(site, 0) - field(site+2, 0) - 
            field(site+1, 0) + field(site-1, 1) - field(site, 1) - field(site+0-1, 1) + 
            field(site+0, 1) + field(site-2, 2) - field(site, 2) - field(site+0-2, 2) + 
            field(site+0, 2);
          break;
        case 1:
          grad = -BxShiftZ + Bx + field(site-0, 0) - field(site-0+1, 0) - field(site, 0) + field(site+1, 0) - 
            field(site-0, 1) - field(site-2, 1) + 4*field(site, 1) - field(site+2, 1) - 
            field(site+0, 1) + field(site-2, 2) - field(site, 2) - field(site+1-2, 2) + 
            field(site+1, 2);
          break;
        case 2:
          grad = BxShiftY - Bx + field(site-0, 0) - field(site-0+2, 0) - field(site, 0) + field(site+2, 0) + 
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