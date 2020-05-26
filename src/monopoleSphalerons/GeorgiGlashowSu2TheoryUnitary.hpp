#ifndef GEORGIGLASHOWSU2THEORYUNITARY_HPP
#define GEORGIGLASHOWSU2THEORYUNITARY_HPP

#include "Theory.hpp"
#include "Su2Tools.hpp"

namespace monsta {
  class GeorgiGlashowSu2Theory: public Theory {
  public:
    GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling);
    GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling, std::vector<int> boundaryConditions);
    GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling, std::vector<int> boundaryConditions, bool tHooftLine);
    double getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    double getSymmetricEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    std::complex<double> getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const;
    monsta::Matrix getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;
    double getVev() const { return vev_; }
    double getGaugeCoupling() const { return gaugeCoupling_; }
    double getSelfCoupling() const { return selfCoupling_; }
    double getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    int getMonopoleNumber(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    void postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;
    void applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const;
    double getHiggsMagnitude(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;

  private:
    bool tHooftLineCheck(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    void applyDirichletBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const;
    void applyCPeriodicBoundaryConditions(LATfield2::Field< std::complex<double> > &field, int dir) const;
    void applyTwistedBoundaryConditions(LATfield2::Field< std::complex<double> > &field, int dir) const;
    void applyTwistedBoundaryConditions2(LATfield2::Field< std::complex<double> > &field, int dir) const;
    void applyTwistedBoundaryConditions3(LATfield2::Field< std::complex<double> > &field, int dir) const;
    monsta::Matrix getPlaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const;
    monsta::Matrix getStaple(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, bool isUp) const;
    std::complex<double> getU1Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const;

    int numFieldMatrices_ = 4;
    int numRows_ = 2;
    int numCols_ = 2;
    int matSize_ = 2;
    double gaugeCoupling_;
    double vev_;
    double selfCoupling_;
    std::vector<int> boundaryConditions_ = {0, 0, 0};
    bool tHooftLine_ = false;
    std::complex<double> fluxQuanta_ = 1;

    bool monopolesPresent_ = false;
    std::vector<int> monopolePos1_;
    std::vector<int> monopolePos2_;
  };

  GeorgiGlashowSu2Theory::GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling)
  : Theory(numFieldMatrices_, matSize_, matSize_), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling)
  {}

  GeorgiGlashowSu2Theory::GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling, std::vector<int> boundaryConditions)
  : Theory(numFieldMatrices_, matSize_, matSize_), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling),
    boundaryConditions_(boundaryConditions)
  {}

  GeorgiGlashowSu2Theory::GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling, std::vector<int> boundaryConditions, bool tHooftLine)
  : Theory(numFieldMatrices_, matSize_, matSize_), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling),
    boundaryConditions_(boundaryConditions), tHooftLine_(tHooftLine)
  {}

  double GeorgiGlashowSu2Theory::getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    double E = 0;

    // Wilson action
    E += 2.0/pow(gaugeCoupling_,2)*(2 - real(trace(getPlaquette(field, site, 1, 2))));
    E += 2.0/pow(gaugeCoupling_,2)*(2 - real(trace(getPlaquette(field, site, 2, 0))));
    E += 2.0/pow(gaugeCoupling_,2)*(2 - real(trace(getPlaquette(field, site, 0, 1))));

    // Covariant Derivative
    for (int ii = 0; ii < 3; ii++)
    {
      LATfield2::Site shiftedSite = site + ii;
      Matrix covDeriv = field(shiftedSite, 3, 0, 0)*Matrix(field, site, ii)*pauli3*conjugateTranspose(Matrix(field, site, ii)) - field(site, 3, 0, 0)*pauli3;
      E += real(trace(covDeriv*covDeriv));
    }

    // Higgs Potential
    E += real(selfCoupling_*pow(2.0*pow(field(site, 3, 0, 0),2) - pow(vev_, 2),2));

    return E;
  }

  double GeorgiGlashowSu2Theory::getSymmetricEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    double E = 0;
    LATfield2::Site tempSite(site);

    // Wilson action
    for (int ii = 0; ii < 3; ii++)
    {
      for (int jj = 0; jj < 3; jj++)
      {
        if (jj >= ii) { continue; }
        E += 0.5/pow(gaugeCoupling_,2)*(2 - real(trace(getPlaquette(field, site, ii, jj))));
        tempSite = site - ii;
        E += 0.5/pow(gaugeCoupling_,2)*(2 - real(trace(getPlaquette(field, tempSite, ii, jj))));
        tempSite = site - jj;
        E += 0.5/pow(gaugeCoupling_,2)*(2 - real(trace(getPlaquette(field, tempSite, ii, jj))));
        tempSite = site - ii - jj;
        E += 0.5/pow(gaugeCoupling_,2)*(2 - real(trace(getPlaquette(field, tempSite, ii, jj))));
      }
    }

    // Covariant Derivative
    for (int ii = 0; ii < 3; ii++)
    {
      tempSite = site + ii;
      Matrix covDeriv = field(tempSite, 3, 0, 0)*Matrix(field, site, ii)*pauli3*conjugateTranspose(Matrix(field, site, ii)) - field(site, 3, 0, 0)*pauli3;
      E += 0.5*real(trace(covDeriv*covDeriv));
      tempSite = site - ii;
      covDeriv = field(site, 3, 0, 0)*Matrix(field, tempSite, ii)*pauli3*conjugateTranspose(Matrix(field, tempSite, ii)) - field(tempSite, 3, 0, 0)*pauli3;
      E += 0.5*real(trace(covDeriv*covDeriv));
    }

    // Higgs Potential
    E += real(selfCoupling_*pow(2.0*pow(field(site, 3, 0, 0),2) - pow(vev_, 2),2));

    return E;
  }

  monsta::Matrix GeorgiGlashowSu2Theory::getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
  {
    monsta::Matrix grad(2);
    if (matIdx < 3)
    {
      int dir1 = (matIdx + 1) % 3;
      int dir2 = (matIdx + 2) % 3;
      int derivIdx = site.index();

      // Derivative of Wilson action
      Matrix plaquetteDerivMat = getStaple(field, site, matIdx, dir1, true);
      plaquetteDerivMat = plaquetteDerivMat + getStaple(field, site, matIdx, dir1, false);
      plaquetteDerivMat = plaquetteDerivMat + getStaple(field, site, matIdx, dir2, true);
      plaquetteDerivMat = plaquetteDerivMat + getStaple(field, site, matIdx, dir2, false);

      grad = grad - 2.0/pow(gaugeCoupling_,2)*plaquetteDerivMat;
      // grad = grad - 2.0/pow(gaugeCoupling_,2)*getTHooftDeriv(field, site, matIdx, fluxQuanta_);

      // Derivative of kinetic term
      LATfield2::Site tempSite(site);
      tempSite = site + matIdx;
      Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
      Matrix scalarMatShifted = field(tempSite, 3, 0, 0)*pauli3;
      Matrix gaugeMat(field, site, matIdx);
      Matrix kineticDerivMat = gaugeMat*scalarMatShifted*conjugateTranspose(gaugeMat)*gaugeMat*scalarMatShifted;
      kineticDerivMat = kineticDerivMat + gaugeMat*conjugateTranspose(scalarMatShifted)*conjugateTranspose(gaugeMat)*gaugeMat*conjugateTranspose(scalarMatShifted);
      kineticDerivMat = kineticDerivMat - conjugateTranspose(scalarMat)*gaugeMat*conjugateTranspose(scalarMatShifted);
      kineticDerivMat = kineticDerivMat - scalarMat*gaugeMat*scalarMatShifted;
      grad = grad + 2.*kineticDerivMat;

      grad = grad - 0.5*trace(grad*conjugateTranspose(Matrix(field, site, matIdx)))*Matrix(field, site, matIdx);

    } else {
      for (int dir = 0; dir < 3; dir++)
      {
        // Deriviative of kinetic term
        LATfield2::Site tempSite = site + dir;
        Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
        Matrix scalarMatShiftedFwd = field(tempSite, 3, 0, 0)*pauli3;
        Matrix gaugeMat(field, site, dir);

        tempSite = site - dir;
        Matrix gaugeMatShiftedBwd(field, tempSite, dir);
        Matrix scalarMatShiftedBwd = field(tempSite, 3, 0, 0)*pauli3;

        Matrix kineticDerivMat = 2*conjugateTranspose(gaugeMatShiftedBwd)*gaugeMatShiftedBwd*conjugateTranspose(scalarMat)*conjugateTranspose(gaugeMatShiftedBwd)*gaugeMatShiftedBwd;
        kineticDerivMat = kineticDerivMat - 2*gaugeMat*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(gaugeMat);
        kineticDerivMat = kineticDerivMat - 2*conjugateTranspose(gaugeMatShiftedBwd)*conjugateTranspose(scalarMatShiftedBwd)*gaugeMatShiftedBwd;
        kineticDerivMat = kineticDerivMat + 2*conjugateTranspose(scalarMat);

        grad(0,0) = grad(0,0) + 2.*kineticDerivMat(0,0);
      }

      // Derivative of Higgs Potential
      grad(0,0) = grad(0,0) + 8.0*selfCoupling_*field(site, 3, 0, 0)*(2.0*pow(field(site, 3, 0, 0),2) - pow(vev_, 2));
    }

    return grad;
  }

  std::complex<double> GeorgiGlashowSu2Theory::getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const
  {
    monsta::Matrix gradMat = getLocalGradient(field, site, matIdx);
    return gradMat(rowIdx, colIdx);
    return 0;
  }

  void GeorgiGlashowSu2Theory::postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
  {
    if (matIdx < 3) // Project to SU(2)
    {
      std::complex<double> determinant = field(site, matIdx, 0, 0)*field(site, matIdx, 1, 1) - field(site, matIdx, 0, 1)*field(site, matIdx, 1, 0);
      for (int rowIdx = 0; rowIdx < numRows_; rowIdx++)
      {
        for (int colIdx = 0; colIdx < numCols_; colIdx++)
        {
          field(site, matIdx, rowIdx, colIdx) = field(site, matIdx, rowIdx, colIdx) / sqrt(determinant);
        }
        field(site, matIdx, 0, 0) = conj(field(site, matIdx, 1, 1));
        field(site, matIdx, 1, 0) = -1.0*conj(field(site, matIdx, 0, 1));
      }
    }
    else 
    {
      field(site, matIdx, 0, 0) = real(field(site, matIdx, 0, 0));
    }
  }

  void GeorgiGlashowSu2Theory::applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const
  {
    field.updateHalo();
    for (int dir = 0; dir < 3; dir++)
    {
      applyCPeriodicBoundaryConditions(field, dir);
    }
  }

  void GeorgiGlashowSu2Theory::applyDirichletBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const
  {
    // This is inefficient but works for now
    LATfield2::Field< std::complex<double> > tempField(field.lattice(), numFieldMatrices_, numRows_, numCols_, 0);

    LATfield2::Site site(field.lattice());
    int xMax = field.lattice().size(0) - 1;
    int yMax = field.lattice().size(1) - 1;
    int zMax = field.lattice().size(2) - 1;

    for (site.haloFirst(); site.haloTest(); site.haloNext())
    {
      int xCoord = site.coord(0);
      int yCoord = site.coord(1);
      int zCoord = site.coord(2);
      if (xCoord < 0 || yCoord < 0 || zCoord < 0 || xCoord > xMax || yCoord > yMax || zCoord > zMax) // Sites on boundary
      {
        for (int matIdx = 0; matIdx < numFieldMatrices_; matIdx++)
        {
          tempField(site, matIdx, 0, 0) = field(site, matIdx, 0, 0);
          tempField(site, matIdx, 0, 1) = field(site, matIdx, 0, 1);
          tempField(site, matIdx, 1, 0) = field(site, matIdx, 1, 0);
          tempField(site, matIdx, 1, 1) = field(site, matIdx, 1, 1);
        }
      }
    }
    // field.updateHalo();
    for (site.haloFirst(); site.haloTest(); site.haloNext())
    {
      int xCoord = site.coord(0);
      int yCoord = site.coord(1);
      int zCoord = site.coord(2);
      if (xCoord < 0 || yCoord < 0 || zCoord < 0 || xCoord > xMax || yCoord > yMax || zCoord > zMax) // Sites on boundary
      {
        for (int matIdx = 0; matIdx < numFieldMatrices_; matIdx++)
        {
          field(site, matIdx, 0, 0) = tempField(site, matIdx, 0, 0);
          field(site, matIdx, 0, 1) = tempField(site, matIdx, 0, 1);
          field(site, matIdx, 1, 0) = tempField(site, matIdx, 1, 0);
          field(site, matIdx, 1, 1) = tempField(site, matIdx, 1, 1);
        }
      }
    }
  }

  void GeorgiGlashowSu2Theory::applyCPeriodicBoundaryConditions(LATfield2::Field< std::complex<double> > &field, int dir) const
  {
    LATfield2::Site site(field.lattice());

    for (site.haloFirst(); site.haloTest(); site.haloNext())
    {
      int coord = site.coord(dir);
      int maxCoord = field.lattice().size(dir) - 1;
      if (coord < 0 || coord > maxCoord)
      {
        for (int matIdx = 0; matIdx < numFieldMatrices_; matIdx++) // Apply C-periodic boundary conditions
        {
          int pauliMatNum = boundaryConditions_[dir];
          Matrix boundaryMat(field, site, matIdx);
          if (matIdx < 3)
          {
            switch(pauliMatNum)
            {
              case 0:
                break;
              case 1:
                boundaryMat = pauli1*boundaryMat*pauli1;
                break;
              case 2:
                boundaryMat = pauli2*boundaryMat*pauli2;
                break;
              case 3:
                boundaryMat = pauli3*boundaryMat*pauli3;
                break;
            }
            if (tHooftLine_ && dir == 2 && site.coord(1) == 0)
            {
              if (matIdx == 1)
              {
                boundaryMat = -1.0*boundaryMat;
              }
            }
          }
          else
          {
            if (pauliMatNum == 0) { break; }
            if (pauliMatNum == 3)
            {
              boundaryMat(0,0) = -1.0*boundaryMat(0,0);
            }
            else
            {
              boundaryMat(0,0) = boundaryMat(0,0);
            }
          }

          field(site, matIdx, 0, 0) = boundaryMat(0, 0);
          field(site, matIdx, 0, 1) = boundaryMat(0, 1);
          field(site, matIdx, 1, 0) = boundaryMat(1, 0);
          field(site, matIdx, 1, 1) = boundaryMat(1, 1);
        }
      }
    }
  }

  monsta::Matrix GeorgiGlashowSu2Theory::getPlaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix plaquette(field, site, dir1);
    tempSite = tempSite+dir1;
    plaquette = plaquette*Matrix(field, tempSite, dir2);
    tempSite = tempSite-dir1+dir2;
    plaquette = plaquette*conjugateTranspose(Matrix(field, tempSite, dir1));
    tempSite = tempSite-dir2;
    plaquette = plaquette*conjugateTranspose(Matrix(field, tempSite, dir2));

    return plaquette;
  }

  monsta::Matrix GeorgiGlashowSu2Theory::getStaple(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, bool isUp) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix staple(field, site, dir2);

    if (isUp)
    {
      tempSite = tempSite + dir2;
      staple = staple*Matrix(field, tempSite, dir1);
      tempSite = tempSite-dir2+dir1;
      staple = staple*conjugateTranspose(Matrix(field, tempSite, dir2));
    }
    else
    {
      tempSite = tempSite - dir2;
      staple = conjugateTranspose(Matrix(field, tempSite, dir2));
      staple = staple*Matrix(field, tempSite, dir1);
      tempSite = tempSite + dir1;
      staple = staple*Matrix(field, tempSite, dir2);
    }

    // cout << "jonnyMorris" << endl;

    return staple;

  }

  std::complex<double> GeorgiGlashowSu2Theory::getU1Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
  {
    if (real(field(site, 3, 0, 0)) > 0)
    {
      return field(site, cpt, 0, 0);
    } else
    {
      return field(site, cpt, 1, 1);
    }
  }

  double GeorgiGlashowSu2Theory::getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
  {
    int dir1 = (cpt + 1) % 3;
    int dir2 = (cpt + 2) % 3;


    LATfield2::Site tempSite(site);
    std::complex<double> u1Plaquette = field(site, dir1, 0, 0);
    tempSite = tempSite+dir1;
    u1Plaquette = u1Plaquette*field(tempSite, dir2, 0, 0);
    tempSite = tempSite-dir1+dir2;
    u1Plaquette = u1Plaquette*conj(field(tempSite, dir1, 0, 0));
    tempSite = tempSite-dir2;
    u1Plaquette = u1Plaquette*conj(field(tempSite, dir2, 0, 0));

    double magneticField = 2./gaugeCoupling_*arg(u1Plaquette);

    return (magneticField);
  }

  double GeorgiGlashowSu2Theory::getHiggsMagnitude(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    return real(field(site, 3, 0, 0));
  }

  int GeorgiGlashowSu2Theory::getMonopoleNumber(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    LATfield2::Site shiftedSite(field.lattice());
    double divB = 0;
    double pi = 4*atan(1);

    for (int dir = 0; dir < 3; dir++)
    {
      shiftedSite = site+dir;
      divB += getMagneticField(field, shiftedSite, dir);
      divB -= getMagneticField(field, site, dir);
    }

    return round(divB*gaugeCoupling_ / (4*pi));
  }

}

#endif