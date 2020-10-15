#ifndef GEORGIGLASHOWSU2THEORY4D_HPP
#define GEORGIGLASHOWSU2THEORY4D_HPP

#include "../Theory.hpp"
#include "../Su2Tools.hpp"
#include "../Matrix.hpp"

namespace monsta {
  class GeorgiGlashowSu2Theory4d: public Theory {
  public:
    GeorgiGlashowSu2Theory4d(double gaugeCoupling, double vev, double selfCoupling);
    GeorgiGlashowSu2Theory4d(double gaugeCoupling, double vev, double selfCoupling, int fluxQuanta);
    double getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    double getAsymmetricEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
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
    Matrix getHiggsField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    Matrix getSu2Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const;
    double getRFromSite(LATfield2::Site &site) const;
    void setSu2Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt, Matrix su2Mat) const;
    void setHiggsMagnitude(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, double higgsVal) const;

  private:
    bool tHooftLineCheck(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    monsta::Matrix getPlaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const;
    monsta::Matrix getStaple(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, bool isUp) const;
    std::complex<double> getU1Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const;
    double getSiteJacobian(LATfield2::Site &site) const;
    double getLinkJacobian(LATfield2::Site &site, int dir) const;
    double getPlaquetteJacobian(LATfield2::Site &site, int dir1, int dir2) const;
    void symmetrise(LATfield2::Field< std::complex<double> > &field) const;



    int numFieldMatrices_ = 4;
    int numRows_ = 2;
    int numCols_ = 2;
    int matSize_ = 2;
    double gaugeCoupling_;
    double vev_;
    double selfCoupling_;
    std::vector<int> boundaryConditions_ = {0, 0, 0};
    bool tHooftLine_ = false;
    int fluxQuanta_ = 0;

    bool monopolesPresent_ = false;
    std::vector<int> monopolePos1_;
    std::vector<int> monopolePos2_;
  };

  GeorgiGlashowSu2Theory4d::GeorgiGlashowSu2Theory4d(double gaugeCoupling, double vev, double selfCoupling)
  : Theory(numFieldMatrices_, matSize_, matSize_), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling)
  {}

  GeorgiGlashowSu2Theory4d::GeorgiGlashowSu2Theory4d(double gaugeCoupling, double vev, double selfCoupling, int fluxQuanta)
  : Theory(numFieldMatrices_, matSize_, matSize_), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling), fluxQuanta_(fluxQuanta)
  {}

  double GeorgiGlashowSu2Theory4d::getRFromSite(LATfield2::Site &site) const
  {
    int xSize = site.lattice().size(0);
    if (site.coord(0) == -1)
    {
      return xSize - xSize / 2 - 0.5;
    }
    if (site.coord(0) == xSize)
    {
      return -xSize / 2 + 0.5;
    }
    return site.coord(0) - xSize / 2 + 0.5;
  }

  double GeorgiGlashowSu2Theory4d::getSiteJacobian(LATfield2::Site &site) const
  {
    return abs(getRFromSite(site));
  }

  double GeorgiGlashowSu2Theory4d::getLinkJacobian(LATfield2::Site &site, int dir) const
  {
    if (dir == 0)
    {
      return abs(getRFromSite(site) + 0.5);
    }
    else
    {
      return abs(getRFromSite(site));
    }
  }

  double GeorgiGlashowSu2Theory4d::getPlaquetteJacobian(LATfield2::Site &site, int dir1, int dir2) const
  {
    if (dir1 == dir2) { return 0; }
    if (dir1 == 0 || dir2 == 0)
    {
      return abs(getRFromSite(site) + 0.5);
    }
    else
    {
      return abs(getRFromSite(site));
    }
  }  

  double GeorgiGlashowSu2Theory4d::getAsymmetricEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    double E = 0;

    // Wilson action
    for (int ii = 0; ii < 3; ii++)
    {
      for (int jj = 0; jj < 3; jj++)
      {
        if (ii >= jj) { continue; }
        double plaquetteJacobian = getPlaquetteJacobian(site, ii, jj);
        E += plaquetteJacobian*2.0/pow(gaugeCoupling_,2)*(2 - real(trace(getPlaquette(field, site, ii, jj))));
      }
    }

    // // Covariant Derivative
    for (int ii = 0; ii < 3; ii++)
    {
      double linkJacobian = getLinkJacobian(site, ii);
      LATfield2::Site shiftedSite = site + ii;
      Matrix covDeriv = field(shiftedSite, 3, 0, 0)*Matrix(field, site, ii)*pauli3*conjugateTranspose(Matrix(field, site, ii)) - field(site, 3, 0, 0)*pauli3;
      E += linkJacobian*real(trace(covDeriv*covDeriv));
    }

    // Higgs Potential
    double siteJacobian = getSiteJacobian(site);
    E += siteJacobian*real(selfCoupling_*pow(2.0*pow(field(site, 3, 0, 0),2) - pow(vev_, 2),2));

    return pi*E;
  }

  double GeorgiGlashowSu2Theory4d::getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    double E = 0;
    double jacobian = getSiteJacobian(site);
    LATfield2::Site tempSite(site);
    // Wilson action
    for (int ii = 0; ii < 3; ii++)
    {
      for (int jj = 0; jj < 3; jj++)
      {
        if (jj >= ii) { continue; }
        E += 0.5/pow(gaugeCoupling_,2)*jacobian*(2 - real(trace(getPlaquette(field, site, ii, jj))));
        tempSite = site - ii;
        E += 0.5/pow(gaugeCoupling_,2)*jacobian*(2 - real(trace(getPlaquette(field, tempSite, ii, jj))));
        tempSite = site - jj;
        E += 0.5/pow(gaugeCoupling_,2)*jacobian*(2 - real(trace(getPlaquette(field, tempSite, ii, jj))));
        tempSite = site - ii - jj;
        E += 0.5/pow(gaugeCoupling_,2)*jacobian*(2 - real(trace(getPlaquette(field, tempSite, ii, jj))));
      }
    }

    // Covariant Derivative
    for (int ii = 0; ii < 3; ii++)
    {
      tempSite = site + ii;
      Matrix covDeriv = field(tempSite, 3, 0, 0)*Matrix(field, site, ii)*pauli3*conjugateTranspose(Matrix(field, site, ii)) - field(site, 3, 0, 0)*pauli3;
      E += 0.5*jacobian*real(trace(covDeriv*covDeriv));
      tempSite = site - ii;
      covDeriv = field(site, 3, 0, 0)*Matrix(field, tempSite, ii)*pauli3*conjugateTranspose(Matrix(field, tempSite, ii)) - field(tempSite, 3, 0, 0)*pauli3;
      E += 0.5*jacobian*real(trace(covDeriv*covDeriv));
    }

    // Higgs Potential
    E += jacobian*real(selfCoupling_*pow(2.0*pow(field(site, 3, 0, 0),2) - pow(vev_, 2),2));

    return pi*E;
  }


  monsta::Matrix GeorgiGlashowSu2Theory4d::getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
  {
    monsta::Matrix grad(2);


    if (matIdx < 3)
    {

      int dir1 = (matIdx + 1) % 3;
      int dir2 = (matIdx + 2) % 3;
      // Derivative of Wilson action
      Matrix plaquetteDerivMat(2);
      for (int dir = 0; dir < 3; dir++)
      {
        if (dir == matIdx) { continue; }
        LATfield2::Site siteShiftedDown = site - dir; 
        plaquetteDerivMat = plaquetteDerivMat + getPlaquetteJacobian(site, matIdx, dir)*getStaple(field, site, matIdx, dir, true);
        plaquetteDerivMat = plaquetteDerivMat + getPlaquetteJacobian(siteShiftedDown, matIdx, dir)*getStaple(field, site, matIdx, dir, false);
      }

      grad = grad - 2.0/pow(gaugeCoupling_,2)*plaquetteDerivMat;

      // Derivative of kinetic term
      LATfield2::Site siteShiftedFwd = site + matIdx;
      Matrix scalarMat = getHiggsField(field, site);
      Matrix scalarMatShiftedFwd = getHiggsField(field, siteShiftedFwd);
      Matrix gaugeMat = getSu2Link(field, site, matIdx);
      Matrix kineticDerivMat = gaugeMat*scalarMatShiftedFwd*scalarMatShiftedFwd;

      double linkJacobian = getLinkJacobian(site, matIdx);
      double linkJacobianShiftedFwd = getLinkJacobian(siteShiftedFwd, matIdx);

      kineticDerivMat = kineticDerivMat + gaugeMat*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(scalarMatShiftedFwd);
      kineticDerivMat = kineticDerivMat - conjugateTranspose(scalarMat)*gaugeMat*conjugateTranspose(scalarMatShiftedFwd);
      kineticDerivMat = kineticDerivMat - scalarMat*gaugeMat*scalarMatShiftedFwd;
      grad = grad + 2*linkJacobian*kineticDerivMat;

      // Project onto SU(2) Lie group: COMMENT OUT IF TESTING GRADIENTS
      grad = grad - 0.5*trace(grad*conjugateTranspose(Matrix(field, site, matIdx)))*Matrix(field, site, matIdx);

    } else {
      for (int dir = 0; dir < 3; dir++)
      {
        // Deriviative of kinetic term
        LATfield2::Site siteShiftedFwd = site + dir;
        Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
        Matrix scalarMatShiftedFwd = field(siteShiftedFwd, 3, 0, 0)*pauli3;
        Matrix gaugeMat(field, site, dir);

        LATfield2::Site siteShiftedBwd = site - dir;
        Matrix gaugeMatShiftedBwd(field, siteShiftedBwd, dir);
        Matrix scalarMatShiftedBwd = field(siteShiftedBwd, 3, 0, 0)*pauli3;

        double linkJacobian = getLinkJacobian(site, dir);
        double linkJacobianShiftedBwd = getLinkJacobian(siteShiftedBwd, dir);

        Matrix kineticDerivMat = 2*linkJacobian*conjugateTranspose(scalarMat);
        kineticDerivMat = kineticDerivMat - 2*linkJacobian*gaugeMat*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(gaugeMat);
        kineticDerivMat = kineticDerivMat - 2*linkJacobianShiftedBwd*conjugateTranspose(gaugeMatShiftedBwd)*conjugateTranspose(scalarMatShiftedBwd)*gaugeMatShiftedBwd;
        kineticDerivMat = kineticDerivMat + 2*linkJacobianShiftedBwd*conjugateTranspose(scalarMat);

        grad(0,0) = grad(0,0) + 2.*kineticDerivMat(0,0);
      }

      // Derivative of Higgs Potential
      double siteJacobian = getSiteJacobian(site);
      grad(0,0) = grad(0,0) + 8.0*siteJacobian*selfCoupling_*field(site, 3, 0, 0)*(2.0*pow(field(site, 3, 0, 0),2) - pow(vev_, 2));
    }

    // Scaling: COMMENT OUT IF TESTING GRADIENTS
    double r = getRFromSite(site);
    if (abs(r) > 1e-10)
    {
      grad = grad/abs(r);
    }

    return pi*grad;
  }

  std::complex<double> GeorgiGlashowSu2Theory4d::getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const
  {
    monsta::Matrix gradMat = getLocalGradient(field, site, matIdx);
    return gradMat(rowIdx, colIdx);
    return 0;
  }

  void GeorgiGlashowSu2Theory4d::postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
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
      }
      field(site, matIdx, 0, 0) = conj(field(site, matIdx, 1, 1));
      field(site, matIdx, 1, 0) = -1.0*conj(field(site, matIdx, 0, 1));
    }
    else 
    {
      field(site, matIdx, 0, 0) = real(field(site, matIdx, 0, 0));
    }
  }

  void GeorgiGlashowSu2Theory4d::applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const
  {
    field.updateHalo();
    symmetrise(field);
  }

  void GeorgiGlashowSu2Theory4d::symmetrise(LATfield2::Field< std::complex<double> > &field) const
  {
    LATfield2::Site site(field.lattice());
    
    for (site.first(); site.test(); site.next())
    {
      double r = getRFromSite(site);
      if (abs(r + 0.5) > 1e-10) { continue; }
      LATfield2::Site siteShiftedFwd = site + 0;

      for (int matIdx = 0; matIdx < numFieldMatrices_; matIdx++)
      {
        if (matIdx < 3)
        {
          if (matIdx == 0) { continue; }
          Matrix su2Mat = getSu2Link(field, site, matIdx);
          su2Mat = pauli3*su2Mat*pauli3;
          setSu2Link(field, siteShiftedFwd, matIdx, su2Mat);
        }
        else
        {
          double higgsVal = getHiggsMagnitude(field, site);
          setHiggsMagnitude(field, siteShiftedFwd, higgsVal);
        }
      }
    }
  }

  Matrix GeorgiGlashowSu2Theory4d::getSu2Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
  {
    return Matrix(field, site, cpt);
  }

  monsta::Matrix GeorgiGlashowSu2Theory4d::getPlaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const
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

  monsta::Matrix GeorgiGlashowSu2Theory4d::getStaple(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, bool isUp) const
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

    return staple;

  }

  std::complex<double> GeorgiGlashowSu2Theory4d::getU1Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
  {
    if (real(field(site, 3, 0, 0)) > 0)
    {
      return field(site, cpt, 0, 0);
    } else
    {
      return field(site, cpt, 1, 1);
    }
  }

  double GeorgiGlashowSu2Theory4d::getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
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

  Matrix GeorgiGlashowSu2Theory4d::getHiggsField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    return real(field(site, 3, 0, 0))*pauli3;
  }

  double GeorgiGlashowSu2Theory4d::getHiggsMagnitude(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    return real(field(site, 3, 0, 0));
  }

  int GeorgiGlashowSu2Theory4d::getMonopoleNumber(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
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

  void GeorgiGlashowSu2Theory4d::setSu2Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, Matrix su2Mat) const
  {
    field(site, matIdx, 0, 0) = su2Mat(0, 0);
    field(site, matIdx, 0, 1) = su2Mat(0, 1);
    field(site, matIdx, 1, 0) = su2Mat(1, 0);
    field(site, matIdx, 1, 1) = su2Mat(1, 1);
  }

  void GeorgiGlashowSu2Theory4d::setHiggsMagnitude(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, double higgsVal) const
  {
    field(site, 3, 0, 0) = higgsVal;
  }

}

#endif
