#ifndef GEORGIGLASHOWSU2THEORYAXISYMMETRIC
#define GEORGIGLASHOWSU2THEORYAXISYMMETRIC

#include "../Theory.hpp"
#include "../Su2Tools.hpp"

namespace monsta {
  class GeorgiGlashowSu2TheoryAxisymmetric: public Theory {
  public:
    GeorgiGlashowSu2TheoryAxisymmetric(double gaugeCoupling, double vev, double selfCoupling, int numThetaPoints);
    double getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    std::complex<double> getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const;
    monsta::Matrix getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;
    double getVev() const { return vev_; }
    double getGaugeCoupling() const { return gaugeCoupling_; }
    double getSelfCoupling() const { return selfCoupling_; }
    double getNumThetaPoints() const { return numThetaPoints_; }
    double getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    double getHiggsMagnitude(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    void applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const;
    void postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;
    monsta::Matrix getSu2Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const;
    std::complex<double> getU1Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const { return 0; };
    void setSu2Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt, monsta::Matrix su2val) const;
    void setHiggsMagnitude(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, std::complex<double> higgsVal) const;
    Matrix getPlaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const;
    std::complex<double> getU1Plaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const;


  private:
    int numFieldMatrices_ = 4;
    int numRows_ = 2;
    int numCols_ = 2;
    int matSize_ = 2;
    double gaugeCoupling_;
    double vev_;
    double selfCoupling_;
    bool tHooftLine_ = false;
    int numThetaPoints_;
    std::vector<Matrix> r0BoundaryVals;
    std::vector< std::vector<Matrix> > rMaxBoundaryVals;

    Matrix getCovDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    Matrix getStaple(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, bool isUp) const;
    double getRFromSite(LATfield2::Site &site) const;
    double getJacobian(LATfield2::Site &site) const;
    double getInverseMetricCpt(LATfield2::Site &site, int cpt) const;
    double getPlaquetteScaleFactor(int cpt1, int cpt2) const;
    double getKineticScaleFactor(int cpt) const;
    LATfield2::Site shiftSite(LATfield2::Site &site, int dir, int sign) const;
    double getPlaquetteAveragedMetricFactors(LATfield2::Site &site, int dir1, int dir2) const;
    double getFourAveragedInverseMetricFactor(LATfield2::Site &site, int dir1, int dir2) const;



  };

  GeorgiGlashowSu2TheoryAxisymmetric::GeorgiGlashowSu2TheoryAxisymmetric(double gaugeCoupling, double vev, double selfCoupling, int numThetaPoints)
  : Theory(numFieldMatrices_, matSize_, matSize_), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling), numThetaPoints_(numThetaPoints)
  {}

  LATfield2::Site GeorgiGlashowSu2TheoryAxisymmetric::shiftSite(LATfield2::Site &site, int dir, int sign) const
  {
    if (dir > 1) { return site; }
    if (sign > 0) { return site + dir; }
    return site - dir;
  }

  double GeorgiGlashowSu2TheoryAxisymmetric::getRFromSite(LATfield2::Site &site) const
  {
    // Currently set up so gradient checker will be satisfied with r-perioic boundary
    int rSize = site.lattice().size(1);
    // if (site.coord(1) == -1) { return 1e-5; }
    // return site.coord(1) + 1;

    if (site.coord(1) == -1) { return rSize; }
    if (site.coord(1) == rSize) { return 1; }
    return site.coord(1) + 1;
  }

  double GeorgiGlashowSu2TheoryAxisymmetric::getJacobian(LATfield2::Site &site) const
  {
    return getRFromSite(site);
  }

  double GeorgiGlashowSu2TheoryAxisymmetric::getInverseMetricCpt(LATfield2::Site &site, int cpt) const
  {
    if (cpt != 2) { return 1; }
    return pow(getRFromSite(site), -2);
  }

  double GeorgiGlashowSu2TheoryAxisymmetric::getPlaquetteScaleFactor(int cpt1, int cpt2) const
  {
    if (cpt1 == cpt2) { return 0; }
    if (cpt1 == 2 || cpt2 == 2 ) { return numThetaPoints_ / (2*pi); }
    return 2*pi / numThetaPoints_;
  }

  double GeorgiGlashowSu2TheoryAxisymmetric::getKineticScaleFactor(int cpt) const
  {
    if (cpt == 2) { return numThetaPoints_ / (2*pi); }
    else { return 2*pi / numThetaPoints_; }
  }


  double GeorgiGlashowSu2TheoryAxisymmetric::getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    double E = 0;
    double jacobian = getJacobian(site);

    for (int iCpt = 0; iCpt < 3; iCpt++)
    {
      for (int jCpt = 0; jCpt < 3; jCpt++)
      {

        double inverseMetricFactor = getInverseMetricCpt(site, iCpt)*getInverseMetricCpt(site, jCpt);
        double scaleFactor = getPlaquetteScaleFactor(iCpt, jCpt);
        if (iCpt >= jCpt) { continue; }

        Matrix avgPlaquette = 0.25*getPlaquette(field, site, iCpt, jCpt);
        LATfield2::Site tempSite = shiftSite(site, iCpt, -1);
        avgPlaquette = avgPlaquette + 0.25*getPlaquette(field, tempSite, iCpt, jCpt);
        tempSite = shiftSite(site, jCpt, -1);
        avgPlaquette = avgPlaquette + 0.25*getPlaquette(field, tempSite, iCpt, jCpt);
        tempSite = shiftSite(tempSite, iCpt, -1);
        avgPlaquette = avgPlaquette + 0.25*getPlaquette(field, tempSite, iCpt, jCpt);

        E += 2.0/pow(gaugeCoupling_,2)*scaleFactor*inverseMetricFactor*(2 - real(trace(avgPlaquette)));
      }
    }

    // Covariant derivative
    for (int dir = 0; dir < 3; dir++)
    {
      LATfield2::Site siteShiftedBwd = shiftSite(site, dir, -1);
      double inverseMetricCpt = getInverseMetricCpt(site, dir);
      double scaleFactor = getKineticScaleFactor(dir);

      Matrix covDeriv = getCovDeriv(field, site, dir);
      E += 0.5*inverseMetricCpt*scaleFactor*real(trace(covDeriv * covDeriv));

      Matrix covDerivShiftedBwd = getCovDeriv(field, siteShiftedBwd, dir);
      E += 0.5*inverseMetricCpt*scaleFactor*real(trace(covDerivShiftedBwd * covDerivShiftedBwd));
    }

    // Higgs potential
    E += real(selfCoupling_*pow(2.0*pow(getHiggsMagnitude(field, site),2) - pow(vev_, 2),2));

    return numThetaPoints_*jacobian*E;
  }

  Matrix GeorgiGlashowSu2TheoryAxisymmetric::getCovDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const
  {
    LATfield2::Site siteShiftedFwd = shiftSite(site, dir, 1);
    Matrix scalarMat = getHiggsMagnitude(field, site)*pauli3;
    Matrix scalarMatShiftedFwd = getHiggsMagnitude(field, siteShiftedFwd)*pauli3;
    Matrix gaugeMat = getSu2Link(field, site, dir);

    return gaugeMat*scalarMatShiftedFwd*conjugateTranspose(gaugeMat) - scalarMat;
  }

  monsta::Matrix GeorgiGlashowSu2TheoryAxisymmetric::getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
  {
    monsta::Matrix grad(2);
    double jacobian = getJacobian(site);

    if (matIdx < 3)
    {
      LATfield2::Site siteShiftedFwd = shiftSite(site, matIdx, 1);
      LATfield2::Site siteShiftedBwd = shiftSite(site, matIdx, -1);
      double inverseMetricCpt = getInverseMetricCpt(site, matIdx); // This only works because \partial_\theta g^{\theta \theta} = 0
      double kineticScaleFactor = getKineticScaleFactor(matIdx);
      double jacobianShiftedFwd = getJacobian(siteShiftedFwd);

      // Derivative of Wilson action
      Matrix plaquetteDerivMat(2);

      for (int dir = 0; dir < 3; dir++)
      {
        double plaquetteScaleFactor = getPlaquetteScaleFactor(matIdx, dir);
        if (dir == matIdx) { continue; }
        LATfield2::Site siteShiftedDown = shiftSite(site, dir, -1); 
        plaquetteDerivMat = plaquetteDerivMat + plaquetteScaleFactor*getPlaquetteAveragedMetricFactors(site, matIdx, dir)*getStaple(field, site, matIdx, dir, true);
        plaquetteDerivMat = plaquetteDerivMat + plaquetteScaleFactor*getPlaquetteAveragedMetricFactors(siteShiftedDown, matIdx, dir)*getStaple(field, site, matIdx, dir, false);
      }

      grad = grad - 2.0/pow(gaugeCoupling_,2)*plaquetteDerivMat;


      // Derivative of kinetic term
      Matrix scalarMat = getHiggsMagnitude(field, site)*pauli3;
      Matrix scalarMatShiftedFwd = getHiggsMagnitude(field, siteShiftedFwd)*pauli3;
      Matrix gaugeMat = getSu2Link(field, site, matIdx);
      Matrix kineticDerivMat = gaugeMat*scalarMatShiftedFwd*scalarMatShiftedFwd;
      kineticDerivMat = kineticDerivMat + gaugeMat*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(scalarMatShiftedFwd);
      kineticDerivMat = kineticDerivMat - conjugateTranspose(scalarMat)*gaugeMat*conjugateTranspose(scalarMatShiftedFwd);
      kineticDerivMat = kineticDerivMat - scalarMat*gaugeMat*scalarMatShiftedFwd;
      grad = grad + kineticScaleFactor*(jacobian + jacobianShiftedFwd)*inverseMetricCpt*kineticDerivMat;

      // Project to SU(2) Lie group: COMMENT OUT IF TESTING GRADIENTS
      // grad = grad - 0.5*trace(grad*conjugateTranspose(getSu2Link(field, site, matIdx)))*getSu2Link(field, site, matIdx);

    }
    else
    {
      for (int dir = 0; dir < 3; dir++)
      {
        // Deriviative of kinetic term
        LATfield2::Site siteShiftedFwd = shiftSite(site, dir, 1);
        LATfield2::Site siteShiftedBwd = shiftSite(site, dir, -1);
        double inverseMetricCpt = getInverseMetricCpt(site, dir); // This only works because \partial_\theta g^{\theta \theta} = 0
        double scaleFactor = getKineticScaleFactor(dir);
        double jacobianShiftedFwd = getJacobian(siteShiftedFwd);
        double jacobianShiftedBwd = getJacobian(siteShiftedBwd);

        Matrix scalarMat = getHiggsMagnitude(field, site)*pauli3;
        Matrix scalarMatShiftedFwd = getHiggsMagnitude(field, siteShiftedFwd)*pauli3;
        Matrix gaugeMat = getSu2Link(field, site, dir);
        Matrix gaugeMatShiftedBwd = getSu2Link(field, siteShiftedBwd, dir);
        Matrix scalarMatShiftedBwd = getHiggsMagnitude(field, siteShiftedBwd)*pauli3;

        Matrix kineticDerivMat = 2*(jacobian + jacobianShiftedFwd)*conjugateTranspose(scalarMat);
        kineticDerivMat = kineticDerivMat - 2*(jacobian + jacobianShiftedFwd)*gaugeMat*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(gaugeMat);
        kineticDerivMat = kineticDerivMat - 2*(jacobian + jacobianShiftedBwd)*conjugateTranspose(gaugeMatShiftedBwd)*conjugateTranspose(scalarMatShiftedBwd)*gaugeMatShiftedBwd;
        kineticDerivMat = kineticDerivMat + 2*(jacobian + jacobianShiftedBwd)*conjugateTranspose(scalarMat);

        grad(0,0) = grad(0,0) + scaleFactor*inverseMetricCpt*kineticDerivMat(0,0);
      }

      // Derivative of Higgs Potential
      grad(0,0) = grad(0,0) + 8.0*jacobian*selfCoupling_*getHiggsMagnitude(field, site)*(2.0*pow(getHiggsMagnitude(field, site),2) - pow(vev_, 2));
    }

    return numThetaPoints_*grad;
  }

  std::complex<double> GeorgiGlashowSu2TheoryAxisymmetric::getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const
  {
    monsta::Matrix gradMat = getLocalGradient(field, site, matIdx);
    return gradMat(rowIdx, colIdx);
  }

  void GeorgiGlashowSu2TheoryAxisymmetric::postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
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
      field(site, matIdx, 0, 1) = 0;
      field(site, matIdx, 1, 0) = 0;
      field(site, matIdx, 1, 1) = 0;
    }
  }

  void GeorgiGlashowSu2TheoryAxisymmetric::applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const
  {
    int rSize = field.lattice().size(1);

    field.updateHalo();
    LATfield2::Site site(field.lattice());
    for (site.haloFirst(); site.haloTest(); site.haloNext())
    {
      // if (r < 0)
      // {
      //   setHiggsMagnitude(field, site, vev_ / sqrt(2));
      //   for (int cpt = 0; cpt < 3; cpt++)
      //   {
      //     setSu2Link(field, site, cpt, identity);         
      //   }
      // }
      // if (site.coord(1) < 0 || site.coord(1) >= rSize)
      // {
      //   double r = getRFromSite(site);
      //   setHiggsMagnitude(field, site, vev_ / sqrt(2));
      //   setSu2Link(field, site, 0, identity);
      //   setSu2Link(field, site, 1, identity);
      //   std::vector<double> su2Vec = {0, 0, 0.01*pow(r,2) / numThetaPoints_};
      //   setSu2Link(field, site, 2, vecToSu2(su2Vec));
      // }
    }
  }

  Matrix GeorgiGlashowSu2TheoryAxisymmetric::getSu2Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
  {
    return Matrix(field, site, cpt);
  }

  void GeorgiGlashowSu2TheoryAxisymmetric::setSu2Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt, monsta::Matrix su2Val) const
  {
    field(site, cpt, 0, 0) = su2Val(0, 0);
    field(site, cpt, 0, 1) = su2Val(0, 1);
    field(site, cpt, 1, 0) = su2Val(1, 0);
    field(site, cpt, 1, 1) = su2Val(1, 1);
  }

  double GeorgiGlashowSu2TheoryAxisymmetric::getHiggsMagnitude(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    return real(field(site, 3, 0, 0));
  }

  void GeorgiGlashowSu2TheoryAxisymmetric::setHiggsMagnitude(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, std::complex<double> higgsVal) const
  {
    field(site, 3, 0, 0) = higgsVal;
  }

  monsta::Matrix GeorgiGlashowSu2TheoryAxisymmetric::getPlaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix plaquette = getSu2Link(field, site, dir1);
    tempSite = shiftSite(site, dir1, +1);
    plaquette = plaquette*getSu2Link(field, tempSite, dir2);
    tempSite = shiftSite(site, dir2, +1);
    plaquette = plaquette*conjugateTranspose(getSu2Link(field, tempSite, dir1));
    tempSite = site;
    plaquette = plaquette*conjugateTranspose(getSu2Link(field, tempSite, dir2));

    return plaquette;
  }

  monsta::Matrix GeorgiGlashowSu2TheoryAxisymmetric::getStaple(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, bool isUp) const
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

  double GeorgiGlashowSu2TheoryAxisymmetric::getPlaquetteAveragedMetricFactors(LATfield2::Site &site, int dir1, int dir2) const
  {
    double avgJacobian = 0.25*getJacobian(site)*getInverseMetricCpt(site, dir1)*getInverseMetricCpt(site, dir2);;
    LATfield2::Site tempSite = shiftSite(site, dir1, +1);
    avgJacobian += 0.25*getJacobian(tempSite)*getInverseMetricCpt(tempSite, dir1)*getInverseMetricCpt(tempSite, dir2);
    tempSite = shiftSite(tempSite, dir2, +1);
    avgJacobian += 0.25*getJacobian(tempSite)*getInverseMetricCpt(tempSite, dir1)*getInverseMetricCpt(tempSite, dir2);
    tempSite = shiftSite(tempSite, dir1, -1);
    avgJacobian += 0.25*getJacobian(tempSite)*getInverseMetricCpt(tempSite, dir1)*getInverseMetricCpt(tempSite, dir2);

    return avgJacobian;
  }

  std::complex<double> GeorgiGlashowSu2TheoryAxisymmetric::getU1Plaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const
  {
    LATfield2::Site tempSite(site);
    std::complex<double> u1Plaquette = getSu2Link(field, site, dir1).operator()(0,0);
    tempSite = shiftSite(site, dir1, +1);
    u1Plaquette = u1Plaquette*getSu2Link(field, tempSite, dir2).operator()(0,0);
    tempSite = shiftSite(site, dir2, +1);
    u1Plaquette = u1Plaquette*conj(getSu2Link(field, tempSite, dir1).operator()(0,0));
    tempSite = site;
    u1Plaquette = u1Plaquette*conj(getSu2Link(field, tempSite, dir2).operator()(0,0));

    return u1Plaquette;
  }

  double GeorgiGlashowSu2TheoryAxisymmetric::getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
  {
    int dir1 = (cpt + 1) % 3;
    int dir2 = (cpt + 2) % 3;

    std::complex<double> avgU1Plaquette = 0.25*getU1Plaquette(field, site, dir1, dir2);
    LATfield2::Site tempSite = shiftSite(site, dir1, -1);
    avgU1Plaquette = avgU1Plaquette + 0.25*getU1Plaquette(field, tempSite, dir1, dir2);
    tempSite = shiftSite(site, dir2, -1);
    avgU1Plaquette = avgU1Plaquette + 0.25*getU1Plaquette(field, tempSite, dir1, dir2);
    tempSite = shiftSite(tempSite, dir1, -1);
    avgU1Plaquette = avgU1Plaquette + 0.25*getU1Plaquette(field, tempSite, dir1, dir2);

    double scaleFactor;
    if (cpt == 2) { scaleFactor = 1; }
    else { scaleFactor = numThetaPoints_ / (2*pi); }
    double jacobian = getJacobian(site);
    double inverseMetricFactor = getInverseMetricCpt(site, dir1)*getInverseMetricCpt(site, dir2);

    return 2./gaugeCoupling_ * jacobian * inverseMetricFactor * scaleFactor * arg(avgU1Plaquette);
  }

}

#endif