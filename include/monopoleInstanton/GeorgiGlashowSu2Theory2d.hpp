#ifndef GEORGIGLASHOWSU2THEORY2D
#define GEORGIGLASHOWSU2THEORY2D

#include "../Theory.hpp"
#include "../Su2Tools.hpp"

namespace monsta {
  class GeorgiGlashowSu2TheoryAxisymmetric2d: public Theory {
  public:
    GeorgiGlashowSu2TheoryAxisymmetric2d(double gaugeCoupling, double vev, double selfCoupling, double numThetaPoints);
    double getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    std::complex<double> getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const;
    monsta::Matrix getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;
    double getVev() const { return vev_; }
    double getGaugeCoupling() const { return gaugeCoupling_; }
    double getSelfCoupling() const { return selfCoupling_; }
    double getNumThetaPoints() const { return numThetaPoints_; }
    double getRFromSite(LATfield2::Site &site) const;
    double getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    double getHiggsMagnitude(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    void applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const;
    void postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;
    monsta::Matrix getSu2Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const;
    std::complex<double> getU1Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const { return 0; };
    void setSu2Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt, monsta::Matrix su2val) const;
    void setHiggsField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, monsta::Matrix higgsVal) const;  
    void setHiggsMagnitude(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, std::complex<double> higgsVal) const;
    Matrix getPlaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const;
    Matrix getAbelianPlaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const;
    std::complex<double> getU1Plaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const;
    Matrix getHiggsField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    Matrix getScalarField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const;


  private:
    int numFieldMatrices_ = 4;
    int numRows_ = 2;
    int numCols_ = 2;
    int matSize_ = 2;
    double gaugeCoupling_;
    double vev_;
    double selfCoupling_;
    bool tHooftLine_ = false;
    double numThetaPoints_;
    std::vector<Matrix> r0BoundaryVals;
    std::vector< std::vector<Matrix> > rMaxBoundaryVals;
    bool zPeriodic = false;

    Matrix getCovDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    Matrix getCovDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int gaugeCpt, int scalarCpt) const;
    Matrix getStaple(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, bool isUp) const;
    double getJacobian(LATfield2::Site &site) const;
    double getInverseMetricCpt(LATfield2::Site &site, int cpt) const;
    double getPlaquetteScaleFactor(int cpt1, int cpt2) const;
    LATfield2::Site shiftSite(LATfield2::Site &site, int dir, int sign) const;
    double getPlaquetteAveragedMetricFactors(LATfield2::Site &site, int dir1, int dir2) const;
    monsta::Matrix getU1Projector(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    monsta::Matrix getAbelianLink(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const;
    void applyConstantFieldBoundary(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    Matrix getKineticDerivWrtScalar(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int scalarCpt) const;

  };

  GeorgiGlashowSu2TheoryAxisymmetric2d::GeorgiGlashowSu2TheoryAxisymmetric2d(double gaugeCoupling, double vev, double selfCoupling, double numThetaPoints)
  : Theory(numFieldMatrices_, matSize_, matSize_), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling), numThetaPoints_(numThetaPoints)
  {}

  LATfield2::Site GeorgiGlashowSu2TheoryAxisymmetric2d::shiftSite(LATfield2::Site &site, int dir, int sign) const
  {
    if (dir > 1) { return site; }
    if (sign > 0) { return site + dir; }
    return site - dir;
  }

  double GeorgiGlashowSu2TheoryAxisymmetric2d::getRFromSite(LATfield2::Site &site) const
  {
    int rSize = site.lattice().size(1);

    // if (site.coord(1) == 0)
    // {
    //   return rSize - rSize / 2 + 0.5;
    // }

    return site.coord(1) - rSize / 2 + 0.5;
  }

  double GeorgiGlashowSu2TheoryAxisymmetric2d::getJacobian(LATfield2::Site &site) const
  {
    return abs(getRFromSite(site));
  }

  double GeorgiGlashowSu2TheoryAxisymmetric2d::getInverseMetricCpt(LATfield2::Site &site, int cpt) const
  {
    if (cpt != 2) { return 1; }
    return pow(getRFromSite(site), -2);
  }

  double GeorgiGlashowSu2TheoryAxisymmetric2d::getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    double E = 0;
    double jacobian = getJacobian(site);

    Matrix plaquette = getPlaquette(field, site, 0, 1);

    double jacobianTimesInvMetric = getPlaquetteAveragedMetricFactors(site, 0, 1);

    E += 2.0/pow(gaugeCoupling_,2)*jacobianTimesInvMetric*(2 - real(trace(plaquette)));

    // Covariant derivative of scalar fields
    for (int dir = 0; dir < 2; dir++)
    {
      for (int scalarCpt = 0; scalarCpt < 2; scalarCpt++)
      {
        LATfield2::Site siteShiftedFwd = site + dir;
        double inverseMetricCpt = getInverseMetricCpt(site, dir);
        double inverseMetricCptShiftedFwd = getInverseMetricCpt(siteShiftedFwd, dir);
        if (scalarCpt == 1) 
        {
          inverseMetricCpt = getInverseMetricCpt(site, 2); 
          inverseMetricCptShiftedFwd = getInverseMetricCpt(siteShiftedFwd, 2);
        }
        double jacobianShiftedFwd = getJacobian(siteShiftedFwd);

        Matrix covDeriv = getCovDeriv(field, site, dir, scalarCpt);
        if (scalarCpt == 1)
        {
          covDeriv = 1.0/gaugeCoupling_ * covDeriv;
        }
        // covDeriv.print();
        E += 0.5*(inverseMetricCpt*jacobian + inverseMetricCptShiftedFwd*jacobianShiftedFwd)*real(trace(covDeriv * covDeriv));
      }
    }

    // A_\theta-Higgs interaction
    double inverseMetricCpt = getInverseMetricCpt(site, 2);
    Matrix aTheta = getScalarField(field, site, 1);
    Matrix higgsMat = getScalarField(field, site, 0);
    Matrix gaugeHiggsCommutator = commutator(aTheta, higgsMat);
    E -= 0.5*jacobian*inverseMetricCpt*real(trace(gaugeHiggsCommutator*gaugeHiggsCommutator));

    // Higgs potential
    E += jacobian*real(selfCoupling_*pow(2.0*pow(getHiggsMagnitude(field, site),2) - pow(vev_, 2),2));

    return pi*E;
  }

  Matrix GeorgiGlashowSu2TheoryAxisymmetric2d::getCovDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const
  {
    LATfield2::Site siteShiftedFwd = site + dir;
    Matrix scalarMat = getHiggsField(field, site);
    Matrix scalarMatShiftedFwd = getHiggsField(field, siteShiftedFwd);
    Matrix gaugeMat = getSu2Link(field, site, dir);

    return gaugeMat*scalarMatShiftedFwd*conjugateTranspose(gaugeMat) - scalarMat;
  }

  Matrix GeorgiGlashowSu2TheoryAxisymmetric2d::getCovDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int gaugeCpt, int scalarCpt) const
  {
    LATfield2::Site siteShiftedFwd = site + gaugeCpt;
    Matrix scalarMat = getScalarField(field, site, scalarCpt);
    Matrix scalarMatShiftedFwd = getScalarField(field, siteShiftedFwd, scalarCpt);
    Matrix gaugeMat = getSu2Link(field, site, gaugeCpt);

    return gaugeMat*scalarMatShiftedFwd*conjugateTranspose(gaugeMat) - scalarMat;
  }

  monsta::Matrix GeorgiGlashowSu2TheoryAxisymmetric2d::getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
  {
    monsta::Matrix grad(2);
    double jacobian = getJacobian(site);

    if (matIdx < 2)
    {
      LATfield2::Site siteShiftedFwd = shiftSite(site, matIdx, 1);
      LATfield2::Site siteShiftedBwd = shiftSite(site, matIdx, -1);
      double inverseMetricCpt = getInverseMetricCpt(site, matIdx);
      double inverseMetricCptShiftedFwd = getInverseMetricCpt(siteShiftedFwd, matIdx);
      double jacobianShiftedFwd = getJacobian(siteShiftedFwd);

      // Derivative of Wilson action
      Matrix plaquetteDerivMat(2);


      int dir = 1 - matIdx;
      LATfield2::Site siteShiftedDown = site - dir; 
      plaquetteDerivMat = plaquetteDerivMat + getPlaquetteAveragedMetricFactors(site, matIdx, dir)*getStaple(field, site, matIdx, dir, true);
      plaquetteDerivMat = plaquetteDerivMat + getPlaquetteAveragedMetricFactors(siteShiftedDown, matIdx, dir)*getStaple(field, site, matIdx, dir, false);

      grad = grad - 2.0/pow(gaugeCoupling_,2)*plaquetteDerivMat;


      // Derivative of kinetic term
      for (int scalarCpt = 0; scalarCpt < 2; scalarCpt++)
      { 
        if (scalarCpt == 1) 
        {
          inverseMetricCpt = getInverseMetricCpt(site, 2); 
          inverseMetricCptShiftedFwd = getInverseMetricCpt(siteShiftedFwd, 2);
        }
        Matrix scalarMat = getScalarField(field, site, scalarCpt);
        Matrix scalarMatShiftedFwd = getScalarField(field, siteShiftedFwd, scalarCpt);
        Matrix gaugeMat = getSu2Link(field, site, matIdx);
        Matrix kineticDerivMat = gaugeMat*scalarMatShiftedFwd*scalarMatShiftedFwd;
        kineticDerivMat = kineticDerivMat + gaugeMat*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(scalarMatShiftedFwd);
        kineticDerivMat = kineticDerivMat - conjugateTranspose(scalarMat)*gaugeMat*conjugateTranspose(scalarMatShiftedFwd);
        kineticDerivMat = kineticDerivMat - scalarMat*gaugeMat*scalarMatShiftedFwd;
        if (scalarCpt == 1)
        {
          kineticDerivMat = 1.0 / gaugeCoupling_ * kineticDerivMat;
        }
        grad = grad + (jacobian*inverseMetricCpt + jacobianShiftedFwd*inverseMetricCptShiftedFwd)*kineticDerivMat;
      }

      // Project to SU(2) Lie group: COMMENT OUT IF TESTING GRADIENTS
      // grad = grad - 0.5*trace(grad*conjugateTranspose(getSu2Link(field, site, matIdx)))*getSu2Link(field, site, matIdx);

    }
    else if (matIdx == 2)
    {
      // Derivative of gauge kinetic term
      double inverseMetricCpt = getInverseMetricCpt(site, matIdx);
      grad = grad + 1.0/ gaugeCoupling_ * getKineticDerivWrtScalar(field, site, 1);

      // Derivative of commutator term
      Matrix aThetaMat = getScalarField(field, site, 1);
      Matrix higgsMat = getScalarField(field, site, 0);

      Matrix gaugeHiggsCommutator = commutator(aThetaMat, higgsMat);
      grad = grad + jacobian*inverseMetricCpt*commutator(gaugeHiggsCommutator, higgsMat);
    }
    else if (matIdx == 3)
    {
      // Derivative of kinetic term
      Matrix kineticDerivMat = getKineticDerivWrtScalar(field, site, 0);
      grad(0,0) = grad(0,0) + 2.0*kineticDerivMat(0,0);

      // Derivative of commutator
      double inverseMetricCpt = getInverseMetricCpt(site, 2);
      Matrix aThetaMat = getScalarField(field, site, 1);
      Matrix higgsMat = getScalarField(field, site, 0);

      Matrix gaugeHiggsCommutator = commutator(aThetaMat, higgsMat);
      Matrix commutatorDerivMat = commutator(aThetaMat, gaugeHiggsCommutator);

      grad(0,0) = grad(0,0) + 2.0*jacobian*inverseMetricCpt*commutatorDerivMat(0,0);

      // Derivative of Higgs Potential
      double trscalarSq = real(trace(getHiggsField(field,site)*getHiggsField(field,site)));
      grad(0,0) = grad(0,0) + 8.0*jacobian*selfCoupling_*getHiggsMagnitude(field, site)*(2.0*pow(getHiggsMagnitude(field, site),2) - pow(vev_, 2));
    }

    return pi*grad;
  }

  std::complex<double> GeorgiGlashowSu2TheoryAxisymmetric2d::getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const
  {
    monsta::Matrix gradMat = getLocalGradient(field, site, matIdx);
    return gradMat(rowIdx, colIdx);
  }

  Matrix GeorgiGlashowSu2TheoryAxisymmetric2d::getKineticDerivWrtScalar(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int scalarCpt) const
  {
    Matrix kineticDerivMat(2);

    for (int gaugeCpt = 0; gaugeCpt < 2; gaugeCpt++)
    {
      LATfield2::Site siteShiftedFwd = site + gaugeCpt;
      LATfield2::Site siteShiftedBwd = site - gaugeCpt;
      double inverseMetricCpt = getInverseMetricCpt(site, gaugeCpt);
      double inverseMetricCptShiftedFwd = getInverseMetricCpt(siteShiftedFwd, gaugeCpt);
      double inverseMetricCptShiftedBwd = getInverseMetricCpt(siteShiftedBwd, gaugeCpt);
      if (scalarCpt == 1) 
      {
        inverseMetricCpt = getInverseMetricCpt(site, 2); 
        inverseMetricCptShiftedFwd = getInverseMetricCpt(siteShiftedFwd, 2);
        inverseMetricCptShiftedBwd = getInverseMetricCpt(siteShiftedBwd, 2);
      }
      double jacobian = getJacobian(site);
      double jacobianShiftedFwd = getJacobian(siteShiftedFwd);
      double jacobianShiftedBwd = getJacobian(siteShiftedBwd);

      Matrix scalarMat = getScalarField(field, site, scalarCpt);
      Matrix scalarMatShiftedFwd = getScalarField(field, siteShiftedFwd, scalarCpt);
      Matrix gaugeMat = getSu2Link(field, site, gaugeCpt);
      Matrix gaugeMatShiftedBwd = getSu2Link(field, siteShiftedBwd, gaugeCpt);
      Matrix scalarMatShiftedBwd = getScalarField(field, siteShiftedBwd, scalarCpt);

      kineticDerivMat = kineticDerivMat + (jacobian*inverseMetricCpt + jacobianShiftedFwd*inverseMetricCptShiftedFwd)*conjugateTranspose(scalarMat);
      kineticDerivMat = kineticDerivMat - (jacobian*inverseMetricCpt + jacobianShiftedFwd*inverseMetricCptShiftedFwd)*gaugeMat*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(gaugeMat);

      kineticDerivMat = kineticDerivMat - (jacobian*inverseMetricCpt + jacobianShiftedBwd*inverseMetricCptShiftedBwd)*conjugateTranspose(gaugeMatShiftedBwd)*conjugateTranspose(scalarMatShiftedBwd)*gaugeMatShiftedBwd;
      kineticDerivMat = kineticDerivMat + (jacobian*inverseMetricCpt + jacobianShiftedBwd*inverseMetricCptShiftedBwd)*conjugateTranspose(scalarMat);
    }

    return kineticDerivMat;
  }

  void GeorgiGlashowSu2TheoryAxisymmetric2d::postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
  {
    if (matIdx < 2) // Project to SU(2)
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
    else if (matIdx == 2)
    {
      field(site, matIdx, 0, 0) = real(field(site, matIdx, 0, 0));
      field(site, matIdx, 1, 1) = -real(field(site, matIdx, 0, 0));
      field(site, matIdx, 0, 1) = conj(field(site, matIdx, 1, 0));
    }
    else
    {
      field(site, matIdx, 0, 0) = real(field(site, matIdx, 0, 0));
      field(site, matIdx, 0, 1) = 0;
      field(site, matIdx, 1, 0) = 0;
      field(site, matIdx, 1, 1) = 0;
    }
  }

  void GeorgiGlashowSu2TheoryAxisymmetric2d::applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const
  {
    int rSize = field.lattice().size(1);
    int zSize = field.lattice().size(0);

    field.updateHalo();
    LATfield2::Site site(field.lattice());
    for (site.haloFirst(); site.haloTest(); site.haloNext())
    {
      int rCoord = site.coord(1);
      int zCoord = site.coord(0);
      std::vector<double> su2Vec2(3);
      double r = getRFromSite(site);

      applyConstantFieldBoundary(field, site);

      if (!zPeriodic)
      {
        if (zCoord < 0 || zCoord > zSize - 1)
        {
          for (int matIdx = 0; matIdx < numFieldMatrices_; matIdx++)
          {
            Matrix boundaryMat(field, site, matIdx);
            if (matIdx < 3)
            {
              boundaryMat = pauli2*boundaryMat*pauli2;
              setSu2Link(field, site, matIdx, boundaryMat);
            }
            else if (matIdx == 3)
            {
              setHiggsField(field, site, boundaryMat);
            }
          }
        }
      }
    }
  }


  void GeorgiGlashowSu2TheoryAxisymmetric2d::applyConstantFieldBoundary(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    int rCoord = site.coord(1);
    int rMax = field.lattice().size(1) - 1;
    double r = getRFromSite(site);
    double zSize = field.lattice().size(0);
    double z = site.coord(0) - zSize/2 + 0.5;

    if (rCoord < 0)
    {
      std::vector<double> su2Vec2 = {0, 0, 0};
      for (int matIdx = 0; matIdx < 2; matIdx++)
      {
        setSu2Link(field, site, matIdx, identity);
      }
      setSu2Link(field, site, 2, vecToSu2LieAlg(su2Vec2));
      setHiggsMagnitude(field, site, vev_ / sqrt(2));
    }

    if (rCoord > rMax)
    {
      std::vector<double> su2Vec2 = {0, 0, 1};
      for (int matIdx = 0; matIdx < 2; matIdx++)
      {
        setSu2Link(field, site, matIdx, identity);
      }
      setSu2Link(field, site, 2, vecToSu2LieAlg(su2Vec2));
      setHiggsMagnitude(field, site, vev_ / sqrt(2));
    }
  }

  Matrix GeorgiGlashowSu2TheoryAxisymmetric2d::getSu2Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
  {
    return Matrix(field, site, cpt);
  }

  void GeorgiGlashowSu2TheoryAxisymmetric2d::setSu2Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt, monsta::Matrix su2Val) const
  {
    field(site, cpt, 0, 0) = su2Val(0, 0);
    field(site, cpt, 0, 1) = su2Val(0, 1);
    field(site, cpt, 1, 0) = su2Val(1, 0);
    field(site, cpt, 1, 1) = su2Val(1, 1);
  }

  Matrix GeorgiGlashowSu2TheoryAxisymmetric2d::getHiggsField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    return field(site, 3, 0, 0)*pauli3;
  }

  Matrix GeorgiGlashowSu2TheoryAxisymmetric2d::getScalarField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
  {
    if (cpt == 0) { return getHiggsField(field, site); }
    else { return getSu2Link(field, site, 2); }
  }

  double GeorgiGlashowSu2TheoryAxisymmetric2d::getHiggsMagnitude(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    return real(field(site, 3, 0, 0));
  }

  void GeorgiGlashowSu2TheoryAxisymmetric2d::setHiggsField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, monsta:: Matrix higgsVal) const
  {
    field(site, 3, 0, 0) = higgsVal(0, 0);
    field(site, 3, 0, 1) = higgsVal(0, 1);
    field(site, 3, 1, 0) = higgsVal(1, 0);
    field(site, 3, 1, 1) = higgsVal(1, 1);
  }

  void GeorgiGlashowSu2TheoryAxisymmetric2d::setHiggsMagnitude(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, std::complex<double> higgsVal) const
  {
    setHiggsField(field, site, higgsVal*pauli3);
  }

  monsta::Matrix GeorgiGlashowSu2TheoryAxisymmetric2d::getPlaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix plaquette = getSu2Link(field, site, dir1);
    tempSite = site + dir1;
    plaquette = plaquette*getSu2Link(field, tempSite, dir2);
    tempSite = site + dir2;
    plaquette = plaquette*conjugateTranspose(getSu2Link(field, tempSite, dir1));
    tempSite = site;
    plaquette = plaquette*conjugateTranspose(getSu2Link(field, tempSite, dir2));

    return plaquette;
  }

  monsta::Matrix GeorgiGlashowSu2TheoryAxisymmetric2d::getStaple(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, bool isUp) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix staple(field, site, dir2);

    if (isUp)
    {
      tempSite = site + dir2;
      staple = staple*Matrix(field, tempSite, dir1);
      tempSite = site + dir1;
      staple = staple*conjugateTranspose(Matrix(field, tempSite, dir2));
    }
    else
    {
      // if (site.coord(1) == 0 && dir1 == 0 && dir2 == 1) { return Matrix(2); }
      // if (site.coord(1) == 0 && dir1 == 2 && dir2 == 1) { return Matrix(2); }
      tempSite = site - dir2;
      staple = conjugateTranspose(Matrix(field, tempSite, dir2));
      staple = staple*Matrix(field, tempSite, dir1);
      tempSite = tempSite + dir1;
      staple = staple*Matrix(field, tempSite, dir2);
    }

    return staple;

  }

  double GeorgiGlashowSu2TheoryAxisymmetric2d::getPlaquetteAveragedMetricFactors(LATfield2::Site &site, int dir1, int dir2) const
  {
    double avgJacobian = 0.25*getJacobian(site)*getInverseMetricCpt(site, dir1)*getInverseMetricCpt(site, dir2);;
    LATfield2::Site tempSite = site + dir1;
    avgJacobian += 0.25*getJacobian(tempSite)*getInverseMetricCpt(tempSite, dir1)*getInverseMetricCpt(tempSite, dir2);
    tempSite = tempSite + dir2;
    avgJacobian += 0.25*getJacobian(tempSite)*getInverseMetricCpt(tempSite, dir1)*getInverseMetricCpt(tempSite, dir2);
    tempSite = tempSite - dir1;
    avgJacobian += 0.25*getJacobian(tempSite)*getInverseMetricCpt(tempSite, dir1)*getInverseMetricCpt(tempSite, dir2);

    return avgJacobian;
  }

  monsta::Matrix GeorgiGlashowSu2TheoryAxisymmetric2d::getU1Projector(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    monsta::Matrix scalarMat= getHiggsField(field, site);
    double trScalarSq = real(monsta::trace(scalarMat*scalarMat));

    double zeroTol = 1e-15;
    if (trScalarSq < zeroTol)
    {
      return 0.5*monsta::identity;
    }

    return 0.5*(monsta::identity + sqrt(2./trScalarSq)*scalarMat);
  }


  std::complex<double> GeorgiGlashowSu2TheoryAxisymmetric2d::getU1Plaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const
  {
    LATfield2::Site tempSite(site);
    std::complex<double> u1Plaquette = getSu2Link(field, site, dir1).operator()(0,0);
    tempSite = site + dir1;
    u1Plaquette = u1Plaquette*getSu2Link(field, tempSite, dir2).operator()(0,0);
    tempSite = site + dir2;
    u1Plaquette = u1Plaquette*conj(getSu2Link(field, tempSite, dir1).operator()(0,0));
    tempSite = site;
    u1Plaquette = u1Plaquette*conj(getSu2Link(field, tempSite, dir2).operator()(0,0));

    return u1Plaquette;
  }


  double GeorgiGlashowSu2TheoryAxisymmetric2d::getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
  {
    double metricFactor = 1;//getPlaquetteAveragedMetricFactors(field, site, dir1, dir2);
    if (cpt == 2)
    {
      int dir1 = (cpt + 1) % 3;
      int dir2 = (cpt + 2) % 3;

      double u1Plaquette = arg(getU1Plaquette(field, site, dir1, dir2));

      return 2./gaugeCoupling_ * metricFactor * u1Plaquette;
    }
    else
    {
      int derivDir = 1 - cpt;
      Matrix covDeriv = getCovDeriv(field, site, derivDir, 1);

      return pow(-1, cpt) * 1./gaugeCoupling_ * metricFactor * real(trace(pauli3 * covDeriv));
    }

  }
}

#endif