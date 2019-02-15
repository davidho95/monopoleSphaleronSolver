#ifndef GEORGIGLASHOWSU2THEORY_HPP
#define GEORGIGLASHOWSU2THEORY_HPP

#include "LATfield2.hpp"
#include "Theory.hpp"
#include "Matrix.hpp"
#include <cmath>
#include <cstdio>

namespace monsta {
  class GeorgiGlashowSu2Theory: public Theory {
  public:
    GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling);
    GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling, bool skipScalar);
    GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling, int *monopolePos1, int *monopolePos2);
    GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling, int* monopolePos1, int*monopolePos2, bool skipScalar);
    int getScalarFieldStart();

    double getLocalEnergyDensity(LATfield2::Field<double> &field, LATfield2::Site &site,
      LATfield2::Field<double> &precomputedValues) const;
    double getLocalGradient(LATfield2::Field<double> &field, LATfield2::Site &site,
      int cpt, LATfield2::Field<double> &precomputedValues) const;
    monsta::Matrix vecToLieAlg(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt) const;
    monsta::Matrix vecToSu2(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt, int twist=1) const;
    monsta::Matrix getSu2Deriv(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt, int derivCpt, int twist=1) const;
    double su2ToVec(monsta::Matrix su2Mat, int cpt) const;
    monsta::Matrix getKineticTerm(LATfield2::Field<double> &field, LATfield2::Site &site, int dir) const;
    monsta::Matrix getPlaquette(LATfield2::Field<double> &field, LATfield2::Site &site, int dir1, int dir2, int twist=1) const;
    monsta::Matrix getPlaquetteDeriv(
    LATfield2::Field<double> &field, LATfield2::Site &site, int dir1, int dir2, int derivIdx, int derivDir, int derivCpt) const;
    monsta::Matrix getU1Projector(LATfield2::Field<double> &field, LATfield2::Site &site) const;
    monsta::Matrix getU1Projection(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt) const;
    monsta::Matrix getU1Plaquette(LATfield2::Field<double> &field, LATfield2::Site &sitem, int dir1, int dir2) const;
    double getMagneticField(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt) const;
    bool stringTest(LATfield2::Site &site) const;

  private:
    double gaugeCoupling_;
    double vev_;
    double selfCoupling_;
    int *monopolePos2_;
    int *monopolePos1_;
    bool noMonopoles_ = true;
    bool skipScalar_ = false;

    int scalarFieldStart_ = 9;

    // U(2) generators
    monsta::Matrix identity_;
    monsta::Matrix pauli1_;
    monsta::Matrix pauli2_;
    monsta::Matrix pauli3_;
  };

    GeorgiGlashowSu2Theory::GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling, int* monopolePos1, int*monopolePos2, bool skipScalar)
  : Theory(), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling), monopolePos1_(monopolePos1), monopolePos2_(monopolePos2),
    identity_(2), pauli1_(2), pauli2_(2), pauli3_(2), noMonopoles_(false), skipScalar_(skipScalar)
  {
    setNumSpatialCpts(12);  // 3 x 3 gauge field, 3 scalar

    identity_(0,0) = 1;
    identity_(1,1) = 1;

    pauli1_(0,1) = 1;
    pauli1_(1,0) = 1;

    pauli2_(0,1) = -1i;
    pauli2_(1,0) = 1i;

    pauli3_(0,0) = 1;
    pauli3_(1,1) = -1;
  }

  GeorgiGlashowSu2Theory::GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling, int* monopolePos1, int*monopolePos2)
  : Theory(), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling), monopolePos1_(monopolePos1), monopolePos2_(monopolePos2),
    identity_(2), pauli1_(2), pauli2_(2), pauli3_(2), noMonopoles_(false)
  {
    setNumSpatialCpts(12);  // 3 x 3 gauge field, 3 scalar

    identity_(0,0) = 1;
    identity_(1,1) = 1;

    pauli1_(0,1) = 1;
    pauli1_(1,0) = 1;

    pauli2_(0,1) = -1i;
    pauli2_(1,0) = 1i;

    pauli3_(0,0) = 1;
    pauli3_(1,1) = -1;
  }

  GeorgiGlashowSu2Theory::GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling)
  : Theory(), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling),
    identity_(2), pauli1_(2), pauli2_(2), pauli3_(2)
  {
    setNumSpatialCpts(12);  // 3 x 3 gauge field, 3 scalar

    identity_(0,0) = 1;
    identity_(1,1) = 1;

    pauli1_(0,1) = 1;
    pauli1_(1,0) = 1;

    pauli2_(0,1) = -1i;
    pauli2_(1,0) = 1i;

    pauli3_(0,0) = 1;
    pauli3_(1,1) = -1;
  }

  GeorgiGlashowSu2Theory::GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling, bool skipScalar)
  : Theory(), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling), skipScalar_(skipScalar),
    identity_(2), pauli1_(2), pauli2_(2), pauli3_(2)
  {
    setNumSpatialCpts(12);  // 3 x 3 gauge field, 3 scalar

    identity_(0,0) = 1;
    identity_(1,1) = 1;

    pauli1_(0,1) = 1;
    pauli1_(1,0) = 1;

    pauli2_(0,1) = -1i;
    pauli2_(1,0) = 1i;

    pauli3_(0,0) = 1;
    pauli3_(1,1) = -1;
  }

  int GeorgiGlashowSu2Theory::getScalarFieldStart()
  {
    return scalarFieldStart_;
  }

  double GeorgiGlashowSu2Theory::getLocalEnergyDensity(LATfield2::Field<double> &field, LATfield2::Site &site,
    LATfield2::Field<double> &precomputedValues) const
  {
    double E = 0;

    int numDimensions = 3;

    int twist = 1;

    // Wilson action
    E += 2/gaugeCoupling_*(2 - real(trace(getPlaquette(field, site, 1, 2))));
    E += 2/gaugeCoupling_*(2 - real(trace(getPlaquette(field, site, 2, 0))));
    E += 2/gaugeCoupling_*(2 - real(trace(getPlaquette(field, site, 0, 1))));

    // Covariant Derivative
    for (int ii = 0; ii < numDimensions; ii++)
    {
      E += real(trace(getKineticTerm(field, site, ii)));
    }

    // Higgs Potential
    double scalarVecNormSq = pow(field(site, scalarFieldStart_), 2)
                   + pow(field(site, scalarFieldStart_ + 1), 2)
                   + pow(field(site, scalarFieldStart_ + 2), 2);
    E += selfCoupling_*pow((2*scalarVecNormSq - pow(vev_,2)),2);

    return E;
  }

  double GeorgiGlashowSu2Theory::getLocalGradient(LATfield2::Field<double> &field, LATfield2::Site &site,
    int cpt, LATfield2::Field<double> &precomputedValues) const
  {
    int numDimensions = 3;

    int vectorCpt = cpt / numDimensions;
    int spatialCpt = cpt % numDimensions;

    monsta::Matrix pauliCpt(2);
    monsta::Matrix scalarMat = vecToLieAlg(field, site, 3);

    double grad = 0;

    if (cpt < scalarFieldStart_)
    {
      monsta::Matrix gaugeMat = vecToSu2(field, site, vectorCpt);
      LATfield2::Site shiftedSite = site+vectorCpt;
      monsta::Matrix scalarMatShiftedFwd = vecToLieAlg(field, shiftedSite, 3);
      // Derivative of Wilson action
      int dir1 = (vectorCpt + 1) % numDimensions;
      int dir2 = (vectorCpt + 2) % numDimensions;

      int derivIdx = site.index();

      grad += real(trace(getPlaquetteDeriv(field, site, vectorCpt, dir1, derivIdx, vectorCpt, spatialCpt)));
      grad += real(trace(getPlaquetteDeriv(field, site, dir2, vectorCpt, derivIdx, vectorCpt, spatialCpt)));

      LATfield2::Site tempSite(site);
      tempSite = tempSite-dir1;
      grad += real(trace(getPlaquetteDeriv(field, tempSite, vectorCpt, dir1, derivIdx, vectorCpt, spatialCpt)));

      tempSite = tempSite+dir1-dir2;
      grad += real(trace(getPlaquetteDeriv(field, tempSite, dir2, vectorCpt, derivIdx, vectorCpt, spatialCpt)));

      grad = -2/gaugeCoupling_*grad;

      // Derivative of kinetic term
      grad -= 4*real(trace(scalarMat*getSu2Deriv(field, site, vectorCpt, spatialCpt)*scalarMatShiftedFwd*conjugateTranspose(gaugeMat)));

    } else {
      if (skipScalar_) { return 0; }
      int xCoord = site.coord(0);
      int yCoord = site.coord(1);
      int zCoord = site.coord(2);

      if (!noMonopoles_)
      {
        if (xCoord == monopolePos1_[0] && yCoord == monopolePos1_[1] && zCoord == monopolePos1_[2])
        {
          return 0;
        } else if (xCoord == monopolePos2_[0] && yCoord == monopolePos2_[1] && zCoord == monopolePos2_[2])
        {
          return 0;
        }
      }
      
      switch (spatialCpt)
      {
          case 0:
            pauliCpt = pauli1_;
            break;
          case 1:
            pauliCpt = pauli2_;
            break;
          case 2:
            pauliCpt = pauli3_;
            break;
        }

      // Derivative of kinetic term
      grad += 12*field(site, cpt);     

      for (int ii = 0; ii < 3; ii++)
      {
        monsta::Matrix gaugeMat = vecToSu2(field, site, ii);
        LATfield2::Site shiftedSite = site+ii;
        monsta::Matrix scalarMatShiftedFwd = vecToLieAlg(field, shiftedSite, 3);
        shiftedSite = shiftedSite-ii-ii;
        monsta::Matrix scalarMatShiftedBwd = vecToLieAlg(field, shiftedSite, 3);
        monsta::Matrix gaugeMatShiftedBwd = vecToSu2(field, shiftedSite, ii);
        grad -= real(trace(pauliCpt*gaugeMat*scalarMatShiftedFwd*conjugateTranspose(gaugeMat)));
        grad -= real(trace(scalarMatShiftedBwd*gaugeMatShiftedBwd*pauliCpt*conjugateTranspose(gaugeMatShiftedBwd)));
          }
      // Derivative of Higgs Potential
      double scalarVecNormSq = pow(field(site, scalarFieldStart_), 2)
                   + pow(field(site, scalarFieldStart_ + 1), 2)
                   + pow(field(site, scalarFieldStart_ + 2), 2);
      grad += 4*selfCoupling_*field(site, cpt)*(2*scalarVecNormSq - pow(vev_, 2));

      grad = 2*grad;
    }

    return grad;

  }

  bool GeorgiGlashowSu2Theory::stringTest(LATfield2::Site &site) const
  {
    int xCoord = site.coord(0);
    int yCoord = site.coord(1);
    int zCoord = site.coord(2);

    if (yCoord == monopolePos1_[1] && zCoord == monopolePos1_[2])
    {
      if (xCoord >= monopolePos1_[0] && xCoord < monopolePos2_[0]) { return true; }
    } else { return false; }
  }

  monsta::Matrix GeorgiGlashowSu2Theory::vecToLieAlg(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt) const
  {
    int su2VecSize = 3;
    int vecStartIdx = su2VecSize*cpt;

    return field(site, vecStartIdx)*pauli1_
          +field(site, vecStartIdx+1)*pauli2_
          +field(site, vecStartIdx+2)*pauli3_;
  }

  monsta::Matrix GeorgiGlashowSu2Theory::vecToSu2(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt, int twist) const
  {
    int su2VecSize = 3;
    int vecStartIdx = su2VecSize*cpt;

    double vecCpt1 = field(site, vecStartIdx);
    double vecCpt2 = field(site, vecStartIdx + 1);
    double vecCpt3 = field(site, vecStartIdx + 2);

    double vecNorm = sqrt(pow(vecCpt1,2) + pow(vecCpt2,2) + pow(vecCpt3, 2));

    double zeroTol = 1e-10;
    if (vecNorm < zeroTol)
    {
      return identity_;
    }

    return cos(vecNorm)*identity_ + twist*1i*(sin(vecNorm)/vecNorm)*(vecCpt1*pauli1_ + vecCpt2*pauli2_ + vecCpt3*pauli3_);
  }

  monsta::Matrix GeorgiGlashowSu2Theory::getSu2Deriv(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt, int derivCpt, int twist) const
  {
    int su2VecSize = 3;
    int vecStartIdx = su2VecSize*cpt;

    double vecCpt1 = field(site, vecStartIdx);
    double vecCpt2 = field(site, vecStartIdx + 1);
    double vecCpt3 = field(site, vecStartIdx + 2);

    double vecCptDerivDir = field(site, vecStartIdx+derivCpt);
    monsta::Matrix pauliCpt(2);
    switch(derivCpt)
    {
      case 0:
        pauliCpt = pauli1_;
        break;
      case 1:
        pauliCpt = pauli2_;
        break;
      case 2:
        pauliCpt = pauli3_;
        break;
    }

    double vecNorm = sqrt(pow(vecCpt1,2) + pow(vecCpt2,2) + pow(vecCpt3, 2));
    double cosVecNorm = cos(vecNorm);
    double sinVecNorm = sin(vecNorm);
    monsta::Matrix pauliPart = vecCpt1*pauli1_ + vecCpt2*pauli2_ + vecCpt3*pauli3_;

    // cout << vecCptDerivDir << " " << vecNorm << " " << cosVecNorm << " " << sinVecNorm << " " << endl;

    if (vecNorm == 0) {
      return 1i * pauliCpt;
    }

    monsta::Matrix deriv = -(vecCptDerivDir*sinVecNorm/vecNorm)*identity_;
    // deriv.print();
    deriv = deriv + twist*1i*sinVecNorm/vecNorm * pauliCpt;
    // (twist*1i*sinVecNorm/vecNorm * pauliCpt).print();
    deriv = deriv + twist*1i*vecCptDerivDir*cosVecNorm*pauliPart/(pow(vecNorm,2));
    // (twist*1i*vecCptDerivDir*cosVecNorm*pauliPart/(pow(vecNorm,2))).print();
    deriv = deriv - twist*1i*vecCptDerivDir*sinVecNorm*pauliPart/pow(vecNorm,3);
    // (- twist*1i*vecCptDerivDir*sinVecNorm*pauliPart/pow(vecNorm,3)).print();
    return deriv;
  }

  double GeorgiGlashowSu2Theory::su2ToVec(monsta::Matrix su2Mat, int cpt) const
  {
    monsta::Matrix pauliCpt(2);
    switch (cpt)
    {
      case 0:
        pauliCpt = pauli1_;
        break;
      case 1:
        pauliCpt = pauli2_;
        break;
      case 2:
        pauliCpt = pauli3_;
        break;
    }

    double cosVecNorm = 0.5*real(trace(su2Mat));
    double vecNorm = acos(cosVecNorm);
    double sinVecNorm = sqrt(1 - pow(cosVecNorm,2));
    double sinVecNormProj = 0.5*imag(trace(su2Mat*pauliCpt));

    return (sinVecNormProj/sinVecNorm)*vecNorm;
  }

  monsta::Matrix GeorgiGlashowSu2Theory::getKineticTerm(LATfield2::Field<double> &field, LATfield2::Site &site, int dir) const
  {
    monsta::Matrix scalarMat = vecToLieAlg(field, site, 3);
    LATfield2::Site shiftedSite = site+dir;
    monsta::Matrix scalarMatShifted = vecToLieAlg(field, shiftedSite, 3);
    monsta::Matrix gaugeMat = vecToSu2(field, site, dir);

    return (gaugeMat*scalarMatShifted*conjugateTranspose(gaugeMat) - scalarMat)*(gaugeMat*scalarMatShifted*conjugateTranspose(gaugeMat) - scalarMat);
    // return 0.5*(scalarMat*scalarMat + scalarMatShifted*scalarMatShifted - 2*scalarMat*gaugeMat*scalarMatShifted*conjugateTranspose(gaugeMat));

  }

  monsta::Matrix GeorgiGlashowSu2Theory::getPlaquette(LATfield2::Field<double> &field, LATfield2::Site &site, int dir1, int dir2, int twist) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix plaquette = vecToSu2(field, site, dir1, twist);
    tempSite = tempSite+dir1;
    plaquette = plaquette*vecToSu2(field, tempSite, dir2, twist);
    tempSite = tempSite-dir1+dir2;
    plaquette = plaquette*conjugateTranspose(vecToSu2(field, tempSite, dir1, twist));
    tempSite = tempSite-dir2;
    plaquette = plaquette*conjugateTranspose(vecToSu2(field, tempSite, dir2, twist));

    return plaquette;
  }

  monsta::Matrix GeorgiGlashowSu2Theory::getPlaquetteDeriv(
    LATfield2::Field<double> &field, LATfield2::Site &site, int dir1, int dir2, int derivIdx, int derivDir, int derivCpt) const
  {

    int twist = 1;
    LATfield2::Site tempSite(site);
    monsta::Matrix plaquetteDeriv = identity_;
    if (tempSite.index() == derivIdx && derivDir == dir1)
    {
      plaquetteDeriv = plaquetteDeriv*getSu2Deriv(field, tempSite, dir1, derivCpt, twist);
    } else {
      plaquetteDeriv = plaquetteDeriv*vecToSu2(field, tempSite, dir1, twist);
    }
    tempSite = tempSite+dir1;
    if (tempSite.index() == derivIdx && derivDir == dir2)
    {
      plaquetteDeriv = plaquetteDeriv*getSu2Deriv(field, tempSite, dir2, derivCpt, twist);
    } else {
      plaquetteDeriv = plaquetteDeriv*vecToSu2(field, tempSite, dir2, twist);
    }
    tempSite = tempSite-dir1;
    tempSite = tempSite+dir2;
    if (tempSite.index() == derivIdx && derivDir == dir1)
    {
      plaquetteDeriv = plaquetteDeriv*conjugateTranspose(getSu2Deriv(field, tempSite, dir1, derivCpt, twist));
    } else {
      plaquetteDeriv = plaquetteDeriv*conjugateTranspose(vecToSu2(field, tempSite, dir1, twist));
    }
    tempSite = tempSite-dir2;
    if (tempSite.index() == derivIdx && derivDir == dir2)
    {
      plaquetteDeriv = plaquetteDeriv*conjugateTranspose(getSu2Deriv(field, tempSite, dir2, derivCpt, twist));
    } else {
      plaquetteDeriv = plaquetteDeriv*conjugateTranspose(vecToSu2(field, tempSite, dir2, twist));
    }
    return plaquetteDeriv;
  }

  monsta::Matrix GeorgiGlashowSu2Theory::getU1Projector(LATfield2::Field<double> &field, LATfield2::Site &site) const
  {
    monsta::Matrix scalarMat = vecToLieAlg(field, site, 3);
    double vecNormSq = real(trace(scalarMat*scalarMat));

    double zeroTol = 1e-10;
    if (vecNormSq < zeroTol) { return 0.5*identity_; }

    return 0.5*(identity_ + sqrt(2./vecNormSq)*scalarMat);
  }

  monsta::Matrix GeorgiGlashowSu2Theory::getU1Projection(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt) const
  {
    monsta::Matrix gaugeMat = vecToSu2(field, site, cpt);
    monsta::Matrix projector = getU1Projector(field, site);
    LATfield2::Site shiftedSite = site + cpt;
    monsta::Matrix projectorShiftedFwd = getU1Projector(field, shiftedSite);

    return projector*gaugeMat*projectorShiftedFwd;
  }

  monsta::Matrix GeorgiGlashowSu2Theory::getU1Plaquette(LATfield2::Field<double> &field, LATfield2::Site &site, int dir1, int dir2) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix u1Plaquette = getU1Projection(field, site, dir1);
    tempSite = tempSite+dir1;
    u1Plaquette = u1Plaquette*getU1Projection(field, tempSite, dir2);
    tempSite = tempSite-dir1+dir2;
    u1Plaquette = u1Plaquette*conjugateTranspose(getU1Projection(field, tempSite, dir1));
    tempSite = tempSite-dir2;
    u1Plaquette = u1Plaquette*conjugateTranspose(getU1Projection(field, tempSite, dir2));
    return u1Plaquette;
  }

  double GeorgiGlashowSu2Theory::getMagneticField(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt) const
  {
    int dir1 = (cpt + 1) % 3;
    int dir2 = (cpt + 2) % 3;

    return 2./gaugeCoupling_*arg(trace(getU1Plaquette(field, site, dir1, dir2)));
  }

}

#endif