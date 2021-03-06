#ifndef GEORGIGLASHOWSU2EOMTHEORY_HPP
#define GEORGIGLASHOWSU2EOMTHEORY_HPP

#include "../Theory.hpp"
#include "../Su2Tools.hpp"

namespace monsta {
  class GeorgiGlashowSu2EomTheory: public Theory {
  public:
    GeorgiGlashowSu2EomTheory(double gaugeCoupling, double vev, double selfCoupling);
    GeorgiGlashowSu2EomTheory(double gaugeCoupling, double vev, double selfCoupling, std::vector<int> boundaryConditions);
    GeorgiGlashowSu2EomTheory(double gaugeCoupling, double vev, double selfCoupling, std::vector<int> boundaryConditions, bool tHooftLine);
    double getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    std::complex<double> getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const;
    monsta::Matrix getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;
    double getVev() const { return vev_; }
    double getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    void postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;
    void applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const;
    monsta::Matrix getCovDerivFwd(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    monsta::Matrix getCovDerivBwd(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    monsta::Matrix getFwdBwdCovDerivSum(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    monsta::Matrix getKineticDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    monsta::Matrix getDirectedKineticDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir, bool isFwd) const;
    monsta::Matrix getInsert2(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    monsta::Matrix getDirectedInsert2(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir, bool isFwd) const;
    monsta::Matrix getPlaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const;
    monsta::Matrix getStaple(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, bool isUp) const;
    monsta::Matrix getHinge1(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int derivDir, int dir1, bool isPos1, int dir2, bool isPos2) const;
    monsta::Matrix getHinge2(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int derivDir, int dir1, bool isPos1, int dir2, bool isPos2) const;
    monsta::Matrix getHinge3(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int derivDir, int dir1, bool isPos1, int dir2, bool isPos2) const;
    monsta::Matrix getDblPlaquetteDeriv1(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int derivDir, int dir1, bool isPos1, int dir2, bool isPos2) const;
    monsta::Matrix getDblPlaquetteDeriv2(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int derivDir, int dir1, bool isPos1, int dir2, bool isPos2) const;
    monsta::Matrix getDblPlaquetteDeriv3(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int derivDir, int dir1, bool isPos1, int dir2, bool isPos2) const;
    monsta::Matrix getDblPlaquetteDeriv4(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int derivDir, int dir1, bool isPos1, int dir2, bool isPos2) const;
    monsta::Matrix getKinDerivInsertStaple(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, bool isUp, int insertIdx) const;
    monsta::Matrix getInsertStaple2(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, bool isUp, int insertIdx) const;
    monsta::Matrix getDirectedLink(LATfield2::Field< std::complex<double> > &field,  LATfield2::Site &site, int dir, bool isFwd) const;
    double getVevDistance(LATfield2::Field< std::complex<double> > &field, LATfield2::Site site) const;
    monsta::Matrix getTHooftHinge(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;

  private:
    bool tHooftLineCheck(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    void applyDirichletBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const;
    void applyCPeriodicBoundaryConditions(LATfield2::Field< std::complex<double> > &field, int dir) const;
    double getTHooftOperator(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir, std::complex<double> fluxQuanta=1) const;
    monsta::Matrix getTHooftDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, std::complex<double> fluxQuanta=1) const;

    int numFieldMatrices_ = 4;
    int numRows_ = 2;
    int numCols_ = 2;
    int matSize_ = 2;
    double gaugeCoupling_;
    double vev_;
    double selfCoupling_;
    std::vector<int> boundaryConditions_ = {1, 1, 1};
    bool tHooftLine_ = false;
    std::complex<double> fluxQuanta_ = 1;

    bool monopolesPresent_ = false;
    std::vector<int> monopolePos1_;
    std::vector<int> monopolePos2_;
  };

  GeorgiGlashowSu2EomTheory::GeorgiGlashowSu2EomTheory(double gaugeCoupling, double vev, double selfCoupling)
  : Theory(numFieldMatrices_, matSize_, matSize_), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling)
  {}

  GeorgiGlashowSu2EomTheory::GeorgiGlashowSu2EomTheory(double gaugeCoupling, double vev, double selfCoupling, std::vector<int> boundaryConditions)
  : Theory(numFieldMatrices_, matSize_, matSize_), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling),
    boundaryConditions_(boundaryConditions)
  {}

  GeorgiGlashowSu2EomTheory::GeorgiGlashowSu2EomTheory(double gaugeCoupling, double vev, double selfCoupling, std::vector<int> boundaryConditions, bool tHooftLine)
  : Theory(numFieldMatrices_, matSize_, matSize_), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling),
    boundaryConditions_(boundaryConditions), tHooftLine_(tHooftLine)
  {}

  double GeorgiGlashowSu2EomTheory::getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    double E = 0;

    // // Wilson action
    // E += 2.0/pow(gaugeCoupling_,2)*(2 - real(trace(getPlaquette(field, site, 1, 2))));
    // E += 2.0/pow(gaugeCoupling_,2)*(2 - real(trace(getPlaquette(field, site, 2, 0))));
    // E += 2.0/pow(gaugeCoupling_,2)*(2 - real(trace(getPlaquette(field, site, 0, 1))));

    // for (int ii = 0; ii < 3; ii++)
    // {
    //   if (tHooftLineCheck(field, site, ii))
    //   {
    //     E += getTHooftOperator(field, site, ii, fluxQuanta_);
    //   }
    // }

    // Covariant Derivative
    // for (int ii = 0; ii < 3; ii++)
    // {
    //   LATfield2::Site shiftedSite = site + ii;
    //   Matrix covDeriv = field(shiftedSite, 3, 0, 0)*Matrix(field, site, ii)*pauli3*conjugateTranspose(Matrix(field, site, ii)) - field(site, 3, 0, 0)*pauli3;
    //   E += real(trace(covDeriv*covDeriv));
    // }

    // Higgs Potential
    // E += real(selfCoupling_*pow(2.0*pow(field(site, 3, 0, 0),2) - pow(vev_, 2),2));

    for (int matIdx = 0; matIdx < 4; matIdx++)
    {
      monsta::Matrix grad(2);

      if (matIdx < 3)
      {
        int dir1 = (matIdx + 1) % 3;
        int dir2 = (matIdx + 2) % 3;

        // Derivative of Wilson action
        Matrix plaquetteDerivMat(2);
        plaquetteDerivMat = plaquetteDerivMat + getStaple(field, site, matIdx, dir1, true);
        plaquetteDerivMat = plaquetteDerivMat + getStaple(field, site, matIdx, dir1, false);
        plaquetteDerivMat = plaquetteDerivMat + getStaple(field, site, matIdx, dir2, true);
        plaquetteDerivMat = plaquetteDerivMat + getStaple(field, site, matIdx, dir2, false);

        plaquetteDerivMat = -2.0/pow(gaugeCoupling_,2)*plaquetteDerivMat;
        grad = grad + plaquetteDerivMat;
        Matrix tHooftDerivMat = 2.0/pow(gaugeCoupling_,2)*getTHooftDeriv(field, site, matIdx, fluxQuanta_);
        grad = grad - tHooftDerivMat;
        // grad.print();

        // Derivative of kinetic term
        Matrix kineticDerivMat = getKineticDeriv(field, site, matIdx);
        grad = grad + kineticDerivMat;
        // E += real(trace())
        E += real(trace(grad*conjugateTranspose(grad)));

        monsta::Matrix gradProj = grad*conjugateTranspose(Matrix(field, site, matIdx));

        E -= real(trace(gradProj*gradProj));
        // E -= 0.5*pow(real(trace(gradProj)),2);
        // COUT << real(trace(grad*conjugateTranspose(grad))) - real(trace(gradProj*gradProj)) << endl;

      } else {
        // if (abs(field(site, 3, 0, 0)) < 1e-15) { return grad; }
        Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
        Matrix kineticDerivMat(2);
        for (int dir = 0; dir < 3; dir++)
        {
          // Deriviative of kinetic term
          LATfield2::Site tempSite = site + dir;
          Matrix scalarMatShiftedFwd = field(tempSite, 3, 0, 0)*pauli3;
          Matrix gaugeMat(field, site, dir);

          tempSite = site - dir;
          Matrix gaugeMatShiftedBwd(field, tempSite, dir);
          Matrix scalarMatShiftedBwd = field(tempSite, 3, 0, 0)*pauli3;

          // Matrix kineticDerivMat(2);
          kineticDerivMat = kineticDerivMat + 2*conjugateTranspose(gaugeMatShiftedBwd)*gaugeMatShiftedBwd*conjugateTranspose(scalarMat)*conjugateTranspose(gaugeMatShiftedBwd)*gaugeMatShiftedBwd;
          kineticDerivMat = kineticDerivMat - 2*gaugeMat*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(gaugeMat);
          kineticDerivMat = kineticDerivMat - 2*conjugateTranspose(gaugeMatShiftedBwd)*conjugateTranspose(scalarMatShiftedBwd)*gaugeMatShiftedBwd;
          kineticDerivMat = kineticDerivMat + 2*conjugateTranspose(scalarMat);
        }
        grad = grad + kineticDerivMat;

        // Derivative of Higgs Potential
        Matrix higgsPotentialDeriv = 8.0*0.5*selfCoupling_*(trace(scalarMat*scalarMat) - pow(vev_, 2))*scalarMat;
        grad = grad + higgsPotentialDeriv;

        E += real(trace(grad*grad));
      }
    }

    return E;
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
  {
    monsta::Matrix grad(2);

    if (matIdx < 3)
    {
      grad = grad + getTHooftHinge(field, site, matIdx);
      for (int a = 0; a < 3; a++)
      {
        if (a == matIdx) { continue; }
        for (bool sgnA : {false, true})
        {
          for (int b = 0; b < 3; b++)
          {
            if (b == matIdx) { continue; }
            for (bool sgnB : {false, true})
            {
              // Pure gauge, magnitude
              grad = grad + 2*pow(2.0/pow(gaugeCoupling_,2),2)*getHinge1(field, site, matIdx, a, sgnA, b, sgnB);

              // Pure gauge, projected
              grad = grad - 2*pow(2.0/pow(gaugeCoupling_,2),2)*getDblPlaquetteDeriv1(field, site, matIdx, a, sgnA, b, sgnB);
              grad = grad - 2*pow(2.0/pow(gaugeCoupling_,2),2)*getDblPlaquetteDeriv4(field, site, matIdx, a, sgnA, b, sgnB);
            }
          }
        }
      }

      for (int a = 0; a < 3; a++)
      {
        for (bool sgnA : {false, true})
        {
          for (int b = 0; b < 3; b++)
          {
            if (b == matIdx) { continue; }
            if (b == a) { continue; }
            for (bool sgnB : {false, true})
            {
              // Pure gauge, magnitude
              grad = grad + 2*pow(2.0/pow(gaugeCoupling_,2),2)*getHinge2(field, site, matIdx, a, sgnA, b, sgnB);
              grad = grad + 2*pow(2.0/pow(gaugeCoupling_,2),2)*getHinge3(field, site, matIdx, a, sgnA, b, sgnB);
              // Pure gauge, projected
              grad = grad - 2*pow(2.0/pow(gaugeCoupling_,2),2)*getDblPlaquetteDeriv2(field, site, matIdx, a, sgnA, b, sgnB);
              grad = grad - 2*pow(2.0/pow(gaugeCoupling_,2),2)*getDblPlaquetteDeriv3(field, site, matIdx, a, sgnA, b, sgnB);
            }
          }
        }
      }

      LATfield2::Site siteShifted = site + matIdx;
      Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
      Matrix scalarMatShifted = field(siteShifted, 3, 0, 0)*pauli3;
      Matrix gaugeMat(field, site, matIdx);
      Matrix covDerivSum = getCovDerivFwd(field, site, matIdx) + conjugateTranspose(getCovDerivFwd(field, site, matIdx));

      for (int ii = 0; ii < 3; ii++)
      {
        if (ii == matIdx) { continue; }
        for (int stapleLinkIdx = 0; stapleLinkIdx < 3; stapleLinkIdx++)
        {
          for (bool stapleSgn : {true, false})
          {
            // Cross terms from squared magnitude
            grad = grad - 2*2.0/pow(gaugeCoupling_,2)*getKinDerivInsertStaple(field, site, matIdx, ii, stapleSgn, stapleLinkIdx);

            // Cross terms from projected square
            grad = grad + 4*2.0/pow(gaugeCoupling_,2)*getInsertStaple2(field, site, matIdx, ii, stapleSgn, stapleLinkIdx);
          }
        }

        Matrix stapleSum = getStaple(field, site, matIdx, ii, true) + getStaple(field, site, matIdx, ii, false);
        Matrix stapleSumPlus = stapleSum*conjugateTranspose(scalarMatShifted)*conjugateTranspose(gaugeMat);
        Matrix stapleSumPlus2 = gaugeMat*scalarMatShifted*conjugateTranspose(gaugeMat)*stapleSum*conjugateTranspose(gaugeMat);

        // Cross terms from squared magnitude
        grad = grad - 4*2.0/pow(gaugeCoupling_,2)*(stapleSumPlus + conjugateTranspose(stapleSumPlus))*gaugeMat*(scalarMatShifted + conjugateTranspose(scalarMatShifted));
        grad = grad - 4*2.0/pow(gaugeCoupling_,2)*covDerivSum*stapleSum*conjugateTranspose(scalarMatShifted);

          // Cross terms from projected square
          grad = grad + 4*2.0/pow(gaugeCoupling_,2)*(stapleSumPlus2 + conjugateTranspose(stapleSumPlus2))*gaugeMat*(scalarMatShifted + conjugateTranspose(scalarMatShifted));
          grad = grad + 4*2.0/pow(gaugeCoupling_,2)*covDerivSum*gaugeMat*conjugateTranspose(stapleSum)*gaugeMat*conjugateTranspose(scalarMatShifted);
          grad = grad + 4*2.0/pow(gaugeCoupling_,2)*stapleSum*conjugateTranspose(gaugeMat)*covDerivSum*gaugeMat*scalarMatShifted;
          grad = grad + 4*2.0/pow(gaugeCoupling_,2)*covDerivSum*gaugeMat*scalarMatShifted*conjugateTranspose(gaugeMat)*stapleSum;
      }

      // Square of kinetic term
      Matrix kineticDerivMat = getKineticDeriv(field, site, matIdx);
      grad = grad + 8*2.0
        *kineticDerivMat*conjugateTranspose(scalarMatShifted)*scalarMatShifted;//*conjugateTranspose(gaugeMat)*gaugeMat
        // *(scalarMatShifted + conjugateTranspose(scalarMatShifted));
      // grad = grad + 4*2.0
        // *gaugeMat*scalarMatShifted*conjugateTranspose(scalarMatShifted)*conjugateTranspose(gaugeMat)*covDerivSum*gaugeMat
        // *(scalarMatShifted + conjugateTranspose(scalarMatShifted));
      grad = grad + 4*2.0
        *covDerivSum*covDerivSum*gaugeMat*scalarMatShifted*conjugateTranspose(scalarMatShifted);

      // Square of projected kinetic term
      grad = grad - 4*2.0
        *gaugeMat*conjugateTranspose(scalarMatShifted)*conjugateTranspose(gaugeMat)*covDerivSum*gaugeMat
        *conjugateTranspose(scalarMatShifted)*(scalarMatShifted + conjugateTranspose(scalarMatShifted));
      grad = grad - 4*2.0
        *gaugeMat*scalarMatShifted*conjugateTranspose(gaugeMat)*covDerivSum*gaugeMat*scalarMatShifted
        *(scalarMatShifted + conjugateTranspose(scalarMatShifted));
      grad = grad - 8*2.0
        *covDerivSum*gaugeMat*scalarMatShifted*conjugateTranspose(gaugeMat)*covDerivSum*gaugeMat*scalarMatShifted;

        // Square of scalar derivative term
        Matrix fwdBwdCovDerivSum = getFwdBwdCovDerivSum(field, site);
        Matrix fwdBwdCovDerivSumShifted = getFwdBwdCovDerivSum(field, siteShifted);
        grad = grad + 16*fwdBwdCovDerivSum*gaugeMat*scalarMatShifted;
        grad = grad - 16*gaugeMat*scalarMatShifted*fwdBwdCovDerivSumShifted;
        grad = grad - 16*gaugeMat*fwdBwdCovDerivSumShifted*scalarMatShifted;
        grad = grad + 16*scalarMat*gaugeMat*fwdBwdCovDerivSumShifted;

        // Cross terms in scalar deriv term
        double vevDistance = getVevDistance(field, site);
        double vevDistanceShifted = getVevDistance(field, siteShifted);

        grad = grad + 128*0.5*selfCoupling_*vevDistanceShifted*gaugeMat*scalarMatShifted*scalarMatShifted;
        grad = grad - 64*0.5*selfCoupling_*vevDistance*scalarMat*gaugeMat*scalarMatShifted;
        grad = grad - 64*0.5*selfCoupling_*vevDistanceShifted*scalarMat*gaugeMat*scalarMatShifted;
    }
    else
    {
      Matrix fwdBwdCovDerivSum = getFwdBwdCovDerivSum(field, site);
      Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
      for (int ii = 0; ii < 3; ii++)
      {
      for (int jj = 0; jj < 3; jj++)
        {
          if (jj == ii) { continue; }
          LATfield2::Site siteShiftedFwd = site + ii;
          LATfield2::Site siteShiftedBwd = site - ii;
          LATfield2::Site siteShiftedUp = site + jj;
          LATfield2::Site siteShiftedDown = site - jj;

          Matrix gaugeMatFwd(field, site, ii);
          Matrix gaugeMatShiftedBwd(field, siteShiftedBwd, ii);

          Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
          Matrix scalarMatShiftedFwd = field(siteShiftedFwd, 3, 0, 0)*pauli3;

          Matrix stapleSumShiftedBwd = getStaple(field, siteShiftedBwd, ii, jj, true) + getStaple(field, siteShiftedBwd, ii, jj, false);

          // Cross terms from squared magnitude
          Matrix gradMat1 = conjugateTranspose(gaugeMatShiftedBwd)*stapleSumShiftedBwd*conjugateTranspose(scalarMat)*conjugateTranspose(gaugeMatShiftedBwd)*gaugeMatShiftedBwd;
          gradMat1 = gradMat1 + conjugateTranspose(gradMat1);
          grad(0,0) = grad(0,0) - 8*2.0/pow(gaugeCoupling_,2)*gradMat1(0,0);

          Matrix covDerivSumShiftedBwd = getCovDerivFwd(field, siteShiftedBwd, ii) + conjugateTranspose(getCovDerivFwd(field, siteShiftedBwd, ii));
          Matrix gradMat2 = conjugateTranspose(gaugeMatShiftedBwd)*covDerivSumShiftedBwd*stapleSumShiftedBwd;
          grad(0,0) = grad(0,0) - 8*2.0/pow(gaugeCoupling_,2)*real(gradMat2(0,0));

          Matrix stapleSum = getStaple(field, site, ii, jj, true) + getStaple(field, site, ii, jj, false);
          Matrix gradMat3 = stapleSum*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(gaugeMatFwd);
          gradMat3 = gradMat3 + conjugateTranspose(gradMat3);
          grad(0,0) = grad(0,0) + 8*2.0/pow(gaugeCoupling_,2)*gradMat3(0,0);

          // Cross terms from projected square
          Matrix gradMat7 = conjugateTranspose(stapleSumShiftedBwd)*gaugeMatShiftedBwd*conjugateTranspose(scalarMat);
          gradMat7 = gradMat7 + conjugateTranspose(gradMat7);
          grad(0,0) = grad(0,0) + 8*2.0/pow(gaugeCoupling_,2)*real(gradMat7(0,0));

          Matrix gradMat8 = conjugateTranspose(gaugeMatShiftedBwd)*covDerivSumShiftedBwd*gaugeMatShiftedBwd
            *conjugateTranspose(stapleSumShiftedBwd)*gaugeMatShiftedBwd;
          grad(0,0) = grad(0,0) + 8*2.0/pow(gaugeCoupling_,2)*real(gradMat8(0,0));

          Matrix gradMat9 = gaugeMatFwd*conjugateTranspose(stapleSum)*gaugeMatFwd*scalarMatShiftedFwd*conjugateTranspose(gaugeMatFwd);
          gradMat9 = gradMat9 + conjugateTranspose(gradMat9);
          grad(0,0) = grad(0,0) - 8*2.0/pow(gaugeCoupling_,2)*real(gradMat9(0,0));          

        }

      LATfield2::Site siteShiftedFwd = site + ii;
      LATfield2::Site siteShiftedBwd = site - ii;
      Matrix gaugeMat(field, site, ii);
      Matrix gaugeMatShiftedBwd(field, siteShiftedBwd, ii);
      Matrix scalarMatShiftedFwd = field(siteShiftedFwd, 3, 0, 0)*pauli3;
      Matrix scalarMatShiftedBwd = field(siteShiftedBwd, 3, 0, 0)*pauli3;

      Matrix covDerivSumShiftedBwd = getCovDerivFwd(field, siteShiftedBwd, ii) + conjugateTranspose(getCovDerivFwd(field, siteShiftedBwd, ii));
      Matrix covDerivSum = getCovDerivFwd(field, site, ii) + conjugateTranspose(getCovDerivFwd(field, site, ii));

      // Square of kinetic term
      Matrix gradMat4 = conjugateTranspose(gaugeMatShiftedBwd)*covDerivSumShiftedBwd*gaugeMatShiftedBwd*scalarMat
        *conjugateTranspose(scalarMat);//*conjugateTranspose(gaugeMatShiftedBwd)*gaugeMatShiftedBwd;
      gradMat4 = gradMat4 + conjugateTranspose(gradMat4);
      grad(0,0) = grad(0,0) + 16.0*gradMat4(0,0);

      Matrix gradMat5 = conjugateTranspose(gaugeMatShiftedBwd)*covDerivSumShiftedBwd*covDerivSumShiftedBwd*gaugeMatShiftedBwd*scalarMat;
      grad(0,0) = grad(0,0) + 16.0*gradMat5(0,0);

      Matrix gradMat6 = getKineticDeriv(field, site, ii)*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(gaugeMat);
      // Matrix gradMat6 = covDerivSum*gaugeMat*scalarMatShiftedFwd*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(gaugeMat);
      gradMat6 = gradMat6 + conjugateTranspose(gradMat6);
      grad(0,0) = grad(0,0) - 8.0*gradMat6(0,0);

      // Square of projected kinetic term
      Matrix gradMat10 = scalarMat*conjugateTranspose(gaugeMatShiftedBwd)*covDerivSumShiftedBwd*gaugeMatShiftedBwd*scalarMat;
      gradMat10 = gradMat10 + conjugateTranspose(gradMat10);
      grad(0,0) = grad(0,0) - 16.0*real(gradMat10(0,0));

      Matrix gradMat11 = conjugateTranspose(gaugeMatShiftedBwd)*covDerivSumShiftedBwd*gaugeMatShiftedBwd
        *scalarMat*conjugateTranspose(gaugeMatShiftedBwd)*covDerivSumShiftedBwd*gaugeMatShiftedBwd;
      grad(0,0) = grad(0,0) - 16.0*gradMat11(0,0);

      Matrix gradMat12 = gaugeMat*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(gaugeMat)
        *covDerivSum*gaugeMat*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(gaugeMat);
      gradMat12 = gradMat12 + conjugateTranspose(gradMat12);
      grad(0,0) = grad(0,0) + 16.0*gradMat12(0,0);

      Matrix fwdBwdCovDerivSumShiftedFwd = getFwdBwdCovDerivSum(field, siteShiftedFwd);
      Matrix fwdBwdCovDerivSumShiftedBwd = getFwdBwdCovDerivSum(field, siteShiftedBwd);

      Matrix gradMat13(2);
      gradMat13 = gradMat13 + 8*conjugateTranspose(gaugeMatShiftedBwd)*fwdBwdCovDerivSumShiftedBwd*gaugeMatShiftedBwd;
      gradMat13 = gradMat13 - 16*fwdBwdCovDerivSum;
      gradMat13 = gradMat13 + 8*gaugeMat*fwdBwdCovDerivSumShiftedFwd*conjugateTranspose(gaugeMat);
      grad(0,0) = grad(0,0) + 2.0*gradMat13(0,0);

      double vevDistanceShiftedFwd = getVevDistance(field, siteShiftedFwd);
      double vevDistanceShiftedBwd = getVevDistance(field, siteShiftedBwd);
      Matrix gradMat14(2);
      gradMat14 = gradMat14 - 32*vevDistanceShiftedBwd*conjugateTranspose(gaugeMatShiftedBwd)*scalarMatShiftedBwd*gaugeMatShiftedBwd;
      gradMat14 = gradMat14 - 32*vevDistanceShiftedFwd*gaugeMat*scalarMatShiftedFwd*conjugateTranspose(gaugeMat);
      grad(0,0) = grad(0,0) + 2.0*0.5*selfCoupling_*gradMat14(0,0);

      }

      double vevDistance = getVevDistance(field, site);

      Matrix gradMat15(2);
      gradMat15 = gradMat15 - 32*real(trace(scalarMat*fwdBwdCovDerivSum))*scalarMat;
      gradMat15 = gradMat15 - 32*real(trace(fwdBwdCovDerivSum*scalarMat))*scalarMat;
      gradMat15 = gradMat15 - 32*vevDistance*fwdBwdCovDerivSum;
      gradMat15 = gradMat15 + 192*vevDistance*scalarMat;
      grad(0,0) = grad(0,0) + 2.0*0.5*selfCoupling_*gradMat15(0,0);

      double trSq = real(trace(scalarMat*scalarMat));
      Matrix gradMat16(2);
      gradMat16 = gradMat16 + 256*pow(0.5*selfCoupling_,2)*vevDistance*trSq*scalarMat;
      gradMat16 = gradMat16 + 128*pow(0.5*selfCoupling_,2)*pow(vevDistance, 2)*scalarMat;
      grad(0,0) = grad(0,0) + 2.0*gradMat16(0,0);

    }


    return grad;
  }

  std::complex<double> GeorgiGlashowSu2EomTheory::getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const
  {
    monsta::Matrix gradMat = getLocalGradient(field, site, matIdx);
    return gradMat(rowIdx, colIdx);
    return 0;
  }

  double GeorgiGlashowSu2EomTheory::getTHooftOperator(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir, std::complex<double> fluxQuanta) const
  {
    double E = 2.0/pow(gaugeCoupling_,2)*real(trace(getPlaquette(field, site, (dir + 1) % 3, (dir + 2) % 3))); // Subtracts unflipped plaquette
    E += 2.0/pow(gaugeCoupling_,2)*real(trace(fluxQuanta*getPlaquette(field, site, (dir + 1) % 3, (dir + 2) % 3))); // Adds flipped plaquette
    return E;
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getTHooftDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, std::complex<double> fluxQuanta) const
  {
    int dir1 = (matIdx + 1) % 3;
    int dir2 = (matIdx + 2) % 3;
    int derivIdx = site.index();

    Matrix tHooftDerivMat(2);

    bool jonnyMorris = false;
    // if (dir2 == 0)
    // {
    //   tHooftDerivMat = tHooftDerivMat - 2*getStaple(field, site, matIdx, dir1, true);
    // }

    if (tHooftLineCheck(field, site, dir2)) {
      tHooftDerivMat = tHooftDerivMat - 2*getStaple(field, site, matIdx, dir1, true);
      // cout << site.coord(0) << " " << site.coord(1) << " " << site.coord(2) << "; " << matIdx << "+" << endl;
      jonnyMorris = true;
    }
    if (tHooftLineCheck(field, site, dir1)) {
      tHooftDerivMat = tHooftDerivMat - 2*getStaple(field, site, matIdx, dir2, true);
      // cout << site.coord(0) << " " << site.coord(1) << " " << site.coord(2) << "; " << matIdx << "+" << endl;
      jonnyMorris = true;
    }

    LATfield2::Site tempSite(site);
    tempSite = tempSite-dir1;
    if (tHooftLineCheck(field, tempSite, dir2)) {
      tHooftDerivMat = tHooftDerivMat - 2*getStaple(field, site, matIdx, dir1, false);
      // cout << site.coord(0) << " " << site.coord(1) << " " << site.coord(2) << "; " << matIdx << "-" << endl;
      jonnyMorris = true;
    }

    tempSite = tempSite+dir1-dir2;
    if (tHooftLineCheck(field, tempSite, dir1)) {
      tHooftDerivMat = tHooftDerivMat - 2*getStaple(field, site, matIdx, dir2, false);
      // cout << site.coord(0) << " " << site.coord(1) << " " << site.coord(2) << "; " << matIdx << "-" << endl;
       jonnyMorris = true;
    }

    return tHooftDerivMat;
  }

    monsta::Matrix GeorgiGlashowSu2EomTheory::getTHooftHinge(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
  {
    int dir1 = (matIdx + 1) % 3;
    int dir2 = (matIdx + 2) % 3;
    int derivIdx = site.index();

    Matrix tHooftHingeMat(2);


    if (tHooftLineCheck(field, site, dir2)) {
      tHooftHingeMat = tHooftHingeMat - 128*getDblPlaquetteDeriv4(field, site, matIdx, dir1, true, dir1, false);
      tHooftHingeMat = tHooftHingeMat + 96*Matrix(field, site, matIdx);
    }
    if (tHooftLineCheck(field, site, dir1)) {
      tHooftHingeMat = tHooftHingeMat - 128*getDblPlaquetteDeriv4(field, site, matIdx, dir2, true, dir2, false);
      tHooftHingeMat = tHooftHingeMat + 96*Matrix(field, site, matIdx);
    }

    LATfield2::Site tempSite(site);
    tempSite = tempSite-dir1;
    if (tHooftLineCheck(field, tempSite, dir2)) {
      tHooftHingeMat = tHooftHingeMat - 128*getDblPlaquetteDeriv4(field, site, matIdx, dir1, false, dir1, true);
      tHooftHingeMat = tHooftHingeMat + 96*Matrix(field, site, matIdx);
    }

    tempSite = tempSite+dir1-dir2;
    if (tHooftLineCheck(field, tempSite, dir1)) {
      tHooftHingeMat = tHooftHingeMat - 128*getDblPlaquetteDeriv4(field, site, matIdx, dir2, false, dir2, true);
      tHooftHingeMat = tHooftHingeMat + 96*Matrix(field, site, matIdx);
    }


    return tHooftHingeMat;
  }

  void GeorgiGlashowSu2EomTheory::postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
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

  void GeorgiGlashowSu2EomTheory::applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const
  {
    // if (boundaryConditions_[0] == 0 && boundaryConditions_[1] == 0 && boundaryConditions_[2] == 0)
    // {
    //   applyDirichletBoundaryConditions(field);
    //   return;
    // }
    field.updateHalo();
    for (int dir = 0; dir < 3; dir++)
    {
      applyCPeriodicBoundaryConditions(field, dir);
    }
  }

  void GeorgiGlashowSu2EomTheory::applyDirichletBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const
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

  void GeorgiGlashowSu2EomTheory::applyCPeriodicBoundaryConditions(LATfield2::Field< std::complex<double> > &field, int dir) const
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

  monsta::Matrix GeorgiGlashowSu2EomTheory::getDirectedLink(LATfield2::Field< std::complex<double> > &field,  LATfield2::Site &site, int dir, bool isFwd) const
  {
    if (isFwd)
    {
      return Matrix(field, site, dir);
    }
    else
    {
      LATfield2::Site reversedSite = site - dir;
      return conjugateTranspose(Matrix(field, reversedSite, dir));
    }
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getCovDerivFwd(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const
  {
    LATfield2::Site tempSite(site);
    tempSite = site + dir;
    Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
    Matrix scalarMatShifted = field(tempSite, 3, 0, 0)*pauli3;
    Matrix gaugeMat(field, site, dir);

    return gaugeMat*scalarMatShifted*conjugateTranspose(gaugeMat) - scalarMat;
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getCovDerivBwd(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const
  {
    LATfield2::Site tempSite(site);
    tempSite = site - dir;
    Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
    Matrix scalarMatShifted = field(tempSite, 3, 0, 0)*pauli3;
    Matrix gaugeMatShifted(field, tempSite, dir);

    return scalarMat - conjugateTranspose(gaugeMatShifted)*scalarMatShifted*gaugeMatShifted;
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getFwdBwdCovDerivSum(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    monsta::Matrix output(2);
    for (int dir = 0; dir < 3; dir++)
    {
      output = output + getCovDerivFwd(field, site, dir);
      output = output - getCovDerivBwd(field, site, dir);
    }
    return output;
  }


  monsta::Matrix  GeorgiGlashowSu2EomTheory::getKineticDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const
  {
    LATfield2::Site tempSite(site);
    tempSite = site + dir;
    Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
    Matrix scalarMatShifted = field(tempSite, 3, 0, 0)*pauli3;
    Matrix gaugeMat(field, site, dir);

    Matrix covDeriv = getCovDerivFwd(field, site, dir);

    return 2*(covDeriv + conjugateTranspose(covDeriv))*gaugeMat*scalarMatShifted;
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getDirectedKineticDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir, bool isFwd) const
  {
    if (isFwd)
    {
      return getKineticDeriv(field, site, dir);
    }
    else
    {
      LATfield2::Site reversedSite = site - dir;
      return conjugateTranspose(getKineticDeriv(field, reversedSite, dir));
    }
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getInsert2(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const
  {
    LATfield2::Site tempSite(site);
    tempSite = site + dir;
    Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
    Matrix scalarMatShifted = field(tempSite, 3, 0, 0)*pauli3;
    Matrix gaugeMat(field, site, dir);

    Matrix covDeriv = getCovDerivFwd(field, site, dir);

    return gaugeMat*conjugateTranspose(scalarMatShifted)*conjugateTranspose(gaugeMat)
      *(covDeriv + conjugateTranspose(covDeriv))*gaugeMat;
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getDirectedInsert2(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir, bool isFwd) const
  {
    if (isFwd)
    {
      return getInsert2(field, site, dir);
    }
    else
    {
      LATfield2::Site reversedSite = site - dir;
      return conjugateTranspose(getInsert2(field, reversedSite, dir));
    }
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getPlaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const
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

  monsta::Matrix GeorgiGlashowSu2EomTheory::getStaple(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, bool isPos2) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix staple = getDirectedLink(field, site, dir2, isPos2);
    // staple.print();

    tempSite = isPos2 ? tempSite + dir2 : tempSite - dir2;
    staple = staple*Matrix(field, tempSite, dir1);
    tempSite = tempSite + dir1;
    tempSite = isPos2 ? tempSite - dir2 : tempSite + dir2;
    staple = staple*conjugateTranspose(getDirectedLink(field, tempSite, dir2, isPos2));

    return staple;

  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getHinge1(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site,
    int derivDir, int dir1, bool isPos1, int dir2, bool isPos2) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix hinge = identity;

    hinge = hinge*getDirectedLink(field, site, dir1, isPos1);
    tempSite = isPos1 ? tempSite + dir1 : tempSite - dir1;
    hinge = hinge*getDirectedLink(field, tempSite, dir2, isPos2);
    tempSite = isPos2 ? tempSite + dir2 : tempSite - dir2;
    hinge = hinge*Matrix(field, tempSite, derivDir);
    tempSite = isPos2 ? tempSite - dir2 + derivDir : tempSite + dir2 + derivDir;
    hinge = hinge*conjugateTranspose(getDirectedLink(field, tempSite, dir2, isPos2));
    tempSite = isPos1 ? tempSite - dir1 : tempSite + dir1;
    hinge = hinge*conjugateTranspose(getDirectedLink(field, tempSite, dir1, isPos1));

    return hinge;
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getHinge2(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site,
    int derivDir, int dir1, bool isPos1, int dir2, bool isPos2) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix hinge = identity;

    hinge = hinge*getDirectedLink(field, site, dir2, isPos2);
    tempSite = isPos2 ? tempSite + dir2 : tempSite - dir2;
    hinge = hinge*Matrix(field, tempSite, derivDir);
    tempSite = tempSite + derivDir;
    hinge = hinge*getDirectedLink(field, tempSite, dir1, isPos1);
    tempSite = isPos1 ? tempSite + dir1 : tempSite - dir1;
    tempSite = isPos2 ? tempSite - dir2 : tempSite + dir2;
    hinge = hinge*conjugateTranspose(getDirectedLink(field, tempSite, dir2, isPos2));
    tempSite = isPos1 ? tempSite - dir1 : tempSite + dir1;
    hinge = hinge*conjugateTranspose(getDirectedLink(field, tempSite, dir1, isPos1));

    return hinge;
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getHinge3(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site,
    int derivDir, int dir1, bool isPos1, int dir2, bool isPos2) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix hinge = identity;

    tempSite = isPos1 ? tempSite - dir1 : tempSite + dir1;
    hinge = hinge*conjugateTranspose(getDirectedLink(field, tempSite, dir1, isPos1));
    hinge = hinge*getDirectedLink(field, tempSite, dir2, isPos2);
    tempSite = isPos2 ? tempSite + dir2 : tempSite - dir2;
    hinge = hinge*getDirectedLink(field, tempSite, dir1, isPos1);
    tempSite = isPos1 ? tempSite + dir1 : tempSite - dir1;
    hinge = hinge*Matrix(field, tempSite, derivDir);
    tempSite = isPos2 ? tempSite + derivDir - dir2 : tempSite + derivDir + dir2;
    hinge = hinge*conjugateTranspose(getDirectedLink(field, tempSite, dir2, isPos2));

    return hinge;
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getDblPlaquetteDeriv1(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site,
    int derivDir, int dir1, bool isPos1, int dir2, bool isPos2) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix derivMat = identity;

    derivMat = derivMat*getDirectedLink(field, tempSite, dir1, isPos1);
    tempSite = isPos1 ? tempSite + dir1 : tempSite - dir1;
    derivMat = derivMat*Matrix(field, tempSite, derivDir);
    tempSite = tempSite + derivDir;
    derivMat = derivMat*getDirectedLink(field, tempSite, dir2, isPos2);
    tempSite = isPos2 ? tempSite + dir2 : tempSite - dir2;
    tempSite = tempSite - derivDir;
    derivMat = derivMat*conjugateTranspose(Matrix(field, tempSite, derivDir));
    tempSite = isPos2 ? tempSite - dir2 : tempSite + dir2;
    derivMat = derivMat*conjugateTranspose(getDirectedLink(field, tempSite, dir2, isPos2));
    derivMat = derivMat*Matrix(field, tempSite, derivDir);
    tempSite = tempSite + derivDir;
    tempSite = isPos1 ? tempSite - dir1 : tempSite + dir1;
    derivMat = derivMat*conjugateTranspose(getDirectedLink(field, tempSite, dir1, isPos1));

    return derivMat;
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getDblPlaquetteDeriv2(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site,
    int derivDir, int dir1, bool isPos1, int dir2, bool isPos2) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix derivMat = identity;

    derivMat = derivMat*getDirectedLink(field, tempSite, dir2, isPos2);
    tempSite = isPos2 ? tempSite + dir2 : tempSite - dir2;
    derivMat = derivMat*Matrix(field, tempSite, derivDir);
    tempSite = tempSite + derivDir;
    tempSite = isPos2 ? tempSite - dir2 : tempSite + dir2;
    derivMat = derivMat*conjugateTranspose(getDirectedLink(field, tempSite, dir2, isPos2));
    derivMat = derivMat*getDirectedLink(field, tempSite, dir1, isPos1);
    tempSite = isPos1 ? tempSite + dir1 : tempSite - dir1;
    derivMat = derivMat*getDirectedLink(field, tempSite, dir2, isPos2);
    tempSite = isPos2 ? tempSite + dir2 : tempSite - dir2;
    tempSite = isPos1 ? tempSite - dir1 : tempSite + dir1;
    derivMat = derivMat*conjugateTranspose(getDirectedLink(field, tempSite, dir1, isPos1));
    tempSite = isPos2 ? tempSite - dir2 : tempSite + dir2;
    derivMat = derivMat*conjugateTranspose(getDirectedLink(field, tempSite, dir2, isPos2));

    return derivMat;
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getDblPlaquetteDeriv3(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site,
    int derivDir, int dir1, bool isPos1, int dir2, bool isPos2) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix derivMat = identity;

    derivMat = derivMat*getDirectedLink(field, tempSite, dir2, isPos2);
    tempSite = isPos2 ? tempSite + dir2 : tempSite - dir2;
    tempSite = isPos1 ? tempSite - dir1 : tempSite + dir1;
    derivMat = derivMat*conjugateTranspose(getDirectedLink(field, tempSite, dir1, isPos1));
    tempSite = isPos2 ? tempSite - dir2 : tempSite + dir2;
    derivMat = derivMat*conjugateTranspose(getDirectedLink(field, tempSite, dir2, isPos2));
    derivMat = derivMat*getDirectedLink(field, tempSite, dir1, isPos1);
    tempSite = isPos1 ? tempSite + dir1 : tempSite - dir1;
    derivMat = derivMat*getDirectedLink(field, tempSite, dir2, isPos2);
    tempSite = isPos2 ? tempSite + dir2 : tempSite - dir2;
    derivMat = derivMat*Matrix(field, tempSite, derivDir);
    tempSite = tempSite + derivDir;
    tempSite = isPos2 ? tempSite - dir2 : tempSite + dir2;
    derivMat = derivMat*conjugateTranspose(getDirectedLink(field, tempSite, dir2, isPos2));

    return derivMat;
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getDblPlaquetteDeriv4(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site,
    int derivDir, int dir1, bool isPos1, int dir2, bool isPos2) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix derivMat = identity;

    derivMat = derivMat*getDirectedLink(field, tempSite, dir1, isPos1);
    tempSite = isPos1 ? tempSite + dir1 : tempSite - dir1;
    derivMat = derivMat*Matrix(field, tempSite, derivDir);
    tempSite = tempSite + derivDir;
    tempSite = isPos1 ? tempSite - dir1 : tempSite + dir1;
    derivMat = derivMat*conjugateTranspose(getDirectedLink(field, tempSite, dir1, isPos1));
    tempSite = tempSite - derivDir;
    derivMat = derivMat*conjugateTranspose(Matrix(field, site, derivDir));
    tempSite = isPos2 ? tempSite - dir2 : tempSite + dir2;
    derivMat = derivMat*conjugateTranspose(getDirectedLink(field, tempSite, dir2, isPos2));
    derivMat = derivMat*Matrix(field, tempSite, derivDir);
    tempSite = tempSite + derivDir;
    derivMat = derivMat*getDirectedLink(field, tempSite, dir2, isPos2);

    return derivMat;
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getKinDerivInsertStaple(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site,
    int dir1, int dir2, bool isPos2, int insertPos) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix insertedStaple = identity;

    if (insertPos == 0)
    {
      insertedStaple = insertedStaple*getDirectedKineticDeriv(field, tempSite, dir2, isPos2);
    }
    else
    {
      insertedStaple = insertedStaple*getDirectedLink(field, tempSite, dir2, isPos2);
    }

    tempSite = isPos2 ? tempSite + dir2 : tempSite - dir2;

    if (insertPos == 1)
    {
      insertedStaple = insertedStaple*getKineticDeriv(field, tempSite, dir1);
    }
    else 
    {
      insertedStaple = insertedStaple*Matrix(field, tempSite, dir1);
    }

    tempSite = tempSite + dir1;
    tempSite = isPos2 ? tempSite - dir2 : tempSite + dir2;

    if (insertPos == 2)
    {
      insertedStaple = insertedStaple*conjugateTranspose(getDirectedKineticDeriv(field, tempSite, dir2, isPos2));
    }
    else
    {
      insertedStaple = insertedStaple*conjugateTranspose(getDirectedLink(field, tempSite, dir2, isPos2));
    }

    return insertedStaple;

  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getInsertStaple2(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site,
    int dir1, int dir2, bool isPos2, int insertPos) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix insertedStaple = identity;

    if (insertPos == 0)
    {
      insertedStaple = insertedStaple*getDirectedInsert2(field, tempSite, dir2, isPos2);
    }
    else
    {
      insertedStaple = insertedStaple*getDirectedLink(field, tempSite, dir2, isPos2);
    }

    tempSite = isPos2 ? tempSite + dir2 : tempSite - dir2;

    if (insertPos == 1)
    {
      insertedStaple = insertedStaple*getInsert2(field, tempSite, dir1);
    }
    else 
    {
      insertedStaple = insertedStaple*Matrix(field, tempSite, dir1);
    }

    tempSite = tempSite + dir1;
    tempSite = isPos2 ? tempSite - dir2 : tempSite + dir2;

    if (insertPos == 2)
    {
      insertedStaple = insertedStaple*conjugateTranspose(getDirectedInsert2(field, tempSite, dir2, isPos2));
    }
    else
    {
      insertedStaple = insertedStaple*conjugateTranspose(getDirectedLink(field, tempSite, dir2, isPos2));
    }

    return insertedStaple;

  }

  double GeorgiGlashowSu2EomTheory::getVevDistance(LATfield2::Field< std::complex<double> > &field, LATfield2::Site site) const
  {
    Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
    return real(trace(scalarMat*scalarMat)) - pow(vev_,2);
  }

  double GeorgiGlashowSu2EomTheory::getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
  {
    int dir1 = (cpt + 1) % 3;
    int dir2 = (cpt + 2) % 3;

    monsta::Matrix plaquette = getPlaquette(field, site, dir1, dir2);
    double magneticField = 2./gaugeCoupling_*arg(plaquette(0,0));
    if (tHooftLineCheck(field, site, cpt))
    {
      magneticField = magneticField > 0 ? magneticField - 2*3.141592654 : magneticField + 2*3.141592654;
    }
    return (magneticField);
  }

  bool GeorgiGlashowSu2EomTheory::tHooftLineCheck(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const
  {
    // if (!tHooftLine_) { return false; }
    // int xCoord = site.coord(0);
    // int yCoord = site.coord(1);
    // int zCoord = site.coord(2);
    // int xSize = field.lattice().size(0);
    // int ySize = field.lattice().size(1);
    // int zSize = field.lattice().size(2);

    // if (dir == 0 && yCoord == ySize/2 && zCoord == zSize/2)
    // {
    //   return true;
    // }

    return false;
  }

}

#endif