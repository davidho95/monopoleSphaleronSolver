#ifndef GEORGIGLASHOWSU2EOMTHEORY_HPP
#define GEORGIGLASHOWSU2EOMTHEORY_HPP

#include "Theory.hpp"
#include "Su2Tools.hpp"

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
    monsta::Matrix getCovDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    monsta::Matrix getKineticDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    monsta::Matrix getDirectedKineticDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
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
    monsta::Matrix getDirectedLink(LATfield2::Field< std::complex<double> > &field,  LATfield2::Site &site, int dir, bool isFwd) const;
    monsta::Matrix getDirectedKineticDeriv(LATfield2::Field< std::complex<double> > &field,  LATfield2::Site &site, int dir, bool isFwd) const;

  private:
    bool tHooftLineCheck(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    void applyDirichletBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const;
    void applyCPeriodicBoundaryConditions(LATfield2::Field< std::complex<double> > &field, int dir) const;
    void applyTwistedBoundaryConditions(LATfield2::Field< std::complex<double> > &field, int dir) const;
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

    for (int matIdx = 0; matIdx < 4; matIdx++)
    {
      monsta::Matrix grad(2);

      if (matIdx < 3)
      {
        int dir1 = (matIdx + 1) % 3;
        int dir2 = (matIdx + 2) % 3;

        // E -= real(trace(plaqMat));

        // Derivative of Wilson action
        Matrix plaquetteDerivMat(2);
        plaquetteDerivMat = plaquetteDerivMat + getStaple(field, site, matIdx, dir1, true);
        plaquetteDerivMat = plaquetteDerivMat + getStaple(field, site, matIdx, dir1, false);
        plaquetteDerivMat = plaquetteDerivMat + getStaple(field, site, matIdx, dir2, true);
        plaquetteDerivMat = plaquetteDerivMat + getStaple(field, site, matIdx, dir2, false);

        plaquetteDerivMat = -2.0/pow(gaugeCoupling_,2)*plaquetteDerivMat;
        grad = grad + plaquetteDerivMat;
        // grad = grad - 2.0/pow(gaugeCoupling_,2)*getTHooftDeriv(field, site, matIdx, fluxQuanta_);

        // Derivative of kinetic term
        Matrix kineticDerivMat = getKineticDeriv(field, site, matIdx);
        grad = grad + kineticDerivMat;

        Matrix plaqGradProj = plaquetteDerivMat*conjugateTranspose(Matrix(field, site, matIdx));
        Matrix kineticDerivProj = kineticDerivMat*conjugateTranspose(Matrix(field, site, matIdx));

        // cout << real(trace(plaqGradProj*kineticDerivProj + kineticDerivProj*plaqGradProj)) << endl;
        E -= real(trace(plaqGradProj*kineticDerivProj + kineticDerivProj*plaqGradProj));
        // E += real(trace(grad*conjugateTranspose(kineticDerivMat) + kineticDerivMat*conjugateTranspose(grad)));
        // E += real(trace(kineticDerivMat*conjugateTranspose(kineticDerivMat)));

        // E += real(trace(grad*conjugateTranspose(grad)));
        // monsta::Matrix gradProj = grad*conjugateTranspose(Matrix(field, site, matIdx));
        // E -= real(trace(gradProj*gradProj));

      } else {
        // if (abs(field(site, 3, 0, 0)) < 1e-15) { return grad; }
        for (int dir = 0; dir < 3; dir++)
        {
          // // Deriviative of kinetic term
          // LATfield2::Site tempSite = site + dir;
          // Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
          // Matrix scalarMatShiftedFwd = field(tempSite, 3, 0, 0)*pauli3;
          // Matrix gaugeMat(field, site, dir);

          // tempSite = site - dir;
          // Matrix gaugeMatShiftedBwd(field, tempSite, dir);
          // Matrix scalarMatShiftedBwd = field(tempSite, 3, 0, 0)*pauli3;

          // Matrix kineticDerivMat = 2*conjugateTranspose(gaugeMatShiftedBwd)*gaugeMatShiftedBwd*conjugateTranspose(scalarMat)*conjugateTranspose(gaugeMatShiftedBwd)*gaugeMatShiftedBwd;
          // kineticDerivMat = kineticDerivMat - 2*gaugeMat*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(gaugeMat);
          // kineticDerivMat = kineticDerivMat - 2*conjugateTranspose(gaugeMatShiftedBwd)*conjugateTranspose(scalarMatShiftedBwd)*gaugeMatShiftedBwd;
          // kineticDerivMat = kineticDerivMat + 2*conjugateTranspose(scalarMat);

          // grad(0,0) = grad(0,0) + 2.*kineticDerivMat(0,0);
        }

        // Derivative of Higgs Potential
        // grad(0,0) = grad(0,0) + 8.0*selfCoupling_*field(site, 3, 0, 0)*(2.0*pow(field(site, 3, 0, 0),2) - pow(vev_, 2));

        // E += pow(real(grad(0,0)),2);
      }
    }

    cout << E << endl;

    return E;
  }

  monsta::Matrix GeorgiGlashowSu2EomTheory::getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
  {
    monsta::Matrix grad(2);

    if (matIdx < 3)
    {
      // for (int a = 0; a < 3; a++)
      // {
      //   if (a == matIdx) { continue; }
      //   for (bool sgnA : {false, true})
      //   {
      //     for (int b = 0; b < 3; b++)
      //     {
      //       if (b == matIdx) { continue; }
      //       for (bool sgnB : {false, true})
      //       {
      //         grad = grad + 4*2.0/pow(gaugeCoupling_,2)*getHinge1(field, site, matIdx, a, sgnA, b, sgnB);
      //         grad = grad - 4*2.0/pow(gaugeCoupling_,2)*getDblPlaquetteDeriv1(field, site, matIdx, a, sgnA, b, sgnB);
      //         grad = grad - 4*2.0/pow(gaugeCoupling_,2)*getDblPlaquetteDeriv4(field, site, matIdx, a, sgnA, b, sgnB);
      //       }
      //     }
      //   }
      // }

      // for (int a = 0; a < 3; a++)
      // {
      //   for (bool sgnA : {false, true})
      //   {
      //     for (int b = 0; b < 3; b++)
      //     {
      //       if (b == matIdx) { continue; }
      //       if (b == a) { continue; }
      //       for (bool sgnB : {false, true})
      //       {
      //         grad = grad + 4*2.0/pow(gaugeCoupling_,2)*getHinge2(field, site, matIdx, a, sgnA, b, sgnB);
      //         grad = grad + 4*2.0/pow(gaugeCoupling_,2)*getHinge3(field, site, matIdx, a, sgnA, b, sgnB);
      //         grad = grad - 4*2.0/pow(gaugeCoupling_,2)*getDblPlaquetteDeriv2(field, site, matIdx, a, sgnA, b, sgnB);
      //         grad = grad - 4*2.0/pow(gaugeCoupling_,2)*getDblPlaquetteDeriv3(field, site, matIdx, a, sgnA, b, sgnB);
      //       }
      //     }
      //   }
      // }

      // LATfield2::Site siteShifted = site + matIdx;
      // Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
      // Matrix scalarMatShifted = field(siteShifted, 3, 0, 0)*pauli3;
      // Matrix gaugeMat(field, site, matIdx);
      // Matrix covDerivSum = getCovDeriv(field, site, matIdx) + conjugateTranspose(getCovDeriv(field, site, matIdx));

      // for (int ii = 0; ii < 3; ii++)
      // {
      //   if (ii == matIdx) { continue; }
      //   for (int stapleLinkIdx = 0; stapleLinkIdx < 3; stapleLinkIdx++)
      //   {
      //     for (bool stapleSgn : {true, false})
      //     {
      //       grad = grad - 2*2.0/pow(gaugeCoupling_,2)*getKinDerivInsertStaple(field, site, matIdx, ii, stapleSgn, stapleLinkIdx);
      //     }
      //   }

      //   Matrix stapleSum = getStaple(field, site, matIdx, ii, true) + getStaple(field, site, matIdx, ii, false);
      //   Matrix stapleSumPlus = stapleSum*conjugateTranspose(scalarMatShifted)*conjugateTranspose(gaugeMat);

      //   grad = grad - 4*2.0/pow(gaugeCoupling_,2)*(stapleSumPlus + conjugateTranspose(stapleSumPlus))*gaugeMat*(scalarMatShifted + conjugateTranspose(scalarMatShifted));

      //   grad = grad - 4*2.0/pow(gaugeCoupling_,2)*covDerivSum*stapleSum*conjugateTranspose(scalarMatShifted);
      // }

      // Matrix kineticDerivMat = getKineticDeriv(field, site, matIdx);
      // grad = grad + 8*2.0/pow(gaugeCoupling_,2)
      //   *kineticDerivMat*conjugateTranspose(scalarMatShifted)*scalarMatShifted;//*conjugateTranspose(gaugeMat)*gaugeMat
      //   // *(scalarMatShifted + conjugateTranspose(scalarMatShifted));
      // // grad = grad + 4*2.0/pow(gaugeCoupling_,2)
      //   // *gaugeMat*scalarMatShifted*conjugateTranspose(scalarMatShifted)*conjugateTranspose(gaugeMat)*covDerivSum*gaugeMat
      //   // *(scalarMatShifted + conjugateTranspose(scalarMatShifted));
      // grad = grad + 4*2.0/pow(gaugeCoupling_,2)
      //   *covDerivSum*covDerivSum*gaugeMat*scalarMatShifted*conjugateTranspose(scalarMatShifted);
    }
    else
    {
      for (int ii = 0; ii < 3; ii++)
      {
      // for (int jj = 0; jj < 3; jj++)
      //   {
      //     if (jj == ii) { continue; }
      //     LATfield2::Site siteShiftedFwd = site + ii;
      //     LATfield2::Site siteShiftedBwd = site - ii;
      //     LATfield2::Site siteShiftedUp = site + jj;
      //     LATfield2::Site siteShiftedDown = site - jj;

      //     Matrix gaugeMatFwd(field, site, ii);
      //     Matrix gaugeMatShiftedBwd(field, siteShiftedBwd, ii);

      //     Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
      //     Matrix scalarMatShiftedFwd = field(siteShiftedFwd, 3, 0, 0)*pauli3;

      //     Matrix stapleSumShiftedBwd = getStaple(field, siteShiftedBwd, ii, jj, true) + getStaple(field, siteShiftedBwd, ii, jj, false);
      //     Matrix gradMat1 = conjugateTranspose(gaugeMatShiftedBwd)*stapleSumShiftedBwd*conjugateTranspose(scalarMat)*conjugateTranspose(gaugeMatShiftedBwd)*gaugeMatShiftedBwd;
      //     gradMat1 = gradMat1 + conjugateTranspose(gradMat1);
      //     grad(0,0) = grad(0,0) - 8*2.0/pow(gaugeCoupling_,2)*gradMat1(0,0);

      //     Matrix kinDerivSumBwd = getCovDeriv(field, siteShiftedBwd, ii) + conjugateTranspose(getCovDeriv(field, siteShiftedBwd, ii));
      //     Matrix gradMat2 = conjugateTranspose(gaugeMatShiftedBwd)*kinDerivSumBwd*stapleSumShiftedBwd;
      //     grad(0,0) = grad(0,0) - 8*2.0/pow(gaugeCoupling_,2)*real(gradMat2(0,0));

      //     Matrix stapleSum = getStaple(field, site, ii, jj, true) + getStaple(field, site, ii, jj, false);
      //     Matrix gradMat3 = stapleSum*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(gaugeMatFwd);
      //     gradMat3 = gradMat3 + conjugateTranspose(gradMat3);
      //     grad(0,0) = grad(0,0) + 8*2.0/pow(gaugeCoupling_,2)*gradMat3(0,0);

      //   }

      // LATfield2::Site siteShiftedFwd = site + ii;
      // LATfield2::Site siteShiftedBwd = site - ii;
      // Matrix gaugeMat(field, site, ii);
      // Matrix gaugeMatShiftedBwd(field, siteShiftedBwd, ii);
      // Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
      // Matrix scalarMatShiftedFwd = field(siteShiftedFwd, 3, 0, 0)*pauli3;
      // Matrix scalarMatShiftedBwd = field(siteShiftedBwd, 3, 0, 0)*pauli3;

      // Matrix covDerivSumShiftedBwd = getCovDeriv(field, siteShiftedBwd, ii) + conjugateTranspose(getCovDeriv(field, siteShiftedBwd, ii));
      // Matrix covDerivSum = getCovDeriv(field, site, ii) + conjugateTranspose(getCovDeriv(field, site, ii));

      // Matrix gradMat4 = conjugateTranspose(gaugeMatShiftedBwd)*covDerivSumShiftedBwd*gaugeMatShiftedBwd*scalarMat
      //   *conjugateTranspose(scalarMat);//*conjugateTranspose(gaugeMatShiftedBwd)*gaugeMatShiftedBwd;
      // gradMat4 = gradMat4 + conjugateTranspose(gradMat4);
      // grad(0,0) = grad(0,0) + 16*gradMat4(0,0);

      // Matrix gradMat5 = conjugateTranspose(gaugeMatShiftedBwd)*covDerivSumShiftedBwd*covDerivSumShiftedBwd*gaugeMatShiftedBwd*scalarMat;
      // grad(0,0) = grad(0,0) + 16*gradMat5(0,0);

      // Matrix gradMat6 = getKineticDeriv(field, site, ii)*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(gaugeMat);
      // // Matrix gradMat6 = covDerivSum*gaugeMat*scalarMatShiftedFwd*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(gaugeMat);
      // gradMat6 = gradMat6 + conjugateTranspose(gradMat6);
      // grad(0,0) = grad(0,0) - 8*gradMat6(0,0);

      }
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

    if (tHooftLineCheck(field, site, dir2)) {
      tHooftDerivMat = tHooftDerivMat - 2*getStaple(field, site, matIdx, dir1, true);
    }
    if (tHooftLineCheck(field, site, dir1)) {
      tHooftDerivMat = tHooftDerivMat - 2*getStaple(field, site, matIdx, dir2, true);
    }

    LATfield2::Site tempSite(site);
    tempSite = tempSite-dir1;
    if (tHooftLineCheck(field, tempSite, dir2)) {
      tHooftDerivMat = tHooftDerivMat - 2*getStaple(field, site, matIdx, dir1, false);
    }

    tempSite = tempSite+dir1-dir2;
    if (tHooftLineCheck(field, tempSite, dir1)) {
      tHooftDerivMat = tHooftDerivMat - 2*getStaple(field, site, matIdx, dir2, false);
    }

    return tHooftDerivMat;
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
    if (boundaryConditions_[0] == 0 && boundaryConditions_[1] == 0 && boundaryConditions_[2] == 0)
    {
      applyDirichletBoundaryConditions(field);
      return;
    }
    field.updateHalo();
    field.updateHalo();
    // if (boundaryConditions_[0] == -1 && boundaryConditions_[1] == -1 && boundaryConditions_[2] == -1)
    // {
    //   for (int dir = 0; dir < 3; dir++)
    //   {
    //     applyTwistedBoundaryConditions(field, dir);
    //     return;
    //   }
    // }
    for (int dir = 0; dir < 3; dir++)
    {
      if (boundaryConditions_[dir] == -1)
      {
        applyCPeriodicBoundaryConditions(field, dir);
      }
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
          Matrix boundaryMat(field, site, matIdx);
          if (matIdx < 3)
          {
            boundaryMat = pauli2*boundaryMat*pauli2;
          }
          field(site, matIdx, 0, 0) = boundaryMat(0, 0);
          field(site, matIdx, 0, 1) = boundaryMat(0, 1);
          field(site, matIdx, 1, 0) = boundaryMat(1, 0);
          field(site, matIdx, 1, 1) = boundaryMat(1, 1);
        }
      }
    }
  }

  void GeorgiGlashowSu2EomTheory::applyTwistedBoundaryConditions(LATfield2::Field< std::complex<double> > &field, int dir) const
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
          Matrix boundaryMat(field, site, matIdx);
          if (matIdx < 3)
          {
            switch(dir)
            {
              case 0:
                boundaryMat = pauli1*boundaryMat*pauli1;
                break;
              case 1:
                boundaryMat = pauli2*boundaryMat*pauli2;
                break;
              case 2:
                boundaryMat = pauli3*boundaryMat*pauli3;
                break;
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

  monsta::Matrix GeorgiGlashowSu2EomTheory::getCovDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const
  {
    LATfield2::Site tempSite(site);
    tempSite = site + dir;
    Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
    Matrix scalarMatShifted = field(tempSite, 3, 0, 0)*pauli3;
    Matrix gaugeMat(field, site, dir);

    return gaugeMat*scalarMatShifted*conjugateTranspose(gaugeMat) - scalarMat;
  }

  monsta::Matrix  GeorgiGlashowSu2EomTheory::getKineticDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const
  {
    LATfield2::Site tempSite(site);
    tempSite = site + dir;
    Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
    Matrix scalarMatShifted = field(tempSite, 3, 0, 0)*pauli3;
    Matrix gaugeMat(field, site, dir);

    Matrix covDeriv = getCovDeriv(field, site, dir);

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
    if (!tHooftLine_) { return false; }
    int xCoord = site.coord(0);
    int yCoord = site.coord(1);
    int zCoord = site.coord(2);
    int xSize = field.lattice().size(0);
    int ySize = field.lattice().size(1);
    int zSize = field.lattice().size(2);

    if (dir == 0 && yCoord == ySize/2 && zCoord == zSize/2)
    {
      return true;
    }
    // if (dir == 1 && xCoord == xSize/2 && zCoord == zSize/2)
    // {
    //   return true;
    // }
    // if (dir == 2 && yCoord == ySize/2 && xCoord == xSize/2)
    // {
    //   return true;
    // }


    return false;
  }

}

#endif