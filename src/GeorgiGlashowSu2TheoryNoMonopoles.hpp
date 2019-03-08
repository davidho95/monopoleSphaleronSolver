#ifndef GEORGIGLASHOWSU2THEORY_HPP
#define GEORGIGLASHOWSU2THEORY_HPP

#include "Theory.hpp"
#include "Su2Tools.hpp"

namespace monsta {
  class GeorgiGlashowSu2Theory: public Theory {
  public:
    GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling);
    GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling, bool twistedBoundaryConditions);
    GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling, bool twistedBoundaryConditions, bool tHooftLine);
    double getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    std::complex<double> getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const;
    monsta::Matrix getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;
    monsta::Matrix getPlaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const;
    monsta::Matrix getTrPlaquetteDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, int derivIdx, int derivDir) const;
    monsta::Matrix getU1Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    monsta::Matrix getU1Projector(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    monsta::Matrix getU1Plaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const;
    double getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    void postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;
    void applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const;
    double getTHooftOperator(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir, std::complex<double> fluxQuanta=1) const;
    monsta::Matrix getTHooftDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, std::complex<double> fluxQuanta=1) const;

  private:
    bool tHooftLineCheck(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;

    int numFieldMatrices_ = 4;
    int numRows_ = 2;
    int numCols_ = 2; 
    int matSize_ = 2;
    double gaugeCoupling_;
    double vev_;
    double selfCoupling_;
    bool twistedBoundaryConditions_ = false;
    bool tHooftLine_ = false;
    std::complex<double> fluxQuanta_ = 1;

    bool monopolesPresent_ = false;
    std::vector<int> monopolePos1_;
    std::vector<int> monopolePos2_;
  };

  GeorgiGlashowSu2Theory::GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling)
  : Theory(numFieldMatrices_, matSize_, matSize_), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling)
  {}

  GeorgiGlashowSu2Theory::GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling, bool twistedBoundaryConditions)
  : Theory(numFieldMatrices_, matSize_, matSize_), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling),
    twistedBoundaryConditions_(twistedBoundaryConditions)
  {}

  GeorgiGlashowSu2Theory::GeorgiGlashowSu2Theory(double gaugeCoupling, double vev, double selfCoupling, bool twistedBoundaryConditions, bool tHooftLine)
  : Theory(numFieldMatrices_, matSize_, matSize_), gaugeCoupling_(gaugeCoupling), vev_(vev), selfCoupling_(selfCoupling),
    twistedBoundaryConditions_(twistedBoundaryConditions), tHooftLine_(tHooftLine)
  {}

  double GeorgiGlashowSu2Theory::getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    double E = 0;

    // int numDimensions = 3;

    // Wilson action
    E += 2/gaugeCoupling_*(2 - real(trace(getPlaquette(field, site, 1, 2))));
    E += 2/gaugeCoupling_*(2 - real(trace(getPlaquette(field, site, 2, 0))));
    E += 2/gaugeCoupling_*(2 - real(trace(getPlaquette(field, site, 0, 1))));

    for (int ii = 0; ii < 3; ii++)
    {
      if (tHooftLineCheck(field, site, ii))
      {
        E += getTHooftOperator(field, site, ii, fluxQuanta_);
      }
    }

    // Covariant Derivative
    for (int ii = 0; ii < 3; ii++)
    {
      LATfield2::Site shiftedSite = site + ii;
      Matrix covDeriv = Matrix(field, site, ii)*Matrix(field, shiftedSite, 3)*conjugateTranspose(Matrix(field, site, ii)) - Matrix(field, site, 3);
      E += 0.5*real(trace(covDeriv*covDeriv));
    }

    // Higgs Potential
    double trScalarSq = real(trace(Matrix(field, site, 3)*Matrix(field, site, 3)));
    E += selfCoupling_*pow((trScalarSq - pow(vev_,2)),2);

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
      Matrix plaquetteDerivMat = getTrPlaquetteDeriv(field, site, matIdx, dir1, derivIdx, matIdx);
      plaquetteDerivMat = plaquetteDerivMat + getTrPlaquetteDeriv(field, site, dir2, matIdx, derivIdx, matIdx);
      LATfield2::Site tempSite(site);
      tempSite = tempSite-dir1;
      plaquetteDerivMat = plaquetteDerivMat + getTrPlaquetteDeriv(field, tempSite, matIdx, dir1, derivIdx, matIdx);
      tempSite = tempSite+dir1-dir2;
      plaquetteDerivMat = plaquetteDerivMat + getTrPlaquetteDeriv(field, tempSite, dir2, matIdx, derivIdx, matIdx);

      grad = grad - 2/gaugeCoupling_*plaquetteDerivMat;
      grad = grad - 2/gaugeCoupling_*getTHooftDeriv(field, site, matIdx, fluxQuanta_);

      // Derivative of kinetic term
      tempSite = site + matIdx;
      Matrix scalarMat(field, site, 3);
      Matrix scalarMatShifted(field, tempSite, 3);
      Matrix gaugeMat(field, site, matIdx);
      Matrix kineticDerivMat = gaugeMat*scalarMatShifted*conjugateTranspose(gaugeMat)*gaugeMat*scalarMatShifted;
      kineticDerivMat = kineticDerivMat + gaugeMat*conjugateTranspose(scalarMatShifted)*conjugateTranspose(gaugeMat)*gaugeMat*conjugateTranspose(scalarMatShifted);
      kineticDerivMat = kineticDerivMat - conjugateTranspose(scalarMat)*gaugeMat*conjugateTranspose(scalarMatShifted);
      kineticDerivMat = kineticDerivMat - scalarMat*gaugeMat*scalarMatShifted;
      grad = grad + kineticDerivMat;


    } else {
      // if (abs(trace(Matrix(field, site, 3)*Matrix(field, site, 3))) < 1e-15)
      // {
      //   return grad;
      // }
      for (int dir = 0; dir < 3; dir++)
      {
        // Deriviative of kinetic term
        LATfield2::Site tempSite = site + dir;
        Matrix scalarMat(field, site, 3);
        Matrix scalarMatShiftedFwd(field, tempSite, 3);
        Matrix gaugeMat(field, site, dir);

        tempSite = site - dir;
        Matrix gaugeMatShiftedBwd(field, tempSite, dir);
        Matrix scalarMatShiftedBwd(field, tempSite, 3);

        Matrix kineticDerivMat = conjugateTranspose(gaugeMatShiftedBwd)*gaugeMatShiftedBwd*conjugateTranspose(scalarMat)*conjugateTranspose(gaugeMatShiftedBwd)*gaugeMatShiftedBwd;
        kineticDerivMat = kineticDerivMat - gaugeMat*conjugateTranspose(scalarMatShiftedFwd)*conjugateTranspose(gaugeMat);
        kineticDerivMat = kineticDerivMat - conjugateTranspose(gaugeMatShiftedBwd)*conjugateTranspose(scalarMatShiftedBwd)*gaugeMatShiftedBwd;
        kineticDerivMat = kineticDerivMat + conjugateTranspose(scalarMat);

        grad = grad + kineticDerivMat;
      }

      // Derivative of Higgs Potential
      double trscalarSq = real(trace(Matrix(field, site, matIdx)*Matrix(field, site, matIdx)));
      grad = grad + 4.0*selfCoupling_*(trscalarSq - pow(vev_,2))*conjugateTranspose(Matrix(field, site, matIdx));
    }

    return grad;
  }

  std::complex<double> GeorgiGlashowSu2Theory::getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const
  {
    monsta::Matrix gradMat = getLocalGradient(field, site, matIdx);
    return gradMat(rowIdx, colIdx);
    return 0;
  }

  double GeorgiGlashowSu2Theory::getTHooftOperator(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir, std::complex<double> fluxQuanta) const
  {
    double E = 2/gaugeCoupling_*real(trace(getPlaquette(field, site, (dir + 1) % 3, (dir + 2) % 3)));
    E += 2/gaugeCoupling_*real(trace(fluxQuanta*getPlaquette(field, site, (dir + 1) % 3, (dir + 2) % 3)));
    return E;
  }

  monsta::Matrix GeorgiGlashowSu2Theory::getTHooftDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, std::complex<double> fluxQuanta) const
  {
    int dir1 = (matIdx + 1) % 3;
    int dir2 = (matIdx + 2) % 3;
    int derivIdx = site.index();

    bool blah = false;

    Matrix tHooftDerivMat(2);
    if (tHooftLineCheck(field, site, dir2)) {
      tHooftDerivMat = tHooftDerivMat - getTrPlaquetteDeriv(field, site, matIdx, dir1, derivIdx, matIdx);
      tHooftDerivMat = tHooftDerivMat - fluxQuanta*getTrPlaquetteDeriv(field, site, matIdx, dir1, derivIdx, matIdx);
      blah = true;
    }
    if (tHooftLineCheck(field, site, dir1)) {
      tHooftDerivMat = tHooftDerivMat - getTrPlaquetteDeriv(field, site, dir2, matIdx, derivIdx, matIdx);
      tHooftDerivMat = tHooftDerivMat - fluxQuanta*getTrPlaquetteDeriv(field, site, dir2, matIdx, derivIdx, matIdx);
      // blah = true;
    }

    LATfield2::Site tempSite(site);
    tempSite = tempSite-dir1;
    if (tHooftLineCheck(field, tempSite, dir2)) {
      tHooftDerivMat = tHooftDerivMat - getTrPlaquetteDeriv(field, tempSite, matIdx, dir1, derivIdx, matIdx);
      tHooftDerivMat = tHooftDerivMat - fluxQuanta*getTrPlaquetteDeriv(field, tempSite, matIdx, dir1, derivIdx, matIdx);
      // blah = true;
    }

    tempSite = tempSite+dir1-dir2;
    if (tHooftLineCheck(field, tempSite, dir1)) {
      tHooftDerivMat = tHooftDerivMat - getTrPlaquetteDeriv(field, tempSite, dir2, matIdx, derivIdx, matIdx);
      tHooftDerivMat = tHooftDerivMat - fluxQuanta*getTrPlaquetteDeriv(field, tempSite, dir2, matIdx, derivIdx, matIdx);
      blah = true;
    }

    if (blah)
    {
      std::complex<double> tmp = tHooftDerivMat(0,0);
      tHooftDerivMat(0,0) = conj(tHooftDerivMat(1,1));
      tHooftDerivMat(1,1) = conj(tmp);
    }

    return tHooftDerivMat;
  }

  void GeorgiGlashowSu2Theory::postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
  {
    if (matIdx < 3)
    {
      std::complex<double> determinant = field(site, matIdx, 0, 0)*field(site, matIdx, 1, 1) - field(site, matIdx, 0, 1)*field(site, matIdx, 1, 0);
      for (int rowIdx = 0; rowIdx < numRows_; rowIdx++)
      {
        for (int colIdx = 0; colIdx < numCols_; colIdx++)
        {
          field(site, matIdx, rowIdx, colIdx) = field(site, matIdx, rowIdx, colIdx) / sqrt(determinant);
        }
        // field(site, matIdx, 0, 0) = 0.5*(field(site, matIdx, 0, 0) + conj(field(site, matIdx, 1, 1)));
        field(site, matIdx, 0, 0) = conj(field(site, matIdx, 1, 1));
        // field(site, matIdx, 0, 1) = 0.5*(field(site, matIdx, 0, 1) - field(site, matIdx, 1, 0));
        field(site, matIdx, 1, 0) = -1.0*conj(field(site, matIdx, 0, 1));
      }
    }
    else
    {
      if (twistedBoundaryConditions_ && site.coord(0) == 7 && site.coord(1) == 7 & site.coord(2) == 7)
      {
        // field(site, 3, 0, 0) = 0;
        // field(site, 3, 0, 1) = 0;
        // field(site, 3, 1, 0) = 0;
        // field(site, 3, 1, 1) = 0;
      }
      field(site, 3, 0, 0) = real(field(site, 3, 0, 0));
      field(site, 3, 1, 1) = real(field(site, 3, 1, 1));
    }
  }

  void GeorgiGlashowSu2Theory::applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const
  {
    field.updateHalo();
    monsta::Matrix xMat = identity;
    monsta::Matrix yMat = identity;
    monsta::Matrix zMat = identity;
    if (!twistedBoundaryConditions_)
    {
      return;
      xMat = identity;
      yMat = pauli2;
      zMat = pauli3;
    } else 
    {
      xMat = pauli1;
      yMat = pauli2;
      zMat = pauli3;
    }

    LATfield2::Site site(field.lattice());

    int count = 0;
    for (site.haloFirst(); site.haloTest(); site.haloNext())
    {
      int xCoord = site.coord(0);
      int yCoord = site.coord(1);
      int zCoord = site.coord(2);
      int xMax = field.lattice().size(0) - 1;
      int yMax = field.lattice().size(1) - 1;
      int zMax = field.lattice().size(2) - 1;

      if (xCoord > 0 && yCoord > 0 && zCoord > 0 && xCoord < xMax && yCoord < yMax && zCoord < zMax) // Skip sites in local halo but not on boundary
      {
        continue;
      }
      if (xCoord < 0 || xCoord > xMax)
      {
        for (int matIdx = 0; matIdx < numFieldMatrices_; matIdx++) // Apply C-periodic boundary conditions
        {
          Matrix boundaryMat(field, site, matIdx);
          boundaryMat = xMat*boundaryMat*xMat;
          if (matIdx == 3 && twistedBoundaryConditions_)
          {
            boundaryMat = -1*boundaryMat;
          }
          field(site, matIdx, 0, 0) = boundaryMat(0, 0);
          field(site, matIdx, 0, 1) = boundaryMat(0, 1);
          field(site, matIdx, 1, 0) = boundaryMat(1, 0);
          field(site, matIdx, 1, 1) = boundaryMat(1, 1);
        }
      }
      if (yCoord < 0 || yCoord > yMax)
      {
        for (int matIdx = 0; matIdx < numFieldMatrices_; matIdx++) // Apply C-periodic boundary conditions
        {
          Matrix boundaryMat(field, site, matIdx);
          boundaryMat = yMat*boundaryMat*yMat;
          if (matIdx == 3)
          {
            boundaryMat = -1*boundaryMat;
          }
          field(site, matIdx, 0, 0) = boundaryMat(0, 0);
          field(site, matIdx, 0, 1) = boundaryMat(0, 1);
          field(site, matIdx, 1, 0) = boundaryMat(1, 0);
          field(site, matIdx, 1, 1) = boundaryMat(1, 1);
        }
      }
      if (zCoord < 0 || zCoord > zMax)
      {
        for (int matIdx = 0; matIdx < numFieldMatrices_; matIdx++) // Apply C-periodic boundary conditions
        {
          Matrix boundaryMat(field, site, matIdx);
          boundaryMat = zMat*boundaryMat*zMat;
          if (matIdx == 3)
          {
            boundaryMat = -1*boundaryMat;
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

  monsta::Matrix GeorgiGlashowSu2Theory::getTrPlaquetteDeriv(
    LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, int derivIdx, int derivDir) const
  {

    LATfield2::Site siteX(site);
    siteX = site + dir1;
    LATfield2::Site siteY(site);
    siteY = site + dir2;

    Matrix plaquetteDeriv(2);

    if (derivIdx == site.index() && derivDir == dir1)
    {
      plaquetteDeriv = Matrix(field, siteX, dir2)*conjugateTranspose(Matrix(field, siteY, dir1))*conjugateTranspose(Matrix(field, site, dir2));
      plaquetteDeriv = conjugateTranspose(plaquetteDeriv);
    }
    if (derivIdx == siteX.index() && derivDir == dir2)
    {
      plaquetteDeriv = conjugateTranspose(Matrix(field, siteY, dir1))*conjugateTranspose(Matrix(field, site, dir2))*Matrix(field, site, dir1);
      plaquetteDeriv = conjugateTranspose(plaquetteDeriv);
    }
    if (derivIdx == siteY.index() && derivDir == dir1)
    {
      plaquetteDeriv = conjugateTranspose(Matrix(field, site, dir2))*Matrix(field, site, dir1)*Matrix(field, siteX, dir2);
    }
    if (derivIdx == site.index() && derivDir == dir2)
    {
      plaquetteDeriv = Matrix(field, site, dir1)*Matrix(field, siteX, dir2)*conjugateTranspose(Matrix(field, siteY, dir1));
    }

    return plaquetteDeriv;
  }

  monsta::Matrix GeorgiGlashowSu2Theory::getU1Projector(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    monsta::Matrix scalarMat(field, site, 3);
    double trScalarSq = real(monsta::trace(scalarMat*scalarMat));

    double zeroTol = 1e-15;
    if (trScalarSq < zeroTol)
    {
      return 0.5*monsta::identity;
    }

    return 0.5*(monsta::identity + sqrt(2./trScalarSq)*scalarMat);
  }

  monsta::Matrix GeorgiGlashowSu2Theory::getU1Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
  {
    monsta::Matrix gaugeMat(field, site, cpt);
    monsta::Matrix projector = getU1Projector(field, site);

    LATfield2::Site shiftedSite = site+cpt;
    monsta::Matrix projectorShifted = getU1Projector(field, shiftedSite);

    return projector*gaugeMat*projectorShifted;
  }

  monsta::Matrix GeorgiGlashowSu2Theory::getU1Plaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix u1Plaquette = getU1Link(field, site, dir1);
    tempSite = tempSite+dir1;
    u1Plaquette = u1Plaquette*getU1Link(field, tempSite, dir2);
    tempSite = tempSite-dir1+dir2;
    u1Plaquette = u1Plaquette*conjugateTranspose(getU1Link(field, tempSite, dir1));
    tempSite = tempSite-dir2;
    u1Plaquette = u1Plaquette*conjugateTranspose(getU1Link(field, tempSite, dir2));

    return u1Plaquette;
  }

  double GeorgiGlashowSu2Theory::getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
  {
    int dir1 = (cpt + 1) % 3;
    int dir2 = (cpt + 2) % 3;

    double magneticField = 2./gaugeCoupling_*arg(trace(getU1Plaquette(field, site, dir1, dir2)));
    // if (tHooftLineCheck(field, site, cpt))
    // {
    //   magneticField = magneticField > 0 ? magneticField - 2*3.141592654 : magneticField + 2*3.141592654;
    // }
    return (magneticField);
  }

  bool GeorgiGlashowSu2Theory::tHooftLineCheck(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const
  {
    if (!tHooftLine_) { return false; }
    int xCoord = site.coord(0);
    int yCoord = site.coord(1);
    int zCoord = site.coord(2);
    int xSize = field.lattice().size(0);
    int ySize = field.lattice().size(1);
    int zSize = field.lattice().size(2);

    if (dir == 0 && xCoord > 2 && xCoord <= 6 && yCoord == ySize/2 - 2 && zCoord == zSize/2)
    {
      return true;
    }
    if (dir == 0 && xCoord > 2 && xCoord <= 6 && yCoord == ySize/2 + 2 && zCoord == zSize/2)
    {
      return true;
    }
    if (dir == 1 && xCoord == 2 && yCoord > ySize/2 - 2 && yCoord <= ySize/2 + 2 && zCoord == zSize/2)
    {
      return true;
    }
    if (dir == 1 && xCoord == 6 && yCoord > ySize/2 - 2 && yCoord <= ySize/2 + 2 && zCoord == zSize/2)
    {
      return true;
    }

    if (dir == 0 && xCoord > 2 && xCoord <= 6 && zCoord == zSize/2 - 2 && yCoord == ySize/2)
    {
      return true;
    }
    if (dir == 0 && xCoord > 2 && xCoord <= 6 && zCoord == zSize/2 + 2 && yCoord == ySize/2)
    {
      return true;
    }
    if (dir == 2 && xCoord == 2 && zCoord > zSize/2 - 2 && zCoord <= zSize/2 + 2 && yCoord == ySize/2)
    {
      return true;
    }
    if (dir == 2 && xCoord == 6 && zCoord > zSize/2 - 2 && zCoord <= zSize/2 + 2 && yCoord == ySize/2)
    {
      return true;
    }

    // if (dir == 0 && xCoord > 2 && xCoord <= 6 && zCoord == zSize/2 && yCoord == ySize/2)
    // {
    //   return true;
    // }

    // if (dir == 1 && xCoord == xSize/2 && zCoord == zSize/2)
    // {
    //   return true;
    // }

    return false;
  }

}

#endif