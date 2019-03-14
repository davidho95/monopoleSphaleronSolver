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
    std::complex<double> getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const;
    monsta::Matrix getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;
    double getVev() const { return vev_; }
    double getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    void postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;
    void applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const;

  private:
    bool tHooftLineCheck(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    void applyDirichletBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const;
    void applyCPeriodicBoundaryConditions(LATfield2::Field< std::complex<double> > &field, int dir) const;
    monsta::Matrix getPlaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const;
    monsta::Matrix getTrPlaquetteDeriv(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, int derivIdx, int derivDir) const;
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
      Matrix covDeriv = field(shiftedSite, 3, 0, 0)*Matrix(field, site, ii)*pauli3*conjugateTranspose(Matrix(field, site, ii)) - field(site, 3, 0, 0)*pauli3;
      E += 0.5*real(trace(covDeriv*covDeriv));
    }

    // Higgs Potential
    E += real(selfCoupling_*pow(pow(field(site, 3, 0, 0),2) - pow(vev_, 2),2));

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
      Matrix scalarMat = field(site, 3, 0, 0)*pauli3;
      Matrix scalarMatShifted = field(tempSite, 3, 0, 0)*pauli3;
      Matrix gaugeMat(field, site, matIdx);
      Matrix kineticDerivMat = gaugeMat*scalarMatShifted*conjugateTranspose(gaugeMat)*gaugeMat*scalarMatShifted;
      kineticDerivMat = kineticDerivMat + gaugeMat*conjugateTranspose(scalarMatShifted)*conjugateTranspose(gaugeMat)*gaugeMat*conjugateTranspose(scalarMatShifted);
      kineticDerivMat = kineticDerivMat - conjugateTranspose(scalarMat)*gaugeMat*conjugateTranspose(scalarMatShifted);
      kineticDerivMat = kineticDerivMat - scalarMat*gaugeMat*scalarMatShifted;
      grad = grad + kineticDerivMat;


    } else {
      if (abs(field(site, 3, 0, 0)) < 1e-15) { return grad; }
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

        grad(0,0) = grad(0,0) + kineticDerivMat(0,0);
      }

      // Derivative of Higgs Potential
      grad(0,0) = grad(0,0) + 4.0*selfCoupling_*field(site, 3, 0, 0)*(pow(field(site, 3, 0, 0),2) - pow(vev_, 2));
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
    double E = 2/gaugeCoupling_*real(trace(getPlaquette(field, site, (dir + 1) % 3, (dir + 2) % 3))); // Subtracts unflipped plaquette
    E += 2/gaugeCoupling_*real(trace(fluxQuanta*getPlaquette(field, site, (dir + 1) % 3, (dir + 2) % 3))); // Adds flipped plaquette
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
    if (boundaryConditions_[0] == 0 && boundaryConditions_[1] == 0 && boundaryConditions_[2] == 0)
    {
      applyDirichletBoundaryConditions(field);
      return;
    }
    field.updateHalo();
    for (int dir = 0; dir < 3; dir++)
    {
      if (boundaryConditions_[dir] == -1)
      {
        applyCPeriodicBoundaryConditions(field, dir);
      }
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

  double GeorgiGlashowSu2Theory::getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
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

  bool GeorgiGlashowSu2Theory::tHooftLineCheck(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const
  {
    if (!tHooftLine_) { return false; }
    int xCoord = site.coord(0);
    int yCoord = site.coord(1);
    int zCoord = site.coord(2);
    int xSize = field.lattice().size(0);
    int ySize = field.lattice().size(1);
    int zSize = field.lattice().size(2);

    // if (dir == 0 && xCoord > 2 && xCoord <= 6 && yCoord == ySize/2 - 2 && zCoord == zSize/2)
    // {
    //   return true;
    // }
    // if (dir == 0 && xCoord > 2 && xCoord <= 6 && yCoord == ySize/2 + 2 && zCoord == zSize/2)
    // {
    //   return true;
    // }
    // if (dir == 1 && xCoord == 2 && yCoord > ySize/2 - 2 && yCoord <= ySize/2 + 2 && zCoord == zSize/2)
    // {
    //   return true;
    // }
    // if (dir == 1 && xCoord == 6 && yCoord > ySize/2 - 2 && yCoord <= ySize/2 + 2 && zCoord == zSize/2)
    // {
    //   return true;
    // }

    // if (dir == 0 && xCoord > 2 && xCoord <= 6 && zCoord == zSize/2 - 2 && yCoord == ySize/2)
    // {
    //   return true;
    // }
    // if (dir == 0 && xCoord > 2 && xCoord <= 6 && zCoord == zSize/2 + 2 && yCoord == ySize/2)
    // {
    //   return true;
    // }
    // if (dir == 2 && xCoord == 2 && zCoord > zSize/2 - 2 && zCoord <= zSize/2 + 2 && yCoord == ySize/2)
    // {
    //   return true;
    // }
    // if (dir == 2 && xCoord == 6 && zCoord > zSize/2 - 2 && zCoord <= zSize/2 + 2 && yCoord == ySize/2)
    // {
    //   return true;
    // }

    // if (dir == 0 && xCoord > 2 && xCoord <= 6 && zCoord == zSize/2 && yCoord == ySize/2)
    // {
    //   return true;
    // }

    // if (dir == 0 && yCoord == ySize/2 && zCoord == zSize/2 && xCoord == xSize/2)
    // {
    //   return false;
    // }
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
    // if (dir == 2 && yCoord == ySize/2 && xCoord == xSize/2+2)
    // {
    //   return true;
    // }

    // if (dir == 1 && xCoord == xSize/2 + 2 && zCoord == zSize/2)
    // {
    //   return true;
    // }

    return false;
  }

}

#endif