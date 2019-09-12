#ifndef ELECTROWEAKTHEORYNONCOMPACT_HPP
#define ELECTROWEAKTHEORYNONCOMPACT_HPP

#include "Theory.hpp"
#include "Su2Tools.hpp"

namespace monsta {
  class ElectroweakTheory: public Theory {
  public:
    ElectroweakTheory(double gaugeCoupling, double mixingAngle, double vev, double selfCoupling);
    ElectroweakTheory(double gaugeCoupling, double mixingAngle, double vev, double selfCoupling, std::vector<int> boundaryConditions);
    ElectroweakTheory(double gaugeCoupling, double mixingAngle, double vev, double selfCoupling, std::vector<int> boundaryConditions, bool tHooftLine);
    double getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    std::complex<double> getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const;
    monsta::Matrix getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;
    double getVev() const { return vev_; }
    double getGaugeCoupling() const { return gaugeCoupling_; }
    double getSelfCoupling() const { return selfCoupling_; }
    monsta::Matrix getSu2Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const;
    std::complex<double> getU1Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const;
    double getHiggsField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    double getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    void postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;
    void applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const;

  private:
    void applyCPeriodicBoundaryConditions(LATfield2::Field< std::complex<double> > &field, int dir) const;
    monsta::Matrix getSu2Plaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const;
    std::complex<double> getU1Plaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const;
    monsta::Matrix getSu2Staple(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, bool isUp) const;
    std::complex<double> getU1Staple(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, bool isUp) const;

    int numFieldMatrices_ = 4;
    int numRows_ = 2;
    int numCols_ = 2;
    int matSize_ = 2;
    double gaugeCoupling_;
    double vev_;
    double selfCoupling_;
    double mixingAngle_;
    std::vector<int> boundaryConditions_ = {0, 0, 0};
    bool tHooftLine_ = false;

    bool monopolesPresent_ = false;
    std::vector<int> monopolePos1_;
    std::vector<int> monopolePos2_;
  };

  ElectroweakTheory::ElectroweakTheory(double gaugeCoupling, double mixingAngle, double vev, double selfCoupling)
  : Theory(numFieldMatrices_, matSize_, matSize_), gaugeCoupling_(gaugeCoupling), mixingAngle_(mixingAngle), vev_(vev), selfCoupling_(selfCoupling)
  {}

  ElectroweakTheory::ElectroweakTheory(double gaugeCoupling, double mixingAngle, double vev, double selfCoupling, std::vector<int> boundaryConditions)
  : Theory(numFieldMatrices_, matSize_, matSize_), gaugeCoupling_(gaugeCoupling), mixingAngle_(mixingAngle), vev_(vev), selfCoupling_(selfCoupling),
    boundaryConditions_(boundaryConditions)
  {}

  ElectroweakTheory::ElectroweakTheory(double gaugeCoupling, double mixingAngle, double vev, double selfCoupling, std::vector<int> boundaryConditions, bool tHooftLine)
  : Theory(numFieldMatrices_, matSize_, matSize_), gaugeCoupling_(gaugeCoupling), mixingAngle_(mixingAngle), vev_(vev), selfCoupling_(selfCoupling),
    boundaryConditions_(boundaryConditions), tHooftLine_(tHooftLine)
  {}

  double ElectroweakTheory::getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    double E = 0;

    // Wilson action
    // E += 2.0/pow(gaugeCoupling_,2)*(2 - real(trace(getSu2Plaquette(field, site, 1, 2))));
    // E += 2.0/pow(gaugeCoupling_,2)*(2 - real(trace(getSu2Plaquette(field, site, 2, 0))));
    // E += 2.0/pow(gaugeCoupling_,2)*(2 - real(trace(getSu2Plaquette(field, site, 0, 1))));

    if (mixingAngle_ > 1e-6)
    {
      E += 1.0/(4*pow(gaugeCoupling_,2)*pow(tan(mixingAngle_),2))*pow(real(getU1Plaquette(field, site, 1, 2)),2);
      E += 1.0/(4*pow(gaugeCoupling_,2)*pow(tan(mixingAngle_),2))*pow(real(getU1Plaquette(field, site, 2, 0)),2);
      E += 1.0/(4*pow(gaugeCoupling_,2)*pow(tan(mixingAngle_),2))*pow(real(getU1Plaquette(field, site, 0, 1)),2);
    }

    // Covariant Derivative
    // for (int ii = 0; ii < 3; ii++)
    // {
    //   LATfield2::Site shiftedSite = site + ii;
    //   E += 2*pow(getHiggsField(field, site),2);
    //   E -= getHiggsField(field, site)*getHiggsField(field, shiftedSite)*real(trace(getSu2Link(field, site, ii)))*cos(real(getU1Link(field, site, ii)));
    //   E += getHiggsField(field, site)*getHiggsField(field, shiftedSite)*imag(trace(getSu2Link(field, site, ii)*pauli3))*sin(real(getU1Link(field, site, ii)));
    // }

    // Higgs Potential
    // E += real(selfCoupling_*pow(2.0*pow(getHiggsField(field, site),2) - pow(vev_, 2),2));

    return E;
  }

  monsta::Matrix ElectroweakTheory::getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
  {
    monsta::Matrix grad(2);
    if (matIdx < 3)
    {
      // int dir1 = (matIdx + 1) % 3;
      // int dir2 = (matIdx + 2) % 3;

      // // Derivative of SU(2) Wilson action
      // Matrix su2PlaquetteDerivMat = getSu2Staple(field, site, matIdx, dir1, true);
      // su2PlaquetteDerivMat = su2PlaquetteDerivMat + getSu2Staple(field, site, matIdx, dir1, false);
      // su2PlaquetteDerivMat = su2PlaquetteDerivMat + getSu2Staple(field, site, matIdx, dir2, true);
      // su2PlaquetteDerivMat = su2PlaquetteDerivMat + getSu2Staple(field, site, matIdx, dir2, false);

      // grad = grad - 2.0/pow(gaugeCoupling_,2)*su2PlaquetteDerivMat;

      // Derivative of kinetic term
      LATfield2::Site tempSite(site);
      tempSite = site + matIdx;
      double scalarField = getHiggsField(field, site);
      double scalarFieldShifted = getHiggsField(field, tempSite);
      Matrix su2GaugeMat = getSu2Link(field, site, matIdx);
      std::complex<double> u1GaugeField = getU1Link(field, site, matIdx);

      Matrix kineticDerivMat = -scalarField*scalarFieldShifted*cos(real(u1GaugeField))*identity;
      kineticDerivMat = kineticDerivMat + 1i*scalarField*scalarFieldShifted*sin(real(u1GaugeField))*pauli3;

      // grad = grad + kineticDerivMat;

      // grad = grad - 0.5*trace(grad*conjugateTranspose(getSu2Link(field, site, matIdx)))*getSu2Link(field, site, matIdx);

    } else {
      for (int ii = 0; ii < 4; ii++)
      {
        if (ii < 3)
        {
          std::complex<double> u1Deriv = 0;
          std::complex<double> u1GaugeField = getU1Link(field, site, ii);
          // Derivative of U(1) Wilson action
          LATfield2::Site tempSite(site);
          int dir1 = (ii + 1) % 3;
          int dir2 = (ii + 2) % 3;
          tempSite = site - dir1;
          std::complex<double> u1PlaquetteDeriv = -2.0*getU1Plaquette(field, site, ii, dir1);
          u1PlaquetteDeriv = u1PlaquetteDeriv - 2.0*getU1Plaquette(field, tempSite, ii, dir1);
          tempSite = site - dir2;
          u1PlaquetteDeriv = u1PlaquetteDeriv = +2.0*getU1Plaquette(field, site, ii, dir2);
          u1PlaquetteDeriv = u1PlaquetteDeriv +2.0*getU1Plaquette(field, tempSite, ii, dir2);

          u1Deriv -= 1.0/(4*(pow(gaugeCoupling_,2)*pow(tan(mixingAngle_),2)))*u1PlaquetteDeriv;

          // Derivative of kinetic term
          tempSite = site + ii;
          double scalarField = getHiggsField(field, site);
          double scalarFieldShifted = getHiggsField(field, tempSite);
          Matrix su2GaugeMat = getSu2Link(field, site, ii);

          std::complex<double> u1KineticDeriv = scalarField*scalarFieldShifted*imag(trace(su2GaugeMat*pauli3))*cos(real(u1GaugeField));
          u1KineticDeriv = u1KineticDeriv + scalarField*scalarFieldShifted*real(trace(su2GaugeMat))*sin(real(u1GaugeField));

          // u1Deriv += u1KineticDeriv;

          // Project perpendicular component:
          // u1Deriv = u1Deriv - real(u1Deriv*conj(u1GaugeField))*u1GaugeField;

          grad((ii + 1) % 2, (ii + 1) / 2) = grad((ii + 1) % 2, (ii + 1) / 2) + u1Deriv;

        } else {
          for (int dir = 0; dir < 3; dir++)
          {
            // Deriviative of kinetic term
            LATfield2::Site tempSite = site + dir;
            double scalarField = getHiggsField(field, site);
            double scalarFieldShiftedFwd = getHiggsField(field, tempSite);
            Matrix su2GaugeMat = getSu2Link(field, site, dir);
            std::complex<double> u1GaugeField = getU1Link(field, site, dir);

            tempSite = site - dir;
            Matrix su2GaugeMatShiftedBwd = getSu2Link(field, tempSite, dir);
            double scalarFieldShiftedBwd = getHiggsField(field, tempSite);
            std::complex<double> u1GaugeFieldShiftedBwd = getU1Link(field, tempSite, dir);

            double higgsKineticDeriv = 0;
            higgsKineticDeriv += 4*scalarField;
            higgsKineticDeriv -= scalarFieldShiftedFwd*real(trace(su2GaugeMat))*cos(real(u1GaugeField));
            higgsKineticDeriv += scalarFieldShiftedFwd*imag(trace(su2GaugeMat*pauli3))*sin(real(u1GaugeField));
            higgsKineticDeriv -= scalarFieldShiftedBwd*real(trace(su2GaugeMatShiftedBwd))*cos(real(u1GaugeFieldShiftedBwd));
            higgsKineticDeriv += scalarFieldShiftedBwd*imag(trace(su2GaugeMatShiftedBwd*pauli3))*sin(real(u1GaugeFieldShiftedBwd));

            // grad(0,0) = grad(0,0) + higgsKineticDeriv;
          }

          // Derivative of Higgs Potential
          // grad(0,0) = grad(0,0) + 8.0*selfCoupling_*getHiggsField(field, site)*(2.0*pow(getHiggsField(field, site),2) - pow(vev_, 2));
        }
      }
    }

    return grad;
  }

  std::complex<double> ElectroweakTheory::getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const
  {
    monsta::Matrix gradMat = getLocalGradient(field, site, matIdx);
    return gradMat(rowIdx, colIdx);
    return 0;
  }

  void ElectroweakTheory::postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
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
      for (int ii = 0; ii < 4; ii++)
      {
        if (ii < 3) // Project to U(1)
        {
          double norm = abs(getU1Link(field, site, ii));
          field(site, matIdx, (ii + 1) % 2, (ii + 1) / 2) = field(site, matIdx, (ii + 1) % 2, (ii + 1) / 2) / norm;
          if (norm == 0)
          {
            field(site, matIdx, (ii + 1) % 2, (ii + 1) / 2) = 1;
          }
        }
        else
        {
          field(site, matIdx, 0, 0) = real(field(site, matIdx, 0, 0));
        }
      }
    }
  }

  void ElectroweakTheory::applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const
  {
    field.updateHalo();
    for (int dir = 0; dir < 3; dir++)
    {
      applyCPeriodicBoundaryConditions(field, dir);
    }
  }

  void ElectroweakTheory::applyCPeriodicBoundaryConditions(LATfield2::Field< std::complex<double> > &field, int dir) const
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
            // if (tHooftLine_ && dir == 2 && site.coord(1) == 0)
            // {
            //   if (matIdx == 1)
            //   {
            //     boundaryMat = -1.0*boundaryMat;
            //   }
            // }
          }
          else
          {
            for (int ii = 0; ii < 4; ii ++)
            {
              if (ii < 3)
              {
                if (pauliMatNum != 0)
                {
                  boundaryMat((ii + 1) % 2, (ii + 1) / 2) = conj(boundaryMat((ii + 1) % 2, (ii + 1) / 2));
                }
              } else {
                continue;
              }
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

  monsta::Matrix ElectroweakTheory::getSu2Plaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix plaquette = getSu2Link(field, site, dir1);
    tempSite = tempSite+dir1;
    plaquette = plaquette*getSu2Link(field, tempSite, dir2);
    tempSite = tempSite-dir1+dir2;
    plaquette = plaquette*conjugateTranspose(getSu2Link(field, tempSite, dir1));
    tempSite = tempSite-dir2;
    plaquette = plaquette*conjugateTranspose(getSu2Link(field, tempSite, dir2));

    return plaquette;
  }

  std::complex<double> ElectroweakTheory::getU1Plaquette(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2) const
  {
    LATfield2::Site tempSite(site);
    std::complex<double> plaquette = getU1Link(field, site, dir1);
    tempSite = tempSite+dir1;
    plaquette = plaquette + getU1Link(field, tempSite, dir2);
    tempSite = tempSite-dir1+dir2;
    plaquette = plaquette - getU1Link(field, tempSite, dir1);
    tempSite = tempSite-dir2;
    plaquette = plaquette - getU1Link(field, tempSite, dir2);

    return plaquette;
  }

  monsta::Matrix ElectroweakTheory::getSu2Staple(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, bool isUp) const
  {
    LATfield2::Site tempSite(site);
    monsta::Matrix staple = getSu2Link(field, site, dir2);

    if (isUp)
    {
      tempSite = tempSite + dir2;
      staple = staple*getSu2Link(field, tempSite, dir1);
      tempSite = tempSite-dir2+dir1;
      staple = staple*conjugateTranspose(getSu2Link(field, tempSite, dir2));
    }
    else
    {
      tempSite = tempSite - dir2;
      staple = conjugateTranspose(getSu2Link(field, tempSite, dir2));
      staple = staple*getSu2Link(field, tempSite, dir1);
      tempSite = tempSite + dir1;
      staple = staple*getSu2Link(field, tempSite, dir2);
    }


    return staple;
  }

  std::complex<double> ElectroweakTheory::getU1Staple(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir1, int dir2, bool isUp) const
  {
    LATfield2::Site tempSite(site);
    std::complex<double> staple = getU1Link(field, site, dir2);

    if (isUp)
    {
      tempSite = tempSite + dir2;
      staple = staple + getU1Link(field, tempSite, dir1);
      tempSite = tempSite-dir2+dir1;
      staple = staple - conj(getU1Link(field, tempSite, dir2));
    }
    else
    {
      tempSite = tempSite - dir2;
      staple = -conj(getU1Link(field, tempSite, dir2));
      staple = staple + getU1Link(field, tempSite, dir1);
      tempSite = tempSite + dir1;
      staple = staple - getU1Link(field, tempSite, dir2);
    }


    return staple;
  }

  double ElectroweakTheory::getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
  {
    int dir1 = (cpt + 1) % 3;
    int dir2 = (cpt + 2) % 3;


    // LATfield2::Site tempSite(site);
    // std::complex<double> zPlaquette = field(site, dir1, 0, 0);
    // tempSite = tempSite+dir1;
    // zPlaquette = zPlaquette*field(tempSite, dir2, 0, 0);
    // tempSite = tempSite-dir1+dir2;
    // zPlaquette = zPlaquette*conj(field(tempSite, dir1, 0, 0));
    // tempSite = tempSite-dir2;
    // zPlaquette = zPlaquette*conj(field(tempSite, dir2, 0, 0));

    std::complex<double> zPlaquette = getU1Plaquette(field, site, dir1, dir2);

    double magneticField = 2./gaugeCoupling_*arg(zPlaquette);

    return (magneticField);
  }

  monsta::Matrix ElectroweakTheory::getSu2Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
  {
    return Matrix(field, site, cpt);
  }

  std::complex<double> ElectroweakTheory::getU1Link(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
  {
    return real(field(site, 3, (cpt+ 1) % 2, (cpt + 1) / 2));
  }

  double ElectroweakTheory::getHiggsField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    return real(field(site, 3, 0, 0));
  }


}

#endif