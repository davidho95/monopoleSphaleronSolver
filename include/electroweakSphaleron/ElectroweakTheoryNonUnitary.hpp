#ifndef ELECTROWEAKTHEORYNONUNITARY_HPP
#define ELECTROWEAKTHEORYNONUNITARY_HPP

#include "../Theory.hpp"
#include "../Su2Tools.hpp"

namespace monsta {
  class ElectroweakTheoryNonUnitary: public ElectroweakTheory {
  public:
    ElectroweakTheoryNonUnitary(double gaugeCoupling, double tanSqMixingAngle, double vev, double selfCoupling);
    ElectroweakTheoryNonUnitary(double gaugeCoupling, double tanSqMixingAngle, double vev, double selfCoupling, std::vector<int> boundaryConditions);
    ElectroweakTheoryNonUnitary(double gaugeCoupling, double tanSqMixingAngle, double vev, double selfCoupling, std::vector<int> boundaryConditions, bool tHooftLine);
    double getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    std::complex<double> getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const;
    monsta::Matrix getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;
    double getVev() const { return vev_; }
    double getGaugeCoupling() const { return gaugeCoupling_; }
    double getSelfCoupling() const { return selfCoupling_; }
    monsta::Matrix getHiggsField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    double getHiggsMagnitude(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const;
    double getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int dir) const;
    void postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const;
    void applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const;

  private:
    void applyCPeriodicBoundaryConditions(LATfield2::Field< std::complex<double> > &field, int dir) const;

  };

  ElectroweakTheoryNonUnitary::ElectroweakTheoryNonUnitary(double gaugeCoupling, double tanSqMixingAngle, double vev, double selfCoupling)
  : ElectroweakTheory(gaugeCoupling, tanSqMixingAngle, vev, selfCoupling)
  {
    numFieldMatrices_ = 5;
  }

  ElectroweakTheoryNonUnitary::ElectroweakTheoryNonUnitary(double gaugeCoupling, double tanSqMixingAngle, double vev, double selfCoupling, std::vector<int> boundaryConditions)
  : ElectroweakTheory(gaugeCoupling, tanSqMixingAngle, vev, selfCoupling, boundaryConditions)
  {
    numFieldMatrices_ = 5;
  }

  ElectroweakTheoryNonUnitary::ElectroweakTheoryNonUnitary(double gaugeCoupling, double tanSqMixingAngle, double vev, double selfCoupling, std::vector<int> boundaryConditions, bool tHooftLine)
  : ElectroweakTheory(gaugeCoupling, tanSqMixingAngle, vev, selfCoupling, boundaryConditions, tHooftLine)
  {
    numFieldMatrices_ = 5;
  }

  double ElectroweakTheoryNonUnitary::getLocalEnergyDensity(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    double E = 0;

    // Wilson action
    E += 2.0/pow(gaugeCoupling_,2)*(2 - real(trace(getSu2Plaquette(field, site, 1, 2))));
    E += 2.0/pow(gaugeCoupling_,2)*(2 - real(trace(getSu2Plaquette(field, site, 2, 0))));
    E += 2.0/pow(gaugeCoupling_,2)*(2 - real(trace(getSu2Plaquette(field, site, 0, 1))));

    if (tanSqMixingAngle_ > 1e-6)
    {
      E += 2.0/(pow(gaugeCoupling_,2)*tanSqMixingAngle_)*(1 - real(getU1Plaquette(field, site, 1, 2)));
      E += 2.0/(pow(gaugeCoupling_,2)*tanSqMixingAngle_)*(1 - real(getU1Plaquette(field, site, 2, 0)));
      E += 2.0/(pow(gaugeCoupling_,2)*tanSqMixingAngle_)*(1 - real(getU1Plaquette(field, site, 0, 1)));
    }

    // Covariant Derivative
    for (int ii = 0; ii < 3; ii++)
    {
      LATfield2::Site shiftedSite = site + ii;

      monsta::Matrix higgsVec = getHiggsField(field, site);
      monsta::Matrix higgsVecShifted = getHiggsField(field, shiftedSite);

      E += pow(getHiggsMagnitude(field, site),2);
      E += pow(getHiggsMagnitude(field, shiftedSite),2);

      monsta::Matrix crossTermMat = conj(getU1Link(field, site, ii))*conjugateTranspose(higgsVecShifted)*conjugateTranspose(getSu2Link(field, site, ii))*higgsVec;
      crossTermMat = crossTermMat + getU1Link(field, site, ii)*conjugateTranspose(higgsVec)*getSu2Link(field, site, ii)*higgsVecShifted;

      E -= real(crossTermMat(0,0));
    }

    // Higgs Potential
    E += real(selfCoupling_*pow(2.0*pow(getHiggsMagnitude(field, site),2) - pow(vev_, 2),2));

    return E;
  }

  monsta::Matrix ElectroweakTheoryNonUnitary::getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
  {
    monsta::Matrix grad(2);
    if (matIdx < 3)
    {
      int dir1 = (matIdx + 1) % 3;
      int dir2 = (matIdx + 2) % 3;

      // Derivative of SU(2) Wilson action
      Matrix su2PlaquetteDerivMat = getSu2Staple(field, site, matIdx, dir1, true);
      su2PlaquetteDerivMat = su2PlaquetteDerivMat + getSu2Staple(field, site, matIdx, dir1, false);
      su2PlaquetteDerivMat = su2PlaquetteDerivMat + getSu2Staple(field, site, matIdx, dir2, true);
      su2PlaquetteDerivMat = su2PlaquetteDerivMat + getSu2Staple(field, site, matIdx, dir2, false);

      grad = grad - 2.0/pow(gaugeCoupling_,2)*su2PlaquetteDerivMat;

      // Derivative of kinetic term
      LATfield2::Site tempSite(site);
      tempSite = site + matIdx;
      Matrix higgsVec = getHiggsField(field, site);
      Matrix higgsVecShifted = getHiggsField(field, tempSite);
      Matrix su2GaugeMat = getSu2Link(field, site, matIdx);
      std::complex<double> u1GaugeField = getU1Link(field, site, matIdx);


      Matrix kineticDerivMat = -conj(u1GaugeField)*higgsVec*conjugateTranspose(higgsVecShifted);
      kineticDerivMat = kineticDerivMat - conjugateTranspose(u1GaugeField*higgsVecShifted*conjugateTranspose(higgsVec));

      // kineticDerivMat = 0.5*(kineticDerivMat + conjugateTranspose(kineticDerivMat));
      // std::complex<double> tr = trace(kineticDerivMat);
      // kineticDerivMat(0,0) -= 0.5*tr;
      // kineticDerivMat(1,1) -= 0.5*tr;


      // grad = grad + kineticDerivMat;

      // grad = 0.5*(grad + conjugateTranspose(grad));
      // std::complex<double> tr = trace(grad);
      // grad(0,0) -= 0.5*tr;
      // grad(1,1) -= 0.5*tr;

      // grad = grad - 0.5*trace(grad*conjugateTranspose(getSu2Link(field, site, matIdx)))*getSu2Link(field, site, matIdx);

    }
    if (matIdx == 3)
    {
      for (int ii = 0; ii < 3; ii++)
      {
        std::complex<double> u1Deriv = 0;
        // Derivative of U(1) Wilson action
        int dir1 = (ii + 1) % 3;
        int dir2 = (ii + 2) % 3;
        std::complex<double> u1PlaquetteDeriv = getU1Staple(field, site, ii, dir1, true);
        u1PlaquetteDeriv = u1PlaquetteDeriv + getU1Staple(field, site, ii, dir1, false);
        u1PlaquetteDeriv = u1PlaquetteDeriv + getU1Staple(field, site, ii, dir2, true);
        u1PlaquetteDeriv = u1PlaquetteDeriv + getU1Staple(field, site, ii, dir2, false);

        if (tanSqMixingAngle_ > 1e-6)
        {
          u1Deriv -= 2.0/(pow(gaugeCoupling_,2)*tanSqMixingAngle_)*u1PlaquetteDeriv;
        }

        // Derivative of kinetic term
        LATfield2::Site tempSite(site);
        tempSite = site + ii;
        Matrix higgsVec = getHiggsField(field, site);
        Matrix higgsVecShifted = getHiggsField(field, tempSite);
        Matrix su2GaugeMat = getSu2Link(field, site, ii);
        std::complex<double> u1GaugeField = getU1Link(field, site, ii);

        monsta::Matrix crossTermMat = conjugateTranspose(higgsVecShifted)*conjugateTranspose(getSu2Link(field, site, ii))*higgsVec;
        std::complex<double> u1KineticDeriv = -2i*imag(crossTermMat(0,0));
        u1KineticDeriv -= 2*real(crossTermMat(0,0));

        u1Deriv += u1KineticDeriv;

        // Project perpendicular component:
        // u1Deriv = u1Deriv - real(u1Deriv*conj(u1GaugeField))*u1GaugeField;

        grad((ii + 1) % 2, (ii + 1) / 2) = grad((ii + 1) % 2, (ii + 1) / 2) + u1Deriv;
      }
    }
    if (matIdx == 4)
    {
      for (int dir = 0; dir < 3; dir++)
      {
        // Deriviative of kinetic term
        LATfield2::Site tempSite = site + dir;
        Matrix higgsField = getHiggsField(field, site);
        Matrix higgsFieldShiftedFwd = getHiggsField(field, tempSite);
        Matrix su2GaugeMat = getSu2Link(field, site, dir);
        std::complex<double> u1GaugeField = getU1Link(field, site, dir);

        tempSite = site - dir;
        Matrix su2GaugeMatShiftedBwd = getSu2Link(field, tempSite, dir);
        Matrix higgsFieldShiftedBwd = getHiggsField(field, tempSite);
        std::complex<double> u1GaugeFieldShiftedBwd = getU1Link(field, tempSite, dir);

        Matrix higgsDerivMat = 4*higgsField;
        higgsDerivMat = higgsDerivMat - u1GaugeField*su2GaugeMat*higgsFieldShiftedFwd;
        higgsDerivMat = higgsDerivMat - u1GaugeField*su2GaugeMat*higgsFieldShiftedFwd;
        higgsDerivMat = higgsDerivMat - conj(u1GaugeFieldShiftedBwd)*conjugateTranspose(su2GaugeMatShiftedBwd)*higgsFieldShiftedBwd;
        higgsDerivMat = higgsDerivMat - conj(u1GaugeFieldShiftedBwd)*conjugateTranspose(su2GaugeMatShiftedBwd)*higgsFieldShiftedBwd;

        grad(0,0) = grad(0,0) + higgsDerivMat(0,0);
        grad(1,0) = grad(1,0) + higgsDerivMat(1,0);

      }

        // Derivative of Higgs Potential
        grad(0,0) = grad(0,0) + 8.0*selfCoupling_*field(site, 4, 0, 0)*(2.0*pow(getHiggsMagnitude(field, site),2) - pow(vev_, 2));
        grad(1,0) = grad(1,0) + 8.0*selfCoupling_*field(site, 4, 1, 0)*(2.0*pow(getHiggsMagnitude(field, site),2) - pow(vev_, 2));

    }
    return grad;
  }

  std::complex<double> ElectroweakTheoryNonUnitary::getLocalGradient(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const
  {
    monsta::Matrix gradMat = getLocalGradient(field, site, matIdx);
    return gradMat(rowIdx, colIdx);
    return 0;
  }

  void ElectroweakTheoryNonUnitary::postProcess(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int matIdx) const
  {
    if (matIdx < 3) // Project to SU(2)
    {
      field(site, matIdx, 0, 0) = conj(field(site, matIdx, 1, 1));
      field(site, matIdx, 1, 0) = -1.0*conj(field(site, matIdx, 0, 1));
      std::complex<double> determinant = field(site, matIdx, 0, 0)*field(site, matIdx, 1, 1) - field(site, matIdx, 0, 1)*field(site, matIdx, 1, 0);
      for (int rowIdx = 0; rowIdx < numRows_; rowIdx++)
      {
        for (int colIdx = 0; colIdx < numCols_; colIdx++)
        {
          field(site, matIdx, rowIdx, colIdx) = field(site, matIdx, rowIdx, colIdx) / sqrt(determinant);
        }
      }
    }
    if (matIdx == 3)
    {
      field(site, matIdx, 0, 0) = 0;
      for (int ii = 0; ii < 3; ii++)
      {
        double norm = abs(getU1Link(field, site, ii));
        field(site, matIdx, (ii + 1) % 2, (ii + 1) / 2) = field(site, matIdx, (ii + 1) % 2, (ii + 1) / 2) / norm;
        if (norm == 0)
        {
          field(site, matIdx, (ii + 1) % 2, (ii + 1) / 2) = 1;
        }
      }
    }
    if (matIdx == 4)
    {
      field(site, matIdx, 0, 1) = 0;
      field(site, matIdx, 1, 1) = 0;
    }
  }

  void ElectroweakTheoryNonUnitary::applyBoundaryConditions(LATfield2::Field< std::complex<double> > &field) const
  {
    field.updateHalo();
    for (int dir = 0; dir < 3; dir++)
    {
      applyCPeriodicBoundaryConditions(field, dir);
    }
  }

  void ElectroweakTheoryNonUnitary::applyCPeriodicBoundaryConditions(LATfield2::Field< std::complex<double> > &field, int dir) const
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
            for (int ii = 0; ii < 4; ii ++)
            {
              if (ii < 3)
              {
                if (tHooftLine_ && dir == 2 && site.coord(1) == 0)
                {
                  if (ii == 1)
                  {
                    boundaryMat((ii + 1) % 2, (ii + 1) / 2) = -boundaryMat((ii + 1) % 2, (ii + 1) / 2);
                  }
                }
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

  double ElectroweakTheoryNonUnitary::getMagneticField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, int cpt) const
  {
    int dir1 = (cpt + 1) % 3;
    int dir2 = (cpt + 2) % 3;


    LATfield2::Site tempSite(site);
    std::complex<double> zPlaquette = field(site, dir1, 0, 0);
    tempSite = tempSite+dir1;
    zPlaquette = zPlaquette*field(tempSite, dir2, 0, 0);
    tempSite = tempSite-dir1+dir2;
    zPlaquette = zPlaquette*conj(field(tempSite, dir1, 0, 0));
    tempSite = tempSite-dir2;
    zPlaquette = zPlaquette*conj(field(tempSite, dir2, 0, 0));

    std::complex<double> u1Plaquette = getU1Plaquette(field, site, dir1, dir2);

    double magneticField =  2./gaugeCoupling_*arg(zPlaquette);
    if (tanSqMixingAngle_ > 1e-6)
    {
      magneticField += 2./(gaugeCoupling_*sqrt(tanSqMixingAngle_)) * arg(u1Plaquette);
    }

    return (magneticField);
  }

  monsta::Matrix ElectroweakTheoryNonUnitary::getHiggsField(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    return Matrix(field, site, 4);
  }

  double ElectroweakTheoryNonUnitary::getHiggsMagnitude(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site) const
  {
    // return real(field(site, 3, 0, 0));
    return sqrt(real(field(site, 4, 0, 0)*conj(field(site, 4, 0, 0)) + field(site, 4, 1, 0)*conj(field(site, 4, 1, 0))));
  }

}

#endif