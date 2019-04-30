#ifndef THEORYCHECKER_HPP
#define THEORYCHECKER_HPP

namespace monsta {
  class TheoryChecker {
  public:
    TheoryChecker(double tol, bool verbose = true)
    : tol_(tol), verbose_(verbose) {
      if (parallel.size() != 1)
      {
        throw std::runtime_error("Gradient checker will not perform correctly with multiple processors.");
      }
    }

    bool checkGradients(monsta::Theory &theory, LATfield2::Field< std::complex<double> > &field) const
    {
      LATfield2::Site site(field.lattice());
      bool checkPassed;
      int numRows = field.rows();
      int numCols = field.cols();
      int numMatrices = field.components() / (numRows*numCols);

      double maxErr = 0;
      int count = 0;
      for (site.first(); site.test(); site.next())
      {
        for (int matIdx = 0; matIdx < numMatrices; matIdx++)
        {
          for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
          {
            for (int colIdx = 0; colIdx < numCols; colIdx++)
            {
              std::complex<double> suppliedGradVal = theory.getLocalGradient(field, site, matIdx, rowIdx, colIdx);
              std::complex<double> fdGradVal = getFdGradient(theory, field, site, matIdx, rowIdx, colIdx);
              double relError;
              if (abs(suppliedGradVal) > tol_)
              {
                relError = abs(suppliedGradVal - fdGradVal) / abs(fdGradVal);
              }
              else
              {
                relError = abs(suppliedGradVal - fdGradVal);
              }
              if (abs(relError) > abs(maxErr))
              {
                maxErr = relError;
              }

              if (verbose_ && abs(relError) > tol_)
              {
                std::cout << "(" << matIdx << " " << rowIdx << " " << colIdx << ";"
                                 << site.coord(0) << ","
                                 << site.coord(1) << ","
                                 << site.coord(2) << ") : ";
                std::cout << suppliedGradVal << " " << fdGradVal << endl;
              }
            }
          }
        }
      }
      checkPassed = (abs(maxErr) < tol_);
      if (checkPassed)
      {
        std::cout << "Gradient check passed: maximum deviation from supplied gradient: " << maxErr << std::endl;
      }
      else
      {
        std::cout << "Gradient check failed: maximum deviation from supplied gradient: " << maxErr << std::endl;
      }
      return checkPassed;
    }

  private:
    double tol_;
    bool verbose_;
    double fdStep_ = 1e-5;

    std::complex<double> getFdGradient(monsta::Theory &theory, LATfield2::Field< std::complex<double> > &field,
      LATfield2::Site &site, int matIdx, int rowIdx, int colIdx) const
    {
      std::complex<double> originalFieldVal = field(site, matIdx, rowIdx, colIdx);

      //Compute gradient numerically
      field(site,  matIdx, rowIdx, colIdx) = originalFieldVal + fdStep_;
      theory.applyBoundaryConditions(field);
      theory.applyBoundaryConditions(field);
      double EPlus = theory.computeEnergy(field);
      field(site,  matIdx, rowIdx, colIdx) = originalFieldVal - fdStep_;
      theory.applyBoundaryConditions(field);
      theory.applyBoundaryConditions(field);
      double EMinus = theory.computeEnergy(field);

      field(site,  matIdx, rowIdx, colIdx) = originalFieldVal;
      theory.applyBoundaryConditions(field);
      theory.applyBoundaryConditions(field);

      double gradReal = (EPlus - EMinus) / (2*fdStep_);

      std::complex<double> imagStep(0, fdStep_);
      field(site,  matIdx, rowIdx, colIdx) = originalFieldVal + imagStep;
      theory.applyBoundaryConditions(field);
      theory.applyBoundaryConditions(field);
      EPlus = theory.computeEnergy(field);
      field(site,  matIdx, rowIdx, colIdx) = originalFieldVal - imagStep;
      theory.applyBoundaryConditions(field);
      theory.applyBoundaryConditions(field);
      EMinus = theory.computeEnergy(field);

      field(site,  matIdx, rowIdx, colIdx) = originalFieldVal;
      theory.applyBoundaryConditions(field);
      theory.applyBoundaryConditions(field);

      double gradImag = (EPlus - EMinus) / (2*fdStep_);

      std::complex<double> gradCplx(gradReal, gradImag);

      return gradCplx;
    }
  };
}

#endif