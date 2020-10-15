#ifndef GRADDESCENTSOLVERBBSTEPNOCHECKERBOARD_HPP
#define GRADDESCENTSOLVERBBSTEPNOCHECKERBOARD_HPP

#include <ctime>
#include <iomanip>


namespace monsta
{
  class GradDescentSolverNoCheckerboard
  {
  public:
   GradDescentSolverNoCheckerboard(double tol, int maxIterations, double initialStepSize, double maxStepSize);
    GradDescentSolverNoCheckerboard(double tol, int maxIterations, double initialStepSize, double maxStepSize, int minSteps);
    GradDescentSolverNoCheckerboard(double tol, int maxIterations, double initialStepSize, double maxStepSize, int minSteps, bool minimiseGrad);

    void setVerbosity(bool isVerbose) { isVerbose_ = isVerbose; };

    void solve(monsta::Theory &theory, LATfield2::Field< std::complex<double> > &field);

    void setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize);
    void setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize, std::vector<int> skipCpts);
    void setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize, int minSteps);

  private:
    double tol_;
    int maxIterations_;
    double stepSize_;
    double maxStepSize_;
    double maxGrad_ = 1e6;
    double relEnergyChange = 1e6;
    double energy;
    double energyOld;
    bool isVerbose_ = true;
    int minSteps_ = 0;
    bool minimiseGrad_ = false;
    std::vector<int> skipCpts_;
    LATfield2::Field< std::complex<double> > oldGrads_;
    LATfield2::Field< std::complex<double> > gradField_;
    int numRows_;
    int numCols_;
    int numMatrices_;


    void iterate(LATfield2::Field< std::complex<double> > &field, monsta::Theory &theory);
  };

  GradDescentSolverNoCheckerboard::GradDescentSolverNoCheckerboard(double tol, int maxIterations, double initialStepSize, double maxStepSize)
  : tol_(tol), maxIterations_(maxIterations), stepSize_(initialStepSize), maxStepSize_(maxStepSize) {}
  GradDescentSolverNoCheckerboard::GradDescentSolverNoCheckerboard(double tol, int maxIterations, double initialStepSize, double maxStepSize, int minSteps)
  : tol_(tol), maxIterations_(maxIterations), stepSize_(initialStepSize), maxStepSize_(maxStepSize), minSteps_(minSteps) {}
  GradDescentSolverNoCheckerboard::GradDescentSolverNoCheckerboard(double tol, int maxIterations, double initialStepSize, double maxStepSize, int minSteps, bool minimiseGrad)
  : tol_(tol), maxIterations_(maxIterations), stepSize_(initialStepSize), maxStepSize_(maxStepSize), minSteps_(minSteps), minimiseGrad_(minimiseGrad) {}

  void GradDescentSolverNoCheckerboard::setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize)
  {
    tol_ = tol;
    maxIterations_ = maxIterations;
    stepSize_ = initialStepSize;
    maxStepSize_ = maxStepSize;
  }

  void GradDescentSolverNoCheckerboard::setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize, std::vector<int> skipCpts)
  {
    tol_ = tol;
    maxIterations_ = maxIterations;
    stepSize_ = initialStepSize;
    maxStepSize_ = maxStepSize;
    skipCpts_ = skipCpts;
  }

  void GradDescentSolverNoCheckerboard::setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize, int minSteps)
  {
    tol_ = tol;
    maxIterations_ = maxIterations;
    stepSize_ = initialStepSize;
    maxStepSize_ = maxStepSize;
    minSteps_ = minSteps;
  }

  void GradDescentSolverNoCheckerboard::solve(monsta::Theory &theory, LATfield2::Field< std::complex<double> > &field)
  {
    numRows_ = field.rows();
    numCols_ = field.cols();
    numMatrices_ = field.components() / (numRows_*numCols_);

    energy = theory.computeEnergy(field);
    int numMatrices_ = field.components() / (field.rows()*field.cols());
    oldGrads_.initialize(field.lattice(), numMatrices_, field.rows(), field.cols(), field.symmetry());
    oldGrads_.alloc();

    gradField_.initialize(field.lattice(), numMatrices_, field.rows(), field.cols(), field.symmetry());
    gradField_.alloc();

    iterate(field, theory);
    energyOld = energy;
    energy = theory.computeEnergy(field);
    relEnergyChange = (energy - energyOld);
    int numIters = 1;

    double maxGradOld = 0;
    if (minimiseGrad_)
    {
      while (maxGradOld < maxGrad_)
      {
        numIters++;
        maxGradOld = maxGrad_;
        iterate(field, theory);
        energyOld = energy;
        energy = theory.computeEnergy(field);
        relEnergyChange = abs(energy - energyOld) / energy;
        if (isVerbose_)
        {
          COUT << std::fixed;
          COUT << std::setprecision(10);
          COUT << energy << std::endl;
          COUT << maxGrad_ << std::endl;
        }
      }
    }
    while (numIters < minSteps_)
    {
      numIters++;
      maxGradOld = maxGrad_;
      iterate(field, theory);
      energyOld = energy;
      energy = theory.computeEnergy(field);
      relEnergyChange = (energy - energyOld);
      if (isVerbose_)
      {
        COUT << std::fixed;
        COUT << std::setprecision(10);
        COUT << energy << std::endl;
        COUT << maxGrad_ << std::endl;
      }
    }

    while (abs(maxGrad_) > tol_ && numIters < maxIterations_)
    {
      if (minimiseGrad_ && maxGradOld < maxGrad_) { break; }
      numIters++;
      maxGradOld = maxGrad_;
      iterate(field, theory);
      energyOld = energy;
      energy = theory.computeEnergy(field);
      relEnergyChange = (energy - energyOld);
      if (isVerbose_)
      {
        COUT << std::fixed;
        COUT << std::setprecision(10);
        COUT << energy << std::endl;
        COUT << maxGrad_ << std::endl;
      }
    }

    double finalEnergy = theory.computeEnergy(field);

    if (numIters < maxIterations_) {
      COUT << "Gradient descent finished in " << numIters << " iterations." << std::endl;
      COUT << "Minimum energy: " << finalEnergy << std::endl;
    } else {
      COUT << "Gradient descent aborted after " << maxIterations_ << " iterations." << std::endl;
      // COUT << "Maximum gradient: " << maxGrad_ << std::endl;
      COUT << "Energy reached: " << finalEnergy << std::endl;
    }
    parallel.barrier();
  }

  void GradDescentSolverNoCheckerboard::iterate(LATfield2::Field< std::complex<double> > &field, monsta::Theory &theory)
  {
    clock_t begin = clock();
    LATfield2::Site site(field.lattice());

    std::vector< std::complex<double> > gradVals(numRows_*numCols_*numMatrices_);
    int vecIdx;
    monsta::Matrix gradMat(numRows_);
    std::complex<double> gradVal;
    std::complex<double> oldGradVal;
    std::complex<double> fieldVal;

    double stepChangeNumerator = 0;
    double stepChangeDenominator = 0;

    double maxGrad = 0;

    for (site.first(); site.test(); site.next())
    {
      for (int matIdx = 0; matIdx < numMatrices_; matIdx++)
      {
        for (int rowIdx = 0; rowIdx < numRows_; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols_; colIdx++)
          {
            gradField_(site, matIdx, rowIdx, colIdx) = theory.getLocalGradient(field, site, matIdx, rowIdx, colIdx);
          }
        }
      }
    }

    for (site.first(); site.test(); site.next())
    {
      for (int matIdx = 0; matIdx < numMatrices_; matIdx++)
      {
        bool skip = false;
        for (int ii = 0; ii < skipCpts_.size(); ii++)
        {
          if (matIdx == skipCpts_[ii])
          { 
            skip = true;
            continue;
          }
        }
        if (skip) { continue; }
        Matrix gradMat(gradField_, site, matIdx);
        Matrix gradProjOld(oldGrads_, site, matIdx);

        Matrix gradProj = gradMat;
        if (matIdx < 3)
        {
          gradProj = gradMat - 0.5*trace(gradMat*conjugateTranspose(Matrix(field, site, matIdx)))*Matrix(field, site, matIdx);
          stepChangeNumerator += 0.5*real(trace(gradProjOld*conjugateTranspose(gradProjOld - gradProj)));
          stepChangeDenominator += 0.5*real(trace((gradProjOld - gradProj)*(gradProjOld - gradProj)));
          for (int rowIdx = 0; rowIdx < numRows_; rowIdx++)
          {
            for (int colIdx = 0; colIdx < numCols_; colIdx++)
            {
              oldGrads_(site, matIdx, rowIdx, colIdx) = gradProj(rowIdx, colIdx);
            }
          }
        }
        else
        {
          for (int rowIdx = 0; rowIdx < numRows_; rowIdx++)
          {
            for (int colIdx = 0; colIdx < numCols_; colIdx++)
            {
              oldGradVal = oldGrads_(site, matIdx, rowIdx, colIdx);
              gradVal = gradProj(rowIdx, colIdx);

              stepChangeNumerator += abs(real(oldGradVal) * (real(oldGradVal - gradVal)));
              stepChangeNumerator += abs(imag(oldGradVal) * (imag(oldGradVal - gradVal)));

              stepChangeDenominator += pow(abs(gradVal - oldGradVal),2);
              oldGrads_(site, matIdx, rowIdx, colIdx) = gradProj(rowIdx, colIdx);
            }
          }
        }
        for (int rowIdx = 0; rowIdx < numRows_; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols_; colIdx++)
          {
            vecIdx = colIdx + numCols_ * (rowIdx + numRows_ * matIdx);
            gradVal = gradProj(rowIdx, colIdx);
            if (abs(gradVal) > abs(maxGrad)) { maxGrad = abs(gradVal); }
            gradVals[vecIdx] = gradVal;
          }
        }
      }
      for (int matIdx = 0; matIdx < numMatrices_; matIdx++)
      {
        for (int rowIdx = 0; rowIdx < numRows_; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols_; colIdx++)
          {
            vecIdx = colIdx + numCols_ * (rowIdx + numRows_ * matIdx);
            fieldVal = field(site, matIdx, rowIdx, colIdx);
            field(site, matIdx, rowIdx, colIdx) = fieldVal - stepSize_*gradVals[vecIdx];
          }
        }
        theory.postProcess(field, site, matIdx);
      }
    }

    parallel.sum(stepChangeNumerator);
    parallel.sum(stepChangeDenominator);

    double stepChange = stepChangeNumerator / stepChangeDenominator;
    if (stepChange > 1e-8) { stepSize_ *= stepChange; }
    if (stepSize_ > maxStepSize_) { stepSize_ = maxStepSize_; }

    theory.applyBoundaryConditions(field);
    // double stepChange = system.updateGradientReturnStepChange();
    // stepSize_ *= stepChange;
    maxGrad_ = maxGrad;
    parallel.max(maxGrad_);

    clock_t end = clock();
    // COUT << double(end - begin) / CLOCKS_PER_SEC << endl;
  }
}

#endif