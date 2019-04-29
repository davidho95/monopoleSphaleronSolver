#ifndef GRADDESCENTSOLVER_HPP
#define GRADDESCENTSOLVER_HPP

#include <ctime>

namespace monsta
{
  class GradDescentSolver
  {
  public:
    GradDescentSolver(double tol, int maxIterations, double initialStepSize, double maxStepSize);
    GradDescentSolver(double tol, int maxIterations, double initialStepSize, double maxStepSize, std::vector<int> skipCpts);

    void setVerbosity(bool isVerbose) { isVerbose_ = isVerbose; };

    void solve(monsta::Theory &theory, LATfield2::Field< std::complex<double> > &field);

    void setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize);
    void setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize, std::vector<int> skipCpts);

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
    std::vector<int> skipCpts_;
    LATfield2::Field< std::complex<double> > oldGrads_;

    void iterate(LATfield2::Field< std::complex<double> > &field, monsta::Theory &theory);
  };

  GradDescentSolver::GradDescentSolver(double tol, int maxIterations, double initialStepSize, double maxStepSize)
  : tol_(tol), maxIterations_(maxIterations), stepSize_(initialStepSize), maxStepSize_(maxStepSize) {}
  GradDescentSolver::GradDescentSolver(double tol, int maxIterations, double initialStepSize, double maxStepSize, std::vector<int> skipCpts)
  : tol_(tol), maxIterations_(maxIterations), stepSize_(initialStepSize), maxStepSize_(maxStepSize), skipCpts_(skipCpts) {}

  void GradDescentSolver::setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize)
  {
    tol_ = tol;
    maxIterations_ = maxIterations;
    stepSize_ = initialStepSize;
    maxStepSize_ = maxStepSize;
  }

  void GradDescentSolver::setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize, std::vector<int> skipCpts)
  {
    tol_ = tol;
    maxIterations_ = maxIterations;
    stepSize_ = initialStepSize;
    maxStepSize_ = maxStepSize;
    skipCpts_ = skipCpts;
  }

  void GradDescentSolver::solve(monsta::Theory &theory, LATfield2::Field< std::complex<double> > &field)
  {
    energy = theory.computeEnergy(field);
    oldGrads_.initialize(field.lattice(), field.rows(), field.cols(), field.symmetry());
    oldGrads_.alloc();

    iterate(field, theory);
    energyOld = energy;
    energy = theory.computeEnergy(field);
    relEnergyChange = (energy - energyOld);
    int numIters = 1;

    while (abs(relEnergyChange) > tol_ && numIters < maxIterations_)
    {
      numIters++;
      iterate(field, theory);
      energyOld = energy;
      energy = theory.computeEnergy(field);
      relEnergyChange = (energy - energyOld);
      if (isVerbose_)
      {
        COUT << energy << std::endl;
        // COUT << maxGrad_ << std::endl;
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

  void GradDescentSolver::iterate(LATfield2::Field< std::complex<double> > &field, monsta::Theory &theory)
  {
    clock_t begin = clock();
    LATfield2::Site site(field.lattice());
    int numRows = field.rows();
    int numCols = field.cols();
    int numMatrices = field.components() / (numRows*numCols);

    std::vector< std::complex<double> > gradVals(numRows*numCols*numMatrices);
    int vecIdx;
    monsta::Matrix gradMat(numRows);
    std::complex<double> gradVal;
    std::complex<double> oldGradVal;
    std::complex<double> fieldVal;

    double stepChangeNumerator = 0;
    double stepChangeDenominator = 0;

    double maxGrad = 0;
    for (site.first(); site.test(); site.next())
    {
      if (site.index() % 2 == 1) { continue; }
      for (int matIdx = 0; matIdx < numMatrices; matIdx++)
      {
        gradMat = theory.getLocalGradient(field, site, matIdx);
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
        for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols; colIdx++)
          {
            vecIdx = colIdx + numCols * (rowIdx + numRows * matIdx);
            gradVal = gradMat(rowIdx, colIdx);
            gradVals[vecIdx] = gradVal;
            // if (abs(gradVal) > abs(maxGrad)) { maxGrad = abs(gradVal); }
            if (matIdx == 3)
            {
              if (abs(gradVal) > abs(maxGrad)) { maxGrad = abs(gradVal); }
              oldGradVal = oldGrads_(site, rowIdx, colIdx);

              stepChangeNumerator += abs(real(oldGradVal) * (real(oldGradVal - gradVal)));
              stepChangeNumerator += abs(imag(oldGradVal) * (imag(oldGradVal - gradVal)));

              stepChangeDenominator += pow(abs(gradVal - oldGradVal),2);
              oldGrads_(site, rowIdx, colIdx) = gradVal;
            }
          }
        }
      }
      for (int matIdx = 0; matIdx < numMatrices; matIdx++)
      {
        for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols; colIdx++)
          {
            vecIdx = colIdx + numCols * (rowIdx + numRows * matIdx);
            fieldVal = field(site, matIdx, rowIdx, colIdx);
            field(site, matIdx, rowIdx, colIdx) = fieldVal - stepSize_*gradVals[vecIdx];
          }
        }
        theory.postProcess(field, site, matIdx);
      }
    }


    for (site.first(); site.test(); site.next())
    {
      if (site.index() % 2 == 0) { continue; }
      for (int matIdx = 0; matIdx < numMatrices; matIdx++)
      {
        gradMat = theory.getLocalGradient(field, site, matIdx);
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
        for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols; colIdx++)
          {
            vecIdx = colIdx + numCols * (rowIdx + numRows * matIdx);
            gradVal = gradMat(rowIdx, colIdx);
            gradVals[vecIdx] = gradVal;
            // if (abs(gradVal) > abs(maxGrad)) { maxGrad = abs(gradVal); }
            if (matIdx == 3)
            {
              if (abs(gradVal) > abs(maxGrad)) { maxGrad = abs(gradVal); }
              oldGradVal = oldGrads_(site, rowIdx, colIdx);

              stepChangeNumerator += abs(real(oldGradVal) * (real(oldGradVal - gradVal)));
              stepChangeNumerator += abs(imag(oldGradVal) * (imag(oldGradVal - gradVal)));

              stepChangeDenominator += pow(abs(gradVal - oldGradVal),2);
              oldGrads_(site, rowIdx, colIdx) = gradVal;
            }
          }
        }
      }
      for (int matIdx = 0; matIdx < numMatrices; matIdx++)
      {
        for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols; colIdx++)
          {
            vecIdx = colIdx + numCols * (rowIdx + numRows * matIdx);
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

    clock_t end = clock();
    // COUT << double(end - begin) / CLOCKS_PER_SEC << endl;
  }
}

#endif