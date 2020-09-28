#ifndef GRADDESCENTSOLVER_HPP
#define GRADDESCENTSOLVER_HPP

#include "Theory.hpp"
#include <ctime>
#include <iomanip>

namespace monsta
{
  class GradDescentSolver
  {
  public:
    GradDescentSolver(double tol, int maxIterations, double initialStepSize, double maxStepSize);
    GradDescentSolver(double tol, int maxIterations, double initialStepSize, double maxStepSize, int minSteps);
    GradDescentSolver(double tol, int maxIterations, double initialStepSize, double maxStepSize, int minSteps, bool minimiseGrad);
    GradDescentSolver(double tol, int maxIterations, double initialStepSize, double maxStepSize, int minSteps, bool minimiseGrad, bool constantStep);

    void setVerbosity(bool isVerbose) { isVerbose_ = isVerbose; };

    bool solve(monsta::Theory &theory, LATfield2::Field< std::complex<double> > &field);

    void setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize);

  private:
    double tol_;
    int maxIterations_;
    double stepSize_;
    double maxStepSize_;
    double maxGrad_ = 1e6;
    bool isVerbose_ = true;
    int minSteps_ = 0;
    bool minimiseGrad_ = false;
    bool constantStep_ = false;
    // LATfield2::Field< std::complex<double> > oldGrads;

    void iterate(LATfield2::Field< std::complex<double> > &field, monsta::Theory &theory, LATfield2::Field< std::complex<double> > &oldGrads);
    void checkerboardUpdate(LATfield2::Field< std::complex<double> > &field, monsta::Theory &theory, LATfield2::Field< std::complex<double> > &oldGrads, bool isEven, double &stepChangeNumerator, double &stepChangeDenominator, double &maxGrad);
  };

  GradDescentSolver::GradDescentSolver(double tol, int maxIterations, double initialStepSize, double maxStepSize)
  : tol_(tol), maxIterations_(maxIterations), stepSize_(initialStepSize), maxStepSize_(maxStepSize) {}
    GradDescentSolver::GradDescentSolver(double tol, int maxIterations, double initialStepSize, double maxStepSize, int minSteps)
  : tol_(tol), maxIterations_(maxIterations), stepSize_(initialStepSize), maxStepSize_(maxStepSize), minSteps_(minSteps) {}
    GradDescentSolver::GradDescentSolver(double tol, int maxIterations, double initialStepSize, double maxStepSize, int minSteps, bool minimiseGrad)
  : tol_(tol), maxIterations_(maxIterations), stepSize_(initialStepSize), maxStepSize_(maxStepSize), minSteps_(minSteps), minimiseGrad_(minimiseGrad) {}
    GradDescentSolver::GradDescentSolver(double tol, int maxIterations, double initialStepSize, double maxStepSize, int minSteps, bool minimiseGrad, bool constantStep)
  : tol_(tol), maxIterations_(maxIterations), stepSize_(initialStepSize), maxStepSize_(maxStepSize), minSteps_(minSteps), minimiseGrad_(minimiseGrad),
    constantStep_(constantStep) {}

  void GradDescentSolver::setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize)
  {
    tol_ = tol;
    maxIterations_ = maxIterations;
    stepSize_ = initialStepSize;
    maxStepSize_ = maxStepSize;
  }

  bool GradDescentSolver::solve(monsta::Theory &theory, LATfield2::Field< std::complex<double> > &field)
  {
    double energy = theory.computeEnergy(field);
    double relEnergyChange = 1e6;
    double energyOld;
    LATfield2::Field< std::complex<double> > oldGrads(field.lattice(), field.rows(), field.cols(), field.symmetry());

    iterate(field, theory, oldGrads);
    energyOld = energy;
    energy = theory.computeEnergy(field);
    relEnergyChange = abs(energy - energyOld) / energy;
    int numIters = 1;

    if (constantStep_) { stepSize_ = maxStepSize_; }

    double maxGradOld = 0;
    if (minimiseGrad_)
    {
      while (maxGradOld < maxGrad_)
      {
        numIters++;
        maxGradOld = maxGrad_;
        iterate(field, theory, oldGrads);
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
      maxGradOld = maxGrad_;
      numIters++;
      iterate(field, theory, oldGrads);
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
    
    while (maxGrad_ > tol_ && numIters < maxIterations_)
    {
      if (minimiseGrad_ && maxGradOld < maxGrad_) { break; }
      numIters++;
      maxGradOld = maxGrad_;
      iterate(field, theory, oldGrads);
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

    double finalEnergy = theory.computeEnergy(field);

    if (maxGrad_ < tol_) {
      COUT << "Gradient descent finished in " << numIters << " iterations." << std::endl;
      COUT << "Minimum energy: " << finalEnergy << std::endl;
      return true;
    } else {
      COUT << "Gradient descent aborted after " << numIters << " iterations." << std::endl;
      COUT << "Energy reached: " << finalEnergy << std::endl;
      return false;
    }
  }

  void GradDescentSolver::iterate(LATfield2::Field< std::complex<double> > &field, monsta::Theory &theory, LATfield2::Field< std::complex<double> > &oldGrads)
  {
    double stepChangeNumerator = 0;
    double stepChangeDenominator = 0;

    double maxGrad = 0;

    bool isEven = true;
    checkerboardUpdate(field, theory, oldGrads, isEven, stepChangeNumerator, stepChangeDenominator, maxGrad);
    isEven = false;
    checkerboardUpdate(field, theory, oldGrads, isEven, stepChangeNumerator, stepChangeDenominator, maxGrad);

    parallel.sum(stepChangeNumerator);
    parallel.sum(stepChangeDenominator);

    double stepChange = stepChangeNumerator / stepChangeDenominator;
    if (stepChange > 1e-8) { stepSize_ *= stepChange; }
    if (stepSize_ > maxStepSize_) { stepSize_ = maxStepSize_; }
    if (constantStep_) { stepSize_ = maxStepSize_; }

    theory.applyBoundaryConditions(field);
    maxGrad_ = maxGrad;
    parallel.max(maxGrad_);
  }

  void GradDescentSolver::checkerboardUpdate(LATfield2::Field< std::complex<double> > &field, monsta::Theory &theory, LATfield2::Field< std::complex<double> > &oldGrads, bool isEven, double &stepChangeNumerator, double &stepChangeDenominator, double &maxGrad)
  {
    LATfield2::Site site(field.lattice());
    int numRows = field.rows();
    int numCols = field.cols();
    int numMatrices = field.components() / (numRows*numCols);

    for (site.first(); site.test(); site.next())
    {
      int vecIdx;
      std::vector< std::complex<double> > gradVals(numRows*numCols*numMatrices);
      if (site.index() % 2 == isEven) { continue; }
      for (int matIdx = 0; matIdx < numMatrices; matIdx++)
      {
        monsta::Matrix fieldMat(field, site, matIdx);
        monsta::Matrix gradMat = theory.getLocalGradient(field, site, matIdx);
        monsta::Matrix gradMatProj = gradMat;
        if (matIdx < 3)
        {
          gradMatProj = gradMat - 0.5*real(trace(gradMat*conjugateTranspose(fieldMat)))*fieldMat;
        }
        for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols; colIdx++)
          {
            vecIdx = colIdx + numCols * (rowIdx + numRows * matIdx);
            std::complex<double> gradVal = gradMat(rowIdx, colIdx);
            gradVals[vecIdx] = gradVal;
            std::complex<double> projGradVal = gradMatProj(rowIdx, colIdx);
            if (matIdx < 3)
            {
              if (abs(projGradVal) > abs(maxGrad)) { maxGrad = abs(projGradVal); }
            }
            if (matIdx == 3)
            {
              if (abs(projGradVal) > abs(maxGrad)) { maxGrad = abs(projGradVal); }
              std::complex<double> oldGradVal = oldGrads(site, rowIdx, colIdx);

              stepChangeNumerator += abs(real(oldGradVal) * (real(oldGradVal - gradVal)));
              stepChangeNumerator += abs(imag(oldGradVal) * (imag(oldGradVal - gradVal)));

              stepChangeDenominator += pow(abs(gradVal - oldGradVal),2);
              oldGrads(site, rowIdx, colIdx) = gradVal;
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
            std::complex<double> fieldVal = field(site, matIdx, rowIdx, colIdx);
            field(site, matIdx, rowIdx, colIdx) = fieldVal - stepSize_*gradVals[vecIdx];
          }
        }
        theory.postProcess(field, site, matIdx);
      }
    }
  }
}

#endif
