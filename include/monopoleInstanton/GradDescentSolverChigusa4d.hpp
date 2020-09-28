#ifndef GRADDESCENTSOLVERCHIGUSA4D_HPP
#define GRADDESCENTSOLVERCHIGUSA4D_HPP

#include <ctime>
#include "../MonopoleFieldTools.hpp"

namespace monsta
{
  class GradDescentSolverChigusa
  {
  public:
    GradDescentSolverChigusa(double tol, int maxIterations, double initialStepSize, double maxStepSize, double correctionCoeff);
    GradDescentSolverChigusa(double tol, int maxIterations, double initialStepSize, double maxStepSize, double correctionCoeff, double abortGrad);
    GradDescentSolverChigusa(double tol, int maxIterations, double initialStepSize, double maxStepSize, double correctionCoeff, double abortGrad, int minSteps);
    GradDescentSolverChigusa(double tol, int maxIterations, double initialStepSize, double maxStepSize, double correctionCoeff, double abortGrad, int minSteps, bool minimiseGrads);

    void setVerbosity(bool isVerbose) { isVerbose_ = isVerbose; };

    bool solve(monsta::Theory &theory, LATfield2::Field< std::complex<double> > &field, LATfield2::Field< std::complex<double> > &referenceField);

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
    double abortGrad_ = 1e6;
    double correctionCoeff_ = 2;
    LATfield2::Field< std::complex<double> > oldGrads_;
    int minSteps_ = 0;
    bool minimiseGrads_ = false;

    void iterate(LATfield2::Field< std::complex<double> > &field, monsta::Theory &theory, LATfield2::Field< std::complex<double> > &referenceField);
  };

  GradDescentSolverChigusa::GradDescentSolverChigusa(double tol, int maxIterations, double initialStepSize, double maxStepSize, double correctionCoeff)
  : tol_(tol), maxIterations_(maxIterations), stepSize_(initialStepSize), maxStepSize_(maxStepSize), correctionCoeff_(correctionCoeff) {}
  GradDescentSolverChigusa::GradDescentSolverChigusa(double tol, int maxIterations, double initialStepSize, double maxStepSize, double correctionCoeff, double abortGrad)
  : tol_(tol), maxIterations_(maxIterations), stepSize_(initialStepSize), maxStepSize_(maxStepSize), correctionCoeff_(correctionCoeff), abortGrad_(abortGrad) {}
  GradDescentSolverChigusa::GradDescentSolverChigusa(double tol, int maxIterations, double initialStepSize, double maxStepSize, double correctionCoeff, double abortGrad, int minSteps)
  : tol_(tol), maxIterations_(maxIterations), stepSize_(initialStepSize), maxStepSize_(maxStepSize), correctionCoeff_(correctionCoeff), abortGrad_(abortGrad), minSteps_(minSteps) {}
  GradDescentSolverChigusa::GradDescentSolverChigusa(double tol, int maxIterations, double initialStepSize, double maxStepSize, double correctionCoeff, double abortGrad, int minSteps, bool minimiseGrads)
  : tol_(tol), maxIterations_(maxIterations), stepSize_(initialStepSize), maxStepSize_(maxStepSize), correctionCoeff_(correctionCoeff), abortGrad_(abortGrad), minSteps_(minSteps), minimiseGrads_(minimiseGrads) {}

  void GradDescentSolverChigusa::setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize)
  {
    tol_ = tol;
    maxIterations_ = maxIterations;
    stepSize_ = initialStepSize;
    maxStepSize_ = maxStepSize;
  }

  bool GradDescentSolverChigusa::solve(monsta::Theory &theory, LATfield2::Field< std::complex<double> > &field, LATfield2::Field< std::complex<double> > &referenceField)
  {
    energy = theory.computeEnergy(field);
    oldGrads_.initialize(field.lattice(), field.rows(), field.cols(), field.symmetry());
    oldGrads_.alloc();

    iterate(field, theory, referenceField);
    energyOld = energy;
    energy = theory.computeEnergy(field);
    relEnergyChange = (energy - energyOld);
    int numIters = 1;

    while (numIters < minSteps_)
    {
      numIters++;
      iterate(field, theory, referenceField);
      energyOld = energy;
      energy = theory.computeEnergy(field);
      relEnergyChange = (energy - energyOld);
      if (isVerbose_)
      {
        COUT << std::fixed;
        COUT << std::setprecision(10);
        COUT << energy << std::endl;
        COUT << maxGrad_ << std::endl;
        // COUT << stepSize_ << std::endl;
      }
    }

    double maxGradOld = 1e6;
    while (maxGrad_ > tol_ && numIters < maxIterations_ && maxGrad_ < abortGrad_)
    {
      if (minimiseGrads_ && maxGrad_ > maxGradOld) { break; }
      maxGradOld = maxGrad_; 
      numIters++;
      iterate(field, theory, referenceField);
      energyOld = energy;
      energy = theory.computeEnergy(field);
      relEnergyChange = (energy - energyOld);
      if (isVerbose_)
      {
        COUT << std::fixed;
        COUT << std::setprecision(10);
        COUT << energy << std::endl;
        COUT << maxGrad_ << std::endl;
        // COUT << stepSize_ << std::endl;
      }
    }

    double finalEnergy = theory.computeEnergy(field);

    if (maxGrad_ > abortGrad_)
    {
      COUT << "Gradient descent aborted after exceeding maximum allowed gradient" << std::endl;
      // COUT << "Maximum gradient: " << maxGrad_ << std::endl;
      COUT << "Energy reached: " << finalEnergy << std::endl;
      return false;
    }
    if (maxGrad_ < tol_)
    {
      COUT << "Gradient descent finished in " << numIters << " iterations." << std::endl;
      COUT << "Minimum energy: " << finalEnergy << std::endl;
      return true;
    }
    else
    {
      COUT << "Gradient descent aborted after " << numIters << " iterations." << std::endl;
      // COUT << "Maximum gradient: " << maxGrad_ << std::endl;
      COUT << "Energy reached: " << finalEnergy << std::endl;
      return false;
    }
  }

  void GradDescentSolverChigusa::iterate(LATfield2::Field< std::complex<double> > &field, monsta::Theory &theory, LATfield2::Field< std::complex<double> > &referenceField)
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

    LATfield2::Field< std::complex<double> > gradField(field.lattice(), numMatrices, numRows, numCols, field.symmetry());
    LATfield2::Field< std::complex<double> > refFieldProj(field.lattice(), numMatrices, numRows, numCols, field.symmetry());


    for (site.first(); site.test(); site.next())
    {
      for (int matIdx = 0; matIdx < numMatrices; matIdx++)
      {
        monsta::Matrix projMat(referenceField, site, matIdx);
        monsta::Matrix gradMat = theory.getLocalGradient(field, site, matIdx);
        if (matIdx < 3)
        {
          monsta::Matrix fieldMat(field, site, matIdx);
          projMat = projMat - 0.5*real(trace(projMat*conjugateTranspose(fieldMat)))*fieldMat;
          gradMat = gradMat - 0.5*real(trace(gradMat*conjugateTranspose(fieldMat)))*fieldMat;
        }
        for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols; colIdx++)
          {
            gradField(site, matIdx, rowIdx, colIdx) = gradMat(rowIdx, colIdx);
            refFieldProj(site, matIdx, rowIdx, colIdx) = projMat(rowIdx, colIdx);
          }
        }
      }
    }

    double gradDotRef = monsta::innerProduct(gradField, refFieldProj);
    site.first();
    // cout << gradDotRef << endl;

    double maxGrad = 0;
    for (site.first(); site.test(); site.next())
    {
      if (site.index() % 2 == 1) { continue; }
      for (int matIdx = 0; matIdx < numMatrices; matIdx++)
      {
        gradMat = Matrix(gradField, site, matIdx);
        monsta::Matrix fieldMat(field, site, matIdx);
        for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols; colIdx++)
          {
            gradVal = gradMat(rowIdx, colIdx);
            if (abs(gradVal) > abs(maxGrad)) { maxGrad = abs(gradVal); }
          }
        }
        gradMat = gradMat - correctionCoeff_*gradDotRef*monsta::Matrix(referenceField, site, matIdx);
        for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols; colIdx++)
          {
            vecIdx = colIdx + numCols * (rowIdx + numRows * matIdx);
            gradVal = gradMat(rowIdx, colIdx);
            gradVals[vecIdx] = gradVal;
            if (matIdx == 3)
            {
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
        gradMat = Matrix(gradField, site, matIdx);
        monsta::Matrix fieldMat(field, site, matIdx);
        for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols; colIdx++)
          {
            gradVal = gradMat(rowIdx, colIdx);
            if (abs(gradVal) > abs(maxGrad)) { maxGrad = abs(gradVal); }
          }
        }
        gradMat = gradMat - correctionCoeff_*gradDotRef*monsta::Matrix(referenceField, site, matIdx);
        for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols; colIdx++)
          {
            vecIdx = colIdx + numCols * (rowIdx + numRows * matIdx);
            gradVal = gradMat(rowIdx, colIdx);
            gradVals[vecIdx] = gradVal;
            if (matIdx == 3)
            {
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
    parallel.max(maxGrad);
    maxGrad_ = maxGrad;

    clock_t end = clock();
    // COUT << double(end - begin) / CLOCKS_PER_SEC << endl;
  }
}

#endif
