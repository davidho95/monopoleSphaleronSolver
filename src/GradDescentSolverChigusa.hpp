#ifndef GRADDESCENTSOLVERCHIGUSA
#define GRADDESCENTSOLVERCHIGUSA

#include <ctime>
#include "MonopoleFieldTools.hpp"

namespace monsta
{
  class GradDescentSolverChigusa
  {
  public:
    GradDescentSolverChigusa(double tol, int maxIterations, double initialStepSize, double maxStepSize);
    GradDescentSolverChigusa(double tol, int maxIterations, double initialStepSize, double maxStepSize, std::vector<int> skipCpts);

    void setVerbosity(bool isVerbose) { isVerbose_ = isVerbose; };

    bool solve(monsta::GeorgiGlashowSu2Theory &theory, LATfield2::Field< std::complex<double> > &field, LATfield2::Field< std::complex<double> > &referenceField);

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
    std::vector<int> coreIndices_;
    std::vector<int> skipCpts_;
    LATfield2::Field< std::complex<double> > oldGrads_;

    void iterate(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2Theory &theory, LATfield2::Field< std::complex<double> > &referenceField);
    std::vector<int> getMonopoleNums(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, monsta::GeorgiGlashowSu2Theory &theory);
    bool coreCheck(LATfield2::Site &site);
  };

  GradDescentSolverChigusa::GradDescentSolverChigusa(double tol, int maxIterations, double initialStepSize, double maxStepSize)
  : tol_(tol), maxIterations_(maxIterations), stepSize_(initialStepSize), maxStepSize_(maxStepSize) {}
  GradDescentSolverChigusa::GradDescentSolverChigusa(double tol, int maxIterations, double initialStepSize, double maxStepSize, std::vector<int> skipCpts)
  : tol_(tol), maxIterations_(maxIterations), stepSize_(initialStepSize), maxStepSize_(maxStepSize), skipCpts_(skipCpts) {}

  void GradDescentSolverChigusa::setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize)
  {
    tol_ = tol;
    maxIterations_ = maxIterations;
    stepSize_ = initialStepSize;
    maxStepSize_ = maxStepSize;
  }

  void GradDescentSolverChigusa::setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize, std::vector<int> skipCpts)
  {
    tol_ = tol;
    maxIterations_ = maxIterations;
    stepSize_ = initialStepSize;
    maxStepSize_ = maxStepSize;
    skipCpts_ = skipCpts;
  }

  bool GradDescentSolverChigusa::solve(monsta::GeorgiGlashowSu2Theory &theory, LATfield2::Field< std::complex<double> > &field, LATfield2::Field< std::complex<double> > &referenceField)
  {
    energy = theory.computeEnergy(field);
    oldGrads_.initialize(field.lattice(), field.rows(), field.cols(), field.symmetry());
    oldGrads_.alloc();

    iterate(field, theory, referenceField);
    energyOld = energy;
    energy = theory.computeEnergy(field);
    relEnergyChange = (energy - energyOld);
    int numIters = 1;

    LATfield2::Site site(field.lattice());
    for (site.first(); site.test(); site.next())
    {
      LATfield2::Site tempSite(field.lattice());
      int monopoleNum = theory.getMonopoleNumber(field, site);
      if (monopoleNum != 0)
      {
        coreIndices_.push_back(site.index());
      }
      tempSite = site + 0;
      coreIndices_.push_back(tempSite.index());
      tempSite = site + 1;
      coreIndices_.push_back(tempSite.index());
      tempSite = site + 2;
      coreIndices_.push_back(tempSite.index());
      tempSite = site + 0 + 1;
      coreIndices_.push_back(tempSite.index());
      tempSite = site + 1 + 2;
      coreIndices_.push_back(tempSite.index());
      tempSite = site + 2 + 0;
      coreIndices_.push_back(tempSite.index());
      tempSite = site + 0 + 1 + 2;
      coreIndices_.push_back(tempSite.index());
    }

    while (maxGrad_ > tol_ && numIters < maxIterations_)
    {
      numIters++;
      iterate(field, theory, referenceField);
      energyOld = energy;
      energy = theory.computeEnergy(field);
      relEnergyChange = (energy - energyOld);
      if (isVerbose_)
      {
        COUT << energy << std::endl;
        COUT << maxGrad_ << std::endl;
        // COUT << stepSize_ << std::endl;
      }
    }

    double finalEnergy = theory.computeEnergy(field);

    if (numIters < maxIterations_) {
      COUT << "Gradient descent finished in " << numIters << " iterations." << std::endl;
      COUT << "Minimum energy: " << finalEnergy << std::endl;
      return true;
    } else {
      COUT << "Gradient descent aborted after " << maxIterations_ << " iterations." << std::endl;
      // COUT << "Maximum gradient: " << maxGrad_ << std::endl;
      COUT << "Energy reached: " << finalEnergy << std::endl;
      return false;
    }
  }

  void GradDescentSolverChigusa::iterate(LATfield2::Field< std::complex<double> > &field, monsta::GeorgiGlashowSu2Theory &theory, LATfield2::Field< std::complex<double> > &referenceField)
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
            if (abs(gradVal) > abs(maxGrad) && !coreCheck(site)) { maxGrad = abs(gradVal); }
          }
        }
        gradMat = gradMat - 2*gradDotRef*monsta::Matrix(referenceField, site, matIdx);
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
            if (matIdx == 3)
            {
              // if (abs(gradVal) > abs(maxGrad) && !coreCheck(site)) { maxGrad = abs(gradVal); }
              oldGradVal = oldGrads_(site, rowIdx, colIdx);

              stepChangeNumerator += abs(real(oldGradVal) * (real(oldGradVal - gradVal)));
              stepChangeNumerator += abs(imag(oldGradVal) * (imag(oldGradVal - gradVal)));

              stepChangeDenominator += pow(abs(gradVal - oldGradVal),2);
              oldGrads_(site, rowIdx, colIdx) = gradVal;
            }
          }
        }
      }
      std::vector<int> monopoleNums = getMonopoleNums(field, site, theory);
      std::vector< std::complex<double> > oldFieldVals(numRows*numCols*numMatrices);
      for (int matIdx = 0; matIdx < numMatrices; matIdx++)
      {
        for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols; colIdx++)
          {
            vecIdx = colIdx + numCols * (rowIdx + numRows * matIdx);
            fieldVal = field(site, matIdx, rowIdx, colIdx);
            field(site, matIdx, rowIdx, colIdx) = fieldVal - stepSize_*gradVals[vecIdx];
            oldFieldVals[vecIdx] = fieldVal;
          }
        }
        theory.postProcess(field, site, matIdx);
      }


      bool monopoleNumChanged = false;
      std::vector<int> newMonopoleNums = getMonopoleNums(field, site, theory);
      for (int ii = 0; ii < monopoleNums.size(); ii++)
      {
        if (monopoleNums[ii] != newMonopoleNums[ii])
        {
          monopoleNumChanged = true;
        }
      }

      // if (monopoleNumChanged)
      // {
      //   for (int matIdx = 0; matIdx < 3; matIdx++)
      //   {
      //     for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
      //     {
      //       for (int colIdx = 0; colIdx < numCols; colIdx++)
      //       {
      //         vecIdx = colIdx + numCols * (rowIdx + numRows * matIdx);
      //         field(site, matIdx, rowIdx, colIdx) = oldFieldVals[vecIdx];// + stepSize_*gradVals[vecIdx];
      //       }
      //     }
      //     theory.postProcess(field, site, matIdx);
      //   }
      // }
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
            if (abs(gradVal) > abs(maxGrad) && !coreCheck(site)) { maxGrad = abs(gradVal); }
          }
        }
        gradMat = gradMat - 2*gradDotRef*monsta::Matrix(referenceField, site, matIdx);
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
            if (matIdx == 3)
            {
              // if (abs(gradVal) > abs(maxGrad) && !coreCheck(site)) { maxGrad = abs(gradVal); }
              oldGradVal = oldGrads_(site, rowIdx, colIdx);

              stepChangeNumerator += abs(real(oldGradVal) * (real(oldGradVal - gradVal)));
              stepChangeNumerator += abs(imag(oldGradVal) * (imag(oldGradVal - gradVal)));

              stepChangeDenominator += pow(abs(gradVal - oldGradVal),2);
              oldGrads_(site, rowIdx, colIdx) = gradVal;
            }
          }
        }
      }
      std::vector<int> monopoleNums = getMonopoleNums(field, site, theory);
      std::vector< std::complex<double> > oldFieldVals(numRows*numCols*numMatrices);
      for (int matIdx = 0; matIdx < numMatrices; matIdx++)
      {
        for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols; colIdx++)
          {
            vecIdx = colIdx + numCols * (rowIdx + numRows * matIdx);
            fieldVal = field(site, matIdx, rowIdx, colIdx);
            field(site, matIdx, rowIdx, colIdx) = fieldVal - stepSize_*gradVals[vecIdx];
            oldFieldVals[vecIdx] = fieldVal;
          }
        }
        theory.postProcess(field, site, matIdx);
      }


      bool monopoleNumChanged = false;
      std::vector<int> newMonopoleNums = getMonopoleNums(field, site, theory);
      for (int ii = 0; ii < monopoleNums.size(); ii++)
      {
        if (monopoleNums[ii] != newMonopoleNums[ii])
        {
          monopoleNumChanged = true;
        }
      }

      // if (monopoleNumChanged)
      // {
      //   for (int matIdx = 0; matIdx < 3; matIdx++)
      //   {
      //     for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
      //     {
      //       for (int colIdx = 0; colIdx < numCols; colIdx++)
      //       {
      //         vecIdx = colIdx + numCols * (rowIdx + numRows * matIdx);
      //         field(site, matIdx, rowIdx, colIdx) = oldFieldVals[vecIdx];// + stepSize_*gradVals[vecIdx];
      //       }
      //     }
      //     theory.postProcess(field, site, matIdx);
      //   }
      // }
    }

    parallel.sum(stepChangeNumerator);
    parallel.sum(stepChangeDenominator);

    double stepChange = stepChangeNumerator / stepChangeDenominator;
    if (stepChange > 1e-8) { stepSize_ *= stepChange; }
    if (stepSize_ > maxStepSize_) { stepSize_ = maxStepSize_; }

    theory.applyBoundaryConditions(field);
    // double stepChange = system.updateGradientReturnStepChange();
    // stepSize_ *= stepChange;
    parallel.max(maxGrad);
    maxGrad_ = maxGrad;

    clock_t end = clock();
    // COUT << double(end - begin) / CLOCKS_PER_SEC << endl;
  }

  std::vector<int> GradDescentSolverChigusa::getMonopoleNums(LATfield2::Field< std::complex<double> > &field, LATfield2::Site &site, monsta::GeorgiGlashowSu2Theory &theory)
  {
    std::vector<int> monopoleNums(7);

    LATfield2::Site shiftedSite(field.lattice());

    monopoleNums[0] = theory.getMonopoleNumber(field, site);
    shiftedSite = site - 0;
    monopoleNums[1] = theory.getMonopoleNumber(field, shiftedSite);
    shiftedSite = shiftedSite - 1;
    monopoleNums[2] = theory.getMonopoleNumber(field, shiftedSite);
    shiftedSite = shiftedSite + 0;
    monopoleNums[3] = theory.getMonopoleNumber(field, shiftedSite);
    shiftedSite = shiftedSite + 1;
    shiftedSite = shiftedSite - 2;
    monopoleNums[4] = theory.getMonopoleNumber(field, shiftedSite);
    shiftedSite = shiftedSite - 0;
    monopoleNums[5] = theory.getMonopoleNumber(field, shiftedSite);
    shiftedSite = shiftedSite + 0;
    shiftedSite = shiftedSite - 1;
    monopoleNums[6] = theory.getMonopoleNumber(field, shiftedSite);

    return monopoleNums;
  }

  bool GradDescentSolverChigusa::coreCheck(LATfield2::Site &site)
  {
    // for (int ii = 0; ii < coreIndices_.size(); ii++)
    // {
    //   if (site.index() == coreIndices_[ii]) { return true; }
    // }
    return false;
  }
}

#endif