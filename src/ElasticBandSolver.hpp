#ifndef ELASTICBANDSOLVER_HPP
#define ELASTICBANDSOLVER_HPP

#include <ctime>

namespace monsta
{
  class ElasticBandSolver
  {
  public:
    ElasticBandSolver(double tol, int maxIterations, double initialStepSize, double maxStepSize, std::vector<double> springConstants);

    void setVerbosity(bool isVerbose) { isVerbose_ = isVerbose; };

    bool solve(monsta::Theory &theory, std::vector< LATfield2::Field< std::complex<double> >* > &fieldVec);

    void setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize);

  private:
    std::vector<double> springConstants_;
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

    double getTotalEnergy(std::vector< LATfield2::Field< std::complex<double> >* > &fieldVec, monsta::Theory &theory);
    void iterateSingleField(std::vector< LATfield2::Field< std::complex<double> >* > &fieldVec, int fieldIdx, monsta::Theory &theory);
  };

  ElasticBandSolver::ElasticBandSolver(double tol, int maxIterations, double initialStepSize, double maxStepSize, std::vector<double> springConstants)
  : tol_(tol), maxIterations_(maxIterations), stepSize_(initialStepSize), maxStepSize_(maxStepSize), springConstants_(springConstants) {}

  void ElasticBandSolver::setParams(double tol, int maxIterations, double initialStepSize, double maxStepSize)
  {
    tol_ = tol;
    maxIterations_ = maxIterations;
    stepSize_ = initialStepSize;
    maxStepSize_ = maxStepSize;
  }

  bool ElasticBandSolver::solve(monsta::Theory &theory, std::vector< LATfield2::Field< std::complex<double> >* > &fieldVec)
  {
    int numFields = fieldVec.size();
    energy = getTotalEnergy(fieldVec, theory);

    for (int ii = 0; ii < numFields; ii++)
    {
      if (ii == 0 || ii == numFields - 1) { continue; }
      if (ii % 2 == 0) { continue; }
      iterateSingleField(fieldVec, ii, theory);
    }

    for (int ii = 0; ii < numFields; ii++)
    {
      if (ii == 0 || ii == numFields - 1) { continue; }
      if (ii % 2 == 1) { continue; }
      iterateSingleField(fieldVec, ii, theory);
    }

    energyOld = energy;
    energy = getTotalEnergy(fieldVec, theory);
    relEnergyChange = (energy - energyOld);
    int numIters = 1;

    while (abs(relEnergyChange) > tol_ && numIters < maxIterations_)
    {
      numIters++;
      for (int ii = 0; ii < numFields; ii++)
      {
        if (ii == 0 || ii == numFields - 1) { continue; }
        if (ii % 2 == 0) { continue; }
        iterateSingleField(fieldVec, ii, theory);
      }

      for (int ii = 0; ii < numFields; ii++)
      {
        if (ii == 0 || ii == numFields - 1) { continue; }
        if (ii % 2 == 1) { continue; }
        iterateSingleField(fieldVec, ii, theory);
      }
      energyOld = energy;
      energy = getTotalEnergy(fieldVec, theory);
      relEnergyChange = (energy - energyOld);
      if (isVerbose_)
      {
        COUT << energy << std::endl;
        // COUT << maxGrad_ << std::endl;
      }
    }

    double finalEnergy = getTotalEnergy(fieldVec, theory);

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

  double ElasticBandSolver::getTotalEnergy(std::vector< LATfield2::Field< std::complex<double> >* > &fieldVec, monsta::Theory &theory)
  {
    double totalEnergy = 0;

    for (int ii = 0; ii < fieldVec.size(); ii++)
    {
      totalEnergy += theory.computeEnergy(*fieldVec[ii]);
    }

    return totalEnergy;
  }

  void ElasticBandSolver::iterateSingleField(std::vector< LATfield2::Field< std::complex<double> >* > &fieldVec, int fieldIdx, monsta::Theory &theory)
  {
    LATfield2::Site site((*fieldVec[fieldIdx]).lattice());
    int numRows = (*fieldVec[fieldIdx]).rows();
    int numCols = (*fieldVec[fieldIdx]).cols();
    int numMatrices = (*fieldVec[fieldIdx]).components() / (numRows*numCols);

    std::vector< std::complex<double> > gradVals(numRows*numCols*numMatrices);
    int vecIdx;
    monsta::Matrix gradMat(numRows);
    std::complex<double> gradVal;
    std::complex<double> fieldVal;

    double maxGrad = 0;
    for (site.first(); site.test(); site.next())
    {
      if (site.index() % 2 == 1) { continue; }
      for (int matIdx = 0; matIdx < numMatrices; matIdx++)
      {
        gradMat = theory.getLocalGradient((*fieldVec[fieldIdx]), site, matIdx);
        for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols; colIdx++)
          {
            vecIdx = colIdx + numCols * (rowIdx + numRows * matIdx);
            gradVal = gradMat(rowIdx, colIdx);
            if (matIdx == 3)
            {
              if (abs(gradVal) > maxGrad)
              {
                maxGrad = abs(gradVal);
              } 
            }
            gradVals[vecIdx] = gradVal;
          }
        }
      }
      for (int matIdx = 0; matIdx < numMatrices; matIdx++)
      {
        for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols; colIdx++)
          {
            // Standard gradient descent update
            vecIdx = colIdx + numCols * (rowIdx + numRows * matIdx);
            fieldVal = (*fieldVec[fieldIdx])(site, matIdx, rowIdx, colIdx);
            (*fieldVec[fieldIdx])(site, matIdx, rowIdx, colIdx) = fieldVal - stepSize_*gradVals[vecIdx];

            // Add spring force
            std::complex<double> springForce = 0;
            springForce += 2*springConstants_[fieldIdx]*fieldVal;
            springForce -= springConstants_[fieldIdx]*(*fieldVec[fieldIdx - 1])(site, matIdx, rowIdx, colIdx);
            springForce -= springConstants_[fieldIdx]*(*fieldVec[fieldIdx + 1])(site, matIdx, rowIdx, colIdx);

            (*fieldVec[fieldIdx])(site, matIdx, rowIdx, colIdx) -= stepSize_*springForce;
          }
        }
        theory.postProcess((*fieldVec[fieldIdx]), site, matIdx);
      }
    }

    for (site.first(); site.test(); site.next())
    {
      if (site.index() % 2 == 0) { continue; }
      for (int matIdx = 0; matIdx < numMatrices; matIdx++)
      {
        gradMat = theory.getLocalGradient((*fieldVec[fieldIdx]), site, matIdx);
        for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols; colIdx++)
          {
            vecIdx = colIdx + numCols * (rowIdx + numRows * matIdx);
            gradVal = gradMat(rowIdx, colIdx);
            if (matIdx == 3)
            {
              if (abs(gradVal) > maxGrad)
              {
                maxGrad = abs(gradVal);
              } 
            }
            gradVals[vecIdx] = gradVal;
          }
        }
      }
      for (int matIdx = 0; matIdx < numMatrices; matIdx++)
      {
        for (int rowIdx = 0; rowIdx < numRows; rowIdx++)
        {
          for (int colIdx = 0; colIdx < numCols; colIdx++)
          {
            // Standard gradient descent update
            vecIdx = colIdx + numCols * (rowIdx + numRows * matIdx);
            fieldVal = (*fieldVec[fieldIdx])(site, matIdx, rowIdx, colIdx);
            (*fieldVec[fieldIdx])(site, matIdx, rowIdx, colIdx) = fieldVal - stepSize_*gradVals[vecIdx];

            // Add spring force
            std::complex<double> springForce = 0;
            springForce += 2*springConstants_[fieldIdx]*fieldVal;
            springForce -= springConstants_[fieldIdx]*(*fieldVec[fieldIdx - 1])(site, matIdx, rowIdx, colIdx);
            springForce -= springConstants_[fieldIdx]*(*fieldVec[fieldIdx + 1])(site, matIdx, rowIdx, colIdx);

            (*fieldVec[fieldIdx])(site, matIdx, rowIdx, colIdx) -= stepSize_*springForce;
          }
        }
        theory.postProcess((*fieldVec[fieldIdx]), site, matIdx);
      }
    }
    theory.applyBoundaryConditions((*fieldVec[fieldIdx]));

    maxGrad_ = maxGrad;
  }
}

#endif