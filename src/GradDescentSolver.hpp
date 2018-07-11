#ifndef GRADDESCENTSOLVER_HPP
#define GRADDESCENTSOLVER_HPP

namespace monsta
{
  class GradDescentSolver
  {
  public:
    GradDescentSolver(double tol, int maxIterations);

    void setVerbosity(bool isVerbose) { isVerbose_ = isVerbose; };

    void solve(monsta::System &system);

  private:
    double tol_;
    int maxIterations_;
    double maxGrad_ = 1e6;
    double stepSize_ = 0.01;
    bool isVerbose_ = true;

    void iterate(monsta::System &system);
  };

  GradDescentSolver::GradDescentSolver(double tol, int maxIterations)
  : tol_(tol), maxIterations_(maxIterations) {}

  void GradDescentSolver::solve(monsta::System &system)
  {
    maxGrad_ = 1e6;
    stepSize_ = 0.1;
    iterate(system);
    int numIters = 1;
    while (abs(maxGrad_) > tol_ && numIters < maxIterations_)
    {
      numIters++;
      iterate(system);
      if (isVerbose_)
      {
        std::cout << maxGrad_ << std::endl;
      }
    }

    double finalEnergy = system.getEnergy();

    if (numIters < maxIterations_) {
      std::cout << "Gradient descent finished in " << numIters << " iterations." << std::endl;
      std::cout << "Minimum energy: " << finalEnergy << std::endl;
    } else {
      std::cout << "Gradient descent aborted after " << maxIterations_ << " iterations." << std::endl;
      // std::cout << "Maximum gradient: " << maxGrad_ << std::endl;
      std::cout << "Energy reached: " << finalEnergy << std::endl;
    }
  }

  void GradDescentSolver::iterate(monsta::System &system)
  {
    LATfield2::Site site(system.getLattice());
    int numSpatialCpts = system.getNumSpatialCpts();

    double maxGrad = 0;
    for (site.first(); site.test(); site.next())
    {
      for (int iCpt = 0; iCpt < numSpatialCpts; iCpt++)
      {
        double fieldVal = system.getFieldElement(site, iCpt);
        double gradVal = system.getGradientElement(site, iCpt);
        if (abs(gradVal) > abs(maxGrad)) { maxGrad = gradVal; }
        system.setFieldElement(site, iCpt, fieldVal - stepSize_*gradVal);
      }
    }
    system.updateFieldHalo();
    double stepChange = system.updateGradientReturnStepChange();
    stepSize_ *= stepChange;
    maxGrad_ = maxGrad;
  }
}

#endif