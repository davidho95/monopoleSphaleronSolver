#ifndef PHI4THEORY_HPP
#define PHI4THEORY_HPP

#include "LATfield2.hpp"
#include "Theory.hpp"
#include <cmath>
#include <cstdio>

namespace monsta {
  class Phi4Theory: public Theory {
  public:
    Phi4Theory(double vev, double selfCoupling);

    bool validateField(LATfield2::Field<double> &field) const;

  private:
    double vev_;
    double selfCoupling_;

    double getLocalEnergyDensity(LATfield2::Field<double> &field, LATfield2::Site &site) const;
    double getLocalGradient(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt) const;
  };

  Phi4Theory::Phi4Theory(double vev, double selfCoupling)
  : Theory(), vev_(vev), selfCoupling_(selfCoupling) {}

  bool Phi4Theory::validateField(LATfield2::Field<double> &field) const
  {
    int numFieldCpts = field.components();
    if (numFieldCpts != 2 * numSpatialCpts_) {
      std::cout << "Incorrect number of spatial components." << std::endl;
      std::cout << "Required: " << numSpatialCpts_ << ", Provided: " << numFieldCpts / 2 << endl;
      return false;
    }
    return true;
  }

  double Phi4Theory::getLocalEnergyDensity(LATfield2::Field<double> &field, LATfield2::Site &site) const
  {
    double fieldVal = field(site, cpt);
    double E = selfCoupling_ * pow(pow(fieldVal, 2) - pow(vev_, 2), 2);
    return E;
  }

  double Phi4Theory::getLocalGradient(LATfield2::Field<double> &field, LATfield2::Site &site, int cpt) const
  {
    double fieldVal = field(site, cpt);
    double grad = 4 * selfCoupling_ * fieldVal*(pow(fieldVal, 2) - pow(vev_, 2));
    return grad;
  }
}

#endif