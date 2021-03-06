#ifndef SU2TOOLS_HPP
#define SU2TOOLS_HPP

#include "Matrix.hpp"

namespace monsta
{
  const monsta::Matrix identity({1, 0, 0, 1});
  const monsta::Matrix pauli1({0, 1, 1, 0});
  const monsta::Matrix pauli2({0, -1i, 1i, 0});
  const monsta::Matrix pauli3({1, 0, 0, -1});

  const double pi = 4*atan(1);

  int sign(double x)
  {
    return (x > 0) - (x < 0);
  }

  bool su2Check(monsta::Matrix mat)
  {
    double zeroTol = 1e-5;
    std::complex<double> determinant = mat(0,0)*mat(1,1) - mat(0,1)*mat(1,0);
    if (abs(abs(determinant) - 1.0) > zeroTol)
    {
      cout << abs(determinant) - 1.0 << endl;
      return false;
    }
    monsta::Matrix hcProduct = mat*monsta::conjugateTranspose(mat);
    monsta::Matrix checkMatrix = hcProduct - identity;
    
    for (int ii = 0; ii < pow(checkMatrix.getSize(),2); ii++)
    {
      if (abs(checkMatrix(ii)) > zeroTol)
      {
        return false;
      }
    }
    return true;
  }

  bool su2LieAlgCheck(monsta::Matrix mat)
  {
    double zeroTol = 1e-15;
    std::complex<double> trace = monsta::trace(mat);
    if (abs(trace) > zeroTol)
    {
      return false;
    }
    monsta::Matrix checkMatrix = monsta::conjugateTranspose(mat) - mat;
    for (int ii = 0; ii < pow(checkMatrix.getSize(),2); ii++)
    {
      if (abs(checkMatrix(ii)) > zeroTol)
      {
        return false;
      }
    }
    return true;
  }

  void makeTraceless(monsta::Matrix &mat)
  {
    std::complex<double> tr = trace(mat);
    mat(0,0) = mat(0,0) - 0.5*tr;
    mat(1,1) = mat(1,1) - 0.5*tr;
  }

  monsta::Matrix vecToSu2LieAlg(std::vector<double> vec)
  {
    return vec[0]*pauli1 + vec[1]*pauli2+vec[2]*pauli3;
  }

  monsta::Matrix vecToSu2(std::vector<double> vec)
  {
    double zeroTol = 1e-15;
    double vecNorm = sqrt(pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2));

    if (vecNorm < zeroTol)
    {
      return identity;
    }

    return cos(vecNorm)*identity + 1i*sin(vecNorm)/vecNorm * vecToSu2LieAlg(vec);
  }

  std::vector<double> su2LieAlgToVec(monsta::Matrix mat)
  {
    // if (!su2LieAlgCheck(mat))
    // { 
    //   mat.print();
    //   throw std::invalid_argument("Matrix is not an element of the SU(2) Lie algebra");
    // }
    std::vector<double> outputVec(3);
    outputVec[0] = 0.5*real(trace(mat*pauli1));
    outputVec[1] = 0.5*real(trace(mat*pauli2));
    outputVec[2] = 0.5*real(trace(mat*pauli3));

    return outputVec;
  }

  std::vector<double> su2ToVec(monsta::Matrix mat)
  {
    // if (!su2Check(mat))
    // {
    //   throw std::invalid_argument("Matrix is not an element of SU(2)");
    // }
    std::vector<double> outputVec(3);

    double zeroTol = 1e-15;

    double cosVecNorm = 0.5*real(trace(mat));
    // cout << cosVecNorm << endl;
    if (abs(cosVecNorm - 1) < zeroTol)
    {
      return outputVec;
    }
    if (abs(cosVecNorm + 1) < zeroTol)
    {
      outputVec[2] += 4*atan(1);
      return outputVec;
    }
    double vecNorm = acos(cosVecNorm);

    outputVec[0] = 0.5 * vecNorm/sin(vecNorm) * imag(trace(mat*pauli1));
    outputVec[1] = 0.5 * vecNorm/sin(vecNorm) * imag(trace(mat*pauli2));
    outputVec[2] = 0.5 * vecNorm/sin(vecNorm) * imag(trace(mat*pauli3));

    return outputVec;

  }

  int epsilonTensor(int ii, int jj, int kk)
  {
    if ((ii == jj) || (jj == kk) || (kk == ii)) { return 0; }
    if ((ii + 1) % 3 == jj && (jj + 1) % 3 == kk) { return 1; }
    else { return -1; } 
  }

  Matrix commutator(Matrix mat1, Matrix mat2)
  {
    return mat1*mat2 - mat2*mat1;
  }
}

#endif