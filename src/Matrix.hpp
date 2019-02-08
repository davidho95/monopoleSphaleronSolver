#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <complex>
#include <vector>
#include <cmath>

namespace monsta
{
  class Matrix
  {
  public:
    Matrix(int size);

    int getSize() const;
    std::complex<double> operator()(int vecIdx) const;
    std::complex<double>& operator()(int vecIdx);
    std::complex<double> operator()(int colIdx, int rowIdx) const;
    std::complex<double>& operator()(int colIdx, int rowIdx);
    Matrix& operator*=(const std::complex<double> scalar);
    Matrix& operator*=(double scalar);

    void print() const;

  private:
    int size_;
    std::vector<std::complex<double> > elements_;

    int computeVecIdx(int colIdx, int rowIdx) const;
  };

  Matrix::Matrix(int size) : size_(size)
  {
    for (int ii = 0; ii < pow(size,2); ii++)
    {
      elements_.push_back(0);
    }
  }

  int Matrix::getSize() const
  {
    return size_;
  }

  int Matrix::computeVecIdx(int colIdx, int rowIdx) const
  {
    return colIdx*size_ + rowIdx;
  }

  std::complex<double> Matrix::operator()(int vecIdx) const
  {
    return elements_[vecIdx];
  }

  std::complex<double>& Matrix::operator()(int vecIdx)
  {
    return elements_[vecIdx];
  }

  std::complex<double> Matrix::operator()(int colIdx, int rowIdx) const
  {
    int vecIdx = computeVecIdx(colIdx, rowIdx);
    return elements_[vecIdx];
  }

  std::complex<double>& Matrix::operator()(int colIdx, int rowIdx)
  {
    int vecIdx = computeVecIdx(colIdx, rowIdx);
    return elements_[vecIdx];
  }

  Matrix operator*(std::complex<double> scalar, const Matrix &mat)
  { 
    int matSize = mat.getSize();
    Matrix outputMat(matSize);
    for (int ii = 0; ii < pow(matSize,2); ii++)
    {
      outputMat(ii)= scalar*mat(ii); 
    }
    return outputMat;
  }

  Matrix operator*(double scalar, const Matrix &mat)
  { 
    int matSize = mat.getSize();
    Matrix outputMat(matSize);
    for (int ii = 0; ii < pow(matSize,2); ii++)
    {
      outputMat(ii)= scalar*mat(ii); 
    }
    return outputMat;
  }

  Matrix operator/(const Matrix &mat, std::complex<double> scalar)
  { 
    int matSize = mat.getSize();
    Matrix outputMat(matSize);
    for (int ii = 0; ii < pow(matSize,2); ii++)
    {
      outputMat(ii)= mat(ii)/scalar; 
    }
    return outputMat;
  }

  Matrix operator/(const Matrix &mat, double scalar)
  { 
    int matSize = mat.getSize();
    Matrix outputMat(matSize);
    for (int ii = 0; ii < pow(matSize,2); ii++)
    {
      outputMat(ii)= mat(ii)/scalar; 
    }
    return outputMat;
  }

  Matrix& Matrix::operator*=(const std::complex<double> scalar)
  { 
    int matSize = this->size_;
    for (int ii = 0; ii < pow(matSize,2); ii++)
    {
      this->operator()(ii) *= scalar; 
    }
    return *this;
  }

  Matrix& Matrix::operator*=(double scalar)
  { 
    int matSize = this->size_;
    for (int ii = 0; ii < pow(matSize,2); ii++)
    {
      this->operator()(ii) *= scalar; 
    }
    return *this;
  }

  void Matrix::print() const
  {
    for (int ii = 0; ii < size_; ii++)
    {
      for (int jj = 0; jj < size_; jj++)
      {
        std::complex<double> element = this->operator()(ii,jj);
        COUT << element.real() << (element.imag() < 0.0 ? "-" : "+") << std::abs(element.imag()) << "i ";
      }
      COUT << endl;
    }
    COUT << endl;
  }

  Matrix operator*(const Matrix &mat1, const Matrix &mat2)
  {
    if (mat1.getSize() != mat2.getSize())
    {
      throw std::invalid_argument("Only matrices of the same size can be multplied");
    }
    int matSize = mat1.getSize();
    Matrix outputMat(matSize);
    for (int ii = 0; ii < matSize; ii++)
    {
      for (int jj = 0; jj < matSize; jj++)
      {
        std::complex<double> outputElement = 0;
        for (int kk = 0; kk < matSize; kk++)
        {
          outputElement += mat1(ii,kk)*mat2(kk,jj);
        }
        outputMat(ii,jj) = outputElement;
      }
    }
    return outputMat; 
  }

  Matrix operator+(const Matrix &mat1, const Matrix &mat2)
  { 
    if (mat1.getSize() != mat2.getSize())
    {
      throw std::invalid_argument("Only matrices of the same size can be added");
    }
    int matSize = mat1.getSize();
    Matrix outputMat(matSize);
    for (int ii = 0; ii < pow(matSize,2); ii++)
    {
      outputMat(ii) = mat1(ii) + mat2(ii); 
    }
    return outputMat;
  }

  Matrix operator-(const Matrix &mat1, const Matrix &mat2)
  { 
    if (mat1.getSize() != mat2.getSize())
    {
      throw std::invalid_argument("Only matrices of the same size can be subtracted");
    }
    int matSize = mat1.getSize();
    Matrix outputMat(matSize);
    for (int ii = 0; ii < pow(matSize,2); ii++)
    {
      outputMat(ii) = mat1(ii) - mat2(ii); 
    }
    return outputMat;
  }

  std::complex<double> trace(const Matrix &mat)
  {
    std::complex<double> output = 0;
    int matSize = mat.getSize();
    for (int ii = 0; ii < matSize; ii++)
    {
      output += mat(ii,ii);
    }
    return output;
  }

  monsta::Matrix conjugateTranspose(const Matrix &mat)
  {
    monsta::Matrix outputMat = mat;
    int matSize = mat.getSize();

    for (int ii = 0; ii < matSize; ii++)
    {
      for (int jj = 0; jj < matSize; jj++)
      {
        outputMat(ii,jj) = conj(mat(jj,ii));
      }
    }

    return outputMat;
  }
}

#endif