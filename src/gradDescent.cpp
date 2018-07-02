#include <cstdio>
// #include "LATfield2.hpp"
#include <cmath>
#include <ctime>
#include "fieldMath.h"

template <class FieldType>
void gradDesc(LATfield2::Field<FieldType> &field, auto iterate, auto energyFcn, double tol, int maxIter) {
  LATfield2::Lattice &lat = field.lattice();
  LATfield2::Site x(lat);

  double bigNumber = 1e6;
  double maxStepSize = 1;
  double stepSize = 0.01;

  double energy = energyFcn(field, x);
  double energyOld;
  double energyChange = bigNumber;
  double maxGrad = bigNumber;

  int numIters = 0;

  while (maxGrad > tol && numIters < maxIter) {
    numIters++;
    energyOld = energy;
    maxGrad = iterate(field, x, stepSize);
    energy = energyFcn(field, x);
    if (energy > energyOld) {
      stepSize /= 2;
      maxStepSize = stepSize;
      continue;
    }
    energyChange = energy - energyOld;
    if (stepSize < maxStepSize) {
      stepSize *= 2;
    }
    x.first();
    // cout << maxGrad << endl;

  }

  cout << "Gradient descent finished in " << numIters << " iterations." << endl;
}

template <class FieldType>
double iterFcn(LATfield2::Field<FieldType> &field, LATfield2::Site &x, double stepSize) {
  int numCpts = field.components();
  double absFieldSqr;
  double grad;
  double maxGrad = 0;
  for (x.first(); x.test(); x.next()) {
    absFieldSqr = 0;
    for (int iCpt = 0; iCpt < numCpts; iCpt++) {
      absFieldSqr += pow(field(x,iCpt), 2);
    }
    for (int iCpt = 0; iCpt < field.components(); iCpt++) {
      grad = 4*field(x,iCpt)*(absFieldSqr - 1);
      field(x,iCpt) -= stepSize*grad;
      if (abs(grad) > maxGrad) {
        maxGrad = abs(grad);
      }
    }
  }
  field.updateHalo();
  return maxGrad;
}

template<class FieldType>
double energy(LATfield2::Field<FieldType> &field, LATfield2::Site &x) {
  int numCpts = field.components();
  double absFieldSqr;
  double E = 0;
  for (x.first(); x.test(); x.next()) {
    absFieldSqr = 0;
    for (int iCpt = 0; iCpt < numCpts; iCpt++) {
      absFieldSqr += pow(field(x,iCpt), 2);
    }
    E += pow(absFieldSqr - 1, 2);
  }
  return E;
}

double magneticEnergy(LATfield2::Field<double> &gaugeField, LATfield2::Site &x) {
  double E = 0;
  for (x.first(); x.test(); x.next()) {
    E += pow(gaugeField(x+1, 2) - gaugeField(x, 2) - gaugeField(x+2, 1) + gaugeField(x,1), 2);
    E += pow(gaugeField(x+2, 0) - gaugeField(x, 0) - gaugeField(x+0, 2) + gaugeField(x,2), 2);
    E += pow(gaugeField(x+0, 1) - gaugeField(x, 1) - gaugeField(x+1, 0) + gaugeField(x,0), 2);
  }
  return 0.5*E;
}

double magneticGradIterCopy(LATfield2::Field<double> &A, LATfield2::Site &x, double stepSize) {
  int numCpts = 3;
  LATfield2::Field<double> grad(A.lattice(), 3);
  double maxGrad = 0;
  // for (x.first(); x.test(); x.next()) {
  //   grad(x,0) = -A(x + 0 - 1, 1) + A(x - 1, 1) + A(x + 0, 1) - A(x, 1)
  //             - A(x + 0 - 2, 2) + A(x - 2, 2) + A(x + 0, 2) - A(x, 2)
  //             + 4*A(x, 0) - A(x + 1, 0) - A(x - 1, 0) - A(x + 2, 0) - A(x - 2, 0);
  //   grad(x,1) = -A(x + 1 - 0, 0) + A(x - 0, 0) + A(x + 1, 0) - A(x, 0)
  //             - A(x + 1 - 2, 2) + A(x - 2, 2) + A(x + 1, 2) - A(x, 2)
  //             + 4*A(x, 1) - A(x + 0, 1) - A(x - 0, 1) - A(x + 2, 1) - A(x - 2, 1);
  //   grad(x,2) = -A(x + 2 - 1, 1) + A(x - 1, 1) + A(x + 2, 1) - A(x, 1)
  //             - A(x + 2 - 0, 0) + A(x - 0, 0) + A(x + 2, 0) - A(x, 0)
  //             + 4*A(x, 2) - A(x + 1, 2) - A(x - 1, 2) - A(x + 0, 2) - A(x - 0, 2);
  // }
  for (x.first(); x.test(); x.next()) {
    for (int iCpt = 0; iCpt < numCpts; iCpt++) {
      grad(x, iCpt) = 0;
      for (int jj = 0; jj < numCpts; jj++) {
        if (iCpt != jj) {
          grad(x, iCpt) -= A(x + iCpt - jj, jj) - A(x - jj, jj) - A(x + iCpt, jj) + A(x, jj);
          grad(x, iCpt) += 2*A(x, iCpt) - A(x + jj, iCpt) - A(x - jj, iCpt);
        }
      }
    }
  }
  double gradCpt;
  for (x.first(); x.test(); x.next()){
    for (int iCpt = 0; iCpt < numCpts; iCpt++) {
      gradCpt = grad(x, iCpt);
      cout << gradCpt << " ";
      A(x, iCpt) -= stepSize*gradCpt;
      if (abs(gradCpt) > maxGrad) {
        maxGrad = abs(gradCpt);
      }
    }
  }
  cout << endl;
  A.updateHalo();
  return maxGrad;
}

double magneticGradIter2(LATfield2::Field<double> &A, LATfield2::Site &x, double stepSize) {
  int numCpts = 3;
  LATfield2::Field<double> curlA = curl(A, true);
  LATfield2::Field<double> grad = curl(curlA, false);
  double maxGrad = 0;
  double gradCpt;
  for (x.first(); x.test(); x.next()){
    for (int iCpt = 0; iCpt < numCpts; iCpt++) {
      gradCpt = grad(x, iCpt);
      cout << gradCpt << " ";
      A(x, iCpt) -= stepSize*gradCpt;
      if (abs(gradCpt) > maxGrad) {
        maxGrad = abs(gradCpt);
      }
    }
  }
  cout << endl;
  A.updateHalo();
  return maxGrad;
}


int main() {
  parallel.initialize(1,1);

  // srand(time(NULL) + parallel.rank());
  srand(1);

  int dim = 3;
  int sz = 2;
  int latSize[dim] = {sz, sz, sz};
  int haloSize = 1;
  int numCpts = 3;

  double tol = 1e-6;
  int maxIter = 1;

  LATfield2::Lattice lat(dim, latSize, haloSize);
  LATfield2::Field<double> phi1(lat, numCpts);
  LATfield2::Field<double> phi2(lat, numCpts);
  LATfield2::Site x(lat);

  int count = 0;
  for (x.first(); x.test(); x.next()) {
    count++;
    for (int iCpt = 0; iCpt < numCpts; ++iCpt) {
      // phi1(x,iCpt) = (std::rand() % 10) / 10.0;
      // phi2(x, iCpt) = (std::rand() % 10) / 10.0;
      phi1(x,iCpt) = sin(count);
      phi2(x,iCpt) = sin(count);
      cout << phi1(x, iCpt) << " ";
    }
  }
  cout << endl;

  // srand(1);
  // for (x.first(); x.test(); x.next()) {
  //   for (int iCpt = 0; iCpt < numCpts; ++iCpt) {
  //     phi2(x,iCpt) = (std::rand() % 10) / 10.0;
  //     cout << phi1(x, iCpt) << endl;
  //     cout << phi2(x, iCpt) << endl;
  //   }
  // }

  phi1.updateHalo();
  phi2.updateHalo();


  std::cout << magneticEnergy(phi1, x) << endl;
  std::cout << magneticEnergy(phi2, x) << std::endl;

  auto gradLambda = [] (LATfield2::Field<double> &field, LATfield2::Site &x, double stepSize)
    { return magneticGradIter2(field, x, stepSize); };
  auto gradLambda2 = [] (LATfield2::Field<double> &field, LATfield2::Site &x, double stepSize)
    { return magneticGradIterCopy(field, x, stepSize); };

  auto energyLambda = [] (LATfield2::Field<double> &field, LATfield2::Site &x)
    { return magneticEnergy(field, x); };

  // cout << gradLambda(phi1, x, 1) << endl;
  // cout << gradLambda2(phi2, x, 1) << endl;

  // for (x.first(); x.test(); x.next()) {
  //   for (int iCpt = 0; iCpt < numCpts; ++iCpt) {
  //     cout << phi1(x, iCpt) << " " << phi2(x, iCpt) << endl;
  //   }
  // }


  clock_t startTime = clock();
  gradDesc(phi1, gradLambda, energyLambda, tol, maxIter);
  std::cout << "Time taken: " << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << std::endl;

  startTime = clock();
  gradDesc(phi2, gradLambda2, energyLambda, tol, maxIter);
  std::cout << "Time taken: " << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << std::endl;

  // double meanAbsSqr = 0;
  // for (x.first(); x.test(); x.next()) {
  //   for (int iCpt = 0; iCpt < numCpts; ++iCpt) {
  //     meanAbsSqr += pow(phi(x,iCpt), 2) / pow(sz,3);
  //   }
  // }

  // std::cout << meanAbsSqr << std::endl;

  std::cout << magneticEnergy(phi1, x) << std::endl;
  std::cout << magneticEnergy(phi2, x) << std::endl;
}