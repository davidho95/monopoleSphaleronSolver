// #include "LATfield2.hpp"
#include <cmath>
#include <complex>
#include <ctime>
#include <iostream>
#include <array>
#include "fieldMath.h"

// using namespace LATfield2;

// template <class FieldType>
// Field<FieldType> operator+(Field<FieldType> &a, Field<FieldType> &b) {
//   Lattice &lat = a.lattice();
//   Field<FieldType> out(lat);
  
//   Site x(lat);
//   for (x.first(); x.test(); x.next()) {
//     out(x) = a(x) + b(x);
//   }
//   out.updateHalo();
  
//   return out;
// }

// template <class FieldType>
// Field<FieldType> operator+=(Field<FieldType> &a, Field<FieldType> &b) {
//   Lattice &lat = a.lattice();
//   int numElem = lat.sitesGross();
//   for (int iElem = 0; iElem < numElem; iElem++) {
//     a.data_[iElem] += b.data_[iElem];
//   }
// }

// template <class FieldType>
// Field<FieldType> operator-(Field<FieldType> &a, Field<FieldType> &b) {
//   Lattice &lat = a.lattice();
//   Field<FieldType> out(lat);
//   int numElem = lat.sitesGross();
//   for (int iElem = 0; iElem < numElem; iElem++) {
//     out.data_[iElem] = a.data_[iElem] - b.data_[iElem];
//   }
  
//   return out;
// }

// template <class FieldType>
// Field<FieldType> operator-=(Field<FieldType> &a, Field<FieldType> &b) {
//   Lattice &lat = a.lattice();
//   int numElem = lat.sitesGross();
//   for (int iElem = 0; iElem < numElem; iElem++) {
//     a.data_[iElem] -= b.data_[iElem];
//   }
// }

// template <class FieldType>
// Field<FieldType> operator*(Field<FieldType> &a, Field<FieldType> &b) {
//   Lattice &lat = a.lattice();
//   Field<FieldType> out(lat);
//   int numElem = lat.sitesGross();
//   for (int iElem = 0; iElem < numElem; iElem++) {
//     out.data_[iElem] = a.data_[iElem] * b.data_[iElem];
//   }
  
//   return out;
// }

// template <class FieldType>
// Field<FieldType> operator/(Field<FieldType> &a, Field<FieldType> &b) {
//   Lattice &lat = a.lattice();
//   Field<FieldType> out(lat);
//   int numElem = lat.sitesGross();
//   for (int iElem = 0; iElem < numElem; iElem++) {
//     out.data_[iElem] = a.data_[iElem] / b.data_[iElem];
//   }
  
//   return out;
// }

// template <class FieldType, class FuncType>
// auto mapLocalFcn(Field<FieldType> &field, FuncType localFcn) {
//   Lattice &lat = field.lattice();
//   Field<FieldType> out(lat);
//   Site x(lat);
//   int numCpts = field.components();

//   for (x.first(); x.test(); x.next()) {
//     for (int iCpt = 0; iCpt < numCpts; iCpt++) {
//       out(x, iCpt) = localFcn(field, x);
//     }
//   }

//   out.updateHalo();

//   return out;
// }

// template <
//   class FieldType1,
//   class FieldType2,
//   class FuncType,
//   class outType = result_of_t<FuncType(FieldType1, FieldType2)>
//   >
// Field<outType> mapBinaryFcn(Field<FieldType1> &field1, Field<FieldType2> &field2, FuncType binaryFcn) {
//   // TODO: size check
//   Lattice &lat = field1.lattice();
//   Site x(lat);

//   Field<outType> out(lat);

//   for (x.first(); x.test(); x.next()) {
//     out(x) = binaryFcn(field1(x), field2(x));
//   }

//   out.updateHalo();

//   return out;
// } 

// template <class FieldType>
// Field<FieldType> diff(Field<FieldType> &field, int dim, bool isForward) {
//   Lattice &lat = field.lattice();
//   Site x(lat);
//   int numCpts = field.components();

//   Field<FieldType> out(lat);

//   if (isForward) {
//     for (x.first(); x.test(); x.next()) {
//       for (int iCpt = 0; iCpt < numCpts; iCpt++) {
//         out(x, iCpt) = field(x + dim, iCpt) - field(x, iCpt);
//       }
//     }
//   } else {
//     for (x.first(); x.test(); x.next()) {
//       for (int iCpt = 0; iCpt < numCpts; iCpt++) {
//         out(x, iCpt) = -(field(x - dim, iCpt) - field(x, iCpt));
//       }
//     }
//   }

//   out.updateHalo();

//   return out;
// }

// template <class FieldType>
// Field<FieldType> laplacian(Field<FieldType> &field) {
//   Lattice &lat = field.lattice();
//   Site x(lat);
//   int numCpts = field.components();
//   int dim = lat.dim();

//   Field<FieldType> out(lat);

//   for (x.first(); x.test(); x.next()) {
//     for (int iCpt = 0; iCpt < numCpts; iCpt++) {
//       for (int iDim = 0; iDim < dim; iDim++) {
//         out(x, iCpt) += field(x + iDim, iCpt) + field(x - iDim, iCpt) - 2*field(x, iCpt);
//       }
//     }
//   }

//   out.updateHalo();

//   return out;
// }

template <class FieldType>
LATfield2::Field<FieldType> curl(LATfield2::Field<FieldType> &field, bool isForward) {
  LATfield2::Lattice &lat = field.lattice();
  LATfield2::Site x(lat);
  int numCpts = field.components();
  int dim = lat.dim();
  if (numCpts != 3 && dim != 3) {
    std::cerr << "Curl is only defined for 3D vector fields with 3 components" << endl;
  }

  LATfield2::Field<FieldType> out(lat, numCpts);

  if (isForward) {
    for (x.first(); x.test(); x.next()) {
      out(x, 0) = field(x + 1, 2) - field(x, 2) - field(x + 2, 1) + field(x,1);
      out(x, 1) = field(x + 2, 0) - field(x, 0) - field(x + 0, 2) + field(x,2);
      out(x, 2) = field(x + 0, 1) - field(x, 1) - field(x + 1, 0) + field(x,0);
    }
  } else {
    for (x.first(); x.test(); x.next()) {
      out(x, 0) = -(field(x - 1, 2) - field(x, 2) - field(x - 2, 1) + field(x,1));
      out(x, 1) = -(field(x - 2, 0) - field(x, 0) - field(x - 0, 2) + field(x,2));
      out(x, 2) = -(field(x - 0, 1) - field(x, 1) - field(x - 1, 0) + field(x,0));
    }
  }
  out.updateHalo();

  return out;
}

// int main() {
//   parallel.initialize(1,1);

//   srand(time(NULL) + parallel.rank());

//   int dim = 3;
//   int latSize[dim] = {2,2,2};
//   int halo = 1;
//   Lattice lat(dim,latSize,halo);

//   // double vev = 1;

//   Field<double> phi(lat,3);

//   Site x(lat);

//   int count = 0;
//   for (x.first(); x.test(); x.next()) {
//     for (int iCpt = 0; iCpt < 3; iCpt++) {
//       phi(x,iCpt) = count;
//     }
//     count++;
//   }

//   phi.updateHalo();
//   Field<double> curlPhi = curl(phi, true);
//   Field<double> curlCurlPhi = curl(curlPhi, false);

//   for (x.first(); x.test(); x.next()) {
//     for (int iCpt = 0; iCpt < 3; iCpt++) {
//       cout << phi(x,iCpt) << " " << curlCurlPhi(x,iCpt) << endl;
//     }
//   }

  // int numTrials = 100;
  // clock_t tStart = clock();
  // for (int ii = 0; ii < numTrials; ++ii) {
  //   laplacian(phi);
  // }
  // printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

  // tStart = clock();
  // for (int ii = 0; ii < numTrials; ++ii) {
  //   laplacian2(phi);
  // }
  // printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

  // Field<complex<double>> eDen = mapLocalFcn(phi, quarticPotentialFcn);
  // complex<double> E = fieldSum(eDen);
  // cout << E << endl;
  // parallel.sum(E);
  // cout << E << endl;

//   return 0;
// }