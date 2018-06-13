#include "LATfield2.hpp"
#include <cmath>
#include <complex>
#include <ctime>
#include <iostream>

using namespace LATfield2;

template <class FieldType, class FuncType>
auto mapLocalFcn(Field<FieldType> &field, FuncType localFcn) {
  Lattice &lat = field.lattice();
  Field<FieldType> out(lat);
  Site x(lat);

  for (x.first(); x.test(); x.next()) {
    out(x) = localFcn(field(x));
  }

  out.updateHalo();

  return out;
}

template <
  class FieldType1,
  class FieldType2,
  class FuncType,
  class outType = result_of_t<FuncType(FieldType1, FieldType2)>
  >
Field<outType> mapBinaryFcn(Field<FieldType1> &field1, Field<FieldType2> &field2, FuncType binaryFcn) {
  // TODO: size check
  Lattice &lat = field1.lattice();
  Site x(lat);

  Field<outType> out(lat);

  for (x.first(); x.test(); x.next()) {
    out(x) = binaryFcn(field1(x), field2(x));
  }

  out.updateHalo();

  return out;
} 

template <class FieldType, class FuncType>
void mapInplace(Field<FieldType> &field, FuncType localFcn) {
  Lattice &lat = field.lattice();
  Site x(lat);

  for (x.first(); x.test(); x.next()) {
    field(x) = localFcn(field(x));
  }

  field.updateHalo();
}

template <class FieldType>
FieldType fieldSum(Field<FieldType> &field) {
  Site x(field.lattice());
  FieldType total = field(x) - field(x); // Sets total to additive identity of FieldType

  for (x.first(); x.test(); x.next()) {
    total += field(x);
  }

  return total;
}

complex<double> quarticPotential(const complex<double> z, const double vev) {
  complex<double> out = pow(vev,2) - pow(abs(z), 2);
  out = pow(out, 2);

  return out;
}

complex<double> quarticPotentialDeriv(const complex<double> z, const double vev) {
  complex<double> out = pow(vev,2) - pow(abs(z), 2);
  out = z*out;

  return out;
}

int main() {
  parallel.initialize(1,1);

  srand(time(NULL) + parallel.rank());

  int dim = 2;
  int latSize[dim] = {2,2};
  int halo = 1;
  Lattice lat(dim,latSize,halo);

  // double vev = 1;

  Field<complex<double>> phi1(lat);
  Field<complex<double>> phi2(lat);

  Site x(lat);

  for (x.first(); x.test(); x.next()) {
    phi1(x) = ((rand() % 10 + 1) + (rand() % 10)*1i) / 10.0;
    phi2(x) = ((rand() % 10 + 1) + (rand() % 10)*1i) / 10.0;
  }

  // auto quarticPotentialFcn = [vev] (complex<double> z) { return quarticPotential(z, vev); };
  // auto quarticPotentialDerivFcn = [vev] (complex<double> z) { return quarticPotentialDeriv(z, vev); };

  // double step = 0.1;
  // int numSteps = 100;

  // auto gradIter = [vev, step] (complex<double> z) { return z + step*quarticPotentialDeriv(z, vev); };

  // clock_t tStart = clock();
  // for (int ii = 0; ii < numSteps; ++ii) {
  //   mapInplace(phi, gradIter);
  // }
  // printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

  // Field<complex<double>> eDen = mapLocalFcn(phi, quarticPotentialFcn);
  // complex<double> E = fieldSum(eDen);
  // cout << E << endl;
  // parallel.sum(E);
  // cout << E << endl;

  auto mult = [] (complex<double> a, complex<double> b) { return a*b; };

  Field<complex<double>> phi3 = mapBinaryFcn(phi1, phi2, mult);

  for (x.first(); x.test(); x.next()) {
    cout << phi1(x) << endl;
    cout << phi2(x) << endl;
    cout << phi3(x) << endl;
  }

  return 0;
}