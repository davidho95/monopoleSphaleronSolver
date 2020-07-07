#ifndef INTSTANTONFILETOOLS_HPP
#define INTSTANTONFILETOOLS_HPP

#include "LATfield2.hpp"
#include <complex>
#include "../Su2Tools.hpp"
#include <cstdio>
#include <iostream>
#include <fstream>
#include "GeorgiGlashowSu2Theory2d.hpp"
#include "GeorgiGlashowSu2TheoryAxisymmetric2.hpp"

namespace monsta
{
  void writeHiggsField(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowSu2TheoryAxisymmetric &theory)
  {
    ofstream fileStream;
    std::string fileName = getFileName(fileBaseName);
    fileStream.open(fileName);

    LATfield2::Site site(field.lattice());

    for (site.first(); site.test(); site.next())
    {
      fileStream << theory.getHiggsMagnitude(field, site) << endl;
    }
    fileStream.close();
    parallel.barrier();
  }

  void writeHiggsField(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowSu2TheoryAxisymmetric2d &theory)
  {
    ofstream fileStream;
    std::string fileName = getFileName(fileBaseName);
    fileStream.open(fileName);

    LATfield2::Site site(field.lattice());

    for (site.first(); site.test(); site.next())
    {
      fileStream << theory.getHiggsMagnitude(field, site) << endl;
    }
    fileStream.close();
    parallel.barrier();
  }

  void writeMagneticField(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowSu2TheoryAxisymmetric &theory)
  {
    ofstream fileStream;
    std::string fileName = getFileName(fileBaseName);
    fileStream.open(fileName);

    LATfield2::Site site(field.lattice());

    std::vector<double> su2Vec;
    for (site.first(); site.test(); site.next())
    {
      for (int ii = 0; ii < 3; ii++)
      {
        fileStream << theory.getMagneticField(field, site, ii) << " ";
      }
      fileStream << endl;
    }
    fileStream.close();
    parallel.barrier();
  }

  void writeMagneticField(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowSu2TheoryAxisymmetric2d &theory)
  {
    ofstream fileStream;
    std::string fileName = getFileName(fileBaseName);
    fileStream.open(fileName);

    LATfield2::Site site(field.lattice());

    std::vector<double> su2Vec;
    for (site.first(); site.test(); site.next())
    {
      for (int ii = 0; ii < 3; ii++)
      {
        fileStream << theory.getMagneticField(field, site, ii) << " ";
      }
      fileStream << endl;
    }
    fileStream.close();
    parallel.barrier();
  }

  // void writeEnergyDensity(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowSu2TheoryAxisymmetric &theory)
  // {
  //   ofstream fileStream;
  //   std::string fileName = getFileName(fileBaseName);
  //   fileStream.open(fileName);

  //   LATfield2::Site site(field.lattice());

  //   std::vector<double> su2Vec;
  //   for (site.first(); site.test(); site.next())
  //   {
  //     fileStream << theory.getNumThetaPoints() / (2*pi) * theory.getLocalEnergyDensity(field, site) / (abs(theory.getRFromSite(site))) << std::endl;
  //   }
  //   fileStream.close();
  //   parallel.barrier();
  // }

  void writeEnergyDensity(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowSu2TheoryAxisymmetric &theory)
  {
    ofstream fileStream;
    std::string fileName = getFileName(fileBaseName);
    fileStream.open(fileName);

    LATfield2::Site site(field.lattice());

    std::vector<double> su2Vec;
    for (site.first(); site.test(); site.next())
    {
      fileStream << theory.getLocalEnergyDensity(field, site) << std::endl;
    }
    fileStream.close();
    parallel.barrier();
  }

  void writeEnergyDensity(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowSu2TheoryAxisymmetric2d &theory)
  {
    ofstream fileStream;
    std::string fileName = getFileName(fileBaseName);
    fileStream.open(fileName);

    LATfield2::Site site(field.lattice());

    std::vector<double> su2Vec;
    for (site.first(); site.test(); site.next())
    {
      fileStream << theory.getLocalEnergyDensity(field, site) << std::endl;
    }
    fileStream.close();
    parallel.barrier();
  }
}

#endif