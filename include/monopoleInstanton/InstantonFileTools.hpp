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
#include "GeorgiGlashowSu2Theory4d.hpp"
#include "GeorgiGlashowRadialTheory.hpp"

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

  void writeHiggsField(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowSu2Theory4d &theory)
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

  void writeMagneticField(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowSu2Theory4d &theory)
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

  void writeMagneticField(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowRadialTheory &theory)
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

  void writeEnergyDensity(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowRadialTheory &theory)
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

  void writeEnergyDensity(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowSu2Theory4d &theory)
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

  void writeGradients(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowSu2Theory4d &theory)
  {
    ofstream fileStream;
    std::string fileName = getFileName(fileBaseName);
    fileStream.open(fileName);

    LATfield2::Site site(field.lattice());

    for (site.first(); site.test(); site.next())
    {
      double gradVal = 0;
      for (int ii = 0; ii < 4; ii++)
      {
        Matrix grad = theory.getLocalGradient(field, site, ii);
        gradVal += real(trace(grad*conjugateTranspose(grad)));
      }
      fileStream << gradVal << endl;
    }
    fileStream.close();
    parallel.barrier();
  }

  void writeGradients(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowSu2Theory &theory)
  {
    ofstream fileStream;
    std::string fileName = getFileName(fileBaseName);
    fileStream.open(fileName);

    LATfield2::Site site(field.lattice());

    for (site.first(); site.test(); site.next())
    {
      double gradVal = 0;
      for (int ii = 0; ii < 4; ii++)
      {
        Matrix grad = theory.getLocalGradient(field, site, ii);
        gradVal += real(trace(grad*conjugateTranspose(grad)));
      }
      fileStream << gradVal << endl;
    }
    fileStream.close();
    parallel.barrier();
  }

  void writeGradients(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowRadialTheory &theory)
  {
    ofstream fileStream;
    std::string fileName = getFileName(fileBaseName);
    fileStream.open(fileName);

    LATfield2::Site site(field.lattice());

    for (site.first(); site.test(); site.next())
    {
      double gradVal = 0;
      for (int ii = 0; ii < 4; ii++)
      {
        Matrix grad = theory.getLocalGradient(field, site, ii);
        gradVal += real(trace(grad*conjugateTranspose(grad)));
      }
      fileStream << gradVal << endl;
    }
    fileStream.close();
    parallel.barrier();
  }
}

#endif