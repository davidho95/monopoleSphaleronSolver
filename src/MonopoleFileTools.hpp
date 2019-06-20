#ifndef MONOPOLEFILETOOLS_HPP
#define MONOPOLEFILETOOLS_HPP

#include "LATfield2.hpp"
#include <complex>
// #include "GeorgiGlashowSu2TheoryNoMonopoles.hpp"
#include "Su2Tools.hpp"
#include <cstdio>
#include <iostream>
#include <fstream>

namespace monsta
{
  void mergeFiles(std::string fileBaseName, int numFiles)
  {
    // if (parallel.size() == 1) { return; }
    parallel.barrier();
    ofstream fileStream;
    fileStream.open(fileBaseName + ".txt");
    for (int ii = 0; ii < numFiles; ii++)
    {
      ifstream fileToMerge;
      std::string fileToMergeName = fileBaseName + std::to_string(ii) + ".txt";
      fileToMerge.open(fileToMergeName);
      fileStream << fileToMerge.rdbuf();
      fileToMerge.close();
      remove((fileToMergeName).c_str());
    }
    fileStream.close();
  }

  std::string getFileName(std::string fileBaseName)
  {
    if (parallel.size() == 1)
    {
      return fileBaseName + ".txt";
    } else {
      return fileBaseName + std::to_string(parallel.rank()) + ".txt";
    }
  }

  void writeRawField(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName)
  {
    ofstream fileStream;
    std::string fileName = getFileName(fileBaseName);
    fileStream.open(fileName);

    LATfield2::Site site(field.lattice());

    for (site.first(); site.test(); site.next())
    {
      for (int cpt = 0; cpt < field.components(); cpt++)
      {
        fileStream << field(site, cpt) << std::endl;
      }
    }
    fileStream.close();
    parallel.barrier();

    // mergeFiles(fileBaseName);
  }

  void readRawField(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName)
  {
    ifstream fileToRead;
    std::string fileName = getFileName(fileBaseName);
    fileToRead.open(fileName);

    LATfield2::Site site(field.lattice());
    std::string line;

    for (site.first(); site.test(); site.next())
    {
      for (int cpt = 0; cpt < field.components(); cpt++)
      {
        std::getline(fileToRead, line);
        std::istringstream lineStream(line);
        std::complex<double> value;
        lineStream >> value;
        field(site, cpt) = value;
      }
    }

    fileToRead.close();
    parallel.barrier();
  }

  void readHiggsField(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName)
  {
    ifstream fileToRead;
    std::string fileName = getFileName(fileBaseName);
    fileToRead.open(fileName);

    LATfield2::Site site(field.lattice());
    std::string line;

    for (site.first(); site.test(); site.next())
    {
      std::getline(fileToRead, line);
      std::istringstream lineStream(line);
      std::complex<double> value;
      lineStream >> value;
      field(site, 3, 0, 0) = value;
    }

    fileToRead.close();
    parallel.barrier();
  }

  void writeCoords(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName)
  {
    ofstream fileStream;
    std::string fileName = getFileName(fileBaseName);
    fileStream.open(fileName);

    LATfield2::Site site(field.lattice());

    for (site.first(); site.test(); site.next())
    {
      fileStream << site.coord(0) << " " << site.coord(1) << " " << site.coord(2) << " " << std::endl;
    }
    fileStream.close();
    parallel.barrier();

    // mergeFiles(fileBaseName);
  }

  void writeHiggsField(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName)
  {
    ofstream fileStream;
    std::string fileName = getFileName(fileBaseName);
    fileStream.open(fileName);

    LATfield2::Site site(field.lattice());

    std::vector<double> su2Vec;
    for (site.first(); site.test(); site.next())
    {
      su2Vec = monsta::su2LieAlgToVec(monsta::Matrix(field, site, 3));
      fileStream << su2Vec[0] << " " << su2Vec[1] << " " << su2Vec[2] << std::endl;
    }
    fileStream.close();
    parallel.barrier();

    // mergeFiles(fileBaseName);
  }

  void writeHiggsFieldUnitary(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName)
  {
    ofstream fileStream;
    std::string fileName = getFileName(fileBaseName);
    fileStream.open(fileName);

    LATfield2::Site site(field.lattice());

    std::vector<double> su2Vec;
    for (site.first(); site.test(); site.next())
    {
      fileStream << real(field(site, 3, 0, 0)) << std::endl;
    }
    fileStream.close();
    parallel.barrier();

    // mergeFiles(fileBaseName);
  }


  void writeMagneticField(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowSu2Theory &theory)
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

    // mergeFiles(fileBaseName);
  }

  void writeUnitaryGaugeField(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName)
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
        fileStream << arg(field(site, ii, 0, 0)) << " ";
      }
      fileStream << endl;
    }
    fileStream.close();
    parallel.barrier();
  }

  void writeEnergyDensity(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::Theory &theory)
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

    // mergeFiles(fileBaseName);
  }
}

#endif