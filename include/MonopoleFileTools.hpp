#ifndef MONOPOLEFILETOOLS_HPP
#define MONOPOLEFILETOOLS_HPP

#include "LATfield2.hpp"
#include <complex>
// #include "GeorgiGlashowSu2TheoryUnitary.hpp"
#include "electroweakSphaleron/ElectroweakTheory.hpp"
#include "electroweakSphaleron/ElectroweakTheoryNonUnitary.hpp"
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
      std::string fileToMergeName = fileBaseName + std::to_string(ii) + ".txt";
      ifstream fileToMerge(fileToMergeName.c_str());
      if (!((bool)fileToMerge)) { COUT << "jonnyMorris" << endl; break; }
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

  void readFileWithCoords(LATfield2::Field< std::complex<double> > &field, std::string coordFile, std::string rawDataFile)
  {
    ifstream coordStream;
    ifstream rawDataStream;
    coordStream.open(coordFile);
    rawDataStream.open(rawDataFile);

    LATfield2::Site site(field.lattice());
    std::string coordLine;
    std::string rawDataLine;

    while(std::getline(coordStream, coordLine))
    {
      int coord;
      std::vector<int> coords;
      std::istringstream coordLineStream(coordLine);
      while (coordLineStream >> coord)
      {
        coords.push_back(coord);
      }
      bool isLocal = site.setCoord(coords[0], coords[1], coords[2]);
      if (!isLocal) { continue; }
      for (int cpt = 0; cpt < field.components(); cpt++)
      {
        std::getline(rawDataStream, rawDataLine);
        std::istringstream rawDataLineStream(rawDataLine);
        std::complex<double> value;
        rawDataLineStream >> value;
        field(site, cpt) = value;
      }
    }


    coordStream.close();
    rawDataStream.close();
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
      for (int ii = 0; ii < field.lattice().dim(); ii++)
      {
        fileStream << site.coord(ii) << " ";
      }
      fileStream << std::endl;
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

  void writeHiggsField(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::ElectroweakTheory &theory)
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

    // mergeFiles(fileBaseName);
  }

  void writeHiggsFieldMagnitude(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName)
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

  void writeHiggsFieldMagnitude(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowSu2Theory &theory)
  {
    ofstream fileStream;
    std::string fileName = getFileName(fileBaseName);
    fileStream.open(fileName);

    LATfield2::Site site(field.lattice());

    std::vector<double> su2Vec;
    for (site.first(); site.test(); site.next())
    {
      fileStream << theory.getHiggsMagnitude(field, site) << std::endl;
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

  void writeMagneticField(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::ElectroweakTheory &theory)
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
    for (site.first(); site.test(); site.next())
    {
      for (int ii = 0; ii < 3; ii++)
      {
        std::vector<double> su2Vec = su2ToVec(Matrix(field, site, ii));
        fileStream << su2Vec[2] << " ";
      }
      fileStream << endl;
    }
    fileStream.close();
    parallel.barrier();
  }

  void writeUnitaryWField(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName)
  {
    ofstream fileStream;
    std::string fileName = getFileName(fileBaseName);
    fileStream.open(fileName);

    LATfield2::Site site(field.lattice());
    for (site.first(); site.test(); site.next())
    {
      for (int ii = 0; ii < 3; ii++)
      {
        std::vector<double> su2Vec = su2ToVec(Matrix(field, site, ii));
        fileStream << su2Vec[0] << " " << su2Vec[1] << " ";
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

  void writeSymmetricEnergyDensity(LATfield2::Field< std::complex<double> > &field, std::string fileBaseName, monsta::GeorgiGlashowSu2Theory &theory)
  {
    ofstream fileStream;
    std::string fileName = getFileName(fileBaseName);
    fileStream.open(fileName);

    LATfield2::Site site(field.lattice());

    std::vector<double> su2Vec;
    for (site.first(); site.test(); site.next())
    {
      fileStream << theory.getSymmetricEnergyDensity(field, site) << std::endl;
    }
    fileStream.close();
    parallel.barrier();

    // mergeFiles(fileBaseName);
  }
}

#endif