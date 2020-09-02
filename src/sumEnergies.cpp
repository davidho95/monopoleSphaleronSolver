#include "LATfield2.hpp"
#include <complex>
#include "../include/monopoleSphaleron/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../include/Matrix.hpp"
#include "../include/TheoryChecker.hpp"
#include "../include/GradDescentSolverBBStep.hpp"
#include "../include/Su2Tools.hpp"
#include "../include/MonopoleFileTools.hpp"
#include "../include/MonopoleFieldTools.hpp"
#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{

  std::string inputPath;
  int numFiles;

  for (int i=1 ; i < argc ; i++ ){
    if ( argv[i][0] != '-' )
      continue;
    switch(argv[i][1]) {
      case 'p':
        inputPath = argv[++i];
        break;
      case 'n':
        numFiles = atoi(argv[++i]);
        break;
    }
  }

  double totalEnergy = 0;

  for (int ii = 0; ii < numFiles; ii++)
  {
    ifstream fileToRead;
    std::string fileName = inputPath + "/energyData" + std::to_string(ii) + ".txt";

    // std::cout << fileName << endl;
    fileToRead.open(fileName);

    std::string line;

    while (getline(fileToRead, line))
    {
      std::istringstream lineStream(line);
      double value;
      lineStream >> value;
      totalEnergy += value;
    }
  }

  std::cout << totalEnergy << endl;
}