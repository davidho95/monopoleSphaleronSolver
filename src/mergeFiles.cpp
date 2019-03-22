#include "LATfield2.hpp"
#include <complex>
#include "../src/GeorgiGlashowSu2TheoryUnitary.hpp"
#include "../src/Matrix.hpp"
#include "../src/TheoryChecker.hpp"
#include "../src/GradDescentSolverBBStep.hpp"
#include "../src/Su2Tools.hpp"
#include "../src/MonopoleFileTools.hpp"
#include "../src/MonopoleFieldTools.hpp"
#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{

  std::string outputPath;
  int numFiles;

  for (int i=1 ; i < argc ; i++ ){
    if ( argv[i][0] != '-' )
      continue;
    switch(argv[i][1]) {
      case 'p':
        outputPath = argv[++i];
        break;
      case 'n':
        numFiles = atoi(argv[++i]);
        break;
    }
  }

  monsta::mergeFiles(outputPath + "/coords", numFiles);
  monsta::mergeFiles(outputPath + "/higgsData", numFiles);
  monsta::mergeFiles(outputPath + "/magneticFieldData", numFiles);
  monsta::mergeFiles(outputPath + "/gaugeData", numFiles);
  monsta::mergeFiles(outputPath + "/energyData", numFiles);
}