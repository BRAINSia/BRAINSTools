//
// Created by Johnson, Hans J on 11/19/16.
//
// This file is used to compare BVec and BVal Files
// to test if conversion to/from DICOM/FSL occurs properly
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>


int main(int argc, char * argv[])
{
  if( argc != 4 )
  {
    std::cout << "USAGE: " << argv[0] << " <BFile> <BFile> <linesToCompare> <tolerance>" << std::endl;
  }
  const std::string refBValFileName(argv[1]);
  const std::string newBValFileName(argv[2]);
  const int linesToCompare = atoi(argv[3]);
  const double tolerance = atof(argv[4]);

  std::ifstream newFile(newBValFileName);
  std::string newStr;
  std::ifstream refFile(refBValFileName);
  std::string refStr;
  std::vector<double> refVector;
  refVector.reserve(200);
  std::vector<double> newVector;
  newVector.reserve(200);
  unsigned int linesReadSuccessfully = 0;

  double refValue;
  double newValue;
  while (std::getline(refFile, refStr) && std::getline(newFile, newStr))
  {
    std::stringstream refSsin(refStr);
    std::stringstream newSsin(newStr);
    while (refSsin.good() && newSsin.good()){
      refSsin >> refValue;
      newSsin >> newValue;
      refVector.push_back(refValue);
      newVector.push_back(newValue);
    }
    linesReadSuccessfully++;
    // Process str
  }

  if( linesReadSuccessfully != linesToCompare )
  {
    std::cout << "ERROR: Incorrect number of lines read. "
              << linesReadSuccessfully << " != " << linesToCompare << std::endl;
    return EXIT_FAILURE;
  }


  return EXIT_SUCCESS;
}
