/*
 * Author: Ali Ghayoor
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2011
 */

// This program adds a new landmark as the midpoint between AC and PC points to the output landmark fcsv file

// This program gets a landmark1.fcsv file as its input and returns landmark2.fcsv file which contains
// the new added landmark that is the midpoint between AC and PC points.
// Therefore, the first argument is the inputlandmark file and the second argument is the name of the outputlandmark
// file
// For instans:
//               ./midACPCpoint {inputlandmarkFile}.fcsv {NameOfoutputlandmarkFile}.fcsv

#include "itkImage.h"
#include "math.h"
#include "Slicer3LandmarkIO.h"
#include "insertMidACPCpointCLP.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  std::map<std::string, PointType> lmksMap;
  PointType                        ACpoint;
  PointType                        PCpoint;
  PointType                        midACPCpoint;

  lmksMap = ReadSlicer3toITKLmk(inputLandmarkFile);

  ACpoint = lmksMap["AC"];
  PCpoint = lmksMap["PC"];

  midACPCpoint[0] = (PCpoint[0] - ACpoint[0]) / 2;
  midACPCpoint[1] = (PCpoint[1] - ACpoint[1]) / 2;
  midACPCpoint[2] = (PCpoint[2] - ACpoint[2]) / 2;

//	std::cout << "PC :" << PCpoint[0] << ", " << PCpoint[1] << ", " << PCpoint[2] << ")" << std::endl;
//	std::cout << "midpoint :" << midACPCpoint[0] << ", " << midACPCpoint[1] << ", " << midACPCpoint[2] << ")" <<
// std::endl;

  lmksMap["midACPC"] = midACPCpoint;

  WriteITKtoSlicer3Lmk(outputLandmarkFile, lmksMap);

  return EXIT_SUCCESS;
}
