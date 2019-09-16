/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
/*
 * Author: Ali Ghayoor
 * at Psychiatry Imaging Lab,
 * University of Iowa Health Care 2011
 */

// This program adds a new landmark as the midpoint between AC and PC points to the output landmark fcsv file

// This program gets a landmark1.fcsv file as its input and returns landmark2.fcsv file which contains
// the new added landmark "midACPC" which is the midpoint between AC and PC points.
//
// For use:
//             .../midACPCpoint --inputLandmarkFile {NameinputlandmarkFile}.fcsv --outputLandmarkFile
// {NameOfoutputlandmarkFile}.fcsv

#include "itkImage.h"
#include <cmath>
#include "Slicer3LandmarkIO.h"
#include "insertMidACPCpointCLP.h"
#include <BRAINSCommonLib.h>

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();


  std::map<std::string, LandmarkPointType> lmksMap = ReadSlicer3toITKLmk(inputLandmarkFile);

  const LandmarkPointType ACpoint = lmksMap["AC"];
  const LandmarkPointType PCpoint = lmksMap["PC"];

  const LandmarkPointType midACPCpoint = (PCpoint - ACpoint) * 0.5;

  //  std::cout << "PC :" << PCpoint[0] << ", " << PCpoint[1] << ", " << PCpoint[2] << ")" << std::endl;
  //  std::cout << "midpoint :" << midACPCpoint[0] << ", " << midACPCpoint[1] << ", " << midACPCpoint[2] << ")" <<
  // std::endl;

  lmksMap["midACPC"] = midACPCpoint;

  WriteITKtoSlicer3Lmk(outputLandmarkFile, lmksMap);

  return EXIT_SUCCESS;
}
