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
#include "Slicer3LandmarkIO.h"


int main()
{
  bool testPassedStatus = true;

  // read in bogus landmarks file
  try
    {
    std::cout << "Attempting to read in non-exitsant landmarks file" <<std::endl;
    //hopefully there is no landmarks file with this path in your system.
    LandmarksMapType bogusLandmarks = ReadSlicer3toITKLmk("/this/file/does/not/exist.fcsv");

    // An exception should be thrown, and this code should never be reached
    std::cout << "!!Error: Succesfully read non-existant file ?!?!?! :(" <<std::endl;
    testPassedStatus &= false;
    }
  catch (itk::ExceptionObject &exception)
    {
    // We want to succesfully catch this exception
    std::cout << "Successfully caught exception for landmark file reading." << std::endl;
    std::cout << "Exception was:";
    exception.Print(std::cout);
    testPassedStatus &= true;
    }
  try
    {
    std::cout << "Attempting to read in non-exitsant landmark weight file" <<std::endl;
    //hopefully there is no landmarks file with this path in your system.
    LandmarkWeightMapType bogusLandmarksWeights = ReadLandmarkWeights("/this/file/does/not/exist/either.fcsv");

    // An exception should be thrown, and this code should never be reached
    std::cout << "!!Error: Succesfully read non-existant file ?!?!?! :(" <<std::endl;
    testPassedStatus &= false;
    }
  catch (itk::ExceptionObject &exception)
    {
    // We want to succesfully catch this exception
    std::cout << "Successfully caught exception for landmark weights file reading." << std::endl;
    std::cout << "Exception was:";
    exception.Print(std::cout);
    testPassedStatus &= true;
    }

  return testPassedStatus ? EXIT_SUCCESS: EXIT_FAILURE;
}