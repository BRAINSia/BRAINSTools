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
 * University of Iowa Health Care 2013
 */

/*
 This program program uses the input transform to propagate
 reference landmarks to the target landmark file.
*/

#include "itkPoint.h"
#include "itkTransformFileReader.h"
#include "itkCompositeTransform.h"
#include "Slicer3LandmarkIO.h"

#include "GenericTransformImage.h"
#include "BRAINSConstellationLandmarksTransformCLP.h"
#include <BRAINSCommonLib.h>


int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  if (inputLandmarksFile.empty() || inputTransformFile.empty() || outputLandmarksFile.empty())
  {
    std::cerr << "Input and output file names should be given by commandline. " << std::endl;
    std::cerr << "Usage:\n"
              << "~/BRAINSConstellationLandmarksTransform\n"
              << "--inputLandmarksFile (-i)\n"
              << "--inputTransformFile (-t)\n"
              << "--outputLandmarksFile (-o)" << std::endl;
    return EXIT_FAILURE;
  }

  LandmarksMapType origLandmarks = ReadSlicer3toITKLmk(inputLandmarksFile);
  LandmarksMapType transformedLandmarks;

  using ReaderType = itk::TransformFileReader;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputTransformFile);
  reader->Update();

  ReaderType::TransformListType * transformList = reader->GetModifiableTransformList();
  if(transformList->size() != 1)
  {
    std::cerr << "ERROR: Only one transform allowed." << std::endl;
    return EXIT_FAILURE;
  }

  auto transformIterator = transformList->begin();
  itk::Transform<double,3,3>::Pointer generic_transform = dynamic_cast<itk::Transform<double,3,3> *>((*transformIterator).GetPointer());
  if( generic_transform.IsNull())
  {
    std::cerr << "The input transform not castable to  Transform<double,3,3> baseclass type." << std::endl;
    return EXIT_FAILURE;
  }

  LandmarksMapType::const_iterator it = origLandmarks.begin();
  for (; it != origLandmarks.end(); it++)
  {
    transformedLandmarks[it->first] = generic_transform->TransformPoint(it->second);
  }

  WriteITKtoSlicer3Lmk(outputLandmarksFile, transformedLandmarks);
  std::cout << "The transformed landmarks file is written." << std::endl;

  return EXIT_SUCCESS;
}
