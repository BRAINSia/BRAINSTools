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
/*=========================================================================

 Program:   GTRACT (Guided Tensor Restore Anatomical Connectivity Tractography)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2006/03/29 14:53:40 $
 Version:   $Revision: 1.9 $

   Copyright (c) University of Iowa Department of Radiology. All rights reserved.
   See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
   for details.

      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include "itkMetaDataObject.h"
#include "itkIOCommon.h"

#include "extractNrrdVectorIndexCLP.h"
#include "BRAINSThreadControl.h"
#include <BRAINSCommonLib.h>
#include "DWIConvertLib.h"

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  bool debug = true;
  if (debug)
  {
    std::cout << "Input Image: " << inputVolume << std::endl;
    std::cout << "Output Image: " << outputVolume << std::endl;
    std::cout << "Vector Index: " << vectorIndex << std::endl;
    std::cout << "Set Image Orientation: " << setImageOrientation << std::endl << std::endl;
  }

  bool violated = false;
  if (inputVolume.empty())
  {
    violated = true;
    std::cout << "  --inputVolume Required! " << std::endl;
  }
  if (outputVolume.empty())
  {
    violated = true;
    std::cout << "  --outputVolume Required! " << std::endl;
  }
  if (violated)
  {
    return EXIT_FAILURE;
  }

  std::string convertedVolume;
  if (convertInputVolumeToNrrdOrNifti(detectOuputVolumeType(outputVolume), inputVolume, convertedVolume))
  {
    inputVolume = convertedVolume;
  }
  else
  {
    std::cout << "Error: DWI Convert can not read inputVolume." << std::endl;
    return -1;
  }

  using PixelType = signed short;
  using NrrdImageType = itk::VectorImage<PixelType, 3>;

  using FileReaderType = itk::ImageFileReader<NrrdImageType, itk::DefaultConvertPixelTraits<PixelType>>;
  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName(inputVolume);
  reader->Update();

  if ((vectorIndex < 0) || (vectorIndex >= static_cast<int>((reader->GetOutput())->GetVectorLength())))
  {
    const int maxIndex = reader->GetOutput()->GetVectorLength() - 1;
    std::cerr << "Invalid vector image index (" << vectorIndex << "), valid indexes are 0-" << maxIndex << std::endl;
    return 1;
  }

  using IndexImageType = itk::Image<PixelType, 3>;
  using VectorSelectFilterType = itk::VectorIndexSelectionCastImageFilter<NrrdImageType, IndexImageType>;
  using VectorSelectFilterPointer = VectorSelectFilterType::Pointer;

  VectorSelectFilterPointer SelectIndexImageFilter = VectorSelectFilterType::New();
  SelectIndexImageFilter->SetIndex(vectorIndex);
  SelectIndexImageFilter->SetInput(reader->GetOutput());
  try
  {
    SelectIndexImageFilter->Update();
  }
  catch (const itk::ExceptionObject & e)
  {
    std::cout << e << std::endl;
  }

  /* Hack Required for Certain Output Image Types */
  itk::MetaDataDictionary       meta;
  IndexImageType::Pointer       indexImage = SelectIndexImageFilter->GetOutput();
  IndexImageType::DirectionType fixImageDir = indexImage->GetDirection();
#define EncapsulateMD(image, flag)                                                                                     \
  {}
  if (setImageOrientation == "Axial" || setImageOrientation == "AXIAL" || setImageOrientation == "axial")
  {
    fixImageDir.Fill(0);
    fixImageDir[0][0] = 1;
    fixImageDir[1][1] = 1;
    fixImageDir[2][2] = 1;
    indexImage->SetDirection(fixImageDir);
    EncapsulateMD(indexImage, itk::SpatialOrientationEnums::ITK_COORDINATE_ORIENTATION_RPI);
  }
  else if (setImageOrientation == "Coronal" || setImageOrientation == "CORONAL" || setImageOrientation == "coronal")
  {
    fixImageDir.Fill(0);
    fixImageDir[0][0] = 1;
    fixImageDir[1][2] = 1;
    fixImageDir[2][1] = 1;
    indexImage->SetDirection(fixImageDir);
    EncapsulateMD(indexImage, itk::AnatomicalOrientation::PositiveEnum::RIP);
  }
  else if (setImageOrientation == "Sagittal" || setImageOrientation == "SAGITTAL" || setImageOrientation == "sagittal")
  {
    fixImageDir.Fill(0);
    fixImageDir[0][2] = 1;
    fixImageDir[1][1] = 1;
    fixImageDir[2][0] = 1;
    indexImage->SetDirection(fixImageDir);
    EncapsulateMD(indexImage, itk::SpatialOrientationEnums::ITK_COORDINATE_ORIENTATION_PIR);
  }
#undef EncapsulateMD // (image, flag)
  // else, leave it AsAcquired.

  using WriterType = itk::ImageFileWriter<IndexImageType>;
  WriterType::Pointer imageWriter = WriterType::New();
  imageWriter->UseCompressionOn();
  imageWriter->SetInput(indexImage);
  imageWriter->SetFileName(outputVolume);

  try
  {
    imageWriter->Update();
  }
  catch (const itk::ExceptionObject & e)
  {
    std::cout << e << std::endl;
  }

  // IndexImageType::Pointer indexImage = SelectIndexImageFilter->GetOutput();
  // std::cout << "Index Image: "<< indexImage  << std::endl;

  return EXIT_SUCCESS;
}
