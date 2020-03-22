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

#include <iostream>
#include <fstream>

#include <itkImage.h>
#include <itkOrientImageFilter.h>

#include "itkGtractImageIO.h"
#include "gtractImageConformityCLP.h"
#include "GenericTransformImage.h"
#include "BRAINSThreadControl.h"

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
    std::cout << "Reference Image: " << inputReferenceVolume << std::endl;
  }

  bool violated = false;
  if (inputVolume.empty())
  {
    violated = true;
    std::cout << "  --inputVolume Required! " << std::endl;
  }
  if (inputReferenceVolume.empty())
  {
    violated = true;
    std::cout << "  --inputReferenceVolume Required! " << std::endl;
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

  using SpecimenImageType = itk::Image<PixelType, 3>;
  using SpecimenImageReaderType = itk::ImageFileReader<SpecimenImageType>;
  SpecimenImageReaderType::Pointer specimenImageReader = SpecimenImageReaderType::New();
  specimenImageReader->SetFileName(inputVolume);

  try
  {
    specimenImageReader->Update();
  }
  catch (itk::ExceptionObject & ex)
  {
    std::cout << ex << std::endl;
    throw;
  }

  using ReferenceImageType = itk::Image<PixelType, 3>;
  using ReferenceImageReaderType = itk::ImageFileReader<ReferenceImageType>;
  ReferenceImageReaderType::Pointer referenceImageReader = ReferenceImageReaderType::New();
  referenceImageReader->SetFileName(inputReferenceVolume);

  try
  {
    referenceImageReader->Update();
  }
  catch (itk::ExceptionObject & ex)
  {
    std::cout << ex << std::endl;
    throw;
  }

  using OrientFilterType = itk::OrientImageFilter<SpecimenImageType, ReferenceImageType>;
  OrientFilterType::Pointer orientImageFilter = OrientFilterType::New();
  orientImageFilter->SetInput(specimenImageReader->GetOutput());
  orientImageFilter->SetDesiredCoordinateDirection(referenceImageReader->GetOutput()->GetDirection());
  orientImageFilter->UseImageDirectionOn();
  try
  {
    orientImageFilter->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cout << e << std::endl;
    throw;
  }

  ReferenceImageType::Pointer reorientedImage = orientImageFilter->GetOutput();
  reorientedImage->SetOrigin(referenceImageReader->GetOutput()->GetOrigin());

  reorientedImage->SetMetaDataDictionary(specimenImageReader->GetOutput()->GetMetaDataDictionary());

  using ImageFileWriterType = itk::ImageFileWriter<ReferenceImageType>;
  ImageFileWriterType::Pointer ImageWriter = ImageFileWriterType::New();
  ImageWriter->UseCompressionOn();
  ImageWriter->SetFileName(outputVolume);
  ImageWriter->SetInput(reorientedImage);
  try
  {
    ImageWriter->Update();
  }
  catch (itk::ExceptionObject & ex)
  {
    std::cout << ex << std::endl;
    throw;
  }
  return EXIT_SUCCESS;
}
