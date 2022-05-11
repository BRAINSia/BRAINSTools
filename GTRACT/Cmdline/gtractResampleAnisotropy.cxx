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
#include <itkResampleImageFilter.h>
#include <itkTransformFactory.h>
#include <itkVersorRigid3DTransform.h>
#include <itkOrientImageFilter.h>

#include "itkGtractImageIO.h"
#include "gtractResampleAnisotropyCLP.h"
#include "GenericTransformImage.h"
#include "BRAINSThreadControl.h"
int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  bool debug = true;
  if (debug)
  {
    std::cout << "==============================================================" << std::endl;
    std::cout << "Input Image: " << inputAnisotropyVolume << std::endl;
    std::cout << "Output Image: " << outputVolume << std::endl;
    std::cout << "Anatomical Image: " << inputAnatomicalVolume << std::endl;
    std::cout << "Transform File: " << inputTransform << std::endl;
    std::cout << "Transform Type: " << transformType << std::endl;
    std::cout << "==============================================================" << std::endl;
  }

  bool violated = false;
  if (inputAnisotropyVolume.empty())
  {
    violated = true;
    std::cout << "  --inputAnisotropyVolume Required! " << std::endl;
  }
  if (inputAnatomicalVolume.empty())
  {
    violated = true;
    std::cout << "  --inputAnatomicalVolume Required! " << std::endl;
  }
  if (inputTransform.empty())
  {
    violated = true;
    std::cout << "  --inputTransform Required! " << std::endl;
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

  using AnisotropyPixelType = float;

  using AnisotropyImageType = itk::Image<AnisotropyPixelType, 3>;
  using AnisotropyImageReaderType = itk::ImageFileReader<AnisotropyImageType>;
  AnisotropyImageReaderType::Pointer anisotropyImageReader = AnisotropyImageReaderType::New();
  anisotropyImageReader->SetFileName(inputAnisotropyVolume);

  try
  {
    anisotropyImageReader->Update();
  }
  catch (const itk::ExceptionObject & ex)
  {
    std::cout << ex << std::endl;
    throw;
  }

  using PixelType = signed short;

  using ImageType = itk::Image<PixelType, 3>;
  using AnatomicalImageReaderType = itk::ImageFileReader<ImageType>;
  AnatomicalImageReaderType::Pointer anatomicalReader = AnatomicalImageReaderType::New();
  anatomicalReader->SetFileName(inputAnatomicalVolume);

  try
  {
    anatomicalReader->Update();
  }
  catch (const itk::ExceptionObject & ex)
  {
    std::cout << ex << std::endl;
    throw;
  }

  // Read the transform
  using GenericTransformType = itk::Transform<double, 3, 3>;
  GenericTransformType::Pointer baseTransform = itk::ReadTransformFromDisk(inputTransform);

  using AnisotropyImageType = itk::Image<float, 3>;
  using ResampleFilterType = itk::ResampleImageFilter<AnisotropyImageType, AnisotropyImageType>;

  ResampleFilterType::Pointer resample = ResampleFilterType::New();
  {
    resample->SetTransform(baseTransform);
  }
  resample->SetInput(anisotropyImageReader->GetOutput());
  resample->SetOutputParametersFromImage(anatomicalReader->GetOutput());
  resample->SetDefaultPixelValue(0);
  try
  {
    resample->Update();
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    throw;
  }

  AnisotropyImageType::Pointer resampledImage = resample->GetOutput();
  resampledImage->SetMetaDataDictionary(anatomicalReader->GetOutput()->GetMetaDataDictionary());

  using ImageFileWriterType = itk::ImageFileWriter<AnisotropyImageType>;
  ImageFileWriterType::Pointer ImageWriter = ImageFileWriterType::New();
  ImageWriter->UseCompressionOn();
  ImageWriter->SetFileName(outputVolume);
  ImageWriter->SetInput(resampledImage);
  try
  {
    ImageWriter->Update();
  }
  catch (const itk::ExceptionObject & ex)
  {
    std::cout << ex << std::endl;
    throw;
  }
  return EXIT_SUCCESS;
}
