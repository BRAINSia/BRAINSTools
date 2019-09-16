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
#include <itkNearestNeighborInterpolateImageFunction.h>
#include "itkGtractImageIO.h"

#include "gtractResampleCodeImageCLP.h"
#include "itkThinPlateR2LogRSplineKernelTransform.h"
#include "GenericTransformImage.h"
#include "BRAINSThreadControl.h"
int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  bool                                                  debug = true;
  if (debug)
  {
    std::cout << "==============================================================" << std::endl;
    std::cout << "Input Image: " << inputCodeVolume << std::endl;
    std::cout << "Output Image: " << outputVolume << std::endl;
    std::cout << "Reference File: " << inputReferenceVolume << std::endl;
    std::cout << "Transform File: " << inputTransform << std::endl;
    std::cout << "Transform Type: " << transformType << std::endl;
    std::cout << "==============================================================" << std::endl;
  }

  bool violated = false;
  if (inputCodeVolume.size() == 0)
  {
    violated = true;
    std::cout << "  --inputCodeVolume Required! " << std::endl;
  }
  if (inputReferenceVolume.size() == 0)
  {
    violated = true;
    std::cout << "  --inputReferenceVolume Required! " << std::endl;
  }
  if (inputTransform.size() == 0)
  {
    violated = true;
    std::cout << "  --inputTransform Required! " << std::endl;
  }
  if (outputVolume.size() == 0)
  {
    violated = true;
    std::cout << "  --outputVolume Required! " << std::endl;
  }
  if (violated)
  {
    return EXIT_FAILURE;
  }

  using CodePixelType = signed short;

  using CodeImageType = itk::Image<CodePixelType, 3>;
  using CodeImageReaderType = itk::ImageFileReader<CodeImageType>;
  CodeImageReaderType::Pointer codeImageReader = CodeImageReaderType::New();
  codeImageReader->SetFileName(inputCodeVolume);

  try
  {
    codeImageReader->Update();
  }
  catch (itk::ExceptionObject & ex)
  {
    std::cout << ex << std::endl;
    throw;
  }

  using PixelType = signed short;

  using ImageType = itk::Image<PixelType, 3>;
  using ReferenceImageReaderType = itk::ImageFileReader<ImageType>;
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

  using OrientFilterType = itk::OrientImageFilter<CodeImageType, ImageType>;
  OrientFilterType::Pointer orientImageFilter = OrientFilterType::New();
  orientImageFilter->SetInput(referenceImageReader->GetOutput());
  orientImageFilter->SetDesiredCoordinateDirection(codeImageReader->GetOutput()->GetDirection());
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

  // Read the transform
  using GenericTransformType = itk::Transform<double, 3, 3>;

  GenericTransformType::Pointer baseTransform = itk::ReadTransformFromDisk(inputTransform);
  using InterpolatorFunctionType = itk::NearestNeighborInterpolateImageFunction<CodeImageType, double>;
  InterpolatorFunctionType::Pointer interpolatorFunction = InterpolatorFunctionType::New();

  using ResampleFilterType = itk::ResampleImageFilter<CodeImageType, CodeImageType>;
  ResampleFilterType::Pointer resample = ResampleFilterType::New();

  {
    resample->SetTransform(baseTransform);
  }
  std::cout << "Code image:  ";
  codeImageReader->GetOutput()->Print(std::cout);
  std::cout << "Reference image:  ";
  referenceImageReader->GetOutput()->Print(std::cout);

  CodeImageType::PointType        p1;
  CodeImageType::PointType        p2;
  itk::ContinuousIndex<double, 3> imageIndex;
  imageIndex[0] = 0;
  imageIndex[1] = 0;
  imageIndex[2] = 0;
  codeImageReader->GetOutput()->TransformContinuousIndexToPhysicalPoint(imageIndex, p1);
  p2 = resample->GetTransform()->TransformPoint(p1);
  std::cout << "Point " << p1 << " mapped to " << p2 << std::endl;

  resample->SetInput(codeImageReader->GetOutput());
  resample->SetOutputParametersFromImage(orientImageFilter->GetOutput());
  resample->SetDefaultPixelValue(0);
  resample->SetInterpolator(interpolatorFunction);
  try
  {
    resample->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    throw;
  }

  using OrientCodeFilterType = itk::OrientImageFilter<CodeImageType, CodeImageType>;
  OrientCodeFilterType::Pointer orientCodeImageFilter = OrientCodeFilterType::New();
  orientCodeImageFilter->SetInput(resample->GetOutput());
  orientCodeImageFilter->SetDesiredCoordinateDirection(codeImageReader->GetOutput()->GetDirection());
  orientCodeImageFilter->UseImageDirectionOn();
  try
  {
    orientCodeImageFilter->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cout << e << std::endl;
    throw;
  }

  CodeImageType::Pointer resampledImage = orientCodeImageFilter->GetOutput();
  resampledImage->SetMetaDataDictionary(codeImageReader->GetOutput()->GetMetaDataDictionary());

  using ImageFileWriterType = itk::ImageFileWriter<CodeImageType>;
  ImageFileWriterType::Pointer ImageWriter = ImageFileWriterType::New();
  ImageWriter->UseCompressionOn();
  ImageWriter->SetFileName(outputVolume);
  ImageWriter->SetInput(resampledImage);
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
