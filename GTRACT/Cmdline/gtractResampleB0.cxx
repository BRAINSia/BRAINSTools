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
#include <itkExtractImageFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkTransformFactory.h>
#include <itkVersorRigid3DTransform.h>
#include <itkOrientImageFilter.h>
#include "BRAINSThreadControl.h"
#include "itkGtractImageIO.h"
#include "gtractResampleB0CLP.h"
#include "GenericTransformImage.h"
#include "DWIConvertLib.h"

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  bool                                                  debug = true;
  if (debug)
  {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Input Image: " << inputVolume << std::endl;
    std::cout << "Output Image: " << outputVolume << std::endl;
    std::cout << "Anatomical Image: " << inputAnatomicalVolume << std::endl;
    std::cout << "Transform File: " << inputTransform << std::endl;
    std::cout << "Transform Type: " << transformType << std::endl;
    std::cout << "Vector Index: " << vectorIndex << std::endl;
    std::cout << "=====================================================" << std::endl;
  }

  bool violated = false;
  if (inputVolume.empty())
  {
    violated = true;
    std::cout << "  --inputVolume Required! " << std::endl;
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

  using VectorImageType = itk::VectorImage<PixelType, 3>;
  using VectorImageReaderType = itk::ImageFileReader<VectorImageType, itk::DefaultConvertPixelTraits<PixelType>>;
  VectorImageReaderType::Pointer vectorImageReader = VectorImageReaderType::New();
  vectorImageReader->SetFileName(inputVolume);

  try
  {
    vectorImageReader->Update();
  }
  catch (itk::ExceptionObject & ex)
  {
    std::cout << ex << std::endl;
    throw;
  }

  using ImageType = itk::Image<PixelType, 3>;
  using AnatomicalImageReaderType = itk::ImageFileReader<ImageType>;
  AnatomicalImageReaderType::Pointer anatomicalReader = AnatomicalImageReaderType::New();
  anatomicalReader->SetFileName(inputAnatomicalVolume);

  try
  {
    anatomicalReader->Update();
  }
  catch (itk::ExceptionObject & ex)
  {
    std::cout << ex << std::endl;
    throw;
  }

  // Read the transform
  using GenericTransformType = itk::Transform<double, 3, 3>;
  GenericTransformType::Pointer baseTransform = itk::ReadTransformFromDisk(inputTransform);

  using VectorImageType = itk::VectorImage<PixelType, 3>;
  VectorImageType::Pointer    dwiImage = vectorImageReader->GetOutput();
  VectorImageType::RegionType fixedRegion = dwiImage->GetLargestPossibleRegion();
  VectorImageType::SizeType   fixedSize = fixedRegion.GetSize();
  // const VectorImageType::SpacingType fixedSpacing = dwiImage->GetSpacing();
  // const VectorImageType::PointType   fixedOrigin = dwiImage->GetOrigin();

  fixedSize[3] = 0;
  fixedRegion.SetSize(fixedSize);

  VectorImageType::IndexType fixedIndex = fixedRegion.GetIndex();
  fixedIndex[0] = 0;
  fixedIndex[1] = 0;
  fixedIndex[2] = 0;
  fixedIndex[3] = 0;
  fixedRegion.SetIndex(fixedIndex);

  /* Extract the Vector Image Index for Registration */
  using VectorSelectFilterType = itk::VectorIndexSelectionCastImageFilter<VectorImageType, ImageType>;
  using VectorSelectFilterPointer = VectorSelectFilterType::Pointer;
  VectorSelectFilterPointer selectIndexImageFilter = VectorSelectFilterType::New();
  selectIndexImageFilter->SetIndex(vectorIndex);
  selectIndexImageFilter->SetInput(dwiImage);
  try
  {
    selectIndexImageFilter->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cout << e << std::endl;
    throw;
  }

  using ResampleFilterType = itk::ResampleImageFilter<ImageType, ImageType>;

  ResampleFilterType::Pointer resample = ResampleFilterType::New();

  resample->SetTransform(baseTransform);

  resample->SetInput(selectIndexImageFilter->GetOutput());
  resample->SetOutputParametersFromImage(anatomicalReader->GetOutput());
  resample->SetDefaultPixelValue(1);
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

  ImageType::Pointer resampledImage = resample->GetOutput();
  resampledImage->SetMetaDataDictionary(anatomicalReader->GetOutput()->GetMetaDataDictionary());

  using ImageFileWriterType = itk::ImageFileWriter<ImageType>;
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
  std::cout << "wrote image" << std::endl;
  return EXIT_SUCCESS;
}
