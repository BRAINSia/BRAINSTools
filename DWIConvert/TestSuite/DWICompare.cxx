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

#include <iostream>
#include <algorithm>
#include <string>
#include <itkMetaDataObject.h>
#include <itkImage.h>
#include <itkVector.h>
#include <itkVectorImage.h>

#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkNrrdImageIO.h>

#include <itkImageRegionConstIterator.h>

#include <itkMath.h>

#include "DWICompareCLP.h"
#include <BRAINSCommonLib.h>
#include "DWIMetaDataDictionaryValidator.h"

#define DIMENSION 3


namespace
{

template <typename ImageType>
bool
TestIfInformationIsDifferent(typename ImageType::ConstPointer first, typename ImageType::ConstPointer second)
{
  bool        failureStatus = false;
  const float spacingTolerance = 1e-3;
  if (!first->GetSpacing().GetVnlVector().is_equal(second->GetSpacing().GetVnlVector(), spacingTolerance))
  {
    std::cout << "The first image Spacing does not match second image Information" << std::endl;
    std::cout << "First Spacing: " << first->GetSpacing() << std::endl;
    std::cout << "Second Spacing: " << second->GetSpacing() << std::endl;
    failureStatus = true;
  }
  const float originTolerance = 1e-2;
  if (!first->GetOrigin().GetVnlVector().is_equal(second->GetOrigin().GetVnlVector(), originTolerance))
  {
    std::cout << "The first image Origin does not match second image Information" << std::endl;
    std::cout << "First Origin: " << first->GetOrigin() << std::endl;
    std::cout << "Second Origin: " << second->GetOrigin() << std::endl;
    failureStatus = true;
  }
  const float directionTolerance = 1e-3;
  if (!first->GetDirection().GetVnlMatrix().as_ref().is_equal(second->GetDirection().GetVnlMatrix().as_matrix(),
                                                              directionTolerance))
  {
    std::cout << "The first image Direction does not match second image Information" << std::endl;
    std::cout << "First Direction: " << first->GetDirection() << std::endl;
    std::cout << "Second Direction: " << second->GetDirection() << std::endl;
    failureStatus = true;
  }
  if (first->GetLargestPossibleRegion() != second->GetLargestPossibleRegion())
  {
    std::cout << "The first image Size does not match second image Information" << std::endl;
    std::cout << "first: " << first->GetLargestPossibleRegion() << "updated: " << second->GetLargestPossibleRegion()
              << std::endl;
    failureStatus = true;
  }

  using firstVectorImageIterator = itk::ImageRegionConstIterator<ImageType>;
  firstVectorImageIterator itr1(first, first->GetLargestPossibleRegion());
  itr1.GoToBegin();

  using secondVectorImageIterator = itk::ImageRegionConstIterator<ImageType>;
  secondVectorImageIterator itr2(second, second->GetLargestPossibleRegion());
  itr2.GoToBegin();

  while (!itr1.IsAtEnd())
  {
    if ((itr2.IsAtEnd()) || (itr1.Get() != itr2.Get()))
    {
      if (!failureStatus)
      {
        std::cout << "ERROR: Pixel values are different " << std::endl;
      }
      failureStatus = true;
    }
    ++itr1;
    ++itr2;
  }

  return failureStatus;
}

template <typename PixelType>
int
DoIt(int argc, char * argv[], PixelType)
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  using DiffusionImageType = itk::VectorImage<PixelType, DIMENSION>;

  using FileReaderType = itk::ImageFileReader<DiffusionImageType>;
  typename FileReaderType::Pointer firstReader = FileReaderType::New();
  typename FileReaderType::Pointer secondReader = FileReaderType::New();
  firstReader->SetFileName(inputVolume1);
  firstReader->Update();
  typename DiffusionImageType::ConstPointer firstImage = firstReader->GetOutput();
  secondReader->SetFileName(inputVolume2);
  secondReader->Update();
  typename DiffusionImageType::ConstPointer secondImage = secondReader->GetOutput();

  bool failure = TestIfInformationIsDifferent<DiffusionImageType>(firstImage, secondImage);

  // If the angle between two gradients differs more than this value they are
  // considered to be non-colinear
  // const double gradientToleranceForSameness = 1;
  constexpr float bValueTolerance = .05;

  using DictionaryType = itk::MetaDataDictionary;
  const DictionaryType &         firstDictionary = firstReader->GetMetaDataDictionary();
  DWIMetaDataDictionaryValidator firstValidator;
  firstValidator.SetMetaDataDictionary(firstDictionary);
  const DictionaryType &         secondDictionary = secondReader->GetMetaDataDictionary();
  DWIMetaDataDictionaryValidator secondValidator;
  secondValidator.SetMetaDataDictionary(secondDictionary);

  {
    if (std::abs(firstValidator.GetBValue() - secondValidator.GetBValue()) / secondValidator.GetBValue() >
        bValueTolerance)
    {
      std::cerr << "firstBValue String != secondBValueString! " << firstValidator.GetBValue()
                << "!=" << secondValidator.GetBValue() << std::endl;
      return EXIT_FAILURE;
    }
  }
  DWIMetaDataDictionaryValidator::RotationMatrixType fMF = firstValidator.GetMeasurementFrame();
  DWIMetaDataDictionaryValidator::RotationMatrixType sMF = secondValidator.GetMeasurementFrame();
  if (!useIdentityMeasurementFrame)
  {
    for (size_t j = 0; j < 3; ++j)
    {
      for (size_t i = 0; i < 3; ++i)
      {
        if (std::abs(fMF[i][j] - sMF[i][j]) > 1e-5)
        {
          std::cerr << "Measurement Frames do not match:\n" << fMF << "\n != \n" << sMF << "\n" << std::endl;
          return EXIT_FAILURE;
        }
      }
    }
  }

  DWIMetaDataDictionaryValidator::GradientTableType fGD = firstValidator.GetGradientTable();
  DWIMetaDataDictionaryValidator::GradientTableType sGD = secondValidator.GetGradientTable();
  for (size_t idx = 0; idx < fGD.size(); ++idx)
  {
    DWIMetaDataDictionaryValidator::GradientTableType::value_type ufVec = fGD[idx];
    DWIMetaDataDictionaryValidator::GradientTableType::value_type usVec = sGD[idx];
    DWIMetaDataDictionaryValidator::GradientTableType::value_type fVec = fGD[idx];
    DWIMetaDataDictionaryValidator::GradientTableType::value_type sVec = sGD[idx];
    if (useIdentityMeasurementFrame)
    {
      fVec = fMF.GetInverse() * ufVec;
      sVec = sMF.GetInverse() * usVec;
    }
    const double mag = (sVec - fVec).squared_magnitude();
    if (mag > 1e-5)
    {
      std::cerr << "Gradient Vector at index " << idx << " mismatch" << fVec << " != " << sVec << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (failure)
  {
    return EXIT_FAILURE;
  }
  else
  {
    return EXIT_SUCCESS;
  }
}

void
GetImageType(const std::string &                 fileName,
             itk::IOPixelEnum &     pixelType,
             itk::ImageIOBase::IOComponentEnum & componentType)
{
  using ImageType = itk::Image<short, 3>;
  itk::ImageFileReader<ImageType>::Pointer imageReader = itk::ImageFileReader<ImageType>::New();
  imageReader->SetFileName(fileName);
  imageReader->UpdateOutputInformation();

  pixelType = imageReader->GetImageIO()->GetPixelType();
  componentType = imageReader->GetImageIO()->GetComponentType();
}
} // namespace

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  itk::IOPixelEnum     pixelType;
  itk::ImageIOBase::IOComponentEnum componentType;

  try
  {
    // itk::GetImageType (inputVolume1, pixelType, componentType);
    GetImageType(inputVolume1, pixelType, componentType);

    // This filter handles all types

    switch (componentType)
    {
      case itk::IOComponentEnum::UCHAR:
      {
        return DoIt(argc, argv, static_cast<unsigned char>(0));
      }
       case itk::IOComponentEnum::CHAR:
      {
        return DoIt(argc, argv, static_cast<char>(0));
      }
      case itk::IOComponentEnum::USHORT:
      {
        return DoIt(argc, argv, static_cast<unsigned short>(0));
      }
      case itk::IOComponentEnum::SHORT:
      {
        return DoIt(argc, argv, static_cast<short>(0));
      }
      case itk::IOComponentEnum::UINT:
      {
        return DoIt(argc, argv, static_cast<unsigned int>(0));
      }
      case itk::IOComponentEnum::INT:
      {
        return DoIt(argc, argv, static_cast<int>(0));
      }
      case itk::IOComponentEnum::ULONG:
      {
        return DoIt(argc, argv, static_cast<unsigned long>(0));
      }
      case itk::IOComponentEnum::LONG:
      {
        return DoIt(argc, argv, static_cast<long>(0));
      }
      case itk::IOComponentEnum::FLOAT:
      {
        return DoIt(argc, argv, 0.0f);
        // std::cout << "FLOAT type not currently supported." << std::endl;
      }
      case itk::IOComponentEnum::DOUBLE:
      {
        std::cout << "DOUBLE type not currently supported." << std::endl;
      }
      break;
      case itk::IOComponentEnum::UNKNOWNCOMPONENTTYPE:
      default:
      {
        std::cout << "unknown component type" << std::endl;
      }
      break;
    }
  }
  catch (const itk::ExceptionObject & excep)
  {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
