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
#include "itksys/SystemTools.hxx"
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkSubtractImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <vcl_compiler.h>
#include <iostream>
#include <algorithm>
#include "DWIConvertUtils.h"
#include "DWISimpleCompareCLP.h"
#include <BRAINSCommonLib.h>

namespace
{
#define DIMENSION 4

template <typename PixelType>
std::vector<std::vector<double>>
RecoverGVector(typename itk::Image<PixelType, DIMENSION>::Pointer & img)
{
  std::vector<std::vector<double>> rval;

  itk::MetaDataDictionary & dict = img->GetMetaDataDictionary();

  for (unsigned curGradientVec = 0;; ++curGradientVec)
  {
    std::stringstream labelSS;
    labelSS << "DWMRI_gradient_" << std::setw(4) << std::setfill('0') << curGradientVec;
    std::string valString;
    // look for gradients in metadata until none by current name exists
    if (!itk::ExposeMetaData<std::string>(dict, labelSS.str(), valString))
    {
      break;
    }
    std::stringstream   valSS(valString);
    std::vector<double> vec;
    for (;;)
    {
      double curVal;
      valSS >> curVal;
      if (!valSS.fail())
      {
        vec.push_back(curVal);
      }
      else
      {
        break;
      }
    }
    rval.push_back(vec);
  }
  return rval;
}

template <typename PixelType>
int
DoIt(const std::string & inputVolume1, const std::string & inputVolume2, PixelType, bool CheckDWIData)
{
  int rval(EXIT_SUCCESS);

  using ImageType = itk::VectorImage<PixelType, DIMENSION>;
  using FileReaderType = itk::ImageFileReader<ImageType>;

  typename FileReaderType::Pointer firstReader = FileReaderType::New();
  typename FileReaderType::Pointer secondReader = FileReaderType::New();

  firstReader->SetFileName(inputVolume1.c_str());
  secondReader->SetFileName(inputVolume2.c_str());

  firstReader->Update();
  secondReader->Update();
  typename ImageType::Pointer firstImage = firstReader->GetOutput();
  typename ImageType::Pointer secondImage = secondReader->GetOutput();
  // check origin -- problem with file conversion causing some drift
  typename ImageType::PointType firstOrigin(firstImage->GetOrigin());
  typename ImageType::PointType secondOrigin(secondImage->GetOrigin());
  double                        distance = std::sqrt(firstOrigin.SquaredEuclideanDistanceTo(secondOrigin));
  if (distance > 1.0E-3)
  {
    std::cerr << "Origins differ " << firstOrigin << " " << secondOrigin << std::endl;
    return EXIT_FAILURE;
  }
  else if (distance > 1.0E-6)
  {
    // if there is a small difference make them the same
    firstImage->SetOrigin(secondOrigin);
  }
  // same deal with spacing, can be slightly off due to numerical error
  typename ImageType::SpacingType firstSpacing(firstImage->GetSpacing());
  typename ImageType::SpacingType secondSpacing(secondImage->GetSpacing());
  for (unsigned int i = 0; i < ImageType::GetImageDimension(); ++i)
  {
    double diff = std::fabs(firstSpacing[i] - secondSpacing[i]);
    if (diff > 1.0e-6 && diff < 1.0e-4)
    {
      firstSpacing[i] = secondSpacing[i];
    }
  }
  firstImage->SetSpacing(firstSpacing);
#if 1
  itk::ImageRegionConstIterator<ImageType> firstIt(firstImage, firstImage->GetLargestPossibleRegion());

  itk::ImageRegionConstIterator<ImageType> secondIt(secondImage, secondImage->GetLargestPossibleRegion());
  unsigned                                 count = 0;
  unsigned                                 errCount = 0;
  for (firstIt.GoToBegin(), secondIt.GoToBegin(); !firstIt.IsAtEnd() && !secondIt.IsAtEnd();
       ++firstIt, ++secondIt, ++count)
  {
    if (!CloseEnough(firstIt.Get(), secondIt.Get()))
    {
      std::cerr << "Images don't match at voxel " << count << std::endl
                << firstIt.Get() << std::endl
                << secondIt.Get() << std::endl;
      rval = EXIT_FAILURE;
      if (++errCount == 5)
      {
        break;
      }
    }
  }
#else
  using SubtractFilterType = itk::SubtractImageFilter<ImageType>;
  typename SubtractFilterType::Pointer subtractFilter = SubtractFilterType::New();
  subtractFilter->SetInput1(firstReader->GetOutput());
  subtractFilter->SetInput2(secondReader->GetOutput());
  subtractFilter->Update();

  typename ImageType::Pointer subtractImage = subtractFilter->GetOutput();

  using StatisticsFilterType = itk::StatisticsImageFilter<ImageType>;
  typename StatisticsFilterType::Pointer statisticsFilter = StatisticsFilterType::New();
  statisticsFilter->SetInput(subtractImage);
  try
  {
    statisticsFilter->Update();
  }
  catch (const itk::ExceptionObject & e)
  {
    std::cerr << "Exception detected while comparing " << inputVolume1 << "  " << inputVolume2 << e.GetDescription();
    return EXIT_FAILURE;
  }
  if (std::fabs(static_cast<float>(statisticsFilter->GetMaximum())) > 0.0001 ||
      std::fabs(static_cast<float>(statisticsFilter->GetMinimum())) > 0.0001)
  {
    std::cerr << "Image Data Differs -- min diff " << statisticsFilter->GetMinimum() << " max diff "
              << statisticsFilter->GetMaximum() << std::endl;
    using ImageWriter = typename itk::ImageFileWriter<ImageType>;
    typename ImageWriter::Pointer writer = ImageWriter::New();
    writer->SetInput(subtractImage);
    std::string filename = itksys::SystemTools::GetFilenameWithoutExtension(inputVolume1);
    filename = itksys::SystemTools::GetFilenameName(filename);
    filename += "-diff.nrrd";
    writer->SetFileName(filename);
    writer->Write();
    rval = EXIT_FAILURE;
  }
#endif
  if (!CheckDWIData)
  {
    return rval;
  }
  double bVal1;

  double bVal2;
  if (RecoverBValue<ImageType>(firstImage, bVal1) != EXIT_SUCCESS)
  {
    std::cerr << "Missing BValue in " << inputVolume1 << std::endl;
    return EXIT_FAILURE;
  }

  if (RecoverBValue<ImageType>(secondImage, bVal2) != EXIT_SUCCESS)
  {
    std::cerr << "Missing BValue in " << inputVolume2 << std::endl;
    return EXIT_FAILURE;
  }

  if (!CloseEnough(bVal1, bVal2))
  {
    std::cerr << "BValue mismatch: " << bVal1 << " " << bVal2 << std::endl;
    rval = EXIT_FAILURE;
  }

  DWIMetaDataDictionaryValidator::GradientTableType firstGVector;
  RecoverBVectors<ImageType>(firstImage, firstGVector);
  DWIMetaDataDictionaryValidator::GradientTableType secondGVector;
  RecoverBVectors<ImageType>(secondImage, secondGVector);
  if (firstGVector.size() != secondGVector.size())
  {
    std::cerr << "First image Gradient Vectors size (" << firstGVector.size()
              << ") doesn't match second image Gradient vectors size (" << secondGVector.size() << ")" << std::endl;
    return EXIT_FAILURE;
  }
  if (!CloseEnough(firstGVector, secondGVector))
  {
    std::cerr << "Gradient vectors don't match" << std::endl;
    std::cerr << "First Vector ";
    PrintVec(firstGVector);
    std::cerr << "Second Vector ";
    PrintVec(secondGVector);

    rval = EXIT_FAILURE;
  }
  return rval;
}

void
GetImageType(const std::string &                 fileName,
             itk::IOPixelEnum &     pixelType,
             itk::ImageIOBase::IOComponentEnum & componentType)
{
  using ImageType = itk::Image<short, 3>;
  itk::ImageFileReader<ImageType>::Pointer imageReader = itk::ImageFileReader<ImageType>::New();
  imageReader->SetFileName(fileName.c_str());
  imageReader->UpdateOutputInformation();

  // const unsigned int componentNum = imageReader->GetImageIO()->GetNumberOfComponents();
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
    GetImageType(inputVolume2, pixelType, componentType);

    // This filter handles all types

    int rval(EXIT_FAILURE);

    switch (componentType)
    {
      case itk::IOComponentEnum::UCHAR:
        rval = DoIt(inputVolume1, inputVolume2, static_cast<unsigned char>(0), CheckDWIData);
        break;
      case itk::IOComponentEnum::CHAR:
        rval = DoIt(inputVolume1, inputVolume2, static_cast<char>(0), CheckDWIData);
        break;
      case itk::IOComponentEnum::USHORT:
        rval = DoIt(inputVolume1, inputVolume2, static_cast<unsigned short>(0), CheckDWIData);
        break;
      case itk::IOComponentEnum::SHORT:
        rval = DoIt(inputVolume1, inputVolume2, static_cast<short>(0), CheckDWIData);
        break;
      case itk::IOComponentEnum::UINT:
        rval = DoIt(inputVolume1, inputVolume2, static_cast<unsigned int>(0), CheckDWIData);
        break;
      case itk::IOComponentEnum::INT:
        rval = DoIt(inputVolume1, inputVolume2, static_cast<int>(0), CheckDWIData);
        break;
      case itk::IOComponentEnum::ULONG:
        rval = DoIt(inputVolume1, inputVolume2, static_cast<unsigned long>(0), CheckDWIData);
        break;
      case itk::IOComponentEnum::LONG:
        rval = DoIt(inputVolume1, inputVolume2, static_cast<long>(0), CheckDWIData);
        break;
      case itk::IOComponentEnum::FLOAT:
        rval = DoIt(inputVolume1, inputVolume2, 0.0f, CheckDWIData);
        // std::cout << "FLOAT type not currently supported." << std::endl;
        break;
      case itk::IOComponentEnum::DOUBLE:
        std::cout << "DOUBLE type not currently supported." << std::endl;
        break;
      case itk::IOComponentEnum::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
    }
    if (rval == EXIT_FAILURE)
    {
      return EXIT_FAILURE;
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
