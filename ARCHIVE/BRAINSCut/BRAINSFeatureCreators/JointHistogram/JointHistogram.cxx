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
// TODO Output Filename

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "JointHistogramCLP.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkHistogram.h"

#include "itkScalarImageToHistogramGenerator.h"
#include "itkHistogramToProbabilityImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkInvertIntensityImageFilter.h"

#include "itkLabelStatisticsImageFilter.h"
#include <map>

#include <fstream>
#include <BRAINSCommonLib.h>

#define HISTOGRAMSIZE 255

/*
 * Author : Eun Young (Regina) Kim
 */

// ///////////////////////////////////////////////////////////////////////////////////////
/*
 * main function
 */
int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  // define image with type of voxel
  using PixelType = double;
  constexpr unsigned long int Dimension = 3;
  using InputImageType = itk::Image<PixelType, 3>;

  // there has to be two input volumes and label volume
  if ((inputVolumeInXAxis.empty()) || (inputVolumeInYAxis.empty()))
  {
    std::cout << " Wrong Argument! " << std::endl
              << " inputVolumeInXAxis, inputVolumeInYAxis, are necessary! " << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "* InputImage Filename " << inputVolumeInXAxis << std::endl
            << "* InputImage Filename " << inputVolumeInYAxis << std::endl;

  // Image Reader for inputVolumes
  using ImageReaderType = itk::ImageFileReader<InputImageType>;

  ImageReaderType::Pointer imageReader1 = ImageReaderType::New();
  imageReader1->SetFileName(inputVolumeInXAxis);

  ImageReaderType::Pointer imageReader2 = ImageReaderType::New();
  imageReader2->SetFileName(inputVolumeInYAxis);

  // Read Binary Images
  using BinaryPixelType = double;
  using BinaryImageType = itk::Image<BinaryPixelType, Dimension>;

  BinaryImageType::Pointer binaryImageInX = BinaryImageType::New();
  BinaryImageType::Pointer binaryImageInY = BinaryImageType::New();

  if ((!inputMaskVolumeInXAxis.empty()) && (!inputMaskVolumeInYAxis.empty()))
  {
    using BinaryImageReaderType = itk::ImageFileReader<BinaryImageType>;

    BinaryImageReaderType::Pointer binaryImageReaderInX = BinaryImageReaderType::New();
    binaryImageReaderInX->SetFileName(inputMaskVolumeInXAxis);

    BinaryImageReaderType::Pointer binaryImageReaderInY = BinaryImageReaderType::New();
    binaryImageReaderInY->SetFileName(inputMaskVolumeInYAxis);

    try
    {
      binaryImageReaderInX->Update();
      binaryImageReaderInY->Update();
    }
    catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception in Reading." << std::endl;
      std::cerr << e.GetDescription() << std::endl;
      std::cerr << e.GetLocation() << std::endl;
      exit(EXIT_FAILURE);
    }
    binaryImageInX = binaryImageReaderInX->GetOutput();
    binaryImageInY = binaryImageReaderInY->GetOutput();
  }

  // Rescale Input Images
  using RescaleFilterType = itk::RescaleIntensityImageFilter<InputImageType, InputImageType>;

  RescaleFilterType::Pointer rescaler1 = RescaleFilterType::New();

  rescaler1->SetInput(imageReader1->GetOutput());

  rescaler1->SetOutputMaximum(HISTOGRAMSIZE);
  rescaler1->SetOutputMinimum(0);

  RescaleFilterType::Pointer rescaler2 = RescaleFilterType::New();

  rescaler2->SetInput(imageReader2->GetOutput());
  rescaler2->SetOutputMaximum(HISTOGRAMSIZE);
  rescaler2->SetOutputMinimum(0);
  try
  {
    rescaler1->Update();
    rescaler2->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cerr << "Exception in Rescaling." << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    exit(EXIT_FAILURE);
  }

  // Images
  const InputImageType::ConstPointer imageInX = rescaler1->GetOutput();
  const InputImageType::ConstPointer imageInY = rescaler2->GetOutput();

  // Start Iterator
  unsigned int HistogramArray[HISTOGRAMSIZE + 1][HISTOGRAMSIZE + 1] = { { 0 } };

  InputImageType::IndexType currentIndexOfX;
  InputImageType::IndexType currentIndexOfY;

  itk::Point<double, Dimension> currentPhysicalPoint;

  // * Iterator For Image

  itk::ImageRegionConstIterator<InputImageType> it(imageInX, imageInX->GetLargestPossibleRegion());

  InputImageType::SizeType imageInYSize = imageInY->GetLargestPossibleRegion().GetSize();
  std::cout << "imageInYSize::" << imageInYSize << std::endl;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    currentIndexOfX = it.GetIndex();

    // get physical points
    imageInX->TransformIndexToPhysicalPoint(currentIndexOfX, currentPhysicalPoint);

    imageInY->TransformPhysicalPointToIndex(currentPhysicalPoint, currentIndexOfY);

    bool imageInBoundary = true;
    for (InputImageType::SizeType::SizeValueType dimIndex = 0; dimIndex < Dimension; dimIndex++)
    {
      if (currentIndexOfY[dimIndex] > static_cast<InputImageType::IndexType::IndexValueType>(imageInYSize[dimIndex]))
      {
        imageInBoundary = false;
      }
    }

    if (imageInBoundary)
    {
      if (verbose)
      {
        std::cout << "Transform IndexOfX ( " << currentIndexOfX << " ) to Physical Point ( " << currentPhysicalPoint
                  << " ) and Get IndexOfY ( " << currentIndexOfY << " ) " << std::endl;
      }

      // taking account binary (mask) images if they are given

      if (inputMaskVolumeInXAxis.empty() ||
          ((!inputMaskVolumeInXAxis.empty()) && (!inputMaskVolumeInYAxis.empty()) &&
           (binaryImageInX->GetPixel(currentIndexOfX) > 0) && (binaryImageInY->GetPixel(currentIndexOfY) > 0)))
      {
        int temp_image1_intensity = imageInX->GetPixel(currentIndexOfX);
        int temp_image2_intensity = imageInY->GetPixel(currentIndexOfY);
        if (verbose)
        {
          std::cout << " ADD to the Bin (" << temp_image1_intensity << " , " << temp_image2_intensity << " ) ";
        }

        HistogramArray[temp_image1_intensity][temp_image2_intensity]++;

        if (verbose)
        {
          std::cout << " = " << HistogramArray[temp_image1_intensity][temp_image2_intensity] << std::endl;
        }
      }
    }
  }
  // - open file stream and write to the file
  std::ofstream outputFileStream;
  std::string   outputHistogramData = outputJointHistogramImage + ".txt";
  outputFileStream.open(outputHistogramData.c_str());

  // - write header line 1: including image names

  outputFileStream << "[Image1]: " << inputVolumeInXAxis << std::endl
                   << "[Image2]: " << inputVolumeInYAxis << std::endl
                   << std::endl;
  // - write header line 2: colume name

  std::string FrequencyName = "Frequence";
  FrequencyName += "OfIntensity";

  outputFileStream << "label, " << inputVolumeInXAxis << ", " << inputVolumeInYAxis << ", " << FrequencyName
                   << std::endl;
  // - Iterate for each bins
  for (int i = 0; i < HISTOGRAMSIZE; i++)
  {
    for (int j = 0; j < HISTOGRAMSIZE; j++)
    {
      // - Write text file

      outputFileStream << i << "," << j << "," << HistogramArray[i][j] << std::endl;
    }
  }
  outputFileStream.close();
  // - Write Histogram Image
  using HistogramImageType = itk::Image<unsigned int, 2>;
  HistogramImageType::Pointer histogramImg = HistogramImageType::New();

  HistogramImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;

  HistogramImageType::SizeType size;
  size[0] = HISTOGRAMSIZE;
  size[1] = HISTOGRAMSIZE;

  HistogramImageType::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);

  histogramImg->SetRegions(region);

  HistogramImageType::SpacingType space;
  space[0] = 1;
  space[1] = 1;
  histogramImg->SetSpacing(space);

  histogramImg->Allocate();

  // * Iterator For Image

  using HistIteratorType = itk::ImageRegionIterator<HistogramImageType>;
  HistIteratorType hit(histogramImg, histogramImg->GetLargestPossibleRegion());
  for (hit.GoToBegin(); !hit.IsAtEnd(); ++hit)
  {
    HistogramImageType::IndexType currentIdx = hit.GetIndex();
    hit.Set(HistogramArray[currentIdx[0]][currentIdx[1]]);
  }

  // Histogram Image Rescale to 0-255

  using HistogramWritingType = itk::Image<unsigned char, 2>;
  using HistogramRescaleFilterType = itk::RescaleIntensityImageFilter<HistogramImageType, HistogramWritingType>;

  HistogramRescaleFilterType::Pointer histogramRescaler = HistogramRescaleFilterType::New();

  histogramRescaler->SetInput(histogramImg);

  histogramRescaler->SetOutputMaximum(255);
  histogramRescaler->SetOutputMinimum(0);

  // Invert Intensity so that background is white
  using InvertIntensityFilterType = itk::InvertIntensityImageFilter<HistogramWritingType, HistogramWritingType>;
  InvertIntensityFilterType::Pointer invertFilter = InvertIntensityFilterType::New();

  invertFilter->SetInput(histogramRescaler->GetOutput());
  invertFilter->SetMaximum(255);
  invertFilter->Update();

  //  Histogram Writer
  using HistogramWriter = itk::ImageFileWriter<HistogramWritingType>;

  HistogramWriter::Pointer histogramWriter = HistogramWriter::New();

  histogramWriter->SetFileName(outputJointHistogramImage);
  histogramWriter->SetInput(invertFilter->GetOutput());

  try
  {
    histogramWriter->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cerr << "Exception in Resampling." << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    exit(EXIT_FAILURE);
  }

  exit(EXIT_SUCCESS);
}
