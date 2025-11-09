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
 *
 *  Program:   Insight Segmentation & Registration Toolkit
 *  Module:    $RCSfile$
 *  Language:  C++
 *  Date:      $Date: 2007-08-31 11:20:20 -0500 (Fri, 31 Aug 2007) $
 *  Version:   $Revision: 10358 $
 *
 *  Copyright (c) Insight Software Consortium. All rights reserved.
 *  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even
 *  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the above copyright notices for more information.
 *
 *  =========================================================================*/

/**
 * Hans J. Johnson @ The University of Iowa
 * This program is a standalone version of a program for masking and clipping
 *images
 * using the ROIAUTO method that seems to work well for brain images.
 */

#include <itkMultiplyImageFilter.h>
#include <itkImageMaskSpatialObject.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkExtractImageFilter.h>
#include <vcl_compiler.h>
#include <iostream>
#include <algorithm>

#include <BRAINSCommonLib.h>
#include "itkIO.h"
#include "itkLargestForegroundFilledMaskImageFilter.h"
#include "itkBRAINSROIAutoImageFilter.h"
#include "BRAINSROIAutoCLP.h"
#include "BRAINSThreadControl.h"
using VolumeImageType = itk::Image<signed int, 3>;
using VolumeMaskType = itk::Image<unsigned char, 3>;
using SOImageMaskType = itk::SpatialObject<3>;

/**
 * This file contains utility functions that are common to a few of the
 *BRAINSFit Programs.
 */

template <typename PixelType>
void
BRAINSROIAUTOWriteOutputVolume(const VolumeImageType::Pointer & image,
                               const VolumeMaskType::Pointer &  mask,
                               std::string &                    fileName,
                               const bool                       MaskImage,
                               const bool                       CropImage)
{
  using WriteOutImageType = typename itk::Image<PixelType, VolumeImageType::ImageDimension>;
  typename WriteOutImageType::Pointer finalOutput;
  {
    using CasterType = itk::CastImageFilter<VolumeImageType, WriteOutImageType>;
    const typename CasterType::Pointer myCaster = CasterType::New();
    myCaster->SetInput(image);
    myCaster->Update();
    finalOutput = myCaster->GetOutput();
  }
  if (MaskImage)
  {
    using MultiplierType = typename itk::MultiplyImageFilter<VolumeMaskType, WriteOutImageType, WriteOutImageType>;

    const typename MultiplierType::Pointer clipper = MultiplierType::New();
    clipper->SetInput1(mask);
    clipper->SetInput2(finalOutput);
    clipper->Update();
    finalOutput = clipper->GetOutput();
  }
  if (CropImage)
  {
    typename VolumeMaskType::IndexType minIndex;
    typename VolumeMaskType::IndexType maxIndex;
    for (VolumeMaskType::IndexType::IndexValueType i = 0; i < VolumeMaskType::ImageDimension; ++i)
    {
      minIndex[i] = std::numeric_limits<VolumeMaskType::IndexType::IndexValueType>::max();
      maxIndex[i] = std::numeric_limits<VolumeMaskType::IndexType::IndexValueType>::min();
    }
    itk::ImageRegionConstIteratorWithIndex<VolumeMaskType> maskIt(mask, mask->GetLargestPossibleRegion());
    while (!maskIt.IsAtEnd())
    {
      if (maskIt.Get() > 0)
      {
        const typename VolumeMaskType::IndexType & currIndex = maskIt.GetIndex();
        for (VolumeMaskType::IndexType::IndexValueType i = 0; i < VolumeMaskType::ImageDimension; ++i)
        {
          minIndex[i] = std::min(minIndex[i], currIndex[i]);
          maxIndex[i] = std::max(maxIndex[i], currIndex[i]);
        }
      }
      ++maskIt;
    }

    VolumeMaskType::SizeType desiredSize;
    for (VolumeMaskType::IndexType::IndexValueType i = 0; i < VolumeMaskType::ImageDimension; ++i)
    {
      desiredSize[i] = maxIndex[i] - minIndex[i] - 1;
    }
    // Now make the desiredSize an even number.
    for (VolumeMaskType::IndexType::IndexValueType i = 0; i < VolumeMaskType::ImageDimension; ++i)
    {
      desiredSize[i] += (desiredSize[i] % 2); // If even number, then add 0, if odd number, then add 1.
    }
    const VolumeMaskType::RegionType desiredRegion(minIndex, desiredSize);

    using ExtractorType = itk::ExtractImageFilter<WriteOutImageType, WriteOutImageType>;
    const typename ExtractorType::Pointer myExtractor = ExtractorType::New();
    myExtractor->SetExtractionRegion(desiredRegion);
    myExtractor->SetInput(finalOutput);
    myExtractor->SetDirectionCollapseToIdentity(); // This is required.
    myExtractor->Update();
    finalOutput = myExtractor->GetOutput();
  }
  itkUtil::WriteImage<WriteOutImageType>(finalOutput, fileName);
}

int
main(int argc, char * argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  if (inputVolume.empty())
  {
    std::cerr << argv[0] << ": Missing required --inputVolume parameter" << std::endl;
    return EXIT_FAILURE;
  }
  const VolumeImageType::Pointer ImageInput = itkUtil::ReadImage<VolumeImageType>(inputVolume);

  using ROIAutoType = itk::BRAINSROIAutoImageFilter<VolumeImageType, VolumeMaskType>;
  const ROIAutoType::Pointer ROIFilter = ROIAutoType::New();
  ROIFilter->SetInput(ImageInput);
  ROIFilter->SetOtsuPercentileThreshold(otsuPercentileThreshold);
  ROIFilter->SetClosingSize(closingSize);
  ROIFilter->SetThresholdCorrectionFactor(thresholdCorrectionFactor);
  ROIFilter->SetDilateSize(ROIAutoDilateSize);
  ROIFilter->Update();
  // const SOImageMaskType::Pointer maskWrapper = ROIFilter->GetSpatialObjectROI();
  const VolumeMaskType::Pointer MaskImage = ROIFilter->GetOutput();

  if (!outputROIMaskVolume.empty())
  {
    itkUtil::WriteImage<VolumeMaskType>(MaskImage, outputROIMaskVolume);
  }

  if (!outputVolume.empty())
  {
    //      std::cout << "=========== resampledImage :\n" <<
    // resampledImage->GetDirection() << std::endl;
    // Set in PARSEARGS const bool scaleOutputValues=false;//INFO: Make this a
    // command line parameter
    if (outputVolumePixelType == "float")
    {
      BRAINSROIAUTOWriteOutputVolume<float>(ImageInput, MaskImage, outputVolume, maskOutput, cropOutput);
    }
    else if (outputVolumePixelType == "short")
    {
      BRAINSROIAUTOWriteOutputVolume<signed short>(ImageInput, MaskImage, outputVolume, maskOutput, cropOutput);
    }
    else if (outputVolumePixelType == "ushort")
    {
      BRAINSROIAUTOWriteOutputVolume<unsigned short>(ImageInput, MaskImage, outputVolume, maskOutput, cropOutput);
    }
    else if (outputVolumePixelType == "int")
    {
      BRAINSROIAUTOWriteOutputVolume<signed int>(ImageInput, MaskImage, outputVolume, maskOutput, cropOutput);
    }
    else if (outputVolumePixelType == "uint")
    {
      BRAINSROIAUTOWriteOutputVolume<unsigned int>(ImageInput, MaskImage, outputVolume, maskOutput, cropOutput);
    }
    else if (outputVolumePixelType == "uchar")
    {
      BRAINSROIAUTOWriteOutputVolume<unsigned char>(ImageInput, MaskImage, outputVolume, maskOutput, cropOutput);
    }
  }
  return EXIT_SUCCESS;
}
