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

#include "itkIO.h"
#include "itkMultiplyImageFilter.h"
#include "BRAINSMultiModeSegmentCLP.h"
#include "BRAINSThreadControl.h"
#include "itkMultiModeHistogramThresholdBinaryImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkBRAINSROIAutoImageFilter.h"
#include <BRAINSCommonLib.h>

int
main( int argc, char * argv[] )
{
  PARSE_ARGS;
  using ImageType = itk::Image< signed short, 3 >;
  using MaskImageType = itk::Image< unsigned char, 3 >;

  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder( numberOfThreads );
  if ( inputVolumes.size() < 1 )
  {
    std::cerr << argv[0] << ": Missing required --inputVolumes parameter" << std::endl;
    return EXIT_FAILURE;
  }
  const unsigned int numberOfModes = inputVolumes.size();
  using ThresholdRegionFinderType = itk::MultiModeHistogramThresholdBinaryImageFilter< ImageType, MaskImageType >;
  ThresholdRegionFinderType::ThresholdArrayType QuantileLowerThreshold( numberOfModes );
  ThresholdRegionFinderType::ThresholdArrayType QuantileUpperThreshold( numberOfModes );

  ThresholdRegionFinderType::Pointer thresholdRegionFinder = ThresholdRegionFinderType::New();

  MaskImageType::Pointer RegionMaskVolume = nullptr;
  if ( inputMaskVolume != "" )
  {
    RegionMaskVolume = itkUtil::ReadImage< MaskImageType >( inputMaskVolume );
  }
  for ( unsigned int modeIndex = 0; modeIndex < numberOfModes; modeIndex++ )
  {
    {
      ImageType::Pointer ImageInput = itkUtil::ReadImage< ImageType >( inputVolumes[modeIndex] );
      thresholdRegionFinder->SetInput( modeIndex, ImageInput );
      if ( RegionMaskVolume.IsNull() ) // USE ROIAUTO if no explicit mask is
                                       // specified on command line.
      {
        using ROIAutoType = itk::BRAINSROIAutoImageFilter< ImageType, MaskImageType >;
        ROIAutoType::Pointer ROIFilter = ROIAutoType::New();
        ROIFilter->SetInput( ImageInput );
        ROIFilter->Update();
        RegionMaskVolume = ROIFilter->GetOutput();
      }
      thresholdRegionFinder->SetBinaryPortionImage( RegionMaskVolume );
    }
    {
      const float lower = lowerThreshold[modeIndex];
      const float upper = upperThreshold[modeIndex];
      QuantileLowerThreshold.SetElement( modeIndex, lower );
      QuantileUpperThreshold.SetElement( modeIndex, upper );
    }
  }
  // std::cout << "============ Starting Threshold for Prior index: " << i <<
  // std::endl;
  // Assume upto (2*0.025)% of intensities are noise that corrupts the image
  // min/max values
  thresholdRegionFinder->SetLinearQuantileThreshold( 0.025 );
  thresholdRegionFinder->SetQuantileLowerThreshold( QuantileLowerThreshold );
  thresholdRegionFinder->SetQuantileUpperThreshold( QuantileUpperThreshold );
  // thresholdRegionFinder->SetInsideValue(1);
  // thresholdRegionFinder->SetOutsideValue(0);//Greatly reduce the value to
  // zero.
  thresholdRegionFinder->Update();
  // INFO:  Investigate if a small gaussian filter is needed here after
  // clipping.
  MaskImageType::Pointer MaskImage = thresholdRegionFinder->GetOutput();

  if ( outputROIMaskVolume != "" )
  {
    itkUtil::WriteImage< MaskImageType >( MaskImage, outputROIMaskVolume );
  }

  if ( outputClippedVolumeROI != "" )
  {
    std::cout << "WARNING:  This feature is not yet implemented! " << std::endl;
    //    itkUtil::WriteImage<unsigned char>(,MaskImage,outputClippedVolumeROI);
  }

  return 0;
}
