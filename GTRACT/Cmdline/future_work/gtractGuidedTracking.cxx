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
#include <itkImageFileReader.h>
#include <itkSpatialOrientationAdapter.h>
#include <itkThresholdImageFilter.h>

#include "gtractGuidedTrackingCLP.h"
#include "itkDtiGuidedTrackingFilter.h"

template <typename TImageType>
void AdaptOriginAndDirection( typename TImageType::Pointer image )
{
  typename TImageType::DirectionType imageDir = image->GetDirection();
  typename TImageType::PointType origin = image->GetOrigin();

  int dominantAxisRL = itk::Function::Max3(imageDir[0][0], imageDir[1][0], imageDir[2][0]);
  int signRL = itk::Function::Sign(imageDir[dominantAxisRL][0]);
  int dominantAxisAP = itk::Function::Max3(imageDir[0][1], imageDir[1][1], imageDir[2][1]);
  int signAP = itk::Function::Sign(imageDir[dominantAxisAP][1]);
  // int dominantAxisSI =
  // itk::Function::Max3(imageDir[0][2],imageDir[1][2],imageDir[2][2]);
  // int signSI = itk::Function::Sign(imageDir[dominantAxisSI][2]);

  /* This current  algorithm needs to be verified.
     I had previously though that it should be
     signRL == 1
     signAP == -1
     signSI == 1
     This appears to be incorrect with the NRRD file format
     at least. Visually this appears to work
     signRL == 1
     signAP == 1
     signSI == -1
  */
  typename TImageType::DirectionType DirectionToRAS;
  DirectionToRAS.SetIdentity();
  if( signRL == 1 )
    {
    DirectionToRAS[dominantAxisRL][dominantAxisRL] = -1.0;
    origin[dominantAxisRL] *= -1.0;
    }
  if( signAP == 1 )
    {
    DirectionToRAS[dominantAxisAP][dominantAxisAP] = -1.0;
    origin[dominantAxisAP] *= -1.0;
    }
  /* This component still seems to be a problem */
  /*
  if (signSI == -1)
    {
    DirectionToRAS[dominantAxisSI][dominantAxisSI] = -1.0;
    origin[dominantAxisSI] *= -1.0;
    }
  */
  imageDir *= DirectionToRAS;
  image->SetDirection( imageDir  );
  image->SetOrigin( origin );
}

int main(int argc, char * *argv)
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  const bool debug = true;
  if( debug )
    {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Tensor Image: " <<  inputTensorVolume << std::endl;
    std::cout << "Anisotropy Image: " <<  inputAnisotropyVolume << std::endl;
    std::cout << "Starting Seeds LabelMap Image: " <<  inputStartingSeedsLabelMapVolume << std::endl;
    std::cout << "Ending Seeds LabelMap Image: " <<  inputEndingSeedsLabelMapVolume << std::endl;
    std::cout << "Input Guide Tract: " <<  inputTract << std::endl;
    std::cout << "Output Tract: " <<  outputTract << std::endl;
    std::cout << "Use XML PolyData Format: " << writeXMLPolyDataFile  << std::endl;
    std::cout << "Seed Threshold: " <<  seedThreshold << std::endl;
    std::cout << "Tracking Threshold: " <<  trackingThreshold << std::endl;
    std::cout << "Curvature Threshold: " <<  curvatureThreshold << std::endl;
    std::cout << "Guided Curvature Threshold: " <<  guidedCurvatureThreshold << std::endl;
    std::cout << "Minimum Length: " <<  minimumLength << std::endl;
    std::cout << "Maximum Length: " <<  maximumLength << std::endl;
    std::cout << "Step Size: " <<  stepSize << std::endl;
    std::cout << "Use Loop Detection: " <<  useLoopDetection << std::endl;
    std::cout << "Use Tensor Deflection: " <<  useTend << std::endl;
    std::cout << "Tend F: " <<  tendF << std::endl;
    std::cout << "Tend G: " <<  tendG << std::endl;
    std::cout << "Starting Label: " <<  startingSeedsLabel << std::endl;
    std::cout << "Ending Label: " <<  endingSeedsLabel << std::endl;
    std::cout << "Guide Distance: " <<  maximumGuideDistance << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  using TensorElementType = double;
  using TensorPixelType = itk::DiffusionTensor3D<TensorElementType>;
  using TensorImageType = itk::Image<TensorPixelType, 3>;
  using TensorImageReaderType = itk::ImageFileReader<TensorImageType>;

  TensorImageReaderType::Pointer tensorImageReader = TensorImageReaderType::New();
  tensorImageReader->SetFileName( inputTensorVolume );

  try
    {
    tensorImageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }

  TensorImageType::Pointer tensorImage = tensorImageReader->GetOutput();
  AdaptOriginAndDirection<TensorImageType>( tensorImage );

  using AnisotropyPixelType = float;
  using AnisotropyImageType = itk::Image<AnisotropyPixelType, 3>;
  using AnisotropyImageReaderType = itk::ImageFileReader<AnisotropyImageType>;
  AnisotropyImageReaderType::Pointer anisotropyImageReader = AnisotropyImageReaderType::New();
  anisotropyImageReader->SetFileName( inputAnisotropyVolume );

  try
    {
    anisotropyImageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }

  AnisotropyImageType::Pointer anisotropyImage = anisotropyImageReader->GetOutput();
  AdaptOriginAndDirection<AnisotropyImageType>( anisotropyImage );

  using MaskPixelType = unsigned char;
  using MaskImageType = itk::Image<MaskPixelType, 3>;
  using MaskImageReaderType = itk::ImageFileReader<MaskImageType>;
  MaskImageReaderType::Pointer startingSeedImageReader = MaskImageReaderType::New();
  startingSeedImageReader->SetFileName( inputStartingSeedsLabelMapVolume );

  try
    {
    startingSeedImageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }

  /* Threshold Starting Label Map */
  using ThresholdFilterType = itk::ThresholdImageFilter<MaskImageType>;
  ThresholdFilterType::Pointer startingThresholdFilter = ThresholdFilterType::New();
  startingThresholdFilter->SetInput( startingSeedImageReader->GetOutput() );
  startingThresholdFilter->SetLower( static_cast<MaskPixelType>( startingSeedsLabel ) );
  startingThresholdFilter->SetUpper( static_cast<MaskPixelType>( startingSeedsLabel ) );
  startingThresholdFilter->Update();

  MaskImageType::Pointer startingSeedMask = startingThresholdFilter->GetOutput();
  AdaptOriginAndDirection<MaskImageType>( startingSeedMask );

  MaskImageReaderType::Pointer endingSeedImageReader = MaskImageReaderType::New();
  endingSeedImageReader->SetFileName( inputEndingSeedsLabelMapVolume );

  try
    {
    endingSeedImageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }

  /* Threshold Ending Label Map */
  ThresholdFilterType::Pointer endingThresholdFilter = ThresholdFilterType::New();
  endingThresholdFilter->SetInput( endingSeedImageReader->GetOutput() );
  endingThresholdFilter->SetLower( static_cast<MaskPixelType>( endingSeedsLabel ) );
  endingThresholdFilter->SetUpper( static_cast<MaskPixelType>( endingSeedsLabel ) );
  endingThresholdFilter->Update();

  MaskImageType::Pointer endingSeedMask = endingThresholdFilter->GetOutput();
  AdaptOriginAndDirection<MaskImageType>( endingSeedMask );

  std::cerr << "Guided Tracking" << std::endl;

  vtkPolyData *guideFiber;
  if( writeXMLPolyDataFile )
    {
    vtkXMLPolyDataReader *guideFiberReader = vtkXMLPolyDataReader::New();
    guideFiberReader->SetFileName( inputTract.c_str() );
    guideFiberReader->Update();
    guideFiber = guideFiberReader->GetOutput();
    }
  else
    {
    vtkPolyDataReader *guideFiberReader = vtkPolyDataReader::New();
    guideFiberReader->SetFileName( inputTract.c_str() );
    guideFiberReader->Update();
    guideFiber = guideFiberReader->GetOutput();
    }

  using GuideTrackingFilterType = itk::DtiGuidedTrackingFilter<TensorImageType, AnisotropyImageType, MaskImageType>;
  GuideTrackingFilterType::Pointer guideTrackingFilter = GuideTrackingFilterType::New();
  guideTrackingFilter->SetAnisotropyImage( anisotropyImage );
  guideTrackingFilter->SetTensorImage( tensorImage );
  guideTrackingFilter->SetStartingRegion( startingSeedMask );
  guideTrackingFilter->SetEndingRegion( endingSeedMask );
  guideTrackingFilter->SetGuideFiber( guideFiber );
  guideTrackingFilter->SetSeedThreshold( seedThreshold );
  guideTrackingFilter->SetAnisotropyThreshold( trackingThreshold );
  guideTrackingFilter->SetCurvatureThreshold( curvatureThreshold );
  guideTrackingFilter->SetGuidedCurvatureThreshold( guidedCurvatureThreshold );
  guideTrackingFilter->SetMaximumGuideDistance( static_cast<double>( maximumGuideDistance ) );
  // Fix this once support for multiple Region tracking is added
  guideTrackingFilter->SetMaximumLength( maximumLength );
  guideTrackingFilter->SetMinimumLength( minimumLength );
  guideTrackingFilter->SetStepSize( stepSize );
  guideTrackingFilter->SetUseLoopDetection( useLoopDetection );
  guideTrackingFilter->SetUseTend( useTend );
  guideTrackingFilter->SetTendG( tendG );
  guideTrackingFilter->SetTendF( tendF );
  // std::cout << guideTrackingFilter;
  guideTrackingFilter->Update();

  vtkPolyData *fibers = guideTrackingFilter->GetOutput();

  if( writeXMLPolyDataFile )
    {
    vtkXMLPolyDataWriter *fiberWriter = vtkXMLPolyDataWriter::New();
    fiberWriter->SetFileName( outputTract.c_str() );
    fiberWriter->SetInput( fibers );
    fiberWriter->Update();
    }
  else
    {
    vtkPolyDataWriter *fiberWriter = vtkPolyDataWriter::New();
    fiberWriter->SetFileName( outputTract.c_str() );
    fiberWriter->SetInput( fibers );
    fiberWriter->Update();
    }
}
