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
#include <map>
#include <string>

#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkMetaDataObject.h>
#include <itkIOCommon.h>
#include <itkThresholdImageFilter.h>
#include <itkSpatialOrientationAdapter.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyDataWriter.h>

#include "itkDtiFastMarchingTrackingFilter.h"
#include "GtractTypes.h"
#include "gtractFastMarchingTrackingCLP.h"
#include "BRAINSThreadControl.h"
#include <BRAINSCommonLib.h>

template < typename TImageType >
void
AdaptOriginAndDirection( typename TImageType::Pointer image )
{
  typename TImageType::DirectionType imageDir = image->GetDirection();
  typename TImageType::PointType     origin = image->GetOrigin();

  int dominantAxisRL = itk::Function::Max3( imageDir[0][0], imageDir[1][0], imageDir[2][0] );
  int signRL = itk::Function::Sign( imageDir[dominantAxisRL][0] );
  int dominantAxisAP = itk::Function::Max3( imageDir[0][1], imageDir[1][1], imageDir[2][1] );
  int signAP = itk::Function::Sign( imageDir[dominantAxisAP][1] );
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
  if ( signRL == 1 )
  {
    DirectionToRAS[dominantAxisRL][dominantAxisRL] = -1.0;
    origin[dominantAxisRL] *= -1.0;
  }
  if ( signAP == 1 )
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
  image->SetDirection( imageDir );
  image->SetOrigin( origin );
}

int
main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder( numberOfThreads );

  const bool debug = true;
  if ( debug )
  {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Cost Image: " << inputCostVolume << std::endl;
    std::cout << "Anisotropy Image: " << inputAnisotropyVolume << std::endl;
    std::cout << "Tensor Image: " << inputTensorVolume << std::endl;
    std::cout << "Output Fiber(s): " << outputTract << std::endl;
    std::cout << "Use XML Polydata Format: " << writeXMLPolyDataFile << std::endl;
    std::cout << "Input Seed LabelMap Image: " << inputStartingSeedsLabelMapVolume << std::endl;
    std::cout << "Input Seed Label: " << startingSeedsLabel << std::endl;
    std::cout << "Number of Iterations: " << numberOfIterations << std::endl;
    std::cout << "Tracking Threshold: " << trackingThreshold << std::endl;
    std::cout << "Start Point Anisotropy Threshold: " << seedThreshold << std::endl;
    std::cout << "Subvoxel Cost Step Size: " << costStepSize << std::endl;
    std::cout << "Maximum Step Size: " << maximumStepSize << std::endl;
    std::cout << "Minimum Step Size: " << minimumStepSize << std::endl;
    std::cout << "=====================================================" << std::endl;
  }

  /* Read Tensor Image */
  using TensorElementType = double;
  using TensorPixelType = itk::DiffusionTensor3D< TensorElementType >;
  using TensorImageType = itk::Image< TensorPixelType, 3 >;
  using TensorImageReaderType = itk::ImageFileReader< TensorImageType >;
  TensorImageReaderType::Pointer tensorImageReader = TensorImageReaderType::New();
  tensorImageReader->SetFileName( inputTensorVolume );

  try
  {
    tensorImageReader->Update();
  }
  catch ( itk::ExceptionObject & ex )
  {
    std::cout << ex << std::endl;
    throw;
  }

  TensorImageType::Pointer tensorImage = tensorImageReader->GetOutput();
  AdaptOriginAndDirection< TensorImageType >( tensorImage );

  /* Read Cost Image */
  using PixelType = float;
  using CostImageType = itk::Image< PixelType, 3 >;
  using CostImageReaderType = itk::ImageFileReader< CostImageType >;
  CostImageReaderType::Pointer costImageReader = CostImageReaderType::New();
  costImageReader->SetFileName( inputCostVolume );

  try
  {
    costImageReader->Update();
  }
  catch ( itk::ExceptionObject & ex )
  {
    std::cout << ex << std::endl;
    throw;
  }

  CostImageType::Pointer costImage = costImageReader->GetOutput();
  AdaptOriginAndDirection< CostImageType >( costImage );

  /* Read Anisotropy Image */
  using AnisotropyImageType = itk::Image< PixelType, 3 >;
  using AnisotropyImageReaderType = itk::ImageFileReader< AnisotropyImageType >;
  AnisotropyImageReaderType::Pointer anisotropyImageReader = AnisotropyImageReaderType::New();
  anisotropyImageReader->SetFileName( inputAnisotropyVolume );

  try
  {
    anisotropyImageReader->Update();
  }
  catch ( itk::ExceptionObject & ex )
  {
    std::cout << ex << std::endl;
    throw;
  }

  AnisotropyImageType::Pointer anisotropyImage = anisotropyImageReader->GetOutput();
  AdaptOriginAndDirection< AnisotropyImageType >( anisotropyImage );

  /* Read the Mask Seed Region */
  using MaskPixelType = signed short;
  using MaskImageType = itk::Image< MaskPixelType, 3 >;
  using MaskImageReaderType = itk::ImageFileReader< MaskImageType >;
  MaskImageReaderType::Pointer startingSeedImageReader = MaskImageReaderType::New();
  startingSeedImageReader->SetFileName( inputStartingSeedsLabelMapVolume );

  try
  {
    startingSeedImageReader->Update();
  }
  catch ( itk::ExceptionObject & ex )
  {
    std::cout << ex << std::endl;
    throw;
  }

  /* Threshold Starting Label Map */
  using ThresholdFilterType = itk::ThresholdImageFilter< MaskImageType >;
  ThresholdFilterType::Pointer startingThresholdFilter = ThresholdFilterType::New();
  startingThresholdFilter->SetInput( startingSeedImageReader->GetOutput() );
  startingThresholdFilter->SetLower( static_cast< MaskPixelType >( startingSeedsLabel ) );
  startingThresholdFilter->SetUpper( static_cast< MaskPixelType >( startingSeedsLabel ) );
  startingThresholdFilter->Update();

  MaskImageType::Pointer startingSeedMask = startingThresholdFilter->GetOutput();
  AdaptOriginAndDirection< MaskImageType >( startingSeedMask );

  /*Set the Parameters and Run the DtiFastMarchingTrackingFilter */
  using TrackingFilterType =
    itk::DtiFastMarchingTrackingFilter< TensorImageType, AnisotropyImageType, CostImageType, MaskImageType >;
  TrackingFilterType::Pointer trackFilter = TrackingFilterType::New();

  trackFilter->SetCostImage( costImage );
  trackFilter->SetAnisotropyImage( anisotropyImage );
  trackFilter->SetTensorImage( tensorImage );
  trackFilter->SetStartingRegion( startingSeedMask );
  trackFilter->SetNumberOfIterations( numberOfIterations );
  trackFilter->SetMaxStepSize( maximumStepSize );
  trackFilter->SetMinStepSize( minimumStepSize );
  trackFilter->SetCostFunctionStepSize( costStepSize );
  trackFilter->SetSeedThreshold( seedThreshold );
  trackFilter->SetAnisotropyThreshold( trackingThreshold );
  trackFilter->Update();

  vtkPolyData * fibers = trackFilter->GetOutput();

  if ( writeXMLPolyDataFile )
  {
    vtkXMLPolyDataWriter * fiberWriter = vtkXMLPolyDataWriter::New();
    fiberWriter->SetFileName( outputTract.c_str() );
#if ( VTK_MAJOR_VERSION < 6 )
    fiberWriter->SetInput( fibers );
#else
    fiberWriter->SetInputData( fibers );
#endif
    fiberWriter->Update();
  }
  else
  {
    vtkPolyDataWriter * fiberWriter = vtkPolyDataWriter::New();
    fiberWriter->SetFileName( outputTract.c_str() );
#if ( VTK_MAJOR_VERSION < 6 )
    fiberWriter->SetInput( fibers );
#else
    fiberWriter->SetInputData( fibers );
#endif
    fiberWriter->Update();
  }
  return EXIT_SUCCESS;
}
