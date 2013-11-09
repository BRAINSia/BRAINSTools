/*=========================================================================

 Program:   GTRACT (Guided Tensor Restore Anatomical Connectivity Tractography)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2010/05/03 09:23:40 $
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

#include <vtkMatrix4x4.h>
#include <vtkMatrixToLinearTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include "itkDtiGuidedTrackingFilter.h"
#include "itkDtiFreeTrackingFilter.h"
#include "itkDtiGraphSearchTrackingFilter.h"
#include "itkDtiStreamlineTrackingFilter.h"

#include "gtractFiberTrackingCLP.h"
#include "BRAINSThreadControl.h"
#include <BRAINSCommonLib.h>

template <class TImageType>
void AdaptOriginAndDirection( typename TImageType::Pointer image )
{
  typename TImageType::DirectionType imageDir = image->GetDirection();
  typename TImageType::PointType origin = image->GetOrigin();
  typename TImageType::SpacingType spacing = image->GetSpacing();

  imageDir.Fill(0);
  imageDir[0][0] = 1.0;
  imageDir[1][1] = 1.0;
  imageDir[2][2] = 1.0;

  origin.Fill(0);
  spacing.Fill(1.0);

  image->SetDirection( imageDir  );
  image->SetOrigin( origin );
  image->SetSpacing( spacing );
}

template <class TImageType>
vtkMatrix4x4 * CreateIjkToRasMatrix( typename TImageType::Pointer image )
{
  double        spacing[3];
  double        origin[3];
  vtkMatrix4x4* IjkToLpsMatrix = vtkMatrix4x4::New();
  vtkMatrix4x4* RasToIjkMatrix = vtkMatrix4x4::New();

  IjkToLpsMatrix->Identity();
  for( unsigned int i = 0; i < 3; i++ )
    {
    spacing[i] = image->GetSpacing()[i];
    origin[i] = image->GetOrigin()[i];
    // Get IJK to LPS direction vector
    for( unsigned int j = 0; j < image->GetImageDimension(); j++ )
      {
      IjkToLpsMatrix->SetElement(j, i, spacing[i] * image->GetDirection()[j][i]);
      }
    }

  vtkMatrix4x4* LpsToRasMatrix = vtkMatrix4x4::New();
  LpsToRasMatrix->Identity();
  LpsToRasMatrix->SetElement(0, 0, -1);
  LpsToRasMatrix->SetElement(1, 1, -1);

  vtkMatrix4x4::Multiply4x4(LpsToRasMatrix, IjkToLpsMatrix, RasToIjkMatrix);

  origin[0] *= -1;   // L -> R
  origin[1] *= -1;   // P -> A
  for( unsigned int j = 0; j < 3; j++ )
    {
    RasToIjkMatrix->SetElement(j, 3, origin[j]);
    }
  RasToIjkMatrix->SetElement(3, 3, 1.0);

  return RasToIjkMatrix;
}

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  const bool debug = true;
  if( debug )
    {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Tensor Image: " <<  inputTensorVolume << std::endl;
    std::cout << "Anisotropy Image: " <<  inputAnisotropyVolume << std::endl;
    std::cout << "Output Tract: " <<  outputTract << std::endl;
    std::cout << "Starting Seeds LabelMap Image: " <<  inputStartingSeedsLabelMapVolume << std::endl;
    std::cout << "Ending Seeds LabelMap Image: " <<  inputEndingSeedsLabelMapVolume << std::endl;
    std::cout << "Input Guide Tract: " <<  inputTract << std::endl;
    std::cout << "Use XML PolyData Format: " << writeXMLPolyDataFile  << std::endl;
    std::cout << "Seed Threshold: " <<  seedThreshold << std::endl;
    std::cout << "Tracking Threshold: " <<  trackingThreshold << std::endl;
    std::cout << "Curvature Threshold: " <<  curvatureThreshold << std::endl;
    std::cout << "Guided Curvature Threshold: " <<  guidedCurvatureThreshold << std::endl;
    std::cout << "Branching Threshold: " <<  branchingThreshold << std::endl;
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
    std::cout << "Maximum Branch Points: " <<  maximumBranchPoints << std::endl;
    std::cout << "Use Random Walk: " <<  useRandomWalk << std::endl;
    std::cout << "Random Seed: " <<  randomSeed << std::endl;
    std::cout << "Branching Angle: " <<  branchingAngle << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  typedef double                                    TensorElementType;
  typedef itk::DiffusionTensor3D<TensorElementType> TensorPixelType;
  typedef itk::Image<TensorPixelType, 3>            TensorImageType;
  typedef itk::ImageFileReader<TensorImageType>     TensorImageReaderType;
  if( inputTensorVolume == "" )
    {
    std::cerr << "Missing Filename for input Tensor Volume (--inputTensorVolume)" << std::endl;
    return EXIT_FAILURE;
    }
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
  vtkMatrix4x4*            IjkToRasMatrix = CreateIjkToRasMatrix<TensorImageType>( tensorImage );
  AdaptOriginAndDirection<TensorImageType>( tensorImage );

  // std::cout <<  "Tensor Image : " << tensorImage << std::endl;
  if( inputAnisotropyVolume == "" )
    {
    std::cerr << "Missing filename for input Anisotropy Volume (--inputAnisotropyVolume)"
              << std::endl;
    return EXIT_FAILURE;
    }
  typedef float                                     AnisotropyPixelType;
  typedef itk::Image<AnisotropyPixelType, 3>        AnisotropyImageType;
  typedef itk::ImageFileReader<AnisotropyImageType> AnisotropyImageReaderType;
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
  // std::cout << "Anisotropy Image: " << anisotropyImage << std::endl;
  AdaptOriginAndDirection<AnisotropyImageType>( anisotropyImage );
  // std::cout << "Anisotropy Image Updated: " << anisotropyImage << std::endl;

  if( inputStartingSeedsLabelMapVolume == "" )
    {
    std::cerr << "Missing filename for input Starting Seeds Label Map Volume (--inputStartingSeedsLabelMapVolume)"
              << std::endl;
    }
  typedef unsigned char                       MaskPixelType;
  typedef itk::Image<MaskPixelType, 3>        MaskImageType;
  typedef itk::ImageFileReader<MaskImageType> MaskImageReaderType;
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
  typedef itk::ThresholdImageFilter<MaskImageType> ThresholdFilterType;
  ThresholdFilterType::Pointer startingThresholdFilter = ThresholdFilterType::New();
  startingThresholdFilter->SetInput( startingSeedImageReader->GetOutput() );
  startingThresholdFilter->SetLower( static_cast<MaskPixelType>( startingSeedsLabel ) );
  startingThresholdFilter->SetUpper( static_cast<MaskPixelType>( startingSeedsLabel ) );
  startingThresholdFilter->Update();

  MaskImageType::Pointer startingSeedMask = startingThresholdFilter->GetOutput();
  AdaptOriginAndDirection<MaskImageType>( startingSeedMask );
  // std::cout <<  "Seed Mask : " << startingSeedMask << std::endl;

  MaskImageType::Pointer endingSeedMask;

  if( trackingMethod != "Free" )
    {
    if( inputEndingSeedsLabelMapVolume == "" )
      {
      std::cerr << "Missing filename for input Ending Seeds Label Map (--inputEndingSeedsLabelMapVolume)"
                << std::endl;
      return EXIT_FAILURE;
      }
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

    endingSeedMask = endingThresholdFilter->GetOutput();
    AdaptOriginAndDirection<MaskImageType>( endingSeedMask );
    }

  vtkPolyData *fibers;
  if( trackingMethod == "Guided" )
    {
    vtkPolyData *guideFiber;
    if( inputTract == "" )
      {
      std::cerr << "Missing Input Guide Tract (--inputTract)" << std::endl;
      return EXIT_FAILURE;
      }
    else
      {
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
      }
    /* Put the guide fiber into IJK space for Tracking */
    vtkMatrix4x4* RasToIjkMatrix = vtkMatrix4x4::New();
    RasToIjkMatrix->DeepCopy(IjkToRasMatrix);
    RasToIjkMatrix->Invert();
    vtkMatrixToLinearTransform *rasToijkTransform = vtkMatrixToLinearTransform::New();
    rasToijkTransform->SetInput(RasToIjkMatrix);

    vtkTransformPolyDataFilter *transformGuideFiber = vtkTransformPolyDataFilter::New();
    transformGuideFiber->SetTransform(rasToijkTransform);
    transformGuideFiber->SetInput( guideFiber );
    transformGuideFiber->Update();

    typedef itk::DtiGuidedTrackingFilter<TensorImageType, AnisotropyImageType, MaskImageType> GuideTrackingFilterType;
    GuideTrackingFilterType::Pointer acturalTrackingFilter = GuideTrackingFilterType::New();
    acturalTrackingFilter->SetEndingRegion( endingSeedMask );
    acturalTrackingFilter->SetGuideFiber( transformGuideFiber->GetOutput() );
    acturalTrackingFilter->SetCurvatureThreshold( curvatureThreshold );
    acturalTrackingFilter->SetGuidedCurvatureThreshold( guidedCurvatureThreshold );
    acturalTrackingFilter->SetMaximumGuideDistance( static_cast<double>( maximumGuideDistance ) );
    // Fix this once support for multiple Region tracking is added
    acturalTrackingFilter->SetAnisotropyImage( anisotropyImage );
    acturalTrackingFilter->SetTensorImage( tensorImage );
    acturalTrackingFilter->SetStartingRegion( startingSeedMask );
    acturalTrackingFilter->SetMaximumLength( maximumLength );
    acturalTrackingFilter->SetMinimumLength( minimumLength );
    acturalTrackingFilter->SetStepSize( stepSize );
    acturalTrackingFilter->SetTendG( tendG );
    acturalTrackingFilter->SetTendF( tendF );
    acturalTrackingFilter->SetUseTend( useTend );
    acturalTrackingFilter->SetUseLoopDetection( useLoopDetection );
    acturalTrackingFilter->SetSeedThreshold( seedThreshold );
    acturalTrackingFilter->SetAnisotropyThreshold( trackingThreshold );
    acturalTrackingFilter->Update();
    fibers = acturalTrackingFilter->GetOutput();
    }
  else if( trackingMethod == "Streamline" )
    {
    typedef itk::DtiStreamlineTrackingFilter<TensorImageType, AnisotropyImageType,
                                             MaskImageType> StreamTrackingFilterType;
    StreamTrackingFilterType::Pointer acturalTrackingFilter = StreamTrackingFilterType::New();
    acturalTrackingFilter->SetEndingRegion( endingSeedMask );
    acturalTrackingFilter->SetCurvatureThreshold( curvatureThreshold );
    acturalTrackingFilter->SetAnisotropyImage( anisotropyImage );
    acturalTrackingFilter->SetTensorImage( tensorImage );
    acturalTrackingFilter->SetStartingRegion( startingSeedMask );
    acturalTrackingFilter->SetMaximumLength( maximumLength );
    acturalTrackingFilter->SetMinimumLength( minimumLength );
    acturalTrackingFilter->SetStepSize( stepSize );
    acturalTrackingFilter->SetTendG( tendG );
    acturalTrackingFilter->SetTendF( tendF );
    acturalTrackingFilter->SetUseTend( useTend );
    acturalTrackingFilter->SetUseLoopDetection( useLoopDetection );
    acturalTrackingFilter->SetSeedThreshold( seedThreshold );
    acturalTrackingFilter->SetAnisotropyThreshold( trackingThreshold );
    acturalTrackingFilter->Update();
    fibers = acturalTrackingFilter->GetOutput();
    }
  else if( trackingMethod == "Free" )
    {
    typedef  itk::DtiFreeTrackingFilter<TensorImageType, AnisotropyImageType, MaskImageType> TrackingType;
    TrackingType::Pointer acturalTrackingFilter = TrackingType::New();
    acturalTrackingFilter->SetCurvatureThreshold( curvatureThreshold ); /*
                                                                          Convert
                                                                          to cos
                                                                          (curvature*pi/180)

                                                                          within
                                                                          method
                                                                          */
    acturalTrackingFilter->SetAnisotropyImage( anisotropyImage );
    acturalTrackingFilter->SetTensorImage( tensorImage );
    acturalTrackingFilter->SetStartingRegion( startingSeedMask );
    acturalTrackingFilter->SetMaximumLength( maximumLength );
    acturalTrackingFilter->SetMinimumLength( minimumLength );
    acturalTrackingFilter->SetStepSize( stepSize );
    acturalTrackingFilter->SetTendG( tendG );
    acturalTrackingFilter->SetTendF( tendF );
    acturalTrackingFilter->SetUseTend( useTend );
    acturalTrackingFilter->SetUseLoopDetection( useLoopDetection );
    acturalTrackingFilter->SetSeedThreshold( seedThreshold );
    acturalTrackingFilter->SetAnisotropyThreshold( trackingThreshold );
    acturalTrackingFilter->Update();
    fibers = acturalTrackingFilter->GetOutput();
    }
  else if( trackingMethod == "GraphSearch" )
    {
    typedef itk::DtiGraphSearchTrackingFilter<TensorImageType, AnisotropyImageType,
                                              MaskImageType> GraphTrackingFilterType;
    GraphTrackingFilterType::Pointer acturalTrackingFilter = GraphTrackingFilterType::New();
    acturalTrackingFilter->SetEndingRegion( endingSeedMask  );
    acturalTrackingFilter->SetAnisotropyBranchingValue( branchingThreshold );
    acturalTrackingFilter->SetCurvatureBranchAngle( curvatureThreshold );
    acturalTrackingFilter->SetMaximumBranches( maximumBranchPoints );
    acturalTrackingFilter->SetUseRandomWalk( useRandomWalk );
    acturalTrackingFilter->SetRandomWalkAngle( branchingAngle );
    acturalTrackingFilter->SetRandomSeed( randomSeed );
    acturalTrackingFilter->SetAnisotropyImage( anisotropyImage );
    acturalTrackingFilter->SetTensorImage( tensorImage );
    acturalTrackingFilter->SetStartingRegion( startingSeedMask );
    acturalTrackingFilter->SetMaximumLength( maximumLength );
    acturalTrackingFilter->SetMinimumLength( minimumLength );
    acturalTrackingFilter->SetStepSize( stepSize );
    acturalTrackingFilter->SetTendG( tendG );
    acturalTrackingFilter->SetTendF( tendF );
    acturalTrackingFilter->SetUseTend( useTend );
    acturalTrackingFilter->SetUseLoopDetection( useLoopDetection );
    acturalTrackingFilter->SetSeedThreshold( seedThreshold );
    acturalTrackingFilter->SetAnisotropyThreshold( trackingThreshold );
    acturalTrackingFilter->Update();
    fibers = acturalTrackingFilter->GetOutput();
    }
  else
    {
    std::cout << "No correct tracking method!" << std::endl;
    return EXIT_FAILURE;
    }

  //   vtkPolyData *fibers = acturalTrackingFilter->GetOutput();

  vtkMatrixToLinearTransform *ijkToRasTransform = vtkMatrixToLinearTransform::New();
  ijkToRasTransform->SetInput(IjkToRasMatrix);

  vtkTransformPolyDataFilter *transformPolyData = vtkTransformPolyDataFilter::New();
  transformPolyData->SetTransform(ijkToRasTransform);
  transformPolyData->SetInput( fibers );
  transformPolyData->Update();

  if( writeXMLPolyDataFile )
    {
    vtkXMLPolyDataWriter *fiberWriter = vtkXMLPolyDataWriter::New();
    fiberWriter->SetFileName( outputTract.c_str() );
    fiberWriter->SetInput( transformPolyData->GetOutput() );
    fiberWriter->Update();
    }
  else
    {
    vtkPolyDataWriter *fiberWriter = vtkPolyDataWriter::New();
    fiberWriter->SetFileName( outputTract.c_str() );
    fiberWriter->SetInput( transformPolyData->GetOutput() );
    fiberWriter->Update();
    }
  return EXIT_SUCCESS;
}
