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
#include <itkBSplineDeformableTransform.h>
#include <itkTransformFactory.h>
#include <itkVersorRigid3DTransform.h>
#include <itkOrientImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include "itkGtractImageIO.h"

#include "gtractResampleCodeImageCLP.h"
#include "itkThinPlateR2LogRSplineKernelTransform.h"
#include "GenericTransformImage.h"
#include "BRAINSThreadControl.h"
int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  typedef double                              BSplineCoordinateRepType;
  typedef itk::VersorRigid3DTransform<double> RigidTransformType;
  typedef itk::BSplineDeformableTransform<BSplineCoordinateRepType, 3, 3>
    BSplineTransformType;

  typedef itk::ThinPlateR2LogRSplineKernelTransform<BSplineCoordinateRepType, 3>
    ThinPlateSplineTransformType;

  bool debug = true;
  if( debug )
    {
    std::cout << "==============================================================" << std::endl;
    std::cout << "Input Image: " <<  inputCodeVolume << std::endl;
    std::cout << "Output Image: " <<  outputVolume << std::endl;
    std::cout << "Reference File: " <<  inputReferenceVolume << std::endl;
    std::cout << "Transform File: " << inputTransform << std::endl;
    std::cout << "Transform Type: " << transformType << std::endl;
    std::cout << "==============================================================" << std::endl;
    }

  bool violated = false;
  if( inputCodeVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputCodeVolume Required! "  << std::endl;
    }
  if( inputReferenceVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputReferenceVolume Required! "
                               << std::endl;
    }
  if( inputTransform.size() == 0 )
    {
    violated = true; std::cout << "  --inputTransform Required! "  << std::endl;
    }
  if( outputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --outputVolume Required! "  << std::endl;
    }
  if( violated )
    {
    return EXIT_FAILURE;
    }

  typedef signed short CodePixelType;

  typedef itk::Image<CodePixelType, 3>        CodeImageType;
  typedef itk::ImageFileReader<CodeImageType> CodeImageReaderType;
  CodeImageReaderType::Pointer codeImageReader = CodeImageReaderType::New();
  codeImageReader->SetFileName( inputCodeVolume );

  try
    {
    codeImageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }

  typedef signed short PixelType;

  typedef itk::Image<PixelType, 3>        ImageType;
  typedef itk::ImageFileReader<ImageType> ReferenceImageReaderType;
  ReferenceImageReaderType::Pointer referenceImageReader = ReferenceImageReaderType::New();
  referenceImageReader->SetFileName( inputReferenceVolume );

  try
    {
    referenceImageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }

  typedef itk::OrientImageFilter<CodeImageType, ImageType> OrientFilterType;
  OrientFilterType::Pointer orientImageFilter = OrientFilterType::New();
  orientImageFilter->SetInput( referenceImageReader->GetOutput() );
  orientImageFilter->SetDesiredCoordinateDirection( codeImageReader->GetOutput()->GetDirection() );
  orientImageFilter->UseImageDirectionOn();
  try
    {
    orientImageFilter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << e << std::endl;
    throw;
    }

  // Read the transform
  GenericTransformType::Pointer baseTransform = itk::ReadTransformFromDisk(inputTransform);
  typedef itk::NearestNeighborInterpolateImageFunction<CodeImageType, double> InterpolatorFunctionType;
  InterpolatorFunctionType::Pointer interpolatorFunction = InterpolatorFunctionType::New();

  typedef itk::ResampleImageFilter<CodeImageType, CodeImageType> ResampleFilterType;
  ResampleFilterType::Pointer resample = ResampleFilterType::New();

    {
    resample->SetTransform( baseTransform );
    }
  std::cout << "Code image:  ";
  codeImageReader->GetOutput()->Print(std::cout);
  std::cout << "Reference image:  ";
  referenceImageReader->GetOutput()->Print(std::cout);

  CodeImageType::PointType        p1;
  CodeImageType::PointType        p2;
  itk::ContinuousIndex<double, 3> imageIndex;
  imageIndex[0] = 0; imageIndex[1] = 0; imageIndex[2] = 0;
  codeImageReader->GetOutput()->TransformContinuousIndexToPhysicalPoint(imageIndex, p1);
  p2 = resample->GetTransform()->TransformPoint(p1);
  std::cout << "Point " << p1 << " mapped to " << p2 << std::endl;

  resample->SetInput( codeImageReader->GetOutput() );
  resample->SetOutputParametersFromImage( orientImageFilter->GetOutput() );
  resample->SetDefaultPixelValue( 0 );
  resample->SetInterpolator( interpolatorFunction );
  try
    {
    resample->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    throw;
    }

  typedef itk::OrientImageFilter<CodeImageType, CodeImageType> OrientCodeFilterType;
  OrientCodeFilterType::Pointer orientCodeImageFilter = OrientCodeFilterType::New();
  orientCodeImageFilter->SetInput( resample->GetOutput() );
  orientCodeImageFilter->SetDesiredCoordinateDirection( codeImageReader->GetOutput()->GetDirection() );
  orientCodeImageFilter->UseImageDirectionOn();
  try
    {
    orientCodeImageFilter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << e << std::endl;
    throw;
    }

  CodeImageType::Pointer resampledImage = orientCodeImageFilter->GetOutput();
  resampledImage->SetMetaDataDictionary( codeImageReader->GetOutput()->GetMetaDataDictionary() );

  typedef itk::ImageFileWriter<CodeImageType> ImageFileWriterType;
  ImageFileWriterType::Pointer ImageWriter =  ImageFileWriterType::New();
  ImageWriter->UseCompressionOn();
  ImageWriter->SetFileName( outputVolume );
  ImageWriter->SetInput( resampledImage );
  try
    {
    ImageWriter->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }
  return EXIT_SUCCESS;
}
