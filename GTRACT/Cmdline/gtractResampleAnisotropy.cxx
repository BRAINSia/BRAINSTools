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

#include "itkGtractImageIO.h"
#include "gtractResampleAnisotropyCLP.h"
#include "GenericTransformImage.h"
#include "BRAINSThreadControl.h"
int main(int argc, char *argv[])
{
  PARSE_ARGS;
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  itk::AddExtraTransformRegister();

  typedef double                                                          BSplineCoordinateRepType;
  typedef itk::VersorRigid3DTransform<double>                             RigidTransformType;
  typedef itk::BSplineDeformableTransform<BSplineCoordinateRepType, 3, 3> BSplineTransformType;

  bool debug = true;
  if( debug )
    {
    std::cout << "==============================================================" << std::endl;
    std::cout << "Input Image: " <<  inputAnisotropyVolume << std::endl;
    std::cout << "Output Image: " <<  outputVolume << std::endl;
    std::cout << "Anatomical Image: " <<  inputAnatomicalVolume << std::endl;
    std::cout << "Transform File: " << inputTransform << std::endl;
    std::cout << "Transform Type: " << transformType << std::endl;
    std::cout << "==============================================================" << std::endl;
    }

  bool violated = false;
  if( inputAnisotropyVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputAnisotropyVolume Required! "
                               << std::endl;
    }
  if( inputAnatomicalVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputAnatomicalVolume Required! "
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

  typedef float AnisotropyPixelType;

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

  typedef signed short PixelType;

  typedef itk::Image<PixelType, 3>        ImageType;
  typedef itk::ImageFileReader<ImageType> AnatomicalImageReaderType;
  AnatomicalImageReaderType::Pointer anatomicalReader = AnatomicalImageReaderType::New();
  anatomicalReader->SetFileName( inputAnatomicalVolume );

  try
    {
    anatomicalReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }

  // Read the transform
  GenericTransformType::Pointer baseTransform = itk::ReadTransformFromDisk(inputTransform);

  typedef itk::Image<float, 3>                                               AnisotropyImageType;
  typedef itk::ResampleImageFilter<AnisotropyImageType, AnisotropyImageType> ResampleFilterType;

  ResampleFilterType::Pointer resample = ResampleFilterType::New();
    {
    resample->SetTransform( baseTransform );
    }
  resample->SetInput( anisotropyImageReader->GetOutput() );
  resample->SetOutputParametersFromImage( anatomicalReader->GetOutput() );
  resample->SetDefaultPixelValue( 0 );
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

  AnisotropyImageType::Pointer resampledImage = resample->GetOutput();
  resampledImage->SetMetaDataDictionary( anatomicalReader->GetOutput()->GetMetaDataDictionary() );

  typedef itk::ImageFileWriter<AnisotropyImageType> ImageFileWriterType;
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
