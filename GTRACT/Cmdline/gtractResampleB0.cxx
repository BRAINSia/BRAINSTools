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
#include <itkExtractImageFilter.h>
#include <itkBSplineDeformableTransform.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkTransformFactory.h>
#include <itkVersorRigid3DTransform.h>
#include <itkOrientImageFilter.h>
#include "BRAINSThreadControl.h"
#include "itkGtractImageIO.h"
#include "gtractResampleB0CLP.h"
#include "GenericTransformImage.h"

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
    std::cout << "=====================================================" << std::endl;
    std::cout << "Input Image: " <<  inputVolume << std::endl;
    std::cout << "Output Image: " <<  outputVolume << std::endl;
    std::cout << "Anatomical Image: " <<  inputAnatomicalVolume << std::endl;
    std::cout << "Transform File: " << inputTransform << std::endl;
    std::cout << "Transform Type: " << transformType << std::endl;
    std::cout << "Vector Index: " << vectorIndex << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  bool violated = false;
  if( inputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputVolume Required! "  << std::endl;
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

  typedef signed short PixelType;

  typedef itk::VectorImage<PixelType, 3> VectorImageType;
  typedef itk::ImageFileReader<VectorImageType,
                               itk::DefaultConvertPixelTraits<PixelType> > VectorImageReaderType;
  VectorImageReaderType::Pointer vectorImageReader = VectorImageReaderType::New();
  vectorImageReader->SetFileName( inputVolume );

  try
    {
    vectorImageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }

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

  typedef itk::VectorImage<PixelType, 3> VectorImageType;
  VectorImageType::Pointer    dwiImage = vectorImageReader->GetOutput();
  VectorImageType::RegionType fixedRegion = dwiImage->GetLargestPossibleRegion();
  VectorImageType::SizeType   fixedSize = fixedRegion.GetSize();
  // const VectorImageType::SpacingType fixedSpacing = dwiImage->GetSpacing();
  // const VectorImageType::PointType   fixedOrigin = dwiImage->GetOrigin();

  fixedSize[3] = 0;
  fixedRegion.SetSize(fixedSize);

  VectorImageType::IndexType fixedIndex = fixedRegion.GetIndex();
  fixedIndex[0] = 0;
  fixedIndex[1] = 0;
  fixedIndex[2] = 0;
  fixedIndex[3] = 0;
  fixedRegion.SetIndex(fixedIndex);

  /* Extract the Vector Image Index for Registration */
  typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ImageType> VectorSelectFilterType;
  typedef VectorSelectFilterType::Pointer                                      VectorSelectFilterPointer;
  VectorSelectFilterPointer selectIndexImageFilter = VectorSelectFilterType::New();
  selectIndexImageFilter->SetIndex( vectorIndex );
  selectIndexImageFilter->SetInput( dwiImage );
  try
    {
    selectIndexImageFilter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << e << std::endl;
    throw;
    }

  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;

  ResampleFilterType::Pointer resample = ResampleFilterType::New();

  resample->SetTransform( baseTransform );

  resample->SetInput( selectIndexImageFilter->GetOutput() );
  resample->SetOutputParametersFromImage( anatomicalReader->GetOutput() );
  resample->SetDefaultPixelValue( 1 );
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

  ImageType::Pointer resampledImage = resample->GetOutput();
  resampledImage->SetMetaDataDictionary( anatomicalReader->GetOutput()->GetMetaDataDictionary() );

  typedef itk::ImageFileWriter<ImageType> ImageFileWriterType;
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
  std::cout << "wrote image" << std::endl;
  return EXIT_SUCCESS;
}
