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
#include <itkOrientImageFilter.h>

#include "itkGtractImageIO.h"
#include "gtractCopyImageOrientationCLP.h"
#include "GenericTransformImage.h"
#include "BRAINSThreadControl.h"
int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  bool debug = true;
  if( debug )
    {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Input Image: " <<  inputVolume << std::endl;
    std::cout << "Output Image: " <<  outputVolume << std::endl;
    std::cout << "Reference Image: " <<  inputReferenceVolume << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  bool violated = false;
  if( inputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputVolume Required! "  << std::endl;
    }
  if( inputReferenceVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputReferenceVolume Required! "
                               << std::endl;
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

  typedef itk::Image<PixelType, 3>                SpecimenImageType;
  typedef itk::ImageFileReader<SpecimenImageType> SpecimenImageReaderType;
  SpecimenImageReaderType::Pointer specimenImageReader = SpecimenImageReaderType::New();
  specimenImageReader->SetFileName( inputVolume );

  try
    {
    specimenImageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }

  typedef itk::Image<PixelType, 3>                 ReferenceImageType;
  typedef itk::ImageFileReader<ReferenceImageType> ReferenceImageReaderType;
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

  typedef itk::OrientImageFilter<SpecimenImageType, ReferenceImageType> OrientFilterType;
  OrientFilterType::Pointer orientImageFilter = OrientFilterType::New();
  orientImageFilter->SetInput( specimenImageReader->GetOutput() );
  orientImageFilter->SetDesiredCoordinateDirection( referenceImageReader->GetOutput()->GetDirection() );
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

  ReferenceImageType::Pointer reorientedImage = orientImageFilter->GetOutput();
  // reorientedImage->SetOrigin(referenceImageReader->GetOutput()->GetOrigin());
  reorientedImage->SetMetaDataDictionary( specimenImageReader->GetOutput()->GetMetaDataDictionary() );

  typedef itk::ImageFileWriter<ReferenceImageType> ImageFileWriterType;
  ImageFileWriterType::Pointer ImageWriter =  ImageFileWriterType::New();
  ImageWriter->UseCompressionOn();
  ImageWriter->SetFileName( outputVolume );
  ImageWriter->SetInput( reorientedImage );
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
