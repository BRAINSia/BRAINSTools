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

#include <itkArray.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <metaCommand.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkExceptionObject.h>
#include <itkMetaDataObject.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkImageToVectorImageFilter.h>

#include "gtractConcatDwiCLP.h"
#include "BRAINSThreadControl.h"
int main(int argc, char *argv[])
{
  PARSE_ARGS;
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  const int numberOfImages = inputVolume.size();
  bool      debug = true;
  if( debug )
    {
    std::cout << "=====================================================" << std::endl;
    for( int i = 0; i < numberOfImages; i++ )
      {
      std::cout << "Input Volume:                      " <<  inputVolume[i] << std::endl;
      }
    std::cout << "Output Image: " <<  outputVolume << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  bool violated = false;
  if( numberOfImages == 0 )
    {
    violated = true; std::cout << " at least one --inputVolume is required! "  << std::endl;
    }
  if( outputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --outputVolume Required! "  << std::endl;
    }
  if( violated )
    {
    return EXIT_FAILURE;
    }

  typedef signed short                   PixelType;
  typedef itk::VectorImage<PixelType, 3> NrrdImageType;
  typedef NrrdImageType::PixelType       VectorImagePixelType;
  typedef itk::Image<PixelType, 3>       IndexImageType;

  typedef itk::ImageFileReader<NrrdImageType,
                               itk::DefaultConvertPixelTraits<PixelType> > FileReaderType;
  itk::MetaDataDictionary resultMetaData;
  typedef itk::ImageToVectorImageFilter<IndexImageType> VectorImageFilterType;
  VectorImageFilterType::Pointer indexImageToVectorImageFilter = VectorImageFilterType::New();
  int                            vectorIndex = 0;
  double                         baselineBvalue = 0.0;
  for( unsigned i = 0; i < inputVolume.size(); i++ )
    {
    std::cout << "Reading volume:              " <<  inputVolume[i] << std::endl;
    FileReaderType::Pointer imageReader = FileReaderType::New();
    imageReader->SetFileName( inputVolume[i] );
    try
      {
      imageReader->Update();
      }
    catch( itk::ExceptionObject & ex )
      {
      std::cout << ex << std::endl;
      throw;
      }
    // InputImages[i] = imageReader->GetOutput();

    itk::MetaDataDictionary currentMetaData = imageReader->GetOutput()->GetMetaDataDictionary();
    std::string             NrrdValue;
    if( i == 0 )
      {
      resultMetaData = currentMetaData;
      itk::ExposeMetaData<std::string>(currentMetaData, "DWMRI_b-value", NrrdValue);
      baselineBvalue = atof( NrrdValue.c_str() );
      }

    itk::ExposeMetaData<std::string>(currentMetaData, "DWMRI_b-value", NrrdValue);
    double currentBvalue = atof( NrrdValue.c_str() );
    double bValueScale = currentBvalue / baselineBvalue;
    for( unsigned int j = 0; j < imageReader->GetOutput()->GetVectorLength(); j++ )
      {
      typedef itk::VectorIndexSelectionCastImageFilter<NrrdImageType, IndexImageType> VectorSelectFilterType;
      typedef VectorSelectFilterType::Pointer                                         VectorSelectFilterPointer;

      VectorSelectFilterPointer selectIndexImageFilter = VectorSelectFilterType::New();
      selectIndexImageFilter->SetIndex( j );
      selectIndexImageFilter->SetInput( imageReader->GetOutput() );
      try
        {
        selectIndexImageFilter->Update();
        }
      catch( itk::ExceptionObject e )
        {
        std::cout << e << std::endl;
        }
      indexImageToVectorImageFilter->SetInput( vectorIndex, selectIndexImageFilter->GetOutput() );

      char tmpStr[64];
      char tmpValue[64];
      char tokStr[64];

      sprintf(tmpStr, "DWMRI_gradient_%04d", j);
      itk::ExposeMetaData<std::string>(currentMetaData, tmpStr, NrrdValue);
      strcpy( tokStr, NrrdValue.c_str() );
      double x = atof( strtok( tokStr, " " ) );
      double y = atof( strtok( NULL, " " ) );
      double z = atof( strtok( NULL, " " ) );
      sprintf(tmpStr, "DWMRI_gradient_%04d", vectorIndex);
      sprintf(tmpValue, " %18.15lf %18.15lf %18.15lf", x * bValueScale, y * bValueScale, z * bValueScale);
      NrrdValue = tmpValue;
      itk::EncapsulateMetaData<std::string>(resultMetaData, tmpStr, NrrdValue);
      vectorIndex++;
      }
    }
  indexImageToVectorImageFilter->Update();
  indexImageToVectorImageFilter->GetOutput()->SetMetaDataDictionary(resultMetaData);

  typedef itk::ImageFileWriter<NrrdImageType> WriterType;
  WriterType::Pointer nrrdWriter = WriterType::New();
  nrrdWriter->UseCompressionOn();
  nrrdWriter->UseInputMetaDataDictionaryOn();
  nrrdWriter->SetInput( indexImageToVectorImageFilter->GetOutput() );
  nrrdWriter->SetFileName( outputVolume );
  try
    {
    nrrdWriter->Update();
    }
  catch( itk::ExceptionObject e )
    {
    std::cout << e << std::endl;
    }
  return EXIT_SUCCESS;
}
