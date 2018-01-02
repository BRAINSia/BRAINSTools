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

#include <cmath>

#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkExceptionObject.h>
#include <itkMetaDataObject.h>
#include <itkImageRegionIteratorWithIndex.h>

#include "gtractClipAnisotropyCLP.h"
#include "BRAINSThreadControl.h"
#include <BRAINSCommonLib.h>
#include "DWIConvertLib.h"

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
    std::cout << "Clip First Slice: " <<  clipFirstSlice << std::endl;
    std::cout << "Clip Last Slice: " <<  clipLastSlice << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  bool violated = false;
  if( inputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputVolume Required! "  << std::endl;
    }
  if( outputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --outputVolume Required! "  << std::endl;
    }
  if( violated )
    {
    return EXIT_FAILURE;
    }

  std::string convertedVolume;
  DWIConvert dwiConvert;
  if (0 == dwiConvert.convertInputVolumeToNrrdOrNifti(dwiConvert.detectOuputVolumeType(outputVolume),
                                                      inputVolume,convertedVolume)){
    inputVolume = convertedVolume;
  }
  else{
    std::cout<<"Error: DWI Convert can not read inputVolume."<<std::endl;
    return -1;
  }

  typedef float                      PixelType;
  typedef itk::Image<PixelType, 3>   ImageType;
  typedef ImageType::RegionType      ImageRegionType;
  typedef ImageRegionType::IndexType IndexType;

  typedef itk::ImageFileReader<ImageType> FileReaderType;
  FileReaderType::Pointer imageReader = FileReaderType::New();
  imageReader->SetFileName( inputVolume );

  try
    {
    imageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }

  ImageType::Pointer originalImage = imageReader->GetOutput();

  typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIteratorType;
  ImageIteratorType it( originalImage, originalImage->GetLargestPossibleRegion() );
  ImageRegionType   region = originalImage->GetLargestPossibleRegion();
  int               lastSlice = region.GetSize(2) - 1;
  std::cout << "Last Slice " <<  lastSlice << std::endl;
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    IndexType index = it.GetIndex();
    if( ( index[2] == 0 ) && clipFirstSlice )
      {
      it.Set( 0.0 );
      }
    if( ( index[2] == lastSlice ) && clipLastSlice )
      {
      it.Set( 0.0 );
      }
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer nrrdWriter = WriterType::New();
  nrrdWriter->UseCompressionOn();
  nrrdWriter->UseInputMetaDataDictionaryOn();
  nrrdWriter->SetInput( originalImage );
  nrrdWriter->SetFileName( outputVolume );
  try
    {
    nrrdWriter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << e << std::endl;
    }
  return EXIT_SUCCESS;
}
