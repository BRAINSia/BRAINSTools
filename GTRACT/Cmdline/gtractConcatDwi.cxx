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
#include <itkComposeImageFilter.h>
#include "itkNumberToString.h"


#include "BRAINSThreadControl.h"
#include "DWIMetaDataDictionaryValidator.h"
#include <BRAINSCommonLib.h>

#include "gtractConcatDwiCLP.h"
#include "DWIConvertLib.h"
int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  const int                                             numberOfImages = inputVolume.size();
  bool                                                  debug = true;
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

  std::vector<std::string> inputVolumeNrrd;
  DWIConvert dwiConvert;
  if (0 == dwiConvert.convertInputVolumeVectorToNrrdOrNifti(dwiConvert.detectOuputVolumeType(outputVolume),
                                                            inputVolume,inputVolumeNrrd)){
    inputVolume = inputVolumeNrrd;
  }
  else{
    std::cout<<"Error: DWI Convert can not read inputVolume."<<std::endl;
    return -1;
  }

  typedef signed short                   PixelType;
  typedef itk::VectorImage<PixelType, 3> NrrdImageType;
  typedef itk::Image<PixelType, 3>       IndexImageType;

  typedef itk::ImageFileReader<NrrdImageType,
                               itk::DefaultConvertPixelTraits<PixelType> > FileReaderType;

  DWIMetaDataDictionaryValidator currentMetaDataValidator;
  DWIMetaDataDictionaryValidator resultMetaDataValidator;
  DWIMetaDataDictionaryValidator::GradientTableType resultGradTable;

  typedef itk::ComposeImageFilter<IndexImageType> VectorImageFilterType;
  VectorImageFilterType::Pointer indexImageToVectorImageFilter = VectorImageFilterType::New();

  int                            vectorIndex = 0;
  double                         baselineBvalue = 0.0;

  NrrdImageType::PointType firstOrigin;
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
      std::cout << ex << std::endl << std::flush;
      throw;
      }

    currentMetaDataValidator.SetMetaDataDictionary(imageReader->GetOutput()->GetMetaDataDictionary());
    NrrdImageType::PointType currentOrigin = imageReader->GetOutput()->GetOrigin();
    if( i == 0 )
      {
      firstOrigin = currentOrigin;

      resultMetaDataValidator.SetMetaDataDictionary(imageReader->GetOutput()->GetMetaDataDictionary());
      resultMetaDataValidator.DeleteGradientTable();
      baselineBvalue = resultMetaDataValidator.GetBValue();
      }
    else
      {
      double distance =
        std::sqrt(firstOrigin.SquaredEuclideanDistanceTo(currentOrigin) );
      if( !ignoreOrigins && distance > 1.0E-3 )
        {
        std::cerr << "Origins differ " << firstOrigin
                  << " " << currentOrigin << std::endl;
        return EXIT_FAILURE;
        }
      else if( distance > 1.0E-6 )
        {
        // if there is a small difference make them the same
        imageReader->GetOutput()->SetOrigin(firstOrigin);
        }
      }
    DWIMetaDataDictionaryValidator::GradientTableType currGradTable = currentMetaDataValidator.GetGradientTable();
    double currentBvalue = currentMetaDataValidator.GetBValue();
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
      catch( itk::ExceptionObject & e )
        {
        std::cout << e << std::endl;
        }
      indexImageToVectorImageFilter->SetInput( vectorIndex, selectIndexImageFilter->GetOutput() );

      // Scale the current gradient and put in the result gradient table
      DWIMetaDataDictionaryValidator::Double3x1ArrayType scaledGradient;
      scaledGradient[0] = currGradTable[j][0] * bValueScale;
      scaledGradient[1] = currGradTable[j][1] * bValueScale;
      scaledGradient[2] = currGradTable[j][2] * bValueScale;
      resultGradTable.push_back(scaledGradient);

      vectorIndex++;
      }
    }
  resultMetaDataValidator.SetGradientTable( resultGradTable );

  try
    {
    indexImageToVectorImageFilter->Update();
    indexImageToVectorImageFilter->GetOutput()->SetMetaDataDictionary( resultMetaDataValidator.GetMetaDataDictionary() );
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl << std::flush;
    throw;
    }

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
  catch( itk::ExceptionObject & e )
    {
    std::cout << e << std::endl;
    }
  return EXIT_SUCCESS;
}
