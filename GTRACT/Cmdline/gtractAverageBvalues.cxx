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

#include <itkArray.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkExceptionObject.h>
#include <itkMetaDataObject.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkVectorIndexSelectionCastImageFilter.h>

#include "gtractAverageBvaluesCLP.h"
#include "BRAINSThreadControl.h"
#include "DWIMetaDataDictionaryValidator.h"
#include <BRAINSCommonLib.h>
#include "DWIConvertLib.h"

int buildDirectionLut(itk::Array<int> & lut, itk::Array<int> & count, itk::MetaDataDictionary meta, int numImages,
                      double directionsTolerance, bool averageB0only);

bool areDirectionsEqual(std::string direction1, std::string direction2, double directionsTolerance, bool averageB0only);

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
    std::cout << "Directions Tolerance: " <<  directionsTolerance << std::endl;
    std::cout << "Average B0 Only: " <<  averageB0only << std::endl;
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
  if (convertInputVolumeToNrrdOrNifti(detectOuputVolumeType(outputVolume),
                                      inputVolume,convertedVolume)){
    inputVolume = convertedVolume;
  }
  else{
    std::cout<<"Error: DWI Convert can not read inputVolume."<<std::endl;
    return -1;
  }

  typedef signed short                   PixelType;
  typedef itk::VectorImage<PixelType, 3> NrrdImageType;
  typedef itk::Image<PixelType, 3>       IndexImageType;

  typedef float                             AvgPixelType;
  typedef itk::VectorImage<AvgPixelType, 3> NrrdAvgImageType;

  typedef itk::ImageFileReader<NrrdImageType,
                               itk::DefaultConvertPixelTraits<PixelType> > FileReaderType;
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

  itk::Array<int> lut;
  lut.SetSize( imageReader->GetOutput()->GetVectorLength() );
  lut.Fill( 0 );

  itk::Array<int> count;
  count.SetSize( imageReader->GetOutput()->GetVectorLength() );
  count.Fill( 0 );

  if( debug )
    {
    std::cout << "Original Number of Directions: " <<  imageReader->GetOutput()->GetVectorLength()
              << std::endl;
    }

  const int vectorLength = imageReader->GetOutput()->GetVectorLength();
  int       numUniqueDirections = buildDirectionLut(lut, count,
                                                    imageReader->GetOutput()->GetMetaDataDictionary(),
                                                    imageReader->GetOutput()->GetVectorLength(),
                                                    directionsTolerance, averageB0only);
  if( debug )
    {
    std::cout << "Avg #Directions: " << numUniqueDirections << std::endl;
    }

  // for (int i=0;i<vectorLength;i++)
  //  std::cout << i << " " << lut[i] << " " << count[i] << std::endl;

  NrrdAvgImageType::Pointer avgImage = NrrdAvgImageType::New();
  avgImage->SetRegions( imageReader->GetOutput()->GetLargestPossibleRegion() );
  avgImage->SetSpacing( imageReader->GetOutput()->GetSpacing() );
  avgImage->SetOrigin( imageReader->GetOutput()->GetOrigin() );
  avgImage->SetDirection( imageReader->GetOutput()->GetDirection() );
  avgImage->SetVectorLength( numUniqueDirections );
  avgImage->Allocate();

  NrrdImageType::Pointer outputImage = NrrdImageType::New();
  outputImage->SetRegions( imageReader->GetOutput()->GetLargestPossibleRegion() );
  outputImage->SetSpacing( imageReader->GetOutput()->GetSpacing() );
  outputImage->SetOrigin( imageReader->GetOutput()->GetOrigin() );
  outputImage->SetDirection( imageReader->GetOutput()->GetDirection() );
  outputImage->SetVectorLength( numUniqueDirections );
  outputImage->Allocate();

  typedef itk::VectorIndexSelectionCastImageFilter<NrrdImageType, IndexImageType> ExtractImageFilterType;
  ExtractImageFilterType::Pointer extractImageFilter = ExtractImageFilterType::New();
  extractImageFilter->SetInput( imageReader->GetOutput() );

  typedef itk::ImageRegionConstIterator<IndexImageType> ConstIndexImageIteratorType;
  typedef itk::ImageRegionIterator<NrrdAvgImageType>    VectorImageIteratorType;
  typedef itk::ImageRegionIterator<NrrdImageType>       VectorOutputImageIteratorType;
  for( int i = 0; i < vectorLength; i++ )
    {
    int currentIndex = lut[i];

    extractImageFilter->SetIndex( i );
    extractImageFilter->Update();

    ConstIndexImageIteratorType it( extractImageFilter->GetOutput(),
                                    extractImageFilter->GetOutput()->GetRequestedRegion() );
    VectorImageIteratorType ot( avgImage, avgImage->GetRequestedRegion() );
    for( ot.GoToBegin(), it.GoToBegin(); !ot.IsAtEnd(); ++ot, ++it )
      {
      NrrdAvgImageType::PixelType vectorImagePixel = ot.Get();
      vectorImagePixel[currentIndex] += static_cast<AvgPixelType>( it.Value() );
      ot.Set( vectorImagePixel );
      }
    }

  VectorOutputImageIteratorType tmpt( outputImage, outputImage->GetRequestedRegion() );
  VectorImageIteratorType       ot( avgImage, avgImage->GetRequestedRegion() );
  for( ot.GoToBegin(), tmpt.GoToBegin(); !ot.IsAtEnd(); ++ot, ++tmpt )
    {
    NrrdAvgImageType::PixelType vectorImagePixel = ot.Get();
    NrrdImageType::PixelType    outputImagePixel = tmpt.Get();
    for( int i = 0; i < numUniqueDirections; i++ )
      {
      outputImagePixel[i] = static_cast<PixelType>( vectorImagePixel[i] / static_cast<AvgPixelType>( count[i] ) );
      }
    tmpt.Set( outputImagePixel );
    }

  /* Update the Meta data Header */
  DWIMetaDataDictionaryValidator metaDataValidator;
  metaDataValidator.SetMetaDataDictionary(imageReader->GetOutput()->GetMetaDataDictionary());

  // Get gradient table and update the gradient vectors
  DWIMetaDataDictionaryValidator::GradientTableType gradTable = metaDataValidator.GetGradientTable();
  // Now delete the gradient table to fill with new gradient values
  metaDataValidator.DeleteGradientTable();
  // Update the validator using a new gradient table
  DWIMetaDataDictionaryValidator::GradientTableType newGradTable( numUniqueDirections );
  int currentIndex = 0;
  for( int i = 0; i < vectorLength; i++ )
    {
    if( lut[i] == currentIndex )
      {
      newGradTable[currentIndex][0] = gradTable[i][0];
      newGradTable[currentIndex][1] = gradTable[i][1];
      newGradTable[currentIndex][2] = gradTable[i][2];
      currentIndex++;
      }
    }
  metaDataValidator.SetGradientTable( newGradTable );

  outputImage->SetMetaDataDictionary( metaDataValidator.GetMetaDataDictionary() );

  typedef itk::ImageFileWriter<NrrdImageType> WriterType;
  WriterType::Pointer nrrdWriter = WriterType::New();
  nrrdWriter->UseCompressionOn();
  nrrdWriter->UseInputMetaDataDictionaryOn();
  nrrdWriter->SetInput( outputImage );
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

int buildDirectionLut(itk::Array<int> & lut,
                      itk::Array<int> & count,
                      itk::MetaDataDictionary meta,
                      int numImages,
                      double directionsTolerance,
                      bool averageB0only)
{
  int numElements = 0;

  for( int i = 0; i < numImages; i++ )
    {
    std::string direction1;
    std::string direction2;
    char        tmpStr[64];
    sprintf(tmpStr, "DWMRI_gradient_%04d", i);

    itk::ExposeMetaData<std::string>(meta, tmpStr, direction1);
    int j;
    for( j = 0; j < i; j++ )
      {
      sprintf(tmpStr, "DWMRI_gradient_%04d", j);
      itk::ExposeMetaData<std::string>(meta, tmpStr, direction2);
      if( areDirectionsEqual(direction1, direction2, directionsTolerance, averageB0only) )
        {
        if( lut[i] == 0 )
          {
          lut[i] = j;
          count[j]++;
          break;
          }
        }
      }
    if( i == j )
      {
      lut[i] = numElements;
      count[numElements]++;
      numElements++;
      }
    }
  // for (int i=0;i<numImages;i++) std::cout << i << " " << lut[i] << " " <<
  // count[i] << std::endl;

  /* Shuffle the count Elements down to the start of the Array */
  int index = 0;
  for( int i = 0; i < numImages; i++ )
    {
    if( count[i] > 0 )
      {
      count[index] = count[i];
      index++;
      }
    }
  for( int i = numElements; i < numImages; i++ )
    {
    count[i] = 0;
    }

  // for (int i=0;i<numImages;i++) std::cout << i << " " << lut[i] << " " <<
  // count[i] << std::endl;

  return numElements;
}

bool areDirectionsEqual(std::string direction1, std::string direction2, double directionsTolerance, bool averageB0only)
{
  const unsigned int MAXSTR = 256;
  char               tmpDir1[MAXSTR];
  char               tmpDir2[MAXSTR];

  strncpy( tmpDir1, direction1.c_str(), MAXSTR - 1 );
  strncpy( tmpDir2, direction2.c_str(), MAXSTR - 1 );

  const double x1 = atof( strtok(tmpDir1, " ") );
  const double y1 = atof( strtok(nullptr, " ") );
  const double z1 = atof( strtok(nullptr, " ") );

  const double x2 = atof( strtok(tmpDir2, " ") );
  const double y2 = atof( strtok(nullptr, " ") );
  const double z2 = atof( strtok(nullptr, " ") );

  if( averageB0only )
    {
    if( ( x1 == 0.0 ) && ( y1 == 0.0 ) && ( z1 == 0.0 ) && ( x2 == 0.0 ) && ( y2 == 0.0 ) && ( z2 == 0.0 ) )
      {
      return true;
      }
    else
      {
      return false;
      }
    }
  else
    {
    const double dist = std::sqrt( ( x1 - x2 ) * ( x1 - x2 ) + ( y1 - y2 ) * ( y1 - y2 ) + ( z1 - z2 ) * ( z1 - z2 ) );
    if( dist > directionsTolerance )
      {
      return false;
      }
    }
  return true;
}
