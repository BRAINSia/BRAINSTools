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

#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>

/* Defines the Additional Anisotropy Metrics */
#include "../Common/gtractDiffusionTensor3D.h"
#include "algo.h"
#include "GtractTypes.h"
#include "gtractAnisotropyMapCLP.h"
#include "BRAINSThreadControl.h"
#include <BRAINSCommonLib.h>

int main(int argc, char *argv[])
{
  typedef double                                            TensorComponentType;
  typedef itk::gtractDiffusionTensor3D<TensorComponentType> TensorPixelType;
  typedef itk::Image<TensorPixelType, 3>                    TensorImageType;
  typedef itk::Image<float, 3>                              AnisotropyImageType;

  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  bool debug = true;
  if( debug )
    {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Input Tensor Image: " <<  inputTensorVolume << std::endl;
    std::cout << "Output Anisotropy Image: " <<  outputVolume << std::endl;
    std::cout << "Anisotropy Type: " <<  anisotropyType << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  bool violated = false;
  if( inputTensorVolume.size() == 0 )
    {
    violated = true; std::cout << "  --inputTensorVolume Required! "  << std::endl;
    }
  if( outputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --outputVolume Required! "  << std::endl;
    }
  if( violated )
    {
    return EXIT_FAILURE;
    }

  typedef itk::ImageFileReader<TensorImageType> TensorImageReaderType;
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

  AnisotropyImageType::Pointer anisotropyImage =  AnisotropyImageType::New();
  anisotropyImage->SetRegions( tensorImage->GetLargestPossibleRegion() );
  anisotropyImage->SetSpacing( tensorImage->GetSpacing() );
  anisotropyImage->SetOrigin( tensorImage->GetOrigin() );
  anisotropyImage->SetDirection( tensorImage->GetDirection() );
  anisotropyImage->Allocate();

  typedef itk::ImageRegionIterator<AnisotropyImageType> IteratorType;
  IteratorType anisoIt( anisotropyImage, anisotropyImage->GetRequestedRegion() );

  typedef itk::ImageRegionConstIterator<TensorImageType> ConstIteratorType;
  ConstIteratorType tensorIt( tensorImage, tensorImage->GetRequestedRegion() );

  float anisotropy = 0.0;
  for( anisoIt.GoToBegin(), tensorIt.GoToBegin(); !anisoIt.IsAtEnd(); ++anisoIt, ++tensorIt )
    {
    TensorPixelType tensorPixel = tensorIt.Get();
    if( anisotropyType == "ADC"  ||  anisotropyType == "adc" )
      {
      anisotropy = static_cast<float>( tensorPixel.GetTrace() / 3.0 );
      }
    else if( anisotropyType == "FA"  ||  anisotropyType == "fa" )
      {
      anisotropy = static_cast<float>( tensorPixel.GetFractionalAnisotropy() );
      }
    else if( anisotropyType == "RA"  ||  anisotropyType == "ra" )
      {
      anisotropy = static_cast<float>( tensorPixel.GetRelativeAnisotropy() );
      }
    else if( anisotropyType == "VR"  ||  anisotropyType == "vr" )
      {
      anisotropy = static_cast<float>( tensorPixel.GetVolumeRatio() );
      }
    else if( anisotropyType == "AD"  ||  anisotropyType == "ad" )
      {
      anisotropy = static_cast<float>( tensorPixel.GetAxialDiffusivity() );
      }
    else if( anisotropyType == "RD"  ||  anisotropyType == "rd" )
      {
      anisotropy = static_cast<float>( tensorPixel.GetRadialDiffusivity() );
      }
    else if( anisotropyType == "LI"  ||  anisotropyType == "li" )
      {
      anisotropy = static_cast<float>( tensorPixel.GetLatticeIndex() );
      }
    anisoIt.Set( anisotropy );
    }

  typedef itk::ImageFileWriter<AnisotropyImageType> WriterType;
  WriterType::Pointer anisotropyWriter = WriterType::New();
  anisotropyWriter->UseCompressionOn();
  anisotropyWriter->SetInput( anisotropyImage );
  anisotropyWriter->SetFileName( outputVolume );
  try
    {
    anisotropyWriter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << e << std::endl;
    }
  return EXIT_SUCCESS;
}
