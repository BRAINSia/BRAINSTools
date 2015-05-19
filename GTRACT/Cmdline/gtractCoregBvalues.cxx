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
#include <itkVectorImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkExceptionObject.h>
#include <itkVectorImageToImageAdaptor.h>
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkTransformFileWriter.h>

#include "BRAINSFitHelper.h"

#include "gtractCoregBvaluesCLP.h"
#include "BRAINSThreadControl.h"
#include "itkOrthogonalize3DRotationMatrix.h"
#include "itkNumberToString.h"

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  itk::NumberToString<double>                                        doubleConvert;
  bool                                                  debug = true;
  if( debug )
    {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Moving Image: " <<  movingVolume << std::endl;
    std::cout << "Fixed Image: " <<  fixedVolume << std::endl;
    std::cout << "Fixed Image Index: " <<  fixedVolumeIndex << std::endl;
    std::cout << "Output Image: " <<  outputVolume << std::endl;
    std::cout << "Output Transform: " <<  outputTransform << std::endl;
    std::cout << "Eddy Current Correction: " << eddyCurrentCorrection << std::endl;
    std::cout << "Iterations: " << numberOfIterations << std::endl;
    std::cout << "Translation Scale: " << spatialScale << std::endl;
    std::cout << "Maximum Step Size: " << maximumStepSize << std::endl;
    std::cout << "Minimum Step Size: " << minimumStepSize << std::endl;
    std::cout << "Relaxation Factor: " << relaxationFactor << std::endl;
    std::cout << "Samples: " << numberOfSpatialSamples << std::endl;
    std::cout << "Register B0 Only: " << registerB0Only << std::endl;
    std::cout << "Debug Level: " << debugLevel << std::endl;
    std::cout << "=====================================================" << std::endl;
    }

  bool violated = false;
  if( movingVolume.size() == 0 )
    {
    violated = true; std::cout << "  --movingVolume Required! "  << std::endl;
    }
  if( fixedVolume.size() == 0 )
    {
    violated = true; std::cout << "  --fixedVolume Required! "  << std::endl;
    }
  if( outputVolume.size() == 0 )
    {
    violated = true; std::cout << "  --outputVolume Required! "  << std::endl;
    }
  if( outputTransform.size() == 0 )
    {
    violated = true; std::cout << "  --outputTransform Required! "  << std::endl;
    }
  if( violated )
    {
    return EXIT_FAILURE;
    }

  typedef signed short                         OutputPixelType;
  typedef float                                PixelType;
  typedef itk::VectorImage<PixelType, 3>       NrrdImageType;
  typedef itk::VectorImage<OutputPixelType, 3> OutputImageType;
  typedef itk::Image<PixelType, 3>             InputIndexImageType;
  typedef itk::Image<OutputPixelType, 3>       OutputIndexImageType;
  typedef itk::AffineTransform<double, 3>      LocalAffineTransformType;
  typedef itk::VersorRigid3DTransform<double>  RigidTransformType;

  typedef itk::ImageFileReader<NrrdImageType,
                               itk::DefaultConvertPixelTraits<PixelType> > FileReaderType;
  FileReaderType::Pointer movingImageReader = FileReaderType::New();
  movingImageReader->SetFileName( movingVolume );
  movingImageReader->Update();

  try
    {
    movingImageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }

  // typedef itk::VectorIndexSelectionCastImageFilter< NrrdImageType,
  // IndexImageType > InputImageType;
  typedef itk::VectorIndexSelectionCastImageFilter<NrrdImageType, InputIndexImageType> ExtractImageFilterType;
  ExtractImageFilterType::Pointer movingImageExtractionFilter = ExtractImageFilterType::New();
  movingImageExtractionFilter->SetInput(movingImageReader->GetOutput() );

  typedef itk::ImageFileReader<NrrdImageType,
                               itk::DefaultConvertPixelTraits<PixelType> > FileReaderType;
  FileReaderType::Pointer fixedImageReader = FileReaderType::New();
  fixedImageReader->SetFileName( fixedVolume );
  fixedImageReader->Update();

  try
    {
    fixedImageReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }

  /* Extract Image Index to be used for Coregistration */
  ExtractImageFilterType::Pointer fixedImageExtractionFilter = ExtractImageFilterType::New();
  fixedImageExtractionFilter->SetIndex( fixedVolumeIndex );
  fixedImageExtractionFilter->SetInput( fixedImageReader->GetOutput() );
  fixedImageExtractionFilter->Update();

  /* Pointer Used to Hold the Resulting Coregistered Image */
  OutputImageType::Pointer RegisteredImage;

  typedef itk::BRAINSFitHelper RegisterFilterType;
  RegisterFilterType::Pointer registerImageFilter = RegisterFilterType::New();

  std::vector<double> minStepLength;
  minStepLength.push_back( (double)minimumStepSize);

  std::vector<std::string> rigidTransformTypes;
  rigidTransformTypes.push_back("ScaleVersor3D");

  std::vector<std::string> affineTransformTypes;
  affineTransformTypes.push_back("Affine");

  std::vector<int> iterations;
  iterations.push_back(numberOfIterations);

  // Allocate output image
  RegisteredImage = OutputImageType::New();
  RegisteredImage->SetRegions( movingImageReader->GetOutput()->GetLargestPossibleRegion() );
  RegisteredImage->SetSpacing( movingImageReader->GetOutput()->GetSpacing() );
  RegisteredImage->SetOrigin( movingImageReader->GetOutput()->GetOrigin() );
  RegisteredImage->SetDirection( movingImageReader->GetOutput()->GetDirection() );
  RegisteredImage->SetVectorLength( movingImageReader->GetOutput()->GetVectorLength() );
  RegisteredImage->SetMetaDataDictionary( movingImageReader->GetOutput()->GetMetaDataDictionary() );
  RegisteredImage->Allocate();

  typedef itk::ImageRegionConstIterator<OutputIndexImageType> ConstIteratorType;
  typedef itk::ImageRegionIterator<OutputImageType>           IteratorType;

  IteratorType               ot( RegisteredImage, RegisteredImage->GetRequestedRegion() );
  OutputImageType::PixelType vectorImagePixel;
  for( unsigned int i = 0; i < movingImageReader->GetOutput()->GetVectorLength(); i++ )
    {
    // Get Current Gradient Direction
    vnl_vector<double> curGradientDirection(3);
    char               tmpStr[64];
    sprintf(tmpStr, "DWMRI_gradient_%04u", i);
    std::string KeyString(tmpStr);
    std::string NrrdValue;

    itk::MetaDataDictionary inputMetaDataDictionary = movingImageReader->GetOutput()->GetMetaDataDictionary();
    itk::ExposeMetaData<std::string>(inputMetaDataDictionary, KeyString, NrrdValue);
    /* %lf is 'long float', i.e., double. */
    // std::cout<<"KeyString: "<<KeyString<<std::endl;
    // std::cout<<"NrrdValue: "<<NrrdValue<<std::endl;
    sscanf(
      NrrdValue.c_str(), " %lf %lf %lf", &curGradientDirection[0], &curGradientDirection[1], &curGradientDirection[2]);
    // std::cout<<" before: "<<std::endl;
    // std::cout<<curGradientDirection<<std::endl;

    if( eddyCurrentCorrection == 0 )
      {
      std::cout << "Rigid Registration: " << std::endl;
      registerImageFilter->SetTransformType(rigidTransformTypes);
      }
    else
      {
      std::cout << "Full Affine Registration: " << std::endl;
      registerImageFilter->SetTransformType(affineTransformTypes);
      }
    movingImageExtractionFilter->SetIndex( i );
    movingImageExtractionFilter->Update();

    registerImageFilter->SetTranslationScale( spatialScale );
    registerImageFilter->SetMaximumStepLength( maximumStepSize );
    registerImageFilter->SetMinimumStepLength( minStepLength );
    registerImageFilter->SetRelaxationFactor( relaxationFactor );
    registerImageFilter->SetNumberOfIterations( iterations );
    if(numberOfSpatialSamples > 0)
      {
        const unsigned long numberOfAllSamples = fixedImageExtractionFilter->GetOutput()->GetBufferedRegion().GetNumberOfPixels();
        samplingPercentage = static_cast<double>( numberOfSpatialSamples )/numberOfAllSamples;
        std::cout << "WARNING --numberOfSpatialSamples is deprecated, please use --samplingPercentage instead " << std::endl;
        std::cout << "WARNING: Replacing command line --samplingPercentage " << samplingPercentage << std::endl;
      }
    registerImageFilter->SetSamplingPercentage( samplingPercentage );
    registerImageFilter->SetMovingVolume( movingImageExtractionFilter->GetOutput() );
    registerImageFilter->SetFixedVolume( fixedImageExtractionFilter->GetOutput() );
    registerImageFilter->SetDebugLevel(debugLevel);
    registerImageFilter->SetInitializeTransformMode("useMomentsAlign" );


    try
      {
      registerImageFilter->Update();
      }
    catch( itk::ExceptionObject & ex )
      {
      std::cout << ex << std::endl;
      throw;
      }
    typedef itk::Transform<double, 3, 3> GenericTransformType;
    // restore writing out transform if specified on command line.
    if( outputTransform.size() != 0 )
      {
      GenericTransformType::Pointer transform =
        registerImageFilter->GetCurrentGenericTransform()->GetNthTransform(0);
      itk::TransformFileWriter::Pointer xfrmWriter =
        itk::TransformFileWriter::New();
      xfrmWriter->SetFileName(outputTransform);
      xfrmWriter->SetInput(transform);
      xfrmWriter->Update();
      }
    typedef itk::ResampleImageFilter<InputIndexImageType, OutputIndexImageType, double> ResampleFilterType;
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetTransform( registerImageFilter->GetCurrentGenericTransform() );
    resampler->SetInput( movingImageExtractionFilter->GetOutput() );
    // Remember:  the Data is Moving's, the shape is Fixed's.
    resampler->SetOutputParametersFromImage(fixedImageExtractionFilter->GetOutput() );
    resampler->SetDefaultPixelValue( 0 );
    resampler->Update();

    if( eddyCurrentCorrection == 0 )
      {
      RigidTransformType::Pointer rigidTransform;
      rigidTransform =
        dynamic_cast<RigidTransformType *>(registerImageFilter->GetCurrentGenericTransform().GetPointer() );
      curGradientDirection = rigidTransform->GetMatrix().GetVnlMatrix() * curGradientDirection;
      }
    else
      {
      LocalAffineTransformType::Pointer affineTransform;
      affineTransform =
        dynamic_cast<LocalAffineTransformType *>(registerImageFilter->GetCurrentGenericTransform().GetPointer() );
      itk::Matrix<double, 3, 3> NonOrthog = affineTransform->GetMatrix();
      itk::Matrix<double, 3, 3> Orthog( itk::Orthogonalize3DRotationMatrix(NonOrthog) );
      curGradientDirection = Orthog.GetVnlMatrix() * curGradientDirection;
      // std::cout<<"i = "<<i<<std::endl;
      // std::cout<<curGradientDirection<<std::endl;
      }

    // Write RegisteredImage
    ConstIteratorType it( resampler->GetOutput(), resampler->GetOutput()->GetRequestedRegion() );
    for( ot.GoToBegin(), it.GoToBegin(); !ot.IsAtEnd(); ++ot, ++it )
      {
      vectorImagePixel = ot.Get();
      vectorImagePixel[i] = it.Value();
      ot.Set( vectorImagePixel );
      }

    // Add the gradient direction to the resulting image
    // sprintf(tmpStr, " %18.15lf %18.15lf %18.15lf", curGradientDirection[0], curGradientDirection[1],
    //         curGradientDirection[2]);
    // NrrdValue = tmpStr;
    NrrdValue = " ";
    for( unsigned dir = 0; dir < 3; ++dir )
      {
      if( i > 0 )
        {
        NrrdValue += " ";
        }
      NrrdValue += doubleConvert(curGradientDirection[dir]);
      }
    itk::EncapsulateMetaData<std::string>(RegisteredImage->GetMetaDataDictionary(), KeyString, NrrdValue);
    }

  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  WriterType::Pointer nrrdWriter = WriterType::New();
  nrrdWriter->UseCompressionOn();
  nrrdWriter->UseInputMetaDataDictionaryOn();
  nrrdWriter->SetInput( RegisteredImage );
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
