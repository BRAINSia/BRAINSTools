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
#include <itkVectorIndexSelectionCastImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkOrientImageFilter.h>

#include "itkGtractImageIO.h"
#include "BRAINSFitHelper.h"
#include "gtractCoRegAnatomyRigidCLP.h"

int main(int argc, char * *argv)
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();

  bool debug = true;
  if( debug )
    {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Input Image: " <<  inputVolume << std::endl;
    std::cout << "Output Rigid Transform: " <<  outputRigidTransform << std::endl;
    std::cout << "Anatomical Image: " <<  inputAnatomicalVolume << std::endl;
    std::cout << "Translation Scale: " << spatialScale << std::endl;
    std::cout << "Maximum Step Length: " << maximumStepSize << std::endl;
    std::cout << "Minimum Step Length: " << minimumStepSize << std::endl;
    std::cout << "Relaxation Factor: " << relaxationFactor << std::endl;
    std::cout << "Iterations: " << numberOfIterations << std::endl;
    std::cout << "Samples: " << numberOfSamples << std::endl;
    std::cout << "Index: " << vectorIndex << std::endl;
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
  if( outputRigidTransform.size() == 0 )
    {
    violated = true; std::cout << "  --outputRigidTransform Required! "
                               << std::endl;
    }
  if( violated )
    {
    return EXIT_FAILURE;
    }

  //  typedef signed short                      PixelType;
  typedef float                          PixelType;
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

  typedef itk::Image<PixelType, 3>                  AnatomicalImageType;
  typedef itk::ImageFileReader<AnatomicalImageType> AnatomicalImageReaderType;
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

  // std::cout << "anatomicalReader produced:  ";
  // itk::Indent tabbed(4);
  // anatomicalReader->GetOutput()->Print(std::cout, tabbed);

  /* Extract the Vector Image Index for Registration */
  typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, AnatomicalImageType> VectorSelectFilterType;
  typedef VectorSelectFilterType::Pointer                                                VectorSelectFilterPointer;

  VectorSelectFilterPointer selectIndexImageFilter = VectorSelectFilterType::New();
  selectIndexImageFilter->SetIndex( vectorIndex );
  selectIndexImageFilter->SetInput( vectorImageReader->GetOutput() );
  try
    {
    selectIndexImageFilter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << e << std::endl;
    throw;
    }

  /* The Threshold Image Filter is used to produce the brain clipping mask. */
  typedef itk::ThresholdImageFilter<AnatomicalImageType> ThresholdFilterType;
  constexpr PixelType              imageThresholdBelow  = 100;
  ThresholdFilterType::Pointer brainOnlyFilter = ThresholdFilterType::New();
  brainOnlyFilter->SetInput( selectIndexImageFilter->GetOutput() );
  brainOnlyFilter->ThresholdBelow( imageThresholdBelow );
  try
    {
    brainOnlyFilter->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << e << std::endl;
    throw;
    }

  typedef itk::BRAINSFitHelper RegisterFilterType;
  RegisterFilterType::Pointer registerImageFilter = RegisterFilterType::New();

  std::vector<double> minStepLength;
  minStepLength.push_back( (double)minimumStepSize );

  std::vector<std::string> transformTypes;
  transformTypes.push_back("ScaleVersor3D");

  std::vector<int> iterations;
  iterations.push_back(numberOfIterations);

  registerImageFilter->SetTranslationScale( spatialScale );
  registerImageFilter->SetMaximumStepLength( maximumStepSize );
  registerImageFilter->SetMinimumStepLength(minStepLength  );
  registerImageFilter->SetRelaxationFactor( relaxationFactor );
  registerImageFilter->SetNumberOfIterations( iterations);
  if(numberOfSamples > 0)
    {
    const unsigned long numberOfAllSamples = extractFixedVolume->GetBufferedRegion().GetNumberOfPixels();
    samplingPercentage = static_cast<double>( numberOfSamples )/numberOfAllSamples;
    std::cout << "WARNING --numberOfSamples is deprecated, please use --samplingPercentage instead " << std::endl;
    std::cout << "WARNING: Replacing command line --samplingPercentage " << samplingPercentage << std::endl;
    }
  registerImageFilter->SetSamplePercentage( samplingPercentage );
  registerImageFilter->SetFixedVolume( anatomicalReader->GetOutput() );
  registerImageFilter->SetMovingVolume( brainOnlyFilter->GetOutput() );
  registerImageFilter->SetTransformType(transformTypes);
  // registerImageFilter->SetInitialRotationAngle( initialRotationAngle );
  // registerImageFilter->SetInitialRotationAxis( initialRotationAxis );
  try
    {
    // registerImageFilter->Update( );
    registerImageFilter->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }
  GenericTransformType::Pointer versor3DTransform = registerImageFilter->GetCurrentGenericTransform();

  itk::WriteTransformToDisk<double>(versor3DTransform, outputRigidTransform);
}
