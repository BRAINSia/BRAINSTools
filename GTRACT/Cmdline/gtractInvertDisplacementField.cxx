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

#include "itkImage.h"
#include "itkVector.h"
#include "itkGtractInverseDisplacementFieldImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "gtractInvertDisplacementFieldCLP.h"
#include "BRAINSThreadControl.h"
#include <BRAINSCommonLib.h>

int main( int argc, char *argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  const unsigned int Dimension = 3;
  typedef   signed short                            ScalarPixelType;
  typedef   itk::Image<ScalarPixelType,  Dimension> ScalarImageType;

  // Read the image defining the space for the inverse field
  typedef itk::ImageFileReader<ScalarImageType> ScalarReaderType;

  ScalarReaderType::Pointer scalarReader = ScalarReaderType::New();
  scalarReader->SetFileName( baseImage );
  try
    {
    scalarReader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown by scalar image reader" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  // Read the deformationfield field
  typedef   float                                       VectorComponentType;
  typedef   itk::Vector<VectorComponentType, Dimension> VectorType;

  typedef itk::Image<VectorType,  Dimension> DisplacementFieldType;

  typedef itk::ImageFileReader<DisplacementFieldType> VectorReaderType;

  VectorReaderType::Pointer vectorReader = VectorReaderType::New();
  vectorReader->SetFileName( deformationImage );
  try
    {
    vectorReader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown by vector image reader" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  // Invert the deformationfield field
  typedef itk::GtractInverseDisplacementFieldImageFilter<
      DisplacementFieldType,
      DisplacementFieldType
      >  FilterType;

  FilterType::Pointer inverseFilter = FilterType::New();
  inverseFilter->SetOutputSpacing( scalarReader->GetOutput()->GetSpacing() );
  inverseFilter->SetOutputOrigin( scalarReader->GetOutput()->GetOrigin() );
  ScalarImageType::RegionType region = scalarReader->GetOutput()->GetLargestPossibleRegion();

  inverseFilter->SetSize( region.GetSize() );
  inverseFilter->SetInput( vectorReader->GetOutput() );
  inverseFilter->SetSubsamplingFactor( subsamplingFactor );

  try
    {
    inverseFilter->UpdateLargestPossibleRegion();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown in Inverse filter" << std::endl;
    std::cerr << excp << std::endl;
    }

  // Write an image for regression testing
  typedef itk::ImageFileWriter<DisplacementFieldType> WriterType;

  WriterType::Pointer writer = WriterType::New();
  writer->UseCompressionOn();
  writer->SetInput( inverseFilter->GetOutput() );
  writer->SetFileName( outputVolume );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown by writer" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
