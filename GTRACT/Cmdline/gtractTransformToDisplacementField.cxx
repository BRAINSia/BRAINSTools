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

Purpose:   Get Displacement Field from the input Transform
Date:      11/13/07
Author:    Madhura A Ingalhalikar


Purpose:  Replace TransformToDisplacementField with GTRACTTransformToDisplacementField
Date:     04/26/10
Author:   Yongqiang Zhao

Copyright (c) University of Iowa Department of Radiology. All rights reserved.
See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include "itkIO.h"

#include "TransformToDisplacementField.h"

#include "gtractTransformToDisplacementFieldCLP.h"
#include "BRAINSThreadControl.h"
#include "GenericTransformImage.h"

int
main( int argc, char * argv[] )
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder( numberOfThreads );

  std::cout << "Input Transform: " << inputTransform << std::endl;
  std::cout << "Reference Image: " << inputReferenceVolume << std::endl;
  std::cout << "Output Displacement Field: " << outputDeformationFieldVolume << std::endl;

  using DeformationPixelType = itk::Vector< float, 3 >;
  using DisplacementFieldType = itk::Image< DeformationPixelType, 3 >;
  using ImageType = DisplacementFieldType;
  using WriterType = itk::ImageFileWriter< ImageType >;
  using ReferenceImageType = itk::Image< signed short, 3 >;

  ReferenceImageType::Pointer image;
  try
  {
    image = itkUtil::ReadImage< ReferenceImageType >( inputReferenceVolume );
  }
  catch ( ... )
  {
    std::cout << "Error while reading Reference file." << std::endl;
    throw;
  }

  // Read the transform
  using GenericTransformType = itk::Transform< double, 3, 3 >;
  GenericTransformType::Pointer baseTransform = itk::ReadTransformFromDisk( inputTransform );

  DisplacementFieldType::Pointer displacementField =
    TransformToDisplacementField< ImageType::Pointer, GenericTransformType::Pointer >( image, baseTransform );

  // Write out Displacement field

  WriterType::Pointer writer = WriterType::New();
  writer->UseCompressionOn();
  writer->SetFileName( outputDeformationFieldVolume );
  writer->SetInput( displacementField );
  try
  {
    writer->Update();
    std::cout << "step 2 " << std::endl;
  }
  catch ( ... )
  {
    std::cout << "Error while writing deformation file." << std::endl;
    // std::cerr<< std::exp <<std::endl;
    throw;
  }
  return EXIT_SUCCESS;
}
