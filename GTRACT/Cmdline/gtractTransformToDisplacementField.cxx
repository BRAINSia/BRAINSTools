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

#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkImage.h>
#include <itkAffineTransform.h>
#include <itkTransformFileReader.h>
#include <itkTransformFactory.h>
#include <itkExceptionObject.h>
#include <itkBSplineDeformableTransform.h>
#include <itkThinPlateR2LogRSplineKernelTransform.h>
#include <itkIdentityTransform.h>
#include <itkVersorRigid3DTransform.h>
// #include <itkInverseDisplacementFieldImageFilter.h>

#include "TransformToDisplacementField.h"

#include "gtractTransformToDisplacementFieldCLP.h"
#include "BRAINSThreadControl.h"
#include "GenericTransformImage.h"

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  std::cout << "Input Transform: " <<  inputTransform << std::endl;
  std::cout << "Reference Image: " <<  inputReferenceVolume << std::endl;
  std::cout << "Output Displacement Field: " << outputDeformationFieldVolume << std::endl;

  typedef itk::Vector<float, 3>                    DeformationPixelType;
  typedef itk::Image<DeformationPixelType, 3>      DisplacementFieldType;
  typedef DisplacementFieldType                    ImageType;
  typedef itk::ImageFileWriter<ImageType>          WriterType;
  typedef itk::Image<signed short, 3>              ReferenceImageType;
  typedef itk::ImageFileReader<ReferenceImageType> ReferenceReaderType;

  typedef itk::ThinPlateR2LogRSplineKernelTransform<double, 3> ThinPlateSplineTransformType;

#if (ITK_VERSION_MAJOR < 4)
  typedef itk::TransformToDeformationFieldSource<DisplacementFieldType, double> DisplacementFieldGeneratorType;
#else
  typedef itk::TransformToDisplacementFieldSource<DisplacementFieldType, double> DisplacementFieldGeneratorType;
#endif
  typedef DisplacementFieldGeneratorType::TransformType TransformType;

  // Read the Reference file
  ReferenceReaderType::Pointer referenceReader = ReferenceReaderType::New();
  referenceReader->SetFileName(inputReferenceVolume);

  try
    {
    referenceReader->Update();
    }
  catch( ... )
    {
    std::cout << "Error while reading Reference file." << std::endl;
    throw;
    }

  ReferenceImageType::Pointer image = referenceReader->GetOutput();

  // Read the transform
  GenericTransformType::Pointer baseTransform = itk::ReadTransformFromDisk(inputTransform);

  DisplacementFieldGeneratorType::Pointer defGenerator = DisplacementFieldGeneratorType::New();
  defGenerator->SetOutputParametersFromImage( image );
  defGenerator->SetTransform(baseTransform);
  try
    {
    defGenerator->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  // Write out Displacement field

  WriterType::Pointer writer = WriterType::New();
  writer->UseCompressionOn();
  writer->SetFileName( outputDeformationFieldVolume );
  writer->SetInput(defGenerator->GetOutput() );
  try
    {
    writer->Update();
    std::cout << "step 2 " << std::endl;
    }
  catch( ... )
    {
    std::cout << "Error while writing deformation file." << std::endl;
    // std::cerr<< vcl_exp <<std::endl;
    throw;
    }
  return EXIT_SUCCESS;
}
