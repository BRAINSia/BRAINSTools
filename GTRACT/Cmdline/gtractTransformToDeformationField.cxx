/*=========================================================================

Program:   GTRACT (Guided Tensor Restore Anatomical Connectivity Tractography)
Module:    $RCSfile: $
Language:  C++
Date:      $Date: 2006/03/29 14:53:40 $
Version:   $Revision: 1.9 $

Purpose:   Get Deformation Field from the input Transform
Date:      11/13/07
Author:    Madhura A Ingalhalikar


Purpose:  Replace itkTransformToDeformationFieldSource with GTRACTTransformToDeformationField
Date:     04/26/10
Author:   Yongqiang Zhao

Copyright (c) University of Iowa Department of Radiology. All rights reserved.
See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

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
#include <itkInverseDeformationFieldImageFilter.h>

#include "itkTransformToDeformationFieldSource.h"

#include "gtractTransformToDeformationFieldCLP.h"
#include "BRAINSThreadControl.h"
#include "GenericTransformImage.h"

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSUtils::SetThreadCount(numberOfThreads);
  itk::AddExtraTransformRegister();

  std::cout << "Input Transform: " <<  inputTransform << std::endl;
  std::cout << "Reference Image: " <<  inputReferenceVolume << std::endl;
  std::cout << "Output Deformation Field: " << outputDeformationFieldVolume << std::endl;

  typedef itk::Vector<float, 3>                    DeformationPixelType;
  typedef itk::Image<DeformationPixelType, 3>      DeformationFieldType;
  typedef DeformationFieldType                     ImageType;
  typedef itk::ImageFileWriter<ImageType>          WriterType;
  typedef itk::Image<signed short, 3>              ReferenceImageType;
  typedef itk::ImageFileReader<ReferenceImageType> ReferenceReaderType;

  // TODO:  May need to add the TPS transform type to
  // "AddExtraTransformRegister"
  typedef itk::ThinPlateR2LogRSplineKernelTransform<double, 3> ThinPlateSplineTransformType;

  typedef itk::TransformToDeformationFieldSource<DeformationFieldType, double> DeformationFieldGeneratorType;
  typedef DeformationFieldGeneratorType::TransformType                         TransformType;

  itk::AddExtraTransformRegister();

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

  DeformationFieldGeneratorType::Pointer defGenerator = DeformationFieldGeneratorType::New();
#if 0
  defGenerator->SetOutputSize( image->GetLargestPossibleRegion().GetSize() );
  defGenerator->SetOutputSpacing( image->GetSpacing() );
  defGenerator->SetOutputOrigin( image->GetOrigin() );
  defGenerator->SetOutputIndex( image->GetLargestPossibleRegion().GetIndex() );
  defGenerator->SetOutputDirection( image->GetDirection() );
#else
  defGenerator->SetOutputParametersFromImage( image );
#endif
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

  // Write out Deformation field

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
