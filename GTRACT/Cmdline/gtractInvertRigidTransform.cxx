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
#include <itkResampleImageFilter.h>
#include <itkBSplineDeformableTransform.h>
#include <itkTransformFactory.h>
#include <itkVersorRigid3DTransform.h>
#include <itkOrientImageFilter.h>

#include "itkGtractImageIO.h"
#include "itkInvertBSplineFilter.h"
#include "gtractInvertRigidTransformCLP.h"
#include "GenericTransformImage.h"
#include "BRAINSThreadControl.h"
int main(int argc, char *argv[])
{
  PARSE_ARGS;
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);
  itk::AddExtraTransformRegister();

  bool debug = true;
  if( debug )
    {
    std::cout << "==============================================================" << std::endl;
    std::cout << "Input Transform: " <<  inputTransform << std::endl;
    std::cout << "Output Transform: " <<  outputTransform << std::endl;
    std::cout << "==============================================================" << std::endl;
    }

  typedef itk::VersorRigid3DTransform<double> RigidTransformType;
  // Read the transform
  GenericTransformType::Pointer forwardTransform = itk::ReadTransformFromDisk(inputTransform);
  RigidTransformType::Pointer   reverseTransform = RigidTransformType::New();
  forwardTransform->GetInverse(reverseTransform);
  itk::WriteTransformToDisk(reverseTransform, outputTransform);
  return EXIT_SUCCESS;
}
