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
#include <itkResampleImageFilter.h>
#include <itkTransformFactory.h>
#include <itkVersorRigid3DTransform.h>
#include <itkOrientImageFilter.h>

#include "itkGtractImageIO.h"
#include "itkInvertBSplineFilter.h"
#include "gtractInvertBSplineTransformCLP.h"
#include "BRAINSThreadControl.h"
#include "GenericTransformImage.h"

int main(int argc, char *argv[])
{
  PARSE_ARGS;
  BRAINSRegisterAlternateIO();
  const BRAINSUtils::StackPushITKDefaultNumberOfThreads TempDefaultNumberOfThreadsHolder(numberOfThreads);

  int xSize = landmarkDensity[0];
  int ySize = landmarkDensity[1];
  int zSize = landmarkDensity[2];

  bool debug = true;
  if( debug )
    {
    std::cout << "==============================================================" << std::endl;
    std::cout << "Input Transform: " <<  inputTransform << std::endl;
    std::cout << "Output Transform: " <<  outputTransform << std::endl;
    std::cout << "Reference Image: " <<  inputReferenceVolume << std::endl;
    std::cout << "Density: " << xSize << " " << ySize << " " << zSize << std::endl;
    std::cout << "==============================================================" << std::endl;
    }

  typedef double                              BSplineCoordinateRepType;
  typedef itk::VersorRigid3DTransform<double> RigidTransformType;
  typedef itk::ThinPlateR2LogRSplineKernelTransform<
      BSplineCoordinateRepType, 3>     ThinPlateSplineTransformType;


  typedef signed short PixelType;

  typedef itk::Image<PixelType, 3>        ImageType;
  typedef itk::ImageFileReader<ImageType> AnatomicalImageReaderType;
  AnatomicalImageReaderType::Pointer anatomicalReader = AnatomicalImageReaderType::New();
  anatomicalReader->SetFileName( inputReferenceVolume );
  try
    {
    anatomicalReader->Update();
    }
  catch( itk::ExceptionObject & ex )
    {
    std::cout << ex << std::endl;
    throw;
    }

  ImageType::Pointer ExampleImage;
  ExampleImage = anatomicalReader->GetOutput();

  std::cout << "Using Image: " << ExampleImage << std::endl;

  // Read the transform
  typedef itk::Transform<double, 3, 3> GenericTransformType;
  GenericTransformType::Pointer baseTransform = itk::ReadTransformFromDisk(inputTransform);

  typedef itk::InvertBSplineFilter InvertFilterType;
  std::cout << "Running Inversion using TPS" << std::endl;
  InvertFilterType::Pointer invertTransformFilter = InvertFilterType::New();
    {
    InvertFilterType::BsplineTransformTypePointer myBSpline
      = dynamic_cast<InvertFilterType::BsplineTransformType *>( baseTransform.GetPointer() );
    invertTransformFilter->SetInput( myBSpline );
    }
  invertTransformFilter->SetXgridSize( xSize );
  invertTransformFilter->SetYgridSize( ySize );
  invertTransformFilter->SetZgridSize( zSize );
  invertTransformFilter->SetExampleImage( ExampleImage );
  invertTransformFilter->Update();
  std::cout << "TPS Inversion Complete" << std::endl;

  itk::WriteTransformToDisk<double>(invertTransformFilter->GetOutput(), outputTransform);
  return EXIT_SUCCESS;
}
