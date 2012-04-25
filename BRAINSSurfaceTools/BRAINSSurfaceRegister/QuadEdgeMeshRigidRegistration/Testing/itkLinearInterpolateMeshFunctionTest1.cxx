/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFreeSurferBinarySurfaceReaderQuadEdgeMeshTest.cxx,v $
  Language:  C++
  Date:      $Date: 2008-06-16 02:07:17 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "itkQuadEdgeMesh.h"
#include "itkFreeSurferBinarySurfaceReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkLinearInterpolateMeshFunction.h"

#include <iostream>

int main(int argc, char* argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage: itkFreeSurferBinarySurfaceReaderTest inputFilename inputDataFilename";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::QuadEdgeMesh<float, 3>                  MeshType;
  typedef itk::FreeSurferBinarySurfaceReader<MeshType> ReaderType;

  ReaderType::Pointer surfaceReader = ReaderType::New();

  typedef ReaderType::PointType  PointType;
  typedef ReaderType::VectorType VectorType;

  surfaceReader->SetFileName(argv[1]);
  surfaceReader->SetDataFileName(argv[2]);

  try
    {
    surfaceReader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error during Update() " << std::endl;
    std::cerr << excp << std::endl;
    }

  MeshType::Pointer mesh = surfaceReader->GetOutput();

  std::cout << "Testing itk::FreeSurferBinarySurfaceReader" << std::endl;

  unsigned int numberOfPoints = mesh->GetNumberOfPoints();
  unsigned int numberOfCells  = mesh->GetNumberOfCells();

  std::cout << "numberOfPoints= " << numberOfPoints << std::endl;
  std::cout << "numberOfCells= " << numberOfCells << std::endl;

  if( !numberOfPoints )
    {
    std::cerr << "ERROR: numberOfPoints= " << numberOfPoints << std::endl;
    return EXIT_FAILURE;
    }

  if( !numberOfCells )
    {
    std::cerr << "ERROR: numberOfCells= " << numberOfCells << std::endl;
    return EXIT_FAILURE;
    }

  //  Exercise the Interpolation
  typedef itk::LinearInterpolateMeshFunction<MeshType> InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  interpolator->SetInputMesh( mesh );

  interpolator->Initialize();

  typedef InterpolatorType::Superclass InterpolatorSuperclassType;

  std::cout << interpolator->InterpolatorSuperclassType::GetNameOfClass() << std::endl;
  interpolator->InterpolatorSuperclassType::Print( std::cout );

  typedef InterpolatorSuperclassType::Superclass InterpolatorSuperSuperclassType;

  std::cout << interpolator->InterpolatorSuperSuperclassType::GetNameOfClass() << std::endl;
  interpolator->InterpolatorSuperSuperclassType::Print( std::cout );

  std::cout << "Interpolator " << std::endl;
  interpolator->Print( std::cout );

  std::cout << interpolator->GetNameOfClass() << std::endl;

  PointType point;

  mesh->GetPoint( 0, &point );

  InterpolatorType::RealType interpolatedValue = interpolator->Evaluate( point );

  std::cout << "Value at " << point << " = " << interpolatedValue << std::endl;

  std::cout << "Test passed" << std::endl;

  return EXIT_SUCCESS;
}
