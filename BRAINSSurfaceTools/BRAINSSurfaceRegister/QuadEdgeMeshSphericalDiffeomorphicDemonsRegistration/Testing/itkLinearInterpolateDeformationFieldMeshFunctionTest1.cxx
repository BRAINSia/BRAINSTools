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
#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkLinearInterpolateDeformationFieldMeshFunction.h"

#include <iostream>

int main(int argc, char* argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: ";
    std::cerr << argv[0] << " ";
    std::cerr << "inputMeshFilename ";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::QuadEdgeMesh<float, 3>                  MeshType;
  typedef itk::QuadEdgeMeshVTKPolyDataReader<MeshType> ReaderType;

  ReaderType::Pointer surfaceReader = ReaderType::New();

  typedef MeshType::PointType PointType;

  typedef MeshType::PointsContainer InputPointsContainer;

  surfaceReader->SetFileName( argv[1] );

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

  std::cout << "Verifying Reader" << std::endl;

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
  typedef itk::LinearInterpolateDeformationFieldMeshFunction<MeshType, InputPointsContainer> InterpolatorType;

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

  InputPointsContainer::Pointer inputPoints = mesh->GetPoints();

  PointType point;

  mesh->GetPoint( 0, &point );

  InterpolatorType::PointType pointToEvaluate;

  pointToEvaluate.CastFrom(point);

  InterpolatorType::PointType evaluatedPoint;

  const bool evaluate = interpolator->Evaluate( inputPoints, point, evaluatedPoint );

  if( !evaluate )
    {
    std::cout << "Couldn't evaluate point " << point << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Value at " << point << " = " << evaluatedPoint << std::endl;

  // PointType point1;
  // point1[0] = 100.0;
  // point1[1] =   0.0;
  // point1[2] =  10.7076;

  // typedef InterpolatorSuperclassType::InstanceIdentifierVectorType  InstanceIdentifierVectorType;

  // InstanceIdentifierVectorType pointIds(3);

  // const bool foundTriangle1 = interpolator->FindTriangle( point1, pointIds );

  // if( !foundTriangle1 )
  //  {
  //  std::cout << "Couldn't find triangle for point " << point1 << std::endl;
  //  return EXIT_FAILURE;
  //  }

  // PointType point2;
  // point2[0] =   0.0830543;
  // point2[1] =   0.0;
  // point2[2] = 100.0;

  // const bool foundTriangle2 = interpolator->FindTriangle( point2, pointIds );

  // if( !foundTriangle2 )
  //  {
  //  std::cout << "Couldn't find triangle for point " << point2 << std::endl;
  //  return EXIT_FAILURE;
  //  }

  std::cout << "Test PASSED !" << std::endl;

  return EXIT_SUCCESS;
}
