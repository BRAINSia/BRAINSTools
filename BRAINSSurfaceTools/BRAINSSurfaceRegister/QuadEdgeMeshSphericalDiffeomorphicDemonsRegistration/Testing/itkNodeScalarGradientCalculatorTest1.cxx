/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNodeScalarGradientCalculatorTest1.cxx,v $
  Language:  C++
  Date:      $Date: 2008-03-10 19:46:31 $
  Version:   $Revision: 1.37 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkNodeScalarGradientCalculator.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkMeshGeneratorHelper.h"
#include "itkCovariantVector.h"
#include "itkTestingMacros.h"
#include "itkTriangleBasisSystem.h"
#include "itkTriangleBasisSystemCalculator.h"
#include "itkTriangleListBasisSystemCalculator.h"

int main( int argc, char *argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedMeshFile  movingMeshFile " << std::endl;
    return EXIT_FAILURE;
    }

  typedef float MeshPixelType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> FixedMeshType;
  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> MovingMeshType;

  typedef MovingMeshType::PointDataContainer MovingPointDataContainerType;

  typedef itk::NodeScalarGradientCalculator<
      FixedMeshType, MovingPointDataContainerType>   GradientCalculatorType;

  GradientCalculatorType::Pointer gradientCalculator = GradientCalculatorType::New();

  std::cout << gradientCalculator->GetNameOfClass() << std::endl;
  gradientCalculator->Print( std::cout );

  typedef itk::QuadEdgeMeshVTKPolyDataReader<FixedMeshType>  FixedReaderType;
  typedef itk::QuadEdgeMeshVTKPolyDataReader<MovingMeshType> MovingReaderType;

  FixedReaderType::Pointer fixedReader = FixedReaderType::New();
  fixedReader->SetFileName( argv[1] );

  MovingReaderType::Pointer movingReader = MovingReaderType::New();
  movingReader->SetFileName( argv[2] );

  try
    {
    fixedReader->Update();
    movingReader->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while reading the input file " << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "test 4 \n";

  GradientCalculatorType::PointType center1;
  center1.Fill( 17.0 );
  gradientCalculator->SetSphereCenter( center1 );
  TEST_SET_GET_VALUE( center1, gradientCalculator->GetSphereCenter() );

  GradientCalculatorType::PointType center2;
  center2.Fill( 31.0 );
  gradientCalculator->SetSphereCenter( center2 );
  TEST_SET_GET_VALUE( center2, gradientCalculator->GetSphereCenter() );

  const double radius1 = 19.0;
  gradientCalculator->SetSphereRadius( radius1 );
  TEST_SET_GET_VALUE( radius1, gradientCalculator->GetSphereRadius() );

  const double radius2 = 29.0;
  gradientCalculator->SetSphereRadius( radius2 );
  TEST_SET_GET_VALUE( radius2, gradientCalculator->GetSphereRadius() );

  GradientCalculatorType::PointType center0;
  center0.Fill( 0.0 );

  gradientCalculator->SetSphereCenter( center0 );
  gradientCalculator->SetSphereRadius( 100.0 );

  // Have not properly initialized gradientCalculator yet...
  TRY_EXPECT_EXCEPTION( gradientCalculator->Initialize(); );

  // First item needed: triangle basis list.
  const unsigned int SurfaceDimension = 2;
  typedef FixedMeshType::PointType                               PointType;
  typedef PointType::VectorType                                  VectorType;
  typedef itk::TriangleBasisSystem<VectorType, SurfaceDimension> TriangleBasisSystemType;
  TriangleBasisSystemType triangleBasisSystem;

  typedef itk::TriangleListBasisSystemCalculator<FixedMeshType, TriangleBasisSystemType>
    TriangleListBasisSystemCalculatorType;

  TriangleListBasisSystemCalculatorType::Pointer triangleListBasisSystemCalculator =
    TriangleListBasisSystemCalculatorType::New();

  triangleListBasisSystemCalculator->SetInputMesh( fixedReader->GetOutput() );

  // Should be able to compute basis list with fixed mesh.
  TRY_EXPECT_NO_EXCEPTION( triangleListBasisSystemCalculator->Calculate() );

  std::cout << "test 3 \n";

  // Still have not properly initialized gradientCalculator yet... 2 to go
  TRY_EXPECT_EXCEPTION( gradientCalculator->Initialize(); );

  gradientCalculator->SetInputMesh( triangleListBasisSystemCalculator->GetInputMesh() );

  if( gradientCalculator->GetInputMesh() != triangleListBasisSystemCalculator->GetInputMesh() )
    {
    std::cerr << "Error in SetFixedMesh()/GetFixedMesh() " << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "test 2 \n";

  // Still have not properly initialized gradientCalculator yet... 1 to go
  TRY_EXPECT_EXCEPTION( gradientCalculator->Initialize(); );

  MovingPointDataContainerType::Pointer movingData = MovingPointDataContainerType::New();
  movingData->Reserve( fixedReader->GetOutput()->GetNumberOfPoints() );

  MovingPointDataContainerType::Iterator dataItr = movingData->Begin();
  MovingPointDataContainerType::Iterator dataEnd = movingData->End();

  while( dataItr != dataEnd )
    {
    dataItr.Value() = 1.0;
    ++dataItr;
    }

  std::cout << "GetPointData from moving mesh " << movingData.GetPointer() << std::endl;

  gradientCalculator->SetDataContainer( movingData );

  if( gradientCalculator->GetDataContainer() != movingData )
    {
    std::cerr << "Error in SetDataContainer()/GetDataContainer() " << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "test 1 gradientCalculator->GetDataContainer() "
            << gradientCalculator->GetDataContainer() << "\n";

  // Have still not properly initialized gradientCalculator yet...
  TRY_EXPECT_EXCEPTION( gradientCalculator->Initialize(); );

  gradientCalculator->SetBasisSystemList( triangleListBasisSystemCalculator->GetBasisSystemList() );

  if( gradientCalculator->GetBasisSystemList() != triangleListBasisSystemCalculator->GetBasisSystemList() )
    {
    std::cerr << "Error in SetBasisSystemList()/GetBasisSystemList() " << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "test 1 gradientCalculator->GetDataContainer() "
            << gradientCalculator->GetDataContainer() << "\n";

  // It is initialized correctly, we expect no exception
  TRY_EXPECT_NO_EXCEPTION( gradientCalculator->Initialize(); );

  // Now that Initialize() has been called successfully, we can call Compute(). */
  TRY_EXPECT_NO_EXCEPTION( gradientCalculator->Compute(); );

  try
    {
    gradientCalculator->Evaluate( 17 );
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
