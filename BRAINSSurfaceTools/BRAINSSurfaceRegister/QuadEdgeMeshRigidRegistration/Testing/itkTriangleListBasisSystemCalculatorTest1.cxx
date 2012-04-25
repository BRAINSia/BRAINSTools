/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTriangleBasisSystemCalculatorTest.cxx,v $
  Language:  C++
  Date:      $Date: 2005-07-26 15:55:12 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkPoint.h"
#include "itkVector.h"
#include "itkQuadEdgeMesh.h"
#include "itkTriangleBasisSystem.h"
#include "itkTriangleBasisSystemCalculator.h"
#include "itkTriangleListBasisSystemCalculator.h"
#include "itkRegularSphereMeshSource.h"
#include "itkMeshGeneratorHelper.h"
#include "itkTestingMacros.h"

int main( int, char * [] )
{
  const unsigned int SurfaceDimension = 2;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<float, Dimension> MovingMeshType;
  typedef itk::QuadEdgeMesh<float, Dimension> MeshType;

  MeshType::Pointer       mesh;
  MovingMeshType::Pointer movingMesh;
  typedef MeshType::PointType   PointType;
  typedef PointType::VectorType VectorType;

  typedef MeshGeneratorHelper<MeshType, MovingMeshType> GeneratorType;

  GeneratorType::GenerateMeshes( mesh, movingMesh );

  typedef itk::TriangleBasisSystem<VectorType, SurfaceDimension> TriangleBasisSystemType;
  TriangleBasisSystemType triangleBasisSystem;

  typedef itk::TriangleListBasisSystemCalculator<MeshType, TriangleBasisSystemType>
    TriangleListBasisSystemCalculatorType;

  TriangleListBasisSystemCalculatorType::Pointer triangleListBasisSystemCalculator =
    TriangleListBasisSystemCalculatorType::New();

  std::cout << triangleListBasisSystemCalculator->GetNameOfClass() << std::endl;
  triangleListBasisSystemCalculator->Print( std::cout );

  std::cout << "test 1 \n";

  // Mesh not set yet. Exception!
  TRY_EXPECT_EXCEPTION( triangleListBasisSystemCalculator->Calculate() );

  std::cout << "test 2 \n";

  TRY_EXPECT_NO_EXCEPTION( triangleListBasisSystemCalculator->SetInputMesh( mesh ) );

  std::cout << "test 3 \n";

  // List should be allocated with constructor now. No exception.
  TRY_EXPECT_NO_EXCEPTION( triangleListBasisSystemCalculator->GetBasisSystemList() );

  std::cout << "test 4 \n";

  TEST_SET_GET( mesh, triangleListBasisSystemCalculator->GetInputMesh() );

  // This should produce a list
  TRY_EXPECT_NO_EXCEPTION( triangleListBasisSystemCalculator->Calculate() );

  std::cout << "after test 4 \n";

  const TriangleListBasisSystemCalculatorType::BasisSystemListType * basisSystemList =
    triangleListBasisSystemCalculator->GetBasisSystemList();

  TriangleListBasisSystemCalculatorType::BasisSystemListIterator basisSystemListIterator;

  std::cout << "Number of elements in the basis system list " << basisSystemList->Size() << std::endl;

  basisSystemListIterator = basisSystemList->Begin();

  std::cout << " basis list first element " << std::endl;

  typedef TriangleListBasisSystemCalculatorType::BasisSystemType BasisSystemType;

  VectorType reference;
  reference[0] = 17.0;
  reference[1] = 19.0;
  reference[2] = 21.0;

  std::cout << "Reference = " << reference << std::endl;

  BasisSystemType basis = basisSystemListIterator->Value();
  VectorType      vector0 = basis.GetVector(0);
  VectorType      vector1 = basis.GetVector(1);

  std::cout << vector0 << std::endl;
  std::cout << vector1 << std::endl;

  std::cout << "Test passed." << std::endl;

  return EXIT_SUCCESS;
}
