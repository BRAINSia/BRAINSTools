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

#include "itkVector.h"
#include "itkQuadEdgeMesh.h"
#include "itkTriangleBasisSystem.h"
#include "itkTriangleBasisSystemCalculator.h"
#include "itkTriangleCell.h"
#include "itkTestingMacros.h"

int main( int, char * [] )
{
  const unsigned int Dimension = 3;
  const unsigned int SurfaceDimension = 2;

  typedef itk::QuadEdgeMesh<double, Dimension> MeshType;

  MeshType::Pointer mesh = MeshType::New();

  typedef MeshType::PointType   PointType;
  typedef PointType::VectorType VectorType;

  typedef itk::TriangleBasisSystem<VectorType, SurfaceDimension> TriangleBasisSystemType;
  TriangleBasisSystemType triangleBasisSystem;

  typedef itk::TriangleBasisSystemCalculator<MeshType, TriangleBasisSystemType> TriangleBasisSystemCalculatorType;
  TriangleBasisSystemCalculatorType::Pointer triangleBasisSystemCalculator = TriangleBasisSystemCalculatorType::New();

  std::cout << triangleBasisSystemCalculator->GetNameOfClass() << std::endl;
  triangleBasisSystemCalculator->Print( std::cout );

  // Define a simple triangular cell
  typedef MeshType::PointType PointType;
  PointType p0;
  PointType p1;
  PointType p2;

  p0[0] = 0.0; p0[1] = 0.0; p0[2] = 0.0;
  p1[0] = 1.0; p1[1] = 0.0; p1[2] = 0.0;
  p2[0] = 0.5; p2[1] = 1.0; p2[2] = 0.0;

  mesh->SetPoint( 0, p0 );
  mesh->SetPoint( 1, p1 );
  mesh->SetPoint( 2, p2 );

  TRY_EXPECT_EXCEPTION( triangleBasisSystemCalculator->CalculateTriangle( 0, triangleBasisSystem ) );

  TRY_EXPECT_NO_EXCEPTION( triangleBasisSystemCalculator->SetInputMesh( mesh ) );

  TEST_SET_GET( mesh, triangleBasisSystemCalculator->GetInputMesh() );

  TRY_EXPECT_EXCEPTION( triangleBasisSystemCalculator->CalculateTriangle( 0, triangleBasisSystem ) );

  typedef MeshType::CellType          CellType;
  typedef itk::TriangleCell<CellType> TriangleType;

  CellType::CellAutoPointer cellpointer;
  cellpointer.TakeOwnership( new TriangleType );
  cellpointer->SetPointId( 0, 0 );
  cellpointer->SetPointId( 1, 1 );
  cellpointer->SetPointId( 2, 2 );

  mesh->SetCell( 0, cellpointer );

  TRY_EXPECT_NO_EXCEPTION( triangleBasisSystemCalculator->CalculateTriangle( 0, triangleBasisSystem ) );
  TRY_EXPECT_EXCEPTION( triangleBasisSystemCalculator->CalculateTriangle( 1, triangleBasisSystem ) );

  VectorType expectedV01;
  VectorType expectedV02;

  std::cout << "Test passed." << std::endl;

  return EXIT_SUCCESS;
}
