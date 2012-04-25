/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTriangleBasisSystemTest.cxx,v $
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
#include "itkTestingMacros.h"

int main( int, char * [] )
{
  const unsigned int Dimension = 3;
  const unsigned int SurfaceDimension = 2;

  typedef itk::QuadEdgeMesh<double, Dimension> MeshType;

  typedef MeshType::PointType   PointType;
  typedef PointType::VectorType VectorType;

  typedef itk::TriangleBasisSystem<VectorType, SurfaceDimension> TriangleBasisSystemType;

  TriangleBasisSystemType triangleBasisSystem;

  /*
   * Define a simple triangular cell
   *
   *     ^      p2
   *     |     /  \
   *     |    /    \
   *     |   /      \
   *     |  /        \
   *     | /          \
   *     |/            \
   *    p0--------------p1-->
   *
   */

  PointType p0;
  PointType p1;
  PointType p2;

  p0[0] = 0.0; p0[1] = 0.0; p0[2] = 0.0;
  p1[0] = 1.0; p1[1] = 0.0; p1[2] = 0.0;
  p2[0] = 0.5; p2[1] = 1.0; p2[2] = 0.0;

  VectorType v01 = p1 - p0;
  VectorType v02 = p2 - p0;

  // Should be SurfaceDimension-1 at most
  TRY_EXPECT_EXCEPTION( triangleBasisSystem.SetVector(SurfaceDimension, v01) );

  // Should be Dimension-1 at most
  TRY_EXPECT_NO_EXCEPTION( triangleBasisSystem.SetVector(0, v01) );

  // Should be Dimension-1 at most
  TRY_EXPECT_NO_EXCEPTION( triangleBasisSystem.SetVector(1, v02) );

  // Should be Dimension-1 at most
  TRY_EXPECT_EXCEPTION( triangleBasisSystem.GetVector(SurfaceDimension) );
  TRY_EXPECT_EXCEPTION( triangleBasisSystem.GetVector(SurfaceDimension + 5) );
  TRY_EXPECT_NO_EXCEPTION( triangleBasisSystem.GetVector(SurfaceDimension - 1) );

  VectorType vb01 = triangleBasisSystem.GetVector(0);
  VectorType vb02 = triangleBasisSystem.GetVector(1);

  std::cout << std::endl;
  std::cout << "Basis vector 0 " << vb01 << std::endl;
  std::cout << "Basis vector 1 " << vb02 << std::endl;
  std::cout << std::endl;

  const double tolerance = 1e-5;
  for( unsigned int k = 0; k < Dimension; k++ )
    {
    if( vnl_math_abs( vb01[k] - v01[k] ) > tolerance )
      {
      std::cerr << "Error, SetVector()/GetVector() failed" << std::endl;
      return EXIT_FAILURE;
      }

    if( vnl_math_abs( vb02[k] - v02[k] ) > tolerance )
      {
      std::cerr << "Error, SetVector()/GetVector() failed" << std::endl;
      return EXIT_FAILURE;
      }
    }

  TriangleBasisSystemType system1( triangleBasisSystem );

  TriangleBasisSystemType system2 = triangleBasisSystem;

  TriangleBasisSystemType system3;

  system3 = triangleBasisSystem;

  VectorType diff1 = system3.GetVector(0) - triangleBasisSystem.GetVector(0);
  VectorType diff2 = system3.GetVector(1) - triangleBasisSystem.GetVector(1);

  if( diff1.GetNorm() > vnl_math::eps || diff2.GetNorm() > vnl_math::eps )
    {
    std::cerr << "Error in operator=() " << std::endl;
    return EXIT_FAILURE;
    }

  TriangleBasisSystemType system4;

  system4 = system3 = system2;

  std::cout << "Test passed." << std::endl;

  return EXIT_SUCCESS;
}
