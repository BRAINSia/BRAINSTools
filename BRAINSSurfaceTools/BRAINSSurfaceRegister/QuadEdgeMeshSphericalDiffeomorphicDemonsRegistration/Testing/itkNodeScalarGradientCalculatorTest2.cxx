/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNodeScalarGradientCalculatorTest2.cxx,v $
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
#include "itkMeshGeneratorHelper.h"
#include "itkCovariantVector.h"
#include "itkTestingMacros.h"
#include "itkTriangleBasisSystem.h"
#include "itkTriangleBasisSystemCalculator.h"
#include "itkTriangleListBasisSystemCalculator.h"
#include "itkMeshGeneratorDerivativeHelper.h"
#include "itkVersorTransform.h"

// Test that creates two spherical meshes, with a PI/4
// rotation applied to moving mesh, ascertains that A) gradient estimated
// at various positions of moving mesh compares numerically to B) gradient
// estimated by 1) assigning moving mesh values to fixed vertices based on
// fixed vertex positions transformed to to moving mesh, 2) estimating
// gradients at fixed faces from vertex values, and 3) averaging
// gradient values at each fixed vertex based on incident faces.
int main( int, char * [] )
{
  const unsigned int Dimension = 3;
  const double       piOver4 = atan( 1.0 );
  const double       pi = piOver4 * 4.0;
  const double       twoPi = pi * 2.0;

  typedef itk::QuadEdgeMesh<float, Dimension> MovingMeshType;
  typedef itk::QuadEdgeMesh<float, Dimension> FixedMeshType;

  FixedMeshType::Pointer  fixedMesh;
  MovingMeshType::Pointer movingMesh;

  typedef MeshGeneratorDerivativeHelper<FixedMeshType, MovingMeshType> GeneratorType;

  GeneratorType::GenerateMeshes( fixedMesh, movingMesh );

  typedef MovingMeshType::PointDataContainer MovingPointDataContainerType;

  typedef itk::NodeScalarGradientCalculator<
      FixedMeshType, MovingPointDataContainerType>   GradientCalculatorType;

  GradientCalculatorType::Pointer gradientCalculator = GradientCalculatorType::New();

  std::cout << gradientCalculator->GetNameOfClass() << std::endl;
  gradientCalculator->Print( std::cout );

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

  triangleListBasisSystemCalculator->SetInputMesh( fixedMesh );

  // Should be able to compute basis list with fixed mesh.
  TRY_EXPECT_NO_EXCEPTION( triangleListBasisSystemCalculator->Calculate() );

  gradientCalculator->SetInputMesh( triangleListBasisSystemCalculator->GetInputMesh() );

  if( gradientCalculator->GetInputMesh() != triangleListBasisSystemCalculator->GetInputMesh() )
    {
    std::cerr << "Error in SetFixedMesh()/GetFixedMesh() " << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::VersorTransform<double> TransformType;
  typedef TransformType::VersorType    VersorType;

  VersorType rotation;
  VectorType axis;
  PointType  myCenter;
  myCenter.Fill( 0.0 );

  axis[0] = 0.0;
  axis[1] = 0.0;
  axis[2] = 1.0;

  const double angle = piOver4;

  rotation.Set(  axis, angle  );

  TransformType::Pointer transform = TransformType::New();
  transform->SetRotation( rotation );

  MovingPointDataContainerType::Pointer movingData = MovingPointDataContainerType::New();
  movingData->Reserve( fixedMesh->GetNumberOfPoints() );

  MovingPointDataContainerType::Iterator dataItr = movingData->Begin();
  MovingPointDataContainerType::Iterator dataEnd = movingData->End();

  typedef FixedMeshType::PointsContainer      FixedPointsContainer;
  typedef FixedPointsContainer::ConstIterator FixedPointIterator;

  FixedPointIterator pointIterator = fixedMesh->GetPoints()->Begin();
  FixedPointIterator pointEnd = fixedMesh->GetPoints()->End();

  typedef MovingMeshType::PointType MovingPointType;
  typedef FixedMeshType::PointType  FixedPointType;
  MovingPointType mappedMovingMeshPoint;

  while( pointIterator != pointEnd )
    {
    FixedPointType fixedPoint;
    fixedPoint.CastFrom( pointIterator.Value() );

    mappedMovingMeshPoint = transform->TransformPoint( fixedPoint );

    const VectorType radial = mappedMovingMeshPoint - myCenter;
    const double     radius = radial.GetNorm(); // assuming radius is not valued 1

    const double theta = vcl_atan2(mappedMovingMeshPoint[1], mappedMovingMeshPoint[0]);
    const double phi = vcl_acos(mappedMovingMeshPoint[2] / radius);

    float thetaOffset = theta - piOver4;
    if( thetaOffset < -pi )
      {
      thetaOffset += twoPi;
      }

    dataItr.Value() = mapSphericalCoordinatesFunction(phi, thetaOffset);
    ++dataItr;
    ++pointIterator;
    }

  gradientCalculator->SetDataContainer( movingData );

  if( gradientCalculator->GetDataContainer() != movingData )
    {
    std::cerr << "Error in SetDataContainer()/GetDataContainer() " << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "test 1 gradientCalculator->GetDataContainer() "
            << gradientCalculator->GetDataContainer() << "\n";

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

  pointIterator = fixedMesh->GetPoints()->Begin();
  pointEnd = fixedMesh->GetPoints()->End();

  typedef GradientCalculatorType::InputType  InputType;
  typedef GradientCalculatorType::OutputType OutputType;
  InputType pointId = 0;
  float     maximumDifferenceMagnitude = 0.0;
  while( pointIterator != pointEnd )
    {
    OutputType     analyticalDerivative;
    OutputType     computedDerivative;
    float          analyticalValue;
    FixedPointType fixedPoint;
    fixedPoint.CastFrom( pointIterator.Value() );

    mappedMovingMeshPoint = transform->TransformPoint( fixedPoint );
    const VectorType radial = mappedMovingMeshPoint - myCenter;
    const double     radius = radial.GetNorm(); // assuming radius is not valued 1

    const double theta = vcl_atan2(mappedMovingMeshPoint[1], mappedMovingMeshPoint[0]);
    const double phi = vcl_acos(mappedMovingMeshPoint[2] / radius);

    float thetaOffset = theta - piOver4;
    if( thetaOffset < -pi )
      {
      thetaOffset += twoPi;
      }

    analyticalValue =
      mapSphericalCoordinatesFunction(phi, thetaOffset);

    analyticalDerivative =
      mapSphericalCoordinatesFunctionGradient(phi, thetaOffset, false);

    computedDerivative = gradientCalculator->Evaluate( pointId );

    std::cout << "  pointId " << pointId
              << "  analytical derivative at moving mesh  "
              << analyticalDerivative << std::endl;

    std::cout << "  computed derivative fixed vertex values"
              << computedDerivative
              << std::endl;

    OutputType differenceOfDerivatives = analyticalDerivative - computedDerivative;
    float      differenceMagnitude = differenceOfDerivatives.GetNorm();

    if( differenceMagnitude >  maximumDifferenceMagnitude )
      {
      maximumDifferenceMagnitude = differenceMagnitude;

      std::cout << " maximum difference magnitude increasing " << maximumDifferenceMagnitude
                << "  pointId " << pointId << std::endl;
      }

    pointIterator++;
    pointId++;
    }

  std::cout << " maximum difference magnitude " << maximumDifferenceMagnitude << std::endl;

  std::cout << "Test End " << std::endl;

  return EXIT_SUCCESS;
}
