/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNodeVectorJacobianCalculatorTest2.cxx,v $
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

#include "itkNodeVectorJacobianCalculator.h"
#include "itkQuadEdgeMesh.h"
#include "itkCovariantVector.h"
#include "itkTestingMacros.h"
#include "itkTriangleBasisSystem.h"
#include "itkTriangleBasisSystemCalculator.h"
#include "itkTriangleListBasisSystemCalculator.h"
#include "itkVersorTransform.h"
#include "itkVectorContainer.h"

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNodeScalarGradientCalculatorTest2.cxx,v $
  Language:  C++
  Date:      $Date: 2007-09-06 17:44:24 $
  Version:   $Revision: 1.3 $

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
#include "itkMatrix.h"
#include "itkRegularSphereMeshSource.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkNodeVectorJacobianCalculator.h"

#define MeshDimension 3
typedef itk::QuadEdgeMesh<float, MeshDimension> FixedMeshType;
typedef FixedMeshType::PointType                FixedPointType;
typedef FixedPointType::VectorType              FixedVectorType;
typedef itk::Matrix<
    float, MeshDimension, MeshDimension>                           JacobianType;

static void mapSphericalCoordinatesVectorFunction(float inPhi, float inTheta, FixedVectorType & result);

static void mapSphericalCoordinatesVectorFunctionJacobian(float inPhi, float inTheta, JacobianType & jacobianResult,
                                                          bool printFlag);

// Really simple example: each component is a sinusoid mapping between
// 0 and 1 as a function of theta, constant in phi
static void
mapSphericalCoordinatesVectorFunction(float itkNotUsed(inPhi), float inTheta, FixedVectorType & result)
{
  // Really simple sinusoidal vector function
  result[0] = vcl_sin(inTheta);
  result[1] = vcl_cos(inTheta);
  result[2] = -vcl_sin(inTheta);
}

static void
mapSphericalCoordinatesVectorFunctionJacobian(float inPhi, float inTheta,
                                              JacobianType & jacobianFunction,
                                              bool printFlag)
{
  JacobianType phiComponent;
  JacobianType thetaComponent;

  float cosTheta = vcl_cos(inTheta);
  float sinTheta = vcl_sin(inTheta);
  float cosPhi = vcl_cos(inPhi);
  float sinPhi = vcl_sin(inPhi);

  // derivative of phiFactor over Theta
  FixedVectorType functionDerivativeOverTheta;

  functionDerivativeOverTheta[0] = vcl_cos(inTheta);
  functionDerivativeOverTheta[1] = -vcl_sin(inTheta);
  functionDerivativeOverTheta[2] = -vcl_cos(inTheta);
  // Need to multiply dF/dtheta by unit vector
  // unit vector= vcl_cos(phi)*vcl_cos(theta) i + vcl_cos(phi)*vcl_sin(theta) j - vcl_sin(phi)k
  for( int i = 0; i < 3; i++ )
    {
    thetaComponent[i][0] = -sinTheta * functionDerivativeOverTheta[i];
    thetaComponent[i][1] = cosTheta * functionDerivativeOverTheta[i];
    thetaComponent[i][2] = 0.0;
    }

  jacobianFunction = thetaComponent;

  if( printFlag )
    {
    std::cout << "  inTheta " << inTheta << "  inPhi " << inPhi
              << "  sinTheta " << sinTheta << "  cosTheta " << cosTheta
              << "  sinPhi " << sinPhi << "  cosPhi " << cosPhi
              << "  dfdphi " << functionDerivativeOverTheta
              << "  thetaComponent " << thetaComponent << " \n";
    }

  return;
}

int main( int, char * [] )
{
  // const double piOver4 = atan( 1.0 );
  // const double pi = piOver4 * 4.0;

  typedef itk::QuadEdgeMesh<float, MeshDimension> FixedMeshType;
  FixedMeshType::Pointer fixedMesh;

  typedef itk::CovariantVector<float, MeshDimension>     ArrayType;
  typedef itk::VectorContainer<unsigned long, ArrayType> FixedPointDataContainerType;

  typedef itk::NodeVectorJacobianCalculator<
      FixedMeshType, FixedPointDataContainerType>             JacobianCalculatorType;

  JacobianCalculatorType::Pointer jacobianCalculator = JacobianCalculatorType::New();

  FixedPointType fixedCenter;
  fixedCenter.Fill( 0.0 );
  FixedVectorType fixedScale;
  fixedScale.Fill( 1.0 );

  typedef itk::RegularSphereMeshSource<FixedMeshType> FixedSphereMeshSourceType;

  FixedSphereMeshSourceType::Pointer fixedSphereMeshSource = FixedSphereMeshSourceType::New();

  fixedSphereMeshSource->SetCenter( fixedCenter );
  fixedSphereMeshSource->SetResolution( 4 );
  fixedSphereMeshSource->SetScale( fixedScale );
  fixedSphereMeshSource->Modified();

  try
    {
    fixedSphereMeshSource->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error during source Update() " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  fixedMesh = fixedSphereMeshSource->GetOutput();

  typedef FixedMeshType::PointIdentifier
    FixedPointIdentifier;
  typedef itk::VectorContainer<FixedPointIdentifier, FixedVectorType>
    VectorContainerType;
  FixedPointDataContainerType::Pointer fixedVectorContainer = FixedPointDataContainerType::New();

  fixedVectorContainer->Reserve( fixedMesh->GetNumberOfPoints() );

  FixedPointType fixedPoint;
  fixedPoint.Fill(0.0f);
  for( unsigned int i = 0; i < fixedMesh->GetNumberOfPoints(); i++ )
    {
    fixedMesh->GetPoint(i, &fixedPoint);

    FixedVectorType fixedVtr = fixedPoint - fixedCenter;

    const float radius = fixedVtr.GetNorm();
    const float theta  = atan2(fixedVtr[1], fixedVtr[0]);
    const float phi    = acos(fixedVtr[2] / radius);

    FixedVectorType fixedVector;
    mapSphericalCoordinatesVectorFunction(phi, theta, fixedVector);

    fixedVectorContainer->ElementAt(i) = fixedVector;
    }

  std::cout << jacobianCalculator->GetNameOfClass() << std::endl;
  jacobianCalculator->Print( std::cout );

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

  jacobianCalculator->SetInputMesh( triangleListBasisSystemCalculator->GetInputMesh() );

  if( jacobianCalculator->GetInputMesh() != triangleListBasisSystemCalculator->GetInputMesh() )
    {
    std::cerr << "Error in SetInputMesh()/GetInputMesh() " << std::endl;
    return EXIT_FAILURE;
    }

  PointType myCenter;
  myCenter.Fill( 0.0 );

  typedef FixedMeshType::PointsContainer      FixedPointsContainer;
  typedef FixedPointsContainer::ConstIterator FixedPointIterator;

  FixedPointIterator pointIterator = fixedMesh->GetPoints()->Begin();
  FixedPointIterator pointEnd = fixedMesh->GetPoints()->End();

  typedef FixedMeshType::PointType FixedPointType;

  jacobianCalculator->SetVectorContainer( fixedVectorContainer );

  if( jacobianCalculator->GetVectorContainer() != fixedVectorContainer )
    {
    std::cerr << "Error in SetVectorContainer()/GetVectorContainer() " << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "test 2 jacobianCalculator->GetVectorContainer() "
            << jacobianCalculator->GetVectorContainer() << "\n";

  jacobianCalculator->SetBasisSystemList( triangleListBasisSystemCalculator->GetBasisSystemList() );

  if( jacobianCalculator->GetBasisSystemList() != triangleListBasisSystemCalculator->GetBasisSystemList() )
    {
    std::cerr << "Error in SetBasisSystemList()/GetBasisSystemList() " << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "test 2 jacobianCalculator->GetVectorContainer() "
            << jacobianCalculator->GetVectorContainer() << "\n";

  // It is not initialized correctly, we expect no exception
  TRY_EXPECT_NO_EXCEPTION( jacobianCalculator->Initialize(); );

  // Now that Initialize() has been called successfully, we can call Compute(). */
  TRY_EXPECT_NO_EXCEPTION( jacobianCalculator->Compute(); );

  pointIterator = fixedMesh->GetPoints()->Begin();
  pointEnd = fixedMesh->GetPoints()->End();

  typedef JacobianCalculatorType::InputType  InputType;
  typedef JacobianCalculatorType::OutputType OutputType;

  FixedPointIdentifier pointId = 0;
  float                maximumDifferenceMagnitude = 0.0;
  while( pointIterator != pointEnd )
    {
    JacobianType analyticalJacobian;
    OutputType   computedJacobian;
    fixedPoint.CastFrom( pointIterator.Value() );

    const VectorType radial = fixedPoint - myCenter;
    const double     radius = radial.GetNorm(); // assuming radius is not valued 1
    const double     theta = vcl_atan2(radial[1], radial[0]);
    const double     phi = vcl_acos(radial[2] / radius);

    mapSphericalCoordinatesVectorFunctionJacobian(phi, theta, analyticalJacobian, false);

    computedJacobian = jacobianCalculator->Evaluate( pointId );

    std::cout << "  pointId " << pointId
              << "  analytical derivative at moving mesh  "
              << analyticalJacobian << std::endl;

    std::cout << "  computed derivative fixed vertex values"
              << computedJacobian
              << std::endl;

    JacobianType differenceOfJacobians;

    float largestDifferenceMatrixElement = 0.0;
    for( int i = 0; i < 3; i++ )
      {
      for( int j = 0; j < 3; j++ )
        {
        differenceOfJacobians[i][j] = analyticalJacobian[i][j] - computedJacobian[i][j];
        if( fabs( differenceOfJacobians[i][j]) > largestDifferenceMatrixElement )
          {
          largestDifferenceMatrixElement = fabs(differenceOfJacobians[i][j]);
          }
        }
      }

    if( largestDifferenceMatrixElement >  maximumDifferenceMagnitude )
      {
      maximumDifferenceMagnitude = largestDifferenceMatrixElement;

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
