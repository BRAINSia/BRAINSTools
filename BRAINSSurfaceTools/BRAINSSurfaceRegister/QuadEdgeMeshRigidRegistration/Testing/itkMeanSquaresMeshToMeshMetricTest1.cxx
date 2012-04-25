/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMeanSquaresMeshToMeshMetricTest1.cxx,v $
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

#include "itkVersorTransform.h"
#include "itkQuadEdgeMesh.h"
#include "itkRegularSphereMeshSource.h"
#include "itkMeanSquaresMeshToMeshMetric.h"
#include "itkNearestNeighborInterpolateMeshFunction.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkCommand.h"
#include "itkTestingMacros.h"
#include "itkMeshGeneratorHelper.h"
#include "itkEllipseSpatialObject.h"

int main( int, char * [] )
{
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<float, Dimension> MovingMeshType;
  typedef itk::QuadEdgeMesh<float, Dimension> FixedMeshType;

  FixedMeshType::Pointer  fixedMesh;
  MovingMeshType::Pointer movingMesh;

  typedef MeshGeneratorHelper<FixedMeshType, MovingMeshType> GeneratorType;

  GeneratorType::GenerateMeshes( fixedMesh, movingMesh );

  // -----------------------------------------------------------
  // Set up  the Metric
  // -----------------------------------------------------------
  typedef itk::MeanSquaresMeshToMeshMetric<
      FixedMeshType, MovingMeshType>   MetricType;

  MetricType::Pointer metric = MetricType::New();

  TRY_EXPECT_EXCEPTION( metric->Initialize() );

  // -----------------------------------------------------------
  // Plug the Meshes into the metric
  // -----------------------------------------------------------
  metric->SetFixedMesh( fixedMesh );

  TRY_EXPECT_EXCEPTION( metric->Initialize() );

  metric->SetMovingMesh( movingMesh );

  TRY_EXPECT_EXCEPTION( metric->Initialize() );

  // -----------------------------------------------------------
  // Set up a Transform
  // -----------------------------------------------------------
  typedef itk::VersorTransform<MetricType::TransformComputationType> TransformType;

  TransformType::Pointer transform = TransformType::New();

  TRY_EXPECT_EXCEPTION( metric->Initialize() );

  metric->SetTransform( transform );

// ------------------------------------------------------------
// Set up an Interpolator
// ------------------------------------------------------------
  typedef itk::NearestNeighborInterpolateMeshFunction<MovingMeshType> InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  interpolator->SetInputMesh( movingMesh );

  metric->SetInterpolator( interpolator );

  interpolator->Initialize();

  std::cout << metric << std::endl;

// ------------------------------------------------------------
// Set up transform parameters
// ------------------------------------------------------------
  const unsigned int numberOfTransformParameters =
    transform->GetNumberOfParameters();

  typedef MetricType::TransformType         TransformBaseType;
  typedef TransformBaseType::ParametersType ParametersType;

  ParametersType parameters( numberOfTransformParameters );

  const unsigned int numberOfSampleValues = 10;

  const float toRadians = 8.0 * atan( 1.0 );

  typedef itk::Versor<double> VersorType;
  VersorType versor;

  VersorType::VectorType axis;
  axis[0] = 0.0;
  axis[1] = 0.0;
  axis[2] = 1.0;

  MetricType::DerivativeType derivative1;
  MetricType::DerivativeType derivative2;

  MetricType::MeasureType value1;
  MetricType::MeasureType value2;

  const double tolerance = 1e-5;
  for( unsigned int p = 0; p < numberOfSampleValues; p++ )
    {
    const float angle = p * toRadians / numberOfSampleValues;

    versor.Set( axis, angle );

    parameters[0] = versor.GetX();
    parameters[1] = versor.GetY();
    parameters[2] = versor.GetZ();

    value1 = metric->GetValue( parameters );
    metric->GetDerivative( parameters, derivative1 );
    metric->GetValueAndDerivative( parameters, value2, derivative2 );

    if( vnl_math_abs( value2 - value1 ) > tolerance )
      {
      std::cerr << "GetValue() does not match GetValueAndDerivative()" << std::endl;
      return EXIT_FAILURE;
      }
    }

  // -------------------------------------------------------
  // exercise Print() method
  // -------------------------------------------------------
  metric->Print( std::cout );

  // -------------------------------------------------------
  // exercise misc member functions
  // -------------------------------------------------------
  std::cout << "FixedMesh: " << metric->GetFixedMesh() << std::endl;
  std::cout << "MovingMesh: " << metric->GetMovingMesh() << std::endl;
  std::cout << "Transform: " << metric->GetTransform() << std::endl;

  std::cout << "Check case when Target is NULL" << std::endl;
  metric->SetFixedMesh( NULL );

  TRY_EXPECT_EXCEPTION( metric->GetValue( parameters ) );
  TRY_EXPECT_EXCEPTION( metric->GetDerivative( parameters, derivative1 ) );
  TRY_EXPECT_EXCEPTION( metric->GetValueAndDerivative( parameters, value1, derivative1  ) );

  // ---------------------------------------------------
  // exercise masks
  // ---------------------------------------------------
  std::cout << "Exercising Masks" << std::endl;
  metric->SetFixedMesh( fixedMesh );

  typedef itk::EllipseSpatialObject<Dimension> MaskType;

  MaskType::Pointer fixedMask  = MaskType::New();
  fixedMask->SetRadius(0.1);

  metric->SetFixedMask( fixedMask );
  TEST_SET_GET( fixedMask, metric->GetFixedMask() );

  // Since the masks are outside, we expect exceptions in the following calls.
  TRY_EXPECT_EXCEPTION( metric->GetValue( parameters ) );
  TRY_EXPECT_EXCEPTION( metric->GetDerivative( parameters, derivative1 ) );
  TRY_EXPECT_EXCEPTION( metric->GetValueAndDerivative( parameters, value1, derivative1  ) );

  fixedMask->SetRadius(5.1); // enlarge radius to accept all points.

  MaskType::Pointer movingMask = MaskType::New();
  movingMask->SetRadius(0.1);

  metric->SetMovingMask( movingMask );
  TEST_SET_GET( movingMask, metric->GetMovingMask() );

  // Since the masks are outside, we expect exceptions in the following calls.
  TRY_EXPECT_EXCEPTION( metric->GetValue( parameters ) );
  TRY_EXPECT_EXCEPTION( metric->GetDerivative( parameters, derivative1 ) );
  TRY_EXPECT_EXCEPTION( metric->GetValueAndDerivative( parameters, value1, derivative1  ) );

  movingMask->SetRadius(5.1);

  // Since the masks are outside, we expect exceptions in the following calls.
  TRY_EXPECT_NO_EXCEPTION( metric->GetValue( parameters ) );
  TRY_EXPECT_NO_EXCEPTION( metric->GetDerivative( parameters, derivative1 ) );
  TRY_EXPECT_NO_EXCEPTION( metric->GetValueAndDerivative( parameters, value1, derivative1  ) );

  std::cout << "Test passed. " << std::endl;

  return EXIT_SUCCESS;
}
