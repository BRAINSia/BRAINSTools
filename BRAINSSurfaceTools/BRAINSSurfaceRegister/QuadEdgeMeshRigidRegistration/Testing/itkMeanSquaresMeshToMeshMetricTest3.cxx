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

#include "itkVersorRigid3DTransform.h"
#include "itkQuadEdgeMesh.h"
#include "itkRegularSphereMeshSource.h"
#include "itkMeanSquaresMeshToMeshMetric.h"
#include "itkMeshToMeshRegistrationMethod.h"
#include "itkNearestNeighborInterpolateMeshFunction.h"
#include "itkAmoebaOptimizer.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkCommand.h"
#include "itkMeshGeneratorHelper.h"

class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate Self;
  typedef  itk::Command           Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate()
  {
    this->IterationCounter = 0;
  }

public:
  typedef itk::AmoebaOptimizer    OptimizerType;
  typedef   const OptimizerType * OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    OptimizerPointer optimizer =
      dynamic_cast<OptimizerPointer>( object );

    if( !itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }
    std::cout << "  Iteration " << IterationCounter++ << "   ";
    std::cout << optimizer->GetCachedValue() << "   ";
    std::cout << optimizer->GetCachedCurrentPosition() << std::endl;
  }

private:
  unsigned int IterationCounter;
};

int main( int, char * [] )
{
  typedef itk::QuadEdgeMesh<float, 3> MovingMeshType;
  typedef itk::QuadEdgeMesh<float, 3> FixedMeshType;

  FixedMeshType::Pointer  fixedMesh;
  MovingMeshType::Pointer movingMesh;

  typedef MeshGeneratorHelper<FixedMeshType, MovingMeshType> GeneratorType;

  GeneratorType::GenerateMeshes( fixedMesh, movingMesh );

  // Registration Method
  typedef itk::MeshToMeshRegistrationMethod<
      FixedMeshType,
      MovingMeshType>    RegistrationType;
  RegistrationType::Pointer registration  = RegistrationType::New();

// -----------------------------------------------------------
// Set up  the Metric
// -----------------------------------------------------------
  typedef itk::MeanSquaresMeshToMeshMetric<FixedMeshType,
                                           MovingMeshType>
    MetricType;

  typedef MetricType::TransformType         TransformBaseType;
  typedef TransformBaseType::ParametersType ParametersType;

  MetricType::Pointer metric = MetricType::New();
  registration->SetMetric( metric );

// -----------------------------------------------------------
// Plug the Meshes into the metric
// -----------------------------------------------------------
  registration->SetFixedMesh( fixedMesh );
  registration->SetMovingMesh( movingMesh );

// -----------------------------------------------------------
// Set up a Transform
// -----------------------------------------------------------

  typedef itk::VersorTransform<MetricType::TransformComputationType> TransformType;

  TransformType::Pointer transform = TransformType::New();

  registration->SetTransform( transform );

// ------------------------------------------------------------
// Set up an Interpolator
// ------------------------------------------------------------
  typedef itk::NearestNeighborInterpolateMeshFunction<MovingMeshType> InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  interpolator->SetInputMesh( movingMesh );

  registration->SetInterpolator( interpolator );

  interpolator->Initialize();

  std::cout << metric << std::endl;

// ------------------------------------------------------------
// Set up transform parameters
// ------------------------------------------------------------
  const unsigned int numberOfTransformParameters =
    transform->GetNumberOfParameters();

  ParametersType parameters( numberOfTransformParameters );
  // initialize the offset/vector part
  for( unsigned int k = 0; k < numberOfTransformParameters; k++ )
    {
    parameters[k] = 0.0f;
    }

  registration->SetInitialTransformParameters( parameters );

  // Optimizer Type
  typedef itk::AmoebaOptimizer OptimizerType;

  OptimizerType::Pointer optimizer     = OptimizerType::New();

  OptimizerType::ParametersType simplexDelta( numberOfTransformParameters );
  simplexDelta.Fill( 0.1 );

  optimizer->AutomaticInitialSimplexOff();
  optimizer->SetInitialSimplexDelta( simplexDelta );
  optimizer->SetParametersConvergenceTolerance( 0.01 ); // radiands ?
  optimizer->SetFunctionConvergenceTolerance(0.001);    // 0.1%
  optimizer->SetMaximumNumberOfIterations( 200 );

  // Create the Command observer and register it with the optimizer.
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  registration->SetOptimizer( optimizer );

  try
    {
    registration->StartRegistration();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Metric initialization failed" << std::endl;
    std::cout << "Reason " << e.GetDescription() << std::endl;
    return EXIT_FAILURE;
    }

// ---------------------------------------------------------
// Print out metric values
// for parameters[1] = {-10,10}  (arbitrary choice...)
// ---------------------------------------------------------

  OptimizerType::ParametersType finalParameters =
    registration->GetLastTransformParameters();

  const double bestValue = optimizer->GetValue();

  std::cout << "final params " << finalParameters << "  final value " << bestValue << std::endl;

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
  try
    {
    std::cout << "Value = " << metric->GetValue( parameters );
    std::cout << "If you are reading this message the Metric " << std::endl;
    std::cout << "is NOT managing exceptions correctly    " << std::endl;
    return EXIT_FAILURE;
    }
  catch( itk::ExceptionObject & e )
    {
    std::cout << "Exception received (as expected) "    << std::endl;
    std::cout << "Description : " << e.GetDescription() << std::endl;
    std::cout << "Location    : " << e.GetLocation()    << std::endl;
    std::cout << "Test for exception throwing... PASSED ! " << std::endl;
    }

  std::cout << "Test passed. " << std::endl;

  return EXIT_SUCCESS;
}
