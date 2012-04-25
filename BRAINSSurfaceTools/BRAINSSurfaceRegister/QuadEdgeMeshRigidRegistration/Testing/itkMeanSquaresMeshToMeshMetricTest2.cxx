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
#include "itkVersorTransformOptimizer.h"
#include "itkQuadEdgeMesh.h"
#include "itkRegularSphereMeshSource.h"
#include "itkMeanSquaresMeshToMeshMetric.h"
#include "itkMeshToMeshRegistrationMethod.h"
#include "itkLinearInterpolateMeshFunction.h"
#include "itkGradientDescentOptimizer.h"
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
  typedef itk::VersorTransformOptimizer OptimizerType;
  typedef   const OptimizerType *       OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    OptimizerPointer optimizer = dynamic_cast<OptimizerPointer>( object );

    if( !itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }

    std::cout << "  Iteration " << IterationCounter++ << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
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

  transform->SetIdentity();

  registration->SetTransform( transform );

// ------------------------------------------------------------
// Set up an Interpolator
// ------------------------------------------------------------
  typedef itk::LinearInterpolateMeshFunction<MovingMeshType> InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  interpolator->SetInputMesh( movingMesh );

  registration->SetInterpolator( interpolator );

  interpolator->Initialize();

  std::cout << metric << std::endl;

// ------------------------------------------------------------
// Set up optimizer
// ------------------------------------------------------------

  typedef itk::VersorTransformOptimizer OptimizerType;

  OptimizerType::Pointer optimizer     = OptimizerType::New();

// ------------------------------------------------------------
// Set up transform parameters
// ------------------------------------------------------------

  typedef OptimizerType::ScalesType ScalesType;
  typedef TransformType::VersorType VersorType;

  // We start with a null rotation
  VersorType::VectorType axis;
  axis[0] =  1.0f;
  axis[1] =  0.0f;
  axis[2] =  0.0f;

  VersorType::ValueType angle = 0.0f;

  VersorType initialRotation;

  initialRotation.Set( axis, angle );

  transform->SetRotation( initialRotation );

  registration->SetInitialTransformParameters( transform->GetParameters() );

  const unsigned int numberOfTransformParameters = transform->GetNumberOfParameters();

  ScalesType parametersScale( numberOfTransformParameters );
  parametersScale[0] = 1.0;
  parametersScale[1] = 1.0;
  parametersScale[2] = 1.0;

  optimizer->MinimizeOn();
  optimizer->SetScales( parametersScale );
  optimizer->SetGradientMagnitudeTolerance( 1e-15 );
  optimizer->SetMaximumStepLength( 0.1745 ); // About 10 degrees
  optimizer->SetMinimumStepLength( 1e-9 );
  optimizer->SetNumberOfIterations( 30 );

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

  std::cout << "Test passed. " << std::endl;

  return EXIT_SUCCESS;
}
