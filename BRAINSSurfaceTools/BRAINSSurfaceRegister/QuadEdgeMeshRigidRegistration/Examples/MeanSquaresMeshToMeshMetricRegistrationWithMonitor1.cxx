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

#include "itkMeanSquaresMeshToMeshMetric.h"
#include "itkMeshToMeshRegistrationMethod.h"
#include "itkLinearInterpolateMeshFunction.h"
#include "itkVersorTransformOptimizer.h"
#include "itkVersorTransform.h"
#include "itkQuadEdgeMesh.h"

#include "itkCommand.h"
#include "itkVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

#include "AffineRegistrationMonitor.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

#include "itkResampleQuadEdgeMeshFilter.h"

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
    iterationCounter = 0;
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

    std::cout << " Iteration " << ++iterationCounter;
    std::cout << "  Value " << optimizer->GetValue() << "   ";
    std::cout << "  Position " << optimizer->GetCurrentPosition() << std::endl;
  }

private:

  unsigned int iterationCounter;
};

int main( int argc, char * argv [] )
{
  if( argc < 8 )
    {
    std::cerr << "Missing arguments" << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << std::endl;
    std::cerr << "inputFixedMesh inputMovingMesh ";
    std::cerr << "outputResampledMesh ";
    std::cerr << "axisX axisY axisZ angle(radians) " << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::QuadEdgeMesh<float, 3> FixedMeshType;
  typedef itk::QuadEdgeMesh<float, 3> MovingMeshType;

  typedef itk::VTKPolyDataReader<FixedMeshType>  FixedReaderType;
  typedef itk::VTKPolyDataReader<MovingMeshType> MovingReaderType;

  FixedReaderType::Pointer fixedMeshReader = FixedReaderType::New();
  fixedMeshReader->SetFileName( argv[1] );

  MovingReaderType::Pointer movingMeshReader = MovingReaderType::New();
  movingMeshReader->SetFileName( argv[2] );

  try
    {
    fixedMeshReader->Update();
    movingMeshReader->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  FixedMeshType::ConstPointer  meshFixed  = fixedMeshReader->GetOutput();
  MovingMeshType::ConstPointer meshMoving = movingMeshReader->GetOutput();

  typedef itk::MeshToMeshRegistrationMethod<
      FixedMeshType,
      MovingMeshType>    RegistrationType;

  RegistrationType::Pointer registration  = RegistrationType::New();

  typedef itk::MeanSquaresMeshToMeshMetric<FixedMeshType,
                                           MovingMeshType>
    MetricType;

  MetricType::Pointer metric = MetricType::New();

  registration->SetMetric( metric );

  registration->SetFixedMesh( meshFixed );
  registration->SetMovingMesh( meshMoving );

  typedef itk::VersorTransform<MetricType::TransformComputationType> TransformType;

  TransformType::Pointer transform = TransformType::New();

  registration->SetTransform( transform );

  typedef itk::LinearInterpolateMeshFunction<MovingMeshType> InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  registration->SetInterpolator( interpolator );

  const unsigned int numberOfTransformParameters =
    transform->GetNumberOfParameters();

  typedef TransformType::ParametersType ParametersType;
  ParametersType parameters( numberOfTransformParameters );

  TransformType::AxisType  axis;
  TransformType::AngleType angle;

  axis[0] = atof( argv[4] );
  axis[1] = atof( argv[5] );
  axis[2] = atof( argv[6] );

  angle = atof( argv[7] );

  transform->SetRotation( axis, angle );

  parameters = transform->GetParameters();

  registration->SetInitialTransformParameters( parameters );

  typedef itk::VersorTransformOptimizer OptimizerType;

  OptimizerType::Pointer optimizer     = OptimizerType::New();

  registration->SetOptimizer( optimizer );

  typedef OptimizerType::ScalesType ScalesType;

  ScalesType parametersScale( numberOfTransformParameters );
  parametersScale[0] = 1.0;
  parametersScale[1] = 1.0;
  parametersScale[2] = 1.0;

  optimizer->SetScales( parametersScale );

  optimizer->MinimizeOn();
  optimizer->SetGradientMagnitudeTolerance( 1e-6 );
  optimizer->SetMaximumStepLength( 0.05 );
  optimizer->SetMinimumStepLength( 1e-9 );
  optimizer->SetRelaxationFactor( 0.9 );
  optimizer->SetNumberOfIterations( 30 );

  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  // Configuration of Visual Registration Monitor
  vtkSmartPointer<vtkPolyDataReader> vtkFixedMeshReader =
    vtkSmartPointer<vtkPolyDataReader>::New();

  vtkSmartPointer<vtkPolyDataReader> vtkMovingMeshReader =
    vtkSmartPointer<vtkPolyDataReader>::New();

  vtkFixedMeshReader->SetFileName( argv[1] );
  vtkMovingMeshReader->SetFileName( argv[2] );

  vtkFixedMeshReader->Update();
  vtkMovingMeshReader->Update();

  AffineRegistrationMonitor visualMonitor;

  visualMonitor.SetFixedSurface( vtkFixedMeshReader->GetOutput() );
  visualMonitor.SetMovingSurface( vtkMovingMeshReader->GetOutput() );

  visualMonitor.SetScalarRange( -1.0, 1.0 );

  visualMonitor.SetNumberOfIterationsPerUpdate( 1 );

  visualMonitor.Observe( optimizer.GetPointer() );
  visualMonitor.ObserveData( transform.GetPointer() );
  visualMonitor.SetScreenShotsBaseFileName("rigidRegistration");

  // visualMonitor.SetVerbose( false );
  visualMonitor.SetVerbose( true );

  try
    {
    registration->StartRegistration();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Registration failed" << std::endl;
    std::cout << "Reason " << e << std::endl;
    return EXIT_FAILURE;
    }

  OptimizerType::ParametersType finalParameters =
    registration->GetLastTransformParameters();

  const double bestValue = optimizer->GetValue();

  std::cout << "final parameters = " << finalParameters << std::endl;
  std::cout << "final value      = " << bestValue << std::endl;

  transform->SetParameters( finalParameters );

  typedef itk::ResampleQuadEdgeMeshFilter<
      FixedMeshType, MovingMeshType>  ResamplingFilterType;

  ResamplingFilterType::Pointer resampler = ResamplingFilterType::New();

  resampler->SetReferenceMesh( meshFixed );
  resampler->SetInput( meshMoving );

  resampler->SetTransform( transform );
  resampler->SetInterpolator( interpolator );

  resampler->Update();

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<FixedMeshType> WriterType;

  WriterType::Pointer writer = WriterType::New();

  writer->SetInput( resampler->GetOutput()  );
  writer->SetFileName( argv[3] );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error during writer Update() " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
