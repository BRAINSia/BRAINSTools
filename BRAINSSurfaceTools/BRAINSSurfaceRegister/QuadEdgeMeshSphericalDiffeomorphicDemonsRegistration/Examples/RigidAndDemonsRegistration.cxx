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

#include "itkMeanSquaresMeshToMeshMetric.h"
#include "itkMeshToMeshRegistrationMethod.h"
#include "itkLinearInterpolateMeshFunction.h"
#include "itkVersorTransformOptimizer.h"
#include "itkVersorTransform.h"
#include "itkQuadEdgeMesh.h"

#include "itkCommand.h"
#include "itkVTKPolyDataReader.h"

#include "itkQuadEdgeMeshSphericalDiffeomorphicDemonsRegistrationConfigure.h"
#include "itkResampleQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshVectorDataVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshSphericalDiffeomorphicDemonsFilter.h"
#include "itkDeformationFieldFromTransformMeshFilter.h"

#include "itkQuadEdgeMeshGenerateDeformationFieldFilter.h"
#include "itkIdentityTransform.h"
#include "itkReplaceDestinationPointsQuadEdgeMeshFilter.h"
#include "itkAssignScalarValuesQuadEdgeMeshFilter.h"

#ifdef USE_VTK
#include "DeformableAndAffineRegistrationMonitor.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#endif

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
  if( argc < 10 )
    {
    std::cerr << "Missing arguments" << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << std::endl;
    std::cerr << "inputFixedMesh inputMovingMesh radius ";
    std::cerr << "sigmaX maximumNumberOfDemonsIterations ";
    std::cerr << "outputResampledRigid outputResampledDemons ";
    std::cerr << "outputDeformedFixedMesh DeformationField ";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  typedef float MeshPixelType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> FixedMeshType;
  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> MovingMeshType;

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
  std::cout << "read fixed and moving mesh!" << std::endl;

  FixedMeshType::Pointer  meshFixed  = fixedMeshReader->GetOutput();
  MovingMeshType::Pointer meshMoving = movingMeshReader->GetOutput();

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

  const unsigned int numberOfTransformParameters = transform->GetNumberOfParameters();

  typedef TransformType::ParametersType ParametersType;
  ParametersType parameters( numberOfTransformParameters );

  transform->SetIdentity();

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
  optimizer->SetMaximumStepLength( 0.01 );
  optimizer->SetMinimumStepLength( 1e-9 );
  optimizer->SetRelaxationFactor( 0.9 );
  optimizer->SetNumberOfIterations( 30 );

  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> RegisteredMeshType;

  typedef itk::QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<
      FixedMeshType, MovingMeshType, RegisteredMeshType>   DemonsFilterType;

  typedef DemonsFilterType::DestinationPointSetType DestinationPointSetType;

#ifdef USE_VTK
  typedef DeformableAndAffineRegistrationMonitor<
      DemonsFilterType, DestinationPointSetType> RegistrationMonitorType;

  RegistrationMonitorType visualMonitor;
  visualMonitor.SetNumberOfIterationsPerUpdate( 1 );

  visualMonitor.SetBaseAnnotationText("Rigid Registration ");

  visualMonitor.Observe( optimizer.GetPointer() );
  visualMonitor.ObserveData( transform.GetPointer() );

  visualMonitor.SetVerbose( false );
  visualMonitor.SetScreenShotsBaseFileName( "rigidAndDemonsRegistration" );

  visualMonitor.SetScalarRange( -1.0, 1.0 );

  vtkSmartPointer<vtkPolyDataReader> vtkFixedMeshReader =
    vtkSmartPointer<vtkPolyDataReader>::New();

  vtkSmartPointer<vtkPolyDataReader> vtkMovingMeshReader =
    vtkSmartPointer<vtkPolyDataReader>::New();

  vtkFixedMeshReader->SetFileName( fixedMeshReader->GetFileName() );
  vtkMovingMeshReader->SetFileName( movingMeshReader->GetFileName() );

  vtkFixedMeshReader->Update();
  vtkMovingMeshReader->Update();

  visualMonitor.SetFixedSurface( vtkFixedMeshReader->GetOutput() );
  visualMonitor.SetMovingSurface( vtkMovingMeshReader->GetOutput() );
#endif

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

  std::cout << "final parameters = " << finalParameters << std::endl;
  std::cout << "final value      = " << optimizer->GetValue() << std::endl;

  transform->SetParameters( finalParameters );

  typedef FixedMeshType::Traits                               MeshTraits;
  typedef itk::PointSet<MeshPixelType, Dimension, MeshTraits> PointSetType;

  typedef itk::DeformationFieldFromTransformMeshFilter<
      FixedMeshType, PointSetType>  DeformationFieldFromTransformFilterType;

  DeformationFieldFromTransformFilterType::Pointer deformationFieldFromTransform =
    DeformationFieldFromTransformFilterType::New();

  deformationFieldFromTransform->SetInput( fixedMeshReader->GetOutput() );
  deformationFieldFromTransform->SetTransform( transform );

  try
    {
    deformationFieldFromTransform->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cout << excp << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Transform Parameters = " << transform->GetParameters() << std::endl;

  PointSetType::ConstPointer destinationPoints = deformationFieldFromTransform->GetOutput();

  typedef FixedMeshType::PointType PointType;
  typedef PointType::VectorType    VectorType;

  typedef itk::QuadEdgeMeshTraits<VectorType, Dimension, bool, bool> VectorPointSetTraits;

  typedef itk::QuadEdgeMesh<VectorType, Dimension, VectorPointSetTraits> MeshWithVectorsType;

  typedef itk::QuadEdgeMeshGenerateDeformationFieldFilter<
      FixedMeshType, PointSetType, MeshWithVectorsType>   DeformationFilterType;

  DeformationFilterType::Pointer deformationFilter = DeformationFilterType::New();

  deformationFilter->SetInputMesh( fixedMeshReader->GetOutput() );
  deformationFilter->SetDestinationPoints( destinationPoints );
  deformationFilter->Update();

  typedef itk::QuadEdgeMeshVectorDataVTKPolyDataWriter<MeshWithVectorsType> VectorMeshWriterType;
  VectorMeshWriterType::Pointer vectorMeshWriter = VectorMeshWriterType::New();
  vectorMeshWriter->SetInput( deformationFilter->GetOutput() );
  vectorMeshWriter->SetFileName("VectorMesh0.vtk");
  vectorMeshWriter->Update();
  std::cout << "Deformation VectorMesh0.vtk  Saved" << std::endl;

  // write resampled moving mesh after rigid
  typedef itk::ResampleQuadEdgeMeshFilter<FixedMeshType, MovingMeshType> ResamplingFilterType;

  ResamplingFilterType::Pointer resampler = ResamplingFilterType::New();

  resampler->SetReferenceMesh( fixedMeshReader->GetOutput() ); // referenceMesh is copied to output
  resampler->SetInput( movingMeshReader->GetOutput() );        // this is where the scalars come from
  resampler->SetTransform( transform );
  resampler->SetInterpolator( interpolator );

  resampler->Update();

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MovingMeshType> RigidWriterType;
  RigidWriterType::Pointer rigidWriter = RigidWriterType::New();
  rigidWriter->SetInput( resampler->GetOutput() );
  rigidWriter->SetFileName(argv[6]);
  rigidWriter->Update();
  std::cout << "resampled mesh after rigid registration Saved" << std::endl;

  DemonsFilterType::Pointer demonsFilter = DemonsFilterType::New();

  typedef itk::ReplaceDestinationPointsQuadEdgeMeshFilter<
      FixedMeshType, DestinationPointSetType> ReplacePointsFilterType;

  ReplacePointsFilterType::Pointer replacePointsFilter = ReplacePointsFilterType::New();

  replacePointsFilter->SetInput( fixedMeshReader->GetOutput() );
  replacePointsFilter->SetDestinationPoints( destinationPoints );

  replacePointsFilter->Update();

  demonsFilter->SetFixedMesh( replacePointsFilter->GetOutput() );
  demonsFilter->SetMovingMesh( movingMeshReader->GetOutput() );

  DemonsFilterType::PointType center;
  center.Fill( 0.0 );

  const double radius = atof( argv[3] );

  demonsFilter->SetSphereCenter( center );
  demonsFilter->SetSphereRadius( radius );

  const double sigmaX =  atof( argv[4] ); // It should be proportional to inter-vertex distance
  const double epsilon = 1.0 / ( sigmaX * sigmaX );

  const double       lambda = 1.0;
  const unsigned int maximumNumberOfSmoothingIterations = 5;
  const unsigned int maximumNumberOfIterations = atoi( argv[5] );

  demonsFilter->SetEpsilon( epsilon );
  demonsFilter->SetSigmaX( sigmaX );

  // Internally refine values of SigmaX and Epsilon.
  demonsFilter->SelfRegulatedModeOn();

  demonsFilter->SetMaximumNumberOfIterations( maximumNumberOfIterations );

  demonsFilter->SetLambda( lambda );
  demonsFilter->SetMaximumNumberOfSmoothingIterations( maximumNumberOfSmoothingIterations );

  // Initialize the deformable registration stage with
  // the results of the Rigid Registration.
  // NOT EQUIVALENT demonsFilter->SetInitialDestinationPoints( destinationPoints );

#ifdef USE_VTK
  visualMonitor.SetBaseAnnotationText("Demons Registration Level 1");
  visualMonitor.Observe( demonsFilter.GetPointer() );
  visualMonitor.ObserveData( demonsFilter->GetFinalDestinationPoints() );
#endif

  try
    {
    demonsFilter->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  // assign resampled moving values to fixedMesh
  // inputMesh: fixedMeshReader->GetOutput()
  // sourceMesh: demonsFilter->GetOutput()
  typedef itk::AssignScalarValuesQuadEdgeMeshFilter<FixedMeshType, RegisteredMeshType, FixedMeshType> AssignFilterType;

  AssignFilterType::Pointer assignFilter = AssignFilterType::New();

  assignFilter->SetInputMesh(fixedMeshReader->GetOutput() );
  assignFilter->SetSourceMesh(demonsFilter->GetOutput() );

  assignFilter->Update();

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MovingMeshType> DemonsWriterType;
  DemonsWriterType::Pointer demonsWriter = DemonsWriterType::New();
  demonsWriter->SetFileName( argv[7] );
  demonsWriter->SetInput( assignFilter->GetOutput() );

  try
    {
    demonsWriter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  DestinationPointSetType::Pointer finalDestinationPoints = demonsFilter->GetFinalDestinationPoints();

  finalDestinationPoints->DisconnectPipeline();

  deformationFilter->SetInputMesh( fixedMeshReader->GetOutput() );
  deformationFilter->SetDestinationPoints( finalDestinationPoints );
  std::cout << "BEFORE Computing deformation field " << std::endl;
  deformationFilter->Update();
  std::cout << "AFTER Computing deformation field " << std::endl;

  vectorMeshWriter->SetInput( deformationFilter->GetOutput() );
  vectorMeshWriter->SetFileName(argv[9]);
  vectorMeshWriter->Update();
  std::cout << "Deformation Field Saved" << std::endl;

  demonsWriter->SetInput( demonsFilter->GetDeformedFixedMesh() );
  demonsWriter->SetFileName(argv[8]);
  demonsWriter->Update();
  std::cout << "Deformed Fixed Mesh  Saved" << std::endl;

  return EXIT_SUCCESS;
}
