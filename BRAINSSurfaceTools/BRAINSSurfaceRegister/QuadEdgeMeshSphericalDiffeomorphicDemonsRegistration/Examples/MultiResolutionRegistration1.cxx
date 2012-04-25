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
//
//
//                  WARNING:
//
//       THIS FILE HAS BEEN DEPRECATED
//
//  MOST OF THE CODE HAS BEEN REIMPLEMENTED IN
//  THE NEW MULTI-RESOLUTION CLASS IN THE SOURCE
//  DIRECTORY. PLEASE USE NOW THE FILE:
//
//     MultiResolutionRegistration2.cxx
//
//
//
//
//
//

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

#include "itkResampleQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshVectorDataVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshSphericalDiffeomorphicDemonsFilter.h"
#include "itkDeformationFieldFromTransformMeshFilter.h"
#include "itkResampleDestinationPointsQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshGenerateDeformationFieldFilter.h"
#include "itkIdentityTransform.h"
#include "itkReplaceDestinationPointsQuadEdgeMeshFilter.h"

#include "itkQuadEdgeMeshSphericalDiffeomorphicDemonsRegistrationConfigure.h"

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
  if( argc < 17 )
    {
    std::cerr << "Missing arguments" << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << std::endl;
    std::cerr << "inputFixedMeshRes1 inputMovingMeshRes1 ";
    std::cerr << "outputResampledMeshRes1 ";
    std::cerr << "inputFixedMeshRes2 inputMovingMeshRes2 ";
    std::cerr << "outputResampledMeshRes2 ";
    std::cerr << "inputFixedMeshRes3 inputMovingMeshRes3 ";
    std::cerr << "outputResampledMeshRes3 ";
    std::cerr << "inputFixedMeshRes4 inputMovingMeshRes4 ";
    std::cerr << "outputResampledMeshRes4 ";
    std::cerr << "inputFixedMeshRes5 inputMovingMeshRes5 ";
    std::cerr << "outputResampledMeshRes5 ";
    std::cerr << "screenshotsTag ";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  typedef float MeshPixelType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> FixedMeshType;
  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> MovingMeshType;

  typedef itk::VTKPolyDataReader<FixedMeshType>  FixedReaderType;
  typedef itk::VTKPolyDataReader<MovingMeshType> MovingReaderType;

  FixedReaderType::Pointer fixedMeshReader1 = FixedReaderType::New();
  fixedMeshReader1->SetFileName( argv[1] );

  MovingReaderType::Pointer movingMeshReader1 = MovingReaderType::New();
  movingMeshReader1->SetFileName( argv[2] );

  try
    {
    fixedMeshReader1->Update();
    movingMeshReader1->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  FixedMeshType::ConstPointer  fixedMesh1  = fixedMeshReader1->GetOutput();
  MovingMeshType::ConstPointer meshMoving = movingMeshReader1->GetOutput();

  typedef itk::MeshToMeshRegistrationMethod<
      FixedMeshType,
      MovingMeshType>    RegistrationType;

  RegistrationType::Pointer registration  = RegistrationType::New();

  typedef itk::MeanSquaresMeshToMeshMetric<FixedMeshType,
                                           MovingMeshType>
    MetricType;

  MetricType::Pointer metric = MetricType::New();

  registration->SetMetric( metric );

  registration->SetFixedMesh( fixedMesh1 );
  registration->SetMovingMesh( meshMoving );

  typedef itk::VersorTransform<MetricType::TransformComputationType> TransformType;

  TransformType::Pointer transform = TransformType::New();

  registration->SetTransform( transform );

  typedef itk::LinearInterpolateMeshFunction<MovingMeshType> InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  registration->SetInterpolator( interpolator );

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> RegisteredMeshType;

  typedef itk::QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<
      FixedMeshType, MovingMeshType, RegisteredMeshType>   DemonsFilterType;

  typedef DemonsFilterType::DestinationPointSetType DestinationPointSetType;

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
  optimizer->SetMaximumStepLength( 1e-2 );
  optimizer->SetMinimumStepLength( 1e-9 );
  optimizer->SetRelaxationFactor( 0.9 );
  optimizer->SetNumberOfIterations( 32 );

  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

#ifdef USE_VTK
  typedef DeformableAndAffineRegistrationMonitor<
      DemonsFilterType, DestinationPointSetType> RegistrationMonitorType;

  RegistrationMonitorType visualMonitor;
  visualMonitor.SetNumberOfIterationsPerUpdate( 1 );

  visualMonitor.SetBaseAnnotationText("Rigid Registration Level 1");

  visualMonitor.Observe( optimizer.GetPointer() );
  visualMonitor.ObserveData( transform.GetPointer() );

  visualMonitor.SetVerbose( false );
  visualMonitor.SetScreenShotsBaseFileName( argv[16] );

  visualMonitor.SetScalarRange( -1.0, 1.0 );

  vtkSmartPointer<vtkPolyDataReader> vtkFixedMeshReader =
    vtkSmartPointer<vtkPolyDataReader>::New();

  vtkSmartPointer<vtkPolyDataReader> vtkMovingMeshReader =
    vtkSmartPointer<vtkPolyDataReader>::New();

  vtkFixedMeshReader->SetFileName( fixedMeshReader1->GetFileName() );
  vtkMovingMeshReader->SetFileName( movingMeshReader1->GetFileName() );

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

  deformationFieldFromTransform->SetInput( fixedMeshReader1->GetOutput() );
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

  deformationFilter->SetInputMesh( fixedMeshReader1->GetOutput() );
  deformationFilter->SetDestinationPoints( destinationPoints );
  deformationFilter->Update();

  typedef itk::QuadEdgeMeshVectorDataVTKPolyDataWriter<MeshWithVectorsType> VectorMeshWriterType;
  VectorMeshWriterType::Pointer vectorMeshWriter = VectorMeshWriterType::New();
  vectorMeshWriter->SetInput( deformationFilter->GetOutput() );
  vectorMeshWriter->SetFileName("VectorMesh0.vtk");
  vectorMeshWriter->Update();
  std::cout << "Deformation VectorMesh0.vtk  Saved" << std::endl;

  DemonsFilterType::Pointer demonsFilter = DemonsFilterType::New();

  typedef itk::ReplaceDestinationPointsQuadEdgeMeshFilter<
      FixedMeshType, DestinationPointSetType> ReplacePointsFilterType;

  ReplacePointsFilterType::Pointer replacePointsFilter = ReplacePointsFilterType::New();

  replacePointsFilter->SetInput( fixedMeshReader1->GetOutput() );
  replacePointsFilter->SetDestinationPoints( destinationPoints );

  replacePointsFilter->Update();

  demonsFilter->SetFixedMesh( replacePointsFilter->GetOutput() );
  demonsFilter->SetMovingMesh( movingMeshReader1->GetOutput() );

  DemonsFilterType::PointType center;
  center.Fill( 0.0 );

  const double radius = 100.0;

  demonsFilter->SetSphereCenter( center );
  demonsFilter->SetSphereRadius( radius );

  const double sigmaXFactor = 2.5;

  // proportional to intervertex spacing of IC4 (at radius 100.0).
  double sigmaX = vcl_sqrt( radius / 10.0 ) * sigmaXFactor;
  double epsilon = 1.0 / (sigmaX * sigmaX );

  const double lambda = 1.0;

  std::cout << "SigmaX = " << sigmaX << " Epsilon = " << epsilon << std::endl;

  demonsFilter->SetEpsilon( epsilon );
  demonsFilter->SetSigmaX( sigmaX );

  // Internally refine values of SigmaX and Epsilon.
  demonsFilter->SelfRegulatedModeOn();

  demonsFilter->SetLambda( lambda );

  demonsFilter->SetMaximumNumberOfIterations( 15 );
  demonsFilter->SetMaximumNumberOfSmoothingIterations( 2 );

  // Initialize the deformable registration stage with
  // the results of the Rigid Registration.
// FIXME   demonsFilter->SetInitialDestinationPoints( destinationPoints );

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

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<FixedMeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetFileName( argv[3] );
  writer->SetInput( demonsFilter->GetOutput() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  DestinationPointSetType::Pointer finalDestinationPoints1 = demonsFilter->GetFinalDestinationPoints();

  finalDestinationPoints1->DisconnectPipeline();

  deformationFilter->SetInputMesh( fixedMeshReader1->GetOutput() );
  deformationFilter->SetDestinationPoints( finalDestinationPoints1 );
  std::cout << "BEFORE Computing deformation field " << std::endl;
  deformationFilter->Update();
  std::cout << "AFTER Computing deformation field " << std::endl;

  vectorMeshWriter->SetInput( deformationFilter->GetOutput() );
  vectorMeshWriter->SetFileName("VectorMesh1.vtk");
  vectorMeshWriter->Update();
  std::cout << "Deformation VectorMesh1.vtk  Saved" << std::endl;

  writer->SetInput( demonsFilter->GetDeformedFixedMesh() );
  writer->SetFileName("DeformedMesh1.vtk");
  writer->Update();
  std::cout << "Deformation DeformedMesh1.vtk  Saved" << std::endl;

  //  Starting process for the second resolution level (IC5).
  FixedReaderType::Pointer fixedMeshReader2 = FixedReaderType::New();
  fixedMeshReader2->SetFileName( argv[4] );

  MovingReaderType::Pointer movingMeshReader2 = MovingReaderType::New();
  movingMeshReader2->SetFileName( argv[5] );

  std::cout << "Fixed  mesh second level = " << argv[4] << std::endl;
  std::cout << "Moving mesh second level = " << argv[5] << std::endl;
  try
    {
    fixedMeshReader2->Update();
    movingMeshReader2->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  FixedMeshType::Pointer fixedMesh2 = fixedMeshReader2->GetOutput();
  fixedMesh2->DisconnectPipeline();

  // Supersample the list of destination points using the mesh at the next resolution level.
  typedef itk::ResampleDestinationPointsQuadEdgeMeshFilter<
      PointSetType, FixedMeshType, FixedMeshType, PointSetType> UpsampleDestinationPointsFilterType;

  UpsampleDestinationPointsFilterType::Pointer upsampleDestinationPoints =
    UpsampleDestinationPointsFilterType::New();

  upsampleDestinationPoints->SetSphereCenter( center );
  upsampleDestinationPoints->SetSphereRadius( radius );
  upsampleDestinationPoints->SetInput( finalDestinationPoints1 );
  upsampleDestinationPoints->SetFixedMesh( fixedMesh1 );
  upsampleDestinationPoints->SetReferenceMesh( fixedMesh2 );
  upsampleDestinationPoints->SetTransform( itk::IdentityTransform<double>::New() );

  try
    {
    std::cout << "BEFORE upsampleDestinationPoints Update()" << std::endl;
    upsampleDestinationPoints->Update();
    std::cout << "AFTER upsampleDestinationPoints Update()" << std::endl;
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  // Here build a Mesh using the upsampled destination points and
  // the scalar values of the fixed IC5 mesh.

  PointSetType::ConstPointer upsampledPointSet = upsampleDestinationPoints->GetOutput();

  const PointSetType::PointsContainer * upsampledPoints = upsampledPointSet->GetPoints();

  PointSetType::PointsContainerConstIterator upsampledPointsItr = upsampledPoints->Begin();
  PointSetType::PointsContainerConstIterator upsampledPointsEnd = upsampledPoints->End();

  FixedMeshType::PointsContainer::Pointer fixedPoints2 = fixedMesh2->GetPoints();

  FixedMeshType::PointsContainerIterator fixedPoint2Itr = fixedPoints2->Begin();
  std::cout << "UPSAMPLED POINTS BEGIN" << std::endl;
  while( upsampledPointsItr != upsampledPointsEnd )
    {
    // Point in the QuadEdgeMesh must also keep their pointer to Edge
    fixedPoint2Itr.Value().SetPoint( upsampledPointsItr.Value() );
    std::cout << upsampledPointsItr.Value() << std::endl;
    ++fixedPoint2Itr;
    ++upsampledPointsItr;
    }

  std::cout << "UPSAMPLED POINTS END" << std::endl;

  // Now feed this mesh into the Rigid registration of the second resolution level.

  registration->SetFixedMesh( fixedMesh2 );
  registration->SetMovingMesh( movingMeshReader2->GetOutput() );

  transform->SetIdentity();
  parameters = transform->GetParameters();

  registration->SetInitialTransformParameters( parameters );

  //   Save deformed fixed mesh, before initiating rigid registration
  writer->SetInput( fixedMesh2 );
  writer->SetFileName("DeformedMesh2BeforeRigidRegistration.vtk");
  writer->Update();

// #if NOT_VERIFYING_INTERMEDIATE_LEVELS

//// FIXME INCORRECT TRANSITION TO NEXT RESOLUTION BEGINS

#ifdef USE_VTK
  //  Reconnecting Registration monitor
  vtkFixedMeshReader->SetFileName("DeformedMesh2BeforeRigidRegistration.vtk");
  vtkMovingMeshReader->SetFileName( movingMeshReader2->GetFileName() );

  vtkFixedMeshReader->Update();
  vtkMovingMeshReader->Update();

  visualMonitor.SetBaseAnnotationText("Rigid Registration Level 2");
  visualMonitor.SetFixedSurface( vtkFixedMeshReader->GetOutput() );
  visualMonitor.SetMovingSurface( vtkMovingMeshReader->GetOutput() );
  visualMonitor.Observe( optimizer.GetPointer() );
  visualMonitor.ObserveData( transform.GetPointer() );
#endif

  //  Running Second Resolution Level Rigid Registration.

  std::cout << "Running Second Resolution Level Rigid Registration." << std::endl;

  optimizer->SetMaximumStepLength( 1e-2 );
  optimizer->SetMinimumStepLength( 1e-9 );
  optimizer->SetRelaxationFactor( 0.8 );
  optimizer->SetNumberOfIterations( 32 );

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

  finalParameters = registration->GetLastTransformParameters();

  std::cout << "final parameters = " << finalParameters << std::endl;
  std::cout << "final value      = " << optimizer->GetValue() << std::endl;

  transform->SetParameters( finalParameters );

  deformationFieldFromTransform->SetInput( fixedMesh2 );
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

//
// FIXME   demonsFilter->SetInitialDestinationPoints( deformationFieldFromTransform->GetOutput() );
//
  replacePointsFilter->SetInput( fixedMesh2 );
  replacePointsFilter->SetDestinationPoints( deformationFieldFromTransform->GetOutput() );

  replacePointsFilter->Update();

  demonsFilter->SetFixedMesh( replacePointsFilter->GetOutput() );
  demonsFilter->SetMovingMesh( movingMeshReader2->GetOutput() );

  demonsFilter->SetMaximumNumberOfIterations( 15 );
  demonsFilter->SetMaximumNumberOfSmoothingIterations( 5 );

#ifdef USE_VTK
  demonsFilter->MakeOutput(2);
  visualMonitor.SetBaseAnnotationText("Demons Registration Level 2");
  visualMonitor.Observe( demonsFilter.GetPointer() );
  visualMonitor.ObserveData( demonsFilter->GetFinalDestinationPoints() );
#endif

  epsilon = demonsFilter->GetEpsilon();
  sigmaX = demonsFilter->GetSigmaX();
  epsilon *= 125.0;
  sigmaX  /= 25.0;
  demonsFilter->SetEpsilon( epsilon );
  demonsFilter->SetSigmaX( sigmaX );

  //  Running Second Resolution Level Demons Registration.
  std::cout << "Running Second Resolution Level Demons Registration." << std::endl;

  try
    {
    demonsFilter->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  writer->SetFileName( argv[6] );
  writer->SetInput( demonsFilter->GetOutput() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  DestinationPointSetType::Pointer finalDestinationPoints2 = demonsFilter->GetFinalDestinationPoints();
  finalDestinationPoints2->DisconnectPipeline();

  deformationFilter->SetInputMesh( fixedMeshReader2->GetOutput() );
  deformationFilter->SetDestinationPoints( finalDestinationPoints2 );
  std::cout << "BEFORE Computing deformation field " << std::endl;
  deformationFilter->Update();
  std::cout << "AFTER Computing deformation field " << std::endl;

  vectorMeshWriter->SetInput( deformationFilter->GetOutput() );
  vectorMeshWriter->SetFileName("VectorMesh2.vtk");
  vectorMeshWriter->Update();
  std::cout << "Deformation VectorMesh2.vtk  Saved" << std::endl;

  writer->SetInput( demonsFilter->GetDeformedFixedMesh() );
  writer->SetFileName("DeformedMesh2.vtk");
  writer->Update();
  std::cout << "Deformation DeformedMesh2.vtk  Saved" << std::endl;

  //  Starting process for the Third resolution level (IC6).
  FixedReaderType::Pointer fixedMeshReader3 = FixedReaderType::New();
  fixedMeshReader3->SetFileName( argv[7] );

  MovingReaderType::Pointer movingMeshReader3 = MovingReaderType::New();
  movingMeshReader3->SetFileName( argv[8] );

  std::cout << "Fixed  mesh third level = " << argv[7] << std::endl;
  std::cout << "Moving mesh third level = " << argv[8] << std::endl;
  try
    {
    fixedMeshReader3->Update();
    movingMeshReader3->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  // Supersample the list of destination points using the mesh at the next resolution level.
  fixedMeshReader2->Update(); // go back to the original..
  upsampleDestinationPoints->SetInput( finalDestinationPoints2 );
  upsampleDestinationPoints->SetFixedMesh( fixedMeshReader2->GetOutput() ); // go back to the original.
  upsampleDestinationPoints->SetReferenceMesh( fixedMeshReader3->GetOutput() );
  upsampleDestinationPoints->SetTransform( itk::IdentityTransform<double>::New() );

  try
    {
    std::cout << "BEFORE upsampleDestinationPoints Update()" << std::endl;
    upsampleDestinationPoints->Update();
    std::cout << "AFTER upsampleDestinationPoints Update()" << std::endl;
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  // Here build a Mesh using the upsampled destination points and
  // the scalar values of the fixed IC6 mesh.

  FixedMeshType::Pointer fixedMesh3 = fixedMeshReader3->GetOutput();
  fixedMesh3->DisconnectPipeline();

  upsampledPointSet = upsampleDestinationPoints->GetOutput();

  upsampledPoints = upsampledPointSet->GetPoints();

  upsampledPointsItr = upsampledPoints->Begin();
  upsampledPointsEnd = upsampledPoints->End();

  FixedMeshType::PointsContainer::Pointer fixedPoints3 = fixedMesh3->GetPoints();

  FixedMeshType::PointsContainerIterator fixedPoint3Itr = fixedPoints3->Begin();

  while( upsampledPointsItr != upsampledPointsEnd )
    {
    fixedPoint3Itr.Value().SetPoint( upsampledPointsItr.Value() );
    ++fixedPoint3Itr;
    ++upsampledPointsItr;
    }

  // Now feed this mesh into the Rigid registration of the third resolution level.

  registration->SetFixedMesh( fixedMesh3 );
  registration->SetMovingMesh( movingMeshReader3->GetOutput() );

  transform->SetIdentity();
  parameters = transform->GetParameters();

  registration->SetInitialTransformParameters( parameters );

  //   Save deformed fixed mesh, before initiating rigid registration
  writer->SetInput( fixedMesh3 );
  writer->SetFileName("DeformedMesh3BeforeRigidRegistration.vtk");
  writer->Update();

//// FIXME INCORRECT TRANSITION TO NEXT RESOLUTION ENDS

#if NOT_VERIFYING_INTERMEDIATE_LEVELS

#ifdef USE_VTK
  //  Reconnecting Registration monitor
  vtkFixedMeshReader->SetFileName("DeformedMesh3BeforeRigidRegistration.vtk");
  vtkMovingMeshReader->SetFileName( movingMeshReader3->GetFileName() );

  vtkFixedMeshReader->Update();
  vtkMovingMeshReader->Update();

  visualMonitor.SetBaseAnnotationText("Rigid Registration Level 3");
  visualMonitor.SetFixedSurface( vtkFixedMeshReader->GetOutput() );
  visualMonitor.SetMovingSurface( vtkMovingMeshReader->GetOutput() );
  visualMonitor.Observe( optimizer.GetPointer() );
  visualMonitor.ObserveData( transform.GetPointer() );
#endif

  //  Running Third Resolution Level Rigid Registration.

  std::cout << "Running Third Resolution Level Rigid Registration." << std::endl;

  optimizer->SetMaximumStepLength( 1e-2 );
  optimizer->SetMinimumStepLength( 1e-9 );
  optimizer->SetRelaxationFactor( 0.6 );
  optimizer->SetNumberOfIterations( 16 );

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

  finalParameters = registration->GetLastTransformParameters();

  std::cout << "final parameters = " << finalParameters << std::endl;
  std::cout << "final value      = " << optimizer->GetValue() << std::endl;

  transform->SetParameters( finalParameters );

  deformationFieldFromTransform->SetInput( fixedMesh3 );
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

// FIXME   demonsFilter->SetInitialDestinationPoints( deformationFieldFromTransform->GetOutput() );

  replacePointsFilter->SetInput( fixedMesh3 );
  replacePointsFilter->SetDestinationPoints( deformationFieldFromTransform->GetOutput() );

  replacePointsFilter->Update();

  demonsFilter->SetFixedMesh( replacePointsFilter->GetOutput() );
  demonsFilter->SetMovingMesh( movingMeshReader3->GetOutput() );

  demonsFilter->SetMaximumNumberOfIterations( 15 );
  demonsFilter->SetMaximumNumberOfSmoothingIterations( 7 );

#ifdef USE_VTK
  demonsFilter->MakeOutput(2);
  visualMonitor.SetBaseAnnotationText("Demons Registration Level 3");
  visualMonitor.Observe( demonsFilter.GetPointer() );
  visualMonitor.ObserveData( demonsFilter->GetFinalDestinationPoints() );
#endif

  epsilon = demonsFilter->GetEpsilon();
  sigmaX = demonsFilter->GetSigmaX();
  epsilon *= 4.0;
  sigmaX  /= 2.0;
  demonsFilter->SetEpsilon( epsilon );
  demonsFilter->SetSigmaX( sigmaX );

  //  Running Third Resolution Level Demons Registration.
  std::cout << "Running Third Resolution Level Demons Registration." << std::endl;

  try
    {
    demonsFilter->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  writer->SetFileName( argv[9] );
  writer->SetInput( demonsFilter->GetOutput() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  DestinationPointSetType::Pointer finalDestinationPoints3 = demonsFilter->GetFinalDestinationPoints();
  finalDestinationPoints3->DisconnectPipeline();

  deformationFilter->SetInputMesh( fixedMeshReader3->GetOutput() );
  deformationFilter->SetDestinationPoints( finalDestinationPoints3 );
  std::cout << "BEFORE Computing deformation field " << std::endl;
  deformationFilter->Update();
  std::cout << "AFTER Computing deformation field " << std::endl;

  vectorMeshWriter->SetInput( deformationFilter->GetOutput() );
  vectorMeshWriter->SetFileName("VectorMesh3.vtk");
  vectorMeshWriter->Update();
  std::cout << "Deformation VectorMesh3.vtk  Saved" << std::endl;

  writer->SetInput( demonsFilter->GetDeformedFixedMesh() );
  writer->SetFileName("DeformedMesh3.vtk");
  writer->Update();
  std::cout << "Deformation DeformedMesh3.vtk  Saved" << std::endl;

  //  Starting process for the Fourth resolution level (IC7).
  FixedReaderType::Pointer fixedMeshReader4 = FixedReaderType::New();
  fixedMeshReader4->SetFileName( argv[10] );

  MovingReaderType::Pointer movingMeshReader4 = MovingReaderType::New();
  movingMeshReader4->SetFileName( argv[11] );

  std::cout << "Fixed  mesh fourth level = " << argv[10] << std::endl;
  std::cout << "Moving mesh fourth level = " << argv[11] << std::endl;
  try
    {
    fixedMeshReader4->Update();
    movingMeshReader4->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  // Supersample the list of destination points using the mesh at the next resolution level.
  upsampleDestinationPoints->SetInput( finalDestinationPoints3 );
  upsampleDestinationPoints->SetFixedMesh( fixedMesh3 );
  upsampleDestinationPoints->SetReferenceMesh( fixedMeshReader4->GetOutput() );
  upsampleDestinationPoints->SetTransform( itk::IdentityTransform<double>::New() );

  try
    {
    std::cout << "BEFORE upsampleDestinationPoints Update()" << std::endl;
    upsampleDestinationPoints->Update();
    std::cout << "AFTER upsampleDestinationPoints Update()" << std::endl;
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  // Here build a Mesh using the upsampled destination points and
  // the scalar values of the fixed IC6 mesh.

  FixedMeshType::Pointer fixedMesh4 = fixedMeshReader4->GetOutput();
  fixedMesh4->DisconnectPipeline();

  upsampledPointSet = upsampleDestinationPoints->GetOutput();

  upsampledPoints = upsampledPointSet->GetPoints();

  upsampledPointsItr = upsampledPoints->Begin();
  upsampledPointsEnd = upsampledPoints->End();

  FixedMeshType::PointsContainer::Pointer fixedPoints4 = fixedMesh4->GetPoints();

  FixedMeshType::PointsContainerIterator fixedPoint4Itr = fixedPoints4->Begin();

  while( upsampledPointsItr != upsampledPointsEnd )
    {
    fixedPoint4Itr.Value().SetPoint( upsampledPointsItr.Value() );
    ++fixedPoint4Itr;
    ++upsampledPointsItr;
    }

  // Now feed this mesh into the Rigid registration of the third resolution level.

  registration->SetFixedMesh( fixedMesh4 );
  registration->SetMovingMesh( movingMeshReader4->GetOutput() );

  transform->SetIdentity();
  parameters = transform->GetParameters();

  registration->SetInitialTransformParameters( parameters );

  //   Save deformed fixed mesh, before initiating rigid registration
  writer->SetInput( fixedMesh4 );
  writer->SetFileName("DeformedMesh4BeforeRigidRegistration.vtk");
  writer->Update();

// #if NOT_VERIFYING_INTERMEDIATE_LEVELS

#ifdef USE_VTK
  //  Reconnecting Registration monitor
  vtkFixedMeshReader->SetFileName("DeformedMesh4BeforeRigidRegistration.vtk");
  vtkMovingMeshReader->SetFileName( movingMeshReader4->GetFileName() );

  vtkFixedMeshReader->Update();
  vtkMovingMeshReader->Update();

  visualMonitor.SetBaseAnnotationText("Rigid Registration Level 4");
  visualMonitor.SetFixedSurface( vtkFixedMeshReader->GetOutput() );
  visualMonitor.SetMovingSurface( vtkMovingMeshReader->GetOutput() );
  visualMonitor.Observe( optimizer.GetPointer() );
  visualMonitor.ObserveData( transform.GetPointer() );
#endif

  //  Running Fourth Resolution Level Rigid Registration.

  std::cout << "Running Fourth Resolution Level Rigid Registration." << std::endl;

  optimizer->SetMaximumStepLength( 1e-2 );
  optimizer->SetMinimumStepLength( 1e-9 );
  optimizer->SetRelaxationFactor( 0.5 );
  optimizer->SetNumberOfIterations( 8 );

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

  finalParameters = registration->GetLastTransformParameters();

  std::cout << "final parameters = " << finalParameters << std::endl;
  std::cout << "final value      = " << optimizer->GetValue() << std::endl;

  transform->SetParameters( finalParameters );

  deformationFieldFromTransform->SetInput( fixedMesh4 );
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

//
// FIXME   demonsFilter->SetInitialDestinationPoints( deformationFieldFromTransform->GetOutput() );
//
  replacePointsFilter->SetInput( fixedMesh4 );
  replacePointsFilter->SetDestinationPoints( deformationFieldFromTransform->GetOutput() );

  replacePointsFilter->Update();

  demonsFilter->SetFixedMesh( replacePointsFilter->GetOutput() );
  demonsFilter->SetMovingMesh( movingMeshReader4->GetOutput() );

  demonsFilter->SetMaximumNumberOfIterations( 15 );
  demonsFilter->SetMaximumNumberOfSmoothingIterations( 10 );

#ifdef USE_VTK
  demonsFilter->MakeOutput(2);
  visualMonitor.SetBaseAnnotationText("Demons Registration Level 4");
  visualMonitor.Observe( demonsFilter.GetPointer() );
  visualMonitor.ObserveData( demonsFilter->GetFinalDestinationPoints() );
#endif

  epsilon = demonsFilter->GetEpsilon();
  sigmaX = demonsFilter->GetSigmaX();
  epsilon *= 4.0;
  sigmaX  /= 2.0;
  demonsFilter->SetEpsilon( epsilon );
  demonsFilter->SetSigmaX( sigmaX );

  //  Running Fourth Resolution Level Demons Registration.
  std::cout << "Running Fourth Resolution Level Demons Registration." << std::endl;

  try
    {
    demonsFilter->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  writer->SetFileName( argv[12] );
  writer->SetInput( demonsFilter->GetOutput() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  DestinationPointSetType::Pointer finalDestinationPoints4 = demonsFilter->GetFinalDestinationPoints();
  finalDestinationPoints4->DisconnectPipeline();

  deformationFilter->SetInputMesh( fixedMeshReader4->GetOutput() );
  deformationFilter->SetDestinationPoints( finalDestinationPoints4 );
  std::cout << "BEFORE Computing deformation field " << std::endl;
  deformationFilter->Update();
  std::cout << "AFTER Computing deformation field " << std::endl;

  vectorMeshWriter->SetInput( deformationFilter->GetOutput() );
  vectorMeshWriter->SetFileName("VectorMesh4.vtk");
  vectorMeshWriter->Update();
  std::cout << "Deformation VectorMesh4.vtk  Saved" << std::endl;

  writer->SetFileName("DeformedMesh4.vtk");
  writer->SetInput( demonsFilter->GetDeformedFixedMesh() );
  writer->Update();
  std::cout << "Deformation DeformedMesh4.vtk  Saved" << std::endl;

#endif // NOT_VERIFYING_INTERMEDIATE_LEVELS

  //  Starting process for the Fifth resolution level (FINAL).
  FixedReaderType::Pointer fixedMeshReader5 = FixedReaderType::New();
  fixedMeshReader5->SetFileName( argv[13] );

  MovingReaderType::Pointer movingMeshReader5 = MovingReaderType::New();
  movingMeshReader5->SetFileName( argv[14] );

  std::cout << "Fixed  mesh fifth level = " << argv[13] << std::endl;
  std::cout << "Moving mesh fifth level = " << argv[14] << std::endl;
  try
    {
    fixedMeshReader5->Update();
    movingMeshReader5->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  // Supersample the list of destination points using the mesh at the next resolution level.
#if NOT_VERIFYING_INTERMEDIATE_LEVELS
  upsampleDestinationPoints->SetInput( finalDestinationPoints4 );
  upsampleDestinationPoints->SetFixedMesh( fixedMesh4 );
#else
  std::cout << "finalDestinationPoints size = " << finalDestinationPoints2->GetNumberOfPoints() << std::endl;
  std::cout << "fixedMesh              size = " << fixedMesh2->GetNumberOfPoints() << std::endl;
  upsampleDestinationPoints->SetInput( finalDestinationPoints2 );
  upsampleDestinationPoints->SetFixedMesh( fixedMesh2 );
#endif
  upsampleDestinationPoints->SetReferenceMesh( fixedMeshReader5->GetOutput() );
  upsampleDestinationPoints->SetTransform( itk::IdentityTransform<double>::New() );

  try
    {
    std::cout << "BEFORE upsampleDestinationPoints Update()" << std::endl;
    upsampleDestinationPoints->Update();
    std::cout << "AFTER upsampleDestinationPoints Update()" << std::endl;
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  // Here build a Mesh using the upsampled destination points and
  // the scalar values of the fixed FINAL mesh.

  FixedMeshType::Pointer fixedMesh5 = fixedMeshReader5->GetOutput();
  fixedMesh5->DisconnectPipeline();

  upsampledPointSet = upsampleDestinationPoints->GetOutput();

  deformationFilter->SetInputMesh( fixedMeshReader5->GetOutput() );
  deformationFilter->SetDestinationPoints( upsampleDestinationPoints->GetOutput() );
  std::cout << "BEFORE Computing deformation field " << std::endl;
  deformationFilter->Update();
  std::cout << "AFTER Computing deformation field " << std::endl;

  vectorMeshWriter->SetInput( deformationFilter->GetOutput() );
  vectorMeshWriter->SetFileName("VectorMesh5.vtk");
  vectorMeshWriter->Update();
  std::cout << "Deformation VectorMesh5.vtk  Saved" << std::endl;

  upsampledPoints = upsampledPointSet->GetPoints();

  upsampledPointsItr = upsampledPoints->Begin();
  upsampledPointsEnd = upsampledPoints->End();

  FixedMeshType::PointsContainer::Pointer fixedPoints5 = fixedMesh5->GetPoints();

  FixedMeshType::PointsContainerIterator fixedPoint5Itr = fixedPoints5->Begin();

  while( upsampledPointsItr != upsampledPointsEnd )
    {
    fixedPoint5Itr.Value().SetPoint( upsampledPointsItr.Value() );
    ++fixedPoint5Itr;
    ++upsampledPointsItr;
    }

  writer->SetFileName( argv[15] );
  writer->SetInput( fixedMesh5 );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  writer->SetFileName("DeformedMesh5.vtk");
  writer->SetInput( fixedMesh5 );
  writer->Update();
  std::cout << "Deformation DeformedMesh5.vtk  Saved" << std::endl;

  return EXIT_SUCCESS;
}
