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

#include "itkQuadEdgeMeshSphericalDiffeomorphicDemonsRegistrationConfigure.h"

#include "itkCommand.h"

#include "itkVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

#include "itkMultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter.h"
#include "itkDeformQuadEdgeMeshFilter.h"

#ifdef USE_VTK
#include "MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#endif

int main( int argc, char * argv [] )
{
  if( argc < 15 )
    {
    std::cerr << "Missing arguments" << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << std::endl;
    std::cerr << "inputFixedMeshRes1 inputMovingMeshRes1 ";
    std::cerr << "inputFixedMeshRes2 inputMovingMeshRes2 ";
    std::cerr << "inputFixedMeshRes3 inputMovingMeshRes3 ";
    std::cerr << "inputFixedMeshRes4 inputMovingMeshRes4 ";
    std::cerr << "outputResampledMeshRes4 ";
    std::cerr << "inputFixedMeshRes5 inputMovingMeshRes5 ";
    std::cerr << "outputResampledMeshRes5 ";
    std::cerr << "screenshotsTag ";
    std::cerr << "targetMesh ";
    std::cerr << "numberOfLevelsToUse ";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  typedef float MeshPixelType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> MeshType;

  typedef itk::MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<
      MeshType>  MultiResolutionDemonsFilterType;

  MultiResolutionDemonsFilterType::Pointer multiResDemonsFilter =
    MultiResolutionDemonsFilterType::New();

  typedef MultiResolutionDemonsFilterType::DestinationPointSetType DestinationPointSetType;

  const unsigned int maximumNumberOfResolutions = 4;
  const unsigned int finestResolution = atoi( argv[15] );

  multiResDemonsFilter->SetNumberOfResolutionLevels( finestResolution );

  typedef MultiResolutionDemonsFilterType::IntegerArrayType IntegerArrayType;
  typedef MultiResolutionDemonsFilterType::DoubleArrayType  DoubleArrayType;

  IntegerArrayType smoothingIterations(maximumNumberOfResolutions);

  smoothingIterations[0] = 10;
  smoothingIterations[1] = 10;
  smoothingIterations[2] = 10;
  smoothingIterations[3] = 10;

  multiResDemonsFilter->SetSmoothingIterations( smoothingIterations );

  IntegerArrayType demonsIterations(maximumNumberOfResolutions);

  demonsIterations[0] = 15;
  demonsIterations[1] = 15;
  demonsIterations[2] = 15;
  demonsIterations[3] = 15;

  multiResDemonsFilter->SetDemonsIterations( demonsIterations );

  IntegerArrayType rigidIterations(maximumNumberOfResolutions);

  rigidIterations[0] = 32;
  rigidIterations[1] = 32;
  rigidIterations[2] = 32;
  rigidIterations[3] = 32;

  multiResDemonsFilter->SetRigidRegistrationIterations( rigidIterations );

  DoubleArrayType rigidStepLength(maximumNumberOfResolutions);

  rigidStepLength[0] = 1e-2;
  rigidStepLength[1] = rigidStepLength[0] / 2.0;
  rigidStepLength[2] = rigidStepLength[1] / 2.0;
  rigidStepLength[3] = rigidStepLength[2] / 2.0;

  multiResDemonsFilter->SetRigidRegistrationStepLength( rigidStepLength );

  const double shortestEdgeLengthAtIC4 = 6.91822;

  DoubleArrayType epsilon(maximumNumberOfResolutions);

  epsilon[0] = 1.0 / ( 8.0 * shortestEdgeLengthAtIC4 );
  epsilon[1] = epsilon[0] * 2.0;
  epsilon[2] = epsilon[1] * 2.0;
  epsilon[3] = epsilon[2] * 2.0;

  multiResDemonsFilter->SetEpsilonValues( epsilon );

  DoubleArrayType sigmaX(maximumNumberOfResolutions);

  sigmaX[0] = 1.0 / vcl_sqrt( epsilon[0] );
  sigmaX[1] = 1.0 / vcl_sqrt( epsilon[1] );
  sigmaX[2] = 1.0 / vcl_sqrt( epsilon[2] );
  sigmaX[3] = 1.0 / vcl_sqrt( epsilon[3] );

  multiResDemonsFilter->SetSigmaXValues( sigmaX );

  typedef itk::VTKPolyDataReader<MeshType> ReaderType;

  ReaderType::Pointer fixedMeshReader[maximumNumberOfResolutions];
  ReaderType::Pointer movingMeshReader[maximumNumberOfResolutions];
  for( unsigned int i = 0; i < maximumNumberOfResolutions; i++ )
    {
    fixedMeshReader[i] = ReaderType::New();
    movingMeshReader[i] = ReaderType::New();

    fixedMeshReader[i]->SetFileName( argv[2 * i + 1] );
    movingMeshReader[i]->SetFileName( argv[2 * i + 2] );

    std::cout << "Fixed   " << argv[2 * i + 1] << std::endl;
    std::cout << "Moving  " << argv[2 * i + 2] << std::endl;

    fixedMeshReader[i]->Update();
    movingMeshReader[i]->Update();

    multiResDemonsFilter->SetFixedMesh( i, fixedMeshReader[i]->GetOutput() );
    multiResDemonsFilter->SetMovingMesh( i, movingMeshReader[i]->GetOutput() );
    }

  ReaderType::Pointer fixedMeshReader5 = ReaderType::New();
  fixedMeshReader5->SetFileName( argv[10] );

  ReaderType::Pointer movingMeshReader5 = ReaderType::New();
  movingMeshReader5->SetFileName( argv[11] );

  ReaderType::Pointer fixedMeshReader6 = ReaderType::New();
  fixedMeshReader6->SetFileName( argv[14] );

  try
    {
    fixedMeshReader5->Update();
    movingMeshReader5->Update();
    fixedMeshReader6->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

#ifdef USE_VTK
  typedef MultiResolutionDeformableAndAffineRegistrationMonitorWithTargetTracking<
      MultiResolutionDemonsFilterType, DestinationPointSetType> RegistrationMonitorType;

  RegistrationMonitorType visualMonitor;

  visualMonitor.SetMultiResolutionDemonsFilter( multiResDemonsFilter );

  visualMonitor.SetNumberOfIterationsPerUpdate( 1 );

  visualMonitor.SetBaseAnnotationText("Rigid Registration Level 1");

  visualMonitor.SetVerbose( true );
  visualMonitor.SetScreenShotsBaseFileName( argv[13] );

  visualMonitor.SetScalarRange( -1.0, 1.0 );

  visualMonitor.SetNumberOfResolutionLevels( maximumNumberOfResolutions );

  vtkSmartPointer<vtkPolyDataReader> vtkFixedMeshReader[maximumNumberOfResolutions];
  vtkSmartPointer<vtkPolyDataReader> vtkMovingMeshReader[maximumNumberOfResolutions];
  for( unsigned int i = 0; i < maximumNumberOfResolutions; i++ )
    {
    vtkFixedMeshReader[i] = vtkSmartPointer<vtkPolyDataReader>::New();
    vtkMovingMeshReader[i] = vtkSmartPointer<vtkPolyDataReader>::New();

    vtkFixedMeshReader[i]->SetFileName( fixedMeshReader[i]->GetFileName() );
    vtkMovingMeshReader[i]->SetFileName( movingMeshReader[i]->GetFileName() );

    vtkFixedMeshReader[i]->Update();
    vtkMovingMeshReader[i]->Update();

    visualMonitor.SetFixedSurface( i, vtkFixedMeshReader[i]->GetOutput() );
    visualMonitor.SetMovingSurface( i, vtkMovingMeshReader[i]->GetOutput() );
    }

  visualMonitor.SetFixedMeshSource( fixedMeshReader5->GetOutput() );
  visualMonitor.SetFixedMeshTarget( fixedMeshReader6->GetOutput() );
  visualMonitor.SetEvaluateDistanceToTarget(true);

#endif

  MultiResolutionDemonsFilterType::PointType center;
  center.Fill( 0.0 );

  const double radius = 100.0;

  multiResDemonsFilter->SetSphereCenter( center );
  multiResDemonsFilter->SetSphereRadius( radius );

  try
    {
    multiResDemonsFilter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::DeformQuadEdgeMeshFilter<
      MeshType, MeshType, DestinationPointSetType>  DeformMeshFilterType;

  DeformMeshFilterType::Pointer deformFilter = DeformMeshFilterType::New();

  deformFilter->SetInput( fixedMeshReader[finestResolution - 1]->GetOutput() );
  deformFilter->SetReferenceMesh( fixedMeshReader[finestResolution - 1]->GetOutput() );

  deformFilter->SetDestinationPoints( multiResDemonsFilter->GetFinalDestinationPoints() );

  deformFilter->SetSphereRadius( radius );
  deformFilter->SetSphereCenter( center );

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetFileName( argv[9] );
  writer->SetInput( deformFilter->GetOutput() );

  try
    {
    deformFilter->Update();
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  deformFilter->SetInput( fixedMeshReader5->GetOutput() );
  writer->SetFileName( argv[12] );

  try
    {
    deformFilter->Update();
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
