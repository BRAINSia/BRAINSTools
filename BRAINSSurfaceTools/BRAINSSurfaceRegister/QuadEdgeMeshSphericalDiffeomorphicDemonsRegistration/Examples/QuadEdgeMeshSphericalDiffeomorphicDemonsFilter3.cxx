/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: QuadEdgeMeshSphericalDiffeomorphicDemonsFilter3.cxx,v $
  Language:  C++
  Date:      $Date: 2010-09-06 11:46:31 $
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

#include "itkQuadEdgeMeshSphericalDiffeomorphicDemonsFilter.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMesh.h"

#include "itkQuadEdgeMeshSphericalDiffeomorphicDemonsRegistrationConfigure.h"

#ifdef USE_VTK
#include "DeformableRegistrationMonitor.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#endif

int main( int argc, char *argv[] )
{
  if( argc < 11 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedMeshFile  movingMeshFile ";
    std::cerr << " outputResampledMovingMeshfile ";
    std::cerr << " outputDeformedFixedMeshfile ";
    std::cerr << " sphereRadius epsilon sigmaX ";
    std::cerr << " lambda smoothingIterations ";
    std::cerr << " numberOfIterations ";
    std::cerr << " MetricSignificance" << std::endl;
    return EXIT_FAILURE;
    }

  typedef float MeshPixelType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> FixedMeshType;
  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> MovingMeshType;
  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> RegisteredMeshType;

  typedef itk::QuadEdgeMeshSphericalDiffeomorphicDemonsFilter<
      FixedMeshType, MovingMeshType, RegisteredMeshType>   DemonsFilterType;

  DemonsFilterType::Pointer demonsFilter = DemonsFilterType::New();

  typedef itk::QuadEdgeMeshVTKPolyDataReader<FixedMeshType>  FixedReaderType;
  typedef itk::QuadEdgeMeshVTKPolyDataReader<MovingMeshType> MovingReaderType;

  FixedReaderType::Pointer fixedReader = FixedReaderType::New();
  fixedReader->SetFileName( argv[1] );

  MovingReaderType::Pointer movingReader = MovingReaderType::New();
  movingReader->SetFileName( argv[2] );

  try
    {
    fixedReader->Update();
    movingReader->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while reading the input file " << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  demonsFilter->SetFixedMesh( fixedReader->GetOutput() );
  demonsFilter->SetMovingMesh( movingReader->GetOutput() );

  DemonsFilterType::PointType center;
  center.Fill( 0.0 );

  const double radius = atof( argv[5] );

  demonsFilter->SetSphereCenter( center );
  demonsFilter->SetSphereRadius( radius );

  const double       epsilon = atof( argv[6] );
  const double       sigmaX = atof( argv[7] );
  const double       lambda = atof( argv[8] );
  const unsigned int maximumNumberOfSmoothingIterations = atoi( argv[9] );
  const unsigned int maximumNumberOfIterations = atoi( argv[10] );
  const double       metricSignificance = atof( argv[11] );

  demonsFilter->SetEpsilon( epsilon );
  demonsFilter->SetSigmaX( sigmaX );

  // Internally refine values of SigmaX and Epsilon.
  demonsFilter->SelfRegulatedModeOn();

  demonsFilter->SetMaximumNumberOfIterations( maximumNumberOfIterations );

  demonsFilter->SetLambda( lambda );
  demonsFilter->SetMaximumNumberOfSmoothingIterations( maximumNumberOfSmoothingIterations );

  demonsFilter->SelfStopModeOn();
  demonsFilter->SetMetricSignificance( metricSignificance );

#ifdef USE_VTK
  typedef DemonsFilterType::DestinationPointSetType DestinationPointSetType;

  typedef DeformableRegistrationMonitor<
      DemonsFilterType, DestinationPointSetType> RegistrationMonitorType;

  RegistrationMonitorType visualMonitor;
  visualMonitor.SetNumberOfIterationsPerUpdate( 1 );

  visualMonitor.SetBaseAnnotationText("Demons Registration ");

  // visualMonitor.SetVerbose( false );
  visualMonitor.SetScreenShotsBaseFileName( "demonsExample3" );

  visualMonitor.SetScalarRange( -1.0, 1.0 );

  visualMonitor.SetCameraZoomFactor( 1.5 );
  visualMonitor.SetCameraAzimuthAngle( 90.0 );
  visualMonitor.SetCameraElevationAngle( 0.0 );

  vtkSmartPointer<vtkPolyDataReader> vtkFixedMeshReader =
    vtkSmartPointer<vtkPolyDataReader>::New();

  vtkSmartPointer<vtkPolyDataReader> vtkMovingMeshReader =
    vtkSmartPointer<vtkPolyDataReader>::New();

  vtkFixedMeshReader->SetFileName( fixedReader->GetFileName() );
  vtkMovingMeshReader->SetFileName( movingReader->GetFileName() );

  vtkFixedMeshReader->Update();
  vtkMovingMeshReader->Update();

  visualMonitor.SetFixedSurface( vtkFixedMeshReader->GetOutput() );
  visualMonitor.SetMovingSurface( vtkMovingMeshReader->GetOutput() );

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

  writer->SetFileName( argv[4] );
  writer->SetInput( demonsFilter->GetDeformedFixedMesh() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
