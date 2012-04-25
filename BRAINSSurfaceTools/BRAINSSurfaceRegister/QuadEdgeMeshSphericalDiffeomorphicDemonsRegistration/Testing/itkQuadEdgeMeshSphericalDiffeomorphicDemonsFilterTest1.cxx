/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ImageRegistration8.cxx,v $
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

#include "itkQuadEdgeMeshSphericalDiffeomorphicDemonsFilter.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkMeshGeneratorHelper.h"
#include "itkMeshWriterHelper1.h"
#include "itkTestingMacros.h"
#include "itkFilterWatcher.h"

int main( int argc, char *argv[] )
{
  if( argc < 9 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedMeshFile  movingMeshFile ";
    std::cerr << " outputMeshfile epsilon sigmaX ";
    std::cerr << " lambda smoothingIterations ";
    std::cerr << " numberOfIterations" << std::endl;
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

  std::cout << demonsFilter->GetNameOfClass() << std::endl;
  demonsFilter->Print( std::cout );

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

  if( demonsFilter->GetFixedMesh() != fixedReader->GetOutput() )
    {
    std::cerr << "Error in SetFixedMesh()/GetFixedMesh() " << std::endl;
    return EXIT_FAILURE;
    }

  if( demonsFilter->GetMovingMesh() != movingReader->GetOutput() )
    {
    std::cerr << "Error in SetMovingMesh()/GetMovingMesh() " << std::endl;
    return EXIT_FAILURE;
    }

  DemonsFilterType::PointType center1;
  center1.Fill( 17.0 );
  demonsFilter->SetSphereCenter( center1 );
  TEST_SET_GET_VALUE( center1, demonsFilter->GetSphereCenter() );

  DemonsFilterType::PointType center2;
  center2.Fill( 31.0 );
  demonsFilter->SetSphereCenter( center2 );
  TEST_SET_GET_VALUE( center2, demonsFilter->GetSphereCenter() );

  const double radius1 = 19.0;
  demonsFilter->SetSphereRadius( radius1 );
  TEST_SET_GET_VALUE( radius1, demonsFilter->GetSphereRadius() );

  const double radius2 = 29.0;
  demonsFilter->SetSphereRadius( radius2 );
  TEST_SET_GET_VALUE( radius2, demonsFilter->GetSphereRadius() );

  DemonsFilterType::PointType center0;
  center0.Fill( 0.0 );

  demonsFilter->SetSphereCenter( center0 );
  demonsFilter->SetSphereRadius( 100.0 );

  const double epsilon1 = 2.0;
  demonsFilter->SetEpsilon( epsilon1 );
  TEST_SET_GET_VALUE( epsilon1, demonsFilter->GetEpsilon() );

  const double epsilon2 = 3.0;
  demonsFilter->SetEpsilon( epsilon2 );
  TEST_SET_GET_VALUE( epsilon2, demonsFilter->GetEpsilon() );

  const double epsilon = atof( argv[4] );
  demonsFilter->SetEpsilon( epsilon );

  const double sigmaX1 = 2.0;
  demonsFilter->SetSigmaX( sigmaX1 );
  TEST_SET_GET_VALUE( sigmaX1, demonsFilter->GetSigmaX() );

  const double sigmaX2 = 3.0;
  demonsFilter->SetSigmaX( sigmaX2 );
  TEST_SET_GET_VALUE( sigmaX2, demonsFilter->GetSigmaX() );

  const double sigmaX = atof( argv[5] );
  demonsFilter->SetSigmaX( sigmaX );

  const double lambda1 = 2.0;
  demonsFilter->SetLambda( lambda1 );
  TEST_SET_GET_VALUE( lambda1, demonsFilter->GetLambda() );

  const double lambda2 = 3.0;
  demonsFilter->SetLambda( lambda2 );
  TEST_SET_GET_VALUE( lambda2, demonsFilter->GetLambda() );

  const double lambda = atof( argv[6] );
  demonsFilter->SetLambda( lambda );

  unsigned int maximumNumberOfSmoothingIterations = 10;
  demonsFilter->SetMaximumNumberOfSmoothingIterations( maximumNumberOfSmoothingIterations );
  TEST_SET_GET_VALUE( maximumNumberOfSmoothingIterations, demonsFilter->GetMaximumNumberOfSmoothingIterations() );

  maximumNumberOfSmoothingIterations = 15;
  demonsFilter->SetMaximumNumberOfSmoothingIterations( maximumNumberOfSmoothingIterations );
  TEST_SET_GET_VALUE( maximumNumberOfSmoothingIterations, demonsFilter->GetMaximumNumberOfSmoothingIterations() );

  maximumNumberOfSmoothingIterations = atoi( argv[7] );
  demonsFilter->SetMaximumNumberOfSmoothingIterations( maximumNumberOfSmoothingIterations );

  unsigned int maximumNumberOfIterations = 10;
  demonsFilter->SetMaximumNumberOfIterations( maximumNumberOfIterations );
  TEST_SET_GET_VALUE( maximumNumberOfIterations, demonsFilter->GetMaximumNumberOfIterations() );

  maximumNumberOfIterations = 15;
  demonsFilter->SetMaximumNumberOfIterations( maximumNumberOfIterations );
  TEST_SET_GET_VALUE( maximumNumberOfIterations, demonsFilter->GetMaximumNumberOfIterations() );

  maximumNumberOfIterations = atoi( argv[8] );
  demonsFilter->SetMaximumNumberOfIterations( maximumNumberOfIterations );

  FilterWatcher watcher( demonsFilter, "Demons Filter");

  try
    {
    demonsFilter->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while running the Demons filter " << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  itk::MeshWriterHelper1<RegisteredMeshType>::WriteMeshToFile( demonsFilter->GetOutput(), argv[3] );

  typedef DemonsFilterType::DestinationPointSetType DestinationPointSetType;
  typedef DestinationPointSetType::ConstPointer     DestinationPointContainerConstPointer;

  DestinationPointContainerConstPointer destinationPoints = demonsFilter->GetFinalDestinationPoints();

  std::cout << "Number of Destination points = " << destinationPoints->GetNumberOfPoints() << std::endl;

  typedef DemonsFilterType::BasisSystemContainerType BasisSystemContainerType;
  typedef BasisSystemContainerType::ConstPointer     BasisSystemContainerPointer;

  BasisSystemContainerPointer basisSystems = demonsFilter->GetBasisSystemAtNode();

  std::cout << "Number of basis systems = " << basisSystems->Size() << std::endl;

  std::cout << std::endl;
  std::cout << "Testing second call to Update(), the filter shouldn't run again" << std::endl;
  demonsFilter->Update();

  return EXIT_SUCCESS;
}
