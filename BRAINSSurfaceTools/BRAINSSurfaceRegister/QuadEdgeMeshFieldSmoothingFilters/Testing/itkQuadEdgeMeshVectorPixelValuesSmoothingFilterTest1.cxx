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

#include "itkQuadEdgeMeshVectorPixelValuesSmoothingFilter.h"
#include "itkQuadEdgeMesh.h"
#include "itkMeshWriterHelper2.h"
#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkTestingMacros.h"
#include "itkFilterWatcher.h"

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputMeshFile  outputMeshfile lambda iterations" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef itk::Vector<float, Dimension> MeshPixelType;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> InputMeshType;
  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> OutputMeshType;

  typedef itk::QuadEdgeMeshVectorPixelValuesSmoothingFilter<
      InputMeshType, OutputMeshType>                       FilterType;

  FilterType::Pointer filter = FilterType::New();

  std::cout << filter->GetNameOfClass() << std::endl;
  filter->Print( std::cout );

  typedef itk::QuadEdgeMeshVTKPolyDataReader<InputMeshType> ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while reading the input file " << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  filter->SetInput( reader->GetOutput() );

  if( filter->GetInput() != reader->GetOutput() )
    {
    std::cerr << "Error in SetFixedMesh()/GetFixedMesh() " << std::endl;
    return EXIT_FAILURE;
    }

  unsigned int maximumNumberOfIterations = 10;
  filter->SetMaximumNumberOfIterations( maximumNumberOfIterations );
  TEST_SET_GET_VALUE( maximumNumberOfIterations, filter->GetMaximumNumberOfIterations() );

  maximumNumberOfIterations = 15;
  filter->SetMaximumNumberOfIterations( maximumNumberOfIterations );
  TEST_SET_GET_VALUE( maximumNumberOfIterations, filter->GetMaximumNumberOfIterations() );

  double lambda = 2.0;
  filter->SetLambda( lambda );
  TEST_SET_GET_VALUE( lambda, filter->GetLambda() );

  lambda = 3.0;
  filter->SetLambda( lambda );
  TEST_SET_GET_VALUE( lambda, filter->GetLambda() );

  FilterType::OutputPointType center1;
  center1.Fill( 17.0 );
  filter->SetSphereCenter( center1 );
  TEST_SET_GET_VALUE( center1, filter->GetSphereCenter() );

  FilterType::OutputPointType center2;
  center2.Fill( 31.0 );
  filter->SetSphereCenter( center2 );
  TEST_SET_GET_VALUE( center2, filter->GetSphereCenter() );

  const double radius1 = 19.0;
  filter->SetSphereRadius( radius1 );
  TEST_SET_GET_VALUE( radius1, filter->GetSphereRadius() );

  const double radius2 = 29.0;
  filter->SetSphereRadius( radius2 );
  TEST_SET_GET_VALUE( radius2, filter->GetSphereRadius() );

  FilterType::OutputPointType center0;
  center0.Fill( 0.0 );

  filter->SetSphereCenter( center0 );
  filter->SetSphereRadius( 100.0 );

  FilterWatcher watcher( filter, "Smoothing Filter");

  filter->SetLambda( atof( argv[3] ) );
  filter->SetMaximumNumberOfIterations( atoi( argv[4] ) );

  TRY_EXPECT_NO_EXCEPTION( filter->Update() );

  itk::MeshWriterHelper2<OutputMeshType>::WriteMeshToFile( filter->GetOutput(), argv[2] );

#ifdef TEST_EXCEPTIONS
  InputMeshType::Pointer inputMesh = reader->GetOutput();
  inputMesh->DisconnectPipeline();

  reader = NULL; // Destroy the reader

  InputMeshType::PointsContainerPointer    points = inputMesh->GetPoints();
  InputMeshType::PointDataContainerPointer pointsData = inputMesh->GetPointData();

  inputMesh->SetPointData( NULL );

  TRY_EXPECT_EXCEPTION( filter->Update() );

  inputMesh->SetPointData( pointsData );

  TRY_EXPECT_NO_EXCEPTION( filter->Update() );

  inputMesh->SetPoints( NULL );

  TRY_EXPECT_EXCEPTION( filter->Update() );

  inputMesh->SetPoints( points );

  TRY_EXPECT_NO_EXCEPTION( filter->Update() );
#endif

  return EXIT_SUCCESS;
}
