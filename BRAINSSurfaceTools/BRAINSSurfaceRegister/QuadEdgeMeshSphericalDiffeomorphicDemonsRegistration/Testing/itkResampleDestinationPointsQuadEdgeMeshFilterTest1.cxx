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
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkMeshGeneratorHelper.h"
#include "itkMeshWriterHelper1.h"
#include "itkTestingMacros.h"
#include "itkFilterWatcher.h"
#include "itkResampleDestinationPointsQuadEdgeMeshFilter.h"

int main( int argc, char *argv[] )
{
  if( argc < 10 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedMeshFile  movingMeshFile outputMeshfile ";
    std::cerr << " referenceMeshFile resampledOutputMeshFile epsilon sigmaX ";
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

  FixedReaderType::Pointer referenceReader = FixedReaderType::New();
  referenceReader->SetFileName( argv[4] );

  try
    {
    fixedReader->Update();
    movingReader->Update();
    referenceReader->Update();
    }
  catch( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while reading the input file " << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  demonsFilter->SetFixedMesh( fixedReader->GetOutput() );
  demonsFilter->SetMovingMesh( movingReader->GetOutput() );

  DemonsFilterType::PointType center0;
  center0.Fill( 0.0 );

  demonsFilter->SetSphereCenter( center0 );
  demonsFilter->SetSphereRadius( 100.0 );

  const double epsilon = atof( argv[6] );
  demonsFilter->SetEpsilon( epsilon );

  const double sigmaX = atof( argv[7] );
  demonsFilter->SetSigmaX( sigmaX );

  const double lambda = atof( argv[8] );
  demonsFilter->SetLambda( lambda );

  unsigned int maximumNumberOfSmoothingIterations = atoi( argv[9] );
  demonsFilter->SetMaximumNumberOfSmoothingIterations( maximumNumberOfSmoothingIterations );

  unsigned int maximumNumberOfIterations = atoi( argv[10] );
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

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<FixedMeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetFileName( argv[3] );
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

  // Supersample the list of destination points using the mesh at the next resolution level.
  typedef FixedMeshType::Traits                               MeshTraits;
  typedef itk::PointSet<MeshPixelType, Dimension, MeshTraits> PointSetType;

  typedef itk::ResampleDestinationPointsQuadEdgeMeshFilter<
      PointSetType, FixedMeshType, FixedMeshType, PointSetType> UpsampleDestinationPointsFilterType;

  UpsampleDestinationPointsFilterType::Pointer upsampleDestinationPoints =
    UpsampleDestinationPointsFilterType::New();

  upsampleDestinationPoints->SetInput( demonsFilter->GetFinalDestinationPoints() );
  upsampleDestinationPoints->SetFixedMesh( fixedReader->GetOutput() );
  upsampleDestinationPoints->SetReferenceMesh( referenceReader->GetOutput() );
  upsampleDestinationPoints->SetTransform( itk::IdentityTransform<double>::New() );

  upsampleDestinationPoints->SetSphereCenter( center0 );
  upsampleDestinationPoints->SetSphereRadius( 100.0 );

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

  FixedMeshType::Pointer referenceMesh = referenceReader->GetOutput();
  referenceMesh->Print(std::cout);

  referenceMesh->DisconnectPipeline();

  PointSetType::ConstPointer upsampledPointSet = upsampleDestinationPoints->GetOutput();

  upsampledPointSet->Print( std::cout );

  const PointSetType::PointsContainer * upsampledPoints = upsampledPointSet->GetPoints();

  PointSetType::PointsContainerConstIterator upsampledPointsItr = upsampledPoints->Begin();
  PointSetType::PointsContainerConstIterator upsampledPointsEnd = upsampledPoints->End();

  FixedMeshType::PointsContainer::Pointer referencePoints = referenceMesh->GetPoints();

  FixedMeshType::PointsContainerIterator referencePointItr = referencePoints->Begin();

  while( upsampledPointsItr != upsampledPointsEnd )
    {
    referencePointItr.Value() = upsampledPointsItr.Value();
    ++referencePointItr;
    ++upsampledPointsItr;
    }

  writer->SetFileName( argv[5] );
  writer->SetInput( referenceMesh );

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
