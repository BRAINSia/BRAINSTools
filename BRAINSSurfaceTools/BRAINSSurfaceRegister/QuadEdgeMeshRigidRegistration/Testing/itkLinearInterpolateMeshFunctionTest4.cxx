/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFreeSurferBinarySurfaceReaderQuadEdgeMeshTest.cxx,v $
  Language:  C++
  Date:      $Date: 2008-06-16 02:07:17 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "itkQuadEdgeMesh.h"
#include "itkVTKPolyDataReader.h"
#include "itkQuadEdgeMeshVectorDataVTKPolyDataWriter.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkLinearInterpolateMeshFunction.h"
#include "itkVersor.h"

#include <iostream>

int main(int argc, char* argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: itkLinearInterpolateMeshFunctionTest inputFilename radius outputDeformationFieldMesh";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<double, Dimension> MeshType;
  typedef itk::VTKPolyDataReader<MeshType>     ReaderType;

  typedef MeshType::PointsContainer PointsContainer;

  ReaderType::Pointer surfaceReader = ReaderType::New();

  typedef ReaderType::PointType PointType;
  typedef PointType::VectorType VectorType;

  typedef itk::Versor<double> VersorType;

  surfaceReader->SetFileName(argv[1]);

  try
    {
    surfaceReader->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error during Update() " << std::endl;
    std::cerr << excp << std::endl;
    }

  MeshType::Pointer mesh = surfaceReader->GetOutput();

  typedef itk::LinearInterpolateMeshFunction<MeshType> InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  interpolator->SetInputMesh( mesh );

  interpolator->Initialize();

  PointsContainer::Pointer points = mesh->GetPoints();

  const unsigned int                             numberOfVerticesInTriangle = 3;
  InterpolatorType::InstanceIdentifierVectorType pointIds(numberOfVerticesInTriangle);

  PointsContainer::ElementIdentifier pointId = 0;

  PointsContainer::ConstIterator pointItr = points->Begin();
  PointsContainer::ConstIterator pointEnd = points->End();

  std::cout << "Checking the original points " << std::endl;

  while( pointItr != pointEnd )
    {
    const PointType & point = pointItr.Value();

    bool triangleFound = interpolator->FindTriangle( point, pointIds );

    if( !triangleFound )
      {
      std::cerr << "Failed to find triangle for point = " << point << std::endl;
      return EXIT_FAILURE;
      }

    if( pointId != pointIds[0]  && pointId != pointIds[1]  && pointId != pointIds[2] )
      {
      std::cerr << "The triangle found doesn't has the point as vertex" << std::endl;
      std::cerr << "pointId " << pointId << std::endl;
      std::cerr << "pointIds[0] " << pointIds[0] << std::endl;
      std::cerr << "pointIds[1] " << pointIds[1] << std::endl;
      std::cerr << "pointIds[2] " << pointIds[2] << std::endl;
      std::cerr << "point = " << point << std::endl;
      std::cerr << "closest point = " << mesh->GetPoint( pointId ) << std::endl;
      return EXIT_FAILURE;
      }

    ++pointItr;
    ++pointId;
    }

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "SUCCESS: Checking the original points " << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "Checked " << pointId << " points " << std::endl;
  std::cout << "from a mesh with " << mesh->GetNumberOfPoints() << " points " << std::endl;

  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "Repeating test with perturbation " << std::endl;

  // Storing the perturbation field in a mesh.
  typedef itk::QuadEdgeMeshTraits<VectorType, Dimension, bool, bool> VectorPointSetTraits;

  typedef itk::QuadEdgeMesh<VectorType, Dimension, VectorPointSetTraits> MeshWithVectorsType;
  MeshWithVectorsType::Pointer vectorMesh = MeshWithVectorsType::New();

  typedef MeshWithVectorsType::PointDataContainer PointDataContainer;
  PointDataContainer::Pointer perturbationVectors = PointDataContainer::New();

  perturbationVectors->Reserve( mesh->GetNumberOfPoints() );

  vectorMesh->SetPoints( const_cast<MeshType::PointsContainer *>( points.GetPointer() ) );
  vectorMesh->SetPointData( perturbationVectors );

  PointDataContainer::Iterator vectorItr = perturbationVectors->Begin();

  const double factor = 0.01;

  const double sphereRadius = atoi( argv[2] );

  const double perturbationDistance = sphereRadius * factor;

  const double angle = vcl_atan( factor );

  VersorType versor;
  PointType  perturbedPoint;

  PointType center;

  center.Fill(0.0);

  pointItr = points->Begin();
  pointEnd = points->End();

  pointId = 0;

  while( pointItr != pointEnd )
    {
    const PointType & point = pointItr.Value();

    unsigned int minimumComponent = 0;

    if( vnl_math_abs( point[minimumComponent] ) > vnl_math_abs( point[1] ) )
      {
      minimumComponent = 1;
      }

    if( vnl_math_abs( point[minimumComponent] ) > vnl_math_abs( point[2] ) )
      {
      minimumComponent = 2;
      }

    VectorType perturbationVector;
    perturbationVector.Fill(0.0);

    perturbationVector[minimumComponent] = 1.0;

    VectorType radialVector = point - center;

    VectorType axis = CrossProduct( perturbationVector, radialVector );
    axis.Normalize();

    versor.Set( axis, angle );

    perturbedPoint.CastFrom( versor.Transform( point ) );

    const double distance = perturbedPoint.EuclideanDistanceTo( point );

    if( distance > perturbationDistance * 2.0 )
      {
      std::cerr << "Perturbed point is too far " << std::endl;
      std::cerr << "original  point = " << point << std::endl;
      std::cerr << "perturbed point = " << perturbedPoint << std::endl;
      std::cerr << "distance = " << distance << std::endl;
      }

    if( distance < perturbationDistance / 2.0 )
      {
      std::cerr << "Perturbed point is too close " << std::endl;
      std::cerr << "original  point = " << point << std::endl;
      std::cerr << "perturbed point = " << perturbedPoint << std::endl;
      std::cerr << "distance = " << distance << std::endl;
      }

    vectorItr.Value() = perturbedPoint - point;

    bool triangleFound = interpolator->FindTriangle( point, pointIds );

    if( !triangleFound )
      {
      std::cerr << "Failed to find triangle for point = " << perturbedPoint << std::endl;
      return EXIT_FAILURE;
      }

    if( pointId != pointIds[0]  && pointId != pointIds[1]  && pointId != pointIds[2] )
      {
      std::cerr << "The triangle found doesn't has the point as vertex" << std::endl;
      std::cerr << "pointId " << pointId << std::endl;
      std::cerr << "pointIds[0] " << pointIds[0] << std::endl;
      std::cerr << "pointIds[1] " << pointIds[1] << std::endl;
      std::cerr << "pointIds[2] " << pointIds[2] << std::endl;
      std::cerr << "point = " << point << std::endl;
      std::cerr << "perturbedPoint = " << perturbedPoint << std::endl;
      std::cout << std::endl;
      std::cout << "SUCCESS: Checking the perturbed points " << std::endl;
      std::cout << std::endl;
      return EXIT_FAILURE;
      }

    ++pointItr;
    ++pointId;
    ++vectorItr;
    }

  typedef itk::QuadEdgeMeshVectorDataVTKPolyDataWriter<MeshWithVectorsType> VectorMeshWriterType;
  VectorMeshWriterType::Pointer vectorMeshWriter = VectorMeshWriterType::New();
  vectorMeshWriter->SetInput( vectorMesh );
  vectorMeshWriter->SetFileName( argv[3] );
  vectorMeshWriter->Update();

  std::cout << "SUCCESS: Checking the perturbed points " << std::endl;
  std::cout << std::endl;

  std::cout << "Test passed" << std::endl;

  return EXIT_SUCCESS;
}
