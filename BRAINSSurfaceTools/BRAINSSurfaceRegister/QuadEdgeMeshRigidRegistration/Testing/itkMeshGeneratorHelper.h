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

#include "itkQuadEdgeMesh.h"
#include "itkRegularSphereMeshSource.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

// Linear mapping between 0 and 1 as a function of phi and theta
static float
sineMapSphericalCoordinatesFunction(float inPhi, float inTheta)
{
  float result;

  const float pi = atan(1.0) * 4.0;

  float phiFactor = (pi - inPhi) / pi;            // inPhi should be in [0,PI]; peak at North Pole: phi=0.
  float thetaFactor = (sin(inTheta) + 1.0) / 2.0; // inTheta should be in [-PI,PI];

  result = phiFactor * thetaFactor;

  return result;
}

template <class TFixedMesh, class TMovingMesh>
class MeshGeneratorHelper
{
public:
  static void GenerateMeshes(
    typename TFixedMesh::Pointer & fixedMesh, typename TMovingMesh::Pointer & movingMesh )
  {
    typedef TMovingMesh MovingMeshType;
    typedef TFixedMesh  FixedMeshType;

    typedef itk::RegularSphereMeshSource<MovingMeshType> MovingSphereMeshSourceType;
    typedef itk::RegularSphereMeshSource<FixedMeshType>  FixedSphereMeshSourceType;

    typename MovingSphereMeshSourceType::Pointer  movingSphereMeshSource = MovingSphereMeshSourceType::New();
    typename FixedSphereMeshSourceType::Pointer   fixedShpereMeshSource = FixedSphereMeshSourceType::New();

    typedef typename MovingSphereMeshSourceType::PointType  MovingPointType;
    typedef typename FixedSphereMeshSourceType::PointType   FixedPointType;
    typedef typename MovingSphereMeshSourceType::VectorType MovingVectorType;
    typedef typename FixedSphereMeshSourceType::VectorType  FixedVectorType;

    // Set up synthetic data. Two spherical meshes, one is rotated theta=pi/4
    // from the other

    MovingPointType movingCenter;
    movingCenter.Fill( 0.0 );
    MovingPointType fixedCenter;
    fixedCenter.Fill( 0.0 );

    MovingVectorType movingScale;
    movingScale.Fill( 1.0 );
    FixedVectorType fixedScale;
    fixedScale.Fill( 1.0 );

    const unsigned int numberOfSubdivisions = 4;

    movingSphereMeshSource->SetCenter( movingCenter );
    movingSphereMeshSource->SetResolution( numberOfSubdivisions );
    movingSphereMeshSource->SetScale( movingScale );
    movingSphereMeshSource->Modified();

    fixedShpereMeshSource->SetCenter( fixedCenter );
    fixedShpereMeshSource->SetResolution( numberOfSubdivisions );
    fixedShpereMeshSource->SetScale( fixedScale );
    fixedShpereMeshSource->Modified();

    try
      {
      movingSphereMeshSource->Update();
      fixedShpereMeshSource->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Error during source Update() " << std::endl;
      std::cerr << excp << std::endl;
      return;
      }

    std::cout << "movingSphereMeshSource: " << movingSphereMeshSource;
    std::cout << "fixedShpereMeshSource: " << fixedShpereMeshSource;

    fixedMesh = fixedShpereMeshSource->GetOutput();
    movingMesh = movingSphereMeshSource->GetOutput();

    MovingPointType movingPt;
    FixedPointType  fixedPt;

    movingPt.Fill(0.0f);
    fixedPt.Fill(0.0f);

    MovingVectorType movingVtr;
    FixedVectorType  fixedVtr;

    movingVtr.Fill(0.0f);
    fixedVtr.Fill(0.0f);

    std::cout << "Testing itk::RegularSphereMeshSource " << std::endl;

    fixedMesh->Print( std::cout );
    for( unsigned int i = 0; i < fixedMesh->GetNumberOfPoints(); i++ )
      {
      fixedMesh->GetPoint(i, &fixedPt);

      fixedVtr = fixedPt - fixedCenter;

      const float radius = fixedVtr.GetNorm();
      const float theta  = atan2(fixedVtr[1], fixedVtr[0]);
      const float phi    = acos(fixedVtr[2] / radius);

      const float fixedValue = sineMapSphericalCoordinatesFunction(phi, theta);

      fixedMesh->SetPointData(i, fixedValue);
      }

    const double pi = 4.0 * atan( 1.0 );
    for( unsigned int i = 0; i < movingMesh->GetNumberOfPoints(); i++ )
      {
      movingMesh->GetPoint(i, &movingPt);

      movingVtr = movingPt - movingCenter;

      const float radius = movingVtr.GetNorm();
      const float theta  = atan2(movingVtr[1], movingVtr[0]);
      const float phi    = acos(movingVtr[2] / radius);

      float movingTheta = theta + pi / 4.0;
      if( movingTheta > pi )
        {
        movingTheta -= 2.0 * pi;
        }

      const float movingValue = sineMapSphericalCoordinatesFunction(phi, movingTheta);

      movingMesh->SetPointData(i, movingValue);
      }

    typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<FixedMeshType> FixedWriterType;
    typename FixedWriterType::Pointer fixedWriter = FixedWriterType::New();
    fixedWriter->SetInput( fixedMesh );
    fixedWriter->SetFileName( "FixedMesh.vtk" );

    try
      {
      fixedWriter->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Error during fixedWriter Update() " << std::endl;
      std::cerr << excp << std::endl;
      return;
      }

    typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MovingMeshType> MovingWriterType;
    typename MovingWriterType::Pointer movingWriter = MovingWriterType::New();
    movingWriter->SetInput( movingMesh );
    movingWriter->SetFileName( "MovingMesh.vtk" );

    try
      {
      movingWriter->Update();
      }
    catch( itk::ExceptionObject & excp )
      {
      std::cerr << "Error during movingWriter Update() " << std::endl;
      std::cerr << excp << std::endl;
      return;
      }
  }
};
