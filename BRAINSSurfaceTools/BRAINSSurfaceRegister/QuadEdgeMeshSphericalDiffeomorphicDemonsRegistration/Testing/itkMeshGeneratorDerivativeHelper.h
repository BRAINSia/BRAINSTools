/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNodeScalarGradientCalculatorTest2.cxx,v $
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
static float mapSphericalCoordinatesFunction(float inPhi, float inTheta);

static itk::CovariantVector<float, 3>
mapSphericalCoordinatesFunctionGradient(float inPhi, float inTheta, bool printFlag);

// Really simple example: a sinusoid mapping between 0 and 1 as a function of theta, constant in phi
static float
mapSphericalCoordinatesFunction(float itkNotUsed(inPhi), float inTheta)
{
  float result;

  float thetaFactor = vcl_sin(inTheta); // simplest non-constant function of theta

  // that is smooth at 0-2pi boundary

  result = thetaFactor;

  return result;
}

static itk::CovariantVector<float, 3>
mapSphericalCoordinatesFunctionGradient(float inPhi, float inTheta, bool printFlag)
{
  itk::CovariantVector<float, 3> phiComponent;
  itk::CovariantVector<float, 3> thetaComponent;
  itk::CovariantVector<float, 3> result;

  float cosTheta = vcl_cos(inTheta);
  float sinTheta = vcl_sin(inTheta);
  float cosPhi = vcl_cos(inPhi);
  float sinPhi = vcl_sin(inPhi);

  // derivative of phiFactor over Theta
  float functionDerivativeOverTheta = vcl_cos(inTheta);

  // Need to multiply dF/dtheta by unit vector
  // unit vector= vcl_cos(phi)*vcl_cos(theta) i + vcl_cos(phi)*vcl_sin(theta) j - vcl_sin(phi)k

  thetaComponent[0] = -sinTheta * functionDerivativeOverTheta;
  thetaComponent[1] = cosTheta * functionDerivativeOverTheta;
  thetaComponent[2] = 0.0;

  result = thetaComponent;

  if( printFlag )
    {
    std::cout << "  inTheta " << inTheta << "  inPhi " << inPhi
              << "  sinTheta " << sinTheta << "  cosTheta " << cosTheta
              << "  sinPhi " << sinPhi << "  cosPhi " << cosPhi
              << "  dfdphi " << functionDerivativeOverTheta
              << "  thetaComponent " << thetaComponent
              << "  phiComponent " << phiComponent << "  result " << result << " \n";
    }

  return result;
}

template <class TFixedMesh, class TMovingMesh>
class MeshGeneratorDerivativeHelper
{
public:
  static void GenerateMeshes(
    typename TFixedMesh::Pointer & fixedMesh, typename TMovingMesh::Pointer & movingMesh )
  {
    typedef itk::QuadEdgeMesh<float, 3> MovingMeshType;
    typedef itk::QuadEdgeMesh<float, 3> FixedMeshType;

    typedef itk::RegularSphereMeshSource<MovingMeshType> MovingSphereMeshSourceType;
    typedef itk::RegularSphereMeshSource<FixedMeshType>  FixedSphereMeshSourceType;

    MovingSphereMeshSourceType::Pointer movingSphereMeshSource = MovingSphereMeshSourceType::New();
    FixedSphereMeshSourceType::Pointer  fixedShpereMeshSource = FixedSphereMeshSourceType::New();

    typedef MovingSphereMeshSourceType::PointType  MovingPointType;
    typedef FixedSphereMeshSourceType::PointType   FixedPointType;
    typedef MovingSphereMeshSourceType::VectorType MovingVectorType;
    typedef FixedSphereMeshSourceType::VectorType  FixedVectorType;

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

    movingSphereMeshSource->SetCenter( movingCenter );
    movingSphereMeshSource->SetResolution( 4 );
    movingSphereMeshSource->SetScale( movingScale );
    movingSphereMeshSource->Modified();

    fixedShpereMeshSource->SetCenter( fixedCenter );
    fixedShpereMeshSource->SetResolution( 4 );
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

    typedef MovingPointType::VectorType MovingVectorType;
    typedef FixedPointType::VectorType  FixedVectorType;

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

      const float fixedValue = mapSphericalCoordinatesFunction(phi, theta);

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

      const float movingValue = mapSphericalCoordinatesFunction(phi, movingTheta);

      movingMesh->SetPointData(i, movingValue);
      }

    typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<FixedMeshType> FixedWriterType;
    FixedWriterType::Pointer fixedWriter = FixedWriterType::New();
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
    MovingWriterType::Pointer movingWriter = MovingWriterType::New();
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
