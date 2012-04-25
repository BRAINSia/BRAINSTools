/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLinearInterpolateMeshFunctionTest3.cxx,v $
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

#include <math.h>
#include "itkVector.h"
#include "itkQuadEdgeMesh.h"
#include "itkRegularSphereMeshSource.h"
#include "itkMeanSquaresMeshToMeshMetric.h"
#include "itkMeshToMeshRegistrationMethod.h"
#include "itkLinearInterpolateMeshFunction.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"
#include "itkMeshGeneratorHelper.h"
#include <iostream>

static float mapSphericalCoordinatesFunction(float inPhi);

static itk::CovariantVector<float, 3> mapSphericalCoordinatesFunctionGradient(float inPhi, float inTheta,
                                                                              bool printFlag);

// Really simple example: a sinusoid mapping between 0 and 1 as a function of theta, constant in phi
static float
mapSphericalCoordinatesFunction(float inTheta)
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

  // Need to multiply dF/dtheta by unit vector theta
  // unit vector= -sinTheta i + cosTheta j

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

int main( int argc, char * argv [] )
{
  if( argc < 2 )
    {
    std::cerr << "Missing Argument " << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " outputDirectory " << std::endl;
    }

  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<float, Dimension> MovingMeshType;
  typedef itk::QuadEdgeMesh<float, Dimension> FixedMeshType;

  FixedMeshType::Pointer  fixedMesh;
  MovingMeshType::Pointer movingMesh;

  typedef MeshGeneratorHelper<FixedMeshType, MovingMeshType> GeneratorType;

  GeneratorType::GenerateMeshes( fixedMesh, movingMesh );

  typedef itk::QuadEdgeMesh<float, 3>            MeshType;
  typedef itk::RegularSphereMeshSource<MeshType> SphereMeshSourceType;

  SphereMeshSourceType::Pointer mySphereMeshSource = SphereMeshSourceType::New();

  typedef SphereMeshSourceType::PointType  PointType;
  typedef SphereMeshSourceType::VectorType VectorType;

  // Set up synthetic data.

  PointType myCenter;
  myCenter.Fill( 0.0 );

  VectorType myScale;
  myScale.Fill( 1.0 );

  const unsigned int numberOfSubdivisions = 4;

  mySphereMeshSource->SetCenter( myCenter );
  mySphereMeshSource->SetResolution( numberOfSubdivisions );
  mySphereMeshSource->SetScale( myScale );
  mySphereMeshSource->Modified();

  try
    {
    mySphereMeshSource->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error during source Update() " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "mySphereMeshSource: " << mySphereMeshSource;

  MeshType::Pointer myMesh = mySphereMeshSource->GetOutput();

  PointType myPt;
  myPt.Fill(0.0f);

  std::cout << "Testing itk::RegularSphereMeshSource " << std::endl;
  for( unsigned int i = 0; i < myMesh->GetNumberOfPoints(); i++ )
    {
    myMesh->GetPoint(i, &myPt);

    const VectorType radial = myPt - myCenter;

    const double radius = radial.GetNorm(); // assuming radius is not valued 1

    const double theta = vcl_atan2( myPt[1], myPt[0] );

    const double phi = vcl_acos( myPt[2] / radius );

    const double myValue = mapSphericalCoordinatesFunction(theta);

    myMesh->SetPointData(i, myValue);

    std::cout << "Point[" << i << "]: " << myPt << " radius "
              << radius << "  theta " << theta << "  phi "
              << phi  << "  myValue " << myValue << std::endl;
    }

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType> WriterType;
  WriterType::Pointer myWriter = WriterType::New();
  myWriter->SetInput( myMesh );

  std::string filename = argv[1];
  filename += "/itkLinearInterpolateMeshFunctionTest3.vtk";

  myWriter->SetFileName( filename );

  try
    {
    myWriter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Error during myWriter Update() " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }

// ------------------------------------------------------------
// Set up an Interpolator
// ------------------------------------------------------------
  typedef itk::LinearInterpolateMeshFunction<
      MeshType> InterpolatorType;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  interpolator->SetInputMesh( myMesh );

  interpolator->Initialize();

// ------------------------------------------------------------
// Simple Test of derivative computed by Interpolator
// ------------------------------------------------------------

// For all cells of mesh, find the middle point, and compare analytical
// value of derivative with value returned by triangle-based estimation.

  typedef MeshType::CellsContainerPointer  CellsContainerPointer;
  typedef MeshType::CellsContainerIterator CellsContainerIterator;
  typedef MeshType::CellType               CellType;
  typedef CellType::PointType              CellPointType;
  typedef CellType::PointIdIterator        PointIdIterator;

  CellsContainerPointer cells = myMesh->GetCells();

  unsigned faceId = 0;

  float maximumDifferenceMagnitude = 0.0;
  for( MeshType::CellsContainerIterator cells_it = cells->Begin();
       cells_it != cells->End();
       ++cells_it, faceId++ )
    {
    CellType* cellPointer = cells_it.Value();

    if( cellPointer->GetType() != 1 )
      {
      std::cout << "Face " << faceId << " has " << cellPointer->GetNumberOfPoints()
                << " points" << std::endl;
      }

    PointIdIterator pointIdIterator = cellPointer->PointIdsBegin();

    if( cellPointer->GetNumberOfPoints() != 3 )  // ignore non triangular cells.
      {
      std::cerr << "Ignoring non-triangular face " << std::endl;
      continue;
      }

    CellPointType points[3];  // three nodes in the triangular cell
    for( unsigned int pointId = 0; pointId < 3; pointId++ )
      {
      points[pointId] = myMesh->GetPoint(*pointIdIterator);
      pointIdIterator++;
      }

    PointType myCellCenter;
    myCellCenter.SetToBarycentricCombination( points[0], points[1], points[2], 0.5, 0.5 );

    VectorType radial = myCellCenter.GetVectorFromOrigin();

    const double cellCenterRadius = radial.GetNorm();
    const double cellCenterTheta  = vcl_atan2( myCellCenter[1], myCellCenter[0] );
    const double cellCenterPhi    = vcl_acos( myCellCenter[2] / cellCenterRadius );

    if( faceId == 1365 )
      {
      std::cout << " problem case \n";
      std::cout << " stop here \n";
      }

    InterpolatorType::DerivativeType computedDerivative;
    interpolator->EvaluateDerivative( myCellCenter, computedDerivative );

    itk::CovariantVector<float, 3> analyticalDerivative = mapSphericalCoordinatesFunctionGradient(cellCenterPhi,
                                                                                                  cellCenterTheta,
                                                                                                  false);

    std::cout << " faceId  " << faceId << "  cell Center " << myCellCenter
              << "  analytical derivative "
              << analyticalDerivative
              << "  computed derivative "
              << computedDerivative
              << std::endl;

    itk::CovariantVector<float, 3> differenceValue = analyticalDerivative - computedDerivative;
    const float                    differenceValueMagnitude = differenceValue.GetNorm();

    if( differenceValueMagnitude >  maximumDifferenceMagnitude )
      {
      maximumDifferenceMagnitude = differenceValueMagnitude;

      std::cout << " maximum difference magnitude increasing " << maximumDifferenceMagnitude << " faceId " << faceId
                << std::endl;
      }
    }

  std::cout << " maximum difference magnitude " << maximumDifferenceMagnitude << std::endl;

  std::cout << "Test End " << std::endl;

  return EXIT_SUCCESS;
}
