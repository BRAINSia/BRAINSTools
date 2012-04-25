/*=========================================================================

 Program:   BRAINS (Brain Research: Analysis of Images, Networks, and Systems)
 Module:    $RCSfile: $
 Language:  C++
 Date:      $Date: 2011/07/09 14:53:40 $
 Version:   $Revision: 1.0 $

 Copyright (c) University of Iowa Department of Radiology. All rights reserved.
 See GTRACT-Copyright.txt or http://mri.radiology.uiowa.edu/copyright/GTRACT-Copyright.txt
 for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
#include <vector>
#include "itkQuadEdgeMeshSphericalDiffeomorphicDemonsRegistrationConfigure.h"

#include "itkCommand.h"

#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshTraits.h"
#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

#include "itkMultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter.h"
#include "itkDeformQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshGenerateDeformationFieldFilter.h"
#include "itkResampleDestinationPointsQuadEdgeMeshFilter.h"
#include "MultiResolutionRegistrationCLP.h"

int main( int argc, char * argv [] )
{
  PARSE_ARGS;
  // check the input fixed mesh list
  if( fixedMeshFileList.size() != 4 )
    {
    std::cerr << "fixed meshes in 4 levels should be specified" << std::endl;
    return 1;
    }
  // check the input moving mesh list
  if( movingMeshFileList.size() != 4 )
    {
    std::cerr << "moving meshes in 4 levels should be specified" << std::endl;
    return 1;
    }
  // check the input fixed mesh -- original
  if( fixedMeshFileName == "" )
    {
    std::cerr << "No fixed original mesh file specified" << std::endl;
    return 1;
    }
  // check the input moving mesh -- original
  if( movingMeshFileName == "" )
    {
    std::cerr << "No moving original mesh file specified" << std::endl;
    return 1;
    }
  // check the output deformed fixed -- highest res
  if( deformedFileNameRes4 == "" )
    {
    std::cerr << "No output deformed res4 file specified" << std::endl;
    return 1;
    }
  // check the output deformed fixed -- original
  if( deformedFileName == "" )
    {
    std::cerr << "No output deformed file specified" << std::endl;
    return 1;
    }
  // check the output deformation field -- original
  if( deformationFieldFileName == "" )
    {
    std::cerr << "No output deformation field file specified" << std::endl;
    return 1;
    }
  // ONLY 4 resolution levels are allowed
  if( resolutionLevels != 4 )
    {
    std::cerr << "ONLY 4 resolution levels are allowed" << std::endl;
    return 1;
    }

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Resolution levels: " << resolutionLevels << std::endl;
  std::cout << "Input Fixed Mesh Res1: " << fixedMeshFileList[0] << std::endl;
  std::cout << "Input Moving Mesh Res1: " << movingMeshFileList[0] << std::endl;
  std::cout << "Input Fixed Mesh Res2: " << fixedMeshFileList[1] << std::endl;
  std::cout << "Input Moving Mesh Res2: " << movingMeshFileList[1] << std::endl;
  std::cout << "Input Fixed Mesh Res3: " << fixedMeshFileList[2] << std::endl;
  std::cout << "Input Moving Mesh Res3: " << movingMeshFileList[2] << std::endl;
  std::cout << "Input Fixed Mesh Res4: " << fixedMeshFileList[3] << std::endl;
  std::cout << "Input Moving Mesh Res4: " << movingMeshFileList[3] << std::endl;
  std::cout << "Input Fixed Mesh Orig: " << fixedMeshFileName << std::endl;
  std::cout << "Input Moving Mesh Orig: " << movingMeshFileName << std::endl;
  std::cout << "Output Deformed Fixed Res4: " << deformedFileNameRes4 << std::endl;
  std::cout << "Output Deformed Fixed Orig: " << deformedFileName << std::endl;
  std::cout << "Output Deformation Field Orig: " << deformationFieldFileName << std::endl;
  std::cout << "------------------" << std::endl;
  std::cout << "Registration Parameters: " << std::endl;
  std::cout << "Rigid Iterations Res1: " << rigidIterations[0] << std::endl;
  std::cout << "Rigid Iterations Res2: " << rigidIterations[1] << std::endl;
  std::cout << "Rigid Iterations Res3: " << rigidIterations[2] << std::endl;
  std::cout << "Rigid Iterations Res4: " << rigidIterations[3] << std::endl;
  std::cout << "Demons Iterations Res1: " << demonsIterations[0] << std::endl;
  std::cout << "Demons Iterations Res2: " << demonsIterations[1] << std::endl;
  std::cout << "Demons Iterations Res3: " << demonsIterations[2] << std::endl;
  std::cout << "Demons Iterations Res4: " << demonsIterations[3] << std::endl;
  std::cout << "Smooth Iterations Res1: " << smoothIterations[0] << std::endl;
  std::cout << "Smooth Iterations Res2: " << smoothIterations[1] << std::endl;
  std::cout << "Smooth Iterations Res3: " << smoothIterations[2] << std::endl;
  std::cout << "Smooth Iterations Res4: " << smoothIterations[3] << std::endl;
  std::cout << "Significance Res1: " << metricSignificance[0] << std::endl;
  std::cout << "Significance Res2: " << metricSignificance[1] << std::endl;
  std::cout << "Significance Res3: " << metricSignificance[2] << std::endl;
  std::cout << "Significance Res4: " << metricSignificance[3] << std::endl;
  std::cout << "---------------------------------------------------" << std::endl;

  typedef float MeshPixelType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> MeshType;

  typedef itk::MultiResolutionQuadEdgeMeshSphericalDiffeomorphicDemonsFilter<
      MeshType>  MultiResolutionDemonsFilterType;

  MultiResolutionDemonsFilterType::Pointer multiResDemonsFilter =
    MultiResolutionDemonsFilterType::New();

  typedef MultiResolutionDemonsFilterType::DestinationPointSetType DestinationPointSetType;

  multiResDemonsFilter->SetNumberOfResolutionLevels( resolutionLevels );

  typedef MultiResolutionDemonsFilterType::IntegerArrayType IntegerArrayType;

  IntegerArrayType smoothingIterations(resolutionLevels);

  smoothingIterations[0] = smoothIterations[0];
  smoothingIterations[1] = smoothIterations[1];
  smoothingIterations[2] = smoothIterations[2];
  smoothingIterations[3] = smoothIterations[3];

  multiResDemonsFilter->SetSmoothingIterations( smoothingIterations );

  IntegerArrayType demonsRegIterations(resolutionLevels);

  demonsRegIterations[0] = demonsIterations[0];
  demonsRegIterations[1] = demonsIterations[1];
  demonsRegIterations[2] = demonsIterations[2];
  demonsRegIterations[3] = demonsIterations[3];

  multiResDemonsFilter->SetDemonsIterations( demonsRegIterations );

  IntegerArrayType rigidRegIterations(resolutionLevels);

  rigidRegIterations[0] = rigidIterations[0];
  rigidRegIterations[1] = rigidIterations[1];
  rigidRegIterations[2] = rigidIterations[2];
  rigidRegIterations[3] = rigidIterations[3];

  multiResDemonsFilter->SetRigidRegistrationIterations( rigidRegIterations );

  typedef MultiResolutionDemonsFilterType::DoubleArrayType DoubleArrayType;
  DoubleArrayType rigidStepLength(resolutionLevels);

  rigidStepLength[0] = 1e-2;
  rigidStepLength[1] = rigidStepLength[0] / 2.0;
  rigidStepLength[2] = rigidStepLength[1] / 2.0;
  rigidStepLength[3] = rigidStepLength[2] / 2.0;

  multiResDemonsFilter->SetRigidRegistrationStepLength( rigidStepLength );

  typedef MultiResolutionDemonsFilterType::DoubleArrayType DoubleArrayType;

  const double shortestEdgeLengthAtIC4 = 6.92;
  const double shortestEdgeLengthAtIC5 = 3.46;
  const double shortestEdgeLengthAtIC6 = 1.73;
  const double shortestEdgeLengthAtIC7 = 0.86;

  DoubleArrayType epsilon(resolutionLevels);

  epsilon[0] = 1.0 / ( 250.0 * shortestEdgeLengthAtIC4 );
  epsilon[1] = 1.0 / ( 250.0 * shortestEdgeLengthAtIC5 );
  epsilon[2] = 1.0 / ( 250.0 * shortestEdgeLengthAtIC6 );
  epsilon[3] = 1.0 / ( 250.0 * shortestEdgeLengthAtIC7 );

  multiResDemonsFilter->SetEpsilonValues( epsilon );

  DoubleArrayType sigmaX(resolutionLevels);

  sigmaX[0] = 1.0 / vcl_sqrt(epsilon[0]);
  sigmaX[1] = 1.0 / vcl_sqrt(epsilon[1]);
  sigmaX[2] = 1.0 / vcl_sqrt(epsilon[2]);
  sigmaX[3] = 1.0 / vcl_sqrt(epsilon[3]);

  multiResDemonsFilter->SetSigmaXValues( sigmaX );

  DoubleArrayType significance(resolutionLevels);
  significance[0] = metricSignificance[0];
  significance[1] = metricSignificance[1];
  significance[2] = metricSignificance[2];
  significance[3] = metricSignificance[3];

  multiResDemonsFilter->SetMetricSignificances(significance);

  typedef itk::QuadEdgeMeshVTKPolyDataReader<MeshType> ReaderType;

  std::vector<ReaderType::Pointer> fixedMeshReader(resolutionLevels);
  std::vector<ReaderType::Pointer> movingMeshReader(resolutionLevels);
  for( int i = 0; i < resolutionLevels; i++ )
    {
    fixedMeshReader[i] = ReaderType::New();
    movingMeshReader[i] = ReaderType::New();

    fixedMeshReader[i]->SetFileName( fixedMeshFileList[i].c_str() );
    movingMeshReader[i]->SetFileName( movingMeshFileList[i].c_str() );

    fixedMeshReader[i]->Update();
    movingMeshReader[i]->Update();

    multiResDemonsFilter->SetFixedMesh( i, fixedMeshReader[i]->GetOutput() );
    multiResDemonsFilter->SetMovingMesh( i, movingMeshReader[i]->GetOutput() );
    }

  ReaderType::Pointer fixedMeshReader5 = ReaderType::New();
  fixedMeshReader5->SetFileName( fixedMeshFileName.c_str() );

  ReaderType::Pointer movingMeshReader5 = ReaderType::New();
  movingMeshReader5->SetFileName( movingMeshFileName.c_str() );

  fixedMeshReader5->Update();
  movingMeshReader5->Update();

  MultiResolutionDemonsFilterType::PointType center;
  center.Fill( 0.0 );

  const double radius = 100.0;

  multiResDemonsFilter->SetSphereCenter( center );
  multiResDemonsFilter->SetSphereRadius( radius );

  multiResDemonsFilter->SelfRegulatedModeOn();
  multiResDemonsFilter->SelfStopModeOn();

  multiResDemonsFilter->Update();

  typedef itk::DeformQuadEdgeMeshFilter<
      MeshType, MeshType, DestinationPointSetType>  DeformMeshFilterType;

  DeformMeshFilterType::Pointer deformFilter = DeformMeshFilterType::New();

  deformFilter->SetInput( fixedMeshReader[resolutionLevels - 1]->GetOutput() );
  deformFilter->SetReferenceMesh( fixedMeshReader[resolutionLevels - 1]->GetOutput() );

  deformFilter->SetDestinationPoints( multiResDemonsFilter->GetFinalDestinationPoints() );

  deformFilter->SetSphereRadius( radius );
  deformFilter->SetSphereCenter( center );

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetFileName( deformedFileNameRes4.c_str() );
  writer->SetInput( deformFilter->GetOutput() );

  deformFilter->Update();
  writer->Update();

  deformFilter->SetInput( fixedMeshReader5->GetOutput() );
  writer->SetFileName( deformedFileName.c_str() );

  deformFilter->Update();
  writer->Update();

  // upsample finalDestinationPoints
  typedef itk::ResampleDestinationPointsQuadEdgeMeshFilter<
      DestinationPointSetType, MeshType, MeshType, DestinationPointSetType> UpsampleDestinationPointsFilterType;

  UpsampleDestinationPointsFilterType::Pointer upsampleDestinationPoints =
    UpsampleDestinationPointsFilterType::New();

  upsampleDestinationPoints->SetInput( multiResDemonsFilter->GetFinalDestinationPoints() );
  upsampleDestinationPoints->SetFixedMesh( fixedMeshReader[resolutionLevels - 1]->GetOutput() );
  upsampleDestinationPoints->SetReferenceMesh( fixedMeshReader5->GetOutput() );
  upsampleDestinationPoints->SetTransform( itk::IdentityTransform<double>::New() );

  upsampleDestinationPoints->SetSphereCenter(center);
  upsampleDestinationPoints->SetSphereRadius(radius);
  upsampleDestinationPoints->Update();

  typedef MeshType::PointType   PointType;
  typedef PointType::VectorType VectorType;

  typedef itk::QuadEdgeMeshTraits<VectorType, Dimension, bool, bool> VectorPointSetTraits;

  typedef itk::QuadEdgeMesh<VectorType, Dimension, VectorPointSetTraits> MeshWithVectorsType;

  typedef itk::QuadEdgeMeshGenerateDeformationFieldFilter<
      MeshType, DestinationPointSetType, MeshWithVectorsType>   DeformationFilterType;

  DeformationFilterType::Pointer deformationFilter = DeformationFilterType::New();

  deformationFilter->SetInputMesh( fixedMeshReader5->GetOutput() );
  deformationFilter->SetDestinationPoints( upsampleDestinationPoints->GetOutput() );
  deformationFilter->Update();

  // write out deformation field
  typedef itk::QuadEdgeMeshVectorDataVTKPolyDataWriter<MeshWithVectorsType> VectorMeshWriterType;
  VectorMeshWriterType::Pointer vectorMeshWriter = VectorMeshWriterType::New();

  vectorMeshWriter->SetInput( deformationFilter->GetOutput() );
  vectorMeshWriter->SetFileName( deformationFieldFileName.c_str() );
  vectorMeshWriter->Update();

  return 0;
}
