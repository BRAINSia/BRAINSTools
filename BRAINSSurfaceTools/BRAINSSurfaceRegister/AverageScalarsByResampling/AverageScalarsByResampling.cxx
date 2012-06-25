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

#include "itkQuadEdgeMeshVTKPolyDataReader.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

#include "itkQuadEdgeMesh.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateMeshFunction.h"
#include "itkResampleQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshAddScalarsFilter.h"
#include "itkAssignScalarValuesQuadEdgeMeshFilter.h"

#include "AverageScalarsByResamplingCLP.h"

int main( int argc, char * argv [] )
{
  PARSE_ARGS;
  // check the deformed template files
  if( int(deformedTemplateMeshList.size() ) != numSubs )
    {
    std::cerr << "number of deformed meshes do not agree with number of subjects" << std::endl;
    return 1;
    }
  // check the sphere files with scalar values
  if( int(sphereWithScalarsList.size() ) != numSubs )
    {
    std::cerr << "number of spheres with scalars do not agree with number of subjects" << std::endl;
    return 1;
    }
  // check template surface
  if( templateSurfaceFile == "" )
    {
    std::cerr << "the template surface file should be specified" << std::endl;
    return 1;
    }
  // check the template sphere with scalars
  if( templateSphereFile == "" )
    {
    std::cerr << "the template sphere with scalar values should be specified" << std::endl;
    return 1;
    }

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Add scalar values from: " << templateSphereFile << std::endl;
  std::cout << "Template sphere: " << templateSphereFile << std::endl;
  std::cout << "With scalar values from " << numSubs << " subjects calculated as: " << std::endl;
  for( int i = 0; i < numSubs; i++ )
    {
    std::cout << "Use deformed template sphere: " << deformedTemplateMeshList[i] << std::endl;
    std::cout << "To resample scalar values of: " << sphereWithScalarsList[i] << std::endl;
    }
  std::cout << "Map the average scalars onto: " << templateSurfaceFile << std::endl;
  std::cout << "Output surface: " << templateSurfaceWithAverageScalars << std::endl;
  std::cout << "---------------------------------------------------" << std::endl;

  typedef double MeshPixelType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh<MeshPixelType, Dimension> MeshType;

  typedef itk::QuadEdgeMeshVTKPolyDataReader<MeshType> InputMeshReaderType;

  // read input fixed mesh
  InputMeshReaderType::Pointer fixedMeshReader = InputMeshReaderType::New();
  fixedMeshReader->SetFileName(templateSurfaceFile.c_str() );
  fixedMeshReader->Update();

  // set up deformed fixed mesh reader
  InputMeshReaderType::Pointer deformedFixedMeshReader = InputMeshReaderType::New();

  // set up reference mesh reader
  InputMeshReaderType::Pointer referenceMeshReader = InputMeshReaderType::New();

  // set up Interpolators Filter
  typedef itk::IdentityTransform<double> TransformType;

  TransformType::Pointer transform = TransformType::New();

  typedef itk::LinearInterpolateMeshFunction<MeshType> LinearInterpolatorType;

  LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();

  // set up Reample filter
  typedef itk::ResampleQuadEdgeMeshFilter<MeshType, MeshType> ResamplingFilterType;

  ResamplingFilterType::Pointer resampleFilter = ResamplingFilterType::New();

  resampleFilter->SetTransform( transform );
  resampleFilter->SetInterpolator( interpolator );

  // set up Add Scalar Filter
  typedef itk::QuadEdgeMeshAddScalarsFilter<
      MeshType, MeshType, MeshType> AddScalarsFilterType;

  AddScalarsFilterType::Pointer addFilter = AddScalarsFilterType::New();

  // set up Assign Scalar Filter
  typedef itk::AssignScalarValuesQuadEdgeMeshFilter<
      MeshType,
      MeshType,
      MeshType>    AssignFilterType;

  AssignFilterType::Pointer assignFilter  = AssignFilterType::New();

  // read one deformed fixed mesh and one reference mesh once at a time
  // perform resampling and adding in one loop
  MeshType::Pointer deformedFixedMesh = MeshType::New();
  MeshType::Pointer referenceMesh = MeshType::New();
  MeshType::Pointer resampledMesh = MeshType::New();

  MeshType::Pointer outputMesh = MeshType::New();
  MeshType::Pointer newMesh = MeshType::New();

  InputMeshReaderType::Pointer inputWithScalarsReader = InputMeshReaderType::New();
  inputWithScalarsReader->SetFileName(templateSphereFile.c_str() );
  inputWithScalarsReader->Update();

  outputMesh = inputWithScalarsReader->GetOutput();
  outputMesh->DisconnectPipeline();
  for( int i = 0; i < numSubs; i++ )
    {
    // read deformed fixed mesh
    deformedFixedMeshReader->SetFileName( deformedTemplateMeshList[i] );
    deformedFixedMeshReader->Update();

    deformedFixedMesh = deformedFixedMeshReader->GetOutput();
    deformedFixedMesh->DisconnectPipeline();

    // read reference mesh (sphere with scalar values)
    referenceMeshReader->SetFileName( sphereWithScalarsList[i] );
    referenceMeshReader->Update();

    referenceMesh = referenceMeshReader->GetOutput();
    referenceMesh->DisconnectPipeline();

    // resample reference mesh by deformed fixed mesh
    // resample the input by the reference
    resampleFilter->SetInput(referenceMesh);
    resampleFilter->SetReferenceMesh(deformedFixedMesh);
    resampleFilter->Update();
    std::cout << "stop here!!!!" << std::endl;
    resampledMesh = resampleFilter->GetOutput();
    resampledMesh->DisconnectPipeline();

    // assign the resampled scalars to input fixed mesh
    assignFilter->SetInputMesh(fixedMeshReader->GetOutput() );
    assignFilter->SetSourceMesh(resampledMesh);
    assignFilter->Update();

    newMesh = assignFilter->GetOutput();
    newMesh->DisconnectPipeline();

    // add point.Value of newDF to point.Value of output
    addFilter->SetInput1(newMesh);
    addFilter->SetInput2(outputMesh);
    addFilter->Update();

    outputMesh = addFilter->GetOutput();

    outputMesh->DisconnectPipeline();
    }

  // divide the sum of scalars with numSubs
  // iterate through pointData of output
  typedef MeshType::PointDataContainer ScalarContainer;

  ScalarContainer * scalars = outputMesh->GetPointData();

  typedef ScalarContainer::Iterator scalarIterator;

  scalarIterator scalarItr = scalars->Begin();
  scalarIterator scalarEnd = scalars->End();

  while( scalarItr != scalarEnd )
    {
    if( numSubs != 0 )
      {
      scalarItr.Value() = scalarItr.Value() / double(numSubs + 1);
      }
    ++scalarItr;
    }

  // write out deformation field
  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter<MeshType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetInput( outputMesh );
  writer->SetFileName( templateSurfaceWithAverageScalars.c_str() );
  writer->Update();

  return EXIT_SUCCESS;
}
